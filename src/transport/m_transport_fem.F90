! FEM discrete-ordinates (SN) transport solver.
! Implements a DG sweep over a Lagrange-FEM mesh with precomputed LU factors.
! Angular flux lives on a per-element DG DOF layout: index = (ee-1)*n_basis + local.
!
! t_fem_transport extends t_solver and owns all solver state as components,
! eliminating module-level save variables.
!
! Public:
!   t_fem_transport                      -- concrete solver type
!   SolveTransport_FEM                   -- high-level entry
!   Transport_Sweep_FEM                  -- angle sweep over all elements and groups
!   Source_DGFEM_FEM                     -- build total source (scatter + fission)
!   Calculate_Total_Production_DGFEM_FEM -- scalar production integral for k-eff
module m_transport_fem
    use m_constants
    use m_types
    use m_quadrature
    use m_material
    use m_power_iteration, only: t_solver, PowerIteration
    use m_output_fem,      only: export_transport_vtk_fem
    implicit none
    public :: t_fem_transport, SolveTransport_FEM
    public :: Transport_Sweep_FEM, Source_DGFEM_FEM, Calculate_Total_Production_DGFEM_FEM

    ! ------------------------------------------------------------------
    ! Concrete FEM-DG transport solver.
    ! ------------------------------------------------------------------
    type, extends(t_solver) :: t_fem_transport
        type(t_mesh_fem),      pointer :: mesh    => null()
        type(t_basis_fem),     pointer :: FE      => null()
        type(t_sn_quadrature), pointer :: sn_quad => null()
        type(t_fem_dg),        pointer :: TD      => null()
        type(t_material),      pointer :: mats(:) => null()
        real(dp), allocatable :: ang_flux(:,:,:)
        integer,  allocatable :: sweep_order(:,:)
        integer,  allocatable :: ref_ids(:)
        integer :: n_groups = 0
    contains
        procedure :: build_source => fem_build_source
        procedure :: do_solve     => fem_do_solve
        procedure :: compute_prod => fem_compute_prod
        procedure :: snapshot     => fem_snapshot_impl
    end type t_fem_transport

    interface
        subroutine dgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            character, intent(in)  :: trans
            integer, intent(in)    :: n, nrhs, lda, ldb
            real(dp), intent(in)   :: a(lda, *)
            integer, intent(in)    :: ipiv(*)
            real(dp), intent(inout):: b(ldb, *)
            integer, intent(out)   :: info
        end subroutine dgetrs
    end interface

contains

    subroutine Transport_Sweep_FEM(mesh, FE, sn_quad, TD, ang_flux, scalar_flux, &
                                    total_source, sweep_order, ref_ID)
        type(t_mesh_fem),      intent(in)    :: mesh
        type(t_basis_fem),     intent(in)    :: FE
        type(t_sn_quadrature), intent(in)    :: sn_quad
        type(t_fem_dg),        intent(in)    :: TD
        real(dp),              intent(inout) :: ang_flux(:,:,:), scalar_flux(:,:)
        real(dp),              intent(in)    :: total_source(:,:)
        integer,               intent(in)    :: sweep_order(:,:), ref_ID(:)

        integer  :: mm, ee, ie, g, f, info, k
        integer  :: neighbor_elem_id, i_face_node, upwind_dof_idx, local_dof_idx, dof_idx, m_ref
        integer  :: ref_ie, lu_ie
        real(dp) :: b(FE%n_basis), o_n, dir(3), w_mm
        real(dp) :: face_term(FE%n_basis, FE%n_basis)
        integer  :: idx_start, idx_end, f_map(FE%n_nodes_per_face)
        logical  :: is_inflow(mesh%n_faces_per_elem)

        scalar_flux = 0.0_dp
        !$OMP PARALLEL DO PRIVATE(mm, dir, w_mm, ee, ie, ref_ie, lu_ie, idx_start, idx_end, &
        !$OMP&   f, o_n, is_inflow, g, b, f_map, face_term, neighbor_elem_id, i_face_node, &
        !$OMP&   upwind_dof_idx, m_ref, local_dof_idx, info, k, dof_idx)
        do mm = 1, sn_quad%n_angles
            dir  = sn_quad%dirs(mm, :)
            w_mm = sn_quad%weights(mm)

            do ee = 1, mesh%n_elems
                ie     = sweep_order(ee, mm)
                ref_ie = TD%elem_ref_id(ie)
                lu_ie  = TD%elem_lu_id(ie)
                idx_start = (ie-1)*FE%n_basis + 1
                idx_end   =  ie   *FE%n_basis

                do f = 1, mesh%n_faces_per_elem
                    o_n = dot_product(dir, TD%face_normals(:,f,ie))
                    is_inflow(f) = (o_n < 0.0_dp)
                end do

                do g = 1, mesh%n_groups
                    b = total_source(idx_start:idx_end, g)

                    do f = 1, mesh%n_faces_per_elem
                        if (.not. is_inflow(f)) cycle
                        f_map = FE%face_node_map(:, f)
                        face_term = dir(1)*TD%face_mass_x(:,:,f,ref_ie) + &
                                    dir(2)*TD%face_mass_y(:,:,f,ref_ie) + &
                                    dir(3)*TD%face_mass_z(:,:,f,ref_ie)
                        neighbor_elem_id = TD%face_connectivity(1, f, ie)

                        if (neighbor_elem_id > 0) then
                            do i_face_node = 1, FE%n_nodes_per_face
                                upwind_dof_idx = TD%upwind_idx(i_face_node, f, ie)
                                b = b - ang_flux(upwind_dof_idx, mm, g) * face_term(:, f_map(i_face_node))
                            end do
                        else if (TD%face_connectivity(4,f,ie) > 0 .and. &
                                 any(TD%face_connectivity(4,f,ie) == ref_ID)) then
                            m_ref = TD%reflect_map(mm, f, ie)
                            do i_face_node = 1, FE%n_nodes_per_face
                                local_dof_idx = idx_start - 1 + f_map(i_face_node)
                                b = b - ang_flux(local_dof_idx, m_ref, g) * face_term(:, f_map(i_face_node))
                            end do
                        end if
                    end do

                    call dgetrs('N', FE%n_basis, 1, &
                                TD%local_lu(:,:,lu_ie,mm,g), FE%n_basis, &
                                TD%local_pivots(:,lu_ie,mm,g), b, FE%n_basis, info)

                    where (b < 0.0_dp) b = 0.0_dp

                    if (any(b /= b)) then
                        write(*,'(A,2I6)') "FATAL: NaN in FEM sweep, elem/angle=", ie, mm
                        stop
                    end if

                    ang_flux(idx_start:idx_end, mm, g) = b
                    do k = 1, FE%n_basis
                        dof_idx = idx_start + k - 1
                        !$OMP ATOMIC
                        scalar_flux(dof_idx, g) = scalar_flux(dof_idx, g) + w_mm * b(k)
                    end do
                end do
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine Transport_Sweep_FEM

    subroutine Source_DGFEM_FEM(total_src, scalar_flux, k_eff, materials, mesh, FE, TD, &
                                 n_groups, is_adjoint, is_eigenvalue)
        real(dp),              intent(inout) :: total_src(:,:)
        real(dp),              intent(in)    :: scalar_flux(:,:)
        real(dp),              intent(in)    :: k_eff
        type(t_material),      intent(in)    :: materials(:)
        type(t_mesh_fem),      intent(in)    :: mesh
        type(t_basis_fem),     intent(in)    :: FE
        type(t_fem_dg),        intent(in)    :: TD
        integer,               intent(in)    :: n_groups
        logical,               intent(in)    :: is_adjoint, is_eigenvalue

        integer  :: g_to, ee, mat_id, idx_start, idx_end, ref_ee
        real(dp) :: M_phi(FE%n_basis, n_groups), fission_rate(FE%n_basis)

        total_src = 0.0_dp
        !$OMP PARALLEL DO PRIVATE(ee, mat_id, ref_ee, idx_start, idx_end, M_phi, fission_rate, g_to)
        do ee = 1, mesh%n_elems
            mat_id    = mesh%material_ids(ee)
            ref_ee    = TD%elem_ref_id(ee)
            idx_start = (ee-1)*FE%n_basis + 1
            idx_end   =  ee   *FE%n_basis

            M_phi = matmul(TD%elem_mass_matrix(:,:,ref_ee), scalar_flux(idx_start:idx_end,:))

            if (is_adjoint) then
                total_src(idx_start:idx_end,:) = matmul(M_phi, transpose(materials(mat_id)%SigmaS))
            else
                total_src(idx_start:idx_end,:) = matmul(M_phi, materials(mat_id)%SigmaS)
            end if

            if (is_eigenvalue) then
                fission_rate = matmul(M_phi, &
                    merge(materials(mat_id)%NuSigF, materials(mat_id)%Chi, .not. is_adjoint)) / k_eff
                do g_to = 1, n_groups
                    total_src(idx_start:idx_end,g_to) = total_src(idx_start:idx_end,g_to) + &
                        fission_rate * merge(materials(mat_id)%Chi(g_to), materials(mat_id)%NuSigF(g_to), .not. is_adjoint)
                end do
            else
                do g_to = 1, n_groups
                    total_src(idx_start:idx_end,g_to) = total_src(idx_start:idx_end,g_to) + &
                        materials(mat_id)%Src(g_to) * TD%basis_integrals_vol(:,ref_ee)
                end do
            end if

            where (total_src(idx_start:idx_end,:) < 0.0_dp) total_src(idx_start:idx_end,:) = 0.0_dp

            if (any(total_src(idx_start:idx_end,:) /= total_src(idx_start:idx_end,:))) then
                write(*,'(A,I6)') "FATAL: NaN in FEM source, elem=", ee; stop
            end if
        end do
        !$OMP END PARALLEL DO
    end subroutine Source_DGFEM_FEM

    subroutine Calculate_Total_Production_DGFEM_FEM(total_prod, scalar_flux, materials, mesh, FE, TD, is_adjoint)
        real(dp),              intent(out) :: total_prod
        real(dp),              intent(in)  :: scalar_flux(:,:)
        type(t_material),      intent(in)  :: materials(:)
        type(t_mesh_fem),      intent(in)  :: mesh
        type(t_basis_fem),     intent(in)  :: FE
        type(t_fem_dg),        intent(in)  :: TD
        logical,               intent(in)  :: is_adjoint

        integer :: g, ee, mat_id, ref_ee, idx_start, idx_end

        total_prod = 0.0_dp
        !$OMP PARALLEL DO PRIVATE(ee, mat_id, ref_ee, idx_start, idx_end, g) REDUCTION(+:total_prod)
        do ee = 1, mesh%n_elems
            mat_id    = mesh%material_ids(ee)
            ref_ee    = TD%elem_ref_id(ee)
            idx_start = (ee-1)*FE%n_basis + 1
            idx_end   =  ee   *FE%n_basis
            do g = 1, mesh%n_groups
                total_prod = total_prod + &
                    merge(materials(mat_id)%NuSigF(g), materials(mat_id)%Chi(g), .not. is_adjoint) * &
                    dot_product(max(scalar_flux(idx_start:idx_end,g), 0.0_dp), TD%basis_integrals_vol(:,ref_ee))
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine Calculate_Total_Production_DGFEM_FEM

    ! ================================================================== !
    !  Type-bound procedure implementations for t_fem_transport            !
    ! ================================================================== !

    subroutine fem_build_source(self, src, k_eff, is_eigenvalue, is_adjoint)
        class(t_fem_transport), intent(inout) :: self
        real(dp), intent(inout) :: src(:,:)
        real(dp), intent(in)    :: k_eff
        logical,  intent(in)    :: is_eigenvalue, is_adjoint
        call Source_DGFEM_FEM(src, self%scalar_flux, k_eff, self%mats, self%mesh, self%FE, self%TD, &
                               self%n_groups, is_adjoint, is_eigenvalue)
    end subroutine fem_build_source

    subroutine fem_do_solve(self, src)
        class(t_fem_transport), intent(inout) :: self
        real(dp), intent(in) :: src(:,:)
        call Transport_Sweep_FEM(self%mesh, self%FE, self%sn_quad, self%TD, &
                                 self%ang_flux, self%scalar_flux, src, self%sweep_order, self%ref_ids)
    end subroutine fem_do_solve

    subroutine fem_compute_prod(self, prod, is_adjoint)
        class(t_fem_transport), intent(inout) :: self
        real(dp), intent(out) :: prod
        logical,  intent(in)  :: is_adjoint
        call Calculate_Total_Production_DGFEM_FEM(prod, self%scalar_flux, self%mats, self%mesh, &
                                                   self%FE, self%TD, is_adjoint)
    end subroutine fem_compute_prod

    subroutine fem_snapshot_impl(self, k_eff, iter)
        class(t_fem_transport), intent(inout) :: self
        real(dp), intent(in) :: k_eff
        integer,  intent(in) :: iter
        character(len=32)  :: lbl
        character(len=512) :: stag
        self%snap_count = self%snap_count + 1
        write(lbl,'(A,I4.4)') "snap", self%snap_count
        stag = trim(self%snap_tag) // "_" // trim(lbl)
        call export_transport_vtk_fem(trim(self%snap_dir), trim(stag), &
            self%mesh, self%FE, self%sn_quad, self%scalar_flux, self%n_groups, self%vtk_refine)
    end subroutine fem_snapshot_impl

    ! ==================================================================
    ! High-level entry point.
    ! Builds a t_fem_transport object, runs PowerIteration, then moves
    ! results to caller's arrays.
    ! ==================================================================
    subroutine SolveTransport_FEM(solver, mesh, materials, FE, sn_quad, TD, k_eff, &
                                   sweep_order, ref_ids, max_outer, tol, &
                                   is_adjoint, is_eigenvalue, &
                                   snap_dir, snap_tag, vtk_ref)
        type(t_fem_transport), intent(out)           :: solver
        type(t_mesh_fem),      intent(in),    target :: mesh
        type(t_material),      intent(in),    target :: materials(:)
        type(t_basis_fem),     intent(in),    target :: FE
        type(t_sn_quadrature), intent(in),    target :: sn_quad
        type(t_fem_dg),        intent(in),    target :: TD
        real(dp),              intent(out)           :: k_eff
        integer,               intent(in)            :: sweep_order(:,:)
        integer,               intent(in)            :: ref_ids(:)
        integer,               intent(in)            :: max_outer
        real(dp),              intent(in)            :: tol
        logical,               intent(in)            :: is_adjoint, is_eigenvalue
        character(len=*), optional, intent(in)       :: snap_dir, snap_tag
        integer,          optional, intent(in)       :: vtk_ref

        integer :: n_dof

        n_dof = mesh%n_elems * FE%n_basis

        solver%mesh    => mesh
        solver%FE      => FE
        solver%sn_quad => sn_quad
        solver%TD      => TD
        solver%mats    => materials
        solver%n_groups = mesh%n_groups

        solver%have_snap = present(snap_dir) .and. present(snap_tag)
        if (solver%have_snap) then
            solver%snap_dir   = trim(snap_dir)
            solver%snap_tag   = trim(snap_tag)
            solver%vtk_refine = merge(vtk_ref, 4, present(vtk_ref))
            solver%snap_count = 0
        end if

        allocate(solver%ang_flux(n_dof, sn_quad%n_angles, mesh%n_groups), source=0.0_dp)
        allocate(solver%sweep_order, source=sweep_order)
        allocate(solver%ref_ids,     source=ref_ids)
        allocate(solver%scalar_flux(n_dof, mesh%n_groups), source=1.0_dp)

        call PowerIteration(solver, k_eff, max_outer, tol, is_eigenvalue, is_adjoint)
    end subroutine SolveTransport_FEM

end module m_transport_fem
