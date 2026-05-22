! IGA discrete-ordinates (SN) transport solver.
! Implements a DG sweep over an IGA mesh with precomputed LU factors.
! Angular flux lives on a per-element DG DOF layout: index = (ee-1)*n_basis + local.
!
! Public:
!   SolveTransport                   -- high-level entry: initialise state + power iteration
!   Transport_Sweep                  -- angle sweep over all elements and groups
!   Source_DGFEM                     -- build total source (scatter + fission)
!   Calculate_Total_Production_DGFEM -- scalar production integral for k-eff
module m_transport_iga
    use m_constants
    use m_types
    use m_types_iga
    use m_quadrature
    use m_material
    use m_power_iteration, only: PowerIteration
    implicit none
    public :: SolveTransport
    public :: Transport_Sweep, Source_DGFEM, Calculate_Total_Production_DGFEM

    ! ------------------------------------------------------------------
    ! Module-level saved state — set by SolveTransport, read by the
    ! private callback wrappers below.
    ! ------------------------------------------------------------------
    type(t_mesh_iga),      pointer, save :: s_mesh    => null()
    type(t_basis_iga),    pointer, save :: s_FE      => null()
    type(t_sn_quadrature), pointer, save :: s_sn_quad => null()
    type(t_fem_dg), pointer, save :: s_TD      => null()
    type(t_material),      pointer, save :: s_mats(:) => null()
    real(dp), allocatable, save          :: s_ang_flux(:,:,:)
    integer,  allocatable, save          :: s_sweep_order(:,:)
    integer,  allocatable, save          :: s_ref_ids(:)
    integer,  save                       :: s_n_groups = 0

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

    subroutine Transport_Sweep(mesh, FE, sn_quad, TD, ang_flux, scalar_flux, &
                                total_source, sweep_order, ref_ID)
        type(t_mesh_iga),      intent(in)    :: mesh
        type(t_basis_iga),    intent(in)    :: FE
        type(t_sn_quadrature), intent(in)    :: sn_quad
        type(t_fem_dg), intent(in)    :: TD
        real(dp),              intent(inout) :: ang_flux(:,:,:), scalar_flux(:,:)
        real(dp),              intent(in)    :: total_source(:,:)
        integer,               intent(in)    :: sweep_order(:,:), ref_ID(:)

        integer  :: mm, ee, ie, g, f, info, k
        integer  :: neighbor_elem_id, i_face_node, upwind_dof_idx, local_dof_idx, dof_idx, m_ref
        real(dp) :: b(FE%n_basis), o_n, dir(3), w_mm
        real(dp) :: face_term(FE%n_basis, FE%n_basis)
        integer  :: idx_start, idx_end, f_map(FE%n_nodes_per_face)
        logical  :: is_inflow(mesh%n_faces_per_elem)

        scalar_flux = 0.0_dp

        !$OMP PARALLEL DO PRIVATE(mm, dir, w_mm, ee, ie, idx_start, idx_end, f, o_n, is_inflow, &
        !$OMP&   g, b, f_map, face_term, neighbor_elem_id, i_face_node, upwind_dof_idx, &
        !$OMP&   m_ref, local_dof_idx, info, k, dof_idx)
        do mm = 1, sn_quad%n_angles
            dir  = sn_quad%dirs(mm, :)
            w_mm = sn_quad%weights(mm)

            do ee = 1, mesh%n_elems
                ie        = sweep_order(ee, mm)
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
                        face_term = dir(1)*TD%face_mass_x(:,:,f,ie) + &
                                    dir(2)*TD%face_mass_y(:,:,f,ie) + &
                                    dir(3)*TD%face_mass_z(:,:,f,ie)
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
                                TD%local_lu(:,:,ie,mm,g), FE%n_basis, &
                                TD%local_pivots(:,ie,mm,g), b, FE%n_basis, info)

                    where (b < 0.0_dp) b = 0.0_dp

                    if (any(b /= b)) then
                        write(*,'(A,2I6)') "FATAL: NaN in sweep, elem/angle=", ie, mm
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
    end subroutine Transport_Sweep

    subroutine Source_DGFEM(total_src, scalar_flux, k_eff, materials, mesh, FE, TD, &
                             n_groups, is_adjoint, is_eigenvalue)
        real(dp),              intent(inout) :: total_src(:,:)
        real(dp),              intent(in)    :: scalar_flux(:,:)
        real(dp),              intent(in)    :: k_eff
        type(t_material),      intent(in)    :: materials(:)
        type(t_mesh_iga),      intent(in)    :: mesh
        type(t_basis_iga),    intent(in)    :: FE
        type(t_fem_dg), intent(in)    :: TD
        integer,               intent(in)    :: n_groups
        logical,               intent(in)    :: is_adjoint, is_eigenvalue

        integer  :: g_to, ee, mat_id, idx_start, idx_end
        real(dp) :: M_phi(FE%n_basis, n_groups), fission_rate(FE%n_basis)

        total_src = 0.0_dp
        !$OMP PARALLEL DO PRIVATE(ee, mat_id, idx_start, idx_end, M_phi, fission_rate, g_to)
        do ee = 1, mesh%n_elems
            mat_id    = mesh%material_ids(ee)
            idx_start = (ee-1)*FE%n_basis + 1
            idx_end   =  ee   *FE%n_basis

            M_phi = matmul(TD%elem_mass_matrix(:,:,ee), scalar_flux(idx_start:idx_end,:))

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
                        materials(mat_id)%Src(g_to) * TD%basis_integrals_vol(:,ee)
                end do
            end if

            where (total_src(idx_start:idx_end,:) < 0.0_dp) total_src(idx_start:idx_end,:) = 0.0_dp

            if (any(total_src(idx_start:idx_end,:) /= total_src(idx_start:idx_end,:))) then
                write(*,'(A,I6)') "FATAL: NaN in source, elem=", ee; stop
            end if
        end do
        !$OMP END PARALLEL DO
    end subroutine Source_DGFEM

    subroutine Calculate_Total_Production_DGFEM(total_prod, scalar_flux, materials, mesh, FE, TD, is_adjoint)
        real(dp),              intent(out) :: total_prod
        real(dp),              intent(in)  :: scalar_flux(:,:)
        type(t_material),      intent(in)  :: materials(:)
        type(t_mesh_iga),      intent(in)  :: mesh
        type(t_basis_iga),    intent(in)  :: FE
        type(t_fem_dg), intent(in)  :: TD
        logical,               intent(in)  :: is_adjoint

        integer :: g, ee, mat_id, idx_start, idx_end

        total_prod = 0.0_dp
        !$OMP PARALLEL DO PRIVATE(ee, mat_id, idx_start, idx_end, g) REDUCTION(+:total_prod)
        do ee = 1, mesh%n_elems
            mat_id    = mesh%material_ids(ee)
            idx_start = (ee-1)*FE%n_basis + 1
            idx_end   =  ee   *FE%n_basis
            do g = 1, mesh%n_groups
                total_prod = total_prod + &
                    merge(materials(mat_id)%NuSigF(g), materials(mat_id)%Chi(g), .not. is_adjoint) * &
                    dot_product(max(scalar_flux(idx_start:idx_end,g), 0.0_dp), TD%basis_integrals_vol(:,ee))
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine Calculate_Total_Production_DGFEM

    ! ==================================================================
    ! High-level entry point.
    ! Sets up module context, allocates working arrays, then delegates
    ! to the shared PowerIteration driver via three callback wrappers.
    ! ==================================================================
    subroutine SolveTransport(mesh, materials, FE, sn_quad, TD, &
                               scalar_flux, ang_flux_out, k_eff, &
                               sweep_order, ref_ids, max_outer, tol, &
                               is_adjoint, is_eigenvalue)
        type(t_mesh_iga),      intent(in),    target :: mesh
        type(t_material),      intent(in),    target :: materials(:)
        type(t_basis_iga),    intent(in),    target :: FE
        type(t_sn_quadrature), intent(in),    target :: sn_quad
        type(t_fem_dg), intent(in),    target :: TD
        real(dp), allocatable, intent(out)           :: scalar_flux(:,:)
        real(dp), allocatable, intent(out)           :: ang_flux_out(:,:,:)
        real(dp),              intent(out)           :: k_eff
        integer,               intent(in)            :: sweep_order(:,:)
        integer,               intent(in)            :: ref_ids(:)
        integer,               intent(in)            :: max_outer
        real(dp),              intent(in)            :: tol
        logical,               intent(in)            :: is_adjoint, is_eigenvalue

        integer :: n_dof

        n_dof      = mesh%n_elems * FE%n_basis
        s_n_groups = mesh%n_groups

        s_mesh    => mesh
        s_FE      => FE
        s_sn_quad => sn_quad
        s_TD      => TD
        s_mats    => materials

        if (allocated(s_ang_flux))    deallocate(s_ang_flux)
        if (allocated(s_sweep_order)) deallocate(s_sweep_order)
        if (allocated(s_ref_ids))     deallocate(s_ref_ids)

        allocate(s_ang_flux(n_dof, sn_quad%n_angles, s_n_groups), source=0.0_dp)
        allocate(s_sweep_order, source=sweep_order)
        allocate(s_ref_ids,     source=ref_ids)
        allocate(scalar_flux(n_dof, s_n_groups), source=1.0_dp)

        call PowerIteration(scalar_flux, k_eff, max_outer, tol, &
                             is_eigenvalue, is_adjoint,           &
                             transport_source, transport_solve, transport_production)

        call move_alloc(s_ang_flux, ang_flux_out)
    end subroutine SolveTransport

    ! ------------------------------------------------------------------
    ! Callback: build scatter + fission source for all groups.
    ! ------------------------------------------------------------------
    subroutine transport_source(src, flux, k_eff, is_eigenvalue, is_adjoint)
        real(dp), intent(inout) :: src(:,:)
        real(dp), intent(in)    :: flux(:,:), k_eff
        logical,  intent(in)    :: is_eigenvalue, is_adjoint
        call Source_DGFEM(src, flux, k_eff, s_mats, s_mesh, s_FE, s_TD, &
                          s_n_groups, is_adjoint, is_eigenvalue)
    end subroutine transport_source

    ! ------------------------------------------------------------------
    ! Callback: perform one transport sweep (updates flux from src).
    ! ------------------------------------------------------------------
    subroutine transport_solve(flux, src)
        real(dp), intent(inout) :: flux(:,:)
        real(dp), intent(in)    :: src(:,:)
        call Transport_Sweep(s_mesh, s_FE, s_sn_quad, s_TD, &
                             s_ang_flux, flux, src,          &
                             s_sweep_order, s_ref_ids)
    end subroutine transport_solve

    ! ------------------------------------------------------------------
    ! Callback: compute total fission production for k-eff update.
    ! ------------------------------------------------------------------
    subroutine transport_production(prod, flux, is_adjoint)
        real(dp), intent(out) :: prod
        real(dp), intent(in)  :: flux(:,:)
        logical,  intent(in)  :: is_adjoint
        call Calculate_Total_Production_DGFEM(prod, flux, s_mats, s_mesh, s_FE, s_TD, is_adjoint)
    end subroutine transport_production

end module m_transport_iga
