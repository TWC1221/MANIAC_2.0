! Patchwise DG-IGA discrete-ordinates transport solver.
! Each NURBS patch is one DG "element"; the angular flux DOF layout is
!   global_dof = (pp-1)*n_basis_patch + local_patch_dof   (pp = 1..n_patches)
!
! The sweep, source, and production routines are structurally identical to
! m_transport_iga but index over patches instead of knot spans.
!
! Public:
!   SolveTransport_IGA              -- high-level entry
!   Transport_Sweep_Patch                -- angle sweep over patches
!   Source_Patch_DGFEM                   -- scatter + fission source
!   Calculate_Production_Patch_DGFEM    -- fission production for k-eff
!   Remap_Patch_To_Elem_Flux            -- utility: map patch flux → elem layout
module m_transport_iga_pdg
    use m_constants
    use m_types
    use m_quadrature
    use m_material
    use m_power_iteration, only: PowerIteration
    implicit none
    public :: SolveTransport_IGA
    public :: Transport_Sweep_Patch, Source_Patch_DGFEM
    public :: Calculate_Production_Patch_DGFEM, Remap_Patch_To_Elem_Flux

    ! Module-level context set by SolveTransport_IGA, read by callbacks.
    type(t_mesh_iga),      pointer, save :: s_mesh    => null()
    type(t_sn_quadrature), pointer, save :: s_sn_quad => null()
    type(t_patch_dg), pointer, save :: s_PD   => null()
    type(t_material),      pointer, save :: s_mats(:) => null()
    real(dp), allocatable, save          :: s_ang_flux(:,:,:)
    integer,  allocatable, save          :: s_sweep_order(:,:)
    integer,  allocatable, save          :: s_ref_ids(:)
    integer,  save                       :: s_n_groups = 0
    integer,  save                       :: s_n_patches = 0

    interface
        subroutine dgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            character, intent(in)  :: trans
            integer,   intent(in)  :: n, nrhs, lda, ldb
            real(dp),  intent(in)  :: a(lda, *)
            integer,   intent(in)  :: ipiv(*)
            real(dp),  intent(inout) :: b(ldb, *)
            integer,   intent(out) :: info
        end subroutine dgetrs
    end interface

contains

    ! ------------------------------------------------------------------
    ! Angle sweep over all patches.
    ! ------------------------------------------------------------------
    subroutine Transport_Sweep_Patch(mesh, sn_quad, PD, ang_flux, scalar_flux, &
                                      total_source, sweep_order, ref_ID)
        type(t_mesh_iga),      intent(in)    :: mesh
        type(t_sn_quadrature), intent(in)    :: sn_quad
        type(t_patch_dg), intent(in) :: PD
        real(dp),              intent(inout) :: ang_flux(:,:,:), scalar_flux(:,:)
        real(dp),              intent(in)    :: total_source(:,:)
        integer,               intent(in)    :: sweep_order(:,:), ref_ID(:)

        integer  :: mm, ee, pp, g, f, info, k
        integer  :: j_f, local_dof, m_ref, n_face
        real(dp) :: w_mm, dir(3), o_n
        integer  :: idx_start, idx_end, nb
        real(dp) :: b(PD%n_basis_patch)
        real(dp) :: upwind_flux(PD%n_face_basis_max)
        real(dp) :: fm(PD%n_basis_patch, PD%n_basis_patch)
        real(dp), allocatable :: ang_flux_snap(:,:,:)

        nb = PD%n_basis_patch
        scalar_flux = 0.0_dp

        ! Snapshot ang_flux before the parallel sweep so reflective BCs can
        ! safely read paired-angle fluxes without racing concurrent writes.
        ! This is the standard lagged-reflective approach; source iteration
        ! converges the lag to zero at each outer level.
        if (size(ref_ID) > 0) allocate(ang_flux_snap, source=ang_flux)

        !$OMP PARALLEL DO &
        !$OMP& PRIVATE(mm, w_mm, dir, o_n, ee, pp, idx_start, idx_end, f, &
        !$OMP&          g, b, n_face, j_f, local_dof, m_ref, info, k, upwind_flux, fm)
        do mm = 1, sn_quad%n_angles
            w_mm = sn_quad%weights(mm)
            if (PD%matrix_free) dir = sn_quad%dirs(mm, :)

            do ee = 1, size(sweep_order, 1)
                pp        = sweep_order(ee, mm)
                idx_start = (pp - 1)*nb + 1
                idx_end   =  pp      *nb

                do g = 1, mesh%n_groups
                    b = total_source(idx_start:idx_end, g)

                    if (PD%matrix_free) then
                        do f = 1, mesh%n_faces_per_elem
                            o_n = dot_product(dir, PD%face_normals(:,f,pp))
                            if (o_n >= 0.0_dp) cycle  ! outflow
                            n_face = PD%n_face_basis_f(f)
                            if (PD%face_connectivity(1,f,pp) > 0) then
                                do j_f = 1, n_face
                                    upwind_flux(j_f) = ang_flux(PD%upwind_idx(j_f,f,pp), mm, g)
                                end do
                            else if (PD%face_connectivity(4,f,pp) > 0 .and. &
                                     any(PD%face_connectivity(4,f,pp) == ref_ID)) then
                                m_ref = PD%reflect_map(mm, f, pp)
                                do j_f = 1, n_face
                                    local_dof = idx_start - 1 + PD%face_node_map_patch(j_f, f)
                                    upwind_flux(j_f) = ang_flux_snap(local_dof, m_ref, g)
                                end do
                            else
                                cycle  ! vacuum
                            end if
                            fm = dir(1)*PD%face_mass_x(:,:,f,pp) + &
                                 dir(2)*PD%face_mass_y(:,:,f,pp) + &
                                 dir(3)*PD%face_mass_z(:,:,f,pp)
                            do j_f = 1, n_face
                                b = b - upwind_flux(j_f) * fm(:, PD%face_node_map_patch(j_f,f))
                            end do
                        end do
                    else
                        do f = 1, mesh%n_faces_per_elem
                            n_face = PD%n_face_basis_f(f)
                            if (PD%face_connectivity(1,f,pp) > 0) then
                                ! Interior face: gather upwind flux from neighbor patch
                                do j_f = 1, n_face
                                    upwind_flux(j_f) = ang_flux(PD%upwind_idx(j_f,f,pp), mm, g)
                                end do
                            else if (PD%face_connectivity(4,f,pp) > 0 .and. &
                                     any(PD%face_connectivity(4,f,pp) == ref_ID)) then
                                ! Reflective BC — read from snapshot to avoid OMP race
                                m_ref = PD%reflect_map(mm, f, pp)
                                do j_f = 1, n_face
                                    local_dof = idx_start - 1 + PD%face_node_map_patch(j_f, f)
                                    upwind_flux(j_f) = ang_flux_snap(local_dof, m_ref, g)
                                end do
                            else
                                cycle  ! vacuum: zero upwind flux, no contribution
                            end if
                            ! face_mass_in is zero for outflow spans, so this is always safe
                            do j_f = 1, n_face
                                b = b - upwind_flux(j_f) * &
                                        PD%face_mass_in(:, PD%face_node_map_patch(j_f,f), f, pp, mm)
                            end do
                        end do
                    end if

                    call dgetrs('N', nb, 1, PD%local_lu(:,:,pp,mm,g), nb, &
                                PD%local_pivots(:,pp,mm,g), b, nb, info)

                    where (b < 0.0_dp) b = 0.0_dp
                    if (any(b /= b)) then
                        write(*,'(A,2I6)') "FATAL: NaN in patch sweep, patch/angle=", pp, mm; stop
                    end if

                    ang_flux(idx_start:idx_end, mm, g) = b
                    do k = 1, nb
                        !$OMP ATOMIC
                        scalar_flux(idx_start + k - 1, g) = scalar_flux(idx_start + k - 1, g) + w_mm * b(k)
                    end do
                end do
            end do
        end do
        !$OMP END PARALLEL DO

        if (allocated(ang_flux_snap)) deallocate(ang_flux_snap)
    end subroutine Transport_Sweep_Patch

    ! ------------------------------------------------------------------
    ! Build total source (scatter + fission) at patch level.
    ! ------------------------------------------------------------------
    subroutine Source_Patch_DGFEM(total_src, scalar_flux, k_eff, materials, mesh, &
                                   PD, n_patches, n_groups, is_adjoint, is_eigenvalue)
        real(dp),              intent(inout) :: total_src(:,:)
        real(dp),              intent(in)    :: scalar_flux(:,:)
        real(dp),              intent(in)    :: k_eff
        type(t_material),      intent(in)    :: materials(:)
        type(t_mesh_iga),      intent(in)    :: mesh
        type(t_patch_dg), intent(in) :: PD
        integer,               intent(in)    :: n_patches, n_groups
        logical,               intent(in)    :: is_adjoint, is_eigenvalue

        integer  :: pp, g_to, mat_id, ee_rep, idx_start, idx_end, nb
        real(dp) :: M_phi(PD%n_basis_patch, n_groups)
        real(dp) :: fission_rate(PD%n_basis_patch)

        nb = PD%n_basis_patch
        total_src = 0.0_dp

        !$OMP PARALLEL DO PRIVATE(pp, mat_id, ee_rep, idx_start, idx_end, M_phi, fission_rate, g_to)
        do pp = 1, n_patches
            ee_rep    = PD%patch_elem_list(PD%patch_elem_start(pp))
            mat_id    = mesh%material_ids(ee_rep)
            idx_start = (pp - 1)*nb + 1
            idx_end   =  pp      *nb

            M_phi = matmul(PD%patch_mass(:,:,pp), scalar_flux(idx_start:idx_end, :))

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
                        fission_rate * merge(materials(mat_id)%Chi(g_to), &
                                             materials(mat_id)%NuSigF(g_to), .not. is_adjoint)
                end do
            else
                do g_to = 1, n_groups
                    total_src(idx_start:idx_end,g_to) = total_src(idx_start:idx_end,g_to) + &
                        materials(mat_id)%Src(g_to) * PD%basis_integrals_vol(:, pp)
                end do
            end if

            where (total_src(idx_start:idx_end,:) < 0.0_dp) total_src(idx_start:idx_end,:) = 0.0_dp
            if (any(total_src(idx_start:idx_end,:) /= total_src(idx_start:idx_end,:))) then
                write(*,'(A,I6)') "FATAL: NaN in patch source, patch=", pp; stop
            end if
        end do
        !$OMP END PARALLEL DO
    end subroutine Source_Patch_DGFEM

    ! ------------------------------------------------------------------
    ! Fission production integral for k-eff.
    ! ------------------------------------------------------------------
    subroutine Calculate_Production_Patch_DGFEM(total_prod, scalar_flux, materials, &
                                                  mesh, PD, n_patches, is_adjoint)
        real(dp),              intent(out) :: total_prod
        real(dp),              intent(in)  :: scalar_flux(:,:)
        type(t_material),      intent(in)  :: materials(:)
        type(t_mesh_iga),      intent(in)  :: mesh
        type(t_patch_dg), intent(in) :: PD
        integer,               intent(in)  :: n_patches
        logical,               intent(in)  :: is_adjoint

        integer :: g, pp, mat_id, ee_rep, idx_start, idx_end, nb

        nb = PD%n_basis_patch
        total_prod = 0.0_dp

        !$OMP PARALLEL DO PRIVATE(pp, mat_id, ee_rep, idx_start, idx_end, g) REDUCTION(+:total_prod)
        do pp = 1, n_patches
            ee_rep    = PD%patch_elem_list(PD%patch_elem_start(pp))
            mat_id    = mesh%material_ids(ee_rep)
            idx_start = (pp - 1)*nb + 1
            idx_end   =  pp      *nb
            do g = 1, mesh%n_groups
                total_prod = total_prod + &
                    merge(materials(mat_id)%NuSigF(g), materials(mat_id)%Chi(g), .not. is_adjoint) * &
                    dot_product(max(scalar_flux(idx_start:idx_end,g), 0.0_dp), PD%basis_integrals_vol(:,pp))
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine Calculate_Production_Patch_DGFEM

    ! ------------------------------------------------------------------
    ! Map patchwise scalar flux to element (knot-span) DOF layout so the
    ! existing VTK export routines can be reused.
    !   flux_patch(pp-1)*nb + pdof, g)  →  elem_flux((ee-1)*nb_elem + a, g)
    ! ------------------------------------------------------------------
    subroutine Remap_Patch_To_Elem_Flux(flux_patch, elem_flux, mesh, FE, PD)
        real(dp),              intent(in)  :: flux_patch(:,:)
        real(dp),              intent(out) :: elem_flux(:,:)
        type(t_mesh_iga),      intent(in)  :: mesh
        type(t_basis_iga),    intent(in)  :: FE
        type(t_patch_dg), intent(in) :: PD

        integer :: ee, pp, a, pdof, n_groups, g, nb

        nb       = PD%n_basis_patch
        n_groups = size(flux_patch, 2)
        do ee = 1, mesh%n_elems
            pp = mesh%elem_patch_id(ee)
            do a = 1, FE%n_basis
                pdof = PD%elem_to_patch_dof(a, ee)
                do g = 1, n_groups
                    elem_flux((ee-1)*FE%n_basis + a, g) = flux_patch((pp-1)*nb + pdof, g)
                end do
            end do
        end do
    end subroutine Remap_Patch_To_Elem_Flux

    ! ==================================================================
    ! High-level entry.  Wires module state, delegates to PowerIteration.
    ! ==================================================================
    subroutine SolveTransport_IGA(mesh, materials, sn_quad, PD, &
                                        scalar_flux, ang_flux_out, k_eff, &
                                        sweep_order, ref_ids, max_outer, tol, &
                                        is_adjoint, is_eigenvalue)
        type(t_mesh_iga),      intent(in),    target :: mesh
        type(t_material),      intent(in),    target :: materials(:)
        type(t_sn_quadrature), intent(in),    target :: sn_quad
        type(t_patch_dg), intent(in), target :: PD
        real(dp), allocatable, intent(out)           :: scalar_flux(:,:)
        real(dp), allocatable, intent(out)           :: ang_flux_out(:,:,:)
        real(dp),              intent(out)           :: k_eff
        integer,               intent(in)            :: sweep_order(:,:)
        integer,               intent(in)            :: ref_ids(:)
        integer,               intent(in)            :: max_outer
        real(dp),              intent(in)            :: tol
        logical,               intent(in)            :: is_adjoint, is_eigenvalue

        integer :: n_dof, n_patches

        n_patches  = size(sweep_order, 1)   ! correct for both FDG (n_elems) and PDG (n_patches)
        n_dof      = n_patches * PD%n_basis_patch
        s_n_groups = mesh%n_groups
        s_n_patches= n_patches

        s_mesh    => mesh
        s_sn_quad => sn_quad
        s_PD      => PD
        s_mats    => materials

        if (allocated(s_ang_flux))    deallocate(s_ang_flux)
        if (allocated(s_sweep_order)) deallocate(s_sweep_order)
        if (allocated(s_ref_ids))     deallocate(s_ref_ids)

        allocate(s_ang_flux(n_dof, sn_quad%n_angles, s_n_groups), source=0.0_dp)
        allocate(s_sweep_order, source=sweep_order)
        allocate(s_ref_ids,     source=ref_ids)
        allocate(scalar_flux(n_dof, s_n_groups), source=1.0_dp)

        call PowerIteration(scalar_flux, k_eff, max_outer, tol, &
                             is_eigenvalue, is_adjoint, &
                             patch_source, patch_solve, patch_production)

        call move_alloc(s_ang_flux, ang_flux_out)
    end subroutine SolveTransport_IGA

    ! ---- PowerIteration callbacks ------------------------------------

    subroutine patch_source(src, flux, k_eff, is_eigenvalue, is_adjoint)
        real(dp), intent(inout) :: src(:,:)
        real(dp), intent(in)    :: flux(:,:), k_eff
        logical,  intent(in)    :: is_eigenvalue, is_adjoint
        call Source_Patch_DGFEM(src, flux, k_eff, s_mats, s_mesh, s_PD, &
                                 s_n_patches, s_n_groups, is_adjoint, is_eigenvalue)
    end subroutine patch_source

    subroutine patch_solve(flux, src)
        real(dp), intent(inout) :: flux(:,:)
        real(dp), intent(in)    :: src(:,:)
        call Transport_Sweep_Patch(s_mesh, s_sn_quad, s_PD, &
                                    s_ang_flux, flux, src, s_sweep_order, s_ref_ids)
    end subroutine patch_solve

    subroutine patch_production(prod, flux, is_adjoint)
        real(dp), intent(out) :: prod
        real(dp), intent(in)  :: flux(:,:)
        logical,  intent(in)  :: is_adjoint
        call Calculate_Production_Patch_DGFEM(prod, flux, s_mats, s_mesh, s_PD, &
                                               s_n_patches, is_adjoint)
    end subroutine patch_production

end module m_transport_iga_pdg
