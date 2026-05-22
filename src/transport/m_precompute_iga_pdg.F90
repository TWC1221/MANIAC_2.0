! Patchwise DG-IGA transport precomputation.
! Each NURBS patch is treated as a single DG "element".  Within a patch the
! full C^(p-1) IGA space is retained; DG upwind coupling occurs only at patch
! interfaces.  Volume/face matrices are assembled from element-level (knot-span)
! contributions using standard FEM assembly.
!
! DOF layout: patch pp → global DOFs  (pp-1)*n_basis_patch + 1 .. pp*n_basis_patch
!
! Public:
!   InitialiseTransport_PatchDG  -- master entry
module m_transport_precompute_patchdg
    use m_constants
    use m_types
    use m_types_iga
    use m_quadrature
    use m_material
    use m_sweep_order,              only: connectivity_and_normals, generate_sweep_order, &
                                          precompute_reflective_map
    use m_transport_precompute_iga, only: precompute_integrals, precompute_patch_lu
    implicit none
    private
    public :: InitialiseTransport_PatchDG

contains

    ! ------------------------------------------------------------------
    ! Master entry.  Builds all patch-level data and per-angle sweep order.
    ! ------------------------------------------------------------------
    subroutine InitialiseTransport_PatchDG(mesh, FE, sn_quad, QuadVol, QuadFace, &
                                            materials, PD, sweep_order_patch)
        type(t_mesh_iga),      intent(in)    :: mesh
        type(t_basis_iga),    intent(in)    :: FE
        type(t_sn_quadrature), intent(in)    :: sn_quad
        type(t_quadrature),    intent(in)    :: QuadVol, QuadFace
        type(t_material),      intent(in)    :: materials(:)
        type(t_patch_dg), intent(inout) :: PD
        integer, allocatable,  intent(out)   :: sweep_order_patch(:,:)

        type(t_fem_dg) :: TD_elem
        integer :: n_patches, mm

        n_patches = size(mesh%patches)
        write(*,'(A,I0,A,I0,A,I0)') "  [PatchDG]  n_patches=", n_patches, &
                                    "  n_spans=", mesh%n_elems, "  dim=", mesh%dim
        ! Warn if spans-per-patch is high: CG-within-patch becomes unstable at
        ! high Peclet number (|Omega|/(SigmaT*h_span) >> 1).  Rule of thumb:
        ! keep SigmaT * L_patch > ~2 (optically thick) for stable convergence.
        if (mesh%n_elems / n_patches > 8) then
            write(*,'(A,I0,A)') "  [PatchDG]  WARNING: ", mesh%n_elems / n_patches, &
                " spans/patch.  PDG may be unstable for optically thin patches."
            write(*,'(A)') "             Consider using FDG (span-DG) for fine meshes."
        end if

        call precompute_integrals(mesh, FE, QuadVol, QuadFace, TD_elem)
        call connectivity_and_normals(mesh, FE, QuadFace, TD_elem)

        call setup_patch_dof_info(mesh, FE, n_patches, PD)
        call assemble_patch_matrices(mesh, FE, TD_elem, n_patches, PD)
        call build_patch_connectivity(mesh, TD_elem, n_patches, PD)
        call precompute_patch_upwind(mesh, n_patches, PD)
        call precompute_reflective_map(n_patches, mesh%n_faces_per_elem, &
                                        sn_quad, PD%face_normals, PD%reflect_map)
        call precompute_patch_lu(mesh, sn_quad, materials, mesh%n_groups, n_patches, PD)

        allocate(sweep_order_patch(n_patches, sn_quad%n_angles))
        do mm = 1, sn_quad%n_angles
            call generate_sweep_order(n_patches, mesh%n_faces_per_elem, &
                                       PD%face_normals, PD%face_connectivity, &
                                       sn_quad%dirs(mm, 1:3), sweep_order_patch(:, mm))
        end do
    end subroutine InitialiseTransport_PatchDG

    ! ------------------------------------------------------------------
    ! Patch DOF metadata, span→patch mapping, face node map.
    ! ------------------------------------------------------------------
    subroutine setup_patch_dof_info(mesh, FE, n_patches, PD)
        type(t_mesh_iga),    intent(in)    :: mesh
        type(t_basis_iga),  intent(in)    :: FE
        integer,             intent(in)    :: n_patches
        type(t_patch_dg), intent(inout) :: PD

        integer :: pp, ee, ii, jj, kk, p, q, r
        integer :: local_idx, cp_xi, cp_eta, cp_zeta, patch_dof
        integer :: span_xi, span_eta, span_zeta, p_idx
        integer :: n_cp_xi, n_cp_eta, n_cp_zeta, nf
        integer :: n_fb_12, n_fb_34, n_fb_56
        integer, allocatable :: patch_count(:)

        p  = FE%p_order;  q  = FE%q_order
        r  = merge(FE%r_order, 0, mesh%dim == 3)
        nf = mesh%n_faces_per_elem

        ! ---- Patch element lists -----------------------------------
        allocate(patch_count(n_patches), source=0)
        do ee = 1, mesh%n_elems
            patch_count(mesh%elem_patch_id(ee)) = patch_count(mesh%elem_patch_id(ee)) + 1
        end do
        allocate(PD%patch_elem_start(n_patches+1))
        PD%patch_elem_start(1) = 1
        do pp = 1, n_patches
            PD%patch_elem_start(pp+1) = PD%patch_elem_start(pp) + patch_count(pp)
        end do
        allocate(PD%patch_elem_list(mesh%n_elems))
        patch_count = 0
        do ee = 1, mesh%n_elems
            pp = mesh%elem_patch_id(ee)
            PD%patch_elem_list(PD%patch_elem_start(pp) + patch_count(pp)) = ee
            patch_count(pp) = patch_count(pp) + 1
        end do
        deallocate(patch_count)

        ! ---- CP counts from patch 1 (uniform assumption) -----------
        p_idx    = mesh%elem_patch_id(PD%patch_elem_list(1))
        n_cp_xi  = size(mesh%patches(p_idx)%knots_xi)  - (p + 1)
        n_cp_eta = size(mesh%patches(p_idx)%knots_eta) - (q + 1)
        n_cp_zeta = 1
        if (mesh%dim == 3) n_cp_zeta = size(mesh%patches(p_idx)%knots_zeta) - (r + 1)

        do pp = 1, n_patches
            p_idx = mesh%elem_patch_id(PD%patch_elem_list(PD%patch_elem_start(pp)))
            if (size(mesh%patches(p_idx)%knots_xi)  - (p+1) /= n_cp_xi  .or. &
                size(mesh%patches(p_idx)%knots_eta) - (q+1) /= n_cp_eta) then
                write(*,'(A)') "FATAL: non-uniform patch CP counts in patchwise DG-IGA"; stop
            end if
            if (mesh%dim == 3) then
                if (size(mesh%patches(p_idx)%knots_zeta) - (r+1) /= n_cp_zeta) then
                    write(*,'(A)') "FATAL: non-uniform patch CP counts (zeta)"; stop
                end if
            end if
        end do

        PD%n_basis_patch = n_cp_xi * n_cp_eta * n_cp_zeta

        ! ---- Face basis counts -------------------------------------
        allocate(PD%n_face_basis_f(nf))
        if (mesh%dim == 2) then
            PD%n_face_basis_f(1) = n_cp_xi;  PD%n_face_basis_f(2) = n_cp_eta
            PD%n_face_basis_f(3) = n_cp_xi;  PD%n_face_basis_f(4) = n_cp_eta
            PD%n_face_basis_max  = max(n_cp_xi, n_cp_eta)
        else
            n_fb_12 = n_cp_xi * n_cp_eta
            n_fb_34 = n_cp_xi * n_cp_zeta
            n_fb_56 = n_cp_eta * n_cp_zeta
            PD%n_face_basis_f(1) = n_fb_12; PD%n_face_basis_f(2) = n_fb_12
            PD%n_face_basis_f(3) = n_fb_34; PD%n_face_basis_f(4) = n_fb_34
            PD%n_face_basis_f(5) = n_fb_56; PD%n_face_basis_f(6) = n_fb_56
            PD%n_face_basis_max  = max(n_fb_12, n_fb_34, n_fb_56)
        end if

        ! ---- Span → patch DOF mapping ------------------------------
        allocate(PD%elem_to_patch_dof(FE%n_basis, mesh%n_elems))
        do ee = 1, mesh%n_elems
            span_xi  = mesh%elem_span_indices(1, ee)
            span_eta = mesh%elem_span_indices(2, ee)
            local_idx = 0
            if (mesh%dim == 2) then
                do jj = 1, q+1
                    do ii = 1, p+1
                        local_idx = local_idx + 1
                        cp_xi  = span_xi  - p + ii - 1
                        cp_eta = span_eta - q + jj - 1
                        PD%elem_to_patch_dof(local_idx, ee) = (cp_eta - 1)*n_cp_xi + cp_xi
                    end do
                end do
            else
                span_zeta = mesh%elem_span_indices(3, ee)
                do kk = 1, r+1
                    do jj = 1, q+1
                        do ii = 1, p+1
                            local_idx = local_idx + 1
                            cp_xi   = span_xi   - p + ii - 1
                            cp_eta  = span_eta  - q + jj - 1
                            cp_zeta = span_zeta - r + kk - 1
                            PD%elem_to_patch_dof(local_idx, ee) = &
                                (cp_zeta-1)*n_cp_xi*n_cp_eta + (cp_eta-1)*n_cp_xi + cp_xi
                        end do
                    end do
                end do
            end if
        end do

        ! Consistency check against global CP IDs
        do ee = 1, mesh%n_elems
            pp = mesh%elem_patch_id(ee)
            do local_idx = 1, FE%n_basis
                patch_dof = PD%elem_to_patch_dof(local_idx, ee)
                if (patch_dof < 1 .or. patch_dof > PD%n_basis_patch) then
                    write(*,'(A,3I8)') "FATAL: patch DOF out of range, ee,a,pdof=", &
                        ee, local_idx, patch_dof; stop
                end if
                if (mesh%elems(ee, local_idx) /= mesh%patches(pp)%cp_ids(patch_dof)) then
                    write(*,'(A)') "FATAL: span-to-patch DOF ordering mismatch — " // &
                        "check ASMG cp_ids ordering matches xi-fastest convention"; stop
                end if
            end do
        end do

        ! ---- Patch face node map -----------------------------------
        allocate(PD%face_node_map_patch(PD%n_face_basis_max, nf), source=0)
        if (mesh%dim == 2) then
            ! Face 1 (eta=min, cp_eta=1): n_cp_xi nodes
            do ii = 1, n_cp_xi
                PD%face_node_map_patch(ii, 1) = ii
            end do
            ! Face 2 (xi=max, cp_xi=n_cp_xi): n_cp_eta nodes
            do jj = 1, n_cp_eta
                PD%face_node_map_patch(jj, 2) = jj * n_cp_xi
            end do
            ! Face 3 (eta=max, cp_eta=n_cp_eta): n_cp_xi nodes
            do ii = 1, n_cp_xi
                PD%face_node_map_patch(ii, 3) = (n_cp_eta - 1)*n_cp_xi + ii
            end do
            ! Face 4 (xi=min, cp_xi=1): n_cp_eta nodes
            do jj = 1, n_cp_eta
                PD%face_node_map_patch(jj, 4) = (jj - 1)*n_cp_xi + 1
            end do
        else
            ! Face 1 (zeta=min): (cp_xi,cp_eta) with cp_zeta=1
            local_idx = 0
            do jj = 1, n_cp_eta
                do ii = 1, n_cp_xi
                    local_idx = local_idx + 1
                    PD%face_node_map_patch(local_idx, 1) = (jj-1)*n_cp_xi + ii
                end do
            end do
            ! Face 2 (zeta=max): cp_zeta=n_cp_zeta
            local_idx = 0
            do jj = 1, n_cp_eta
                do ii = 1, n_cp_xi
                    local_idx = local_idx + 1
                    PD%face_node_map_patch(local_idx, 2) = &
                        (n_cp_zeta-1)*n_cp_xi*n_cp_eta + (jj-1)*n_cp_xi + ii
                end do
            end do
            ! Face 3 (eta=min): (cp_xi,cp_zeta) with cp_eta=1
            local_idx = 0
            do kk = 1, n_cp_zeta
                do ii = 1, n_cp_xi
                    local_idx = local_idx + 1
                    PD%face_node_map_patch(local_idx, 3) = (kk-1)*n_cp_xi*n_cp_eta + ii
                end do
            end do
            ! Face 4 (eta=max): cp_eta=n_cp_eta
            local_idx = 0
            do kk = 1, n_cp_zeta
                do ii = 1, n_cp_xi
                    local_idx = local_idx + 1
                    PD%face_node_map_patch(local_idx, 4) = &
                        (kk-1)*n_cp_xi*n_cp_eta + (n_cp_eta-1)*n_cp_xi + ii
                end do
            end do
            ! Face 5 (xi=min): (cp_eta,cp_zeta) with cp_xi=1
            local_idx = 0
            do kk = 1, n_cp_zeta
                do jj = 1, n_cp_eta
                    local_idx = local_idx + 1
                    PD%face_node_map_patch(local_idx, 5) = &
                        (kk-1)*n_cp_xi*n_cp_eta + (jj-1)*n_cp_xi + 1
                end do
            end do
            ! Face 6 (xi=max): cp_xi=n_cp_xi
            local_idx = 0
            do kk = 1, n_cp_zeta
                do jj = 1, n_cp_eta
                    local_idx = local_idx + 1
                    PD%face_node_map_patch(local_idx, 6) = &
                        (kk-1)*n_cp_xi*n_cp_eta + (jj-1)*n_cp_xi + n_cp_xi
                end do
            end do
        end if
    end subroutine setup_patch_dof_info

    ! ------------------------------------------------------------------
    ! Assemble patch-level volume and boundary-face matrices from TD.
    ! ------------------------------------------------------------------
    subroutine assemble_patch_matrices(mesh, FE, TD, n_patches, PD)
        type(t_mesh_iga),      intent(in)    :: mesh
        type(t_basis_iga),    intent(in)    :: FE
        type(t_fem_dg), intent(in)   :: TD
        integer,               intent(in)    :: n_patches
        type(t_patch_dg), intent(inout) :: PD

        integer :: ee, pp, a, b, ia, ib, nf
        integer :: sx, sy, sz
        integer, allocatable :: mn_sx(:), mx_sx(:), mn_sy(:), mx_sy(:)
        integer, allocatable :: mn_sz(:), mx_sz(:)

        nf = mesh%n_faces_per_elem

        allocate(PD%patch_mass      (PD%n_basis_patch, PD%n_basis_patch, n_patches), source=0.0_dp)
        allocate(PD%patch_stiff_x   (PD%n_basis_patch, PD%n_basis_patch, n_patches), source=0.0_dp)
        allocate(PD%patch_stiff_y   (PD%n_basis_patch, PD%n_basis_patch, n_patches), source=0.0_dp)
        allocate(PD%patch_stiff_z   (PD%n_basis_patch, PD%n_basis_patch, n_patches), source=0.0_dp)
        allocate(PD%basis_integrals_vol(PD%n_basis_patch, n_patches),                source=0.0_dp)
        allocate(PD%face_mass_x(PD%n_basis_patch, PD%n_basis_patch, nf, n_patches),  source=0.0_dp)
        allocate(PD%face_mass_y(PD%n_basis_patch, PD%n_basis_patch, nf, n_patches),  source=0.0_dp)
        allocate(PD%face_mass_z(PD%n_basis_patch, PD%n_basis_patch, nf, n_patches),  source=0.0_dp)

        ! Per-patch span bounds
        allocate(mn_sx(n_patches), source= huge(0)); allocate(mx_sx(n_patches), source=-huge(0))
        allocate(mn_sy(n_patches), source= huge(0)); allocate(mx_sy(n_patches), source=-huge(0))
        if (mesh%dim == 3) then
            allocate(mn_sz(n_patches), source= huge(0)); allocate(mx_sz(n_patches), source=-huge(0))
        end if
        do ee = 1, mesh%n_elems
            pp = mesh%elem_patch_id(ee)
            sx = mesh%elem_span_indices(1,ee); sy = mesh%elem_span_indices(2,ee)
            mn_sx(pp) = min(mn_sx(pp),sx); mx_sx(pp) = max(mx_sx(pp),sx)
            mn_sy(pp) = min(mn_sy(pp),sy); mx_sy(pp) = max(mx_sy(pp),sy)
            if (mesh%dim == 3) then
                sz = mesh%elem_span_indices(3,ee)
                mn_sz(pp) = min(mn_sz(pp),sz); mx_sz(pp) = max(mx_sz(pp),sz)
            end if
        end do

        ! Volume assembly
        do ee = 1, mesh%n_elems
            pp = mesh%elem_patch_id(ee)
            do b = 1, FE%n_basis
                ib = PD%elem_to_patch_dof(b, ee)
                PD%basis_integrals_vol(ib,pp) = PD%basis_integrals_vol(ib,pp) + TD%basis_integrals_vol(b,ee)
                do a = 1, FE%n_basis
                    ia = PD%elem_to_patch_dof(a, ee)
                    PD%patch_mass(ia,ib,pp)    = PD%patch_mass(ia,ib,pp)    + TD%elem_mass_matrix(a,b,ee)
                    PD%patch_stiff_x(ia,ib,pp) = PD%patch_stiff_x(ia,ib,pp) + TD%elem_stiffness_x(a,b,ee)
                    PD%patch_stiff_y(ia,ib,pp) = PD%patch_stiff_y(ia,ib,pp) + TD%elem_stiffness_y(a,b,ee)
                    PD%patch_stiff_z(ia,ib,pp) = PD%patch_stiff_z(ia,ib,pp) + TD%elem_stiffness_z(a,b,ee)
                end do
            end do
        end do

        ! Patch-boundary face assembly
        do ee = 1, mesh%n_elems
            pp = mesh%elem_patch_id(ee)
            sx = mesh%elem_span_indices(1,ee); sy = mesh%elem_span_indices(2,ee)

            if (mesh%dim == 2) then
                if (sy == mn_sy(pp)) call accum_face(TD, FE, PD, ee, pp, 1)
                if (sx == mx_sx(pp)) call accum_face(TD, FE, PD, ee, pp, 2)
                if (sy == mx_sy(pp)) call accum_face(TD, FE, PD, ee, pp, 3)
                if (sx == mn_sx(pp)) call accum_face(TD, FE, PD, ee, pp, 4)
            else
                sz = mesh%elem_span_indices(3,ee)
                if (sz == mn_sz(pp)) call accum_face(TD, FE, PD, ee, pp, 1)
                if (sz == mx_sz(pp)) call accum_face(TD, FE, PD, ee, pp, 2)
                if (sy == mn_sy(pp)) call accum_face(TD, FE, PD, ee, pp, 3)
                if (sy == mx_sy(pp)) call accum_face(TD, FE, PD, ee, pp, 4)
                if (sx == mn_sx(pp)) call accum_face(TD, FE, PD, ee, pp, 5)
                if (sx == mx_sx(pp)) call accum_face(TD, FE, PD, ee, pp, 6)
            end if
        end do

        deallocate(mn_sx, mx_sx, mn_sy, mx_sy)
        if (mesh%dim == 3) deallocate(mn_sz, mx_sz)
    end subroutine assemble_patch_matrices

    ! Accumulate element face mass matrices into patch face mass.
    subroutine accum_face(TD, FE, PD, ee, pp, f)
        type(t_fem_dg), intent(in)         :: TD
        type(t_basis_iga),     intent(in)         :: FE
        type(t_patch_dg), intent(inout) :: PD
        integer, intent(in) :: ee, pp, f
        integer :: a, b, ia, ib
        do b = 1, FE%n_basis
            ib = PD%elem_to_patch_dof(b, ee)
            do a = 1, FE%n_basis
                ia = PD%elem_to_patch_dof(a, ee)
                PD%face_mass_x(ia,ib,f,pp) = PD%face_mass_x(ia,ib,f,pp) + TD%face_mass_x(a,b,f,ee)
                PD%face_mass_y(ia,ib,f,pp) = PD%face_mass_y(ia,ib,f,pp) + TD%face_mass_y(a,b,f,ee)
                PD%face_mass_z(ia,ib,f,pp) = PD%face_mass_z(ia,ib,f,pp) + TD%face_mass_z(a,b,f,ee)
            end do
        end do
    end subroutine accum_face

    ! ------------------------------------------------------------------
    ! Derive patch face connectivity and normals from element-level TD.
    ! ------------------------------------------------------------------
    subroutine build_patch_connectivity(mesh, TD, n_patches, PD)
        type(t_mesh_iga),      intent(in)    :: mesh
        type(t_fem_dg), intent(in)   :: TD
        integer,               intent(in)    :: n_patches
        type(t_patch_dg), intent(inout) :: PD

        integer :: ee, pp, pp2, ee2, nf, f
        integer :: sx, sy, sz
        integer, allocatable :: mn_sx(:), mx_sx(:), mn_sy(:), mx_sy(:)
        integer, allocatable :: mn_sz(:), mx_sz(:)
        real(dp) :: n3(3)

        nf = mesh%n_faces_per_elem

        allocate(PD%face_connectivity(4, nf, n_patches))
        allocate(PD%face_normals(3, nf, n_patches), source=0.0_dp)
        PD%face_connectivity(1,:,:) = -1
        PD%face_connectivity(2,:,:) = -1
        PD%face_connectivity(3,:,:) =  1
        PD%face_connectivity(4,:,:) = BC_VACUUM

        allocate(mn_sx(n_patches), source= huge(0)); allocate(mx_sx(n_patches), source=-huge(0))
        allocate(mn_sy(n_patches), source= huge(0)); allocate(mx_sy(n_patches), source=-huge(0))
        if (mesh%dim == 3) then
            allocate(mn_sz(n_patches), source= huge(0)); allocate(mx_sz(n_patches), source=-huge(0))
        end if
        do ee = 1, mesh%n_elems
            pp = mesh%elem_patch_id(ee)
            sx = mesh%elem_span_indices(1,ee); sy = mesh%elem_span_indices(2,ee)
            mn_sx(pp) = min(mn_sx(pp),sx); mx_sx(pp) = max(mx_sx(pp),sx)
            mn_sy(pp) = min(mn_sy(pp),sy); mx_sy(pp) = max(mx_sy(pp),sy)
            if (mesh%dim == 3) then
                sz = mesh%elem_span_indices(3,ee)
                mn_sz(pp) = min(mn_sz(pp),sz); mx_sz(pp) = max(mx_sz(pp),sz)
            end if
        end do

        do ee = 1, mesh%n_elems
            pp = mesh%elem_patch_id(ee)
            sx = mesh%elem_span_indices(1,ee); sy = mesh%elem_span_indices(2,ee)

            if (mesh%dim == 2) then
                if (sy == mn_sy(pp)) call proc_face(mesh, TD, PD, ee, pp, 1)
                if (sx == mx_sx(pp)) call proc_face(mesh, TD, PD, ee, pp, 2)
                if (sy == mx_sy(pp)) call proc_face(mesh, TD, PD, ee, pp, 3)
                if (sx == mn_sx(pp)) call proc_face(mesh, TD, PD, ee, pp, 4)
            else
                sz = mesh%elem_span_indices(3,ee)
                if (sz == mn_sz(pp)) call proc_face(mesh, TD, PD, ee, pp, 1)
                if (sz == mx_sz(pp)) call proc_face(mesh, TD, PD, ee, pp, 2)
                if (sy == mn_sy(pp)) call proc_face(mesh, TD, PD, ee, pp, 3)
                if (sy == mx_sy(pp)) call proc_face(mesh, TD, PD, ee, pp, 4)
                if (sx == mn_sx(pp)) call proc_face(mesh, TD, PD, ee, pp, 5)
                if (sx == mx_sx(pp)) call proc_face(mesh, TD, PD, ee, pp, 6)
            end if
        end do

        ! Normalise averaged normals
        do pp = 1, n_patches
            do f = 1, nf
                n3 = PD%face_normals(:,f,pp)
                if (norm2(n3) > dp_EPSILON) PD%face_normals(:,f,pp) = n3 / norm2(n3)
            end do
        end do

        ! Symmetry check
        do pp = 1, n_patches
            do f = 1, nf
                pp2 = PD%face_connectivity(1,f,pp)
                if (pp2 > 0) then
                    if (PD%face_connectivity(1, PD%face_connectivity(2,f,pp), pp2) /= pp) &
                        write(*,'(A,3I5)') "[WARN] Asymmetric patch connectivity: pp,f,pp2=", pp,f,pp2
                end if
            end do
        end do

        deallocate(mn_sx, mx_sx, mn_sy, mx_sy)
        if (mesh%dim == 3) deallocate(mn_sz, mx_sz)
    end subroutine build_patch_connectivity

    ! Accumulate normal and set connectivity for one patch boundary (element,face).
    subroutine proc_face(mesh, TD, PD, ee, pp, f)
        type(t_mesh_iga),      intent(in)    :: mesh
        type(t_fem_dg), intent(in)   :: TD
        type(t_patch_dg), intent(inout) :: PD
        integer, intent(in) :: ee, pp, f
        integer :: ee2

        ! Accumulate normal (average over all boundary elements on this patch face)
        PD%face_normals(:,f,pp) = PD%face_normals(:,f,pp) + TD%face_normals(:,f,ee)

        ! Derive connectivity once (use first contributing element)
        if (PD%face_connectivity(1,f,pp) /= -1) return

        ee2 = TD%face_connectivity(1,f,ee)
        if (ee2 > 0) then
            PD%face_connectivity(1,f,pp) = mesh%elem_patch_id(ee2)
            PD%face_connectivity(2,f,pp) = TD%face_connectivity(2,f,ee)
            PD%face_connectivity(3,f,pp) = TD%face_connectivity(3,f,ee)
            PD%face_connectivity(4,f,pp) = 0
        else
            ! Physical boundary: inherit BC from element
            PD%face_connectivity(1,f,pp) = -1
            PD%face_connectivity(4,f,pp) = TD%face_connectivity(4,f,ee)
        end if
    end subroutine proc_face

    ! ------------------------------------------------------------------
    ! Patch upwind DOF indices.
    ! Interior: match shared CP IDs into neighbor patch's DOF space.
    ! Boundary: own face DOFs (for reflective BC; vacuum needs no upwind).
    ! ------------------------------------------------------------------
    subroutine precompute_patch_upwind(mesh, n_patches, PD)
        type(t_mesh_iga), intent(in)    :: mesh
        integer,          intent(in)    :: n_patches
        type(t_patch_dg), intent(inout) :: PD

        integer :: pp, pp2, f, j_f, j_nb, gcp, n_face, nb

        nb = PD%n_basis_patch
        allocate(PD%upwind_idx(PD%n_face_basis_max, mesh%n_faces_per_elem, n_patches), source=0)

        do pp = 1, n_patches
            do f = 1, mesh%n_faces_per_elem
                n_face = PD%n_face_basis_f(f)
                pp2    = PD%face_connectivity(1,f,pp)

                if (pp2 > 0) then
                    do j_f = 1, n_face
                        gcp = mesh%patches(pp)%cp_ids(PD%face_node_map_patch(j_f, f))
                        do j_nb = 1, nb
                            if (mesh%patches(pp2)%cp_ids(j_nb) == gcp) then
                                PD%upwind_idx(j_f,f,pp) = (pp2-1)*nb + j_nb
                                exit
                            end if
                        end do
                        if (PD%upwind_idx(j_f,f,pp) == 0) then
                            write(*,'(A,3I5)') "FATAL: upwind CP not found, pp,f,j_f=",pp,f,j_f; stop
                        end if
                    end do
                else
                    ! Boundary: own face DOF (only used for reflective BC)
                    do j_f = 1, n_face
                        PD%upwind_idx(j_f,f,pp) = (pp-1)*nb + PD%face_node_map_patch(j_f,f)
                    end do
                end if
            end do
        end do
    end subroutine precompute_patch_upwind

end module m_transport_precompute_patchdg
