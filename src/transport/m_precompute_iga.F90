! IGA transport precomputation — span-DG and patch-DG modes.
!
! Both modes share the same t_patch_dg solver (m_transport_iga_pdg).
! SpanDG promotes each knot span to a trivial single-span "patch" via
! promote_spans_to_patches, so the sweep kernel sees an identical layout.
! PatchDG assembles span contributions into multi-span patches with C^(p-1)
! continuity within each patch and DG upwind coupling at patch interfaces.
!
! Shared private helpers (no duplication):
!   precompute_integrals  -- NURBS volume/face quadrature
!   precompute_patch_lu   -- per-patch LU factorisation
!
! Public:
!   InitialiseTransport_SpanDG   -- span-DG entry
!   InitialiseTransport_PatchDG  -- patch-DG entry
module m_transport_precompute_iga
    use m_constants
    use m_types
    use m_quadrature
    use m_basis_iga,   only: GetMapping2D, GetMapping3D
    use m_material
    use m_sweep_order, only: connectivity_and_normals, precompute_upwind_indices, &
                              precompute_reflective_map, generate_sweep_order
    implicit none
    private
    public :: InitialiseTransport_SpanDG, InitialiseTransport_PatchDG

    interface
        subroutine dgetrf(m, n, a, lda, ipiv, info)
            import :: dp
            integer, intent(in)    :: m, n, lda
            real(dp), intent(inout):: a(lda, *)
            integer, intent(out)   :: ipiv(*), info
        end subroutine dgetrf
    end interface

contains

    ! ==================================================================
    ! SPAN-DG (FDG): each knot span is its own DG element.
    ! Promotes span data to t_patch_dg with identity DOF map so the
    ! unified patch solver handles both modes.
    ! ==================================================================
    subroutine InitialiseTransport_SpanDG(mesh, FE, sn_quad, QuadVol, QuadFace, &
                                           materials, PD, sweep_order)
        type(t_mesh_iga),      intent(in)    :: mesh
        type(t_basis_iga),     intent(in)    :: FE
        type(t_sn_quadrature), intent(in)    :: sn_quad
        type(t_quadrature),    intent(in)    :: QuadVol, QuadFace
        type(t_material),      intent(in)    :: materials(:)
        type(t_patch_dg),      intent(inout) :: PD
        integer, allocatable,  intent(out)   :: sweep_order(:,:)

        type(t_fem_dg) :: TD
        integer  :: mm
        real(dp) :: dir_tmp(3)

        write(*,'(A,I0,A,I0)') "  [SpanDG]  n_spans=", mesh%n_elems, "  dim=", mesh%dim

        call precompute_integrals(mesh, FE, QuadVol, QuadFace, TD)
        call connectivity_and_normals(mesh, FE, QuadFace, TD)
        call precompute_upwind_indices(mesh%n_elems, mesh%n_faces_per_elem, &
                                       FE%n_nodes_per_face, FE%n_basis, &
                                       FE%face_node_map, TD%face_connectivity, TD%upwind_idx)
        call precompute_reflective_map(mesh%n_elems, mesh%n_faces_per_elem, &
                                       sn_quad, TD%face_normals, TD%reflect_map)

        call promote_spans_to_patches(mesh, FE, sn_quad, TD, PD)
        call precompute_patch_lu(mesh, sn_quad, materials, mesh%n_groups, mesh%n_elems, PD)

        allocate(sweep_order(mesh%n_elems, sn_quad%n_angles))
        do mm = 1, sn_quad%n_angles
            dir_tmp = sn_quad%dirs(mm, :)
            call generate_sweep_order(mesh%n_elems, mesh%n_faces_per_elem, &
                                      PD%face_normals, PD%face_connectivity, &
                                      dir_tmp, sweep_order(:, mm))
        end do

        write(*,'(A)') "  [SpanDG]  Precompute complete."
    end subroutine InitialiseTransport_SpanDG

    ! ==================================================================
    ! PATCH-DG (PDG): each NURBS patch is one DG element.
    ! Span contributions are assembled into patch matrices; C^(p-1)
    ! continuity is retained within each patch.
    ! ==================================================================
    subroutine InitialiseTransport_PatchDG(mesh, FE, sn_quad, QuadVol, QuadFace, &
                                            materials, PD, sweep_order_patch)
        type(t_mesh_iga),      intent(in)    :: mesh
        type(t_basis_iga),     intent(in)    :: FE
        type(t_sn_quadrature), intent(in)    :: sn_quad
        type(t_quadrature),    intent(in)    :: QuadVol, QuadFace
        type(t_material),      intent(in)    :: materials(:)
        type(t_patch_dg),      intent(inout) :: PD
        integer, allocatable,  intent(out)   :: sweep_order_patch(:,:)

        type(t_fem_dg) :: TD_elem
        integer :: n_patches, mm
        real(dp) :: dir_tmp(3)

        n_patches = size(mesh%patches)
        write(*,'(A,I0,A,I0,A,I0)') "  [PatchDG]  n_patches=", n_patches, &
                                    "  n_spans=", mesh%n_elems, "  dim=", mesh%dim
        if (mesh%n_elems / n_patches > 8) then
            write(*,'(A,I0,A)') "  [PatchDG]  WARNING: ", mesh%n_elems / n_patches, &
                " spans/patch.  PDG may be unstable for optically thin patches."
            write(*,'(A)') "             Consider using FDG (span-DG) for fine meshes."
        end if

        call precompute_integrals(mesh, FE, QuadVol, QuadFace, TD_elem)
        call connectivity_and_normals(mesh, FE, QuadFace, TD_elem)

        call setup_patch_dof_info(mesh, FE, n_patches, PD)
        call assemble_patch_matrices(mesh, FE, sn_quad, TD_elem, n_patches, PD)
        call build_patch_connectivity(mesh, TD_elem, n_patches, PD)
        call precompute_patch_upwind(mesh, n_patches, PD)
        call precompute_reflective_map(n_patches, mesh%n_faces_per_elem, &
                                        sn_quad, PD%face_normals, PD%reflect_map)
        call precompute_patch_lu(mesh, sn_quad, materials, mesh%n_groups, n_patches, PD)

        allocate(sweep_order_patch(n_patches, sn_quad%n_angles))
        do mm = 1, sn_quad%n_angles
            dir_tmp = sn_quad%dirs(mm, 1:3)
            call generate_sweep_order(n_patches, mesh%n_faces_per_elem, &
                                       PD%face_normals, PD%face_connectivity, &
                                       dir_tmp, sweep_order_patch(:, mm))
        end do

        write(*,'(A)') "  [PatchDG]  Precompute complete."
    end subroutine InitialiseTransport_PatchDG

    ! ==================================================================
    ! SHARED PRIVATE HELPERS
    ! ==================================================================

    ! ------------------------------------------------------------------
    ! Volume mass/stiffness matrices and face mass matrices (both modes).
    ! ------------------------------------------------------------------
    subroutine precompute_integrals(mesh, FE, Quad, QuadFace, TD)
        type(t_mesh_iga),   intent(in)    :: mesh
        type(t_basis_iga),  intent(in)    :: FE
        type(t_quadrature), intent(in)    :: Quad, QuadFace
        type(t_fem_dg),     intent(inout) :: TD

        integer  :: ee, q, f, nf
        real(dp) :: nodes(FE%n_basis, 3)
        real(dp) :: dN_dx(FE%n_basis), dN_dy(FE%n_basis), dN_dz(FE%n_basis)
        real(dp) :: detJ, dV, R(FE%n_basis)
        real(dp) :: xi_f, eta_f, zeta_f, J(3,3), J2(2,2), dA(3), s1, s2
        real(dp) :: u1, u2, v1, v2, w1, w2

        nf = mesh%n_faces_per_elem

        allocate(TD%elem_mass_matrix(FE%n_basis, FE%n_basis, mesh%n_elems), &
                 TD%elem_stiffness_x(FE%n_basis, FE%n_basis, mesh%n_elems), &
                 TD%elem_stiffness_y(FE%n_basis, FE%n_basis, mesh%n_elems), &
                 TD%elem_stiffness_z(FE%n_basis, FE%n_basis, mesh%n_elems), &
                 TD%face_mass_x(FE%n_basis, FE%n_basis, nf, mesh%n_elems),  &
                 TD%face_mass_y(FE%n_basis, FE%n_basis, nf, mesh%n_elems),  &
                 TD%face_mass_z(FE%n_basis, FE%n_basis, nf, mesh%n_elems),  &
                 TD%basis_integrals_vol(FE%n_basis, mesh%n_elems))

        TD%elem_mass_matrix = 0.0_dp
        TD%elem_stiffness_x = 0.0_dp; TD%elem_stiffness_y = 0.0_dp; TD%elem_stiffness_z = 0.0_dp
        TD%face_mass_x = 0.0_dp; TD%face_mass_y = 0.0_dp; TD%face_mass_z = 0.0_dp
        TD%basis_integrals_vol = 0.0_dp

        !$OMP PARALLEL DO PRIVATE(ee, nodes, u1, u2, v1, v2, w1, w2, q, &
        !$OMP&   dN_dx, dN_dy, dN_dz, detJ, R, dV, f, xi_f, eta_f, zeta_f, J, J2, dA, s1, s2)
        do ee = 1, mesh%n_elems
            nodes = mesh%nodes(mesh%elems(ee, 1:FE%n_basis), :)
            u1=mesh%elem_u_min(ee); u2=mesh%elem_u_max(ee)
            v1=mesh%elem_v_min(ee); v2=mesh%elem_v_max(ee)
            if (mesh%dim == 3) then
                w1=mesh%elem_w_min(ee); w2=mesh%elem_w_max(ee)
            else
                w1 = 0.0_dp; w2 = 1.0_dp
            end if

            do q = 1, Quad%n_points
                if (mesh%dim == 3) then
                    call GetMapping3D(FE, ee, mesh, q, Quad, u1, u2, v1, v2, w1, w2, nodes, &
                                      dN_dx, dN_dy, dN_dz, detJ, R)
                else
                    call GetMapping2D(FE, ee, mesh, q, Quad, u1, u2, v1, v2, nodes(:,1:2), &
                                      dN_dx, dN_dy, detJ, R)
                    dN_dz = 0.0_dp
                end if
                dV = detJ * Quad%weights(q)
                TD%elem_mass_matrix(:,:,ee) = TD%elem_mass_matrix(:,:,ee) + &
                    spread(R,2,FE%n_basis) * spread(R,1,FE%n_basis) * dV
                TD%elem_stiffness_x(:,:,ee) = TD%elem_stiffness_x(:,:,ee) + &
                    spread(dN_dx,2,FE%n_basis) * spread(R,1,FE%n_basis) * dV
                TD%elem_stiffness_y(:,:,ee) = TD%elem_stiffness_y(:,:,ee) + &
                    spread(dN_dy,2,FE%n_basis) * spread(R,1,FE%n_basis) * dV
                TD%elem_stiffness_z(:,:,ee) = TD%elem_stiffness_z(:,:,ee) + &
                    spread(dN_dz,2,FE%n_basis) * spread(R,1,FE%n_basis) * dV
                TD%basis_integrals_vol(:,ee) = TD%basis_integrals_vol(:,ee) + R * dV
            end do

            do f = 1, nf
                do q = 1, QuadFace%n_points
                    if (mesh%dim == 3) then
                        select case(f)
                            case(1); xi_f=QuadFace%xi(q); eta_f=QuadFace%eta(q); zeta_f=-1.0_dp
                            case(2); xi_f=QuadFace%xi(q); eta_f=QuadFace%eta(q); zeta_f= 1.0_dp
                            case(3); xi_f=QuadFace%xi(q); eta_f=-1.0_dp;         zeta_f=QuadFace%eta(q)
                            case(4); xi_f=QuadFace%xi(q); eta_f= 1.0_dp;         zeta_f=QuadFace%eta(q)
                            case(5); xi_f=-1.0_dp;        eta_f=QuadFace%xi(q);  zeta_f=QuadFace%eta(q)
                            case(6); xi_f= 1.0_dp;        eta_f=QuadFace%xi(q);  zeta_f=QuadFace%eta(q)
                        end select
                        call GetMapping3D(FE, ee, mesh, q, Quad, u1, u2, v1, v2, w1, w2, nodes, &
                                          dN_dx, dN_dy, dN_dz, detJ, R, &
                                          xi_custom=xi_f, eta_custom=eta_f, zeta_custom=zeta_f, J_out=J)
                        select case(f)
                            case(1,2)
                                s1=0.5_dp*(u2-u1); s2=0.5_dp*(v2-v1)
                                dA(1)=(J(1,2)*J(2,3)-J(1,3)*J(2,2))*s1*s2
                                dA(2)=(J(1,3)*J(2,1)-J(1,1)*J(2,3))*s1*s2
                                dA(3)=(J(1,1)*J(2,2)-J(1,2)*J(2,1))*s1*s2
                                if (f==1) dA=-dA
                            case(3,4)
                                s1=0.5_dp*(w2-w1); s2=0.5_dp*(u2-u1)
                                dA(1)=(J(3,2)*J(1,3)-J(3,3)*J(1,2))*s1*s2
                                dA(2)=(J(3,3)*J(1,1)-J(3,1)*J(1,3))*s1*s2
                                dA(3)=(J(3,1)*J(1,2)-J(3,2)*J(1,1))*s1*s2
                                if (f==3) dA=-dA
                            case(5,6)
                                s1=0.5_dp*(v2-v1); s2=0.5_dp*(w2-w1)
                                dA(1)=(J(2,2)*J(3,3)-J(2,3)*J(3,2))*s1*s2
                                dA(2)=(J(2,3)*J(3,1)-J(2,1)*J(3,3))*s1*s2
                                dA(3)=(J(2,1)*J(3,2)-J(2,2)*J(3,1))*s1*s2
                                if (f==5) dA=-dA
                        end select
                    else
                        select case(f)
                            case(1); xi_f=QuadFace%xi(q); eta_f=-1.0_dp
                            case(2); xi_f= 1.0_dp;        eta_f=QuadFace%xi(q)
                            case(3); xi_f=QuadFace%xi(q); eta_f= 1.0_dp
                            case(4); xi_f=-1.0_dp;        eta_f=QuadFace%xi(q)
                        end select
                        call GetMapping2D(FE, ee, mesh, q, Quad, u1, u2, v1, v2, nodes(:,1:2), &
                                          dN_dx, dN_dy, detJ, R, &
                                          xi_custom=xi_f, eta_custom=eta_f, J_out=J2)
                        dA(2) = 0.0_dp
                        select case(f)
                            case(1,3)
                                s1 = 0.5_dp*(u2-u1)
                                dA(1) =  J2(1,2) * s1
                                dA(2) = -J2(1,1) * s1
                                if (f==3) dA = -dA
                            case(2,4)
                                s1 = 0.5_dp*(v2-v1)
                                dA(1) =  J2(2,2) * s1
                                dA(2) = -J2(2,1) * s1
                                if (f==4) dA = -dA
                        end select
                    end if
                    dA = dA * QuadFace%weights(q)
                    TD%face_mass_x(:,:,f,ee) = TD%face_mass_x(:,:,f,ee) + &
                        spread(R,2,FE%n_basis)*spread(R,1,FE%n_basis)*dA(1)
                    TD%face_mass_y(:,:,f,ee) = TD%face_mass_y(:,:,f,ee) + &
                        spread(R,2,FE%n_basis)*spread(R,1,FE%n_basis)*dA(2)
                    TD%face_mass_z(:,:,f,ee) = TD%face_mass_z(:,:,f,ee) + &
                        spread(R,2,FE%n_basis)*spread(R,1,FE%n_basis)*dA(3)
                end do
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine precompute_integrals

    ! ------------------------------------------------------------------
    ! Per-patch LU factorisation.  Works for both SpanDG (n_patches=n_elems)
    ! and PatchDG (n_patches=size(mesh%patches)).
    ! ------------------------------------------------------------------
    subroutine precompute_patch_lu(mesh, sn_quad, materials, n_groups, n_patches, PD)
        type(t_mesh_iga),      intent(in)    :: mesh
        type(t_sn_quadrature), intent(in)    :: sn_quad
        type(t_material),      intent(in)    :: materials(:)
        integer,               intent(in)    :: n_groups, n_patches
        type(t_patch_dg),      intent(inout) :: PD

        integer  :: pp, mm, g, f, info, mat_id, ee_rep, nb, nf
        real(dp) :: dir(3), o_n
        real(dp) :: Stiff(PD%n_basis_patch, PD%n_basis_patch)
        real(dp) :: A    (PD%n_basis_patch, PD%n_basis_patch)

        nb = PD%n_basis_patch;  nf = mesh%n_faces_per_elem

        do pp = 1, n_patches
            ee_rep = PD%patch_elem_list(PD%patch_elem_start(pp))
            mat_id = mesh%material_ids(ee_rep)
            if (.not. allocated(materials(mat_id)%SigmaT)) then
                write(*,'(A,I0)') "FATAL: patch material has no SigmaT, mat_id=", mat_id; stop
            end if
        end do

        allocate(PD%local_lu    (nb, nb, n_patches, sn_quad%n_angles, n_groups))
        allocate(PD%local_pivots(nb,     n_patches, sn_quad%n_angles, n_groups))

        !$OMP PARALLEL DO PRIVATE(mm, dir, pp, Stiff, ee_rep, mat_id, f, o_n, g, A, info)
        do mm = 1, sn_quad%n_angles
            dir = sn_quad%dirs(mm, 1:3)
            do pp = 1, n_patches
                Stiff = -(dir(1)*PD%patch_stiff_x(:,:,pp) + &
                           dir(2)*PD%patch_stiff_y(:,:,pp) + &
                           dir(3)*PD%patch_stiff_z(:,:,pp))
                if (PD%matrix_free) then
                    do f = 1, nf
                        o_n = dot_product(dir, PD%face_normals(:,f,pp))
                        if (o_n > 0.0_dp) &
                            Stiff = Stiff + dir(1)*PD%face_mass_x(:,:,f,pp) + &
                                            dir(2)*PD%face_mass_y(:,:,f,pp) + &
                                            dir(3)*PD%face_mass_z(:,:,f,pp)
                    end do
                else
                    do f = 1, nf
                        Stiff = Stiff + PD%face_mass_out(:,:,f,pp,mm)
                    end do
                end if
                ee_rep = PD%patch_elem_list(PD%patch_elem_start(pp))
                mat_id = mesh%material_ids(ee_rep)
                do g = 1, n_groups
                    A = materials(mat_id)%SigmaT(g) * PD%patch_mass(:,:,pp) + Stiff
                    call dgetrf(nb, nb, A, nb, PD%local_pivots(:,pp,mm,g), info)
                    if (info /= 0) then
                        write(*,'(A,2I5)') "FATAL: LU failed, patch/span,angle=", pp, mm; stop
                    end if
                    PD%local_lu(:,:,pp,mm,g) = A
                end do
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine precompute_patch_lu

    ! ==================================================================
    ! SPAN-DG PRIVATE: promotes span-level t_fem_dg to t_patch_dg with
    ! identity DOF map (one span = one patch).
    ! ==================================================================
    subroutine promote_spans_to_patches(mesh, FE, sn_quad, TD, PD)
        type(t_mesh_iga),      intent(in)    :: mesh
        type(t_basis_iga),     intent(in)    :: FE
        type(t_sn_quadrature), intent(in)    :: sn_quad
        type(t_fem_dg),        intent(in)    :: TD
        type(t_patch_dg),      intent(inout) :: PD

        integer  :: ee, a, f, mm, nf, nb
        real(dp) :: dir(3), o_n

        nb = FE%n_basis
        nf = mesh%n_faces_per_elem

        PD%n_basis_patch    = nb
        PD%n_face_basis_max = FE%n_nodes_per_face

        allocate(PD%n_face_basis_f(nf),        source=FE%n_nodes_per_face)
        allocate(PD%face_node_map_patch(FE%n_nodes_per_face, nf), source=FE%face_node_map)

        allocate(PD%elem_to_patch_dof(nb, mesh%n_elems))
        do ee = 1, mesh%n_elems
            PD%elem_to_patch_dof(:, ee) = [(a, a=1,nb)]
        end do

        allocate(PD%patch_elem_start(mesh%n_elems + 1))
        allocate(PD%patch_elem_list (mesh%n_elems))
        do ee = 1, mesh%n_elems
            PD%patch_elem_start(ee) = ee
            PD%patch_elem_list(ee)  = ee
        end do
        PD%patch_elem_start(mesh%n_elems + 1) = mesh%n_elems + 1

        allocate(PD%patch_mass      (nb, nb, mesh%n_elems), source=TD%elem_mass_matrix)
        allocate(PD%patch_stiff_x   (nb, nb, mesh%n_elems), source=TD%elem_stiffness_x)
        allocate(PD%patch_stiff_y   (nb, nb, mesh%n_elems), source=TD%elem_stiffness_y)
        allocate(PD%patch_stiff_z   (nb, nb, mesh%n_elems), source=TD%elem_stiffness_z)
        allocate(PD%basis_integrals_vol(nb, mesh%n_elems),  source=TD%basis_integrals_vol)

        allocate(PD%face_connectivity(4, nf, mesh%n_elems), source=TD%face_connectivity)
        allocate(PD%face_normals     (3, nf, mesh%n_elems), source=TD%face_normals)
        allocate(PD%upwind_idx(FE%n_nodes_per_face, nf, mesh%n_elems), source=TD%upwind_idx)
        allocate(PD%reflect_map(sn_quad%n_angles, nf, mesh%n_elems), source=TD%reflect_map)

        if (PD%matrix_free) then
            ! Mode 2: keep angle-independent components; compute dir·face_mass inline in sweep
            allocate(PD%face_mass_x(nb, nb, nf, mesh%n_elems), source=TD%face_mass_x)
            allocate(PD%face_mass_y(nb, nb, nf, mesh%n_elems), source=TD%face_mass_y)
            allocate(PD%face_mass_z(nb, nb, nf, mesh%n_elems), source=TD%face_mass_z)
        else
            ! Full precompute: per-angle per-span inflow/outflow split
            allocate(PD%face_mass_out(nb, nb, nf, mesh%n_elems, sn_quad%n_angles), source=0.0_dp)
            allocate(PD%face_mass_in (nb, nb, nf, mesh%n_elems, sn_quad%n_angles), source=0.0_dp)
            !$OMP PARALLEL DO PRIVATE(mm, dir, ee, f, o_n)
            do mm = 1, sn_quad%n_angles
                dir = sn_quad%dirs(mm,:)
                do ee = 1, mesh%n_elems
                    do f = 1, nf
                        o_n = dot_product(dir, TD%face_normals(:,f,ee))
                        if (o_n >= 0.0_dp) then
                            PD%face_mass_out(:,:,f,ee,mm) = dir(1)*TD%face_mass_x(:,:,f,ee) + &
                                                             dir(2)*TD%face_mass_y(:,:,f,ee) + &
                                                             dir(3)*TD%face_mass_z(:,:,f,ee)
                        else
                            PD%face_mass_in(:,:,f,ee,mm)  = dir(1)*TD%face_mass_x(:,:,f,ee) + &
                                                             dir(2)*TD%face_mass_y(:,:,f,ee) + &
                                                             dir(3)*TD%face_mass_z(:,:,f,ee)
                        end if
                    end do
                end do
            end do
            !$OMP END PARALLEL DO
        end if
    end subroutine promote_spans_to_patches

    ! ==================================================================
    ! PATCH-DG PRIVATE HELPERS
    ! ==================================================================

    subroutine setup_patch_dof_info(mesh, FE, n_patches, PD)
        type(t_mesh_iga),  intent(in)    :: mesh
        type(t_basis_iga), intent(in)    :: FE
        integer,           intent(in)    :: n_patches
        type(t_patch_dg),  intent(inout) :: PD

        integer :: pp, ee, ii, jj, kk, p, q, r
        integer :: local_idx, cp_xi, cp_eta, cp_zeta, patch_dof
        integer :: span_xi, span_eta, span_zeta, p_idx
        integer :: n_cp_xi, n_cp_eta, n_cp_zeta, nf
        integer :: n_fb_12, n_fb_34, n_fb_56
        integer, allocatable :: patch_count(:)

        p  = FE%p_order;  q  = FE%q_order
        r  = merge(FE%r_order, 0, mesh%dim == 3)
        nf = mesh%n_faces_per_elem

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

        allocate(PD%face_node_map_patch(PD%n_face_basis_max, nf), source=0)
        if (mesh%dim == 2) then
            do ii = 1, n_cp_xi;  PD%face_node_map_patch(ii, 1) = ii; end do
            do jj = 1, n_cp_eta; PD%face_node_map_patch(jj, 2) = jj * n_cp_xi; end do
            do ii = 1, n_cp_xi;  PD%face_node_map_patch(ii, 3) = (n_cp_eta - 1)*n_cp_xi + ii; end do
            do jj = 1, n_cp_eta; PD%face_node_map_patch(jj, 4) = (jj - 1)*n_cp_xi + 1; end do
        else
            local_idx = 0
            do jj = 1, n_cp_eta; do ii = 1, n_cp_xi
                local_idx = local_idx + 1
                PD%face_node_map_patch(local_idx, 1) = (jj-1)*n_cp_xi + ii
            end do; end do
            local_idx = 0
            do jj = 1, n_cp_eta; do ii = 1, n_cp_xi
                local_idx = local_idx + 1
                PD%face_node_map_patch(local_idx, 2) = &
                    (n_cp_zeta-1)*n_cp_xi*n_cp_eta + (jj-1)*n_cp_xi + ii
            end do; end do
            local_idx = 0
            do kk = 1, n_cp_zeta; do ii = 1, n_cp_xi
                local_idx = local_idx + 1
                PD%face_node_map_patch(local_idx, 3) = (kk-1)*n_cp_xi*n_cp_eta + ii
            end do; end do
            local_idx = 0
            do kk = 1, n_cp_zeta; do ii = 1, n_cp_xi
                local_idx = local_idx + 1
                PD%face_node_map_patch(local_idx, 4) = &
                    (kk-1)*n_cp_xi*n_cp_eta + (n_cp_eta-1)*n_cp_xi + ii
            end do; end do
            local_idx = 0
            do kk = 1, n_cp_zeta; do jj = 1, n_cp_eta
                local_idx = local_idx + 1
                PD%face_node_map_patch(local_idx, 5) = &
                    (kk-1)*n_cp_xi*n_cp_eta + (jj-1)*n_cp_xi + 1
            end do; end do
            local_idx = 0
            do kk = 1, n_cp_zeta; do jj = 1, n_cp_eta
                local_idx = local_idx + 1
                PD%face_node_map_patch(local_idx, 6) = &
                    (kk-1)*n_cp_xi*n_cp_eta + (jj-1)*n_cp_xi + n_cp_xi
            end do; end do
        end if
    end subroutine setup_patch_dof_info

    subroutine assemble_patch_matrices(mesh, FE, sn_quad, TD, n_patches, PD)
        type(t_mesh_iga),      intent(in)    :: mesh
        type(t_basis_iga),     intent(in)    :: FE
        type(t_sn_quadrature), intent(in)    :: sn_quad
        type(t_fem_dg),        intent(in)    :: TD
        integer,               intent(in)    :: n_patches
        type(t_patch_dg),      intent(inout) :: PD

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
        allocate(PD%face_mass_out(PD%n_basis_patch, PD%n_basis_patch, nf, n_patches, sn_quad%n_angles), source=0.0_dp)
        allocate(PD%face_mass_in (PD%n_basis_patch, PD%n_basis_patch, nf, n_patches, sn_quad%n_angles), source=0.0_dp)

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

        do ee = 1, mesh%n_elems
            pp = mesh%elem_patch_id(ee)
            sx = mesh%elem_span_indices(1,ee); sy = mesh%elem_span_indices(2,ee)
            if (mesh%dim == 2) then
                if (sy == mn_sy(pp)) call accum_face_split(TD, FE, sn_quad, PD, ee, pp, 1)
                if (sx == mx_sx(pp)) call accum_face_split(TD, FE, sn_quad, PD, ee, pp, 2)
                if (sy == mx_sy(pp)) call accum_face_split(TD, FE, sn_quad, PD, ee, pp, 3)
                if (sx == mn_sx(pp)) call accum_face_split(TD, FE, sn_quad, PD, ee, pp, 4)
            else
                sz = mesh%elem_span_indices(3,ee)
                if (sz == mn_sz(pp)) call accum_face_split(TD, FE, sn_quad, PD, ee, pp, 1)
                if (sz == mx_sz(pp)) call accum_face_split(TD, FE, sn_quad, PD, ee, pp, 2)
                if (sy == mn_sy(pp)) call accum_face_split(TD, FE, sn_quad, PD, ee, pp, 3)
                if (sy == mx_sy(pp)) call accum_face_split(TD, FE, sn_quad, PD, ee, pp, 4)
                if (sx == mn_sx(pp)) call accum_face_split(TD, FE, sn_quad, PD, ee, pp, 5)
                if (sx == mx_sx(pp)) call accum_face_split(TD, FE, sn_quad, PD, ee, pp, 6)
            end if
        end do

        deallocate(mn_sx, mx_sx, mn_sy, mx_sy)
        if (mesh%dim == 3) deallocate(mn_sz, mx_sz)
    end subroutine assemble_patch_matrices

    subroutine accum_face_split(TD, FE, sn_quad, PD, ee, pp, f)
        type(t_fem_dg),        intent(in)    :: TD
        type(t_basis_iga),     intent(in)    :: FE
        type(t_sn_quadrature), intent(in)    :: sn_quad
        type(t_patch_dg),      intent(inout) :: PD
        integer, intent(in) :: ee, pp, f
        integer  :: a, b, ia, ib, mm
        real(dp) :: dir(3), o_n, fm_ab
        do mm = 1, sn_quad%n_angles
            dir = sn_quad%dirs(mm,:)
            o_n = dot_product(dir, TD%face_normals(:,f,ee))
            do b = 1, FE%n_basis
                ib = PD%elem_to_patch_dof(b, ee)
                do a = 1, FE%n_basis
                    ia = PD%elem_to_patch_dof(a, ee)
                    fm_ab = dir(1)*TD%face_mass_x(a,b,f,ee) + &
                            dir(2)*TD%face_mass_y(a,b,f,ee) + &
                            dir(3)*TD%face_mass_z(a,b,f,ee)
                    if (o_n >= 0.0_dp) then
                        PD%face_mass_out(ia,ib,f,pp,mm) = PD%face_mass_out(ia,ib,f,pp,mm) + fm_ab
                    else
                        PD%face_mass_in (ia,ib,f,pp,mm) = PD%face_mass_in (ia,ib,f,pp,mm) + fm_ab
                    end if
                end do
            end do
        end do
    end subroutine accum_face_split

    subroutine build_patch_connectivity(mesh, TD, n_patches, PD)
        type(t_mesh_iga),  intent(in)    :: mesh
        type(t_fem_dg),    intent(in)    :: TD
        integer,           intent(in)    :: n_patches
        type(t_patch_dg),  intent(inout) :: PD

        integer :: ee, pp, nf, f
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

        do pp = 1, n_patches
            do f = 1, nf
                n3 = PD%face_normals(:,f,pp)
                if (norm2(n3) > dp_EPSILON) PD%face_normals(:,f,pp) = n3 / norm2(n3)
            end do
        end do

        do pp = 1, n_patches
            do f = 1, nf
                if (PD%face_connectivity(1,f,pp) > 0) then
                    if (PD%face_connectivity(1, PD%face_connectivity(2,f,pp), &
                                             PD%face_connectivity(1,f,pp)) /= pp) &
                        write(*,'(A,3I5)') "[WARN] Asymmetric patch connectivity: pp,f,pp2=", &
                            pp, f, PD%face_connectivity(1,f,pp)
                end if
            end do
        end do

        deallocate(mn_sx, mx_sx, mn_sy, mx_sy)
        if (mesh%dim == 3) deallocate(mn_sz, mx_sz)
    end subroutine build_patch_connectivity

    subroutine proc_face(mesh, TD, PD, ee, pp, f)
        type(t_mesh_iga),  intent(in)    :: mesh
        type(t_fem_dg),    intent(in)    :: TD
        type(t_patch_dg),  intent(inout) :: PD
        integer, intent(in) :: ee, pp, f
        integer :: ee2

        PD%face_normals(:,f,pp) = PD%face_normals(:,f,pp) + TD%face_normals(:,f,ee)
        if (PD%face_connectivity(1,f,pp) /= -1) return

        ee2 = TD%face_connectivity(1,f,ee)
        if (ee2 > 0) then
            PD%face_connectivity(1,f,pp) = mesh%elem_patch_id(ee2)
            PD%face_connectivity(2,f,pp) = TD%face_connectivity(2,f,ee)
            PD%face_connectivity(3,f,pp) = TD%face_connectivity(3,f,ee)
            PD%face_connectivity(4,f,pp) = 0
        else
            PD%face_connectivity(1,f,pp) = -1
            PD%face_connectivity(4,f,pp) = TD%face_connectivity(4,f,ee)
        end if
    end subroutine proc_face

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
                    do j_f = 1, n_face
                        PD%upwind_idx(j_f,f,pp) = (pp-1)*nb + PD%face_node_map_patch(j_f,f)
                    end do
                end if
            end do
        end do
    end subroutine precompute_patch_upwind

end module m_transport_precompute_iga
