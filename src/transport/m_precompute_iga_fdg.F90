! IGA transport precomputation — span (knot-element) DG mode.
!
! Computes NURBS volume/face integrals then promotes each knot span to a
! trivial "patch" (identity DOF map, n_patches = n_elems) so both span-DG
! and patch-DG share the same solver (m_transport_iga_pdg / m_transport_iga).
!
! precompute_patch_lu is also kept here as a public utility because it
! operates purely on t_patch_dg and is reused by m_precompute_iga_pdg.
!
! Public:
!   InitialiseTransport_SpanDG  -- master entry; fills t_patch_dg + sweep order
!   precompute_integrals        -- NURBS volume/face quadrature (also used by pdg)
!   precompute_patch_lu         -- per-patch LU factorisation    (also used by pdg)
module m_transport_precompute_iga
    use m_constants
    use m_types
    use m_types_iga
    use m_quadrature
    use m_basis_iga,   only: GetMapping2D, GetMapping3D
    use m_material
    use m_sweep_order, only: connectivity_and_normals, precompute_upwind_indices, &
                              precompute_reflective_map, generate_sweep_order
    implicit none
    public :: InitialiseTransport_SpanDG, precompute_integrals, precompute_patch_lu

    interface
        subroutine dgetrf(m, n, a, lda, ipiv, info)
            import :: dp
            integer, intent(in)    :: m, n, lda
            real(dp), intent(inout):: a(lda, *)
            integer, intent(out)   :: ipiv(*), info
        end subroutine dgetrf
    end interface

contains

    ! ------------------------------------------------------------------
    ! Master entry for span-DG.  Wraps the full precompute chain and
    ! returns a t_patch_dg where each knot span is a degenerate "patch".
    ! ------------------------------------------------------------------
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

        call promote_spans_to_patches(mesh, FE, sn_quad%n_angles, TD, PD)
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

    ! ------------------------------------------------------------------
    ! Promote span-level t_fem_dg data to a t_patch_dg with one span
    ! per patch and an identity DOF map.  The two layouts are numerically
    ! identical; this lets the unified IGA solver handle both modes.
    ! ------------------------------------------------------------------
    subroutine promote_spans_to_patches(mesh, FE, n_angles, TD, PD)
        type(t_mesh_iga),  intent(in)    :: mesh
        type(t_basis_iga), intent(in)    :: FE
        integer,           intent(in)    :: n_angles
        type(t_fem_dg),    intent(in)    :: TD
        type(t_patch_dg),  intent(inout) :: PD

        integer :: ee, a, nf, nb

        nb = FE%n_basis
        nf = mesh%n_faces_per_elem

        PD%n_basis_patch    = nb
        PD%n_face_basis_max = FE%n_nodes_per_face

        allocate(PD%n_face_basis_f(nf),        source=FE%n_nodes_per_face)
        allocate(PD%face_node_map_patch(FE%n_nodes_per_face, nf), source=FE%face_node_map)

        ! Identity DOF map: each span ee is its own "patch", DOF a → patch DOF a
        allocate(PD%elem_to_patch_dof(nb, mesh%n_elems))
        do ee = 1, mesh%n_elems
            PD%elem_to_patch_dof(:, ee) = [(a, a=1,nb)]
        end do

        ! Trivial element lists: patch pp contains exactly span pp
        allocate(PD%patch_elem_start(mesh%n_elems + 1))
        allocate(PD%patch_elem_list (mesh%n_elems))
        do ee = 1, mesh%n_elems
            PD%patch_elem_start(ee) = ee
            PD%patch_elem_list(ee)  = ee
        end do
        PD%patch_elem_start(mesh%n_elems + 1) = mesh%n_elems + 1

        ! Volume matrices — direct copy
        allocate(PD%patch_mass      (nb, nb, mesh%n_elems), source=TD%elem_mass_matrix)
        allocate(PD%patch_stiff_x   (nb, nb, mesh%n_elems), source=TD%elem_stiffness_x)
        allocate(PD%patch_stiff_y   (nb, nb, mesh%n_elems), source=TD%elem_stiffness_y)
        allocate(PD%patch_stiff_z   (nb, nb, mesh%n_elems), source=TD%elem_stiffness_z)
        allocate(PD%basis_integrals_vol(nb, mesh%n_elems),  source=TD%basis_integrals_vol)

        ! Face mass matrices — direct copy
        allocate(PD%face_mass_x(nb, nb, nf, mesh%n_elems), source=TD%face_mass_x)
        allocate(PD%face_mass_y(nb, nb, nf, mesh%n_elems), source=TD%face_mass_y)
        allocate(PD%face_mass_z(nb, nb, nf, mesh%n_elems), source=TD%face_mass_z)

        ! Topology — direct copy (patch index = span index throughout)
        allocate(PD%face_connectivity(4, nf, mesh%n_elems), source=TD%face_connectivity)
        allocate(PD%face_normals     (3, nf, mesh%n_elems), source=TD%face_normals)
        allocate(PD%upwind_idx(FE%n_nodes_per_face, nf, mesh%n_elems), source=TD%upwind_idx)
        allocate(PD%reflect_map(n_angles, nf, mesh%n_elems), source=TD%reflect_map)
    end subroutine promote_spans_to_patches

    ! ------------------------------------------------------------------
    ! Volume mass/stiffness matrices and face mass matrices.
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
                                dA(1) =  J2(1,2) * s1;  dA(2) = -J2(1,1) * s1
                                if (f==3) dA = -dA
                            case(2,4)
                                s1 = 0.5_dp*(v2-v1)
                                dA(1) =  J2(2,2) * s1;  dA(2) = -J2(2,1) * s1
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
    ! Per-patch LU factorisation.  Generic over t_patch_dg: works for
    ! span-DG (n_patches = n_elems) and multi-span patch-DG alike.
    ! Automatic (non-allocatable) local arrays so OMP PRIVATE is valid.
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
                do f = 1, nf
                    o_n = dot_product(dir, PD%face_normals(:,f,pp))
                    if (o_n > 0.0_dp) &
                        Stiff = Stiff + dir(1)*PD%face_mass_x(:,:,f,pp) + &
                                        dir(2)*PD%face_mass_y(:,:,f,pp) + &
                                        dir(3)*PD%face_mass_z(:,:,f,pp)
                end do
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

end module m_transport_precompute_iga
