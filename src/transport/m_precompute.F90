! Transport integral precomputation and LU factorisation for IGA meshes.
! All output goes into t_transport_iga (TD); the mesh itself is read-only.
! Supports both 2D (quad elements, 4 faces) and 3D (hex elements, 6 faces).
!
! Public:
!   InitialiseTransport  -- integrals + reflective angle map + LU factors
module m_transport_precompute
    use m_constants
    use m_types
    use m_types_iga
    use m_quadrature
    use m_basis_iga, only: GetMapping2D, GetMapping3D
    use m_material
    implicit none
    public :: InitialiseTransport

    interface
        subroutine dgetrf(m, n, a, lda, ipiv, info)
            import :: dp
            integer, intent(in)    :: m, n, lda
            real(dp), intent(inout):: a(lda, *)
            integer, intent(out)   :: ipiv(*)
            integer, intent(out)   :: info
        end subroutine dgetrf

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

    subroutine InitialiseTransport(mesh, FE, sn_quad, Quad2D, Quad1D, materials, TD)
        type(t_mesh_iga),      intent(in)    :: mesh
        type(t_finite_iga),    intent(in)    :: FE
        type(t_sn_quadrature), intent(in)    :: sn_quad
        type(t_quadrature),    intent(in)    :: Quad2D, Quad1D
        type(t_material),      intent(in)    :: materials(:)
        type(t_transport_iga), intent(inout) :: TD

        call precompute_integrals(mesh, FE, Quad2D, Quad1D, TD)
        call precompute_reflective_map(mesh, sn_quad, TD)
        call precompute_lu(mesh, FE, sn_quad, materials, mesh%n_groups, TD)
    end subroutine InitialiseTransport

    ! ------------------------------------------------------------------
    ! Volume mass/stiffness matrices and face mass matrices.
    ! Quad    = volume quadrature (3D hex or 2D quad)
    ! QuadFace = face/edge quadrature (2D quad for 3D faces, 1D linear for 2D edges)
    ! ------------------------------------------------------------------
    subroutine precompute_integrals(mesh, FE, Quad, QuadFace, TD)
        type(t_mesh_iga),      intent(in)    :: mesh
        type(t_finite_iga),    intent(in)    :: FE
        type(t_quadrature),    intent(in)    :: Quad, QuadFace
        type(t_transport_iga), intent(inout) :: TD

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

            ! ---- Volume integrals ----
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

            ! ---- Face/edge integrals ----
            do f = 1, nf
                do q = 1, QuadFace%n_points
                    if (mesh%dim == 3) then
                        ! 3D: QuadFace is 2D (xi + eta)
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
                        ! 2D: QuadFace is 1D (xi only)
                        ! Face 1: eta=-1 (bottom),  Face 2: xi=+1 (right)
                        ! Face 3: eta=+1 (top),     Face 4: xi=-1 (left)
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
                            case(1,3)  ! tangent along xi: dA = (dy/dxi, -dx/dxi, 0)
                                s1 = 0.5_dp*(u2-u1)
                                dA(1) =  J2(1,2) * s1
                                dA(2) = -J2(1,1) * s1
                                if (f==3) dA = -dA
                            case(2,4)  ! tangent along eta: dA = (dy/deta, -dx/deta, 0)
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

    subroutine precompute_reflective_map(mesh, sn_quad, TD)
        type(t_mesh_iga),      intent(in)    :: mesh
        type(t_sn_quadrature), intent(in)    :: sn_quad
        type(t_transport_iga), intent(inout) :: TD

        integer :: ee, f, mm, m_iter
        real(dp) :: normal(3), dir(3), ref_dir(3), max_dot, dprod

        allocate(TD%reflect_map(sn_quad%n_angles, mesh%n_faces_per_elem, mesh%n_elems))
        TD%reflect_map = 0

        !$OMP PARALLEL DO PRIVATE(ee, f, normal, mm, dir, ref_dir, max_dot, m_iter, dprod)
        do ee = 1, mesh%n_elems
            do f = 1, mesh%n_faces_per_elem
                normal = TD%face_normals(:,f,ee)
                do mm = 1, sn_quad%n_angles
                    dir = sn_quad%dirs(mm, :)
                    ref_dir = dir - 2.0_dp * dot_product(dir, normal) * normal
                    max_dot = -2.0_dp
                    do m_iter = 1, sn_quad%n_angles
                        if (abs(ref_dir(3) - sn_quad%dirs(m_iter,3)) > SMALL_NUMBER) cycle
                        dprod = dot_product(ref_dir, sn_quad%dirs(m_iter,:))
                        if (dprod > max_dot) then
                            max_dot = dprod
                            TD%reflect_map(mm,f,ee) = m_iter
                        end if
                    end do
                end do
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine precompute_reflective_map

    subroutine precompute_lu(mesh, FE, sn_quad, materials, n_groups, TD)
        type(t_mesh_iga),      intent(in)    :: mesh
        type(t_finite_iga),    intent(in)    :: FE
        type(t_sn_quadrature), intent(in)    :: sn_quad
        type(t_material),      intent(in)    :: materials(:)
        integer,               intent(in)    :: n_groups
        type(t_transport_iga), intent(inout) :: TD

        integer  :: ee, mm, g, f, info
        real(dp) :: A(FE%n_basis, FE%n_basis), dir(3), o_n
        real(dp) :: StiffOut(FE%n_basis, FE%n_basis)

        do ee = 1, mesh%n_elems
            if (.not. allocated(materials(mesh%material_ids(ee))%SigmaT)) then
                write(*,'(A,I0,A)') "FATAL: Material ID ", mesh%material_ids(ee), &
                    " has no SigmaT. Check mat_file and material IDs."
                stop
            end if
        end do

        allocate(TD%local_lu(FE%n_basis, FE%n_basis, mesh%n_elems, sn_quad%n_angles, n_groups), &
                 TD%local_pivots(FE%n_basis, mesh%n_elems, sn_quad%n_angles, n_groups))

        !$OMP PARALLEL DO PRIVATE(mm, dir, ee, StiffOut, f, o_n, g, A, info)
        do mm = 1, sn_quad%n_angles
            dir = sn_quad%dirs(mm, 1:3)
            do ee = 1, mesh%n_elems
                StiffOut = -(dir(1)*TD%elem_stiffness_x(:,:,ee) + &
                              dir(2)*TD%elem_stiffness_y(:,:,ee) + &
                              dir(3)*TD%elem_stiffness_z(:,:,ee))
                do f = 1, mesh%n_faces_per_elem
                    o_n = dot_product(dir, TD%face_normals(:,f,ee))
                    if (o_n > 0.0_dp) &
                        StiffOut = StiffOut + (dir(1)*TD%face_mass_x(:,:,f,ee) + &
                                               dir(2)*TD%face_mass_y(:,:,f,ee) + &
                                               dir(3)*TD%face_mass_z(:,:,f,ee))
                end do

                do g = 1, n_groups
                    A = materials(mesh%material_ids(ee))%SigmaT(g) * TD%elem_mass_matrix(:,:,ee) + StiffOut
                    call dgetrf(FE%n_basis, FE%n_basis, A, FE%n_basis, TD%local_pivots(:,ee,mm,g), info)
                    if (info /= 0) then
                        write(*,'(A,2I6)') "FATAL: LU failed for elem,angle=", ee, mm; stop
                    end if
                    TD%local_lu(:,:,ee,mm,g) = A
                end do
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine precompute_lu

end module m_transport_precompute
