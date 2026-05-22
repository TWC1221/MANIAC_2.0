! Transport integral precomputation, geometry, and LU factorisation for FEM meshes.
! Supports 2D (quad, 4 faces) and 3D (hex, 6 faces), all polynomial orders.
! Lagrange basis arrays are precomputed in t_basis_fem by InitialiseBasisFEM;
! this module assembles the element/face integrals, connectivity, normals,
! and per-angle LU factors into t_fem_dg (TD).
!
! FEM meshes are read via read_asmg_mesh into t_mesh_iga.  Each ASMG patch
! is a single element (trivial knot vectors).  There is no intra-patch span
! connectivity -- all neighbour matching uses the geometric search.
!
! Public:
!   InitialiseTransport_FEM   -- integrals + reflective map + LU factors
!   InitialiseGeometry_FEM    -- connectivity, normals, upwind indices, sweep order
module m_transport_precompute_fem
    use m_constants
    use m_types
    use m_types_fem
    use m_quadrature
    use m_basis_fem,   only: GetMapping2D_FEM, GetMapping3D_FEM, EvalLagrange2D, EvalLagrange3D
    use m_material
    use m_sweep_order, only: generate_sweep_order, precompute_reflective_map, &
                             precompute_upwind_indices, all_nodes_in_list
    implicit none
    public :: InitialiseTransport_FEM, InitialiseGeometry_FEM

    interface
        subroutine dgetrf(m, n, a, lda, ipiv, info)
            import :: dp
            integer,  intent(in)    :: m, n, lda
            real(dp), intent(inout) :: a(lda, *)
            integer,  intent(out)   :: ipiv(*)
            integer,  intent(out)   :: info
        end subroutine dgetrf

        subroutine dgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            character, intent(in)   :: trans
            integer,   intent(in)   :: n, nrhs, lda, ldb
            real(dp),  intent(in)   :: a(lda, *)
            integer,   intent(in)   :: ipiv(*)
            real(dp),  intent(inout):: b(ldb, *)
            integer,   intent(out)  :: info
        end subroutine dgetrs
    end interface

contains

    ! ------------------------------------------------------------------
    ! Master entry: volume/face integrals, reflective angle map, LU factors.
    ! ------------------------------------------------------------------
    subroutine InitialiseTransport_FEM(mesh, FE, sn_quad, QuadVol, QuadFace, materials, TD)
        type(t_mesh_fem),      intent(in)    :: mesh
        type(t_basis_fem),    intent(in)    :: FE
        type(t_sn_quadrature), intent(in)    :: sn_quad
        type(t_quadrature),    intent(in)    :: QuadVol, QuadFace
        type(t_material),      intent(in)    :: materials(:)
        type(t_fem_dg), intent(inout) :: TD

        call precompute_integrals_fem(mesh, FE, QuadVol, QuadFace, TD)
        call precompute_reflective_map(mesh%n_elems, mesh%n_faces_per_elem, &
                                       sn_quad, TD%face_normals, TD%reflect_map)
        call precompute_lu_fem(mesh, FE, sn_quad, materials, mesh%n_groups, TD)
    end subroutine InitialiseTransport_FEM

    ! ------------------------------------------------------------------
    ! Master geometry entry: connectivity, normals, upwind indices, sweep order.
    ! ------------------------------------------------------------------
    subroutine InitialiseGeometry_FEM(mesh, FE, QuadSn, QuadFace, TD, sweep_order)
        type(t_mesh_fem),      intent(in)    :: mesh
        type(t_basis_fem),    intent(in)    :: FE
        type(t_sn_quadrature), intent(in)    :: QuadSn
        type(t_quadrature),    intent(in)    :: QuadFace
        type(t_fem_dg), intent(inout) :: TD
        integer, allocatable,  intent(out)   :: sweep_order(:,:)

        integer :: mm
        real(dp) :: dir_tmp(3)

        call connectivity_and_normals_fem(mesh, FE, QuadFace, TD)
        call precompute_upwind_indices(mesh%n_elems, mesh%n_faces_per_elem, &
                                       FE%n_nodes_per_face, FE%n_basis, &
                                       FE%face_node_map, TD%face_connectivity, TD%upwind_idx)

        allocate(sweep_order(mesh%n_elems, QuadSn%n_angles))
        do mm = 1, QuadSn%n_angles
            dir_tmp = QuadSn%dirs(mm, :)
            call generate_sweep_order(mesh%n_elems, mesh%n_faces_per_elem, &
                                      TD%face_normals, TD%face_connectivity, &
                                      dir_tmp, sweep_order(:, mm))
        end do
    end subroutine InitialiseGeometry_FEM

    ! ------------------------------------------------------------------
    ! Face connectivity, outward normals, and boundary BC assignment.
    ! FEM: no intra-patch span lookup -- each patch is one element,
    ! so all connectivity is found via geometric face matching.
    ! ------------------------------------------------------------------
    subroutine connectivity_and_normals_fem(mesh, FE, QuadFace, TD)
        type(t_mesh_fem),      intent(in)    :: mesh
        type(t_basis_fem),    intent(in)    :: FE
        type(t_quadrature),    intent(in)    :: QuadFace
        type(t_fem_dg), intent(inout) :: TD

        integer  :: ee, e1, e2, f, f1, f2, q, s_idx, n_pts, orient_val
        real(dp) :: nodes(FE%n_basis, 3)
        real(dp) :: dN_dx(FE%n_basis), dN_dy(FE%n_basis), dN_dz(FE%n_basis)
        real(dp) :: N(FE%n_basis)
        real(dp) :: xi_f, eta_f, zeta_f
        real(dp) :: J(3,3), J2(2,2), dA(3)
        real(dp), allocatable :: centroids(:,:,:)
        integer  :: nf
        logical  :: found_neighbor

        nf = mesh%n_faces_per_elem

        allocate(TD%face_connectivity(4, nf, mesh%n_elems))
        allocate(TD%face_normals(3, nf, mesh%n_elems))
        allocate(centroids(3, nf, mesh%n_elems))
        TD%face_normals = 0.0_dp; centroids = 0.0_dp

        ! ---- Face normals and centroids ----
        do ee = 1, mesh%n_elems
            nodes = mesh%nodes(mesh%elems(ee, 1:FE%n_basis), :)
            n_pts = QuadFace%n_points

            if (mesh%dim == 3) then
                do f = 1, 6
                    dA = 0.0_dp
                    do q = 1, n_pts
                        select case (f)
                            case(1); xi_f=QuadFace%xi(q); eta_f=QuadFace%eta(q); zeta_f=-1.0_dp
                            case(2); xi_f=QuadFace%xi(q); eta_f=QuadFace%eta(q); zeta_f= 1.0_dp
                            case(3); xi_f=QuadFace%xi(q); eta_f=-1.0_dp;         zeta_f=QuadFace%eta(q)
                            case(4); xi_f=QuadFace%xi(q); eta_f= 1.0_dp;         zeta_f=QuadFace%eta(q)
                            case(5); xi_f=-1.0_dp;        eta_f=QuadFace%xi(q);  zeta_f=QuadFace%eta(q)
                            case(6); xi_f= 1.0_dp;        eta_f=QuadFace%xi(q);  zeta_f=QuadFace%eta(q)
                        end select
                        call eval_fem_3d_at(FE, xi_f, eta_f, zeta_f, nodes, N, dN_dx, dN_dy, dN_dz, J)
                        centroids(:,f,ee) = centroids(:,f,ee) + &
                            matmul(transpose(nodes), N) / real(n_pts, dp)
                        ! Area element from surface tangent cross product
                        select case (f)
                            case(1,2)
                                dA(1)=dA(1)+(J(1,2)*J(2,3)-J(1,3)*J(2,2))*QuadFace%weights(q)
                                dA(2)=dA(2)+(J(1,3)*J(2,1)-J(1,1)*J(2,3))*QuadFace%weights(q)
                                dA(3)=dA(3)+(J(1,1)*J(2,2)-J(1,2)*J(2,1))*QuadFace%weights(q)
                            case(3,4)
                                dA(1)=dA(1)+(J(3,2)*J(1,3)-J(3,3)*J(1,2))*QuadFace%weights(q)
                                dA(2)=dA(2)+(J(3,3)*J(1,1)-J(3,1)*J(1,3))*QuadFace%weights(q)
                                dA(3)=dA(3)+(J(3,1)*J(1,2)-J(3,2)*J(1,1))*QuadFace%weights(q)
                            case(5,6)
                                dA(1)=dA(1)+(J(2,2)*J(3,3)-J(2,3)*J(3,2))*QuadFace%weights(q)
                                dA(2)=dA(2)+(J(2,3)*J(3,1)-J(2,1)*J(3,3))*QuadFace%weights(q)
                                dA(3)=dA(3)+(J(2,1)*J(3,2)-J(2,2)*J(3,1))*QuadFace%weights(q)
                        end select
                    end do
                    if (f==1 .or. f==3 .or. f==5) dA = -dA
                    if (norm2(dA) > dp_EPSILON) TD%face_normals(:,f,ee) = dA / norm2(dA)
                end do

            else  ! dim == 2
                do f = 1, 4
                    dA = 0.0_dp
                    do q = 1, n_pts
                        select case (f)
                            case(1); xi_f=QuadFace%xi(q); eta_f=-1.0_dp
                            case(2); xi_f= 1.0_dp;        eta_f=QuadFace%xi(q)
                            case(3); xi_f=QuadFace%xi(q); eta_f= 1.0_dp
                            case(4); xi_f=-1.0_dp;        eta_f=QuadFace%xi(q)
                        end select
                        call eval_fem_2d_at(FE, xi_f, eta_f, nodes(:,1:2), N, J2)
                        centroids(1,f,ee) = centroids(1,f,ee) + dot_product(N, nodes(:,1)) / real(n_pts,dp)
                        centroids(2,f,ee) = centroids(2,f,ee) + dot_product(N, nodes(:,2)) / real(n_pts,dp)
                        select case (f)
                            case(1,3)
                                dA(1) = dA(1) + ( J2(1,2))*QuadFace%weights(q)
                                dA(2) = dA(2) + (-J2(1,1))*QuadFace%weights(q)
                            case(2,4)
                                dA(1) = dA(1) + ( J2(2,2))*QuadFace%weights(q)
                                dA(2) = dA(2) + (-J2(2,1))*QuadFace%weights(q)
                        end select
                    end do
                    dA(3) = 0.0_dp
                    if (f==3 .or. f==4) dA = -dA
                    if (norm2(dA) > dp_EPSILON) TD%face_normals(:,f,ee) = dA / norm2(dA)
                end do
            end if
        end do

        ! ---- Initialise connectivity to unset/vacuum ----
        TD%face_connectivity(1,:,:) = -1
        TD%face_connectivity(2,:,:) = -1
        TD%face_connectivity(3,:,:) =  0
        TD%face_connectivity(4,:,:) = BC_VACUUM

        ! ---- Geometric face matching (all connectivity for FEM) ----
        do e1 = 1, mesh%n_elems
            if (.not. any(TD%face_connectivity(1,:,e1) == -1)) cycle
            do f1 = 1, nf
                if (TD%face_connectivity(1,f1,e1) /= -1) cycle
                do e2 = e1+1, mesh%n_elems
                    if (.not. any(TD%face_connectivity(1,:,e2) == -1)) cycle
                    do f2 = 1, nf
                        if (TD%face_connectivity(1,f2,e2) /= -1) cycle
                        found_neighbor = .false.
                        if (all_nodes_in_list(mesh%elems(e1, FE%face_node_map(:,f1)), &
                                              mesh%elems(e2, FE%face_node_map(:,f2)))) then
                            found_neighbor = .true.
                        else if (norm2(centroids(:,f1,e1) - centroids(:,f2,e2)) < 1.0e-4_dp) then
                            if (dot_product(TD%face_normals(:,f1,e1), &
                                            TD%face_normals(:,f2,e2)) < -0.9_dp) found_neighbor = .true.
                        end if
                        if (found_neighbor) then
                            orient_val = face_orient_fem(mesh, FE, e1, f1, e2, f2)
                            TD%face_connectivity(1,f1,e1)=e2; TD%face_connectivity(2,f1,e1)=f2
                            TD%face_connectivity(1,f2,e2)=e1; TD%face_connectivity(2,f2,e2)=f1
                            TD%face_connectivity(3,f1,e1)=orient_val
                            TD%face_connectivity(3,f2,e2)=orient_val
                            TD%face_connectivity(4,f1,e1)=0; TD%face_connectivity(4,f2,e2)=0
                            exit
                        end if
                    end do
                    if (TD%face_connectivity(1,f1,e1) /= -1) exit
                end do
            end do
        end do
        deallocate(centroids)

        ! ---- BC assignment from surfaces ----
        do e1 = 1, mesh%n_elems
            do f1 = 1, nf
                if (TD%face_connectivity(1,f1,e1) == -1) then
                    found_neighbor = .false.
                    do s_idx = 1, size(mesh%surfaces)
                        if (all_nodes_in_list(mesh%elems(e1, FE%face_node_map(:,f1)), &
                                              mesh%surfaces(s_idx)%cp_ids)) then
                            TD%face_connectivity(4,f1,e1) = mesh%surfaces(s_idx)%bc_id
                            found_neighbor = .true.
                            exit
                        end if
                    end do
                    if (.not. found_neighbor) TD%face_connectivity(4,f1,e1) = BC_VACUUM
                end if
            end do
        end do
    end subroutine connectivity_and_normals_fem

    ! ------------------------------------------------------------------
    ! Volume mass/stiffness and face mass matrices.
    ! FEM: basis precomputed in FE arrays; no knot-span Jacobian scaling.
    ! ------------------------------------------------------------------
    subroutine precompute_integrals_fem(mesh, FE, QuadVol, QuadFace, TD)
        type(t_mesh_fem),      intent(in)    :: mesh
        type(t_basis_fem),    intent(in)    :: FE
        type(t_quadrature),    intent(in)    :: QuadVol, QuadFace
        type(t_fem_dg), intent(inout) :: TD

        integer  :: ee, q, f, nf
        real(dp) :: nodes(FE%n_basis, 3)
        real(dp) :: dN_dx(FE%n_basis), dN_dy(FE%n_basis), dN_dz(FE%n_basis)
        real(dp) :: detJ, dV, N(FE%n_basis)
        real(dp) :: xi_f, eta_f, zeta_f, J(3,3), J2(2,2), dA(3)

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

        !$OMP PARALLEL DO PRIVATE(ee, nodes, q, dN_dx, dN_dy, dN_dz, detJ, N, dV, &
        !$OMP&   f, xi_f, eta_f, zeta_f, J, J2, dA)
        do ee = 1, mesh%n_elems
            nodes = mesh%nodes(mesh%elems(ee, 1:FE%n_basis), :)

            ! ---- Volume integrals ----
            do q = 1, QuadVol%n_points
                if (mesh%dim == 3) then
                    call GetMapping3D_FEM(FE, q, nodes, dN_dx, dN_dy, dN_dz, detJ, N)
                else
                    call GetMapping2D_FEM(FE, q, nodes(:,1:2), dN_dx, dN_dy, detJ, N)
                    dN_dz = 0.0_dp
                end if
                dV = detJ * QuadVol%weights(q)
                TD%elem_mass_matrix(:,:,ee) = TD%elem_mass_matrix(:,:,ee) + &
                    spread(N,2,FE%n_basis) * spread(N,1,FE%n_basis) * dV
                TD%elem_stiffness_x(:,:,ee) = TD%elem_stiffness_x(:,:,ee) + &
                    spread(dN_dx,2,FE%n_basis) * spread(N,1,FE%n_basis) * dV
                TD%elem_stiffness_y(:,:,ee) = TD%elem_stiffness_y(:,:,ee) + &
                    spread(dN_dy,2,FE%n_basis) * spread(N,1,FE%n_basis) * dV
                TD%elem_stiffness_z(:,:,ee) = TD%elem_stiffness_z(:,:,ee) + &
                    spread(dN_dz,2,FE%n_basis) * spread(N,1,FE%n_basis) * dV
                TD%basis_integrals_vol(:,ee) = TD%basis_integrals_vol(:,ee) + N * dV
            end do

            ! ---- Face integrals ----
            do f = 1, nf
                do q = 1, QuadFace%n_points
                    if (mesh%dim == 3) then
                        select case (f)
                            case(1); xi_f=QuadFace%xi(q); eta_f=QuadFace%eta(q); zeta_f=-1.0_dp
                            case(2); xi_f=QuadFace%xi(q); eta_f=QuadFace%eta(q); zeta_f= 1.0_dp
                            case(3); xi_f=QuadFace%xi(q); eta_f=-1.0_dp;         zeta_f=QuadFace%eta(q)
                            case(4); xi_f=QuadFace%xi(q); eta_f= 1.0_dp;         zeta_f=QuadFace%eta(q)
                            case(5); xi_f=-1.0_dp;        eta_f=QuadFace%xi(q);  zeta_f=QuadFace%eta(q)
                            case(6); xi_f= 1.0_dp;        eta_f=QuadFace%xi(q);  zeta_f=QuadFace%eta(q)
                        end select
                        call eval_fem_3d_at(FE, xi_f, eta_f, zeta_f, nodes, N, dN_dx, dN_dy, dN_dz, J)
                        select case (f)
                            case(1,2)
                                dA(1)=(J(1,2)*J(2,3)-J(1,3)*J(2,2))
                                dA(2)=(J(1,3)*J(2,1)-J(1,1)*J(2,3))
                                dA(3)=(J(1,1)*J(2,2)-J(1,2)*J(2,1))
                                if (f==1) dA=-dA
                            case(3,4)
                                dA(1)=(J(3,2)*J(1,3)-J(3,3)*J(1,2))
                                dA(2)=(J(3,3)*J(1,1)-J(3,1)*J(1,3))
                                dA(3)=(J(3,1)*J(1,2)-J(3,2)*J(1,1))
                                if (f==3) dA=-dA
                            case(5,6)
                                dA(1)=(J(2,2)*J(3,3)-J(2,3)*J(3,2))
                                dA(2)=(J(2,3)*J(3,1)-J(2,1)*J(3,3))
                                dA(3)=(J(2,1)*J(3,2)-J(2,2)*J(3,1))
                                if (f==5) dA=-dA
                        end select
                    else
                        select case (f)
                            case(1); xi_f=QuadFace%xi(q); eta_f=-1.0_dp
                            case(2); xi_f= 1.0_dp;        eta_f=QuadFace%xi(q)
                            case(3); xi_f=QuadFace%xi(q); eta_f= 1.0_dp
                            case(4); xi_f=-1.0_dp;        eta_f=QuadFace%xi(q)
                        end select
                        call eval_fem_2d_at(FE, xi_f, eta_f, nodes(:,1:2), N, J2)
                        dA(3) = 0.0_dp
                        select case (f)
                            case(1,3)
                                dA(1) =  J2(1,2)
                                dA(2) = -J2(1,1)
                                if (f==3) dA = -dA
                            case(2,4)
                                dA(1) =  J2(2,2)
                                dA(2) = -J2(2,1)
                                if (f==4) dA = -dA
                        end select
                    end if

                    dA = dA * QuadFace%weights(q)
                    TD%face_mass_x(:,:,f,ee) = TD%face_mass_x(:,:,f,ee) + &
                        spread(N,2,FE%n_basis)*spread(N,1,FE%n_basis)*dA(1)
                    TD%face_mass_y(:,:,f,ee) = TD%face_mass_y(:,:,f,ee) + &
                        spread(N,2,FE%n_basis)*spread(N,1,FE%n_basis)*dA(2)
                    TD%face_mass_z(:,:,f,ee) = TD%face_mass_z(:,:,f,ee) + &
                        spread(N,2,FE%n_basis)*spread(N,1,FE%n_basis)*dA(3)
                end do
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine precompute_integrals_fem

    ! ------------------------------------------------------------------
    ! Per-element, per-angle, per-group LU factorisation.
    ! ------------------------------------------------------------------
    subroutine precompute_lu_fem(mesh, FE, sn_quad, materials, n_groups, TD)
        type(t_mesh_fem),      intent(in)    :: mesh
        type(t_basis_fem),    intent(in)    :: FE
        type(t_sn_quadrature), intent(in)    :: sn_quad
        type(t_material),      intent(in)    :: materials(:)
        integer,               intent(in)    :: n_groups
        type(t_fem_dg), intent(inout) :: TD

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
                    A = materials(mesh%material_ids(ee))%SigmaT(g) * &
                        TD%elem_mass_matrix(:,:,ee) + StiffOut
 
                    call dgetrf(FE%n_basis, FE%n_basis, A, FE%n_basis, &
                                TD%local_pivots(:,ee,mm,g), info)
                    if (info /= 0) then
                        write(*,'(A,2I6)') "FATAL: LU failed for elem,angle=", ee, mm; stop
                    end if
                    TD%local_lu(:,:,ee,mm,g) = A
                end do
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine precompute_lu_fem

    ! ------------------------------------------------------------------
    ! Evaluate Lagrange basis + Jacobian at an arbitrary reference point.
    ! Used for face normal and face integral computations.
    ! ------------------------------------------------------------------
    subroutine eval_fem_2d_at(FE, xi, eta, coords2d, N, J)
        type(t_basis_fem), intent(in)  :: FE
        real(dp),           intent(in)  :: xi, eta, coords2d(:,:)
        real(dp),           intent(out) :: N(:), J(2,2)
        real(dp) :: dN_dxi(FE%n_basis), dN_deta(FE%n_basis)
        call EvalLagrange2D(FE, xi, eta, N, dN_dxi, dN_deta)
        J(1,:) = matmul(dN_dxi,  coords2d)
        J(2,:) = matmul(dN_deta, coords2d)
    end subroutine eval_fem_2d_at

    subroutine eval_fem_3d_at(FE, xi, eta, zeta, coords3d, N, dN_dx, dN_dy, dN_dz, J)
        type(t_basis_fem), intent(in)  :: FE
        real(dp),           intent(in)  :: xi, eta, zeta, coords3d(:,:)
        real(dp),           intent(out) :: N(:), dN_dx(:), dN_dy(:), dN_dz(:), J(3,3)
        real(dp) :: dN_dxi(FE%n_basis), dN_deta(FE%n_basis), dN_dzeta(FE%n_basis)
        real(dp) :: invJ(3,3), detJ
        call EvalLagrange3D(FE, xi, eta, zeta, N, dN_dxi, dN_deta, dN_dzeta)
        J(1,:) = matmul(dN_dxi,   coords3d)
        J(2,:) = matmul(dN_deta,  coords3d)
        J(3,:) = matmul(dN_dzeta, coords3d)
        detJ = J(1,1)*(J(2,2)*J(3,3)-J(2,3)*J(3,2)) &
             - J(1,2)*(J(2,1)*J(3,3)-J(2,3)*J(3,1)) &
             + J(1,3)*(J(2,1)*J(3,2)-J(2,2)*J(3,1))
        if (abs(detJ) < dp_EPSILON) detJ = sign(dp_EPSILON, detJ)
        invJ(1,1)=(J(2,2)*J(3,3)-J(2,3)*J(3,2))/detJ; invJ(1,2)=(J(1,3)*J(3,2)-J(1,2)*J(3,3))/detJ
        invJ(1,3)=(J(1,2)*J(2,3)-J(1,3)*J(2,2))/detJ; invJ(2,1)=(J(2,3)*J(3,1)-J(2,1)*J(3,3))/detJ
        invJ(2,2)=(J(1,1)*J(3,3)-J(1,3)*J(3,1))/detJ; invJ(2,3)=(J(1,3)*J(2,1)-J(1,1)*J(2,3))/detJ
        invJ(3,1)=(J(2,1)*J(3,2)-J(2,2)*J(3,1))/detJ; invJ(3,2)=(J(1,2)*J(3,1)-J(1,1)*J(3,2))/detJ
        invJ(3,3)=(J(1,1)*J(2,2)-J(1,2)*J(2,1))/detJ
        dN_dx = invJ(1,1)*dN_dxi + invJ(1,2)*dN_deta + invJ(1,3)*dN_dzeta
        dN_dy = invJ(2,1)*dN_dxi + invJ(2,2)*dN_deta + invJ(2,3)*dN_dzeta
        dN_dz = invJ(3,1)*dN_dxi + invJ(3,2)*dN_deta + invJ(3,3)*dN_dzeta
    end subroutine eval_fem_3d_at

    ! ------------------------------------------------------------------
    ! Face orientation: 1 if shared face nodes start with the same node,
    ! -1 otherwise (reversed). Used for upwind DOF alignment.
    ! ------------------------------------------------------------------
    integer function face_orient_fem(mesh, FE, e1, f1, e2, f2) result(orient)
        type(t_mesh_fem),   intent(in) :: mesh
        type(t_basis_fem), intent(in) :: FE
        integer,            intent(in) :: e1, f1, e2, f2
        integer :: first_node_e1
        first_node_e1 = mesh%elems(e1, FE%face_node_map(1, f1))
        orient = merge(1, -1, mesh%elems(e2, FE%face_node_map(1, f2)) == first_node_e1)
    end function face_orient_fem

    subroutine print_matrix(M, n)
        real(dp), intent(in) :: M(:,:)
        integer,  intent(in) :: n
        integer :: i
        do i = 1, n
            write(*,'(*(F10.5,1X))') M(i, 1:n)
        end do
    end subroutine print_matrix

end module m_transport_precompute_fem
