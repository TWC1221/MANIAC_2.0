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
    use m_quadrature
    use m_basis_fem,   only: GetMapping2D_FEM, GetMapping3D_FEM, &
                             EvalAtFace2D_FEM, EvalAtFace3D_FEM
    use m_material
    use m_sweep_order, only: generate_sweep_order, precompute_reflective_map, &
                             precompute_upwind_indices, all_nodes_in_list
    implicit none
    public :: InitialiseTransport_FEM, InitialiseGeometry_FEM, identify_and_compact_fem

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
    ! FEM: no intra-patch span lookup -- each patch is one element.
    ! Normals/centroids are computed in parallel; connectivity uses a
    ! spatial hash (O(n log n)) instead of the O(n^2) brute-force search.
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

        ! Spatial-hash variables
        integer(kind=8), allocatable :: face_keys(:)
        integer,         allocatable :: face_elem(:), face_face(:), sort_idx(:)
        integer :: n_faces_total, flat_i, grp_start, grp_end, ii, jj
        integer(kind=8) :: kx, ky, kz
        real(dp) :: cmin(3), cmax(3), crange, htol
        integer(kind=8), parameter :: NBINS = 1000003_8

        nf = mesh%n_faces_per_elem

        allocate(TD%face_connectivity(4, nf, mesh%n_elems))
        allocate(TD%face_normals(3, nf, mesh%n_elems))
        allocate(centroids(3, nf, mesh%n_elems))
        TD%face_normals = 0.0_dp; centroids = 0.0_dp

        ! ---- Parallel face normals and centroids ----
        !$OMP PARALLEL DO PRIVATE(ee, nodes, n_pts, f, q, xi_f, eta_f, zeta_f, &
        !$OMP&                    J, J2, dA, N, dN_dx, dN_dy, dN_dz)
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
                        call EvalAtFace3D_FEM(FE, xi_f, eta_f, zeta_f, nodes, N, dN_dx, dN_dy, dN_dz, J)
                        centroids(:,f,ee) = centroids(:,f,ee) + &
                            matmul(transpose(nodes), N) / real(n_pts, dp)
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
                        call EvalAtFace2D_FEM(FE, xi_f, eta_f, nodes(:,1:2), N, J2)
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
        !$OMP END PARALLEL DO

        ! ---- Initialise connectivity to unset/vacuum ----
        TD%face_connectivity(1,:,:) = -1
        TD%face_connectivity(2,:,:) = -1
        TD%face_connectivity(3,:,:) =  0
        TD%face_connectivity(4,:,:) = BC_VACUUM

        ! ---- Spatial-hash O(n log n) face matching ----
        ! Each face centroid is bucketed by rounded physical coordinate.
        ! After sorting by bucket key, only faces in the same bucket are
        ! compared -- O(1) candidates per face for a conforming mesh.
        n_faces_total = nf * mesh%n_elems
        allocate(face_keys(n_faces_total), face_elem(n_faces_total), &
                 face_face(n_faces_total), sort_idx(n_faces_total))

        ! Bounding box → uniform grid cell size (1 ppm of domain extent)
        cmin =  huge(0.0_dp); cmax = -huge(0.0_dp)
        do ee = 1, mesh%n_elems
            do f = 1, nf
                cmin(1) = min(cmin(1), centroids(1,f,ee))
                cmin(2) = min(cmin(2), centroids(2,f,ee))
                cmin(3) = min(cmin(3), centroids(3,f,ee))
                cmax(1) = max(cmax(1), centroids(1,f,ee))
                cmax(2) = max(cmax(2), centroids(2,f,ee))
                cmax(3) = max(cmax(3), centroids(3,f,ee))
            end do
        end do
        crange = max(maxval(cmax - cmin), 1.0e-12_dp)
        htol   = crange / real(NBINS - 1_8, dp)

        flat_i = 0
        do ee = 1, mesh%n_elems
            do f = 1, nf
                flat_i = flat_i + 1
                face_elem(flat_i) = ee
                face_face(flat_i) = f
                sort_idx(flat_i)  = flat_i
                kx = int(nint((centroids(1,f,ee) - cmin(1)) / htol), kind=8)
                ky = int(nint((centroids(2,f,ee) - cmin(2)) / htol), kind=8)
                kz = int(nint((centroids(3,f,ee) - cmin(3)) / htol), kind=8)
                face_keys(flat_i) = kx + NBINS * ky + NBINS * NBINS * kz
            end do
        end do

        call qsort_int64(face_keys, sort_idx, n_faces_total)

        ! Scan sorted list: faces sharing a bucket key are neighbour candidates
        flat_i = 1
        do while (flat_i <= n_faces_total)
            grp_start = flat_i
            do while (flat_i <= n_faces_total)
                if (face_keys(sort_idx(flat_i)) /= face_keys(sort_idx(grp_start))) exit
                flat_i = flat_i + 1
            end do
            grp_end = flat_i - 1

            do ii = grp_start, grp_end - 1
                e1 = face_elem(sort_idx(ii)); f1 = face_face(sort_idx(ii))
                if (TD%face_connectivity(1,f1,e1) /= -1) cycle
                do jj = ii + 1, grp_end
                    e2 = face_elem(sort_idx(jj)); f2 = face_face(sort_idx(jj))
                    if (TD%face_connectivity(1,f2,e2) /= -1) cycle
                    found_neighbor = .false.
                    if (all_nodes_in_list(mesh%elems(e1, FE%face_node_map(:,f1)), &
                                          mesh%elems(e2, FE%face_node_map(:,f2)))) then
                        found_neighbor = .true.
                    else if (dot_product(TD%face_normals(:,f1,e1), &
                                         TD%face_normals(:,f2,e2)) < -0.9_dp) then
                        found_neighbor = .true.
                    end if
                    if (found_neighbor) then
                        orient_val = face_orient_fem(mesh, FE, e1, f1, e2, f2)
                        TD%face_connectivity(1,f1,e1) = e2; TD%face_connectivity(2,f1,e1) = f2
                        TD%face_connectivity(1,f2,e2) = e1; TD%face_connectivity(2,f2,e2) = f1
                        TD%face_connectivity(3,f1,e1) = orient_val
                        TD%face_connectivity(3,f2,e2) = orient_val
                        TD%face_connectivity(4,f1,e1) = 0; TD%face_connectivity(4,f2,e2) = 0
                        exit
                    end if
                end do
            end do
        end do

        deallocate(face_keys, face_elem, face_face, sort_idx, centroids)

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
                        call EvalAtFace3D_FEM(FE, xi_f, eta_f, zeta_f, nodes, N, dN_dx, dN_dy, dN_dz, J)
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
                        call EvalAtFace2D_FEM(FE, xi_f, eta_f, nodes(:,1:2), N, J2)
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

        call identify_and_compact_fem(mesh%n_elems, mesh%n_faces_per_elem, FE%n_basis, TD)
    end subroutine precompute_integrals_fem

    ! ------------------------------------------------------------------
    ! Identify geometrically identical elements and compact shared matrices.
    ! Fingerprint: (element_volume, ||Kx||_F^2, ||Ky||_F^2, ||Kz||_F^2).
    ! Rotated elements are distinct (different stiffness matrices → different fp).
    ! Sets TD%elem_ref_id, TD%n_ref_elems; reallocates all compacted arrays.
    ! ------------------------------------------------------------------
    subroutine identify_and_compact_fem(n_elems, n_faces_per_elem, n_basis, TD)
        integer,        intent(in)    :: n_elems, n_faces_per_elem, n_basis
        type(t_fem_dg), intent(inout) :: TD

        real(dp), parameter :: RTOL = 1.0e-4_dp, ATOL = 1.0e-12_dp

        integer  :: ee, rr, ref_count, nb, nf
        integer,  allocatable :: rep(:)
        real(dp) :: tol_x, tol_y, tol_z, mem_before, mem_after

        real(dp), allocatable :: tmp_mass(:,:,:), tmp_Kx(:,:,:)
        real(dp), allocatable :: tmp_Ky(:,:,:),   tmp_Kz(:,:,:)
        real(dp), allocatable :: tmp_fx(:,:,:,:), tmp_fy(:,:,:,:), tmp_fz(:,:,:,:)
        real(dp), allocatable :: tmp_vol(:,:)

        nb = n_basis
        nf = n_faces_per_elem

        allocate(rep(n_elems))
        allocate(TD%elem_ref_id(n_elems))
        ref_count = 0

        ! Matrix-max-norm comparison: tolerance = RTOL * max(||K||_inf, ATOL).
        ! Per-entry scaling would give ~1e-36 for near-zero entries, which is
        ! tighter than floating-point noise and causes false mismatches at p>1.
        outer: do ee = 1, n_elems
            do rr = 1, ref_count
                tol_x = RTOL * max(maxval(abs(TD%elem_stiffness_x(:,:,rep(rr)))), ATOL)
                tol_y = RTOL * max(maxval(abs(TD%elem_stiffness_y(:,:,rep(rr)))), ATOL)
                tol_z = RTOL * max(maxval(abs(TD%elem_stiffness_z(:,:,rep(rr)))), ATOL)
                if (any(abs(TD%elem_stiffness_x(:,:,ee) - TD%elem_stiffness_x(:,:,rep(rr))) > tol_x)) cycle
                if (any(abs(TD%elem_stiffness_y(:,:,ee) - TD%elem_stiffness_y(:,:,rep(rr))) > tol_y)) cycle
                if (any(abs(TD%elem_stiffness_z(:,:,ee) - TD%elem_stiffness_z(:,:,rep(rr))) > tol_z)) cycle
                TD%elem_ref_id(ee) = rr
                cycle outer
            end do
            ref_count      = ref_count + 1
            rep(ref_count) = ee
            TD%elem_ref_id(ee) = ref_count
        end do outer
        TD%n_ref_elems = ref_count

        mem_before = real((nb**2*4 + nb**2*nf*3 + nb) * n_elems,    dp) * 8.0e-6_dp
        mem_after  = real((nb**2*4 + nb**2*nf*3 + nb) * ref_count,  dp) * 8.0e-6_dp
        write(*,'(A,I0,A,I0,A,F7.2,A,F7.2,A)') &
            "  FEM ref. elems: ", ref_count, " / ", n_elems, &
            "  (geom. matrices: ", mem_before, " -> ", mem_after, " MB)"

        allocate(tmp_mass(nb,nb,ref_count), tmp_Kx(nb,nb,ref_count), &
                 tmp_Ky(nb,nb,ref_count),   tmp_Kz(nb,nb,ref_count), &
                 tmp_fx(nb,nb,nf,ref_count), tmp_fy(nb,nb,nf,ref_count), &
                 tmp_fz(nb,nb,nf,ref_count), tmp_vol(nb,ref_count))

        do rr = 1, ref_count
            ee = rep(rr)
            tmp_mass(:,:,rr) = TD%elem_mass_matrix(:,:,ee)
            tmp_Kx(:,:,rr)   = TD%elem_stiffness_x(:,:,ee)
            tmp_Ky(:,:,rr)   = TD%elem_stiffness_y(:,:,ee)
            tmp_Kz(:,:,rr)   = TD%elem_stiffness_z(:,:,ee)
            tmp_fx(:,:,:,rr) = TD%face_mass_x(:,:,:,ee)
            tmp_fy(:,:,:,rr) = TD%face_mass_y(:,:,:,ee)
            tmp_fz(:,:,:,rr) = TD%face_mass_z(:,:,:,ee)
            tmp_vol(:,rr)    = TD%basis_integrals_vol(:,ee)
        end do

        call move_alloc(tmp_mass, TD%elem_mass_matrix)
        call move_alloc(tmp_Kx,   TD%elem_stiffness_x)
        call move_alloc(tmp_Ky,   TD%elem_stiffness_y)
        call move_alloc(tmp_Kz,   TD%elem_stiffness_z)
        call move_alloc(tmp_fx,   TD%face_mass_x)
        call move_alloc(tmp_fy,   TD%face_mass_y)
        call move_alloc(tmp_fz,   TD%face_mass_z)
        call move_alloc(tmp_vol,  TD%basis_integrals_vol)

    end subroutine identify_and_compact_fem

    ! ------------------------------------------------------------------
    ! Per-LU-class, per-angle, per-group LU factorisation.
    ! LU class = unique (ref_geom_id × material_id) pair.
    ! Face normals come from the representative element of each LU class;
    ! within a reference class all elements have identical normals.
    ! ------------------------------------------------------------------
    subroutine precompute_lu_fem(mesh, FE, sn_quad, materials, n_groups, TD)
        type(t_mesh_fem),      intent(in)    :: mesh
        type(t_basis_fem),     intent(in)    :: FE
        type(t_sn_quadrature), intent(in)    :: sn_quad
        type(t_material),      intent(in)    :: materials(:)
        integer,               intent(in)    :: n_groups
        type(t_fem_dg),        intent(inout) :: TD

        integer  :: ee, mm, g, f, info, n_lu, ll, rr, mid, n_mat_max, ee_rep, ref_id
        integer, allocatable :: lu_class_map(:,:), lu_rep_elem(:)
        real(dp) :: A(FE%n_basis, FE%n_basis), dir(3), o_n
        real(dp) :: StiffOut(FE%n_basis, FE%n_basis)

        do ee = 1, mesh%n_elems
            if (.not. allocated(materials(mesh%material_ids(ee))%SigmaT)) then
                write(*,'(A,I0,A)') "FATAL: Material ID ", mesh%material_ids(ee), &
                    " has no SigmaT. Check mat_file and material IDs."
                stop
            end if
        end do

        n_mat_max = maxval(mesh%material_ids)
        allocate(lu_class_map(TD%n_ref_elems, n_mat_max), source=0)
        allocate(lu_rep_elem(mesh%n_elems))
        allocate(TD%elem_lu_id(mesh%n_elems))

        n_lu = 0
        do ee = 1, mesh%n_elems
            rr  = TD%elem_ref_id(ee)
            mid = mesh%material_ids(ee)
            if (lu_class_map(rr, mid) == 0) then
                n_lu = n_lu + 1
                lu_class_map(rr, mid) = n_lu
                lu_rep_elem(n_lu) = ee
            end if
            TD%elem_lu_id(ee) = lu_class_map(rr, mid)
        end do
        TD%n_lu_classes = n_lu

        write(*,'(A,I0,A,I0,A)') "  FEM LU classes: ", n_lu, " (was ", &
            mesh%n_elems * sn_quad%n_angles * n_groups, " factorizations)"

        allocate(TD%local_lu(FE%n_basis, FE%n_basis, n_lu, sn_quad%n_angles, n_groups), &
                 TD%local_pivots(FE%n_basis, n_lu, sn_quad%n_angles, n_groups))

        !$OMP PARALLEL DO PRIVATE(mm, dir, ll, ee_rep, ref_id, mid, StiffOut, f, o_n, g, A, info)
        do mm = 1, sn_quad%n_angles
            dir = sn_quad%dirs(mm, 1:3)
            do ll = 1, n_lu
                ee_rep = lu_rep_elem(ll)
                ref_id = TD%elem_ref_id(ee_rep)
                mid    = mesh%material_ids(ee_rep)

                StiffOut = -(dir(1)*TD%elem_stiffness_x(:,:,ref_id) + &
                              dir(2)*TD%elem_stiffness_y(:,:,ref_id) + &
                              dir(3)*TD%elem_stiffness_z(:,:,ref_id))
                do f = 1, mesh%n_faces_per_elem
                    o_n = dot_product(dir, TD%face_normals(:,f,ee_rep))
                    if (o_n > 0.0_dp) &
                        StiffOut = StiffOut + (dir(1)*TD%face_mass_x(:,:,f,ref_id) + &
                                               dir(2)*TD%face_mass_y(:,:,f,ref_id) + &
                                               dir(3)*TD%face_mass_z(:,:,f,ref_id))
                end do
                do g = 1, n_groups
                    A = materials(mid)%SigmaT(g) * TD%elem_mass_matrix(:,:,ref_id) + StiffOut
                    call dgetrf(FE%n_basis, FE%n_basis, A, FE%n_basis, &
                                TD%local_pivots(:,ll,mm,g), info)
                    if (info /= 0) then
                        write(*,'(A,2I6)') "FATAL: LU failed for class,angle=", ll, mm; stop
                    end if
                    TD%local_lu(:,:,ll,mm,g) = A
                end do
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine precompute_lu_fem

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

    ! ------------------------------------------------------------------
    ! Indirect quicksort: reorders idx(1:n) so that keys(idx) is
    ! non-decreasing.  Iterative Hoare partition; O(log n) stack depth.
    ! ------------------------------------------------------------------
    subroutine qsort_int64(keys, idx, n)
        integer(kind=8), intent(in)    :: keys(:)
        integer,         intent(inout) :: idx(:)
        integer,         intent(in)    :: n

        integer :: lo, hi, i, j, tmp
        integer(kind=8) :: pv
        integer :: stk(128)   ! 64 levels deep — sufficient for n up to 2^64
        integer :: sp

        if (n <= 1) return
        sp = 0
        sp = sp + 1; stk(sp) = 1
        sp = sp + 1; stk(sp) = n

        do while (sp > 0)
            hi = stk(sp); sp = sp - 1
            lo = stk(sp); sp = sp - 1
            if (lo >= hi) cycle

            ! Hoare partition around median element
            pv = keys(idx((lo + hi) / 2))
            i = lo - 1; j = hi + 1
            do
                do; i = i + 1; if (keys(idx(i)) >= pv) exit; end do
                do; j = j - 1; if (keys(idx(j)) <= pv) exit; end do
                if (i >= j) exit
                tmp = idx(i); idx(i) = idx(j); idx(j) = tmp
            end do

            ! Push larger partition first so smaller is processed next,
            ! bounding stack depth to O(log n).
            if ((j - lo) >= (hi - j - 1)) then
                if (lo  < j  ) then; sp=sp+1; stk(sp)=lo;  sp=sp+1; stk(sp)=j;  end if
                if (j+1 < hi ) then; sp=sp+1; stk(sp)=j+1; sp=sp+1; stk(sp)=hi; end if
            else
                if (j+1 < hi ) then; sp=sp+1; stk(sp)=j+1; sp=sp+1; stk(sp)=hi; end if
                if (lo  < j  ) then; sp=sp+1; stk(sp)=lo;  sp=sp+1; stk(sp)=j;  end if
            end if
        end do
    end subroutine qsort_int64

end module m_transport_precompute_fem
