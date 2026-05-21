! Transport sweep-order and geometry initialisation for IGA meshes.
! Builds element-to-element connectivity, outward face normals, upwind
! index tables, and angle-dependent sweep orderings.
! All output populates t_transport_data (TD), not the mesh itself.
! Supports 2D (quad, 4 faces) and 3D (hex, 6 faces).
!
! Public:
!   InitialiseGeometry        -- master entry: connectivity + normals + sweep order
!   connectivity_and_normals  -- also public for unit testing
module m_sweep_order
    use m_constants
    use m_types
    use m_types_iga
    use m_quadrature
    use m_basis_iga, only: GetMapping2D, GetMapping3D
    implicit none
    private
    public :: InitialiseGeometry, connectivity_and_normals
    public :: generate_sweep_order, precompute_reflective_map, precompute_upwind_indices
    public :: all_nodes_in_list

contains

    subroutine InitialiseGeometry(mesh, FE, QuadSn, QuadFace, TD, sweep_order)
        type(t_mesh_iga),      intent(in)    :: mesh
        type(t_finite_iga),    intent(in)    :: FE
        type(t_sn_quadrature), intent(in)    :: QuadSn
        type(t_quadrature),    intent(in)    :: QuadFace
        type(t_transport_data), intent(inout) :: TD
        integer, allocatable,  intent(out)   :: sweep_order(:,:)

        integer :: mm
        real(dp) :: dir_tmp(3)

        call connectivity_and_normals(mesh, FE, QuadFace, TD)
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
    end subroutine InitialiseGeometry

    subroutine connectivity_and_normals(mesh, FE, QuadFace, TD)
        type(t_mesh_iga),      intent(in)    :: mesh
        type(t_finite_iga),    intent(in)    :: FE
        type(t_quadrature),    intent(in)    :: QuadFace
        type(t_transport_data), intent(inout) :: TD

        integer  :: ee, e1, e2, f, f1, f2, q, s_idx, u, v, w, nid, k_iter, n_pts, p, orient_val
        real(dp) :: nodes(FE%n_basis, 3)
        real(dp) :: dN_dx(FE%n_basis), dN_dy(FE%n_basis), dN_dz(FE%n_basis)
        real(dp) :: detJ, R(FE%n_basis)
        real(dp) :: u1, u2, v1, v2, w1, w2, xi_f, eta_f, zeta_f
        real(dp) :: J(3,3), J2(2,2), dA(3), s1, s2, pos(3)
        real(dp), allocatable :: centroids(:,:,:)
        integer  :: nf
        logical  :: found_neighbor

        nf = mesh%n_faces_per_elem

        allocate(TD%face_connectivity(4, nf, mesh%n_elems))
        allocate(TD%face_normals(3, nf, mesh%n_elems))
        allocate(centroids(3, nf, mesh%n_elems))
        TD%face_normals = 0.0_dp; centroids = 0.0_dp

        do ee = 1, mesh%n_elems
            nodes = mesh%nodes(mesh%elems(ee, 1:FE%n_basis), :)
            u1 = mesh%elem_u_min(ee); u2 = mesh%elem_u_max(ee)
            v1 = mesh%elem_v_min(ee); v2 = mesh%elem_v_max(ee)

            if (mesh%dim == 3) then
                w1 = mesh%elem_w_min(ee); w2 = mesh%elem_w_max(ee)

                do f = 1, 6
                    dA = 0.0_dp
                    n_pts = QuadFace%n_points
                    do q = 1, n_pts
                        select case(f)
                            case(1); xi_f=QuadFace%xi(q); eta_f=QuadFace%eta(q); zeta_f=-1.0_dp
                            case(2); xi_f=QuadFace%xi(q); eta_f=QuadFace%eta(q); zeta_f= 1.0_dp
                            case(3); xi_f=QuadFace%xi(q); eta_f=-1.0_dp;         zeta_f=QuadFace%eta(q)
                            case(4); xi_f=QuadFace%xi(q); eta_f= 1.0_dp;         zeta_f=QuadFace%eta(q)
                            case(5); xi_f=-1.0_dp;        eta_f=QuadFace%xi(q);  zeta_f=QuadFace%eta(q)
                            case(6); xi_f= 1.0_dp;        eta_f=QuadFace%xi(q);  zeta_f=QuadFace%eta(q)
                        end select
                        call GetMapping3D(FE, ee, mesh, q, QuadFace, u1, u2, v1, v2, w1, w2, nodes, &
                                          dN_dx, dN_dy, dN_dz, detJ, R, &
                                          xi_custom=xi_f, eta_custom=eta_f, zeta_custom=zeta_f, &
                                          J_out=J, R_out=pos)
                        centroids(:,f,ee) = centroids(:,f,ee) + pos / real(n_pts, dp)
                        select case(f)
                            case(1,2)
                                s1=0.5_dp*(u2-u1); s2=0.5_dp*(v2-v1)
                                dA(1)=dA(1)+(J(1,2)*J(2,3)-J(1,3)*J(2,2))*s1*s2*QuadFace%weights(q)
                                dA(2)=dA(2)+(J(1,3)*J(2,1)-J(1,1)*J(2,3))*s1*s2*QuadFace%weights(q)
                                dA(3)=dA(3)+(J(1,1)*J(2,2)-J(1,2)*J(2,1))*s1*s2*QuadFace%weights(q)
                            case(3,4)
                                s1=0.5_dp*(w2-w1); s2=0.5_dp*(u2-u1)
                                dA(1)=dA(1)+(J(3,2)*J(1,3)-J(3,3)*J(1,2))*s1*s2*QuadFace%weights(q)
                                dA(2)=dA(2)+(J(3,3)*J(1,1)-J(3,1)*J(1,3))*s1*s2*QuadFace%weights(q)
                                dA(3)=dA(3)+(J(3,1)*J(1,2)-J(3,2)*J(1,1))*s1*s2*QuadFace%weights(q)
                            case(5,6)
                                s1=0.5_dp*(v2-v1); s2=0.5_dp*(w2-w1)
                                dA(1)=dA(1)+(J(2,2)*J(3,3)-J(2,3)*J(3,2))*s1*s2*QuadFace%weights(q)
                                dA(2)=dA(2)+(J(2,3)*J(3,1)-J(2,1)*J(3,3))*s1*s2*QuadFace%weights(q)
                                dA(3)=dA(3)+(J(2,1)*J(3,2)-J(2,2)*J(3,1))*s1*s2*QuadFace%weights(q)
                        end select
                    end do
                    if (f==1 .or. f==3 .or. f==5) dA = -dA
                    if (norm2(dA) > dp_EPSILON) TD%face_normals(:,f,ee) = dA / norm2(dA)
                end do

            else  ! dim == 2
                ! QuadFace is 1D (xi only); face numbering:
                ! 1=bottom(eta=-1), 2=right(xi=+1), 3=top(eta=+1), 4=left(xi=-1)
                do f = 1, 4
                    dA = 0.0_dp
                    n_pts = QuadFace%n_points
                    do q = 1, n_pts
                        select case(f)
                            case(1); xi_f=QuadFace%xi(q); eta_f=-1.0_dp
                            case(2); xi_f= 1.0_dp;        eta_f=QuadFace%xi(q)
                            case(3); xi_f=QuadFace%xi(q); eta_f= 1.0_dp
                            case(4); xi_f=-1.0_dp;        eta_f=QuadFace%xi(q)
                        end select
                        call GetMapping2D(FE, ee, mesh, 1, QuadFace, u1, u2, v1, v2, nodes(:,1:2), &
                                          dN_dx, dN_dy, detJ, R, &
                                          xi_custom=xi_f, eta_custom=eta_f, J_out=J2)
                        centroids(1,f,ee) = centroids(1,f,ee) + dot_product(R, nodes(:,1)) / real(n_pts,dp)
                        centroids(2,f,ee) = centroids(2,f,ee) + dot_product(R, nodes(:,2)) / real(n_pts,dp)
                        select case(f)
                            case(1,3)
                                s1 = 0.5_dp*(u2-u1)
                                dA(1) = dA(1) + ( J2(1,2)*s1)*QuadFace%weights(q)
                                dA(2) = dA(2) + (-J2(1,1)*s1)*QuadFace%weights(q)
                            case(2,4)
                                s1 = 0.5_dp*(v2-v1)
                                dA(1) = dA(1) + ( J2(2,2)*s1)*QuadFace%weights(q)
                                dA(2) = dA(2) + (-J2(2,1)*s1)*QuadFace%weights(q)
                        end select
                    end do
                    dA(3) = 0.0_dp
                    if (f==3 .or. f==4) dA = -dA
                    if (norm2(dA) > dp_EPSILON) TD%face_normals(:,f,ee) = dA / norm2(dA)
                end do
            end if
        end do

        TD%face_connectivity(1,:,:) = -1
        TD%face_connectivity(2,:,:) = -1
        TD%face_connectivity(3,:,:) =  0
        TD%face_connectivity(4,:,:) = BC_VACUUM

        ! Intra-patch connectivity via parametric span indexing
        do ee = 1, mesh%n_elems
            p = mesh%elem_patch_id(ee)
            u = mesh%elem_span_indices(1, ee)
            v = mesh%elem_span_indices(2, ee)

            if (mesh%dim == 2) then
                ! Face 1 (eta=-1): lower v neighbour's face 3
                do k_iter = v-1, 1, -1
                    nid = mesh%elem_map_2d(p, u, k_iter)
                    if (nid > 0) then
                        TD%face_connectivity(1,1,ee)=nid; TD%face_connectivity(2,1,ee)=3
                        TD%face_connectivity(3,1,ee)=-1;  TD%face_connectivity(4,1,ee)=0
                        TD%face_connectivity(1,3,nid)=ee; TD%face_connectivity(2,3,nid)=1
                        TD%face_connectivity(3,3,nid)=-1; TD%face_connectivity(4,3,nid)=0
                        exit
                    end if
                end do
                ! Face 4 (xi=-1): lower u neighbour's face 2
                do k_iter = u-1, 1, -1
                    nid = mesh%elem_map_2d(p, k_iter, v)
                    if (nid > 0) then
                        TD%face_connectivity(1,4,ee)=nid; TD%face_connectivity(2,4,ee)=2
                        TD%face_connectivity(3,4,ee)=-1;  TD%face_connectivity(4,4,ee)=0
                        TD%face_connectivity(1,2,nid)=ee; TD%face_connectivity(2,2,nid)=4
                        TD%face_connectivity(3,2,nid)=-1; TD%face_connectivity(4,2,nid)=0
                        exit
                    end if
                end do

            else  ! dim == 3
                w = mesh%elem_span_indices(3, ee)
                ! Face 1 (z=-1): lower w
                do k_iter = w-1, 1, -1
                    nid = mesh%elem_map_3d(p, u, v, k_iter)
                    if (nid > 0) then
                        TD%face_connectivity(1,1,ee)=nid; TD%face_connectivity(2,1,ee)=2
                        TD%face_connectivity(3,1,ee)=1;   TD%face_connectivity(4,1,ee)=0
                        TD%face_connectivity(1,2,nid)=ee; TD%face_connectivity(2,2,nid)=1
                        TD%face_connectivity(3,2,nid)=1;  TD%face_connectivity(4,2,nid)=0
                        exit
                    end if
                end do
                ! Face 3 (y=-1): lower v
                do k_iter = v-1, 1, -1
                    nid = mesh%elem_map_3d(p, u, k_iter, w)
                    if (nid > 0) then
                        TD%face_connectivity(1,3,ee)=nid; TD%face_connectivity(2,3,ee)=4
                        TD%face_connectivity(3,3,ee)=1;   TD%face_connectivity(4,3,ee)=0
                        TD%face_connectivity(1,4,nid)=ee; TD%face_connectivity(2,4,nid)=3
                        TD%face_connectivity(3,4,nid)=1;  TD%face_connectivity(4,4,nid)=0
                        exit
                    end if
                end do
                ! Face 5 (x=-1): lower u
                do k_iter = u-1, 1, -1
                    nid = mesh%elem_map_3d(p, k_iter, v, w)
                    if (nid > 0) then
                        TD%face_connectivity(1,5,ee)=nid; TD%face_connectivity(2,5,ee)=6
                        TD%face_connectivity(3,5,ee)=1;   TD%face_connectivity(4,5,ee)=0
                        TD%face_connectivity(1,6,nid)=ee; TD%face_connectivity(2,6,nid)=5
                        TD%face_connectivity(3,6,nid)=1;  TD%face_connectivity(4,6,nid)=0
                        exit
                    end if
                end do
            end if
        end do

        ! Inter-patch connectivity via geometric matching
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
                        else if (norm2(centroids(:,f1,e1)-centroids(:,f2,e2)) < 1e-4_dp) then
                            if (dot_product(TD%face_normals(:,f1,e1), &
                                            TD%face_normals(:,f2,e2)) < -0.9_dp) found_neighbor = .true.
                        end if
                        if (found_neighbor) then
                            orient_val = face_orient_from_nodes(mesh, FE, e1, f1, e2, f2)
                            TD%face_connectivity(1,f1,e1)=e2; TD%face_connectivity(2,f1,e1)=f2
                            TD%face_connectivity(1,f2,e2)=e1; TD%face_connectivity(2,f2,e2)=f1
                            TD%face_connectivity(3,f1,e1)=orient_val
                            TD%face_connectivity(3,f2,e2)=orient_val
                            TD%face_connectivity(4,f1,e1)=0;  TD%face_connectivity(4,f2,e2)=0
                            exit
                        end if
                    end do
                    if (TD%face_connectivity(1,f1,e1) /= -1) exit
                end do
            end do
        end do
        deallocate(centroids)

        ! Connectivity symmetry check
        do e1 = 1, mesh%n_elems
            do f1 = 1, nf
                e2 = TD%face_connectivity(1,f1,e1)
                if (e2 > 0) then
                    f2 = TD%face_connectivity(2,f1,e1)
                    if (TD%face_connectivity(1,f2,e2) /= e1) &
                        write(*,'(A,4I5)') "[CRITICAL] Asymmetric connectivity: e1,f1,e2,f2=",e1,f1,e2,f2
                end if
            end do
        end do

        ! Assign boundary BC from iga_surfaces
        do e1 = 1, mesh%n_elems
            do f1 = 1, nf
                if (TD%face_connectivity(1,f1,e1) == -1) then
                    found_neighbor = .false.
                    do s_idx = 1, size(mesh%iga_surfaces)
                        if (all_nodes_in_list(mesh%elems(e1, FE%face_node_map(:,f1)), &
                                              mesh%iga_surfaces(s_idx)%cp_ids)) then
                            TD%face_connectivity(4,f1,e1) = mesh%iga_surfaces(s_idx)%bc_id
                            found_neighbor = .true.
                            exit
                        end if
                    end do
                    if (.not. found_neighbor) TD%face_connectivity(4,f1,e1) = BC_VACUUM
                end if
            end do
        end do
    end subroutine connectivity_and_normals

    subroutine precompute_upwind_indices(n_elems, n_faces_per_elem, n_nodes_per_face, n_basis, &
                                         face_node_map, face_connectivity, upwind_idx)
        integer, intent(in)                       :: n_elems, n_faces_per_elem
        integer, intent(in)                       :: n_nodes_per_face, n_basis
        integer, intent(in)                       :: face_node_map(:,:)
        integer, intent(in)                       :: face_connectivity(:,:,:)
        integer, allocatable, intent(out)         :: upwind_idx(:,:,:)

        integer :: ee, f, j_f, j_nf, nid, n_fac, orient

        allocate(upwind_idx(n_nodes_per_face, n_faces_per_elem, n_elems))
        upwind_idx = 0

        do ee = 1, n_elems
            do f = 1, n_faces_per_elem
                nid = face_connectivity(1, f, ee)
                if (nid > 0) then
                    n_fac  = face_connectivity(2, f, ee)
                    orient = face_connectivity(3, f, ee)
                    do j_f = 1, n_nodes_per_face
                        j_nf = merge(n_nodes_per_face - j_f + 1, j_f, orient /= 1)
                        upwind_idx(j_f,f,ee) = (nid-1)*n_basis + face_node_map(j_nf, n_fac)
                    end do
                else
                    do j_f = 1, n_nodes_per_face
                        upwind_idx(j_f,f,ee) = (ee-1)*n_basis + face_node_map(j_f, f)
                    end do
                end if
            end do
        end do
    end subroutine precompute_upwind_indices

    subroutine generate_sweep_order(n_elems, n_faces_per_elem, face_normals, &
                                    face_connectivity, direction, sweep_order)
        integer,             intent(in)            :: n_elems, n_faces_per_elem
        real(dp),            intent(in)            :: face_normals(:,:,:)
        integer,             intent(in)            :: face_connectivity(:,:,:)
        real(dp),            intent(in)            :: direction(3)
        integer, contiguous, intent(out)           :: sweep_order(:)

        integer :: e1, e2, f1, nid
        integer, allocatable :: queue(:), incoming(:)
        integer :: head, tail, sweep_idx, level_end

        allocate(queue(n_elems), incoming(n_elems))
        incoming = 0; head = 1; tail = 0

        do e1 = 1, n_elems
            do f1 = 1, n_faces_per_elem
                if (dot_product(face_normals(:,f1,e1), direction) < -1e-12_dp) then
                    if (face_connectivity(1,f1,e1) > 0) incoming(e1) = incoming(e1) + 1
                end if
            end do
            if (incoming(e1) == 0) then
                tail = tail + 1; queue(tail) = e1
            end if
        end do

        sweep_idx = 0; level_end = tail
        do while (head <= tail)
            do while (head <= level_end)
                e1 = queue(head); head = head + 1
                sweep_idx = sweep_idx + 1; sweep_order(sweep_idx) = e1
                do f1 = 1, n_faces_per_elem
                    if (dot_product(face_normals(:,f1,e1), direction) > 1e-12_dp) then
                        e2 = face_connectivity(1,f1,e1)
                        if (e2 > 0) then
                            incoming(e2) = incoming(e2) - 1
                            if (incoming(e2) == 0) then
                                tail = tail + 1; queue(tail) = e2
                            end if
                        end if
                    end if
                end do
            end do
            level_end = tail
        end do

        if (sweep_idx /= n_elems) then
            write(*,'(A,I0,A,I0)') "Sweep Error: processed ", sweep_idx, " / ", n_elems
            do e1 = 1, n_elems
                if (incoming(e1) > 0) then
                    write(*,'(A,I0,A,I0,A)') "Element ", e1, " stalled (waiting ", incoming(e1), ")"
                    do f1 = 1, n_faces_per_elem
                        if (dot_product(face_normals(:,f1,e1), direction) < -1e-12_dp) then
                            nid = face_connectivity(1,f1,e1)
                            if (nid > 0) write(*,'(A,2I5)') "  face -> elem", f1, nid
                        end if
                    end do
                end if
            end do
            stop "STOP: cycle detected in sweep graph"
        end if
        deallocate(queue, incoming)
    end subroutine generate_sweep_order

    subroutine precompute_reflective_map(n_elems, n_faces_per_elem, sn_quad, face_normals, reflect_map)
        integer,               intent(in)    :: n_elems, n_faces_per_elem
        type(t_sn_quadrature), intent(in)    :: sn_quad
        real(dp),              intent(in)    :: face_normals(:,:,:)
        integer, allocatable,  intent(out)   :: reflect_map(:,:,:)

        integer :: ee, f, mm, m_iter
        real(dp) :: normal(3), dir(3), ref_dir(3), max_dot, dprod

        allocate(reflect_map(sn_quad%n_angles, n_faces_per_elem, n_elems))
        reflect_map = 0

        !$OMP PARALLEL DO PRIVATE(ee, f, normal, mm, dir, ref_dir, max_dot, m_iter, dprod)
        do ee = 1, n_elems
            do f = 1, n_faces_per_elem
                normal = face_normals(:,f,ee)
                do mm = 1, sn_quad%n_angles
                    dir = sn_quad%dirs(mm, :)
                    ref_dir = dir - 2.0_dp * dot_product(dir, normal) * normal
                    max_dot = -2.0_dp
                    do m_iter = 1, sn_quad%n_angles
                        if (abs(ref_dir(3) - sn_quad%dirs(m_iter,3)) > SMALL_NUMBER) cycle
                        dprod = dot_product(ref_dir, sn_quad%dirs(m_iter,:))
                        if (dprod > max_dot) then
                            max_dot = dprod
                            reflect_map(mm,f,ee) = m_iter
                        end if
                    end do
                end do
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine precompute_reflective_map

    integer function face_orient_from_nodes(mesh, FE, e1, f1, e2, f2) result(orient)
        type(t_mesh_iga),   intent(in) :: mesh
        type(t_finite_iga), intent(in) :: FE
        integer,            intent(in) :: e1, f1, e2, f2
        integer :: first_node_e1
        first_node_e1 = mesh%elems(e1, FE%face_node_map(1, f1))
        if (mesh%elems(e2, FE%face_node_map(1, f2)) == first_node_e1) then
            orient = 1
        else
            orient = -1
        end if
    end function face_orient_from_nodes

    function all_nodes_in_list(subset, superset) result(is_subset)
        integer, intent(in) :: subset(:), superset(:)
        logical :: is_subset
        integer :: ii
        is_subset = .true.
        do ii = 1, size(subset)
            if (.not. any(superset == subset(ii))) then
                is_subset = .false.; return
            end if
        end do
    end function all_nodes_in_list

end module m_sweep_order
