! ASMG format mesh reader for IGA (2D and 3D).
! Unified mesh reader for IGA and FEM (FEM uses trivial open knot vectors).
!
! ASMG format summary:
!   dim=2:  $2D_Patch_Description_Start/End  -- volume patches
!           $1D_Patch_Description_Start/End  -- boundary edges
!   dim=3:  $3D_Patch_Description_Start/End  -- volume patches
!           $2D_Patch_Description_Start/End  -- boundary surface faces
!   Header keywords: PolyOrder, Groups, Dims
!   POINTS block: x y z weight per control point (one per line)
!
! Public: read_asmg_mesh, write_mesh_to_files`
module m_asmg
    use m_constants
    use m_types_iga
    implicit none
    private
    public :: read_asmg_mesh, write_mesh_to_files

contains

    ! ------------------------------------------------------------------
    ! Top-level ASMG mesh reader.
    ! Performs two passes over the file:
    !   Pass 1: count patches/surfaces, read header fields, find max knot size.
    !   Pass 2: parse each block into mesh%patches(:) and mesh%surfaces(:).
    ! Auto-detects dim from presence of $3D_Patch blocks.
    ! ------------------------------------------------------------------
    subroutine read_asmg_mesh(filepath, mesh)
        character(len=*),       intent(in)    :: filepath
        type(t_mesh_iga), allocatable, intent(inout) :: mesh

        integer :: unit, ios, ii, jj, kk, k, p_idx
        integer :: s_u, s_v, s_w, cp_idx, n_tot
        character(len=1024) :: line
        integer :: n_vol, n_bnd, max_k, p, q, r, nxi, neta, pos
        integer, allocatable :: tmp_cp(:)
        real(dp) :: xi_t(2), eta_t(2), cross_z

        if (.not. allocated(mesh)) allocate(mesh)

        ! --- Pass 1: scan for counts and header values ---
        open(newunit=unit, file=filepath, status='old', action='read')
        mesh%order = 1; mesh%n_nodes = 0; n_vol = 0; n_bnd = 0; max_k = 0
        mesh%dim = 2; mesh%n_groups = 1

        do
            read(unit,'(A)',iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)

            if (index(line,'$3D_Patch_Description_Start') > 0) then
                n_vol = n_vol + 1; mesh%dim = 3
            else if (index(line,'$2D_Patch_Description_Start') > 0) then
                if (mesh%dim == 3) then
                    n_bnd = n_bnd + 1   ! 3D boundary surface
                else
                    n_vol = n_vol + 1   ! 2D volume patch
                end if
            else if (index(line,'$1D_Patch_Description_Start') > 0) then
                n_bnd = n_bnd + 1       ! 2D boundary edge
            end if

            if (index(line,'POINTS') == 1) read(line(7:),*,iostat=ios) mesh%n_nodes

            pos = scan(line,': ')
            if (pos > 0) then
                if (index(line,'PolyOrder') > 0) read(line(pos+1:),*,iostat=ios) mesh%order
                if (index(line,'Groups')    > 0) read(line(pos+1:),*,iostat=ios) mesh%n_groups
                if (index(line,'Dims')      > 0) read(line(pos+1:),*,iostat=ios) mesh%dim
                if (index(line,'KnotVector') > 0) then
                    read(line(pos+1:),*,iostat=ios) k
                    if (ios == 0 .and. k > max_k) max_k = k
                end if
            end if
        end do

        ! A second scan is needed if $3D patches were found after $2D ones.
        ! Re-count with correct dim.
        if (mesh%dim == 3 .and. n_vol == 0) then
            rewind(unit)
            n_vol = 0; n_bnd = 0
            do
                read(unit,'(A)',iostat=ios) line
                if (ios /= 0) exit
                line = adjustl(line)
                if (index(line,'$3D_Patch_Description_Start') > 0) n_vol = n_vol + 1
                if (index(line,'$2D_Patch_Description_Start') > 0) n_bnd = n_bnd + 1
            end do
        end if

        ! For FEM ASMG files (no KnotVector lines) max_k stays 0; ensure the
        ! elem_map arrays are large enough to hold the span index p+1.
        max_k = max(max_k + 5, 2*(mesh%order + 1))

        p = mesh%order; q = p; r = p
        mesh%n_faces_per_elem = merge(6, 4, mesh%dim == 3)

        allocate(mesh%nodes(mesh%n_nodes, 3), mesh%weights(mesh%n_nodes))
        allocate(mesh%patches(n_vol), mesh%surfaces(n_bnd))
        do k = 1, n_vol; mesh%patches(k)%face_to_surface = 0; end do

        ! --- Pass 2: parse blocks ---
        rewind(unit)
        ii = 0; jj = 0

        do
            read(unit,'(A)',iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0 .or. line(1:1) == '!') cycle

            if (index(line,'POINTS') == 1) then
                do k = 1, mesh%n_nodes
                    read(unit,*) mesh%nodes(k,1:3), mesh%weights(k)
                end do

            else if (index(line,'$3D_Patch_Description_Start') > 0) then
                jj = jj + 1
                call parse_block(unit, mesh%patches(jj)%cp_ids, mesh%patches(jj)%material_id, &
                                 mesh%patches(jj)%knots_xi, mesh%patches(jj)%knots_eta, &
                                 mesh%patches(jj)%knots_zeta)
                if (.not. allocated(mesh%patches(jj)%knots_xi))   call gen_trivial_knots(p, mesh%patches(jj)%knots_xi)
                if (.not. allocated(mesh%patches(jj)%knots_eta))  call gen_trivial_knots(q, mesh%patches(jj)%knots_eta)
                if (.not. allocated(mesh%patches(jj)%knots_zeta)) call gen_trivial_knots(r, mesh%patches(jj)%knots_zeta)

            else if (index(line,'$2D_Patch_Description_Start') > 0) then
                if (mesh%dim == 3) then
                    ! 3D boundary surface
                    ii = ii + 1
                    call parse_block(unit, mesh%surfaces(ii)%cp_ids, mesh%surfaces(ii)%bc_id, &
                                     mesh%surfaces(ii)%knots_xi, mesh%surfaces(ii)%knots_eta)
                    if (.not. allocated(mesh%surfaces(ii)%knots_xi))  call gen_trivial_knots(p, mesh%surfaces(ii)%knots_xi)
                    if (.not. allocated(mesh%surfaces(ii)%knots_eta)) call gen_trivial_knots(q, mesh%surfaces(ii)%knots_eta)
                else
                    ! 2D volume patch
                    jj = jj + 1
                    call parse_block(unit, mesh%patches(jj)%cp_ids, mesh%patches(jj)%material_id, &
                                     mesh%patches(jj)%knots_xi, mesh%patches(jj)%knots_eta)
                    if (.not. allocated(mesh%patches(jj)%knots_xi))  call gen_trivial_knots(p, mesh%patches(jj)%knots_xi)
                    if (.not. allocated(mesh%patches(jj)%knots_eta)) call gen_trivial_knots(q, mesh%patches(jj)%knots_eta)
                end if

            else if (index(line,'$1D_Patch_Description_Start') > 0) then
                ! 2D boundary edge
                ii = ii + 1
                call parse_block(unit, mesh%surfaces(ii)%cp_ids, mesh%surfaces(ii)%bc_id, &
                                 mesh%surfaces(ii)%knots_xi)
                if (.not. allocated(mesh%surfaces(ii)%knots_xi)) call gen_trivial_knots(p, mesh%surfaces(ii)%knots_xi)
            end if
        end do
        close(unit)

        ! --- Precompute 3D surface elements (span-indexed, mirroring volume) ---
        if (mesh%dim == 3) then
            do k = 1, n_bnd
                call build_surface_elems(mesh%surfaces(k), p)
            end do
        end if

        ! --- Ensure consistent eta orientation for 2D patches ---
        ! For each patch, compute the cross product of the xi and eta tangent
        ! vectors at the first corner.  If the z-component is negative the patch
        ! is parameterised with eta decreasing in y (inverted orientation), which
        ! makes face normals point inward and causes sweep deadlocks at higher SN.
        ! Fix: reverse the row ordering of cp_ids (knot vector unchanged).
        if (mesh%dim == 2) then
            do p_idx = 1, n_vol
                associate(ptch => mesh%patches(p_idx))
                if (.not. allocated(ptch%cp_ids)) cycle
                nxi  = size(ptch%knots_xi)  - p - 1
                neta = size(ptch%knots_eta) - q - 1
                if (nxi < 2 .or. neta < 2) cycle
                if (size(ptch%cp_ids) < nxi + 1) cycle
                xi_t  = mesh%nodes(ptch%cp_ids(2),     1:2) - mesh%nodes(ptch%cp_ids(1), 1:2)
                eta_t = mesh%nodes(ptch%cp_ids(nxi+1), 1:2) - mesh%nodes(ptch%cp_ids(1), 1:2)
                cross_z = xi_t(1)*eta_t(2) - xi_t(2)*eta_t(1)
                if (cross_z < -dp_EPSILON) then
                    allocate(tmp_cp(nxi))
                    do ii = 1, neta/2
                        jj = neta + 1 - ii
                        tmp_cp = ptch%cp_ids((ii-1)*nxi+1 : ii*nxi)
                        ptch%cp_ids((ii-1)*nxi+1 : ii*nxi) = ptch%cp_ids((jj-1)*nxi+1 : jj*nxi)
                        ptch%cp_ids((jj-1)*nxi+1 : jj*nxi) = tmp_cp
                    end do
                    deallocate(tmp_cp)
                    write(*,'(A,I6,A)') "  [ASMG] Patch ", p_idx, " eta reversed (was inverted)."
                end if
                end associate
            end do
        end if

        ! --- Subdivide patches into knot-span elements ---
        n_tot = 0
        do p_idx = 1, n_vol
            associate(ptch => mesh%patches(p_idx))
            do s_u = p+1, size(ptch%knots_xi) - p - 1
                if (ptch%knots_xi(s_u+1) <= ptch%knots_xi(s_u) + dp_EPSILON) cycle
                do s_v = q+1, size(ptch%knots_eta) - q - 1
                    if (ptch%knots_eta(s_v+1) <= ptch%knots_eta(s_v) + dp_EPSILON) cycle
                    if (mesh%dim == 3) then
                        do s_w = r+1, size(ptch%knots_zeta) - r - 1
                            if (ptch%knots_zeta(s_w+1) > ptch%knots_zeta(s_w) + dp_EPSILON) n_tot = n_tot + 1
                        end do
                    else
                        n_tot = n_tot + 1
                    end if
                end do
            end do
            end associate
        end do

        mesh%n_elems = n_tot
        k = merge((p+1)**3, (p+1)**2, mesh%dim == 3)

        allocate(mesh%elems(n_tot, k), mesh%material_ids(n_tot))
        allocate(mesh%elem_patch_id(n_tot))
        allocate(mesh%elem_span_indices(mesh%dim, n_tot))
        allocate(mesh%elem_u_min(n_tot), mesh%elem_u_max(n_tot))
        allocate(mesh%elem_v_min(n_tot), mesh%elem_v_max(n_tot))

        if (mesh%dim == 3) then
            allocate(mesh%elem_w_min(n_tot), mesh%elem_w_max(n_tot))
            allocate(mesh%elem_map_3d(n_vol, max_k, max_k, max_k))
            mesh%elem_map_3d = 0
        else
            allocate(mesh%elem_map_2d(n_vol, max_k, max_k))
            mesh%elem_map_2d = 0
        end if

        k = 0
        do p_idx = 1, n_vol
            associate(ptch => mesh%patches(p_idx))
            nxi  = size(ptch%knots_xi) - p - 1
            neta = size(ptch%knots_eta) - q - 1

            do s_u = p+1, size(ptch%knots_xi) - p - 1
                if (ptch%knots_xi(s_u+1) <= ptch%knots_xi(s_u) + dp_EPSILON) cycle
                do s_v = q+1, size(ptch%knots_eta) - q - 1
                    if (ptch%knots_eta(s_v+1) <= ptch%knots_eta(s_v) + dp_EPSILON) cycle

                    if (mesh%dim == 3) then
                        do s_w = r+1, size(ptch%knots_zeta) - r - 1
                            if (ptch%knots_zeta(s_w+1) <= ptch%knots_zeta(s_w) + dp_EPSILON) cycle

                            k = k + 1
                            mesh%elem_patch_id(k)       = p_idx
                            mesh%elem_span_indices(1,k) = s_u
                            mesh%elem_span_indices(2,k) = s_v
                            mesh%elem_span_indices(3,k) = s_w
                            mesh%elem_map_3d(p_idx, s_u, s_v, s_w) = k
                            mesh%material_ids(k)  = ptch%material_id
                            mesh%elem_u_min(k)    = ptch%knots_xi(s_u)
                            mesh%elem_u_max(k)    = ptch%knots_xi(s_u+1)
                            mesh%elem_v_min(k)    = ptch%knots_eta(s_v)
                            mesh%elem_v_max(k)    = ptch%knots_eta(s_v+1)
                            mesh%elem_w_min(k)    = ptch%knots_zeta(s_w)
                            mesh%elem_w_max(k)    = ptch%knots_zeta(s_w+1)

                            cp_idx = 0
                            do kk = s_w - r, s_w
                                do jj = s_v - q, s_v
                                    do ii = s_u - p, s_u
                                        cp_idx = cp_idx + 1
                                        mesh%elems(k, cp_idx) = ptch%cp_ids( &
                                            (kk-1)*nxi*neta + (jj-1)*nxi + ii)
                                    end do
                                end do
                            end do
                        end do

                    else  ! dim == 2
                        k = k + 1
                        mesh%elem_patch_id(k)       = p_idx
                        mesh%elem_span_indices(1,k) = s_u
                        mesh%elem_span_indices(2,k) = s_v
                        mesh%elem_map_2d(p_idx, s_u, s_v) = k
                        mesh%material_ids(k)  = ptch%material_id
                        mesh%elem_u_min(k)    = ptch%knots_xi(s_u)
                        mesh%elem_u_max(k)    = ptch%knots_xi(s_u+1)
                        mesh%elem_v_min(k)    = ptch%knots_eta(s_v)
                        mesh%elem_v_max(k)    = ptch%knots_eta(s_v+1)

                        cp_idx = 0
                        do jj = s_v - q, s_v
                            do ii = s_u - p, s_u
                                cp_idx = cp_idx + 1
                                mesh%elems(k, cp_idx) = ptch%cp_ids((jj-1)*nxi + ii)
                            end do
                        end do
                    end if
                end do
            end do
            end associate
        end do

    end subroutine read_asmg_mesh

    ! ------------------------------------------------------------------
    ! Parse a single $xD_Patch block from the open file unit u.
    ! k1 (xi) is always required; k2 (eta) is for 2D+ patches;
    ! k3 (zeta) is optional and used only for 3D volume patches.
    ! ------------------------------------------------------------------
    subroutine parse_block(u, cp, id, k1, k2, k3)
        integer,                       intent(in)           :: u
        integer, allocatable,          intent(out)          :: cp(:)
        integer,                       intent(out)          :: id
        real(dp), allocatable,         intent(out)          :: k1(:)
        real(dp), allocatable,         intent(out), optional :: k2(:)
        real(dp), allocatable,         intent(out), optional :: k3(:)

        character(len=1024) :: l, key
        integer :: n, ios
        real(dp), allocatable :: temp_k(:)

        id = -1

        do
            read(u,'(A)',iostat=ios) l
            if (ios /= 0 .or. index(l,'_End') > 0) exit
            l = adjustl(l)
            if (l == '' .or. l(1:1) == '!') cycle

            read(l,*,iostat=ios) key, n

            if (index(l,'control_points') > 0) then
                if (ios /= 0) read(u,*) n
                allocate(cp(n)); read(u,*) cp; cp = cp + 1

            else if (index(l,'KnotVector_Xi') > 0 .or. index(l,'KnotVectorXi') > 0) then
                if (ios /= 0) read(u,*) n
                if (allocated(k1)) deallocate(k1)
                allocate(k1(n)); read(u,*) k1

            else if (index(l,'KnotVector_Eta') > 0 .or. index(l,'KnotVectorEta') > 0) then
                if (.not. present(k2)) cycle
                if (ios /= 0) read(u,*) n
                if (allocated(k2)) deallocate(k2)
                allocate(k2(n)); read(u,*) k2

            else if (index(l,'KnotVector_Zeta') > 0) then
                if (.not. present(k3)) cycle
                if (ios /= 0) read(u,*) n
                if (allocated(k3)) deallocate(k3)
                allocate(k3(n)); read(u,*) k3

            else if (index(l,'KnotVector') > 0) then
                ! Generic single knot vector — copy to all directions
                if (ios /= 0) read(u,*) n
                allocate(temp_k(n)); read(u,*) temp_k
                if (.not. allocated(k1)) allocate(k1(n), source=temp_k)
                if (present(k2)) then
                    if (.not. allocated(k2)) allocate(k2(n), source=temp_k)
                end if
                if (present(k3)) then
                    if (.not. allocated(k3)) allocate(k3(n), source=temp_k)
                end if
                deallocate(temp_k)

            else if (index(l,'Material_ID') > 0 .or. index(l,'BC') > 0) then
                if (ios == 0) then
                    id = n
                else
                    read(u,*,iostat=ios) id
                end if
            end if
        end do
    end subroutine parse_block

    ! ------------------------------------------------------------------
    ! Write mesh summary files to outdir for diagnostics.
    ! ------------------------------------------------------------------
    subroutine write_mesh_to_files(mesh, outdir)
        type(t_mesh_iga),  intent(in) :: mesh
        character(len=*),  intent(in) :: outdir
        integer :: u, ii

        open(newunit=u, file=trim(outdir)//'/nodes.dat', status='replace')
        write(u,'(A)') "# ID | X | Y | Z | Weight"
        do ii = 1, mesh%n_nodes
            write(u,'(I8,4F15.8)') ii, mesh%nodes(ii,1:3), mesh%weights(ii)
        end do
        close(u)

        open(newunit=u, file=trim(outdir)//'/elements.dat', status='replace')
        write(u,'(A)') "# Element ID | Control Point IDs..."
        do ii = 1, mesh%n_elems
            write(u,'(I8,A,500I8)') ii, " : ", mesh%elems(ii,:)
        end do
        close(u)

        open(newunit=u, file=trim(outdir)//'/edges.dat', status='replace')
        write(u,'(A)') "# Surface ID | BC_ID | Control Point IDs..."
        do ii = 1, size(mesh%surfaces)
            write(u,'(I8,A,I4,A,500I8)') ii, " BC:", mesh%surfaces(ii)%bc_id, " : ", mesh%surfaces(ii)%cp_ids
        end do
        close(u)

        open(newunit=u, file=trim(outdir)//'/materials.dat', status='replace')
        write(u,'(A)') "# Elem ID | Material ID"
        do ii = 1, mesh%n_elems
            write(u,'(2I10)') ii, mesh%material_ids(ii)
        end do
        close(u)

        open(newunit=u, file=trim(outdir)//'/knotspans.dat', status='replace')
        if (mesh%dim == 3) then
            write(u,'(A)') "# Elem ID | U_range | V_range | W_range"
            do ii = 1, mesh%n_elems
                write(u,'(I10,6F12.6)') ii, &
                    mesh%elem_u_min(ii), mesh%elem_u_max(ii), &
                    mesh%elem_v_min(ii), mesh%elem_v_max(ii), &
                    mesh%elem_w_min(ii), mesh%elem_w_max(ii)
            end do
        else
            write(u,'(A)') "# Elem ID | U_range | V_range"
            do ii = 1, mesh%n_elems
                write(u,'(I10,4F12.6)') ii, &
                    mesh%elem_u_min(ii), mesh%elem_u_max(ii), &
                    mesh%elem_v_min(ii), mesh%elem_v_max(ii)
            end do
        end if
        close(u)
    end subroutine write_mesh_to_files

    ! ------------------------------------------------------------------
    ! Build surface element connectivity for a single 3D boundary surface.
    ! Each non-zero-measure (xi,eta) span pair becomes one surface element.
    ! surf%elems(ee,:) holds global CP IDs in xi-fastest order.
    ! ------------------------------------------------------------------
    subroutine build_surface_elems(surf, p)
        type(t_surface_iga), intent(inout) :: surf
        integer,             intent(in)    :: p

        integer :: n_cp_xi, s_u, s_v, ee, cp_idx, ii, jj
        integer :: n_elem

        n_cp_xi = size(surf%knots_xi) - p - 1

        ! Count non-zero-measure span pairs
        n_elem = 0
        do s_u = p+1, size(surf%knots_xi) - p - 1
            if (surf%knots_xi(s_u+1) <= surf%knots_xi(s_u) + dp_EPSILON) cycle
            do s_v = p+1, size(surf%knots_eta) - p - 1
                if (surf%knots_eta(s_v+1) > surf%knots_eta(s_v) + dp_EPSILON) n_elem = n_elem + 1
            end do
        end do

        surf%n_elements = n_elem
        allocate(surf%elem_span_indices(2, n_elem))
        allocate(surf%elems(n_elem, (p+1)*(p+1)))

        ee = 0
        do s_u = p+1, size(surf%knots_xi) - p - 1
            if (surf%knots_xi(s_u+1) <= surf%knots_xi(s_u) + dp_EPSILON) cycle
            do s_v = p+1, size(surf%knots_eta) - p - 1
                if (surf%knots_eta(s_v+1) <= surf%knots_eta(s_v) + dp_EPSILON) cycle
                ee = ee + 1
                surf%elem_span_indices(1, ee) = s_u
                surf%elem_span_indices(2, ee) = s_v
                cp_idx = 0
                do jj = s_v - p, s_v
                    do ii = s_u - p, s_u
                        cp_idx = cp_idx + 1
                        surf%elems(ee, cp_idx) = surf%cp_ids((jj-1)*n_cp_xi + ii)
                    end do
                end do
            end do
        end do
    end subroutine build_surface_elems

    ! ------------------------------------------------------------------
    ! Generate a trivial open knot vector for a single-span FEM patch.
    ! Result: [0,0,...,0, 1,1,...,1] with p+1 zeros and p+1 ones.
    ! ------------------------------------------------------------------
    subroutine gen_trivial_knots(p, kv)
        integer,              intent(in)  :: p
        real(dp), allocatable, intent(out) :: kv(:)
        allocate(kv(2*(p+1)))
        kv(1:p+1)        = 0.0_dp
        kv(p+2:2*(p+1)) = 1.0_dp
    end subroutine gen_trivial_knots

end module m_asmg
