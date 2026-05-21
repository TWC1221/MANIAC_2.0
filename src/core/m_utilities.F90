! Core utility routines: string helpers, splash output, mesh quality diagnostics.
module m_utilities
    use m_constants
    use m_types
    implicit none

    public :: print_splash, derive_case_nametag, int_to_str, check_mesh_quality

    ! 74-character wide report box  (2 + 71 content + 1)
    character(len=74), parameter :: &
        BOX_H = " |=======================================================================|", &
        BOX_S = " |-----------------------------------------------------------------------|"

contains

    ! ------------------------------------------------------------------
    pure function int_to_str(ii) result(res)
        integer, intent(in) :: ii
        character(len=12) :: res
        write(res, '(I0)') ii
        res = adjustl(res)
    end function int_to_str

    pure function derive_case_nametag(filename) result(tag)
        character(len=*), intent(in) :: filename
        character(len=:), allocatable :: tag
        integer :: pos, dot_pos
        pos = index(filename, '/', back=.true.)
        if (pos == 0) pos = index(filename, '\', back=.true.)
        tag = filename(pos+1:)
        dot_pos = index(tag, '.', back=.true.)
        if (dot_pos > 0) then
            tag = tag(1:dot_pos-1) // ".vtk"
        else
            tag = tag // ".vtk"
        end if
    end function derive_case_nametag

    ! ------------------------------------------------------------------
    ! Print ASCII logo then a run-configuration summary box.
    ! All key config values are passed in; solver/PC are integers mapped
    ! to human-readable names.  sn_order/solver_type/preconditioner are
    ! optional so diffusion and transport callers can omit irrelevant fields.
    ! ------------------------------------------------------------------
    subroutine print_splash(solver, method, mesh_file, mat_file, output_dir, &
                             max_outer, tol, is_eigenvalue, is_adjoint, &
                             sn_order, solver_type, preconditioner)
        character(len=*), intent(in)           :: solver, method, mesh_file, mat_file, output_dir
        integer,          intent(in)           :: max_outer
        real(dp),         intent(in)           :: tol
        logical,          intent(in)           :: is_eigenvalue, is_adjoint
        integer, optional, intent(in)          :: sn_order, solver_type, preconditioner

        integer :: i
        character(len=148), dimension(9) :: logo
        character(len=8)  :: date_str
        character(len=10) :: time_str
        character(len=16) :: ksp_name, pc_name

        logo(1) = "__/\\\\\_____/\\\________/\\\\\\\\\\\________/\\\_______/\\\_____________/\\\\\_____________/\\\\\_____/\\\_        "
        logo(2) = " _\/\\\\\\___\/\\\_______\/////\\\///________\///\\\___/\\\/____________/\\\///\\\__________\/\\\\\\___\/\\\_       "
        logo(3) = "  _\/\\\/\\\__\/\\\___________\/\\\_____________\///\\\\\\/____________/\\\/__\///\\\________\/\\\/\\\__\/\\\_      "
        logo(4) = "   _\/\\\//\\\_\/\\\___________\/\\\_______________\//\\\\_____________/\\\______\//\\\_______\/\\\//\\\_\/\\\_     "
        logo(5) = "    _\/\\\\//\\\\/\\\___________\/\\\________________\/\\\\____________\/\\\_______\/\\\_______\/\\\\//\\\\/\\\_    "
        logo(6) = "     _\/\\\_\//\\\/\\\___________\/\\\________________/\\\\\\___________\//\\\______/\\\________\/\\\_\//\\\/\\\_   "
        logo(7) = "      _\/\\\__\//\\\\\\___________\/\\\______________/\\\////\\\__________\///\\\__/\\\__________\/\\\__\//\\\\\\_  "
        logo(8) = "       _\/\\\___\//\\\\\__/\\\__/\\\\\\\\\\\__/\\\__/\\\/___\///\\\__/\\\____\///\\\\\/______/\\\_\/\\\___\//\\\\\_ "
        logo(9) = "        _\///_____\/////__\///__\///////////__\///__\///_______\///__\///_______\/////_______\///__\///_____\/////__"

        call date_and_time(date_str, time_str)

        write(*,*)
        do i = 1, 9; write(*,'(A)') logo(i); end do
        write(*,*)
        write(*,'(A)') "            N O D A L    I G A - F E M    E X A S C A L E    O R D I N A T E S - D I F F U S I O N    N E U T R O N I C S"

        ! Solver/PC strings
        if (present(solver_type)) then
            select case (solver_type)
            case (SOLVER_KSP_CG);    ksp_name = "CG"
            case (SOLVER_KSP_GMRES); ksp_name = "GMRES"
            case (SOLVER_KSP_BCGS);  ksp_name = "BiCGStab"
            case default;            ksp_name = "UNKNOWN"
            end select
        else
            ksp_name = "N/A"
        end if
        if (present(preconditioner)) then
            select case (preconditioner)
            case (PRECON_NONE);    pc_name = "NONE"
            case (PRECON_JACOBI);  pc_name = "JACOBI"
            case (PRECON_ILU);     pc_name = "ILU"
            case (PRECON_CHOLESKY);pc_name = "ICC"
            case (PRECON_GAMG);    pc_name = "GAMG"
            case default;          pc_name = "UNKNOWN"
            end select
        else
            pc_name = "N/A"
        end if

        write(*,*)
        write(*,'(A)') BOX_H
        write(*,'(A, A4, "-", A2, "-", A2, " @ ", A2, ":", A2, ":", A2, T74, A)') &
            " |  Timestamp      :  ", date_str(1:4), date_str(5:6), date_str(7:8), &
            time_str(1:2), time_str(3:4), time_str(5:6), "|"
        write(*,'(A)') BOX_S
        write(*,'(A, A, T74, A)') " |  Solver         :  ", trim(solver),    "|"
        write(*,'(A, A, T74, A)') " |  Discretisation :  ", trim(method),    "|"
        write(*,'(A, A, T74, A)') " |  Mesh File      :  ", trim(mesh_file), "|"
        write(*,'(A, A, T74, A)') " |  Material File  :  ", trim(mat_file),  "|"

        if (present(sn_order) .and. trim(solver) == "transport") &
            write(*,'(A, I0, T74, A)') " |  Sn Order       :  ", sn_order, "|"
        if (present(solver_type) .and. trim(solver) == "diffusion") &
            write(*,'(A, A, " / ", A, T74, A)') " |  Linear Solver  :  ", trim(ksp_name), trim(pc_name), "|"
        write(*,'(A)') BOX_H
        write(*,*)
    end subroutine print_splash

    ! ==================================================================
    ! Mesh quality checker — works for both IGA and FEM via class(t_mesh).
    !
    ! Element quality is assessed using bilinear (2D) or trilinear (3D)
    ! interpolation at corner control points.  For IGA meshes this is an
    ! approximation over the control polygon; the NURBS geometry can be
    ! smoother but the control polygon reveals gross quality issues.
    !
    ! Optional arguments:
    !   out_unit  — write report to this Fortran unit (default: stdout)
    !   log_path  — also write to this file path
    ! ==================================================================
    subroutine check_mesh_quality(mesh, out_unit, log_path)
        class(t_mesh),    intent(in)           :: mesh
        integer,          intent(in), optional :: out_unit
        character(len=*), intent(in), optional :: log_path

        integer :: u_arr(2), n_u, u_log, u, ios
        integer :: dup_count, bad_jac_count, degenerate_count, high_skew_count
        real(dp) :: g_min_jac, g_avg_jac, g_max_jac
        real(dp) :: g_min_skew, g_max_skew, g_avg_skew
        real(dp) :: g_min_ar, g_max_ar, g_avg_ar
        real(dp) :: total_meas

        u_arr(1) = 6
        if (present(out_unit)) u_arr(1) = out_unit
        n_u = 1
        if (present(log_path)) then
            open(newunit=u_log, file=trim(log_path), status='replace', action='write', iostat=ios)
            if (ios == 0) then; n_u = 2; u_arr(2) = u_log; end if
        end if

        call compute_elem_quality(mesh, &
            g_min_jac, g_avg_jac, g_max_jac, total_meas,     &
            g_min_skew, g_avg_skew, g_max_skew,               &
            g_min_ar,   g_avg_ar,   g_max_ar,                 &
            bad_jac_count, degenerate_count, high_skew_count)
        call count_coincident_nodes(mesh, dup_count)

        do u = 1, n_u
            call print_quality_report(mesh, u_arr(u),                           &
                g_min_jac, g_avg_jac, g_max_jac, total_meas,   &
                g_min_skew, g_avg_skew, g_max_skew,             &
                g_min_ar,   g_avg_ar,   g_max_ar,               &
                bad_jac_count, degenerate_count, high_skew_count, dup_count)
        end do

        if (n_u == 2) close(u_log)
    end subroutine check_mesh_quality

    ! ------------------------------------------------------------------
    ! Loop over all elements, extract corner nodes, compute quality.
    ! ------------------------------------------------------------------
    subroutine compute_elem_quality(mesh, &
            g_min_jac, g_avg_jac, g_max_jac, total_meas,     &
            g_min_skew, g_avg_skew, g_max_skew,               &
            g_min_ar,   g_avg_ar,   g_max_ar,                 &
            bad_jac, degenerate, high_skew)
        class(t_mesh), intent(in) :: mesh
        real(dp), intent(out) :: g_min_jac, g_avg_jac, g_max_jac, total_meas
        real(dp), intent(out) :: g_min_skew, g_avg_skew, g_max_skew
        real(dp), intent(out) :: g_min_ar,   g_avg_ar,   g_max_ar
        integer,  intent(out) :: bad_jac, degenerate, high_skew

        integer  :: ee, p, c(8)
        real(dp) :: xc(8,3), detJ, skew, ar
        real(dp) :: inv_ne

        p = mesh%order

        g_min_jac  = huge(1.0_dp);  g_max_jac  = -huge(1.0_dp); g_avg_jac  = 0.0_dp
        g_min_skew = huge(1.0_dp);  g_max_skew = 0.0_dp;         g_avg_skew = 0.0_dp
        g_min_ar   = huge(1.0_dp);  g_max_ar   = 0.0_dp;         g_avg_ar   = 0.0_dp
        total_meas = 0.0_dp
        bad_jac = 0; degenerate = 0; high_skew = 0

        ! Corner local indices for order-p elements (x-fastest, 1-based)
        if (mesh%dim == 2) then
            c(1) = 1;         c(2) = p+1
            c(3) = p*(p+1)+1; c(4) = (p+1)**2
        else
            c(1) = 1;                         c(2) = p+1
            c(3) = p*(p+1)+1;                 c(4) = (p+1)**2
            c(5) = p*(p+1)**2+1;              c(6) = p*(p+1)**2+p+1
            c(7) = p*(p+1)**2+p*(p+1)+1;      c(8) = (p+1)**3
        end if

        do ee = 1, mesh%n_elems
            if (mesh%dim == 2) then
                xc(1,:) = mesh%nodes(mesh%elems(ee, c(1)), :)
                xc(2,:) = mesh%nodes(mesh%elems(ee, c(2)), :)
                xc(3,:) = mesh%nodes(mesh%elems(ee, c(3)), :)
                xc(4,:) = mesh%nodes(mesh%elems(ee, c(4)), :)
                call quad_center_quality(xc(1:4, :), detJ, skew, ar)
                total_meas = total_meas + max(detJ, 0.0_dp) * 4.0_dp
            else
                xc(1,:) = mesh%nodes(mesh%elems(ee, c(1)), :)
                xc(2,:) = mesh%nodes(mesh%elems(ee, c(2)), :)
                xc(3,:) = mesh%nodes(mesh%elems(ee, c(3)), :)
                xc(4,:) = mesh%nodes(mesh%elems(ee, c(4)), :)
                xc(5,:) = mesh%nodes(mesh%elems(ee, c(5)), :)
                xc(6,:) = mesh%nodes(mesh%elems(ee, c(6)), :)
                xc(7,:) = mesh%nodes(mesh%elems(ee, c(7)), :)
                xc(8,:) = mesh%nodes(mesh%elems(ee, c(8)), :)
                call hex_center_quality(xc, detJ, skew, ar)
                total_meas = total_meas + max(detJ, 0.0_dp) * 8.0_dp
            end if

            g_min_jac  = min(g_min_jac,  detJ)
            g_max_jac  = max(g_max_jac,  detJ)
            g_avg_jac  = g_avg_jac  + detJ
            g_min_skew = min(g_min_skew, skew)
            g_max_skew = max(g_max_skew, skew)
            g_avg_skew = g_avg_skew + skew
            g_min_ar   = min(g_min_ar,   ar)
            g_max_ar   = max(g_max_ar,   ar)
            g_avg_ar   = g_avg_ar   + ar

            if (detJ < 0.0_dp)              bad_jac    = bad_jac    + 1
            if (abs(detJ) < dp_EPSILON)     degenerate = degenerate + 1
            if (skew > 0.85_dp)             high_skew  = high_skew  + 1
        end do

        inv_ne = 1.0_dp / real(mesh%n_elems, dp)
        g_avg_jac  = g_avg_jac  * inv_ne
        g_avg_skew = g_avg_skew * inv_ne
        g_avg_ar   = g_avg_ar   * inv_ne
    end subroutine compute_elem_quality

    ! ------------------------------------------------------------------
    ! Bilinear quad quality at element center (xi=eta=0).
    ! Corners in FEM x-fastest order: BL(1), BR(2), TL(3), TR(4).
    ! Jacobian uses bilinear shape gradients: dN/dxi=[-1,1,-1,1]/4,
    !                                          dN/deta=[-1,-1,1,1]/4.
    ! Skewness: equiangular (max|angle - 90|/90).
    ! Aspect ratio: max_edge / min_edge.
    ! ------------------------------------------------------------------
    subroutine quad_center_quality(x, detJ, skew, ar)
        real(dp), intent(in)  :: x(4,3)
        real(dp), intent(out) :: detJ, skew, ar

        real(dp) :: J11, J12, J21, J22
        real(dp) :: e1(2), e2(2), e3(2), e4(2), l(4)
        real(dp) :: ang(4), max_dev
        real(dp), parameter :: PI_2 = 0.5_dp * (4.0_dp * atan(1.0_dp))
        integer :: k

        ! Bilinear Jacobian at center
        J11 = 0.25_dp * (-x(1,1) + x(2,1) - x(3,1) + x(4,1))
        J12 = 0.25_dp * (-x(1,2) + x(2,2) - x(3,2) + x(4,2))
        J21 = 0.25_dp * (-x(1,1) - x(2,1) + x(3,1) + x(4,1))
        J22 = 0.25_dp * (-x(1,2) - x(2,2) + x(3,2) + x(4,2))
        detJ = J11*J22 - J12*J21

        ! Edges: BL→BR, BR→TR, TR→TL, TL→BL
        e1 = x(2,1:2) - x(1,1:2)
        e2 = x(4,1:2) - x(2,1:2)
        e3 = x(3,1:2) - x(4,1:2)
        e4 = x(1,1:2) - x(3,1:2)
        l(1) = max(norm2(e1), dp_EPSILON)
        l(2) = max(norm2(e2), dp_EPSILON)
        l(3) = max(norm2(e3), dp_EPSILON)
        l(4) = max(norm2(e4), dp_EPSILON)

        ar = maxval(l) / minval(l)

        ! Interior angles at each corner (angle between incoming and outgoing edge)
        ang(1) = safe_angle2d(-e4, e1)
        ang(2) = safe_angle2d(-e1, e2)
        ang(3) = safe_angle2d(-e2, e3)
        ang(4) = safe_angle2d(-e3, e4)
        max_dev = 0.0_dp
        do k = 1, 4; max_dev = max(max_dev, abs(ang(k) - PI_2)); end do
        skew = max_dev / PI_2
    end subroutine quad_center_quality

    ! ------------------------------------------------------------------
    ! Trilinear hex quality at element center (xi=eta=zeta=0).
    ! Corners in FEM x-fastest order (see m_basis_fem for layout).
    ! dN/dxi=[-1,1,-1,1,-1,1,-1,1]/8, dN/deta=[-1,-1,1,1,-1,-1,1,1]/8,
    ! dN/dzeta=[-1,-1,-1,-1,1,1,1,1]/8.
    ! Skewness from Jacobian column orthogonality.
    ! ------------------------------------------------------------------
    subroutine hex_center_quality(x, detJ, skew, ar)
        real(dp), intent(in)  :: x(8,3)
        real(dp), intent(out) :: detJ, skew, ar

        real(dp) :: J(3,3), e1(3), e2(3), edges(12)
        real(dp) :: max_dev, ang
        real(dp), parameter :: &
            W8(8,3) = reshape([ &
                -1.0_dp, 1.0_dp,-1.0_dp, 1.0_dp,-1.0_dp, 1.0_dp,-1.0_dp, 1.0_dp, &  ! dNdxi
                -1.0_dp,-1.0_dp, 1.0_dp, 1.0_dp,-1.0_dp,-1.0_dp, 1.0_dp, 1.0_dp, &  ! dNdeta
                -1.0_dp,-1.0_dp,-1.0_dp,-1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp  &  ! dNdzeta
            ], [8,3]) * 0.125_dp
        real(dp), parameter :: PI_2 = 0.5_dp * (4.0_dp * atan(1.0_dp))
        integer :: d

        do d = 1, 3; J(d,:) = matmul(W8(:,d), x); end do

        detJ = J(1,1)*(J(2,2)*J(3,3) - J(2,3)*J(3,2)) &
             - J(1,2)*(J(2,1)*J(3,3) - J(2,3)*J(3,1)) &
             + J(1,3)*(J(2,1)*J(3,2) - J(2,2)*J(3,1))

        ! 12 hex edges (corner indexing: see hex layout in subroutine header)
        edges(1)  = norm2(x(2,:)-x(1,:)); edges(2)  = norm2(x(4,:)-x(2,:))
        edges(3)  = norm2(x(3,:)-x(4,:)); edges(4)  = norm2(x(1,:)-x(3,:))
        edges(5)  = norm2(x(6,:)-x(5,:)); edges(6)  = norm2(x(8,:)-x(6,:))
        edges(7)  = norm2(x(7,:)-x(8,:)); edges(8)  = norm2(x(5,:)-x(7,:))
        edges(9)  = norm2(x(5,:)-x(1,:)); edges(10) = norm2(x(6,:)-x(2,:))
        edges(11) = norm2(x(7,:)-x(3,:)); edges(12) = norm2(x(8,:)-x(4,:))
        ar = maxval(edges) / max(minval(edges), dp_EPSILON)

        ! Skewness: max deviation from orthogonality of J column directions
        max_dev = 0.0_dp
        e1 = J(1,:) / max(norm2(J(1,:)), dp_EPSILON)
        e2 = J(2,:) / max(norm2(J(2,:)), dp_EPSILON)
        ang = acos(max(-1.0_dp, min(1.0_dp, dot_product(e1, e2))))
        max_dev = max(max_dev, abs(ang - PI_2))

        e2 = J(3,:) / max(norm2(J(3,:)), dp_EPSILON)
        ang = acos(max(-1.0_dp, min(1.0_dp, dot_product(e1, e2))))
        max_dev = max(max_dev, abs(ang - PI_2))

        e1 = J(2,:) / max(norm2(J(2,:)), dp_EPSILON)
        ang = acos(max(-1.0_dp, min(1.0_dp, dot_product(e1, e2))))
        max_dev = max(max_dev, abs(ang - PI_2))

        skew = max_dev / PI_2
    end subroutine hex_center_quality

    pure real(dp) function safe_angle2d(a, b)
        real(dp), intent(in) :: a(2), b(2)
        real(dp) :: denom
        real(dp), parameter :: PI_2 = 0.5_dp * (4.0_dp * atan(1.0_dp))
        denom = norm2(a) * norm2(b)
        if (denom < dp_EPSILON) then
            safe_angle2d = PI_2
        else
            safe_angle2d = acos(max(-1.0_dp, min(1.0_dp, dot_product(a, b) / denom)))
        end if
    end function safe_angle2d

    ! ------------------------------------------------------------------
    ! Count coincident nodes using a sorted-x scan.
    ! Tolerance = 1e-8 * bounding-box diagonal.
    ! ------------------------------------------------------------------
    subroutine count_coincident_nodes(mesh, dup_count)
        class(t_mesh), intent(in) :: mesh
        integer,       intent(out) :: dup_count

        real(dp) :: tol, extent
        real(dp), allocatable :: xkey(:)
        integer,  allocatable :: idx(:)
        logical,  allocatable :: flagged(:)
        integer :: i, j, n

        n = mesh%n_nodes
        dup_count = 0
        if (n <= 1) return

        extent = norm2(maxval(mesh%nodes(:, 1:mesh%dim), 1) - &
                       minval(mesh%nodes(:, 1:mesh%dim), 1))
        tol = 1.0e-8_dp * max(extent, dp_EPSILON)

        allocate(xkey(n), idx(n), flagged(n))
        xkey    = mesh%nodes(:, 1)
        idx     = [(i, i=1,n)]
        flagged = .false.
        call qsort_idx(xkey, idx, 1, n)

        do i = 1, n-1
            if (flagged(idx(i))) cycle
            do j = i+1, n
                if (xkey(j) - xkey(i) > tol) exit
                if (flagged(idx(j))) cycle
                if (norm2(mesh%nodes(idx(j), 1:mesh%dim) - &
                          mesh%nodes(idx(i), 1:mesh%dim)) < tol) then
                    flagged(idx(j)) = .true.
                    dup_count = dup_count + 1
                end if
            end do
        end do
        deallocate(xkey, idx, flagged)
    end subroutine count_coincident_nodes

    recursive subroutine qsort_idx(a, idx, lo, hi)
        real(dp), intent(inout) :: a(:)
        integer,  intent(inout) :: idx(:)
        integer,  intent(in)    :: lo, hi
        integer  :: i, j, ti
        real(dp) :: pivot, tr
        if (lo >= hi) return
        pivot = a((lo+hi)/2)
        i = lo; j = hi
        do
            do while (a(i) < pivot); i = i + 1; end do
            do while (a(j) > pivot); j = j - 1; end do
            if (i >= j) exit
            tr = a(i); a(i) = a(j); a(j) = tr
            ti = idx(i); idx(i) = idx(j); idx(j) = ti
            i = i + 1; j = j - 1
        end do
        call qsort_idx(a, idx, lo, j)
        call qsort_idx(a, idx, j+1, hi)
    end subroutine qsort_idx

    ! ------------------------------------------------------------------
    ! Format and write the quality report to unit u.
    ! ------------------------------------------------------------------
    subroutine print_quality_report(mesh, u,                             &
            g_min_jac, g_avg_jac, g_max_jac, total_meas,                &
            g_min_skew, g_avg_skew, g_max_skew,                          &
            g_min_ar,   g_avg_ar,   g_max_ar,                            &
            bad_jac, degenerate, high_skew, dup_count)
        class(t_mesh), intent(in) :: mesh
        integer,       intent(in) :: u
        real(dp),      intent(in) :: g_min_jac, g_avg_jac, g_max_jac, total_meas
        real(dp),      intent(in) :: g_min_skew, g_avg_skew, g_max_skew
        real(dp),      intent(in) :: g_min_ar, g_avg_ar, g_max_ar
        integer,       intent(in) :: bad_jac, degenerate, high_skew, dup_count

        integer :: s, m, b, n_surf, max_mat, max_bc, n_bcs
        integer, allocatable :: mat_cnt(:), bc_cnt(:)
        character(len=12) :: elem_str, iga_note

        n_surf = size(mesh%surfaces)

        select case (mesh%dim)
        case (2); elem_str = "Quad-4 (2D)"
        case (3); elem_str = "Hex-8  (3D)"
        case default; elem_str = "Unknown"
        end select

        ! IGA note: quality based on control polygon approximation
        select type (mesh)
        type is (t_mesh)
            iga_note = ""
        class default
            iga_note = " [approx]"
        end select

        write(u,*)
        write(u,'(A)') BOX_H
        write(u,'(A)') " |                         MESH QUALITY REPORT                           |"
        write(u,'(A)') BOX_H

        write(u,'(A)') " |  [1] GEOMETRY & TOPOLOGY INVENTORY                                    |"
        write(u,'(A, I15, T74, A)')   " |      - Total Nodes        : ", mesh%n_nodes,  "|"
        write(u,'(A, I15, T74, A)')   " |      - Total Elements     : ", mesh%n_elems,  "|"
        write(u,'(A, A15, T74, A)')   " |      - Element Type       : ", trim(elem_str),"|"
        write(u,'(A, I15, T74, A)')   " |      - Polynomial Order   : ", mesh%order,    "|"
        write(u,'(A, I15, T74, A)')   " |      - Energy Groups      : ", mesh%n_groups, "|"
        write(u,'(A, I15, T74, A)')   " |      - Boundary Surfaces  : ", n_surf,        "|"
        write(u,'(A, I15, T74, A)')   " |      - Coincident Nodes   : ", dup_count,     "|"

        write(u,'(A)') BOX_S
        write(u,'(A, A, A)') " |  [2] ELEMENT QUALITY METRICS", trim(iga_note), &
            repeat(' ', 42 - len_trim(iga_note)) // "|"
        write(u,'(A, ES15.4, T74, A)') " |      - Min Jacobian Det   : ", g_min_jac,   "|"
        write(u,'(A, ES15.4, T74, A)') " |      - Avg Jacobian Det   : ", g_avg_jac,   "|"
        write(u,'(A, ES15.4, T74, A)') " |      - Max Jacobian Det   : ", g_max_jac,   "|"
        if (mesh%dim == 2) then
            write(u,'(A, ES15.4, T74, A)') " |      - Total Area         : ", total_meas, "|"
        else
            write(u,'(A, ES15.4, T74, A)') " |      - Total Volume       : ", total_meas, "|"
        end if
        write(u,'(A, F15.4, T74, A)') " |      - Min Skewness       : ", g_min_skew,  "|"
        write(u,'(A, F15.4, T74, A)') " |      - Avg Skewness       : ", g_avg_skew,  "|"
        write(u,'(A, F15.4, T74, A)') " |      - Max Skewness       : ", g_max_skew,  "|"
        write(u,'(A, F15.4, T74, A)') " |      - Min Aspect Ratio   : ", g_min_ar,    "|"
        write(u,'(A, F15.4, T74, A)') " |      - Avg Aspect Ratio   : ", g_avg_ar,    "|"
        write(u,'(A, F15.4, T74, A)') " |      - Max Aspect Ratio   : ", g_max_ar,    "|"
        write(u,'(A, I15, T74, A)')   " |      - Inverted (J<0)     : ", bad_jac,     "|"
        write(u,'(A, I15, T74, A)')   " |      - Degenerate (|J|~0) : ", degenerate,  "|"
        write(u,'(A, I15, T74, A)')   " |      - High Skew  (>0.85) : ", high_skew,   "|"

        write(u,'(A)') BOX_S
        write(u,'(A)') " |  [3] MATERIAL ID INVENTORY                                            |"
        if (mesh%n_elems > 0) then
            max_mat = maxval(mesh%material_ids)
            if (max_mat > 0) then
                allocate(mat_cnt(max_mat)); mat_cnt = 0
                do m = 1, mesh%n_elems
                    mat_cnt(mesh%material_ids(m)) = mat_cnt(mesh%material_ids(m)) + 1
                end do
                do m = 1, max_mat
                    if (mat_cnt(m) > 0) &
                        write(u,'(A, I4, A, I10, T74, A)') &
                            " |      - Material ID ", m, "  : ", mat_cnt(m), "|"
                end do
                deallocate(mat_cnt)
            end if
        end if

        write(u,'(A)') BOX_S
        write(u,'(A)') " |  [4] BOUNDARY ID INVENTORY                                            |"
        if (n_surf > 0) then
            max_bc = 0
            do s = 1, n_surf
                if (mesh%surfaces(s)%bc_id > max_bc) max_bc = mesh%surfaces(s)%bc_id
            end do
            if (max_bc > 0) then
                allocate(bc_cnt(max_bc)); bc_cnt = 0
                do s = 1, n_surf
                    b = mesh%surfaces(s)%bc_id
                    if (b > 0) bc_cnt(b) = bc_cnt(b) + 1
                end do
                n_bcs = count(bc_cnt > 0)
                if (n_bcs > 0) then
                    do b = 1, max_bc
                        if (bc_cnt(b) > 0) &
                            write(u,'(A, I4, A, I10, T74, A)') &
                                " |      - BC ID       ", b, "  : ", bc_cnt(b), "|"
                    end do
                else
                    write(u,'(A)') " |      (all surfaces have bc_id <= 0)                                   |"
                end if
                deallocate(bc_cnt)
            else
                write(u,'(A)') " |      (no named boundary surfaces)                                     |"
            end if
        else
            write(u,'(A)') " |      (no boundary surfaces registered)                                |"
        end if

        write(u,'(A)') BOX_S
        write(u,'(A)') " |  [5] HEALTH SUMMARY                                                   |"
        if (degenerate == 0 .and. bad_jac == 0 .and. dup_count == 0 .and. high_skew == 0) then
            write(u,'(A)') " |      >> STATUS: [ PASS ]  No quality issues detected.                 |"
        else
            write(u,'(A)') " |      >> STATUS: [ WARNING ]  Mesh issues detected:                    |"
            if (degenerate > 0) &
                write(u,'(A, I8, T74, A)') " |         -- DEGENERATE ELEMENTS       : ", degenerate, "|"
            if (bad_jac > 0) &
                write(u,'(A, I8, T74, A)') " |         -- INVERTED ELEMENTS (J<0)   : ", bad_jac,    "|"
            if (high_skew > 0) &
                write(u,'(A, I8, T74, A)') " |         -- HIGH SKEWNESS (>0.85)     : ", high_skew,  "|"
            if (dup_count > 0) &
                write(u,'(A, I8, T74, A)') " |         -- COINCIDENT NODES          : ", dup_count,  "|"
        end if
        write(u,'(A)') BOX_H
        write(u,*)
    end subroutine print_quality_report

end module m_utilities
