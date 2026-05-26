! Lagrange basis functions for FEM (2D and 3D).
! Basis is precomputed once on the reference element [-1,1]^d and stored
! in t_basis_fem.  GetMapping routines then index into these arrays by
! quadrature point, avoiding per-element basis evaluation.
!
! Face numbering and face_node_map ordering match m_basis_iga exactly so
! that transport/diffusion solvers can be written against a common convention.
!   2D:  Face 1 y=-1,  Face 2 x=+1,  Face 3 y=+1,  Face 4 x=-1
!   3D:  Face 1 z=-1,  Face 2 z=+1,  Face 3 y=-1,  Face 4 y=+1,  Face 5 x=-1,  Face 6 x=+1
!
! Public:
!   InitialiseBasisFEM   -- set up t_basis_fem and precompute all arrays
!   GetMapping2D_FEM     -- Lagrange Jacobian + physical gradients at precomputed quad point
!   GetMapping3D_FEM     -- same for 3D
!   EvalAtFace2D_FEM     -- Jacobian at arbitrary face (xi,eta) coords — face normals/centroids
!   EvalAtFace3D_FEM     -- same for 3D
!   EvalLagrange1D       -- evaluate 1D basis at arbitrary xi
!   EvalLagrange2D       -- evaluate 2D basis at arbitrary (xi,eta)
!   EvalLagrange3D       -- evaluate 3D basis at arbitrary (xi,eta,zeta)
module m_basis_fem
    use m_constants
    use m_types
    implicit none
    private
    public :: InitialiseLagrangeBasis
    public :: GetMapping2D_FEM, GetMapping3D_FEM
    public :: EvalAtFace2D_FEM, EvalAtFace3D_FEM
    public :: EvalLagrange1D, EvalLagrange2D, EvalLagrange3D

contains

    ! ------------------------------------------------------------------
    ! Initialise t_basis_fem for a mesh of given dim and polynomial order.
    ! Builds equispaced Lagrange nodes, face_node_map, then precomputes
    ! basis and gradient arrays at all volume and face quadrature points.
    ! ------------------------------------------------------------------
    subroutine InitialiseLagrangeBasis(FE, dim, order, QuadVol, QuadFace)
        type(t_basis_fem), intent(inout) :: FE
        integer,            intent(in)    :: dim, order
        type(t_quadrature), intent(in)    :: QuadVol, QuadFace

        integer  :: i, j, k, ii, jj, kk, q, idx, p, n
        real(dp) :: Li, Lj, Lk, dLi, dLj, dLk

        FE%dim   = dim
        FE%order = order
        p        = order
        FE%p_order = p; FE%q_order = p; FE%r_order = merge(p, 0, dim == 3)

        ! Equispaced Lagrange nodes in [-1, 1]: -1, -1+2/p, ..., 1
        if (allocated(FE%node_roots)) deallocate(FE%node_roots)
        allocate(FE%node_roots(p+1))
        if (p == 0) then
            FE%node_roots(1) = 0.0_dp
        else
            do i = 1, p+1
                FE%node_roots(i) = -1.0_dp + (i-1) * 2.0_dp / real(p, dp)
            end do
        end if

        if (dim == 3) then

            FE%n_basis          = (p+1)**3
            FE%n_nodes_per_face = (p+1)**2
            if (allocated(FE%face_node_map)) deallocate(FE%face_node_map)
            allocate(FE%face_node_map(FE%n_nodes_per_face, 6))

            ! Face 1 z=-1 (k=1), Face 2 z=+1 (k=p+1) — vary x,y
            FE%face_node_map(:,1) = [(ii, ii=1, FE%n_nodes_per_face)]
            FE%face_node_map(:,2) = FE%face_node_map(:,1) + p*(p+1)*(p+1)

            ! Face 3 y=-1 (j=1), Face 4 y=+1 (j=p+1) — vary x,z
            n = 0
            do kk = 0, p
                do ii = 1, p+1
                    n = n + 1
                    FE%face_node_map(n, 3) = kk*(p+1)*(p+1) + ii
                    FE%face_node_map(n, 4) = kk*(p+1)*(p+1) + p*(p+1) + ii
                end do
            end do

            ! Face 5 x=-1 (i=1), Face 6 x=+1 (i=p+1) — vary y,z
            n = 0
            do kk = 0, p
                do jj = 0, p
                    n = n + 1
                    FE%face_node_map(n, 5) = kk*(p+1)*(p+1) + jj*(p+1) + 1
                    FE%face_node_map(n, 6) = kk*(p+1)*(p+1) + jj*(p+1) + (p+1)
                end do
            end do

            ! Volume basis (x-fastest: idx = i + (j-1)*(p+1) + (k-1)*(p+1)^2)
            if (allocated(FE%basis_at_quad))  deallocate(FE%basis_at_quad)
            if (allocated(FE%dbasis_dxi))     deallocate(FE%dbasis_dxi)
            if (allocated(FE%dbasis_deta))    deallocate(FE%dbasis_deta)
            if (allocated(FE%dbasis_dzeta))   deallocate(FE%dbasis_dzeta)
            allocate(FE%basis_at_quad (QuadVol%n_points, FE%n_basis))
            allocate(FE%dbasis_dxi    (QuadVol%n_points, FE%n_basis))
            allocate(FE%dbasis_deta   (QuadVol%n_points, FE%n_basis))
            allocate(FE%dbasis_dzeta  (QuadVol%n_points, FE%n_basis))

            do q = 1, QuadVol%n_points
                idx = 0
                do k = 1, p+1
                    Lk  = Lagrange1D(FE, k, QuadVol%zeta(q))
                    dLk = dLagrange1D(FE, k, QuadVol%zeta(q))
                    do j = 1, p+1
                        Lj  = Lagrange1D(FE, j, QuadVol%eta(q))
                        dLj = dLagrange1D(FE, j, QuadVol%eta(q))
                        do i = 1, p+1
                            Li  = Lagrange1D(FE, i, QuadVol%xi(q))
                            dLi = dLagrange1D(FE, i, QuadVol%xi(q))
                            idx = idx + 1
                            FE%basis_at_quad(q, idx) = Li * Lj * Lk
                            FE%dbasis_dxi   (q, idx) = dLi * Lj * Lk
                            FE%dbasis_deta  (q, idx) = Li * dLj * Lk
                            FE%dbasis_dzeta (q, idx) = Li * Lj * dLk
                        end do
                    end do
                end do
            end do

            ! Face basis (2D Lagrange on the face reference quad [-1,1]^2)
            if (allocated(FE%basis_at_face_quad)) deallocate(FE%basis_at_face_quad)
            if (allocated(FE%dbasis_face_dxi))    deallocate(FE%dbasis_face_dxi)
            if (allocated(FE%dbasis_face_deta))   deallocate(FE%dbasis_face_deta)
            allocate(FE%basis_at_face_quad(QuadFace%n_points, FE%n_nodes_per_face))
            allocate(FE%dbasis_face_dxi   (QuadFace%n_points, FE%n_nodes_per_face))
            allocate(FE%dbasis_face_deta  (QuadFace%n_points, FE%n_nodes_per_face))

            do q = 1, QuadFace%n_points
                idx = 0
                do j = 1, p+1
                    Lj  = Lagrange1D(FE, j, QuadFace%eta(q))
                    dLj = dLagrange1D(FE, j, QuadFace%eta(q))
                    do i = 1, p+1
                        Li  = Lagrange1D(FE, i, QuadFace%xi(q))
                        dLi = dLagrange1D(FE, i, QuadFace%xi(q))
                        idx = idx + 1
                        FE%basis_at_face_quad(q, idx) = Li * Lj
                        FE%dbasis_face_dxi   (q, idx) = dLi * Lj
                        FE%dbasis_face_deta  (q, idx) = Li * dLj
                    end do
                end do
            end do

        else  ! dim == 2

            FE%n_basis          = (p+1)**2
            FE%n_nodes_per_face = p + 1
            if (allocated(FE%face_node_map)) deallocate(FE%face_node_map)
            allocate(FE%face_node_map(FE%n_nodes_per_face, 4))

            ! Face 1 y=-1 (j=1), Face 2 x=+1 (i=p+1), Face 3 y=+1 (j=p+1), Face 4 x=-1 (i=1)
            ! Faces 3 and 4 are reversed for consistent counterclockwise boundary traversal.
            do ii = 1, p+1
                FE%face_node_map(ii, 1) = ii
                FE%face_node_map(ii, 2) = (ii-1)*(p+1) + (p+1)
                FE%face_node_map(ii, 3) = (p+1)*(p+1) - (ii-1)
                FE%face_node_map(ii, 4) = (p+1)*(p+1 - ii) + 1
            end do

            ! Volume basis (x-fastest: idx = i + (j-1)*(p+1))
            if (allocated(FE%basis_at_quad))  deallocate(FE%basis_at_quad)
            if (allocated(FE%dbasis_dxi))     deallocate(FE%dbasis_dxi)
            if (allocated(FE%dbasis_deta))    deallocate(FE%dbasis_deta)
            allocate(FE%basis_at_quad(QuadVol%n_points, FE%n_basis))
            allocate(FE%dbasis_dxi   (QuadVol%n_points, FE%n_basis))
            allocate(FE%dbasis_deta  (QuadVol%n_points, FE%n_basis))

            do q = 1, QuadVol%n_points
                idx = 0
                do j = 1, p+1
                    Lj  = Lagrange1D(FE, j, QuadVol%eta(q))
                    dLj = dLagrange1D(FE, j, QuadVol%eta(q))
                    do i = 1, p+1
                        Li  = Lagrange1D(FE, i, QuadVol%xi(q))
                        dLi = dLagrange1D(FE, i, QuadVol%xi(q))
                        idx = idx + 1
                        FE%basis_at_quad(q, idx) = Li * Lj
                        FE%dbasis_dxi   (q, idx) = dLi * Lj
                        FE%dbasis_deta  (q, idx) = Li * dLj
                    end do
                end do
            end do

            ! Face basis (1D Lagrange on edge reference interval [-1,1])
            if (allocated(FE%basis_at_face_quad)) deallocate(FE%basis_at_face_quad)
            if (allocated(FE%dbasis_face_dxi))    deallocate(FE%dbasis_face_dxi)
            allocate(FE%basis_at_face_quad(QuadFace%n_points, FE%n_nodes_per_face))
            allocate(FE%dbasis_face_dxi   (QuadFace%n_points, FE%n_nodes_per_face))

            do q = 1, QuadFace%n_points
                do i = 1, p+1
                    FE%basis_at_face_quad(q, i) = Lagrange1D (FE, i, QuadFace%xi(q))
                    FE%dbasis_face_dxi   (q, i) = dLagrange1D(FE, i, QuadFace%xi(q))
                end do
            end do

        end if
    end subroutine InitialiseLagrangeBasis

    ! ------------------------------------------------------------------
    ! 2D Lagrange Jacobian mapping for a quad element.
    ! Uses precomputed basis arrays — no per-element knot evaluation.
    !
    !   q           -- quadrature point index (indexes into FE%basis_at_quad etc.)
    !   elem_coords -- (n_basis, 2) physical node coordinates
    ! ------------------------------------------------------------------
    subroutine GetMapping2D_FEM(FE, q, elem_coords, dN_dx, dN_dy, detJ, N_basis, N_mat, J_out)
        type(t_basis_fem), intent(in)  :: FE
        integer,            intent(in)  :: q
        real(dp),           intent(in)  :: elem_coords(:,:)
        real(dp),           intent(out) :: dN_dx(:), dN_dy(:), detJ, N_basis(:)
        real(dp), optional, intent(out) :: N_mat(:,:)
        real(dp), optional, intent(out) :: J_out(2,2)

        real(dp) :: J(2,2), invJ(2,2), detJ_raw
        real(dp) :: dNdXiEta(2, FE%n_basis)

        N_basis(:)       = FE%basis_at_quad(q,:)
        dNdXiEta(1,:)    = FE%dbasis_dxi(q,:)
        dNdXiEta(2,:)    = FE%dbasis_deta(q,:)
        J                = matmul(dNdXiEta, elem_coords)
        if (present(J_out)) J_out = J

        detJ_raw = J(1,1)*J(2,2) - J(1,2)*J(2,1)
        detJ     = detJ_raw

        if (detJ < 0.0_dp) write(*,*) " WARNING: inverted 2D FEM Jacobian at q ", q, " detJ=", detJ
        if (abs(detJ_raw) < dp_EPSILON) detJ_raw = sign(dp_EPSILON, detJ_raw)

        invJ(1,1) =  J(2,2) / detJ_raw;  invJ(1,2) = -J(1,2) / detJ_raw
        invJ(2,1) = -J(2,1) / detJ_raw;  invJ(2,2) =  J(1,1) / detJ_raw

        dN_dx = invJ(1,1)*FE%dbasis_dxi(q,:) + invJ(1,2)*FE%dbasis_deta(q,:)
        dN_dy = invJ(2,1)*FE%dbasis_dxi(q,:) + invJ(2,2)*FE%dbasis_deta(q,:)

        if (present(N_mat)) &
            N_mat = spread(N_basis, dim=2, ncopies=FE%n_basis) * &
                    spread(N_basis, dim=1, ncopies=FE%n_basis)
    end subroutine GetMapping2D_FEM

    ! ------------------------------------------------------------------
    ! 3D Lagrange Jacobian mapping for a hex element.
    ! ------------------------------------------------------------------
    subroutine GetMapping3D_FEM(FE, q, elem_coords, dN_dx, dN_dy, dN_dz, detJ, N_basis, J_out)
        type(t_basis_fem), intent(in)  :: FE
        integer,            intent(in)  :: q
        real(dp),           intent(in)  :: elem_coords(:,:)
        real(dp),           intent(out) :: dN_dx(:), dN_dy(:), dN_dz(:), detJ, N_basis(:)
        real(dp), optional, intent(out) :: J_out(3,3)

        real(dp) :: J(3,3), invJ(3,3), detJ_raw
        real(dp) :: dNdXiEtaZeta(3, FE%n_basis)

        N_basis(:)          = FE%basis_at_quad(q,:)
        dNdXiEtaZeta(1,:)   = FE%dbasis_dxi(q,:)
        dNdXiEtaZeta(2,:)   = FE%dbasis_deta(q,:)
        dNdXiEtaZeta(3,:)   = FE%dbasis_dzeta(q,:)
        J                   = matmul(dNdXiEtaZeta, elem_coords)
        if (present(J_out)) J_out = J

        detJ_raw = J(1,1)*(J(2,2)*J(3,3) - J(2,3)*J(3,2)) &
                 - J(1,2)*(J(2,1)*J(3,3) - J(2,3)*J(3,1)) &
                 + J(1,3)*(J(2,1)*J(3,2) - J(2,2)*J(3,1))
        detJ = detJ_raw

        if (detJ < 0.0_dp) write(*,*) " WARNING: inverted 3D FEM Jacobian at q ", q, " detJ=", detJ
        if (abs(detJ_raw) < dp_EPSILON) detJ_raw = sign(dp_EPSILON, detJ_raw)

        invJ(1,1) = (J(2,2)*J(3,3) - J(2,3)*J(3,2)) / detJ_raw
        invJ(1,2) = (J(1,3)*J(3,2) - J(1,2)*J(3,3)) / detJ_raw
        invJ(1,3) = (J(1,2)*J(2,3) - J(1,3)*J(2,2)) / detJ_raw
        invJ(2,1) = (J(2,3)*J(3,1) - J(2,1)*J(3,3)) / detJ_raw
        invJ(2,2) = (J(1,1)*J(3,3) - J(1,3)*J(3,1)) / detJ_raw
        invJ(2,3) = (J(1,3)*J(2,1) - J(1,1)*J(2,3)) / detJ_raw
        invJ(3,1) = (J(2,1)*J(3,2) - J(2,2)*J(3,1)) / detJ_raw
        invJ(3,2) = (J(1,2)*J(3,1) - J(1,1)*J(3,2)) / detJ_raw
        invJ(3,3) = (J(1,1)*J(2,2) - J(1,2)*J(2,1)) / detJ_raw

        dN_dx = invJ(1,1)*FE%dbasis_dxi(q,:) + invJ(1,2)*FE%dbasis_deta(q,:) + invJ(1,3)*FE%dbasis_dzeta(q,:)
        dN_dy = invJ(2,1)*FE%dbasis_dxi(q,:) + invJ(2,2)*FE%dbasis_deta(q,:) + invJ(2,3)*FE%dbasis_dzeta(q,:)
        dN_dz = invJ(3,1)*FE%dbasis_dxi(q,:) + invJ(3,2)*FE%dbasis_deta(q,:) + invJ(3,3)*FE%dbasis_dzeta(q,:)
    end subroutine GetMapping3D_FEM

    ! ------------------------------------------------------------------
    ! Evaluate 1D Lagrange basis at arbitrary xi.
    ! Used for 2D boundary edge integration.
    ! ------------------------------------------------------------------
    subroutine EvalLagrange1D(FE, xi, N, dN_dxi)
        type(t_basis_fem), intent(in)  :: FE
        real(dp),           intent(in)  :: xi
        real(dp),           intent(out) :: N(:), dN_dxi(:)
        integer :: i
        do i = 1, FE%order + 1
            N(i)      = Lagrange1D (FE, i, xi)
            dN_dxi(i) = dLagrange1D(FE, i, xi)
        end do
    end subroutine EvalLagrange1D

    ! ------------------------------------------------------------------
    ! Evaluate 2D Lagrange basis at arbitrary (xi, eta).
    ! Used for postprocessing and 3D boundary face integration.
    ! ------------------------------------------------------------------
    subroutine EvalLagrange2D(FE, xi, eta, N, dN_dxi, dN_deta)
        type(t_basis_fem), intent(in)  :: FE
        real(dp),           intent(in)  :: xi, eta
        real(dp),           intent(out) :: N(:), dN_dxi(:), dN_deta(:)
        integer  :: i, j, idx
        real(dp) :: Li, Lj, dLi, dLj

        idx = 0
        do j = 1, FE%order + 1
            Lj  = Lagrange1D (FE, j, eta)
            dLj = dLagrange1D(FE, j, eta)
            do i = 1, FE%order + 1
                Li  = Lagrange1D (FE, i, xi)
                dLi = dLagrange1D(FE, i, xi)
                idx = idx + 1
                N(idx)       = Li * Lj
                dN_dxi(idx)  = dLi * Lj
                dN_deta(idx) = Li * dLj
            end do
        end do
    end subroutine EvalLagrange2D

    ! ------------------------------------------------------------------
    ! Evaluate 3D Lagrange basis at arbitrary (xi, eta, zeta).
    ! Used for VTK postprocessing (arbitrary reference coordinates).
    ! ------------------------------------------------------------------
    subroutine EvalLagrange3D(FE, xi, eta, zeta, N, dN_dxi, dN_deta, dN_dzeta)
        type(t_basis_fem), intent(in)  :: FE
        real(dp),           intent(in)  :: xi, eta, zeta
        real(dp),           intent(out) :: N(:), dN_dxi(:), dN_deta(:), dN_dzeta(:)
        integer  :: i, j, k, idx
        real(dp) :: Li, Lj, Lk, dLi, dLj, dLk

        idx = 0
        do k = 1, FE%order + 1
            Lk  = Lagrange1D (FE, k, zeta)
            dLk = dLagrange1D(FE, k, zeta)
            do j = 1, FE%order + 1
                Lj  = Lagrange1D (FE, j, eta)
                dLj = dLagrange1D(FE, j, eta)
                do i = 1, FE%order + 1
                    Li  = Lagrange1D (FE, i, xi)
                    dLi = dLagrange1D(FE, i, xi)
                    idx = idx + 1
                    N(idx)        = Li * Lj * Lk
                    dN_dxi(idx)   = dLi * Lj * Lk
                    dN_deta(idx)  = Li * dLj * Lk
                    dN_dzeta(idx) = Li * Lj * dLk
                end do
            end do
        end do
    end subroutine EvalLagrange3D

    ! ------------------------------------------------------------------
    ! i-th 1D Lagrange basis function at xi, using node_roots.
    ! ------------------------------------------------------------------
    pure real(dp) function Lagrange1D(FE, i, xi)
        type(t_basis_fem), intent(in) :: FE
        integer,            intent(in) :: i
        real(dp),           intent(in) :: xi
        integer  :: p
        real(dp) :: val
        val = 1.0_dp
        do p = 1, FE%order + 1
            if (p /= i) val = val * (xi - FE%node_roots(p)) / (FE%node_roots(i) - FE%node_roots(p))
        end do
        Lagrange1D = val
    end function Lagrange1D

    ! ------------------------------------------------------------------
    ! Derivative of the i-th 1D Lagrange basis function at xi.
    ! ------------------------------------------------------------------
    pure real(dp) function dLagrange1D(FE, i, xi)
        type(t_basis_fem), intent(in) :: FE
        integer,            intent(in) :: i
        real(dp),           intent(in) :: xi
        integer  :: p, q
        real(dp) :: sum_val, prod_val
        sum_val = 0.0_dp
        do p = 1, FE%order + 1
            if (p /= i) then
                prod_val = 1.0_dp / (FE%node_roots(i) - FE%node_roots(p))
                do q = 1, FE%order + 1
                    if (q /= i .and. q /= p) &
                        prod_val = prod_val * (xi - FE%node_roots(q)) / (FE%node_roots(i) - FE%node_roots(q))
                end do
                sum_val = sum_val + prod_val
            end if
        end do
        dLagrange1D = sum_val
    end function dLagrange1D

    ! ------------------------------------------------------------------
    ! Evaluate basis + 2D Jacobian at an arbitrary face (xi,eta) point.
    ! Used for face normals and centroids in connectivity_and_normals_fem.
    ! ------------------------------------------------------------------
    subroutine EvalAtFace2D_FEM(FE, xi, eta, coords2d, N, J)
        type(t_basis_fem), intent(in)  :: FE
        real(dp),          intent(in)  :: xi, eta, coords2d(:,:)
        real(dp),          intent(out) :: N(:), J(2,2)
        real(dp) :: dN_dxi(FE%n_basis), dN_deta(FE%n_basis)
        call EvalLagrange2D(FE, xi, eta, N, dN_dxi, dN_deta)
        J(1,:) = matmul(dN_dxi,  coords2d)
        J(2,:) = matmul(dN_deta, coords2d)
    end subroutine EvalAtFace2D_FEM

    ! ------------------------------------------------------------------
    ! Evaluate basis + 3D Jacobian + physical gradients at an arbitrary
    ! face (xi,eta,zeta) point.
    ! ------------------------------------------------------------------
    subroutine EvalAtFace3D_FEM(FE, xi, eta, zeta, coords3d, N, dN_dx, dN_dy, dN_dz, J)
        type(t_basis_fem), intent(in)  :: FE
        real(dp),          intent(in)  :: xi, eta, zeta, coords3d(:,:)
        real(dp),          intent(out) :: N(:), dN_dx(:), dN_dy(:), dN_dz(:), J(3,3)
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
    end subroutine EvalAtFace3D_FEM

end module m_basis_fem
