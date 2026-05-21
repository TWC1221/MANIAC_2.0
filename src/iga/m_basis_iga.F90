! NURBS/B-spline basis functions for IGA (2D and 3D).
! Basis functions are evaluated on-the-fly per element (not precomputed).
! The knot-span parametric coordinate is mapped from reference [-1,1]
! to physical parametric space [u_min, u_max] inside GetMapping.
!
! Public:
!   InitialiseBasis  -- set FE from mesh, build face_node_map
!   GetMapping2D     -- NURBS Jacobian + gradients for a 2D patch element
!   GetMapping3D     -- NURBS Jacobian + gradients for a 3D patch element
!   EvalNURBS2D      -- rational B-spline basis (2D element)
!   EvalNURBS3D      -- rational B-spline basis (3D element)
!   EvalNURBS1D      -- rational B-spline basis (1D edge, for 2D boundaries)
!   FindSpan         -- knot-span binary search (NURBS Book Algorithm A2.1)
!   DersBasisFuns    -- B-spline basis + 1st derivative (Algorithm A2.3)
module m_basis_iga
    use m_constants
    use m_types
    use m_types_iga
    implicit none
    private
    public :: InitialiseNurbsBasis
    public :: GetMapping2D, GetMapping3D
    public :: EvalNURBS1D, EvalNURBS2D, EvalNURBS3D
    public :: DersNurbsBasis

    interface EvalNURBS2D
        module procedure EvalNURBS2D_vol
        module procedure EvalNURBS2D_surf
    end interface EvalNURBS2D

contains

    ! ------------------------------------------------------------------
    ! Populate FE from the mesh order and build the face node map.
    ! The face_node_map gives the local basis index of the n-th node
    ! on each face, in the same x-fastest tensor ordering used by elems.
    ! ------------------------------------------------------------------
    subroutine InitialiseNurbsBasis(FE, mesh)
        type(t_finite_iga), intent(inout) :: FE
        type(t_mesh_iga),   intent(in)    :: mesh
        integer :: ii, jj, kk, p, q, r, n

        FE%dim    = mesh%dim
        FE%order  = mesh%order
        p = mesh%order; q = p; r = p
        FE%p_order = p; FE%q_order = q
        FE%r_order = merge(r, 0, mesh%dim == 3)

        if (mesh%dim == 3) then
            FE%n_basis        = (p+1) * (q+1) * (r+1)
            FE%n_nodes_per_face = (p+1) * (q+1)
            if (allocated(FE%face_node_map)) deallocate(FE%face_node_map)
            allocate(FE%face_node_map(FE%n_nodes_per_face, 6))

            ! 3D hex face map (x-fastest: n = i + (j-1)*(p+1) + (k-1)*(p+1)^2)
            ! Face 1: z=-1 (k=1),   Face 2: z=+1 (k=r+1)
            ! Face 3: y=-1 (j=1),   Face 4: y=+1 (j=q+1)
            ! Face 5: x=-1 (i=1),   Face 6: x=+1 (i=p+1)
            FE%face_node_map(:,1) = [(ii, ii=1, FE%n_nodes_per_face)]
            FE%face_node_map(:,2) = FE%face_node_map(:,1) + r*(p+1)*(q+1)

            n = 0
            do kk = 0, r
                do ii = 1, p+1
                    n = n + 1
                    FE%face_node_map(n, 3) = kk*(p+1)*(q+1) + ii
                    FE%face_node_map(n, 4) = kk*(p+1)*(q+1) + q*(p+1) + ii
                end do
            end do

            n = 0
            do kk = 0, r
                do jj = 0, q
                    n = n + 1
                    FE%face_node_map(n, 5) = kk*(p+1)*(q+1) + jj*(p+1) + 1
                    FE%face_node_map(n, 6) = kk*(p+1)*(q+1) + jj*(p+1) + (p+1)
                end do
            end do

        else  ! dim == 2
            FE%n_basis          = (p+1) * (q+1)
            FE%n_nodes_per_face = p + 1
            if (allocated(FE%face_node_map)) deallocate(FE%face_node_map)
            allocate(FE%face_node_map(FE%n_nodes_per_face, 4))

            ! 2D quad face map (x-fastest: n = i + (j-1)*(p+1))
            ! Face 1: y=-1 (j=1),   Face 2: x=+1 (i=p+1)
            ! Face 3: y=+1 (j=q+1), Face 4: x=-1 (i=1)
            do ii = 1, p+1
                FE%face_node_map(ii,1) = ii
                FE%face_node_map(ii,2) = (ii-1)*(p+1) + (p+1)
                FE%face_node_map(ii,3) = (q+1)*(p+1) - (ii-1)
                FE%face_node_map(ii,4) = (p+1)*(q+1-ii) + 1
            end do
        end if
    end subroutine InitialiseNurbsBasis

    ! ------------------------------------------------------------------
    ! 2D NURBS Jacobian mapping for a quad knot-span element.
    ! Maps reference coords [xi_ref, eta_ref] in [-1,1]^2 through the
    ! parametric map to physical space, returns NURBS basis R_basis and
    ! physical-space gradients dN_dx, dN_dy.
    !
    !   ee          -- element index (for span info and patch lookup)
    !   mesh        -- full IGA mesh
    !   q           -- quadrature point index (for Quad%xi/eta)
    !   Quad        -- 2D Gauss quadrature rule
    !   u1,u2,v1,v2 -- knot span bounds in parametric space
    !   elem_coords -- (n_basis, 2) physical coordinates of CPs
    ! ------------------------------------------------------------------
    subroutine GetMapping2D(FE, ee, mesh, q, Quad, u1, u2, v1, v2, &
                            elem_coords, dN_dx, dN_dy, detJ, R_basis, &
                            R_mat, xi_custom, eta_custom, J_out)
        type(t_finite_iga), intent(in)  :: FE
        integer,            intent(in)  :: ee, q
        type(t_mesh_iga),   intent(in)  :: mesh
        type(t_quadrature), intent(in)  :: Quad
        real(dp),           intent(in)  :: u1, u2, v1, v2
        real(dp),           intent(in)  :: elem_coords(:,:)
        real(dp),           intent(out) :: dN_dx(:), dN_dy(:), detJ, R_basis(:)
        real(dp), optional, intent(out) :: R_mat(:,:)
        real(dp), optional, intent(in)  :: xi_custom, eta_custom
        real(dp), optional, intent(out) :: J_out(2,2)

        real(dp) :: J(2,2), invJ(2,2)
        real(dp) :: dRdXiEta(2, size(R_basis))
        real(dp) :: dR_dxi(size(R_basis)), dR_deta(size(R_basis))
        real(dp) :: xi, eta, xi_ref, eta_ref, detJ_param, detJ_raw

        if (present(xi_custom) .and. present(eta_custom)) then
            xi_ref = xi_custom; eta_ref = eta_custom
        else
            xi_ref = Quad%xi(q); eta_ref = Quad%eta(q)
        end if

        xi         = 0.5_dp * ((u2 - u1)*xi_ref  + (u2 + u1))
        eta        = 0.5_dp * ((v2 - v1)*eta_ref + (v2 + v1))
        detJ_param = 0.25_dp * (u2 - u1) * (v2 - v1)

        call EvalNURBS2D(FE, ee, mesh, xi, eta, R_basis, dR_dxi, dR_deta)

        dRdXiEta(1,:) = dR_dxi
        dRdXiEta(2,:) = dR_deta
        J = matmul(dRdXiEta, elem_coords)
        if (present(J_out)) J_out = J

        detJ_raw = J(1,1)*J(2,2) - J(1,2)*J(2,1)
        detJ     = detJ_raw * detJ_param

        if (detJ < 0.0_dp) &
            write(*,*) " WARNING: inverted 2D Jacobian at elem ", ee, " q ", q, " detJ=", detJ
        if (abs(detJ_raw) < dp_EPSILON) detJ_raw = sign(dp_EPSILON, detJ_raw)

        invJ(1,1) =  J(2,2) / detJ_raw
        invJ(1,2) = -J(1,2) / detJ_raw
        invJ(2,1) = -J(2,1) / detJ_raw
        invJ(2,2) =  J(1,1) / detJ_raw

        dN_dx = invJ(1,1)*dR_dxi + invJ(1,2)*dR_deta
        dN_dy = invJ(2,1)*dR_dxi + invJ(2,2)*dR_deta

        if (present(R_mat)) &
            R_mat = spread(R_basis, dim=2, ncopies=FE%n_basis) * &
                    spread(R_basis, dim=1, ncopies=FE%n_basis)
    end subroutine GetMapping2D

    ! ------------------------------------------------------------------
    ! 3D NURBS Jacobian mapping for a hex knot-span element.
    ! ------------------------------------------------------------------
    subroutine GetMapping3D(FE, ee, mesh, q, Quad, u1, u2, v1, v2, w1, w2, &
                            elem_coords, dN_dx, dN_dy, dN_dz, detJ, R_basis, &
                            xi_custom, eta_custom, zeta_custom, J_out, R_out)
        type(t_finite_iga), intent(in)  :: FE
        integer,            intent(in)  :: ee, q
        type(t_mesh_iga),   intent(in)  :: mesh
        type(t_quadrature), intent(in)  :: Quad
        real(dp),           intent(in)  :: u1, u2, v1, v2, w1, w2
        real(dp),           intent(in)  :: elem_coords(:,:)
        real(dp),           intent(out) :: dN_dx(:), dN_dy(:), dN_dz(:), detJ, R_basis(:)
        real(dp), optional, intent(in)  :: xi_custom, eta_custom, zeta_custom
        real(dp), optional, intent(out) :: J_out(3,3), R_out(3)

        real(dp) :: J(3,3), invJ(3,3)
        real(dp) :: dRdXiEtaZeta(3, size(R_basis))
        real(dp) :: dR_dxi(size(R_basis)), dR_deta(size(R_basis)), dR_dzeta(size(R_basis))
        real(dp) :: xi, eta, zeta, xi_ref, eta_ref, zeta_ref
        real(dp) :: detJ_param, detJ_raw

        if (present(xi_custom) .and. present(eta_custom) .and. present(zeta_custom)) then
            xi_ref = xi_custom; eta_ref = eta_custom; zeta_ref = zeta_custom
        else
            xi_ref = Quad%xi(q); eta_ref = Quad%eta(q); zeta_ref = Quad%zeta(q)
        end if

        xi         = 0.5_dp * ((u2 - u1)*xi_ref   + (u2 + u1))
        eta        = 0.5_dp * ((v2 - v1)*eta_ref  + (v2 + v1))
        zeta       = 0.5_dp * ((w2 - w1)*zeta_ref + (w2 + w1))
        detJ_param = 0.125_dp * (u2 - u1) * (v2 - v1) * (w2 - w1)

        call EvalNURBS3D(FE, ee, mesh, xi, eta, zeta, R_basis, dR_dxi, dR_deta, dR_dzeta)

        dRdXiEtaZeta(1,:) = dR_dxi
        dRdXiEtaZeta(2,:) = dR_deta
        dRdXiEtaZeta(3,:) = dR_dzeta
        J = matmul(dRdXiEtaZeta, elem_coords)
        if (present(J_out)) J_out = J

        if (present(R_out)) then
            R_out(1) = dot_product(R_basis, elem_coords(:,1))
            R_out(2) = dot_product(R_basis, elem_coords(:,2))
            R_out(3) = dot_product(R_basis, elem_coords(:,3))
        end if

        detJ_raw = J(1,1)*(J(2,2)*J(3,3) - J(2,3)*J(3,2)) &
                 - J(1,2)*(J(2,1)*J(3,3) - J(2,3)*J(3,1)) &
                 + J(1,3)*(J(2,1)*J(3,2) - J(2,2)*J(3,1))
        detJ = detJ_raw * detJ_param

        if (detJ < 0.0_dp) &
            write(*,*) " WARNING: inverted 3D Jacobian at elem ", ee, " q ", q, " detJ=", detJ
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

        dN_dx = invJ(1,1)*dR_dxi + invJ(1,2)*dR_deta + invJ(1,3)*dR_dzeta
        dN_dy = invJ(2,1)*dR_dxi + invJ(2,2)*dR_deta + invJ(2,3)*dR_dzeta
        dN_dz = invJ(3,1)*dR_dxi + invJ(3,2)*dR_deta + invJ(3,3)*dR_dzeta
    end subroutine GetMapping3D

    ! ------------------------------------------------------------------
    ! Evaluate 2D NURBS basis and parametric derivatives at (xi, eta).
    ! Volume element variant: uses precomputed span indices from mesh.
    ! ------------------------------------------------------------------
    subroutine EvalNURBS2D_vol(FE, ee, mesh, xi, eta, R_basis, dR_dxi, dR_deta)
        type(t_finite_iga), intent(in)  :: FE
        integer,            intent(in)  :: ee
        type(t_mesh_iga),   intent(in)  :: mesh
        real(dp),           intent(in)  :: xi, eta
        real(dp),           intent(out) :: R_basis(:), dR_dxi(:), dR_deta(:)

        integer :: p, q, span_xi, span_eta, ii, jj, idx, p_idx
        real(dp) :: dN_xi(2, FE%p_order+1), dN_eta(2, FE%q_order+1)
        real(dp) :: W, dW_dxi, dW_deta, w_ij, invW, invW2

        p = FE%p_order; q = FE%q_order
        p_idx    = mesh%elem_patch_id(ee)
        span_xi  = mesh%elem_span_indices(1, ee)
        span_eta = mesh%elem_span_indices(2, ee)

        R_basis = 0.0_dp; dR_dxi = 0.0_dp; dR_deta = 0.0_dp

        call DersNurbsBasis(span_xi,  xi,  p, 1, mesh%patches(p_idx)%knots_xi,  dN_xi)
        call DersNurbsBasis(span_eta, eta, q, 1, mesh%patches(p_idx)%knots_eta, dN_eta)

        W = 0.0_dp; dW_dxi = 0.0_dp; dW_deta = 0.0_dp
        idx = 0
        do jj = 1, q+1
            do ii = 1, p+1
                idx = idx + 1
                w_ij = mesh%weights(mesh%elems(ee, idx))
                W       = W       + dN_xi(1,ii) * dN_eta(1,jj) * w_ij
                dW_dxi  = dW_dxi  + dN_xi(2,ii) * dN_eta(1,jj) * w_ij
                dW_deta = dW_deta + dN_xi(1,ii) * dN_eta(2,jj) * w_ij
            end do
        end do

        invW = 1.0_dp / W; invW2 = invW * invW
        idx = 0
        do jj = 1, q+1
            do ii = 1, p+1
                idx = idx + 1
                w_ij = mesh%weights(mesh%elems(ee, idx))
                R_basis(idx) = (dN_xi(1,ii) * dN_eta(1,jj) * w_ij) * invW
                dR_dxi(idx)  = ((dN_xi(2,ii) * dN_eta(1,jj) * w_ij) * W  &
                              - (dN_xi(1,ii) * dN_eta(1,jj) * w_ij) * dW_dxi)  * invW2
                dR_deta(idx) = ((dN_xi(1,ii) * dN_eta(2,jj) * w_ij) * W  &
                              - (dN_xi(1,ii) * dN_eta(1,jj) * w_ij) * dW_deta) * invW2
            end do
        end do
    end subroutine EvalNURBS2D_vol

    ! ------------------------------------------------------------------
    ! Surface element variant: uses precomputed span indices from surf.
    ! mesh_weights(:) is mesh%weights from the parent t_mesh_iga.
    ! ------------------------------------------------------------------
    subroutine EvalNURBS2D_surf(FE, ee, surf, mesh_weights, xi, eta, R_basis, dR_dxi, dR_deta)
        type(t_finite_iga),  intent(in)  :: FE
        integer,             intent(in)  :: ee
        type(t_surface_iga), intent(in)  :: surf
        real(dp),            intent(in)  :: mesh_weights(:)
        real(dp),            intent(in)  :: xi, eta
        real(dp),            intent(out) :: R_basis(:), dR_dxi(:), dR_deta(:)

        integer :: p, q, span_xi, span_eta, ii, jj, idx
        real(dp) :: dN_xi(2, FE%p_order+1), dN_eta(2, FE%q_order+1)
        real(dp) :: W, dW_dxi, dW_deta, w_ij, invW, invW2

        p = FE%p_order; q = FE%q_order
        span_xi  = surf%elem_span_indices(1, ee)
        span_eta = surf%elem_span_indices(2, ee)

        R_basis = 0.0_dp; dR_dxi = 0.0_dp; dR_deta = 0.0_dp

        call DersNurbsBasis(span_xi,  xi,  p, 1, surf%knots_xi,  dN_xi)
        call DersNurbsBasis(span_eta, eta, q, 1, surf%knots_eta, dN_eta)

        W = 0.0_dp; dW_dxi = 0.0_dp; dW_deta = 0.0_dp
        idx = 0
        do jj = 1, q+1
            do ii = 1, p+1
                idx = idx + 1
                w_ij    = mesh_weights(surf%elems(ee, idx))
                W       = W       + dN_xi(1,ii) * dN_eta(1,jj) * w_ij
                dW_dxi  = dW_dxi  + dN_xi(2,ii) * dN_eta(1,jj) * w_ij
                dW_deta = dW_deta + dN_xi(1,ii) * dN_eta(2,jj) * w_ij
            end do
        end do

        invW = 1.0_dp / W; invW2 = invW * invW
        idx = 0
        do jj = 1, q+1
            do ii = 1, p+1
                idx = idx + 1
                w_ij = mesh_weights(surf%elems(ee, idx))
                R_basis(idx) = (dN_xi(1,ii) * dN_eta(1,jj) * w_ij) * invW
                dR_dxi(idx)  = ((dN_xi(2,ii) * dN_eta(1,jj) * w_ij) * W  &
                              - (dN_xi(1,ii) * dN_eta(1,jj) * w_ij) * dW_dxi)  * invW2
                dR_deta(idx) = ((dN_xi(1,ii) * dN_eta(2,jj) * w_ij) * W  &
                              - (dN_xi(1,ii) * dN_eta(1,jj) * w_ij) * dW_deta) * invW2
            end do
        end do
    end subroutine EvalNURBS2D_surf

    ! ------------------------------------------------------------------
    ! Evaluate 3D NURBS basis and parametric derivatives at (xi, eta, zeta).
    ! ------------------------------------------------------------------
    subroutine EvalNURBS3D(FE, ee, mesh, xi, eta, zeta, R_basis, dR_dxi, dR_deta, dR_dzeta)
        type(t_finite_iga), intent(in)  :: FE
        integer,            intent(in)  :: ee
        type(t_mesh_iga),   intent(in)  :: mesh
        real(dp),           intent(in)  :: xi, eta, zeta
        real(dp),           intent(out) :: R_basis(:), dR_dxi(:), dR_deta(:), dR_dzeta(:)

        integer :: p, q, r, span_xi, span_eta, span_zeta, ii, jj, kk, idx, p_idx
        real(dp) :: dN_xi(2,FE%p_order+1), dN_eta(2,FE%q_order+1), dN_zeta(2,FE%r_order+1)
        real(dp) :: W, dW_dxi, dW_deta, dW_dzeta, w_ijk, invW, invW2

        p = FE%p_order; q = FE%q_order; r = FE%r_order
        p_idx     = mesh%elem_patch_id(ee)
        span_xi   = mesh%elem_span_indices(1, ee)
        span_eta  = mesh%elem_span_indices(2, ee)
        span_zeta = mesh%elem_span_indices(3, ee)

        R_basis = 0.0_dp; dR_dxi = 0.0_dp; dR_deta = 0.0_dp; dR_dzeta = 0.0_dp

        call DersNurbsBasis(span_xi,   xi,   p, 1, mesh%patches(p_idx)%knots_xi,   dN_xi)
        call DersNurbsBasis(span_eta,  eta,  q, 1, mesh%patches(p_idx)%knots_eta,  dN_eta)
        call DersNurbsBasis(span_zeta, zeta, r, 1, mesh%patches(p_idx)%knots_zeta, dN_zeta)

        W = 0.0_dp; dW_dxi = 0.0_dp; dW_deta = 0.0_dp; dW_dzeta = 0.0_dp
        idx = 0
        do kk = 1, r+1
            do jj = 1, q+1
                do ii = 1, p+1
                    idx = idx + 1
                    w_ijk    = mesh%weights(mesh%elems(ee, idx))
                    W        = W        + dN_xi(1,ii)*dN_eta(1,jj)*dN_zeta(1,kk)*w_ijk
                    dW_dxi   = dW_dxi   + dN_xi(2,ii)*dN_eta(1,jj)*dN_zeta(1,kk)*w_ijk
                    dW_deta  = dW_deta  + dN_xi(1,ii)*dN_eta(2,jj)*dN_zeta(1,kk)*w_ijk
                    dW_dzeta = dW_dzeta + dN_xi(1,ii)*dN_eta(1,jj)*dN_zeta(2,kk)*w_ijk
                end do
            end do
        end do

        invW = 1.0_dp / W; invW2 = invW * invW
        idx = 0
        do kk = 1, r+1
            do jj = 1, q+1
                do ii = 1, p+1
                    idx = idx + 1
                    w_ijk = mesh%weights(mesh%elems(ee, idx))
                    R_basis(idx)  = (dN_xi(1,ii)*dN_eta(1,jj)*dN_zeta(1,kk)*w_ijk) * invW
                    dR_dxi(idx)   = ((dN_xi(2,ii)*dN_eta(1,jj)*dN_zeta(1,kk)*w_ijk)*W  &
                                   - (dN_xi(1,ii)*dN_eta(1,jj)*dN_zeta(1,kk)*w_ijk)*dW_dxi)  * invW2
                    dR_deta(idx)  = ((dN_xi(1,ii)*dN_eta(2,jj)*dN_zeta(1,kk)*w_ijk)*W  &
                                   - (dN_xi(1,ii)*dN_eta(1,jj)*dN_zeta(1,kk)*w_ijk)*dW_deta) * invW2
                    dR_dzeta(idx) = ((dN_xi(1,ii)*dN_eta(1,jj)*dN_zeta(2,kk)*w_ijk)*W  &
                                   - (dN_xi(1,ii)*dN_eta(1,jj)*dN_zeta(1,kk)*w_ijk)*dW_dzeta)* invW2
                end do
            end do
        end do
    end subroutine EvalNURBS3D

    ! ------------------------------------------------------------------
    ! Evaluate 1D NURBS basis at xi for a standalone edge/surface.
    ! Used for 2D boundary condition integration.
    ! span     -- precomputed knot span index (knots(span) <= xi < knots(span+1))
    ! Output R and dR_dxi have size n_cp; only the p+1 active entries
    ! (at indices cp_start:cp_start+p) are filled — rest remain zero.
    ! ------------------------------------------------------------------
    subroutine EvalNURBS1D(p, span, knots, weights, xi, R, dR_dxi)
        integer,  intent(in)  :: p, span
        real(dp), intent(in)  :: knots(:), weights(:), xi
        real(dp), intent(out) :: R(:), dR_dxi(:)

        integer  :: ii, cp_start
        real(dp) :: dN(2, p+1), W, dW_dxi, w_loc, invW, invW2

        R = 0.0_dp; dR_dxi = 0.0_dp
        call DersNurbsBasis(span, xi, p, 1, knots, dN)
        cp_start = span - p

        W = 0.0_dp; dW_dxi = 0.0_dp
        do ii = 1, p+1
            w_loc   = weights(cp_start + ii - 1)
            W       = W       + dN(1,ii) * w_loc
            dW_dxi  = dW_dxi  + dN(2,ii) * w_loc
        end do

        invW = 1.0_dp / W; invW2 = invW * invW
        do ii = 1, p+1
            w_loc              = weights(cp_start + ii - 1)
            R(cp_start+ii-1)      = (dN(1,ii)*w_loc) * invW
            dR_dxi(cp_start+ii-1) = ((dN(2,ii)*w_loc)*W - (dN(1,ii)*w_loc)*dW_dxi) * invW2
        end do
    end subroutine EvalNURBS1D

    ! ------------------------------------------------------------------
    ! NURBS Book Algorithm A2.3: B-spline basis and 1st derivatives.
    ! ders(1, 1:p+1) = N_i values
    ! ders(2, 1:p+1) = dN_i/du values
    ! ------------------------------------------------------------------
    subroutine DersNurbsBasis(ii, u, p, n, UU, ders)
        integer,  intent(in)  :: ii, p, n
        real(dp), intent(in)  :: u, UU(:)
        real(dp), intent(out) :: ders(n+1, p+1)

        real(dp) :: ndu(p+1,p+1), left(p+1), right(p+1), saved, temp, d, a(2,p+1)
        integer  :: jj, r, k, s1, s2, rk, pk, j1, j2

        ndu(1,1) = 1.0_dp
        do jj = 1, p
            left(jj+1)  = u - UU(ii+1-jj)
            right(jj+1) = UU(ii+jj) - u
            saved = 0.0_dp
            do r = 1, jj
                ndu(jj+1,r) = right(r+1) + left(jj-r+2)
                temp = ndu(r,jj) / ndu(jj+1,r)
                ndu(r,jj+1) = saved + right(r+1)*temp
                saved = left(jj-r+2) * temp
            end do
            ndu(jj+1,jj+1) = saved
        end do
        ders(1,:) = ndu(:,p+1)

        do r = 0, p
            s1 = 0; s2 = 1; a(1,1) = 1.0_dp
            do k = 1, n
                d = 0.0_dp; rk = r-k; pk = p-k
                if (r >= k) then
                    a(s2+1,1) = a(s1+1,1) / ndu(pk+2,rk+1)
                    d = a(s2+1,1) * ndu(rk+1,pk+1)
                end if
                j1 = merge(1, -rk+1, rk >= -1)
                j2 = merge(k-1, p-r, r <= p-k)
                do jj = j1, j2
                    a(s2+1,jj+1) = (a(s1+1,jj+1) - a(s1+1,jj)) / ndu(pk+2,rk+jj+1)
                    d = d + a(s2+1,jj+1) * ndu(rk+jj+1,pk+1)
                end do
                if (r <= pk) then
                    a(s2+1,k+1) = -a(s1+1,k) / ndu(pk+2,r+1)
                    d = d + a(s2+1,k+1) * ndu(r+1,pk+1)
                end if
                ders(k+1,r+1) = d
                jj = s1; s1 = s2; s2 = jj
            end do
        end do

        r = p
        do k = 1, n
            ders(k+1,:) = ders(k+1,:) * r
            r = r * (p-k)
        end do
    end subroutine DersNurbsBasis

end module m_basis_iga
