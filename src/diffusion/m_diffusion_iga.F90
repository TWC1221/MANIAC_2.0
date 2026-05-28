#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>

! IGA diffusion solver: element assembly, BC application, and power iteration.
!
! t_iga_diffusion extends t_solver and owns all solver state (PETSc handles,
! flux arrays) as components, eliminating module-level save variables.
!
! Public:
!   t_iga_diffusion     -- concrete solver type
!   SolveDiffusion      -- high-level entry: assemble + power iteration
!   assemble_petsc_iga  -- PETSc matrices (A, F, S, ProdVec, FixedSrc)
!   apply_bcs_iga       -- Robin (vacuum/albedo) and Dirichlet BCs
module m_diffusion_iga
    use m_constants
    use m_types
    use m_material
    use m_quadrature
    use m_basis_iga
    use m_petsc,           only: setup_ksp, petsc_build_sparsity, petsc_create_diff_mats, &
                                  petsc_assemble_diff_mats, petsc_setup_ksp_group, &
                                  petsc_destroy_diff_state
    use m_power_iteration, only: t_solver, PowerIteration
    use m_output_iga,      only: export_diffusion_vtk_iga
    use petscsys
    use petscvec
    use petscmat
    use petscksp
    implicit none
    public :: t_iga_diffusion, SolveDiffusion
    public :: assemble_petsc_iga, apply_bcs_iga

    ! ------------------------------------------------------------------
    ! Concrete IGA diffusion solver.
    ! All state previously held in module-level save variables lives here.
    ! ------------------------------------------------------------------
    type, extends(t_solver) :: t_iga_diffusion
        type(t_mesh_iga),   pointer :: mesh    => null()
        type(t_basis_iga),  pointer :: FE      => null()
        type(t_material),   pointer :: mats(:) => null()
        KSP, allocatable :: KSPs(:)
        Vec, allocatable :: X_petsc(:)
        Mat, allocatable :: MAT_F(:,:), MAT_S(:,:)
        Vec, allocatable :: FixedSrc(:)
        Vec  :: tmp_b, tmp_x
        logical :: tmp_valid = .false.
        real(dp), allocatable :: prod_dense(:,:)
        integer :: n_groups = 0
        integer :: n_nodes  = 0
    contains
        procedure :: build_source => iga_diff_build_source
        procedure :: do_solve     => iga_diff_do_solve
        procedure :: compute_prod => iga_diff_compute_prod
        procedure :: snapshot     => iga_diff_snapshot
    end type t_iga_diffusion

contains

    ! ------------------------------------------------------------------
    ! Build PETSc multigroup diffusion matrices for an IGA mesh.
    ! ------------------------------------------------------------------
    subroutine assemble_petsc_iga(A_MAT, MAT_F, MAT_S, PROD_VEC, FixedSrc, &
                                   mesh, FE, Quad, mats, n_groups, is_adjoint)
        Mat, allocatable, intent(out) :: A_MAT(:), MAT_F(:,:), MAT_S(:,:)
        Vec, allocatable, intent(out) :: PROD_VEC(:), FixedSrc(:)
        type(t_mesh_iga),   intent(in) :: mesh
        type(t_basis_iga),  intent(in) :: FE
        type(t_quadrature), intent(in) :: Quad
        type(t_material),   intent(in) :: mats(:)
        integer,            intent(in) :: n_groups
        logical,            intent(in) :: is_adjoint

        PetscInt, allocatable :: nnz(:)
        integer :: ee, count

        call petsc_build_sparsity(mesh%n_nodes, mesh%n_elems, FE%n_basis, mesh%elems, nnz)
        call petsc_create_diff_mats(mesh%n_nodes, n_groups, nnz, A_MAT, MAT_F, MAT_S, PROD_VEC, FixedSrc)
        deallocate(nnz)

        count = 0
        write(*,'(A)') " [ IGA MATRIX ] :: Starting element-wise assembly..."
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ee)
        do ee = 1, mesh%n_elems
            call assemble_iga_elem_petsc(ee, mesh, FE, Quad, mats, n_groups, is_adjoint, &
                                         A_MAT, MAT_F, MAT_S, PROD_VEC, FixedSrc, count)
        end do
        !$OMP END PARALLEL DO

        call petsc_assemble_diff_mats(A_MAT, MAT_F, MAT_S, PROD_VEC, FixedSrc, n_groups)
    end subroutine assemble_petsc_iga

    ! ------------------------------------------------------------------
    ! Internal: assemble one IGA element into PETSc matrices.
    ! ------------------------------------------------------------------
    subroutine assemble_iga_elem_petsc(ee, mesh, FE, Quad, mats, n_groups, is_adjoint, &
                                       A_MAT, MAT_F, MAT_S, PROD_VEC, FixedSrc, count)
        integer,            intent(in)    :: ee, n_groups
        type(t_mesh_iga),   intent(in)    :: mesh
        type(t_basis_iga),  intent(in)    :: FE
        type(t_quadrature), intent(in)    :: Quad
        type(t_material),   intent(in)    :: mats(:)
        logical,            intent(in)    :: is_adjoint
        Mat,                intent(inout) :: A_MAT(:), MAT_F(:,:), MAT_S(:,:)
        Vec,                intent(inout) :: PROD_VEC(:), FixedSrc(:)
        integer,            intent(inout) :: count

        integer  :: i, j, q, mat_id, g_to, g_from
        real(dp) :: dN_dx(FE%n_basis), dN_dy(FE%n_basis), dN_dz(FE%n_basis)
        real(dp) :: R_basis(FE%n_basis), detJ, dV
        real(dp) :: ec(FE%n_basis, 3)
        real(dp) :: u1, u2, v1, v2, w1, w2
        real(dp) :: N_i, N_j, stiff, nsf, chi, sgs
        real(dp) :: loc_A(FE%n_basis, FE%n_basis, n_groups)
        real(dp) :: loc_F(FE%n_basis, FE%n_basis, n_groups, n_groups)
        real(dp) :: loc_S(FE%n_basis, FE%n_basis, n_groups, n_groups)
        real(dp) :: loc_Prod(FE%n_basis, n_groups), loc_Src(FE%n_basis, n_groups)
        PetscInt :: idx(FE%n_basis)
        PetscErrorCode :: ierr

        loc_A = 0.0_dp; loc_F = 0.0_dp; loc_S = 0.0_dp
        loc_Prod = 0.0_dp; loc_Src = 0.0_dp

        mat_id = mesh%material_ids(ee)
        idx    = mesh%elems(ee,:) - 1
        do i = 1, FE%n_basis; ec(i,:) = mesh%nodes(mesh%elems(ee,i),:); end do

        u1 = mesh%elem_u_min(ee); u2 = mesh%elem_u_max(ee)
        v1 = mesh%elem_v_min(ee); v2 = mesh%elem_v_max(ee)

        do q = 1, Quad%n_points
            if (mesh%dim == 2) then
                call GetMapping2D(FE, ee, mesh, q, Quad, u1, u2, v1, v2, &
                                  ec(:,1:2), dN_dx, dN_dy, detJ, R_basis)
                dN_dz = 0.0_dp
            else
                w1 = mesh%elem_w_min(ee); w2 = mesh%elem_w_max(ee)
                call GetMapping3D(FE, ee, mesh, q, Quad, u1, u2, v1, v2, w1, w2, &
                                  ec, dN_dx, dN_dy, dN_dz, detJ, R_basis)
            end if
            dV = detJ * Quad%weights(q)

            do g_to = 1, n_groups
                do i = 1, FE%n_basis
                    N_i = R_basis(i)
                    nsf = merge(mats(mat_id)%NuSigF(g_to), mats(mat_id)%Chi(g_to), .not. is_adjoint)
                    loc_Prod(i,g_to) = loc_Prod(i,g_to) + nsf * N_i * dV
                    loc_Src(i,g_to)  = loc_Src(i,g_to)  + mats(mat_id)%Src(g_to) * N_i * dV
                    do j = 1, FE%n_basis
                        N_j   = R_basis(j)
                        stiff = mats(mat_id)%D(g_to) * &
                                (dN_dx(i)*dN_dx(j) + dN_dy(i)*dN_dy(j) + dN_dz(i)*dN_dz(j)) &
                              + mats(mat_id)%SigmaR(g_to) * N_i * N_j
                        loc_A(i,j,g_to) = loc_A(i,j,g_to) + stiff * dV
                        do g_from = 1, n_groups
                            nsf = merge(mats(mat_id)%NuSigF(g_from),mats(mat_id)%Chi(g_from),.not.is_adjoint)
                            chi = merge(mats(mat_id)%Chi(g_to),mats(mat_id)%NuSigF(g_to),.not.is_adjoint)
                            loc_F(i,j,g_to,g_from) = loc_F(i,j,g_to,g_from) + chi*nsf*N_i*N_j*dV
                            if (g_from /= g_to) then
                                sgs = merge(mats(mat_id)%SigmaS(g_from,g_to), &
                                            mats(mat_id)%SigmaS(g_to,g_from), .not.is_adjoint)
                                loc_S(i,j,g_to,g_from) = loc_S(i,j,g_to,g_from) + sgs*N_i*N_j*dV
                            end if
                        end do
                    end do
                end do
            end do
        end do

        !$OMP CRITICAL
        do g_to = 1, n_groups
            call VecSetValues(PROD_VEC(g_to), FE%n_basis, idx, loc_Prod(:,g_to), ADD_VALUES, ierr)
            call VecSetValues(FixedSrc(g_to), FE%n_basis, idx, loc_Src(:,g_to),  ADD_VALUES, ierr)
            call MatSetValues(A_MAT(g_to), FE%n_basis, idx, FE%n_basis, idx, &
                              reshape(loc_A(:,:,g_to), [FE%n_basis**2]), ADD_VALUES, ierr)
            do g_from = 1, n_groups
                call MatSetValues(MAT_F(g_to,g_from), FE%n_basis, idx, FE%n_basis, idx, &
                                  reshape(loc_F(:,:,g_to,g_from), [FE%n_basis**2]), ADD_VALUES, ierr)
                if (g_from /= g_to) &
                    call MatSetValues(MAT_S(g_to,g_from), FE%n_basis, idx, FE%n_basis, idx, &
                                      reshape(loc_S(:,:,g_to,g_from), [FE%n_basis**2]), ADD_VALUES, ierr)
            end do
        end do
        count = count + 1
        if (mod(count, 7500) == 0) write(*,'(A,I0,A,I0,A,F6.2,A)') &
            ">>> Assembled ", count, " / ", mesh%n_elems, " (", &
            real(count,dp)/real(mesh%n_elems,dp)*100.0_dp, "%)"
        !$OMP END CRITICAL
    end subroutine assemble_iga_elem_petsc

    ! ------------------------------------------------------------------
    ! Apply diffusion boundary conditions to a PETSc matrix.
    ! ------------------------------------------------------------------
    subroutine apply_bcs_iga(mesh, QuadBound, bc_cfg, A)
        type(t_mesh_iga),   intent(in)    :: mesh
        type(t_quadrature), intent(in)    :: QuadBound
        type(t_bc_config),  intent(in)    :: bc_cfg
        Mat,                intent(inout) :: A

        integer, allocatable :: bdy_nodes(:)
        logical, allocatable :: mask(:)
        integer :: s, k, node_id
        real(dp) :: eff_alpha
        PetscErrorCode :: ierr

        select case (bc_cfg%bc_type)
        case (BC_DIRICHLET)
            allocate(mask(mesh%n_nodes)); mask = .false.
            do s = 1, size(mesh%iga_surfaces)
                if (mesh%iga_surfaces(s)%bc_id /= bc_cfg%mat_id) cycle
                do k = 1, size(mesh%iga_surfaces(s)%cp_ids)
                    node_id = mesh%iga_surfaces(s)%cp_ids(k)
                    if (node_id > 0) mask(node_id) = .true.
                end do
            end do
            if (any(mask)) then
                bdy_nodes = pack([(k, k=1, mesh%n_nodes)], mask)
                do k = 1, size(bdy_nodes)
                    call MatSetValue(A, bdy_nodes(k)-1, bdy_nodes(k)-1, PENALTY, ADD_VALUES, ierr)
                end do
                deallocate(bdy_nodes)
            end if
            deallocate(mask)

        case (BC_VACUUM, BC_ALBEDO)
            eff_alpha = merge(bc_cfg%value, 0.0_dp, bc_cfg%bc_type == BC_ALBEDO)
            call iga_robin_petsc(mesh, QuadBound, bc_cfg%mat_id, eff_alpha, A)
        end select
    end subroutine apply_bcs_iga

    ! ------------------------------------------------------------------
    ! Robin BC — 1D edge integral (2D mesh) or 2D surface integral (3D).
    ! ------------------------------------------------------------------
    subroutine iga_robin_petsc(mesh, QuadBound, target_id, alpha, A)
        type(t_mesh_iga),   intent(in)    :: mesh
        type(t_quadrature), intent(in)    :: QuadBound
        integer,            intent(in)    :: target_id
        real(dp),           intent(in)    :: alpha
        Mat,                intent(inout) :: A

        integer :: s, gp, i_span, row_b, col_b
        integer :: n_knots_xi, ncp_global
        real(dp) :: beta, u1, u2, xi, det_param, dV
        real(dp) :: dx_du, dy_du
        PetscErrorCode :: ierr

        type(t_basis_iga) :: surf_FE
        integer :: p, q, ncp_local, ee_surf, span_xi, span_eta
        real(dp) :: v1, v2, eta, det_2d
        real(dp) :: dx_dv, dy_dv, dz_dv, dz_du
        real(dp) :: nx, ny, nz

        real(dp), allocatable :: R(:), dR_du(:), dR_dv(:)
        real(dp), allocatable :: surf_w(:)
        PetscInt, allocatable  :: enodes(:)
        PetscScalar, allocatable :: vals(:)

        beta = 0.5_dp * (1.0_dp - alpha) / (1.0_dp + alpha)

        do s = 1, size(mesh%iga_surfaces)
            if (mesh%iga_surfaces(s)%bc_id /= target_id) cycle

            n_knots_xi = size(mesh%iga_surfaces(s)%knots_xi)

            if (mesh%dim == 2) then
                ncp_global = size(mesh%iga_surfaces(s)%cp_ids)

                allocate(enodes(ncp_global), vals(ncp_global*ncp_global))
                allocate(R(ncp_global), dR_du(ncp_global))
                allocate(surf_w(ncp_global))

                surf_w = mesh%weights(mesh%iga_surfaces(s)%cp_ids)
                enodes = mesh%iga_surfaces(s)%cp_ids - 1
                vals   = 0.0_dp

                do i_span = 1, n_knots_xi - 1
                    u1 = mesh%iga_surfaces(s)%knots_xi(i_span)
                    u2 = mesh%iga_surfaces(s)%knots_xi(i_span+1)
                    if (abs(u2-u1) < 1.0e-10_dp) cycle

                    do gp = 1, QuadBound%n_points
                        xi        = 0.5_dp * ((u2-u1)*QuadBound%xi(gp) + (u2+u1))
                        det_param = 0.5_dp * (u2 - u1)

                        R = 0.0_dp; dR_du = 0.0_dp
                        call EvalNURBS1D(mesh%order, i_span, mesh%iga_surfaces(s)%knots_xi, surf_w, xi, R, dR_du)

                        dx_du = dot_product(mesh%nodes(mesh%iga_surfaces(s)%cp_ids, 1), dR_du)
                        dy_du = dot_product(mesh%nodes(mesh%iga_surfaces(s)%cp_ids, 2), dR_du)

                        dV = sqrt(dx_du**2 + dy_du**2) * det_param * QuadBound%weights(gp) * beta

                        do row_b = 1, ncp_global
                            do col_b = 1, ncp_global
                                vals((row_b-1)*ncp_global+col_b) = vals((row_b-1)*ncp_global+col_b) + &
                                    R(row_b) * R(col_b) * dV
                            end do
                        end do
                    end do
                end do

                call MatSetValues(A, ncp_global, enodes, ncp_global, enodes, vals, ADD_VALUES, ierr)
                deallocate(enodes, vals, R, dR_du, surf_w)

            else
                p = mesh%order
                q = mesh%order
                ncp_local = (p + 1) * (q + 1)

                allocate(enodes(ncp_local), vals(ncp_local * ncp_local))
                allocate(R(ncp_local), dR_du(ncp_local), dR_dv(ncp_local))

                surf_FE%p_order = p
                surf_FE%q_order = q

                do ee_surf = 1, mesh%iga_surfaces(s)%n_elements

                    span_xi  = mesh%iga_surfaces(s)%elem_span_indices(1, ee_surf)
                    span_eta = mesh%iga_surfaces(s)%elem_span_indices(2, ee_surf)

                    u1 = mesh%iga_surfaces(s)%knots_xi(span_xi)
                    u2 = mesh%iga_surfaces(s)%knots_xi(span_xi + 1)
                    if (abs(u2 - u1) < 1.0e-10_dp) cycle

                    v1 = mesh%iga_surfaces(s)%knots_eta(span_eta)
                    v2 = mesh%iga_surfaces(s)%knots_eta(span_eta + 1)
                    if (abs(v2 - v1) < 1.0e-10_dp) cycle

                    enodes = mesh%iga_surfaces(s)%elems(ee_surf, 1:ncp_local)
                    vals = 0.0_dp

                    do gp = 1, QuadBound%n_points
                        xi        = 0.5_dp * ((u2 - u1) * QuadBound%xi(gp)  + (u2 + u1))
                        eta       = 0.5_dp * ((v2 - v1) * QuadBound%eta(gp) + (v2 + v1))
                        det_param = 0.25_dp * (u2 - u1) * (v2 - v1)

                        call EvalNURBS2D(surf_FE, ee_surf, mesh%iga_surfaces(s), mesh%weights, xi, eta, R, dR_du, dR_dv)

                        dx_du = dot_product(mesh%nodes(enodes, 1), dR_du)
                        dy_du = dot_product(mesh%nodes(enodes, 2), dR_du)
                        dz_du = dot_product(mesh%nodes(enodes, 3), dR_du)

                        dx_dv = dot_product(mesh%nodes(enodes, 1), dR_dv)
                        dy_dv = dot_product(mesh%nodes(enodes, 2), dR_dv)
                        dz_dv = dot_product(mesh%nodes(enodes, 3), dR_dv)

                        nx = dy_du*dz_dv - dz_du*dy_dv
                        ny = dz_du*dx_dv - dx_du*dz_dv
                        nz = dx_du*dy_dv - dy_du*dx_dv
                        det_2d = sqrt(nx**2 + ny**2 + nz**2)

                        dV = det_2d * det_param * QuadBound%weights(gp) * beta

                        do row_b = 1, ncp_local
                            do col_b = 1, ncp_local
                                vals((row_b-1)*ncp_local + col_b) = vals((row_b-1)*ncp_local + col_b) + &
                                                                    R(row_b) * R(col_b) * dV
                            end do
                        end do
                    end do

                    call MatSetValues(A, ncp_local, enodes - 1, ncp_local, enodes - 1, vals, ADD_VALUES, ierr)
                end do

                deallocate(enodes, vals, R, dR_du, dR_dv)
            end if
        end do
    end subroutine iga_robin_petsc

    ! ================================================================== !
    !  Type-bound procedure implementations for t_iga_diffusion            !
    ! ================================================================== !

    subroutine iga_diff_build_source(self, src, k_eff, is_eigenvalue, is_adjoint)
        class(t_iga_diffusion), intent(inout) :: self
        real(dp), intent(inout) :: src(:,:)
        real(dp), intent(in)    :: k_eff
        logical,  intent(in)    :: is_eigenvalue, is_adjoint
        integer        :: g, gp
        PetscErrorCode :: ierr
        PetscScalar, pointer :: parr(:)
        src = 0.0_dp
        do g = 1, self%n_groups
            if (.not. is_eigenvalue) then
                call VecGetArrayRead(self%FixedSrc(g), parr, ierr)
                src(:,g) = parr(:)
                call VecRestoreArrayRead(self%FixedSrc(g), parr, ierr)
            end if
            do gp = 1, self%n_groups
                call VecGetArray(self%tmp_x, parr, ierr)
                parr(:) = self%scalar_flux(:,gp)
                call VecRestoreArray(self%tmp_x, parr, ierr)
                if (gp /= g) then
                    call MatMult(self%MAT_S(g,gp), self%tmp_x, self%tmp_b, ierr)
                    call VecGetArrayRead(self%tmp_b, parr, ierr)
                    src(:,g) = src(:,g) + parr(:)
                    call VecRestoreArrayRead(self%tmp_b, parr, ierr)
                end if
                if (is_eigenvalue) then
                    call MatMult(self%MAT_F(g,gp), self%tmp_x, self%tmp_b, ierr)
                    call VecGetArrayRead(self%tmp_b, parr, ierr)
                    src(:,g) = src(:,g) + parr(:) / k_eff
                    call VecRestoreArrayRead(self%tmp_b, parr, ierr)
                end if
            end do
        end do
    end subroutine iga_diff_build_source

    subroutine iga_diff_do_solve(self, src)
        class(t_iga_diffusion), intent(inout) :: self
        real(dp), intent(in) :: src(:,:)
        integer        :: g
        PetscErrorCode :: ierr
        PetscScalar, pointer :: parr(:)
        do g = 1, self%n_groups
            call VecGetArray(self%tmp_b, parr, ierr)
            parr(:) = src(:,g)
            call VecRestoreArray(self%tmp_b, parr, ierr)
            call KSPSolve(self%KSPs(g), self%tmp_b, self%X_petsc(g), ierr)
            call VecGetArrayRead(self%X_petsc(g), parr, ierr)
            self%scalar_flux(:,g) = parr(:)
            call VecRestoreArrayRead(self%X_petsc(g), parr, ierr)
        end do
    end subroutine iga_diff_do_solve

    subroutine iga_diff_compute_prod(self, prod, is_adjoint)
        class(t_iga_diffusion), intent(inout) :: self
        real(dp), intent(out) :: prod
        logical,  intent(in)  :: is_adjoint
        integer :: g
        prod = 0.0_dp
        do g = 1, self%n_groups
            prod = prod + dot_product(self%prod_dense(:,g), self%scalar_flux(:,g))
        end do
    end subroutine iga_diff_compute_prod

    subroutine iga_diff_snapshot(self, k_eff, iter)
        class(t_iga_diffusion), intent(inout) :: self
        real(dp), intent(in) :: k_eff
        integer,  intent(in) :: iter
        character(len=32)  :: lbl
        character(len=512) :: stag
        self%snap_count = self%snap_count + 1
        write(lbl,'(A,I4.4)') "snap", self%snap_count
        stag = trim(self%snap_tag) // "_" // trim(lbl)
        call export_diffusion_vtk_iga(trim(self%snap_dir), trim(stag), &
            self%mesh, self%FE, self%scalar_flux, self%n_groups, self%vtk_refine)
    end subroutine iga_diff_snapshot

    ! ==================================================================
    ! High-level entry point.
    ! ==================================================================
    subroutine SolveDiffusion(solver, mesh, FE, Quad, QuadBound, mats,   &
                               solver_type, preconditioner, ref_ids, &
                               max_outer, tol, is_eigenvalue, is_adjoint, &
                               k_eff_out, snap_dir, snap_tag, vtk_ref)
        type(t_iga_diffusion), intent(out)     :: solver
        type(t_mesh_iga),   intent(in), target :: mesh
        type(t_basis_iga),  intent(in), target :: FE
        type(t_quadrature), intent(in)         :: Quad, QuadBound
        type(t_material),   intent(in), target :: mats(:)
        integer,            intent(in)         :: solver_type, preconditioner
        integer,            intent(in)         :: ref_ids(:)
        integer,            intent(in)         :: max_outer
        real(dp),           intent(in)         :: tol
        logical,            intent(in)         :: is_eigenvalue, is_adjoint
        real(dp),              intent(out)     :: k_eff_out
        character(len=*), optional, intent(in) :: snap_dir, snap_tag
        integer,          optional, intent(in) :: vtk_ref

        integer :: g, ss, n_uniq
        type(t_bc_config)    :: bc_vac
        integer, allocatable :: uniq_bc_ids(:)
        PetscErrorCode       :: ierr
        Mat, allocatable :: loc_A(:), loc_MF(:,:), loc_MS(:,:)
        Vec, allocatable :: loc_PV(:), loc_FS(:)

        solver%mesh     => mesh
        solver%FE       => FE
        solver%mats     => mats
        solver%n_groups  = mesh%n_groups
        solver%n_nodes   = mesh%n_nodes
        solver%tmp_valid = .false.

        solver%have_snap = present(snap_dir) .and. present(snap_tag)
        if (solver%have_snap) then
            solver%snap_dir   = trim(snap_dir)
            solver%snap_tag   = trim(snap_tag)
            solver%vtk_refine = merge(vtk_ref, 4, present(vtk_ref))
            solver%snap_count = 0
        end if

        call assemble_petsc_iga(loc_A, loc_MF, loc_MS, loc_PV, loc_FS, &
                                 mesh, FE, Quad, mats, mesh%n_groups, is_adjoint)

        call collect_vacuum_ids(mesh, ref_ids, uniq_bc_ids, n_uniq)
        bc_vac%bc_type = BC_VACUUM; bc_vac%value = 0.0_dp
        do ss = 1, n_uniq
            bc_vac%mat_id = uniq_bc_ids(ss)
            do g = 1, mesh%n_groups
                call apply_bcs_iga(mesh, QuadBound, bc_vac, loc_A(g))
            end do
        end do
        deallocate(uniq_bc_ids)

        allocate(solver%MAT_F(mesh%n_groups, mesh%n_groups), &
                 solver%MAT_S(mesh%n_groups, mesh%n_groups), &
                 solver%FixedSrc(mesh%n_groups))
        solver%MAT_F    = loc_MF
        solver%MAT_S    = loc_MS
        solver%FixedSrc = loc_FS
        deallocate(loc_MF, loc_MS, loc_FS)

        call petsc_setup_ksp_group(loc_A, mesh%n_nodes, mesh%n_groups, &
                                    solver_type, preconditioner, solver%KSPs, solver%X_petsc)
        deallocate(loc_A)

        call extract_prod_dense(loc_PV, mesh%n_groups, mesh%n_nodes, solver%prod_dense)
        do g = 1, mesh%n_groups; call VecDestroy(loc_PV(g), ierr); end do
        deallocate(loc_PV)

        call VecCreateSeq(PETSC_COMM_SELF, mesh%n_nodes, solver%tmp_b, ierr)
        call VecCreateSeq(PETSC_COMM_SELF, mesh%n_nodes, solver%tmp_x, ierr)
        solver%tmp_valid = .true.

        allocate(solver%scalar_flux(mesh%n_nodes, mesh%n_groups), source=1.0_dp)

        call PowerIteration(solver, k_eff_out, max_outer, tol, is_eigenvalue, is_adjoint)

        call petsc_destroy_diff_state(solver%KSPs, solver%X_petsc, solver%MAT_F, solver%MAT_S, &
                                       solver%FixedSrc, solver%tmp_b, solver%tmp_x, solver%tmp_valid)
    end subroutine SolveDiffusion

    subroutine extract_prod_dense(prod_vecs, n_groups, n_nodes, prod_dense_out)
        Vec,     intent(in)  :: prod_vecs(:)
        integer, intent(in)  :: n_groups, n_nodes
        real(dp), allocatable, intent(out) :: prod_dense_out(:,:)
        integer :: g
        PetscErrorCode :: ierr
        PetscScalar, pointer :: parr(:)

        allocate(prod_dense_out(n_nodes, n_groups))
        do g = 1, n_groups
            call VecGetArrayRead(prod_vecs(g), parr, ierr)
            prod_dense_out(:,g) = parr(:)
            call VecRestoreArrayRead(prod_vecs(g), parr, ierr)
        end do
    end subroutine extract_prod_dense

    subroutine collect_vacuum_ids(mesh, ref_ids, uniq_ids, n_uniq)
        type(t_mesh_iga), intent(in)           :: mesh
        integer,          intent(in)           :: ref_ids(:)
        integer, allocatable, intent(out)      :: uniq_ids(:)
        integer,              intent(out)      :: n_uniq
        integer :: s, bc_id, j
        logical :: found

        allocate(uniq_ids(size(mesh%iga_surfaces)))
        n_uniq = 0
        do s = 1, size(mesh%iga_surfaces)
            bc_id = mesh%iga_surfaces(s)%bc_id
            if (any(bc_id == ref_ids)) cycle
            found = .false.
            do j = 1, n_uniq
                if (uniq_ids(j) == bc_id) then; found = .true.; exit; end if
            end do
            if (.not. found) then
                n_uniq = n_uniq + 1
                uniq_ids(n_uniq) = bc_id
            end if
        end do
    end subroutine collect_vacuum_ids

end module m_diffusion_iga
