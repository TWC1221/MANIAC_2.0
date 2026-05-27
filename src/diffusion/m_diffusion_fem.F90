#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>

! FEM CG multigroup diffusion solver: element assembly, BC application, and power iteration.
! Supports 2D (quad elements) and 3D (hex elements) Lagrange meshes.
!
! Public:
!   SolveDiffusion_FEM     -- high-level entry: assemble + power iteration
!   assemble_petsc_fem     -- PETSc matrices (A, F, S, ProdVec, FixedSrc)
!   apply_bcs_fem          -- Robin (vacuum/albedo) and Dirichlet BCs
module m_diffusion_fem
    use m_constants
    use m_types
    use m_material
    use m_quadrature
    use m_basis_fem,       only: GetMapping2D_FEM, GetMapping3D_FEM
    use m_petsc,           only: setup_ksp, petsc_build_sparsity, petsc_create_diff_mats, &
                                  petsc_assemble_diff_mats, petsc_setup_ksp_group, &
                                  petsc_destroy_diff_state
    use m_power_iteration, only: PowerIteration
    use m_output_fem,      only: export_diffusion_vtk_fem
    use petscsys
    use petscvec
    use petscmat
    use petscksp
    implicit none
    public :: SolveDiffusion_FEM
    public :: assemble_petsc_fem, apply_bcs_fem

    ! ------------------------------------------------------------------
    ! Module-level saved state — set by SolveDiffusion_FEM, read by callbacks.
    ! ------------------------------------------------------------------
    type(t_mesh_fem),   pointer, save :: s_d_mesh    => null()
    type(t_basis_fem), pointer, save :: s_d_FE      => null()
    type(t_material),   pointer, save :: s_d_mats(:) => null()

    KSP, allocatable, save :: s_d_KSPs(:)
    Vec, allocatable, save :: s_d_X_petsc(:)
    Mat, allocatable, save :: s_d_MAT_F(:,:)
    Mat, allocatable, save :: s_d_MAT_S(:,:)
    Vec, allocatable, save :: s_d_FixedSrc(:)
    Vec, save              :: s_d_tmp_b, s_d_tmp_x
    logical, save          :: s_d_tmp_valid = .false.
    real(dp), allocatable, save :: s_d_prod_dense(:,:)
    integer, save :: s_d_n_groups = 0
    integer, save :: s_d_n_nodes  = 0

    ! Snapshot state
    logical,            save :: s_d_have_snap  = .false.
    character(len=256), save :: s_d_snap_dir   = ""
    character(len=256), save :: s_d_snap_tag   = ""
    integer,            save :: s_d_vtk_refine = 4
    integer,            save :: s_d_snap_count = 0

contains

    ! ------------------------------------------------------------------
    ! Build PETSc multigroup diffusion matrices for a FEM mesh.
    ! ------------------------------------------------------------------
    subroutine assemble_petsc_fem(A_MAT, MAT_F, MAT_S, PROD_VEC, FixedSrc, &
                                   mesh, FE, Quad, mats, n_groups, is_adjoint)
        Mat, allocatable, intent(out) :: A_MAT(:), MAT_F(:,:), MAT_S(:,:)
        Vec, allocatable, intent(out) :: PROD_VEC(:), FixedSrc(:)
        type(t_mesh_fem),   intent(in) :: mesh
        type(t_basis_fem), intent(in) :: FE
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
        write(*,'(A)') " [ FEM MATRIX ] :: Starting element-wise assembly..."
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ee)
        do ee = 1, mesh%n_elems
            call assemble_fem_elem_petsc(ee, mesh, FE, Quad, mats, n_groups, is_adjoint, &
                                         A_MAT, MAT_F, MAT_S, PROD_VEC, FixedSrc, count)
        end do
        !$OMP END PARALLEL DO

        call petsc_assemble_diff_mats(A_MAT, MAT_F, MAT_S, PROD_VEC, FixedSrc, n_groups)
    end subroutine assemble_petsc_fem

    ! ------------------------------------------------------------------
    ! Internal: assemble one FEM element into PETSc matrices.
    ! ------------------------------------------------------------------
    subroutine assemble_fem_elem_petsc(ee, mesh, FE, Quad, mats, n_groups, is_adjoint, &
                                       A_MAT, MAT_F, MAT_S, PROD_VEC, FixedSrc, count)
        integer,            intent(in)    :: ee, n_groups
        type(t_mesh_fem),   intent(in)    :: mesh
        type(t_basis_fem), intent(in)    :: FE
        type(t_quadrature), intent(in)    :: Quad
        type(t_material),   intent(in)    :: mats(:)
        logical,            intent(in)    :: is_adjoint
        Mat,                intent(inout) :: A_MAT(:), MAT_F(:,:), MAT_S(:,:)
        Vec,                intent(inout) :: PROD_VEC(:), FixedSrc(:)
        integer,            intent(inout) :: count

        integer  :: i, j, q, mat_id, g_to, g_from
        real(dp) :: dN_dx(FE%n_basis), dN_dy(FE%n_basis), dN_dz(FE%n_basis)
        real(dp) :: N_basis(FE%n_basis), detJ, dV
        real(dp) :: ec(FE%n_basis, 3)
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

        do q = 1, Quad%n_points
            if (mesh%dim == 2) then
                call GetMapping2D_FEM(FE, q, ec(:,1:2), dN_dx, dN_dy, detJ, N_basis)
                dN_dz = 0.0_dp
            else
                call GetMapping3D_FEM(FE, q, ec, dN_dx, dN_dy, dN_dz, detJ, N_basis)
            end if
            dV = abs(detJ) * Quad%weights(q)

            do g_to = 1, n_groups
                do i = 1, FE%n_basis
                    N_i = N_basis(i)
                    nsf = merge(mats(mat_id)%NuSigF(g_to), mats(mat_id)%Chi(g_to), .not. is_adjoint)
                    loc_Prod(i,g_to) = loc_Prod(i,g_to) + nsf * N_i * dV
                    loc_Src(i,g_to)  = loc_Src(i,g_to)  + mats(mat_id)%Src(g_to) * N_i * dV
                    do j = 1, FE%n_basis
                        N_j   = N_basis(j)
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
    end subroutine assemble_fem_elem_petsc

    ! ------------------------------------------------------------------
    ! Apply diffusion boundary conditions to a PETSc matrix.
    ! ------------------------------------------------------------------
    subroutine apply_bcs_fem(mesh, FE, QuadFace, bc_cfg, A)
        type(t_mesh_fem),   intent(in)    :: mesh
        type(t_basis_fem), intent(in)    :: FE
        type(t_quadrature), intent(in)    :: QuadFace
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
            do s = 1, size(mesh%surfaces)
                if (mesh%surfaces(s)%bc_id /= bc_cfg%mat_id) cycle
                do k = 1, size(mesh%surfaces(s)%cp_ids)
                    node_id = mesh%surfaces(s)%cp_ids(k)
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
            call fem_robin_petsc(mesh, FE, QuadFace, bc_cfg%mat_id, eff_alpha, A)
        end select
    end subroutine apply_bcs_fem

    ! ------------------------------------------------------------------
    ! Robin BC — scan all element faces for boundary faces with target_id.
    ! 2D: 1D edge integrals.  3D: 2D face integrals.
    ! ------------------------------------------------------------------
    subroutine fem_robin_petsc(mesh, FE, QuadFace, target_id, alpha, A)
        type(t_mesh_fem),   intent(in)    :: mesh
        type(t_basis_fem), intent(in)    :: FE
        type(t_quadrature), intent(in)    :: QuadFace
        integer,            intent(in)    :: target_id
        real(dp),           intent(in)    :: alpha
        Mat,                intent(inout) :: A

        integer, allocatable :: node_bc_id(:)
        integer :: s, k, ee, f, q, i, j, nf
        integer :: face_global(FE%n_nodes_per_face)
        real(dp) :: beta, dS, dV
        real(dp) :: x_face(FE%n_nodes_per_face, 3)
        real(dp) :: N_face(FE%n_nodes_per_face)
        real(dp) :: dN_dxi(FE%n_nodes_per_face), dN_deta(FE%n_nodes_per_face)
        real(dp) :: dx_du, dy_du, dz_du, dx_dv, dy_dv, dz_dv
        real(dp) :: nx, ny, nz
        real(dp) :: loc_B(FE%n_nodes_per_face, FE%n_nodes_per_face)
        PetscInt :: idx(FE%n_nodes_per_face)
        PetscErrorCode :: ierr

        beta = 0.5_dp * (1.0_dp - alpha) / (1.0_dp + alpha)
        nf   = FE%n_nodes_per_face

        ! Build per-node bc_id lookup (last surface wins if shared, which shouldn't happen)
        allocate(node_bc_id(mesh%n_nodes)); node_bc_id = 0
        do s = 1, size(mesh%surfaces)
            do k = 1, size(mesh%surfaces(s)%cp_ids)
                node_bc_id(mesh%surfaces(s)%cp_ids(k)) = mesh%surfaces(s)%bc_id
            end do
        end do

        do ee = 1, mesh%n_elems
            do f = 1, mesh%n_faces_per_elem
                ! Collect global node IDs for this face
                do k = 1, nf
                    face_global(k) = mesh%elems(ee, FE%face_node_map(k,f))
                end do

                ! Skip if not all face nodes belong to the target BC surface
                if (.not. all(node_bc_id(face_global) == target_id)) cycle

                ! Physical coordinates of face nodes
                do k = 1, nf
                    x_face(k,:) = mesh%nodes(face_global(k),:)
                end do

                loc_B = 0.0_dp
                idx   = face_global - 1

                if (mesh%dim == 2) then
                    ! 1D edge integral
                    do q = 1, QuadFace%n_points
                        N_face  = FE%basis_at_face_quad(q,:)
                        dN_dxi  = FE%dbasis_face_dxi(q,:)
                        dx_du = dot_product(dN_dxi, x_face(:,1))
                        dy_du = dot_product(dN_dxi, x_face(:,2))
                        dS    = sqrt(dx_du**2 + dy_du**2)
                        dV    = dS * QuadFace%weights(q) * beta
                        do i = 1, nf
                            do j = 1, nf
                                loc_B(i,j) = loc_B(i,j) + N_face(i) * N_face(j) * dV
                            end do
                        end do
                    end do
                else
                    ! 2D face integral (3D mesh)
                    do q = 1, QuadFace%n_points
                        N_face  = FE%basis_at_face_quad(q,:)
                        dN_dxi  = FE%dbasis_face_dxi(q,:)
                        dN_deta = FE%dbasis_face_deta(q,:)
                        dx_du = dot_product(dN_dxi,  x_face(:,1))
                        dy_du = dot_product(dN_dxi,  x_face(:,2))
                        dz_du = dot_product(dN_dxi,  x_face(:,3))
                        dx_dv = dot_product(dN_deta, x_face(:,1))
                        dy_dv = dot_product(dN_deta, x_face(:,2))
                        dz_dv = dot_product(dN_deta, x_face(:,3))
                        nx = dy_du*dz_dv - dz_du*dy_dv
                        ny = dz_du*dx_dv - dx_du*dz_dv
                        nz = dx_du*dy_dv - dy_du*dx_dv
                        dS = sqrt(nx**2 + ny**2 + nz**2)
                        dV = dS * QuadFace%weights(q) * beta
                        do i = 1, nf
                            do j = 1, nf
                                loc_B(i,j) = loc_B(i,j) + N_face(i) * N_face(j) * dV
                            end do
                        end do
                    end do
                end if

                call MatSetValues(A, nf, idx, nf, idx, reshape(loc_B, [nf**2]), ADD_VALUES, ierr)
            end do
        end do

        deallocate(node_bc_id)
    end subroutine fem_robin_petsc

    ! ==================================================================
    ! High-level entry point.
    ! ==================================================================
    subroutine SolveDiffusion_FEM(mesh, FE, Quad, QuadFace, mats,   &
                                   solver_type, preconditioner, ref_ids, &
                                   max_outer, tol, is_eigenvalue, is_adjoint, &
                                   phi_out, k_eff_out, snap_dir, snap_tag, vtk_ref)
        type(t_mesh_fem),   intent(in), target :: mesh
        type(t_basis_fem), intent(in), target :: FE
        type(t_quadrature), intent(in)         :: Quad, QuadFace
        type(t_material),   intent(in), target :: mats(:)
        integer,            intent(in)         :: solver_type, preconditioner
        integer,            intent(in)         :: ref_ids(:)
        integer,            intent(in)         :: max_outer
        real(dp),           intent(in)         :: tol
        logical,            intent(in)         :: is_eigenvalue, is_adjoint
        real(dp), allocatable, intent(out)     :: phi_out(:,:)
        real(dp),              intent(out)     :: k_eff_out
        character(len=*), optional, intent(in) :: snap_dir, snap_tag
        integer,          optional, intent(in) :: vtk_ref

        integer :: g, ss, n_uniq
        type(t_bc_config)    :: bc_vac
        integer, allocatable :: uniq_bc_ids(:)
        PetscErrorCode       :: ierr
        Mat, allocatable :: loc_A(:), loc_MF(:,:), loc_MS(:,:)
        Vec, allocatable :: loc_PV(:), loc_FS(:)

        s_d_mesh     => mesh
        s_d_FE       => FE
        s_d_mats     => mats
        s_d_n_groups = mesh%n_groups
        s_d_n_nodes  = mesh%n_nodes

        s_d_have_snap = present(snap_dir) .and. present(snap_tag)
        if (s_d_have_snap) then
            s_d_snap_dir   = trim(snap_dir)
            s_d_snap_tag   = trim(snap_tag)
            s_d_vtk_refine = merge(vtk_ref, 4, present(vtk_ref))
            s_d_snap_count = 0
        end if

        call petsc_destroy_diff_state(s_d_KSPs, s_d_X_petsc, s_d_MAT_F, s_d_MAT_S, &
                                       s_d_FixedSrc, s_d_tmp_b, s_d_tmp_x, s_d_tmp_valid)

        call assemble_petsc_fem(loc_A, loc_MF, loc_MS, loc_PV, loc_FS, &
                                 mesh, FE, Quad, mats, mesh%n_groups, is_adjoint)

        call collect_vacuum_ids_fem(mesh, ref_ids, uniq_bc_ids, n_uniq)
        bc_vac%bc_type = BC_VACUUM; bc_vac%value = 0.0_dp
        do ss = 1, n_uniq
            bc_vac%mat_id = uniq_bc_ids(ss)
            do g = 1, mesh%n_groups
                call apply_bcs_fem(mesh, FE, QuadFace, bc_vac, loc_A(g))
            end do
        end do
        deallocate(uniq_bc_ids)

        allocate(s_d_MAT_F(mesh%n_groups, mesh%n_groups), &
                 s_d_MAT_S(mesh%n_groups, mesh%n_groups), &
                 s_d_FixedSrc(mesh%n_groups))
        s_d_MAT_F    = loc_MF
        s_d_MAT_S    = loc_MS
        s_d_FixedSrc = loc_FS
        deallocate(loc_MF, loc_MS, loc_FS)

        call petsc_setup_ksp_group(loc_A, mesh%n_nodes, mesh%n_groups, &
                                    solver_type, preconditioner, s_d_KSPs, s_d_X_petsc)
        deallocate(loc_A)

        call extract_prod_dense_fem(loc_PV, mesh%n_groups, mesh%n_nodes)
        do g = 1, mesh%n_groups; call VecDestroy(loc_PV(g), ierr); end do
        deallocate(loc_PV)

        call VecCreateSeq(PETSC_COMM_SELF, mesh%n_nodes, s_d_tmp_b, ierr)
        call VecCreateSeq(PETSC_COMM_SELF, mesh%n_nodes, s_d_tmp_x, ierr)
        s_d_tmp_valid = .true.

        allocate(phi_out(mesh%n_nodes, mesh%n_groups), source=1.0_dp)

        if (s_d_have_snap) then
            call PowerIteration(phi_out, k_eff_out, max_outer, tol, &
                                 is_eigenvalue, is_adjoint, &
                                 diffusion_source_fem, diffusion_solve_fem, diffusion_production_fem, &
                                 snapshot_export=diff_fem_snapshot)
        else
            call PowerIteration(phi_out, k_eff_out, max_outer, tol, &
                                 is_eigenvalue, is_adjoint, &
                                 diffusion_source_fem, diffusion_solve_fem, diffusion_production_fem)
        end if
    end subroutine SolveDiffusion_FEM

    subroutine diffusion_source_fem(src, flux, k_eff, is_eigenvalue, is_adjoint)
        real(dp), intent(inout) :: src(:,:)
        real(dp), intent(in)    :: flux(:,:), k_eff
        logical,  intent(in)    :: is_eigenvalue, is_adjoint
        integer        :: g, gp
        PetscErrorCode :: ierr
        PetscScalar, pointer :: parr(:)
        src = 0.0_dp
        do g = 1, s_d_n_groups
            if (.not. is_eigenvalue) then
                call VecGetArrayRead(s_d_FixedSrc(g), parr, ierr)
                src(:,g) = parr(:)
                call VecRestoreArrayRead(s_d_FixedSrc(g), parr, ierr)
            end if
            do gp = 1, s_d_n_groups
                call VecGetArray(s_d_tmp_x, parr, ierr)
                parr(:) = flux(:,gp)
                call VecRestoreArray(s_d_tmp_x, parr, ierr)
                if (gp /= g) then
                    call MatMult(s_d_MAT_S(g,gp), s_d_tmp_x, s_d_tmp_b, ierr)
                    call VecGetArrayRead(s_d_tmp_b, parr, ierr)
                    src(:,g) = src(:,g) + parr(:)
                    call VecRestoreArrayRead(s_d_tmp_b, parr, ierr)
                end if
                if (is_eigenvalue) then
                    call MatMult(s_d_MAT_F(g,gp), s_d_tmp_x, s_d_tmp_b, ierr)
                    call VecGetArrayRead(s_d_tmp_b, parr, ierr)
                    src(:,g) = src(:,g) + parr(:) / k_eff
                    call VecRestoreArrayRead(s_d_tmp_b, parr, ierr)
                end if
            end do
        end do
    end subroutine diffusion_source_fem

    subroutine diffusion_solve_fem(flux, src)
        real(dp), intent(inout) :: flux(:,:)
        real(dp), intent(in)    :: src(:,:)
        integer        :: g
        PetscErrorCode :: ierr
        PetscScalar, pointer :: parr(:)

        do g = 1, s_d_n_groups
            call VecGetArray(s_d_tmp_b, parr, ierr)
            parr(:) = src(:,g)
            call VecRestoreArray(s_d_tmp_b, parr, ierr)
            call KSPSolve(s_d_KSPs(g), s_d_tmp_b, s_d_X_petsc(g), ierr)
            call VecGetArrayRead(s_d_X_petsc(g), parr, ierr)
            flux(:,g) = parr(:)
            call VecRestoreArrayRead(s_d_X_petsc(g), parr, ierr)
        end do
    end subroutine diffusion_solve_fem

    subroutine diffusion_production_fem(prod, flux, is_adjoint)
        real(dp), intent(out) :: prod
        real(dp), intent(in)  :: flux(:,:)
        logical,  intent(in)  :: is_adjoint
        integer :: g

        prod = 0.0_dp
        do g = 1, s_d_n_groups
            prod = prod + dot_product(s_d_prod_dense(:,g), flux(:,g))
        end do
    end subroutine diffusion_production_fem

    subroutine extract_prod_dense_fem(prod_vecs, n_groups, n_nodes)
        Vec,     intent(in) :: prod_vecs(:)
        integer, intent(in) :: n_groups, n_nodes
        integer :: g
        PetscErrorCode :: ierr
        PetscScalar, pointer :: parr(:)

        if (allocated(s_d_prod_dense)) deallocate(s_d_prod_dense)
        allocate(s_d_prod_dense(n_nodes, n_groups))
        do g = 1, n_groups
            call VecGetArrayRead(prod_vecs(g), parr, ierr)
            s_d_prod_dense(:,g) = parr(:)
            call VecRestoreArrayRead(prod_vecs(g), parr, ierr)
        end do
    end subroutine extract_prod_dense_fem

    subroutine collect_vacuum_ids_fem(mesh, ref_ids, uniq_ids, n_uniq)
        type(t_mesh_fem), intent(in)           :: mesh
        integer,          intent(in)           :: ref_ids(:)
        integer, allocatable, intent(out)      :: uniq_ids(:)
        integer,              intent(out)      :: n_uniq
        integer :: s, bc_id, j
        logical :: found

        allocate(uniq_ids(size(mesh%surfaces)))
        n_uniq = 0
        do s = 1, size(mesh%surfaces)
            bc_id = mesh%surfaces(s)%bc_id
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
    end subroutine collect_vacuum_ids_fem

    subroutine diff_fem_snapshot(flux, k_eff, iter)
        real(dp), intent(in) :: flux(:,:)
        real(dp), intent(in) :: k_eff
        integer,  intent(in) :: iter
        character(len=32)  :: lbl
        character(len=512) :: stag
        s_d_snap_count = s_d_snap_count + 1
        write(lbl,'(A,I4.4)') "snap", s_d_snap_count
        stag = trim(s_d_snap_tag) // "_" // trim(lbl)
        call export_diffusion_vtk_fem(trim(s_d_snap_dir), trim(stag), &
            s_d_mesh, s_d_FE, flux, s_d_n_groups, s_d_vtk_refine)
    end subroutine diff_fem_snapshot

end module m_diffusion_fem
