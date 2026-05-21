#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>

! PETSc helpers for the multigroup diffusion solver.
!
! Public:
!   setup_ksp               -- configure a KSP solver from integer selectors
!   petsc_build_sparsity    -- compute per-row nnz from mesh topology
!   petsc_create_diff_mats  -- allocate and create multigroup Mat/Vec objects
!   petsc_assemble_diff_mats-- call MatAssembly/VecAssembly for all objects
!   petsc_setup_ksp_group   -- create per-group KSP and solution Vec
!   petsc_destroy_diff_state-- destroy all solver PETSc objects
module m_petsc
    use m_constants
    use petscsys
    use petscvec
    use petscmat
    use petscksp
    implicit none
    public :: setup_ksp
    public :: petsc_build_sparsity, petsc_create_diff_mats, petsc_assemble_diff_mats
    public :: petsc_setup_ksp_group, petsc_destroy_diff_state

contains

    subroutine setup_ksp(ksp_solver, A_mat, ksp_choice, pc_choice)
        KSP, intent(inout) :: ksp_solver
        Mat, intent(inout) :: A_mat
        integer, intent(in) :: ksp_choice, pc_choice

        PetscErrorCode :: ierr
        PC             :: pc_ctx

        PetscCall(MatAssemblyBegin(A_mat, MAT_FINAL_ASSEMBLY, ierr))
        PetscCall(MatAssemblyEnd(A_mat, MAT_FINAL_ASSEMBLY, ierr))
        PetscCall(MatSetFromOptions(A_mat, ierr))
        PetscCall(KSPSetOperators(ksp_solver, A_mat, A_mat, ierr))

        select case (ksp_choice)
        case (SOLVER_KSP_CG)
            PetscCall(KSPSetType(ksp_solver, KSPCG, ierr))
        case (SOLVER_KSP_GMRES)
            PetscCall(KSPSetType(ksp_solver, KSPGMRES, ierr))
        case (SOLVER_KSP_BCGS)
            PetscCall(KSPSetType(ksp_solver, KSPBCGS, ierr))
        end select

        PetscCall(KSPGetPC(ksp_solver, pc_ctx, ierr))
        select case (pc_choice)
        case (PRECON_NONE)
            PetscCall(PCSetType(pc_ctx, PCNONE, ierr))
        case (PRECON_JACOBI)
            PetscCall(PCSetType(pc_ctx, PCJACOBI, ierr))
        case (PRECON_ILU)
            PetscCall(PCSetType(pc_ctx, PCILU, ierr))
        case (PRECON_CHOLESKY)
            PetscCall(PCSetType(pc_ctx, PCICC, ierr))
        case (PRECON_GAMG)
            PetscCall(PCSetType(pc_ctx, PCGAMG, ierr))
        end select

        PetscCall(KSPSetFromOptions(ksp_solver, ierr))
    end subroutine setup_ksp

    ! ------------------------------------------------------------------
    ! Compute per-row non-zero counts from mesh topology.
    ! Each node receives contributions from every element it belongs to,
    ! and each element contributes n_basis columns.
    ! ------------------------------------------------------------------
    subroutine petsc_build_sparsity(n_nodes, n_elems, n_basis, elems, nnz)
        integer,   intent(in)  :: n_nodes, n_elems, n_basis
        integer,   intent(in)  :: elems(:,:)
        PetscInt, allocatable, intent(out) :: nnz(:)
        integer :: ee, i

        allocate(nnz(n_nodes)); nnz = 0
        do ee = 1, n_elems
            do i = 1, n_basis
                nnz(elems(ee,i)) = nnz(elems(ee,i)) + n_basis
            end do
        end do
        do i = 1, n_nodes
            nnz(i) = min(nnz(i), n_nodes)
        end do
    end subroutine petsc_build_sparsity

    ! ------------------------------------------------------------------
    ! Allocate and create all PETSc Mat/Vec objects for multigroup diffusion.
    ! A_MAT(g)        -- group stiffness matrix
    ! MAT_F(g,gp)     -- fission coupling matrix g <- gp
    ! MAT_S(g,gp)     -- scatter coupling matrix g <- gp
    ! PROD_VEC(g)     -- fission production vector
    ! FixedSrc(g)     -- fixed external source vector
    ! ------------------------------------------------------------------
    subroutine petsc_create_diff_mats(n_nodes, n_groups, nnz, &
                                       A_MAT, MAT_F, MAT_S, PROD_VEC, FixedSrc)
        integer,    intent(in)  :: n_nodes, n_groups
        PetscInt,   intent(in)  :: nnz(:)
        Mat, allocatable, intent(out) :: A_MAT(:), MAT_F(:,:), MAT_S(:,:)
        Vec, allocatable, intent(out) :: PROD_VEC(:), FixedSrc(:)

        integer :: g, gp
        PetscErrorCode :: ierr

        allocate(A_MAT(n_groups), MAT_F(n_groups,n_groups), &
                 MAT_S(n_groups,n_groups), PROD_VEC(n_groups), FixedSrc(n_groups))

        do g = 1, n_groups
            call VecCreateSeq(PETSC_COMM_SELF, n_nodes, FixedSrc(g), ierr)
            call VecCreateSeq(PETSC_COMM_SELF, n_nodes, PROD_VEC(g), ierr)
            call MatCreateSeqAIJ(PETSC_COMM_SELF, n_nodes, n_nodes, 0, nnz, A_MAT(g), ierr)
            call MatSetOption(A_MAT(g), MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
            do gp = 1, n_groups
                call MatCreateSeqAIJ(PETSC_COMM_SELF, n_nodes, n_nodes, 0, nnz, MAT_F(g,gp), ierr)
                call MatSetOption(MAT_F(g,gp), MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
                call MatCreateSeqAIJ(PETSC_COMM_SELF, n_nodes, n_nodes, 0, nnz, MAT_S(g,gp), ierr)
                call MatSetOption(MAT_S(g,gp), MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
            end do
        end do
    end subroutine petsc_create_diff_mats

    ! ------------------------------------------------------------------
    ! Call assembly begin/end on all multigroup diffusion objects.
    ! ------------------------------------------------------------------
    subroutine petsc_assemble_diff_mats(A_MAT, MAT_F, MAT_S, PROD_VEC, FixedSrc, n_groups)
        Mat, intent(inout) :: A_MAT(:), MAT_F(:,:), MAT_S(:,:)
        Vec, intent(inout) :: PROD_VEC(:), FixedSrc(:)
        integer, intent(in) :: n_groups

        integer :: g, gp
        PetscErrorCode :: ierr

        do g = 1, n_groups
            call VecAssemblyBegin(FixedSrc(g), ierr); call VecAssemblyEnd(FixedSrc(g), ierr)
            call VecAssemblyBegin(PROD_VEC(g), ierr); call VecAssemblyEnd(PROD_VEC(g), ierr)
            call MatAssemblyBegin(A_MAT(g), MAT_FINAL_ASSEMBLY, ierr)
            call MatAssemblyEnd  (A_MAT(g), MAT_FINAL_ASSEMBLY, ierr)
            do gp = 1, n_groups
                call MatAssemblyBegin(MAT_F(g,gp), MAT_FINAL_ASSEMBLY, ierr)
                call MatAssemblyEnd  (MAT_F(g,gp), MAT_FINAL_ASSEMBLY, ierr)
                call MatAssemblyBegin(MAT_S(g,gp), MAT_FINAL_ASSEMBLY, ierr)
                call MatAssemblyEnd  (MAT_S(g,gp), MAT_FINAL_ASSEMBLY, ierr)
            end do
        end do
    end subroutine petsc_assemble_diff_mats

    ! ------------------------------------------------------------------
    ! For each group: create a KSP from A_MAT(g), create solution Vec.
    ! A_MAT is consumed (destroyed) after KSP setup.
    ! ------------------------------------------------------------------
    subroutine petsc_setup_ksp_group(A_MAT, n_nodes, n_groups, &
                                      solver_type, preconditioner, KSPs, X_petsc)
        Mat, intent(inout) :: A_MAT(:)
        integer, intent(in) :: n_nodes, n_groups, solver_type, preconditioner
        KSP, allocatable, intent(out) :: KSPs(:)
        Vec, allocatable, intent(out) :: X_petsc(:)

        integer :: g
        PetscErrorCode :: ierr

        allocate(KSPs(n_groups), X_petsc(n_groups))
        do g = 1, n_groups
            call KSPCreate(PETSC_COMM_SELF, KSPs(g), ierr)
            call setup_ksp(KSPs(g), A_MAT(g), solver_type, preconditioner)
            call VecCreateSeq(PETSC_COMM_SELF, n_nodes, X_petsc(g), ierr)
            call VecSet(X_petsc(g), 1.0_dp/real(n_nodes,dp), ierr)
            call MatDestroy(A_MAT(g), ierr)
        end do
    end subroutine petsc_setup_ksp_group

    ! ------------------------------------------------------------------
    ! Destroy all PETSc solver state from a previous SolveDiffusion call.
    ! Guards against unallocated arrays.
    ! ------------------------------------------------------------------
    subroutine petsc_destroy_diff_state(KSPs, X_petsc, MAT_F, MAT_S, FixedSrc, &
                                         tmp_b, tmp_x, tmp_valid)
        KSP, allocatable, intent(inout) :: KSPs(:)
        Vec, allocatable, intent(inout) :: X_petsc(:), FixedSrc(:)
        Mat, allocatable, intent(inout) :: MAT_F(:,:), MAT_S(:,:)
        Vec,              intent(inout) :: tmp_b, tmp_x
        logical,          intent(inout) :: tmp_valid

        integer :: g, gp
        PetscErrorCode :: ierr

        if (allocated(KSPs)) then
            do g = 1, size(KSPs); call KSPDestroy(KSPs(g), ierr); end do
            deallocate(KSPs)
        end if
        if (allocated(X_petsc)) then
            do g = 1, size(X_petsc); call VecDestroy(X_petsc(g), ierr); end do
            deallocate(X_petsc)
        end if
        if (allocated(MAT_F)) then
            do g = 1, size(MAT_F,1); do gp = 1, size(MAT_F,2)
                call MatDestroy(MAT_F(g,gp), ierr)
                call MatDestroy(MAT_S(g,gp), ierr)
            end do; end do
            deallocate(MAT_F, MAT_S)
        end if
        if (allocated(FixedSrc)) then
            do g = 1, size(FixedSrc); call VecDestroy(FixedSrc(g), ierr); end do
            deallocate(FixedSrc)
        end if
        if (tmp_valid) then
            call VecDestroy(tmp_b, ierr); call VecDestroy(tmp_x, ierr)
            tmp_valid = .false.
        end if
    end subroutine petsc_destroy_diff_state

end module m_petsc
