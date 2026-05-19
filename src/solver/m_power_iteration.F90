! Shared outer power-iteration driver for eigenvalue and fixed-source problems.
!
! Physics modules (transport, diffusion) supply three callback procedures that
! conform to the abstract interfaces defined here:
!   build_source  -- construct the distributed source from the current flux
!   do_solve      -- advance the flux by one linear solve / sweep
!   compute_prod  -- evaluate total fission production (for k-eff update)
!
! The outer loop is identical for transport and diffusion, for IGA and FEM.
!
! Public: PowerIteration
module m_power_iteration
    use m_constants
    implicit none
    public :: PowerIteration

    ! ------------------------------------------------------------------
    ! Abstract interfaces for the three physics callbacks.
    ! ------------------------------------------------------------------
    abstract interface

        ! Build the distributed source for all groups from the current flux.
        ! For transport : wraps Source_DGFEM
        ! For diffusion : assembles scatter + fission RHS via CSR matvec
        subroutine source_fn(src, flux, k_eff, is_eigenvalue, is_adjoint)
            import :: dp
            real(dp), intent(inout) :: src(:,:)          ! (n_dof, n_groups) output
            real(dp), intent(in)    :: flux(:,:)         ! (n_dof, n_groups) current flux
            real(dp), intent(in)    :: k_eff
            logical,  intent(in)    :: is_eigenvalue, is_adjoint
        end subroutine source_fn

        ! Advance the flux by one step (sweep or linear solve).
        ! For transport : wraps Transport_Sweep  (rewrites flux from scratch)
        ! For diffusion : solves A x = src group by group (PCG or PETSc KSP)
        subroutine solve_fn(flux, src)
            import :: dp
            real(dp), intent(inout) :: flux(:,:)         ! in: initial guess / out: new flux
            real(dp), intent(in)    :: src(:,:)
        end subroutine solve_fn

        ! Compute total fission production  P = sum_g integral(nuSigF_g * phi_g dV).
        ! Used to update k-eff between outer iterations.
        subroutine production_fn(prod, flux, is_adjoint)
            import :: dp
            real(dp), intent(out) :: prod
            real(dp), intent(in)  :: flux(:,:)
            logical,  intent(in)  :: is_adjoint
        end subroutine production_fn

    end interface

contains

    ! ------------------------------------------------------------------
    ! Shared outer power-iteration loop.
    !
    ! flux        -- (n_dof, n_groups) initial guess on entry, converged
    !                solution on exit; caller allocates and initialises.
    ! k_eff       -- eigenvalue (always 1 for fixed-source problems).
    ! build_source, do_solve, compute_prod -- physics callbacks.
    ! ------------------------------------------------------------------
    subroutine PowerIteration(flux, k_eff, max_outer, tol, &
                               is_eigenvalue, is_adjoint,   &
                               build_source, do_solve, compute_prod)
        real(dp), intent(inout)      :: flux(:,:)
        real(dp), intent(out)        :: k_eff
        integer,  intent(in)         :: max_outer
        real(dp), intent(in)         :: tol
        logical,  intent(in)         :: is_eigenvalue, is_adjoint
        procedure(source_fn)         :: build_source
        procedure(solve_fn)          :: do_solve
        procedure(production_fn)     :: compute_prod

        real(dp), allocatable :: src(:,:), flux_old(:,:)
        real(dp) :: prod_new, prod_old, k_eff_old, err_phi, err_k, norm_phi
        integer  :: outer

        allocate(src    (size(flux,1), size(flux,2)), source=0.0_dp)
        allocate(flux_old(size(flux,1), size(flux,2)))

        k_eff = 1.0_dp
        call compute_prod(prod_old, flux, is_adjoint)
        prod_old = max(prod_old, 1.0e-20_dp)

        write(*,'(/,A)') " ========== POWER ITERATION =========="
        write(*,'(A5,A15,A12,A12)') "Iter", "k_eff", "err_k", "err_phi"

        do outer = 1, max_outer
            flux_old  = flux
            k_eff_old = k_eff

            call build_source(src, flux, k_eff, is_eigenvalue, is_adjoint)
            call do_solve(flux, src)

            if (is_eigenvalue) then
                call compute_prod(prod_new, flux, is_adjoint)
                k_eff    = k_eff_old * prod_new / max(prod_old, 1.0e-20_dp)
                prod_old = prod_new
            end if

            norm_phi = maxval(abs(flux))
            err_phi  = merge(maxval(abs(flux - flux_old)) / norm_phi, 0.0_dp, norm_phi > 0.0_dp)
            err_k    = abs(k_eff - k_eff_old)

            if (mod(outer, 10) == 0 .or. outer == 1) &
                write(*,'(I5,F15.8,E12.3,E12.3)') outer, k_eff, err_k, err_phi

            if (err_phi < tol .and. (err_k < tol .or. .not. is_eigenvalue)) then
                write(*,'(A,I5,A,F12.8)') " Converged in ", outer, " iterations.  k_eff = ", k_eff
                exit
            end if
        end do

        if (outer > max_outer) &
            write(*,'(A)') " WARNING: Maximum outer iterations reached without convergence."

    end subroutine PowerIteration

end module m_power_iteration
