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

        ! Optional on-demand snapshot: called when SNAPSHOT trigger file is detected.
        ! Caller exports the current flux to VTK and returns; iteration continues.
        subroutine snapshot_fn(flux, k_eff, iter)
            import :: dp
            real(dp), intent(in) :: flux(:,:)
            real(dp), intent(in) :: k_eff
            integer,  intent(in) :: iter
        end subroutine snapshot_fn

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
                               build_source, do_solve, compute_prod, &
                               snapshot_export)
        real(dp), intent(inout)           :: flux(:,:)
        real(dp), intent(out)             :: k_eff
        integer,  intent(in)              :: max_outer
        real(dp), intent(in)              :: tol
        logical,  intent(in)              :: is_eigenvalue, is_adjoint
        procedure(source_fn)              :: build_source
        procedure(solve_fn)               :: do_solve
        procedure(production_fn)          :: compute_prod
        procedure(snapshot_fn), optional  :: snapshot_export

        real(dp), allocatable :: src(:,:), flux_old(:,:)
        real(dp) :: prod_new, prod_old, k_eff_old, err_phi, err_k, norm_phi, norm_old
        integer  :: outer, u_snap
        logical  :: snap_exists
        integer, parameter :: PRINT_STRIDE_EIG = 10, PRINT_STRIDE_FS = 5

        allocate(src    (size(flux,1), size(flux,2)), source=0.0_dp)
        allocate(flux_old(size(flux,1), size(flux,2)))

        k_eff    = 1.0_dp
        norm_old = 0.0_dp
        call compute_prod(prod_old, flux, is_adjoint)
        prod_old = max(prod_old, 1.0e-20_dp)

        write(*,*)
        write(*,'(A)') " |=======================================================================|"
        if (is_eigenvalue) then
            write(*,'(A)') " |                            POWER ITERATION                            |"
            write(*,'(A)') " |=======================================================================|"
            write(*,'(A)') " |   Iter        k_eff       dk_eff (pcm)                               |"
        else
            write(*,'(A)') " |                       FIXED-SOURCE ITERATION                          |"
            write(*,'(A)') " |=======================================================================|"
            write(*,'(A)') " |   Iter    Flux residual    Flux norm                                  |"
        end if
        write(*,'(A)') " |-----------------------------------------------------------------------|"

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

            if (is_eigenvalue) then
                if (mod(outer, PRINT_STRIDE_EIG) == 0 .or. outer == 1) &
                    write(*,'(A, I5, F15.8, F12.2, T74, A)') " |  ", outer, k_eff, err_k * 1.0e5_dp, "|"
            else
                if (mod(outer, PRINT_STRIDE_FS) == 0 .or. outer == 1) &
                    write(*,'(A, I5, 2X, ES15.6, 2X, ES15.6, T74, A)') " |  ", outer, err_phi, norm_phi, "|"
                ! Divergence guard: warn if flux has doubled since last print stride
                if (norm_old > 0.0_dp .and. norm_phi > 10.0_dp * norm_old) then
                    write(*,'(A)') " |-----------------------------------------------------------------------|"
                    write(*,'(A)') " |  WARNING: Flux is growing - problem may be super-critical.            |"
                    write(*,'(A)') " |  Check scatter ratio c = sum(SigmaS)/SigmaT and boundary conditions.  |"
                    write(*,'(A)') " |  For vacuum-only BC use n_ref_ids = 0 in config.nml.                  |"
                    write(*,'(A)') " |=======================================================================|"
                    write(*,*)
                    return
                end if
                if (mod(outer, PRINT_STRIDE_FS) == 0) norm_old = norm_phi
            end if

            ! ---- On-demand VTK snapshot trigger --------------------------------
            inquire(file="SNAPSHOT", exist=snap_exists)
            if (snap_exists) then
                if (present(snapshot_export)) then
                    write(*,'(A)') " |-----------------------------------------------------------------------|"
                    write(*,'(A, I5, A, F12.8, T74, A)') &
                        " |  [SNAP] iter", outer, "  k_eff = ", k_eff, "|"
                    call snapshot_export(flux, k_eff, outer)
                end if
                open(newunit=u_snap, file="SNAPSHOT")
                close(u_snap, status='delete')
            end if

            if (err_phi < tol .and. (err_k < tol .or. .not. is_eigenvalue)) then
                write(*,'(A)') " |-----------------------------------------------------------------------|"
                if (is_eigenvalue) then
                    write(*,'(A, I5, A, F12.8, T74, A)') " |  Converged in", outer, " iterations.  k_eff = ", k_eff, "|"
                else
                    write(*,'(A, I5, A, ES11.4, T74, A)') " |  Converged in", outer, &
                        " iterations.  flux residual = ", err_phi, "|"
                end if
                write(*,'(A)') " |=======================================================================|"
                write(*,*)
                exit
            end if
        end do

        if (outer > max_outer) then
            write(*,'(A)') " |-----------------------------------------------------------------------|"
            write(*,'(A)') " |  WARNING: Maximum outer iterations reached without convergence.       |"
            write(*,'(A)') " |=======================================================================|"
            write(*,*)
        end if

    end subroutine PowerIteration

end module m_power_iteration
