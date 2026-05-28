! Shared outer power-iteration driver for eigenvalue and fixed-source problems.
!
! All physics backends extend t_solver and implement three deferred TBPs.
! The solver owns its flux (scalar_flux component); no flux array is threaded
! through the callback chain.
!
!   build_source(self, src, k_eff, is_eigenvalue, is_adjoint)
!       reads  self%scalar_flux → writes src
!   do_solve(self, src)
!       reads  src              → writes self%scalar_flux (+ self%ang_flux for transport)
!   compute_prod(self, prod, is_adjoint)
!       reads  self%scalar_flux → writes prod
!
! On-demand VTK snapshots are triggered by a "SNAPSHOT" file in the run
! directory.  Override snapshot() in the concrete type and set have_snap=.true.
!
! Public: PowerIteration, t_solver
module m_power_iteration
    use m_constants
    implicit none
    public :: PowerIteration, t_solver

    ! ------------------------------------------------------------------
    ! Abstract solver base — all physics backends extend this.
    ! ------------------------------------------------------------------
    type, abstract :: t_solver
        real(dp), allocatable :: scalar_flux(:,:)   ! (n_dof, n_groups) — owned by solver
        logical            :: have_snap  = .false.
        character(len=256) :: snap_dir   = ""
        character(len=256) :: snap_tag   = ""
        integer            :: vtk_refine = 4
        integer            :: snap_count = 0
    contains
        procedure(i_build_source), deferred :: build_source
        procedure(i_do_solve),     deferred :: do_solve
        procedure(i_compute_prod), deferred :: compute_prod
        procedure :: snapshot => noop_snapshot
    end type t_solver

    abstract interface
        ! Build total source from self%scalar_flux into src.
        subroutine i_build_source(self, src, k_eff, is_eigenvalue, is_adjoint)
            import :: t_solver, dp
            class(t_solver), intent(inout) :: self
            real(dp), intent(inout) :: src(:,:)
            real(dp), intent(in)    :: k_eff
            logical,  intent(in)    :: is_eigenvalue, is_adjoint
        end subroutine i_build_source

        ! Advance self%scalar_flux by one sweep / linear solve using src.
        subroutine i_do_solve(self, src)
            import :: t_solver, dp
            class(t_solver), intent(inout) :: self
            real(dp), intent(in) :: src(:,:)
        end subroutine i_do_solve

        ! Integrate fission production from self%scalar_flux into prod.
        subroutine i_compute_prod(self, prod, is_adjoint)
            import :: t_solver, dp
            class(t_solver), intent(inout) :: self
            real(dp), intent(out) :: prod
            logical,  intent(in)  :: is_adjoint
        end subroutine i_compute_prod
    end interface

contains

    subroutine noop_snapshot(self, k_eff, iter)
        class(t_solver), intent(inout) :: self
        real(dp), intent(in) :: k_eff
        integer,  intent(in) :: iter
    end subroutine noop_snapshot

    ! ------------------------------------------------------------------
    ! Outer power-iteration loop.
    !
    ! solver%scalar_flux  -- initial guess on entry, converged solution on exit;
    !                        caller allocates and initialises before calling.
    ! k_eff               -- eigenvalue (always 1 for fixed-source problems).
    ! ------------------------------------------------------------------
    subroutine PowerIteration(solver, k_eff, max_outer, tol, is_eigenvalue, is_adjoint)
        class(t_solver), intent(inout) :: solver
        real(dp), intent(out)          :: k_eff
        integer,  intent(in)           :: max_outer
        real(dp), intent(in)           :: tol
        logical,  intent(in)           :: is_eigenvalue, is_adjoint

        real(dp), allocatable :: src(:,:), flux_old(:,:)
        real(dp) :: prod_new, prod_old, k_eff_old, err_phi, err_k, norm_phi, norm_old
        integer  :: outer, u_snap
        logical  :: snap_exists
        integer, parameter :: PRINT_STRIDE_EIG = 10, PRINT_STRIDE_FS = 5

        allocate(src    (size(solver%scalar_flux,1), size(solver%scalar_flux,2)), source=0.0_dp)
        allocate(flux_old(size(solver%scalar_flux,1), size(solver%scalar_flux,2)))

        k_eff    = 1.0_dp
        norm_old = 0.0_dp
        call solver%compute_prod(prod_old, is_adjoint)
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
            flux_old  = solver%scalar_flux
            k_eff_old = k_eff

            call solver%build_source(src, k_eff, is_eigenvalue, is_adjoint)
            call solver%do_solve(src)

            if (is_eigenvalue) then
                call solver%compute_prod(prod_new, is_adjoint)
                k_eff    = k_eff_old * prod_new / max(prod_old, 1.0e-20_dp)
                prod_old = prod_new
            end if

            norm_phi = maxval(abs(solver%scalar_flux))
            err_phi  = merge(maxval(abs(solver%scalar_flux - flux_old)) / norm_phi, &
                             0.0_dp, norm_phi > 0.0_dp)
            err_k    = abs(k_eff - k_eff_old)

            if (is_eigenvalue) then
                if (mod(outer, PRINT_STRIDE_EIG) == 0 .or. outer == 1) &
                    write(*,'(A, I5, F15.8, F12.2, T74, A)') " |  ", outer, k_eff, err_k * 1.0e5_dp, "|"
            else
                if (mod(outer, PRINT_STRIDE_FS) == 0 .or. outer == 1) &
                    write(*,'(A, I5, 2X, ES15.6, 2X, ES15.6, T74, A)') " |  ", outer, err_phi, norm_phi, "|"
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
                if (solver%have_snap) then
                    write(*,'(A)') " |-----------------------------------------------------------------------|"
                    write(*,'(A, I5, A, F12.8, T74, A)') &
                        " |  [SNAP] iter", outer, "  k_eff = ", k_eff, "|"
                    call solver%snapshot(k_eff, outer)
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
