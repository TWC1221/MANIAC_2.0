module m_constants
    implicit none

    ! Precision
    integer, parameter :: sp = selected_real_kind(4)
    integer, parameter :: dp = selected_real_kind(8)
    integer, parameter :: qp = selected_real_kind(16)

    ! Mathematical constants
    real(dp), parameter :: PI            = 4.0_dp * datan(1.0_dp)
    real(dp), parameter :: dp_EPSILON   = 1.0e-12_dp
    real(dp), parameter :: VSMALL_NUMBER = 1.0e-9_dp
    real(dp), parameter :: SMALL_NUMBER  = 1.0e-6_dp
    real(dp), parameter :: LARGE_NUMBER  = 1.0e+6_dp
    real(dp), parameter :: VLARGE_NUMBER = 1.0e+9_dp

    integer,  parameter :: MAX_ITERATIONS = 10000
    real(dp), parameter :: ADJUSTED_NEUTRON_MASS = 1.04625e-8_dp

    ! Formatting
    character,          parameter :: COMMENT_CHAR = '!'
    character(len=2),   parameter :: tab   = "  "
    character,          parameter :: space = " "

    ! PETSc KSP solver choices (diffusion linear solves)
    integer, parameter :: SOLVER_KSP_CG    = 1
    integer, parameter :: SOLVER_KSP_GMRES = 2
    integer, parameter :: SOLVER_KSP_BCGS  = 3

    ! Preconditioner choices
    integer, parameter :: PRECON_NONE      = 0
    integer, parameter :: PRECON_JACOBI    = 1
    integer, parameter :: PRECON_ILU       = 2
    integer, parameter :: PRECON_CHOLESKY  = 3
    integer, parameter :: PRECON_GAMG      = 4

    ! Boundary condition types
    integer, parameter :: BC_VACUUM     = 1
    integer, parameter :: BC_REFLECTIVE = 2
    integer, parameter :: BC_DIRICHLET  = 3
    integer, parameter :: BC_ALBEDO     = 4

    real(dp), parameter :: PENALTY = 1.0e10_dp

end module m_constants
