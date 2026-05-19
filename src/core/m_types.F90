! Shared derived types used across all solvers (diffusion/transport, FEM/IGA).
module m_types
    use m_constants
    implicit none

    ! ------------------------------------------------------------------
    ! Nuclear cross-section data for one material, all energy groups.
    ! ------------------------------------------------------------------
    type t_material
        character(len=32)     :: name
        real(dp), allocatable :: D(:)        ! diffusion coefficient [cm]
        real(dp), allocatable :: SigmaT(:)   ! total XS [cm^-1]  (transport)
        real(dp), allocatable :: SigmaR(:)   ! removal XS [cm^-1] (diffusion)
        real(dp), allocatable :: SigA(:)     ! absorption XS [cm^-1]
        real(dp), allocatable :: SigF(:)     ! fission XS [cm^-1]
        real(dp), allocatable :: Nu(:)       ! average fission production [-]
        real(dp), allocatable :: NuSigF(:)   ! nu*fission XS [cm^-1]
        real(dp), allocatable :: Chi(:)      ! fission spectrum [-]
        real(dp), allocatable :: Src(:)      ! fixed source [cm^-3 s^-1]
        real(dp), allocatable :: SigmaS(:,:) ! scatter matrix SigmaS(from_g, to_g)
    end type t_material

    ! ------------------------------------------------------------------
    ! Spatial quadrature rule.
    ! ------------------------------------------------------------------
    type t_quadrature
        integer :: n_points
        real(dp), allocatable :: xi(:), eta(:), zeta(:), weights(:)
    end type t_quadrature

    ! ------------------------------------------------------------------
    ! Discrete-ordinates (Sn) angular quadrature set.
    ! ------------------------------------------------------------------
    type t_sn_quadrature
        integer :: order
        integer :: n_angles
        real(dp), allocatable :: dirs(:,:)    ! (n_angles, 3)
        real(dp), allocatable :: weights(:)   ! (n_angles)
    end type t_sn_quadrature

    ! ------------------------------------------------------------------
    ! Boundary condition for one physical surface ID.
    ! ------------------------------------------------------------------
    type t_bc_config
        integer  :: mat_id
        integer  :: bc_type   ! BC_VACUUM, BC_REFLECTIVE, BC_DIRICHLET, BC_ALBEDO
        real(dp) :: value
    end type t_bc_config

end module m_types
