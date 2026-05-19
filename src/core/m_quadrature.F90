!-----------------------------------------------------------------------
! Spatial and angular quadrature generators.
!
! Exported routines (all mesh/FE-type agnostic):
!   LinearQuadrature        -- 1-D Gauss-Legendre nodes via LAPACK DSTEQR
!   QuadrilateralQuadrature -- 2-D tensor-product Gauss-Legendre
!   HexahedralQuadrature    -- 3-D tensor-product Gauss-Legendre
!   Spectral1DQuadrature    -- 1-D GLL (Gauss-Lobatto-Legendre) nodes
!   Spectral2DQuadrature    -- 2-D tensor-product GLL
!   AngularQuadrature       -- Level-symmetric Sn set (2-D or 3-D)
!
! Method-specific InitialiseQuadrature wrappers belong in src/fem/ and
! src/iga/ where the FE / mesh types are known.
!-----------------------------------------------------------------------
module m_quadrature
    use m_constants
    use m_types, only: t_quadrature, t_sn_quadrature
    implicit none
    private
    public :: LinearQuadrature, QuadrilateralQuadrature, HexahedralQuadrature, &
              Spectral1DQuadrature, Spectral2DQuadrature, AngularQuadrature

contains

    ! ------------------------------------------------------------------
    ! 1-D Gauss-Legendre quadrature (n points) via symmetric tridiagonal
    ! eigenproblem solved by LAPACK DSTEQR.
    ! ------------------------------------------------------------------
    subroutine LinearQuadrature(Quad, n)
        type(t_quadrature), intent(out) :: Quad
        integer,            intent(in)  :: n

        real(dp), allocatable :: d(:), e(:), z(:,:), work(:)
        integer :: info, i

        Quad%n_points = n
        allocate(Quad%xi(n), Quad%weights(n))
        allocate(d(n), e(n-1), z(n,n), work(max(1, 2*n-2)))

        d = 0.0_dp
        do i = 1, n-1
            e(i) = real(i, dp) / sqrt(4.0_dp * real(i,dp)**2 - 1.0_dp)
        end do
        z = 0.0_dp
        do i = 1, n; z(i,i) = 1.0_dp; end do

        call DSTEQR('V', n, d, e, z, n, work, info)
        if (info /= 0) stop "LinearQuadrature: DSTEQR failed."

        Quad%xi = d
        do i = 1, n
            Quad%weights(i) = 2.0_dp * z(1,i)**2
        end do
    end subroutine LinearQuadrature

    ! ------------------------------------------------------------------
    ! 2-D tensor-product Gauss-Legendre over [-1,1]^2.
    ! ------------------------------------------------------------------
    subroutine QuadrilateralQuadrature(Quad, n)
        type(t_quadrature), intent(out) :: Quad
        integer,            intent(in)  :: n

        type(t_quadrature) :: Q1
        integer :: i, j, k

        call LinearQuadrature(Q1, n)
        Quad%n_points = n**2
        allocate(Quad%xi(n**2), Quad%eta(n**2), Quad%weights(n**2))

        k = 0
        do i = 1, n
            do j = 1, n
                k = k + 1
                Quad%xi(k)      = Q1%xi(i)
                Quad%eta(k)     = Q1%xi(j)
                Quad%weights(k) = Q1%weights(i) * Q1%weights(j)
            end do
        end do
    end subroutine QuadrilateralQuadrature

    ! ------------------------------------------------------------------
    ! 3-D tensor-product Gauss-Legendre over [-1,1]^3.
    ! ------------------------------------------------------------------
    subroutine HexahedralQuadrature(Quad, n)
        type(t_quadrature), intent(out) :: Quad
        integer,            intent(in)  :: n

        type(t_quadrature) :: Q1
        integer :: ii, jj, kk, gp

        call LinearQuadrature(Q1, n)
        Quad%n_points = n**3
        allocate(Quad%xi(n**3), Quad%eta(n**3), Quad%zeta(n**3), Quad%weights(n**3))

        gp = 0
        do kk = 1, n
            do jj = 1, n
                do ii = 1, n
                    gp = gp + 1
                    Quad%xi(gp)      = Q1%xi(ii)
                    Quad%eta(gp)     = Q1%xi(jj)
                    Quad%zeta(gp)    = Q1%xi(kk)
                    Quad%weights(gp) = Q1%weights(ii) * Q1%weights(jj) * Q1%weights(kk)
                end do
            end do
        end do
    end subroutine HexahedralQuadrature

    ! ------------------------------------------------------------------
    ! 1-D Gauss-Lobatto-Legendre (GLL / spectral-element) nodes.
    ! Points include the endpoints -1 and +1.
    ! ------------------------------------------------------------------
    subroutine Spectral1DQuadrature(Quad, n)
        type(t_quadrature), intent(out) :: Quad
        integer,            intent(in)  :: n   ! number of points = order + 1

        Quad%n_points = n
        allocate(Quad%xi(n), Quad%weights(n))

        Quad%xi(1) = -1.0_dp; Quad%xi(n) = 1.0_dp
        Quad%weights(1) = 2.0_dp / (real(n-1,dp) * real(n,dp))
        Quad%weights(n) = Quad%weights(1)

        select case (n)
        case (3)
            Quad%xi(2)      = 0.0_dp
            Quad%weights(2) = 4.0_dp / 3.0_dp
        case (4)
            Quad%xi(2)      = -1.0_dp / sqrt(5.0_dp)
            Quad%xi(3)      = -Quad%xi(2)
            Quad%weights(2) = 5.0_dp / 6.0_dp
            Quad%weights(3) = Quad%weights(2)
        case (5)
            Quad%xi(2)      = -sqrt(3.0_dp/7.0_dp)
            Quad%xi(3)      = 0.0_dp
            Quad%xi(4)      = -Quad%xi(2)
            Quad%weights(2) = 49.0_dp / 90.0_dp
            Quad%weights(3) = 32.0_dp / 45.0_dp
            Quad%weights(4) = Quad%weights(2)
        case (6)
            Quad%xi(2)      = -0.765055323929465_dp
            Quad%xi(3)      = -0.285231516480645_dp
            Quad%xi(4)      = -Quad%xi(3)
            Quad%xi(5)      = -Quad%xi(2)
            Quad%weights(2) = 0.37847495629785_dp
            Quad%weights(3) = 0.55485837703549_dp
            Quad%weights(4) = Quad%weights(3)
            Quad%weights(5) = Quad%weights(2)
        case (7)
            Quad%xi(2)      = -0.830223896278567_dp
            Quad%xi(3)      = -0.468848793470714_dp
            Quad%xi(4)      = 0.0_dp
            Quad%xi(5)      = -Quad%xi(3)
            Quad%xi(6)      = -Quad%xi(2)
            Quad%weights(2) = 0.27682604736157_dp
            Quad%weights(3) = 0.43174538120986_dp
            Quad%weights(4) = 0.48761904761905_dp
            Quad%weights(5) = Quad%weights(3)
            Quad%weights(6) = Quad%weights(2)
        case (8)
            Quad%xi(2)      = -0.871740148509606_dp
            Quad%xi(3)      = -0.591700181433142_dp
            Quad%xi(4)      = -0.209299217902478_dp
            Quad%xi(5)      = -Quad%xi(4)
            Quad%xi(6)      = -Quad%xi(3)
            Quad%xi(7)      = -Quad%xi(2)
            Quad%weights(2) = 0.21070422714350_dp
            Quad%weights(3) = 0.34112269248350_dp
            Quad%weights(4) = 0.41245879465870_dp
            Quad%weights(5) = Quad%weights(4)
            Quad%weights(6) = Quad%weights(3)
            Quad%weights(7) = Quad%weights(2)
        case (9)
            Quad%xi(2)      = -0.899757995411460_dp
            Quad%xi(3)      = -0.677186279510737_dp
            Quad%xi(4)      = -0.363117463826178_dp
            Quad%xi(5)      = 0.0_dp
            Quad%xi(6)      = -Quad%xi(4)
            Quad%xi(7)      = -Quad%xi(3)
            Quad%xi(8)      = -Quad%xi(2)
            Quad%weights(2) = 0.16549536156080688_dp
            Quad%weights(3) = 0.274538712500162_dp
            Quad%weights(4) = 0.3464285109730465_dp
            Quad%weights(5) = 0.3715192743764172_dp
            Quad%weights(6) = Quad%weights(4)
            Quad%weights(7) = Quad%weights(3)
            Quad%weights(8) = Quad%weights(2)
        case (10)
            Quad%xi(2)      = -0.919533908166459_dp
            Quad%xi(3)      = -0.738773865105505_dp
            Quad%xi(4)      = -0.477924949810444_dp
            Quad%xi(5)      = -0.165278957666387_dp
            Quad%xi(6)      = -Quad%xi(5)
            Quad%xi(7)      = -Quad%xi(4)
            Quad%xi(8)      = -Quad%xi(3)
            Quad%xi(9)      = -Quad%xi(2)
            Quad%weights(2) = 0.13330599085107228_dp
            Quad%weights(3) = 0.2248893420631255_dp
            Quad%weights(4) = 0.2920426836796838_dp
            Quad%weights(5) = 0.32753976118389755_dp
            Quad%weights(6) = Quad%weights(5)
            Quad%weights(7) = Quad%weights(4)
            Quad%weights(8) = Quad%weights(3)
            Quad%weights(9) = Quad%weights(2)
        case (11)
            Quad%xi(2)      = -0.934001430408059_dp
            Quad%xi(3)      = -0.784483473663144_dp
            Quad%xi(4)      = -0.565235326996205_dp
            Quad%xi(5)      = -0.295758135586939_dp
            Quad%xi(6)      = 0.0_dp
            Quad%xi(7)      = -Quad%xi(5)
            Quad%xi(8)      = -Quad%xi(4)
            Quad%xi(9)      = -Quad%xi(3)
            Quad%xi(10)     = -Quad%xi(2)
            Quad%weights(2) = 0.10961227326699513_dp
            Quad%weights(3) = 0.18716988178030833_dp
            Quad%weights(4) = 0.24804810426402857_dp
            Quad%weights(5) = 0.2868791247790081_dp
            Quad%weights(6) = 0.3002175954556907_dp
            Quad%weights(7) = Quad%weights(5)
            Quad%weights(8) = Quad%weights(4)
            Quad%weights(9) = Quad%weights(3)
            Quad%weights(10) = Quad%weights(2)
        case (12)
            Quad%xi(2)      = -0.9448992722296681_dp
            Quad%xi(3)      = -0.8192793216440067_dp
            Quad%xi(4)      = -0.6328761530318606_dp
            Quad%xi(5)      = -0.3995309409653489_dp
            Quad%xi(6)      = -0.1365529328549276_dp
            Quad%xi(7)      = -Quad%xi(6)
            Quad%xi(8)      = -Quad%xi(5)
            Quad%xi(9)      = -Quad%xi(4)
            Quad%xi(10)     = -Quad%xi(3)
            Quad%xi(11)     = -Quad%xi(2)
            Quad%weights(2) = 0.09168451741320352_dp
            Quad%weights(3) = 0.15797470556437104_dp
            Quad%weights(4) = 0.21250841776102014_dp
            Quad%weights(5) = 0.25127560319920128_dp
            Quad%weights(6) = 0.2714052409106962_dp
            Quad%weights(7) = Quad%weights(6)
            Quad%weights(8) = Quad%weights(5)
            Quad%weights(9) = Quad%weights(4)
            Quad%weights(10) = Quad%weights(3)
            Quad%weights(11) = Quad%weights(2)
        case (13)
            Quad%xi(2)      = -0.9533098466421639_dp
            Quad%xi(3)      = -0.8463475646518723_dp
            Quad%xi(4)      = -0.6861884690817574_dp
            Quad%xi(5)      = -0.4829098210913362_dp
            Quad%xi(6)      = -0.249286930106240_dp
            Quad%xi(7)      = 0.0_dp
            Quad%xi(8)      = -Quad%xi(6)
            Quad%xi(9)      = -Quad%xi(5)
            Quad%xi(10)     = -Quad%xi(4)
            Quad%xi(11)     = -Quad%xi(3)
            Quad%xi(12)     = -Quad%xi(2)
            Quad%weights(2) = 0.07780168674682487_dp
            Quad%weights(3) = 0.13498192668960732_dp
            Quad%weights(4) = 0.1836468652035501_dp
            Quad%weights(5) = 0.2207677935661101_dp
            Quad%weights(6) = 0.2440157903066763_dp
            Quad%weights(7) = 0.2519308493334467_dp
            Quad%weights(8) = Quad%weights(6)
            Quad%weights(9) = Quad%weights(5)
            Quad%weights(10) = Quad%weights(4)
            Quad%weights(11) = Quad%weights(3)
            Quad%weights(12) = Quad%weights(2)
        case default
            stop "Spectral1DQuadrature: unsupported order."
        end select
    end subroutine Spectral1DQuadrature

    ! ------------------------------------------------------------------
    ! 2-D tensor-product GLL quadrature from a pre-computed 1-D GLL set.
    ! ------------------------------------------------------------------
    subroutine Spectral2DQuadrature(Quad, Quad1D)
        type(t_quadrature), intent(out) :: Quad
        type(t_quadrature), intent(in)  :: Quad1D

        integer :: ii, jj, k, np

        np = Quad1D%n_points
        Quad%n_points = np**2
        allocate(Quad%xi(np**2), Quad%eta(np**2), Quad%weights(np**2))

        k = 0
        do ii = 1, np
            do jj = 1, np
                k = k + 1
                Quad%xi(k)      = Quad1D%xi(ii)
                Quad%eta(k)     = Quad1D%xi(jj)
                Quad%weights(k) = Quad1D%weights(ii) * Quad1D%weights(jj)
            end do
        end do
    end subroutine Spectral2DQuadrature

    ! ------------------------------------------------------------------
    ! Level-symmetric Sn angular quadrature set.
    !
    !   dim      – 2 (four quadrants, mu > 0 hemisphere) or 3 (eight octants)
    !   sn_order – must be 2, 4, 6, 8, 12, or 16
    !   adjoint  – if .true., reverse all direction cosines
    !
    ! Weights are normalised so sum(w) = 1.
    ! Verification line printed: M0 should equal 1, M2*3 should equal 1.
    ! ------------------------------------------------------------------
    subroutine AngularQuadrature(dim, sn_order, QuadSn, adjoint)
        integer,               intent(in)    :: dim, sn_order
        type(t_sn_quadrature), intent(inout) :: QuadSn
        logical,               intent(in)    :: adjoint

        integer  :: i, j, k, m, n_levels, n_octant, q, i_sign
        integer  :: ids(3)
        real(dp), allocatable :: levels(:)
        real(dp) :: signs(8,3), m0, m2
        integer  :: n_sectors

        n_sectors = merge(4, 8, dim == 2)
        n_levels  = sn_order / 2
        n_octant  = (n_levels * (n_levels + 1)) / 2
        QuadSn%order    = sn_order
        QuadSn%n_angles = n_octant * n_sectors

        if (allocated(QuadSn%dirs))    deallocate(QuadSn%dirs)
        if (allocated(QuadSn%weights)) deallocate(QuadSn%weights)
        allocate(QuadSn%dirs(QuadSn%n_angles, 3), QuadSn%weights(QuadSn%n_angles))
        allocate(levels(n_levels))

        select case (sn_order)
        case (2);  levels = [0.57735027_dp]
        case (4);  levels = [0.35002120_dp, 0.86889030_dp]
        case (6);  levels = [0.26663550_dp, 0.68150760_dp, 0.92618080_dp]
        case (8);  levels = [0.21821790_dp, 0.57735030_dp, 0.78679580_dp, 0.95118970_dp]
        case (12); levels = [0.16721260_dp, 0.45954760_dp, 0.62802360_dp, 0.76002100_dp, 0.87227060_dp, 0.97163770_dp]
        case (16); levels = [0.13895680_dp, 0.39228930_dp, 0.53709660_dp, 0.65042640_dp, 0.74675060_dp, 0.83199660_dp, 0.90928550_dp, 0.98050090_dp]
        case default; stop "AngularQuadrature: sn_order not supported (use 2,4,6,8,12,16)."
        end select

        ! First octant
        m = 0
        do i = 1, n_levels
            do j = 1, n_levels - i + 1
                k = n_levels - i - j + 2
                m = m + 1
                QuadSn%dirs(m, :) = [levels(i), levels(j), levels(k)]

                ids = [i, j, k]
                if (ids(1) > ids(2)) call swap(ids(1), ids(2))
                if (ids(2) > ids(3)) call swap(ids(2), ids(3))
                if (ids(1) > ids(2)) call swap(ids(1), ids(2))

                select case (sn_order)
                case (2, 4)
                    QuadSn%weights(m) = 1.0_dp / 3.0_dp
                case (6)
                    if     (ids(1)==1 .and. ids(2)==1) then; QuadSn%weights(m) = 0.1761263_dp
                    else;                                    QuadSn%weights(m) = 0.1572071_dp; end if
                case (8)
                    if     (ids(1)==1 .and. ids(2)==1) then; QuadSn%weights(m) = 0.1209877_dp
                    else if(ids(1)==2 .and. ids(2)==2) then; QuadSn%weights(m) = 0.0925926_dp
                    else;                                    QuadSn%weights(m) = 0.0907407_dp; end if
                case (12)
                    if     (ids(1)==1 .and. ids(2)==1) then; QuadSn%weights(m) = 0.0707626_dp
                    else if(ids(1)==1 .and. ids(2)==2) then; QuadSn%weights(m) = 0.0558811_dp
                    else if(ids(1)==1 .and. ids(2)==3) then; QuadSn%weights(m) = 0.0373377_dp
                    else if(ids(1)==2 .and. ids(2)==2) then; QuadSn%weights(m) = 0.0502819_dp
                    else;                                    QuadSn%weights(m) = 0.0258513_dp; end if
                case (16)
                    if     (ids(1)==1 .and. ids(2)==1) then; QuadSn%weights(m) = 0.0489872_dp
                    else if(ids(1)==1 .and. ids(2)==2) then; QuadSn%weights(m) = 0.0413296_dp
                    else if(ids(1)==1 .and. ids(2)==3) then; QuadSn%weights(m) = 0.0212326_dp
                    else if(ids(1)==1 .and. ids(2)==4) then; QuadSn%weights(m) = 0.0256207_dp
                    else if(ids(1)==2 .and. ids(2)==2) then; QuadSn%weights(m) = 0.0360486_dp
                    else if(ids(1)==2 .and. ids(2)==3) then; QuadSn%weights(m) = 0.0144589_dp
                    else if(ids(1)==2 .and. ids(2)==4) then; QuadSn%weights(m) = 0.0344958_dp
                    else if(ids(1)==3 .and. ids(2)==3) then; QuadSn%weights(m) = 0.0085179_dp
                    else if(ids(1)==3 .and. ids(2)==4) then; QuadSn%weights(m) = 0.0144589_dp
                    else;                                    QuadSn%weights(m) = 0.0256207_dp; end if
                end select
            end do
        end do

        ! Remaining sectors / octants
        signs(1,:) = [ 1.0_dp,  1.0_dp,  1.0_dp]
        signs(2,:) = [-1.0_dp,  1.0_dp,  1.0_dp]
        signs(3,:) = [-1.0_dp, -1.0_dp,  1.0_dp]
        signs(4,:) = [ 1.0_dp, -1.0_dp,  1.0_dp]
        signs(5,:) = [ 1.0_dp,  1.0_dp, -1.0_dp]
        signs(6,:) = [-1.0_dp,  1.0_dp, -1.0_dp]
        signs(7,:) = [-1.0_dp, -1.0_dp, -1.0_dp]
        signs(8,:) = [ 1.0_dp, -1.0_dp, -1.0_dp]

        do q = 2, n_sectors
            do i_sign = 1, n_octant
                m = m + 1
                QuadSn%dirs(m, :)  = QuadSn%dirs(i_sign, :) * signs(q, :)
                QuadSn%weights(m)  = QuadSn%weights(i_sign)
            end do
        end do

        if (adjoint) QuadSn%dirs = -QuadSn%dirs

        QuadSn%weights = QuadSn%weights / sum(QuadSn%weights)

        m0 = sum(QuadSn%weights)
        m2 = sum(QuadSn%weights * QuadSn%dirs(:,1)**2)
        write(*,'(A,I2,A,F10.8,A,F10.8)') " S", sn_order, " Verification: M0=", m0, "  M2*3=", m2*3.0_dp

    contains
        subroutine swap(a, b)
            integer, intent(inout) :: a, b
            integer :: tmp
            tmp = a; a = b; b = tmp
        end subroutine swap
    end subroutine AngularQuadrature

end module m_quadrature
