module sph_m
  !
  ! sph.f03: A library for spherical harmonics expansion.
  !
  !   by Akira Kageyama (kage@port.kobe-u.ac.jp)
  !   on 2016.11.08
  !   at Kobe, Japan.
  ! 
  !      input : a function on a spherical surface
  !                 f(lat,lon)
  !     output : sphereical harmonics coefficients of f
  !                 ff(lm2lm(l,m)) 
  !                 where 0<=l<=NN-1, -l<=m<=l
  !                 Note: l is upto NN-1.
  !   
  !
  ! Method: 
  !   This is basically Dilts' method for the spherical 
  !   harmonics expansion:
  !     Reference: "Computation of Spherical Harmonic Expansion
  !                 Coefficients via FFT's", Gary A. Dilts,
  !                 J. Comput. Phys., Vol.57, p.439, 1985
  !
  !   A slight difference is that Fourier coeeficients of
  !   associated Legendre function, B_j^{lm} in the above
  !   paper, are numerically calculated by FFT in this program.
  !   (In the the above paper, they are analytically derived
  !   by recurrence formula.) This FFT approach is simple and
  !   accurate since associated Legendre functions are Fourier 
  !   series. I did not try to save memory for zero values of
  !   B_j^{lm}. 
  !
  !   For applying FFT in latitudinal direction to the
  !   input function f and Legendre functions Plm, their
  !   domain [0,pi] is expanded to [0,2*pi].
  !
  ! Normalization:
  !    
  !     ff_{lm} = \int_S dS f(\theta,\phi) Y_{lm}^\dagger
  !   where
  !     \dagger = complex conjugate
  !     \int_S dS = \int_0^{2\pi} dphi \int_0^\pi d\cos\theta 
  !     \int_S dS Y_{lm} Y_{ab}^\dagger = \delta_{la} \delta{mb}
  !     Y_{lm}(\theta,\phi) = C_{lm} P_l^m e^{im\phi}
  !     C_{lm} = \sqrt{(2l+1)/(4\pi) (l-m)!/(l+m)!
  !   for negative m, 
  !     P_l^{-m} = P_l^{m}
  !     C_l^{-m} = C_l^m
  !  
  ! Usage:
  !   (0) set input function f to be expanded.
  !   (1) call sph%initialize(nn,lat,lon)
  !   (2) ff = sph%expand(nn,lat,lon,f)
  ! 
  ! Dependence and a note on the size of NN: 
  !   This module uses "fft_m", which is an FFT library.
  !   Since fft_m is a simple FFT program, the data size 
  !   is limited to be power of two. Therefore, a key
  !   number NN in this program is always power of two.
  !   But this is not essential for the basic algorithm
  !   of the Dilts' method. You can replace fft_m with
  !   other FFT library so that NN can be arbitrary number.
  !
  ! History:
  !   This program is converted from my old Fortran code 
  !   "dilts.f" that was developed in July 1993.
  !   
  !   The followings are the comments:
  !!----------------
  !!     This subroutine calculates spherical harmonic coefficient.
  !!
  !!                             by KAGEYAMA, Akira
  !!                                Theory and Computer Simulation Center,
  !!                                National Institute for Fusion Science 
  !!                                <93.06.26>
  !!                      completed <93.07.09>
  !!
  !!     Reference: "Computation of Spherical Harmonic Expansion
  !!                 Coefficients via FFT's", Gary A. Dilts,
  !!                 J. Comput. Phys., Vol.57, p.439, 1985
  !!
  !!     Before you call this routine, you must calculate 2D Fourier 
  !!     coefficients f(n,m) via FFT. 
  !!     Fourier coefficients of spherical harmonic function 
  !!     bjlm are also required.
  !!                           
  !!      input: nn,f,bjlm
  !!     output: ff
  !!                           
  !!             nn: (integer) maximum mode number 
  !!                       (written as N in Dilts' paper).
  !!             In the finite difference grid mesh,
  !!                 longitudinal grid number = 2*nn
  !!                  latitudinal grid number = nn
  !!             So, longitudinal mode number = -nn:nn
  !!                  latitudinal mode number = -nn:nn (This is the trick.)
  !!----------------
  
  use fft_m
  implicit none
  private
  public :: sph__lm2lm

  integer,  parameter :: DR = selected_real_kind(15)
  real(DR), parameter :: PI = 3.1415926535897932_DR
  real(DR), parameter :: TWOPI = PI*2
  real(DR), parameter :: FOURPI = PI*4
  complex(DR), parameter :: ZERO = cmplx(0,kind=DR)
  complex(DR), parameter :: IPI  = cmplx((0,PI),kind=DR)

  type, public :: sph_t
    integer, private :: lat ! latitude
    integer, private :: lon ! longitude
    integer, private :: nn  ! max fourier mode
    integer, private :: max_mode_l  ! max Legendre mode
    logical :: initialize_done = .false.
    complex(DR), dimension(:,:), allocatable :: bjlm
      ! bjlm : 2-D Fourier coefficients of P_l^m.
    complex(DR), dimension(:,:,:), allocatable :: ylm
      ! y lm : spherical harmonics Y_l^m.
  contains
    procedure, public, pass   :: initialize => sph_t__initialize
    procedure, public, pass   :: expand => sph_t__expand_by_dilts
    procedure, public, pass   :: get_int => sph_t__get_int
    procedure, public, nopass :: plm_norm => normalized_legendre_base
  end type sph_t

contains

!
!---private---
!

  subroutine assert(condition, last_will)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: last_will
    if (.not.condition) then
      call fatal(last_will)
    end if
  end subroutine assert


  subroutine calc_fourier_transform_of_legendre(this)
    class(sph_t), intent(in out) :: this

    complex(DR), allocatable :: cwork(:)
    integer :: j, l, m, n, lm, nprim, stat, nn, nn2
    real(DR) :: theta, dtht, plm

! In order to calculate latitudinal Fourier modes of Legendre 
! function Plm, we temporally set Plm(cos(theta)), where
!          0<=l<=nn,   
!         -l<=m<=l,
!          0<=theta_j<=twopi, (j=1,2,3,...,nn2)
! The data are set on equi-spaced grid points, starting from
! the north pole.
!
    nn = this%nn
    nn2 = nn*2
    call assert(is_power_of_two(nn2), 'nn2 must be power of 2.')

    allocate(cwork(nn2),stat=stat)
    call assert(stat==0, 'alloc cwork failed.')

! --- Fourier transformation of the associated Legendre function ---
    dtht = TWOPI / nn2
    do l = 0 , nn
      do m = -l , l
        lm = sph__lm2lm(l,m)
        do j = 1 , nn2
          theta = dtht * (j-1)
          plm = extended_normalized_legendre(l,m,theta)
          cwork(j) = cmplx(plm,0.0_DR,kind=DR)
        end do  
        call fft__1d(nn2,cwork,isign=1)
        cwork(:) = cwork(:) / nn2 ! normalize
        do n = -nn+1 , nn-1
          j = merge(n+1,n+nn2+1,n>=0)
          this%bjlm(n,lm) = cwork(j)
        end do
        ! split the max mode into plus and minus nn modes. 
        this%bjlm( nn,lm) = cwork(nn+1)/2
        this%bjlm(-nn,lm) = conjg(cwork(nn+1))/2
      end do
    end do
  end subroutine calc_fourier_transform_of_legendre


  function extended_normalized_legendre(l,m,theta) result(plm)
    integer, intent(in) :: l, m 
    real(DR), intent(in) :: theta
    real(DR) :: plm

!  This function calculates "extended" (normalized) legendre function.
!  It is "extended" in two senses: 
!      (1) to accept theta>PI.
!      (2) to accept negative m.
!  
    logical  :: prerequisite

    prerequisite = l>=0 .and.  &
                   abs(m)<=l .and.  &
                   theta>=0.0_DR .and.  &
                   theta<=TWOPI
    call assert( prerequisite,   &
                 'arg(s) out of range in extended_normalized_legendre.')

    if ( theta <= PI ) then
      plm = normalized_legendre_base(l,abs(m),theta)
    else 
      plm = normalized_legendre_base(l,abs(m),TWOPI-theta)
      if (mod(m,2)/=0) then ! Note that mod(-3,2) = -1
        plm = -plm
      end if
    end if

!   if (mod(m,2)==1) then
!     plm = -plm
!   end if
  end function extended_normalized_legendre


  subroutine fatal(last_will)
    character(len=*), intent(in) :: last_will
    print *, "***<sph> fatal error: "//trim(last_will)
    stop
  end subroutine fatal


  function is_power_of_two(n) result(ans)
    integer, intent(in) :: n
    logical :: ans
    ! n = 2^k,  k = log_2(n) = log(n)/log(2)
    ans = n==2**nint(log(real(n))/log(2.0)) 
  end function is_power_of_two


  function legendre(l,m,x) result(y)
    integer,  intent(in)  :: l, m
    real(DR), intent(in)  :: x
    real(DR) :: y
!
!  This function calculates associated Legendre function ${P_l}^m$.
!  Reference: Numerical Recipe in C, Cambridge Univ. Chapter 6 p254
!        0 <= m <= l,  -1. <= x <= 1.
!
    real(DR) :: pmm, fact, somx2, pmmp1, pll
    integer  :: i, ll

    call assert( m>=0 .and. m<=l .and. abs(x)<=1.0_DR, &
               ' arg(s) out of range in legendre.')

    pmm = 1.0_DR
    if (m>0) then
      somx2 = sqrt((1.0_DR-x)*(1.0_DR+x))
      fact = 1.0_DR
      do i = 1 , m
! Note the sign before the fact. This definition follows "Recipe".
        pmm = -fact * somx2 * pmm
        fact = fact + 2.0_DR
      end do
    end if

    if (l==m) then
      y = pmm
      return   
    else
      pmmp1 = x * (2.0_DR*m+1.0_DR) * pmm
      if (l==m+1) then
        y = pmmp1
        return   
      else
        do ll = m+2 , l
          pll = ( x*(2.0_DR*ll-1.0_DR)*pmmp1  &
                 - (ll+m-1.0_DR)*pmm ) / real(ll-m,DR) 
          pmm = pmmp1
          pmmp1 = pll
        end do  
        y = pll
        return
      endif
    endif
  end function legendre


  pure function make_it_complex_and_double_in_lat(lat,lon,f) result(ff)
    integer,  intent(in) :: lat, lon
    real(DR), intent(in) :: f(lat,lon) 
    complex(DR) :: ff(2*lat-2,lon)

! north pole             equator                south pole
!   |                       |                       |
! theta=0                   pi                    twopi
!   |                       |                       |
! j=1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16  _
!   |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
!                     .   lat=9   .
!                     .   (nn=8)  .
!                     j           j_mirror

    integer :: k, j, lonh, kop, j_mirror

    do k = 1 , lon
      do j = 1 , lat
        ff(j,k) = cmplx(f(j,k),0,kind=DR)
      end do
    end do

    lonh = lon/2
    do k = 1 , lon
      kop = merge(k+lonh,k-lonh,k<=lonh) ! opposite longitude
      do j = lat+1 , lat*2-2
        j_mirror = 2*lat - j  ! (j+j_mirror)/2 = lat
        ff(j,k) = cmplx(f(j_mirror,kop),0,kind=DR)
      end do
    end do
  end function make_it_complex_and_double_in_lat


  function normalized_legendre_base(l,m,theta) result(plm)
    integer, intent(in) :: l, m ! P_l^m(cos(theta))
    real(DR), intent(in) :: theta
    real(DR) :: plm

    real(DR) :: clm
    logical  :: prerequisite

    prerequisite = m>=0 .and.  &
                   l>=0 .and.  &
                   m<=l .and.  &
                   theta>=0.0_DR .and.  &
                   theta<=PI

    call assert( prerequisite, &
                'arg(s) out of range in normalized_legendre_base.')

    clm = sqrt((2*l+1)/FOURPI)*sqrt_fact_lmm_over_fact_lpm(l,m)
    plm = clm * legendre(l,m,cos(theta))
  end function normalized_legendre_base


  function rearrange_indeces(jmax,kmax,f) result(ans)
    integer, intent(in) :: jmax, kmax
    complex(DR), intent(in) :: f(jmax,kmax)
    complex(DR) :: ans(-jmax/2:jmax/2,-kmax/2:kmax/2)
!          k   1   2   3   4   5   6   7   8=kmax
!          m   0   1   2   3   4   5   6   7
!         mm   0   1   2   3   4  -3  -2  -1
!-------------------------------------------------------
!          j   1   2   3   4   5   6   7   8=jmax
!          l   0   1   2   3   4   5   6   7
!         ll   0   1   2   3   4  -3  -2  -1
    integer :: j, k, l, m, ll, mm
    call assert( mod(jmax,2)==0 .and. mod(kmax,2)==0, &
                'jmax and kmax must be even in rearrange_indeces')
    do k = 1 , kmax
      m = k-1
      mm = merge(m,m-kmax,m<=kmax/2)
      do j = 1 , jmax
        l = j-1
        ll = merge(l,l-jmax,l<=jmax/2)
        ans(ll,mm) = f(j,k)
      end do
    end do
      
      ! Split f(max) to ans(max/2) & ans(-max/2)
    ans(jmax/2,:) = ans(jmax/2,:) / 2
    ans(:,kmax/2) = ans(:,kmax/2) / 2

    ans(-jmax/2,:) = conjg(ans(jmax/2,:))
    ans(:,-kmax/2) = conjg(ans(:,kmax/2))
  end function rearrange_indeces


  subroutine set_spherical_harmonics_ylm(this)
    class(sph_t), intent(in out) :: this

    real(DR) :: dtht, dphi, theta, phi, plm
    integer  :: j, k, l, m, lm


    associate (  nn=>this%nn, &
                lat=>this%lat, &
                lon=>this%lon)

      dtht = PI / (lat-1)
      dphi = TWOPI / lon

      do j = 1 , lat
        do k = 1 , lon
          do l = 0 , nn
            do m = -l , l
              theta = dtht * (j-1)
              phi = dphi * (k-1)
              plm = extended_normalized_legendre(l,m,theta)
              lm = sph__lm2lm(l,m)
              this%ylm(lm,j,k) = plm*exp(cmplx(0.0_DR,m*phi,kind=DR))
            end do
          end do
        end do
      end do
    end associate
  end subroutine set_spherical_harmonics_ylm


  function sqrt_fact_lmm_over_fact_lpm(l,m) result(ans)
    integer, intent(in) :: l, m
    real(DR) :: ans
! This function calculates sqrt((l-m)!/(l+m)!)  ( 0<=m<=l).
!
! let
!   ans = sqrt(g)
! where
!   g(l,m) = (l-m)!/(l+m)! 
!   ..
!   g(5,0) = (5-0)!/(5+0) = 5*4*3*2*1 / 5*4*3*2*1 = 1
!   g(5,1) = (5-1)!/(5+1) = 4*3*2*1 / 6*5*4*3*2*1 = 1 / 6*5
!   g(5,5) = (5-5)!/(5+5) = 1 / 10*9*8*7*6*5*4*3*2*1
!   ..
!   g(l,m) = 1 / (l+m)*(l+m-1)*...*(l-m+1)

    integer :: i
    call assert(l>=0 .and. m>=0 .and. m<=l, &
                'arg(s) out of rage in sqrt_fact_lmm_over_fact_lpm')
    ans = 1.0_DR
    do i = l-m+1 , l+m
      ans = ans / sqrt(real(i,DR))
    end do
  end function sqrt_fact_lmm_over_fact_lpm

!
!---public---
!

  pure function sph__lm2lm(l,m)
    integer, intent(in) :: l, m
    integer :: sph__lm2lm
    !
    ! [l]   minus<---[m]--->plus
    !  |              0
    !  |             123
    !  |            45678
    !  |           9012345
    ! \|/         678901234
    !  v
      sph__lm2lm = l*(l+1) + m
  end function sph__lm2lm


  function sph_t__expand_by_dilts(this,nn,lat,lon,input_real) result(ff)
    class(sph_t), intent(in out) :: this
    integer,  intent(in) :: nn, lat, lon
    real(DR), intent(in) :: input_real(lat,lon)
    complex(DR) :: ff(0:(nn-1)*(nn+1)) ! =ff(0:lm2lm(nn-1,nn-1))
!
! This is the main routine of this module. See comments at the
! top of this module.
!
    integer     :: j, l, m, a, amj, lm
    real(DR)    :: shift_tht, shift_phi, m1_to_amj, amjr
    complex(DR) :: add, term1, term2, term3
    complex(DR) :: fc(2*lat-2,lon)  ! Fourier modes of input
    complex(DR) ::  f(-nn:nn,-nn:nn) ! reordered-and-shifted fc

    call assert( lon==2*(lat-1), &
                'lon and lat inconsistent. <sph_t__expand_by_dilts>')
    call assert( is_power_of_two(nn), &
                'nn must be power of two. <sph_t__expand_by_dilts>')
    call assert( nn*2==lon, &
                'nn and lon inconsistent. <sph_t__expand_by_dilts>')
    call assert( this%initialize_done, &
                'You forgot initialization <sph_t__expand_by_dilts>')

    fc = make_it_complex_and_double_in_lat(lat,lon,input_real) 

      ! complex FFT in both (lat & lon) directions.
    call fft__2d(2*lat-2,lon,fc,switch=1)
      ! normalization of FFT
    fc(:,:) = fc(:,:) / lon / (2*lat-2)  

    f(:,:) = rearrange_indeces(2*lat-2,lon,fc)

    associate ( nn=>this%nn )
      do l = 0 , this%max_mode_l
        do m = -l , l
          lm = sph__lm2lm(l,m)
          add = ZERO
          do j = -l , l
            !--- 1st term ---
            term1 = ZERO
            do a = -nn , j-2
              amj = a-j
              amjr = real(amj,kind=DR) 
              m1_to_amj = real((-1)**amj,kind=DR)
              term1 = ((m1_to_amj+1)/(amjr**2-1)) * f(a,m) + term1
            end do  
            !--- 2nd term ---
            if (j==nn) then
              term2 = IPI*( f(j-1,m) -   ZERO   ) - 4*f(j,m)
            else if (j==-nn) then
              term2 = IPI*(    ZERO  - f(j+1,m) ) - 4*f(j,m)
            else
              term2 = IPI*( f(j-1,m) - f(j+1,m) ) - 4*f(j,m)
            endif
            !--- 3rd term ---
            term3 = ZERO
            do a = j+2 , nn
              amj = a-j
              amjr = real(amj,kind=DR) 
              m1_to_amj = real((-1)**amj,kind=DR)
              term3 = ((m1_to_amj+1)/(amjr**2-1)) * f(a,m) + term3
            end do  
            !--- sum of the three terms ---
            add = add  &
                + (2*term1+term2+2*term3)*conjg(this%bjlm(j,lm)) 
          end do  
          ff(lm) = -PI * add
        end do
      end do
    end associate
  end function sph_t__expand_by_dilts


  function sph_t__get_int(this,which) result(ans)
    class(sph_t), intent(in out) :: this
    character(len=*), intent(in) :: which
    integer :: ans
    select case (which)
      case ('nn')
        ans = this%nn
      case ('lat')
        ans = this%lat
      case ('lon') 
        ans = this%lon
      case ('max_mode_l')
        ans = this%max_mode_l
    end select
  end function sph_t__get_int


  subroutine sph_t__initialize(this,nn,lat,lon)
    class(sph_t), intent(in out) :: this
    integer, intent(in) :: nn, lat, lon 

    integer :: lmmax, stat
    logical :: prerequisite

    prerequisite = lat-1==nn .and. &
                   lon==2*nn .and. &
                   is_power_of_two(nn)
    call assert(prerequisite, "bad params of nn/lat/lon.") 

    this%nn  = nn
    this%lat = lat
    this%lon = lon
    this%max_mode_l = nn-1 ! Limit inherent in the Dilts' method.

    lmmax = sph__lm2lm(nn,nn)
    allocate(this%bjlm(-nn:nn,0:lmmax),stat=stat)
    call assert(stat==0, &
               'alloc this%bjlm failed. in sph_t__initialize.')
    allocate(this%ylm(0:lmmax,lat,lon),stat=stat)
    call assert(stat==0, &
               'alloc this%ylm failed. in sph_t__initialize.')

    call calc_fourier_transform_of_legendre(this)
    call set_spherical_harmonics_ylm(this)

    this%initialize_done = .true.
  end subroutine sph_t__initialize

end module sph_m


program test
  use sph_m
  implicit none

  integer,  parameter :: DR = selected_real_kind(15)
  integer,  parameter :: SP = kind(1.0)
  real(DR), parameter :: PI = 3.1415926535897932_DR
  real(DR), parameter :: TWOPI = PI*2

  integer, parameter :: NR = 124
  integer, parameter :: NT = 162
  integer, parameter :: NP = 322
  !integer, parameter :: NR = 2**6
  !integer, parameter :: NT = NR+1 
  !integer, parameter :: NP = NR*2
  real(SP) :: sfield(NR,NT,NP)
  integer, parameter :: R_INDEX = nint(real(NR)*0.5_DR)
  
  type(sph_t) :: sph
  integer, parameter :: NN=2**6 ! power of 2.
  integer, parameter :: LAT=NN+1
  integer, parameter :: LON=NN*2
  integer :: j, k, l, m, lm, max_l
  !real(DR) :: dtht, dphi, theta, phi, sum
  real(DR) :: dtht, dphi, theta, phi, sum, cull_rate_lat, cull_rate_lon
  !real(DR) :: theta(LAT_RAW), phi(LON_RAW)
  real(DR) :: error, max_f, max_s, max_e
  real(DR) :: f_raw(NT,NP), f(LAT, LON)
  real(DR) :: legendre 
  complex(DR), allocatable :: ff(:)
  character(2) s

  dtht = PI / (LAT-1)
  dphi = TWOPI / LON

  print *, "### sph%initialize started"
  call sph%initialize(nn=NN,lat=LAT,lon=LON)
  
  ! --- read file ---
  print *, "### read_file started"
  !open(11, file="./input/rad0.9@eMar03c.005.sfield.vr.n001575375.t00190", status="old")
  open(11, file="input/eMar03c.005.sfield.br.n001575375.t00190", form="unformatted", status="old")
  read(11)   sfield
  close(11)
  print *, "R_INDEX = ", R_INDEX

  do k = 1 , NP
    do j = 1 , NT
!     f_raw(j,k) = sfield(R_INDEX,j,k)
      theta = dtht*(j-1)
      phi   = dphi*(k-1)
!     call ut__legendren(10,3,cos(theta),legendre)
!     f_raw(j,k) = legendre*cos(3*phi)
      call ut__legendren(1,0,cos(theta),legendre)
      f_raw(j,k) = legendre*cos(0*phi) 
!     call ut__legendren(5,3,cos(theta),legendre)
!     f(j,k) = f(j,k) + legendre*cos(3*phi) * 0.3_DR
!     f(j,k) = cos(theta)*cos(m*phi)
!     f(j,k) = cos(theta)*cos(3*phi)
!     f_raw(j,k) = 1.0_DR
      !print *, f_raw(j,k)
    end do
  end do
  print *, "### read_file ended"

  print *, "### interpolate_data start"
  call liner_interpolation(f_raw,f)
  print *, "### interpolate_data ended"

  ! --- cull data to apply NN=2**x  ---
!  print *, "### cull data started"
!  cull_rate_lat = real(NT)/real(LAT)   ! NT > LAT
!  cull_rate_lon = real(NP)/real(LON)   ! NP > LON
!  print *, "cull_rate_lat, cull_rate_lon", cull_rate_lat, cull_rate_lon
!  do k = 1, LON
!    do j = 1, LAT
!      f(j,k) = f_raw(nint(cull_rate_lat*j),nint(cull_rate_lon*k))
!      !print *, f(j,k)
!    end do
!  end do
!  print *, "### cull data ended"

  max_l = sph%get_int("max_mode_l")
  allocate(ff(0:sph__lm2lm(max_l,max_l)))

  print *, "### expand_data start"
  ff = sph%expand(NN,LAT,LON,f)
  print *, "### expand_data ended"

  open(12, file="./output/expanded3D.dat", status="replace")
  open(13, file="./output/expandedM0.dat", status="replace")
  open(14, file="./output/expandedM1.dat", status="replace")
  open(15, file="./output/expandedM2.dat", status="replace")
  open(16, file="./output/expandedM3.dat", status="replace")
  open(17, file="./output/expandedM4.dat", status="replace")
  open(18, file="./output/expandedM5.dat", status="replace")
  do l = 0 , max_l
    do m = -l , l
      !write(s,'(i2)') m
      !open(12, file="./output/expandedM=" // s // ".dat", status="replace")
      lm = sph__lm2lm(l,m)
      !if ( abs(ff(lm))>1.e-5 ) then
        print *, l,m,ff(lm),abs(ff(lm))
        if(m==0) then
         write(12,*) l, m, abs(ff(lm))
        end if
        if(m==1) then
         write(13,*) l, abs(ff(lm))
        end if
        if(m==2) then
         write(14,*) l, abs(ff(lm))
        end if
        if(m==3) then
         write(15,*) l, abs(ff(lm))
        end if
        if(m==4) then
         write(16,*) l, abs(ff(lm))
        end if
        if(m==5) then
         write(17,*) l, abs(ff(lm))
        end if
        if(m==6) then
         write(18,*) l, abs(ff(lm))
        end if
      !end if
      !close(12)
!      if ( abs(ff(lm))>1.e-5 ) then
!       print *,l,m,ff(lm)
!      end if
    end do
    write(12,*)
  end do
  close(12)

!   ! --- test 02 ---
!   l = 9
!   m = 2
!   do k = 1 , LON
!     do j = 1 , LAT
!       theta = dtht*(j-1)
!       phi   = dphi*(k-1)
!       f(j,k) = sph%plm_norm(l,m,theta)*sin(m*phi)
!     end do
!   end do
! 
!   ff = sph%expand(NN,LAT,LON,f)
! 
!   do l = 0 , max_l
!     do m = -l , l
!       lm = sph__lm2lm(l,m)
!       if ( abs(ff(lm))>1.e-5 ) then
!         print *,l,m,ff(lm)
!       end if
!     end do
!   end do

  ! Check if input is recovered by inverse spherical exp.
  max_f = 0.0_DR
  max_s = 0.0_DR
  max_e = 0.0_DR
  do j = 1 , LAT
    do k = 1 , lon
      sum = 0.0_DR
      do l = 0 , max_l
        do m = -l , l
          lm = sph__lm2lm(l,m)
          sum = sum + ff(lm)*sph%ylm(lm,j,k)
        end do
      end do
      error = f(j,k) - sum
      max_f = max(max_f, abs(f(j,k)))
      max_s = max(max_s, abs(sum))
      max_e = max(max_e, abs(error))
    end do
  end do
  print *, 'max of f, s, error = ', max_f, max_s, max_e


contains
  subroutine liner_interpolation(f_,f)
    real(DR), intent(in)  :: f_(:,:)
    real(DR), intent(out) :: f(:,:)
    !real, parameter :: PI=3.14159265
    integer :: j_,k_,j,k,count_j,count_k
    integer :: north,south,west,east
    real(DR) :: dt_,dp_,dt,dp
    real(DR) :: theta, phi
    real(DR) :: theta_1, phi_1, theta_2, phi_2
    real(DR) :: alpha, beta

    dt_ = PI/(NT-1)
    dp_ = TWOPI/NP
    dt  = PI/(LAT-1)
    dp  = TWOPI/LON
    count_j = 1
    count_k = 1

    do j=1, LAT
      do k=1, LON
        theta = dt*real(j-1)
        phi   = dp*real(k-1)
        north=-1
        west=-1
        ! find neighbors
        do j_ = 1, NT
          theta_1 = dt_*real(j_-1)
          theta_2 = dt_*real(j_)
          if(theta_1 <= theta .and. &
             theta <= theta_2) then
             alpha = (theta-theta_1)/dt_
             north=j_
             south=j_+1
             exit
          end if
        end do
        do k_ = 1, NP
          phi_1 = dp_*real(k_-1)
          phi_2 = dp_*real(k_)
          if(phi_1 <= phi .and. &
             phi <= phi_2) then
             beta = (phi-phi_1)/dp_
             west=k_
             east=k_+1
             exit
          end if
        end do
        if(north==-1 .or. west==-1) then
          print *,"error: Not find neighbors"
        end if
        !print *, "LAT=", j, "LON=", k
        !print *, "alpha=", alpha, "beta=", beta
        !print *, "north=", north
        !print *, "south=", south
        !print *, "west=", west
        !print *, "east=", east
        ! caluculate f
        f(j,k)=f_(north, west)*(1-alpha)*(1-beta) &
              +f_(south, west)*alpha*(1-beta)     &
              +f_(north, east)*(1-alpha)*beta     &
              +f_(south, east)*alpha*beta
      end do
    end do
  end subroutine liner_interpolation
!_______________________________________________________public__
!
  subroutine ut__legendren(l,m,x,y)
    integer, intent(in) :: l, m
    real(DR), intent(in) :: x
    real(DR), intent(out) :: y
!_______________________________________________________________
!
! - "legendren" stands for Legendr-normalized.
!
! - Calculates the normalized associated Legendre function
!                                          $C_{lm}{P_l}^m$.
! - Reference: Numerical Recipe in C, Chapter 6 p254
!
!       0 <= m <= l,  -1. <= x <= 1.
!_______________________________________________________________/
!

    real(DR) :: pmm, fact, fact1, fact2, pll, pmmp1, somx2
    integer :: ll, i

    pmm = 1.0_DR

    if (m>0) then
       somx2 = sqrt((1.0_DR-x)*(1.0_DR+x))
       do i = 1, m
          ! Note the sign before the fact. This definition follows "Recipe".
          fact = sqrt(real(2*i-1,DR)/real(2*i,DR))
          pmm = -fact * somx2 * pmm
       end do
       pmm = pmm * sqrt(real(2*m+1,DR)/(4*PI))
    else
       pmm = sqrt(1.0_DR/(4*PI))
    end if

    if (l==m) then
       y = pmm
    else
       pmmp1 = x * sqrt(real(2*m+3,DR)) * pmm
       if (l==m+1) then
          y = pmmp1
       else
          do ll = m+2, l
             fact1 = sqrt( (real(2*ll+1,DR)*real(ll-m,DR)) &
                         / (real(2*ll-1,DR)*real(ll+m,DR)) &
                         )
             fact2 = sqrt( (real(2*ll+1,DR)*real(ll-m,DR)*real(ll-m-1,DR)) &
                         / (real(2*ll-3,DR)*real(ll+m,DR)*real(ll+m-1,DR)) &
                         )
             pll = ( x*real(2*ll-1,DR)*fact1*pmmp1  &
                     - real(ll+m-1,DR)*fact2*pmm ) / real(ll-m,DR)
             pmm = pmmp1
             pmmp1 = pll
          end do
          y = pll
       endif
    endif

  end subroutine ut__legendren

  function test_plm(x,l,m) result(y)
    real(DR), intent(in) :: x
    real(DR) :: y, fact
    integer  :: i, l, m
    ! P_l^m = P_10^3
     y = 6435.0_DR/16*Sqrt(1-x**2)*(-1+x**2)*(-7*x+105*x**3-357*x**5+323*x**7)

    fact = 1.0_DR
    do i = l-m+1 , l+m
      fact = fact * real(i,DR)
    end do
    fact = sqrt((2*l+1)/(PI*4)/fact)
    y = fact*y
  end function test_plm
end program test
