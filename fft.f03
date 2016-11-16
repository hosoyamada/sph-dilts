module fft_m
  !
  ! FFT module; an FFT library.
  !   - Developed by Akira Kageyama (kage@port.kobe-u.ac.jp)
  !   - on 2016.11.07
  !   - Converted from my old library dcft2.f.
  ! Note:
  !   - Normalization is not applied. Apply 1/N for fft__1d
  !     or 1/(N1*N2) for fft__2d, if necessary.
  !
  ! The followings are the old comments for dcft2.f
  !
  !!<>c<>c<>c<>c<>c<>c<>c<>c<>c<>c<>c<>c<>c<>c<>c<>c<>c<>c<>c<>c<>c<>c<>c<>
  !!
  !! This is dcft2.f, used in my old dynamo code.
  !!            2008.04.07: Akira Kageyama (kage@jamstec.go.jp)
  !!
  !!<>c<>c<>c<>c<>c<>c<>c<>c<>c<>c<>c<>c<>c<>c<>c<>c<>c<>c<>c<>c<>c<>c<>c<>
  !!
  !! this is a (1-d) fft routine which transforms, simultaneously, ll-set
  !!  of nn-dimensional double precision complex data 'data'.
  !!
  !! nn must be a power of 2. this is not checked for.
  !!
  !! vectorization is taken in the ll-loop: so this routine is effective
  !!  when you have many same (nn) dimensional data to be transformed,
  !!  i.e., ll>nn.
  !!
  !!        data( 1,real(1))
  !!        data( 2,real(1))
  !!              !
  !!        data(ll,real(1))
  !!        data( 1,imag(1))          do 1  n = 1 , nn
  !!              !                     do 2  ndata = 1 , ll
  !!        data(ll,imag(1))              data(ndata,2*n-1) = re(n-th data)
  !!        data( 1,real(2))              data(ndata,2*n)   = im(n-th data)
  !!              !                 2   continue
  !!        data(ll,real(2))        1 continue
  !!        data( 1,imag(2))
  !!              !
  !!        data(ll,imag(2))
  !!        data( 1,real(3))
  !!              !
  !!              !
  !!              !
  !!        data( 1,real(nn))
  !!              !
  !!        data(ll,real(nn))
  !!        data( 1,imag(nn))
  !!              !
  !!        data(ll,imag(nn))
  !!
  !! data is replaced by nn times of its fourier transform (isign=1)
  !!  or is replaced by its inverse fourier transform     (isign=-1).
  !! so you must devide data by nn, if necessary.
  !!
  !! 1-d fft used in this routine is based on danielson-lanczos algorithm.
  !! referece:  "numerical recipe in c", press et al.
  !!  cambridge university press, 1988, p412
  !!
  !!  by kageyama, a    1992.04.04
  !!
  
  implicit none
  private
  public :: fft__1d, fft__2d

  interface fft__1d
    module procedure fft_1d_complex_multi, &
                     fft_1d_complex_single, &
                     fft_1d_real_multi, &
                     fft_1d_real_single 
  end interface fft__1d

  interface fft__2d
    module procedure fft_2d_complex_multi
  end interface fft__2d

  integer,  parameter :: DR = selected_real_kind(15)
  real(DR), parameter :: PI = 3.1415926535897932_DR
  real(DR), parameter :: TWOPI = PI*2

contains

! Private

  subroutine fft_1d_complex_multi(ll,nn,fc,isign)
    !  n=1  2  3  4  5  6  7  8    (nn=8)
    !        . . 
    !     .        .  <--- f(1,n)
    !    +--+--+--+--.--+--+--+--+
    !                  .        .
    !                      . .
    !    
    integer, intent(in) :: ll    ! number of data
    integer, intent(in) :: nn    ! size of data
    complex(DR), dimension(ll,nn), intent(in out) :: fc
    integer, intent(in) :: isign ! forward (+1) or inverse (-1)

    real(DR), dimension(ll,2*nn) :: fr

    integer :: n, l
   
    ! complex ==> real
    do n = 1 , nn
      do l = 1 , ll
        fr(l,2*n-1) =  real(fc(l,n))
        fr(l,2*n  ) = aimag(fc(l,n))
      end do
    end do

    call fft_1d_real_multi(ll,2*nn,fr,isign)

    ! real ==> complex
    do n = 1 , nn
      do l = 1 , ll
        fc(l,n) = cmplx(fr(l,2*n-1),  &  ! real part
                        fr(l,2*n  ),  &  ! imag part
                        DR)
      end do
    end do
  end subroutine fft_1d_complex_multi


  subroutine fft_1d_complex_single(nn,fc,isign)
    !  n=1  2  3  4  5  6  7  8    (nn=8)
    !        . . 
    !     .        .  <--- f(n)
    !    +--+--+--+--.--+--+--+--+
    !                  .        .
    !                      . .
    !    
    integer, intent(in) :: nn    ! size of data
    complex(DR), dimension(nn), intent(in out) :: fc
    integer, intent(in) :: isign ! forward (+1) or inverse (-1)

    real(DR), dimension(2*nn) :: fr

    integer :: n
   
    ! complex ==> real
    do n = 1 , nn
      fr(2*n-1) = real(fc(n))
      fr(2*n  ) = aimag(fc(n))
    end do

    call fft_1d_real_single(2*nn,fr,isign)

    ! real ==> complex
    do n = 1 , nn
      fc(n) = cmplx(fr(2*n-1),  &  ! real part
                    fr(2*n  ),  &  ! imag part
                    DR) 
    end do
  end subroutine fft_1d_complex_single
  

  subroutine fft_1d_real_multi(ll,two_n,data,isign)
    integer, intent(in) :: ll, two_n
    real(DR), dimension(ll,two_n), intent(in out) :: data
    integer, intent(in) :: isign
   
    integer  :: i, j, l, m, mmax, istep
    real(DR) :: theta, wpr, wpi, wr, wi
    real(DR) :: tempr, tempi, wtemp

    j = 1
    do i = 1 , two_n , 2
      if (j>i) then
        do l = 1 , ll
                tempr = data(l,j  )
                tempi = data(l,j+1)
          data(l,j  ) = data(l,i  )
          data(l,j+1) = data(l,i+1)
          data(l,i  ) = tempr
          data(l,i+1) = tempi
        end do
      endif
      m = two_n/2
      do while ( (m>=2).and.(j>m) )
        j = j-m
        m = m/2
      end do
      j = j + m
    end do
  
    mmax = 2
    do while ( two_n>mmax )
      istep = 2*mmax
      theta = -TWOPI/real(isign*mmax,DR)
      wpr   = -2.0_DR*sin(0.5_DR*theta)**2
      wpi   = sin(theta)
      wr    = 1.0_DR
      wi    = 0.0_DR
      do m = 1 , mmax , 2
        do i = m , two_n , istep
          j = i + mmax
          do l = 1 , ll
                  tempr = wr*data(l,j  )-wi*data(l,j+1)
                  tempi = wr*data(l,j+1)+wi*data(l,j  )
            data(l,j  ) = data(l,i  ) - tempr
            data(l,j+1) = data(l,i+1) - tempi
            data(l,i  ) = data(l,i  ) + tempr
            data(l,i+1) = data(l,i+1) + tempi
          end do
        end do
        wtemp = wr
        wr = wr*wpr -    wi*wpi + wr
        wi = wi*wpr + wtemp*wpi + wi
      end do
      mmax = istep
    end do
  end subroutine fft_1d_real_multi
  

  subroutine fft_1d_real_single(two_n,data,isign)
    integer, intent(in) :: two_n
    real(DR), dimension(two_n), intent(inout) :: data
    integer, intent(in) :: isign
  !
  ! this is a (1-d) fft routine which transforms, simultaneously, ll-set
  !  of nn-dimensional double precision complex data 'data'.
  !
  ! nn must be a power of 2. this is not checked for.
  !
  ! vectorization is taken in the ll-loop: so this routine is effective
  !  when you have many same (nn) dimensional data to be transformed,
  !  i.e., ll>nn.
  !
  !        data(real(1))
  !        data(real(1))
  !              !
  !        data(real(1))
  !        data(imag(1))          
  !              !                  do n = 1 , nn
  !        data(imag(1))              data(ndata,2*n-1) = re(n-th data)
  !        data(real(2))              data(ndata,2*n)   = im(n-th data)
  !              !                  end do
  !        data(real(2))        
  !        data(imag(2))
  !              !
  !        data(imag(2))
  !        data(real(3))
  !              !
  !              !
  !              !
  !        data(real(nn))
  !              !
  !        data(real(nn))
  !        data(imag(nn))
  !              !
  !        data(imag(nn))
  !
  ! data is replaced by nn times of its fourier transform (isign=1)
  !  or is replaced by its inverse fourier transform     (isign=-1).
  ! so you must devide data by nn, if necessary.
  !
  ! 1-d fft by Danielson-Lanczos algorithm.
  ! referece: "numerical recipe in c", press et al.
  !            cambridge university press, 1988, p412
  !
  !  by kageyama, a    1992.04.04
  !
    integer  :: i, j, m, mmax, istep
    real(DR) :: theta, wpr, wpi, wr, wi
    real(DR) :: tempr, tempi, wtemp

    j = 1
    do i = 1 , two_n , 2
      if (j>i) then
              tempr = data(j)
              tempi = data(j+1)
          data(j)   = data(i)
          data(j+1) = data(i+1)
          data(i)   = tempr
          data(i+1) = tempi
      endif
      m = two_n/2
      do while ( (m>=2).and.(j>m) )
        j = j-m
        m = m/2
      end do
      j = j + m
    end do
  
    mmax = 2
    do while ( two_n>mmax )
      istep = 2*mmax
      theta = -TWOPI/real(isign*mmax,DR)
      wpr   = -2.0_DR*sin(0.5_DR*theta)**2
      wpi   = sin(theta)
      wr    = 1.0_DR
      wi    = 0.0_DR
      do m = 1 , mmax , 2
        do i = m , two_n , istep
          j = i + mmax
              tempr = wr*data(j  )-wi*data(j+1)
              tempi = wr*data(j+1)+wi*data(j  )
          data(j  ) = data(i  ) - tempr
          data(j+1) = data(i+1) - tempi
          data(i  ) = data(i  ) + tempr
          data(i+1) = data(i+1) + tempi
        end do
        wtemp = wr
        wr = wr*wpr -    wi*wpi + wr
        wi = wi*wpr + wtemp*wpi + wi
      end do
      mmax = istep
    end do
  end subroutine fft_1d_real_single


  subroutine fft_2d_complex_multi(n1,n2,f,switch)
    integer, intent(in) :: n1, n2
    complex(DR), intent(in out) :: f(n1,n2)
    integer, intent(in) :: switch
    complex(DR) :: work(n2,n1)

    integer :: i, j

    ! --- FFT in n2 direction ---
    call fft_1d_complex_multi(n1,n2,f,switch)

    do j = 1 , n2
      do i = 1 , n1
        work(j,i) = f(i,j)
      end do
    end do 

    ! --- FFT in n1 direction ---
    call fft_1d_complex_multi(n2,n1,work,switch)

    do j = 1 , n2
      do i = 1 , n1
        f(i,j) = work(j,i)
      end do
    end do 
  end subroutine fft_2d_complex_multi

end module fft_m


!program test
!  use fft_m
!  implicit none
!
!  integer,  parameter :: DR = selected_real_kind(15)
!  real(DR), parameter :: PI = 3.1415926535897932_DR
!  real(DR), parameter :: TWOPI = PI*2
!
!  integer :: isign, n, l
!  integer, parameter :: nn = 32
!  integer, parameter :: ll = 64
!  real(DR) :: tht, dtht = TWOPI / ll
!  real(DR) :: phi, dphi = TWOPI / nn
!  real(DR), dimension(ll,2*nn) :: data
!  complex(DR), dimension(ll,nn) :: c
!  integer :: mode, mode_l, mode_m
!
!  if (.true.) then
!    mode_l = 8
!    do n = 1 , nn
!      phi = dphi*(n-1)
!      do l = 1 , ll
!        c(l,n) = cmplx(cos(mode_l*phi), sin(mode_l*phi), DR)
!      end do
!    end do
!
!    call iTest00
!  end if
!
!  if (.false.) then
!    do n = 1 , nn
!      do l = 1 , ll
!  !     data(l,2*n-1) = re(n-th data)
!  !     data(l,2*n)   = im(n-th data)
!        mode = l
!        phi = dphi*(n-1)
!        data(l,2*n-1) = cos(mode*phi)
!        data(l,2*n)   = sin(mode*phi)
!      end do
!    end do
!
!    call iTest01
!  end if
!
!  if (.true.) then
!    mode_l = -9
!    mode_m = -21
!    do n = 1 , nn
!      do l = 1 , ll
!        tht = dtht*(l-1)
!        phi = dphi*(n-1)
!        c(l,n) = exp(cmplx(0.0_DR,mode_l*tht,DR)) &
!               * exp(cmplx(0.0_DR,mode_m*phi,DR)) 
!      end do
!    end do
!    call iTest02
!  end if
!
!contains
!
!  subroutine iTest00
!    isign = 1
!    call fft__1d(ll,nn,c,isign)
!    c(:,:) = c(:,:) / nn
!
!    print *, '------ iTest00 ------'
!    print *, ' l, mode_n, c'
!    do l = 1 , ll
!      do n = 1 , nn
!        if ( abs(c(l,n)) > 1.e-10 ) then
!          print *, l, n-1, c(l,n)
!        end if
!      end do
!    end do
!  end subroutine iTest00
!
!  subroutine iTest01
!    isign = 1
!    call fft__1d(ll,nn,data,isign)
!    data(:,:) = data(:,:) / nn
!
!    do n = 1 , nn
!      do l = 1 , ll
!        if ( abs(data(l,n)) > 1.e-10 ) then
!          print *, l, n, data(l,n)
!        end if
!      end do
!    end do
!  end subroutine iTest01
!
!  subroutine iTest02
!    isign = 1
!    call fft__2d(ll,nn,c,isign)
!    c(:,:) = c(:,:) / (nn*ll)
!
!    print *,' ---- '
!    print *,' mode_l, mode_m, c'
!    do n = 1 , nn
!      do l = 1 , ll
!        if ( abs(c(l,n)) > 1.e-10 ) then
!          print *, l-1, n-1, c(l,n)
!        end if
!      end do
!    end do
!  end subroutine iTest02
!
!end program test

