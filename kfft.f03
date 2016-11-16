module kfft_m
  !
  ! Kobe FFT module; an FFT library.
  !   - Developed by Akira Kageyama (kage@port.kobe-u.ac.jp)
  !   - on 2016.11.07
  !   - Converted from my old library dcft2.f.
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
  
  implicit none
  private
  public :: kfft

  interface kfft
    module procedure fft_complex_multi, &
                     fft_complex_single, &
                     fft_multi, &
                     fft_single 
  end interface kfft

  integer,  parameter :: DR = selected_real_kind(15)
  real(DR), parameter :: PI = 3.1415926535897932_DR
  real(DR), parameter :: TWOPI = PI*2

contains

! Private

  subroutine normalize_multi(ll,nn,data)
    integer, intent(in) :: ll, nn
    real(DR), dimension(ll,2*nn), intent(inout) :: data

    real(DR) :: norm
    norm = 1.0_DR / nn
    data(:,:) = data(:,:) * norm
  end subroutine normalize_multi


  subroutine normalize_single(nn,data)
    integer, intent(in) :: nn
    real(DR), dimension(2*nn), intent(inout) :: data

    real(DR) :: norm
    norm = 1.0_DR / nn
    data(:) = data(:) * norm
  end subroutine normalize_single


  subroutine fft_complex_multi(ll,nn,fc,isign)
    !  n=1  2  3  4  5  6  7  8    (nn=8)
    !        . . 
    !     .        .  <--- f(1,n)
    !    +--+--+--+--.--+--+--+--+
    !                  .        .
    !                      . .
    !    
    integer, intent(in) :: ll    ! number of data
    integer, intent(in) :: nn    ! size of data
    complex(DR), dimension(ll,nn), intent(inout) :: fc
    integer, intent(in) :: isign ! forward (+1) or inverse (-1)

    real(DR), dimension(ll,2*nn) :: fr
    real(DR) :: norm

    integer :: n, l
   
    ! complex ==> real
    do n = 1 , nn
      do l = 1 , ll
        fr(l,2*n-1) = real(fc(l,n))
        fr(l,2*n  ) = aimag(fc(l,n))
      end do
    end do

    call fft_multi(ll,nn,fr,isign)

    if ( isign == 1 ) then
      call normalize_multi(ll,nn,fr)
    end if 

    ! real ==> complex
    do n = 1 , nn
      do l = 1 , ll
        fc(l,n) = cmplx(fr(l,2*n-1),  &  ! real part
                        fr(l,2*n  ),  &  ! imag part
                        DR)
      end do
    end do

  end subroutine fft_complex_multi


  subroutine fft_complex_single(nn,fc,isign)
    !  n=1  2  3  4  5  6  7  8    (nn=8)
    !        . . 
    !     .        .  <--- f(n)
    !    +--+--+--+--.--+--+--+--+
    !                  .        .
    !                      . .
    !    
    integer, intent(in) :: nn    ! size of data
    complex(DR), dimension(nn), intent(inout) :: fc
    integer, intent(in) :: isign ! forward (+1) or inverse (-1)

    real(DR), dimension(2*nn) :: fr
    real(DR) :: norm

    integer :: n
   
    ! complex ==> real
    do n = 1 , nn
      fr(2*n-1) = real(fc(n))
      fr(2*n  ) = aimag(fc(n))
    end do

    call fft_single(nn,fr,isign)

    if ( isign == 1 ) then
      call normalize_single(nn,fr)
    end if 

    ! real ==> complex
    do n = 1 , nn
        fc(n) = cmplx(fr(2*n-1),  &  ! real part
                      fr(2*n  ),  &  ! imag part
                      DR)
    end do

  end subroutine fft_complex_single
  

  subroutine fft_multi(ll,nn,data,isign)
    integer, intent(in) :: ll, nn
    real(DR), dimension(ll,2*nn), intent(inout) :: data
    integer, intent(in) :: isign
   
    integer  :: n2, i, j, l, m, mmax, istep
    real(DR) :: theta, wpr, wpi, wr, wi
    real(DR) :: tempr, tempi, wtemp

    n2 = 2*nn
    j = 1
  
    do i = 1 , n2 , 2
      if (j>i) then
        do l = 1 , ll
              tempr   = data(l,j)
              tempi   = data(l,j+1)
          data(l,j)   = data(l,i)
          data(l,j+1) = data(l,i+1)
          data(l,i)   = tempr
          data(l,i+1) = tempi
        end do
      endif
      m = nn
      do while ( (m>=2).and.(j>m) )
        j = j-m
        m = m/2
      end do
      j = j + m
    end do
  
    mmax = 2
  
    do while ( n2>mmax )
      istep = 2*mmax
      theta = -TWOPI/real(isign*mmax,DR)
      wpr   = -2.0_DR*sin(0.5_DR*theta)**2
      wpi   = sin(theta)
      wr    = 1.0_DR
      wi    = 0.0_DR
      do m = 1 , mmax , 2
        do i = m , n2 , istep
          j = i + mmax
          do l = 1 , ll
            tempr = wr*data(l,j)-wi*data(l,j+1)
            tempi = wr*data(l,j+1)+wi*data(l,j)
            data(l,j)   = data(l,i) - tempr
            data(l,j+1) = data(l,i+1) - tempi
            data(l,i)   = data(l,i) + tempr
            data(l,i+1) = data(l,i+1) + tempi
          end do
        end do
        wtemp = wr
        wr = wr*wpr -    wi*wpi + wr
        wi = wi*wpr + wtemp*wpi + wi
      end do
      mmax = istep
    end do

  end subroutine fft_multi
  

  subroutine fft_single(nn,data,isign)
    integer, intent(in) :: nn
    real(DR), dimension(2*nn), intent(inout) :: data
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
    integer  :: n2, i, j, l, m, mmax, istep
    real(DR) :: theta, wpr, wpi, wr, wi
    real(DR) :: tempr, tempi, wtemp

    n2 = 2*nn
    j = 1
  
    do i = 1 , n2 , 2
      if (j>i) then
              tempr   = data(j)
              tempi   = data(j+1)
          data(j)   = data(i)
          data(j+1) = data(i+1)
          data(i)   = tempr
          data(i+1) = tempi
      endif
      m = nn
      do while ( (m>=2).and.(j>m) )
        j = j-m
        m = m/2
      end do
      j = j + m
    end do
  
    mmax = 2
  
    do while ( n2>mmax )
      istep = 2*mmax
      theta = -TWOPI/real(isign*mmax,DR)
      wpr   = -2.0_DR*sin(0.5_DR*theta)**2
      wpi   = sin(theta)
      wr    = 1.0_DR
      wi    = 0.0_DR
      do m = 1 , mmax , 2
        do i = m , n2 , istep
          j = i + mmax
          tempr = wr*data(j)-wi*data(j+1)
          tempi = wr*data(j+1)+wi*data(j)
          data(j)   = data(i) - tempr
          data(j+1) = data(i+1) - tempi
          data(i)   = data(i) + tempr
          data(i+1) = data(i+1) + tempi
        end do
        wtemp = wr
        wr = wr*wpr -    wi*wpi + wr
        wi = wi*wpr + wtemp*wpi + wi
      end do
      mmax = istep
    end do
  end subroutine fft_single
  
end module kfft_m


!program test
!  use kfft_m
!  implicit none
!
!  integer,  parameter :: DR = selected_real_kind(15)
!  real(DR), parameter :: PI = 3.1415926535897932_DR
!  real(DR), parameter :: TWOPI = PI*2
!
!  integer :: isign, n, ndata
!  integer, parameter :: nn = 256
!  integer, parameter :: ll = 10
!  real(DR) :: phi, dphi = TWOPI / nn
!  real(DR), dimension(ll,2*nn) :: data
!  complex(DR), dimension(ll,nn) :: c
!  integer :: mode
!
!  do n = 1 , nn
!    do ndata = 1 , ll
!      mode = ndata
!      phi = dphi*(n-1)
!      c(ndata,n) = cmplx(cos(mode*phi), sin(mode*phi), DR)
!    end do
!  end do
!
!
!  call iTest00
!
!  do n = 1 , nn
!    do ndata = 1 , ll
!!     data(ndata,2*n-1) = re(n-th data)
!!     data(ndata,2*n)   = im(n-th data)
!      mode = ndata
!      phi = dphi*(n-1)
!      data(ndata,2*n-1) = cos(mode*phi)
!      data(ndata,2*n)   = sin(mode*phi)
!    end do
!  end do
!
!  call iTest01
!
!contains
!
!  subroutine iTest00
!    isign = 1
!    call kfft(ll,nn,c,isign)
!
!    do n = 1 , nn
!      do ndata = 1 , ll
!        if ( abs(c(ndata,n)) > 1.e-10 ) then
!          print *, ndata, n, c(ndata,n)
!        end if
!      end do
!    end do
!  end subroutine iTest00
!
!  subroutine iTest01
!    isign = 1
!    call kfft(ll,nn,data,isign)
!
!    data(:,:) = data(:,:) / nn
!
!    do n = 1 , nn
!      do ndata = 1 , ll
!        if ( abs(data(ndata,n)) > 1.e-10 ) then
!          print *, ndata, n, data(ndata,n)
!        end if
!      end do
!    end do
!  end subroutine iTest01
!
!end program test

