# sph-dilts
Spherical harmonics expansion of data on a spherical surface

## fft.f03 ##
A FFT module

## sph.f03 ##
 A library for spherical harmonics expansion.
 
>      input : a function on a spherical surface
>                 f(lat,lon)
>     output : sphereical harmonics coefficients of f
>                 ff(lm2lm(l,m)) 
>                 where 0<=l<=NN-1, -l<=m<=l
>                 Note: l is upto NN-1.

### Method ###

This is basically Dilts' method for the spherical 
harmonics expansion:
  Reference: "Computation of Spherical Harmonic Expansion
              Coefficients via FFT's", Gary A. Dilts,
              J. Comput. Phys., Vol.57, p.439, 1985

A slight difference is that Fourier coeeficients of
associated Legendre function, B_j^{lm} in the above
paper, are numerically calculated by FFT in this program.
(In the the above paper, they are analytically derived
by recurrence formula.) This FFT approach is simple and
accurate since associated Legendre functions are Fourier 
series. I did not try to save memory for zero values of
B_j^{lm}. 

For applying FFT in latitudinal direction to the
input function f and Legendre functions Plm, their
domain [0,pi] is expanded to [0,2*pi].

### Normalization ###


>     ff_{lm} = \int_S dS f(\theta,\phi) Y_{lm}^\dagger
>   where
>     \dagger = complex conjugate
>     \int_S dS = \int_0^{2\pi} dphi \int_0^\pi d\cos\theta 
>     \int_S dS Y_{lm} Y_{ab}^\dagger = \delta_{la} \delta{mb}
>     Y_{lm}(\theta,\phi) = C_{lm} P_l^m e^{im\phi}
>     C_{lm} = \sqrt{(2l+1)/(4\pi) (l-m)!/(l+m)!
>   for negative m, 
>     P_l^{-m} = P_l^{m}
>     C_l^{-m} = C_l^m

  
### Usage ###
   
>   (0) set input function f to be expanded.
>   (1) call sph%initialize(nn,lat,lon)
>   (2) ff = sph%expand(nn,lat,lon,f)
 

 
### Note on the size of N ###
This module uses "fft_m", which is an FFT library.
Since fft_m is a simple FFT program, the data size 
is limited to be power of two. Therefore, a key
number NN in this program is always power of two.
But this is not essential for the basic algorithm
of the Dilts' method. You can replace fft_m with
other FFT library so that NN can be arbitrary number.

### Hisotory ###
This program is converted from my old Fortran code 
"dilts.f" that was developed in July 1993.
