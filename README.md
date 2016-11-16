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
