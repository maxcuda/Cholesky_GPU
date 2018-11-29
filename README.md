# Cholesky_GPU

Generic Fortran interfaces for :

   Cholesky factorization: 
   dpotrf(double precision)  and zpotrf( double complex matrix)

Inverse of a triangular matrix: 
   dtrtri(double precision)  and ztrtri( double complex matrix)

For the TRTRI routines, only the path used by Quantum-Espresso ( Lower, non-unit)
 is implemented on GPU.
