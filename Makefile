all: test_cholesky test_zttri

test_cholesky:  kind.f90 cholesky.f90 main_cholesky.f90 
	pgf90 -Mcuda=cc70 -Mpreprocess -D__CUDA kind.f90 cholesky.f90 main_cholesky.f90 -o test_cholesky  -lblas -lcusolver

test_zttri:  kind.f90 nvtx.f90 cholesky.f90 main_trtri.f90 
	pgf90 -Mcuda=cc70 -Mpreprocess -D__CUDA kind.f90 nvtx.f90   cholesky.f90 main_trtri.f90 -o test_trtri  -lblas -lcusolver -L/usr/local/cuda/lib -lnvToolsExt

clean:
	rm *.o test_cholesky test_trtri  *.mod
