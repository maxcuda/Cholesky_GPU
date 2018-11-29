#if defined(__CUDA)
module cholesky
use cudafor
use cublas
use cusolverDn


logical, save         :: cusolver_init=.false.
type(cusolverDnHandle):: handle_cusolver
integer, device       :: devInfo_d
integer               :: istat
complex(8), parameter::  cone=cmplx( 1.0d+0, 0.0d+0 ), czero=cmplx( 0.0D+0, 0.0D+0 )
real(8)   , parameter::  one=1.0d+0, zero=0.0D+0




! Generic interfaces for compute:
!  Cholesky factorization: 
!     dpotrf(double precision)  and zpotrf( double complex matrix)
!  Inverse of a triangular matrix: 
!    dtrtri(double precision)  and ztrtri( double complex matrix)
! For the TRTRI routines, only the path used by Quantum-Espresso ( Lower, non-unit)
! is implemented on GPU.
! The CPU routines will be in the BLAS library, we are just specifying the interface
! to dispatch to the CPU or GPU based on the properties of the matrix A

interface dpotrf
    subroutine dpotrf( uplo, n, A, lda, info )
      use kinds, only : dp
      character:: uplo
      integer:: n, lda,info
      real(dp):: A(lda,*)
    end subroutine dpotrf

    module procedure dpotrf_gpu
end interface dpotrf

interface zpotrf
    subroutine zpotrf( uplo, n, A, lda, info )
      use kinds, only : dp
      character:: uplo
      integer :: n, lda,info
      complex(dp):: A(lda,*)
    end subroutine zpotrf

    module procedure zpotrf_gpu
end interface zpotrf

interface dtrtri
    subroutine dtrtri( uplo, diag, n, A, lda, info )
      use kinds, only : dp
      character:: uplo, diag
      integer:: n, lda, info
      real(dp):: A(lda,*)
    end subroutine dtrtri

    module procedure dtrtri_gpu
end interface dtrtri

interface ztrtri
    subroutine ztrtri( uplo, diag, n, A, lda, info )
      use kinds, only : dp
      character:: uplo, diag
      integer :: n, lda,info
      complex(dp):: A(lda,*)
    end subroutine ztrtri

    module procedure ztrtri_gpu
end interface ztrtri

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE dpotrf_gpu( uplo, n, A, lda, info )
      use kinds, only : dp
      implicit none
      character::          UPLO
      integer :: n,lda,info
      real(dp), intent(inout), device:: A(lda,*)
      !
      ! given a matrix A, returns the Cholesky decomposition of A
      ! for cholesky matrices
      !
      real(dp),  device:: work(1)
      integer :: lwork=1
      logical:: upper
      logical,external:: lsame

      if (.not. cusolver_init) then
        istat = cusolverDnCreate(handle_cusolver)
        if (istat /= CUSOLVER_STATUS_SUCCESS) write(*,*) 'cusolver handle creation failed'
      end if

      info = -1
      upper = lsame( uplo, 'U' )
      if(upper) then
       istat= cusolverDnDpotrf(handle_cusolver, CUBLAS_FILL_MODE_UPPER, n, A, lda, work, lwork, devInfo_d)
      else
       istat= cusolverDnDpotrf(handle_cusolver, CUBLAS_FILL_MODE_LOWER, n, A, lda, work, lwork, devInfo_d)
      endif
      info=devInfo_d


    END SUBROUTINE dpotrf_gpu

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     SUBROUTINE zpotrf_gpu( uplo, n, A, lda, info )
      use kinds, only : dp
      implicit none
      character::          UPLO
      integer :: n,lda,info
      complex(dp), intent(inout), device:: A(lda,*)
      !
      ! given a matrix A, returns the Cholesky decomposition of A
      ! for cholesky matrices
      !
      complex(dp),  device:: work(1)
      integer :: lwork=1
      logical:: upper
      logical,external:: lsame

      if (.not. cusolver_init) then
        istat = cusolverDnCreate(handle_cusolver)
        if (istat /= CUSOLVER_STATUS_SUCCESS) write(*,*) 'cusolver handle creation failed' 
      end if

      info = -1
      upper = lsame( uplo, 'U' )
      if(upper) then
       istat= cusolverDnZpotrf(handle_cusolver, CUBLAS_FILL_MODE_UPPER, n, A, lda, work, lwork, devInfo_d)
      else
       istat= cusolverDnZpotrf(handle_cusolver, CUBLAS_FILL_MODE_LOWER, n, A, lda, work, lwork, devInfo_d)
      endif
      info=devInfo_d


    END SUBROUTINE zpotrf_gpu





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The two following routines are the GPU implementation of the ZTRTRI LAPACK routine. 
! This routine computes the inverse of a complex lower non-unit triangular matrix.
! This is the unblocked code that operates on a block 32x32
   attributes(global) subroutine ztrti2_gpu(n,a,lda)
      use cudafor
      implicit none
      integer, value    :: n,lda
      complex(8),device :: a(lda,*)
      complex(8),shared :: a_s(32,32)
      complex(8)        :: cv
      integer           :: tx,ty,tl,i,j,ii


     if (n == 0 ) return

      tx=threadIdx%x
      ty=threadIdx%y
      ! Linear id of the thread (tx,ty)
      tl=tx+ blockDim%x*(ty-1)

      ! Load a_d in shared memory
      if (tx <= n .and. ty <= n) then
         a_s(tx,ty)=a(tx,ty)
      endif

      call syncthreads()

     ! Compute all the diagonal elements

     if (tl <=n ) a_s(tl,tl)=cone/a_s(tl,tl)

      call syncthreads()


     ! For each column working backward
      do i=n-1,1,-1

         if ( tl >i .and. tl <= n) then
            cv=czero
             do j= i+1,tl
               cv = cv + a_s(j , i ) * a_s(tl ,j)
             end do
         end if
!        call syncthreads()
         if ( tl >i .and. tl <= n) then
           a_s(tl,i)=-cv*a_s(i,i)
         end if

!        call syncthreads()

       end do

        call syncthreads()

! Write back to global memory
      if (tx <= n .and. ty <= n) then
         a(tx,ty)=a_s(tx,ty)
      endif

     end subroutine ztrti2_gpu

! This routine computes the inverse of a complex lower non-unit triangular matrix.
! This is the blocked code that operates on a generic size matrix
 subroutine ztrtri_gpu(uplo, diag, n, a, lda, info )
      use cublas
      implicit none

      character :: uplo, diag
      integer   :: lda, n, info,sizeb
      complex(8),device:: a( lda, * )

      logical:: nounit, upper, supported
      logical, external:: lsame
      integer:: nb,nn,j,jb
      type(dim3):: threads


      if (n == 0 ) return

      upper = lsame( uplo, 'U' )
      nounit = lsame( diag, 'N' )
      IF( .not.upper .and. .not.lsame( uplo, 'L' ) ) THEN
         info = -1
      else if( .not.nounit .and. .not.lsame( diag, 'u' ) ) then
         info = -2
      else if( n.lt.0 ) theN
         info = -3
      else if( lda.lt.max( 1, n ) ) theN
         info = -5
      end iF

      ! Check if this is a supported configuratio, lower non unit triangular
      supported=(.not.upper) .and. nounit
      if ( .not.supported) then
         print *," Only L and  no-unit  implemented on GPU", upper,nounit
         info=-6
         return
      endif 
      
      ! Block size for zttri2 on GPU is 32
       nb=32
       threads=dim3(32,32,1)

      if( nb .ge. n ) then
       ! If the problem size is smaller than block size, call directly unblocked code
           call ztrti2_gpu<<<1,threads>>>( n, a, lda )
      else
       ! Use blocked code
       !
       ! Compute inverse of lower triangular matrix
       !
            nn = ( ( n-1 ) / nb )*nb + 1
            do j = nn, 1, -nb
               jb = min( nb, n-j+1 )
               if( j+jb.le.n ) then
                !
                ! Compute rows j+jb:n of current block column
                !
                  call ztrmm( 'Left', 'Lower', 'No transpose', diag, &
                              n-j-jb+1, jb, cone, a( j+jb, j+jb ), lda, &
                              a( j+jb, j ), lda )
                  call ztrsm( 'Right', 'Lower', 'No transpose', diag, &
                              n-j-jb+1, jb, -cone, a( j, j ), lda, &
                              a( j+jb, j ), lda )
               end if
               !
               ! Compute inverse of current diagonal block
               !
               call ztrti2_gpu<<<1,threads>>>( jb, a( j, j ), lda )
            end do
      endif

     end subroutine ztrtri_gpu

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The two following routines are the GPU implementation of the ZTRTRI LAPACK routine. 
! This routine computes the inverse of a double precision lower non-unit triangular matrix.
! This is the unblocked code that operates on a block 32x32
   attributes(global) subroutine dtrti2_gpu(n,a,lda)
      use cudafor
      implicit none
      integer, value    :: n,lda
      real(8),device :: a(lda,*)
      real(8),shared :: a_s(32,32)
      real(8)        :: cv
      integer           :: tx,ty,tl,i,j,ii


     if (n == 0 ) return

      tx=threadIdx%x
      ty=threadIdx%y
      ! Linear id of the thread (tx,ty)
      tl=tx+ blockDim%x*(ty-1)

      ! Load a_d in shared memory
      if (tx <= n .and. ty <= n) then
         a_s(tx,ty)=a(tx,ty)
      endif

      call syncthreads()

     ! Compute all the diagonal elements

     if (tl <=n ) a_s(tl,tl)=one/a_s(tl,tl)

      call syncthreads()


     ! For each column working backward
      do i=n-1,1,-1

         if ( tl >i .and. tl <= n) then
            cv=zero
             do j= i+1,tl
               cv = cv + a_s(j , i ) * a_s(tl ,j)
             end do
         end if
!        call syncthreads()
         if ( tl >i .and. tl <= n) then
           a_s(tl,i)=-cv*a_s(i,i)
         end if

!        call syncthreads()

       end do

        call syncthreads()

! Write back to global memory
      if (tx <= n .and. ty <= n) then
         a(tx,ty)=a_s(tx,ty)
      endif

     end subroutine dtrti2_gpu

! This routine computes the inverse of a complex lower non-unit triangular matrix.
! This is the blocked code that operates on a generic size matrix
 subroutine dtrtri_gpu(uplo, diag, n, a, lda, info )
      use cublas
      implicit none

      character :: uplo, diag
      integer   :: lda, n, info,sizeb
      real(8),device:: a( lda, * )

      logical:: nounit, upper
      logical, external:: lsame
      integer:: nb,nn,j,jb
      type(dim3):: threads


      if (n == 0 ) return

      upper = lsame( uplo, 'U' )
      nounit = lsame( diag, 'N' )
      IF( .not.upper .and. .not.lsame( uplo, 'L' ) ) THEN
         info = -1
      else if( .not.nounit .and. .not.lsame( diag, 'u' ) ) then
         info = -2
      else if( n.lt.0 ) theN
         info = -3
      else if( lda.lt.max( 1, n ) ) theN
         info = -5
      end iF

      if (.not. upper .or. .not.nounit ) then
         print *," Only L no unit  implemented on GPU"
      endif 
      
      ! Block size for dttri2 on GPU is 32
       nb=32
       threads=dim3(32,32,1)

      if( nb .ge. n ) then
       ! If the problem size is smaller than block size, call directly unblocked code
           call dtrti2_gpu<<<1,threads>>>( n, a, lda )
      else
       ! Use blocked code
       !
       ! Compute inverse of lower triangular matrix
       !
            nn = ( ( n-1 ) / nb )*nb + 1
            do j = nn, 1, -nb
               jb = min( nb, n-j+1 )
               if( j+jb.le.n ) then
                !
                ! Compute rows j+jb:n of current block column
                !
                  call dtrmm( 'Left', 'Lower', 'No transpose', diag, &
                              n-j-jb+1, jb, one, a( j+jb, j+jb ), lda, &
                              a( j+jb, j ), lda )
                  call dtrsm( 'Right', 'Lower', 'No transpose', diag, &
                              n-j-jb+1, jb, -one, a( j, j ), lda, &
                              a( j+jb, j ), lda )
               end if
               !
               ! Compute inverse of current diagonal block
               !
               call dtrti2_gpu<<<1,threads>>>( jb, a( j, j ), lda )
            end do
      endif

     end subroutine dtrtri_gpu

end module

#endif
