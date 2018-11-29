program inverse
use cholesky
use cudafor
use cublas
integer, parameter :: n=1024
complex(8):: A_h(n,n),B_h(n,n)
complex(8),device :: B(n,n)

real*8:: x,y,threshold=1.d-14
integer:: info
type(dim3):: threads
type(cudaEvent):: startEvent,stopEvent
type(cusolverDnHandle):: h
real:: time
logical:: print=.false.

istat=cudaEventCreate(startEvent)
istat=cudaEventCreate(stopEvent)

! Generate hermitian positive definite matrix

do j=1,n
 do i=j,n
  call random_number(x)
  call random_number(y)
  A_h(i,j)=cmplx(x,y)
  A_h(j,i)=conjg(A_h(i,j))
  if (i==j) A_h(i,j)=real(A_h(i,j))
 end do
end do

call zgemm('N', 'C', N, N, N, cmplx(1.0, 0.0, 8), A_h, N, A_h, N, cmplx(0.0, 0.0, 8), B_h, N)
A_h=B_h

!A_h is an hermitian positive definite matrix

! Copy the matrix to the device
B=A_h

print *,"calling CPU zpotrf and ztrtri"
! Call Cholesky
call zpotrf('L', n, A_h, n, INFO )
! Compute the inverse 
CALL ztrtri( 'L', 'N', n, A_h, n, INFO )

if (info .lt.  0) print *,"Error in CPU path !!!!!!!!!!!!!",info

print *,"calling GPU zpotrf and ztrtri"
! Call Cholesky 
call zpotrf('L', n, B, n, INFO )
! Compute the inverse 
CALL ZTRTRI( 'L', 'N', n, B, n, INFO )
if (info .lt.  0) print *,"Error in GPU path !!!!!!!!!!!!!",info


! Copy result back to host to compare it to CPU solution
B_h=B

print *,"Difference "
   print *,(B_h(1,1)-A_h(1,1))
   print *,(B_h(n/2,n/2)-A_h(n/2,n/2))
   print *,(B_h(n,n)-A_h(n,n))
end program inverse
