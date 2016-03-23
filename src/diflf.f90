subroutine dflrtz()
! The general form of the differential of the lorentz force, as a function
! of parameter beta 

!    Input
!  a	n-by-n real general matrix
!  n	dimension of matrix a
!  isv	1:compute also the right eigenvector, 0: only eigenvalue

!	Output
!  wr,wi the real and imaginary part the eigenvalues
!  vr    right eigenvector

implicit none
integer, parameter :: dp = kind(1.d0), lwmax = 1000

integer, intent(in) :: n, isv
real(kind=dp), intent(in) :: a(n,n)
real(kind=dp), intent(out) :: wr(n), wi(n), vr(n,n)

character:: jobvl, jobvr

integer  ::  lda, ldvl, ldvr, info, lwork
!
