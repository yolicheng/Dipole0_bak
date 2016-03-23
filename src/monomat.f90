subroutine monomat(yi,tp, mmat, deriv, gr_cj)
! compute the Monodramy matrix of the p.o. of period tp
! 	Input Varialbles
!   yi		initial state(x0,y0,z0,vx0,vy0,vz0)
!   tp 		period of the P.O.      

! 	Output Variables
!   mmat	Monodramy matrix, which is also the state transition matrix after one period

!  Routine used: dfcr, fctn, champ, adams, gr_cjlf

 
! 20150408 --- the vertical lyapunov checked  
!     imax is set to 1 to get the first intersection, pay attention to
!     the computation of g for new PO

! Finally revised : 20150525
!***********************************************************************
 
implicit none
  
!  both the state and the variational equations are desired, so n=42
integer, parameter:: n = 42, dp = kind(1.d0) 
  
real(kind=dp), intent(in)  :: yi(6), tp 
real(kind=dp), intent(out) :: mmat(6,6)
external ::  deriv, gr_cj ! here should be gr_lf (the vector field for lorentz force problem)

integer :: i, fout 
real(kind=dp) :: y(n), yf(n), t, h, hmin, hmax, e, r(13,42),b(42),f(42), &
   		 cj, cj2, dcj


h    = 1.d-2

! this is the same... should we use the value from the module pomod? no...
hmin = 1.d-10
hmax = 1.d0
e    = 1.d-13

t    = 0.d0 ! start time is 0
! initialize vector field and variational matrix(identity)
y(1:6)    = yi
y(7:42)   = 0.d0
y(7:42:7) = 1.d0


! check if the energy is conversative 
call gr_cj(yi, cj)
 
do while (dabs(t+h) .le. dabs(tp))

!SUBROUTINE GR_RK78 (X,Y,N,H,HMI,HMAX,E1,R,B,F,DERIV)
  call gr_rk78(t,y, n,h,hmin,hmax,e,r,b,f, deriv)
enddo  


if (dabs(t-tp) > 1.d-9 ) h = tp - t 

! h for sure will be smaller than the previous h, so 1 step is enough
call gr_rk78(t,y, n,h,hmin,hmax,e,r,b,f, deriv )

! check if the energy integral is conservative
call gr_cj(y(1:6),cj2) 
dcj = cj2 - cj 

!print*,'check  if the energy integral is conservative -- monomat'  
!print*, 'dcj, cj, cj2', dcj, cj, cj2
!read*

mmat = reshape(y(7:42), (/6,6/)) 
  
! check the final state, should be the same with the initial state 
!print*,'check  the final state -- monomat'  
!print*, 'yf, y0', y(1:6), yi(1:6)
 
!do i = 1, 6
!  write(*,'(6d20.10)') mmat(i,:)
!enddo  
 
  

return
end subroutine monomat
	 
	 
