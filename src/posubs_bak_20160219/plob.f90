subroutine plob(y0,t0,tf, tdir, ftag, deriv, gr_cj, y) 
!  this subroutine is to plot a segment of trajectory, with the intial state y0(6), 
!  and time span for integration is [0, tf], always save the energy or jacobi constant in the eight-th column to check the numerical error.
!  And save all the date in file ftag, with 2 blank lines between 2 different orbits

! Note: for backward in time integration, the stepsize h should be negative, and because h=1.d-3 is assigned inside the subroutine
!       we introduce  tdir to control the integration direction    :  1:unstable; -1:stable  
!       h = tdir * h   
   
!  better to do the write loop outside the subroutine--- 
!  in case to change the name of the data file, and write during the do loop

!  Input variables  
!    y0(6) :  inital state of trajectory 
!    t0    :  start time for the integration
!    tf    :  end time for the integration 
!    tdir  :  control the integration direction :  1: forward;  -1: backward  
!    ftag  :  data file to save the result of integration, 1:t, 2-7:y, 8:cj 
!    deriv :  external subroutine to compute the vector field
!   gr_cj  :  external subroutine to compute the conservative quantity: energy or Jacobi Constant

! Output Variables
!    y(6)  :  finial state of trajectory 
  
!  ROUTINE USED: GR_RK78     GR_RTBP(A,B,N,F)
!  revised by Yu  -- 20160219

implicit none 
integer, parameter:: dp = kind(1.d0) 

integer, intent(in):: ftag, tdir   
real(kind=dp), intent(in):: y0(6),t0, tf
real(kind=dp), intent(out):: y(6) ! As the final state
external :: deriv, gr_cj ! vector field  && energy or Jacobi constant

! Local Variables
real(kind=dp) :: r(13,6),b(6),f(6),  hmin, hmax, e, t, h, & ! rk78
		 cj 

! Initial value for Earth-Moon system
 
hmin = 1.d-10
hmax = 1.d-1
e    = 1.d-13
 
y = y0      
 
! Integrate the periodic orbits by gr_rk78, write the x, y in txt file with ftag
t =  t0 !0.d0
h =  tdir*1.d-3

! if the evaluation condition is (t+h .lt. tf), only consider the positive step size, integate forward
! ( dabs(t+h-t0) .lt. dabs(tf-t0) ), both positive and negative cases are considered.
!write(*,*) 'before plob, t0, tf, y0', t0, tf, y0  !ck

do while( dabs(t+h-t0) .lt. dabs(tf-t0) ) 
  call gr_cj(y(1:6), cj)
  write(ftag,'(7e20.10)') t, y, cj
  call gr_rk78(t,y,6,h,hmin,hmax,e,r,b,f,  deriv)
enddo
     
!if you want the strict periodic orbit, a control of tf must be made
if(dabs(t-tf) .gt. 1.d-9) then
  h =  tf - t
  call gr_rk78(t,y,6,h,hmin,hmax,e,r,b,f, deriv)
endif 
        
write(ftag,'(7e20.10)') t, y

!write(*,*) 'Finish plob, t0, t, tf, yf', t0, t, tf, y !ck
!read(*,*)   !ck

write(ftag,*)  ! better to save as a block than index
  
end subroutine plob


