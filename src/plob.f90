!***********************************************************************    
subroutine plob(y0,t0,tf, n, tdir, ftag, deriv, gr_cj,  y) 
!  this subroutine is to plot a segment of trajectory, with the intial state y0(6), 
!  and time span for integration is [0, tf], always save the energy or jacobi constant in the eight-th column to check the numerical error.
!  And save all the date in file ftag, with 2 blank lines between 2 different orbits

! Note: for backward in time integration, the stepsize h should be negative, and because h=1.d-3 is assigned inside the subroutine
!       we introduce  tdir to control the integration direction    :  1:unstable; -1:stable  
!       h = tdir * h   
   
!  better to do the write loop outside the subroutine--- 
!  in case to change the name of the data file, and write during the do loop

!  Input variables  
!    y0(*) :  inital state of trajectory 
!    t0    :  start time for the integration
!    tf    :  end time for the integration 
!    n     :  the dimension of state   
!    tdir  :  control the integration direction :  1: forward;  -1: backward
!    ftag  :  data file to save the result of integration, 1:t, 2-7:y, 8:cj 
!    deriv :  external subroutine to compute the vector field
!   gr_cj  :  external subroutine to compute the conservative quantity: energy or Jacobi Constant

! Output Variables
!    y(*)  :  finial state of trajectory 
  
!  ROUTINE USED: GR_RK78     GR_RTBP(A,B,N,F)
!  revised by Yu  -- 20160219
! --------------------------------------------------------------------------

  implicit none 
  integer, parameter :: dp=kind(1.d0)
    
  integer, intent(in) :: n, ftag, tdir    
  real(kind=dp), intent(in) :: y0(n),t0, tf
  real(kind=dp), intent(out) :: y(n) ! As the final state
  external :: deriv, gr_cj ! vector field  && energy or Jacobi constant

! Local Variables
  real(kind=dp) :: r(13,n),b(n),f(n), t, h, hmin, hmax, e,  & ! rk78 
		   cj  ! energy 

! copy the initial condition
  y = y0    

  t =  t0   ! initial time
  
  ! maybe small value is better, for a p.o. with small period, the orbit could be quite coarse 
  h =  tdir*1.d-3 

! specify the error control 
  hmin = 1.d-10
  hmax = 1.d0
  e    = 1.d-14

! if the evaluation condition is (t+h .lt. tf), only consider the positive step size, integate forward
! ( dabs(t+h-t0) .lt. dabs(tf-t0) ), both positive and negative cases are considered.

!  write(*,*) 'before plob, t0, tf, y0', t0, tf, y0  !ck

  do while( dabs(t+h-t0) .lt. dabs(tf-t0) ) 
  
    call gr_cj(y(1:6), cj)
    write(ftag,'(8e20.10)')  t, y(1:6), cj
    
!    write(*,'(8e20.10)')  t, y(1:6), cj!ck
    call gr_rk78(t,y, n, h,hmin,hmax,e,r,b,f, deriv )
    
  enddo
     
! if you want the strict periodic orbit, we need an additional step  
  if(dabs(t-tf) .gt. 1.d-9) then
  
    h =  tdir*tf - t ! for the stable orbit, the final time should be -tf
    
 ! to make sure the next step can be executed within allowable interval [hmin hmax]
!    hmin should be dabs(h)
    call gr_rk78(t,y, n, h, dabs(h), hmax, e,r,b,f, deriv)
    
  endif 
  
  call gr_cj(y(1:6), cj)      
  write(ftag,'(8e20.10)') t, y(1:6), cj

!  write(*,*) 'Finish plob, t0, t, tf, yf', t0, t, tf, y(1:6) !ck
  
  write(ftag,*)  ! better to save as a block than index
  
  end subroutine plob
