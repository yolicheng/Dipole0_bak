subroutine plob_fxd(y0,t0, h0, nsm, smpl, xmax, ftag, yall,  deriv) 
!********************************************************************
!  this subroutine is to plot a segment of trajectory, the intial state y0(6), 
! time span for integration is [0, tf] 

!  And save all the date in file ftag, with 2 blank lines between 2 different orbits

! 20160128
!  add the escape constraint, stop integration once we are outside this domain. 
   
!  better to do the write loop outside the subroutine--- 
!  in case to change the name of the data file, and write during the do loop

!  Insmut variables  
!    y0(6) :  inital state of trajectory 
!    t0    :  start time for the integration
!    h     :  the fixed step size
!    nsm   :  number of sampling points  
!    smpl  :  the sampling increment,  to save every smpl points
!   xmax   :  the maximum modulu of x component, to detect escape
! 	  	we still have the flaw that it's hard to detect which component is the one to moniter the escape 
!               try the maximum?

!    ftag  :  data file to save the result of integration 
!   deriv  :  the vector field to be integrated 

! Output Variables
!yall(nsm)  :  the state of all the sampling points
  
!  ROUTINE USED: GR_RK78     DERIV(A,B,N,F)
!********************************************************* ***********

implicit none	
integer, parameter:: dp = kind(1.d0) 

integer, intent(in):: nsm, smpl, ftag   
real(kind=dp), intent(in):: y0(6),t0, h0, xmax
real(kind=dp), intent(out):: yall(nsm, 6) ! As the final state
external deriv

! 	Local variables
integer :: ni, ismp, nsmp 
real(kind=dp) :: r(13,6),b(6),f(6), h, hmin, hmax, e1, t, y(6)  !rk78  


! Initial value for Earth-Moon system
h = h0

! use the same value to make it as fixed step size
hmin = h  
hmax = h  
e1   = 1.d-13
 
y = y0      
t =  t0 !0.d0

ni = 1 ! number of integration step
! Integrate the periodic orbits by gr_rk78, write the x, y in txt file with ftag

nsmp = 0 
do while( ni .lt. nsm*smpl) 
  
  ismp = mod(ni-1, smpl)
  
! save the integration result every smpl steps, in this way, we could use a small h, just save the necessary smpl points to do fft
  if ( ismp == 0) then  
    nsmp = nsmp + 1  ! record the number of sampling points
    yall(nsmp,:) = y
    write(ftag,'(7e20.10)') t+t0, y
!  write(*,'(1I5, 8e20.10)') ni, t, h, y
  endif
  
  call gr_rk78(t,y,6,h,hmin,hmax,e1,r,b,f, deriv)
  
  if( maxval(dabs(y(1:3))) > xmax  ) then ! once escape, return
!  if(dabs(y(1)) > xmax .or. dabs(y(3)) > xmax .or. dabs(y(2)) > xmax) then ! once we escape, return
    write(ftag,*) 
    print*, 'Escape!'
    read(*,*)
!   close(ftag) ! doesn't need to close the file 
    return
  endif  
  ni = ni+1
enddo

! save the last point, only needed when smpl == 1
if (smpl == 1) then 
  yall(ni,:) = y
  write(ftag,'(7e20.10)') t+t0, y    
endif 
 
!write(*,*) 'after plob_fxd, t0, t, tf, yf', t0, t, tf, y !ck
!read(*,*)  aaa !ck

write(ftag,*)  ! add one blank line to seperate orbit
!write(ftag,*) 
	 
end subroutine plob_fxd


