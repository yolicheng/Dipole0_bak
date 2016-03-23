subroutine poinc_n(y0, imax,  xmax,  fpc, deriv, gr_cj,  tf, yf, ncrs)  
! This routine is to do a Poincare map for imax revolutions,  -- only works for lf problem, 
!  the maximum time and domain for one revolution is tmax and xmax, respectively.
!  the poincare section is defined by ind-th component of the state vector 
!    y(ind) = z0, with velocity in that direction y(ind+3) of positive value 


! 20160202 -- Add domain constraint, if x,y and z component is beyong some value, xmax, we say the orbit escape, stop the integratin  

! 	Input
!  y0(6)	the initial state 
!  imax		maximum number of crossing through the poincare section  y(ind) = sec 
!  tmax		the maximum integration time (allowed to be updated)
!  xmax		the maximum domain for the orbit, beyond that the orbit is treated as escape 
!  ind 		the index of component in state vector for poincare map y(ind) = z0
!  fpc		file tag to save tpc,ypc at every crossing


! 	Output 
!  yf(6)	the final state at ncrs-th crossing
!  tf		the final time
!  ncrs		real number of crossing through the poincare section (may in)

!use lfmod

!   ----- pomod ----- 
use pomod, only: ind, sec,  hmin, hmax, e, tmax &  ! parameters 
                 sect, poinc  !subroutines
implicit none
integer, parameter :: dp=kind(1.d0)


! input and output declarations
integer, intent(in)  :: fpc, imax
integer, intent(out) :: ncrs
real(kind=dp), intent(in)  :: y0(6), xmax !tmax- possibly be updated, can not specified as input here, (not able to change the value then)
real(kind=dp), intent(out) :: yf(6), tf
external :: deriv, gr_cj! the differential of the vector field and the variational matrix



! local variables
integer :: i,ispc  

real(kind=dp) ::  yi(42), ypci(42), t,  tpci, tpc, &
                  cj, hminim 
                  

! Initial the input data for poinc, both vector field and variational matrix(identity matrix)
yi         = 0.d0
yi(1:6)    = y0  
yi(7:42:7) = 1.d0  

!print*,'poinc_n, x0', y0
!read(*,*) aaa
! use poinc to obtain one crossing, and then use this point as the new initial point to get the nect poincare map
! save the time and state at each crossing to file fpc(eqpoinc.dat)

tpci = 0.d0
ncrs = 0

do  
 
! subroutine poinc(yi,imax, tf,yf, hminim, ispc, deriv, gr_cj)  
  call poinc(yi,  1,  tpc, ypci, hminim, ispc, deriv, gr_cj)  
! !  call poinc_z0(yi,1, tmax, ind, z0, tpc, ypci, ispc, deriv)  
  
  if(ispc == 0 ) exit ! tpc, ypci not updated if there is no intersecion

! 20160202 -- Add domain constraint, if x,y and z component is beyong some value, xmax, we say the orbit escape, stop the integratin  
  if ( maxval( dabs( ypci(1:3)) ) > xmax )  then 
    ispc = 0
    exit
  endif
  
  if (ypci( 3+ind) .le. 0.d0) cycle ! only keep the plus velocity when crossing the poincare section
  
  ncrs = ncrs + 1
  tpci = tpci + tpc
  yi = ypci
  
  write(fpc, '(7f16.8, 1I5)') tpci, ypci(1:6),  ncrs
!  write(*, '(7f16.8)') tpci, ypci(1:6), ncrs
!  read(*,*) aaa
  
  if(5*tpc > tmax) tmax = 5*tpc  ! update maximum integration time, which is used as the stop condition
  
  if(ncrs .ge. imax) exit ! use the maximum intersecion as the stop condition
enddo 

! return the state at the last crossing
tf = tpci
yf = ypci(1:6)

return
end subroutine poinc_n



















  
