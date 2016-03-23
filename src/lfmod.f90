module lfmod
!  Put all the commom variables into this module for lorentz force problem, to avoid passing a lot of shared arguments down to all
!  subroutines that use them

!  for the safe side, declare the variables that are only used in this module private attribute, and use subroutine init_lfmod to do the assignment

! Todo - 20151223 -finished for planar symmetry, and x-axis symmetry --20160302
! define the type of symmetry and the corresponding index of the componts to 
!   be as control and target variables 
! at this moment, we only deal with symmetric wrt xz plane, and initial point is 
!   on the y=0 plane
!
! common variables for all 3 cases
! beta 	 = n/w_c, the most important parameter,  where
! 	 n   is angular velocity of the mean orbital rate of the chief’s circular reference orbit
!        w_c is the magnitude of the dipole’s rotational angular rate
! eq 	 the chosen equilibrium 

implicit none
integer, parameter :: dp = kind(1.d0), n = 6 ! 6 dimensional state
integer, parameter, private :: debug = 0 

! Public -- the most common used ones, do not need to declare in the main routine 
real(kind=dp) ::  eq(n) 
integer :: sgn

! Private -- accessible only by the module subroutines and functions
integer, parameter, private  :: neq1 = 3, neq2 = 3,  neq3 = 3
integer, private, save       :: cs, sgn1(neq1), sgn2(neq2), sgn3(neq3) 
real(kind=dp), private, save :: beta, x0, eq1(neq1,n), eq2(neq2,n), eq3(neq3,n) 

contains

! system based subroutines -- in seperate files to avoid to long of a single file
!  include 'x2g.f90'
!  include 'invx2g.f90'
!  include 'dflrtz_g3.f90'

  include 'dflrtz.f90' ! put the first three into dlfrtz 
  include 'gr_lf.f90'
  include 'lfunit.f90' ! the function to compute the unit of distance, time and velocity
  
! Subroutines   
!1.	subroutine init_lfmod   -- private constants initialization
!2.	subroutine init_lf(beta0, cs0, ieq)
!3.	subroutine gr_cjlf(y, cj)

!subroutine x2g(x0, xg)
!subroutine invx2g(x0, xg)
!subroutine dflrtz_g3(x0, dlf)
!subroutine dflrtz(x0, dlf)
!subroutine gr_lf(t, y, neq, f)

  subroutine init_lfmod  
! finish only case 1, with eq1, sgn1 -20160302

 ! cs : case ---  1: N=[ 1 0 0]; 2: N=[ 0 1 0]; 3: N=[ 0 0 1]
  eq1 = 0.d0 
    
  ! the first equilibrium point, q/m > 0
  sgn1(1) = 1
  ! x=0,y=0,z= \pm 1
  eq1(1,3) = 1.d0
    
  ! the second equilibrium point, q/m < 0
  ! x = \pm (2 sqrt(9))^(1/3), y² = 2x²,z = 0
  sgn1(2) = -1
  x0 = ( 2.d0 / (9.d0*dsqrt(3.d0)) ) ** (1.d0/3.d0)
 
  eq1(2,1) = x0
  eq1(2,2) = dsqrt(2.d0) * x0

  ! the third equilibrium point, q/m < 0, this is what we choose to study, at beta=10, has 3 famililes of center manifold
  ! we are taking z = x in this case. 
  sgn1(3) = -1
  ! x = \pm (1/4/sqrt(2)^(1/3), y = 0,z² = x²
  x0 = ( 1.d0/ (4.d0*dsqrt(2.d0)) ) ** (1.d0/3.d0)
  eq1(3,1) = x0
  eq1(3,3) = x0 !z0 

  end subroutine init_lfmod


  subroutine init_lf(beta0, cs0, ieq)
  !to initialize the specified case to study, parameters to be initialized: beta,cs,sgn,ieq
  
    real(kind=dp),intent(in) :: beta0  
    integer, intent(in) :: cs0, ieq 
    cs = cs0
    beta = beta0

! the state of the equilibrium point, positon+velocity 
    if (cs0 == 1) then 
      eq = eq1(ieq,:)
      sgn = sgn1(ieq)
!      print*, 'check init_lf: ieq, sgn, sgn1', ieq, sgn, sgn1 !ckd
!      read(*,*) aaa
     
    elseif (cs0 == 2) then 
      eq = eq2(ieq,:)
      sgn = sgn2(ieq)
     
    elseif (cs0 == 3) then 
      eq = eq3(ieq,:)
      sgn = sgn3(ieq)
   
    endif  
  end subroutine init_lf


!***********************************************************************
subroutine gr_cjlf(y, cj)
! use include to split this module into different files, with gr_rk78, deriv, gr_lf in external fiels
!***********************************************************************
!  for the case of lorentz force, compute the energy integral
!  c = 3x²-z² - sgn* 2(y²+z²)/r³ - (dotx² + doty² + dotz²), sgn:1 if q/m>0; -1 if q/m<0

!	input
!  y(6) 	the initial state
!  sgn   	1 if q/m>0; -1 if q/m<0
!------------------------------------------------------------------------------

implicit none 

real(kind=dp), intent(in)  :: y(6)
real(kind=dp), intent(out) :: cj 
 
! Local Variables 
real(kind=dp) :: y1, y2, y3, r, r3, dv2   

y1=y(1)**2
y2=y(2)**2
y3=y(3)**2
        
r = dsqrt(y1 + y2 + y3) 
r3 = r**3
dv2 = y(4)**2 + y(5)**2 + y(6)**2 

 cj = 3*y1 - y3 - sgn * 2 * (y2+y3)/r3 - dv2 

return 
end subroutine gr_cjlf
 

end module lfmod

