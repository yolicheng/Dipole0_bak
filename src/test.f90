program test
! this is to compute the eigenvalues of the first case, to see if it matches the result in Peng's paper

!  consider the 1st equilibrium point  ! 1:  q/m > 0,  x=0,y=0,z= \pm 1
!  take beta = 1, we have 1 pair of pure imaginary eigenvalues,  and 2 pairs of real eigenvalues()

use lfmod, only : dp, cs, sgn, eq, n, init_lfmod, init_lf

implicit none

!np = 2^21, lt = 2.d3 *2*pi is the best match, we have x-z clearly extracted.
! fs = nsm / lt /2/pi  !  sampling frequency: per cycle


! number of sampling points on the orbit: 2**32, 2**24, 2**16 points should be good
integer, parameter ::   nsm = 2**18,  &  !number of sampling points 
                        smpl = 2**5, np = nsm*smpl,  & ! steps of integration  
			neq = 42, & ! number of dimensions for monodrony matrix computation 
			nfq = 10, & ! keep the first 4 frequencies with the biggest amplitude
			nob = 10 !number of initial points to explore
			 
real(kind=dp), parameter :: pi = 4*datan(1.d0), lt = 5.d2 * 2 * pi, & !! time length of sampling 
 			    xmax = 1.d0 !  when x > xmax, escape,  stop the integration


! specify the case, equilibrium point and parameter beta that we are going to study
integer ::  cs0, ieq 
real(kind=dp) ::  beta0,  dlf(n,n), beta

! local variables
integer :: i, j, ix, iy, icj, iscj,  imf, imax, ncrs, feqmf, feqpc, feqcj,  feq1, &
	   isign, f_fft1, f_fft2, f_fft3, f_fft4, f_fft, fmx, obi, obf, cji, cjf!  --fft

real(kind=dp) ::  wr(n),wi(n), vr(n,n), y0(neq), &   ! stability of equilibrium points
                  x0(6), vr_usd(6,3), dx0(6), epsl, epsl_dx(3,nob), dx, cj_tar, & !manifold of eq 
                  xf(6), tf, xpc(6), tmax, tfpl, & !  poinc_n of eq
                  cj_eq, cj, cj2, aaa, & ! debug auxillary variables
                  xl(nsm,6), x(nsm),y(nsm), z(nsm), a(nsm), h, t, & !, h fft  
                  f(nsm/2), amp1, amp2,  amp3, amp4, amp_min, iamp, ramp, &! write to file  
                  tst, tend   ! evaluate time  
      
DOUBLE COMPLEX fft1(nsm),fft2(nsm), fft3(nsm),fft4(nsm) 
                 
external :: gr_lf ! the differential of the vector field and the variational matrix

! use this to initialize all the equilibrium points and the respective sign of (q/m), 
! which are the private variable in module lfmod
call init_lfmod


! the parameter beta, ratio between the angular velocity of mean motion of the chief around the earth and 
!                     the angular velocity of the rotaion of the deputy

! case 1: N=[1 0 0]; 2:  N=[0 1 0]; 3:  N=[0 0 1]; 
 cs0 = 1 ! use 2 to test the swap rule is coded correctly--ckd!
ieq = 1  !  q/m > 0,  x=0,y=0,z= \pm 1
beta0 = 1.d0

feq1 = 20;  
open(feq1,file='./dat/eq1.dat',access ='append',status='replace')
 
! beta, cs, sgn, eq  are initialized by subroutine init_lf in  module lfmod
! provided the case, beta, and ieq 
call init_lf(beta0, cs0, ieq) 

print*, 'check the assignment with module' !--ckd
print*, 'beta, cs, sgn, eq', beta, cs, sgn, eq

! the differential matrix of the lorentz force with respective to state 
do i = 1, 1000 
  beta = -10 + (i-1) * 0.02
  call dflrtz( eq, beta, cs, sgn, dlf)

!do i = 1, n
! write(*,'(6f8.4)') dlf(i,:) 
!enddo

! check the energy 
!C1 = -3   
!C2 = 2.2894  
!C3 = 1.8899 

call gr_cjlf(eq, cj_eq)
print*,'check energy of equilibrium point!, cj, ieq,', cj_eq, ieq, eq
!read(*,*)  

! set the target energy to be cj_eq
 cj_tar = cj_eq  


! compute the eigenvalues and eigenvectors 
!subroutine eigrg(a,n,isv, wr,wi,vr)
  call eigrg(dlf,n,0, wr,wi,vr)
  write(feq1, '(1f14.8)', advance='no') beta
  
  do j = 1, 6
    if (dabs(wr(j)) < 1.d-8) then 
       write(feq1, '(1f14.1)', advance='no') 0.d0
    else 
       write(feq1, '(1f14.8)', advance='no') wr(j)
    endif
    
    if (dabs(wi(j)) < 1.d-8) then 
       write(feq1, '(1f14.1)', advance='no') 0.d0
    else 
       write(feq1, '(1f14.8)', advance='no') wi(j)
    endif
         
!    read(*,*) 
  enddo  
  write(feq1,*)
!  read(*,*)  
enddo
 
 close(feq1);
 
stop
end program test



















  
