! to test if there is something wrong with the computation of the variational matrix 
use lfmod

implicit none 
!  integer, parameter :: dp=kind(1.d0), n = 42
integer, parameter ::  neq = 42
    
integer ::   ftag, tdir, i, cs, ieq  
real(kind=dp) :: y0(neq),t0, tf
real(kind=dp) :: y(neq) ! As the final state
  

! Local Variables
real(kind=dp) :: r(13,neq), b(neq), f(neq), t, h, hmin, hmax, e,  & ! rk78 
		 cj , po0(6), phi(6,6),  & 
		 beta, &
		 dy(6), phi1(6), yck(neq), vf(neq)

!  external :: gr_lf, gr_cjlf
 cs = 1   
ieq = 3  
beta = 10.d0; 

call init_lfmod
call init_lf(beta, cs, ieq) 

! 0.56123418697415695        0.0000000000000000       0.56122987432829607        0.0000000000000000        9.9943356318814896E-005   0.0000000000000000 

! 0.42413555  0.56123418  0.00000000  0.56122988  0.00000000  0.00009994  0.00000000  1.88988156

tdir = 1
tf = 0.42413555/2 ! the period
!  tf = 1.d-0

po0 = (/0.56123418d0,  0.d0,   0.56122988d0,  0.d0,  0.00009994d0, 0.d0/)  ! the first p.o. 
 
!po0 = (/0.56123102415468651d0, 0.d0, 0.56123102415468651d0,  0.d0,  0.d0, 0.d0/)  ! eq2
!po0 = (/0.d0, 0.d0, 1.d0,  0.d0,  0.d0, 0.d0/)  ! eq1
 
ftag = 22 
open(ftag,file='./dat/testphi.dat', access ='append',status='replace')
 
! Integrate the periodic orbits by gr_rk78, write the x, y in txt file with ftag
t =  0.d0
t0 = 0.d0
h =  1.d-3! maybe small value is better, for a p.o. with small period, the orbit could be quite corse 

! specify the error control 
hmin = 1.d-10
hmax = 1.d0
e    = 1.d-14

y = 0.d0
y(1:6) = po0 
y(7:42:7) = 1.d0  
  
 
! Integrate the orbit in interval [t0, tf]
do while( dabs(t+h-t0) .lt. dabs(tf-t0) )
  
! check the vector field
  call  gr_lf(t, y, neq, vf)
  print*, 'check the variatioal matrix'
  phi =  reshape(vf(7:42), (/6,6/))
  do i = 1, 6
    write(*,'(6f12.8)') phi(i,:) 
  enddo
  print*
! --------------------------------------


  call gr_cjlf(y(1:6), cj)
    
  write(ftag,'(8e20.10)')  t, y(1:6), cj
  write(*,'(8e20.10)')  t, y(1:6), cj!ck
  call gr_rk78(t,y, neq ,h,hmin,hmax,e,r,b,f, gr_lf )
  
    
  ! check the state transition matrixf
  print*; print*, 'Phi- t=',  t
  phi = reshape(y(7:42), (/6,6/)) ! better than equivalence declaration.... 
  do i = 1, 6
    write(*,'(6e20.10)')  phi(i,:)
  enddo 
!    read*
  
    
enddo
     
!if we want a periodic orbit with exact 1 period, a control of tf must be made
if(dabs(t-tf) .gt. 1.d-9) then
  
  h =  tdir*tf - t ! for the stable orbit, the final time should be -tf
! to make sure the next step can be executed within allowable interval [hmin hmax]
!    hmin should be dabs(h)
  call gr_rk78(t,y, neq,h, dabs(h), hmax, e,r,b,f, gr_lf)
    
 endif 

! save the last point  
call gr_cjlf(y(1:6), cj)      
write(ftag,'(8e20.10)') t, y(1:6), cj
  
write(ftag,*)  ! better to save as a block than index
  

  
! test with difference  
yck = 0.d0
yck(1:6) = po0 
yck(7:42:7) = 1.d0
yck(1) = yck(1)+1.d-5 
!  yck(2) = yck(2)+1.d-5 
  
  
t = 0.d0
h = 1.d-3
do while( dabs(t+h-t0) .lt. dabs(tf-t0) )
  write(ftag,'(7e20.10)')  t, y(1:6) 
  write(*,'(7e20.10)')  t, y(1:6) 
  call gr_rk78(t,yck, neq ,h,hmin,hmax,e,r,b,f, gr_lf )
enddo
  
!if you want the strict periodic orbit, a control of tf must be made
if(dabs(t-tf) .gt. 1.d-9) then
  h =  tdir*tf - t 
  call gr_rk78(t,yck, neq,h, dabs(h), hmax, e,r,b,f, gr_lf)
endif 
  
! by integration
print*
print*, 'By integration, phi'
do i = 1, 6
  write(*,'(6e20.10)') phi(i,:)
enddo 
  
print*
  
! the difference between yck and y
dy = yck(1:6) - y(1:6)
phi1 = dy / 1.d-5
print*, 'By difference, difference in state -- first column of phi'
!print*, dy
print*, phi1
  
read*
  
  
end  
