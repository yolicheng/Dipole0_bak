	subroutine deltx(f,g, ipo,dx)
c
C++++++++++++++++++++++++++  
c********************************************************************** 
c   Nov.26,2014
c Modify the definement with 3 free variables(control), 2 target equations
c we have Xk = X(k-1) + deltX(k-1)
c where we require G( deltX(k-1) ) = - F( X(k-1) )

c minimize the norm
c     (deltX)T * Q * deltX
c      Q is a given diagonal positive definite weight matrix
c        here,just set it to identity
c so, the solution of this problem of minima is given by
c     (deltXk-1)= -Q-1 (G)T * (G Q-1 (G)T)-1 * F( X(k-1) ) 
c     (deltXk-1)= -(G)T* (G*(G)T)-1 * F( X(k-1) ) 
c     the unknowns: x,z,vy
c     the equations: f1=vx = 0, f2=vz = 0 for halo and planar lyapunov
c                    f1= z = 0, f2=vx = 0 for vertical lyapunov
 
c    Input variables
c  f(*)	target variables to be set 0, (here vx,vz )
c  g(*,*) 	Jacobian matrix of f(*) 
c    Output variables
c  ipo  1: norm of dx is below the precision, 0: default
c  dx   the difference for the initial state to meet the final target
c
c subroutine used: none
c **********************************************************************	

	implicit real*8(a-h,o-z)
	dimension  f(2),g(2,3),dx(3), g2(2,2),b(2)
 
	prsc = 1.d-13
	ipo = 0
	
 	do 2 i = 1,2
  	do 2 j = 1,2
2 	g2(i,j) = g(i,1)*g(j,1)+ g(i,2)*g(j,2)+ g(i,3)*g(j,3)
 
	det = g2(1,1)*g2(2,2)-g2(1,2)*g2(2,1)
  	
  	b(1) = g2(2,2)*f(1)-g2(1,2)*f(2)
  	b(2) = -g2(2,1)*f(1)+g2(1,1)*f(2)
  	
  	do 5 i = 1,3
 5	dx(i) = -(g(1,i)*b(1)+g(2,i)*b(2))/det
  
	dxm = dx(1)*dx(1) + dx(2)*dx(2)+dx(3)*dx(3)
	if (dsqrt(dxm).lt. prsc)  ipo = 1 
      
      return      
      end 
