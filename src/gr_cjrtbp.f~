        SUBROUTINE GR_CJRTBP(X,CJA) !mu
C******************************************************
C COMPUTATION OF THE JACOBI CONSTANT, CJA, OF A POINT X
C IN THE SPATIAL RTBP.
C IN THIS RTBP THE BIG PRIMARY IS AT (MU,0,0) WITH
C MASS 1-MU, AND THE SMALL ONE AT (MU-1,0,0) WITH
C MASS MU.
C******************************************************
	use pomod, only : mu
        IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION X(6)
        X1=X(1)*X(1)
        X2=X(2)*X(2)
        X3=X(3)*X(3)
        R1=DSQRT((X(1)-MU)**2+X2+X3)
        R2=DSQRT((X(1)-MU+1)**2+X2+X3)
        OME=(X1+X2)*.5+(1-MU)/R1+MU/R2+.5*MU*(1-MU)
        CJA=2*OME-X(4)*X(4)-X(5)*X(5)-X(6)*X(6)
        RETURN
        END

