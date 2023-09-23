        IMPLICIT REAL *8 (A-H,O-Z)
        DIMENSION FF(1000),FF2(1000)
        CALL PRINI(6,13)
C 
C 
        DONE=1
        HALF=DONE/2
        x=2
        gam=dgamma(x)
        call prin2('x=*',x,1)
        call prin2('gam=*',gam,1)
        call prin2('gam-1=*',gam-1,1)
        STOP
        END
C 
C 
C 
C 
C 
        function dGAMMA(X)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        data i99/-7/
c 
c       if the need be - initialize
c 
        if(i99 .eq. 11) goto 1100
        DONE=1
        PI=DATAN(DONE)*4
        COEF=DSQRT(2*PI)
        E=DEXP(DONE)
        i99=11
 1100 continue
C 
C       EVALUATE THE GAMMA FUNCTION F OF X
C 
C 
C       . . .  DETERMINE THE SHIFT NECESSARY TO GET STIRLING
C              FORMULA TO WORK
C 
        LEN=130
        ISHFT=LEN-X+1
        IF(ISHFT .LT. 0) ISHFT=0
        Z=X+ISHFT
C 
C       EVALUATE GAMMA FUNCTION AT THE SHIFTED POINT VIA
C       STIRLING FORMULA
C 
        D=(Z/E) **Z / DSQRT(Z) * COEF
        F=1  +  1 / (12*Z)  +  1 / (288*Z**2)
     1    -  139 / (51840 * Z**3)  -  571 / (2488320 * Z**4 )
        F=F*D
C 
C       SHIFT THE VALUE OF THE GAMMA FUNCTION
C 
        dgamma=f
        IF (ISHFT .EQ. 0) RETURN
        DO 1200 I=1,ISHFT
        F=F/(ISHFT-I+X   )
 1200 CONTINUE
        dgamma=f
        RETURN
        end
