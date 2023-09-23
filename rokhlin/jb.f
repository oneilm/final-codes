        IMPLICIT REAL *8 (A-H,O-Z)
        COMPLEX *16 FUNS(40000),Z,Z1,Z2,Z3,IMA,dz,fun1,df
        real *8 dfuns(40000),ts(40000)
        DATA IMA/(0.0D0,1.0D0)/
        CALL PRINI(6,13)
        Z=9
c 
         PRINT *, 'ENTER X '
         READ *,X
         call prin2('x=*',x,1)
c 
         PRINT *, 'ENTER eps2'
         READ *,eps2
         call prin2('eps2=*',eps2,1)
        Z=10 * ima +x
         z=x * ima
         z=x
        EPS=1.0D-20
        CALL JBFUN(Z,FUNS(2),N,EPS)
        CALL PRINF ('N= *',N,1)
        CALL PRIN2('Z=*',Z,2)
CCC     IF(2.NE.3) STOP
        CALL PRIN2('FUNS=*',FUNS,N*2)
 1200 FORMAT('AND MORE ACCURATELY',4(2X,E22.15) )
cccc         WRITE(6,1200) (FUNS(I),I=1,N)
cccc         WRITE(13,1200) (FUNS(I),I=1,N)
c 
c       now, find the number after which all functions are
c       less than eps2
c 
        jj=0
        do 2200 i=1,n-1
        d=cdabs(funs(i+1))
        if (d.gt. eps2) jj=jj+1
 2200 continue
        call prinf('number of js greater than eps2 is*',
     1    jj,1)
c 
c       scale the bessel functions and plot them
c 
        imax=1
        dmax=cdabs(funs(1))
        do 2400 i=1,n
        if(dmax .gt. cdabs(funs(i)) ) goto 2400
        imax=i
        dmax=cdabs(funs(i))
 2400 continue
         call prinf('imax=*',imax,1)
         call prin2('dmax=*',dmax,1)
c 
         done=1
        do 2600 i=1,n
ccc        funs(i)=funs(i)/dexp(x)
        funs(i)=funs(i)/dmax
        ts(i)=i
        dfuns(i)=cdabs(funs(i))
        if(dfuns(i) .lt. 1.0d-15) dfuns(I)=0
 2600 continue
        call prin2('and rescaled funs are*',funs,n*2)
c 
c        plot them things
c 
        iw=21
c 
        call RSGRAF(ts,dfuns,N,IW,'scaled i-functions*')
  
  
        STOP
        END
C 
C 
C 
C 
C 
       SUBROUTINE JBFUN(Z,J,N,EPS)
C 
C      THIS SUBROUTINE COMPUTES ARRAY OF VALUES OF BESSEL
C      FUNCTIONS "J" J<0>(Z),...,J<N>(Z) FOR A COMPLEX ARGUMENT
C      "Z".
C      J<I> I=N+1,... ARE CONSIDERED TO BE 0 (N .GE. 3).
C      FOR I .LE. N THE SMALLER "I" THE MORE PRECISE THE VALUE
C      OF J<I>(Z) IS.
C 
C      INPUT
C      Z - COMPLEX*16 ARGUMENT;
C      EPS - VALUE BASED ON WHICH WE DETERMINE "N" SUCH THAT
C            J<I>(Z) .LT. EPS FOR I=N+1,... .
C            ACTUALLY J<N>(Z) IS MUCH SMALLER THAN EPS;
C            THE SMALLER ABS(Z) THE MORE SIGNIFICANT THE
C            DIFFERENCE IS.
C 
C      OUTPUT
C      N - INTEGER NUMBER SUCH THAT J<I>(Z) .LT. EPS
C          FOR I=N+1,... (N .GE. 3)
C      J - ARRAY OF J<0>(Z),...,J<N>(Z);
C 
C      REMARK
C      ARRAY "J" MUST HAVE 0 SUBSCRIPT.
C 
        save
       INTEGER     K,N,I0
       REAL*8      EPS,BIG,REA(2)
       COMPLEX*16  Z,J(1),SUM,SC,IM,CR1,CTWO,COM
       COMPLEX*8   Z8
        EQUIVALENCE(REA(1),COM)
       DATA        I0 /0/,IM /(0.0D0,1.0D0)/,CTWO /(2.0D0,0.0D0)/
C 
C      IF (CABS(Z) .GT. 1.0D-30) GO TO 11
        COM=Z
        IF(DABS(REA(1))+DABS(REA(2))  .GT. 1.0D-30) GOTO 11
       J(I0) = 1
        J(1) = 0
        J(2) = 0
        J(3) = 0
       N = 3
C 
       GO TO 999
C 
11     CONTINUE
C 
       BIG = 1/EPS
C 
       CALL FINDK(K,Z,J,BIG)
C 
       CALL JREC(K,Z,J)
C       CALL PRIN2('AFTER JREC J IS*',J,K*2)
C 
       Z8 = Z
       N = K + 1
       N0 = N/4
       SUM = J(I0)
C 
       NR = N - 4*N0
       IF (NR .LE. 0) GO TO 44
       DO 444 I=1,NR
       SUM = SUM + J(4*N0+I)
444    CONTINUE
C 
44      CONTINUE
C 
       IF (AIMAG(Z8) .GE. 0) GO TO 22
C 
       CR1 =CDEXP(Z*IM)
C 
       DO 111 I=1,N0
C      SUM = SUM +  J(I)*(T**I + ((-1)/T)**I)
       SUM = SUM + J(4*I-3)*CTWO*IM
       SUM = SUM - J(4*I-2)*CTWO
       SUM = SUM - J(4*I-1)*CTWO*IM
       SUM = SUM + J(4*I)*CTWO
111    CONTINUE
C 
       GO TO 33
C 
22     CONTINUE
       CR1 =CDEXP(-Z*IM)
C 
       DO 122 I=1,N0
       SUM = SUM - J(4*I-3)*CTWO*IM
       SUM = SUM - J(4*I-2)*CTWO
       SUM = SUM + J(4*I-1)*CTWO*IM
       SUM = SUM + J(4*I)*CTWO
122    CONTINUE
C 
33     CONTINUE
       SC = SUM/CR1
C 
C         WRITE (8,1122) (SUM,CR1,J(I0),(I,J(I),I=1,N))
C         WRITE (13,1122) (SUM,CR1,J(I0),(I,J(I),I=1,N))
C1122     FORMAT (//' ','SUM=',D12.5,1X,D12.5/
C    *             ' ','EXP=',D12.5,1X,D12.5/
C    *              ' ','J(0)=',D12.5,1X,D12.5/
C    *             ' ',('I=',I4,2X,'J(I)=',D12.5,1X,D12.5))
C 
C 
       J(I0) = J(I0)/SC
C 
       DO 222 I=1,N
       J(I) = J(I)/SC
222    CONTINUE
C 
C 
C       REDUCE THE NUMBER OF J-FUNCTIONS TO THE CORRECT ONE
C 
        K=1
        N1=N-1
        DO 1200 I=1,N1
        IF(CDABS(J(I)).GE.EPS) K=I
 1200 CONTINUE
        N = MAX0(K,3)
C 
999     CONTINUE
       RETURN
       END
C 
C 
C 
C 
C 
       SUBROUTINE JREC(K,Z,J)
C 
C      THIS SUBROUTINE COMPUTES "UNSCALED" BESSEL FUNCTIONS J<K> OF
C      COMPLEX ARGUMENT "Z".
C      THE FOLLOWING RECURSION RELATIONS ARE USED:
C      J<K-1> = (2*K/Z)*J<K> - J<K+1>
C 
        save
       INTEGER     K,I,K1
       COMPLEX*16  Z,J(1),NLNJK
       REAL*8      PI,E
        DATA PI/3.1415 92653 58979 32384 D0/,E/2.7182 81828 45904 52353
     1       D0/
C 
CCC    PI = DATAN(1.0D0)*4
CCC    E =  DEXP(1.0D0)
       K1 = K + 1
C 
CCC    J(K+1) = ((E*Z/(2*K1))**K1)/DSQRT(2*PI*K1)
CCC    J(K) =((E*Z/(2*K))**K)/DSQRT(2*PI*K)
        J(K+1)=NLNJK(Z,K+1)
        J(K)=NLNJK(Z,K)
C       CALL PRIN2('J(K+1) IS*',J(K+1),2)
C       CALL PRIN2('AND J(K) IS*',J(K),2)
C 
       DO 222 I=1,K
       J(K-I) = (2*(K-I+1)/Z)*J(K-I+1) - J(K-I+2)
222    CONTINUE
C 
C 
       RETURN
       END
C 
C 
C 
C 
C 
       SUBROUTINE FINDK(K,Z,J,BIG)
C 
C      THE SUBROUTINE DETERMINES "K" SUCH THAT BESSEL FUNCTIONS Y<K>
C      OF COMPLEX ARGUMENT "Z" IS .GE. "BIG" (K .GE. 3).
C 
        save
       INTEGER     K
       REAL*8      J(1),R,BIG
       COMPLEX*16  Z
C 
       R =CDABS(Z)
       I0 = 0
       J(I0) = 1
       J(1) = 0
       I = 1
C 
11     CONTINUE
       J(I+1) = (2*I/R)*J(I) - J(I-1)
       IF (DABS(J(I+1)) .GE. BIG) GO TO 22
       I = I + 1
       GO TO 11
C 
22     CONTINUE
       K = I + 1
        K = MAX0(K,3)
C 
C 
       RETURN
       END
C 
C 
C 
C 
C 
        FUNCTION NLNJK(W,N)
        save
        COMPLEX *16 W,Z,CNUM,CDEN,NLNJK,CDS1,CDS2
        REAL *8 PI2,DONE
        DATA PI2/6.2831853D0/,DONE/1.0D0/
C       CALL PRIN2(' INSIDE NLNJK W IS*',W,2)
C       CALL PRINF('AND N IS*',N,1)
C 
C        COMPUTE Z
C 
        Z=W/N
C 
C        EVALUATE THE APPROXIMATION TO J  (W)=J  (N*Z)
C                                       N      N
C 
C 
        CDS1=CDSQRT(DONE-Z**2)
        CDS2=CDLOG(DONE+CDS1)
        CNUM=CDEXP(N*(CDLOG(Z)+CDS1-CDS2))
        CDEN=CDSQRT(PI2*N*CDS1)
C         CNUM=CEXP(N*(CLOG(Z)+CSQRT(DONE-Z**2)))
C       CALL PRIN2('INSIDE NLNJK CNUM IS*',CNUM,2)
C        CDEN=CSQRT(PI2*N*CSQRT(DONE-Z**2)) *(DONE
C    1    +CSQRT(DONE-Z**2))**N
C       CALL PRIN2('CDEN IS*',CDEN,2)
       NLNJK=CNUM/CDEN
        RETURN
        END
