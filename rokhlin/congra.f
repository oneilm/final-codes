C 
C 
C 
C          DATA SET CONGRA     AT LEVEL 001 AS OF 10/25/83
        implicit real *8 (a-h,o-z)
        REAL *8 A(10000),Y(200),XK(200),PK(200),PKM1(200),
     1         RK(200),err(300)
        real *8 w(10000),WORK(10000)
        EXTERNAL MULT
C 
C        CREATE THE SYSTEM
C 
        CALL PRINI(6,13)
        N=20
C 
        PRINT *,'ENTER N'
        READ *, N
        CALL PRINF('N=*',N,1)
C 
        PRINT *,'ENTER eps'
        READ *, eps
        CALL PRIN2('eps=*',eps,1)
        NUMIT=30
        CALL CREMAT(w,A,Y,N,eps)
        CALL PRIN2('A IS*',A,100)
        CALL PRIN2('AND Y IS*',Y,10)
C 
C        SOLVE THE SYSTEM
C 
        EPS2=1.0d-16
        NUMIT=500
         CALL CONGRA(IER,MULT,Y,N,XK,NUMIT,EPS2,NITER,WORK,ERR,a)
        CALL PRIN2('AND X IS*',XK,N)
         CALL PRINF('IER=*',IER,1)
         CALL PRIN2('ERR IS*',ERR,niter)
        CALL PRINF('NITER IS*',NITER,1)
        STOP
        END
C 
C 
C 
C 
        SUBROUTINE MULT(a,X,Y,n)
        implicit real *8 (a-h,o-z)
        save
        REAL *8 A(N,N),X(1),Y(1),D
        DO 1400 I=1,N
        D=0
        DO 1200 J=1,N
        D=D+A(I,J)*X(J)
 1200 CONTINUE
        Y(I)=D
 1400 CONTINUE
        RETURN
        END
C 
C 
C 
C 
C 
        SUBROUTINE CREMAT(A,b,Y,N,eps)
        implicit real *8 (a-h,o-z)
        save
        REAL *8 A(N,N),Y(1)
ccc     EPS=0.001
        DO 1400 I=1,N
        DO 1200 J=1,N
        A(I,J)=EPS *(i-j)**2
 1200 CONTINUE
        A(I,I)=1.0/I
        Y(I)=20*I
 1400 CONTINUE
C       A(10,10)=0.1
c 
        call matmua(a,a,b,n)
c 
        RETURN
        END
c 
c 
c 
c 
c 
        subroutine matmua(a,b,c,n)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),b(n,n),c(n,n)
c 
        do 2000 i=1,n
        do 1800 j=1,n
        d=0
        do 1600 k=1,n
        d=d+ a(i,k)*b(j,k)
 1600 continue
        c(i,j)=d
 1800 continue
 2000 continue
        return
        end
C 
C 
C 
C 
C 
        SUBROUTINE CONGRA(IER,MULT,Y,N,X,NUMIT,EPS,NITER,
     1          W,ERR,a)
        implicit real *8 (a-h,o-z)
        save
        REAL *8 Y(1),X(1),W(1),a(1),err(1)
         EXTERNAL MULT
c 
c       this subroutine solves a linear system with a
c       symmetric positive definite matrix via the standard
c       conjugate gradient algorithm.
c 
C                      INPUT PARAMETERS:
C 
C   MULT - THE USER-SUPPLIED MATRIX-VECTOR MULTIPLICATION ROUTINE.
C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      DESCRIPTION OF THE CALLING SYSTEM OF SUBROUTINE MULT
C 
C       THE CALLING SEQUENCE OF MULT IS AS FOLLOWS:
C 
C          MULT(A,X,Y,N)
C 
C             WITH THE FOLLOWING PARAMETERS:
C      A - THE MATRIX OF THE SYSTEM (OR WHATEVER OTHER PARAMETER) (INPUT)
C      X - THE VECTOr TO WHICH A IS TO BE APPLIED (INPUT)
C      Y -=A(X) (OUTPUT)
C      N - THE DIMENSIONALITY OF A,X,Y (INPUT)
C 
C      END OF DESCRIPTION OF THE CALLING SYSTEM OF SUBROUTINE MULT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
C   Y - THE RIGHT-HAND SIDE OF THE SYSTEM TO BE SOLVED
C   N - THE DIMENSIONALITY OF THE SYSTEM
C   NUMIT - THE MAXIMUM NUMBER OF ITERATIONS PERMITTED
C   EPS - THED ACCURACY TO WHICH THE SYSTEM IS TO BE SOLVED
C   A - THE MATRIX OF THE SYSTEM (OR WHATEVER OTHER PARAMETER)
C 
C                       OUTPUT PARAMETERS:
C 
C   IER - ERROR RETURN CODE
C     IER=0 MEANS SUCCESSFUL CONCLUSION
C     IER=2 MEANS THAT no element of the right-hand side is greater
c           than 1.0d-60. in this case, x is set to zero and the
c           execurion is terminated
C     IER=8 MEANS THAT THE MAXIMUM PERTMITTED NUMBER OF ITERATIONS
C           (NUMIT) HAS BEEN PERFORMED WITHOUT THE DESIRED PRECISION
C           BEING ACHIEVED
C     IER=16 MEANS THAT THE QUADRATIC FORM FAILED TO DECREASE BEFORE
C           NUMIT ITERATIONS HAVE BEEN PERFORMED OR THE ACCURACY
C           EPS HAS BEEN ACHIEVED
C   X - THE SOLUTION OF THE SYSTEM
C   NITER - THE NUMBER OF ITERATIONS ACTUALLY PERFORMED
C   ERR - THE ARRAY OF ERRORS PRODUCED BY THE ALGORITHM ON VARIOUS STEPS,
C           THAT IS ERR(I) = ||A(X_I)-Y||, WHERE X_I IS THE APPROXIMATION
C           TO THE SOLUTION OBTAINED ON THE I-TH STEP OF THE CONJUGATE
C           GRADIENT PROCESS, AND || * || IS THE L^2 NORM
C 
C                         WORK ARRAY :
C 
C   W - MUST BE AT LEAST 3*N+10 real *8 ELEMENTS LONG.
C 
C        . . . ALLOCATE MEMORY FOR THE CONJUGATE GRADIENT SCHEME
C 
C 
C        . . . ALLOCATE MEMORY FOR WORK ARRAYS
C 
        IPK=1
        IRK=IPK+N
        IW=IRK+N
ccc      call prinf('in congra, n=*',n,1)
C 
C        APPLY THE CG ALGORITHM
C 
        NREC=5
        CALL CONGR1(Y,N,  MULT,X,W(IRK),W(IPK),W(IW),NUMIT,
     1       IER,NREC,NITER,EPS,ERR,a)
         RETURN
         END
C 
C 
C 
        SUBROUTINE CONGR1(Y,N,MULT,XK,RK,PK,W,NUMIT,IER,NREC,
     1      NITER,EPS, err,a)
        implicit real *8 (a-h,o-z)
        save
        REAL *8 Y(1),XK(1),PK(1),RK(1),W(1),EPS,D1,D2,D,D7,D8,
     1        RKMOD2,D9,D29,a(1),err(1)
        EXTERNAL MULT
C 
C        INITIALIZE THE CG PROCESS
C 
         done=1
         zero=0
         D9= 1.0d60
         D29=1.0d-60
C 
C        CHECK IF Y IS ZERO
C 
        IER=0
        DO 800 I=1,N
        IF(dABS(Y(I)).GT. D29) GOTO 1000
  800 CONTINUE
        DO 900 I=1,N
        XK(I)=zero
  900 CONTINUE
        IER=2
        RETURN
 1000 CONTINUE
        DO 1200 I=1,N
        XK(I)=Y(I)
        RK(I)=zero
        PK(I)=zero
 1200 CONTINUE
        BETKM1=zero
ccc      call prinf('before first mult, n=*',n,1)
        CALL MULT(a,XK,W,n)
ccc      call prinf('after first mult, n=*',n,1)
        RKMOD2=zero
        D8=0
        DO 1400 I=1,N
        D8=D8+Y(I)**2
        D=W(I)-Y(I)
        RKMOD2=RKMOD2+D**2
        RK(I)=D
 1400 CONTINUE
        D=dSQRT(RKMOD2)
        D8=done/D8
ccc     CALL PRIN2('AFTER INITIALIZATION THE ERROR IS*',D,1)
C 
C        PERFORM THE ITERATIONS
C 
        DO 4000 K=1,NUMIT
        NITER=K
ccc     CALL PRINF('K=*',K,1)
C 
C        COMPUTE THE CURRENT PK
C 
        DO 2200 I=1,N
        PK(I)=BETKM1*PK(I)-RK(I)
 2200 CONTINUE
C 
C        COMPUTE A(PK)
C 
        CALL MULT(a,PK,W,n)
ccc     call prin2('after second mult, w=*',w,n)
C 
C        COMPUTE  (RK,RK) AND (PK,A(PK))
C 
        D1=0
        D2=0
        DO 2400 I=1,N
        D1=D1+RK(I)**2
        D2=D2+PK(I)*W(I)
 2400 CONTINUE
C 
C        CHECK IF A(XK)-YX IS ZERO
C 
        IF(D1.GT.D29)   GOTO 2500
        IER=1
        RETURN
 2500 CONTINUE
C 
C        COMPUTE XKP1,RKP1
C 
        ALFA=D1/D2
        DO 2600 I=1,N
        XK(I)=XK(I)+ALFA*PK(I)
        RK(I)=RK(I)+ALFA*W(I)
 2600 CONTINUE
C 
C        IF THE TIME HAS COME TO RECOMPUTE RK - DO SO
C 
C       GOTO 2700
        J=K/NREC
        J=K-J*NREC
        IF(J.NE.0) GOTO 2700
        CALL MULT(a,XK,W,n)
        D3=0
        D2=0
        DO 2650 I=1,N
        D3=D3+W(I)*XK(I)
        D2=D2+XK(I)*Y(I)
        RK(I)=W(I)-Y(I)
 2650 CONTINUE
        D3=D3/2-D2
C       CALL PRIN2('THE JUNK BEING MINIMIZED IS*',D3,1)
        D11=D9-D3
C       CALL PRIN2('THE DIFFERENCE IS*',D11,1)
        IF(D3.LT.D9) GOTO 2670
        IER=16
        GOTO 4200
 2670 CONTINUE
        D9=D3
 2700 CONTINUE
C 
C        COMPUTE BETA
C 
        D=0
        DO 2800 I=1,N
        D=D+RK(I)**2
 2800 CONTINUE
        D7=dSQRT(D)*D8
        err(k)=d7
        IF(D7.LT.EPS) GOTO 4200
C       CALL PRIN2('THE ERROR IS*',D7,1)
        BETKM1=D/D1
 4000 CONTINUE
        IER=8
 4200 CONTINUE
        RETURN
        END
