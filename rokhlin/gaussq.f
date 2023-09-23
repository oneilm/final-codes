      IMPLICIT REAL *8 (A-H,O-Z)
        DIMENSION T(2000),W(2000),B(2000),ENDPTS(2),F(2000)
        CALL PRINI(6,13)
C 
C       CREATE THE INPUT PARAMETERS
C 
        N=128
       KIND=5
c 
         PRINT *, 'ENTER n'
         READ *,n
         CALL PRINF('n=*',n,1 )
c 
         PRINT *, 'ENTER k'
         READ *,k
         CALL PRINF('k=*',k,1 )
c 
        ALPHA=2
        BETA=2
        KPTS=0
C 
C       CONSTRUCT GAUSSIAN QUADRATURES
C 
        call prinf('before gaussq, n=*',n,1)
        CALL GAUSSQ(KIND, N, ALPHA, BETA, KPTS, ENDPTS, B, T, W)
        CALL PRIN2('T =*',T,N)
        CALL PRIN2('W =*',W,N)
C 
C       CONSTRUCT THE TEST FUNCTION AND INTEGRATE IT
C 
        do 1200 i=1,n
        call funt(t(i),k,f(i))
 1200 continue
c 
        rint=0
        do 1400 i=1,n
        rint=rint+w(i)*f(i)
 1400 continue
        call prin2('rint as calculated*',rint,1)
c 
c      . . . and analytically
c 
        done=1
        call fint(done,k,finp)
        call fint(-done,k,finm)
         call prin2('finp=*',finp,1)
         call prin2('finm=*',finm,1)
c 
        call prin2('and analytically,*',finp-finm,1)
        call prin2('and the difference is*',finp-finm-rint,1)
        STOP
        END
c 
c 
c 
c 
c 
       subroutine funt(x,k,f)
       implicit real *8 (a-h,o-x)
c 
        save
       f=x**k
       return
c 
c 
c 
c 
        entry fint(x,k,fin)
        done=1
        fin=done/(k+1)*x**(k+1)-2*done/(k+3)*x**(k+3)
     1      +done/(k+5)*x**(k+5)
        return
        end
C 
C 
C 
C 
C 
        FUNCTION FUN(T)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        FUN=dsin(PI*T*199)
        FUN=t**4
        RETURN
C 
C 
C 
C 
        ENTRY FUNI(IER)
        DONE=1
        PI=DATAN(DONE)*4
        FUNI=PI
        RETURN
        END
C 
C 
C 
C 
C 
C 
C          DATA SET GAUSSQ     AT LEVEL 001 AS OF 07/23/82
C                                1/20/75
      SUBROUTINE GAUSSQ(KIND, N, ALPHA, BETA, KPTS, ENDPTS, B, T, W)
C 
C           THIS SET OF ROUTINES COMPUTES THE NODES T(J) AND WEIGHTS
C        W(J) FOR GAUSSIAN-TYPE QUADRATURE RULES WITH PRE-ASSIGNED
C        NODES.  THESE ARE USED WHEN ONE WISHES TO APPROXIMATE
C 
C                 INTEGRAL (FROM A TO B)  F(X) W(X) DX
C 
C                              N
C        BY                   SUM W  F(T )
C                             J=1  J    J
C 
C        (NOTE W(X) AND W(J) HAVE NO CONNECTION WITH EACH OTHER.)
C        HERE W(X) IS ONE OF SIX POSSIBLE NON-NEGATIVE WEIGHT
C        FUNCTIONS (LISTED BELOW), AND F(X) IS THE
C        FUNCTION TO BE INTEGRATED.  GAUSSIAN QUADRATURE IS PARTICULARLY
C        USEFUL ON INFINITE INTERVALS (WITH APPROPRIATE WEIGHT
C        FUNCTIONS), SINCE THEN OTHER TECHNIQUES OFTEN FAIL.
C 
C           ASSOCIATED WITH EACH WEIGHT FUNCTION W(X) IS A SET OF
C        ORTHOGONAL POLYNOMIALS.  THE NODES T(J) ARE JUST THE ZEROES
C        OF THE PROPER N-TH DEGREE POLYNOMIAL.
C 
C     INPUT PARAMETERS (ALL REAL NUMBERS ARE IN DOUBLE PRECISION)
C 
C        KIND     AN INTEGER BETWEEN 1 AND 6 GIVING THE TYPE OF
C                 QUADRATURE RULE:
C 
C        KIND = 1:  LEGENDRE QUADRATURE, W(X) = 1 ON (-1, 1)
C        KIND = 2:  CHEBYSHEV QUADRATURE OF THE FIRST KIND
C                   W(X) = 1/SQRT(1 - X*X) ON (-1, +1)
C        KIND = 3:  CHEBYSHEV QUADRATURE OF THE SECOND KIND
C                   W(X) = SQRT(1 - X*X) ON (-1, 1)
C        KIND = 4:  HERMITE QUADRATURE, W(X) = EXP(-X*X) ON
C                   (-INFINITY, +INFINITY)
C        KIND = 5:  JACOBI QUADRATURE, W(X) = (1-X)**ALPHA * (1+X)**
C                   BETA ON (-1, 1), ALPHA, BETA .GT. -1.
C                   NOTE: KIND=2 AND 3 ARE A SPECIAL CASE OF THIS.
C        KIND = 6:  GENERALIZED LAGUERRE QUADRATURE, W(X) = EXP(-X)*
C                   X**ALPHA ON (0, +INFINITY), ALPHA .GT. -1
C 
C        N        THE NUMBER OF POINTS USED FOR THE QUADRATURE RULE
C        ALPHA    REAL PARAMETER USED ONLY FOR GAUSS-JACOBI AND GAUSS-
C                 LAGUERRE QUADRATURE (OTHERWISE USE 0.D0).
C        BETA     REAL PARAMETER USED ONLY FOR GAUSS-JACOBI QUADRATURE--
C                 (OTHERWISE USE 0.D0)
C        KPTS     (INTEGER) NORMALLY 0, UNLESS THE LEFT OR RIGHT END-
C                 POINT (OR BOTH) OF THE INTERVAL IS REQUIRED TO BE A
C                 NODE (THIS IS CALLED GAUSS-RADAU OR GAUSS-LOBATTO
C                 QUADRATURE).  THEN KPTS IS THE NUMBER OF FIXED
C                 ENDPOINTS (1 OR 2).
C        ENDPTS   REAL ARRAY OF LENGTH 2.  CONTAINS THE VALUES OF
C                 ANY FIXED ENDPOINTS, IF KPTS = 1 OR 2.
C        B        REAL SCRATCH ARRAY OF LENGTH N
C 
C     OUTPUT PARAMETERS (BOTH DOUBLE PRECISION ARRAYS OF LENGTH N)
C 
C        T        WILL CONTAIN THE DESIRED NODES.
C        W        WILL CONTAIN THE DESIRED WEIGHTS W(J).
C 
C     SUBROUTINES REQUIRED
C 
C        SOLVE, CLASS, AND IMTQL2 ARE PROVIDED.  UNDERFLOW MAY SOMETIMES
C        OCCUR, BUT IT IS HARMLESS IF THE UNDERFLOW INTERRUPTS ARE
C        TURNED OFF.  TO DO THIS, THE FIRST CALL OF THE MAIN PROGRAM
C        SHOULD BE
C                  CALL TRAPS (0, 0, 10000, 0, 0)    IN WATFIV
C        OR
C                  CALL INIT                         IN FORTRAN G OR H.
C 
C     ACCURACY
C 
C        THE ROUTINE WAS TESTED UP TO N = 512 FOR LEGENDRE QUADRATURE,
C        UP TO N = 136 FOR HERMITE, UP TO N = 68 FOR LAGUERRE, AND UP
C        TO N = 10 OR 20 IN OTHER CASES.  IN ALL BUT TWO INSTANCES,
C        COMPARISON WITH TABLES IN REF. 3 SHOWED 12 OR MORE SIGNIFICANT
C        DIGITS OF ACCURACY.  THE TWO EXCEPTIONS WERE THE WEIGHTS FOR
C        HERMITE AND LAGUERRE QUADRATURE, WHERE UNDERFLOW CAUSED SOME
C        VERY SMALL WEIGHTS TO BE SET TO ZERO.  THIS IS, OF COURSE,
C        COMPLETELY HARMLESS.
C 
C     METHOD
C 
C           THE COEFFICIENTS OF THE THREE-TERM RECURRENCE RELATION
C        FOR THE CORRESPONDING SET OF ORTHOGONAL POLYNOMIALS ARE
C        USED TO FORM A SYMMETRIC TRIDIAGONAL MATRIX, WHOSE
C        EIGENVALUES (DETERMINED BY THE IMPLICIT QL-METHOD WITH
C        SHIFTS) ARE JUST THE DESIRED NODES.  THE FIRST COMPONENTS OF
C        THE ORTHONORMALIZED EIGENVECTORS, WHEN PROPERLY SCALED,
C        YIELD THE WEIGHTS.  THIS TECHNIQUE IS MUCH FASTER THAN USING A
C        ROOT-FINDER TO LOCATE THE ZEROES OF THE ORTHOGONAL POLYNOMIAL.
C        FOR FURTHER DETAILS, SEE REF. 1.  REF. 2 CONTAINS DETAILS OF
C        GAUSS-RADAU AND GAUSS-LOBATTO QUADRATURE ONLY.
C 
C     REFERENCES
C 
C        1.  GOLUB, G. H., AND WELSCH, J. H., "CALCULATION OF GAUSSIAN
C            QUADRATURE RULES," MATHEMATICS OF COMPUTATION 23 (APRIL,
C            1969), PP. 221-230.
C        2.  GOLUB, G. H., "SOME MODIFIED MATRIX EIGENVALUE PROBLEMS,"
C            SIAM REVIEW 15 (APRIL, 1973), PP. 318-334 (SECTION 7).
C        3.  STROUD AND SECREST, GAUSSIAN QUADRATURE FORMULAS, PRENTICE-
C            HALL, ENGLEWOOD CLIFFS, N.J., 1966.
C 
C     ..................................................................
C 
        save
      DOUBLE PRECISION B(N), T(N), W(N), ENDPTS(2), MUZERO, T1,
     X GAM, SOLVE, DSQRT, ALPHA, BETA
C 
      CALL CLASS (KIND, N, ALPHA, BETA, B, T, MUZERO)
C 
C           THE MATRIX OF COEFFICIENTS IS ASSUMED TO BE SYMMETRIC.
C           THE ARRAY T CONTAINS THE DIAGONAL ELEMENTS, THE ARRAY
C           B THE OFF-DIAGONAL ELEMENTS.
C           MAKE APPROPRIATE CHANGES IN THE LOWER RIGHT 2 BY 2
C           SUBMATRIX.
C 
      IF (KPTS.EQ.0)  GO TO 100
      IF (KPTS.EQ.2)  GO TO  50
C 
C           IF KPTS=1, ONLY T(N) MUST BE CHANGED
C 
      T(N) = SOLVE(ENDPTS(1), N, T, B)*B(N-1)**2 + ENDPTS(1)
      GO TO 100
C 
C           IF KPTS=2, T(N) AND B(N-1) MUST BE RECOMPUTED
C 
   50 GAM = SOLVE(ENDPTS(1), N, T, B)
      T1 = ((ENDPTS(1) - ENDPTS(2))/(SOLVE(ENDPTS(2), N, T, B) - GAM))
      B(N-1) = DSQRT(T1)
      T(N) = ENDPTS(1) + GAM*T1
C 
C           NOTE THAT THE INDICES OF THE ELEMENTS OF B RUN FROM 1 TO N-1
C           AND THUS THE VALUE OF B(N) IS ARBITRARY.
C           NOW COMPUTE THE EIGENVALUES OF THE SYMMETRIC TRIDIAGONAL
C           MATRIX, WHICH HAS BEEN MODIFIED AS NECESSARY.
C           THE METHOD USED IS A QL-TYPE METHOD WITH ORIGIN SHIFTING
C 
  100 W(1) = 1.0D0
      DO 105 I = 2, N
  105    W(I) = 0.0D0
C 
      CALL IMTQL2 (N, T, B, W, IERR)
      DO 110 I = 1, N
  110    W(I) = MUZERO * W(I) * W(I)
C 
      RETURN
      END
C 
C 
C 
      DOUBLE PRECISION FUNCTION SOLVE(SHIFT, N, A, B)
C 
C       THIS PROCEDURE PERFORMS ELIMINATION TO SOLVE FOR THE
C       N-TH COMPONENT OF THE SOLUTION DELTA TO THE EQUATION
C 
C             (JN - SHIFT*IDENTITY) * DELTA  = EN,
C 
C       WHERE EN IS THE VECTOR OF ALL ZEROES EXCEPT FOR 1 IN
C       THE N-TH POSITION.
C 
C       THE MATRIX JN IS SYMMETRIC TRIDIAGONAL, WITH DIAGONAL
C       ELEMENTS A(I), OFF-DIAGONAL ELEMENTS B(I).  THIS EQUATION
C       MUST BE SOLVED TO OBTAIN THE APPROPRIATE CHANGES IN THE LOWER
C       2 BY 2 SUBMATRIX OF COEFFICIENTS FOR ORTHOGONAL POLYNOMIALS.
C 
C 
      DOUBLE PRECISION SHIFT, A(N), B(N), ALPHA
C 
      ALPHA = A(1) - SHIFT
      NM1 = N - 1
      DO 10 I = 2, NM1
   10    ALPHA = A(I) - SHIFT - B(I-1)**2/ALPHA
      SOLVE = 1.0D0/ALPHA
      RETURN
      END
C 
C 
C 
      SUBROUTINE CLASS(KIND, N, ALPHA, BETA, B, A, MUZERO)
C 
C           THIS PROCEDURE SUPPLIES THE COEFFICIENTS A(J), B(J) OF THE
C        RECURRENCE RELATION
C 
C             B P (X) = (X - A ) P   (X) - B   P   (X)
C              J J            J   J-1       J-1 J-2
C 
C        FOR THE VARIOUS CLASSICAL (NORMALIZED) ORTHOGONAL POLYNOMIALS,
C        AND THE ZERO-TH MOMENT
C 
C             MUZERO = INTEGRAL W(X) DX
C 
C        OF THE GIVEN POLYNOMIAL'S WEIGHT FUNCTION W(X).  SINCE THE
C        POLYNOMIALS ARE ORTHONORMALIZED, THE TRIDIAGONAL MATRIX IS
C        GUARANTEED TO BE SYMMETRIC.
C 
C           THE INPUT PARAMETER ALPHA IS USED ONLY FOR LAGUERRE AND
C        JACOBI POLYNOMIALS, AND THE PARAMETER BETA IS USED ONLY FOR
C        JACOBI POLYNOMIALS.  THE LAGUERRE AND JACOBI POLYNOMIALS
C        REQUIRE THE GAMMA FUNCTION.
C 
C     ..................................................................
C 
        save
      DOUBLE PRECISION A(N), B(N), MUZERO, ALPHA, BETA
      DOUBLE PRECISION ABI, A2B2, DGAMMA, PI, DSQRT, AB
      DATA PI / 3.141592653589793D0/
C 
      NM1 = N - 1
      GO TO (10, 20, 30, 40, 50, 60), KIND
C 
C              KIND = 1:  LEGENDRE POLYNOMIALS P(X)
C              ON (-1, +1), W(X) = 1.
C 
   10 MUZERO = 2.0D0
      DO 11 I = 1, NM1
         A(I) = 0.0D0
         ABI = I
   11    B(I) = ABI/DSQRT(4*ABI*ABI - 1.0D0)
      A(N) = 0.0D0
      RETURN
C 
C              KIND = 2:  CHEBYSHEV POLYNOMIALS OF THE FIRST KIND T(X)
C              ON (-1, +1), W(X) = 1 / SQRT(1 - X*X)
C 
   20 MUZERO = PI
      DO 21 I = 1, NM1
         A(I) = 0.0D0
   21    B(I) = 0.5D0
      B(1) = DSQRT(0.5D0)
      A(N) = 0.0D0
      RETURN
C 
C              KIND = 3:  CHEBYSHEV POLYNOMIALS OF THE SECOND KIND U(X)
C              ON (-1, +1), W(X) = SQRT(1 - X*X)
C 
   30 MUZERO = PI/2.0D0
      DO 31 I = 1, NM1
         A(I) = 0.0D0
   31    B(I) = 0.5D0
      A(N) = 0.0D0
      RETURN
C 
C              KIND = 4:  HERMITE POLYNOMIALS H(X) ON (-INFINITY,
C              +INFINITY), W(X) = EXP(-X**2)
C 
   40 MUZERO = DSQRT(PI)
      DO 41 I = 1, NM1
         A(I) = 0.0D0
   41    B(I) = DSQRT(I/2.0D0)
      A(N) = 0.0D0
      RETURN
C 
C              KIND = 5:  JACOBI POLYNOMIALS P(ALPHA, BETA)(X) ON
C              (-1, +1), W(X) = (1-X)**ALPHA + (1+X)**BETA, ALPHA AND
C              BETA GREATER THAN -1
C 
   50 AB = ALPHA + BETA
      ABI = 2.0D0 + AB
      MUZERO = 2.0D0 ** (AB + 1.0D0) * DGAMMA(ALPHA + 1.0D0) * DGAMMA(
     X BETA + 1.0D0) / DGAMMA(ABI)
      A(1) = (BETA - ALPHA)/ABI
      B(1) = DSQRT(4.0D0*(1.0D0 + ALPHA)*(1.0D0 + BETA)/((ABI + 1.0D0)*
     1  ABI*ABI))
      A2B2 = BETA*BETA - ALPHA*ALPHA
      DO 51 I = 2, NM1
         ABI = 2.0D0*I + AB
         A(I) = A2B2/((ABI - 2.0D0)*ABI)
   51    B(I) = DSQRT (4.0D0*I*(I + ALPHA)*(I + BETA)*(I + AB)/
     1   ((ABI*ABI - 1)*ABI*ABI))
      ABI = 2.0D0*N + AB
      A(N) = A2B2/((ABI - 2.0D0)*ABI)
      RETURN
C 
C              KIND = 6:  LAGUERRE POLYNOMIALS L(ALPHA)(X) ON
C              (0, +INFINITY), W(X) = EXP(-X) * X**ALPHA, ALPHA GREATER
C              THAN -1.
C 
   60 MUZERO = DGAMMA(ALPHA + 1.0D0)
      DO 61 I = 1, NM1
         A(I) = 2.0D0*I - 1.0D0 + ALPHA
   61    B(I) = DSQRT(I*(I + ALPHA))
      A(N) = 2.0D0*N - 1 + ALPHA
      RETURN
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
c 
c 
c 
C 
C 
      SUBROUTINE IMTQL2(N, D, E, Z, IERR)
C 
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE IMTQL2,
C     NUM. MATH. 12, 377-383(1968) BY MARTIN AND WILKINSON,
C     AS MODIFIED IN NUM. MATH. 15, 450(1970) BY DUBRULLE.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).
C     THIS IS A MODIFIED VERSION OF THE 'EISPACK' ROUTINE IMTQL2.
C 
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND FIRST COMPONENTS OF THE
C     EIGENVECTORS OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE IMPLICIT QL
C     METHOD.
C 
C     ON INPUT:
C 
C        N IS THE ORDER OF THE MATRIX;
C 
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX;
C 
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS FIRST N-1 POSITIONS.  E(N) IS ARBITRARY;
C 
C        Z CONTAINS THE FIRST ROW OF THE IDENTITY MATRIX.
C 
C      ON OUTPUT:
C 
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
C          UNORDERED FOR INDICES 1, 2, ..., IERR-1;
C 
C        E HAS BEEN DESTROYED;
C 
C        Z CONTAINS THE FIRST COMPONENTS OF THE ORTHONORMAL EIGENVECTORS
C          OF THE SYMMETRIC TRIDIAGONAL MATRIX.  IF AN ERROR EXIT IS
C          MADE, Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
C          EIGENVALUES;
C 
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C 
C     ------------------------------------------------------------------
C 
        save
      INTEGER I, J, K, L, M, N, II, MML, IERR
      REAL*8 D(N), E(N), Z(N), B, C, F, G, P, R, S, MACHEP
      REAL*8 DSQRT, DABS, DSIGN
C 
C     :::::::::: MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C                MACHEP = 16.0D0**(-13) FOR LONG FORM ARITHMETIC
C                ON S360 ::::::::::
CCCC  DATA MACHEP/Z3410000000000000/
      DATA MACHEP/1.0D-14/
C 
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C 
      E(N) = 0.0D0
      DO 240 L = 1, N
         J = 0
C     :::::::::: LOOK FOR SMALL SUB-DIAGONAL ELEMENT ::::::::::
  105    DO 110 M = L, N
            IF (M .EQ. N) GO TO 120
            IF (DABS(E(M)) .LE. MACHEP * (DABS(D(M)) + DABS(D(M+1))))
     X         GO TO 120
  110    CONTINUE
C 
  120    P = D(L)
         IF (M .EQ. L) GO TO 240
         IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     :::::::::: FORM SHIFT ::::::::::
         G = (D(L+1) - P) / (2.0D0 * E(L))
         R = DSQRT(G*G+1.0D0)
         G = D(M) - P + E(L) / (G + DSIGN(R, G))
         S = 1.0D0
         C = 1.0D0
         P = 0.0D0
         MML = M - L
C 
C     :::::::::: FOR I=M-1 STEP -1 UNTIL L DO -- ::::::::::
         DO 200 II = 1, MML
            I = M - II
            F = S * E(I)
            B = C * E(I)
            IF (DABS(F) .LT. DABS(G)) GO TO 150
            C = G / F
            R = DSQRT(C*C+1.0D0)
            E(I+1) = F * R
            S = 1.0D0 / R
            C = C * S
            GO TO 160
  150       S = F / G
            R = DSQRT(S*S+1.0D0)
            E(I+1) = G * R
            C = 1.0D0 / R
            S = S * C
  160       G = D(I+1) - P
            R = (D(I) - G) * S + 2.0D0 * C * B
            P = S * R
            D(I+1) = G + P
            G = C * R - B
C     :::::::::: FORM FIRST COMPONENT OF VECTOR ::::::::::
            F = Z(I+1)
            Z(I+1) = S * Z(I) + C * F
  200       Z(I) = C * Z(I) - S * F
C 
         D(L) = D(L) - P
         E(L) = G
         E(M) = 0.0D0
         GO TO 105
  240 CONTINUE
C 
C     :::::::::: ORDER EIGENVALUES AND EIGENVECTORS ::::::::::
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
C 
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
C 
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
         P = Z(I)
         Z(I) = Z(K)
         Z(K) = P
  300 CONTINUE
C 
      GO TO 1001
C     :::::::::: SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ::::::::::
 1000 IERR = L
 1001 RETURN
C     :::::::::: LAST CARD OF IMTQL2 ::::::::::
      END
