        IMPLICIT REAL *8 (A-H,O-Z)
        DIMENSION X0(1000),Y0(1000),
     1   w(200 000),ts(10 000),xs(10 000),ys(10 000),
     2      curv(10 000),tangs(2,10 000),zs(2,10 000),
     3      tangs2(2,10 000),diffs(2,10 000),ders2(2,10 000),
     4      diffc(10 000),curv2(10 000)
C 
        CALL PRINI(6,13)
c 
C      INPUT ALL PARAMETERS
C 
         PRINT *, 'ENTER N2'
         READ *,N2
         CALL PRINF('N2=*',N2,1 )
C 
         PRINT *, 'ENTER NMIN'
         READ *,NMIN
         CALL PRINF('NMIN=*',NMIN,1 )
  
  
        nsub=1000
C 
C 
C       CREATE THE TEST CURVE
C 
        CALL CREELI(X0,Y0,N0)
        CALL CRSHIP(X0,Y0,N0)
cccc       CALL CRETRI(X0,Y0,N0)
ccc        CALL STAR (X0,Y0,N0)
cccc        CALL  CRESENT (X0, Y0, N0)
ccc        CALL  STAR2 (X0, Y0, N0)
  
cccc        call cavitget(x0,y0,n0)
  
  
        CALL PRINF('AND N0=*',N0,1)
  
        iw=21
        call lotaplot(iw,x0,y0,n0,'curve as entered*')
  
  
        nlarge=n2
c 
c        resample the user-specified curve in an equispaced manner
c 
        lenw=100 000
        call rsrespni(ier,x0,y0,n0,nmin,nlarge,
     1      curvelen,err,w,lenw,lsave,lused)
  
        call prinf('after rsrespni, ier=*',ier,1)
        call prin2('after rsrespni, err=*',err,1)
  
c 
c       construct equispaced nodes on the interval [0,curvelen]
c 
        h=curvelen/nsub
        do 1200 i=1,nsub
c 
        ts(i)=(i-1)*h
 1200 continue
  
  
        call prin2('after rsrespni, curvelen=*',curvelen,1)
        call prin2('ts as created*',ts,nsub)
c 
c        create the equispaced discretization of the curve
c 
        do 1400 i=1,nsub
c 
        call rsrespnt(ts(i),zs(1,i),tangs(1,i),curv(i),
     1      w,par2)
 1400 continue
c 
c       plot the resulting curve
c 
        call arrsepa(zs,xs,ys,nsub)
c 
        iw=22
        call lotaplot(iw,xs,ys,nsub,'curve as resampled*')
c 
c       test the tangent vectors produced by the
c       subroutine vs. the numerically obtained ones
c 
        do 2200 i=2,nsub-1
c 
        tangs2(1,i)= (zs(1,i+1)-zs(1,i-1))/h/2
        tangs2(2,i)= (zs(2,i+1)-zs(2,i-1))/h/2
c 
 2200 continue
c 
        tangs2(1,1)= (zs(1,2)-zs(1,nsub))/h/2
        tangs2(2,1)= (zs(2,2)-zs(2,nsub))/h/2
c 
        tangs2(1,nsub)= (zs(1,1)-zs(1,nsub-1))/h/2
        tangs2(2,nsub)= (zs(2,1)-zs(2,nsub-1))/h/2
  
  
  
        call prin2('tangs from rsrespnt are*',tangs,nsub*2)
        call prin2('tangs2 numerically are*',tangs2,nsub*2)
  
        do 2400 i=1,nsub
c 
        diffs(1,i)=tangs2(1,i)-tangs(1,i)
        diffs(2,i)=tangs2(2,i)-tangs(2,i)
c 
 2400 continue
c 
        call prin2('and differences are*',diffs,nsub*2)
c 
c        test the calculation of the curvature
c 
        do 2600 i=1,nsub
c 
c 
        ders2(1,i)= (tangs(1,i+1)-tangs(1,i-1))/h/2
        ders2(2,i)= (tangs(2,i+1)-tangs(2,i-1))/h/2
c 
 2600 continue
c 
        ders2(1,1)= (tangs(1,2)-tangs(1,nsub))/h/2
        ders2(2,1)= (tangs(2,2)-tangs(2,nsub))/h/2
c 
        ders2(1,nsub)= (tangs(1,1)-tangs(1,nsub-1))/h/2
        ders2(2,nsub)= (tangs(2,1)-tangs(2,nsub-1))/h/2
  
        do 2800 i=1,nsub
c 
        curv2(i)=ders2(2,i)*tangs(1,i)-ders2(1,i)*tangs(2,i)
c 
        diffc(i)=curv2(i)-curv(i)
 2800 continue
c 
        call prin2('and curv=*',curv,nsub)
        call prin2('and differences are*',diffc,nsub)
  
  
  
  
  
c 
        STOP
        END
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c          Here starts debugging and testing junk
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c 
c 
        subroutine arrsepa(z,x,y,n)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1),z(2,1)
c 
c 
        do 1200 i=1,n
        x(i)=z(1,i)
        y(i)=z(2,i)
 1200 continue
        return
        end
c 
C 
C 
C 
C 
        SUBROUTINE CRETRI(X,Y,N)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        DIMENSION X(1),Y(1),IX(4),IY(4),REA(2)
        COMPLEX *16 SHIP(14),COM
        EQUIVALENCE (REA(1),COM)
CCCC    DATA IX/0,5,4,0/,IY/0,1,3,0/
CCCC    DATA IX/0,1,1,0/,IY/0,0,1,0/
        DATA IX/0,1,1,0/,IY/0,0,1,0/
        DATA SHIP /(18.0,0.0), (28.0,0.0), (32.0,4.0),
     9    (21.0,4.0),
     1    (21.0,5.5),(29.0,6.5), (29.0,7.0), (21.0,6.0),
     2    (21.0,7.0),(15.0,7.0), (15.0,4.0), ( 4.0,4.0),
     3    ( 8.0,0.0),(18.0,0.0)/
        DO 1200 I=1,4
        X(I)=IX(I)
        Y(I)=IY(I)
 1200 CONTINUE
        X(3)=X(3)+0.00001
        N=4
        RETURN
C 
C 
C 
C 
        ENTRY CREELI(X,Y,N)
        A=1
        B=4
        B=1
        N=200
        DONE=1
        PI=DATAN(DONE)*4
        H=2*PI/(N-1)
        DO 2200 I=1,N
        D=(I-1)*H
        X(I)=A*DCOS(D)
        Y(I)=B*DSIN(D)
 2200 CONTINUE
        RETURN
C 
C 
C 
C 
        ENTRY CRSHIP(X,Y,N)
        N=14
        DO 3200 I=1,N
        COM=SHIP(I)
        X(I)=REA(1)
        Y(I)=REA(2)
 3200 CONTINUE
CCCC     IF(2.NE.3) RETURN
        DO 3400 I=1,N
           J=I/2
        J=I-J*2
CCCCC   CALL PRINF('I=*',I,1)
CCCCC   IF( (J.EQ.1) .AND. (I.NE.1) .AND. (I.NE.N) ) X(I)=X(I)+0.001D-5
CCCC    IF(J.NE.0) X(I)=X(I)+0.0000001D0
 3400 CONTINUE
CCC      X(1)=X(1)+0.000001
CCCC     X(1)=X(1)+0.000001D0
CCC       CALL PRIN2('IN CRSHIP, X=*',X,N)
        RETURN
        END
C 
C 
C 
C 
C 
      SUBROUTINE  CRESENT (XX, YY, N)
        save
      REAL*8      XX(1), YY(1)
      REAL*8      R, X0, Y0, THETA, T, DT
      REAL*8      PI
      INTEGER     N, M
      N=101
      PI = DATAN(1.0D0)*4.0D0
C 
C---- N MUST BE ODD (INPUT)
C 
      M = (N-1)/2
      IF (2*M+1 .NE. N) THEN
         PRINT *, ' ERROR. N IS NOT ODD.'
         STOP
      ENDIF
  
      R = 1.0D0
      X0 = 1.0D0
      Y0 = 0.0D0
      THETA = 60.0D0 * PI / 180.0D0
  
      DT = 2.0D0 * (PI - THETA) / FLOAT(M)
  
      INDEX = 1
      DO 10 I = 0, M
         T = THETA + I * DT
         XX(INDEX) = X0 + R * COS(T)
         YY(INDEX) = Y0 + R * SIN(T)
         INDEX = INDEX+1
 10   CONTINUE
  
      R = DSQRT(13.0D0) / 4.0D0
      X0 = 1.250D0
      Y0 = 0.0D0
      THETA = DATAN(2.0D0 * DSQRT(3.0D0))
  
      DT = 2.0D0 * (PI - THETA) / FLOAT(M+1)
  
      DO 20 I = M-1, 0, -1
         T = THETA + I * DT
         XX(INDEX) = X0 + R * COS(T)
         YY(INDEX) = Y0 + R * SIN(T)
         INDEX = INDEX+1
 20   CONTINUE
  
      RETURN
      END
  
  
  
      SUBROUTINE  STAR2(XX, YY, N)
        save
      REAL*8      XX(1), YY(1)
      REAL*8      R, RBIG, X0, Y0, T
      REAL*8      PI
      INTEGER     N
  
      PI = DATAN(1.0D0)*4.0D0
  
C 
C---- N IS OUTPUT.
C 
      N = 11
  
      RBIG =  0.3D0
      X0   =  2.0D0
      Y0   =  0.0D0
  
      R  = RBIG * DSIN(18.0D0 * PI /180.0D0) / SIN(126. * PI /180.0D0)
  
      T = 90. * PI / 180.
  
      XX(1) = X0
      YY(1) = Y0 + RBIG
  
      T = T + 36. * PI / 180.
      XX(2) = X0 + R * COS(T)
      YY(2) = Y0 + R * SIN(T)
  
      T = T + 36.  * PI / 180.
      XX(3) = X0 + RBIG* COS(T)
      YY(3) = Y0 + RBIG* SIN(T)
  
      T = T + 36.  * PI / 180.
      XX(4) = X0 + R * COS(T)
      YY(4) = Y0 + R * SIN(T)
  
      T = T + 36.  * PI / 180.
      XX(5) = X0 + RBIG* COS(T)
      YY(5) = Y0 + RBIG* SIN(T)
  
      T = T + 36.  * PI / 180.
      XX(6) = X0 + R * COS(T)
      YY(6) = Y0 + R * SIN(T)
  
      T = T + 36.  * PI / 180.
      XX(7) = X0 + RBIG* COS(T)
      YY(7) = Y0 + RBIG* SIN(T)
  
      T = T + 36.  * PI / 180.
      XX(8) = X0 + R * COS(T)
      YY(8) = Y0 + R * SIN(T)
  
      T = T + 36.  * PI / 180.
      XX(9) = X0 + RBIG* COS(T)
      YY(9) = Y0 + RBIG* SIN(T)
  
      T = T + 36.  * PI / 180.
      XX(10) = X0 + R * COS(T)
      YY(10) = Y0 + R * SIN(T)
  
C     LAST POINT = FIRST.
  
      XX(11) = X0
      YY(11) = Y0 + RBIG
  
      RETURN
      END
  
  
      SUBROUTINE STAR (XX,YY,N)
        save
      REAL*8      XX(1), YY(1)
      INTEGER     N
  
      REAL*8      PI, R, T, X0, Y0, RBIG
  
      PI = 4.0D0 * DATAN(1.0D0)
  
      RBIG = 10.0D0
      X0   =  0.0D0
      Y0   =  0.0D0
  
      R  = RBIG * SIN(18. * PI /180.) / SIN(126. * PI /180.)
  
      T = 90. * PI / 180.
  
      N = 11
  
      XX(1) = X0
      YY(1) = Y0 + RBIG
  
      T = T + 36. * PI / 180.
      XX(2) = X0 + R * COS(T)
      YY(2) = Y0 + R * SIN(T)
  
      T = T + 36.  * PI / 180.
      XX(3) = X0 + RBIG* COS(T)
      YY(3) = Y0 + RBIG* SIN(T)
  
      T = T + 36.  * PI / 180.
      XX(4) = X0 + R * COS(T)
      YY(4) = Y0 + R * SIN(T)
  
      T = T + 36.  * PI / 180.
      XX(5) = X0 + RBIG* COS(T)
      YY(5) = Y0 + RBIG* SIN(T)
  
      T = T + 36.  * PI / 180.
      XX(6) = X0 + R * COS(T)
      YY(6) = Y0 + R * SIN(T)
  
      T = T + 36.  * PI / 180.
      XX(7) = X0 + RBIG* COS(T)
      YY(7) = Y0 + RBIG* SIN(T)
  
      T = T + 36.  * PI / 180.
      XX(8) = X0 + R * COS(T)
      YY(8) = Y0 + R * SIN(T)
  
      T = T + 36.  * PI / 180.
      XX(9) = X0 + RBIG* COS(T)
      YY(9) = Y0 + RBIG* SIN(T)
  
      T = T + 36.  * PI / 180.
      XX(10) = X0 + R * COS(T)
      YY(10) = Y0 + R * SIN(T)
  
C     LAST POINT = FIRST.
  
      XX(11) = X0
      YY(11) = Y0 + RBIG
  
      RETURN
      END
c 
c 
c 
c 
c 
        subroutine indint(z,n,h,zint,zint0)
        implicit real *8 (a-h,o-z)
        save
        complex *16 z(1),zint(1),zint0
c 
        zint(1)=zint0
ccc        zint(1)=(z(1)+z(2))/2
        do 1200 i=1,n-1
        zint(i+1)=zint(i)+(z(i)+z(i+1))*h/2
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine curtes(der1,der2,w,n)
        implicit real *8 (a-h,o-z)
        save
        dimension der1(2,1),der2(2,1),w(1)
c 
c 
        do 1200 i=1,n
        w(i)=der1(1,i)*der2(2,i)-der1(2,i)*der2(1,i)
 1200 continue
        call prin2('and the curvature is*',w,n)
        return
        end
c 
c 
c 
c 
c 
        subroutine dertes(z,der,w,h,n,der2,der3)
        implicit real *8 (a-h)
        save
        complex *16 z(1),der(1),w(1),der2(1),der3(1)
c 
c       check the first derivative
c 
        do 1200 i=2,n-1
        w(i)=(z(i+1)-z(i-1))/(2*h)
 1200 continue
        call prin2('in dertes, der=*',der,n*2)
        call prin2('in dertes, w=*',w,n*2)
c 
        do 1400 i=1,n
        w(i)=w(i)-der(i)
 1400 continue
        call prin2('and the differences are*',w,n*2)
c 
c        check the second derivative
c 
        do 1600 i=2,n-1
        w(i)=(der(i+1)-der(i-1))/(2*h)
 1600 continue
c 
        call prin2('numerical second derivative is*',
     1   w,n*2)
        call prin2('analytical second derivative is*',
     1   der2,n*2)
c 
        do 1800 i=1,n
        w(i)=w(i)-der2(i)
 1800 continue
        call prin2('and differences are*',w,n*2)
c 
c       check the third derivative
c 
        do 2600 i=2,n-1
        w(i)=(der2(i+1)-der2(i-1))/(2*h)
 2600 continue
c 
        call prin2('numerical third derivative is*',
     1   w,n*2)
        call prin2('analytical third derivative is*',
     1   der3,n*2)
c 
        do 2800 i=1,n
        w(i)=w(i)-der3(i)
 2800 continue
        call prin2('and differences are*',w,n*2)
        return
        end
  
  
  
  
c 
c 
c 
c 
c 
        subroutine cavitget(x0,y0,n0)
        implicit real *8 (a-h,o-z)
        save
        real *8 x0(1),y0(1)
c 
        toclose=-1.9
        toclose=2.3
        toclose=2
        toclose=1.5
  
cccc        toclose=-2.1
  
  
c 
        stretch=2
cccc        stretch=1
  
        x0(1)=0
        y0(1)=0
c 
        x0(2)=7
        y0(2)=0
c 
        x0(3)=7
        y0(3)=12
c 
        x0(4)=3
        y0(4)=12
c 
        x0(5)=3
        y0(5)=10
c 
        x0(6)=5
        y0(6)=10
c 
        x0(7)=5
        y0(7)=2
c 
c       doing the left side
c 
        x0(8)=-5
        y0(8)=2
c 
        x0(9)=-5
        y0(9)=10
c 
        x0(10)=-3
        y0(10)=10
c 
        x0(11)=-3
        y0(11)=12
c 
        x0(12)=-7
        y0(12)=12
c 
        x0(13)=-7
        y0(13)=0
c 
        x0(14)=0
        y0(14)=0
c 
c        more or less, enclose the cavity
c 
        x0(4)=x0(4)-toclose
        x0(5)=x0(5)-toclose
c 
        x0(10)=x0(10)+toclose
        x0(11)=x0(11)+toclose
  
        n0=14
  
c 
c        stretch the cavity vertically
c 
        do 2200 i=1,n0
        y0(i)=y0(i)*stretch
 2200 continue
  
        return
        end
C 
C 
C 
C 
C 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c      THIS IS THE END OF THE TESTING & DEBUGGING CODE AND
C      THE START OF THE ACTUAL RESAMPLING ROUTINES
C 
C 
C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c 
c 
c 
c 
c 
        subroutine rsrespni(ier,x0,y0,n0,nmin,nlarge7,
     1      curvelen,err,w,lenw,lsave,lused)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension z(2),tang(2),w(1),x0(1),y0(1),
     1      z2(2),tang2(2),der22(2)
c 
c 
c       This entry constructs an equispaced discretization of the
c       user-supplied curve, to be used to obtain arbitrary points
c       on that curve, specified by their arc-length distance from
c       the origin. The user provides the curve in the form of the
c       nodes (x0,y0) on the curve; the subroutine smoothes and
c       resamples the curve in an equispaced manner; in reality, it
c       feeds  the subroutine rsresa (see elsewhere) that
c       resamples the curve. Please note that the actual evaluations
c       are performed by the entry rsrespnt of this subroutine (see
c       below). This entry has no use as a stand-alone device.
c 
c                   Input parameters:
c 
c  x0,y0 - user-supplied curve to be resampled
c  n0 - the number of nodes in (x0,y0). note that
c        the input curve must be closed, i.e.
c        x0(n0)=x0(1)
c        y0(n0)=y0(1)
c        otherwise, ier is set to 16 and the execution of the
c        subroutine terminated.
C   NMIN - THE HIGHEST ORDER OF A FOURIER NODE OF THE
C        CURVATURE OF THE USER'S CURVE LEFT UNCHANGED BY THE
C        FILTER.
c  nlarge7  - the number of nodes into which the curve is to be resampled;
c        must be large
c 
c  lenw - the length (in real *8 words) of the array w provided by
c        the user
C 
C                          OUTPUT PARAMETERS:
C 
C   IER -ERROR RETURN CODE.
C        IER=0 MEANS SUCCCESSFUL CONCLUSION.
C        IER=8 means that the newton process closing
c              the curve failed to converge to the prescribed
c              precision. more specifically, the discrepancy
c              failed to go down after an iteration. this is
c              not a fatal error, but a careful examination
c              of the vector errs (see below) is recommended.
C        IER=12 means that the newton process closing
c              the curve failed to converge to the prescribed
c              precision. more specifically, the discrepancy
c              failed to decay below eps after 30 iterations,
c              while it kept going down on every step. this is
c              not a fatal error, but a careful examination
c              of the vector errs (see below) is recommended.
c              furhthermore, this is an extremely suspicious
c              situation that has been never encountered in
c              our experiments.
C        IER=16 MEANS THAT THE USER-SUPPLIED CURVE (X0,Y0)
C              IS NOT A CLOSED ONE (I.E. ITS FIRST AND LAST
C              POINTS DO NOT COINCIDE). THIS IS A FATAL ERROR.
C        IER=32 MEANS THAT n is less than or equal to nmin*4.
c              this is a fatal error
C        IER=64 MEANS THAT n is less than or equal to n0.
c              this is a fatal error
C        IER=96 MEANS THAT both of the preceeding conditions
c               have occured.
c              this is a fatal error
c        ier=16000 means that the length of the user-provided array
c              is insufficient. this is a fatal error
c  curvelen - the length of the filtered and resampled curve
c  err - the accuracy to which (in the opinion of the subroutine) the
c        curve could be resampled by nlarge7 nodes with the parameter
c        nmin as supplied by the user.
c  w - the array to be used by the entry rsrespnt (see below);
c         the first lsave elements of w should not be altered between
c         the call to this entry and the subsequent calls to rsrespnt
c  lsave - the number of elements of w that should not be changed
c         between the call to this entry and the subsequent calls
c         to rsrespnt
c 
c       . . . allocate memory for the resampling of the user-supplied
c             curve at an ungodly number of points
c 
        nlarge=nlarge7
c 
        its=1
        lts=nlarge+4
c 
        izs=its+lts
        lzs=2*nlarge+10
c 
        ider1=izs+lzs
        lder1=2*nlarge+10
c 
        ider2=ider1+lder1
        lder2=2*nlarge+10
c 
        ifis=ider2+lder2
        lfis=nlarge+4
c 
        iw=ifis+lfis
c 
        lsave=ifis-1
c 
        ltot=8*nlarge+n0+100
        lused=ltot
        if(ltot .le. lenw) goto 1100
c 
        ier=16000
        return
 1100 continue
c 
c       resample the user-specified curve at an ungodly
c       number of points
c 
        eps=1.0d-14
        nders=2
c 
        call rsresa(ier,x0,y0,n0,nmin,nders,
     1     w(izs),w(ider1),w(ider2),der3,nlarge,w(ifis),h,
     2     acc,err,w(iw),ltot)
c 
cccc        call prin2('after rsresa, ac=*',acc,1)
cccc        call prin2('after rsresa, err=*',err,1)
c 
        err=err+acc
c 
        do 1200 i=1,nlarge+1
        w(i)=(i-1)*h
c 
 1200 continue
c 
        ifinit=1
        n=20
        xinit=(w(1)+w(2))/2
c 
        call rsresper(w(its),w(izs),w(ider1),w(ider2),
     1      nlarge,xinit,n,z2,tang2,der22,ifinit)
        curvelen=nlarge*h
  
        return
c 
c 
c 
c 
        entry rsrespnt(t,z,tang,curv,w,par2)
c 
c       this entry finds the location of a point on the curve,
c       given the user-supplied distance of this point from the
c       beginning of the curve, along the arc-length. In order
c       for this entry to work, the user must have called the
c       entry analypni (see above).
c 
c                    Input parameters:
c 
c  t - the distance of the point from the beginning of the curve,
c       along the arc-length
c 
c  w - array produced by a prior call to the entry rsrespni (see above)
c 
c                     Output parameters:
c 
c  z - the location of the point in R^2
c  tang - the tangent vectorat the point z (not normalized)
c  curv - the curvature of the curve at the point t. Please note that
c          curv can be either positive or negative; it is positive when
c          the curve is convex
c 
c                     Unused parameters
c 
c  par1, par2 - unused
c 
c        use interpolation to find the value of the original parameter
c        corresponding to the arclength value t
c 
        ifinit=0
        call rsresper(w(its),w(izs),w(ider1),w(ider2),
     1      nlarge,t,n,z,tang,der22,ifinit)
  
c 
        curv=der22(2)*tang(1)-der22(1)*tang(2)
  
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine rsresper(ts,zs,tangs,der2,m,x,n,zout,tangout,
     1      der2out,ifinit)
        implicit real *8 (a-h,o-z)
        save
        real *8 ts(1),coefs(100),coefsx(100),zs(2,1),tangs(2,1),
     1      tsnew(100),tangsnew(2,100),zout(2),tangout(2),
     2      zsnew(2,100),der2(2,1),der2out(2),der2new(2,100)
c 
c       this subroutine interpolates (via standard Lagrange
c       interpolation) a user-specified perioudic function
c       from the user-provided equispaced grid to the user-supplied
c       point x on the line. The point is permitted to be outside the
c       interval where the discretization is located, but by no more
c       than the length of the interval on which the function is
c       originally defined (i.e. by no more than ts(m)-ts(1)). Please
c       note that the subpoutine assumes that ts(m) .neq. ts(1),
c       just like standard FFT routines do; conceptually,
c       ts(m+1)=ts(1).
c 
c 
c                    Input parameters:
c 
c  ts - the equispaced nodes at which the user-supplied function is
c       tabulated
c  fs - the values of the interpolant at the nodes ts
c  m - the number of nodes in the array ts
c  x - the point where the function is to be interpolated
c  ifinit - an integer parameter telling the subroutine whether the
c       Lagrange interpolation routine has to be initialized.
c          ifinit=1 will cause the initialization to be performed
c          ifinit=0 will cause the initialization to be skipped.
c       Please note that ifinit has to be set to one each time the
c       sampling distance in the array ts changes
c 
c                     Output parameters:
c 
c  fout - the interpolated value of the function at the point x
c 
c 
c 
c       if the user so requested - initialize the interpolation routine
c 
        if(ifinit .eq. 1) call rslagrin(ts,n,coefs)
c 
c       use interpolation to find the value of the function
c       for the user-specified x
c 
c       . . . find on which subinterval x lives
c 
        ier=0
c 
        h=ts(2)-ts(1)
        iint=(x-ts(1))/h
        istart=iint-n/2+2
        iend=istart+n-1
c 
       do 2300 i=1,n
c 
       tsnew(i)=ts(1)+(istart+i-2)*h
cccc       fsnew(i)=fs(irsrwrap(m,istart+i-1))
  
       zsnew(1,i)=zs(1,irsrwrap(m,istart+i-1))
       zsnew(2,i)=zs(2,irsrwrap(m,istart+i-1))
c 
       tangsnew(1,i)=tangs(1,irsrwrap(m,istart+i-1))
       tangsnew(2,i)=tangs(2,irsrwrap(m,istart+i-1))
c 
       der2new(1,i)=der2(1,irsrwrap(m,istart+i-1))
       der2new(2,i)=der2(2,irsrwrap(m,istart+i-1))
c 
 2300 continue
c 
cccc        call prin2('in rsresper, x=*',x,1)
c 
       call rslagrco(tsnew,coefs,n,x,coefsx)
c 
cccc        call prin2('in rsresper, coefsx=*',coefsx,n)
cccc        call prin2('in rsresper, coefs=*',coefs,n)
c 
        zout(1)=0
        zout(2)=0
        tangout(1)=0
        tangout(2)=0
  
        der2out(1)=0
        der2out(2)=0
c 
        do 2400 i=1,n
        zout(1)=zout(1)+coefsx(i)*zsnew(1,i)
        zout(2)=zout(2)+coefsx(i)*zsnew(2,i)
c 
        tangout(1)=tangout(1)+coefsx(i)*tangsnew(1,i)
        tangout(2)=tangout(2)+coefsx(i)*tangsnew(2,i)
c 
        der2out(1)=der2out(1)+coefsx(i)*der2new(1,i)
        der2out(2)=der2out(2)+coefsx(i)*der2new(2,i)
c 
 2400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        function irsrwrap(m,i)
c 
        save
        if( (i .le. m) .and. (i .ge. 1)) then
           irsrwrap=i
           return
        endif
c 
        if(i .lt. 1) irsrwrap=m+i
        if(i .gt. m) irsrwrap=i-m
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine rslagrin(t,n,coefs)
        implicit real *8 (a-h,o-z)
        save
        real *8 t(1),coefs(1),coefsx(1),cd,
     1      coesum,snear,x
c 
c       this subroutine performs initialization for the entry
c       rslagrco (see below). its output is the array coefs
c       to be used by rslagrco to construct interpolation
c       coefficients for the interpolation on the real line.
c 
c                         input parameters:
c 
c  t - the nodes on which the interpolation is based.
c  n - the number of elements in array t
c 
c                         output parameters:
c 
c  coefs -  an array of length n to be used by the entry lagreva
c 
c        . . . initialize the interpolation routine based on
c              user-specified nodes t(i)
c 
        done=1
        do 1400 i=1,n
        cd=1
        do 1200 j=1,n
        if(j.eq.i) goto 1200
        cd=cd*(t(i)-t(j))
 1200 continue
        coefs(i)=done/cd
 1400 continue
        return
c 
c 
c 
c 
       entry rslagrco(t,coefs,n,x,coefsx)
c 
c       this entry constructs interpolation coiefficients
c       connecting the values of a function f tabulated at the
c       nodes t(i) on the real line  with the value of f at
c       the real polint x
c 
c                       input parameters:
c 
c  t - the nodes on which the interpolation is based.
c  coefs -  an array of length n produced the entry rslagrin (see)
c  n - the number of elements in arrays t, coefs
c  x - the point in the complex plane where f is to be evaluated
c 
c                         output parameters:
c 
c  coefsx - the array of interpolation coefficients
c 
c       . . . find the node nearest to this point
c 
        knear=1
        cd=x-t(1)
        dd=cd**2
        dist=dd
        do 2400 i=2,n
        cd=x-t(i)
        dd=cd**2
        if(dd .gt. dist) goto 2400
        dist=dd
        knear=i
 2400 continue
c 
c       calculate the coefficient in front of the sum
c 
        coesum=1
        do 4200 i=1,n
        coesum=coesum*(x-t(i))
 4200 continue
c 
c       calculate the sum excluding the term corresponding
c       to knear
c 
        do 4400 i=1,n
        if(i .eq. knear) goto 4400
c 
        coefsx(i)=coefs(i)*coesum/(x-t(i))
 4400 continue
c 
c       account for the term knear
c 
        snear=1
        do 4600 i=1,n
        if(i .eq. knear) goto 4600
        snear=snear*(x-t(i))
 4600 continue
        coefsx(knear)=snear*coefs(knear)
c 
        return
        end
c 
c 
c 
c 
c 
  
        subroutine rsresa(ier,x0,y0,n0,nmin,nders,
     1     z,der1,der2,der3,n,fis,h,
     2     acc,err,w,ltot)
        implicit real *8 (a-h,o-z)
        save
        dimension x0(1),y0(1),z(2,1),der1(2,1),
     1     der2(2,1),der3(2,1),fis(1),w(1)
c 
c        this subroutine resamples a user-specified curve in an
c        equispaced manner with respect to the wavelength.
c        in addition, it produces several derivatives (1, 2, or 3)
c        of the curve with respect to the arc-length, tabulated
c        at the nodes of the discretization. as a matter of fact,
c        this is a memory management routine. all actual
c        processing is performed by the routine rsres2 (see).
c 
c          important note:
c 
c      the curve to be resampled must be closed and specified
c      in a counterclockwise fashion. the subroutine does not
c      check the latter condition, and if it is violated, will
c      produce garbage.
c 
c                      input parameters:
c 
c  x0,y0 - user-supplied curve to be resampled
c  n0 - the number of nodes in (x0,y0). note that
c        the input curve must be closed, i.e.
c        x0(n0)=x0(1)
c        y0(n0)=y0(1)
c        otherwise, ier is set to 16 and the execution of the
c        subroutine terminated.
C   NMIN - THE HIGHEST ORDER OF A FOURIER NODE OF THE
C        CURVATURE OF THE USER'S CURVE LEFT UNCHANGED BY THE
C        FILTER.
C   NDERS - THE NUMBER OF DERIVATIVES (WITH RESPECT TO
C        ARC LENGTH) OF THE RESAMPLED CURVE TO BE PRODUCED
C        BY THE SUBROUTINE. PERMITTED VALUES: 1, 2, 3.
C        THE VALUES OF THE PARAMETERS THAT WERE NOT REQUESTED
C        ARE NOT RETURNED, AND THESE PARAMETERS CAN BE
C        REPLACED BY DUMMIES.
c  n - the number of nodes into which the curve is to be resampled
C 
C                          OUTPUT PARAMETERS:
C 
C   IER -ERROR RETURN CODE.
C        IER=0 MEANS SUCCCESSFUL CONCLUSION.
C        IER=8 means that the newton process closing
c              the curve failed to converge to the prescribed
c              precision. more specifically, the discrepancy
c              failed to go down after in iteration. this is
c              not a fatal error, but a careful examination
c              of the vector errs (see below) is recommended.
C        IER=12 means that the newton process closing
c              the curve failed to converge to the prescribed
c              precision. more specifically, the discrepancy
c              failed to decay below eps after 30 iterations,
c              while it kept going down on every step. this is
c              not a fatal error, but a careful examination
c              of the vector errs (see below) is recommended.
c              furhthermore, this is an extremely suspicious
c              situation that has been never encountered in
c              our experiments.
C        IER=16 MEANS THAT THE USER-SUPPLIED CURVE (X0,Y0)
C              IS NOT A CLOSED ONE (I.E. ITS FIRST AND LAST
C              POINTS DO NOT COINCIDE). THIS IS A FATAL ERROR.
C        IER=32 MEANS THAT n is less than or equal to nmin*4.
c              this is a fatal error
C        IER=64 MEANS THAT n is less than or equal to n0.
c              this is a fatal error
C        IER=96 MEANS THAT both of the preceeding conditions
c               have occured.
c              this is a fatal error
c 
C   Z - RESAMPLED CURVE
C   DER1 - THE FIRST DERIVATIVE OF THE RESAMPLED CURVE
C   DER2 - THE SECOND DERIVATIVE OF THE RESAMPLED CURVE
C   DER3 - THE THIRD DERIVATIVE OF THE RESAMPLED CURVE
c   fis -tangent angles for the smoothed curve at the
c           resampled nodes
C   h - THE sampling interval on THE RESAMPLED CURVE
C   acc -  THE precision with which the curve has been resampled,
c           in the sense that acc is the estimated value of the
c           fourier series of z (as a function of arclength
c           of the curve) near n/2
c 
c   err - the accuracy with which the newton process
c           managed to close the curve
c 
c                       work arrays:
c 
c   w - must be at least 8*n+n0+100 real *8 elements long
c 
c       . . . verify that the user-supplied integer parameters
c             are compatible with each other
c 
        ier=0
        if(n .le. n0) ier=64
        if(n .le. nmin*4) ier=ier+32
        if(ier .ne. 0) return
c 
c        allocate memory for the resampling routine
c 
        is=1
        ls=n0+3
c 
        iw=is+ls
        lw=n*2+2
c 
        iwsave=iw+lw
        lwsave=n*4+30
c 
        ierrs=iwsave+lwsave
        lerrs=26
c 
        iz2=ierrs+lerrs
        lz2=2*n+2
c 
        ltot=iz2+lz2
ccc        call prinf('ltot=*',ltot,1)
c 
c       set the appropriate parameters
c 
        delta=1.0d-5
        eps=1.0d-9
c 
c       . . . resample
c 
ccc        call rsres2(ier,x0,y0,n0,s,fis,
ccc     1    n,w,wsave,nmin,z,z2,eps,niter,errs,delta,acc,
ccc     2    der1,der2,der3,nders,h)
c 
        call rsres2(ier,x0,y0,n0,w(is),fis,
     1    n,w(iw),w(iwsave),nmin,z,w(iz2),eps,niter,
     2    w(ierrs),delta,acc,der1,der2,der3,nders,h)
         err=w(ierrs+niter-1)
ccc         call prinf('after rsres2, ier=*',ier,1)
ccc         call prin2('after rsres2, acc=*',acc,1)
ccc         call prin2('after rsres2, h=*',h,1)
ccc         call prin2('after rsres2, errs=*',w(ierrs),niter)
         return
         end
c 
c 
c 
c 
c 
        subroutine rsres2(ier,x0,y0,n0,s,fis,
     1    n,w,wsave,nmin,z,z2,eps,niter,errs,delta,acc,
     2    der1,der2,der3,nders,h)
        implicit real *8 (a-h,o-z)
        save
        dimension x0(1),y0(1),s(1),fis(1),w(1),
     1     wsave(1),z(2,1),z2(2,1),errs(1),
     2     der1(2,1),der2(2,1),der3(2,1)
c 
c        this subroutine resamples a user-specified curve in an
c        equispaced manner with respect to the wavelength.
c        in addition, it produces several derivatives (1, 2, or 3)
c        of the curve with respect to the arc-length, tabulated
c        at the nodes of the discretization.
c 
c                      input parameters:
c 
c  x0,y0 - user-supplied curve to be resampled
c  n0 - the number of nodes in (x0,y0). note that
c        the input curve must be closed, i.e.
c        x0(n0)=x0(1)
c        y0(n0)=y0(1)
c        otherwise, ier is set to 16 and the execution of the
c        subroutine terminated.
c  n - the number of nodes into which the curve is to be resampled
C   NMIN - THE HIGHEST ORDER OF A FOURIER NODE OF THE
C        CURVATURE OF THE USER'S CURVE LEFT UNCHANGED BY THE
C        FILTER.
C   NDERS - THE NUMBER OF DERIVATIVES (WITH RESPECT TO
C        ARC LENGTH) OF THE RESAMPLED CURVE TO BE PRODUCED
C        BY THE SUBROUTINE. PERMITTED VALUES: 1, 2, 3.
C        THE VALUES OF THE PARAMETERS THAT WERE NOT REQUESTED
C        ARE NOT RETURNED, AND THESE PARAMETERS CAN BE
C        REPLACED BY DUMMIES.
c  eps - precision to which the newton iterations will be
c        conducted in order to close the curve
c  delta - the step used in the evaluation of derivatives
c        via the finite differences during the newton
C 
C                          OUTPUT PARAMETERS:
C 
C   IER -ERROR RETURN CODE.
C        IER=0 MEANS SUCCCESSFUL CONCLUSION.
C        IER=8 means that the newton process closing
c              the curve failed to converge to the prescribed
c              precision. more specifically, the discrepancy
c              failed to go down after in iteration. this is
c              not a fatal error, but a careful examination
c              of the vector errs (see below) is recommended.
C        IER=12 means that the newton process closing
c              the curve failed to converge to the prescribed
c              precision. more specifically, the discrepancy
c              failed to decay below eps after 30 iterations,
c              while it kept going down on every step. this is
c              not a fatal error, but a careful examination
c              of the vector errs (see below) is recommended.
c              furhthermore, this is an extremely suspicious
c              situation that has been never encountered in
c              our experiments.
C        IER=16 MEANS THAT THE USER-SUPPLIED CURVE (X0,Y0)
C              IS NOT A CLOSED ONE (I.E. ITS FIRST AND LAST
C              POINTS DO NOT COINCIDE). THIS IS A FATAL ERROR.
c   fis -tangent angles for the smoothed curve at the
c           resampled nodes
C   Z - RESAMPLED CURVE
c   niter - the number of iterations taken by the newton
c        process to close the curve
c   errs - the array of errors produced by newton on its
c        niter iterations
C   DER1 - THE FIRST DERIVATIVE OF THE RESAMPLED CURVE
C   DER2 - THE SECOND DERIVATIVE OF THE RESAMPLED CURVE
C   DER3 - THE THIRD DERIVATIVE OF THE RESAMPLED CURVE
C   acc -  THE precision with which the curve has been resampled
C   h - THE sampling interval on THE RESAMPLED CURVE
c 
c                       work arrays:
c 
c   s - must be at least n0+3 real *8 elements long
c   w - must be at least 2*n+2 real *8 elements long
c   wsave - must be at least 4*n+20 real *8 elements long
c   z2 - must be at least 4*n+20 real *8 elements long
c 
        ier=0
        numit=25
c 
c        construct the smoothed curve (not properly closed)
c 
         ifz2=0
         ifdtr=0
        call rsres1(ier,x0,y0,n0,s,fis,n,w,
ccc     1   wsave,nmin,z,z2,ifz2,ifdtr,dertrans)
     1   wsave,nmin,z)
         if(ier.ne.0) return
c 
c       conduct newton iterations to close the curve
c 
        ifjump=0
        err2=1.0d20
        do 2000 i=1,numit
ccc        call prinf('i=*',i,1)
        call oneiter(fis,n,z,wsave,nmin,
     1        delta,w,err)
        niter=i
        errs(i)=err
ccc        call prin2('err=*',err,1)
ccc        call prin2('err2=*',err2,1)
  
        if(ifjump .eq. 1) goto 2200
        if(err .lt. err2) goto 1900
        ier=8
        goto 2200
 1900 continue
         err2=err
         if(err2 .lt. eps) ifjump=1
 2000 continue
        ier=12
 2200 continue
ccc        call prinf('after the loop, ier=*',ier,1)
c 
c       from the obtained tangent angles, reconstruct the
c       curve
c 
        ifdtr=0
        if(nders .ge. 2) ifdtr=1
        call rsrec0(fis,z,n,wsave,z2,acc,ifdtr,w)
ccc        call prin2('in rsres2 after rsrec0, w=*',
ccc     1  w,n*2)
c 
        do 2400 i=1,n
        z(1,i)=z2(1,i)
        z(2,i)=z2(2,i)
 2400 continue
c 
c       construct the tangent vectors to the curve
c 
        do 2600 i=1,n
        der1(1,i)=dcos(fis(i))
        der1(2,i)=dsin(fis(i))
 2600 continue
c 
c       scale and shift the resampled curve
c 
        call rsscale(z,n,x0,y0,n0,scale)
c 
c       calculate the sampling interval on the resampled
c       curve
c 
        done=1
        pi=datan(done)*4
        h=2*pi*scale/n
ccc        call prin2('in rsres2, h=*',h,1)
c 
c        if the user so requested - generate higher order
c        derivatives
c 
        if(nders. le. 1) return
        rlen=n*h
c 
ccc        call RSDIFF(Z,Z2,N2,WSAVE,DER,M,RLEN,IDER)
        ider=1
        call RSDIFF(w,Z2,N,WSAVE,DER2,n,RLEN,IDER)
c 
        if(nders. le. 2) return
c 
        ider=2
        call RSDIFF(w,Z2,N,WSAVE,DER3,n,RLEN,IDER)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine rsscale(z,n,x0,y0,n0,scale)
        implicit real *8 (a-h,o-z)
        save
        dimension z(2,1),x0(1),y0(1)
c 
c        this subroutine scales and shifts the resampled
c        curve so that its length is the same as that for
c        the user-specified polyginal one, and that their
c        center of mass coincide.
c 
c                     input parameters
c 
c  z - the curve to be rescaled
c  n - the number of elements in z
c  x0,y0 - user-supplied curve
c  n0 - the number of elements in each of arrays x0, y0
c 
c                     output parameters:
c 
c  z - the rescaled and shifted curve
c  scale - the coefficient by which z has been scaled
c 
c       . . . determine the length of the user-supplied polygon
c             and its center of mass
c 
        uslen=0
        cxu=0
        cyu=0
        do 1200 i=1,n0-1
        d=(x0(i+1)-x0(i))**2+(y0(i+1)-y0(i))**2
        d=dsqrt(d)
        uslen=uslen+d
c 
        dx=(x0(i+1)+x0(i))/2
        dy=(y0(i+1)+y0(i))/2
c 
        cxu=cxu+dx*d
        cyu=cyu+dy*d
 1200 continue
        cxu=cxu /uslen
        cyu=cyu /uslen
ccc        call prin2('uslen=*',uslen,1)
ccc        call prin2('cxu=*',cxu,1)
ccc        call prin2('cyu=*',cyu,1)
c 
c       determine the polygonal length of the resampled curve
c 
        reslen=0
        do 1400 i=1,n
        d=(z(1,i+1)-z(1,i))**2+(z(2,i+1)-z(2,i))**2
        d=dsqrt(d)
        reslen=reslen+d
 1400 continue
ccc        call prin2('reslen=*',reslen,1)
c 
c        rescale the resampled curve
c 
        d=uslen/reslen
        scale=d
        do 1600 i=1,n+1
        z(1,i)=z(1,i)*d
        z(2,i)=z(2,i)*d
 1600 continue
c 
c       find the center of mass of the scaled resampled curve
c 
        cxr=0
        cyr=0
        do 1800 i=1,n
        cxr=cxr+z(1,i)
        cyr=cyr+z(2,i)
 1800 continue
        cxr=cxr/n
        cyr=cyr/n
c 
c       shift the scaled resampled curve
c 
        xshift=cxu-cxr
        yshift=cyu-cyr
        do 2000 i=1,n
        z(1,i)=z(1,i)+xshift
        z(2,i)=z(2,i)+yshift
 2000    continue
        return
        end
C 
C 
C 
C 
C 
        SUBROUTINE RSDIFF(Z,Z2,N2,WSAVE,DER,M,RLEN,IDER)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        COMPLEX *16 Z(1),Z2(1),IMA,DER(1),WSAVE(1)
        DATA IMA/(0.0D0,1.0D0)/
C 
C       THIS SUBROUTINE DIFFERENTIATES (IN THE FREQUENCY DOMAIN)
C       A USER-SPECIFIED CURVE AND CONVERTS IT INTO SPACE DOMAIN.
C       IT ASSUMES THAT THE CURVE IS GIVEN TO IT IN THE
C       FREQUENCY DOMAIN. THE NUMBER M OF POINTS AT WHICH THE
C       DERIVATIVE IS RETURNED DOES NOT HAVE TO BE EQUAL
C       TO THE NUMBER N OF POINTS WHERE IT HAS BEEN EVALUATED,
C       BUT M MUST DIVIDE N.
C 
C                        INPUT PARAMETERS:
C 
C   Z - THE CURVE TO BE DIFFERENTIATED (IN THE FREQUENCY DOMAIN)
C       ATTENTION!!! THIS PARAMETER IS CHANGED BY THE SUBROUTINE!
C       SEE IT IN THE LIST OF OUTPUT PARAMETERS.
C   N - THE NUMBER OF ELEMENTS IN ARRAY Z
C   WSAVE - THE ARRAY TO BE USED BY THE FFT ROUTINE DCFFTB.
C       THE USER MUST HAVE INITIALIZED IT BY CALLING THE
C       ROUTINE DCFFTI (SEE).
C   M - THE NUMBER OF ELEMENTS IN THE OUTPUT ARRAY DER.
C   RLEN - THE LENGTH OF THE CURVE BEING RESAMPLED.
C   IDER - THE ORDER OF DERIVAIVE BEING COMPUTED
C 
C                         OUTPUT PARAMETERS:
C 
C   Z - THE CURVE DIFFERENTIATED IN THE FREQUENCY DOMAIN,
C       I.E. MULTIPLIED BY I*K
C   DER - THE DERIVATIVE OF THE INPUT CURVE IN THE SPACE DOMAIN.
C 
        DONE=1
        PI=DATAN(DONE)*4
        NMULT=N2/M
        DO 1400 I=1,N2/2
        Z2(I)=Z(I)
        Z2(N2-I+1)=Z(N2-I+1)
        Z2(I)=Z(I)*IMA*(I-1)
        Z2(N2-I+1)=-Z(N2-I+1)*IMA*I
        Z(I)=Z2(I)
        Z(N2-I+1)=Z2(N2-I+1)
 1400 CONTINUE
        CALL DCFFTB(N2,Z2,WSAVE)
        D=DONE/N2/RLEN*PI*2
        D=DONE/N2 *(PI*2/RLEN)**IDER
        DO 1600 I=1,M
        J=(I-1)*NMULT+1
        DER(I)=Z2(J)*D
 1600 CONTINUE
        RETURN
        END
c 
c 
c 
c 
c 
        subroutine oneiter(thetas,n,z,wsave,nmin,
     1        delta,w,err)
        implicit real *8 (a-h,o-z)
        save
        dimension thetas(1),z(2,1),wsave(1),w(1)
c 
c       this subroutine performs one iteration of the
c       finite difference newton process for the closing of
c       a curve in two dimensions
c 
c                     input parameters:
c 
c  thetas - the tangent angles of the curve to be closed
c            note: this array will be modified by this subroutine
c  n - the number of nodes in the discretization of the curve
c        (the number of elements in the array thetas)
c  WSAVE - THE ARRAY CONTAINING PARAMETERS CREATED BY THE
C        ENTRY DCFFTI AND USED BY THE ENTRIES DCFFTF, DCFFTB.
c  nmin - the parameter of the filtering subroutine rsfilt (see)
c  delta - the small shift used to calculate the matrix of
c        derivatives in the newton via finite differences
c 
c                     output parameters:
c 
c  thetas - the tangents of the curve to be closed (hopefully,
c        closer to being closed than on entry)
c  err - the error on this iteration
c 
c                     work arrays:
c 
c  z - must be at least n*2+4 real *8 locations
c  w -  must be at least n+2 real *8 locations
c 
  
c        . . . calculate the discrepancy of the user-supplied
c              tangent angle
c 
ccc         call prinf('in oneiter, nmin=*',nmin,1)
        cosint=0
        sinint=0
        do 1200 i=1,n
        z(1,i)=dcos(thetas(i))
        z(2,i)=dsin(thetas(i))
        cosint=cosint+z(1,i)
        sinint=sinint+z(2,i)
 1200 continue
          err=dsqrt(sinint**2+cosint**2)
ccc         call prin2('in oneiter, cosint=*',cosint,1)
ccc         call prin2('in oneiter, sinint=*',sinint,1)
c 
c       obtain the matrix connecting the integrals with the
c       perturbations
c 
c       . . . filter the perturbed tangent angles
c 
        call dcfftf(n,z,wsave)
ccc         call prin2('in oneiter before rsfilt, z=*',z,n*2)
c 
ccc        call RSFILG(W,N,NMIN)
        call RSFILt(z,N,NMIN)
c 
ccc         call prin2('in oneiter, after rsfilg, z=*',z,n*2)
         d=1
         d=d/n
         call dcfftb(n,z,wsave)
         do 1400 i=1,n
         z(1,i)=z(1,i)*d
         z(2,i)=z(2,i)*d
 1400 continue
ccc         call prin2('filtered z are*',z,n*2)
c 
c        construct the perturbed tangent vectors
c        and the corresponding discrepancies
c 
        do 1600 i=1,n
        w(i)=thetas(i)+delta*z(1,i)
 1600 continue
c 
        call comdev(w,n,a11p,a12p)
c 
        do 1800 i=1,n
        w(i)=thetas(i)-delta*z(1,i)
 1800 continue
c 
        call comdev(w,n,a11m,a12m)
c 
        do 2000 i=1,n
        w(i)=thetas(i)+delta*z(2,i)
 2000 continue
c 
        call comdev(w,n,a21p,a22p)
c 
        do 2200 i=1,n
        w(i)=thetas(i)-delta*z(2,i)
 2200 continue
c 
        call comdev(w,n,a21m,a22m)
c 
c         calculate the matrix for the newton
c         via finite differences
c 
        a11=(a11p-a11m)/(2*delta)
        a12=(a12p-a12m)/(2*delta)
        a21=(a21p-a21m)/(2*delta)
        a22=(a22p-a22m)/(2*delta)
c 
        dd=a12
        a12=a21
        a21=dd
c 
ccc         call prin2('a11=*',a11,1)
ccc         call prin2('a21=*',a21,1)
ccc         call prin2('a12=*',a12,1)
ccc         call prin2('a22=*',a22,1)
c 
c        solve the linear system
c 
        y1=cosint
        y2=sinint
c 
        det=a22*a11-a12*a21
        det1=y1*a22-y2*a12
        det2=y2*a11-y1*a21
c 
        delx=det1/det
        dely=det2/det
ccc         call prin2('delx=*',delx,1)
ccc         call prin2('dely=*',dely,1)
c 
c       adjust thetas
c 
        do 2400 i=1,n
        thetas(i)=thetas(i)-(delx*z(1,i)+dely*z(2,i))
 2400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine comdev(theta,n,dx,dy)
        implicit real *8 (a-h,o-z)
        save
        dimension theta(1)
c 
c       this subroutine determines the deviation of a curve
c       specified by its tangent angles from being a closed
c       one
c 
c                        input parameters:
c 
c  thetas - the tangent angles
c  n - the number of nodes in the discretization of the boundary
c       (the number of elements in the array thetas)
c 
c                        output parameters:
c 
c  dx - distance between the ends of the ends of the
c       curve in the x direction
c  dy - distance between the ends of the ends of the
c       curve in the y direction
c 
        dx=0
        dy=0
        do 1200 i=1,n
        ddx=dcos(theta(i))
        ddy=dsin(theta(i))
        dx=dx+ddx
        dy=dy+ddy
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine rsres1(ier,x0,y0,n0,s,fis,n,w,
ccc     1   wsave,nmin,z,z2,ifz2,ifdtr,dertrans)
     1   wsave,nmin,z)
        implicit real *8 (a-h,o-z)
        save
        dimension x0(1),y0(1),s(1),fis(1),wsave(1)
cccc        complex *16 w(1),z(1),z2(1),dertrans(1)
        complex *16 w(1),z(1)
c 
c       this subroutine performs a preliminary resampling,
c       taking a user-specified curve (given by a collection
c       of points in the plane, and producing a smooth curve
c       resampled at equispaced nodes, but not a closed one.
c       the curve will have to be closed by subsequent
c       processing.
c 
C                  INPUT PARAMETERS:
C 
C   X0,Y0 - COORDINATS OF THE USER-SPECIFIED POLYGON.
C   N0 - THE NUMBER OF ELEMENTS IN ARRAYS X0,Y0
C   n - THE NUMBER OF NODES IN THE EQUISPACED RESAMPLING
C       OF THE CURVE
c   nmin - the parameter of the filtering subroutine rsfilt (see)
C 
C                  OUTPUT PARAMETERS:
C 
C   IER - ERROR RETURN CODE.
C      IER=0 MEANS SUCCESSFULL EXECUTION
C      IER=16 MEANS THAT THE USER-SUPPLIED POLYGON IS NOT
C         CLOSED. THIS IS A TERMINAL ERROR.
C   S - THE ARRAY OF CUMULATIVE POLYGONAL ARC-LENGTHS ON THE
C       USER-SUPPLIED POLYGON (NOTE THAT S(1)=0).
C   H - THE SAMPLING INTERVAL ON THE EQUISPACED-RESAMPLED
C       CURVE
c   fis - the tangent angles of the resampled curve
c   wsave - the array initialized y the subroutine dcffti,
c       that can be fed into dcfftf, dcfftb
c 
c                  work arrays:
c 
c   w - must be at least n*2+4 real *8 locations long
c   z - must be at least n*2+4 real *8 locations long. this array
c       is only used if ifz2=1 has been specified (see above)
c       otherwise, it is not used.
c 
c       . . . construct the polygonal resumpling of the curve
c 
        call REPOLY(IER,X0,Y0,n0,S,n,H,fis,w)
ccc        call prin2('in rsres1 after repoly, s=*',s,n0+1)
ccc        call prin2('in rsres1 after repoly, fis=*',fis,n+1)
        if(ier.ne.0) return
c 
c       filter the tangent angles of the curve in the frequency
c       domain using the gaussian filter
c 
        call dcffti(n,wsave)
c 
        done=1
        pi=datan(done)*4
        htet=2*pi/n
        do 1200 i=1,n
        w(i)=fis(i) - (i-1)*htet
 1200 continue
cccc         call prin2('before dcfftf, w=*',w,n*2)
c 
        call dcfftf(n,w,wsave)
ccc         call prin2('before rsfilt, w=*',w,n*2)
c 
ccc        call RSFILG(W,N,NMIN)
        call RSFILt(W,N,NMIN)
c 
ccc         call prin2('after rsfilg, w=*',w,n*2)
         d=1
         d=d/n
        call dcfftb(n,w,wsave)
         do 1300 i=1,n
         w(i)=w(i)*d
 1300     continue
ccc         call prin2('filtered fis are*',w,n+1)
c 
c        construct the tangent vectors for the filtered
c        curve and reconstruct the curve from them
c 
        do 1400 i=1,n
        fis(i)=w(i)+(i-1)*htet
 1400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine rsrec0(thetas,z,n,wsave,z2,acc,ifdtr,
     1    dertrans)
        implicit real *8 (a-h,o-z)
        save
        dimension wsave(1),thetas(1),z(2,1),z2(2,1),
     1    dertrans(1)
c 
c        this subroutine reconstructs a closed curve in R^2
c        from its tangent angles
c 
c                     input parameters:
c 
c  thetas - the tangent angles
c  n - the number of nodes in the discretization of the boundary
c       (the number of elements in the array thetas)
c  wsave - the information array to be used y the subroutines
c       dcfftf, dcfftb. this array must have been precomputed
c       by a prior call to the suroutine dcffti
c 
c                     output parameters:
c 
c  z2 - the reconstructed curve
c  acc - the accuracy of the reconstruction, in the sense
c       that acc is the estimated value of the fourier series
c       of z (as a function of arclength of the curve) near n/2
c 
c                     work array:
c 
c  z - must be 2*n+2 real *8 locations long
c 
c 
c       . . . construct the tangents to the curve from the
c             tangent angles
c 
        do 1400 i=1,n
        z(1,i)=dcos(thetas(i))
        z(2,i)=dsin(thetas(i))
 1400 continue
c 
c        integrate the tangents
c 
        call RSINTC(N,z2,z,WSAVE,acc,ifdtr,dertrans)
ccc         call prin2('after rsintc, acc=*',acc,1)
        return
        end
c 
C 
C 
C 
C 
       SUBROUTINE RSINTC(N,XOUT,W,WSAVE,acc,ifdtr,dertrans)
       IMPLICIT REAL *8 (A-H,O-Z)
        save
       COMPLEX *16 W(1),WSAVE(1),C,IMA,XOUT(1),D,H,dertrans(1)
        DATA IMA/(0.0D0,1.0D0)/
C 
C        THIS SUBROUTINE CALCULATES THE VALUES OF AN INDEFINITE
C        INTEGRAL OF A COMPLEX PERIODIC FUNCTION, GIVEN BY ITS
C        FOURIER TRANSFORM, AT A BUNCH OF EQUISPACED NODES.
C 
C 
C                      INPUT PARAMETERS:
C 
C  N - THE NUMBER OF ELEMENTS IN ARRAYS W, XOUT (SEE BELOW).
C        IT IS IMPORTANT THAT N BE A POWER OF TWO, OR AT LEAST
C        BE A PRODUCT OF SMALL PRIMES. OTHERWIZE, THE SUBROUTINE
C        CAN BE     V E R Y    SLOW.
C  W - ARRAY CONTAINING THE (COMPLEX)  FUNCTION
C        WHOSE INTEGRAL IS BEING EVALUATED. THE
C        CONTENTS OF W ARE DESTROYED BY THIS SUBROUTINE.
C  WSAVE - THE ARRAY CONTAINING PARAMETERS CREATED BY THE
C        ENTRY DCFFTI AND USED BY THE ENTRIES DCFFTF, DCFFTB.
C 
C                       OUTPUT PARAMETERS:
C 
C  XOUT - THE TABLE OF VALUES OF THE INDEFINITE INTEGRAL
c  acc - the accuracy of the integration, in the sense
c       that acc is the estimated value of the fourier series
c       of z (as a function of arclength of the curve) near n/2
c 
C 
C       . . .  FOURIER-TRANSFORM THE INPUT FUNCTION
C 
         IFIN=0
         IFTRAN=1
        IF(IFTRAN .EQ. 0) GOTO 1300
         IF(IFIN .NE.0) CALL DCFFTI(N,WSAVE)
         CALL DCFFTF(N,W,WSAVE)
c 
c        if the user so requested - return to the user
c        the fourier transform of the input vector
c 
         if(ifdtr .eq. 0) goto 1300
        do 1200 i=1,n
        dertrans(i)=w(i)
 1200 continue
c 
 1300 CONTINUE
cccc         call prin2('in rsintc, ffted w is*',w,n*2)
C 
C       INTEGRATE THE FUNCTION IN THE DUAL DOMAIN
C 
        C=W(1)
ccc          call prin2('in rsintc, c=*',c,2)
        W(1)=0
        W(N)=-W(N)/IMA
        DO 1400 I=2,N/2
        W(I)=W(I)/(IMA*(I-1))
        W(N-I+1)=-W(N-I+1)/(IMA*I)
 1400 CONTINUE
c 
c        determine the accuracy of the discretization
c 
        acc=0
        do 1500 i=-3,3
        dd=cdabs(w(n/2+i))
        if(dd .gt. acc) acc=dd
 1500 continue
ccc        call prin2('in rsintc, acc=*',acc,1)
c 
C       FOURIER-TRANSFORM THE THING BACK
C 
        CALL DCFFTB(N,W,WSAVE)
C 
C       SCALE AND ADD IN THE LINEAR TERM
C 
        DONE=1
        PI=DATAN(DONE)*4
        COEF=1
        COEF=COEF/ N
        H=2*PI*C/N**2
        DO 1600 I=1,N
        XOUT(I)=W(I)*COEF+(I-1)*H
 1600 CONTINUE
        XOUT(N+1)=2*PI*C/N + XOUT(1)
C 
C        NOW, THE PARTICULAR REPRESENTATION OF THE INDEFINITE
C        INTEGRAL IS NON-ZERO AT ZERO. SHIFT IT SO IT WILL
C        BE THAT WAY
C 
         D=XOUT(1)
        DO 1800 I=1,N+1
        XOUT(I)=XOUT(I)-D
 1800 CONTINUE
        RETURN
        END
C 
C 
C 
C 
C 
        SUBROUTINE RSFILT(W,N,NMIN)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        COMPLEX *16 ZERO,W(1)
C 
C       THIS SUBROUTINE FILTERS A FUNCTION IN THE FREQUENCY
C       DOMAIN BY APPLYING THE STANDARD COSINE FILTER. THE FILER
C       STARTS AT NMIN-TH FREQUENCY AND IS NMIN FREQUENCIES WIDE.
C 
C                        INPUT PARAMETERS:
C 
C  W - THE FOURIER TRANSFORM OF THE FUNCTION BEING FILTERED
C  N - THE NUMBER OF ELEMENTS IN W
C  NMIN - THE LOWEST FREQUENCY AT WHICH THE FILTER IS .NEQ. 1
C 
C                        OUTPUT PARAMETERS:
C 
C  W - THE FILTERED FOURIER TRANSFORM OF THE FUNCTOIN
C 
C 
C       . . . SET UP THE PARAMETERS OF THE BELL
C 
        DONE=1
        HALF=DONE/2
        PI=DATAN(DONE)*4
        H=PI/(NMIN-1)
C 
C       IMPOSE THE BELL
C 
        DO 1200 I=1,NMIN
        D=I*H
        DD=(DCOS(D)+DONE)*HALF
        W(NMIN+I)=W(NMIN+I)*DD
        W(N-NMIN-I)=W(N-NMIN-I)*DD
 1200 CONTINUE
C 
C        SET TO ZERO EVERYTHING OUTSIDE THE BELL
C 
        ZERO=0
        DO 1400 I=2*NMIN, N-2*NMIN
        W(I)=ZERO
 1400 CONTINUE
        RETURN
        END
C 
C 
C 
c 
c 
        SUBROUTINE RSFILG(W,N,NMIN)
        IMPLICIT REAL *8 (A-H,O-Z)
cccc        COMPLEX *16 ZERO,W(1)
        save
        COMPLEX *16 W(1)
  
C       THIS SUBROUTINE FILTERS A FUNCTION IN THE FREQUENCY
C       DOMAIN BY APPLYING THE STANDARD GAUSSIAN FILTER. THE FILER
C       STARTS AT NMIN-TH FREQUENCY AND IS NMIN FREQUENCIES WIDE.
C 
C                        INPUT PARAMETERS:
C 
C  W - THE FOURIER TRANSFORM OF THE FUNCTION BEING FILTERED
C  N - THE NUMBER OF ELEMENTS IN W
C  NMIN - THE LOWEST FREQUENCY AT WHICH THE FILTER IS .NEQ. 1
C 
C                        OUTPUT PARAMETERS:
C 
C  W - THE FILTERED FOURIER TRANSFORM OF THE FUNCTion
C 
C 
C       . . . SET UP THE PARAMETERS OF THE BELL
C 
        DONE=1
        D=10
        A=DLOG(D)/NMIN**2
C 
C       IMPOSE THE BELL
C 
        DD0=1
        DO 1200 I=1,N/2
        D=I
        DD=DEXP(-A*D**2)
        W(I+1)=W(I+1)*DD
        W(N-I+1)=W(N-I+1)*DD
        DD0=DD0+DD
 1200 CONTINUE
        IF(2.NE.3) RETURN
C 
        DD0=1/DD0
        DO 1400 I=1,N
        W(I)=W(I)*DD0
 1400 CONTINUE
        RETURN
        END
C 
C 
C 
C 
C 
        SUBROUTINE REPOLY(IER,X,Y,N,S,M,H,fis,thetas)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        DIMENSION X(1),Y(1),S(1),fis(1),thetas(1)
C 
C       THIS SUBROUTINE FOR A USER-SPECIFIED POLYGONAL CURVE
C       (X,Y) CONSTRUCTS THE ARRAY S OF ITS CUMULATIVE ARC-LENGTHS,
c       snd the array fis  of its tangent angles
C       THE USER-PROVIDED POLYGON MUST BE CLOSED, I.E.
C       X(N)=X(1),
C       Y(N)=Y(1).
c       note also that m must be greater than n
C 
C                  INPUT PARAMETERS:
C 
C   X,Y - COORDINATS OF THE USER-SPECIFIED POLYGON.
C   N - THE NUMBER OF ELEMENTS IN ARRAYS X,Y
C   M - THE NUMBER OF NODES IN THE EQUISPACED RESAMPLING
C       OF THE CURVE
C 
C                  OUTPUT PARAMETERS:
C 
C   IER - ERROR RETURN CODE.
C      IER=0 MEANS SUCCESSFULL EXECUTION
C      IER=16 MEANS THAT THE USER-SUPPLIED POLYGON IS NOT
C         CLOSED. THIS IS A TERMINAL ERROR.
C   S - THE ARRAY OF CUMULATIVE POLYGONAL ARC-LENGTHS ON THE
C       USER-SUPPLIED POLYGON (NOTE THAT S(1)=0).
C   H - THE SAMPLING INTERVAL ON THE EQUISPACED-RESAMPLED
C       CURVE
c   fis - the tangent angles of the rfesampled curve
c 
c                  work arrays:
c 
c   thetas - must be at least m+2 real *8 locations long
c 
        IER=0
        EPS0=1.0D-13
cccc        CALL PRIN2('IN REPOLY, X=*',X,N)
cccc        CALL PRIN2('IN REPOLY, Y=*',Y,N)
CCC     CALL PRINF('IN REPOLY, N=*',N,1)
CCC     CALL PRINF('IN REPOLY, M=*',M,1)
C 
C       CONSTRUCT THE ARRAY OF USER-PROVIDED LENGTHS FOR ALL
C       NODES OF THE USER-PROVIDED POLYGON
C 
        S(1)=0
        DO 1200 I=2,N
        D= (X(I)-X(I-1))**2 + (Y(I)-Y(I-1))**2
        S(I)=S(I-1)+DSQRT(D)
 1200 CONTINUE
cccc         call prin2('inside repoly, s=*',s,n)
        EPS=S(N)*EPS0
C 
C        IF THIS IS NOT A CLOSED POLYGON - BOMB OUT
C 
        D=DABS( X(N)-X(1) ) +  DABS( Y(N)-Y(1) )
        IF (D.LT.EPS) GOTO 1300
        IER=16
        RETURN
 1300 CONTINUE
C 
C       DETERMINE THE SAMPLING INTERVAL ON THE RESAMPLED CURVE
C 
ccc        H=S(N)/(M-1)
        H=S(N)/m
c 
c        construct the angles between the x axis and the
c        various sides of the polygon
c 
c       . . . construct the angles between the consecutive
c             sides of the polygon
c 
        call REANCR(X,Y,N,FIS)
ccc         call prin2('angles between sides of polygon are*',fis,n)
c 
c       find the angle between the first side of the polygon and
c       the x axis
c 
        wx0=x(1)-1
        wy0=y(1)
c 
        call REANGL(wX0,wY0,X(1),Y(1),X(2),Y(2),thetas(1))
c 
cccc        SUBROUTINE REANGL(X1,Y1,X2,Y2,X3,Y3,FI)
c 
c        construct the angles
c 
         do 2200 i=2,n
         thetas(i)=thetas(i-1)+fis(i)
 2200 continue
ccc         call prin2('tangent angles of original polygon are*',
ccc     1    thetas,n)
c 
c       construct the angles between the x axis and the various
c       sides of the polygon, but tabulated at the nodes of
c       the finer discretization
c 
        j=1
ccc         call prin2('in repoly, h=*',h,1)
        do 2400 i=1,m
        d=(i-1)*h
        if(d .gt.s(j+1)) j=j+1
        fis(i)=thetas(j)
 2400 continue
c 
        RETURN
        END
C 
C 
C 
C 
C 
        SUBROUTINE REANCR(X,Y,N,FIS)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        DIMENSION X(1),Y(1),FIS(1)
C 
C       FOR THE USER-SPECIFIED POLYGON, THIS SUBROUTINE
C       CONSTRUCTS THE ANGLES MADE BY THE CONSECUTIVE SIDES
C       OF THE POLYGON. THE POLYGON IS ASSUMED TO BE CLOSED.
C 
C                  INPUT PARAMETERS:
C 
C   X,Y - THE COORDINATES OF THE VERTICES OF THE USER-SUPPLIED
C         POLYGON
C   N - THE NUMBER OF ELEMENTS IN ARRAYS X,Y
C 
C                  OUTPUT PARAMETERS:
C 
C   FIS - THE ANGLES MADE BY CONSECUTIVE SIDES OF THE POLYGON
C 
C      REMARK:
C       IT IS ASSUMED THAT THE POLYGON IS CLOSED, I.E.
C       X(N)=X(1),
C       Y(N)=Y(1).
C 
C       . . . INITIALIZE THE ANGLE-COMPUTING ROUTINE AND
C             CALCULATE THE ANGLE BETWEEN THE LAST AND THE FIRST
C             SIDES. THIS ANGLE IS DECLARED TO BE AT THE
C             FIRST VERTEX.
C 
        CALL REANGI(IER)
        CALL REANGL(X(N-1),Y(N-1),X(N),Y(N),X(2),Y(2),FIS(1))
C 
C       CALCULATE THE REST OF THE ANGLES
C 
        DO 1200 I=2,N-1
        CALL REANGL(X(I-1),Y(I-1),X(I),Y(I),X(I+1),Y(I+1),FIS(I))
CCCC    CALL PRIN2('AFTER REANGL, FIS(I)  IS*',FIS(I),1)
 1200 CONTINUE
        FIS(N)=FIS(1)
        RETURN
        END
C 
C 
C 
C 
C 
        SUBROUTINE REANGL(X1,Y1,X2,Y2,X3,Y3,FI)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        DATA EPS/1.0D-30/
C 
C 
C       THIS SUBROUTINE FINDS THE ANGLE FI BETWEEN THE
C       LINE SEGMENTS  [(X1,Y1),(X2,Y2)] AND  [(X2,Y2),(X3,Y3)]
C       INITIALIZATION ENTRY POINT REANGI (SEE BELOW) HAS TO BE
C       CALLED PRIOR TO THE INVOCATION OF THIS MAIN ENTRY
C 
C        . . . FIND THE (UNSCALED) SIN AND COS OF THE ANGLE
C 
        UX=X2-X1
        UY=Y2-Y1
C 
        VX=X3-X2
        VY=Y3-Y2
C 
        SINFI=UX*VY-VX*UY
        COSFI=UX*VX+UY*VY
C 
C       IF THE ANGLE IS NOT RIGHT - FIND IT AND EXIT
C 
          IF(DABS(COSFI) .LT. EPS) GOTO 2000
        FI=DATAN(SINFI/COSFI)
CCCC    IF(COSFI .LE. 0) FI=FI+PI
        IF( (COSFI .LE. 0) .AND. (SINFI.GE.0) ) FI=FI+PI
        IF( (COSFI .LE. 0) .AND. (SINFI.LE.0) ) FI=FI-PI
        RETURN
 2000 CONTINUE
C 
C       THE ANGLE IS RIGHT. FIND IT
C 
CCCC    CALL PRIN2('IN REANGL, THE ANGLE IS RIGHT. COSFI IS*',COSFI,1)
        FI=PI2
        IF(SINFI .LT. 0.0D0) FI=-FI
CCC     IF(SINFI .GT. 0.0D0) FI=-FI
        RETURN
C 
C 
C 
C 
        ENTRY REANGI(IER)
C 
C       THIS IS INITIALIZATION ENTRY POINT. IT HAS TO BE CALLED
C       PRIOR TO ANY CALLS TO THE MAIN ENTRY REANGL (ABOVE).
c 
        DONE=1
        PI=DATAN(DONE)*4
        PI2=PI/2
        RETURN
        END
