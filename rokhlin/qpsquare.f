        implicit real *16 (a-h,o-z)
        complex *32  pol(300000),work(41000)
c 
c 
        call prini(6,13)
c 
C 
C       SET ALL PARAMETERS
C 
C 
        PRINT *, 'ENTER n, # of Gaussian nodes on each side'
        READ *,n
        CALL PRINF('n=*',n,1 )
C 
        PRINT *, 'ENTER m, highest order of orth-polynomials'
        READ *,m
        CALL PRINF('m=*',m,1 )
c 
c       construct the polynomials orthogonal on the square
c 
        call orthotes(n,m,pol,work)
        stop
        end
  
c 
c 
c 
c 
c 
        subroutine orthotes(n,m,pol,work)
        implicit real *16 (a-h,o-z)
        save
        complex *32 pol(n*4,m)
        dimension t(2000),w(2000),x(4000),y(4000),rea(2),
     1   weights(4000),rlens(4000),corners(2,4)
        complex *32 com,work(m,m),
     1  ima,zs(4000),ww(2000)
        data ima/(0.0d0,1.0q0)/
        equivalence (rea(1),com)
C 
c       construct the square
c 
c       . . . the four corners
c 
        x1=-1
        y1=-1
c 
        x2=1
        y2=-1
c 
        x3=1
        y3=1
c 
        x4=-1
        y4=1
c 
        corners(1,1)=x1
        corners(2,1)=y1
c 
        corners(1,2)=x2
        corners(2,2)=y2
c 
        corners(1,3)=x3
        corners(2,3)=y3
c 
        corners(1,4)=x4
        corners(2,4)=y4
c 
c       construct the square
c 
        call squacre(n,zs,weights,t,w,rlens,corners)
c 
c       plot the square as created
c 
        do 1400 i=1,n*4
        com=zs(i)
        x(i)=rea(1)
        y(i)=rea(2)
 1400 continue
c 
        iw=21
        call RSPLOT(X,Y,N*4,IW,'square as created*')
cccc        if(2 .ne. 3) stop
  
c 
c        construct the orthogonal polynomials on the square
c 
        do 1700 i=1,4*n
        call qpsquare(zs(i),m,ww)
        do 1600 j=1,m
        pol(i,j)=ww(j)
 1600 continue
 1700 continue
  
         call prinq('after qpsquare, pol(m)=*',pol(1,m-1),8*n)
c 
c        plot the obtained polynomials
c 
        do 2400 i=1,m
        iw1=100+i
        iw2=200+i
c 
        call polplot(pol,n*4,i,rlens,work,iw1,iw2)
         close(iw1)
         close(iw2)
 2400 continue
c 
c       check if the obtained polynomials are really orthonormal
c 
        call prods(pol,n*4,m,weights,work)
        call prinq('and the matrix of inner products*',work,m*m*2)
c 
        d=0
        done=1
        do 2800 i=1,m
        do 2600 j=1,m
        if(i .eq. j) goto 2600
        dd=work(j,i)*qconjg(work(j,i))
        if(dd .gt. d) d=dd
 2600 continue
        dd=(work(i,i)-done)*qconjg(work(i,i)-done)
        if(dd .gt. d) d=dd
 2800 continue
        call prinq('the non-orthonormality is*',qsqrt(d),1)
 3000 continue
        stop
        end
  
c 
c 
c 
c 
c 
        subroutine squacre(n,zs,weights,t,w,rlens,corners)
        implicit real *16 (a-h,o-z)
        save
        complex *32 zs(1),ima
        dimension weights(1),t(1),w(1),rlens(1),corners(2,1)
        data ima/(0.0d0,1.0d0)/
c 
c        this subroutine constructs a quadrilateral polygon
c        from four user-provided corners. it constructs
c        discretizations of the four sides at gaussian
c        points on each side, the corresponding weights
c        for the integration of points on the boundary of the
c        polygon, and the cumulative lengths along the
c        boundary of the polygon of the gaussian nodes on the
c        sides.
c 
c                     input parameters:
c 
c  n - the number of gaussian nodes on each side of the polygon
c  corners - the user-specified corners of the polygon
c 
c    explanation: corners(1,i) is the x-coordinate of the i-th
c    corner. corners(2,i) is the y-coordinate of the i-th corner.
c    the total number of them is, of couse, four.
c 
c                     output parameters:
c 
c  zs - the discretization of the polygon (zs(i) is the i-th
c       point in the complex plane. the total number of them: 4*n
c  weights - quadrature weights at the points zs
c  t - the n gaussian points on the interval [-1,1]
c  w - the n gaussian weights on the interval [-1,1]
c  rlens - the cumulative lengths along the polygon corresponding
c       to the points zs
c 
c        . . . construct the points on the square
c 
        x1=corners(1,1)
        y1=corners(2,1)
c 
        x2=corners(1,2)
        y2=corners(2,2)
c 
        x3=corners(1,3)
        y3=corners(2,3)
c 
        x4=corners(1,4)
        y4=corners(2,4)
c 
        call squacon(n,t,w,
     1    zs(1),zs(n+1),zs(2*n+1),zs(3*n+1),
     2    weights,weights(n+1),x1,y1,x2,y2,x3,y3,x4,y4)
c 
c        construct the weights
c 
        coe1=qsqrt((x2-x1)**2+(y2-y1)**2) /2
        coe2=qsqrt((x3-x2)**2+(y3-y2)**2) /2
        coe3=qsqrt((x4-x3)**2+(y4-y3)**2) /2
        coe4=qsqrt((x4-x1)**2+(y4-y1)**2) /2
c 
cccc        coe2=coe2*10
c 
        do 1200 i=1,n
        weights(i)=w(i) * coe1
        weights(n+i)=w(i)  * coe2
        weights(2*n+i)=w(i)  * coe3
        weights(3*n+i)=w(i)  * coe4
 1200 continue
c 
c       construct the cumulative lenghts along the square
c 
        rlens(1)=cqabs(zs(1)-x1-ima*y1)
        do 2200 i=2,n
        rlens(i)=rlens(i-1)+cqabs(zs(i)-zs(i-1))
 2200 continue
c 
        rlens(n+1)=rlens(n)+cqabs(zs(n)-x2-ima*y2)+
     1   +cqabs(zs(n+1)-x2-ima*y2)
        do 2400 i=n+2,2*n
        rlens(i)=rlens(i-1)+cqabs(zs(i)-zs(i-1))
 2400 continue
  
        rlens(2*n+1)=rlens(2*n)+cqabs(zs(2*n)-x3-ima*y3)+
     1   +cqabs(zs(2*n+1)-x3-ima*y3)
        do 2600 i=2*n+2,3*n
        rlens(i)=rlens(i-1)+cqabs(zs(i)-zs(i-1))
 2600 continue
c 
        rlens(3*n+1)=rlens(3*n)+cqabs(zs(3*n)-x4-ima*y4)+
     1   +cqabs(zs(3*n+1)-x4-ima*y4)
        do 2800 i=3*n+2,4*n
        rlens(i)=rlens(i-1)+cqabs(zs(i)-zs(i-1))
 2800 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine squacon(n,t,w,
     1    z12,z23,z34,z41,wwx,wwy,x1,y1,x2,y2,x3,y3,x4,y4)
        implicit real *16 (a-h,o-z)
        save
        dimension t(1),w(1),z12(2,1),z23(2,1),
     1   z34(2,1),z41(2,1),wwx(1),wwy(1)
c 
c        this subroutine constructs a quadrilateral polygon
c        from four user-provided corners. it constructs
c        discretizations of the four sides at gaussian
c        points on each side. it also constructs the
c        nodes and weights of the n-point gaussian quadrature,
c        to be used for the integration of functions on the
c        sides of the polygon.
c 
c                     input parameters:
c 
c  n - the number of gaussian nodes on each side of the polygon
c  (x1,y1) - coordinates of the first vertex of the polygon.
c  (x2,y2) - coordinates of the second vertex of the polygon.
c  (x3,y3) - coordinates of the third vertex of the polygon.
c  (x4,y4) - coordinates of the fourth vertex of the polygon.
c 
c                     output parameters:
c 
c  z12 - the discretization of first leg of the the polygon
c  z23 - the discretization of second leg of the the polygon
c  z34 - the discretization of third leg of the the polygon
c  z41 - the discretization of fourth leg of the the polygon
c  t - the n gaussian points on the interval [-1,1]
c  w - the n gaussian weights on the interval [-1,1]
c 
c                     work arrays:
c 
c  wwx,wwx - each must be at least n real *16 locations long
c 
c        . . . construct the points on the square
c 
c        construct the gaussian nodes
c 
        kind=1
        alpha=0
        beta=0
        kpts=0
  
  
        itype=1
        call legeexps(itype,n,t,u,v,w)
c 
cccc        CALL GAUSSQ(KIND,N,ALPHA,BETA,KPTS,ENDPTS,wwx,T,W)
c 
c       . . . the sides
c 
        call creside(t,n,x1,x2,wwx)
        call creside(t,n,y1,y2,wwy)
        call arrmove(wwx,wwy,n,z12)
  
c 
        call creside(t,n,x2,x3,wwx)
        call creside(t,n,y2,y3,wwy)
        call arrmove(wwx,wwy,n,z23)
c 
        call creside(t,n,x3,x4,wwx)
        call creside(t,n,y3,y4,wwy)
        call arrmove(wwx,wwy,n,z34)
c 
        call creside(t,n,x4,x1,wwx)
        call creside(t,n,y4,y1,wwy)
        call arrmove(wwx,wwy,n,z41)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine arrmove(x,y,n,z)
        implicit real *16 (a-h,o-z)
        save
        dimension x(1),y(1),z(2,1)
c 
        do 1200 i=1,n
        z(1,i)=x(i)
        z(2,i)=y(i)
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine creside(t,n,x1,x2,x)
        implicit real *16 (a-h,o-z)
        save
        dimension t(1),x(1)
c 
c       construct the linear mapping transforming the interval
c       (-1,1) into the interval (x1,x2)
c 
        beta=(x2+x1)/2
        alpha=(x2-x1)/2
c 
c       map the gaussian nodes from the interval
c       [-2,2] onto the interval [x1,x2]
c 
        do 1200 i=1,n
        x(i)=alpha*t(i)+beta
 1200    continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine scapro(w,n,a,b,prod)
        implicit real *16 (a-h,o-z)
        save
        dimension w(1)
        complex *32 a(1),b(1),prod
c 
c       calculate inner product of two vectors
c 
        prod=0
        do 1200 i=1,n
        prod=prod+qconjg(a(i))*b(i)*w(i)
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine polplot(pol,n,m,rlens,w,iw1,iw2)
        implicit real *16 (a-h,o-z)
        save
        complex *32 pol(n,1),com
        dimension rlen(1),w(1),rea(2),rlens(1)
        equivalence (rea(1),com)
c 
c       plot the real and imaginary parts of the m-th
c       orthogonal polynomial
c 
cccc        call prinq('in polplot, pol=*',pol,n*m*2)
        do 1200 i=1,n
        com=pol(i,m)
        w(i)=rea(1)
 1200 continue
        call RSGRAF(rlens,w,N,IW1,'real part of p_m is*')
c 
        do 1400 i=1,n
        com=pol(i,m)
        w(i)=rea(2)
 1400 continue
        call RSGRAF(rlens,w,N,IW2,'imaginary part of p_m is*')
        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c 
c 
        subroutine legeexps(itype,n,x,u,v,whts)
        implicit real *16 (a-h,o-z)
        save
        dimension x(1),whts(1),u(n,n),v(n,n)
  
c 
c         this subroutine constructs the gaussiaqn nodes
c         on the interval [-1,1], and the weights for the
c         corresponding order n quadrature. it also constructs
c         the matrix v converting the coefficients
c         of a legendre expansion into its values at the n
c         gaussian nodes, and its inverse u, converting the
c         values of a function at n gaussian nodes into the
c         coefficients of the corresponding legendre series.
c         no attempt has been made to make this code efficient,
c         but its speed is normally sufficient, and it is
c         mercifully short.
c 
c                 input parameters:
c 
c  itype - the type of the calculation to be performed
c          itype=0 means that only the gaussian nodes are
c                  to be constructed.
c          itype=1 means that only the nodes and the weights
c                  are to be constructed
c          itype=2 means that the nodes, the weights, and
c                  the matrices u, v are to be constructed
c  n - the number of gaussian nodes and weights to be generated
c 
c                 output parameters:
c 
c  x - the order n gaussian nodes - computed independently
c          of the value of itype.
c  u - the n*n matrix converting the  values at of a polynomial of order
c         n-1 at n legendre nodes into the coefficients of its
c         legendre expansion - computed only in itype=2
c  v - the n*n matrix converting the coefficients
c         of an n-term legendre expansion into its values at
c         n legendre nodes (note that v is the inverse of u)
c          - computed only in itype=2
c  whts - the corresponding quadrature weights - computed only
c         if itype .ge. 1
c 
c       . . . construct the nodes and the weights of the n-point gaussian
c             quadrature
c 
        ifwhts=0
        if(itype. gt. 0) ifwhts=1
        call legewhts(n,x,whts,ifwhts)
c 
c       construct the matrix of values of the legendre polynomials
c       at these nodes
c 
        if(itype .ne. 2) return
        do 1400 i=1,n
c 
        call legepols(x(i),n-1,u(1,i) )
 1400 continue
c 
        do 1800 i=1,n
        do 1600 j=1,n
        v(i,j)=u(j,i)
 1600 continue
 1800 continue
c 
c       now, v converts coefficients of a legendre expansion
c       into its values at the gaussian nodes. construct its
c       inverse u, converting the values of a function at
c       gaussian nodes into the coefficients of a legendre
c       expansion of that function
c 
        do 2800 i=1,n
        d=1
        d=d*(2*i-1)/2
        do 2600 j=1,n
        u(i,j)=v(j,i)*whts(j)*d
 2600 continue
 2800 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine legewhts(n,ts,whts,ifwhts)
        implicit real *16 (a-h,o-z)
        save
        dimension ts(1),whts(1)
c 
c        this subroutine constructs the nodes and the
c        weights of the n-point gaussian quadrature on
c        the interval [-1,1]
c 
c                input parameters:
c 
c  n - the number of nodes in the quadrature
c 
c                output parameters:
c 
c  ts - the nodes of the n-point gaussian quadrature
c  w - the weights of the n-point gaussian quadrature
c 
c       . . . construct the array of initial approximations
c             to the roots of the n-th legendre polynomial
c 
        eps=1.0d-30
        ZERO=0
        DONE=1
        pi=qatan(done)*4
        h=pi/(2*n)
        do 1200 i=1,n
        t=(2*i-1)*h
        ts(n-i+1)=qcos(t)
1200  CONTINUE
c 
c         use newton to find all roots of the legendre polynomial
c 
        ts(n/2+1)=0
        do 2000 i=1,n/2
c 
        xk=ts(i)
        deltold=1
        do 1400 k=1,10
        call legepol(xk,n,pol,der)
        delta=-pol/der
cccc         call prinq('delta=*',delta,1)
        xk=xk+delta
        if(deltold .le. eps) goto 1600
        deltold=qabs(delta)
 1400 continue
 1600 continue
        ts(i)=xk
        ts(n-i+1)=-xk
 2000 continue
c 
c       now, use the explicit integral formulae
c       to obtain the weights
c 
        if(ifwhts .eq. 0) return
        a=-1
        b=1
        do 2200 i=1,n/2+1
        call prodend(a,ts,n,i,fm)
        call prodend(b,ts,n,i,fp)
        whts(i)=fp-fm
        whts(n-i+1)=whts(i)
 2200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine legepol(x,n,pol,der)
        implicit real *16 (a-h,o-z)
c 
        save
        pkm1=1
        pk=x
c 
        pk=1
        pkp1=x
c 
c        if n=0 or n=1 - exit
c 
        if(n .ge. 2) goto 1200
        pol=1
        der=0
        if(n .eq. 0) return
c 
        pol=x
        der=1
        return
 1200 continue
c 
c       n is greater than 1. conduct recursion
c 
        do 2000 k=1,n-1
        pkm1=pk
        pk=pkp1
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
 2000 continue
c 
c        calculate the derivative
c 
        pol=pkp1
        der=n*(x*pkp1-pk)/(x**2-1)
        return
        end
c 
c 
c 
c 
c 
        subroutine prodend(x,xs,n,i,f)
        implicit real *16 (a-h,o-z)
        save
        dimension xs(1)
c 
c      evaluate the product
c 
        f=1
        dlarge=1.0d20
        dsmall=f/dlarge
c 
        large=0
        do 2000 j=1,n
        dd=qabs(f)
        if( dd .gt. dsmall) goto 1200
         f=f*10000
         large=large-1
 1200 continue
c 
        if( dd .lt. dlarge) goto 1400
        f=f/10000
        large=large+1
 1400 continue
        if(j .eq. i) goto 2000
        f=f*(x-xs(j))/(xs(i)-xs(j))
 2000 continue
        d10000=10000
        f=f*d10000**large
        f=f**2*(x-xs(i))
        return
        end
c 
c 
c 
c 
c 
        subroutine legepols(x,n,pols)
        implicit real *16 (a-h,o-z)
        save
        dimension pols(1)
c 
        pkm1=1
        pk=x
c 
        pk=1
        pkp1=x
c 
c 
c        if n=0 or n=1 - exit
c 
        if(n .ge. 2) goto 1200
        pols(1)=1
        if(n .eq. 0) return
c 
        pols(2)=x
        return
 1200 continue
c 
        pols(1)=1
        pols(2)=x
c 
c       n is greater than 2. conduct recursion
c 
        do 2000 k=1,n-1
        pkm1=pk
        pk=pkp1
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
        pols(k+2)=pkp1
 2000 continue
c 
        return
        end
c 
c 
c 
c 
c 
      SUBROUTINE legeFDER(X,VAL,der,PEXP,N)
      IMPLICIT REAL *16 (A-H,O-Z)
        save
      REAL *16 PEXP(1)
C 
C     This subroutine computes the value and the derivative
c     of a gaussian expansion with coefficients PEXP
C     at point X in interval [-1,1]
C 
c                input parameters:
c 
C     X = evaluation point
C     PEXP = expansion coefficients
C     N  = order of expansion (0 ,..., N-1)
c 
c                output parameters:
c 
C     VAL = computed value
C     der = computed value of the derivative
C 
C 
        done=1
        pjm2=1
        pjm1=x
        derjm2=0
        derjm1=1
c 
        val=pexp(1)*pjm2+pexp(2)*pjm1
        der=pexp(2)
c 
        dd=x**1-1
        dd=done/dd
        DO 600 J = 2,N
c 
        pj= ( (2*j-1)*x*pjm1-(j-1)*pjm2 ) / j
        val=val+pexp(j+1)*pj
c 
        derj=j*(x*pj-pjm1)*dd
        der=der+pexp(j+1)*derj
c 
        pjm2=pjm1
        pjm1=pj
        derjm2=derjm1
        derjm1=derj
 600   CONTINUE
c 
      RETURN
      END
c 
c 
c 
c 
c 
c 
c 
c 
c 
cccc        SUBROUTINE RSplot(Z,N,IW)
        SUBROUTINE RSPLOT(X,Y,N,IW,mes)
        IMPLICIT REAL *16 (A-H,O-Z)
        save
        DIMENSION X(1),Y(1),Z(2,1)
        character *1 mes(1)
C 
C        THIS SUBROUTINE PLOTS THE CURVE SPECIFIED BY ITS NODES
C        X,Y (CONCEPTUALLY, THE CURVE CONNECTS THE POINT (X(1),Y(1))
C        WITH THE POINT (X(N),Y(N)). IF THE OUTPUT FILE IS 'INCLUDE' D
C        INTO TD OR TD3 STREAM, A FILE NAMED "PLOT.IMP" IS CREATED.
C        THIS IS A TRUE, UNADULTERATED 'IMPRINT' FILE, WHICH CAN AND
C        SHOULD BE "IMPRINTED".
C     NOTE: NATURALLY, CERTAIN DECISIONS WERE MADE IN WRITING THIS
C        ROUTINE (FRAME, SIZE, AND OTHER PARAMETERS OF THE PLOT.
C        IF YOU DON'T LIKE THEM - YOU ARE FREE TO WRITE YOUR OWN
C        SUBROUTINE!!!
C 
C  INPUT PARAMETERS:
C   X,Y - THE COORDINATES OF THE CURVE TO PLOT
C   N - THE NUMBER OF ELEMENTS IN ARAYS X,Y.
C   IW - THE FORTRAN UNIT NUMBER ON WHICH THE OUTPUT DATA SET IS WRITTEN
C 
C  OUTPUT PARAMETERS: NONE
C 
c 
c       . . . print the label
c 
        IF(IW.EQ.0) RETURN
        call titlpr(MES,Iw)
C 
C 
C . . . FIND THE LIMITS TO BE SET
C 
        XMIN=1.0D20
        YMIN=1.0D20
        XMAX=-1.0D20
        YMAX=-1.0D20
        DO 1000 I=1,N
       IF(XMIN.GT.X(I))XMIN=X(I)
       IF(XMAX.LT.X(I))XMAX=X(I)
       IF(YMIN.GT.Y(I))YMIN=Y(I)
       IF(YMAX.LT.Y(I))YMAX=Y(I)
ccccc        IF(XMIN.GT.Z(1,I))XMIN=Z(1,I)
ccccc        IF(XMAX.LT.Z(1,I))XMAX=Z(1,I)
ccccc        IF(YMIN.GT.Z(2,I))YMIN=Z(2,I)
ccccc        IF(YMAX.LT.Z(2,I))YMAX=Z(2,I)
 1000 CONTINUE
        DX=(XMAX-XMIN)
        DY=YMAX-YMIN
        DD=DX
        IF(DY.GT.DX) DD=DY
        DMIN=XMIN
        IF(YMIN.LT.XMIN) DMIN=YMIN
        DMAX=XMAX
        IF(YMAX.GT.XMAX) DMAX=YMAX
        DMIN=DMIN-DD/10
        DMAX=DMAX+DD/10
 1200 FORMAT(1X,37H#SET DEVICE postscript FILE 'plot.ps',/,
c1200 FORMAT(1X,34H#SET DEVICE IMAGEN FILE 'PLOT.IMP',/,
C1200 FORMAT(1X,26HSET IMAGEN FILE 'PLOT.IMP',/,
ccc     A       1X,'SET DEVICE SUN',/,
C    1  1X,'SET DEV MC',/,
     2  1X,'SET AXIS TOP OFF',/,
     3  1X,'SET AXIS RIGHT OFF',/,
     4  1X,'SET LIMITS X ', 2(1X,E11.5),' Y ',2(1X,E11.5),/,
CCC  5  1X,'SET LIMITS Y ', 2(2X,E11.5),/,
     5  1X,'SET SCALE EQUAL')
          IF(IW.EQ.0) RETURN
          WRITE(IW,1200) DMIN,DMAX,DMIN,DMAX
 1400 FORMAT(2(2X,E15.7))
ccccc       WRITE(IW,1400) (Z(1,I),Z(2,I),I=1,N)
       WRITE(IW,1400) (X(I),Y(I),I=1,N)
cccc 1500 FORMAT(1X,'plot ''*'' ')
 1500 FORMAT(1X,'PLOT')
C1600 FORMAT(1X,'FRAME')
C1800 FORMAT(1X,'EXIT')
 1450 FORMAT(1X,'JOIN')
c1500 FORMAT(1X,'PLOT')
 1600 FORMAT(1X,'#FRAME')
 1800 FORMAT(1X,'#EXIT')
        WRITE(IW,1500)
        WRITE(IW,1600)
        WRITE(IW,1800)
        RETURN
C 
C 
C 
C 
        ENTRY RSGRAF(X,Y,N,IW,mes)
C 
C        THIS SUBROUTINE PLOTS THE Y AS A FUNCTION OF X.
C        IF THE OUTPUT FILE IS 'INCLUDE' D
C        INTO TD OR TD3 STREAM, A FILE NAMED "PLOT.IMP" IS CREATED.
C        THIS IS A TRUE, UNADULTERATED 'IMPRINT' FILE, WHICH CAN AND
C        SHOULD BE "IMPRINTED".
C     NOTE: NATURALLY, CERTAIN DECISIONS WERE MADE IN WRITING THIS
C        ROUTINE (FRAME, SIZE, AND OTHER PARAMETERS OF THE PLOT).
C        IF YOU DON'T LIKE THEM - YOU ARE FREE TO WRITE YOUR OWN
C        SUBROUTINE!!!
C 
C  INPUT PARAMETERS:
C   X,Y - THE COORDINATES OF THE CURVE TO PLOT
C   N - THE NUMBER OF ELEMENTS IN ARAYS X,Y.
C   IW - THE FORTRAN UNIT NUMBER ON WHICH THE OUTPUT DATA SET IS WRITTEN
C 
C  OUTPUT PARAMETERS: NONE
C 
        IF(IW.EQ.0) RETURN
        call titlpr(MES,Iw)
 2200 FORMAT(1X,37H#SET DEVICE postscript FILE 'plot.ps',/,
c2200 FORMAT(1X,26HSET IMAGEN FILE 'PLOT.IMP',/,
ccc     A       1X,'SET DEVICE SUN',/,
ccc  1  1X,'SET DEV MC',/,
     2  1X,'SET AXIS TOP OFF',/,
     3  1X,'SET AXIS RIGHT OFF')
CC   3  1X,'SET AXIS RIGHT OFF',/,
CC   4  1X,'SET SCALE EQUAL')
          WRITE(IW,2200)
          WRITE(IW,1400) (X(I),Y(I),I=1,N)
          WRITE(IW,1450)
          WRITE(IW,1600)
          WRITE(IW,1800)
          RETURN
          END
c 
c 
c 
c 
c 
        SUBROUTINE titlpr(MES,IP)
        save
        CHARACTER *1 MES(1),AST,card(80),blank,title(5),quo
        DATA AST/'*'/,quo/''''/,blank/' '/
     1      title/'t','i','t','l','e'/
C 
C         DETERMINE THE LENGTH OF THE MESSAGE
C 
        I=0
        DO 1400 I=1,10000
        IF(MES(I).EQ.AST) GOTO 1600
        I1=I
 1400 CONTINUE
 1600 CONTINUE
 1800 FORMAT(1X, 80A1)
         do 2000 i=1,80
         card(i)=blank
 2000 continue
         do 2200 i=1,5
         card(i+1)=title(i)
 2200 continue
         card(8)=quo
c 
         do 2400 i=1,i1
         card(i+8)=mes(i)
 2400 continue
         card(8+i1+1)=quo
c 
         write(ip,1800) card
         return
         end
c 
c 
c 
c 
c 
        subroutine prods(pol,n,m,w,out)
        implicit real *16 (a-h,o-z)
        save
        dimension w(1)
        complex *32 pol(n,1),out(m,1)
c 
c       calculate the inner products between m vectors of length n
c 
        do 1400 i=1,m
        do 1200 j=1,m
        call scapro(w,n,pol(1,i),pol(1,j),out(i,j))
 1200 continue
 1400 continue
        return
        end
c 
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code and the beginning
c        of the orthogonal polynomial subroutine itself.
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
  
        subroutine qpsquare(z,n,pols)
        implicit real *16 (a-h,o-z)
        save
        complex *32 pols(1),w(50),z
c 
c       this subroutine evaluates at the point z \in C
c       a set of polynomials orthogonal on the square
c       with the corners at the points (1,-1),(1,2),
c       (-1,1), (-1,-1).
c 
c                   input parameters:
c 
c  z - the point where the polynomials are to be evaluarted
c  n - the order of the last polynomial to be evaluated;
c       must be .leq. 100
c 
c                   output parameters:
c 
c  pols - the array of values at z of polynomials orthogonal
c       on the square.
c 
         real *16
     1  c0o1001(16),c0o1002(16),c0o1003(16),c0o1004(16),
     2  c0o1005(16),c0o1006(16),c0o1007(16),c0o1008(16),
     3  c0o1009(16),c0o1010(16),c0o1011(16),c0o1012(16),
     4  c0o1013(16),c0o1014(16),c0o1015(16),c0o1016(16)
         real *16
     1  c0o1017(16),c0o1018(16),c0o1019(12)
         real *16
     1  c0o2001(2),c0o2002(2),c0o2003(2),c0o2004(2),
     2  c0o2005(2),c0o2006(2),c0o2007(2),c0o2008(2),
     3  c0o2009(2),c0o2010(2),c0o2011(2),c0o2012(2),
     4  c0o2013(2),c0o2014(2),c0o2015(2),c0o2016(2)
         real *16
     1  c0o2017(2),c0o2018(2)
         equivalence
     1  (c0o1001(16),c0o2001(1)),(c0o1002(1),c0o2001(2)),
     2  (c0o1002(16),c0o2002(1)),(c0o1003(1),c0o2002(2)),
     3  (c0o1003(16),c0o2003(1)),(c0o1004(1),c0o2003(2)),
     4  (c0o1004(16),c0o2004(1)),(c0o1005(1),c0o2004(2)),
     5  (c0o1005(16),c0o2005(1)),(c0o1006(1),c0o2005(2)),
     6  (c0o1006(16),c0o2006(1)),(c0o1007(1),c0o2006(2)),
     7  (c0o1007(16),c0o2007(1)),(c0o1008(1),c0o2007(2)),
     8  (c0o1008(16),c0o2008(1)),(c0o1009(1),c0o2008(2)),
     9  (c0o1009(16),c0o2009(1)),(c0o1010(1),c0o2009(2)),
     a  (c0o1010(16),c0o2010(1)),(c0o1011(1),c0o2010(2))
         equivalence
     1  (c0o1011(16),c0o2011(1)),(c0o1012(1),c0o2011(2)),
     2  (c0o1012(16),c0o2012(1)),(c0o1013(1),c0o2012(2)),
     3  (c0o1013(16),c0o2013(1)),(c0o1014(1),c0o2013(2)),
     4  (c0o1014(16),c0o2014(1)),(c0o1015(1),c0o2014(2)),
     5  (c0o1015(16),c0o2015(1)),(c0o1016(1),c0o2015(2)),
     6  (c0o1016(16),c0o2016(1)),(c0o1017(1),c0o2016(2)),
     7  (c0o1017(16),c0o2017(1)),(c0o1018(1),c0o2017(2)),
     8  (c0o1018(16),c0o2018(1)),(c0o1019(1),c0o2018(2))
         data c0o1001/
     1      0.80000000000000000000000000000000Q+00,
     2      -.60168114647228395617615718571876Q+00,
     3      0.13360457724094087730451366815003Q+01,
     4      0.29238350307129225611664122706641Q+00,
     5      -.51424935915482851072034282933979Q+00,
     6      0.12884215877396131000227444189043Q+01,
     7      -.14840544494263769576752771365307Q+00,
     8      0.19232020848719077699948929701190Q+00,
     9      -.46488295131879336836988859709781Q+00,
     a      0.12875421891869047604387000595992Q+01,
     1      0.86847066717594456214183340624058Q-01,
     2      -.92025549806677452434214246310576Q-01,
     3      0.15631437207356396046852919157395Q+00,
     4      -.45986469752206439618247028699794Q+00,
     5      0.12897738485781815213430720944679Q+01,
     6      -.56613701274149330029310176997694Q-01/
         data c0o1002/
     1      0.54463366436480590081721579479806Q-01,
     2      -.66827579278540837104010622853574Q-01,
     3      0.14998377116535163786839553536482Q+00,
     4      -.45994409871857183596019048740205Q+00,
     5      0.12911343628431089384544215924099Q+01,
     6      0.39736165017987716816267523442718Q-01,
     7      -.36470614249711759287770079263964Q-01,
     8      0.36591961219969463142834223578419Q-01,
     9      -.60952355907110944772467432975393Q-01,
     a      0.14885477482439715540158630981998Q+00,
     1      -.46039999337268902212934789161326Q+00,
     2      0.12919505751056088120475407954285Q+01,
     3      -.29399004124558953756172469397314Q-01,
     4      0.26315874333343321310594830254071Q-01,
     5      -.23430900357263456185639296886362Q-01,
     6      0.31613543384247571455420823504773Q-01/
         data c0o1003/
     1      -.59427786751094774167457257321133Q-01,
     2      0.14867368701315501036607395786490Q+00,
     3      -.46078277091333090512436011101123Q+00,
     4      0.12924686746172530660570747621321Q+01,
     5      0.22621004508177874873220841809630Q-01,
     6      -.19958905162346392224868139484840Q-01,
     7      0.16513151175633845570986580091801Q-01,
     8      -.19325645741469654454192881893538Q-01,
     9      0.30058280737368355669952106396848Q-01,
     a      -.58950409656854481888328837032891Q-01,
     1      0.14870791151326211074200401811914Q+00,
     2      -.46106798418340934665917418461823Q+00,
     3      0.12928154965198537557679674469584Q+01,
     4      -.17940384968220442963862944110250Q-01,
     5      0.15689685501457592344726273636290Q-01,
     6      -.12385590152959703849371440859009Q-01/
         data c0o1004/
     1      0.13142738415775158969611644489467Q-01,
     2      -.17879217326885905884181878295942Q-01,
     3      0.29468429659067605657000735251318Q-01,
     4      -.58797231017818669842152516056957Q-01,
     5      0.14878834511584183407702696335018Q+00,
     6      -.46127815574287766341956949990759Q+00,
     7      0.12930581800246665862720250453854Q+01,
     8      0.14574246230399726271039641191457Q-01,
     9      -.12673177016103895676610143625983Q-01,
     a      0.96972116794043660040552260428246Q-02,
     1      -.96055943487158044598590316696790Q-02,
     2      0.11848128359338131770348046554727Q-01,
     3      -.17267077933106481882636164818095Q-01,
     4      0.29223038372269294031057856725262Q-01,
     5      -.58758474422441288272949448226917Q-01,
     6      0.14887027438400444574927346890002Q+00/
         data c0o1005/
     1      -.46143498946576165637100034515655Q+00,
     2      0.12932342829156132910432165141414Q+01,
     3      -.12073188367644651366354282987465Q-01,
     4      0.10457975108766684383659474874318Q-01,
     5      -.78326836419596865285684076934545Q-02,
     6      0.73852417124932169876820221655047Q-02,
     7      -.84658388742663601303254570881895Q-02,
     8      0.11256674540404581035700611291350Q-01,
     9      -.16981647811867185228295169171308Q-01,
     a      0.29118054050298191086640334442718Q-01,
     1      -.58762101604633896255240401910665Q-01,
     2      0.14894240031365929416799186920229Q+00,
     3      -.46155416456075881955040777559555Q+00,
     4      0.12933659680912759962350034802334Q+01,
     5      0.10164491130217015246816401458909Q-01,
     6      -.87809780561093690653420354801521Q-02/
         data c0o1006/
     1      0.64778409337365775673011378174612Q-02,
     2      -.58916714114497206583103777597551Q-02,
     3      0.63880319090109584619727118783220Q-02,
     4      -.79139990445404848560451824669755Q-02,
     5      0.10960112172158354773876091594463Q-01,
     6      -.16841182889719401590500410857822Q-01,
     7      0.29074983827488105428872713448541Q-01,
     8      -.58780919201798188794764410323275Q-01,
     9      0.14900310487279548661612056426344Q+00,
     a      -.46164641989330473198273899543635Q+00,
     1      0.12934669409358940707196502232684Q+01,
     2      -.86749414008400806717093522375129Q-02,
     3      0.74796951008349928744067798475966Q-02,
     4      -.54575780664087646376517865741347Q-02,
     5      0.48325384469559182553747026887451Q-02,
     6      -.50200089235020644745574484830474Q-02/
         data c0o1007/
     1      0.58824621709547997390001439423964Q-02,
     2      -.76221381097529657526156706517834Q-02,
     3      0.10802636061230746889887296930210Q-01,
     4      -.16770302012383377768700554080490Q-01,
     5      0.29060543569820750354583517853305Q-01,
     6      -.58804102830658510892522013919263Q-01,
     7      0.14905349820001861418378490671206Q+00,
     8      -.46171908871897854367483851196483Q+00,
     9      0.12935460244789805937969833607811Q+01,
     a      0.74902715036151197360108305321499Q-02,
     1      -.64490513630143052976742429930831Q-02,
     2      0.46673550480289187904778318123230Q-02,
     3      -.40499596434601206228788557430526Q-02,
     4      0.40692154371249434307301575737978Q-02,
     5      -.45612180514536029468351129606634Q-02,
     6      0.56036171823321866629331310744542Q-02/
         data c0o1008/
     1      -.74588700812583176646708878305435Q-02,
     2      0.10715846608099088691497116881408Q-01,
     3      -.16734729444710027199877453541376Q-01,
     4      0.29059604556242526103089121559448Q-01,
     5      -.58827292090988120541470729200331Q-01,
     6      0.14909526471647272214813493443187Q+00,
     7      -.46177724054294988844249413704667Q+00,
     8      0.12936090979885431449130372497157Q+01,
     9      -.65326414675975187219563430124452Q-02,
     a      0.56185063185194258798410697988258Q-02,
     1      -.40412458807227142939852777105795Q-02,
     2      0.34525548270618485791909691368417Q-02,
     3      -.33792884011528924122421551926564Q-02,
     4      0.36547816934112330151479779516825Q-02,
     5      -.42993534947458049538683156823141Q-02,
     6      0.54413142070355599507979415673383Q-02/
         data c0o1009/
     1      -.73638647998257183663169795806619Q-02,
     2      0.10667004874176233878407151153693Q-01,
     3      -.16717782107080856472302628229339Q-01,
     4      0.29064942217889146097424983015978Q-01,
     5      -.58848814203327863113421765445710Q-01,
     6      0.14913000669070064343073006087496Q+00,
     7      -.46182444019278050562764457577164Q+00,
     8      0.12936601974130969711943757379871Q+01,
     9      0.57475482984872037137470687805245Q-02,
     a      -.49392036167697393977696488795481Q-02,
     1      0.35357885226273988360073774927419Q-02,
     2      -.29843478079318308917453438226893Q-02,
     3      0.28608843854399988314184039996014Q-02,
     4      -.30055474214092689896890285274427Q-02,
     5      0.34113601708158873232594948293751Q-02,
     6      -.41419781819487075454901618361138Q-02/
         data c0o1010/
     1      0.53430911037347686149355120434398Q-02,
     2      -.73070451113758711033228505683553Q-02,
     3      0.10639396649590646319732418577701Q-01,
     4      -.16710904482809669782805114014547Q-01,
     5      0.29072991895346873675997786420875Q-01,
     6      -.58868143594581882119695720137751Q-01,
     7      0.14915907620963071794582735163666Q+00,
     8      -.46186324007322944527583774412301Q+00,
     9      0.12937021658138896750311414721343Q+01,
     a      -.50959098377566674985960372007431Q-02,
     1      0.43763997223238933779943901556552Q-02,
     2      -.31212645976859084784807992030140Q-02,
     3      0.26093953376435748981872102159314Q-02,
     4      -.24600650513341882890857525009020Q-02,
     5      0.25238032755514964918777344236982Q-02,
     6      -.27806188545698104078500932158501Q-02/
         data c0o1011/
     1      0.32611152373227999925757131990924Q-02,
     2      -.40437616572181660898894902811066Q-02,
     3      0.52819039190890439532174492104578Q-02,
     4      -.72724719118379423344835689980826Q-02,
     5      0.10624046689403993502970982314707Q-01,
     6      -.16709527948139433681752033620853Q-01,
     7      0.29081965787648204944729211519660Q-01,
     8      -.58885242339124497561137814601723Q-01,
     9      0.14918356507163166843195975819151Q+00,
     a      -.46189550010948456124093627457368Q+00,
     1      0.12937370517232967710135817532469Q+01,
     2      0.45491120080857171358105658316016Q-02,
     3      -.39048200678752877680692536876692Q-02,
     4      0.27767245022191663997979528255225Q-02,
     5      -.23036859877223840614442264474271Q-02,
     6      0.21427272665144565211364874493019Q-02/
         data c0o1012/
     1      -.21556805518932114754417100013744Q-02,
     2      0.23166495190630719507858639848417Q-02,
     3      -.26386003176911266085960640128705Q-02,
     4      0.31649475481138470411049243662266Q-02,
     5      -.39806755659805497526898528671141Q-02,
     6      0.52429744169149114501442937514458Q-02,
     7      -.72512823377942507585008390558677Q-02,
     8      0.10615935930326936239434837632631Q-01,
     9      -.16711136278194157555453823163324Q-01,
     a      0.29090970070127460593679751738442Q-01,
     1      -.58900268821578275753312880901364Q-01,
     2      0.14920434029790090742524133833399Q+00,
     3      -.46192259824058322523572805512615Q+00,
     4      0.12937663608838511174113835487206Q+01,
     5      -.40858166825882612902099937335302Q-02,
     6      0.35057118524701687956137710506916Q-02/
         data c0o1013/
     1      -.24870146000884484477472061978921Q-02,
     2      0.20506319231350780875818145166298Q-02,
     3      -.18864447409794162644939904483485Q-02,
     4      0.18673503793252164289824661551971Q-02,
     5      -.19651964033641111205972153725554Q-02,
     6      0.21832631662581653412588808728687Q-02,
     7      -.25457318379173592965177616799606Q-02,
     8      0.31016404492369450351357456265118Q-02,
     9      -.39392455753646202343105463908240Q-02,
     a      0.52178485059959730031219314724214Q-02,
     1      -.72383517157714123086153332366627Q-02,
     2      0.10612163732284442738145930557792Q-01,
     3      -.16714312118700789326102338342907Q-01,
     4      0.29099572139357860720694070832975Q-01,
     5      -.58913447434563213029447288928329Q-01,
     6      0.14922208722694911432399939640645Q+00/
         data c0o1014/
     1      -.46194557111183243308041387570071Q+00,
     2      0.12937912198054096692593857399020Q+01,
     3      0.36898497984721834970560680467346Q-02,
     4      -.31649177555082691923474328637082Q-02,
     5      0.22409293434849477226032094147441Q-02,
     6      -.18384403145493496235839761777875Q-02,
     7      0.16759631369209716394482117134675Q-02,
     8      -.16367413132553024779426481120045Q-02,
     9      0.16922670339912940636130222009701Q-02,
     a      -.18404247166706536428949388373711Q-02,
     1      0.20944137066895414069667854668723Q-02,
     2      -.24833288945431662900132752297705Q-02,
     3      0.30590307605259178686027168404832Q-02,
     4      -.39115765904027387067100195798574Q-02,
     5      0.52015056965152305750944789008814Q-02,
     6      -.72306193252935978973571334405047Q-02/
         data c0o1015/
     1      0.10611008406346857442658004837867Q-01,
     2      -.16718245602843817675804124955578Q-01,
     3      0.29107580809085561598559931627149Q-01,
     4      -.58925011023701482046654777598026Q-01,
     5      0.14923734825594246855124750582966Q+00,
     6      -.46196520970968655083976252081850Q+00,
     7      0.12938124846840159238541742395686Q+01,
     8      -.33487713638025784819316743057070Q-02,
     9      0.28715842041499433703536025797017Q-02,
     a      -.20300230945609549835899260086804Q-02,
     1      0.16585175218944624195412694691185Q-02,
     2      -.15005980389160638016683530084096Q-02,
     3      0.14489700030129977299803274413680Q-02,
     4      -.14757551825193059068019904804327Q-02,
     5      0.15758405941723314456887744622340Q-02,
     6      -.17559600176592594332813285755276Q-02/
         data c0o1016/
     1      0.20336520221670886627807938876738Q-02,
     2      -.24404748137079696261169584727948Q-02,
     3      0.30298432975547367080382008916187Q-02,
     4      -.38928740717837484632037732838513Q-02,
     5      0.51908699987107665733705635770646Q-02,
     6      -.72262069584911569346569215603634Q-02,
     7      0.10611425999456250793113502382679Q-01,
     8      -.16722471491258488651529737633206Q-01,
     9      0.29114931515451714749775965078618Q-01,
     a      -.58935176845377029994737327903135Q-01,
     1      0.14925055446820745981355460743911Q+00,
     2      -.46198212546402415225289837219349Q+00,
     3      0.12938308155533857244846640576055Q+01,
     4      0.30528869165577830848261616209622Q-02,
     5      -.26172753725325831358432063642278Q-02,
     6      0.18478234097464735000766660950139Q-02/
         data c0o1017/
     1      -.15044654988000735764336341353383Q-02,
     2      0.13526720404849492165552213437397Q-02,
     3      -.12937072019285124463722780003345Q-02,
     4      0.13008145924857354388014986131480Q-02,
     5      -.13672625873860371519531826541405Q-02,
     6      0.14958921478044567397765891922356Q-02,
     7      -.16973040003332980331481160278106Q-02,
     8      0.19912090561751990536795010876882Q-02,
     9      -.24105220572506901887890202857713Q-02,
     a      0.30095741675150710309288934462617Q-02,
     1      -.38801378252327432076313823132829Q-02,
     2      0.51840062754731414946802471127524Q-02,
     3      -.72239355734657592177518311083441Q-02,
     4      0.10612772993534984953147186062232Q-01,
     5      -.16726724057171906410822083495391Q-01,
     6      0.29121625117972721179923301407940Q-01/
         data c0o1018/
     1      -.58944138031030136832072391368020Q-01,
     2      0.14926205044558105989417497019820Q+00,
     3      -.46199679667673572878279406055438Q+00,
     4      0.12938467277783528009562208208482Q+01,
     5      -.27945513143723037461895609817310Q-02,
     6      0.23953531594805291540684162958545Q-02,
     7      -.16892966444894859858954642488222Q-02,
     8      0.13714307745846062971995421008047Q-02,
     9      -.12265431548594881652970259142542Q-02,
     a      0.11636005019079074989377038680103Q-02,
     1      -.11571908071618943874877909300991Q-02,
     2      0.11997725759882420560370906083121Q-02,
     3      -.12918092187430840189255467723686Q-02,
     4      0.14396149443299077807861193040626Q-02,
     5      -.16557238159285543605690521947865Q-02,
     6      0.19610423002577441152583017928545Q-02/
         data c0o1019/
     1      -.23892855686273518848486573242126Q-02,
     2      0.29953536609466807797039483324281Q-02,
     3      -.38714419166986477693050732092098Q-02,
     4      0.51796692813961793462339439117592Q-02,
     5      -.72230508460234890351623829396864Q-02,
     6      0.10614648014949243315623400173687Q-01,
     7      -.16730854816819207831136642719786Q-01,
     8      0.29127694842195923741504947093461Q-01,
     9      -.58952062104687447246064538885437Q-01,
     a      0.14927211340629088783983004467772Q+00,
     1      -.46200960163050027079298984035060Q+00,
     2      0.12938606284603542796627054212576Q+01/
         real *16
     1  c0n1001(16),c0n1002( 9)
         real *16
     1  c0n2001(2)
         equivalence
     1  (c0n1001(16),c0n2001(1)),(c0n1002(1),c0n2001(2))
         data c0n1001/
     1      0.35355339059327376220042218105236Q+00,
     2      0.52882132014165582085795065151063Q+00,
     3      0.51647995295471645849350956944281Q+00,
     4      0.51659725219953458695448641487232Q+00,
     5      0.51599100859060043358204646628871Q+00,
     6      0.51569104537220442711236613655000Q+00,
     7      0.51553171362430247283499295669508Q+00,
     8      0.51543791936466335082699951218393Q+00,
     9      0.51537826598289081579539851885136Q+00,
     a      0.51533803899495889530775642344854Q+00,
     1      0.51530965121643953170902224645058Q+00,
     2      0.51528888064445106656576029982587Q+00,
     3      0.51527322973216200641175443323523Q+00,
     4      0.51526114554879346150544695700894Q+00,
     5      0.51525162180694704421783924651656Q+00,
     6      0.51524398339914210323255823829549Q+00/
         data c0n1002/
     1      0.51523776379906155623016330940541Q+00,
     2      0.51523263228734505355227293289863Q+00,
     3      0.51522834915050475148513550349916Q+00,
     4      0.51522473721314205045072346746284Q+00,
     5      0.51522166323916688717243273962739Q+00,
     6      0.51521902547849836877914121028240Q+00,
     7      0.51521674514419686285592363935884Q+00,
     8      0.51521476046403341718596161983149Q+00,
     9      0.51521302245483244144601212618395Q+00/
         real *16
     1  c1o1001(16),c1o1002(16),c1o1003(16),c1o1004(16),
     2  c1o1005(16),c1o1006(16),c1o1007(16),c1o1008(16),
     3  c1o1009(16),c1o1010(16),c1o1011(16),c1o1012(16),
     4  c1o1013(16),c1o1014(16),c1o1015(16),c1o1016(16)
         real *16
     1  c1o1017(16),c1o1018(16),c1o1019(12)
         real *16
     1  c1o2001(2),c1o2002(2),c1o2003(2),c1o2004(2),
     2  c1o2005(2),c1o2006(2),c1o2007(2),c1o2008(2),
     3  c1o2009(2),c1o2010(2),c1o2011(2),c1o2012(2),
     4  c1o2013(2),c1o2014(2),c1o2015(2),c1o2016(2)
         real *16
     1  c1o2017(2),c1o2018(2)
         equivalence
     1  (c1o1001(16),c1o2001(1)),(c1o1002(1),c1o2001(2)),
     2  (c1o1002(16),c1o2002(1)),(c1o1003(1),c1o2002(2)),
     3  (c1o1003(16),c1o2003(1)),(c1o1004(1),c1o2003(2)),
     4  (c1o1004(16),c1o2004(1)),(c1o1005(1),c1o2004(2)),
     5  (c1o1005(16),c1o2005(1)),(c1o1006(1),c1o2005(2)),
     6  (c1o1006(16),c1o2006(1)),(c1o1007(1),c1o2006(2)),
     7  (c1o1007(16),c1o2007(1)),(c1o1008(1),c1o2007(2)),
     8  (c1o1008(16),c1o2008(1)),(c1o1009(1),c1o2008(2)),
     9  (c1o1009(16),c1o2009(1)),(c1o1010(1),c1o2009(2)),
     a  (c1o1010(16),c1o2010(1)),(c1o1011(1),c1o2010(2))
         equivalence
     1  (c1o1011(16),c1o2011(1)),(c1o1012(1),c1o2011(2)),
     2  (c1o1012(16),c1o2012(1)),(c1o1013(1),c1o2012(2)),
     3  (c1o1013(16),c1o2013(1)),(c1o1014(1),c1o2013(2)),
     4  (c1o1014(16),c1o2014(1)),(c1o1015(1),c1o2014(2)),
     5  (c1o1015(16),c1o2015(1)),(c1o1016(1),c1o2015(2)),
     6  (c1o1016(16),c1o2016(1)),(c1o1017(1),c1o2016(2)),
     7  (c1o1017(16),c1o2017(1)),(c1o1018(1),c1o2017(2)),
     8  (c1o1018(16),c1o2018(1)),(c1o1019(1),c1o2018(2))
         data c1o1001/
     1      0.11428571428571428571428571428571Q+01,
     2      -.67753057814182123935424546537519Q+00,
     3      0.13190060599031575284873437907739Q+01,
     4      0.29194591934208887414962676110343Q+00,
     5      -.49345590232374692503112160498434Q+00,
     6      0.12846122696083339472822324849799Q+01,
     7      -.14633068665309633867546120723514Q+00,
     8      0.17866185916994398883966386846435Q+00,
     9      -.46141852214686758823109400137326Q+00,
     a      0.12879036607595322907091315311371Q+01,
     1      0.86900243159408782580295514297060Q-01,
     2      -.83131116369903696639955182756150Q-01,
     3      0.15315035330266609913684702675171Q+00,
     4      -.45962391855374544297312679613435Q+00,
     5      0.12901251573320572229011479302981Q+01,
     6      -.57486081591059187129861037049021Q-01/
         data c1o1002/
     1      0.48472571566865523396832947123549Q-01,
     2      -.64224360552210897202610403496682Q-01,
     3      0.14943195231432834530875984231259Q+00,
     4      -.46002316261937310892316876776698Q+00,
     5      0.12913568812473500077676185772970Q+01,
     6      0.40827703065125946075232227203317Q-01,
     7      -.32265678717150179177917618602566Q-01,
     8      0.34525375044894838970294342394787Q-01,
     9      -.60332272322314282772740405659883Q-01,
     a      0.14874920063873143412936020480455Q+00,
     1      -.46049326723620483405617085225699Q+00,
     2      0.12920933309963045070173874635345Q+01,
     3      -.30488290569622505865837931782581Q-01,
     4      0.23247429889703267195977371692081Q-01,
     5      -.21796680647520970419083068998783Q-01,
     6      0.31026066293900842501856112790833Q-01/
         data c1o1003/
     1      -.59242503492829275854440178959362Q-01,
     2      0.14866542678028002125349594105053Q+00,
     3      -.46085720827356528386514107738634Q+00,
     4      0.12925641717544856067346667450715Q+01,
     5      0.23632748813915624876751363899396Q-01,
     6      -.17642839291090817819225748925710Q-01,
     7      0.15210029691090045701829040266084Q-01,
     8      -.18802559237481291423995941692527Q-01,
     9      0.29848543047196743971707277986449Q-01,
     a      -.58889275641420946894147131123285Q-01,
     1      0.14872326569310690523102917175798Q+00,
     2      -.46112402101074063856034481013610Q+00,
     3      0.12928820448142101457786365389326Q+01,
     4      -.18854711693531909656328814485782Q-01,
     5      0.13890735510381523190405393834687Q-01,
     6      -.11333322107435648672958955027403Q-01/
         data c1o1004/
     1      0.12688351681484053962969897183095Q-01,
     2      -.17670846110161020018964422357005Q-01,
     3      0.29382936880961771474998935688485Q-01,
     4      -.58778762624568960046819773017143Q-01,
     5      0.14880781935180175365902510829816Q+00,
     6      -.46132025329016070260002688451027Q+00,
     7      0.12931062270838588346802241583833Q+01,
     8      0.15392104569373400744195861084447Q-01,
     9      -.11241723173197596011782598607457Q-01,
     a      0.88358844591399214932988207475600Q-02,
     1      -.92143232615951502427101276926789Q-02,
     2      0.11652413523449924149779133371688Q-01,
     3      -.17172539768946889735918297820367Q-01,
     4      0.29186370859620579531326047791951Q-01,
     5      -.58755829440175744861679316097025Q-01,
     6      0.14888849418776972438303544032648Q+00/
         data c1o1005/
     1      -.46146701360013271976057720832267Q+00,
     2      0.12932700279883383094867412133187Q+01,
     3      -.12802811065666065715346047917044Q-01,
     4      0.92954170137557795506903514900896Q-02,
     5      -.71182566681433978004296909984720Q-02,
     6      0.70487150014798276109417242272329Q-02,
     7      -.82869649997843863903016721178282Q-02,
     8      0.11161278420647200770378909007719Q-01,
     9      -.16935777134203964788231232348664Q-01,
     a      0.29102458717751803244036561299490Q-01,
     1      -.58765372338320056652818295529731Q-01,
     2      0.14895810911555491883423773070379Q+00,
     3      -.46157892322255976885588915115041Q+00,
     4      0.12933932439478359165424730745473Q+01,
     5      0.10816008095396386806653865587599Q-01,
     6      -.78202140351643504125089424255580Q-02/
         data c1o1006/
     1      0.58778931033888432070639957138536Q-02,
     2      -.56014804795477159370495098780001Q-02,
     3      0.62267647402601468123057945495256Q-02,
     4      -.78219973552459354876425505495775Q-02,
     5      0.10910342425990932812698798282775Q-01,
     6      -.16818107269944553794266489878664Q-01,
     7      0.29069033597947210709487937967733Q-01,
     8      -.58786187531876037554901086210240Q-01,
     9      0.14901628600021955770673404267710Q+00,
     a      -.46166587737986496349118146382612Q+00,
     1      0.12934882084002606657029027824376Q+01,
     2      -.92582745150575706384557756394980Q-02,
     3      0.66737359117291231521686049863703Q-02,
     4      -.49480344506799801794237971475966Q-02,
     5      0.45811885889074860268799426622751Q-02,
     6      -.48755435853941173472814123562336Q-02/
         data c1o1007/
     1      0.57959050387174571454880225607443Q-02,
     2      -.75716051803491098347305015433445Q-02,
     3      0.10775492741482935289911199215446Q-01,
     4      -.16758602009782952332844702079005Q-01,
     5      0.29059147906832214453427910600376Q-01,
     6      -.58809770206328635943737077618095Q-01,
     7      0.14906447519949251580413667365228Q+00,
     8      -.46173461733859659947548919429881Q+00,
     9      0.12935629172591889998348976687990Q+01,
     a      0.80143842050048650711724382131740Q-02,
     1      -.57641498595593936329786768992575Q-02,
     2      0.42301292085414734166920233948710Q-02,
     3      -.38310938379313195829520389596313Q-02,
     4      0.39401020641291065914119393733230Q-02,
     5      -.44809259047295154062215939368991Q-02,
     6      0.55541415104374135911299387837280Q-02/
         data c1o1008/
     1      -.74298294957891865701181786275426Q-02,
     2      0.10700636628522927310909065015407Q-01,
     3      -.16728966433073622895083737070546Q-01,
     4      0.29060363971427797145746118850268Q-01,
     5      -.58832715621876670343352865708455Q-01,
     6      0.14910441229839364951214243733449Q+00,
     7      -.46178980940592999530079679901988Q+00,
     8      0.12936227329998017120729761594514Q+01,
     9      -.70053501283908085037319816727807Q-02,
     a      0.50298779910368252004419321561577Q-02,
     1      -.36625775502967237754446478294082Q-02,
     2      0.32609000209519836481701335905620Q-02,
     3      -.32638891894889816991000814504265Q-02,
     4      0.35808991509346592602560871566689Q-02,
     5      -.42519479257817470630516595429668Q-02,
     6      0.54117569453388821779372042127799Q-02/
         data c1o1009/
     1      -.73466369632996073357828102396412Q-02,
     2      0.10658385434285701902806442975189Q-01,
     3      -.16715212212270765761241725249132Q-01,
     4      0.29066680553733139425160140739239Q-01,
     5      -.58853759813155809072635783954069Q-01,
     6      0.14913766284793801607801460462024Q+00,
     7      -.46183474393170427111774718630417Q+00,
     8      0.12936713581927511166998825812208Q+01,
     9      0.61755591859859908595775522095141Q-02,
     a      -.44282820295185156878067551963071Q-02,
     1      0.32050752050490513883729722170660Q-02,
     2      -.28155712028632718083133870543830Q-02,
     3      0.27575912860523635538849553156717Q-02,
     4      -.29378605358047356445701915313324Q-02,
     5      0.33665410599666670259856310480093Q-02,
     6      -.41127710812272587557753446921823Q-02/
         data c1o1010/
     1      0.53248497384696121851795337400694Q-02,
     2      -.72966070631232439767624131948707Q-02,
     3      0.10634546253450558510251219573590Q-01,
     4      -.16710080281771355203654067817526Q-01,
     5      0.29075118996675204565925592635016Q-01,
     6      -.58872554199004126916463474007399Q-01,
     7      0.14916552227618912955490399843953Q+00,
     8      -.46187178456083366314190129081792Q+00,
     9      0.12937114146513052306651026216167Q+01,
     a      -.54849370337334146842038575793716Q-02,
     1      0.39290290544600189144537995328779Q-02,
     2      -.28302360560637625411811987438370Q-02,
     3      0.24599464156980519517416234910076Q-02,
     4      -.23673953305832498650703849056476Q-02,
     5      0.24619211190623644750721674619290Q-02,
     6      -.27385992642460616310910449281034Q-02/
         data c1o1011/
     1      0.32327878233310865343846539787791Q-02,
     2      -.40251819017052494989081032584608Q-02,
     3      0.52703730492684356880442209941933Q-02,
     4      -.72660760860126248028828209444335Q-02,
     5      0.10621408257948001469361866360043Q-01,
     6      -.16709659203474534788277400164086Q-01,
     7      0.29084183665204819034100231540387Q-01,
     8      -.58889135693391450905341477188667Q-01,
     9      0.14918902797574665827410644821904Q+00,
     a      -.46190265956101615429027627607282Q+00,
     1      0.12937448003261036973376736949523Q+01,
     2      0.49040145271792943251301803334605Q-02,
     3      -.35100358899168390899535769588158Q-02,
     4      0.25188541696089905728221753270865Q-02,
     5      -.21706523992942817349293805836881Q-02,
     6      0.20593574515425075824948260699438Q-02/
         data c1o1012/
     1      -.20991383501056249787512227905652Q-02,
     2      0.22774600257444464466012505583482Q-02,
     3      -.26114601159932563056275256744091Q-02,
     4      0.31464788047681185367567803685450Q-02,
     5      -.39685622494313510762232697504173Q-02,
     6      0.52355612279051211308567654021219Q-02,
     7      -.72473609021093696524683834008433Q-02,
     8      0.10614619025638981899890648508765Q-01,
     9      -.16711780310726972730608505961807Q-01,
     a      0.29093130742388503120611837293296Q-01,
     1      -.58903691509732109698194719048714Q-01,
     2      0.14920900093621803927820381113601Q+00,
     3      -.46192865362450358643689727014156Q+00,
     4      0.12937729161443235254117272550771Q+01,
     5      -.44107295720252855023914252797274Q-02,
     6      0.31549006462445079696915757259166Q-02/
         data c1o1013/
     1      -.22570953184817624399192551258877Q-02,
     2      0.19316173973121642375562325985822Q-02,
     3      -.18112143564211923648811912108394Q-02,
     4      0.18156670373272067095664543283807Q-02,
     5      -.19287608486843548552286592350936Q-02,
     6      0.21574719651212676866592610731460Q-02,
     7      -.25276675638583843000316550552388Q-02,
     8      0.30893029140489620217141931204709Q-02,
     9      -.39311984067992128899238436681924Q-02,
     a      0.52130330402526085741559270962222Q-02,
     1      -.72359783796958326979510024078865Q-02,
     2      0.10611644206570040947990551538425Q-01,
     3      -.16715216528076568264097922265183Q-01,
     4      0.29101606135504974100117947114817Q-01,
     5      -.58916454103487095076421664526719Q-01,
     6      0.14922608942135356769819433588783Q+00/
         data c1o1014/
     1      -.46195073633666493743380964903804Q+00,
     2      0.12937968141751088572399801324874Q+01,
     3      0.39883001934368117686773270329948Q-02,
     4      -.28512260106466292139113567220281Q-02,
     5      0.20347639948391806395390367456033Q-02,
     6      -.17314647686413567211600560476865Q-02,
     7      0.16078645245723461118583784817259Q-02,
     8      -.15894503837002013863541482409589Q-02,
     9      0.16584499834573420668602833613341Q-02,
     a      -.18160493471195638232364129780714Q-02,
     1      0.20769393457157498756019053653242Q-02,
     2      -.24710166359850882937942342355098Q-02,
     3      0.30506287734535578909828971385521Q-02,
     4      -.39061561859939639662825955967685Q-02,
     5      0.51983666313557966958665436488138Q-02,
     6      -.72292307456098919443602520201375Q-02/
         data c1o1015/
     1      0.10610970451163737969765696703812Q-01,
     2      -.16719265130755265598850835065168Q-01,
     3      0.29109460084380785181088362470080Q-01,
     4      -.58927655266057045129378904595848Q-01,
     5      0.14924080657753637377915605592379Q+00,
     6      -.46196964974158965987147835545310Q+00,
     7      0.12938172967820013735767884147218Q+01,
     8      -.36237775288245953332920139211716Q-02,
     9      0.25894954893173436856677063787653Q-02,
     a      -.18441963690376565672922958509851Q-02,
     1      0.15619367184059793706818921585594Q-02,
     2      -.14387606345225443262967555888429Q-02,
     3      0.14056352327534274121087237902670Q-02,
     4      -.14443915081530306242533783706037Q-02,
     5      0.15528868684697102350824455713892Q-02,
     6      -.17391862592717003530329922621630Q-02/
         data c1o1016/
     1      0.20215363637184751115530320533289Q-02,
     2      -.24319208619980346919930969207420Q-02,
     3      0.30240342107263905230513625199216Q-02,
     4      -.38891891285476432406113173678830Q-02,
     5      0.51888329258980587259443726168704Q-02,
     6      -.72254521071497474115738436186479Q-02,
     7      0.10611675814987330193789908340385Q-01,
     8      -.16723522576514348873122288441486Q-01,
     9      0.29116649756464386032579198695540Q-01,
     a      -.58937507527291613299979400696570Q-01,
     1      0.14925356058568507742129153014046Q+00,
     2      -.46198596904894631828062915984160Q+00,
     3      0.12938349844289873611471585146341Q+01,
     4      0.33070412115031400309045763196897Q-02,
     5      -.23623026642958890942343719605859Q-02,
     6      0.16795315019797775567653628815454Q-02/
         data c1o1017/
     1      -.14169069298070273251203808388062Q-02,
     2      0.12963453133127965917443026860181Q-02,
     3      -.12539292061076800754741262962392Q-02,
     4      0.12717277140247812806228400207507Q-02,
     5      -.13456980131150788638007948833192Q-02,
     6      0.14798781501681814975485612397585Q-02,
     7      -.16854999228236171447185996056643Q-02,
     8      0.19826496211963465529538197966059Q-02,
     9      -.24044863016688517242717022272626Q-02,
     a      0.30055112004860560064490091132330Q-02,
     1      -.38776212411949672871507248877853Q-02,
     2      0.51827043835955376738506008954552Q-02,
     3      -.72235912261994597656239168507769Q-02,
     4      0.10613190018076643740263550751732Q-01,
     5      -.16727759280117925008871771912020Q-01,
     6      0.29123187108973321352580022719303Q-01/
         data c1o1018/
     1      -.58946198155400836622268053172042Q-01,
     2      0.14926467811473300206006621405298Q+00,
     3      -.46200014540374766544634932139130Q+00,
     4      0.12938503629617535392106786355903Q+01,
     5      -.30300890538368838364996003890516Q-02,
     6      0.21638116660818334611398040161979Q-02,
     7      -.15362187722011275776362697550130Q-02,
     8      0.12917426306490516726218293047469Q-02,
     9      -.11750802210832796249190320790628Q-02,
     a      0.11270184295645840050195973811013Q-02,
     1      -.11302037227755660817990938780257Q-02,
     2      0.11795415438787115763299513632121Q-02,
     3      -.12765788872832520074853894784788Q-02,
     4      0.14281968169440493649599816144041Q-02,
     5      -.16472635289262294795711884663122Q-02,
     6      0.19549010024992354167177554698027Q-02/
         data c1o1019/
     1      -.23849732622308134618747524551828Q-02,
     2      0.29924886246748791371210064549106Q-02,
     3      -.38697243668552445022765193585288Q-02,
     4      0.51788634665317137820188937925371Q-02,
     5      -.72229726329319226217250755699078Q-02,
     6      0.10615156575622950247445440558841Q-01,
     7      -.16731848181511751498811489186130Q-01,
     8      0.29129110672718365356773783519338Q-01,
     9      -.58953888821401200281632994567222Q-01,
     a      0.14927442233793829976404213571206Q+00,
     1      -.46201253647809594149422491396936Q+00,
     2      0.12938638171619488685476203405141Q+01/
         real *16
     1  c1n1001(16),c1n1002( 9)
         real *16
     1  c1n2001(2)
         equivalence
     1  (c1n1001(16),c1n2001(1)),(c1n1002(1),c1n2001(2))
         data c1n1001/
     1      0.30618621784789726227466050933818Q+00,
     2      0.51638713464146160635172271753002Q+00,
     3      0.51774675527400439355164307846490Q+00,
     4      0.51649685096288276403108817991942Q+00,
     5      0.51590993537539871811322317538614Q+00,
     6      0.51564648222552252757468483874929Q+00,
     7      0.51550543128756880467920191494527Q+00,
     8      0.51542129573136982440539126646301Q+00,
     9      0.51536713518636286545298246411516Q+00,
     a      0.51533023916933444370386606238850Q+00,
     1      0.51530398125266153391523899016724Q+00,
     2      0.51528463319077700017180409672233Q+00,
     3      0.51526996736151384980396266977682Q+00,
     4      0.51525858640626330101693800484386Q+00,
     5      0.51524957783428522483015915716504Q+00,
     6      0.51524232531411918157135907496065Q+00/
         data c1n1002/
     1      0.51523640039755925759918108058719Q+00,
     2      0.51523149775946248213193378059908Q+00,
     3      0.51522739505758907151662072620597Q+00,
     4      0.51522392726939381402670145192214Q+00,
     5      0.51522096983237791086107548236914Q+00,
     6      0.51521842729851615787551193841352Q+00,
     7      0.51521622553407149223869725314174Q+00,
     8      0.51521430625293056072979572935190Q+00,
     9      0.51521262311844071230517040492421Q+00/
         real *16
     1  c2o1001(16),c2o1002(16),c2o1003(16),c2o1004(16),
     2  c2o1005(16),c2o1006(16),c2o1007(16),c2o1008(16),
     3  c2o1009(16),c2o1010(16),c2o1011(16),c2o1012(16),
     4  c2o1013(16),c2o1014(16),c2o1015(16),c2o1016(16)
         real *16
     1  c2o1017(16),c2o1018(16),c2o1019(12)
         real *16
     1  c2o2001(2),c2o2002(2),c2o2003(2),c2o2004(2),
     2  c2o2005(2),c2o2006(2),c2o2007(2),c2o2008(2),
     3  c2o2009(2),c2o2010(2),c2o2011(2),c2o2012(2),
     4  c2o2013(2),c2o2014(2),c2o2015(2),c2o2016(2)
         real *16
     1  c2o2017(2),c2o2018(2)
         equivalence
     1  (c2o1001(16),c2o2001(1)),(c2o1002(1),c2o2001(2)),
     2  (c2o1002(16),c2o2002(1)),(c2o1003(1),c2o2002(2)),
     3  (c2o1003(16),c2o2003(1)),(c2o1004(1),c2o2003(2)),
     4  (c2o1004(16),c2o2004(1)),(c2o1005(1),c2o2004(2)),
     5  (c2o1005(16),c2o2005(1)),(c2o1006(1),c2o2005(2)),
     6  (c2o1006(16),c2o2006(1)),(c2o1007(1),c2o2006(2)),
     7  (c2o1007(16),c2o2007(1)),(c2o1008(1),c2o2007(2)),
     8  (c2o1008(16),c2o2008(1)),(c2o1009(1),c2o2008(2)),
     9  (c2o1009(16),c2o2009(1)),(c2o1010(1),c2o2009(2)),
     a  (c2o1010(16),c2o2010(1)),(c2o1011(1),c2o2010(2))
         equivalence
     1  (c2o1011(16),c2o2011(1)),(c2o1012(1),c2o2011(2)),
     2  (c2o1012(16),c2o2012(1)),(c2o1013(1),c2o2012(2)),
     3  (c2o1013(16),c2o2013(1)),(c2o1014(1),c2o2013(2)),
     4  (c2o1014(16),c2o2014(1)),(c2o1015(1),c2o2014(2)),
     5  (c2o1015(16),c2o2015(1)),(c2o1016(1),c2o2015(2)),
     6  (c2o1016(16),c2o2016(1)),(c2o1017(1),c2o2016(2)),
     7  (c2o1017(16),c2o2017(1)),(c2o1018(1),c2o2017(2)),
     8  (c2o1018(16),c2o2018(1)),(c2o1019(1),c2o2018(2))
         data c2o1001/
     1      0.14965986394557823129251700680272Q+01,
     2      -.70229042437678827554907437078451Q+00,
     3      0.12877499134644976996127047621891Q+01,
     4      0.28398804276438171918721401096329Q+00,
     5      -.47474352997359069550751784323569Q+00,
     6      0.12840324854247836153298597175608Q+01,
     7      -.14330809276578105731084238677209Q+00,
     8      0.16627197908402170577747056266899Q+00,
     9      -.45971365333809431484756394910271Q+00,
     a      0.12884454235579846591754842478211Q+01,
     1      0.86225978969959947309440333389520Q-01,
     2      -.74756418771357353371504192037049Q-01,
     3      0.15120288810139867760212908473314Q+00,
     4      -.45960369924600842561534754373252Q+00,
     5      0.12904569491240628076107206455390Q+01,
     6      -.57733322420031560232578174884221Q-01/
         data c2o1002/
     1      0.42628297361309000774964456693428Q-01,
     2      -.62474304800739393022017143949499Q-01,
     3      0.14910244112371045631498246030700Q+00,
     4      -.46012332828865092131492316160638Q+00,
     5      0.12915598951652865364003403946228Q+01,
     6      0.41414040491829972507802011411482Q-01,
     7      -.28037186086003273211356644178808Q-01,
     8      0.33065184176052492967550466344049Q-01,
     9      -.59905710731172997940363220026111Q-01,
     a      0.14868689566730452067662306711527Q+00,
     1      -.46058515611352708049741113665098Q+00,
     2      0.12922234790671084292004710304299Q+01,
     3      -.31176355965640060295958137283993Q-01,
     4      0.20083014870954604571762970095499Q-01,
     5      -.20603582841021072147982432486618Q-01,
     6      0.30597269004098035940697251540037Q-01/
         data c2o1003/
     1      -.59105600864005479549855801127874Q-01,
     2      0.14866805031990386562719322338680Q+00,
     3      -.46092762136696688103529499541594Q+00,
     4      0.12926517456348597240171415478906Q+01,
     5      0.24324760611126751852162653313335Q-01,
     6      -.15203868821546114893192513951578Q-01,
     7      0.14236413910283358494613665016655Q-01,
     8      -.18407077292157821302593214716261Q-01,
     9      0.29684573196300267453452089704218Q-01,
     a      -.58844116797925612880527102926852Q-01,
     1      0.14874117307377243987800133723700Q+00,
     2      -.46117655573520446110302991889458Q+00,
     3      0.12929434925996950215842509300244Q+01,
     4      -.19511434606162224343716465326308Q-01,
     5      0.11962905479941608294187906063029Q-01,
     6      -.10533644389961121770321551172981Q-01/
         data c2o1004/
     1      0.12336420786112605159364688628027Q-01,
     2      -.17503003264797128745049575250451Q-01,
     3      0.29314615709296288250070070543395Q-01,
     4      -.58766151001797617001894349270731Q-01,
     5      0.14882747858703819315680433472669Q+00,
     6      -.46135970417309459028306083676074Q+00,
     7      0.12931508813055606960141302321641Q+01,
     8      0.15999665931787718980711145387158Q-01,
     9      -.96848645098491169053782200181507Q-02,
     a      0.81728402619602700274797331640402Q-02,
     1      -.89058515163116359822713885670094Q-02,
     2      0.11491576512638789450921650795272Q-01,
     3      -.17094545263269818063331687562319Q-01,
     4      0.29156916766184038486211468104761Q-01,
     5      -.58755416297276653446158053247335Q-01,
     6      0.14890624634073718881782947869812Q+00/
         data c2o1005/
     1      -.46149709494579332431632505411038Q+00,
     2      0.12933034428247231247942894272678Q+01,
     3      -.13358469547156277221571774248442Q-01,
     4      0.80149165841215493744379746145666Q-02,
     5      -.65628091811387660395381816718021Q-02,
     6      0.67797724275315564042848238283006Q-02,
     7      -.81377635305340687934928000742426Q-02,
     8      0.11081039990789551875430965545651Q-01,
     9      -.16897434008528418752227279758949Q-01,
     a      0.29090116454292730068564915991007Q-01,
     1      -.58769461125755367761186966730949Q-01,
     2      0.14897322747502550798004756125906Q+00,
     3      -.46160225443711234150318947957875Q+00,
     4      0.12934188724138968064479355547337Q+01,
     5      0.11321808268568420687792556056180Q-01,
     6      -.67503448436359933252095069301725Q-02/
         data c2o1006/
     1      0.54078047471773538704829924790208Q-02,
     2      -.53670763044926674128007388751509Q-02,
     3      0.60906698628672466834296777217476Q-02,
     4      -.77435380706342413320745157734511Q-02,
     5      0.10867876818908473443479244404326Q-01,
     6      -.16798715337317507255170504342057Q-01,
     7      0.29064613422648122848071891957693Q-01,
     8      -.58791695740386209185482324585746Q-01,
     9      0.14902891576740716768023578109167Q+00,
     a      -.46168427465434858866239707307995Q+00,
     1      0.12935082810529306677012173406404Q+01,
     2      -.97181434583997863249192097350848Q-02,
     3      0.57676440616223033843034503531381Q-02,
     4      -.45463045629587083915891746375945Q-02,
     5      0.43764129415439955524536112414173Q-02,
     6      -.47524702830353384396130359304630Q-02/
         data c2o1007/
     1      0.57212957553834912319595037503687Q-02,
     2      -.75278967405090686916242908421433Q-02,
     3      0.10752124523784556264431867850631Q-01,
     4      -.16748799697519734483872233358242Q-01,
     5      0.29058482921702236983288185225066Q-01,
     6      -.58815438790711145591086066744286Q-01,
     7      0.14907498098331687007089644701357Q+00,
     8      -.46174934741669158389105133450654Q+00,
     9      0.12935789237616824530324842382155Q+01,
     a      0.84327920308872168481850370890071Q-02,
     1      -.49876578023222894852747711969532Q-02,
     2      0.38836973918412932181124833913972Q-02,
     3      -.36515387578814447847681362268725Q-02,
     4      0.38292448231202727207346389089669Q-02,
     5      -.44111155673899015512492901387148Q-02,
     6      0.55109095690845044640643180129765Q-02/
         data c2o1008/
     1      -.74044657340421331431309348678627Q-02,
     2      0.10687480905358419757165823928543Q-01,
     3      -.16724216446602412750939193280836Q-01,
     4      0.29061465895676255311556880027818Q-01,
     5      -.58838046375542002694638221941623Q-01,
     6      0.14911317020943808141406308045731Q+00,
     7      -.46180176798286267251532150484320Q+00,
     8      0.12936356974488529177877094869873Q+01,
     9      -.73867029764517731688135505640030Q-02,
     a      0.43575388040062873471744280221054Q-02,
     1      -.33613304818634355847534916177098Q-02,
     2      0.31027681308168322895177451712729Q-02,
     3      -.31641555567503012378162194992831Q-02,
     4      0.35161977446180268940188651793497Q-02,
     5      -.42101875165142041710998940300899Q-02,
     6      0.53856784004379838118233218843888Q-02/
         data c2o1009/
     1      -.73314910091715006030253624501628Q-02,
     2      0.10650928148177560074578991080605Q-01,
     3      -.16713195068972082089227846243974Q-01,
     4      0.29068567194593504005036141418479Q-01,
     5      -.58858581850570787493050848744028Q-01,
     6      0.14914500053743694109043207913948Q+00,
     7      -.46184457461640577971236988067845Q+00,
     8      0.12936820026672948298580068127392Q+01,
     9      0.65239418203291713990708859346573Q-02,
     a      -.38407997356182127071916433502209Q-02,
     1      0.29411086888913972358029855926615Q-02,
     2      -.26756555071396351605609005410069Q-02,
     3      0.26678227348060047855576201295842Q-02,
     4      -.28782223527946498194957610853009Q-02,
     5      0.33267932981548808220700084084060Q-02,
     6      -.40867984340772159246410897406531Q-02/
         data c2o1010/
     1      0.53086410490755178469410291340408Q-02,
     2      -.72873941468026353650845212532000Q-02,
     3      0.10630372793613105791729670552710Q-01,
     4      -.16709560074434071864189153756518Q-01,
     5      0.29077295854703347220663248406350Q-01,
     6      -.58876837944944503056657720463936Q-01,
     7      0.14917170866266866619027283470380Q+00,
     8      -.46187995745715869058452929673320Q+00,
     9      0.12937202597323223720008823032978Q+01,
     a      -.58040230984115510622780737704786Q-02,
     1      0.34115306836433394086996959787982Q-02,
     2      -.25973122080254173492386451930118Q-02,
     3      0.23355633178497597553759666264981Q-02,
     4      -.22864751674148336749558260127510Q-02,
     5      0.24071112902201242798188840403872Q-02,
     6      -.27011217080577385549939750900177Q-02/
         data c2o1011/
     1      0.32074359461532113346679933522916Q-02,
     2      -.40085420364200610543571401787310Q-02,
     3      0.52600750362655603133575714858575Q-02,
     4      -.72604235056795631374356096293524Q-02,
     5      0.10619172472479876840423492072436Q-01,
     6      -.16709956056784869511714832318802Q-01,
     7      0.29086401532575930552080168801509Q-01,
     8      -.58892910316308492125688608684777Q-01,
     9      0.14919427864780196356749754484048Q+00,
     a      -.46190952349882236885085199182027Q+00,
     1      0.12937522287803963840899636951376Q+01,
     2      0.51970559593618589303706693263601Q-02,
     3      -.30508881950428307338440840084232Q-02,
     4      0.23120023068765516787743157299550Q-02,
     5      -.20595638455398637427111836924438Q-02,
     6      0.19862589065402186353200169229973Q-02/
         data c2o1012/
     1      -.20488304251626505932889836378905Q-02,
     2      0.22423349193833991170271630005344Q-02,
     3      -.25870398132793765781242059353865Q-02,
     4      0.31298347026340533138336555872966Q-02,
     5      -.39576548808363978823507719025424Q-02,
     6      0.52289185898714071894976270958615Q-02,
     7      -.72439013304321297183020238909538Q-02,
     8      0.10613544782115917688679038560957Q-01,
     9      -.16712511320284569005452542517316Q-01,
     a      0.29095267106502099477484468418112Q-01,
     1      -.58907007652352656821391571111534Q-01,
     2      0.14921348732449195499997198674375Q+00,
     3      -.46193447131008125745102572120958Q+00,
     4      0.12937792143762883588454518487947Q+01,
     5      -.46805815173861784318192750259920Q-02,
     6      0.27448833211161273519294301136574Q-02/
         data c2o1013/
     1      -.20723172001067995180379904242411Q-02,
     2      0.18319569504658486657789131118053Q-02,
     3      -.17450174017230408058098532819586Q-02,
     4      0.17694995651476086174457999556062Q-02,
     5      -.18959641827017376318870740496861Q-02,
     6      0.21341579780855140855323381175344Q-02,
     7      -.25113034617963820868720357011913Q-02,
     8      0.30781235834857027699080464205596Q-02,
     9      -.39239230773771836533263233090360Q-02,
     a      0.52087111202847593843258824259454Q-02,
     1      -.72338974266221333783128679575941Q-02,
     2      0.10611271632243499957088487120373Q-01,
     3      -.16716162670479393994637634164341Q-01,
     4      0.29103605001159843067825417965408Q-01,
     5      -.58919367142457050239759621356890Q-01,
     6      0.14922994773506821315603758042791Q+00/
         data c2o1014/
     1      -.46195570835955052691088729796873Q+00,
     2      0.12938021998676135302706911149348Q+01,
     3      0.42374598334073587356200466805523Q-02,
     4      -.24829433115607629149769405319572Q-02,
     5      0.18688146074588515413117203111945Q-02,
     6      -.16416728549009706041796194651504Q-02,
     7      0.15477568180945156542173807658291Q-02,
     8      -.15470590231411956144487757474320Q-02,
     9      0.16278955943348948005947352292909Q-02,
     a      -.17939261604078824271952061749878Q-02,
     1      0.20610396019366696486899908655786Q-02,
     2      -.24598032008282677224531863280958Q-02,
     3      0.30429829177162718187814221060100Q-02,
     4      -.39012418292973347365592393398141Q-02,
     5      0.51955500708245058426232990656211Q-02,
     6      -.72280297836767117662623231291597Q-02/
         data c2o1015/
     1      0.10611021021209334372860009181985Q-01,
     2      -.16720300248205612043259622150345Q-01,
     3      0.29111300546588519423721414573593Q-01,
     4      -.58930218085457285608948261039095Q-01,
     5      0.14924414536191332807604117224743Q+00,
     6      -.46197393122393882123238524897032Q+00,
     7      0.12938219377320989577597198051349Q+01,
     8      -.38544262856609612357233750878309Q-02,
     9      0.22569504515964219211290186387384Q-02,
     a      -.16944195494473292998025346350039Q-02,
     1      0.14807060842498711848599534560800Q-02,
     2      -.13840308836230311932764786132091Q-02,
     3      0.13666702555271759862051630702557Q-02,
     4      -.14159591674258993715719237034763Q-02,
     5      0.15319797041379045960455942258963Q-02,
     6      -.17238653080578536829772559430213Q-02/
         data c2o1016/
     1      0.20104545817426113873199323111285Q-02,
     2      -.24240966120753985076698657118604Q-02,
     3      0.30187304576395445083088848102182Q-02,
     4      -.38858427170084567563432264715593Q-02,
     5      0.51870100203731402722768210192794Q-02,
     6      -.72248186880314448257905738670206Q-02,
     7      0.10611977932818573980976308134098Q-01,
     8      -.16724574250618126778356369802992Q-01,
     9      0.29118329228107947155629878094294Q-01,
     a      -.58939767739604352346616504818792Q-01,
     1      0.14925646676616165410816240286523Q+00,
     2      -.46198968136816637540124743737743Q+00,
     3      0.12938390116369393665865621616049Q+01,
     4      0.35210853966404463000986910558252Q-02,
     5      -.20605842319678772089943691120280Q-02,
     6      0.15437355924281575540615863094668Q-02/
         data c2o1017/
     1      -.13431376220872464818213127173355Q-02,
     2      0.12463734870952688116943569141948Q-02,
     3      -.12180643145442938185279738664411Q-02,
     4      0.12452808977922172110269446145341Q-02,
     5      -.13259938490122982096863618087613Q-02,
     6      0.14652016011107465240536658278556Q-02,
     7      -.16746631585746125806167710777724Q-02,
     8      0.19747870848573734349165030790203Q-02,
     9      -.23989462304316303567142136160118Q-02,
     a      0.30017927038285398179715067582539Q-02,
     1      -.38753349727892555707622859914711Q-02,
     2      0.51815465275142030435409882879742Q-02,
     3      -.72233257322817918505252153186898Q-02,
     4      0.10613636641994856663395299581425Q-01,
     5      -.16728786683834524364332124997250Q-01,
     6      0.29124712238489999396630479383969Q-01/
         data c2o1018/
     1      -.58948197381610866180765055982116Q-01,
     2      0.14926722171429158297043917153030Q+00,
     3      -.46200338456694829777333900797538Q+00,
     4      0.12938538798891815523662488881967Q+01,
     5      -.32291968647174516951979716721257Q-02,
     6      0.18888615455915796881311290080971Q-02,
     7      -.14125817157612760142117863588906Q-02,
     8      0.12245052233434792539070480290308Q-02,
     9      -.11293269260895078307299691991766Q-02,
     a      0.10939544639043562373284605152575Q-02,
     1      -.11056007091653259173041562293333Q-02,
     2      0.11610033699515430398161954297122Q-02,
     3      -.12625784525696814131154601771110Q-02,
     4      0.14176803331783680568708936599831Q-02,
     5      -.16394639552485645661573104861140Q-02,
     6      0.19492398773348735085327154278571Q-02/
         data c2o1019/
     1      -.23810042207297770752847213070636Q-02,
     2      0.29898623912435877915123044014807Q-02,
     3      -.38681657203548988569636238598469Q-02,
     4      0.51781554321103723653385784878560Q-02,
     5      -.72229454687159537565236455504832Q-02,
     6      0.10615680420748593996299559310249Q-01,
     7      -.16732829253567073958218022509117Q-01,
     8      0.29130492393430689873174071174514Q-01,
     9      -.58955662887333557701259503259579Q-01,
     a      0.14927666012180255370143741269026Q+00,
     1      -.46201537918750296567598422105125Q+00,
     2      0.12938669063866312128142879993034Q+01/
         real *16
     1  c2n1001(16),c2n1002( 9)
         real *16
     1  c2n2001(2)
         equivalence
     1  (c2n1001(16),c2n2001(1)),(c2n1002(1),c2n2001(2))
         data c2n1001/
     1      0.25877458475338283167354757423513Q+00,
     2      0.51768571340438849046436368683072Q+00,
     3      0.51806928592740971129449783795901Q+00,
     4      0.51633448710822207877824704058422Q+00,
     5      0.51583581635639711552926963847645Q+00,
     6      0.51560664061663314746562657031711Q+00,
     7      0.51548180028369765369860373151703Q+00,
     8      0.51540620348515997308574265607635Q+00,
     9      0.51535693488138964673236727353047Q+00,
     a      0.51532303265390605175135135535039Q+00,
     1      0.51529870577943426033593717485629Q+00,
     2      0.51528065756268329807576444601157Q+00,
     3      0.51526689807859064160344987080904Q+00,
     4      0.51525616804555370568470477936011Q+00,
     5      0.51524763884647756531097082791665Q+00,
     6      0.51524074707126588761131382033091Q+00/
         data c2n1002/
     1      0.51523509877222317451093956871189Q+00,
     2      0.51523041176315417645280817332360Q+00,
     3      0.51522647961148699913680545124870Q+00,
     4      0.51522314847663577895035146235306Q+00,
     5      0.51522030181031514833431507596788Q+00,
     6      0.51521785000972770451403774870932Q+00,
     7      0.51521572327208443874138328571757Q+00,
     8      0.51521386656631144545342999253052Q+00,
     9      0.51521223603409126362336110984173Q+00/
         real *16
     1  c3o1001(16),c3o1002(16),c3o1003(16),c3o1004(16),
     2  c3o1005(16),c3o1006(16),c3o1007(16),c3o1008(16),
     3  c3o1009(16),c3o1010(16),c3o1011(16),c3o1012(16),
     4  c3o1013(16),c3o1014(16),c3o1015(16),c3o1016(16)
         real *16
     1  c3o1017(16),c3o1018(16),c3o1019(12)
         real *16
     1  c3o2001(2),c3o2002(2),c3o2003(2),c3o2004(2),
     2  c3o2005(2),c3o2006(2),c3o2007(2),c3o2008(2),
     3  c3o2009(2),c3o2010(2),c3o2011(2),c3o2012(2),
     4  c3o2013(2),c3o2014(2),c3o2015(2),c3o2016(2)
         real *16
     1  c3o2017(2),c3o2018(2)
         equivalence
     1  (c3o1001(16),c3o2001(1)),(c3o1002(1),c3o2001(2)),
     2  (c3o1002(16),c3o2002(1)),(c3o1003(1),c3o2002(2)),
     3  (c3o1003(16),c3o2003(1)),(c3o1004(1),c3o2003(2)),
     4  (c3o1004(16),c3o2004(1)),(c3o1005(1),c3o2004(2)),
     5  (c3o1005(16),c3o2005(1)),(c3o1006(1),c3o2005(2)),
     6  (c3o1006(16),c3o2006(1)),(c3o1007(1),c3o2006(2)),
     7  (c3o1007(16),c3o2007(1)),(c3o1008(1),c3o2007(2)),
     8  (c3o1008(16),c3o2008(1)),(c3o1009(1),c3o2008(2)),
     9  (c3o1009(16),c3o2009(1)),(c3o1010(1),c3o2009(2)),
     a  (c3o1010(16),c3o2010(1)),(c3o1011(1),c3o2010(2))
         equivalence
     1  (c3o1011(16),c3o2011(1)),(c3o1012(1),c3o2011(2)),
     2  (c3o1012(16),c3o2012(1)),(c3o1013(1),c3o2012(2)),
     3  (c3o1013(16),c3o2013(1)),(c3o1014(1),c3o2013(2)),
     4  (c3o1014(16),c3o2014(1)),(c3o1015(1),c3o2014(2)),
     5  (c3o1015(16),c3o2015(1)),(c3o1016(1),c3o2015(2)),
     6  (c3o1016(16),c3o2016(1)),(c3o1017(1),c3o2016(2)),
     7  (c3o1017(16),c3o2017(1)),(c3o1018(1),c3o2017(2)),
     8  (c3o1018(16),c3o2018(1)),(c3o1019(1),c3o2018(2))
         data c3o1001/
     1      0.18383838383838383838383838383838Q+01,
     2      -.67645950378691790573535563677811Q+00,
     3      0.12667714034321556284601880979144Q+01,
     4      0.27386959273647190091084089038425Q+00,
     5      -.45998124613020651065121024931895Q+00,
     6      0.12852239391136641850223233531961Q+01,
     7      -.13929818750345940229625517252884Q+00,
     8      0.15578227376068860031460192376812Q+00,
     9      -.45934381396079482784278865178286Q+00,
     a      0.12889772965614738737165535495815Q+01,
     1      0.84709423170628182004711078721237Q-01,
     2      -.67336709772711305756876696728615Q-01,
     3      0.15026344245462888666836613851295Q+00,
     4      -.45966841183910822451768135282858Q+00,
     5      0.12907605321919497095346278527735Q+01,
     6      -.57287860177357382725782241747261Q-01/
         data c3o1002/
     1      0.37266809627776809415741846203879Q-01,
     2      -.61467241741089142260946065639480Q-01,
     3      0.14889236640903169581090205240026Q+00,
     4      -.46023141670304239321410099457771Q+00,
     5      0.12917440044270380472103847627099Q+01,
     6      0.41445352238911723360519667631536Q-01,
     7      -.24050282932560346444864388116748Q-01,
     8      0.32154193685117469990994045967781Q-01,
     9      -.59595543097435666378032149720219Q-01,
     a      0.14865329133928899120926748963552Q+00,
     1      -.46067362068733961610193196195246Q+00,
     2      0.12923420582100297156284239706800Q+01,
     3      -.31420655662347412587067424914391Q-01,
     4      0.17033641581379285828274026617091Q-01,
     5      -.19823272061525034254330757152298Q-01,
     6      0.30269425906554564395388917689181Q-01/
         data c3o1003/
     1      -.59003740156607639647383915437365Q-01,
     2      0.14867801542200400195248407013232Q+00,
     3      -.46099383189495850867682385326825Q+00,
     4      0.12927321194519337830490662008686Q+01,
     5      0.24658961337005246234408257603939Q-01,
     6      -.12811828860651724303289887929744Q-01,
     7      0.13579747234221650695856220307314Q-01,
     8      -.18095508090488011816260856185540Q-01,
     9      0.29555126140823145780975793386366Q-01,
     a      -.58811058761242182333595748660118Q-01,
     1      0.14876062801170991565772660369726Q+00,
     2      -.46122573080019607818284892774416Q+00,
     3      0.12930002985106114526394762596097Q+01,
     4      -.19876126152107314668999882051965Q-01,
     5      0.10044684177168526344952048892933Q-01,
     6      -.99826688437933309947914881200144Q-02/
         data c3o1004/
     1      0.12053283890747886106505361678394Q-01,
     2      -.17366521971846202866801116384616Q-01,
     3      0.29259827490218123752401165533608Q-01,
     4      -.58758059901574254640331317367078Q-01,
     5      0.14884702610151841504029278374030Q+00,
     6      -.46139667171024714338160407168416Q+00,
     7      0.12931924332934704368142038155307Q+01,
     8      0.16365782323428255100153128229039Q-01,
     9      -.81170507316503827506699128007369Q-02,
     a      0.77089529067271922466235203797512Q-02,
     1      -.86537263907778351599748771838171Q-02,
     2      0.11358223261692716179357443445551Q-01,
     3      -.17029893093113213189121190059984Q-01,
     4      0.29133287609170416457310805425812Q-01,
     5      -.58756724565088690971048382187489Q-01,
     6      0.14892345339740468400908411297098Q+00/
         data c3o1005/
     1      -.46152536428953268090243278107669Q+00,
     2      0.12933347154204908932660733471786Q+01,
     3      -.13711997763627086727701142441091Q-01,
     4      0.67122796102005082000877925921628Q-02,
     5      -.61697966665420893012048436809981Q-02,
     6      0.65572170130222029339190888573099Q-02,
     7      -.80122653446114644700367844966387Q-02,
     8      0.11013226072089601162298242809318Q-01,
     9      -.16865308670946627376351377797702Q-01,
     a      0.29080439762602249094798335604220Q-01,
     1      -.58774160423373152748873021354820Q-01,
     2      0.14898774821774304718642021003883Q+00,
     3      -.46162425480708852656575360937530Q+00,
     4      0.12934429782765029466192160161524Q+01,
     5      0.11656416780121992309510314238572Q-01,
     6      -.56525298096806021208381047122312Q-02/
         data c3o1006/
     1      0.50723770869475619730756681488851Q-02,
     2      -.51711622219780612155089349840709Q-02,
     3      0.59748966974714634175804033691693Q-02,
     4      -.76763177266160424343232889524483Q-02,
     5      0.10831536219675312056933557843109Q-01,
     6      -.16782415592247308209347568836550Q-01,
     7      0.29061455848136528739984758620701Q-01,
     8      -.58797357083935183397698012885906Q-01,
     9      0.14904100820863710683716024617566Q+00,
     a      -.46170168168390364987634701161533Q+00,
     1      0.12935272438794040123289791344692Q+01,
     2      -.10031481096851909604023026834432Q-01,
     3      0.48309211610577331805140256110616Q-02,
     4      -.42578262865875461775337970410484Q-02,
     5      0.42038639164178644663201032669441Q-02,
     6      -.46468103533586361703409130698190Q-02/
         data c3o1007/
     1      0.56567004151350238246103820805402Q-02,
     2      -.74899764653423445607128327730690Q-02,
     3      0.10731972317146109834213593189585Q-01,
     4      -.16740608101744201432124820787803Q-01,
     5      0.29058420306217618396552611120921Q-01,
     6      -.58821072868248127904901857823046Q-01,
     7      0.14908503410221282848367039395362Q+00,
     8      -.46176332967365243308182357302610Q+00,
     9      0.12935941032729901592694355802996Q+01,
     a      0.87245689702781993693561275637969Q-02,
     1      -.41796797116294841752410745164536Q-02,
     2      0.36337249372056512468470701644304Q-02,
     3      -.34992178430325741799134561954556Q-02,
     4      0.37333457866982826296550799870905Q-02,
     5      -.43501613070972390249482126762566Q-02,
     6      0.54730214663504554660876746109127Q-02/
         data c3o1008/
     1      -.73822679426464739062906813390279Q-02,
     2      0.10676095540654555514273768382466Q-01,
     3      -.16720332481801276080903961240986Q-01,
     4      0.29062846511930111492711265829167Q-01,
     5      -.58843270961360491297787937464782Q-01,
     6      0.14912155600082302943048198874878Q+00,
     7      -.46181315336544885221835707021904Q+00,
     8      0.12936480336121069157696180492638Q+01,
     9      -.76576717545695945404903962069119Q-02,
     a      0.36539276015689783664338002985306Q-02,
     1      -.31431621029759000577221589845421Q-02,
     2      0.29678619317613145576432249519019Q-02,
     3      -.30773242967461050068383582951648Q-02,
     4      0.34593050768748468008964349564887Q-02,
     5      -.41732948976703192070672222801070Q-02,
     6      0.53626205879391156178672705008558Q-02/
         data c3o1009/
     1      -.73181580374049942482914114958294Q-02,
     2      0.10644481058137272584848791714380Q-01,
     3      -.16711652142658123219647995637085Q-01,
     4      0.29070569724019356345550970215963Q-01,
     5      -.58863276805476168788114189980958Q-01,
     6      0.14915203477405035016287265232842Q+00,
     7      -.46185395969357318918707679135735Q+00,
     8      0.12936921615881741127116967468568Q+01,
     9      0.67753480304830591693522920007420Q-02,
     a      -.32228727347885836898727333407570Q-02,
     1      0.27494054625856444928586193276995Q-02,
     2      -.25557200137130625540515040444440Q-02,
     3      0.25892391148149078812393662297590Q-02,
     4      -.28254676086327623325599474526785Q-02,
     5      0.32914460721686284311411456450836Q-02,
     6      -.40636539658905600578917599917538Q-02/
         data c3o1010/
     1      0.52942167078242741165936219900810Q-02,
     2      -.72792574927805423859730609965085Q-02,
     3      0.10626791451796051341843606321007Q-01,
     4      -.16709300525674944410702344017546Q-01,
     5      0.29079505929135964457728419865258Q-01,
     6      -.58880995605135604606750036615860Q-01,
     7      0.14917764771024522152727527208451Q+00,
     8      -.46188777933182392488141010694677Q+00,
     9      0.12937287238298709691280745296673Q+01,
     a      -.60373077505603244850129654548182Q-02,
     1      0.28647519399924042659900763079803Q-02,
     2      -.24277965203644419881946453739829Q-02,
     3      0.22285095703599176117446786054687Q-02,
     4      -.22153037279943554242016459772332Q-02,
     5      0.23583778362876303297413609925021Q-02,
     6      -.26676064655804328803170334751885Q-02/
         data c3o1011/
     1      0.31847011437773795456619087888854Q-02,
     2      -.39936161937723442239420293546756Q-02,
     3      0.52508685487704232298599815580299Q-02,
     4      -.72554283445816069066668694014606Q-02,
     5      0.10617290357438564342525486986093Q-01,
     6      -.16710394013881158373716424502966Q-01,
     7      0.29088611101211690536920555254974Q-01,
     8      -.58896568691008312974894250531244Q-01,
     9      0.14919932706572872189968171461124Q+00,
     a      -.46191610753551159792870809138639Q+00,
     1      0.12937593542259177633198296568461Q+01,
     2      0.54136954895555896701582390797804Q-02,
     3      -.25637927912997110695848694086399Q-02,
     4      0.21612238092390467224145824556426Q-02,
     5      -.19636204437302365351418425542112Q-02,
     6      0.19217044732450020163378598206219Q-02/
         data c3o1012/
     1      -.20038995997154944410558048632517Q-02,
     2      0.22107716705004630076344749502747Q-02,
     3      -.25650238030158754263126608497409Q-02,
     4      0.31148116611178103682994488178310Q-02,
     5      -.39478218738952486908214504355805Q-02,
     6      0.52229628740205351795774102203710Q-02,
     7      -.72408524654775886849740031820433Q-02,
     8      0.10612684131596114645343280289520Q-01,
     9      -.16713315235409514229498726322930Q-01,
     a      0.29097375240266663068415109791224Q-01,
     1      -.58910220256785066256310814670624Q-01,
     2      0.14921780748514935913281773457737Q+00,
     3      -.46194006329333774239353999230712Q+00,
     4      0.12937852686701771115375538076743Q+01,
     5      -.48820041724614439894352575937593Q-02,
     6      0.23083161793375362479206321130182Q-02/
         data c3o1013/
     1      -.19374729228370414314624450842274Q-02,
     2      0.17456265552454264780134780716644Q-02,
     3      -.16863500150247703411918040208844Q-02,
     4      0.17281050512405811114805208229166Q-02,
     5      -.18663690397220203213788963030232Q-02,
     6      0.21130433249363994993472500345432Q-02,
     7      -.24964569838054602181316441306705Q-02,
     8      0.30679814485577195205122622674806Q-02,
     9      -.39173400013785852130981944474582Q-02,
     a      0.52048316414865592661834622497786Q-02,
     1      -.72320775052022484640385128272384Q-02,
     2      0.10611028370755463258333746054325Q-01,
     3      -.16717142396694021598434874018921Q-01,
     4      0.29105567141408443570330456016516Q-01,
     5      -.58922189548331428270968683548720Q-01,
     6      0.14923366861801645958237534652271Q+00/
         data c3o1014/
     1      -.46196049650506124942191564211761Q+00,
     2      0.12938073870145617545653859042040Q+01,
     3      0.44250047559304678343088310938618Q-02,
     4      -.20895113695151421231985201609584Q-02,
     5      0.17476110186622013770505531986164Q-02,
     6      -.15636902373223771536153462954930Q-02,
     7      0.14943201241313982349690920988890Q-02,
     8      -.15089189335486011715468814995369Q-02,
     9      0.16002214765150296827916957804824Q-02,
     a      -.17738101424378539038627492999116Q-02,
     1      0.20465512044394139099555172520507Q-02,
     2      -.24495781727657108765539287398531Q-02,
     3      0.30360186113398424614292628807154Q-02,
     4      -.38967838690566130691344392745765Q-02,
     5      0.51930239265411572061843186679590Q-02,
     6      -.72269967736602815286462780132100Q-02/
         data c3o1015/
     1      0.10611149226541829155586312931081Q-01,
     2      -.16721346243505443410240160647469Q-01,
     3      0.29113101833600501236209548434454Q-01,
     4      -.58932702243920724865123130665317Q-01,
     5      0.14924736980979529174187675239044Q+00,
     6      -.46197806148395408990057873122581Q+00,
     7      0.12938264154714467736423994619218Q+01,
     8      -.40293264553989753425090309571491Q-02,
     9      0.19006222161627406365131300988618Q-02,
     a      -.15849665855294393774012444926075Q-02,
     1      0.14100006711751630994980138748783Q-02,
     2      -.13352418760941775927735247380642Q-02,
     3      0.13315055997940523001277594811568Q-02,
     4      -.13901222146215751348398169529905Q-02,
     5      0.15129025552818509771970540416688Q-02,
     6      -.17098512476840575292673160157268Q-02/
         data c3o1016/
     1      0.20003064051178388761103604331217Q-02,
     2      -.24169328798961700970191588149978Q-02,
     3      0.30138846743898251651197001505516Q-02,
     4      -.38828030374448145247380556962501Q-02,
     5      0.51853806629060177320640135921916Q-02,
     6      -.72242941227024186424212792355311Q-02,
     7      0.10612325543356194368671658131720Q-01,
     8      -.16725623829757220349431945314735Q-01,
     9      0.29119970204270096679601088605322Q-01,
     a      -.58941959932815537049912133835792Q-01,
     1      0.14925927722192887964044714547258Q+00,
     2      -.46199326823809982618961951663024Q+00,
     3      0.12938429034643518741096531431396Q+01,
     4      0.36844637407240563485170101417771Q-02,
     5      -.17363928459488320943175645075616Q-02,
     6      0.14444649116530621006293276281801Q-02/
         data c3o1017/
     1      -.12788013273707378847496661294213Q-02,
     2      0.12017179947957290978916171182101Q-02,
     3      -.11856091536302434195359268666114Q-02,
     4      0.12211775407236536528271384735349Q-02,
     5      -.13079581417658744449675677582846Q-02,
     6      0.14517320353473624866244715915868Q-02,
     7      -.16647028931810578977931166843183Q-02,
     8      0.19675577014198125860500041631337Q-02,
     9      -.23938572192432584715653905469247Q-02,
     a      0.29983878267543118069944748876200Q-02,
     1      -.38732581895903977242542934314512Q-02,
     2      0.51805191928373327469808873531325Q-02,
     3      -.72231309155083554852221884832073Q-02,
     4      0.10614108572626802412685633799007Q-01,
     5      -.16729804794118564329340437161011Q-01,
     6      0.29126201088203895973984460436532Q-01/
         data c3o1018/
     1      -.58950137840647219190924997553403Q-01,
     2      0.14926968467411069344167566336025Q+00,
     3      -.46200651882731013160896277090778Q+00,
     4      0.12938572835914578963206696325199Q+01,
     5      -.33820694219007675705236533256981Q-02,
     6      0.15926816921038250039030895920005Q-02,
     7      -.13221832928614885448661060017569Q-02,
     8      0.11657653910152867593691217780152Q-02,
     9      -.10883532715768353920790419170364Q-02,
     a      0.10639610175213178227273705662804Q-02,
     1      -.10831186078261624707545972543929Q-02,
     2      0.11439873316797782215595418464289Q-02,
     3      -.12496911122513041239951439374709Q-02,
     4      0.14079833424662415458695575567717Q-02,
     5      -.16322666150053805263690769545226Q-02,
     6      0.19440172388001489205963525328384Q-02/
         data c3o1019/
     1      -.23773489469002927077264147661229Q-02,
     2      0.29874544225661926812403944680221Q-02,
     3      -.38667520634042880973039746158190Q-02,
     4      0.51775361686623191871708167402792Q-02,
     5      -.72229639686162892830826126308169Q-02,
     6      0.10616216836360770987494665694565Q-01,
     7      -.16733797284269460610944098135240Q-01,
     8      0.29131840713380671676151857027393Q-01,
     9      -.58957386135958154162721588147512Q-01,
     a      0.14927882956656159014782047279138Q+00,
     1      -.46201813352679782709665086287276Q+00,
     2      0.12938699001976308196397589574494Q+01/
         real *16
     1  c3n1001(16),c3n1002( 9)
         real *16
     1  c3n2001(2)
         equivalence
     1  (c3n1001(16),c3n2001(1)),(c3n1002(1),c3n2001(2))
         data c3n1001/
     1      0.21347814095749163298377163917748Q+00,
     2      0.53267987661616876667733101622215Q+00,
     3      0.51739402675583304700148753242507Q+00,
     4      0.51619038846808976944095769378383Q+00,
     5      0.51577022397236114683415275166125Q+00,
     6      0.51557120841841984539596174332987Q+00,
     7      0.51546055170819587790169798147708Q+00,
     8      0.51539248308252102618336525602104Q+00,
     9      0.51534757252888772780909220325251Q+00,
     a      0.51531636425157446044269184700452Q+00,
     1      0.51529379064365519613178298390303Q+00,
     2      0.51527693185775585500371130298358Q+00,
     3      0.51526400737227306158734198303803Q+00,
     4      0.51525388056967614645389872329993Q+00,
     5      0.51524579791804329525418733633793Q+00,
     6      0.51523924371306492002919463317422Q+00/
         data c3n1002/
     1      0.51523385530193741329869474762194Q+00,
     2      0.51522937160541520803885199539575Q+00,
     3      0.51522560077686566533801832188383Q+00,
     4      0.51522239927410176729705847251785Q+00,
     5      0.51521965796031566539420950759781Q+00,
     6      0.51521729265863931897221688925687Q+00,
     7      0.51521523760033258368526256882067Q+00,
     8      0.51521344079572315620881331272257Q+00,
     9      0.51521186070883185645840245506287Q+00/
c 
c       construct the orthogonal polynomials of orders 0, 4, 8, ...
c 
        n4=n/4+1
        index=0
        call orteva(z,n4,c0o1001,c0n1001,index,w)
c 
        do 1400 i=1,n4
        j=(i-1)*4+1
        if(j .gt. n) goto 1400
        pols(j)=w(i)
 1400 continue
c 
c       construct the orthogonal polynomials of orders 1, 5, 9, ...
c 
        index=1
        call orteva(z,n4,c1o1001,c1n1001,index,w)
c 
        do 1600 i=1,n4
        j=(i-1)*4+2
        if(j .gt. n) goto 1600
        pols(j)=w(i)
 1600 continue
c 
c       construct the orthogonal polynomials of orders 2, 6, 10, ...
c 
        index=2
        call orteva(z,n4,c2o1001,c2n1001,index,w)
c 
        do 1800 i=1,n4
        j=(i-1)*4+3
        if(j .gt. n) goto 1800
        pols(j)=w(i)
 1800 continue
c 
c       construct the orthogonal polynomials of orders 3, 7, 11, . . .
c 
        index=3
        call orteva(z,n4,c3o1001,c3n1001,index,w)
c 
        do 2000 i=1,n4
        j=(i-1)*4+4
        if(j .gt. n) goto 2000
        pols(j)=w(i)
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine orteva(z,m,coefs,coenorm,index,pols)
        implicit real *16 (a-h,o-z)
        save
        complex *32 z,pols(1),z4
        dimension coefs(1),coenorm(1)
c 
c        start the gram-schmidt process
c 
        z4=z**4
        pols(1)=z**index
        pols(1)=pols(1)*coenorm(1)
c 
c       conduct the ghost of the gram-schmidt process
c 
        icoefs=0
        do 2000 i=2,m
c 
        pols(i)=pols(i-1)*z4
        do 1400 j=1,i-1
        icoefs=icoefs+1
        pols(i)=pols(i)+coefs(icoefs)*pols(j)
 1400 continue
c 
        pols(i)=pols(i)*coenorm(i)
 2000 continue
        return
        end
  
  
  
