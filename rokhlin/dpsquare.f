        implicit real *8 (a-h,o-z)
        complex *16  pol(300000),work(41000)
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
        implicit real *8 (a-h,o-z)
        save
        complex *16 pol(n*4,m)
        dimension t(2000),w(2000),x(4000),y(4000),rea(2),
     1   weights(4000),rlens(4000),corners(2,4)
        complex *16 com,work(m,m),
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
        call dpsquare(zs(i),m,ww)
        do 1600 j=1,m
        pol(i,j)=ww(j)
 1600 continue
 1700 continue
  
         call prin2('after dpsquare, pol(m)=*',pol(1,m-1),8*n)
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
        call prin2('and the matrix of inner products*',work,m*m*2)
c 
        d=0
        done=1
        do 2800 i=1,m
        do 2600 j=1,m
        if(i .eq. j) goto 2600
        dd=work(j,i)*dconjg(work(j,i))
        if(dd .gt. d) d=dd
 2600 continue
        dd=(work(i,i)-done)*dconjg(work(i,i)-done)
        if(dd .gt. d) d=dd
 2800 continue
        call prin2('the non-orthonormality is*',dsqrt(d),1)
 3000 continue
        stop
        end
  
c 
c 
c 
c 
c 
        subroutine squacre(n,zs,weights,t,w,rlens,corners)
        implicit real *8 (a-h,o-z)
        save
        complex *16 zs(1),ima
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
        coe1=dsqrt((x2-x1)**2+(y2-y1)**2) /2
        coe2=dsqrt((x3-x2)**2+(y3-y2)**2) /2
        coe3=dsqrt((x4-x3)**2+(y4-y3)**2) /2
        coe4=dsqrt((x4-x1)**2+(y4-y1)**2) /2
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
        rlens(1)=cdabs(zs(1)-x1-ima*y1)
        do 2200 i=2,n
        rlens(i)=rlens(i-1)+cdabs(zs(i)-zs(i-1))
 2200 continue
c 
        rlens(n+1)=rlens(n)+cdabs(zs(n)-x2-ima*y2)+
     1   +cdabs(zs(n+1)-x2-ima*y2)
        do 2400 i=n+2,2*n
        rlens(i)=rlens(i-1)+cdabs(zs(i)-zs(i-1))
 2400 continue
  
        rlens(2*n+1)=rlens(2*n)+cdabs(zs(2*n)-x3-ima*y3)+
     1   +cdabs(zs(2*n+1)-x3-ima*y3)
        do 2600 i=2*n+2,3*n
        rlens(i)=rlens(i-1)+cdabs(zs(i)-zs(i-1))
 2600 continue
c 
        rlens(3*n+1)=rlens(3*n)+cdabs(zs(3*n)-x4-ima*y4)+
     1   +cdabs(zs(3*n+1)-x4-ima*y4)
        do 2800 i=3*n+2,4*n
        rlens(i)=rlens(i-1)+cdabs(zs(i)-zs(i-1))
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
        implicit real *8 (a-h,o-z)
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
c  wwx,wwx - each must be at least n real *8 locations long
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
        implicit real *8 (a-h,o-z)
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
        implicit real *8 (a-h,o-z)
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
        implicit real *8 (a-h,o-z)
        save
        dimension w(1)
        complex *16 a(1),b(1),prod
c 
c       calculate inner product of two vectors
c 
        prod=0
        do 1200 i=1,n
        prod=prod+dconjg(a(i))*b(i)*w(i)
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine polplot(pol,n,m,rlens,w,iw1,iw2)
        implicit real *8 (a-h,o-z)
        save
        complex *16 pol(n,1),com
        dimension rlen(1),w(1),rea(2),rlens(1)
        equivalence (rea(1),com)
c 
c       plot the real and imaginary parts of the m-th
c       orthogonal polynomial
c 
cccc        call prin2('in polplot, pol=*',pol,n*m*2)
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
        implicit real *8 (a-h,o-z)
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
        implicit real *8 (a-h,o-z)
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
        pi=datan(done)*4
        h=pi/(2*n)
        do 1200 i=1,n
        t=(2*i-1)*h
        ts(n-i+1)=dcos(t)
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
cccc         call prin2('delta=*',delta,1)
        xk=xk+delta
        if(deltold .le. eps) goto 1600
        deltold=dabs(delta)
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
        implicit real *8 (a-h,o-z)
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
        implicit real *8 (a-h,o-z)
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
        dd=dabs(f)
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
        implicit real *8 (a-h,o-z)
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
      IMPLICIT REAL *8 (A-H,O-Z)
        save
      REAL *8 PEXP(1)
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
        IMPLICIT REAL *8 (A-H,O-Z)
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
        implicit real *8 (a-h,o-z)
        save
        dimension w(1)
        complex *16 pol(n,1),out(m,1)
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
  
        subroutine dpsquare(z,n,pols)
        implicit real *8 (a-h,o-z)
        save
        complex *16 pols(1),w(50),z
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
         real *8
     1  c0o1001(16),c0o1002(16),c0o1003(16),c0o1004(16),
     2  c0o1005(16),c0o1006(16),c0o1007(16),c0o1008(16),
     3  c0o1009(16),c0o1010(16),c0o1011(16),c0o1012(16),
     4  c0o1013(16),c0o1014(16),c0o1015(16),c0o1016(16)
         real *8
     1  c0o1017(16),c0o1018(16),c0o1019(12)
         real *8
     1  c0o2001(2),c0o2002(2),c0o2003(2),c0o2004(2),
     2  c0o2005(2),c0o2006(2),c0o2007(2),c0o2008(2),
     3  c0o2009(2),c0o2010(2),c0o2011(2),c0o2012(2),
     4  c0o2013(2),c0o2014(2),c0o2015(2),c0o2016(2)
         real *8
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
     1      0.8000000000000000D+00,-.6016811464722840D+00,
     2      0.1336045772409409D+01,0.2923835030712923D+00,
     3      -.5142493591548285D+00,0.1288421587739613D+01,
     4      -.1484054449426377D+00,0.1923202084871908D+00,
     5      -.4648829513187934D+00,0.1287542189186905D+01,
     6      0.8684706671759446E-01,-.9202554980667745E-01,
     7      0.1563143720735640D+00,-.4598646975220644D+00,
     8      0.1289773848578182D+01,-.5661370127414933E-01/
         data c0o1002/
     1      0.5446336643648059E-01,-.6682757927854084E-01,
     2      0.1499837711653516D+00,-.4599440987185718D+00,
     3      0.1291134362843109D+01,0.3973616501798772E-01,
     4      -.3647061424971176E-01,0.3659196121996946E-01,
     5      -.6095235590711094E-01,0.1488547748243972D+00,
     6      -.4603999933726890D+00,0.1291950575105609D+01,
     7      -.2939900412455895E-01,0.2631587433334332E-01,
     8      -.2343090035726346E-01,0.3161354338424757E-01/
         data c0o1003/
     1      -.5942778675109477E-01,0.1486736870131550D+00,
     2      -.4607827709133309D+00,0.1292468674617253D+01,
     3      0.2262100450817787E-01,-.1995890516234639E-01,
     4      0.1651315117563385E-01,-.1932564574146965E-01,
     5      0.3005828073736836E-01,-.5895040965685448E-01,
     6      0.1487079115132621D+00,-.4610679841834093D+00,
     7      0.1292815496519854D+01,-.1794038496822044E-01,
     8      0.1568968550145759E-01,-.1238559015295970E-01/
         data c0o1004/
     1      0.1314273841577516E-01,-.1787921732688591E-01,
     2      0.2946842965906761E-01,-.5879723101781867E-01,
     3      0.1487883451158418D+00,-.4612781557428777D+00,
     4      0.1293058180024667D+01,0.1457424623039973E-01,
     5      -.1267317701610390E-01,0.9697211679404366E-02,
     6      -.9605594348715804E-02,0.1184812835933813E-01,
     7      -.1726707793310648E-01,0.2922303837226929E-01,
     8      -.5875847442244129E-01,0.1488702743840044D+00/
         data c0o1005/
     1      -.4614349894657617D+00,0.1293234282915613D+01,
     2      -.1207318836764465E-01,0.1045797510876668E-01,
     3      -.7832683641959687E-02,0.7385241712493217E-02,
     4      -.8465838874266360E-02,0.1125667454040458E-01,
     5      -.1698164781186719E-01,0.2911805405029819E-01,
     6      -.5876210160463390E-01,0.1489424003136593D+00,
     7      -.4615541645607588D+00,0.1293365968091276D+01,
     8      0.1016449113021702E-01,-.8780978056109369E-02/
         data c0o1006/
     1      0.6477840933736578E-02,-.5891671411449721E-02,
     2      0.6388031909010958E-02,-.7913999044540485E-02,
     3      0.1096011217215835E-01,-.1684118288971940E-01,
     4      0.2907498382748811E-01,-.5878091920179819E-01,
     5      0.1490031048727955D+00,-.4616464198933047D+00,
     6      0.1293466940935894D+01,-.8674941400840081E-02,
     7      0.7479695100834993E-02,-.5457578066408765E-02,
     8      0.4832538446955918E-02,-.5020008923502064E-02/
         data c0o1007/
     1      0.5882462170954800E-02,-.7622138109752966E-02,
     2      0.1080263606123075E-01,-.1677030201238338E-01,
     3      0.2906054356982075E-01,-.5880410283065851E-01,
     4      0.1490534982000186D+00,-.4617190887189785D+00,
     5      0.1293546024478981D+01,0.7490271503615120E-02,
     6      -.6449051363014305E-02,0.4667355048028919E-02,
     7      -.4049959643460121E-02,0.4069215437124943E-02,
     8      -.4561218051453603E-02,0.5603617182332187E-02/
         data c0o1008/
     1      -.7458870081258318E-02,0.1071584660809909E-01,
     2      -.1673472944471003E-01,0.2905960455624253E-01,
     3      -.5882729209098812E-01,0.1490952647164727D+00,
     4      -.4617772405429499D+00,0.1293609097988543D+01,
     5      -.6532641467597519E-02,0.5618506318519426E-02,
     6      -.4041245880722714E-02,0.3452554827061849E-02,
     7      -.3379288401152892E-02,0.3654781693411233E-02,
     8      -.4299353494745805E-02,0.5441314207035560E-02/
         data c0o1009/
     1      -.7363864799825718E-02,0.1066700487417623E-01,
     2      -.1671778210708086E-01,0.2906494221788915E-01,
     3      -.5884881420332786E-01,0.1491300066907006D+00,
     4      -.4618244401927805D+00,0.1293660197413097D+01,
     5      0.5747548298487204E-02,-.4939203616769739E-02,
     6      0.3535788522627399E-02,-.2984347807931831E-02,
     7      0.2860884385439999E-02,-.3005547421409269E-02,
     8      0.3411360170815887E-02,-.4141978181948708E-02/
         data c0o1010/
     1      0.5343091103734769E-02,-.7307045111375871E-02,
     2      0.1063939664959065E-01,-.1671090448280967E-01,
     3      0.2907299189534687E-01,-.5886814359458188E-01,
     4      0.1491590762096307D+00,-.4618632400732294D+00,
     5      0.1293702165813890D+01,-.5095909837756667E-02,
     6      0.4376399722323893E-02,-.3121264597685908E-02,
     7      0.2609395337643575E-02,-.2460065051334188E-02,
     8      0.2523803275551496E-02,-.2780618854569810E-02/
         data c0o1011/
     1      0.3261115237322800E-02,-.4043761657218166E-02,
     2      0.5281903919089044E-02,-.7272471911837942E-02,
     3      0.1062404668940399E-01,-.1670952794813943E-01,
     4      0.2908196578764820E-01,-.5888524233912450E-01,
     5      0.1491835650716317D+00,-.4618955001094846D+00,
     6      0.1293737051723297D+01,0.4549112008085717E-02,
     7      -.3904820067875288E-02,0.2776724502219166E-02,
     8      -.2303685987722384E-02,0.2142727266514457E-02/
         data c0o1012/
     1      -.2155680551893211E-02,0.2316649519063072E-02,
     2      -.2638600317691127E-02,0.3164947548113847E-02,
     3      -.3980675565980550E-02,0.5242974416914911E-02,
     4      -.7251282337794251E-02,0.1061593593032694E-01,
     5      -.1671113627819416E-01,0.2909097007012746E-01,
     6      -.5890026882157828E-01,0.1492043402979009D+00,
     7      -.4619225982405832D+00,0.1293766360883851D+01,
     8      -.4085816682588261E-02,0.3505711852470169E-02/
         data c0o1013/
     1      -.2487014600088448E-02,0.2050631923135078E-02,
     2      -.1886444740979416E-02,0.1867350379325216E-02,
     3      -.1965196403364111E-02,0.2183263166258165E-02,
     4      -.2545731837917359E-02,0.3101640449236945E-02,
     5      -.3939245575364620E-02,0.5217848505995973E-02,
     6      -.7238351715771412E-02,0.1061216373228444E-01,
     7      -.1671431211870079E-01,0.2909957213935786E-01,
     8      -.5891344743456321E-01,0.1492220872269491D+00/
         data c0o1014/
     1      -.4619455711118324D+00,0.1293791219805410D+01,
     2      0.3689849798472183E-02,-.3164917755508269E-02,
     3      0.2240929343484948E-02,-.1838440314549350E-02,
     4      0.1675963136920972E-02,-.1636741313255302E-02,
     5      0.1692267033991294E-02,-.1840424716670654E-02,
     6      0.2094413706689541E-02,-.2483328894543166E-02,
     7      0.3059030760525918E-02,-.3911576590402739E-02,
     8      0.5201505696515231E-02,-.7230619325293598E-02/
         data c0o1015/
     1      0.1061100840634686E-01,-.1671824560284382E-01,
     2      0.2910758080908556E-01,-.5892501102370148E-01,
     3      0.1492373482559425D+00,-.4619652097096866D+00,
     4      0.1293812484684016D+01,-.3348771363802578E-02,
     5      0.2871584204149943E-02,-.2030023094560955E-02,
     6      0.1658517521894462E-02,-.1500598038916064E-02,
     7      0.1448970003012998E-02,-.1475755182519306E-02,
     8      0.1575840594172331E-02,-.1755960017659259E-02/
         data c0o1016/
     1      0.2033652022167089E-02,-.2440474813707970E-02,
     2      0.3029843297554737E-02,-.3892874071783748E-02,
     3      0.5190869998710767E-02,-.7226206958491157E-02,
     4      0.1061142599945625E-01,-.1672247149125849E-01,
     5      0.2911493151545171E-01,-.5893517684537703E-01,
     6      0.1492505544682075D+00,-.4619821254640242D+00,
     7      0.1293830815553386D+01,0.3052886916557783E-02,
     8      -.2617275372532583E-02,0.1847823409746474E-02/
         data c0o1017/
     1      -.1504465498800074E-02,0.1352672040484949E-02,
     2      -.1293707201928512E-02,0.1300814592485735E-02,
     3      -.1367262587386037E-02,0.1495892147804457E-02,
     4      -.1697304000333298E-02,0.1991209056175199E-02,
     5      -.2410522057250690E-02,0.3009574167515071E-02,
     6      -.3880137825232743E-02,0.5184006275473141E-02,
     7      -.7223935573465759E-02,0.1061277299353498E-01,
     8      -.1672672405717191E-01,0.2912162511797272E-01/
         data c0o1018/
     1      -.5894413803103014E-01,0.1492620504455811D+00,
     2      -.4619967966767357D+00,0.1293846727778353D+01,
     3      -.2794551314372304E-02,0.2395353159480529E-02,
     4      -.1689296644489486E-02,0.1371430774584606E-02,
     5      -.1226543154859488E-02,0.1163600501907907E-02,
     6      -.1157190807161894E-02,0.1199772575988242E-02,
     7      -.1291809218743084E-02,0.1439614944329908E-02,
     8      -.1655723815928554E-02,0.1961042300257744E-02/
         data c0o1019/
     1      -.2389285568627352E-02,0.2995353660946681E-02,
     2      -.3871441916698648E-02,0.5179669281396179E-02,
     3      -.7223050846023489E-02,0.1061464801494924E-01,
     4      -.1673085481681921E-01,0.2912769484219592E-01,
     5      -.5895206210468745E-01,0.1492721134062909D+00,
     6      -.4620096016305003D+00,0.1293860628460354D+01/
         real *8
     1  c0n1001(16),c0n1002( 9)
         real *8
     1  c0n2001(2)
         equivalence
     1  (c0n1001(16),c0n2001(1)),(c0n1002(1),c0n2001(2))
         data c0n1001/
     1      0.3535533905932738D+00,0.5288213201416558D+00,
     2      0.5164799529547165D+00,0.5165972521995346D+00,
     3      0.5159910085906004D+00,0.5156910453722044D+00,
     4      0.5155317136243025D+00,0.5154379193646634D+00,
     5      0.5153782659828908D+00,0.5153380389949589D+00,
     6      0.5153096512164395D+00,0.5152888806444511D+00,
     7      0.5152732297321620D+00,0.5152611455487935D+00,
     8      0.5152516218069470D+00,0.5152439833991421D+00/
         data c0n1002/
     1      0.5152377637990616D+00,0.5152326322873451D+00,
     2      0.5152283491505048D+00,0.5152247372131421D+00,
     3      0.5152216632391669D+00,0.5152190254784984D+00,
     4      0.5152167451441969D+00,0.5152147604640334D+00,
     5      0.5152130224548324D+00/
         real *8
     1  c1o1001(16),c1o1002(16),c1o1003(16),c1o1004(16),
     2  c1o1005(16),c1o1006(16),c1o1007(16),c1o1008(16),
     3  c1o1009(16),c1o1010(16),c1o1011(16),c1o1012(16),
     4  c1o1013(16),c1o1014(16),c1o1015(16),c1o1016(16)
         real *8
     1  c1o1017(16),c1o1018(16),c1o1019(12)
         real *8
     1  c1o2001(2),c1o2002(2),c1o2003(2),c1o2004(2),
     2  c1o2005(2),c1o2006(2),c1o2007(2),c1o2008(2),
     3  c1o2009(2),c1o2010(2),c1o2011(2),c1o2012(2),
     4  c1o2013(2),c1o2014(2),c1o2015(2),c1o2016(2)
         real *8
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
     1      0.1142857142857143D+01,-.6775305781418212D+00,
     2      0.1319006059903158D+01,0.2919459193420889D+00,
     3      -.4934559023237469D+00,0.1284612269608334D+01,
     4      -.1463306866530963D+00,0.1786618591699440D+00,
     5      -.4614185221468676D+00,0.1287903660759532D+01,
     6      0.8690024315940878E-01,-.8313111636990370E-01,
     7      0.1531503533026661D+00,-.4596239185537454D+00,
     8      0.1290125157332057D+01,-.5748608159105919E-01/
         data c1o1002/
     1      0.4847257156686552E-01,-.6422436055221090E-01,
     2      0.1494319523143283D+00,-.4600231626193731D+00,
     3      0.1291356881247350D+01,0.4082770306512595E-01,
     4      -.3226567871715018E-01,0.3452537504489484E-01,
     5      -.6033227232231428E-01,0.1487492006387314D+00,
     6      -.4604932672362048D+00,0.1292093330996305D+01,
     7      -.3048829056962251E-01,0.2324742988970327E-01,
     8      -.2179668064752097E-01,0.3102606629390084E-01/
         data c1o1003/
     1      -.5924250349282928E-01,0.1486654267802800D+00,
     2      -.4608572082735653D+00,0.1292564171754486D+01,
     3      0.2363274881391562E-01,-.1764283929109082E-01,
     4      0.1521002969109005E-01,-.1880255923748129E-01,
     5      0.2984854304719674E-01,-.5888927564142095E-01,
     6      0.1487232656931069D+00,-.4611240210107406D+00,
     7      0.1292882044814210D+01,-.1885471169353191E-01,
     8      0.1389073551038152E-01,-.1133332210743565E-01/
         data c1o1004/
     1      0.1268835168148405E-01,-.1767084611016102E-01,
     2      0.2938293688096177E-01,-.5877876262456896E-01,
     3      0.1488078193518018D+00,-.4613202532901607D+00,
     4      0.1293106227083859D+01,0.1539210456937340E-01,
     5      -.1124172317319760E-01,0.8835884459139921E-02,
     6      -.9214323261595150E-02,0.1165241352344992E-01,
     7      -.1717253976894689E-01,0.2918637085962058E-01,
     8      -.5875582944017574E-01,0.1488884941877697D+00/
         data c1o1005/
     1      -.4614670136001327D+00,0.1293270027988338D+01,
     2      -.1280281106566607E-01,0.9295417013755780E-02,
     3      -.7118256668143398E-02,0.7048715001479828E-02,
     4      -.8286964999784386E-02,0.1116127842064720E-01,
     5      -.1693577713420396E-01,0.2910245871775180E-01,
     6      -.5876537233832006E-01,0.1489581091155549D+00,
     7      -.4615789232225598D+00,0.1293393243947836D+01,
     8      0.1081600809539639E-01,-.7820214035164350E-02/
         data c1o1006/
     1      0.5877893103388843E-02,-.5601480479547716E-02,
     2      0.6226764740260147E-02,-.7821997355245935E-02,
     3      0.1091034242599093E-01,-.1681810726994455E-01,
     4      0.2906903359794721E-01,-.5878618753187604E-01,
     5      0.1490162860002196D+00,-.4616658773798650D+00,
     6      0.1293488208400261D+01,-.9258274515057571E-02,
     7      0.6673735911729123E-02,-.4948034450679980E-02,
     8      0.4581188588907486E-02,-.4875543585394117E-02/
         data c1o1007/
     1      0.5795905038717457E-02,-.7571605180349110E-02,
     2      0.1077549274148294E-01,-.1675860200978295E-01,
     3      0.2905914790683221E-01,-.5880977020632864E-01,
     4      0.1490644751994925D+00,-.4617346173385966D+00,
     5      0.1293562917259189D+01,0.8014384205004865E-02,
     6      -.5764149859559394E-02,0.4230129208541473E-02,
     7      -.3831093837931320E-02,0.3940102064129107E-02,
     8      -.4480925904729515E-02,0.5554141510437414E-02/
         data c1o1008/
     1      -.7429829495789187E-02,0.1070063662852293E-01,
     2      -.1672896643307362E-01,0.2906036397142780E-01,
     3      -.5883271562187667E-01,0.1491044122983936D+00,
     4      -.4617898094059300D+00,0.1293622732999802D+01,
     5      -.7005350128390809E-02,0.5029877991036825E-02,
     6      -.3662577550296724E-02,0.3260900020951984E-02,
     7      -.3263889189488982E-02,0.3580899150934659E-02,
     8      -.4251947925781747E-02,0.5411756945338882E-02/
         data c1o1009/
     1      -.7346636963299607E-02,0.1065838543428570E-01,
     2      -.1671521221227077E-01,0.2906668055373314E-01,
     3      -.5885375981315581E-01,0.1491376628479380D+00,
     4      -.4618347439317043D+00,0.1293671358192751D+01,
     5      0.6175559185985991E-02,-.4428282029518516E-02,
     6      0.3205075205049051E-02,-.2815571202863272E-02,
     7      0.2757591286052364E-02,-.2937860535804736E-02,
     8      0.3366541059966667E-02,-.4112771081227259E-02/
         data c1o1010/
     1      0.5324849738469612E-02,-.7296607063123244E-02,
     2      0.1063454625345056E-01,-.1671008028177136E-01,
     3      0.2907511899667520E-01,-.5887255419900413E-01,
     4      0.1491655222761891D+00,-.4618717845608337D+00,
     5      0.1293711414651305D+01,-.5484937033733415E-02,
     6      0.3929029054460019E-02,-.2830236056063763E-02,
     7      0.2459946415698052E-02,-.2367395330583250E-02,
     8      0.2461921119062364E-02,-.2738599264246062E-02/
         data c1o1011/
     1      0.3232787823331087E-02,-.4025181901705249E-02,
     2      0.5270373049268436E-02,-.7266076086012625E-02,
     3      0.1062140825794800E-01,-.1670965920347453E-01,
     4      0.2908418366520482E-01,-.5888913569339145E-01,
     5      0.1491890279757467D+00,-.4619026595610162D+00,
     6      0.1293744800326104D+01,0.4904014527179294E-02,
     7      -.3510035889916839E-02,0.2518854169608991E-02,
     8      -.2170652399294282E-02,0.2059357451542508E-02/
         data c1o1012/
     1      -.2099138350105625E-02,0.2277460025744446E-02,
     2      -.2611460115993256E-02,0.3146478804768119E-02,
     3      -.3968562249431351E-02,0.5235561227905121E-02,
     4      -.7247360902109370E-02,0.1061461902563898E-01,
     5      -.1671178031072697E-01,0.2909313074238850E-01,
     6      -.5890369150973211E-01,0.1492090009362180D+00,
     7      -.4619286536245036D+00,0.1293772916144324D+01,
     8      -.4410729572025286E-02,0.3154900646244508E-02/
         data c1o1013/
     1      -.2257095318481762E-02,0.1931617397312164E-02,
     2      -.1811214356421192E-02,0.1815667037327207E-02,
     3      -.1928760848684355E-02,0.2157471965121268E-02,
     4      -.2527667563858384E-02,0.3089302914048962E-02,
     5      -.3931198406799213E-02,0.5213033040252609E-02,
     6      -.7235978379695833E-02,0.1061164420657004E-01,
     7      -.1671521652807657E-01,0.2910160613550497E-01,
     8      -.5891645410348710E-01,0.1492260894213536D+00/
         data c1o1014/
     1      -.4619507363366649D+00,0.1293796814175109D+01,
     2      0.3988300193436812E-02,-.2851226010646629E-02,
     3      0.2034763994839181E-02,-.1731464768641357E-02,
     4      0.1607864524572346E-02,-.1589450383700201E-02,
     5      0.1658449983457342E-02,-.1816049347119564E-02,
     6      0.2076939345715750E-02,-.2471016635985088E-02,
     7      0.3050628773453558E-02,-.3906156185993964E-02,
     8      0.5198366631355797E-02,-.7229230745609892E-02/
         data c1o1015/
     1      0.1061097045116374E-01,-.1671926513075527E-01,
     2      0.2910946008438079E-01,-.5892765526605705E-01,
     3      0.1492408065775364D+00,-.4619696497415897D+00,
     4      0.1293817296782001D+01,-.3623777528824595E-02,
     5      0.2589495489317344E-02,-.1844196369037657E-02,
     6      0.1561936718405979E-02,-.1438760634522544E-02,
     7      0.1405635232753427E-02,-.1444391508153031E-02,
     8      0.1552886868469710E-02,-.1739186259271700E-02/
         data c1o1016/
     1      0.2021536363718475E-02,-.2431920861998035E-02,
     2      0.3024034210726391E-02,-.3889189128547643E-02,
     3      0.5188832925898059E-02,-.7225452107149747E-02,
     4      0.1061167581498733E-01,-.1672352257651435E-01,
     5      0.2911664975646439E-01,-.5893750752729161E-01,
     6      0.1492535605856851D+00,-.4619859690489463D+00,
     7      0.1293834984428987D+01,0.3307041211503140E-02,
     8      -.2362302664295889E-02,0.1679531501979778E-02/
         data c1o1017/
     1      -.1416906929807027E-02,0.1296345313312797E-02,
     2      -.1253929206107680E-02,0.1271727714024781E-02,
     3      -.1345698013115079E-02,0.1479878150168181E-02,
     4      -.1685499922823617E-02,0.1982649621196347E-02,
     5      -.2404486301668852E-02,0.3005511200486056E-02,
     6      -.3877621241194967E-02,0.5182704383595538E-02,
     7      -.7223591226199460E-02,0.1061319001807664E-01,
     8      -.1672775928011793E-01,0.2912318710897332E-01/
         data c1o1018/
     1      -.5894619815540084E-01,0.1492646781147330D+00,
     2      -.4620001454037477D+00,0.1293850362961754D+01,
     3      -.3030089053836884E-02,0.2163811666081833E-02,
     4      -.1536218772201128E-02,0.1291742630649052E-02,
     5      -.1175080221083280E-02,0.1127018429564584E-02,
     6      -.1130203722775566E-02,0.1179541543878712E-02,
     7      -.1276578887283252E-02,0.1428196816944049E-02,
     8      -.1647263528926229E-02,0.1954901002499235E-02/
         data c1o1019/
     1      -.2384973262230813E-02,0.2992488624674879E-02,
     2      -.3869724366855245E-02,0.5178863466531714E-02,
     3      -.7222972632931923E-02,0.1061515657562295E-01,
     4      -.1673184818151175E-01,0.2912911067271837E-01,
     5      -.5895388882140120E-01,0.1492744223379383D+00,
     6      -.4620125364780959D+00,0.1293863817161949D+01/
         real *8
     1  c1n1001(16),c1n1002( 9)
         real *8
     1  c1n2001(2)
         equivalence
     1  (c1n1001(16),c1n2001(1)),(c1n1002(1),c1n2001(2))
         data c1n1001/
     1      0.3061862178478973D+00,0.5163871346414616D+00,
     2      0.5177467552740044D+00,0.5164968509628828D+00,
     3      0.5159099353753987D+00,0.5156464822255225D+00,
     4      0.5155054312875688D+00,0.5154212957313698D+00,
     5      0.5153671351863629D+00,0.5153302391693344D+00,
     6      0.5153039812526615D+00,0.5152846331907770D+00,
     7      0.5152699673615138D+00,0.5152585864062633D+00,
     8      0.5152495778342852D+00,0.5152423253141192D+00/
         data c1n1002/
     1      0.5152364003975593D+00,0.5152314977594625D+00,
     2      0.5152273950575891D+00,0.5152239272693938D+00,
     3      0.5152209698323779D+00,0.5152184272985162D+00,
     4      0.5152162255340715D+00,0.5152143062529306D+00,
     5      0.5152126231184407D+00/
         real *8
     1  c2o1001(16),c2o1002(16),c2o1003(16),c2o1004(16),
     2  c2o1005(16),c2o1006(16),c2o1007(16),c2o1008(16),
     3  c2o1009(16),c2o1010(16),c2o1011(16),c2o1012(16),
     4  c2o1013(16),c2o1014(16),c2o1015(16),c2o1016(16)
         real *8
     1  c2o1017(16),c2o1018(16),c2o1019(12)
         real *8
     1  c2o2001(2),c2o2002(2),c2o2003(2),c2o2004(2),
     2  c2o2005(2),c2o2006(2),c2o2007(2),c2o2008(2),
     3  c2o2009(2),c2o2010(2),c2o2011(2),c2o2012(2),
     4  c2o2013(2),c2o2014(2),c2o2015(2),c2o2016(2)
         real *8
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
     1      0.1496598639455782D+01,-.7022904243767883D+00,
     2      0.1287749913464498D+01,0.2839880427643817D+00,
     3      -.4747435299735907D+00,0.1284032485424784D+01,
     4      -.1433080927657811D+00,0.1662719790840217D+00,
     5      -.4597136533380943D+00,0.1288445423557985D+01,
     6      0.8622597896995995E-01,-.7475641877135735E-01,
     7      0.1512028881013987D+00,-.4596036992460084D+00,
     8      0.1290456949124063D+01,-.5773332242003156E-01/
         data c2o1002/
     1      0.4262829736130900E-01,-.6247430480073939E-01,
     2      0.1491024411237105D+00,-.4601233282886509D+00,
     3      0.1291559895165287D+01,0.4141404049182997E-01,
     4      -.2803718608600327E-01,0.3306518417605249E-01,
     5      -.5990571073117300E-01,0.1486868956673045D+00,
     6      -.4605851561135271D+00,0.1292223479067108D+01,
     7      -.3117635596564006E-01,0.2008301487095460E-01,
     8      -.2060358284102107E-01,0.3059726900409804E-01/
         data c2o1003/
     1      -.5910560086400548E-01,0.1486680503199039D+00,
     2      -.4609276213669669D+00,0.1292651745634860D+01,
     3      0.2432476061112675E-01,-.1520386882154611E-01,
     4      0.1423641391028336E-01,-.1840707729215782E-01,
     5      0.2968457319630027E-01,-.5884411679792561E-01,
     6      0.1487411730737724D+00,-.4611765557352045D+00,
     7      0.1292943492599695D+01,-.1951143460616222E-01,
     8      0.1196290547994161E-01,-.1053364438996112E-01/
         data c2o1004/
     1      0.1233642078611261E-01,-.1750300326479713E-01,
     2      0.2931461570929629E-01,-.5876615100179762E-01,
     3      0.1488274785870382D+00,-.4613597041730946D+00,
     4      0.1293150881305561D+01,0.1599966593178772E-01,
     5      -.9684864509849117E-02,0.8172840261960270E-02,
     6      -.8905851516311636E-02,0.1149157651263879E-01,
     7      -.1709454526326982E-01,0.2915691676618404E-01,
     8      -.5875541629727665E-01,0.1489062463407372D+00/
         data c2o1005/
     1      -.4614970949457933D+00,0.1293303442824723D+01,
     2      -.1335846954715628E-01,0.8014916584121549E-02,
     3      -.6562809181138766E-02,0.6779772427531556E-02,
     4      -.8137763530534069E-02,0.1108103999078955E-01,
     5      -.1689743400852842E-01,0.2909011645429273E-01,
     6      -.5876946112575537E-01,0.1489732274750255D+00,
     7      -.4616022544371123D+00,0.1293418872413897D+01,
     8      0.1132180826856842E-01,-.6750344843635993E-02/
         data c2o1006/
     1      0.5407804747177354E-02,-.5367076304492667E-02,
     2      0.6090669862867247E-02,-.7743538070634241E-02,
     3      0.1086787681890847E-01,-.1679871533731751E-01,
     4      0.2906461342264812E-01,-.5879169574038621E-01,
     5      0.1490289157674072D+00,-.4616842746543486D+00,
     6      0.1293508281052931D+01,-.9718143458399786E-02,
     7      0.5767644061622303E-02,-.4546304562958708E-02,
     8      0.4376412941543996E-02,-.4752470283035338E-02/
         data c2o1007/
     1      0.5721295755383491E-02,-.7527896740509069E-02,
     2      0.1075212452378456E-01,-.1674879969751973E-01,
     3      0.2905848292170224E-01,-.5881543879071115E-01,
     4      0.1490749809833169D+00,-.4617493474166916D+00,
     5      0.1293578923761682D+01,0.8432792030887217E-02,
     6      -.4987657802322289E-02,0.3883697391841293E-02,
     7      -.3651538757881445E-02,0.3829244823120273E-02,
     8      -.4411115567389902E-02,0.5510909569084504E-02/
         data c2o1008/
     1      -.7404465734042133E-02,0.1068748090535842E-01,
     2      -.1672421644660241E-01,0.2906146589567626E-01,
     3      -.5883804637554200E-01,0.1491131702094381D+00,
     4      -.4618017679828627D+00,0.1293635697448853D+01,
     5      -.7386702976451773E-02,0.4357538804006287E-02,
     6      -.3361330481863436E-02,0.3102768130816832E-02,
     7      -.3164155556750301E-02,0.3516197744618027E-02,
     8      -.4210187516514204E-02,0.5385678400437984E-02/
         data c2o1009/
     1      -.7331491009171501E-02,0.1065092814817756E-01,
     2      -.1671319506897208E-01,0.2906856719459350E-01,
     3      -.5885858185057079E-01,0.1491450005374369D+00,
     4      -.4618445746164058D+00,0.1293682002667295D+01,
     5      0.6523941820329171E-02,-.3840799735618213E-02,
     6      0.2941108688891397E-02,-.2675655507139635E-02,
     7      0.2667822734806005E-02,-.2878222352794650E-02,
     8      0.3326793298154881E-02,-.4086798434077216E-02/
         data c2o1010/
     1      0.5308641049075518E-02,-.7287394146802635E-02,
     2      0.1063037279361311E-01,-.1670956007443407E-01,
     3      0.2907729585470335E-01,-.5887683794494450E-01,
     4      0.1491717086626687D+00,-.4618799574571587D+00,
     5      0.1293720259732322D+01,-.5804023098411551E-02,
     6      0.3411530683643339E-02,-.2597312208025417E-02,
     7      0.2335563317849760E-02,-.2286475167414834E-02,
     8      0.2407111290220124E-02,-.2701121708057739E-02/
         data c2o1011/
     1      0.3207435946153211E-02,-.4008542036420061E-02,
     2      0.5260075036265560E-02,-.7260423505679563E-02,
     3      0.1061917247247988E-01,-.1670995605678487E-01,
     4      0.2908640153257593E-01,-.5889291031630849E-01,
     5      0.1491942786478020D+00,-.4619095234988224D+00,
     6      0.1293752228780396D+01,0.5197055959361859E-02,
     7      -.3050888195042831E-02,0.2312002306876552E-02,
     8      -.2059563845539864E-02,0.1986258906540219E-02/
         data c2o1012/
     1      -.2048830425162651E-02,0.2242334919383399E-02,
     2      -.2587039813279377E-02,0.3129834702634053E-02,
     3      -.3957654880836398E-02,0.5228918589871407E-02,
     4      -.7243901330432130E-02,0.1061354478211592E-01,
     5      -.1671251132028457E-01,0.2909526710650210E-01,
     6      -.5890700765235266E-01,0.1492134873244920D+00,
     7      -.4619344713100813D+00,0.1293779214376288D+01,
     8      -.4680581517386178E-02,0.2744883321116127E-02/
         data c2o1013/
     1      -.2072317200106800E-02,0.1831956950465849E-02,
     2      -.1745017401723041E-02,0.1769499565147609E-02,
     3      -.1895964182701738E-02,0.2134157978085514E-02,
     4      -.2511303461796382E-02,0.3078123583485703E-02,
     5      -.3923923077377184E-02,0.5208711120284759E-02,
     6      -.7233897426622133E-02,0.1061127163224350E-01,
     7      -.1671616267047939E-01,0.2910360500115984E-01,
     8      -.5891936714245705E-01,0.1492299477350682D+00/
         data c2o1014/
     1      -.4619557083595505D+00,0.1293802199867614D+01,
     2      0.4237459833407359E-02,-.2482943311560763E-02,
     3      0.1868814607458852E-02,-.1641672854900971E-02,
     4      0.1547756818094516E-02,-.1547059023141196E-02,
     5      0.1627895594334895E-02,-.1793926160407882E-02,
     6      0.2061039601936670E-02,-.2459803200828268E-02,
     7      0.3042982917716272E-02,-.3901241829297335E-02,
     8      0.5195550070824506E-02,-.7228029783676712E-02/
         data c2o1015/
     1      0.1061102102120933E-01,-.1672030024820561E-01,
     2      0.2911130054658852E-01,-.5893021808545729E-01,
     3      0.1492441453619133D+00,-.4619739312239388D+00,
     4      0.1293821937732099D+01,-.3854426285660961E-02,
     5      0.2256950451596422E-02,-.1694419549447329E-02,
     6      0.1480706084249871E-02,-.1384030883623031E-02,
     7      0.1366670255527176E-02,-.1415959167425899E-02,
     8      0.1531979704137905E-02,-.1723865308057854E-02/
         data c2o1016/
     1      0.2010454581742611E-02,-.2424096612075399E-02,
     2      0.3018730457639545E-02,-.3885842717008457E-02,
     3      0.5187010020373140E-02,-.7224818688031445E-02,
     4      0.1061197793281857E-01,-.1672457425061813E-01,
     5      0.2911832922810795E-01,-.5893976773960435E-01,
     6      0.1492564667661617D+00,-.4619896813681664D+00,
     7      0.1293839011636939D+01,0.3521085396640446E-02,
     8      -.2060584231967877E-02,0.1543735592428158E-02/
         data c2o1017/
     1      -.1343137622087246E-02,0.1246373487095269E-02,
     2      -.1218064314544294E-02,0.1245280897792217E-02,
     3      -.1325993849012298E-02,0.1465201601110747E-02,
     4      -.1674663158574613E-02,0.1974787084857373E-02,
     5      -.2398946230431630E-02,0.3001792703828540E-02,
     6      -.3875334972789256E-02,0.5181546527514203E-02,
     7      -.7223325732281792E-02,0.1061363664199486E-01,
     8      -.1672878668383452E-01,0.2912471223849000E-01/
         data c2o1018/
     1      -.5894819738161087E-01,0.1492672217142916D+00,
     2      -.4620033845669483D+00,0.1293853879889182D+01,
     3      -.3229196864717452E-02,0.1888861545591580E-02,
     4      -.1412581715761276E-02,0.1224505223343479E-02,
     5      -.1129326926089508E-02,0.1093954463904356E-02,
     6      -.1105600709165326E-02,0.1161003369951543E-02,
     7      -.1262578452569681E-02,0.1417680333178368E-02,
     8      -.1639463955248565E-02,0.1949239877334874E-02/
         data c2o1019/
     1      -.2381004220729777E-02,0.2989862391243588E-02,
     2      -.3868165720354899E-02,0.5178155432110372E-02,
     3      -.7222945468715954E-02,0.1061568042074859E-01,
     4      -.1673282925356707E-01,0.2913049239343069E-01,
     5      -.5895566288733356E-01,0.1492766601218026D+00,
     6      -.4620153791875030D+00,0.1293866906386631D+01/
         real *8
     1  c2n1001(16),c2n1002( 9)
         real *8
     1  c2n2001(2)
         equivalence
     1  (c2n1001(16),c2n2001(1)),(c2n1002(1),c2n2001(2))
         data c2n1001/
     1      0.2587745847533828D+00,0.5176857134043885D+00,
     2      0.5180692859274097D+00,0.5163344871082221D+00,
     3      0.5158358163563971D+00,0.5156066406166331D+00,
     4      0.5154818002836977D+00,0.5154062034851600D+00,
     5      0.5153569348813896D+00,0.5153230326539061D+00,
     6      0.5152987057794343D+00,0.5152806575626833D+00,
     7      0.5152668980785906D+00,0.5152561680455537D+00,
     8      0.5152476388464776D+00,0.5152407470712659D+00/
         data c2n1002/
     1      0.5152350987722232D+00,0.5152304117631542D+00,
     2      0.5152264796114870D+00,0.5152231484766358D+00,
     3      0.5152203018103151D+00,0.5152178500097277D+00,
     4      0.5152157232720844D+00,0.5152138665663114D+00,
     5      0.5152122360340913D+00/
         real *8
     1  c3o1001(16),c3o1002(16),c3o1003(16),c3o1004(16),
     2  c3o1005(16),c3o1006(16),c3o1007(16),c3o1008(16),
     3  c3o1009(16),c3o1010(16),c3o1011(16),c3o1012(16),
     4  c3o1013(16),c3o1014(16),c3o1015(16),c3o1016(16)
         real *8
     1  c3o1017(16),c3o1018(16),c3o1019(12)
         real *8
     1  c3o2001(2),c3o2002(2),c3o2003(2),c3o2004(2),
     2  c3o2005(2),c3o2006(2),c3o2007(2),c3o2008(2),
     3  c3o2009(2),c3o2010(2),c3o2011(2),c3o2012(2),
     4  c3o2013(2),c3o2014(2),c3o2015(2),c3o2016(2)
         real *8
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
     1      0.1838383838383838D+01,-.6764595037869179D+00,
     2      0.1266771403432156D+01,0.2738695927364719D+00,
     3      -.4599812461302065D+00,0.1285223939113664D+01,
     4      -.1392981875034594D+00,0.1557822737606886D+00,
     5      -.4593438139607948D+00,0.1288977296561474D+01,
     6      0.8470942317062818E-01,-.6733670977271131E-01,
     7      0.1502634424546289D+00,-.4596684118391082D+00,
     8      0.1290760532191950D+01,-.5728786017735738E-01/
         data c3o1002/
     1      0.3726680962777681E-01,-.6146724174108914E-01,
     2      0.1488923664090317D+00,-.4602314167030424D+00,
     3      0.1291744004427038D+01,0.4144535223891172E-01,
     4      -.2405028293256035E-01,0.3215419368511747E-01,
     5      -.5959554309743567E-01,0.1486532913392890D+00,
     6      -.4606736206873396D+00,0.1292342058210030D+01,
     7      -.3142065566234741E-01,0.1703364158137929E-01,
     8      -.1982327206152503E-01,0.3026942590655456E-01/
         data c3o1003/
     1      -.5900374015660764E-01,0.1486780154220040D+00,
     2      -.4609938318949585D+00,0.1292732119451934D+01,
     3      0.2465896133700525E-01,-.1281182886065172E-01,
     4      0.1357974723422165E-01,-.1809550809048801E-01,
     5      0.2955512614082315E-01,-.5881105876124218E-01,
     6      0.1487606280117099D+00,-.4612257308001961D+00,
     7      0.1293000298510611D+01,-.1987612615210731E-01,
     8      0.1004468417716853E-01,-.9982668843793331E-02/
         data c3o1004/
     1      0.1205328389074789E-01,-.1736652197184620E-01,
     2      0.2925982749021812E-01,-.5875805990157425E-01,
     3      0.1488470261015184D+00,-.4613966717102471D+00,
     4      0.1293192433293470D+01,0.1636578232342826E-01,
     5      -.8117050731650383E-02,0.7708952906727192E-02,
     6      -.8653726390777835E-02,0.1135822326169272E-01,
     7      -.1702989309311321E-01,0.2913328760917042E-01,
     8      -.5875672456508869E-01,0.1489234533974047D+00/
         data c3o1005/
     1      -.4615253642895327D+00,0.1293334715420491D+01,
     2      -.1371199776362709E-01,0.6712279610200508E-02,
     3      -.6169796666542089E-02,0.6557217013022203E-02,
     4      -.8012265344611464E-02,0.1101322607208960E-01,
     5      -.1686530867094663E-01,0.2908043976260225E-01,
     6      -.5877416042337315E-01,0.1489877482177430D+00,
     7      -.4616242548070885D+00,0.1293442978276503D+01,
     8      0.1165641678012199E-01,-.5652529809680602E-02/
         data c3o1006/
     1      0.5072377086947562E-02,-.5171162221978061E-02,
     2      0.5974896697471463E-02,-.7676317726616042E-02,
     3      0.1083153621967531E-01,-.1678241559224731E-01,
     4      0.2906145584813653E-01,-.5879735708393518E-01,
     5      0.1490410082086371D+00,-.4617016816839036D+00,
     6      0.1293527243879404D+01,-.1003148109685191E-01,
     7      0.4830921161057733E-02,-.4257826286587546E-02,
     8      0.4203863916417864E-02,-.4646810353358636E-02/
         data c3o1007/
     1      0.5656700415135024E-02,-.7489976465342345E-02,
     2      0.1073197231714611E-01,-.1674060810174420E-01,
     3      0.2905842030621762E-01,-.5882107286824813E-01,
     4      0.1490850341022128D+00,-.4617633296736524D+00,
     5      0.1293594103272990D+01,0.8724568970278199E-02,
     6      -.4179679711629484E-02,0.3633724937205651E-02,
     7      -.3499217843032574E-02,0.3733345786698283E-02,
     8      -.4350161307097239E-02,0.5473021466350455E-02/
         data c3o1008/
     1      -.7382267942646474E-02,0.1067609554065456E-01,
     2      -.1672033248180128E-01,0.2906284651193011E-01,
     3      -.5884327096136049E-01,0.1491215560008230D+00,
     4      -.4618131533654489D+00,0.1293648033612107D+01,
     5      -.7657671754569595E-02,0.3653927601568978E-02,
     6      -.3143162102975900E-02,0.2967861931761315E-02,
     7      -.3077324296746105E-02,0.3459305076874847E-02,
     8      -.4173294897670319E-02,0.5362620587939116E-02/
         data c3o1009/
     1      -.7318158037404994E-02,0.1064448105813727E-01,
     2      -.1671165214265812E-01,0.2907056972401936E-01,
     3      -.5886327680547617E-01,0.1491520347740504D+00,
     4      -.4618539596935732D+00,0.1293692161588174D+01,
     5      0.6775348030483059E-02,-.3222872734788584E-02,
     6      0.2749405462585644E-02,-.2555720013713063E-02,
     7      0.2589239114814908E-02,-.2825467608632762E-02,
     8      0.3291446072168628E-02,-.4063653965890560E-02/
         data c3o1010/
     1      0.5294216707824274E-02,-.7279257492780542E-02,
     2      0.1062679145179605E-01,-.1670930052567494E-01,
     3      0.2907950592913596E-01,-.5888099560513560E-01,
     4      0.1491776477102452D+00,-.4618877793318239D+00,
     5      0.1293728723829871D+01,-.6037307750560324E-02,
     6      0.2864751939992404E-02,-.2427796520364442E-02,
     7      0.2228509570359918E-02,-.2215303727994355E-02,
     8      0.2358377836287630E-02,-.2667606465580433E-02/
         data c3o1011/
     1      0.3184701143777380E-02,-.3993616193772344E-02,
     2      0.5250868548770423E-02,-.7255428344581607E-02,
     3      0.1061729035743856E-01,-.1671039401388116E-01,
     4      0.2908861110121169E-01,-.5889656869100831E-01,
     5      0.1491993270657287D+00,-.4619161075355116D+00,
     6      0.1293759354225918D+01,0.5413695489555590E-02,
     7      -.2563792791299711E-02,0.2161223809239047E-02,
     8      -.1963620443730237E-02,0.1921704473245002E-02/
         data c3o1012/
     1      -.2003899599715494E-02,0.2210771670500463E-02,
     2      -.2565023803015875E-02,0.3114811661117810E-02,
     3      -.3947821873895249E-02,0.5222962874020535E-02,
     4      -.7240852465477589E-02,0.1061268413159611E-01,
     5      -.1671331523540951E-01,0.2909737524026666E-01,
     6      -.5891022025678507E-01,0.1492178074851494D+00,
     7      -.4619400632933377D+00,0.1293785268670177D+01,
     8      -.4882004172461444E-02,0.2308316179337536E-02/
         data c3o1013/
     1      -.1937472922837041E-02,0.1745626555245426E-02,
     2      -.1686350015024770E-02,0.1728105051240581E-02,
     3      -.1866369039722020E-02,0.2113043324936399E-02,
     4      -.2496456983805460E-02,0.3067981448557720E-02,
     5      -.3917340001378585E-02,0.5204831641486559E-02,
     6      -.7232077505202248E-02,0.1061102837075546E-01,
     7      -.1671714239669402E-01,0.2910556714140844E-01,
     8      -.5892218954833143E-01,0.1492336686180165D+00/
         data c3o1014/
     1      -.4619604965050612D+00,0.1293807387014562D+01,
     2      0.4425004755930468E-02,-.2089511369515142E-02,
     3      0.1747611018662201E-02,-.1563690237322377E-02,
     4      0.1494320124131398E-02,-.1508918933548601E-02,
     5      0.1600221476515030E-02,-.1773810142437854E-02,
     6      0.2046551204439414E-02,-.2449578172765711E-02,
     7      0.3036018611339842E-02,-.3896783869056613E-02,
     8      0.5193023926541157E-02,-.7226996773660282E-02/
         data c3o1015/
     1      0.1061114922654183E-01,-.1672134624350544E-01,
     2      0.2911310183360050E-01,-.5893270224392072E-01,
     3      0.1492473698097953D+00,-.4619780614839541D+00,
     4      0.1293826415471447D+01,-.4029326455398975E-02,
     5      0.1900622216162741E-02,-.1584966585529439E-02,
     6      0.1410000671175163E-02,-.1335241876094178E-02,
     7      0.1331505599794052E-02,-.1390122214621575E-02,
     8      0.1512902555281851E-02,-.1709851247684058E-02/
         data c3o1016/
     1      0.2000306405117839E-02,-.2416932879896170E-02,
     2      0.3013884674389825E-02,-.3882803037444815E-02,
     3      0.5185380662906018E-02,-.7224294122702419E-02,
     4      0.1061232554335619E-01,-.1672562382975722E-01,
     5      0.2911997020427010E-01,-.5894195993281554E-01,
     6      0.1492592772219289D+00,-.4619932682380998D+00,
     7      0.1293842903464352D+01,0.3684463740724056E-02,
     8      -.1736392845948832E-02,0.1444464911653062E-02/
         data c3o1017/
     1      -.1278801327370738E-02,0.1201717994795729E-02,
     2      -.1185609153630243E-02,0.1221177540723654E-02,
     3      -.1307958141765874E-02,0.1451732035347362E-02,
     4      -.1664702893181058E-02,0.1967557701419813E-02,
     5      -.2393857219243258E-02,0.2998387826754312E-02,
     6      -.3873258189590398E-02,0.5180519192837333E-02,
     7      -.7223130915508355E-02,0.1061410857262680E-01,
     8      -.1672980479411856E-01,0.2912620108820390E-01/
         data c3o1018/
     1      -.5895013784064722E-01,0.1492696846741107D+00,
     2      -.4620065188273101D+00,0.1293857283591458D+01,
     3      -.3382069421900768E-02,0.1592681692103825E-02,
     4      -.1322183292861489E-02,0.1165765391015287E-02,
     5      -.1088353271576835E-02,0.1063961017521318E-02,
     6      -.1083118607826162E-02,0.1143987331679778E-02,
     7      -.1249691112251304E-02,0.1407983342466242E-02,
     8      -.1632266615005381E-02,0.1944017238800149E-02/
         data c3o1019/
     1      -.2377348946900293E-02,0.2987454422566193E-02,
     2      -.3866752063404288E-02,0.5177536168662319E-02,
     3      -.7222963968616289E-02,0.1061621683636077E-01,
     4      -.1673379728426946E-01,0.2913184071338067E-01,
     5      -.5895738613595815E-01,0.1492788295665616D+00,
     6      -.4620181335267978D+00,0.1293869900197631D+01/
         real *8
     1  c3n1001(16),c3n1002( 9)
         real *8
     1  c3n2001(2)
         equivalence
     1  (c3n1001(16),c3n2001(1)),(c3n1002(1),c3n2001(2))
         data c3n1001/
     1      0.2134781409574916D+00,0.5326798766161688D+00,
     2      0.5173940267558330D+00,0.5161903884680898D+00,
     3      0.5157702239723611D+00,0.5155712084184198D+00,
     4      0.5154605517081959D+00,0.5153924830825210D+00,
     5      0.5153475725288877D+00,0.5153163642515745D+00,
     6      0.5152937906436552D+00,0.5152769318577559D+00,
     7      0.5152640073722731D+00,0.5152538805696761D+00,
     8      0.5152457979180433D+00,0.5152392437130649D+00/
         data c3n1002/
     1      0.5152338553019374D+00,0.5152293716054152D+00,
     2      0.5152256007768657D+00,0.5152223992741018D+00,
     3      0.5152196579603157D+00,0.5152172926586393D+00,
     4      0.5152152376003326D+00,0.5152134407957232D+00,
     5      0.5152118607088319D+00/
c 
c       construct the orthogonal polynomials of orders 0, 4, 8, ...
c 
        n4=n/4+1
        index=0
        call dorteva(z,n4,c0o1001,c0n1001,index,w)
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
        call dorteva(z,n4,c1o1001,c1n1001,index,w)
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
        call dorteva(z,n4,c2o1001,c2n1001,index,w)
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
        call dorteva(z,n4,c3o1001,c3n1001,index,w)
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
        subroutine dorteva(z,m,coefs,coenorm,index,pols)
        implicit real *8 (a-h,o-z)
        save
        complex *16 z,pols(1),z4
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
  
  
  
