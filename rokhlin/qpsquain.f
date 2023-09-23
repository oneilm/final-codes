        implicit real *16 (a-h,o-z)
        complex *32  pol(50 000),work(41000),
     1      amatr(850 )
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
     1   weights(4000),rlens(4000),corners(2,4),rnorms(2000),
     2      promax(2000)
        complex *32 com,work(m,m),
     1  ima,zs(4000),ww(2000),zs2(4000)
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
c        construct the outer frame
c 
        three=3
        done=1
        three=done/3
        do 1200 i=1,n*4
        zs2(i)=zs(i)*three
        zs(i+n*4)=zs2(i)
 1200 continue
cccc        call prinq('zs2=*',zs2,n*8)
c 
c       plot the square as created
c 
        do 1400 i=1,n*4*2
        com=zs(i)
        x(i)=rea(1)
        y(i)=rea(2)
 1400 continue
  
        iw=21
        call RSPLOT(X,Y,N*4*2,IW,'square as created*')
cccc        if(2 .ne. 3) stop
  
c        construct the orthogonal polynomials in 1/z on the square
c 
        do 1700 i=1,4*n
        call qpsquain(zs(i),m,ww)
        do 1600 j=1,m
        pol(i,j)=ww(j)
 1600 continue
 1700 continue
  
         call prinq('after qpsquain, pol(m)=*',pol(1,m-1),8*n)
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
cccc        call prinq('and the matrix of inner products*',work,m*m*2)
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
  
c 
c       test the projections
c 
        do 3200 i=1,n
         call prinf('i=*',i,1)
        call projtest(pol,n*4,m,zs2(i),zs,weights,promax)
 3200 continue
        call prinq('and promax=*',promax,m)
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine projtest(pol,n4,m,z0,zs,weights,promax)
        implicit real *16 (a-h,o-z)
        save
        complex *32 pol(n4,1),z0,zs(1),proj(4000),
     1      f(4000),cd
        dimension weights(1),promax(1)
c 
c       construct the potential of the charge on the
c       boundary
c 
        done=1
        do 1400 i=1,n4
        f(i)=done/(zs(i)-z0)
 1400 continue
cccc        call prinq('f as constructed*',f,n4*2)
c 
c       evaluate the projections of the potentials on all
c       orthonormal polynomials
c 
        do 2400 i=1,m
        cd=0
        do 2200 j=1,n4
ccc        cd=cd+weights(j)*pol(j,i)*f(j)
        cd=cd+weights(j)*pol(j,i)*qconjg(f(j) )
 2200 continue
        proj(i)=cd
 2400 continue
cccc        call prinq('and proj=*',proj,m*2)
c 
c       adjust the list of maximal projections
c 
        do 2600 i=1,m
        dd=cqabs(proj(i))
        if(dd .gt. promax(i)) promax(i)=dd
 2600 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine matvec(a,x,y,n)
        implicit complex *32 (a-h,o-z)
        save
        dimension a(n,n),x(1),y(1)
c 
        do 1400 i=1,n
        cd=0
        do 1200 j=1,n
        cd=cd+a(i,j)*qconjg(x(j))
 1200 continue
        y(i)=cd
 1400 continue
        return
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
  
        subroutine qpsquain(z,n,pols)
        implicit real *16 (a-h,o-z)
        save
        complex *32 pols(1),w(50),z
c 
c       this subroutine evaluates at the point z \in C
c       a set of polynomials in 1/z orthogonal on the square
c       with the corners at the points (1,-1),(1,1),
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
c  pols - the array of values at z of polynomials in 1/z
c       orthogonal on the square.
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
     1      -.83333333333333333333333333333334Q-01,
     2      0.23211047127176688168822634929888Q-01,
     3      -.31712273215463247884482265258382Q+00,
     4      -.16924880728180495898224974401788Q-01,
     5      0.60547503353298738410148851480740Q-01,
     6      -.29800486623510989048207167833939Q+00,
     7      0.13901226165827048063685994820129Q-01,
     8      -.31427825230991273212595653174339Q-01,
     9      0.50977559868508883666206395658249Q-01,
     a      -.29642266110660928796872801008052Q+00,
     1      -.12035700952569176616769545618850Q-01,
     2      0.20038047331368438237427929853113Q-01,
     3      -.25552507506684712062029257440938Q-01,
     4      0.49995641810775932368278177133326Q-01,
     5      -.29594743592109428915225837819034Q+00,
     6      0.10738027095494403941846040944161Q-01/
         data c0o1002/
     1      -.14117355374175716267202404389722Q-01,
     2      0.16103872237864699653736899676635Q-01,
     3      -.24862346502013665991493021373410Q-01,
     4      0.49662361256528451564955546785386Q-01,
     5      -.29574504014002031862296104394689Q+00,
     6      -.97685465482991548323350443871223Q-02,
     7      0.10551429415656707385079314522808Q-01,
     8      -.11360713366272683629158745201651Q-01,
     9      0.15591880060694139247989313768072Q-01,
     a      -.24607240246882516732363395175069Q-01,
     1      0.49509235208671219350685702318167Q-01,
     2      -.29564092440710747799130865266313Q+00,
     3      0.90089688421741049763555152884699Q-02,
     4      -.81986047935206547167135367697945Q-02,
     5      0.85743597178176474331224602484202Q-02,
     6      -.10968229272078954283428688956739Q-01/
         data c0o1003/
     1      0.15389694779618102283639427084423Q-01,
     2      -.24483193225042879255537831670119Q-01,
     3      0.49426296436895717608904956293605Q-01,
     4      -.29558053616739057393516687022836Q+00,
     5      -.83932424006243699681675410341893Q-02,
     6      0.65463126490021303005534908453051Q-02,
     7      -.67691392031139103815094177844056Q-02,
     8      0.82669923313292121103685951239329Q-02,
     9      -.10804324185514665879325531713068Q-01,
     a      0.15286829378115039629794068025611Q-01,
     1      -.24413233933932167673249897206434Q-01,
     2      0.49376364252935889998829944989506Q-01,
     3      -.29554246510536227933136123281253Q+00,
     4      0.78811851106207148896447945998471Q-02,
     5      -.53320713826143423157780333180602Q-02,
     6      0.55183622409327852719732697657200Q-02/
         data c0o1004/
     1      -.65250423760923753347468892846332Q-02,
     2      0.81319705033433574291049079966276Q-02,
     3      -.10717680993839262102679427524311Q-01,
     4      0.15226863987130199424277136526044Q-01,
     5      -.24369823318342138372044535168473Q-01,
     6      0.49343986207402818020340708331849Q-01,
     7      -.29551694454003947059800874927973Q+00,
     8      -.74467564265714556970965251552310Q-02,
     9      0.44083514806303537810274893121438Q-02,
     a      -.46082423103734827801320877579986Q-02,
     1      0.53228394516433718505102235019208Q-02,
     2      -.64124995415492171366596395949727Q-02,
     3      0.80581325028421015385861243267910Q-02,
     4      -.10665725254523520899118659118818Q-01,
     5      0.15188701662287123717146732035414Q-01,
     6      -.24341006657810172014562693484379Q-01/
         data c0o1005/
     1      0.49321799361553455852123870742541Q-01,
     2      -.29549901237886886573676583952286Q+00,
     3      0.70722380965527609726000183338314Q-02,
     4      -.36862504588946458231736636454476Q-02,
     5      0.39208100333285244557773384891993Q-02,
     6      -.44509662049362676188551385871550Q-02,
     7      0.52282188625276203107772231121041Q-02,
     8      -.63489995772257839376876685277987Q-02,
     9      0.80127362924985336509701670613545Q-02,
     a      -.10631930634968572451705077198893Q-01,
     1      0.15162858185464141040442675376832Q-01,
     2      -.24320887888151211078585908825885Q-01,
     3      0.49305932389235887056958391293671Q-01,
     4      -.29548593389906395361107371166565Q+00,
     5      -.67450985122787941124757093122268Q-02,
     6      0.31092099914833596872668739164728Q-02/
         data c0o1006/
     1      -.33860824690604161385756026247625Q-02,
     2      0.37942692804141995092115416248686Q-02,
     3      -.43709156742567591663579636395478Q-02,
     4      0.51732167879103654815403049521580Q-02,
     5      -.63090623525264297190812310174778Q-02,
     6      0.79826314385652593640117198106211Q-02,
     7      -.10608642928887976852445058133200Q-01,
     8      0.15144521453517533346209316190731Q-01,
     9      -.24306279167208752005461451426044Q-01,
     a      0.49294191851740593147856235430148Q-01,
     1      -.29547610272582396334441097755885Q+00,
     2      0.64561903856048830921325745910235Q-02,
     3      -.26396502114034758802465012976576Q-02,
     4      0.29601340418733766530111380175197Q-02,
     5      -.32846590402527821637377795299456Q-02,
     6      0.37262555087923563174792074385733Q-02/
         data c0o1007/
     1      -.43230037802853952365725837792236Q-02,
     2      0.51378830072343186908574471376230Q-02,
     3      -.62821088083891394601519269461509Q-02,
     4      0.79615619013285490517753672126989Q-02,
     5      -.10591882788304063173354965075948Q-01,
     6      0.15131028071896080046751231425448Q-01,
     7      -.24295332182796399742311391997521Q-01,
     8      0.49285259618245138823042927485335Q-01,
     9      -.29546852554002158393923206562622Q+00,
     a      -.61986594418943854983239485471685Q-02,
     1      0.22516950849160740088243623483820Q-02,
     2      -.26141174261421880197186796100552Q-02,
     3      0.28795274635724840122258833084369Q-02,
     4      -.32267244522874486789644724550679Q-02,
     5      0.36843342627634123756548826152765Q-02,
     6      -.42915958637128926814659587984899Q-02/
         data c0o1008/
     1      0.51136459550422457244253035480505Q-02,
     2      -.62629772031346306296169544397219Q-02,
     3      0.79462035159647295519516293590746Q-02,
     4      -.10579402559490680026807380081563Q-01,
     5      0.15120803175583389702115016705690Q-01,
     6      -.24286914509762598041465728874740Q-01,
     7      0.49278304893872080981959238823532Q-01,
     8      -.29546256156905323653302018457006Q+00,
     9      0.59672546567319433929635829543803Q-02,
     a      -.19269885746554817347024627869373Q-02,
     1      0.23283702697675625575675909071780Q-02,
     2      -.25509708096976348638514874903380Q-02,
     3      0.28301306700658786693900819097134Q-02,
     4      -.31899198212757800580631335080925Q-02,
     5      0.36563071530431351957473123623575Q-02,
     6      -.42697197577008611302629832928931Q-02/
         data c0o1009/
     1      0.50962175126388474013818394000124Q-02,
     2      -.62488686219443832071876213754886Q-02,
     3      0.79346442386543876790662031037347Q-02,
     4      -.10569850165579247357929734526057Q-01,
     5      0.15112865624985808943678854016332Q-01,
     6      -.24280300560295160166645016903067Q-01,
     7      0.49272782980971606609702420516064Q-01,
     8      -.29545778253974042504286098558987Q+00,
     9      -.57578762967118982792072187097997Q-02,
     a      0.16521758452943285714930388772568Q-02,
     1      -.20890676879700864727455860864571Q-02,
     2      0.22800191115868562255093908599257Q-02,
     3      -.25088784361389463090191867375371Q-02,
     4      0.27977373173933269159766248013520Q-02,
     5      -.31648291270820845530189185020644Q-02,
     6      0.36364981185829494680445547564792Q-02/
         data c0o1010/
     1      -.42537965390243067850281579320942Q-02,
     2      0.50832265694183554458732757237939Q-02,
     3      -.62381458801962524751514013115686Q-02,
     4      0.79257159863455746058788241088702Q-02,
     5      -.10562370757337339257103462886355Q-01,
     6      0.15106577668950440697768764103247Q-01,
     7      -.24275007831750049938681627769015Q-01,
     8      0.49268324704921011707246760997638Q-01,
     9      -.29545389344085025926855584606475Q+00,
     a      0.55672705485038630183424930588245Q-02,
     1      -.14173272623838742804658647375239Q-02,
     2      0.18862293874200307921851859001910Q-02,
     3      -.20533662536634116341159927138899Q-02,
     4      0.22442303871471442714888388998165Q-02,
     5      -.24803196438894100601005664663106Q-02,
     6      0.27752164152268844569070131873362Q-02/
         data c0o1011/
     1      -.31468419012431001329222127503320Q-02,
     2      0.36219122290765571166420942838879Q-02,
     3      -.42418083243566944614086534151087Q-02,
     4      0.50732638733435010367915708474986Q-02,
     5      -.62297944945099323512605665741869Q-02,
     6      0.79186703624592703069350718224852Q-02,
     7      -.10556401528893439447296210455848Q-01,
     8      0.15101509886373210562953149073718Q-01,
     9      -.24270705254261884894263907150032Q-01,
     a      0.49264672662401838442325235961648Q-01,
     1      -.29545068570930239034291568105276Q+00,
     2      -.53928176565019299051862082890518Q-02,
     3      0.12149192249245639936825682567032Q-02,
     4      -.17124852746182299892525508427780Q-02,
     5      0.18614277343275721453128499949231Q-02,
     6      -.20230591684034780422949525787861Q-02/
         data c0o1012/
     1      0.22190291324940109106851133067088Q-02,
     2      -.24600628804672841259495554056698Q-02,
     3      0.27588444941676966412882673341718Q-02,
     4      -.31334500759402550905164927568528Q-02,
     5      0.36108267478049185751568601608660Q-02,
     6      -.42325369797977668646504818180819Q-02,
     7      0.50654444077389102916018489570985Q-02,
     8      -.62231564868359436283872667134643Q-02,
     9      0.79130089151660762384551455846053Q-02,
     a      -.10551559262246147820991459858595Q-01,
     1      0.15097364386567629166277009188124Q-01,
     2      -.24267159488713784685939066470364Q-01,
     3      0.49261642933628351878563202612390Q-01,
     4      -.29544800852707469768984117333189Q+00,
     5      0.52323815037123981482425893997972Q-02,
     6      -.10391560006539772238775603846998Q-02/
         data c0o1013/
     1      0.15622844867974154885499527790772Q-02,
     2      -.16971402601894221217263578211860Q-02,
     3      0.18359197888404480357096770133479Q-02,
     4      -.20008168353096079499815204986697Q-02,
     5      0.22007792845534236602704876112127Q-02,
     6      -.24451309875384323580972431256704Q-02,
     7      0.27465237617457304786582856773867Q-02,
     8      -.31231801282288150965986461229409Q-02,
     9      0.36021853787450358063138760889534Q-02,
     a      -.42252075786748623385024249233401Q-02,
     1      0.50591875451783064529316399895937Q-02,
     2      -.62177889353665829903221744792019Q-02,
     3      0.79083887792584468566286175342609Q-02,
     4      -.10547575483230664339496859176365Q-01,
     5      0.15093929175048553487467303062623Q-01,
     6      -.24264202212454175777668370644142Q-01/
         data c0o1014/
     1      0.49259101266009476022432263453812Q-01,
     2      -.29544575064679830458677652616543Q+00,
     3      -.50842006259052543825776653784375Q-02,
     4      0.88550705881965307253105941367839Q-03,
     5      -.14313734451664485692815607159854Q-02,
     6      0.15551941485476597037155370765565Q-02,
     7      -.16758593254692811723457594439719Q-02,
     8      0.18163002262569265586000444009375Q-02,
     9      -.19843556019731767418572130896057Q-02,
     a      0.21871375449981408311228087850884Q-02,
     1      -.24337751902115025713166527987602Q-02,
     2      0.27369930850485775777394388932177Q-02,
     3      -.31151140221009960491750508249298Q-02,
     4      0.35953075665036978531918145892694Q-02,
     5      -.42193061184304903764262366927370Q-02,
     6      0.50540985182400763745044762754278Q-02/
         data c0o1015/
     1      -.62133843033962996696849024649908Q-02,
     2      0.79045675791324769331391861704370Q-02,
     3      -.10544257450130848005001911272580Q-01,
     4      0.15091050005174151931412515579925Q-01,
     5      -.24261709511200180220695123155760Q-01,
     6      0.49256947800268446941196054286284Q-01,
     7      -.29544382855493235140013747076497Q+00,
     8      0.49468078114960353740803168269750Q-02,
     9      -.75038439302994618634437038426587Q-03,
     a      0.13164424642156568033245966309530Q-02,
     1      -.14315278013451421452005305452721Q-02,
     2      0.15376564786692938044932507616952Q-02,
     3      -.16585767244668376475859515165071Q-02,
     4      0.18014409202409405365900617622916Q-02,
     5      -.19718749889671966449125727028851Q-02,
     6      0.21766545358905912426327344553994Q-02/
         data c0o1016/
     1      -.24249168633120590447969684631802Q-02,
     2      0.27294534877059244436402325012386Q-02,
     3      -.31086526985796364680187722545591Q-02,
     4      0.35897370677382504608319023860340Q-02,
     5      -.42144798027276772070024051746785Q-02,
     6      0.50499008129065545462226162014193Q-02,
     7      -.62097233737276931257676311016655Q-02,
     8      0.79013699069906338743695145372599Q-02,
     9      -.10541463821393172896288999390724Q-01,
     a      0.15088612447219123127028243142339Q-01,
     1      -.24259588512658418954604469017643Q-01,
     2      0.49255107008862043689105692374039Q-01,
     3      -.29544217856812960489127225746115Q+00,
     4      -.48189698280727311181196874930560Q-02,
     5      0.63091302742340002942878789224275Q-03,
     6      -.12148809956306480860540370048476Q-02/
         data c0o1017/
     1      0.13229860734145653387660100887671Q-02,
     2      -.14173212014814928489335286104133Q-02,
     3      0.15224660259654503277651202569981Q-02,
     4      -.16451585830261710518911342166347Q-02,
     5      0.17900096590418907538839250978677Q-02,
     6      -.19621843651115463594221867311112Q-02,
     7      0.21684096242146289454719833219591Q-02,
     8      -.24178602938920434741386436184197Q-02,
     9      0.27233767659911228748591193403040Q-02,
     a      -.31033902933617792107680443719003Q-02,
     1      0.35851578533364547255394673213311Q-02,
     2      -.42104794548511077571645106336889Q-02,
     3      0.50463957861644735530168140406267Q-02,
     4      -.62066463142584228761944059094008Q-02,
     5      0.78986662023540912452718580778046Q-02,
     6      -.10539088999568987091434621493128Q-01/
         data c0o1018/
     1      0.15086530140913223698904349081885Q-01,
     2      -.24257768491740254400440196433797Q-01,
     3      0.49253520908778599888333429042558Q-01,
     4      -.29544075143756593832565053308091Q+00,
     5      0.46996415797616763422764515761143Q-02,
     6      -.52476491126469625118533667167526Q-03,
     7      0.11246046589308004246576211057745Q-02,
     8      -.12270840817029484027040389729818Q-02,
     9      0.13117568457656709223910531878284Q-02,
     a      -.14040117836050778708741150253075Q-02,
     1      0.15103497321060570724633744580088Q-02,
     2      -.16346792849827631192754700933728Q-02,
     3      0.17810409182415476903162064561671Q-02,
     4      -.19545007777544811481684900186360Q-02,
     5      0.21617972431389003267746297854889Q-02,
     6      -.24121392640736430788347350659595Q-02/
         data c0o1019/
     1      0.27184012481074146251014235753990Q-02,
     2      -.30990431448443157735134161332972Q-02,
     3      0.35813449274472966369604545904531Q-02,
     4      -.42071247479071896630545945130159Q-02,
     5      0.50434375951254531853388268344699Q-02,
     6      -.62040342766786386264222922463350Q-02,
     7      0.78963590379707394469478205038611Q-02,
     8      -.10537052780822447966513184623548Q-01,
     9      0.15084736906569943092759558440208Q-01,
     a      -.24256194815214947222387586925835Q-01,
     1      0.49252144385180893409131607647196Q-01,
     2      -.29543950859153243769181147680858Q+00/
         real *16
     1  c0n1001(16),c0n1002( 9)
         real *16
     1  c0n2001(2)
         equivalence
     1  (c0n1001(16),c0n2001(1)),(c0n1002(1),c0n2001(2))
         data c0n1001/
     1      0.35355339059327376220042218105244Q+00,
     2      0.14622959690121313546358260005829Q+01,
     3      0.13583210685174274898977475975570Q+01,
     4      0.13553248581034735071861054682323Q+01,
     5      0.13546116748830643584595354111168Q+01,
     6      0.13543419348999097264273313299400Q+01,
     7      0.13542129260943665108811647388907Q+01,
     8      0.13541417047309853214339511364456Q+01,
     9      0.13540983777605473038305313316318Q+01,
     a      0.13540701067819275587996569530213Q+01,
     1      0.13540506565726352737044188942057Q+01,
     2      0.13540367091042304391754586506701Q+01,
     3      0.13540263692876963362896703740486Q+01,
     4      0.13540184918925274007573488488699Q+01,
     5      0.13540123521757648503011448958122Q+01,
     6      0.13540074735407691349897176286921Q+01/
         data c0n1002/
     1      0.13540035322291192962493940916848Q+01,
     2      0.13540003021018390626586876486882Q+01,
     3      0.13539976213508984507823742232212Q+01,
     4      0.13539953717309198771956661499930Q+01,
     5      0.13539934652083175075883148480713Q+01,
     6      0.13539918351494827249187989890553Q+01,
     7      0.13539904303671879943319708865048Q+01,
     8      0.13539892110131685821928717652236Q+01,
     9      0.13539881456906609949175652803994Q+01/
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
     1      -.16860329539459689051258917558168Q+00,
     2      0.36448908691553899251605683413764Q-01,
     3      -.30693996358347489672803634455535Q+00,
     4      -.25368574586991859628031640396680Q-01,
     5      0.54866288948492238807881492804882Q-01,
     6      -.29726589323903768541131931178576Q+00,
     7      0.20272048675003634148498336731864Q-01,
     8      -.27558559727411542567810691271005Q-01,
     9      0.50469948913772180685747867331449Q-01,
     a      -.29617936595342178097921460992643Q+00,
     1      -.17220434527850774750080327659310Q-01,
     2      0.17156260658097846828259200063029Q-01,
     3      -.25157957533784251168122846984630Q-01,
     4      0.49813722847980518181272448982227Q-01,
     5      -.29583847433463393881133896123719Q+00,
     6      0.15144483148948031403051998178852Q-01/
         data c1o1002/
     1      -.11855832359389535619662727808190Q-01,
     2      0.15780452529766116626341472124458Q-01,
     3      -.24713677325685627939420158636008Q-01,
     4      0.49576539370872947231230398174938Q-01,
     5      -.29568776803018336401444206043890Q+00,
     6      -.13620584780605124764310544215859Q-01,
     7      0.87144307801042282957702380131038Q-02,
     8      -.11086614811440364889480162219581Q-01,
     9      0.15466020307996995130648258829840Q-01,
     a      -.24534692135987615151940101256616Q-01,
     1      0.49462516359229368698314879044204Q-01,
     2      -.29560768796923153788510924004567Q+00,
     3      0.12443717949226454996642180711984Q-01,
     4      -.66695100400784258433652711527251Q-02,
     5      0.83365885764441085361189749186838Q-02,
     6      -.10859178246856505601082189739913Q-01/
         data c1o1003/
     1      0.15326839984541555427329998929180Q-01,
     2      -.24442731588961260587692287597818Q-01,
     3      0.49398498217082639400773902910311Q-01,
     4      -.29555991330218355848454684694933Q+00,
     5      -.11501256575009601238973700255592Q-01,
     6      0.52501090689522049815845868394438Q-02,
     7      -.65592823496474080065376011290347Q-02,
     8      0.81708835924187811111221636589173Q-02,
     9      -.10748957803706093808471521246697Q-01,
     a      0.15251162497387337175072698065416Q-01,
     1      -.24388722286015401342184698470797Q-01,
     2      0.49358795096238936288672774275157Q-01,
     3      -.29552904507249907575786040621680Q+00,
     4      0.10725634409598315007401585611766Q-01,
     5      -.42177202717049742175451307036776Q-02,
     6      0.53306445679732415893156910974598Q-02/
         data c1o1004/
     1      -.64392181659145962093852961253704Q-02,
     2      0.80825806350743622661770046125762Q-02,
     3      -.10685862439299945371124346630102Q-01,
     4      0.15204970461837074370150337111824Q-01,
     5      -.24354119721961294332075263727238Q-01,
     6      0.49332394449605704519610117212020Q-01,
     7      -.29550790246373785817016627965340Q+00,
     8      -.10073592767916600120406396152216Q-01,
     9      0.34395967064158808778375871046218Q-02,
     a      -.44385216733460532349565789487354Q-02,
     1      0.52453883965026508560396797021002Q-02,
     2      -.63679935026504114908028208979118Q-02,
     3      0.80294784843846486030803022819925Q-02,
     4      -.10646001194933391436094721179128Q-01,
     5      0.15174533170296142860775867782362Q-01,
     6      -.24330538740079206110315352436476Q-01/
         data c1o1005/
     1      0.49313907955638127936421165673500Q-01,
     2      -.29549276191101668043480205430435Q+00,
     3      0.95160120090258372551813911454481Q-02,
     4      -.28364420771874819575950328118646Q-02,
     5      0.37660115994281529034786096020698Q-02,
     6      -.43804655791693133001094642854736Q-02,
     7      0.51877785803185251992834749268191Q-02,
     8      -.63229937556881727268769103583338Q-02,
     9      0.79948422080631367219350569766812Q-02,
     a      -.10619069369183036017327884919051Q-01,
     1      0.15153341104941315032348838660317Q-01,
     2      -.24313706759759520285442231835161Q-01,
     3      0.49300435863336771472823087932724Q-01,
     4      -.29548153216415021299117070938474Q+00,
     5      -.90324931074751634441158358595799Q-02,
     6      0.23582076234132339046600939635260Q-02/
         data c1o1006/
     1      -.32438594927376254937616535854571Q-02,
     2      0.37296306212585391388640091752294Q-02,
     3      -.43339119793755908520992207689526Q-02,
     4      0.51494584747423292530877536926553Q-02,
     5      -.62927314497961645526132254886837Q-02,
     6      0.79708973985438950925086320769613Q-02,
     7      -.10599954966080632357902895203956Q-01,
     8      0.15137956033438099228733068046176Q-01,
     9      -.24301250034326664913855620367809Q-01,
     a      0.49290300850495866944469820994352Q-01,
     1      -.29547296322466959932430771309789Q+00,
     2      0.86082742561821553003344259330751Q-02,
     3      -.19718784682809894086298897561592Q-02,
     4      0.28286539768145945043264764952458Q-02,
     5      -.32250306379690733534369530740698Q-02,
     6      0.36921938843619767968285525867353Q-02/
         data c1o1007/
     1      -.43011755790503604622808232449235Q-02,
     2      0.51229013163584825783673816305152Q-02,
     3      -.62713551023059747000573298513999Q-02,
     4      0.79536028014437572548943612209989Q-02,
     5      -.10585865536001264974327554262381Q-01,
     6      0.15126412885751517397777518613204Q-01,
     7      -.24291759787010549757810981164007Q-01,
     8      0.49282476006338557723594260772191Q-01,
     9      -.29546626960301328098755225098317Q+00,
     a      -.82323807395014576152433650659833Q-02,
     1      0.16548727133000638491536055456045Q-02,
     2      -.24919214556421094058638653629941Q-02,
     3      0.28242301517623148157259174238440Q-02,
     4      -.31952091013876867071642423433023Q-02,
     5      0.36641804934182089335928570225709Q-02,
     6      -.42777892130189692778149772665648Q-02/
         data c1o1008/
     1      0.51037510561966407380870625914494Q-02,
     2      -.62556620840343366116872662906266Q-02,
     3      0.79406763819472080993718763433833Q-02,
     4      -.10575162702808461585731797542930Q-01,
     5      0.15117518370923809280452484273655Q-01,
     6      -.24284355103642803168207913716583Q-01,
     7      0.49276303082919902258767997369959Q-01,
     8      -.29546093700294505093605830910477Q+00,
     9      0.78964652104781015292310619417043Q-02,
     a      -.13912696574602049026629122150336Q-02,
     1      0.22142778744082727706525649168444Q-02,
     2      -.24994542175534219067714443793219Q-02,
     3      0.28008398736466978703127062364470Q-02,
     4      -.31712315542569439603344876488540Q-02,
     5      0.36435320050647928813343386618906Q-02,
     6      -.42605821749298565876789280553180Q-02/
         data c1o1009/
     1      0.50894740114103884387755837341899Q-02,
     2      -.62437806094062874705722727123579Q-02,
     3      0.79307449953075785049125035952373Q-02,
     4      -.10566830312095390001712504769241Q-01,
     5      0.15110512027349901781119821418490Q-01,
     6      -.24278461248092145934781982204755Q-01,
     7      0.49271343830754026319641993468356Q-01,
     8      -.29545661692266011873742407851564Q+00,
     9      -.75940535662466084953926057660132Q-02,
     a      0.11695527089055288211614490344538Q-02,
     1      -.19821099245858536487282123109296Q-02,
     2      0.22318309933241294379277700564493Q-02,
     3      -.24815471689068151917730328359519Q-02,
     4      0.27803416924606816192059387097607Q-02,
     5      -.31529659059115185336525208563147Q-02,
     6      0.36280324022273618568072877181233Q-02/
         data c1o1010/
     1      -.42475627237616789976194769253217Q-02,
     2      0.50785329039646652212968710096942Q-02,
     3      -.62345555792599190863170054634142Q-02,
     4      0.79229397152397229168780928982412Q-02,
     5      -.10560209382883182762099223705696Q-01,
     6      0.15104889722071643243887546876981Q-01,
     7      -.24273689785647872106982338728437Q-01,
     8      0.49267297015908486935969892823931Q-01,
     9      -.29545306622045206254216406046707Q+00,
     a      0.73200389936713837489142303735460Q-02,
     1      -.98120503905143354231402459621674Q-03,
     2      0.17856010790712875764676170014166Q-02,
     3      -.20081304549640106540584440168972Q-02,
     4      0.22186377582339572169584311249147Q-02,
     5      -.24640720717501688742140147328168Q-02,
     6      0.27641645411494236669313556067671Q-02/
         data c1o1011/
     1      -.31389755653695797766562348442381Q-02,
     2      0.36161348857550794211376299100553Q-02,
     3      -.42374697674563397979831472800098Q-02,
     4      0.50699538849111180609264316655182Q-02,
     5      -.62272415390517465722352566477984Q-02,
     6      0.79166876850641107749729264492264Q-02,
     7      -.10554856588231682867757852460395Q-01,
     8      0.15100306064381228342855125168516Q-01,
     9      -.24269770257717920075287438077848Q-01,
     a      0.49263949873958477247538184592650Q-01,
     1      -.29545011095676464821425419118176Q+00,
     2      -.70703332350568378606214289842643Q-02,
     3      0.81980618328560481570061374995767Q-03,
     4      -.16175100114785663307770123144036Q-02,
     5      0.18188280255879797085468612742494Q-02,
     6      -.19990191600082996994434500745530Q-02/
         data c1o1012/
     1      0.22038074911187474543756651458453Q-02,
     2      -.24497370131851145381428251984618Q-02,
     3      0.27515156149816093428246394394945Q-02,
     4      -.31280832834343072853742131614681Q-02,
     5      0.36068089223342714170027188247279Q-02,
     6      -.42294816303850105793820247787443Q-02,
     7      0.50630959098642885991880754177934Q-02,
     8      -.62213391580070340038159578912367Q-02,
     9      0.79115981609066874585468672511172Q-02,
     a      -.10550464357065706623542747602286Q-01,
     1      0.15096517720795340259286006727352Q-01,
     2      -.24266509485818020293884808620571Q-01,
     3      0.49261148557084559341766032774007Q-01,
     4      -.29544762395864283963106769040307Q+00,
     5      0.68416204993256629015239983149012Q-02,
     6      -.68043387411723853806701977776926Q-03/
         data c1o1013/
     1      0.14723885204458104769478206078740Q-02,
     2      -.16569082529451088222949542519272Q-02,
     3      0.18132743035517751154019892442914Q-02,
     4      -.19865170783453909761684620288416Q-02,
     5      0.21911064406382791810667071358408Q-02,
     6      -.24382863132858953491488442706519Q-02,
     7      0.27415277113062901260941359827099Q-02,
     8      -.31194528813637444249542608695476Q-02,
     9      0.35993617529012774320075992770403Q-02,
     a      -.42230462477907033774599494829589Q-02,
     1      0.50575227881793063394879163500318Q-02,
     2      -.62165033124465966722104176712050Q-02,
     3      0.79073968416998025707377919006704Q-02,
     4      -.10546813632038647191645248889952Q-01,
     5      0.15093349094757430862683047734887Q-01,
     6      -.24263766407157526788308691940846Q-01/
         data c1o1014/
     1      0.49258779517704165513491402595225Q-01,
     2      -.29544551045279669397668268201655Q+00,
     3      -.66311801755270690014099922556938Q-02,
     4      0.55925785090690876774125438020344Q-03,
     5      -.13460658909026116507228406388496Q-02,
     6      0.15171001149029666293878200090330Q-02,
     7      -.16544728836044391602536499924083Q-02,
     8      0.18028330107698977720331665668844Q-02,
     9      -.19752729153985486299729557470240Q-02,
     a      0.21807309970029484265195023517149Q-02,
     1      -.24291152115328136702510213638434Q-02,
     2      0.27335299126057475639250990579408Q-02,
     3      -.31125017071000640729183769251305Q-02,
     4      0.35933176935239942493054847958320Q-02,
     5      -.42177819680657636004296887073366Q-02,
     6      0.50529291239154277421800687354133Q-02/
         data c1o1015/
     1      -.62124889938544604687226756709062Q-02,
     2      0.79038863658779936135323400423867Q-02,
     3      -.10543744825213193865610362407567Q-01,
     4      0.15090670799093363924724130845456Q-01,
     5      -.24261436003781303807358420934339Q-01,
     6      0.49256757427745749778734851482527Q-01,
     7      -.29544369860823453099994688214117Q+00,
     8      0.64367567472957988518192497345040Q-02,
     9      -.45325782023151965204185725652010Q-03,
     a      0.12353000636891743241212566991349Q-02,
     1      -.13953735966256016524695596712431Q-02,
     2      0.15174121166582667188955399521568Q-02,
     3      -.16458647871535502682359350804030Q-02,
     4      0.17928938827753876604247000422420Q-02,
     5      -.19658664639903176006443607153514Q-02,
     6      0.21723002948970973505451292205442Q-02/
         data c1o1016/
     1      -.24216944044594115669812969484035Q-02,
     2      0.27270343044403357278670223284761Q-02,
     3      -.31068200846206696749264343748403Q-02,
     4      0.35883424710478243343745058436263Q-02,
     5      -.42134181349634249372522228225782Q-02,
     6      0.50490957463064537869904090143324Q-02,
     7      -.62091182039681250362490105111347Q-02,
     8      0.79009216649456131312778950798424Q-02,
     9      -.10541139355016319426937843830275Q-01,
     a      0.15088385733608186396240780671877Q-01,
     1      -.24259438746405294396882828938289Q-01,
     2      0.49255017132046061345770009076384Q-01,
     3      -.29544213315737722452443273713672Q+00,
     4      -.62564627642440646174018012320189Q-02,
     5      0.36002370420473312618195150464856Q-03,
     6      -.11375362328720315995429275864750Q-02/
         data c1o1017/
     1      0.12885996084400350872201094090524Q-02,
     2      -.13981172588576522621084839837560Q-02,
     3      0.15104421173579740049124595703373Q-02,
     4      -.16370996472898555708602003482537Q-02,
     5      0.17843640653543463693255883201323Q-02,
     6      -.19581092027357550780156314553486Q-02,
     7      0.21654072254554077908938191072689Q-02,
     8      -.24156180378554128474549021811482Q-02,
     9      0.27216886097188800469489404372001Q-02,
     a      -.31021151237085610863830112505471Q-02,
     1      0.35841959382046395548432830737686Q-02,
     2      -.42097584273719688388784823655224Q-02,
     3      0.50458619398627532781009096364722Q-02,
     4      -.62062589978094646595808197022910Q-02,
     5      0.78983940949185318287912050245923Q-02,
     6      -.10538907564491786661944127922041Q-01/
         data c1o1018/
     1      0.15086419837391977535934600899458Q-01,
     2      -.24257713486891718897551333707892Q-01,
     3      0.49253508163698320166961581061769Q-01,
     4      -.29544077099948220111665349501692Q+00,
     5      0.60887053964189518614844169019448Q-02,
     6      -.27761163885516877661007538976760Q-03,
     7      0.10507363802866913818603476735181Q-02,
     8      -.11943149541754814974394199502994Q-02,
     9      0.12935044193043821000949185318420Q-02,
     a      -.13926170239992127084925611375617Q-02,
     1      0.15027371919231408630836548102775Q-02,
     2      -.16293657609597319212741986303876Q-02,
     3      0.17772213041678465150958341451221Q-02,
     4      -.19517001202322972409157532901839Q-02,
     5      0.21597174466557615141845858971511Q-02,
     6      -.24105840339784276299011054186189Q-02/
         data c1o1019/
     1      0.27172362714016510074331787083171Q-02,
     2      -.30981735863542710717042437475149Q-02,
     3      0.35807020428039803596636958444460Q-02,
     4      -.42066575739453875448190436793387Q-02,
     5      0.50431076027593439157235372939379Q-02,
     6      -.62038118210610603973896968840157Q-02,
     7      0.78962209354306767169505025307034Q-02,
     8      -.10536980611994420009613339543296Q-01,
     9      0.15084715799828654317288806047185Q-01,
     a      -.24256212566544995922366179566323Q-01,
     1      0.49252190928010126971306575949435Q-01,
     2      -.29543957811111536723324088095258Q+00/
         real *16
     1  c1n1001(16),c1n1002( 9)
         real *16
     1  c1n2001(2)
         equivalence
     1  (c1n1001(16),c1n2001(1)),(c1n1002(1),c1n2001(2))
         data c1n1001/
     1      0.39894228040143267793994605993440Q+00,
     2      0.13999113969194818736245921161687Q+01,
     3      0.13570441972843824753695817220635Q+01,
     4      0.13549794836524748256569817096164Q+01,
     5      0.13544710888929506391801599558685Q+01,
     6      0.13542721438886048350051663505889Q+01,
     7      0.13541739305448133870101059587806Q+01,
     8      0.13541181454912999058423221172216Q+01,
     9      0.13540833423912524985183975076296Q+01,
     a      0.13540601214384060900742898369773Q+01,
     1      0.13540438272972654718607541248863Q+01,
     2      0.13540319364481947692896151058610Q+01,
     3      0.13540229824154468111405209011357Q+01,
     4      0.13540160645692229897452355380050Q+01,
     5      0.13540106043405447423602748763262Q+01,
     6      0.13540062159138947960658881233533Q+01/
         data c1n1002/
     1      0.13540026337790742624233226680443Q+01,
     2      0.13539996702313477673437128226840Q+01,
     3      0.13539971894453495935171010251888Q+01,
     4      0.13539950911170805595074984418910Q+01,
     5      0.13539932998352394817224388538775Q+01,
     6      0.13539917579945915078311526169261Q+01,
     7      0.13539904209616437346537805601980Q+01,
     8      0.13539892537087794819172334624556Q+01,
     9      0.13539882284274530312166347988910Q+01/
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
     1      -.25070998606785569047565324945748Q+00,
     2      0.57512045737301469496887323463574Q-01,
     3      -.29966091357091194746036402311300Q+00,
     4      -.38075166250556302152708092576867Q-01,
     5      0.51180561980448622864172282582733Q-01,
     6      -.29654536614120004177680765355072Q+00,
     7      0.29521345867248174316942070337458Q-01,
     8      -.25245640675711909454951696771571Q-01,
     9      0.50019290378078518657954493919732Q-01,
     a      -.29598223249908628357260227597078Q+00,
     1      -.24553321892902567360077567858587Q-01,
     2      0.15564552791084551251735316516925Q-01,
     3      -.24837632647479883172748677200110Q-01,
     4      0.49679411731572120629048255798306Q-01,
     5      -.29576625792317534127763395320888Q+00,
     6      0.21250438170054574659447317679624Q-01/
         data c2o1002/
     1      -.10702828773914904313383647383310Q-01,
     2      0.15539548457395201395047648340047Q-01,
     3      -.24613665400806050928951153140474Q-01,
     4      0.49525361478379116161075907378633Q-01,
     5      -.29565788808087234488523313544884Q+00,
     6      -.18869801556506937626326087366491Q-01,
     7      0.78530286906855986253736026413730Q-02,
     8      -.10898977005291286619147346372123Q-01,
     9      0.15388911250513444087206933684118Q-01,
     a      -.24495957688816936201228518909272Q-01,
     1      0.49441201487129258328506985721077Q-01,
     2      -.29559505614174383744358239483824Q+00,
     3      0.17058811648166671297338635499160Q-01,
     4      -.60139084151962929417547917413420Q-02,
     5      0.81869271600925051109716487415520Q-02,
     6      -.10798433305087110962763609940897Q-01/
         data c2o1003/
     1      0.15296970264983866266710499906030Q-01,
     2      -.24426815388319371530403365898083Q-01,
     3      0.49389808375965656476673300062069Q-01,
     4      -.29555511926287253331375328714695Q+00,
     5      -.15626924783911591326595653468065Q-01,
     6      0.47461928084031604941134595602552Q-02,
     7      -.64378990273593543219836913770933Q-02,
     8      0.81223505525813977465176564894655Q-02,
     9      -.10725725334021131863346580516152Q-01,
     a      0.15239284796109844560608081159798Q-01,
     1      -.24382674973403876282414875805766Q-01,
     2      0.49355973047138471883806980282255Q-01,
     3      -.29552805067429026024384568162472Q+00,
     4      0.14461433172892855271830958437551Q-01,
     5      -.38293412342137203506017880846949Q-02,
     6      0.52310187820525450719090251372172Q-02/
         data c2o1004/
     1      -.64000943221631212375639379023883Q-02,
     2      0.80644757912280466051031386941143Q-02,
     3      -.10677118482782346454088534743450Q-01,
     4      0.15200976964133169489845534936826Q-01,
     5      -.24352715413795018755951836906480Q-01,
     6      0.49332446724830186004900659648805Q-01,
     7      -.29550880784101047025112577614099Q+00,
     8      -.13491062625513339550448317800546Q-01,
     9      0.31415201628643109443333778408591Q-02,
     a      -.43560713178685161513271337770677Q-02,
     1      0.52136967174382060571303686800530Q-02,
     2      -.63539444873645079888045382915080Q-02,
     3      0.80232210931081411543218659278078Q-02,
     4      -.10643645651832556817340832121611Q-01,
     5      0.15174264844299174760058766539101Q-01,
     6      -.24331412164797108390967270763523Q-01/
         data c2o1005/
     1      0.49315393303559195695621210590795Q-01,
     2      -.29549461353030524413095414823542Q+00,
     3      0.12668357601464749898051046803778Q-01,
     4      -.26104252122969693748782226058162Q-02,
     5      0.36974025719993597219480052257332Q-02,
     6      -.43547629764832539999034194806344Q-02,
     7      0.51769986463606340538666635226866Q-02,
     8      -.63187435712958433057386781364427Q-02,
     9      0.79938130722518029686495804056027Q-02,
     a      -.10619725724486491674799110842932Q-01,
     1      0.15154888605161657282450876975979Q-01,
     2      -.24315701452949671619823438309231Q-01,
     3      0.49302619133142359061619753840045Q-01,
     4      -.29548382952021061469306976759749Q+00,
     5      -.11960399253324328097081434983421Q-01,
     6      0.21906820797622132899452806789992Q-02/
         data c2o1006/
     1      -.31865962578080049994617593849488Q-02,
     2      0.37088343139407781994697259521743Q-02,
     3      -.43258082237806341462298395328184Q-02,
     4      0.51468515807994993670204963994341Q-02,
     5      -.62927895452962441598184135145178Q-02,
     6      0.79723136744096924112837921672388Q-02,
     7      -.10602059017398512934165847479019Q-01,
     8      0.15140376755987876191132997574970Q-01,
     9      -.24303770141166395128156246127867Q-01,
     a      0.49292792260314714138588995688210Q-01,
     1      -.29547543600533194760108327245751Q+00,
     2      0.11343576497026311034676545001935Q-01,
     3      -.18525118825486727597815498487512Q-02,
     4      0.27808244196674673202236819519036Q-02,
     5      -.32083094591041748680372603964288Q-02,
     6      0.36863096997241700014693944255196Q-02/
         data c2o1007/
     1      -.42999301589177507514822001178240Q-02,
     2      0.51238598311954114945315366087894Q-02,
     3      -.62734009828327201660426198515785Q-02,
     4      0.79561686026563307130755520661425Q-02,
     5      -.10588640995717806666712691361164Q-01,
     6      0.15129218647943709522158587407151Q-01,
     7      -.24294487873609523100955952679407Q-01,
     8      0.49285064227868875881669593210709Q-01,
     9      -.29546876992425296788656920970128Q+00,
     a      -.10800477180968630344547417153922Q-01,
     1      0.15756390213810354841877396956315Q-02,
     2      -.24520313009260034425640603408529Q-02,
     3      0.28109337022152414630727943511003Q-02,
     4      -.31911863549941898507935285320266Q-02,
     5      0.36640746153578211314111005463513Q-02,
     6      -.42795001119203409663576518196981Q-02/
         data c2o1008/
     1      0.51063223458695811046213211872927Q-02,
     2      -.62586128615548974727133811667134Q-02,
     3      0.79437474853972851557831174184946Q-02,
     4      -.10578207083308057027859313160984Q-01,
     5      0.15120449066782547843898483716668Q-01,
     6      -.24287119301975378690248271571302Q-01,
     7      0.49278872916948613806329436556029Q-01,
     8      -.29546338438682802828866310154903Q+00,
     9      0.10317955690937434313736342174460Q-01,
     a      -.13458275138677170072375513485276Q-02,
     1      0.21811395256441041201415495829691Q-02,
     2      -.24890651182309551769682806972349Q-02,
     3      0.27983932289087490290739543052347Q-02,
     4      -.31720878938517538291464196813417Q-02,
     5      0.36458764828509365680930514044756Q-02,
     6      -.42635947302251800531871955221221Q-02/
         data c2o1009/
     1      0.50927472537358017074539469020246Q-02,
     2      -.62470986809311429933251066033721Q-02,
     3      0.79339882252795198921071315330013Q-02,
     4      -.10569930052121374729681066140315Q-01,
     5      0.15113427735114292869256676930856Q-01,
     6      -.24281168781441248105317865752166Q-01,
     7      0.49273833490853314204373218381695Q-01,
     8      -.29545896895014412288056189731962Q+00,
     9      -.98858860923732126856605227595853Q-02,
     a      0.11528177492604036069711603798647Q-02,
     1      -.19547656602754354462657195555350Q-02,
     2      0.22239317824488179857006789984612Q-02,
     3      -.24804461702942427384056919935729Q-02,
     4      0.27820167576087646329481100138677Q-02,
     5      -.31558475313911394531528344464644Q-02,
     6      0.36314176360321325792869070360879Q-02/
         data c2o1010/
     1      -.42511072051829392083735541240352Q-02,
     2      0.50820578036174578315486381418204Q-02,
     3      -.62379645618129686774441329461653Q-02,
     4      0.79261797626854479856197808187331Q-02,
     5      -.10563250945149290548247192339473Q-01,
     6      0.15107716158842005990829688110848Q-01,
     7      -.24276292065206184687443432245123Q-01,
     8      0.49269675042632025719792826673293Q-01,
     9      -.29545530231616656547553261719552Q+00,
     a      0.94963315826431238747689544020452Q-02,
     1      -.98904812877758878316982420886227Q-03,
     2      0.17632686566611160160097938269043Q-02,
     3      -.20023800744972146180355418125630Q-02,
     4      0.22186940382833363455671189057051Q-02,
     5      -.24664484392617754900879038874498Q-02,
     6      0.27675042453214967336103919717642Q-02/
         data c2o1011/
     1      -.31426770828520094399406506276497Q-02,
     2      0.36199083453409996045041578756615Q-02,
     3      -.42411682648893107945854018320828Q-02,
     4      0.50735011269461186616189089007287Q-02,
     5      -.62305978740702819815647251893297Q-02,
     6      0.79198329959872647210394312518811Q-02,
     7      -.10557781186993397153563255456637Q-01,
     8      0.15103005703001018340065240431849Q-01,
     9      -.24272243837109777058494196900494Q-01,
     a      0.49266202243387543910053776454295Q-01,
     1      -.29545222317219651788305151620903Q+00,
     2      -.91429758430838716768968919786603Q-02,
     3      0.84883588657010082372455709801968Q-03,
     4      -.15995433026190698734194705092630Q-02,
     5      0.18149447625937951740370637400018Q-02,
     6      -.20000771171573768015052475772656Q-02/
         data c2o1012/
     1      0.22067881081017666739266223335910Q-02,
     2      -.24534693680751412999183921432508Q-02,
     3      0.27554866717154839565583072259663Q-02,
     4      -.31320506229206012967391431942297Q-02,
     5      0.36106533707677339304072597622229Q-02,
     6      -.42331442375077129179857353613590Q-02,
     7      0.50665485280380925372332658365937Q-02,
     8      -.62245697423097867322363468322681Q-02,
     9      0.79146030120712533487674853667458Q-02,
     a      -.10553243846374309789588105639105Q-01,
     1      0.15099073923795022389469566506650Q-01,
     2      -.24268845345500978632538230835532Q-01,
     3      0.49263271292717260147268093711122Q-01,
     4      -.29544961161540670218858655863035Q+00,
     5      0.88207242651927416946375169342741Q-02,
     6      -.72783662432271278393503146722615Q-03/
         data c2o1013/
     1      0.14582484796048271397091206443357Q-02,
     2      -.16546571321181950026707800421242Q-02,
     3      0.18152041870330678392047563849456Q-02,
     4      -.19900210096951354718864859150135Q-02,
     5      0.21951768571083009655878856951023Q-02,
     6      -.24424878756553266110047174260867Q-02,
     7      0.27456595939657646050944319826834Q-02,
     8      -.31234201579008965647459746218120Q-02,
     9      0.36031205839342754946672504146355Q-02,
     a      -.42265784562479855996848468406870Q-02,
     1      0.50608232733989493558202470466833Q-02,
     2      -.62195734894214192611050213120121Q-02,
     3      0.79102410697667615020937076768941Q-02,
     4      -.10549437204036321756877307048788Q-01,
     5      0.15095757143667362675244066143426Q-01,
     6      -.24265963779526013994166989067280Q-01/
         data c2o1014/
     1      0.49260774371008300879575627391877Q-01,
     2      -.29544737697721965340423731069452Q+00,
     3      -.85254183460852500103792358836188Q-02,
     4      0.62267851972305899704284016561199Q-03,
     5      -.13352987827628732702698450608530Q-02,
     6      0.15162834263133561096654210782021Q-02,
     7      -.16571655820547595057630495551955Q-02,
     8      0.18067922186492958901576540405504Q-02,
     9      -.19796355055365971983891886334626Q-02,
     a      0.21851302475426375256272630939283Q-02,
     1      -.24333869679111707820285270230108Q-02,
     2      0.27376005764876229505480930098633Q-02,
     3      -.31163406868081009024665069430848Q-02,
     4      0.35969155106130312959518703495656Q-02,
     5      -.42211395364001081345792982878604Q-02,
     6      0.50560522440920887152170105347622Q-02/
         data c2o1015/
     1      -.62153854145080106858360307566367Q-02,
     2      0.79065641597639842687742533100691Q-02,
     3      -.10546211469285824847030737519406Q-01,
     4      0.15092932709265176859988535271615Q-01,
     5      -.24263498750876259216734855944673Q-01,
     6      0.49258629316644404448724363041212Q-01,
     7      -.29544544964018395388104316410271Q+00,
     8      0.82536272759282890643586270407776Q-02,
     9      -.53070912109770371547839520168450Q-03,
     a      0.12275208882512416068611922542861Q-02,
     1      -.13958236737581230497610665445745Q-02,
     2      0.15207751746720034948902873053193Q-02,
     3      -.16502216593979737471667976075229Q-02,
     4      0.17975098097795230877435674444617Q-02,
     5      -.19704356384647978333119078339152Q-02,
     6      0.21766910630039534604024007790827Q-02/
         data c2o1016/
     1      -.24258520282493485443286985292943Q-02,
     2      0.27309398709854753196242065629899Q-02,
     3      -.31104717497289407361924503370412Q-02,
     4      0.35917463922012975775810348631298Q-02,
     5      -.42165838991161367300200094867974Q-02,
     6      0.50520339909708902523681654343981Q-02,
     7      -.62118393470761925118291127333060Q-02,
     8      0.79034352079561255211379301185236Q-02,
     9      -.10543453495367364879418757950052Q-01,
     a      0.15090507209563366384614978282356Q-01,
     1      -.24261373209723461498864603230473Q-01,
     2      0.49256772609252096872921919929289Q-01,
     3      -.29544377541488128333789450800611Q+00,
     4      -.80024932931278943378647892773677Q-02,
     5      0.44981653337355524164029066083877Q-03,
     6      -.11324159367490077544495261424448Q-02/
         data c2o1017/
     1      0.12901732488039377749530796534884Q-02,
     2      -.14020717930596169498308278326469Q-02,
     3      0.15151475476665757823106733894211Q-02,
     4      -.16419358376772884874454842191403Q-02,
     5      0.17890795395433726440418166659362Q-02,
     6      -.19626012478482209757664883562626Q-02,
     7      0.21696378618583013741563771903355Q-02,
     8      -.24195786994473656568425278082627Q-02,
     9      0.27253841881609261717463132866111Q-02,
     a      -.31055563753367597241951132397192Q-02,
     1      0.35873957413730952567295983213132Q-02,
     2      -.42127298991634030896567067321224Q-02,
     3      0.50486174796016711255509790853179Q-02,
     4      -.62088097987533859322754785653383Q-02,
     5      0.79007498911872464872708718252793Q-02,
     6      -.10541076501981059161237222216556Q-01/
         data c2o1018/
     1      0.15088408468300051983563144876567Q-01,
     2      -.24259527204130971882105650144323Q-01,
     3      0.49255154499450033837417693447414Q-01,
     4      -.29544231159881722916046141900401Q+00,
     5      0.77696152010045088672807549009522Q-02,
     6      -.37830085748887489388931612193451Q-03,
     7      0.10479918688634179020956165664406Q-02,
     8      -.11968890572039049534539931766876Q-02,
     9      0.12979827030851238428606326358980Q-02,
     a      -.13976289042953412089853426323042Q-02,
     1      0.15077653203368652198560252384507Q-02,
     2      -.16342073272036201936524487663500Q-02,
     3      0.17817994769260208400025017811007Q-02,
     4      -.19559918749954750309463200623760Q-02,
     5      0.21637234250948246416739614572194Q-02,
     6      -.24143150991982960407895636539179Q-02/
         data c2o1019/
     1      0.27207072345308353307364710217108Q-02,
     2      -.31014002047234571894692476220713Q-02,
     3      0.35836995773906341796762017738309Q-02,
     4      -.42094401102408249661561948909245Q-02,
     5      0.50456877474696210913037079190236Q-02,
     6      -.62062005785781098690710211829677Q-02,
     7      0.78984276635215535091609471519525Q-02,
     8      -.10539013014192048382786375716531Q-01,
     9      0.15086579987726264386737939139702Q-01,
     a      -.24257913507522827446149168672146Q-01,
     1      0.49253735564627235710216179411966Q-01,
     2      -.29544102416471463430715199336309Q+00/
         real *16
     1  c2n1001(16),c2n1002( 9)
         real *16
     1  c2n2001(2)
         equivalence
     1  (c2n1001(16),c2n2001(1)),(c2n1002(1),c2n2001(2))
         data c2n1001/
     1      0.44101277172455148218908093217173Q+00,
     2      0.13678375475512571695021963113044Q+01,
     3      0.13556744552377051951970767752428Q+01,
     4      0.13546724120646797895047007644671Q+01,
     5      0.13543679291229456362679808729000Q+01,
     6      0.13542309448069017582192585232623Q+01,
     7      0.13541563992446176794653508823721Q+01,
     8      0.13541109494738953008795200758165Q+01,
     9      0.13540810429324468379480140776831Q+01,
     a      0.13540602497085407826147008847001Q+01,
     1      0.13540451753385363185206415835213Q+01,
     2      0.13540338811435611816239580243495Q+01,
     3      0.13540251906723475410109195682646Q+01,
     4      0.13540183546817441001975311456410Q+01,
     5      0.13540128767862101358886689873369Q+01,
     6      0.13540084171138025791474099220596Q+01/
         data c2n1002/
     1      0.13540047363697705638346593653999Q+01,
     2      0.13540016619836837017367177905028Q+01,
     3      0.13539990668830587745054453289147Q+01,
     4      0.13539968557705658039211119503815Q+01,
     5      0.13539949560120545135497667326778Q+01,
     6      0.13539933114420046100327708174745Q+01,
     7      0.13539918780632689748886449179464Q+01,
     8      0.13539906210054152898064971754202Q+01,
     9      0.13539895123367322817281142661907Q+01/
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
     1      -.32610039095003236362680845552716Q+00,
     2      0.83525956091866291472108033017387Q-01,
     3      -.29761405535875304046908911865794Q+00,
     4      -.53964371739008475257088312523792Q-01,
     5      0.50592237433966623998060836882519Q-01,
     6      -.29652003171832477082556514884883Q+00,
     7      0.41024268401731005258654685230476Q-01,
     8      -.25184733128379806989624177950062Q-01,
     9      0.50085473424096489530914430404508Q-01,
     a      -.29606085899559042744795767369759Q+00,
     1      -.33597786752511836817237098034158Q-01,
     2      0.15755732393737850622440665374710Q-01,
     3      -.24946186070838164189236092624550Q-01,
     4      0.49766106423249197763697464883794Q-01,
     5      -.29583876929652110809598702195064Q+00,
     6      0.28718644868595452502488833803639Q-01/
         data c3o1002/
     1      -.11031082545801275079883173737433Q-01,
     2      0.15670090059428654943385955283319Q-01,
     3      -.24705336173552110664707330505604Q-01,
     4      0.49595450639155539268992135375539Q-01,
     5      -.29571549373500087778407083964791Q+00,
     6      -.25239607389285055654588535032567Q-01,
     7      0.82614240652639456382258814423130Q-02,
     8      -.11041430521852461574886412870585Q-01,
     9      0.15482696513220455688258674513218Q-01,
     a      -.24565506309061420369829962884621Q-01,
     1      0.49495460658082982149565324561781Q-01,
     2      -.29564015863451279014251159625497Q+00,
     3      0.22618452848441220083657437665718Q-01,
     4      -.64710929960297297993940776367897Q-02,
     5      0.83358063718111831503852678507299Q-02,
     6      -.10892707340616588419758041493021Q-01/
         data c3o1003/
     1      0.15365626582144274005367654504265Q-01,
     2      -.24479845457075027468622375416636Q-01,
     3      0.49432125622536578362025045690304Q-01,
     4      -.29559080522705707577582892361990Q+00,
     5      -.20563799955945751772796457288165Q-01,
     6      0.52335892320125146332463793980780Q-02,
     7      -.65900170522175585067700440222384Q-02,
     8      0.82161904412647058025280130656575Q-02,
     9      -.10793221336222519202234366128982Q-01,
     a      0.15291141375352548399907072081755Q-01,
     1      -.24423942625620078628692933343576Q-01,
     2      0.49389539569004170788063109583164Q-01,
     3      -.29555673378257307542767798734045Q+00,
     4      0.18904321444991723780659718215880Q-01,
     5      -.43353435917149778841650381627846Q-02,
     6      0.53844393167684380737330326606593Q-02/
         data c3o1004/
     1      -.64929716786521259950396344539154Q-02,
     2      0.81306677640422240856575938565672Q-02,
     3      -.10727781527900395496771028908357Q-01,
     4      0.15241283447444744379739099825709Q-01,
     5      -.24385518895713969127463369077657Q-01,
     6      0.49359547894792857073858295752218Q-01,
     7      -.29553223174822924283959604256425Q+00,
     8      -.17532317478629469933578885669316Q-01,
     9      0.36585866804834165833411379705997Q-02,
     a      -.45095811340257713964616736725229Q-02,
     1      0.53053112658929051943615605439049Q-02,
     2      -.64187714676827351406959132711365Q-02,
     3      0.80726919596286821444310504778314Q-02,
     4      -.10683001000787959536345655415623Q-01,
     5      0.15206363776416039601619889003344Q-01,
     6      -.24357993165262763310215727294041Q-01/
         data c3o1005/
     1      0.49337635428247568511101143933651Q-01,
     2      -.29551402530797923390637268419390Q+00,
     3      0.16376493257666488126207951859985Q-01,
     4      -.31335074803255292938410815090580Q-02,
     5      0.38502266269969259154349255258029Q-02,
     6      -.44449502655940714636493253325012Q-02,
     7      0.52404498267032613124441857915386Q-02,
     8      -.63670452929705130170425480246696Q-02,
     9      0.80322329352558598643554984835868Q-02,
     a      -.10651118716584354438045191230361Q-01,
     1      0.15180976777308250122186728697562Q-01,
     2      -.24337604852941824601624252153571Q-01,
     3      0.49321141455714343737759938539520Q-01,
     4      -.29550012936507634160039253654555Q+00,
     5      -.15387666182556627839478138700043Q-01,
     6      0.27163246388959069470463052827444Q-02/
         data c3o1006/
     1      -.33382349996314064938293166801199Q-02,
     2      0.37975126678015958725973727945952Q-02,
     3      -.43879022741078017438841852473179Q-02,
     4      0.51940209029920522922572755985959Q-02,
     5      -.63303000920177450734721444274607Q-02,
     6      0.80030052991928914196513470056855Q-02,
     7      -.10627638318898476771150356824373Q-01,
     8      0.15161946409743796638552802808705Q-01,
     9      -.24322083828291313351617457927852Q-01,
     a      0.49308416391164388739993674391667Q-01,
     1      -.29548928361012790323267232731478Q+00,
     2      0.14530745216508297250314526878368Q-01,
     3      -.23783108934495655116899294118847Q-02,
     4      0.29309566061677314764804782971190Q-02,
     5      -.32954484513263067702346030948341Q-02,
     6      0.37470821208745715517949918473180Q-02/
         data c3o1007/
     1      -.43460111974508739681191671124381Q-02,
     2      0.51604942036851899461121725266861Q-02,
     3      -.63034053948949582197635425084886Q-02,
     4      0.79812311205529188933284566321535Q-02,
     5      -.10609850322691969310717557333610Q-01,
     6      0.15147314412301601297178754616610Q-01,
     7      -.24309995221528752620979551273514Q-01,
     8      0.49298393663514952249481271603650Q-01,
     9      -.29548065674245449120462007705901Q+00,
     a      -.13779981157925717669909718822091Q-01,
     1      0.20999038219225633183867457477380Q-02,
     2      -.26004537771658527314457773111423Q-02,
     3      0.28965348810779254975332485597322Q-02,
     4      -.32506817736151947074201796021318Q-02,
     5      0.37091147091773250909412650250310Q-02,
     6      -.43152950113012021575963720031153Q-02/
         data c3o1008/
     1      0.51356602796157076880536017687360Q-02,
     2      -.62831608475626906415830662930770Q-02,
     3      0.79645794174875986707864469795258Q-02,
     4      -.10596052151610251550559826997014Q-01,
     5      0.15135821888378846857500640195883Q-01,
     6      -.24300396390562478938110888854357Q-01,
     7      0.49290358924358545423680384070986Q-01,
     8      -.29547368247568910973222667448855Q+00,
     9      0.13116025985881882926052175253679Q-01,
     a      -.18673619894887004116827966002092Q-02,
     1      0.23277287476045442331288810857231Q-02,
     2      -.25731499425606835001472713014711Q-02,
     3      0.28566606917173289271283094595763Q-02,
     4      -.32161350543386630804588571666113Q-02,
     5      0.36808698302824269870334177048328Q-02,
     6      -.42922906308961891812555365715077Q-02/
         data c3o1009/
     1      0.51167901079940307990952799225241Q-02,
     2      -.62675466338089337542873516243816Q-02,
     3      0.79515598082070486065975296695085Q-02,
     4      -.10585132613311680234075451609485Q-01,
     5      0.15126630123268388064797748633768Q-01,
     6      -.24292647449028100302603647584619Q-01,
     7      0.49283819173305154591529768433514Q-01,
     8      -.29546796427030158102309481974882Q+00,
     9      -.12524042838735688627052059039938Q-01,
     a      0.16707751387920365172166885504813Q-02,
     1      -.20994524761767696304827262833486Q-02,
     2      0.23065339276971056226810588236706Q-02,
     3      -.25375362694825355364626359884921Q-02,
     4      0.28251181562497781015396528203843Q-02,
     5      -.31900770104838236784917792029678Q-02,
     6      0.36594976620193196782550950213108Q-02/
         data c3o1010/
     1      -.42746587367134235665127897537346Q-02,
     2      0.51021227055537137716489370109385Q-02,
     3      -.62552502174409310665824222924959Q-02,
     4      0.79411861633051629369423741610580Q-02,
     5      -.10576342082759853673032710378358Q-01,
     6      0.15119162976898100993691500031391Q-01,
     7      -.24286301668154583492571435825557Q-01,
     8      0.49278425281011693913233368587613Q-01,
     9      -.29546321774844293614867812384947Q+00,
     a      0.11992451847781898963673069674939Q-01,
     1      -.15028322912877192847929188403469Q-02,
     2      0.19060214184928327063153510125264Q-02,
     3      -.20835404551852854633852278025650Q-02,
     4      0.22746571668183347890184834084910Q-02,
     5      -.25086495481096401440603252034396Q-02,
     6      0.28010064426127802314469188028790Q-02/
         data c3o1011/
     1      -.31701678265839477153589621413915Q-02,
     2      0.36429848476971920400004508006399Q-02,
     3      -.42608559418802535102039685463159Q-02,
     4      0.50904959478837526098721246544928Q-02,
     5      -.62453920976150299714575894714648Q-02,
     6      0.79327856640743809112269676022517Q-02,
     7      -.10569160125472913527621015828014Q-01,
     8      0.15113014087740765520970306341362Q-01,
     9      -.24281039504190240315158354847976Q-01,
     a      0.49273924395852496571278539795693Q-01,
     1      -.29545923467348918177945004596597Q+00,
     2      -.11512074634807736289912364251915Q-01,
     3      0.13580327788793446614171228016040Q-02,
     4      -.17403563493920234674490169832556Q-02,
     5      0.18947082800128371918750871681602Q-02,
     6      -.20549624286210568997925098670978Q-02/
         data c3o1012/
     1      0.22481322454387840319486393328896Q-02,
     2      -.24862792970792590623833365415423Q-02,
     3      0.27824143737722292589606877585653Q-02,
     4      -.31546696079874965822733270620599Q-02,
     5      0.36299726920109596663815616018581Q-02,
     6      -.42498484101608825981475559320515Q-02,
     7      0.50811218563238171918447475935575Q-02,
     8      -.62373659603132785177848158138447Q-02,
     9      0.79258866994399092387061160008329Q-02,
     a      -.10563216289171786760996137907982Q-01,
     1      0.15107890178832307764809908328559Q-01,
     2      -.24276627489095152040031313803689Q-01,
     3      0.49270129754337218047617507514962Q-01,
     4      -.29545585973624721090789576421641Q+00,
     5      0.11075536531560948118737151129015Q-01,
     6      -.12321658964728256240086932659354Q-02/
         data c3o1013/
     1      0.15971341632205425767398853248761Q-02,
     2      -.17330705380388098298890861705523Q-02,
     3      0.18690589984025340183041614737236Q-02,
     4      -.20305491351894771576440462455360Q-02,
     5      0.22273277227172577225136053005373Q-02,
     6      -.24688779602788666608127617463289Q-02,
     7      0.27678389800115579006149082230144Q-02,
     8      -.31423817877836963368764270591437Q-02,
     9      0.36195376393935189128531990652768Q-02,
     a      -.42409272349365975214026818985349Q-02,
     1      0.50734519675111886898448649920031Q-02,
     2      -.62307430714130236318620758977668Q-02,
     3      0.79201507742534966142346749712372Q-02,
     4      -.10558241057216633495768711585220Q-01,
     5      0.15103575278154433575324324823899Q-01,
     6      -.24272891828131934162133309748078Q-01/
         data c3o1014/
     1      0.49266900985186234771526812222931Q-01,
     2      -.29545297519688947042736953006124Q+00,
     3      -.10676839932368929767469511093746Q-01,
     4      0.11219587905975500620500220890869Q-02,
     5      -.14722818763310794467592985261860Q-02,
     6      0.15933939197403273957013674029600Q-02,
     7      -.17100350562690338496233851403237Q-02,
     8      0.18465428804836019373793428646022Q-02,
     9      -.20111586167396997974513417814787Q-02,
     a      0.22110070637024646309196886703205Q-02,
     1      -.24551445681425676652803204442501Q-02,
     2      0.27562161482154715299322421896793Q-02,
     3      -.31324763020419363076546545575680Q-02,
     4      0.36110395311852769026599313950344Q-02,
     5      -.42335945208088099696482833310168Q-02,
     6      0.50670953044384386184107747067316Q-02/
         data c3o1015/
     1      -.62252133705262255860595346515751Q-02,
     2      0.79153298221178495306176591794544Q-02,
     3      -.10554034538357562765040358024356Q-01,
     4      0.15099907483409322826744506580996Q-01,
     5      -.24269700985566336989908801347898Q-01,
     6      0.49264130936587150663402227587848Q-01,
     7      -.29545049048516104244667163035315Q+00,
     8      0.10311054003952757118869790774723Q-01,
     9      -.10248327299673213009201644341976Q-02,
     a      0.13626346242223143753405485955947Q-02,
     1      -.14716779852471350886406303510401Q-02,
     2      0.15727021799285897477005853893373Q-02,
     3      -.16892310322225354374087398738324Q-02,
     4      0.18284345822970426327972605668254Q-02,
     5      -.19958223311718279305189577953188Q-02,
     6      0.21980442705390330850548218596904Q-02/
         data c3o1016/
     1      -.24441336003791751404526750311085Q-02,
     2      0.27468010559287287416910714559633Q-02,
     3      -.31243731235862953097407001605925Q-02,
     4      0.36040248926722415110112122430676Q-02,
     5      -.42274926504196404107946250038998Q-02,
     6      0.50617671318457828960505186317443Q-02,
     7      -.62205481124653276810383918013284Q-02,
     8      0.79112386675319845351247558399560Q-02,
     9      -.10550445930065041402342727453280Q-01,
     a      0.15096763486382091305932624008586Q-01,
     1      -.24266953924366886218659648385759Q-01,
     2      0.49261736669503656782381348215882Q-01,
     3      -.29544833501043570406221061053790Q+00,
     4      -.99740851921482763021473047799947Q-02,
     5      0.93873096378135627157198989374502Q-03,
     6      -.12656991767766573546093884655179Q-02/
         data c3o1017/
     1      0.13648170522650435054418949724696Q-02,
     2      -.14530968730371272181731545547612Q-02,
     3      0.15534495229657192044996478968163Q-02,
     4      -.16722898495854961687954069401703Q-02,
     5      0.18139980029544696497108791635828Q-02,
     6      -.19835668487388261966562157088024Q-02,
     7      0.21875975632334705221550290740857Q-02,
     8      -.24351732463542115923644400113423Q-02,
     9      0.27390665689965958886420017608500Q-02,
     a      -.31176579886296146324614743988523Q-02,
     1      0.35981656896693057837961791562448Q-02,
     2      -.42223594876093959640861185372305Q-02,
     3      0.50572560962860743730351863021821Q-02,
     4      -.62165755311565149894174449653389Q-02,
     5      0.79077367985006830237734047550509Q-02,
     6      -.10547359719317899969334297152405Q-01/
         data c3o1018/
     1      0.15094048030607674185214048968500Q-01,
     2      -.24264571995118740421900939038494Q-01,
     3      0.49259653176522757666462758512843Q-01,
     4      -.29544645309739904421893573006489Q+00,
     5      0.96625049011437392590083784026265Q-02,
     6      -.86199492870422133998027801019200Q-03,
     7      0.11794870699067452648446807390756Q-02,
     8      -.12703666032642350829572651748870Q-02,
     9      0.13481441157365231228942982085134Q-02,
     a      -.14352552090543054626591322880039Q-02,
     1      0.15375743974001999084482188477437Q-02,
     2      -.16586782083841036226615073141540Q-02,
     3      0.18023935466136784988506428729668Q-02,
     4      -.19736416753602498499276772997064Q-02,
     5      0.21790595352895706672931246275653Q-02,
     6      -.24277831447953913510578243828220Q-02/
         data c3o1019/
     1      0.27326331990442090216639335442685Q-02,
     2      -.31120291979019769850123931758173Q-02,
     3      0.35932200038031733279291819473146Q-02,
     4      -.42179993435308175686103367731683Q-02,
     5      0.50534026233277414845061635965820Q-02,
     6      -.62131645956439701557763893734458Q-02,
     7      0.79047160204162529414742620428600Q-02,
     8      -.10544686208199767330259983472687Q-01,
     9      0.15091686583333803370178653872389Q-01,
     a      -.24262493235006732678538431906539Q-01,
     1      0.49257828922792046092488590114820Q-01,
     2      -.29544480035839374910240454153471Q+00/
         real *16
     1  c3n1001(16),c3n1002( 9)
         real *16
     1  c3n2001(2)
         equivalence
     1  (c3n1001(16),c3n2001(1)),(c3n1002(1),c3n2001(2))
         data c3n1001/
     1      0.47912228281421502510918408085686Q+00,
     2      0.13579033917534436382682716702425Q+01,
     3      0.13553948358591440929811547108701Q+01,
     4      0.13547251920816213488089328831358Q+01,
     5      0.13544384959683960629589261665113Q+01,
     6      0.13542896652957407909907389065094Q+01,
     7      0.13542026311181964497291807222554Q+01,
     8      0.13541473840379199476604960349200Q+01,
     9      0.13541101419524155342335321487167Q+01,
     a      0.13540838554214131751846118570777Q+01,
     1      0.13540646148687859834493308592667Q+01,
     2      0.13540501106618174428032128160893Q+01,
     3      0.13540389072770602552119570237658Q+01,
     4      0.13540300744481664439466799426910Q+01,
     5      0.13540229878937679191033282889822Q+01,
     6      0.13540172160466105845469452348296Q+01/
         data c3n1002/
     1      0.13540124528427575360649716780782Q+01,
     2      0.13540084763470657436068071529793Q+01,
     3      0.13540051224584637966879483128890Q+01,
     4      0.13540022677261654445107967202063Q+01,
     5      0.13539998178399636773024527021420Q+01,
     6      0.13539976997502905017529647244580Q+01,
     7      0.13539958561664734054674222078038Q+01,
     8      0.13539942416469645079882664575353Q+01,
     9      0.13539928197760943959369171001150Q+01/
        complex *32 zinv
c 
c       construct the orthogonal polynomials of orders 0, 4, 8, ...
c 
        done=1
        zinv=done/z
        n4=n/4+1
        index=0
        call orteva(zinv,n4,c0o1001,c0n1001,index,w)
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
        call orteva(zinv,n4,c1o1001,c1n1001,index,w)
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
        call orteva(zinv,n4,c2o1001,c2n1001,index,w)
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
        call orteva(zinv,n4,c3o1001,c3n1001,index,w)
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
  
