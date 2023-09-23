        implicit real *8 (a-h,o-z)
        complex *16  pol(500 000),work(410000),
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
        implicit real *8 (a-h,o-z)
        save
        complex *16 pol(n*4,m)
        dimension t(8000),w(8000),x(8000),y(8000),rea(2),
     1   weights(8000),rlens(8000),corners(2,4),rnorms(2000),
     2      promax(2000)
        complex *16 com,work(m,m),
     1  ima,zs(8000),ww(8000),zs2(8000)
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
cccc        call prin2('zs2=*',zs2,n*8)
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
        call dpsquain(zs(i),m,ww)
        do 1600 j=1,m
        pol(i,j)=ww(j)
 1600 continue
 1700 continue
  
         call prin2('after dpsquain, pol(m)=*',pol(1,m-1),8*n)
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
cccc        call prin2('and the matrix of inner products*',work,m*m*2)
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
  
c 
c       test the projections
c 
        do 3200 i=1,n
         call prinf('i=*',i,1)
        call projtest(pol,n*4,m,zs2(i),zs,weights,promax)
 3200 continue
        call prin2('and promax=*',promax,m)
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine projtest(pol,n4,m,z0,zs,weights,promax)
        implicit real *8 (a-h,o-z)
        save
        complex *16 pol(n4,1),z0,zs(1),proj(8000),
     1      f(8000),cd
        dimension weights(1),promax(1)
c 
c       construct the potential of the charge on the
c       boundary
c 
        done=1
        do 1400 i=1,n4
        f(i)=done/(zs(i)-z0)
 1400 continue
cccc        call prin2('f as constructed*',f,n4*2)
c 
c       evaluate the projections of the potentials on all
c       orthonormal polynomials
c 
        do 2400 i=1,m
        cd=0
        do 2200 j=1,n4
ccc        cd=cd+weights(j)*pol(j,i)*f(j)
        cd=cd+weights(j)*pol(j,i)*dconjg(f(j) )
 2200 continue
        proj(i)=cd
 2400 continue
cccc        call prin2('and proj=*',proj,m*2)
c 
c       adjust the list of maximal projections
c 
        do 2600 i=1,m
        dd=cdabs(proj(i))
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
        implicit complex *16 (a-h,o-z)
        save
        dimension a(n,n),x(1),y(1)
c 
        do 1400 i=1,n
        cd=0
        do 1200 j=1,n
        cd=cd+a(i,j)*dconjg(x(j))
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
  
        subroutine dpsquain(z,n,pols)
        implicit real *8 (a-h,o-z)
        save
        complex *16 pols(1),w(50),z
c 
c       this subroutine evaluates at the point z \in C
c       a set of polynomials in 1/z orthogonal on the square
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
c  pols - the array of values at z of polynomials in 1/z orthogonal
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
     1      -.8333333333333333D-01,0.2321104712717669D-01,
     2      -.3171227321546325D+00,-.1692488072818050D-01,
     3      0.6054750335329874D-01,-.2980048662351099D+00,
     4      0.1390122616582705D-01,-.3142782523099127D-01,
     5      0.5097755986850888D-01,-.2964226611066093D+00,
     6      -.1203570095256918D-01,0.2003804733136844D-01,
     7      -.2555250750668471D-01,0.4999564181077593D-01,
     8      -.2959474359210943D+00,0.1073802709549440D-01/
         data c0o1002/
     1      -.1411735537417572D-01,0.1610387223786470D-01,
     2      -.2486234650201367D-01,0.4966236125652845D-01,
     3      -.2957450401400203D+00,-.9768546548299155D-02,
     4      0.1055142941565671D-01,-.1136071336627268D-01,
     5      0.1559188006069414D-01,-.2460724024688252D-01,
     6      0.4950923520867122D-01,-.2956409244071075D+00,
     7      0.9008968842174105D-02,-.8198604793520655D-02,
     8      0.8574359717817647D-02,-.1096822927207895D-01/
         data c0o1003/
     1      0.1538969477961810D-01,-.2448319322504288D-01,
     2      0.4942629643689572D-01,-.2955805361673906D+00,
     3      -.8393242400624370D-02,0.6546312649002130D-02,
     4      -.6769139203113910D-02,0.8266992331329212D-02,
     5      -.1080432418551467D-01,0.1528682937811504D-01,
     6      -.2441323393393217D-01,0.4937636425293589D-01,
     7      -.2955424651053623D+00,0.7881185110620715D-02,
     8      -.5332071382614342D-02,0.5518362240932785D-02/
         data c0o1004/
     1      -.6525042376092375D-02,0.8131970503343357D-02,
     2      -.1071768099383926D-01,0.1522686398713020D-01,
     3      -.2436982331834214D-01,0.4934398620740282D-01,
     4      -.2955169445400395D+00,-.7446756426571456D-02,
     5      0.4408351480630354D-02,-.4608242310373483D-02,
     6      0.5322839451643372D-02,-.6412499541549217D-02,
     7      0.8058132502842102D-02,-.1066572525452352D-01,
     8      0.1518870166228712D-01,-.2434100665781017D-01/
         data c0o1005/
     1      0.4932179936155346D-01,-.2954990123788689D+00,
     2      0.7072238096552761D-02,-.3686250458894646D-02,
     3      0.3920810033328524D-02,-.4450966204936268D-02,
     4      0.5228218862527620D-02,-.6348999577225784D-02,
     5      0.8012736292498534D-02,-.1063193063496857D-01,
     6      0.1516285818546414D-01,-.2432088788815121D-01,
     7      0.4930593238923589D-01,-.2954859338990640D+00,
     8      -.6745098512278794D-02,0.3109209991483360D-02/
         data c0o1006/
     1      -.3386082469060416D-02,0.3794269280414200D-02,
     2      -.4370915674256759D-02,0.5173216787910365D-02,
     3      -.6309062352526430D-02,0.7982631438565259D-02,
     4      -.1060864292888798D-01,0.1514452145351753D-01,
     5      -.2430627916720875D-01,0.4929419185174059D-01,
     6      -.2954761027258240D+00,0.6456190385604883D-02,
     7      -.2639650211403476D-02,0.2960134041873377D-02,
     8      -.3284659040252782D-02,0.3726255508792356D-02/
         data c0o1007/
     1      -.4323003780285395D-02,0.5137883007234319D-02,
     2      -.6282108808389139D-02,0.7961561901328549D-02,
     3      -.1059188278830406D-01,0.1513102807189608D-01,
     4      -.2429533218279640D-01,0.4928525961824514D-01,
     5      -.2954685255400216D+00,-.6198659441894385D-02,
     6      0.2251695084916074D-02,-.2614117426142188D-02,
     7      0.2879527463572484D-02,-.3226724452287449D-02,
     8      0.3684334262763412D-02,-.4291595863712893D-02/
         data c0o1008/
     1      0.5113645955042246D-02,-.6262977203134631D-02,
     2      0.7946203515964730D-02,-.1057940255949068D-01,
     3      0.1512080317558339D-01,-.2428691450976260D-01,
     4      0.4927830489387208D-01,-.2954625615690532D+00,
     5      0.5967254656731943D-02,-.1926988574655482D-02,
     6      0.2328370269767563D-02,-.2550970809697635D-02,
     7      0.2830130670065879D-02,-.3189919821275780D-02,
     8      0.3656307153043135D-02,-.4269719757700861D-02/
         data c0o1009/
     1      0.5096217512638847D-02,-.6248868621944383D-02,
     2      0.7934644238654388D-02,-.1056985016557925D-01,
     3      0.1511286562498581D-01,-.2428030056029516D-01,
     4      0.4927278298097161D-01,-.2954577825397404D+00,
     5      -.5757876296711898D-02,0.1652175845294329D-02,
     6      -.2089067687970086D-02,0.2280019111586856D-02,
     7      -.2508878436138946D-02,0.2797737317393327D-02,
     8      -.3164829127082085D-02,0.3636498118582949D-02/
         data c0o1010/
     1      -.4253796539024307D-02,0.5083226569418355D-02,
     2      -.6238145880196252D-02,0.7925715986345575D-02,
     3      -.1056237075733734D-01,0.1510657766895044D-01,
     4      -.2427500783175005D-01,0.4926832470492101D-01,
     5      -.2954538934408503D+00,0.5567270548503863D-02,
     6      -.1417327262383874D-02,0.1886229387420031D-02,
     7      -.2053366253663412D-02,0.2244230387147144D-02,
     8      -.2480319643889410D-02,0.2775216415226884D-02/
         data c0o1011/
     1      -.3146841901243100D-02,0.3621912229076557D-02,
     2      -.4241808324356694D-02,0.5073263873343501D-02,
     3      -.6229794494509932D-02,0.7918670362459270D-02,
     4      -.1055640152889344D-01,0.1510150988637321D-01,
     5      -.2427070525426188D-01,0.4926467266240184D-01,
     6      -.2954506857093024D+00,-.5392817656501930D-02,
     7      0.1214919224924564D-02,-.1712485274618230D-02,
     8      0.1861427734327572D-02,-.2023059168403478D-02/
         data c0o1012/
     1      0.2219029132494011D-02,-.2460062880467284D-02,
     2      0.2758844494167697D-02,-.3133450075940255D-02,
     3      0.3610826747804919D-02,-.4232536979797767D-02,
     4      0.5065444407738910D-02,-.6223156486835944D-02,
     5      0.7913008915166076D-02,-.1055155926224615D-01,
     6      0.1509736438656763D-01,-.2426715948871378D-01,
     7      0.4926164293362835D-01,-.2954480085270747D+00,
     8      0.5232381503712398D-02,-.1039156000653977D-02/
         data c0o1013/
     1      0.1562284486797415D-02,-.1697140260189422D-02,
     2      0.1835919788840448D-02,-.2000816835309608D-02,
     3      0.2200779284553424D-02,-.2445130987538432D-02,
     4      0.2746523761745730D-02,-.3123180128228815D-02,
     5      0.3602185378745036D-02,-.4225207578674862D-02,
     6      0.5059187545178306D-02,-.6217788935366583D-02,
     7      0.7908388779258447D-02,-.1054757548323066D-01,
     8      0.1509392917504855D-01,-.2426420221245418D-01/
         data c0o1014/
     1      0.4925910126600948D-01,-.2954457506467983D+00,
     2      -.5084200625905260D-02,0.8855070588196525D-03,
     3      -.1431373445166449D-02,0.1555194148547660D-02,
     4      -.1675859325469281D-02,0.1816300226256927D-02,
     5      -.1984355601973177D-02,0.2187137544998141D-02,
     6      -.2433775190211503D-02,0.2736993085048578D-02,
     7      -.3115114022100996D-02,0.3595307566503698D-02,
     8      -.4219306118430490D-02,0.5054098518240076D-02/
         data c0o1015/
     1      -.6213384303396300D-02,0.7904567579132477D-02,
     2      -.1054425745013085D-01,0.1509105000517415D-01,
     3      -.2426170951120018D-01,0.4925694780026845D-01,
     4      -.2954438285549324D+00,0.4946807811495880D-02,
     5      -.7503843930299636D-03,0.1316442464215657D-02,
     6      -.1431527801345140D-02,0.1537656478669292D-02,
     7      -.1658576724466835D-02,0.1801440920240938D-02,
     8      -.1971874988967194D-02,0.2176654535890589D-02/
         data c0o1016/
     1      -.2424916863312057D-02,0.2729453487705922D-02,
     2      -.3108652698579634D-02,0.3589737067738248D-02,
     3      -.4214479802727675D-02,0.5049900812906553D-02,
     4      -.6209723373727691D-02,0.7901369906990632D-02,
     5      -.1054146382139317D-01,0.1508861244721912D-01,
     6      -.2425958851265842D-01,0.4925510700886204D-01,
     7      -.2954421785681296D+00,-.4818969828076344D-02,
     8      0.6309130274229664D-03,-.1214880995630637D-02/
         data c0o1017/
     1      0.1322986073414602D-02,-.1417321201481542D-02,
     2      0.1522466025965504D-02,-.1645158583026225D-02,
     3      0.1790009659041945D-02,-.1962184365111599D-02,
     4      0.2168409624214681D-02,-.2417860293892094D-02,
     5      0.2723376765991172D-02,-.3103390293361827D-02,
     6      0.3585157853336501D-02,-.4210479454851153D-02,
     7      0.5046395786164518D-02,-.6206646314258466D-02,
     8      0.7898666202354133D-02,-.1053908899956903D-01/
         data c0o1018/
     1      0.1508653014091326D-01,-.2425776849174029D-01,
     2      0.4925352090877864D-01,-.2954407514375660D+00,
     3      0.4699641579689604D-02,-.5247649112738940D-03,
     4      0.1124604658931169D-02,-.1227084081702318D-02,
     5      0.1311756845764754D-02,-.1404011783604066D-02,
     6      0.1510349732105019D-02,-.1634679284981728D-02,
     7      0.1781040918240528D-02,-.1954500777753483D-02,
     8      0.2161797243137926D-02,-.2412139264072693D-02/
         data c0o1019/
     1      0.2718401248106489D-02,-.3099043144843413D-02,
     2      0.3581344927446416D-02,-.4207124747906330D-02,
     3      0.5043437595124614D-02,-.6204034276677819D-02,
     4      0.7896359037969937D-02,-.1053705278082166D-01,
     5      0.1508473690656917D-01,-.2425619481521419D-01,
     6      0.4925214438518015D-01,-.2954395085915317D+00/
         real *8
     1  c0n1001(16),c0n1002( 9)
         real *8
     1  c0n2001(2)
         equivalence
     1  (c0n1001(16),c0n2001(1)),(c0n1002(1),c0n2001(2))
         data c0n1001/
     1      0.3535533905932738D+00,0.1462295969012131D+01,
     2      0.1358321068517427D+01,0.1355324858103474D+01,
     3      0.1354611674883064D+01,0.1354341934899910D+01,
     4      0.1354212926094367D+01,0.1354141704730985D+01,
     5      0.1354098377760547D+01,0.1354070106781928D+01,
     6      0.1354050656572635D+01,0.1354036709104230D+01,
     7      0.1354026369287696D+01,0.1354018491892527D+01,
     8      0.1354012352175765D+01,0.1354007473540769D+01/
         data c0n1002/
     1      0.1354003532229119D+01,0.1354000302101839D+01,
     2      0.1353997621350898D+01,0.1353995371730920D+01,
     3      0.1353993465208318D+01,0.1353991835149483D+01,
     4      0.1353990430367188D+01,0.1353989211013169D+01,
     5      0.1353988145690660D+01/
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
     1      -.1686032953945969D+00,0.3644890869155390D-01,
     2      -.3069399635834749D+00,-.2536857458699186D-01,
     3      0.5486628894849224D-01,-.2972658932390377D+00,
     4      0.2027204867500363D-01,-.2755855972741154D-01,
     5      0.5046994891377218D-01,-.2961793659534218D+00,
     6      -.1722043452785077D-01,0.1715626065809785D-01,
     7      -.2515795753378425D-01,0.4981372284798052D-01,
     8      -.2958384743346339D+00,0.1514448314894803D-01/
         data c1o1002/
     1      -.1185583235938954D-01,0.1578045252976612D-01,
     2      -.2471367732568563D-01,0.4957653937087295D-01,
     3      -.2956877680301834D+00,-.1362058478060512D-01,
     4      0.8714430780104228D-02,-.1108661481144036D-01,
     5      0.1546602030799700D-01,-.2453469213598762D-01,
     6      0.4946251635922937D-01,-.2956076879692315D+00,
     7      0.1244371794922645D-01,-.6669510040078426D-02,
     8      0.8336588576444109D-02,-.1085917824685651D-01/
         data c1o1003/
     1      0.1532683998454156D-01,-.2444273158896126D-01,
     2      0.4939849821708264D-01,-.2955599133021836D+00,
     3      -.1150125657500960D-01,0.5250109068952205D-02,
     4      -.6559282349647408D-02,0.8170883592418781D-02,
     5      -.1074895780370609D-01,0.1525116249738734D-01,
     6      -.2438872228601540D-01,0.4935879509623894D-01,
     7      -.2955290450724991D+00,0.1072563440959832D-01,
     8      -.4217720271704974D-02,0.5330644567973242D-02/
         data c1o1004/
     1      -.6439218165914596D-02,0.8082580635074362D-02,
     2      -.1068586243929995D-01,0.1520497046183707D-01,
     3      -.2435411972196129D-01,0.4933239444960570D-01,
     4      -.2955079024637379D+00,-.1007359276791660D-01,
     5      0.3439596706415881D-02,-.4438521673346053D-02,
     6      0.5245388396502651D-02,-.6367993502650411D-02,
     7      0.8029478484384649D-02,-.1064600119493339D-01,
     8      0.1517453317029614D-01,-.2433053874007921D-01/
         data c1o1005/
     1      0.4931390795563813D-01,-.2954927619110167D+00,
     2      0.9516012009025837D-02,-.2836442077187482D-02,
     3      0.3766011599428153D-02,-.4380465579169313D-02,
     4      0.5187778580318525D-02,-.6322993755688173D-02,
     5      0.7994842208063137D-02,-.1061906936918304D-01,
     6      0.1515334110494132D-01,-.2431370675975952D-01,
     7      0.4930043586333677D-01,-.2954815321641502D+00,
     8      -.9032493107475163D-02,0.2358207623413234D-02/
         data c1o1006/
     1      -.3243859492737625D-02,0.3729630621258539D-02,
     2      -.4333911979375591D-02,0.5149458474742329D-02,
     3      -.6292731449796165D-02,0.7970897398543895D-02,
     4      -.1059995496608063D-01,0.1513795603343810D-01,
     5      -.2430125003432666D-01,0.4929030085049587D-01,
     6      -.2954729632246696D+00,0.8608274256182155D-02,
     7      -.1971878468280989D-02,0.2828653976814595D-02,
     8      -.3225030637969073D-02,0.3692193884361977D-02/
         data c1o1007/
     1      -.4301175579050360D-02,0.5122901316358483D-02,
     2      -.6271355102305975D-02,0.7953602801443757D-02,
     3      -.1058586553600126D-01,0.1512641288575152D-01,
     4      -.2429175978701055D-01,0.4928247600633856D-01,
     5      -.2954662696030133D+00,-.8232380739501458D-02,
     6      0.1654872713300064D-02,-.2491921455642109D-02,
     7      0.2824230151762315D-02,-.3195209101387687D-02,
     8      0.3664180493418209D-02,-.4277789213018969D-02/
         data c1o1008/
     1      0.5103751056196641D-02,-.6255662084034337D-02,
     2      0.7940676381947208D-02,-.1057516270280846D-01,
     3      0.1511751837092381D-01,-.2428435510364280D-01,
     4      0.4927630308291990D-01,-.2954609370029451D+00,
     5      0.7896465210478102D-02,-.1391269657460205D-02,
     6      0.2214277874408273D-02,-.2499454217553422D-02,
     7      0.2800839873646698D-02,-.3171231554256944D-02,
     8      0.3643532005064793D-02,-.4260582174929857D-02/
         data c1o1009/
     1      0.5089474011410388D-02,-.6243780609406287D-02,
     2      0.7930744995307579D-02,-.1056683031209539D-01,
     3      0.1511051202734990D-01,-.2427846124809215D-01,
     4      0.4927134383075403D-01,-.2954566169226601D+00,
     5      -.7594053566246608D-02,0.1169552708905529D-02,
     6      -.1982109924585854D-02,0.2231830993324129D-02,
     7      -.2481547168906815D-02,0.2780341692460682D-02,
     8      -.3152965905911519D-02,0.3628032402227362D-02/
         data c1o1010/
     1      -.4247562723761679D-02,0.5078532903964665D-02,
     2      -.6234555579259919D-02,0.7922939715239723D-02,
     3      -.1056020938288318D-01,0.1510488972207164D-01,
     4      -.2427368978564787D-01,0.4926729701590849D-01,
     5      -.2954530662204521D+00,0.7320038993671384D-02,
     6      -.9812050390514335D-03,0.1785601079071288D-02,
     7      -.2008130454964011D-02,0.2218637758233957D-02,
     8      -.2464072071750169D-02,0.2764164541149424D-02/
         data c1o1011/
     1      -.3138975565369580D-02,0.3616134885755079D-02,
     2      -.4237469767456340D-02,0.5069953884911118D-02,
     3      -.6227241539051747D-02,0.7916687685064111D-02,
     4      -.1055485658823168D-01,0.1510030606438123D-01,
     5      -.2426977025771792D-01,0.4926394987395848D-01,
     6      -.2954501109567646D+00,-.7070333235056838D-02,
     7      0.8198061832856048D-03,-.1617510011478566D-02,
     8      0.1818828025587980D-02,-.1999019160008300D-02/
         data c1o1012/
     1      0.2203807491118747D-02,-.2449737013185115D-02,
     2      0.2751515614981609D-02,-.3128083283434307D-02,
     3      0.3606808922334271D-02,-.4229481630385011D-02,
     4      0.5063095909864289D-02,-.6221339158007034D-02,
     5      0.7911598160906687D-02,-.1055046435706571D-01,
     6      0.1509651772079534D-01,-.2426650948581802D-01,
     7      0.4926114855708456D-01,-.2954476239586428D+00,
     8      0.6841620499325663D-02,-.6804338741172385D-03/
         data c1o1013/
     1      0.1472388520445810D-02,-.1656908252945109D-02,
     2      0.1813274303551775D-02,-.1986517078345391D-02,
     3      0.2191106440638279D-02,-.2438286313285895D-02,
     4      0.2741527711306290D-02,-.3119452881363744D-02,
     5      0.3599361752901277D-02,-.4223046247790703D-02,
     6      0.5057522788179306D-02,-.6216503312446597D-02,
     7      0.7907396841699803D-02,-.1054681363203865D-01,
     8      0.1509334909475743D-01,-.2426376640715753D-01/
         data c1o1014/
     1      0.4925877951770417D-01,-.2954455104527967D+00,
     2      -.6631180175527079D-02,0.5592578509069089D-03,
     3      -.1346065890902612D-02,0.1517100114902967D-02,
     4      -.1654472883604440D-02,0.1802833010769898D-02,
     5      -.1975272915398549D-02,0.2180730997002949D-02,
     6      -.2429115211532814D-02,0.2733529912605748D-02,
     7      -.3112501707100064D-02,0.3593317693523994D-02,
     8      -.4217781968065764D-02,0.5052929123915428D-02/
         data c1o1015/
     1      -.6212488993854461D-02,0.7903886365877994D-02,
     2      -.1054374482521319D-01,0.1509067079909336D-01,
     3      -.2426143600378130D-01,0.4925675742774575D-01,
     4      -.2954436986082345D+00,0.6436756747295543D-02,
     5      -.4532578202315166D-03,0.1235300063689163D-02,
     6      -.1395373596625591D-02,0.1517412116658257D-02,
     7      -.1645864787153542D-02,0.1792893882775380D-02,
     8      -.1965866463990310D-02,0.2172300294897091D-02/
         data c1o1016/
     1      -.2421694404459405D-02,0.2727034304440330D-02,
     2      -.3106820084620664D-02,0.3588342471047819D-02,
     3      -.4213418134963420D-02,0.5049095746306449D-02,
     4      -.6209118203968120D-02,0.7900921664945609D-02,
     5      -.1054113935501632D-01,0.1508838573360818D-01,
     6      -.2425943874640529D-01,0.4925501713204606D-01,
     7      -.2954421331573772D+00,-.6256462764249859D-02,
     8      0.3600237042047616D-03,-.1137536232872269D-02/
         data c1o1017/
     1      0.1288599608440260D-02,-.1398117258857856D-02,
     2      0.1510442117358160D-02,-.1637099647290026D-02,
     3      0.1784364065354505D-02,-.1958109202735903D-02,
     4      0.2165407225455547D-02,-.2415618037855544D-02,
     5      0.2721688609719005D-02,-.3102115123708680D-02,
     6      0.3584195938204753D-02,-.4209758427372078D-02,
     7      0.5045861939862858D-02,-.6206258997809566D-02,
     8      0.7898394094918630D-02,-.1053890756449188D-01/
         data c1o1018/
     1      0.1508641983739207D-01,-.2425771348689181D-01,
     2      0.4925350816369841D-01,-.2954407709994823D+00,
     3      0.6088705396306393D-02,-.2776116388554089D-03,
     4      0.1050736380282367D-02,-.1194314954171309D-02,
     5      0.1293504419300568D-02,-.1392617023995724D-02,
     6      0.1502737191919925D-02,-.1629365760956744D-02,
     7      0.1777221304165051D-02,-.1951700120229665D-02,
     8      0.2159717446653271D-02,-.2410584033976061D-02/
         data c1o1019/
     1      0.2717236271399393D-02,-.3098173586352111D-02,
     2      0.3580702042801907D-02,-.4206657573943392D-02,
     3      0.5043107602757420D-02,-.6203811821059201D-02,
     4      0.7896220935428876D-02,-.1053698061199267D-01,
     5      0.1508471579982696D-01,-.2425621256654335D-01,
     6      0.4925219092800851D-01,-.2954395781111139D+00/
         real *8
     1  c1n1001(16),c1n1002( 9)
         real *8
     1  c1n2001(2)
         equivalence
     1  (c1n1001(16),c1n2001(1)),(c1n1002(1),c1n2001(2))
         data c1n1001/
     1      0.3989422804014327D+00,0.1399911396919482D+01,
     2      0.1357044197284382D+01,0.1354979483652475D+01,
     3      0.1354471088892951D+01,0.1354272143888605D+01,
     4      0.1354173930544813D+01,0.1354118145491300D+01,
     5      0.1354083342391252D+01,0.1354060121438406D+01,
     6      0.1354043827297265D+01,0.1354031936448195D+01,
     7      0.1354022982415447D+01,0.1354016064569223D+01,
     8      0.1354010604340545D+01,0.1354006215913895D+01/
         data c1n1002/
     1      0.1354002633779074D+01,0.1353999670231348D+01,
     2      0.1353997189445350D+01,0.1353995091117081D+01,
     3      0.1353993299835239D+01,0.1353991757994592D+01,
     4      0.1353990420961644D+01,0.1353989253708780D+01,
     5      0.1353988228427450D+01/
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
     1      -.2507099860678557D+00,0.5751204573730147D-01,
     2      -.2996609135709119D+00,-.3807516625055630D-01,
     3      0.5118056198044862D-01,-.2965453661412000D+00,
     4      0.2952134586724817D-01,-.2524564067571191D-01,
     5      0.5001929037807852D-01,-.2959822324990863D+00,
     6      -.2455332189290257D-01,0.1556455279108455D-01,
     7      -.2483763264747988D-01,0.4967941173157212D-01,
     8      -.2957662579231753D+00,0.2125043817005457D-01/
         data c2o1002/
     1      -.1070282877391490D-01,0.1553954845739520D-01,
     2      -.2461366540080605D-01,0.4952536147837912D-01,
     3      -.2956578880808723D+00,-.1886980155650694D-01,
     4      0.7853028690685599D-02,-.1089897700529129D-01,
     5      0.1538891125051344D-01,-.2449595768881694D-01,
     6      0.4944120148712926D-01,-.2955950561417438D+00,
     7      0.1705881164816667D-01,-.6013908415196293D-02,
     8      0.8186927160092505D-02,-.1079843330508711D-01/
         data c2o1003/
     1      0.1529697026498387D-01,-.2442681538831937D-01,
     2      0.4938980837596566D-01,-.2955551192628725D+00,
     3      -.1562692478391159D-01,0.4746192808403160D-02,
     4      -.6437899027359354D-02,0.8122350552581398D-02,
     5      -.1072572533402113D-01,0.1523928479610984D-01,
     6      -.2438267497340388D-01,0.4935597304713847D-01,
     7      -.2955280506742903D+00,0.1446143317289286D-01,
     8      -.3829341234213720D-02,0.5231018782052545D-02/
         data c2o1004/
     1      -.6400094322163121D-02,0.8064475791228047D-02,
     2      -.1067711848278235D-01,0.1520097696413317D-01,
     3      -.2435271541379502D-01,0.4933244672483019D-01,
     4      -.2955088078410105D+00,-.1349106262551334D-01,
     5      0.3141520162864311D-02,-.4356071317868516D-02,
     6      0.5213696717438206D-02,-.6353944487364508D-02,
     7      0.8023221093108141D-02,-.1064364565183256D-01,
     8      0.1517426484429917D-01,-.2433141216479711D-01/
         data c2o1005/
     1      0.4931539330355920D-01,-.2954946135303052D+00,
     2      0.1266835760146475D-01,-.2610425212296969D-02,
     3      0.3697402571999360D-02,-.4354762976483254D-02,
     4      0.5176998646360634D-02,-.6318743571295843D-02,
     5      0.7993813072251803D-02,-.1061972572448649D-01,
     6      0.1515488860516166D-01,-.2431570145294967D-01,
     7      0.4930261913314236D-01,-.2954838295202106D+00,
     8      -.1196039925332433D-01,0.2190682079762213D-02/
         data c2o1006/
     1      -.3186596257808005D-02,0.3708834313940778D-02,
     2      -.4325808223780634D-02,0.5146851580799499D-02,
     3      -.6292789545296244D-02,0.7972313674409692D-02,
     4      -.1060205901739851D-01,0.1514037675598788D-01,
     5      -.2430377014116640D-01,0.4929279226031471D-01,
     6      -.2954754360053319D+00,0.1134357649702631D-01,
     7      -.1852511882548673D-02,0.2780824419667467D-02,
     8      -.3208309459104175D-02,0.3686309699724170D-02/
         data c2o1007/
     1      -.4299930158917751D-02,0.5123859831195411D-02,
     2      -.6273400982832720D-02,0.7956168602656331D-02,
     3      -.1058864099571781D-01,0.1512921864794371D-01,
     4      -.2429448787360952D-01,0.4928506422786888D-01,
     5      -.2954687699242530D+00,-.1080047718096863D-01,
     6      0.1575639021381035D-02,-.2452031300926003D-02,
     7      0.2810933702215241D-02,-.3191186354994190D-02,
     8      0.3664074615357821D-02,-.4279500111920341D-02/
         data c2o1008/
     1      0.5106322345869581D-02,-.6258612861554897D-02,
     2      0.7943747485397285D-02,-.1057820708330806D-01,
     3      0.1512044906678255D-01,-.2428711930197538D-01,
     4      0.4927887291694861D-01,-.2954633843868280D+00,
     5      0.1031795569093743D-01,-.1345827513867717D-02,
     6      0.2181139525644104D-02,-.2489065118230955D-02,
     7      0.2798393228908749D-02,-.3172087893851754D-02,
     8      0.3645876482850937D-02,-.4263594730225180D-02/
         data c2o1009/
     1      0.5092747253735802D-02,-.6247098680931143D-02,
     2      0.7933988225279520D-02,-.1056993005212137D-01,
     3      0.1511342773511429D-01,-.2428116878144125D-01,
     4      0.4927383349085331D-01,-.2954589689501441D+00,
     5      -.9885886092373213D-02,0.1152817749260404D-02,
     6      -.1954765660275435D-02,0.2223931782448818D-02,
     7      -.2480446170294243D-02,0.2782016757608765D-02,
     8      -.3155847531391139D-02,0.3631417636032133D-02/
         data c2o1010/
     1      -.4251107205182939D-02,0.5082057803617458D-02,
     2      -.6237964561812969D-02,0.7926179762685448D-02,
     3      -.1056325094514929D-01,0.1510771615884201D-01,
     4      -.2427629206520618D-01,0.4926967504263203D-01,
     5      -.2954553023161666D+00,0.9496331582643124D-02,
     6      -.9890481287775888D-03,0.1763268656661116D-02,
     7      -.2002380074497215D-02,0.2218694038283336D-02,
     8      -.2466448439261775D-02,0.2767504245321497D-02/
         data c2o1011/
     1      -.3142677082852009D-02,0.3619908345341000D-02,
     2      -.4241168264889311D-02,0.5073501126946119D-02,
     3      -.6230597874070282D-02,0.7919832995987265D-02,
     4      -.1055778118699340D-01,0.1510300570300102D-01,
     5      -.2427224383710978D-01,0.4926620224338754D-01,
     6      -.2954522231721965D+00,-.9142975843083872D-02,
     7      0.8488358865701008D-03,-.1599543302619070D-02,
     8      0.1814944762593795D-02,-.2000077117157377D-02/
         data c2o1012/
     1      0.2206788108101767D-02,-.2453469368075141D-02,
     2      0.2755486671715484D-02,-.3132050622920601D-02,
     3      0.3610653370767734D-02,-.4233144237507713D-02,
     4      0.5066548528038093D-02,-.6224569742309787D-02,
     5      0.7914603012071253D-02,-.1055324384637431D-01,
     6      0.1509907392379502D-01,-.2426884534550098D-01,
     7      0.4926327129271726D-01,-.2954496116154067D+00,
     8      0.8820724265192741D-02,-.7278366243227127D-03/
         data c2o1013/
     1      0.1458248479604827D-02,-.1654657132118195D-02,
     2      0.1815204187033068D-02,-.1990021009695135D-02,
     3      0.2195176857108301D-02,-.2442487875655327D-02,
     4      0.2745659593965765D-02,-.3123420157900897D-02,
     5      0.3603120583934275D-02,-.4226578456247986D-02,
     6      0.5060823273398949D-02,-.6219573489421419D-02,
     7      0.7910241069766761D-02,-.1054943720403632D-01,
     8      0.1509575714366736D-01,-.2426596377952601D-01/
         data c2o1014/
     1      0.4926077437100830D-01,-.2954473769772197D+00,
     2      -.8525418346085266D-02,0.6226785197230611D-03,
     3      -.1335298782762875D-02,0.1516283426313357D-02,
     4      -.1657165582054761D-02,0.1806792218649297D-02,
     5      -.1979635505536598D-02,0.2185130247542638D-02,
     6      -.2433386967911171D-02,0.2737600576487624D-02,
     7      -.3116340686808101D-02,0.3596915510613032D-02,
     8      -.4221139536400109D-02,0.5056052244092089D-02/
         data c2o1015/
     1      -.6215385414508011D-02,0.7906564159763985D-02,
     2      -.1054621146928583D-01,0.1509293270926518D-01,
     3      -.2426349875087626D-01,0.4925862931664440D-01,
     4      -.2954454496401840D+00,0.8253627275927879D-02,
     5      -.5307091210976528D-03,0.1227520888251203D-02,
     6      -.1395823673758093D-02,0.1520775174671978D-02,
     7      -.1650221659397952D-02,0.1797509809779504D-02,
     8      -.1970435638464780D-02,0.2176691063003937D-02/
         data c2o1016/
     1      -.2425852028249334D-02,0.2730939870985462D-02,
     2      -.3110471749728928D-02,0.3591746392201285D-02,
     3      -.4216583899116125D-02,0.5052033990970879D-02,
     4      -.6211839347076182D-02,0.7903435207956115D-02,
     5      -.1054345349536736D-01,0.1509050720956336D-01,
     6      -.2426137320972345D-01,0.4925677260925209D-01,
     7      -.2954437754148813D+00,-.8002493293136952D-02,
     8      0.4498165333746199D-03,-.1132415936749826D-02/
         data c2o1017/
     1      0.1290173248804587D-02,-.1402071793060161D-02,
     2      0.1515147547667047D-02,-.1641935837677707D-02,
     3      0.1789079539543751D-02,-.1962601247848567D-02,
     4      0.2169637861858622D-02,-.2419578699447664D-02,
     5      0.2725384188161206D-02,-.3105556375337024D-02,
     6      0.3587395741373345D-02,-.4212729899163641D-02,
     7      0.5048617479601898D-02,-.6208809798753603D-02,
     8      0.7900749891187455D-02,-.1054107650198126D-01/
         data c2o1018/
     1      0.1508840846830025D-01,-.2425952720413116D-01,
     2      0.4925515449945022D-01,-.2954423115988174D+00,
     3      0.7769615200833096D-02,-.3783008574698660D-03,
     4      0.1047991868848467D-02,-.1196889057191978D-02,
     5      0.1297982703075110D-02,-.1397628904286643D-02,
     6      0.1507765320329131D-02,-.1634207327196627D-02,
     7      0.1781799476919616D-02,-.1955991874989551D-02,
     8      0.2163723425089302D-02,-.2414315099193114D-02/
         data c2o1019/
     1      0.2720707234525947D-02,-.3101400204718825D-02,
     2      0.3583699577386227D-02,-.4209440110236617D-02,
     3      0.5045687747465592D-02,-.6206200578574242D-02,
     4      0.7898427663517831D-02,-.1053901301418846D-01,
     5      0.1508657998772279D-01,-.2425791350751947D-01,
     6      0.4925373556462397D-01,-.2954410241647116D+00/
         real *8
     1  c2n1001(16),c2n1002( 9)
         real *8
     1  c2n2001(2)
         equivalence
     1  (c2n1001(16),c2n2001(1)),(c2n1002(1),c2n2001(2))
         data c2n1001/
     1      0.4410127717245515D+00,0.1367837547551257D+01,
     2      0.1355674455237705D+01,0.1354672412064680D+01,
     3      0.1354367929122946D+01,0.1354230944806902D+01,
     4      0.1354156399244618D+01,0.1354110949473895D+01,
     5      0.1354081042932447D+01,0.1354060249708541D+01,
     6      0.1354045175338536D+01,0.1354033881143561D+01,
     7      0.1354025190672348D+01,0.1354018354681744D+01,
     8      0.1354012876786210D+01,0.1354008417113803D+01/
         data c2n1002/
     1      0.1354004736369771D+01,0.1354001661983684D+01,
     2      0.1353999066883059D+01,0.1353996855770566D+01,
     3      0.1353994956012055D+01,0.1353993311442005D+01,
     4      0.1353991878063269D+01,0.1353990621005416D+01,
     5      0.1353989512336727D+01/
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
     1      -.3261003909500324D+00,0.8352595609186629D-01,
     2      -.2976140553587530D+00,-.5396437173900848D-01,
     3      0.5059223743396662D-01,-.2965200317183248D+00,
     4      0.4102426840173101D-01,-.2518473312837981D-01,
     5      0.5008547342409649D-01,-.2960608589955904D+00,
     6      -.3359778675251184D-01,0.1575573239373785D-01,
     7      -.2494618607083816D-01,0.4976610642324920D-01,
     8      -.2958387692965211D+00,0.2871864486859545D-01/
         data c3o1002/
     1      -.1103108254580128D-01,0.1567009005942865D-01,
     2      -.2470533617355211D-01,0.4959545063915554D-01,
     3      -.2957154937350009D+00,-.2523960738928506D-01,
     4      0.8261424065263946D-02,-.1104143052185246D-01,
     5      0.1548269651322046D-01,-.2456550630906142D-01,
     6      0.4949546065808298D-01,-.2956401586345128D+00,
     7      0.2261845284844122D-01,-.6471092996029730D-02,
     8      0.8335806371811183D-02,-.1089270734061659D-01/
         data c3o1003/
     1      0.1536562658214427D-01,-.2447984545707503D-01,
     2      0.4943212562253658D-01,-.2955908052270571D+00,
     3      -.2056379995594575D-01,0.5233589232012515D-02,
     4      -.6590017052217559D-02,0.8216190441264706D-02,
     5      -.1079322133622252D-01,0.1529114137535255D-01,
     6      -.2442394262562008D-01,0.4938953956900417D-01,
     7      -.2955567337825731D+00,0.1890432144499172D-01,
     8      -.4335343591714978D-02,0.5384439316768438D-02/
         data c3o1004/
     1      -.6492971678652126D-02,0.8130667764042224D-02,
     2      -.1072778152790040D-01,0.1524128344744474D-01,
     3      -.2438551889571397D-01,0.4935954789479286D-01,
     4      -.2955322317482292D+00,-.1753231747862947D-01,
     5      0.3658586680483417D-02,-.4509581134025771D-02,
     6      0.5305311265892905D-02,-.6418771467682735D-02,
     7      0.8072691959628682D-02,-.1068300100078796D-01,
     8      0.1520636377641604D-01,-.2435799316526276D-01/
         data c3o1005/
     1      0.4933763542824757D-01,-.2955140253079792D+00,
     2      0.1637649325766649D-01,-.3133507480325529D-02,
     3      0.3850226626996926D-02,-.4444950265594071D-02,
     4      0.5240449826703261D-02,-.6367045292970513D-02,
     5      0.8032232935255860D-02,-.1065111871658435D-01,
     6      0.1518097677730825D-01,-.2433760485294182D-01,
     7      0.4932114145571434D-01,-.2955001293650763D+00,
     8      -.1538766618255663D-01,0.2716324638895907D-02/
         data c3o1006/
     1      -.3338234999631406D-02,0.3797512667801596D-02,
     2      -.4387902274107802D-02,0.5194020902992052D-02,
     3      -.6330300092017745D-02,0.8003005299192891D-02,
     4      -.1062763831889848D-01,0.1516194640974380D-01,
     5      -.2432208382829131D-01,0.4930841639116439D-01,
     6      -.2954892836101279D+00,0.1453074521650830D-01,
     7      -.2378310893449566D-02,0.2930956606167731D-02,
     8      -.3295448451326307D-02,0.3747082120874572D-02/
         data c3o1007/
     1      -.4346011197450874D-02,0.5160494203685190D-02,
     2      -.6303405394894958D-02,0.7981231120552919D-02,
     3      -.1060985032269197D-01,0.1514731441230160D-01,
     4      -.2430999522152875D-01,0.4929839366351495D-01,
     5      -.2954806567424545D+00,-.1377998115792572D-01,
     6      0.2099903821922563D-02,-.2600453777165853D-02,
     7      0.2896534881077925D-02,-.3250681773615195D-02,
     8      0.3709114709177325D-02,-.4315295011301202D-02/
         data c3o1008/
     1      0.5135660279615708D-02,-.6283160847562691D-02,
     2      0.7964579417487599D-02,-.1059605215161025D-01,
     3      0.1513582188837885D-01,-.2430039639056248D-01,
     4      0.4929035892435855D-01,-.2954736824756891D+00,
     5      0.1311602598588188D-01,-.1867361989488700D-02,
     6      0.2327728747604544D-02,-.2573149942560684D-02,
     7      0.2856660691717329D-02,-.3216135054338663D-02,
     8      0.3680869830282427D-02,-.4292290630896189D-02/
         data c3o1009/
     1      0.5116790107994031D-02,-.6267546633808934D-02,
     2      0.7951559808207049D-02,-.1058513261331168D-01,
     3      0.1512663012326839D-01,-.2429264744902810D-01,
     4      0.4928381917330515D-01,-.2954679642703016D+00,
     5      -.1252404283873569D-01,0.1670775138792037D-02,
     6      -.2099452476176770D-02,0.2306533927697106D-02,
     7      -.2537536269482536D-02,0.2825118156249778D-02,
     8      -.3190077010483824D-02,0.3659497662019320D-02/
         data c3o1010/
     1      -.4274658736713424D-02,0.5102122705553714D-02,
     2      -.6255250217440931D-02,0.7941186163305163D-02,
     3      -.1057634208275985D-01,0.1511916297689810D-01,
     4      -.2428630166815458D-01,0.4927842528101169D-01,
     5      -.2954632177484429D+00,0.1199245184778190D-01,
     6      -.1502832291287719D-02,0.1906021418492833D-02,
     7      -.2083540455185285D-02,0.2274657166818335D-02,
     8      -.2508649548109640D-02,0.2801006442612780D-02/
         data c3o1011/
     1      -.3170167826583948D-02,0.3642984847697192D-02,
     2      -.4260855941880254D-02,0.5090495947883753D-02,
     3      -.6245392097615030D-02,0.7932785664074381D-02,
     4      -.1056916012547291D-01,0.1511301408774077D-01,
     5      -.2428103950419024D-01,0.4927392439585250D-01,
     6      -.2954592346734892D+00,-.1151207463480774D-01,
     7      0.1358032778879345D-02,-.1740356349392023D-02,
     8      0.1894708280012837D-02,-.2054962428621057D-02/
         data c3o1012/
     1      0.2248132245438784D-02,-.2486279297079259D-02,
     2      0.2782414373772229D-02,-.3154669607987497D-02,
     3      0.3629972692010960D-02,-.4249848410160883D-02,
     4      0.5081121856323817D-02,-.6237365960313279D-02,
     5      0.7925886699439909D-02,-.1056321628917179D-01,
     6      0.1510789017883231D-01,-.2427662748909515D-01,
     7      0.4927012975433722D-01,-.2954558597362472D+00,
     8      0.1107553653156095D-01,-.1232165896472825D-02/
         data c3o1013/
     1      0.1597134163220542D-02,-.1733070538038810D-02,
     2      0.1869058998402534D-02,-.2030549135189477D-02,
     3      0.2227327722717258D-02,-.2468877960278867D-02,
     4      0.2767838980011558D-02,-.3142381787783696D-02,
     5      0.3619537639393519D-02,-.4240927234936597D-02,
     6      0.5073451967511189D-02,-.6230743071413024D-02,
     7      0.7920150774253497D-02,-.1055824105721663D-01,
     8      0.1510357527815443D-01,-.2427289182813193D-01/
         data c3o1014/
     1      0.4926690098518623D-01,-.2954529751968895D+00,
     2      -.1067683993236896D-01,0.1121958790597556D-02,
     3      -.1472281876331083D-02,0.1593393919740330D-02,
     4      -.1710035056269036D-02,0.1846542880483604D-02,
     5      -.2011158616739702D-02,0.2211007063702466D-02,
     6      -.2455144568142569D-02,0.2756216148215473D-02,
     7      -.3132476302041938D-02,0.3611039531185278D-02,
     8      -.4233594520808811D-02,0.5067095304438440D-02/
         data c3o1015/
     1      -.6225213370526227D-02,0.7915329822117850D-02,
     2      -.1055403453835756D-01,0.1509990748340932D-01,
     3      -.2426970098556634D-01,0.4926413093658715D-01,
     4      -.2954504904851610D+00,0.1031105400395211D-01,
     5      -.1024832729967176D-02,0.1362634624222219D-02,
     6      -.1471677985247063D-02,0.1572702179928531D-02,
     7      -.1689231032222486D-02,0.1828434582296999D-02,
     8      -.1995822331171789D-02,0.2198044270538998D-02/
         data c3o1016/
     1      -.2444133600379143D-02,0.2746801055928699D-02,
     2      -.3124373123586268D-02,0.3604024892672216D-02,
     3      -.4227492650419616D-02,0.5061767131845760D-02,
     4      -.6220548112465306D-02,0.7911238667531964D-02,
     5      -.1055044593006502D-01,0.1509676348638207D-01,
     6      -.2426695392436687D-01,0.4926173666950364D-01,
     7      -.2954483350104357D+00,-.9974085192162123D-02,
     8      0.9387309637843959D-03,-.1265699176778649D-02/
         data c3o1017/
     1      0.1364817052266551D-02,-.1453096873038355D-02,
     2      0.1553449522966763D-02,-.1672289849586409D-02,
     3      0.1813998002955284D-02,-.1983566848739564D-02,
     4      0.2187597563234146D-02,-.2435173246354836D-02,
     5      0.2739066568997177D-02,-.3117657988630159D-02,
     6      0.3598165689669819D-02,-.4222359487609882D-02,
     7      0.5057256096286536D-02,-.6216575531156955D-02,
     8      0.7907736798501103D-02,-.1054735971931830D-01/
         data c3o1018/
     1      0.1509404803060806D-01,-.2426457199511911D-01,
     2      0.4925965317652312D-01,-.2954464530973994D+00,
     3      0.9662504900888378D-02,-.8619949286497932D-03,
     4      0.1179487069870940D-02,-.1270366603237089D-02,
     5      0.1348144115714400D-02,-.1435255209035489D-02,
     6      0.1537574397383738D-02,-.1658678208369412D-02,
     7      0.1802393546600370D-02,-.1973641675348057D-02,
     8      0.2179059535278297D-02,-.2427783144784890D-02/
         data c3o1019/
     1      0.2732633199034366D-02,-.3112029197892704D-02,
     2      0.3593220003794398D-02,-.4217999343522481D-02,
     3      0.5053402623319795D-02,-.6213164595636372D-02,
     4      0.7904716020408970D-02,-.1054468620819277D-01,
     5      0.1509168658332706D-01,-.2426249323500023D-01,
     6      0.4925782892278573D-01,-.2954448003583880D+00/
         real *8
     1  c3n1001(16),c3n1002( 9)
         real *8
     1  c3n2001(2)
         equivalence
     1  (c3n1001(16),c3n2001(1)),(c3n1002(1),c3n2001(2))
         data c3n1001/
     1      0.4791222828142150D+00,0.1357903391753444D+01,
     2      0.1355394835859144D+01,0.1354725192081621D+01,
     3      0.1354438495968396D+01,0.1354289665295741D+01,
     4      0.1354202631118196D+01,0.1354147384037920D+01,
     5      0.1354110141952416D+01,0.1354083855421413D+01,
     6      0.1354064614868786D+01,0.1354050110661817D+01,
     7      0.1354038907277060D+01,0.1354030074448166D+01,
     8      0.1354022987893768D+01,0.1354017216046611D+01/
         data c3n1002/
     1      0.1354012452842758D+01,0.1354008476347066D+01,
     2      0.1354005122458464D+01,0.1354002267726165D+01,
     3      0.1353999817839964D+01,0.1353997699750291D+01,
     4      0.1353995856166473D+01,0.1353994241646965D+01,
     5      0.1353992819776084D+01/
        complex *16 zinv
c 
c       construct the orthogonal polynomials of orders 0, 4, 8, ...
c 
        done=1
        zinv=done/z
        n4=n/4+1
        index=0
        call dorteva(zinv,n4,c0o1001,c0n1001,index,w)
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
        call dorteva(zinv,n4,c1o1001,c1n1001,index,w)
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
        call dorteva(zinv,n4,c2o1001,c2n1001,index,w)
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
        call dorteva(zinv,n4,c3o1001,c3n1001,index,w)
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
  
