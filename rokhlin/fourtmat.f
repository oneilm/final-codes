        implicit real *8 (a-h,o-z)
        complex *16 w(100 000),a(10 000),
     2      f1(1000),f2(1000),f3(1000),f4(1000)
        dimension tfour(1000),tgauss(1000),whts(1000)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER nfour'
        READ *,nfour
        CALL PRINF('nfour=*',nfour,1 )
C 
        PRINT *, 'ENTER ngauss'
        READ *,ngauss
        CALL PRINF('ngauss=*',ngauss,1 )
c 
c        construct the matrix interpolating trigonometric
c        polynomials from the equispaced to the Gaussian nodes
c 
         call fourtint(tfour,nfour,tgauss,whts,ngauss,a,w)
c 
            call prin2('after fourtint, a=*',a,ngauss*nfour*2)
c 
c         now, evaluate a test function at equispaced
c         nodes on the interval [0, 2*pi], interpolate it to the
c         gaussian nodes, and compare the result with the test
c         function evaluated at the gaussian nodes directly
c 
        do 2200 i=1,nfour
        call funfool(tfour(i),f1(i))
 2200 continue
c 
        call prin2('f1 as created*',f1,nfour*2)
cccc        call matvec(a,f1,f2,nfour)
        call matvec3(a,ngauss,nfour,f1,f2)
cccc        call matvec3(a,nfour,ngauss,f1,f2)
  
  
        call prin2('f2 as interpolated*',f2,ngauss*2)
  
c 
        do 2400 i=1,ngauss
        call funfool(tgauss(i),f3(i))
         f4(i)=f3(i)-f2(i)
 2400 continue
c 
         call prin2('and f3 as computed directly is*',f3,ngauss*2)
         call prin2('and f3 -f2=*',f4,ngauss*2)
  
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine funfool(t,f)
        implicit real *8 (a-h,o-z)
        save
        complex *16 f
c 
        f=dsin(t*2)+dcos(t)**4
        return
        end
c 
c 
c 
c 
c 
        subroutine matvec3(a,n,m,x,y)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,m),x(m),y(n),cd
c 
c        apply the matrix a to the vector x obtaining y
c 
        do 1400 i=1,n
        cd=0
        do 1200 j=1,m
        cd=cd+a(i,j)*x(j)
 1200 continue
        y(i)=cd
 1400 continue
        return
        end
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c       this is the end of the debugging code and the beginning
c       of the actual expansion-interpolation routines
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine fourtint(tfour,nfour,tgauss,whts,ngauss,c,w)
        implicit real *8 (a-h,o-z)
        save
        complex *16 c(ngauss,nfour),ima,w(1)
        dimension tgauss(1),whts(1),tfour(1)
        data ima/(0.0d0,1.0d0)/
c 
c        this subroutine constructs the matrix interpolating a
c        trigonometric polynomial from an equispaced grid
c        onto a Guassian one
c 
c 
c                       input parameters:
c 
c  nfour - the number of equispaced nodes on the interval [0,2*pi]
c        from which the interpolation is to be performed
c  ngauss - the number of Gaussian nodes on the interval [0,2*pi]
c        onto which the interpolation is to be performed
c 
c                       output parameters:
c 
c  tfour - nfour equispaced nodes on the interval [0, 2*pi]
c  tgauss -  ngauss gaussian nodes on the interval [0, 2*pi]
c  whts -  ngauss gaussian weights corresponding to the nodes tgauss
c  c - the (ngauss * nfour) matrix (complex) interpolating from the
c        nfour equispaced nodes onto the ngauss gaussian nodes on the
c        interval [0, 2*pi]
c 
c                       work arrays:
c 
c  w - must be at least nfour**2 +ngauss*nfour+20 complex *16 elements
c 
c        . . . construct the matrix converting the values
c              of a function on the equispaced grid into
c              coefficients of its fourier series
c 
        iv=1
        lv=nfour**2+10
c 
        iu=iv+lv
c 
         call fourtmat(nfour,tfour,w(iu),w(iv) )
cccc         call prin2('in fourtint, w(iu)=*',w(iu),nfour**2*2)
cccc         call prin2('in fourtint, w(iv)=*',w(iv),nfour**2*2)
c 
c        construct the matrix converting the coefficients of a
c        fourier series of a function into its values on the
c        Gaussian grid
c 
        ia=iu
c 
        call fourexeva(ngauss,tgauss,whts,nfour,w(ia) )
cccc         call prin2('in fourtint after fourexeva, w(ia)=*',
cccc     1       w(ia),nfour*ngauss*2)
c 
c       calculate the product of these two matrices
c 
        call fourmatp(w(ia),w(iv),nfour,ngauss,nfour,c)
        return
        end
c 
c 
c 
c 
c 
        subroutine fourmatp(a,b,n,m,k,c)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(m,n),b(n,k),c(m,k),cd
c 
c        multiply matrix a by matrix b getting matrix c
c 
        do 1400 i=1,m
        do 1200 j=1,k
c 
        cd=0
        do 1100 l=1,n
        cd=cd+a(i,l)*b(l,j)
 1100 continue
        c(i,j)=cd
 1200 continue
 1400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine fourexeva(ngauss,tgauss,whts,nfour,a)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(ngauss,nfour),ima
        dimension tgauss(1),whts(1)
        data ima/(0.0d0,1.0d0)/
c 
c        this subroutine constructs the table of values of
c        complex exponentials at the gaussian nodes on the
c        interval [0,2*pi]
c 
c                       input parameters:
c 
c  ngauss - the number of gaussian nodes on the interval [0,2*pi]
c        at which the table of values of exponentials is to be
c        constructed
c  nfour - the number of complex exponentials for which the tables
c        of values are to be constructed.
c 
c                       output parameters:
c 
c  tgauss -  ngauss gaussian nodes on the interval [0, 2*pi]
c  whts -  ngauss gaussian weights corresponding to the nodes tgauss
c  a - the (ngauss*nfour) matrix containing the values of various
c        exponentials at the nodes tgauss
c 
c    EXPLANATION: for j .le. nfour/2, the j-th column of a
c                 contains values of the function cdexp(ima*(j-1)*x)
c                 at the gaussian nodes on the interval [1,2*pi].
c 
c                 for j .gt. nfour/2, the j-th column of a contains
c                 values of the function cdexp(ima*(nfour-j+1)*x)
c                 at the gaussian nodes on the interval [1,2*pi].
c 
c          . . . construct the Gaussian nodes and corresponding
c                weights on the interval [-1,1]
c 
        call gauswhts(ngauss,tgauss,whts)
c 
c        convert these nodes and weights into those on the
c        interval [0,2*pi]
c 
        done=1
        pi=datan(done)*4
        u=pi
        v=pi
        do 1200 i=1,ngauss
        tgauss(i)=u*tgauss(i)+v
        whts(i)=whts(i)*u
 1200 continue
c 
        do 1600 i=1,nfour
        do 1400 j=1,ngauss
c 
        a(j,i)=0
 1400 continue
 1600 continue
c 
c       construct the matrix
c 
        do 2000 i=1,ngauss
        do 1800 j=1,nfour/2
        a(i,j)=cdexp(ima*tgauss(i)*(j-1))
c 
          if(j .eq. nfour/2) goto 1800
        a(i,nfour-j+1)=cdexp(-ima*tgauss(i)*(j) )
 1800 continue
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine fourtmat(n,t,u,v)
        implicit real *8 (a-h,o-z)
        save
        complex *16 u(n,n),ima,v(n,n)
        dimension t(1)
        data ima/(0.0d0,1.0d0)/
c 
c        this subroutine constructs explicit matrices for the DFT and
c        its inverse.
c 
c                       input parameters:
c 
c  n - the number of nodes on the interval [0,2*pi] for which
c        the DFT is to be constructed
c 
c                       output parameters:
c 
c  t - the n equispaced nodes on the interval [0, 2*pi]
c  u - the matrix of values of the complex exponentials at equispaced
c        nodes on the interval [0, 2*pi]
c  v - the matrix converting the values of a function
c        at the nodes into the coefficients of its fourier series
c 
c          NOTE: obviously, u and v are inverses of each other
c 
c 
c        . . . construct the n-point discrtetization of the
c              interval [0,2*pi]
c 
        done=1
        pi=datan(done)*4
        h=2*pi/n
        do 1200 i=1,n
        t(i)=(i-1)*h
 1200 continue
c 
c       construct the matrix of values of various complex exponentials
c       at the nodes
c 
        do 1600 i=1,n
        do 1400 j=1,n
        u(j,i)=cdexp(ima*(t(i)*(j-1)) )
 1400 continue
 1600 continue
c 
c        construct the matrix converting the values of a function
c        at the nodes into the coefficients of its fourier series
c 
        d=n
        d=done/d
        do 2400 i=1,n
        do 2200 j=1,n
        v(j,i)=dconjg(u(i,j))*d
 2200 continue
 2400 continue
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine gauswhts(n,ts,whts)
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
        eps=1.0d-14
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
        call legpol(xk,n,pol,der)
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
        subroutine legpol(x,n,pol,der)
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
