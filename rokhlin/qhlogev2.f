        implicit real *8 (a-h,o-z)
        dimension hilcoefs(100),tlege(100)
        dimension coefslog(100),whts(100),
     1      coefsqua(100),w(100 000)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER n'
        READ *,n
        CALL PRInf('n=*',n,1 )
C 
        PRINT *, 'ENTER x'
        READ *,x
        CALL PRIn2('x=*',x,1 )
c 
c        initialize the array w to be used in the construction,
c        for a user-specified real x and integer n, the
c        coefficients hilcoefs, coefsqua, coefslog of the linear forms
c        converting the values of a function f at n Gaussian nodes on
c        the interval [-1,1] into the values at the point x of the
c        Hilbert, quadrupole, and logarithm transforms of f, respectively.
c 
        call qhlogini(n,tlege,whts,w,ltot,lsave)
c 
c        construct the coefficients hilcoefs, coefsqua, coefslog of the
c        linear forms converting the values of a function f at n Gaussian
c        nodes on the interval [-1,1] into the values at the point x of
c        the Hilbert, quadrupole, and logarithm transforms of f,
c        respectively.
c 
cccc        call qhlogeva(x,hilcoefs,coefsqua,coefslog,w)
  
        a=-1.4
        b=2
  
cccc        a=-1
cccc        b=1
c 
        uu=(b-a)/2
        vv=b-uu
c 
        call qhlogev2(a,b,x,w,hilcoefs,coefslog,coefsqua)
c 
        call prin2('after qhlogeva, hilcoefs=*',hilcoefs,n)
        call prin2('after qhlogeva, coefslog=*',coefslog,n)
        call prin2('after qhlogeva, coefsqua=*',coefsqua,n)
c 
c        evaluate the Hilbert transform by direct integration and
c        compare it to the one obtained vya hilcoefs
c 
        d=0
        do 1200 i=1,n
c 
        dd=funtest(uu*tlege(i)+vv)
c 
cccc        dd=funtest(tlege(i))
c 
        d=d+dd*hilcoefs(i)
 1200 continue
c 
        call prin2('and hilbert transform via quahillo*',d,1)
c 
c 
        call numinhil(x,rint,a,b)
c 
        call prin2('and integral numerically is*',rint,1)
c 
        call prin2('and difference with numerical hilbert transform*',
     1      d-rint,1)
c 
        call prin2('and the ratio*',d/rint,1)
c 
c        evaluate the logarithm transform by direct integration and
c        compare it to the one obtained vya coefslog
c 
        d=0
        ddd=0
        do 2200 i=1,n
c 
        dd=funtest(uu*tlege(i)+vv)
c 
        d=d+dd*coefslog(i)
 2200 continue
c 
        call prin2('and log transform via coefslog*',d,1)
c 
        call numinlog(x,rint,a,b)
  
        call prin2('and integral numerically is*',rint,1)
c 
        call prin2('and difference with numerical log transform*',
     1      d-rint,1)
c 
        call prin2('and the ratio*',d/rint,1)
c 
c        evaluate the quadrupole transform by direct integration and
c        compare it to the one obtained via coefsqua
c 
        d=0
        ddd=0
        do 3200 i=1,n
c 
        dd=funtest(uu*tlege(i)+vv)
c 
        d=d+dd*coefsqua(i)
 3200 continue
c 
        call prin2('and quadrupole transform via coefsqua*',d,1)
c 
        call numinqua(x,rint,a,b)
  
        call prin2('and integral numerically is*',rint,1)
c 
        call prin2('difference with numerical quadrupole transform*',
     1      d-rint,1)
c 
        call prin2('and the ratio*',d/rint,1)
  
  
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine numinqua(y,rint,a,b)
        implicit real *8 (a-h,o-z)
        save
        external funqua
c 
c       evaluate Hilbert transform by adaptive integration
c 
        eps=1.0d-10
        m=10
c 
        call adapgaus(ier,a,b,funqua,y,par2,m,eps,
     1      rint,maxrec,numint)
        return
        end
c 
c 
c 
c 
c 
        function funqua(x,y,par2)
        implicit real *8 (a-h,o-z)
c 
        save
        delta=1.0d-12
        h=1.0d-6
        d=funtest(x)
        funqua=d* ( (x-y)**2-h**2) / ((x-y)**2+h**2)**2
        return
        end
c 
c 
c 
c 
c 
        subroutine numinlog(y,rint,a,b)
        implicit real *8 (a-h,o-z)
        save
        external funlog
c 
c       evaluate Hilbert transform by adaptive integration
c 
        eps=1.0d-10
        m=10
c 
        call adapgaus(ier,a,b,funlog,y,par2,m,eps,
     1      rint,maxrec,numint)
        return
        end
c 
c 
c 
c 
c 
        function funlog(x,y,par2)
        implicit real *8 (a-h,o-z)
c 
        save
  
  
        delta=1.0d-12
        d=funtest(x)
        funlog=d*dlog((y-x)**2+delta) /2
        return
        end
c 
c 
c 
c 
c 
        subroutine numinhil(y,rint,a,b)
        implicit real *8 (a-h,o-z)
        save
        external funhilb
c 
c       evaluate Hilbert transform by adaptive integration
c 
        eps=1.0d-10
        m=10
c 
        call adapgaus(ier,a,b,funhilb,y,par2,m,eps,
     1      rint,maxrec,numint)
        return
        end
c 
c 
c 
c 
c 
        function funhilb(x,y,par2)
        implicit real *8 (a-h,o-z)
c 
        save
        delta=1.0d-12
        d=funtest(x)
        funhilb=d*(y-x)/((y-x)**2+delta)
        return
        end
c 
c 
c 
c 
c 
        function funtest(x)
        implicit real *8 (a-h,o-z)
c 
        save
        done=1
        pi=datan(done)*4
  
        f=sin(3*x)+cos(2*x)
cccc        f=sin(3*x)
cccc        f=sin(x*pi)
  
        f=1+f
c 
        funtest=f
        return
        end
c 
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c          This is the end of the debugging code and the beginning of the
c          actual code for the Hilbert, logarithm, and quadrupole transforms
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        subroutine qhlogev2(a,b,x,w,hilcoefs,coefslog,coefsqua)
        implicit real *8 (a-h,o-z)
        save
        dimension coefsqua(1),w(1),coefslog(1),hilcoefs(1)
c 
c        For the user-supplied real x and integer n, this subroutine
c        constructe the coefficients hilcoefs, coefsqua, coefslog
c        of the linear forms converting the values of a function f
c        at n Gaussian nodes on the interval [a,b] into the values
c        at the point x of the Hilbert, quadrupole, and logarithm
c        transforms of f, respectively. This subroutine uses as its
c        input tha array w, constructed by the initialization
c        subroutine qhlogini (see). This subroutine has no function
c        as a stand-alone device. Furthermore, this is only a scaling
c        and post-processing routine; all actual work is done by the
c        subroutine qhlogeva (see).
c 
c                 Input parameters:
c 
c  a - the left end of the interval on which the function's
c        Hilbert, log, and quadrupule transforms will be evaluated
c  b - the right end of the interval on which the function's
c        Hilbert, log, and quadrupule transforms will be evaluated
c  x - the points where the quadrupole, logarithm, and Hilbert
c        transforms will be evaluated
c  w - array created by the subroutine qhlogini (see). Note that
c        the first lsave elements of the array w must be unchanged
c        between the call to this subroutine and the preceding call
c        to the subroutine qhlogini.
c 
c                 Output paramaters:
c 
c  hilcoefs - coefficients of the linear form converting the values
c        of a function f at n Gaussian nodes into the value at the
c        point x of the Hilbert transform of f.
c  coefsqua - coefficients of the linear form converting the values
c        of a function f at n Gaussian nodes into the value at the
c        point x of the quadrupole transform of f.
c  coefslog - coefficients of the linear form converting the values
c        of a function f at n Gaussian nodes into the value at the
c        point x of the logarithm transform of f.
c 
c 
c        . . . construct the linear mapping converting the
c              interval [a,b] into the interval [-1,1], and
c              apply that mapping to x
c 
        u=2/(b-a)
        v=1-u*b
c 
        y=u*x+v
c 
cccc        call prin2('u=*',u,1)
cccc        call prin2('v=*',v,1)
cccc        call prin2('y=*',y,1)
c 
c       construct the arrays hilcoefs, coefslog, coefsqua for the
c       transformed configuration
c 
        call qhlogeva(y,hilcoefs,coefsqua,coefslog,w)
c 
c        scale the coefficients coefsqua to account for the length
c        of the interval of integration that is not equal to 2
c 
        n=w(4)
        iwhts=w(7)
c 
        d=log(u)
c 
        do 2200 i=1,n
c 
        coefsqua(i)=coefsqua(i)*u
c 
        coefslog(i)=coefslog(i)-w(iwhts+i-1)*d
c 
        coefslog(i)=coefslog(i)/u
 2200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine qhlogeva(x,hilcoefs,coefsqua,coefslog,w)
        implicit real *8 (a-h,o-z)
        save
        dimension coefsqua(1),w(1),coefslog(1),hilcoefs(1)
c 
c        For the user-supplied real x and integer n, this subroutine
c        constructe the coefficients hilcoefs, coefsqua, coefslog
c        of the linear forms converting the values of a function f
c        at n Gaussian nodes on the interval [-1,1] into the values
c        at the point x of the Hilbert, quadrupole, and logarithm
c        transforms of f, respectively. This subroutine uses as its
c        input tha array w, constructed by the initialization
c        subroutine qhlogini (see). This subroutine has no function
c        as a stand-alone device. Furthermore, this is only a memory
c        management routine; all actual work is done by the
c        subroutine quahillo (see).
c 
c                 Input parameters:
c 
c  x - the points where the quadrupole, logarithm, and Hilbert
c        transforms will be evaluated
c  w - array created by the subroutine qhlogini (see). Note that
c        the first lsave elements of the array w must be unchanged
c        between the call to this subroutine and the preceding call
c        to the subroutine qhlogini.
c 
c                 Output paramaters:
c 
c  hilcoefs - coefficients of the linear form converting the values
c        of a function f at n Gaussian nodes into the value at the
c        point x of the Hilbert transform of f.
c  coefsqua - coefficients of the linear form converting the values
c        of a function f at n Gaussian nodes into the value at the
c        point x of the quadrupole transform of f.
c  coefslog - coefficients of the linear form converting the values
c        of a function f at n Gaussian nodes into the value at the
c        point x of the logarithm transform of f.
c 
c        . . . remember the locations in the memory of the arrays
c              to be used
c 
        ibinte=w(1)
        ibdiff=w(2)
        iendinte=w(3)
        n=w(4)
        iqfuns=w(5)
        iu=w(6)
        iwhts=w(7)
c 
cccc        call prinf('in qhlogeva, ibinte=*',ibinte,1)
cccc        call prinf('in qhlogeva, ibdiff=*',ibdiff,1)
cccc        call prinf('in qhlogeva, iendinte=*',iendinte,1)
cccc        call prinf('in qhlogeva, n=*',n,1)
cccc        call prinf('in qhlogeva, iqfuns=*',iqfuns,1)
c 
c        construct the linear forms converting values of
c        a function at Gaussian nodes into values of its
c        Hilbert, log, and quadrupole transforms at the
c        point x
c 
c 
        ifinit=0
        ifeval=1
c 
        call quahillo(n,coefsqua,w(iu),x,tlege,w(ibinte),w(iwhts),
     1      w(ibdiff),coefslog,hilcoefs,ainte,adiff,
     2      v,w(iendinte),w(iqfuns),w,ifinit,ifeval)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine qhlogini(n,tlege,whts,w,ltot,lsave)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),tlege(1),whts(1)
c 
c        This is the initialization subroutine for the subroutine
c        qhlogeva (see). Tha latter, for user-supplied real x and
c        integer n, constructs the coefficients hilcoefs, coefsqua,
c        coefslog of the linear forms converting the values of a
c        function f at n Gaussian nodes on the interval [-1,1] into
c        the values at the point x of the Hilbert, quadrupole, and
c        logarithm transforms of f, respectively. This subroutine's
c        output is the array w, to be used by the subroutine qhlogeva
c        (see). This subroutine has no function as a stand-alone
c        device. Furthermore, this is only a memory management routine;
c        all actual work is done by the subroutine quahillo (see).
  
c 
c                 Input parameters:
c 
c  n - the number of Gaussian nodes on the interval [-1,1] where
c        the function is tabulated whose quadrupole, log, and Hilbert
c        transforms will be evaluated
c 
c                 Output paramaters:
c  tlege - the n Gaussian nodes on the interval [-1,1]
c  whts - the n Gaussian weights on the interval [-1,1]
c  w - array containing various types of information to be used by
c        the subroutine qhlogeva (see). note that the first lsave
c        elements of the array w must not be changed between the call
c        to this subroutine and the subsequent calls to qhlogeva.
c  ltot - the number of elements of array w used by this subroutine
c  lused - the number of elements ov w that will be used by the
c        subroutine qhlogeva (see). these elements should not be
c        changed between the call to this subroutine and the
c        subsequent calls to qhlogeva.
c 
c        ... allocate memory for the construction of the various
c            arrays to be used for the evaluation of the Hilbert,
c             log, and quadrupole transforms
c 
        ibinte=10
        lbinte=n**2+2
c 
        ibdiff=ibinte+lbinte
        lbdiff=n**2+2
c 
        iendinte=ibdiff+lbdiff
        lendinte=n+2
c 
        iu=iendinte+lendinte
        lu=n**2+2
c 
        iwhts=iu+lu
        lwhts=n+1
c 
        iv=iwhts+lwhts
        lv=n**2+2
c 
        iainte=iv+lv
        lainte=n**2+2
c 
        iadiff=iainte+lainte
        ladiff=n**2+2
c 
        iw=iadiff+ladiff
        lw=n**2*4+10
c 
        ltot=iw+lw
c 
c        construct all the necessary arrays
c 
        ifinit=1
        ifeval=0
c 
        call quahillo(n,coefsqua,w(iu),y,tlege,w(ibinte),whts,
     1      w(ibdiff),coefslog,hilcoefs,w(iainte),w(iadiff),
     2      w(iv),w(iendinte),qfuns,w(iw),ifinit,ifeval)
c 
c       store in the first elements of array w the locations
c       in w of various arrays to be used later by the
c       subroutine ghlogeva
c 
        do 2200 i=1,n
        w(iwhts+i-1)=whts(i)
 2200 continue
c 
        iqfuns=iwhts+lwhts
        lqfuns=n+2
c 
        lsave=iqfuns+lqfuns
c 
        w(1)=ibinte+0.1
        w(2)=ibdiff+0.1
        w(3)=iendinte+0.1
        w(4)=n+0.1
        w(5)=iqfuns+0.1
        w(6)=iu+0.1
        w(7)=iwhts+0.1
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine quahillo(n,coefsqua,u,y,tlege,binte,whts,
     1      bdiff,coefslog,hilcoefs,ainte,adiff,v,endinter,
     2      qfuns,w,ifinit,ifeval)
        implicit real *8 (a-h,o-z)
        save
        dimension coefsqua(1),u(n,n),v(1),qfuns(1),
     1      tlege(1),whts(1),ainte(1),adiff(1),
     2      endinter(1),w(1),binte(n,n),bdiff(n,n),
     3      coefslog(1),hilcoefs(1)
c 
c        construct the Gaussian nodes and weights on
c        the interval [-1,1] and the Legendre expansion matrices
c 
        if(ifinit .eq. 0) goto 1200
        itype=2
        call legeexps(itype,n,tlege,u,v,whts)
c 
c        construct the matrices of spectral integration and
c        differentiation on the n Legendre nodes
c 
        itype=3
        call legeinmt(n,ainte,adiff,tlege,whts,endinter,
     1      itype,w)
c 
cccc        call prin2('in quahillo, ainte=*',ainte,n**2)
cccc        call prin2('in quahillo, adiff=*',adiff,n**2)
cccc        call prin2('while u=*',u,n**2)
c 
c        multiply the matrix converting the values of a function
c        at the Gaussian nodes into its Legendre coefficients
c        from the left by the differentiating matrix
c 
        call matmat(u,ainte,n,binte)
        call matmat(u,adiff,n,bdiff)
c 
 1200 continue
c 
c       evaluate the functions Q_k at the point y
c 
        if(ifeval .eq. 0) return
c 
        call qneval(y,n,qfuns(2) )
c 
c       construct the coefficients of the linear function
c       converting the values of a function at Gaussian nodes
c       into the value of its quadrupole transform at the
c       user-specified point y
c 
        call vecmat(qfuns,bdiff,n,coefsqua)
c 
c        correct  for the non-zero values of the function at the ends
c 
        dleft=(y-1)/(y-1)**2
        dright=(y+1)/(y+1)**2
c 
        do 1800 i=1,n
c 
        coefsqua(i)=-coefsqua(i)*2-endinter(n-i+1)*dright
        coefsqua(i)=coefsqua(i)+endinter(i)*dleft
c 
 1800 continue
c 
c       construct the coefficients of the linear function
c       converting the values of a function at Gaussian nodes
c       into the value of its log transform at the user-specified
c       point y
c 
        call vecmat(qfuns,binte,n,coefslog)
c 
c        correct  for the non-zero integral of the function
c 
        done=1
        dd=dlog((y-1)**2) /4
c 
        do 2800 i=1,n
c 
        coefslog(i)=coefslog(i)+whts(i)*dd
c 
        coefslog(i)=coefslog(i)*2
 2800 continue
c 
c       construct the coefficients of the linear function
c       converting the values of a function at Gaussian nodes
c       into the value of its Hilbert transform at the
c       user-specified point y
c 
        call vecmat(qfuns,u,n,hilcoefs)
c 
        do 3800 i=1,n
        hilcoefs(i)=hilcoefs(i)*2
 3800 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine matmat(a,b,n,c)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),b(n,n),c(n,n)
c 
        do 2000 i=1,n
        do 1400 j=1,n
        d=0
        do 1200 k=1,n
        d=d+a(i,k)*b(k,j)
 1200 continue
c 
        c(i,j)=d
 1400 continue
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine qneval(x,n,qfuns)
        implicit real *8 (a-h,o-z)
        save
        dimension qfuns(1)
c 
c        for user-specified real x and integer n, this subroutine
c        evaluates the Legendre functions Q_0, Q_1, . . . ,Q_n at
c        the point x anywhere on the real line (except |x|=1).
c 
c                        Input parameters:
c 
c  x - the point where the functions Q_i will be evaluated
c  n - maximum order of Q to be evaluated
c 
c                        Output parameters:
c 
c  qfuns - the values of the Legendre functions at the
c      point x (n+1 of them things)
c 
c           IMPORTANT NOTE:
c 
c      ARRAY QFUNS MUST HAVE ELEMENT NUMBER 0 !!!!!!
c 
c        . . . if x is inside the interval [-1,1], or sufficiently
c              close to 1 or -1, evaluate the Legenedre functions
c              via the simple recursion
c 
        if(dabs(x) .lt. 1) goto 1100
c 
        delta=dabs(x)-1
        b=1+delta+sqrt(2*delta+delta**2)
c 
        if( (n+1)*log(b) .gt. 2.3) goto 2000
c 
 1100 continue
c 
c        construct Q_0(x),Q_1(x)
c 
        done=1
        i0=0
        d=log( abs((done+x)/(done-x)) )
c 
        qfuns(i0)=d/2
        qfuns(1)=x/2*d-1
c 
c       recurse up
c 
        do 1200 i=1,n-1
c 
        qfuns(i+1)=( (2*i+1)*x*qfuns(i)-i*qfuns(i-1) ) /(i+1)
 1200 continue
c 
        return
c 
 2000 continue
c 
c       the point x is outside the interval [-1,1].
c       determine how far we will have to start to get it by
c       the combination of recursing up and scaling
c 
        eps=1.0d-20
c 
        delta=dabs(x)-1
        b=1+delta+sqrt(2*delta+delta**2)
        nn=-log(eps)/log(b)+1
c 
        call prinf('nn as obtained in qninte is*',nn,1)
c 
c       use recursion down coupled with scaling to get the
c       Legendre functions Q_i
c 
        done=1
        i0=0
        d=log( abs((done+x)/(done-x)) )
c 
        q0=d/2
c 
        do 2050 i=i0,n
        qfuns(i)=0
 2050 continue
c 
c        if nn is greater than n+1, do the preliminary recursion
c        to get the starting values for qfuns(n+1),qfuns(n)
c 
        fip1=0
        fi  =1
c 
        if(nn .le. n+1) goto 2150
c 
        do 2100 i=nn,n+1,-1
c 
         fim1=( (2*i+1)*x*fi-(i+1)*fip1 ) / i
c 
        fip1=fi
        fi=fim1
 2100 continue
c 
        nn=n
 2150 continue
c 
        qfuns(nn+1)=fip1
        qfuns(nn)=fi
c 
c       recurse down starting with i=n
c 
        do 2200 i=nn,1,-1
c 
         qfuns(i-1)=( (2*i+1)*x*qfuns(i)-(i+1)*qfuns(i+1) ) / i
 2200 continue
c 
c        scale the values of the recursed (unscaled) Q_i to
c        obtain the correct one
c 
        rat=q0/qfuns(i0)
c 
        call prin2('and rat=*',rat,1)
c 
        do 2400 i=i0,n+1
c 
        qfuns(i)=rat*qfuns(i)
 2400 continue
 3000 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine vecmat(x,a,n,y)
        implicit real *8 (a-h,o-z)
        save
        dimension x(n),a(n,n),y(n)
c 
        do 1400 i=1,n
        d=0
        do 1200 j=1,n
        d=d+x(j)*a(j,i)
 1200 continue
        y(i)=d
 1400 continue
c 
        return
        end
