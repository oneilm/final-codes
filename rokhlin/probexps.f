        implicit real *8 (a-h,o-z)
c 
        external fun1
c 
        dimension ts(10 000),xs(10 000),whts(10 000),
     1      whts2(10 000),fs(10 000),dxdt(10 000),
     2      u(500 000),v(500 000),coefs(10 000),
     3      xstest(10 000),ftest(10 000),ftest2(10 000),
     4      diffs(10 000),ders(10 000),ders2(10 000),
     5      valsp(10 000),valsm(10 000),diffs2(10 000)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER n'
        READ *,n
        CALL PRINf('n=*',n,1 )
c 
c       construct the sl0 nodes and weights
c 
        itype=2
        ndigits=10
  
        call probexps(ndigits,itype,n,xs,u,v,whts2)
  
        call prin2('after probexps, xs=*',xs,n)
        call prin2('after probexps, whts2=*',whts2,n)
c 
c       construct the values of the test function at the
c       sl0 nodes, and evaluate the integral
c 
        rint=0
        do 1200 i=1,n
c 
        fs(i)=fun1(xs(i),par1,par2)
c 
        rint=rint+whts2(i)*fs(i)
c 
 1200 continue
c 
        call prin2('rint as calculated*',rint,1)
cccc        call prin2('and rint-2=*',rint-2,1)
  
c 
c        evaluate the same integral via adapgaus
c 
        a=-1
        b=1
        m=20
        eps=1.0d-12
  
  
c 
        call adapgaus(ier,a,b,fun1,x,n,m,eps,
     1      rint2,maxrec,numint)
  
        call prinf('after adapgaus, ier=*',ier,1)
        call prin2('after adapgaus, eps=*',eps,1)
  
        call prin2('after adapgaus, rint2=*',rint2,1)
        call prinf('after adapgaus, numint=*',numint,1)
  
        call prin2('and rint2-rint=*',rint2-rint,1)
  
  
cccc        stop
  
c 
c       calculate the cosine expansion of the test function
c       in the resampled variable
c 
        call matvec(u,fs,coefs,n)
  
  
        call prin2('and coefficients of Legendre expansion are*',
     1      coefs,n)
        call prin2('while fs=*',fs,n)
  
  
  
cccc        stop
c 
c        now, construct a bunch of test points on the interval
c        [-1,1], and evaluate the function at these points
c 
        ntest=300
        two=2
        htest=two/ntest
        do 2200 i=1,ntest
c 
        xstest(i)=-1+(i-1)*htest+htest/2
  
        ftest(i)=fun1(xstest(i),par1,par2)
  
        call protFDER(ndigits,Xstest(i),VAL,der,coefs,N)
  
        ftest2(i)=val
  
        diffs(i)=ftest(i)-ftest2(i)
  
  
        ders(i)=der
 2200 continue
  
  
        call prin2('anf xstest=*',xstest,ntest)
        call prin2('and ftest=*',ftest,ntest)
        call prin2('and ftest2=*',ftest2,ntest)
        call prin2('and diffs=*',diffs,ntest)
  
        call prin2('and ders=*',ders,ntest)
  
c 
c       test the derivatives
c 
        h=0.0001
        do 2400 i=1,ntest
c 
        call protFDER(ndigits,Xstest(i)+h,VALsp(i),der,coefs,N)
        call protFDER(ndigits,Xstest(i)-h,VALsm(i),der,coefs,N)
c 
        ders2(i)=(valsp(i)-valsm(i))/h/2
  
        diffs2(i)=ders2(i)-ders(i)
 2400 continue
  
        call prin2('and ders2=*',ders2,ntest)
        call prin2('and diffs2=*',diffs2,ntest)
  
  
        stop
        end
  
  
c 
c 
c 
c 
c 
        subroutine matvec(a,x,y,n)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),x(n),y(n)
c 
        do 1400 i=1,n
        d=0
        do 1200 j=1,n
        d=d+a(i,j)*x(j)
 1200 continue
        y(i)=d
 1400 continue
        return
        end
c 
c 
c 
c 
c 
        function fun1(x,par1,par2)
        implicit real *8 (a-h,o-z)
c 
        save
        bb=1+1.0d-14
        bb=1
  
cccc        bb=0
cccc        fun1=(1-x**2)*log(1-x**2)
        fun1=sqrt(bb-x**2)*log(bb-x**2)
        fun1=log(bb-x**2) *(bb-x**2)**2
cccc        fun1=1/sqrt(sqrt((bb-x**2)))
c 
        return
  
  
  
  
cccc        fun1=(bb-x**2)*log(bb-x**2) *cos(120*x)
  
ccccc        fun1=(bb-x**2)*log(bb-x**2)
  
  
        nnn=6
        call legepol(x,nnn,pol,der)
  
        fun1=fun1*pol
  
        fun1=x**2
  
        return
        end
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c        This is the ens of the debugging code, and the beginning of the
c        "universal" diecretization code proper
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine probexps(ndigits,itype,n,xs,u,v,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),whts(1),u(1),v(1)
c 
c         this subroutine constructs the two-sided "universal"
c         quadrature nodes on the interval [-1,1], and the weights
c         for the corresponding n-point quadrature. It also constructs
c         the matrix v converting the coefficients of a "distorted
c         cosine" expansion into its values at the n "universal" nodes,
c         and its inverse u, converting the values of a function
c         at n "universal" nodes into the coefficients of the
c         corresponding "distorted cosine" series. No attempt has
c         been made to make this code efficient, but its speed is
c         normally sufficient, and it is mercifully short.
c 
c                 input parameters:
c 
c  ndigits - the accuracy (the number of digits) to which the
c        calculations will be conducted; must be between 1 and 16
c  itype - the type of the calculation to be performed
c     itype=1 means that only the nodes and the weights
c                  are to be constructed
c     itype=2 means that the nodes, the weights, and
c                  the matrices u, v are to be constructed
c  n - the number of "universal" nodes and weights to be generated
c 
c                 output parameters:
c 
c  ts - the n "universal" nodes on the interval [-1,1] - computed
c          independently of the value of itype.
c  u - the n*n matrix converting the  values at of a cosine series
c         of order n-1 at n equispaced nodes into the coefficients
c         of its cosine expansion - computed only if itype=2
c  v - the n*n matrix converting the coefficients of an n-term
c         cosine expansion into its values at n equispaced nodes
c         (note that v is the inverse of u)  - computed only if
c         itype=2
c  whts - the quadrature weights corresponding to the "universal"
c         nodes; computed independently of the value of itype.
c 
c 
c        . . . construct the n equispaced nodes on the interval [-1,1],
c              and the corresponding nodes on the interval [0,pi]
c              parametrizing space
c 
c 
c       construct the equispaced nodes and weights
c 
  
        call cosexps(itype,n,xs,whts,u,v)
c 
c       construct the images of the Gaussian nodes and the
c       corresponding weights after prolate resampling
c 
        do 1200 i=1,n
c 
        call proresf(ndigits,xs(i),d,dd)
c 
        xs(i)=d
        whts(i)=whts(i)*dd
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        SUBROUTINE protFDER(ndigits,X,VAL,der,PEXP,N)
        implicit real *8 (a-h,o-z)
        save
        dimension pexp(1)
c 
c        This subroutine evaluates the user-supplied "distorted
c        cosine" expansion (together with its derivative) at the
c        user-specified point x on the interval [-1,1]
c 
c 
c                input parameters:
c 
c  ndigits - the accuracy (the number of digits) to which the
c        calculations will be conducted; must be between 1 and 16
C  X = evaluation point
C     PEXP = expansion coefficients
C  N  = order of expansion
c         IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c 
c                output parameters:
c 
C  VAL = computed value
C  der = computed value of the derivative
C 
c 
c        . . . find the location of the point x in the
c              parametrizing space
c 
        call proresb(ndigits,x,t,dtdx)
c 
c       evaluate the function to be evaluated
c 
        call cosFDER(t,VAL,der,PEXP,N)
c 
        der=der*dtdx
c 
        return
        end
c 
c 
c 
c 
c 
        SUBROUTINE protFDEm(ndigits,X,VALs,ders,PEXPs,N,m)
        implicit real *8 (a-h,o-z)
        save
        dimension pexps(n,m),vals(1),ders(1)
c 
c        This subroutine evaluates n user-supplied "distorted
c        cosine" expansions (together with their derivatives) at
c        the user-specified point x on the interval [-1,1]
c 
c     NOTE: THIS SUBROUTINE IS A MULTI-EXPANSION VERSION OF THE
C     SUBROUTINE PROTFDER (SEE). IT EXISTS BECAUSE OF THE
c     unreasonable cpu time requirements of the subroutine
C     PRORESB (SEE) USED BY BOTH PROTFDER and this subroutine.
c     HOWEVER, PROTFDER CALLS PRORESB ONCE FOR EACH EXPANSION
C     IT EVALUATES, WHILE THIS SUBROUTINE CALLS PROTFDER ONCE
C     FOR M EXPANSIONS IT EVALUATES. ALSO, PLEASE NOTE THAT THE
C     NOTATION USED BY THIS SUBROUTINE IS SOMEWHAT INCONSISTENT
C     WITH THE NOTATION USED PBY PROTFDER: THE PARAMETER N MEANS
C     SLIGHTLY (BY 1) DIFFERENT THINGS.
c 
c                input parameters:
c 
c  ndigits - the accuracy (the number of digits) to which the
c        calculations will be conducted; must be between 1 and 16
C  X - evaluation point
C  PEXPs - coefficients of all m expansions
C  N - number of terms in each of the expansions
c         IMPORTANT NOTE: n is {\bf number of terms} in each of the
c         expansions, which is one greater than the order!
c  m - the number pf expansions to be evaluated
c 
c                output parameters:
c 
C  VALs - computed values
C  ders - computed value of the derivative
C 
c 
c        . . . find the location of the point x in the
c              parametrizing space
c 
        call proresb(ndigits,x,t,dtdx)
c 
c       evaluate the functions to be evaluated
c 
cccc        call prinf('in protfdem, n=*',n,1)
cccc        call prinf('in protfdem, m=*',m,1)
cccc        call prin2('in protfdem, x=*',x,1)
cccc        call prin2('in protfdem, t=*',t,1)
  
        do 1400 i=1,m
c 
        call cosFDER(t,VALs(i),ders(i),PEXPs(1,i),N-1)
c 
        ders(i)=ders(i)*dtdx
 1400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine pro1sexp(ndigits,itype,n,xs,u,v,
     1      whts,w)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),whts(1),u(1),v(1),w(1)
c 
c         this subroutine constructs the one-sided "universal"
c         quadrature nodes on the interval [0,1], and the weights
c         for the corresponding n-point quadrature. It also constructs
c         the matrix v converting the coefficients of a "distorted
c         Legendre" expansion into its values at the n "universal" nodes,
c         and its inverse u, converting the values of a function
c         at n "universal" nodes into the coefficients of the
c         corresponding "distorted Legendre" series. No attempt has
c         been made to make this code efficient, but its speed is
c         normally sufficient, and it is mercifully short.
c 
c                 input parameters:
c 
c  ndigits - the accuracy (the number of digits) to which the
c        calculations will be conducted; must be between 1 and 16
c  itype - the type of the calculation to be performed
c     itype=1 means that only the nodes and the weights
c                  are to be constructed
c     itype=2 means that the nodes, the weights, and
c                  the matrices u, v are to be constructed
c  n - the number of "universal" nodes and weights to be generated
c 
c                 output parameters:
c 
c  ts - the n "universal" nodes on the interval [-1,1] - computed
c          independently of the value of itype.
c  u - the n*n matrix converting the  values at of a cosine series
c         of order n-1 at n equispaced nodes into the coefficients
c         of its cosine expansion - computed only if itype=2
c  v - the n*n matrix converting the coefficients of an n-term
c         cosine expansion into its values at n equispaced nodes
c         (note that v is the inverse of u)  - computed only if
c         itype=2
c  whts - the quadrature weights corresponding to the "universal"
c         nodes; computed independently of the value of itype.
c 
c                 Work arrays:
c 
c  w - must be at least 8*n**2+4*n+100 real *8 locations long
c 
c 
c       . . . allocate memory for the construction of the discretization
c             and the corresponding expansion-evaluation matrices
c 
        ixsgauss=1
        lxsgauss=n*2+2
c 
        iwsgauss=ixsgauss+lxsgauss
        lwsgauss=n*2+2
c 
        iu=iwsgauss+lwsgauss
        lu=n*n*4+10
c 
        iv=iu+lu
        lv=n*n*4+10
c 
c       conatruct the discretization and the corresponding
c       expansion-evaluation matrices (if the latter have been requested)
c 
        if(itype .eq. 2)
     1      call pr01quad(ndigits,n,xs,whts,itype,
     2      w(ixsgauss),w(iwsgauss),w(iu),w(iv),u,v)
c 
        if(itype .eq. 1)
     1      call pr01quad(ndigits,n,xs,whts,itype,
     2      w(ixsgauss),w(iwsgauss),udummy,vdummy,udummy,vdummy)
c 
        return
        end
c 
c 
c 
c 
c 
        SUBROUTINE pro1tFDE(ndigits,X,VAL,der,PEXP,N,w)
        implicit real *8 (a-h,o-z)
        save
        dimension pexp(1),w(1000)
c 
c        This subroutine evaluates the user-supplied "distorted
c        Legendre" expansion (together with its derivative) at the
c        user-specified point x on the interval [0,1]
c 
c 
c                input parameters:
c 
c  ndigits - the accuracy (the number of digits) to which the
c        calculations will be conducted; must be between 1 and 16
C  X = evaluation point
C     PEXP = expansion coefficients
C  N  = order of expansion
c         IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c 
c                output parameters:
c 
C  VAL = computed value
C  der = computed value of the derivative
c 
c                work arrays:
c 
c  w - must be at least n*2+2 real *8 locations long
C 
C 
c        . . . find the location of the point x in the
c              parametrizing space
c 
        call proresb(ndigits,x,t,dtdx)
c 
        t=t-1
c 
c       evaluate the function to be evaluated
c 
        do 1400 i=1,n*2
        w(i)=0
 1400 continue
c 
        do 1500 i=1,n
c 
        w(2*i-1)=pexp(i)
 1500 continue
c 
        call legeFDER(t,VAL,der,w,N*2-2)
        der=der*dtdx
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine pr01quad(ndigits,n,xs,whts,itype,
     1      xsgauss,wsgauss,u,v,uout,vout)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),whts(1),xsgauss(1),wsgauss(1),
     1      uout(n,n),vout(n,n),u(2*n,2*n),v(2*n,2*n)
c 
c       . . . construct the Gaussian nodes and weights
c 
        call legeexps(itype,n*2,xsgauss,u,v,wsgauss)
c 
c       construct the images of the Gaussian nodes and the
c       corresponding weights after prolate resampling
c 
        do 1200 i=1,n
c 
        call proresf(ndigits,xsgauss(i)+1,d,dd)
c 
        xs(i)=d
        whts(i)=wsgauss(i)*dd
 1200 continue
c 
        if(itype .ne. 2) return
c 
c        extract from the matrices u,v the n \times n submatrices that
c        convert the values of the function at the first n nodes
c        into the EVEN-numbered coefficients of the Legendre expansion
c        (assuming theat the function is even)
c 
        do 1600 i=1,n
        do 1400 j=1,n
c 
        uout(j,i)=u(j*2-1,i)*2
c 
        vout(j,i)=v(j,i*2-1)
 1400 continue
 1600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine proresb(ndigits,x,t,dtdx)
        implicit real *8 (a-h,o-z)
c 
c       This subroutine implements the mapping [-1,1] \to [-1,1]
c       that is the inverse of the mapping implemented by the
c       function proresf (see)
c 
c 
c        . . . use bisection to find the crude approximation
c              to t to which the value of x corresponds
c 
        save
        tk=x
        numit=65
c 
        d=1-abs(x)
        if(d .gt. 0.001) numit=10
  
        t1=-1
        t2=1
        x1=-1
        x2=1
  
        do 1600 i=1,numit
c 
        t3=(t1+t2)/2
c 
        call proresf(ndigits,t3,x3,dxdt)
  
        if(x3 .gt. x) goto 1200
c 
        t1=t3
        x1=x3
        goto 1600
c 
 1200 continue
c 
        t2=t3
        x2=x3
c 
 1600 continue
c 
c       use Newton to tighten up the approximate value of t
c 
        tk=(t1+t2)/2
        numit=6
        do 2200 i=1,numit
c 
        call proresf(ndigits,tk,f,dfdt)
  
        f=f-x
        tkp1=tk-f/dfdt
        tk=tkp1
 2200 continue
c 
        t=tk
        dtdx=1/dfdt
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine probquad(n,xs,whts,w)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),whts(1),w(1)
c 
c       This subroutine constructs a "left-unversal" quadrature
c       on the interval [-1,0]. The claim to fame of this subroutine
c       is spectral convergence for functions with a singularity at
c       (or near) the left end. Please note that in order to provide
c       double precision accuracy, the obtained quadrature has to
c       have at least 14 (or so) nodes. Furthermore, for smooth
c       functions, the quadratures produced by this subroutine are
c       far inferior to the Gaussian quadratures; the former should
c       not be used in this environment.
c    Explanation: The functions this quadrature is supposed to
c       integrate have a (hopefully, integrable) singularity at
c       (or near) the point -1, and are smooth on the rest of the
c       interval [-1,0], including the point 0.
c 
c                     Input parameters:
c 
c  n - the number of nodes on the interval [-1,0] to be constructed;
c       should not be less than 14
c 
c                     Output parameters:
c 
c  xs - the n "left-universal" nodes on the interval [-1,0]
c  whts - the quadrature weights corresponding to the nodes xs
c 
c                     Work arrays:
c 
c  w - must be at least 4*n+60 real *8 locations long
c 
c 
c       . . . construct the Gaussian nodes and weights
c 
        ixs=1
        lxs=2*n+6
c 
        iwhts=ixs+lxs
c 
        itype=1
        call legeexps(itype,n*2,w(ixs),u,v,w(iwhts))
c 
c       construct the images of the Gaussian nodes and the
c       corresponding weights after prolate resampling
c 
        ndigits=15
        do 1200 i=1,n
c 
        call proresf(ndigits,w(ixs+i-1)+1,d,dd)
c 
        xs(n-i+1)=-d
        whts(n-i+1)=w(iwhts+i-1)*dd
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        SUBROUTINE cosexps(itype,n,ts,whts,u,v)
        implicit real *8 (a-h,o-z)
        save
        dimension ts(1),u(n,n),v(n,n),whts(1)
c 
c         this subroutine constructs the equispaced nodes
c         on the interval [-1,1], and the weights for the
c         corresponding order n quadrature (obviously, all
c         these weights are equal to 2/n). It also constructs
c         the matrix v converting the coefficients of a cosine
c         expansion into its values at the n equispaced nodes,
c         and its inverse u, converting the values of a function
c         at n equispaced nodes into the coefficients of the
c         corresponding cosine series. No attempt has been made
c         to make this code efficient, but its speed is normally
c         sufficient, and it is mercifully short.
c 
c                 input parameters:
c 
c  itype - the type of the calculation to be performed
c          itype=1 means that only the nodes and the weights
c                  are to be constructed
c          itype=2 means that the nodes, the weights, and
c                  the matrices u, v are to be constructed
c  n - the number of equispaced nodes and weights to be generated
c 
c                 output parameters:
c 
c  ts - the n equispaced nodes on the interval [-1,1] - computed
c          independently of the value of itype.
c  u - the n*n matrix converting the  values at of a cosine series
c         of order n-1 at n equispaced nodes into the coefficients
c         of its cosine expansion - computed only if itype=2
c  v - the n*n matrix converting the coefficients of an n-term
c         cosine expansion into its values at n equispaced nodes
c         (note that v is the inverse of u)  - computed only if
c         itype=2
c  whts - the corresponding quadrature weights - all equal to 2/n;
c         computed independently of the value of itype.
c 
c        . . . construct the n equispaced nodes on the interval [-1,1],
c              and the corresponding nodes on the interval [0,pi]
c              parametrizing space
c 
        done=1
        pi=atan(done)*4
        h=2*done/n
        do 1200 i=1,n
        ts(i)=-1+(i-1)*h+h/2
        whts(i)=h
 1200 continue
c 
        if(itype .eq. 1) return
c 
c       create the matrix of values of cosines at the equispaced nodes
c 
        do 1600 i=1,n
        do 1400 j=1,n
c 
        v(j,i)=cos((i-1)*(pi/2*ts(j)+pi/2) )
        u(i,j)=v(j,i)/n*2
c 
 1400 continue
 1600 continue
c 
        do 1800 j=1,n
        u(1,j)=u(1,j)/2
 1800 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cosfder(t,VAL,der,PEXP,N)
        implicit real *8 (a-h,o-z)
        save
        dimension pexp(1)
c 
c        This subroutine evaluates the user-supplied cosine
c        expansion (together with its derivative) at the
c        user-specified point x on the interval [-1,1]
c 
c 
c                input parameters:
c 
C     X = evaluation point
C     PEXP = expansion coefficients
C     N  = order of expansion
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c 
c                output parameters:
c 
C     VAL = computed value
C     der = computed value of the derivative
C 
C 
c 
c       . . . evaluate the cosine expansion at the point t
c 
        done=1
        pi=atan(done)*4
        x=pi/2*t+pi/2
c 
        val=0
        der=0
  
        sinmx=0
        cosmx=1
        sinx=sin(x)
        cosx=cos(x)
        do 1200 i=1,n+1
c 
  
cccc        val=val+pexp(i)*cos((i-1)*x)
cccc        der=der+pexp(i)*sin((i-1)*x)*(i-1)
  
        val=val+pexp(i)*cosmx
        der=der+pexp(i)*sinmx*(i-1)
c 
        sinnew=sinmx*cosx+cosmx*sinx
        cosnew=cosmx*cosx-sinmx*sinx
c 
        sinmx=sinnew
        cosmx=cosnew
 1200 continue
  
        der=-der*pi/2
c 
        return
        end
