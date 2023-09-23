        implicit real *8 (a-h,o-z)
        external whtfun1,funfunfu
        dimension w(100 000),roots(42 000),whts(42 000),
     1      info(10)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER maxnpts'
        READ *,maxnpts
        CALL PRINF('maxnpts=*',maxnpts,1 )
c 
        a=-1
        b=1
  
        a=0
        b=1
        eps=1.0d-14
c 
c        construct the orthonormal polynomials for the weight whtfun
c 
        lenw=100 000
        call arbgauss(ier,a,b,whtfun1,par1,par2,maxnpts,
     1    eps,roots,whts,w,lenw,info)
c 
        call prinf('after arbgauss, ier=*',ier,1)
        call prinf('after arbgauss, ltot=*',ltot,1)
        call prinf('after arbgauss, info=*',info,3)
  
c 
c        test the gaussian nodes and weights
c 
        call whtstest(roots,whts,maxnpts,a,b)
c 
        call prin2('incidentally, roots(maxnpts^2+1)=*',
     1      roots(maxnpts**2+1),200)
c 
        call prin2('incidentally, whts(maxnpts^2+1)=*',
     1      whts(maxnpts**2+1),200)
  
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine whtstest(roots,whts,maxnpts,a,b)
        implicit real *8 (a-h,o-z)
cccc        dimension roots(maxnpts+2,1),whts(maxnpts+2,1),
        save
        dimension roots(maxnpts,1),whts(maxnpts,1),
     1      rints(1000),rints2(1000),errs(1000)
        external funfunfu
c 
c        print the gaussian roots and corresponding weights
c 
        maxord=maxnpts+2
cccc        do 2000 i=2,maxord-1
         do 2000 i=1,maxord-2
c 
cccc        call prinf('number of roots is =*',i-1,1)
cccc        call prin2('and roots are*',roots(1,i),i-1)
cccc        call prin2('and weights are*',whts(1,i),i-1)
        call prinf('number of roots is =*',i,1)
        call prin2('and roots are*',roots(1,i),i)
        call prin2('and weights are*',whts(1,i),i)
 2000 continue
  
cccc              if(2 .ne. 3) stop
c 
c       now, check how well the obtained quadratures integrate
c       various powers of x
c 
        eps=1.0d-14
        m=20
c 
cccc        do 3000 i=2,maxord-1
        do 3000 i=1,maxord-2
c 
          call prinf('checking the i-point quadrature, i=*',i,1)
c 
c       integrate all  powers up to 2*i-th by using the designed
c       quadrature rule
c 
        do 2400 j=1,i*2
        d=0
cccc        do 2200 l=1,i-1
        do 2200 l=1,i
c 
        d=d+roots(l,i)**(j-1) * whts(l,i)
 2200 continue
        rints(j)=d
c 
c        now, evaluate the same integral via adapgaus
c 
        call adapgaus(ier,a,b,funfunfu,j-1,par2,m,eps,
     1      rints2(j),maxrec,numint)
c 
 2400 continue
  
c 
          call prin2('rints=*',rints,i*2)
          call prin2('rints2=*',rints2,i*2)
  
        do 2600 j=1,i*2
        errs(j)=rints(j)-rints2(j)
 2600 continue
        call prin2('and errs=*',errs,i*2)
  
  
  
 3000 continue
  
        return
        end
c 
c 
c 
c 
c 
        function funfunfu(x,i,par2)
        implicit real *8 (a-h,o-z)
c 
        save
        call whtfun1(x,par1,par2,f)
        funfunfu=f*x**i
        return
        end
c 
c 
c 
c 
c 
        subroutine whtfun1(x,par1,par2,f)
        implicit real *8 (a-h,o-z)
c 
cccc        f=dlog(x**2+1.0d-12)
        save
        f=dlog(x**2)
        f=f**2
cccc        f=1+dsin(x)
cccc        f=dexp(-x)
cccc        f=1
        return
        end
c 
c 
c 
c 
c 
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        this is the end of the debugging code and the
c        start of the actual routines for the construction
c        of Gaussian quadratures for reasonably general
c        weight functions.
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c 
c 
c 
c 
        subroutine arbgauss(ier,a,b,whtfun,par1,par2,maxnpts,
     1    eps,roots,whts,w,lenw,info)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),roots(1),whts(1),par1(1),par2(1),info(1)
c 
c       this subroutine constructs gaussian quadratures on the
c       interval [a,b] with the weight given by the (user-specified)
c       function whtfun. The function whtfun is fairly arbitrary;
c 
c                       input parameters:
c 
c  a,b - the ends of the interval on which the quadratures are to be
c       constructed
c  whtfun - the user-supplied subroutine evaluating the weight function
c       the calling sequence of whtfun must be
c 
c        call whtfun(x,par1,par2,f)                            (1)
c 
c        in (1), x is a point on the interval [a,b] where
c        the weight function is to be evaluated, and par1, par2
c        are two parameters to be used by whtfun; they can be
c        variables or arrays, real or integer, as desired.
c  par1, par2 - parameters to be used by the user-supplied
c        function fun (see above)
c  maxnpts  - the maximum number of nodes for which the quadratures
c        will be constructed; the subroutine constructs quadrature
c        nodes and weights for all m=1,2,...,maxnpts
c  eps - the accuracy (absolute) to which the integrals will be
c        evaluated; when whtsfun is substantially singular, care should
c        be exercised - a very small eps will cause a failure of the
c        subroutine
c  lenw - the length of the work array w in real *8 words
  
c                       output parameters:
c 
c  ier - error return code.
c          ier=0 means normal conclusion
c          ier=4 means that the amount of space allocated by the user
c                (given by the parameter lenw above) is insufficient.
c                it might also mean that the parameter eps is too
c                small.
c          ier=8 means that the amount of space allocated by the user
c                (given by the parameter lenw above) is insufficient.
c                it is unlikely to have other interpretations.
c          ier=32 means that the inverse power method for the calculation
c                of the gaussian weights has demonstrated slower
c                convergence than the subroutine is comfortable with.
c                something is seriously wrong with the function whtfun.
c                this is a fatal error.
c  roots - roots of the Gaussian quadrature, stored in the columns
c        of the two-dimensional array roots. Note that the array
c        roots is dimensioned roots(maxpnts,maxpnts), and only the
c        last column is full. The first column contains only one
c        element (the node of the oen-point quadrature), the second
c        column contains two elements, etc.
c  whts - weights of the Gaussian quadrature, stored in the columns
c        of the two-dimensional array whts. Note that the array
c        whts is dimensioned whts(maxpnts,maxpnts), and only the
c        last column is full. The first column contains only one
c        element (the weight of the one-point quadrature), the second
c        column contains two elements, etc.
c  info - a three-element integer array containing the following
c        information; info(1) - the number of subintervals into
c        which the subroutine subdivided  the interval [a,b] had
c        in oreder to integrate the weight function. info(2) is the
c        total number of nodes in the composite quadrature formula
c        constructed by the subroutine in order to calculate the
c        inner products required by the three-term recursion. info(3)
c        is the total amount of storage in array w that has been used
c        by the subroutine.
c 
c                       work arrays:
c 
c  w - must be long; the amount of storage available in w is
c        communicated to the subroutine by the input parameter lenw
c        (see above); if it is insufficient, the error return code
c        ier is set to eithe 4 or 8 (see above), and the execution
c        of the subroutine is terminated.
c 
c       . . . allocate memory for the subroutine arbincon
c 
        maxord=maxnpts+2
        k=20
        iw=1
        lw=10*k+8*k**2+300
c 
        iab=iw+lw
        lab=lenw-iab-2
c 
c       construct the discretization of the interval [a,b]
c       that will integrate the weight function properly
c 
        ifinit=1
c 
         call arbincon(ier,a,b,k,whtfun,par1,par2,eps,
     1    w(iab),lab,nn,labused,w(iw),ifinit)
c 
        info(1)=nn
c 
        call prinf('after arbincon, ier=*',ier,1)
c 
        if(ier .ne. 0) return
c 
        call prinf('and nn=*',nn,1)
c 
        if(ier .ne. 0) return
c 
c       . . . collect garbage
c 
        do 1200 i=1,nn*2
        w(i)=w(iab+i-1)
 1200 continue
        iab=1
        lab=nn*2+2
c 
c        allocate memory for the construction of orthogonal
c        polynomials - and the rest
c 
        kord=k+maxord
        n=kord*nn
c 
        info(2)=n
c 
        ix=iab+lab
        lx=n+2
c 
        iwhts=ix+lx
        lwhts=n+2
c 
        it=iwhts+lwhts
        lt=kord+2
c 
        iwww=it+lt
        lwww=kord+2
c 
        ifunwhts=iwww+lwww
        lfunwhts=n+2
c 
        ipkm1=ifunwhts+lfunwhts
        lpkm1=n+2
c 
        ipk=ipkm1+lpkm1
        lpk=n+2
c 
        ipkp1=ipk+lpk
        lpkp1=n+2
c 
        ialphas=ipkp1+lpkp1
        lalphas=maxord+4
c 
        ibetas=ialphas+lalphas
        lbetas=maxord+4
c 
        igammas=ibetas+lbetas
        lgammas=maxord+4
c 
        irgamma=igammas+lgammas
        lrgamma=maxord+4
c 
        irdelta=irgamma+lrgamma
        lrdelta=maxord+4
c 
        iaa=irdelta+lrdelta
        laa=maxord+4
c 
        ibb=iaa+laa
        lbb=maxord+4
c 
        icc=ibb+lbb
        lcc=maxord+4
c 
        iuu=icc+lcc
        luu=maxord+4
c 
        ivv=iuu+luu
        lvv=maxord+4
c 
        iww=ivv+lvv
        lww=maxord+4
c 
        ixkxk=iww+lww
        lxkxk=maxord+4
c 
        ltot=ixkxk+lxkxk
        info(3)=ltot
c 
        if(ltot .le. lenw) goto 2200
        ier=8
        return
 2200 continue
c 
        call arbgaus0(ier,a,b,whtfun,maxord,eps,k,n,
     1       w(iab),w(ix),w(iwhts),w(ifunwhts),w(ipkm1),w(ipk),
     2       w(ipkp1),w(ialphas),w(ibetas),w(igammas),roots,
     3       w(irgamma),w(irdelta),w(iaa),w(ibb),w(icc),
     4       w(iuu),w(ivv),w(iww),whts,w(ixkxk),w(it),
     5       w(iwww),nn,par1,par2)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine arbgaus0(ier,a,b,whtfun,maxord,eps,k,n,
     1    ab,x,whts,funwhts,pkm1,pk,pkp1,alphas,betas,gammas,
     2    roots,rgamma,rdelta,aa,bb,cc,uu,vv,ww,whtswhts,xkxk,
     3    t,www,nn,par1,par2)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension ab(2,1),x(1),whts(1),funwhts(1),pkm1(1),
     1    pkp1(1),pk(1),alphas(1),betas(1),gammas(1),
     2    rgamma(1),rdelta(1),par1(1),par2(1),
     3      aa(1),bb(1),cc(1),uu(1),vv(1),ww(1),t(1),www(1),
cccc     4      whtswhts(maxord,1),roots(maxord,1),xkxk(1)
     4      whtswhts(maxord-2,1),roots(maxord-2,1),xkxk(1)
c 
c        construct the actual nodes on the interval [a,b]
c 
        kord=maxord+k
         call prinf('in arbgauss, kord=*',kord,1)
        call arbnodes(ab,nn,kord,x,whts,t,www)
c 
        do 1400 i=1,n
c 
        call whtfun(x(i),par1,par2,f)
        funwhts(i)=whts(i)*f
 1400 continue
c 
c        construct the orthogonal polynomials with the weight whtfun
c 
        ifvect=0
        call arbortho(x,funwhts,n,maxord,
     1      alphas,betas,gammas,pkm1,pk,pkp1,vect,ifvect,
     2      rgamma,rdelta,rint)
c 
c       one order after another, find the roots
c 
cccc        do 3200 m=2,maxord-1
        do 3200 m=1,maxord-2
        mm1=m-1
        if(m .eq. 1) mm1=1
c 
        call arbfndrt(alphas,betas,gammas,
cccc     1    m,roots(1,m-1),roots(1,m),a,b,pk,pkm1,pkp1)
     1    m+1,roots(1,mm1),roots(1,m),a,b,pk,pkm1,pkp1)
c 
 3200 continue
c 
c       now, for each order, find the quadrature weights
c 
cccc        do 3600 i=2,maxord-1
        do 3600 i=1,maxord-2
c 
cccc        call arbwhts(ier,rgamma,rdelta,roots(1,i),i-1,
        call arbwhts(ier,rgamma,rdelta,roots(1,i),i,
     1    whtswhts(1,i),aa,bb,cc,xkxk,uu,vv,ww,rint)
c 
        if(ier .ne. 0) return
 3600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine arbwhts(ier,rgamma,rdelta,roots,n,
     1    whts,a,b,c,xk,u,v,w,rint)
        implicit real *8 (a-h,o-z)
        save
        dimension rgamma(1),rdelta(1),roots(1),
     1      a(1),b(1),c(1),u(1),v(1),w(1),xk(1),whts(1),ds(20)
c 
c        this subroutine constructs the weights of the Gaussian
c        quadrature from the ceofficients of the three-term
c        recursion and the roots of the orthogonal polynomial
c 
c                        input parameters:
c 
c  rgamma, rdelta - coefficients of the three-term recursion as
c        produced by the subroutine arbortho
c  roots - the roots of the n-the orthogonal polynomial
c  rint - the integral of the weight function over the interval
c        where the polynomials are orthogonal
c 
c                        output parameters:
c 
c  whts - the Gaussian weights corresponding to the Gaussian
c        nodes roots
c 
c                        work arrays:
c 
c  a,b,c,xk,u,v,w - each must be at least n+1 real *8 locations long
c 
c 
c        . . . from the coefficients of the recursion (supplied
c              by the user) construct the elements of the
c              tridiagonal matrix whose eigenvalues are the
c              gaussian nodes
c 
        ier=0
        done=1
        numit=6
        eps=1.0d-10
c 
        do 3000 k=1,n
c 
        do 1200 i=1,n
        a(i)=rdelta(i) - roots(k)*(1+eps)
        b(i)=rgamma(i)
        c(i+1)=b(i)
 1200 continue
c 
c       l-u factor this matrix
c 
        call arbthrdf(a,b,c,n,u,v,w)
c 
c       repeatedly apply the inverse to a vector
c 
        do 1400 i=1,n
        xk(i)=1
 1400 continue
c 
        do 2000 i=1,numit
c 
        call arbthrds(u,v,w,n,xk)
c 
        d=0
        do 1600 j=1,n
        d=d+xk(j)**2
 1600 continue
c 
        d=done/dsqrt(d)
        ds(i)=d
        do 1800 j=1,n
        xk(j)=xk(j)*d
 1800 continue
c 
        dlook=d
 2000 continue
c 
        if(dlook .gt. eps*1.0d4) ier=32
        if(ier .eq. 32) return
        whts(k)=xk(1)**2*rint
c 
 3000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine arbthrdf(a,b,c,n,u,v,w)
        implicit real *8 (a-h,o-z)
        save
        dimension a(1),b(1),c(1),u(1),v(1),w(1),rhs(1)
c 
c        eliminate down
c 
        do 1200 i=1,n-1
        d=c(i+1)/a(i)
        a(i+1)=a(i+1)-b(i)*d
        u(i)=d
 1200 continue
c 
c        eliminate up
c 
        do 1400 i=n,2,-1
        d=b(i-1)/a(i)
        v(i)=d
 1400 continue
c 
c       scale the diagonal
c 
        done=1
        do 1600 i=1,n
        w(i)=done/a(i)
 1600 continue
c 
        return
c 
c 
c 
c 
        entry arbthrds(u,v,w,n,rhs)
c 
c        eliminate down
c 
        do 2400 i=1,n-1
        rhs(i+1)=rhs(i+1)-u(i)*rhs(i)
 2400 continue
c 
c        eliminate up
c 
        do 2600 i=n,2,-1
        rhs(i-1)=rhs(i-1)-rhs(i)*v(i)
 2600 continue
c 
c       scale
c 
        do 2800 i=1,n
        rhs(i)=rhs(i)*w(i)
 2800 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine arbfndrt(alphas,betas,gammas,
     1    m,rootsold,roots,a,b,sepa,p,ders)
        implicit real *8 (a-h,o-z)
        save
        dimension alphas(1),betas(1),gammas(1),
     1    roots(1),rootsold(1),sepa(1),p(1),ders(1)
c 
c        this subroutine finds all roots of the orthogonal polynomial
c        of order m-1, given by the three-term recursion that has
c        generated this polynomial.
c 
c                          input parameters:
c 
c  alphas,betas,gammas - the coefficients of the three-term relation
c        as produced by the subroutine arbortho (see).
c  m - the order of the orthogonasl polynomial whose roors are to
c        be found is m-1
c  rootsold - the roots of the orthogonal polynomial of order m-2
c        (not used if m=2), to be used to separate the roots to
c        be found
c  a,b - the ends of the interval on which the polynomials are
c        orthogonal
c 
c                          output parameters:
c 
c  roots - the m-1 roots of the orthogonal polynomial of order m-1
c 
c                          work arrays:
c 
c  sepa - must be at least m+2 real *8 locations long
c 
c       . . . construct the array of separators
c 
        sepa(1)=a
        sepa(m)=b
c 
        if(m .le. 2) goto 1400
c 
        do 1200 i=1,m-2
        sepa(i+1)=rootsold(i)
 1200 continue
c 
 1400 continue
c 
c       use bisection to find roots on each of the intervals
c 
        do 2000 i=1,m-1
c 
        call arbnewt(alphas,betas,gammas,sepa(i),
     1      sepa(i+1),m,roots(i),p,ders)
c 
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine arbnewt(alphas,betas,gammas,a,b,m,x,
     1      p,ders)
        implicit real *8 (a-h,o-z)
        save
        dimension alphas(1),betas(1),gammas(1),p(1),ders(1)
c 
c        use ten bisection steps to better to localize the root
c 
        ak=a
        bk=b
c 
        call arbpeva(ak,alphas,betas,gammas,
     1      p,m+1,ders)
c 
        fak=p(m)
c 
        call arbpeva(bk,alphas,betas,gammas,
     1      p,m+1,ders)
c 
        fbk=p(m)
c 
        do 2000 i=1,9
c 
        ck=(ak+bk)/2
c 
        call arbpeva(ck,alphas,betas,gammas,
     1      p,m+1,ders)
c 
        fck=p(m)
c 
        d=fck*fak
        if(d .lt. 0) goto 1400
        ak=ck
        fak=fck
        goto 2000
c 
 1400 continue
c 
        bk=ck
        fbk=fck
 2000 continue
        x=(ak+bk)/2
c 
c       now, tighten the result by Newton
c 
        xk=x
        do 2400 i=1,5
c 
        call arbpeva(xk,alphas,betas,gammas,
     1      p,m+1,ders)
c 
        fxk=p(m)
        derfxk=ders(m)
        xkp1=xk-fxk/derfxk
        xk=xkp1
 2400 continue
c 
        x=xk
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine arbpeva(x,alphas,betas,gammas,
     1      pols,m,ders)
        implicit real *8 (a-h,o-z)
        save
        dimension alphas(1),betas(1),gammas(1),pols(1),
     1      ders(1)
c 
c       this subroutine evaluates the first m orthogonal polynomials
c       and their derivatives at the point x, using the coefficients
c       (in a somewhat goofey form) of the three-term recursion
c 
c                            input parameters:
c 
c  x - the point where the polynomials are to be evaluated
c  alphas,betas,gammas - the recursiuon coefficients as produced
c       by the subroutine arbortho (see)
c  m - the number of polynomials to be evaluated
c 
c                            output parameters:
c 
c  pols - the orthogonal polynomials evalueted at the point x
c  ders - the derivative of the orthogonal polynomials evaluated
c       at the point x
c 
        pols(1)=gammas(1)
        ders(1)=0
c 
        pols(2)=x*pols(1)-pols(1)*betas(2)
        pols(2)=pols(2)*gammas(2)
        ders(2)=pols(1)*gammas(2)
c 
        do 1200 i=2,m
        pols(i+1)=x*pols(i)-alphas(i+1)*pols(i-1)-
     1     betas(i+1)*pols(i)
        pols(i+1)=pols(i+1)*gammas(i+1)
c 
        ders(i+1)=x*ders(i)+pols(i)-alphas(i+1)*ders(i-1)-
     1     betas(i+1)*ders(i)
        ders(i+1)=ders(i+1)*gammas(i+1)
c 
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine arbortho(x,whts,n,maxord,
     1      alphas,betas,gammas,pkm1,pk,pkp1,vectors,
     2      ifvect,rgamma,rdelta,rint)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),whts(1),alphas(1),betas(1),
     1      gammas(1),pkm1(1),pk(1),pkp1(1),vectors(n,1),
     2      rgamma(1),rdelta(1)
c 
c        this subroutine produces the three-term recursion generating
c        the first maxord orthogonal polynomials corresponding to
c        a user-provided weight function. actually, the polynomials
c        are orthogonal on a discrete set of nodes x with respect
c        to the weights in array funwhts; creating the appropriate arrays
c        x, w is the user's responsibility, and is normally performrd
c        by the subroutines arbincon, arbnodes (see). the subroutine
c        produces three types of output:
c 
c  1. The somewhat goofey version of the three-term recursion, given
c        by the arrays alphas, betas, gammas; these are used by the
c        subroutine arbpeva to evaluate the constructed orthogonal
c        polynomials and their derivatives at arbitrary points
c  2. The coefficients of the standard three-term recursion, as
c        described in Stoer and Bulirsch, for example. these are used
c        to construct the weights corresponding to the roots of the
c        polynomials produced by this subroutine
c  3. (optional) If the user sets the parameter ifvect (see below)
c        to 1, the subroutine will return the values at the nodes x
c        of all the orthogonal polynomials it has created, in the
c        array vect (see). this is principally a debugging option.
c 
c                             input parameters:
c 
c  x - the nodes on which the polynomials will be orthogonal
c  funwhts - the weights with respect to which the polynomials
c        will be orthogonal
c  n - the number of points in arrays x, funwhts
c  maxord - the maximum order of the orthogonal polynomials to be created
c  ifvect - tells the subroutine if the values of the orthononal
c        polynomials are to be returned in array vect.
c 
c                             output parameters:
c 
c  alphas,betas,gammas - arras to be used by the subroutine
c        arbpeva to evaluate the constructed orthogonal polynomials and
c        their derivatives at arbitrary points. each is maxord+2
c        real *8 elements long
c  vectors - the values of the orthonormal polynomials (maxord+1
c        of them) at the points x
c  rgamma,rdelta - the coefficients of the standard thee-term recursion
c        (see, for example, Stoer and Bulirsch), to be used by the
c        subroutine arbwhts to find the weights of the gGaussian
c        quadrature
c  rint - the integral of the weight function over the interval where
c        the polynomials are orthogonal
c 
c                             work arrays:
c 
c  pkm1,pk,pkmp1 - must be at least n+2 rteal *8 elements each
c 
c 
c       . . . initialize the orthogonalization procedure
c 
        done=1
        do 1400 i=1,n
        pk(i)=1
        pkm1(i)=0
 1400 continue
c 
        call arbscap(pk,pk,whts,n,d)
        rint=d
c 
        d=done/dsqrt(d)
        do 1600 i=1,n
        pk(i)=pk(i)*d
 1600 continue
        gammas(1)=d
c 
c        conduct the three-term recursion
c 
        do 3000 i=1,maxord
c 
c        . . . construct the new vector
c 
        do 2200 j=1,n
        pkp1(j)=pk(j)*x(j)
 2200 continue
c 
c       ortogonalize it to the preceding two
c 
        call arbscap(pkp1,pkm1,whts,n,alpha)
        call arbscap(pkp1,pk,whts,n,beta)
c 
        rdelta(i)=beta
c 
        do 2400 j=1,n
        pkp1(j)=pkp1(j)-(alpha*pkm1(j)+beta*pk(j))
 2400 continue
c 
c       normalize the obtained vector
c 
        call arbscap(pkp1,pkp1,whts,n,gamma)
c 
        rgamma(i)=dsqrt(gamma)
c 
        gamma=done/dsqrt(gamma)
        do 2600 j=1,n
        pkp1(j)=pkp1(j)*gamma
 2600 continue
c 
c       memorize alpha, beta, and gamma, and shift pk, pkp1
c 
        alphas(i+1)=alpha
        betas(i+1)=beta
        gammas(i+1)=gamma
c 
        do 2800 j=1,n
        pkm1(j)=pk(j)
        pk(j)=pkp1(j)
c 
 2800 continue
c 
c       if the user so requested, store the current orthogonal
c       polynomial in the array vectors
c 
        if(ifvect .eq. 0) goto 3000
        do 2900 j=1,n
        vectors(j,i)=pkm1(j)
 2900 continue
c 
 3000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine arbscap(x,y,w,n,d)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1),w(1)
c 
        d=0
        do 2000 i=1,n
        d=d+x(i)*y(i)*w(i)
 2000 continue
        return
        end
  
c 
c 
c 
c 
c 
        subroutine arbnodes(ab,nn,k,x,whts,t,w)
        implicit real *8 (a-h,o-z)
        save
        dimension ab(2,1),x(1),whts(1),t(1),w(1)
c 
c       this subroutine constructs the nodes and weights on
c       corresponding to the adptive structure boilt by the
c       subroutine arbincon (see), and given by the array ab
c       of left and right ends of the intervals
c 
c                        input parameters:
c 
c  ab - the array of left and right ends of all subintervals in the
c       structure, as constructed by subroutine arbincon (see)
c  nn - the number of intervals in the structure
c  k - the number of gaussian nodes on each subinterval in the structure
c 
c                        output parameters:
c 
c  x - the array of all gaussian nodes on all subintervals (nn*k of them)
c  whts - the weights corresponding to the nodes x (nn*k of them)
c 
c                        work arrays:
c 
c  t,w - must be at least k+2 real *8 elements each
c 
c       construct the nodes and weights of the k-point Gaussian
c       quadrature on the interval [-1,1]
c 
        itype=1
        iwht=k+2
c 
        call arbexps(itype,k,t,u,v,w)
c 
c       one subinterval after another, construct nodes and weights
c       of the quadrature on the big interval
c 
        do 2000 i=1,nn
        u=(ab(2,i)-ab(1,i))/2
        v=(ab(2,i)+ab(1,i))/2
        ik=(i-1)*k
c 
        do 1200 j=1,k
        x(ik+j)=t(j)*u+v
        whts(ik+j)=w(j)*u
 1200 continue
 2000 continue
c 
        return
        end
c 
c 
c 
c 
c 
         subroutine arbincon(ier,a,b,k,fun,par1,par2,eps,
     1    ab,lab,nn,labused,w,ifinit)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),par1(1),par2(1),ab(2,1)
c 
c       this subroutine finds a subdivision of the interval [a,b],
c       such that on every subinterval, the user-provided function
c       fun is integrated to precision eps by a k-point Gaussian
c       quadrature. Actually, this is a memory management routine;
c       all the work is done by the routine arbinter (see).
c 
c                          input parameters:
c 
c  a,b - the ends of the interval to be subdivided
c  k - the number of Gaussian nodes on each subinterval
c  fun - the user-supplied function to be integrated
c      The calling sequence of funget must be:
c 
c      funget(x,par1,par2,f),
c 
c      where x is the input argument where the function is to be
c      evaluated, f  (output) is the function, and par1,par2 are
c      both input parameters, carrying whatever information funget
c      might need.
c 
c  par1,par2 - whatever the user-supplied subroutine funget needs
c  eps - the accuracy to which the legendre expansion must represent
c     the user-supplied function for the subroutine to decide that
c     the representation is accurate
c  lab - the length of the user-provided array ab (see below)
c  w - the array containing the data to be used by the subroutine,
c     created by a preceding call to it. It is only an input parameter
c     if the user has set the parameter ifinit (see below) to 0.
c     Otherwise, it is an output parameter, and must be at
c     least 10*k+8*k**2+30 real *8 locations long.
c  ifinit - the index telling the subroutine whether it should construct
c     the parameter w (see above), or view it as an input parameter.
c     In the latter case, this  must have been created by a prior call
c     to this subroutine with the same parameter k.
c     ifinit=1 means that the parameter w will be created.
c     ifinit=0 means that the parameter w will be viewed
c     as an input parameter.
c 
c                     output parameters:
c 
c  ier - error return code.
c         ier=0 means successful execution
c         ier=4 means that the space allocated in the array ab is
c         insufficient
c  ab - the array of ends of subintervals created. For each i=1,2,...,nn,
c     ab(1,i) is the left end of the i-th subinterval; ab(2,i) is
c     the right end of the i-th subinterval
c  nn - the number of subintervals created on the interval [a,b];
c     also the number of entries in the array ab (see above)
c  labused - the total amount of space used by the subroutine in the
c     array ab; note that normally, it  is somewhat greater than
c     the expected nn*2; The difference is not large, but should
c     be taken into account when allocating memory.
c  w - the array containing the data used by the subroutine,
c     It is only an output parameter if the user has set the
c     parameter ifinit (see above) to 1. Otherwise, it is an input
c     parameter, and must have been created  by a prior call to
c     this subroutine. Note that it should not be altered between
c     consecutive calls to this subroutine; the piece
c     used by this subroutine is 10*k+8*k**2+30 real *8 locations long.
c 
c       . . . allocate memory for the subroutine arbinter
c 
        ier=0
c 
        it=1
        lt=k*2+2
c 
        ix=it+lt
        lx=k*2+2
c 
        if=ix+lx
        lf=k*2+2
c 
        icoefs=if+lf
        lcoefs=k*2+2
c 
        iwhts=icoefs+lcoefs
        lwhts=k*2+1
c 
        iu=iwhts+lwhts
        lu=k**2*4 +10
c 
        iv=iu+lu
        lv=k**2*4 +10
c 
c       recursively subdivide the interval [a,b] into subintervals
c       such that on each of them, the k-point Gaussian rule
c       will obtain the accuracy epc
c 
         call arbinter(ier,a,b,fun,par1,par2,k,ab,lab,
     1    eps,w(it),w(iu),w(iv),w(iwhts),w(ix),w(if),
     2       w(icoefs),nn,labused,ifinit)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine arbbubbl(ab,nn)
        implicit real *8 (a-h,o-z)
        save
        dimension ab(2,1)
c 
c       sort the array ab with respect to the first coordinate
c 
        do 2000 i=1,nn
        do 1200 j=1,i-1
        if(ab(1,i) .ge. ab(1,j) ) goto 1200
        d=ab(1,i)
        ab(1,i)=ab(1,j)
        ab(1,j)=d
c 
        d=ab(2,i)
        ab(2,i)=ab(2,j)
        ab(2,j)=d
 1200 continue
 2000 continue
        return
        end
c 
c 
c 
c 
c 
         subroutine arbinter(ier,a,b,fun,par1,par2,k,ab,lab,
     1    eps,t,u,v,whts,x,f,coefs,nn,lused,ifinit7)
        implicit real *8 (a-h,o-z)
        save
        dimension t(1),u(2*k,2*k),v(2*k,2*k),whts(1),x(1),
     1      f(1),coefs(1),par1(1),par2(1),ab(2,1)
        equivalence (rea,intint)
c 
c       this subroutine finds a subdivision of the interval [a,b],
c       such that on every subinterval, the user-provided function
c       fun is integrated to precision eps by a k-point Gaussian
c       quadrature.
c 
c                          input parameters:
c 
c  a,b - the ends of the interval to be subdivided
c  fun - the user-supplied function to be integrated
c      The calling sequence of funget must be:
c 
c      funget(x,par1,par2,f),
c 
c      where x is the input argument where the function is to be
c      evaluated, f  (output) is the function, and par1,par2 are
c      both input parameters, carrying whatever information funget
c      might need.
c 
c  par1,par2 - whatever the user-supplied subroutine funget needs
c  k - the number of Gaussian nodes on each subinterval
c  lab - the length of the user-provided array ab (see below)
c  eps - the accuracy to which the legendre expansion must represent
c     the user-supplied function for the subroutine to decide that
c     the representation is accurate
c  t - the Gaussian nodes on the interval [-1,1]. This is only an
c     input parameter is ifinit7 (see below) has been set to 0;
c     otherwise, it is an output parameter.
c  u - the 2*k * 2*k  matrix converting the  values at of a polynomial
c     of order 2*k-1 at 2*k legendre nodes into the coefficients of its
c     legendre expansion. This is only an input parameter is ifinit7
c     (see below) has been set to 0; otherwise, it is an output parameter.
c  v - the 2*k * 2*k matrix converting the coefficients
c     of an 2*k-term legendre expansion into its values at 2*k legendre
c     nodes (note that v is the inverse of u). This is only an input
c     parameter is ifinit7 (see below) has been set to 0; otherwise, it
c     is an output parameter.
c  whts - the Gaussian quadrature weights (2*k of them).  This is only
c     an input parameter is ifinit7 (see below) has been set to 0;
c     otherwise, its an output parameter.
c  ifinit7 - the index telling the subroutine whether it should construct
c     the parameters t,u,v,whts (by calling the subroutine arbexps),
c     or view them as input parameters. In the latter case, these
c     must have been created by a prior call to this subroutine (or to
c     subroutine arbexps, for an adventurous and creative user).
c     ifinit7=1 means that the parameters t,u,v,whts will be created.
c     ifinit7=0 means that the parameters t,u,v,whts will be viewed
c     as input parameters.
c 
c                     output parameters:
c 
c  ier - error return code.
c         ier=0 means successful execution
c         ier=4 means that the space allocated in the array ab is
c         insufficient
c  ab - the array of ends of subintervals created. For each i=1,2,...,nn,
c     ab(1,i) is the left end of the i-th subinterval; ab(2,i) is
c     the right end of the i-th subinterval
c  t,u,v,whts - are only output is the parameter ifinit (see above)
c     has been set by the user to 1.
c  x - the Gaussian nodes on the interval [a,b]
c  f - the values of the user-specified function fun at
c     the nodes x
c  coefs - the coefficients of the Legendre expansion of the
c     function fun on the interval [a,b]
c  nn - the number of subintervals created on the interval [a,b];
c     also the number of entries in the array ab (see above)
c  lused - the total amount of space used by the subroutine in the
c     array ab; note that normally, it  is somewhat greater than
c     the expected nn*2; The difference is not large, but should
c     be taken into account when allocating memory.
c 
c       . . . start the recursive subdivision process
c 
        ab(1,1)=a
        ab(2,1)=b
        nn=0
        ifinit=ifinit7
        istore=1
        ier=0
c 
c       recursively subdivide the intervals till none need
c       to be subdivided
c 
        lused=0
        do 2000 i=1,1 000 000
c 
        nn=nn+1
c 
c       check if any more intervals need to be subdivided
c 
        if(nn .gt. istore) goto 2200
c 
c       check if this interval needs to be subdivided
c 
        d=ab(2,nn)-ab(1,nn)
c 
        call arbsplit(ier,ab(1,nn),ab(2,nn),k,fun,par1,par2,
     1    eps/d,t,u,v,whts,x,ifinit,f,coefs)
c 
        ifinit=0
c 
        if(ier .eq. 0) goto 2000
c 
c       this interval needs to be subdivided. act accordingly
c 
        if(lused .lt. istore*2+4) lused=istore*2+4
        if(istore*2 +4 .lt. lab) goto 1800
        ier=4
        return
 1800 continue
c 
        istore=istore+1
        ab(1,istore)=ab(1,nn)
        ab(2,istore)=(ab(1,nn)+ab(2,nn))/2
c 
        istore=istore+1
        ab(1,istore)=ab(2,istore-1)
        ab(2,istore)=ab(2,nn)
c 
        intint=9654236
        ab(1,nn)=rea
        ab(2,nn)=1.0d35
 2000 continue
 2200 continue
c 
        nn=nn-1
c 
c       scan the intervals we have created, discarding those that have kids
c 
        ii=0
        do 2400 i=1,nn
        rea=ab(1,i)
c 
        if( (ab(2,i) .gt. 1.0d34 ) .and.
     1      (intint .eq. 9654236) ) goto 2400
        ii=ii+1
        ab(1,ii)=ab(1,i)
        ab(2,ii)=ab(2,i)
 2400 continue
        nn=ii
c 
c        bubble-sort the array ab
c 
        call arbbubbl(ab,nn)
        return
        end
c 
c 
c 
c 
c 
        subroutine arbsplit(ier,a,b,k,funget,par1,par2,eps,
     1    t,u,v,whts,x,ifinit,f,coefs)
        implicit real *8 (a-h,o-z)
        save
        dimension t(1),u(2*k,2*k),v(2*k,2*k),whts(1),x(1),
     1      f(1),coefs(1),par1(1),par2(1)
c 
c       this subroutine determines whether a k-point legendre
c       expansion of the user-supplied function f on the interval
c       [a,b] represents it with accuracy eps.
c 
c                          input parameters:
c 
c  a,b - the ends of the interval on which the approximation is checked
c  k - the order of the approximation for which the accuracy is checked
c  funget - the user-supplied function for which the accuracy is checked
c      The calling sequence of funget must be:
c 
c      funget(x,par1,par2,f),
c 
c      where x is the input argument where the function is to be
c      evaluated, f  (output) is the function, and par1,par2 are
c      both input parameters, carrying whatever information funget
c      might need.
c 
c  par1,par2 - whatever the user-supplied subroutine funget needs
c  eps - the accuracy to which the legendre expansion must represent
c     the user-supplied function for the subroutine to decide that
c     the representation is accurate
c  t - the Gaussian nodes on the interval [-1,1]. This is only an
c     input parameter is ifinit (see below) has been set to 0;
c     otherwise, it is an output parameter.
c  u - the 2*k * 2*k  matrix converting the  values at of a polynomial
c     of order 2*k-1 at 2*k legendre nodes into the coefficients of its
c     legendre expansion. This is only an input parameter is ifinit
c     (see below) has been set to 0; otherwise, it is an output parameter.
c  v - the 2*k * 2*k matrix converting the coefficients
c     of an 2*k-term legendre expansion into its values at 2*k legendre
c     nodes (note that v is the inverse of u). This is only an input
c     parameter is ifinit (see below) has been set to 0; otherwise, it
c     is an output parameter.
c  whts - the Gaussian quadrature weights (2*k of them).  This is only
c     an input parameter is ifinit (see below) has been set to 0;
c     otherwise, its an output parameter.
c  ifinit - the index telling the subroutine whether it sholud construct
c     the parameters t,u,v,whts (by calling the subroutine arbexps),
c     or view them as input parameters. In the latter case, these
c     must have been created by a prior call to this subroutine (or to
c     subroutine arbexps, for an adventurous and creative user).
c     ifinit=1 means that the parameters t,u,v,whts will be created.
c     ifinit=0 means that the parameters t,u,v,whts will be viewed
c     as input parameters.
c 
c                     output parameters:
c 
c  t,u,v,whts - are only output is the parameter ifinit (see above)
c     has been set by the user to 1.
c  x - the Gaussian nodes on the interval [a,b]
c  f - the values of the user-specified function funget at
c     the nodes x
c  coefs - the coefficients of the Legendre expansion of the
c     function f on the interval [a,b]
c 
c       . . . if the need be, construct the gaussian nodes on the
c             interval [-1,1], and the matrix conevrting the values
c             of the function at these nodes into coefficients of
c             its legendre expansion
c 
        if(ifinit .eq. 0) goto 1200
c 
        itype=2
        call arbexps(itype,k*2,t,u,v,whts)
c 
 1200 continue
c 
c        discretize the interval and evaluate the
c        user-supplied function at its nodes
c 
        alpha=(b-a)/2
        beta=(b+a)/2
c 
        do 1400 i=1,k*2
        x(i)=alpha*t(i)+beta
        call funget(x(i),par1,par2,f(i) )
 1400 continue
c 
c       construct the legendre expansion of the function on
c       the interval
c 
        call arbmatv(u,f,coefs,k*2)
c 
c       scan the coefficients and see if the last k
c       of them are small
c 
        ier=0
        do 1800 i=1,k
        d=dabs(coefs(k+i))
        if(d .gt. eps) ier=4
 1800 continue
c 
        return
c 
        end
c 
c 
c 
c 
        subroutine arbmatv(a,x,y,n)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),x(1),y(1)
c 
        do 2000 i=1,n
        d=0
        do 1400 j=1,n
        d=d+a(i,j)*x(j)
 1400 continue
        y(i)=d
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine arbexps(itype,n,x,u,v,whts)
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
        call arblwhts(n,x,whts,ifwhts)
c 
c       construct the matrix of values of the legendre polynomials
c       at these nodes
c 
        if(itype .ne. 2) return
        do 1400 i=1,n
c 
        call arblpols(x(i),n-1,u(1,i) )
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
        subroutine arblwhts(n,ts,whts,ifwhts)
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
        call arblpol(xk,n,pol,der)
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
        call arbprode(a,ts,n,i,fm)
        call arbprode(b,ts,n,i,fp)
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
        subroutine arblpol(x,n,pol,der)
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
        subroutine arbprode(x,xs,n,i,f)
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
        subroutine arblpols(x,n,pols)
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
  
