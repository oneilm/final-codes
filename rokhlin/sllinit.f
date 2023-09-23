        implicit real *8 (a-h,o-z)
        real *8 ts(5000),whts(1000),xs(10 000),whtfun(10 000)
c 
        real *8 w(2 000 000),rlams(10000),rlams2(10000),
     1      fs(10000),whts0(10000)
c 
        dimension coefsout(10 000),apprrhs(10 000)
        dimension srccoe(10000),pattcoe(10000)
c 
        dimension pattern(10 000),ts2(10 000),
     1      pjcoefs1(10 000),pjcoefs2(10 000),srcrea(10000),
     2      srcima(10 000),whts2(10 000),rfinpatt(10 000),
     3      apprrhs1(10000),apprrhs2(10000)
  
c 
        complex *16 cdd,ima,srcvals(10 000),
     1      finpatt(10 000)
c 
        data ima/(0.0d0,1.0d0)/
c 
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
c 
        PRINT *, 'ENTER c'
        READ *,c
        CALL PRIN2('c=*',c,1 )
c 
        nn=300
cccc        nn=20
        b=0.5
        a=-b
c 
c       construct the Gaussian nodes and weights on the interval [-1,1]
c 
        itype=1
        call legeexps(itype,nn,ts,u,v,whts)
c 
        aa=-1
        call intedisc(ts,whts,nn,aa,a,xs,whts0)
        call intedisc(ts,whts,nn,a,b,xs(nn+1),whts0(nn+1))
        bb=1
        call intedisc(ts,whts,nn,b,bb,xs(2*nn+1),whts0(2*nn+1))
c 
        n=nn*3
c 
        do 1400 i=1,nn
        fs(i)=0
        fs(nn+i)=2 * xs(nn+i)
        fs(2*nn+i)=0
 1400 continue
c 
        d=0
        do 1600 i=1,n
c 
        whtfun(i)=1
c 
        if(i .le. nn-10) whtfun(i)=1000
        if(i .gt. 2*nn+10) whtfun(i)=1000
 1600 continue
  
        call prin2('whts0 as created*',whts0,n)
        call prin2('whtfun as created*',whtfun,n)
c 
c        initialize the least squares procedure
c 
        lenw=2 000 000
        call sllinit(ier,c,xs,whtfun,whts0,n,w,lenw,m,rlams,
     1      rlams2,keep,nhigh,lused)
c 
        call prinf('after sllinit, keep=*',keep,1)
        call prinf('after sllinit, nhigh=*',nhigh,1)
        call prinf('after sllinit, lused=*',lused,1)
c 
c       use the prepared data to construct the least squares
c       approximation
c 
        sgain=20
        ifsecant=1
c 
        call slleval(ier,w,fs,sgain,ifsecant,coefsout,
     1      sgainout,apprrhs,err,errnosu,srccoe,pattcoe,nlege)
  
        call prinf('after slleval, ier=*',ier,1)
c 
        call prin2('after slleval, pattcoe=*',pattcoe,nlege)
c 
c       plot the right-hand side and the approximation together
c 
        iw=21
        call lotagraph2(iw,xs,fs,n,xs,apprrhs,n,
     1      'approximation together with fs*')
c 
c       use the Legendre expansions of the pattern and the
c       source that creates it; plot the results
c 
        itype=1
        call legeexps(itype,n,ts2,u,v,whts2)
c 
        ninit=1000
        do 2200 i=1,n
c 
        call legeexe2(ts2(i),pattern(i),pattcoe,Nlege,
     1      pjcoefs1,pjcoefs2,ninit)
  
        d=0
        do 2100 j=1,m
c 
        d=d+val*coefsout(j)
 2100 continue
c 
        ninit=0
c 
 2200 continue
c 
c       plot the right-hand side and the approximation together
c 
        iw=22
        call lotagraph2(iw,xs,fs,n,ts2,pattern,n,
     1      'approximation from pattcoe together with fs*')
c 
c       calculate and plot the real and imaginary parts
c       of the source
c 
        do 2400 i=1,n
c 
        call sgainleg(ts2(i),cdd,srccoe,Nlege)
  
        srcrea(i)=cdd
        srcima(i)=cdd*ima
  
        srcvals(i)=cdd
  
 2400 continue
  
c 
c       plot the right-hand side and the approximation together
c 
        iw=23
        call lotagraph2(iw,ts2,srcrea,n,ts2,srcima,n,
     1      'real and imaginary parts of source*')
c 
c        use numerical integration to convert the source into
c        the pattern - the last test
c 
        do 3200 i=1,n
c 
        call prointev(c,ts2,whts2,srcvals,ts2(i),n,finpatt(i) )
  
        rfinpatt(i)=finpatt(i)
 3200 continue
  
        rfinpatt(1)=2
c 
        iw=24
        call lotagraph2(iw,ts2,pattern,n,ts2,rfinpatt,n,
     1      'pattern computed two ways*')
c 
c       finally, design patterns with several different
c       supergains and plot them things together
c 
  
        sgain=2
        sgain=1.0001
cccc        sgain=1.01
        ifsecant=1
c 
        call slleval(ier,w,fs,sgain,ifsecant,coefsout,
     1      sgainout,apprrhs1,err,errnosu,srccoe,pattcoe,nlege)
c 
        call prinf('after slleval, ier=*',ier,1)
        call prin2('and remember, sgain=*',sgain,1)
c 
        sgain=200
        ifsecant=1
c 
        call slleval(ier,w,fs,sgain,ifsecant,coefsout,
     1      sgainout,apprrhs2,err,errnosu,srccoe,pattcoe,nlege)
c 
        call prinf('after slleval, ier=*',ier,1)
c 
        iw=25
        call lotagraph3(iw,xs,fs,n,xs,apprrhs1,n,
     1      xs, apprrhs2,n, 'patterns for different supergains*')
c 
        stop
        end
c 
c 
c 
c 
c 
        subroutine prointev(c,ts,whts,fs,x,n,rint)
        implicit real *8 (a-h,o-z)
        save
        dimension ts(1),whts(1)
        complex *16 rint,fs(1),ima
c 
        data ima/(0.0d0,1.0d0)/
c 
c        use the integral formula to evaluate the prolate function
c 
        rint=0
        do 1600 i=1,n
c 
        rint=rint+whts(i)*fs(i)*exp(ima*c*x*ts(i))
 1600 continue
c 
        return
        end
c 
c 
c 
c 
c 
      SUBROUTINE sgainleg(X,VAL,PEXP,N)
      IMPLICIT REAL *8 (A-H,O-Z)
        save
      complex *16 PEXP(1),val
C 
C     This subroutine computes the value o a Legendre
c     expansion with coefficients PEXP at point X in interval [-1,1]
C 
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
C 
        done=1
        pjm2=1
        pjm1=x
c 
        val=pexp(1)*pjm2+pexp(2)*pjm1
c 
        DO 600 J = 2,N
c 
        pj= ( (2*j-1)*x*pjm1-(j-1)*pjm2 ) / j
        val=val+pexp(j+1)*pj
c 
        pjm2=pjm1
        pjm1=pj
 600   CONTINUE
c 
        RETURN
        END
c 
c 
c 
c 
c 
        subroutine intedisc(ts,whts,n,a,b,xs,whtsx)
        implicit real *8 (a-h,o-z)
        save
        dimension ts(1),xs(1),whtsx(1),whts(1)
c 
        u=(b-a)/2
        v=(b+a)/2
        do 1200 i=1,n
        xs(i)=ts(i)*u+v
c 
        whtsx(i)=whts(i)*u
 1200 continue
c 
        return
        end
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c       This is the end of the debugging code and the beginning of the
c       band-limited least squares approximation code proper
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine slleval(ier,w,fs,sgain,ifsecant,coefsout,
     1      sgainout,apprrhs,err,errnosu,srccoe,pattcoe,nlege)
        implicit real *8 (a-h,o-z)
        save
        dimension fs(1),coefsout(1),apprrhs(1),w(1)
c 
        dimension srccoe(1),pattcoe(1)
c 
c       this subroutine approximates the user-supplied function (to be
c       referred to as the "pattern") on the interval [-1,1] by a
c       band-limited function, having the user-prescribed amount of
c       supergain. Depending on  the value of the input parameter
c       ifsecant (see below) the actual supergain of the constructed
c       function will either be (more or less) exactly equal to the
c       user-specified value, or be only a reasonable approximation to
c       the requested supergain (this way, the scheme is noticeably
c       simpler, and who needs the supergain exactly, anyway!). This
c       subroutine produces four principal pieces of output:
c 
c    1. The array alphas of the coefficients of a prolate expansion
c       (hopefully) approximating the user-supplied input function ff,
c       i.e.
c 
c     \phi(x) = \sum_{j=0}^n alphas(i) \cdot \psi_i^c (x)  \sim ff(x)    (1)
c 
c       on the interval [-1,1].
c 
c    2. The array srccoe of Legendre (nlege of them) coefficients of the
c       function \sigma : [-1,1] \to C^1, such that
c 
c       \phi(x) = \int_{-1}^1 e^{i \cdot c \cdot x \cdot t }
c                  \cdot \sigma(t) dt,                                   (2)
c 
c       with
c 
c       \sigma(x)=sum_0^{nlege} srccoe(j) \cdot P_j(x)                   (3)
c 
c       (had we been electrical engineers, we would say that \sigma
c       is the source distribution generating the pattern \phi that in
c       turn approximates the user-specified function ff).
c 
c    3. The array patt of Legendre (nlege of them) coefficients of the
c       function \phi : [-1,1] \to C^1, such that \phi is close to ff,
c       i.e.
c 
c       \phi(x)=sum_0^{nlege} pattcoe(j) \cdot P_j(x) \simm ff(x)        (4)
c 
c       on the interval [-1,1].
c 
c    4. The array pattexp of coefficients of an exponential expansion
c       of the function \phi on the interval [-1,1], i.e.
c 
c       \phi(x)=\sum_{j=1}^{nexpons} pattexp(j) \cdot
c               e^{i \cdot c \cdot ttexp(j) \cdot x},                    (5)
c 
c       where the points ttexp and their number nexpons have been
c       (hopefully) produced by a prior call to the subroutine
c       slepapin (see)
c 
c                 Input parameters:
c 
c  sgain - the amount of supergain the user would like the appreoximation
c       to possess
c  ifsecant - tells the subroutine whether to enforce the supergain
c       condition exactly, or to settle for an approximate result;
c     ifsecant=0 will cause the subroutine to settle for the approximate
c       result
c     ifsecant=1 will cause the subroutine to enforce the supergain
c       condition exactly
c  fs - the function (pattern) to be approximated
c  w - array of data of varous kinds produced by a prior call to
c       sllinit (see)
c 
c                  Output parameters:
c 
c  ier - error return code;
c         ier=0 means successful execution
c         ier=2 means that the subroutine failed to match the
c               user-specified supergain to its satisfaction. Normally,
c               it means that the supergain specified by the user is
c               too large (the user has been too humble); in any event,
c               the user should check the parameter sgainout (see below)
c  coefsout - the (real) coefficients of the prolate expansion of the
c       pattern approximating the user-supplied function fs
c  sgainout - the actual supergain of the approximating signal (normally,
c        it is pretty close to the user-supplied sgain)
c  apprrhs - the approximation to the rhs defined by the coefficients
c  err - the error of the obtained approximation. Please note that this
c       is the L^2 error in the norm with the weight provided by the
c       user via the array whts, in the input to the subroutine sllinit
c  errnosu - the discrepancy of the expansion of rhs into the basis
c       functions obtained without any supergain control; no
c       supergain-controlled approximation can have a smaller discrepancy
c       in array coefsout (i.e. the actual approximant obtained)
c  srccoe - the (complex) Legendre expansion on the interval [-1,1] of
c       the function \sigma such that
c 
c        fff=F(\sigma),
c 
c        where fff is the pattern approximating the user-supplied function
c        fs, and the operator F is defined by (1) in the description of
c        the subroutine slllinit (see). Please note that this array should
c        have nhigh elements reserved in it (nhigh is returned by the
c        subroutine sllinit)
c 
c  pattcoe - the Legendre expansion of the pattern approximating the
c        user-supplied function fs. Please note that this array should
c        have nhigh elements reserved in it (nhigh is returned by the
c        subroutine sllinit)
c  nlege - the order of each of the expansions srccoe, pattcoe (in other
c        words, each of the arrays pattcoe, srccoe contains nlege+1
c        elements
c 
c       . . . construct the least squares approximation to the
c             user-specified right-hand side fs
c 
        iamatr=w(21)
        irlams2=w(22)
        iwhts0=w(23)
c 
        n=w(24)
        m=w(25)
c 
        iww=w(26)
        irlams=w(27)
c 
        iwork=w(28)
  
        nhigh=w(4)
c 
        call slleaev(ier,w(iamatr),w(iww),fs,sgain,ifsecant,
     1       coefsout,sgainout,apprrhs,err,errnosu,w(iwhts0),
     2       srccoe)
c 
c        scale the coefficients of the prolate expansion to their
c        natural size, and use them to recalculate the approximation
c        to the RHS
c 
        d=0
        do 3400 i=1,m
c 
        d=d+coefsout(i)**2
        coefsout(i)=coefsout(i)*w(irlams2-1+i)
 3400 continue
c 
        call lesqsca(apprrhs,apprrhs,n,d2,w(iwhts0) )
c 
cccc        call prin2('and the supergain recomputed manually is*',
cccc     1      d/d2,1)
  
cccc        call prin2('and sgainout=*',sgainout,1)
cccc        call prin2('and supergain discrepancy is*',d/d2-sgainout,1)
cccc        call prin2('and d=*',d,1)
cccc        call prin2('and d2=*',d2,1)
c 
c       construct the Legendre coefficients of the obtained
c       pattern and of the source that generated it
c 
        call slllecoe(w,w(iwork),nhigh,coefsout,m,
     1      w(irlams),srccoe,pattcoe,nlege)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine sllinit(ier,c,xs,whtfun,whts,n,w,lenw,m,rlams,
     1      rlams2,keep,nhigh,lused)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),w(1),rlams(1),whts(1),whtfun(1)
c 
c         This subroutine prepares data to be used by the
c         subroutine slleval for the construction of band-limited
c         approximations to functions on the interval [-1,1].
c         Please note that this subroutine has no uses as a
c         stand-alone device. Also note that this subroutine
c         should only be used when the approximation has to be
c         with a weight function that is not a constant; for
c         approximation with constant weight function, the much
c         more efficient slepapin (see) should be used.
c 
c     EXPLANATION: The subroutines sllinit (this one) and slleval
c         (its companion) construct a band-limited approximation to
c         the user-supplied right-hand side with the band limit c.
c         The approximation is in the least squares sense, subject
c         to the user-supplied weight function, given by its values at
c         the points xs; these values are provided in the array whtfun.
c         Most importantly, the approximation produced by slleval will
c         have the user-specified supergain - this feature is the claim
c         to fame of these subroutines.
c 
c                   Input parameters:
c 
c  c - the band-limit of the function to be constructed
c  xs - the nodes on the interval [-1,1] where the functions
c       approximated will be tabulated
c  whtfun - the values of the weight function at the nodes xs
c  whts - the quadrature weights at the nodes in the array xs
c  n - the number of elements in each of the arrays xs, whts, whtfun
c  lenw - the amount of memory (in real *8 words) provided in the
c       array w; this parameter is used by the subroutine to decide
c       whether it has been provided enough memory; if the amount of
c       memory is insufficient, the subroutine bombs with an appropriate
c       error return code ier (see below). Generally, this parameter
c       should be reasonably large
c 
c 
c                  Output parameters:
c 
c  ier - error return code;                       TO BE FIXED UP
c         ier=0 means successful execution
c         2 < ier < 2048 means that the subroutine prolcrea has returned
c               a non-zero error code; ier is equal to that code.
c               Normally, the reason is insufficient memory provided by
c               the user in the array w. However, ier=1024 means some
c               sort of serious trouble in prolcrea; it has never been
c               encountered, and probably means a memory screw-up by
c               the user. This is a fatal error.
c         ier=2048 means that the user-provided memory is grossly
c               insufficient. This is a fatal error.
c  w - the first keep (see below) elements of the array w will be used
c       by the subroutine slepapev. Obviously, the frst keep elements of w
c       should not be changed between a call to this subroutine and the
c       subsequent calls to slepapev.
c  m - the number prolate functions used in the approximation
c  rlams - the first n (complex) eigenvalues of the Fourier integral operator
c 
c 
c            F(\sigma) (x) = \int_{-1}^1 e^{i \cdot x \cdot t}            (1)
c                            \cdot \sigma(t) dt
c 
c  rlams2 - the first n (real) eigenvalues of the Fourier integral operator
c 
c 
c            Q(\sigma) (x) = {1 \over \pi}  \int_{-1}^1
c                            { sin(x \cdot (x-t)) \over x-t}              (2)
c                            \cdot \sigma(t) dt
c 
c  keep - the number of elements in array w (see above) that should be
c        unchanged between the call to this subroutine and subsequent calls
c        to slepapev.
c  nhigh - the maximum order of the Legendre approximation of any of the
c        prolate functions encountered by the subroutine prolcrea when the
c        latter was called by this subroutine. It is provided to the user
c        because that is the number of elements that should be reserved in
c        each of the arrays srccoe, pattcoe when the subroutine slleval is
c        called (please note that the array srccoe is COMPLEX *16 while the
c        array pattcoe is REAL * 8.
c  lused - the total amount of space (in real *8 elements) used by the
c        subroutine in the array w
c 
c        . . . construct the prolate functions
c 
        ikhi=1
        lkhi=c+1000
c 
        iw=ikhi+lkhi
c 
        call prolcrea(jer,c,w(iw),lenw,nvects,nhigh,
     1    w(ikhi),rlams,rlams2,keep,lused)
c 
        call prinf('after prolcrea, jer=*',jer,1)
        call prinf('nhigh=*',nhigh,1)
        call prinf('keep=*',keep,1)
        call prinf('lused=*',lused,1)
        call prinf('nvects=*',nvects,1)
c 
        if(jer .eq. 0) goto 2100
c 
        ier=2048
        return
 2100 continue
c 
c       determine the number of prolate functions to be used
c       in the sdesign of the pattern
c 
        eps7=1.0d-29
        do 2200 i=1,nvects
        if(w(irlams22+i-1) .gt. eps7) m=i
 2200 continue
c 
c       compress the array w
c 
        call lesqarrm(w(iw),w,keep+2)
c 
        irlams=keep+6
        lrlams=nvects*2+6
c 
        irlams2=irlams+lrlams
        lrlams2=nvects+4
c 
        call lesqarrm(rlams,w(irlams),lrlams)
        call lesqarrm(rlams2,w(irlams2),lrlams2)
c 
c       using the obtained information, allocate memory in array
c       w for the various types of data to be used by slleval
c 
        iamatr=irlams2+lrlams2
        lamatr=n*m+10
c 
        iwhts=iamatr+lamatr
        lwhts=n+4
c 
        iwork=iwhts+lwhts
        lwork=nhigh*2+10
c 
        if(lwork .le. n+2) lwork=n+2
c 
        call lesqarrm(whts,w(iwhts),n+2)
c 
        iww=iwork+lwork
c 
        w(21)=iamatr+0.1
        w(22)=irlams2+0.1
        w(23)=iwhts+0.1
c 
        w(24)=n+0.1
        w(25)=m+0.1
c 
        w(26)=iww+0.1
        w(27)=irlams+0.1
        w(28)=iwork
c 
        lenww=lenw-iww
c 
        if(lenww .gt. 100) goto 2400
        ier=512
        return
 2400 continue
c 
c       construct the matrix of values of prolate functions at the
c       nodes of the discretization
c 
        do 2600 i=1,n
c 
        w(iwork-1+i)=whts(i)*whtfun(i)
 2600 continue
c 
        call sllinit0(ier,xs,w(iwork),n,m,w(iamatr),w,w(irlams2),
     1      w(iww),lenww,keep2,lused2)
c 
        keep=iww+keep2
        lused=iww+lused2
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine sllinit0(ier,xs,whts,n,m,amatr,w,rlams2,
     1      ww,lenww,keep,lused)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),w(1),amatr(n,m),whts(1)
c 
        real *8 rlams2(1)
c 
c       construct the matrix of values of prolate functions at the
c       nodes of the discretization
c 
        do 1600 i=1,m
        d=0
        do 1400 j=1,n
c 
        call prolevv(i-1,xs(j),w,val)
c 
        amatr(j,i)=val*rlams2(i)
 1400 continue
c 
 1600 continue
c 
c       prepare the data for the construction of the
c       least squares approximation
c 
        call lesqinit(ier,amatr,n,m,whts,ww,lenww,lused,keep)
c 
cccc        call prinf('after lesqinit, ier=*',ier,1)
cccc        call prinf('after lesqinit, keep=*',keep,1)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine slllecoe(w,coefs,nhigh,coefsout,m,
     1      rlams,srccoe,pattcoe,nlege)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),coefsout(1),coefs(1),pattcoe(1)
c 
        complex *16 srccoe(1),rlams(1)
c 
c       construct the Legendre coefficients of the obtained
c       pattern and of the source that generated it
c 
        call prinf('in slllecoe, nhigh=*',nhigh,1)
  
        do 4200 i=1,nhigh
        pattcoe(i)=0
        srccoe(i)=0
 4200 continue
c 
        half=0.5d0
        do 4600 i=1,m
c 
        call prolunpk(i-1,w,coefs,ndum)
c 
        do 4400 j=1,ndum
c 
        pattcoe(j)=pattcoe(j)+coefsout(i)*coefs(j)
c 
        srccoe(j)=srccoe(j)+coefsout(i)*coefs(j)/rlams(i)
c 
 4400 continue
 4600 continue
c 
c        if the need be, truncate the arrays pattcoe, srccoe
c 
        eps=1.0d-15
        nlege=0
c 
        do 4800 i=1,nhigh
c 
        if(abs(srccoe(i)) .gt. eps) nlege=i
 4800 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine slleaev(ier,a,w,rhs,sgain,ifsecant,
     1      coefsout,sgainout,apprrhs,err,errnosu,whts0,work2)
c 
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),rhs(1),coefsout(1),apprrhs(1),errssec(20),
     1      a(1),work2(1)
c 
c       for the user-specified function (given via the array rhs,
c       which is the table of values of the input function at the
c       user-supplied nodes), this subroutine constructs a least
c       squares approximation of rhs with prolate spheroidal wave
c       functions; the approximation is in the norm defined by the
c       user-supplied array whtfun. The least squares approximation
c       constructed by this subroutine is subject to the supergain
c       condition if the input parameter ifsecant has been set to 1,
c       and to the total power condition if ifsecant has been set to 0.
c 
c                     Input parameters:
c 
c  a - the matrix of values of prolate functions at user-specified
c       nodes, scaled by the corresponding elements of array rlams2;
c       (hopefully) obtained by a pripr call to sllinit (see)
c  w - array containing data constructed by the subroutine lesqinit
c       during the preceding call to sllinit
c  rhs - the function (pattern) to be approximated
c  whts - the weights (user-supplied) of the quadrature formula at the
c  sgain - the amount of supergain the user would like the appreoximation
c       to possess
c       user-supplied nodes
c  ifsecant - tells the subroutine whether to enforce the supergain
c       condition exactly, or to settle for an approximate result;
c     ifsecant=0 will cause the subroutine to settle for the approximate
c       result
c     ifsecant=1 will cause the subroutine to enforce the supergain
c       condition exactly
c  whts0 - the weights (user-supplied) of the quadrature formula at the
c       user-supplied nodes (referred to as whts elsewhere)
c 
c                     Output parameters:
c 
c  ier - error return code; please note that if ifsecant has been set
c       to 0, then the only possible values of ier on exit are
c       ier=0 and ier=32; all other values (see below) pertain to the
c       case when ifsecant=1
c 
c     ier=0 means successful execution
c     ier=2 means that the supergain condition has not been activated,
c       i.e. that simply projecting the rhs on the basis functions
c       (represented by the columns of the matrix u) results in a
c       value of the supergain that is less than sgain. Obviously,
c       this is a very mild condition: the user has been too modest.
c     ier=6 means that the secant method used by the subroutine to
c       satisfy the supergain condition failed to converge to 1.0e-6
c       in 50 secant iterations; please note that this condition
c       has never been encountered
c     ier=7, ier=17  mean that the secant process wandered outside
c       any reasonable range; probably, the supergain specified by
c       the user is too small (or, as usual, there might be a memory
c       problem). Please note that neither of these conditions has been
c       encountered by the author
c     ier=8 means that the error failed to decrese after one step of
c       secant algorithm; probably, the user-supplied supergain is too
c       small.
c     ier=32 means that the supergain specified by the user is so small
c       that setting the Lagrange multiplier to 10^10 (which
c       the subroutine lesqbis viewes as infinity) still results in
c       the supergain that is too large. The user has been way too
c       ambitious. As usual, this is a fatal condition
c 
c     ATTENTION: Please note that whenever the subroutine is forced
c       to generate a condition ier .ne. 0, it begins to behave as
c       if ifsecant had been set to 0. In other words, it returns to
c       the user an approximation to the function rhs whose supergain
c       is VERY APPROXIMATELY equal to the user-specified supergain
c  coefsout - the array of coefficients of the expansion into the basis
c       function into prolate functions (hopefully) possessing the
c       user-requested supergain
c  sgainout - the supergain of the obtained approximation
c  apprrhs - the approximation to the rhs defined by the coefficients
c       in array coefsout (i.e. the actual approximant obtained)
c  err - the error of the obtained approximation. Please note that this
c       is the L^2 error in the norm with the weight provided by the
c       user via the array whts, in the input to the subroutine sllinit
c  errnosu - the discrepancy of the expansion of rhs into the basis
c       functions obtained without any supergain control; no
c       supergain-controlled approximation can have a smaller discrepancy
c 
c 
c       . . . unpack the first 10 elements of array w to obtain
c             various forms of discrete information to be used
c             by the subroutine
c 
        ier=0
c 
        n=w(1)
        m=w(2)
        ncols=w(3)
c 
        iu=w(6)
        iwhts=w(7)
        icm=w(8)
        is=w(9)
        iv=w(10)
c 
c        construct the least squares approximation with the supergain
c        specified by the user
c 
        call sllea2ev(ier,w(iu),n,ncols,rhs,w(iwhts),w(icm),sgain,
     1      coefsout,sgainout,apprrhs,ifsecant,errnosu,
     2      errssec,niter,whts0,work2)
c 
        call prinf('in slleaev after sllea2ev, ier=*',ier,1)
        call prin2('in slleaev after errssec=*',erssec,niter)
c 
c       construct the linear combination of original functions
c       that is equal (approximately) to the rhs
c 
        do 2400 i=1,ncols
c 
        coefsout(i)=coefsout(i)/w(is+i-1)
 2400 continue
c 
        call lesqmatv(w(iv),m,ncols,coefsout,apprrhs)
        call lesqarrm(apprrhs,coefsout,m)
c 
c        recompute the rhs from the final coefficients and
c        compare it to the original rhs
c 
        call lesqmatv(a,n,m,coefsout,apprrhs)
c 
c        evaluate the difference between the recomputed rhs and the
c        original one
c 
        err2=0
        do 2600 i=1,n
        err2=(apprrhs(i)-rhs(i))**2*w(iwhts+i-1) +err2
 2600 continue
c 
        err=sqrt(err2)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine sllea2ev(ier,u,n,m,rhs,whts,cm,sgain,
     1      prods,sgainout,work,ifsecant,errnosu,
     2      errssec,niter,whts0,work2)
        implicit real *8 (a-h,o-z)
        save
        dimension u(n,m),rhs(1),whts(1),cm(1),prods(1),
     1      work(1),errssec(1),whts0(1),work2(1)
c 
c       for the user-specified function (given via the array rhs,
c       which is the table of values of the input function at the
c       user-supplied nodes), this subroutine constructs a least
c       squares approximation of rhs with the orthonormal functions
c       that are combinations of prolate functions; these have
c       (hopefully) been obtained by SVD-ing the prolate functions
c       in the norm defined by the user-supplied array whtfun, during
c       a preceding call to the subroutine sllinit). The least squares
c       approximation constructed by this subroutine is subject to the
c       supergain condition if the input parameter ifsecant has been
c       set to 1, and to the total power condition if ifsecant has been
c       set to 0.
c 
c                     Input parameters:
c 
c  u - the right matrix in the SVD of the matrix of values of prolate
c       functions, (hopefully) obtained by a pripr call to sllinit (see)
c  n - the number of points in the discretizations of all functions
c       (to be approximated, to be approximated with, and all others)
c  m - the number of approximating functions (also the number of elements
c       in each of the arrays prods, cm
c  rhs - the function (pattern) to be approximated
c  whts - the weights (user-supplied) of the quadrature formula at the
c       user-supplied nodes
c  cm - the array of penalty coefficients corresponding to the basis
c       functions
c  sgain - the amount of supergain the user would like the appreoximation
c       to possess
c  ifsecant - tells the subroutine whether to enforce the supergain
c       condition exactly, or to settle for an approximate result;
c     ifsecant=0 will cause the subroutine to settle for the approximate
c       result
c     ifsecant=1 will cause the subroutine to enforce the supergain
c       condition exactly
c  whts0 - the weights (user-supplied) of the quadrature formula at the
c       user-supplied nodes (referred to as whts elsewhere)
c 
c                     Output parameters:
c 
c  ier - error return code; please note that if ifsecant has been set
c       to 0, then the only possible values of ier on exit are
c       ier=0 and ier=32; all other values (see below) pertain to the
c       case when ifsecant=1
c 
c     ier=0 means successful execution
c     ier=2 means that the supergain condition has not been activated,
c       i.e. that simply projecting the rhs on the basis functions
c       (represented by the columns of the matrix u) results in a
c       value of the supergain that is less than sgain. Obviously,
c       this is a very mild condition: the user has been too modest.
c     ier=6 means that the secant method used by the subroutine to
c       satisfy the supergain condition failed to converge to 1.0e-6
c       in 50 secant iterations; please note that this condition
c       has never been encountered
c     ier=7, ier=17  mean that the secant process wandered outside
c       any reasonable range; probably, the supergain specified by
c       the user is too small (or, as usual, there might be a memory
c       problem). Please note that neither of these conditions has been
c       encountered by the author
c     ier=8 means that the error failed to decrese after one step of
c       secant algorithm; probably, the user-supplied supergain is too
c       small.
c     ier=32 means that the supergain specified by the user is so small
c       that setting the Lagrange multiplier to 10^10 (which
c       the subroutine lesqbis viewes as infinity) still results in
c       the supergain that is too large. The user has been way too
c       ambitious. As usual, this is a fatal condition
c 
c     ATTENTION: Please note that whenever the subroutine is forced
c       to generate a condition ier .ne. 0, it begins to behave as
c       if ifsecant had been set to 0. In other words, it returns to
c       the user an approximation to the function rhs whose supergain
c       is VERY APPROXIMATELY equal to the user-specified supergain
c  prods - the array of coefficients of the expansion into the basis
c       function (specified by the columns of u) possessing the
c       user-requested supergain
c  sgainout - the supergain of the obtained approximation
c  errnosu - the discrepancy of the expansion of rhs into the basis
c       functions obtained without any supergain control; no
c       supergain-controlled approximation can have a smaller discrepancy
c  errsec - the errors produced by the secant process during all
c       niter iterations
c  niter - the number of iterations performed by the secant process
c 
c                      work arrays:
c 
c  work, work2 - must be at least m real *8 elements each
c 
c        . . . calculate the inner products of all of the basis
c              elements with the right-hand side
c 
        ier=0
c 
        eps=1.0d-6
c 
        errnosu=0
        do 1200 i=1,m
c 
        call lesqsca(u(1,i),rhs,n,prods(i),whts)
c 
        work(i)=prods(i)
c 
 1200 continue
c 
c       evaluate the error of the approximation without supergain
c       control
c 
        do 1500 i=1,n
c 
        d=0
        do 1400 j=1,m
c 
        d=d+prods(j)*u(i,j)
 1400 continue
        errnosu=errnosu+(rhs(i)-d)**2*whts(i)
 1500 continue
c 
        errnosu=sqrt(errnosu)
c 
c       evaluate the norm of the rhs
c 
        call lesqsca(rhs,rhs,n,d,whts0)
c 
        totpowr1=d*sgain
c 
c       construct the approximation to the right-hand side
c       subject to the approximate supergain condition
c 
        call sllearl(jer,cm,prods,n,totpowr1,rlout,denom,
     1      sgainout1,m,whts0,u)
c 
        call lesqarrm(prods,work2,m)
c 
c        if the user so requested - do not enforce the supergain
c        condition exactly; use the approximate one obtained
c        on the first attempt
c 
         if(ifsecant .eq. 1) goto 1530
c 
         if(jer .eq. 16) ier=32
         sgainout=sgainout1
         return
 1530 continue
c 
c        if the supergain condition has been enforced automatically
c        or something has gone terribly wrong - get out
c 
         sgainout = sgainout1
  
         if(jer .eq. 0) goto 1550
c 
         if(jer .eq. 4) ier=2
         if(jer .eq. 16) ier=32
c 
 1550 continue
c 
c        the user has requested that supergain condition be enforced
c        exactly - prepare the initial data for the secant process
c 
        call lesqarrm(work,prods,m)
c 
        totpowr2=denom*sgain
        call sllearl(jer,cm,prods,n,totpowr2,rlout,denom,
     1      sgainout2,m,whts0,u)
c 
        if(jer .eq. 0) goto 1800
c 
        ier=7
        call lesqarrm(work2,prods,m)
        return
 1800 continue
c 
        dd=(sgainout2-sgain)/sgain
c 
        ddold=1.0d20
        do 2500 ijk=1,50
c 
        niter=ijk
c 
        alpha=(sgainout2-sgainout1)/(totpowr2-totpowr1)
        beta=sgainout1-alpha*totpowr1
        totpowr3=(sgain-beta)/alpha
c 
        if(totpowr3 .le. 0) totpowr3=totpowr2/2
        call lesqarrm(work,prods,m)
c 
        call sllearl(jer,cm,prods,n,totpowr3,rlout,denom,
     1       sgainout3,m,whts0,u)
c 
        if(jer .eq. 0) goto 2000
c 
        ier=17
        call lesqarrm(work2,prods,m)
        return
 2000 continue
c 
        if(abs(ddold/sgain) .lt. eps) goto 2550
c 
        ddold=dd
        dd=(sgainout3-sgain)/sgain
c 
        errssec(ijk)=ddold
c 
c       if the error failed to decrease on this step - stop iterations
c 
        if  (ijk .eq. 1) goto 2200
        if (abs(ddold) .gt. abs(dd)) goto 2200
c 
        ier=8
        call lesqarrm(work2,prods,m)
        return
 2200 continue
c 
        totpowr1=totpowr2
        sgainout1=sgainout2
c 
        totpowr2=totpowr3
        sgainout2=sgainout3
c 
 2500 continue
c 
c       precision eps has not been achieved after 12 iterations.
c       set error return code to 6 and reset the coefficients
c       of the expansions to what they were before we began to iterate
c 
        ier=6
        call lesqarrm(work2,prods,m)
 2550 continue
c 
        errssec(niter)=dd/sgain
        sgainout=sgainout3
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine sllearl(ier,cm,prods,n,totpowr,rlout,denom,
     1      sgainout,m,whts,u)
        implicit real *8 (a-h,o-z)
        save
        dimension cm(1),prods(1),whts(1),u(n,m)
c 
c       for the weighted products prods between the basis functions
c       (obtained by SVD-ing the prolate functions in the norm defined
c       by the user-supplied array whtfun), this subroutibe constructs
c       the coefficients of the least squares approximation possessing
c       the user-specified total power totpowr. It also evaluates the
c       actual supergain of the obtained expansion
c 
c                     Input parameters:
c 
c  cm - the array of penalty coefficients corresponding to the basis
c       functions
c  prods - the array of inner products between the function to be
c       approximated and the functions into which it is to be expanded
c       (in the norm defined by the user-supplied weight function
c       whtfun)
c  n - the number of points in the discretizations of all functions
c       (to be approximated, to be approximated with, and all others)
c  totpowr - the total power of the expansion with which the approximation
c       is to be constructed
c  m - the number of approximating functions (also the number of elements
c       in each of the arrays prods, cm
c  whts - the weights (user-supplied) of the quadrature formula at the
c       user-supplied nodes
c  u - the right matrix in the SVD of the matrix of values of prolate
c       functions, (hopefully) obtained by a pripr call to sllinit (see)
c 
c                     Output parameters:
c 
c  ier - error return code
c 
c     ier=0 means successful execution
c     ier=4 means that the supergain condition has not been activated.
c       In other words, the standard least squares approximation to the
c       right-hand side satisfies the supergain condition automatically.
c       Obviously, this is a very mild condition: the user has been
c       too modest.
c     ier=16 means that setting the Lagrange multiplier to 10^10 (which
c       the subroutine lesqbis viewes as infinity) still results in
c       the supergain that is too large. The user has been too ambitious.
c       As ususal, this is a fatal condition
c  rlout - the value of the Lagrange multiplier that produces the
c       supergain totpowr; if rlout=10^10 still produces the supergain
c       that is too large, rlout is set to 10^10
c  denom - the L^2 norm of the obtained approximation
c  sgainout - the supergain of the obtained approximation
c 
c 
c       . . . find the value of the Lagrange multiplier
c 
        call lesqbis(ier,cm,prods,m,totpowr,rlout)
c 
        if(ier .eq. 16) rlout=1.0d10
        if(ier .eq. 4) rlout=0
        if(ier .ne. 0) return
c 
c       construct the projections of the resulting least squares
c       approximation on the orthogonal functions
c 
        sumtes=0
        do 2400 i=1,m
c 
        prods(i)=prods(i)/(1+rlout*cm(i))
c 
        sumtes=sumtes+cm(i)*prods(i)**2
 2400 continue
c 
c       evaluate the denominator in the calculation of the supergain
c 
        denom=0
        do 3600 i=1,n
c 
        d=0
        do 3200 j=1,m
        d=d+u(i,j)*prods(j)
 3200 continue
c 
        denom=denom+whts(i)*d**2
 3600 continue
c 
        sgainout=sumtes/denom
c 
        return
        end
  
