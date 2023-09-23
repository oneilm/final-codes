        implicit real *8 (a-h,o-z)
        real *8 w(8 000 000),rlams(10000),
     1     rlams2(10000),tt(10 000),tt0(10 000),ww(10 000),
     2      ww0(10 000),
     3      pattrea(10 000),pattima(10 000),
     4      srccoe(10 000),pattcoe(100  000),
     5      sourcerea(10 000),sourceima(10 000),ffima(10 000)
c 
        complex *16 source(10 000),cd,ima,
     1      ff(10 000),patt(10 000),alphas(10 000),ff2(10 000),
     2      diffs(10 000),derout,ff1(10 000),dersout(10 000),
     3      dersout2(10 000),pattexp(10 000)
c 
        dimension tt1(10 000),
     1      ff1rea(10 000),ff1ima(10 000),ww1(10 000),ttexp(10 000)
c 
        data ima/(0.0d0,1.0d0)/
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER c'
        READ *,c
        CALL PRIN2('c=*',c,1 )
C 
        PRINT *, 'ENTER b'
        READ *,b
        CALL PRIN2('b=*',b,1 )
C 
        PRINT *, 'ENTER supergain'
        READ *,sgain
        CALL PRIN2('sgain=*',sgain,1 )
C 
        oversamp=1.5
        nplot=500
  
        nn=600
cccc        nn=2000
  
c 
c       discretize the interval [-b,b]
c 
        itype=1
        call legeexps(itype,nn,tt,u,v,ww)
c 
        do 1200 i=1,nn
c 
        tt0(i)=tt(i)
        ww0(i)=ww(i)
  
        tt(i)=tt(i)*b
  
        ww(i)=ww(i)*b
c 
        ff(i)=1 +ima*tt(i)**2
        ff(i)=1  
  
cccc        if(i .eq. 20) ff(i)=ff(i)+0.01
  
cccc        ff(i)=1 +ima*sin(tt(i)*10)
        ffima(i)=ff(i)/ima
 1200 continue
  
  
  
c 
cccc        call prin2('tt as created*',tt,nn)
cccc        call prin2('ww as created*',ww,nn)
cccc        call prin2('ff as created*',ff,nn*2)
c 
c       prepare the data to be used in the construction of the optimal
c       beam corresponding to the user-specified geometry and supergain
c 
        lenw=8 000 000
cccc        lenw=40 000
c 
        call slepapin(ier,c,tt,ww,nn,w,lenw,rlams,
     1      rlams2,nvects,nhigh,oversamp,ttexp,nexpons,keep)
  
  
        call prinf('after slepapin, nhigh=*',nhigh,1)
        call prinf('after slepapin, keep=*',keep,1)
        call prinf('after slepapin, ier=*',ier,1)
  
         if(ier .gt. 2) stop
  
c 
c       construct the optimal beam corresponding to the user-specified
c       geometry and total supergain
c 
        call slepapev(ier,sgain,ff,w,
     1      alphas,n,srccoe,pattcoe,nlege,pattexp,errl2,sgainout)
  
        call prinf('after slepapev, ier=*',ier,1)
        call prinf('after slepapev, nlege=*',nlege,1)
  
        call prin2('l2-error of approximation of original pattern*',
     1      errl2,1)
  
        call prin2('actual supergain of obtained pattern*',sgainout,1)
c 
c       construct the values of the obtained pattern
c       on the interval [-1,1] and plot it
c 
        itype=1
        call legeexps(itype,nplot,tt0,u,v,ww0)
c 
        do 2900 i=1,nplot
c 
        call sgainleg(tt0(i),patt(i),pattcoe,nlege)
  
        pattrea(i)=patt(i)
        pattima(i)=-patt(i)*ima
c 
 2900 continue
c 
        iw=21
        call quagraph(iw,tt0,pattrea,nplot,3,
     1      'real part of pattern *')
c 
        iw=22
        call quagraph(iw,tt0,pattima,nplot,3,
     1   'imaginary part of pattern *')
c 
        iw=27
        call quagraph2(iw,tt0,pattima,nplot,3,tt,ffima,nn,3,
     1   'imaginary part of pattern vs original one*')
c 
c       construct and plot the source distribution that
c       generates the obtained pattern
c 
        do 3600 i=1,nplot
c 
        call sgainleg(tt0(i),cd,srccoe,Nlege)
c 
        source(i)=cd
        sourcerea(i)=cd
        sourceima(i)=-cd*ima
 3600 continue
c 
c       plot the real and imaginary parts of the obtained source
c       distribution
c 
        iw=23
        call quagraph(iw,tt0,sourcerea,nplot,3,
     1      'real part of source generating the pattern*')
c 
        iw=24
        call quagraph(iw,tt0,sourceima,nplot,3,
     1      'imaginary part of source generating the pattern*')
c 
c       evaluate the obtained exponential expansion at the nodes of
c       the discretization and compare it to the original function ff
c 
        diffmax=0
        do 5200 i=1,nn
c 
        call slepextr(c,ttexp,pattexp,nexpons,tt(i),ff2(i),derout )
c 
        diffs(i)=ff2(i)-ff(i)
  
        d=abs(diffs(i))
  
        if(d .gt. diffmax) diffmax=d
 5200 continue
  
        call prin2('maximum error of exponential expansion is*',
     1       diffmax,1)
c 
c       discretize the interval [-2,2] and plot the real and imaginary parts
c 
        nplot=500
        h=1
        h=2*h/(nplot-1)
c 
        do 5400 i=1,nplot
c 
        tt1(i)=-1+(i-1)*h
        ww1(i)=h
c 
        tt1(i)=tt1(i)*3
        ww1(i)=ww1(i)*3
c 
        call slepextr(c,ttexp,pattexp,nexpons,tt1(i),ff1(i),dersout(i))
c 
        ff1rea(i)=ff1(i)
        ff1ima(i)=-ff1(i)*ima
 5400 continue
c 
c 
        iw=25
        call quagraph(iw,tt1,ff1rea,nplot,3,
     1      'extended real part*')
c 
c 
        iw=26
        call quagraph(iw,tt1,ff1ima,nplot,3,
     1      'extended imaginary part*')
c 
c       test the derivative evaluated via exponential expansion via
c       the derivative evaluated numerically
c 
        dererr=0
        do 5600 i=2,nplot-1
c 
        dersout2(i)=(ff1(i+1)-ff1(i-1))/(tt1(i+1)-tt1(i-1))
c 
        d=abs(dersout2(i)-dersout(i))
        if(dererr .lt. d) dererr=d
 5600 continue
c 
cccc        call prin2('and dersout2=*',dersout2,nplot*2)
cccc        call prin2('while dersout=*',dersout,nplot*2)
        call prin2('and maximum derivative error is*',dererr,1)
        stop
        end
  
c 
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c         This is the end of the debugging code and the beginning of the
c         actual band-limited approximatipn routines
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        subroutine slepextr(c,ttt,vals,nnn,x,fout,derout)
        implicit real *8 (a-h,o-z)
        save
        dimension ttt(1)
        complex *16 vals(1),cd,ima,fout,derout,xcima
c 
        data ima/(0.0d0,1.0d0)/
c 
c       Given a collection of user-provided point sources of
c       (complex) intensities vals(j) at points ttt(j) on R^1,
c       this subroutine evaluates the far field generated by
c       these sources at the point x \in R^1, via the
c       formula
c 
c       fout = sum_{j=1}^nnn
c              vals(i) \cdot e^{i \cdot c \cdot ttt(j) \cdot x};      (1)
c 
c       the subroutine also returns the derivative derout of fout
c       with respect to x.
c 
c                       Input parameters:
c 
c  c - the parameter in (1) above
c  ttt - the locations of soyurces on the line
c  vals - the intensities of the sources
c  nnn - the number of sources
c  x - the point where the far field is to be evaluated
c 
c                       Output parameters:
c 
c  fout - the far field generated at the point x by the sources of
c        intensities vals located at the points ttt
c  derout  - the derivative of fout with respect to x
c 
c       . . . evaluate the function and the derivative at the user-specified
c             point x (x can be anywhere on the interval [-2,2])
c 
        xcima=ima*c*x
        fout=0
        derout=0
c 
        do 1200 i=1,nnn
c 
        cd=exp(xcima*ttt(i))*vals(i)
c 
        fout=fout+cd
        derout=derout+cd*ttt(i)
 1200 continue
c 
        derout=derout*ima*c
        return
        end
c 
c 
c 
c 
c 
c 
c 
c 
c 
c 
        subroutine slepapin(ier,c,tt,ww,nn,w,lenw,rlams,
     1      rlams2,nvects,nhigh,oversamp,ttexp,nexpons,keep)
        implicit real *8 (a-h,o-z)
        save
        dimension tt(1),w(1),rlams2(1),ww(1),rlams(1),ttexp(1)
c 
c       this subroutine constructs various types of data to be used
c       by the subroutine slepapev (see) to approximates user-supplied
c       functions (to be referred to as the "patterns") on the interval
c       [-1,1] by band-limited functions, having user-prescribed amounts
c       of supergain. Please note that the actual supergains of the
c       constructed functions are only reasonable approximations to the
c       requested supergains (this way, the scheme is noticeably simpler,
c       and who needs the supergain exactly, anyway!). This subroutine
c       has no uses as a stand-alone device.
c 
c         NOTE:
c 
c     Normally, the output array w produced by this subroutine is to be
c       used by the subroutine slepapev (see). However, at a certain
c       point, this subroutine calls the subroutine prolcrea (see); the
c       output from prolcrea is stored in array w starting with position
c       21. Thus, the chunk of w starting at position 21 can be fed into
c       the subroutines proleva, prolunpk, to be used for the evaluation
c       of the prolate functions, etc.
c 
c 
c                 Input parameters:
c 
c  c - the band-limit of the function to be constructed
c  tt - the nodes on the interval [-1,1] where the function ff to be
c       approximated is tabulated
c  ww - the quadrature weights at the nodes in the array tt
c  nn - the number of elements in each of the arrays tt, ff, ww
c  lenw - the amount of memory (in real *8 words) provided in the
c       array w; this parameter is used by the subroutine to decide
c       whether it has been provided enough memory; if the amount of
c       memory is insufficient, the subroutine bombs with an appropriate
c       error return code ier (see below). Generally, this parameter
c       should not be less than 40 000
c  oversamp - the coefficient telling the subroutine how many nodes on
c       the interval [-1,1] to use in the construction of exponential
c       expansions of user-supplied functions (see subroutines slepapev,
c       slepextr).
c 
c     The value oversamp=1 will guarantee that the expansion is valid on
c       the interval [-1,1];
c     The value oversamp=1.5 will guarantee that the expansion is valid on
c       the interval [-2,2];
c     The value oversamp=2 will guarantee that the expansion is valid on
c       the interval [-3,3], and so on.
c 
c     On the other hand, who cares if the expansion is "valid" outside
c       the interval [-1,1] ? It is correct on the interval [-1,1], and
c       appropriately smooth and bounded elsewhere - and who cares about
c       the rest!
c 
c                  Output parameters:
c 
c  ier - error return code;
c         ier=0 means successful execution
c         2 < ier < 2048 means that the subroutine prolcrea has returned
c               a non-zero error code; ier is equal to that code.
c               Normally, it reason is insufficient memory provided by
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
c  nvects - the number of elements in each of the arrays rlams, rlams2
c  nhigh - the highest order of any of the Legendre expansions of the
c         nvects prolate functions constructed by prolcrea when the latter
c         was called by this subroutine.
c 
c  ttexp - nexpons the Gaussian nodes on the interval [-1,1], to be used
c         as a part of input into the subroutine slepextr (see)
c  nexpons - the number of elements in the array ttexp, to be used
c         as a part of input into the subroutine slepextr (see)
c  keep - the number of elements in array w (see below) that should be
c        unchanged between the call to this subroutine and subsequent calls
c        to slepapev.
c 
c        allocate memory for the construction of prolate functions
c 
         ikhi=1
         lkhi=c*2+20 000
c 
         iw=ikhi+lkhi
         lenw2=lenw-iw
c 
c        construct all the necessary prolate functions
c 
        ier=0
c 
        call prolcrea(jer,c,w(iw),lenw2,nvects,nhigh,
     1    w(ikhi),rlams,rlams2,keep,lused)
c 
        if(jer .eq. 0) goto 1100
c 
        call prinf(
     1  'bombing from sgainpow, since after prolcrea, jer=*',
     2      jer,1)
        ier=16
        if(jer .eq. 1024) ier=1024
        return
 1100 continue
c 
        keep=lused
        n=nvects-2
c 
        if(2 .ne. 3) goto 1150
c 
        call prin2('after prolcrea in sgainpow w(ikhi) are*',
     1      w(ikhi),nvects)
        call prin2('after prolcrea, in sgainpow rlams2=*',
     1      rlams2,nvects)
        call prin2('after prolcrea, in sgainpow rlams=*',
     1      rlams,nvects*2)
c 
        call prinf('ier=*',ier,1)
        call prinf('nhigh=*',nhigh,1)
        call prinf('keep=*',keep,1)
        call prinf('lused=*',lused,1)
        call prinf('nvects=*',nvects,1)
c 
 1150 continue
c 
c       compress the work array, eliminating the space for khi
c 
        iw2=21
        do 1200 i=1,keep
        w(iw2+i-1)=w(iw+i-1)
 1200 continue
c 
        iw=21
        lw=keep
c 
c       allocate memory for the subroutine sgainpow
c 
        ier=0
c 
        icm=iw+lw
        lcm=n+2
c 
        istore=icm+lcm
        lstore=nn*n+20
c 
        ifuns=istore+lstore
        lfuns=n+2
c 
        ifders=ifuns+lfuns
        lfders=n+2
c 
        iwwww=ifders+lfders
        lwwww=lenw-iwwww
c 
c       construct the tables to be used in the design of approximations
c       with prescribed supergain
c 
        call sgainpow(jer,tt,nn,w(iw),rlams2,n,w(icm),w(istore),
     1      w(ifuns),w(ifders),w(iwwww),lwwww)
c 
        if(jer .ne. 0) then
            ier=8
            return
        endif
c 
c       store in the array w the remainder of the information to be
c       used in the design of approximations
c       with prescribed supergain
c 
        irlams=istore+lstore
        lrlams=n*2+4
c 
        irlams2=irlams+lrlams
        lrlams2=n+4
c 
        itt=irlams2+lrlams2
        ltt=nn+4
c 
        iww=itt+ltt
        lww=nn+4
c 
c        construct and store in array w the Gaussian nodes and weights
c        to be used in the design of the exponential expansions
c        of approximations with prescribed supergain
c 
cccc        nexpons=nhigh+20
        nexpons=nhigh*oversamp+20
c 
        ittexp=iww+lww
        lttexp=nexpons+2
c 
        iwwexp=ittexp+lttexp
        lwwexp=nexpons+2
c 
        itype=1
        call legeexps(itype,nexpons,w(ittexp),u,v,w(iwwexp))
c 
        do 1450 i=1,nexpons
c 
        ttexp(i)=w(ittexp+i-1)
 1450 continue
c 
        icoefs=iwwexp+lwwexp
        lcoefs=nhigh
c 
        keep=icoefs+lcoefs
c 
        if(keep .lt. lenw) goto 1500
        ier=4
        return
 1500 continue
c 
        do 1600 i=1,n
c 
        w(irlams2+i-1)=rlams2(i)
        w(irlams+i-1)=rlams(i)
        w(irlams+i+n-1)=rlams(i+n)
 1600 continue
c 
        do 1800 i=1,nn
c 
        w(itt+i-1)=tt(i)
        w(iww+i-1)=ww(i)
 1800 continue
c 
        w(1)=iw+0.1
        w(2)=icm+0.1
        w(3)=istore+0.1
        w(4)=irlams+0.1
        w(5)=irlams2+0.1
c 
        w(6)=n+0.1
        w(7)=nhigh+0.1
        w(8)=nn+0.2
c 
        w(9)=itt+0.1
        w(10)=iww+0.1
c 
        w(11)=c
c 
cccc        w(12)=nlege+0.1
        w(13)=icoefs+0.1
c 
        w(14)=nexpons+0.1
        w(15)=ittexp+0.1
        w(16)=iwwexp+0.1
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine sgainpow(ier,tt,nn,w,rlams2,n,cm,www,funs,ders,
     1      wwww,lww)
        implicit real *8 (a-h,o-z)
        save
        dimension tt(1),w(1),rlams2(1),cm(1),www(1),
     1      wwww(1),funs(1),ders(1)
c 
c       calculate the supergains corresponding to the prolate functions
c 
        ier=0
c 
        do 1200 i=1,n
c 
        cm(i)=0
        if(rlams2(i) .lt. 1.0d-36) goto 1200
c 
        cm(i)=1/rlams2(i)
 1200 continue
c 
c       construct the values of the first n prolate functions
c       at the user-specified nodes
c 
        ii=w(4)
c 
        nhigh=w(4)
        nvects=n
        lww=100 000
  
        ifunpack=1
        do 2800 j=1,nn
c 
        call prolderfast(jer,tt(j),w,ifunpack,nhigh,
     1      funs,ders,nvects,wwww,lww,keep)
c 
        if(jer .ne. 0) then
            ier=64
            return
        endif
c 
        ifunpack=0
c 
        do 2600 i=1,n
c 
        www(nn*(i-1)+j)=funs(i)
 2600 continue
  
 2800 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine slepapev(ier,sgain,ff,w,
     1      alphas,n,srccoe,pattcoe,nlege,pattexp,errl2,sgainout)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1)
c 
        complex *16 pattcoe(1),srccoe(1),ff(1),alphas(1),pattexp(1),cd
c 
c       this subroutine approximates the user-supplied function (to be
c       referred to as the "pattern") on the interval [-1,1] by a
c       band-limited function, having the user-prescribed amount of
c       supergain. Please note that the actual supergain of the
c       constructed function is only a reasonable approximation to the
c       requested supergain (this way, the scheme is noticeably simpler,
c       and who needs the supergain exactly, anyway!). This subroutine
c       produces four principal pieces of output:
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
c  ff - the function (pattern) to be approximated
c  w - array of data of varous kinds produced by a prior call to
c       slepapin (see)
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
c  alphas - the (complex) coefficients of the prolate expansion of the
c       pattern approximating the user-supplied function ff
c  n - the number of elements in each of the array alphas
c  srccoe - the Legendre expansion on the interval [-1,1] of the function
c       \sigma such that
c 
c        fff=F(\sigma),
c 
c        where fff is the pattern approximating the user-supplied function
c        ff, and the operator F is defined by (1) in the description of
c        the subroutine slepapin (see)
c  pattcoe - the Legendre expansion of the pattern approximating the
c        user-supplied function ff
c  nlege - the order of each of the expansions srccoe, pattcoe (in other
c        words, each of the arrays pattcoe, srccoe contains nlege+1
c        elements
c  errl2 - the error (in L^2 sense) with which the pattern fff approximates
c        the user-supplied function ff; please note that the error is
c        estimated on the user-provided nodes tt with user-provided weights
c        ww - caviat user!
c  sgainout - the actual supergain of the approximating signal (normally,
c        it is pretty clslose to the user-supplied sgain)
c 
c 
c       . . . unpack the beginning of the array w to obtain the needed
c             integer data
c 
        iw=w(1)
        icm=w(2)
        istore=w(3)
        irlams=w(4)
        irlams2=w(5)
c 
        n=w(6)
        nhigh=w(7)
        nn=w(8)
c 
        itt=w(9)
        iww=w(10)
c 
        c=w(11)
cccc        nlege=w(12)
        icoefs=w(13)
c 
c       construct the required approximation
c 
        call sgainpow2(ier,sgain,w(itt),ff,w(iww),nn,w(iw),
     1      w(irlams),alphas,n,srccoe,pattcoe,nlege,
     2      errl2,w(icm),w(icoefs),sgainout,w(istore),nhigh)
c 
c        evaluate the obtained source distribution at the nodes
c        of the appropriate Gaussian discretization
c 
  
        call prin2('after sgainpow2, alphas=*',alphas,n*2)
        call prin2('after sgainpow2, srccoe=*',srccoe,nlege*2)
        call prin2('after sgainpow2, pattcoe=*',pattcoe,nlege*2)
  
  
  
        ittexp=w(15)
        iwwexp=w(16)
        nexpons=w(14)
c 
        do 2200 i=1,nexpons
c 
c 
        call sgainleg(w(ittexp+i-1),cd,srccoe,nlege)
        pattexp(i)=cd*w(iwwexp+i-1)
 2200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine sgainpow2(ier,sgain,tt,ff,ww,nn,w,
     1      rlams,alphas,n,srccoe,pattcoe,nlege,errl2,
     2      cm,coefs,sgainout,www,nhigh)
        implicit real *8 (a-h,o-z)
        save
        dimension tt(1),w(1),cm(1),ww(1)
c 
        complex *16 rlams(1),pattcoe(1),
     1      srccoe(1),ff(1),cd,alphas(1)
c 
        real *8 coefs(1),www(nn,1)
c 
c       find the L^2 - norm on the user supplied right-hand
c       side, in order to evaluate the total (visible+invisible)
c       energy of the approximation to be constructed
c 
        ier=0
        eps=1.0d-15
        d=0
        do 1600 i=1,nn
c 
        d=d+ff(i)*conjg(ff(i))*ww(i)
 1600 continue
c 
        totpowr=d*sgain
c 
c       construct the projections of the function ff on all
c       the relevant prolate functions
c 
        do 1800 i=0,n-1
c 
        call sgainsca(ff,www(1,i+1),nn,cd,ww)
c 
        alphas(i+1)=cd
c 
 1800 continue
  
cccc        if(2 .ne. 3) return
  
  
c 
c       find the Lagrange multiplier corresponding to the
c       user-supplied total power
c 
        call sgainbis(jer,cm,alphas,n,totpowr,rlout)
  
        call prinf('after sgainbis, jer=*',jer,1)
c 
        if(jer .ne. 0) ier=2
c 
c       evaluate all the prolate coefficients
c 
        sumtes=0
        denom=0
        do 2400 i=1,n
c 
        alphas(i)=alphas(i)/(1+rlout*cm(i))
  
        if(cm(i) .lt. 0.1) alphas(i)=0
c 
        sumtes=sumtes+cm(i)*abs(alphas(i))**2
c 
        denom=denom+abs(alphas(i))**2
 2400 continue
c 
        sgainout=sumtes/denom
c 
        call prin2('in sgainpow2, cm=*',cm,n+2)
c 
c       construct the Legendre coefficients of the obtained
c       pattern
c 
        do 2450 i=1,nhigh
        pattcoe(i)=0
        srccoe(i)=0
 2450 continue
c 
        do 2550 i=1,n
c 
        if(cm(i) .lt. 0.1) goto 2550
c 
        call prolunpk(i-1,w,coefs,ndum)
c 
        do 2500 j=1,ndum
c 
        pattcoe(j)=pattcoe(j)+alphas(i)*coefs(j)
c 
        srccoe(j)=srccoe(j)+alphas(i)*coefs(j)/rlams(i)
c 
 2500 continue
 2550 continue
c 
        call prin2('in sgainpow2, srccoe=*',srccoe,nhigh)
        call prin2('in sgainpow2, pattcoe=*',pattcoe,nhigh)
c 
c        if the need be, truncate the arrays pattcoe, srccoe
c 
        nlege=0
c 
        do 3200 i=1,nhigh
c 
        if(abs(srccoe(i)) .gt. eps) nlege=i
 3200 continue
c 
c       evaluate the L^2-error of the obtained approximation
c       to the user-supplied pattern
c 
        d=0
        do 3400 i=1,nn
c 
        call sgainleg(tt(i),cd,pattcoe,Nlege)
c 
        cd=cd-ff(i)
        d=d+cd*conjg(cd)*ww(i)
c 
 3400 continue
c 
        errl2=sqrt(d)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine sgainbis(ier,cm,prods,n,totpowr,rlout)
        implicit real *8 (a-h,o-z)
        save
        dimension cm(1),prods(1)
c 
c       initialize the bisection process
c 
        ier=0
        rl1=0
        rl2=1.0d10
        rlout=0
c 
        call sgainsum(cm,prods,n,rl1,f1)
        call sgainsum(cm,prods,n,rl2,f2)
c 
        call prin2('in sgainbis, f1=*',f1,1)
        call prin2('in sgainbis, f2=*',f2,1)
c 
        if( (f1 .gt. totpowr) .and. (f2 .lt. totpowr) ) goto 1200
c 
        if(f1 .lt. totpowr) ier=4
        if(f2 .gt. totpowr) ier=16
        return
c 
 1200 continue
c 
c        reduce rl2 by decimal scaling as far as it will go
c 
        do 1300 i=1,100
c 
        rl22=rl2/10
  
        call sgainsum(cm,prods,n,rl22,f22)
c 
        if(f22 .ge. totpowr) goto 1350
c 
        rl2=rl22
        f2=f22
 1300 continue
 1350 continue
c 
        do 2000 i=1,15
c 
        rl3=(rl1+rl2)/2
c 
        call sgainsum(cm,prods,n,rl3,f3)
c 
        if(f3 .lt. totpowr) goto 1400
c 
        f1=f3
        rl1=rl3
        goto 2000
 1400 continue
c 
        f2=f3
        rl2=rl3
c 
 2000 continue
c 
        rlout=(rl1+rl2)/2
c 
c       use Newton to tighten the obtained value of rlout
c 
        do 1600 i=1,5
c 
        call sgainsum(cm,prods,n,rlout,sum)
        call sgainder(cm,prods,n,rlout,der)
c 
        ff=sum-totpowr
c 
        call prin2('in newton, ff=*',ff,1)
c 
        rlout=rlout-ff/der
 1600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine sgainder(cm,prods,n,rlagr,der)
        implicit real *8 (a-h,o-z)
        save
        dimension cm(1)
        complex *16 prods(1)
c 
c       evaluate the sum that is supposed to be equal to the
c       total power
c 
        der=0
        do 1200 i=1,n
c 
        der=der + cm(i)**2*abs(prods(i))**2/(1+rlagr*cm(i))**3
 1200 continue
c 
        der=-der*2
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine sgainsum(cm,prods,n,rlagr,sum)
        implicit real *8 (a-h,o-z)
        save
        dimension cm(1)
        complex *16 prods(1)
c 
c       evaluate the sum that is supposed to be equal to the
c       total power
c 
        sum=0
        do 1200 i=1,n
c 
        sum=sum + cm(i)*abs(prods(i))**2/(1+rlagr*cm(i))**2
 1200 continue
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine sgainsca(cx,y,n,prod,ww)
        implicit real *8 (a-h,o-z)
        save
        dimension y(1),ww(1)
        complex *16 cx(1),prod
c 
        prod=0
        do 1800 i=1,n
        prod=prod+cx(i)*(y(i)*ww(i))
 1800 continue
  
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
