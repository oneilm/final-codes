        implicit real *8 (a-h,o-z)
        real *8 corrs(1000 000)
        dimension xs(1000),whts(1000),w(1000 000),
     1       ainterp(1000 000),xsaux(4000),whtsaux(4000)
  
  
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER ndigits'
        READ *,ndigits
        CALL PRINf('ndigits=*',ndigits,1 )
c 
        PRINT *, 'ENTER n'
        READ *,n
        CALL PRINf('n=*',n,1 )
  
        naux=n*6
  
        iw=17
  
        call unidiscr(ndigits,n,naux,iw,
     1      xs,whts,xsaux,whtsaux,corrs,ainterp,w)
  
  
        stop
        end
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code, and the beginning
c        of the actual "universal" discretization cpde
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine unidiscr(ndigits,n,naux,iw,
     1      xs,whts,xsaux,whtsaux,corrs,ainterp,w)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension whtsaux(1),xsaux(1),ainterp(1),w(1),
     1      whts2aux(1000),corrs(naux,n),whts(1),xs(1)
c 
c        This subroutine produces a universal discretization of the
c        interval [-1,1], together with certain types of data to be
c        used in the construction of discretization of boundary
c        integral equations in situations where the boundary to be
c        discretized has corners. Please note that the principal
c        result of the operation of this subroutine is the output on
c        the FORTRAN file is (the number iw is provided by the user).
c        It is expected that these data will be read from the file iw
c        by the subroutine uniretr (see) or unipnts (see); all of
c        these data are also returned to the user in arrays corrs,
c        ainterp, xs, whts, xsaux, whtsaux (see below)
c 
c                            Input parameters:
c 
c  ndigits - the number of digits of accuracy the subroutine will attempt
c        to produce. PLEASE NOTE THAT THIS PARAMETER IS NOT USED BY THE
c        SUBROUTINE CORRECTLY. THE ACCURACY PRODUCED GROSSLY EXCEEDS THE
c        REQUESTED ACCURACY.
c  n - the number of nodes in the principal discretization of the
c        interval [-1,1]
c  naux - the number of nodes in the auxiliary discretization of the
c        interval [-1,1]
c  iw - the FORTRAN unit number on which the produced data are to be written
c 
c                            Output parameters:
c 
c  xs - nodes of the "universal" n-point discretization on the
c        interval [-1,1]
c  whts - the weights of the "universal" quadrature corresponding to the
c        nodes xs
c  xsaux - nodes of the "universal" naux-point discretization on the
c        interval [-1,1]
c  whtsaux - the weights of the "universal" quadrature corresponding to the
c        nodes xsaux
c  corrs - the weights of the quadrature connecting the values of a
c        function at the nodes xsaux with the integrals of that function
c        (the i-th column of the matrix corrs corresponds to the singularity
c        of that function sitting the the node xs(i)
c  ainterp - the matrix of interpolation coefficients, connecting the
c        values of a function at the nodes xs with the values of that
c        function at the nodes xsout
c 
        call unidisc0(ndigits,n,naux,corrs,ainterp,
     1      xs,whts,iw,w,xsaux,whts2aux)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine unidisc0(ndigits,n,naux,corrs,ainterp,
     1      xs,whts2,iw,w,xsaux,whts2aux)
        implicit real *8 (a-h,o-z)
        save
        dimension ts(1000),whts(1000),xs(1),
     1      dxdt(1000),whts2(1)
c 
        dimension tsaux(4000),whtsaux(4000),xsaux(4000),
     1      dxdtaux(4000),whts2aux(1),corrs(naux,n)
c 
        real *8 ainterp(1),w(1),ts22(1000),tsaux22(4000)
c 
c       construct the nodes of the discretization
c 
        call uninodes(ndigits,n,ts,whts,xs,dxdt,whts2)
c 
c       construct the auxiliary nodes
c 
        call uninodes(ndigits,naux,tsaux,whtsaux,
     1      xsaux,dxdtaux,whts2aux)
  
        call prin2('after second uninodes, xsaux=*',xsaux,naux)
        call prin2('after second uninodes, dxdtaux=*',dxdtaux,naux)
c 
c        construct the quadratures connecting the values of the
c        charge at auxiliary nodes with the resulting values of
c        the potential at the discretization nodes
c 
        iamatr=1
        lamatr=n*naux*2+2
c 
        iwww=iamatr+lamatr
c 
        call unimatr(n,naux,corrs,w(iamatr),w(iwww) )
c 
c       scale the matrix of corrections to accaount for the
c       "universal" change of variables
c 
        do 2400 i=1,n
c 
        do 2200 j=1,naux
c 
        corrs(j,i)=corrs(j,i)*dxdtaux(j)
 2200 continue
 2400 continue
c 
c        construct the matrix interpolating the values of the charge
c        at the support nodes to the values of the charge at the
c        auxiliary nodes
c 
        call uniinter(n,naux,ainterp,ts22,tsaux22,w)
c 
c       store the obtained data on disk
c 
        ifnew=0
        call unistore(ifnew,iw,ndigits,n,naux,xs,xsaux,whts2,
     1      whts2aux,corrs,ainterp)
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine unistore(ifnew,iw,ndigits,n,naux,xs,xsaux,
     1      whts,whtsaux,corrs,ainterp)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(n),xsaux(naux),whts(n),corrs(naux,n),
     1      ainterp(n,naux),whtsaux(naux)
c 
c       if the user so requested, skip whatever has already
c       been written on this file
c 
 1100 format(10a1)
        if(ifnew .eq. 1) goto 1160
        do 1130 i=1,10000000
        read(iw,1100,end=1150)
 1130 continue
 1150 continue
c 
 1160 continue
c 
        if(ifnew .eq. 1) rewind (iw)
c 
c       store the header of this record
c 
 1200 format('      In this run, ndigits = ',i2,', n=',i4,
     1    ', and naux= ',i4)
c 
        write(iw,1400)
        write(iw,1400)
        write(iw,1200) ndigits,n,naux
c 
 1400 format('          ')
        write(iw,1400)
c 
 1600 format('      support nodes are ')
 1800 format(1x,e21.15,1x,e21.15,1x,e21.15)
c 
        write(iw,1600)
        write(iw,1800) xs
c 
 2200 format('      auxiliary nodes are ')
c 
        write(iw,1400)
        write(iw,2200)
        write(iw,1800) xsaux
c 
 2400 format('      weights at support nodes are ')
c 
        write(iw,1400)
        write(iw,2400)
        write(iw,1800) whts
c 
 2600 format('      weights at auxiliary nodes are ')
c 
        write(iw,1400)
        write(iw,2600)
        write(iw,1800) whtsaux
c 
 2800 format('      matrix interpolating from support to ',
     1       'auxiliary nodes is ')
c 
        write(iw,1400)
        write(iw,2800)
        write(iw,1800) ainterp
c 
 3200 format('      matrix of weights based on auxiliary nodes is')
c 
        write(iw,1400)
        write(iw,3200)
        write(iw,1800) corrs
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine unimatr(n,ncorr,corrs,amatr,w)
        implicit real *8 (a-h,o-z)
        save
        dimension ts(1000),tscorr(4000),amatr(1),
     1      corrs(ncorr,n),w(1)
c 
c        construct the n Gaussian nodes on the interval [-1,1]
c 
        itype=0
        call legeexps(itype,n,ts,u,v,whts)
c 
c       allocate memory for the construction of the quadratures
c 
        iuu=1
        luu=n*ncorr*2+2
c 
        ivv=iuu+luu
        lvv=n*ncorr*2+2
c 
        iww=ivv+lvv
        lww=n*ncorr*2+2
c 
        itt=iww+lww
        ltt=n*ncorr*2+2
c 
        iw=itt+ltt
c 
c       one Gaussian node after another, construct the quadratures
c       based on the nodes tcorr for the evaluation of the
c       four types of integrals
c 
        errsmax=0
        do 1600 i=1,n/2+1
c 
        call prinf('in unimatr, i=*',i,1)
  
        call prini(0,0)
  
        call uniquad(n,ts(i),ncorr,amatr,tscorr,corrs(1,i),errmax,
     1        w(iuu),w(ivv),w(iww),w(itt),w(iw) )
  
        call prini(6,13)
  
        if(errsmax .lt. errmax) errsmax=errmax
  
 1600 continue
c 
c       reflect the first n/2 columns of the matrix corrs,
c       in  order to obtain the second half of the matrix
c 
        do 2400 i=1,n/2
        do 2200 j=1,ncorr
c 
        corrs(ncorr-j+1,n-i+1)=corrs(j,i)
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
        subroutine uniquad(n,x,ncorr,amatr,tscorr,sol,errmax,
     1        uu,vv,ww,tt,w)
        implicit real *8 (a-h,o-z)
        save
        dimension whts(1000),tscorr(1),
     1      tlarge(1000),whtlarge(1000),amatr(n*2,ncorr)
c 
        dimension w(1),hilcoefs(1000),coefsqua(1000),
     1      coefslog(1000),valslog(1000)
c 
        dimension uu(1),vv(1),ww(1),
     1      tt(1),rnorms(1000),rhs(1000),sol(1)
c 
        data nlarge/0/
c 
c 
c       decide if the subroutine for the evaluation of the
c       log, hilbert, and quadrupole transforms should be
c       initialized
c 
        ifinit=0
        if(nlarge .lt. n+2) ifinit=1
        if(nlarge .lt. n+2) nlarge=n+2
c 
c       initialize the subroutine for the evaluation
c       of Hilbert, log, and quadrupole transforms
c 
        if(ifinit .eq. 1)
     1      call qhlogini(nlarge,tlarge,whtlarge,w,ltot,lsave)
c 
c        for the user-specifeid point x, evaluate the
c        coefficients converting the values of a function at
c        the Gaussian nodes into the values of its three
c        transforms at x
c 
        call qhlogeva(x,hilcoefs,coefsqua,coefslog,w)
c 
c       construct the values at x of the three transforms
c       of the first n Legendre polynomials
c 
        do 1400 i=1,n
c 
        dlog=0
        do 1200 j=1,nlarge
  
        call legepol(tlarge(j),i-1,pol,der)
c 
        dlog=dlog+pol*coefslog(j)
 1200 continue
c 
        valslog(i)=dlog
 1400 continue
c 
c        construct the correction nodes
c 
        itype=0
        call legeexps(itype,ncorr,tscorr,u,v,whts)
c 
c       construct the values of the Legendre polynomials
c       multiplied by the integrand at the correction nodes
c 
        ii=0
        do 3000 iftype=1,2
        do 2400 i=1,n
c 
        ii=ii+1
c 
        do 2200 j=1,ncorr
c 
        call legepol(tscorr(j),i-1,pol,der)
c 
        if(iftype .eq. 1) amatr(ii,j)=pol
        if(iftype .eq. 2) amatr(ii,j)=pol*log( (tscorr(j)-x)**2 )/2
 2200 continue
 2400 continue
 3000 continue
c 
c        construct the right-hand side
c 
        do 3200 i=1,n
c 
        rhs(i)=0
        rhs(n+i)=valslog(i)
 3200 continue
c 
        rhs(1)=2
c 
c       use least squares to find the quadrature weights
c 
         eps=1.0d-13
         call leastsq(amatr,uu,ww,tt,n*2,ncorr,ncols,rnorms,eps,vv)
  
        call leastsq2(uu,ww,tt,n*2,ncorr,ncols,rhs,sol,vv)
  
        call prin2('after leastsq2, sol=*',sol,ncorr)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine uninodes(ndigits,n,ts,whts,xs,dxdt,whts2)
        implicit real *8 (a-h,o-z)
        save
        dimension ts(1),whts(1),xs(1),dxdt(1),whts2(1)
c 
c       construct the Gaussian nodes and weights
c 
        itype=1
        call legeexps(itype,n,ts,u,v,whts)
c 
c       construct the images of the Gaussian nodes and the
c       corresponding weights after prolate resampling
c 
        do 1200 i=1,n
c 
        call proresf(ndigits,ts(i),xs(i),dxdt(i))
c 
        whts2(i)=dxdt(i)*whts(i)
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine uniinter(n,ncorr,ainterp,ts,tscorr,w)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),ts(1),tscorr(1),ainterp(1)
c 
c         This subroutine constructs the matrix interpolating
c         a function from n Gaussian nodes on the interval [-1,1]
c         to ncorr Gaussian nodes. Please note that this is a
c         memory management routine; all work is done by the
c         subroutine uniinte0 (see)
c 
c                       Input parameters:
c 
c  n - the number of nodes from which the functions are to be
c       interpolated
c  ncorr - the number of nodes to which the functions are to be
c       interpolated
c 
c                       Output parameters:
c 
c  ainterp - the ncorr \times n matrix iterpolating from n Gaussian
c       nodes to ncorr of them
c  ts - the n Gaussian nodes on the interval [-1,1]
c  tscorr - the ncorr Gaussian nodes on the interval [-1,1]
c 
c                       Work arrays:
c 
c  w - must be at sufficiently long
c 
c        allocate memory for the main subroutine
c 
        iu=1
        lu=n**2+10
c 
        iuembed=iu+lu
        luembed=n*ncorr+10
c 
        iucorr=iuembed+luembed
        lucorr=ncorr**2+10
c 
        ivcorr=iucorr+lucorr
        lvcorr=ncorr**2+10
c 
        iwhtcorr=ivcorr+lvcorr
        lwhtcorr=ncorr+10
c 
        iv=iwhtcorr+lwhtcorr
        lv=n**2+10
c 
        iwhts=iv+lv
        lwhts=n+2
c 
c       construct the interpolation matrix
c 
        call uniinte0(n,ncorr,w(iu),w(iuembed),ainterp,
     1      ts,tscorr,w(iucorr),w(ivcorr),w(iwhtcorr),w(iv),w(iwhts))
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine uniinte0(n,ncorr,u,uembed,ainterp,
     1      ts,tscorr,ucorr,vcorr,whtscorr,v,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension u(n,n),uembed(ncorr,n),ainterp(ncorr,n)
        dimension v(1),ucorr(1),vcorr(1),
     1      ts(1),whts(1),tscorr(1),whtscorr(1)
c 
c         This subroutine constructs the matrix interpolating
c         a function from n Gaussian nodes on the interval [-1,1]
c         to ncorr Gaussian nodes
c 
c                       Input parameters:
c 
c  n - the number of nodes from which the functions are to be
c       interpolated
c  ncorr - the number of nodes to which the functions are to be
c       interpolated
c 
c                       Output parameters:
c 
c  ainterp - the ncorr \times n matrix iterpolating from n Gaussian
c       nodes to ncorr of them
c  ts - the n Gaussian nodes on the interval [-1,1]
c  tscorr - the ncorr Gaussian nodes on the interval [-1,1]
c 
c                       Work arrays:
c 
c  uembed - must be at least n*ncorr real *8 locations long
c  u - must be at least n*n real *8 elements
c  ucorr - must be at least ncorr*ncorr real *8 elements
c  ucorr - must be at least ncorr*ncorr real *8 elements
c  whtscorr - must be at least ncorr real *8 elements
c  v - must be at least n*n real *8 elements
c  whts - must be at least n*n real *8 elements
c 
c       . . . construct the n Gausian nodes on the interval [-1,1];
c             also the corresponding Gaussian expansion-evaluation
c             matrices
c 
        itype=2
        call legeexps(itype,n,ts,u,v,whts)
c 
c       construct the ncorr Gausian nodes on the interval [-1,1];
c       also the corresponding Gaussian expansion-evaluation
c       matrices
c 
        itype=2
        call legeexps(itype,ncorr,tscorr,ucorr,vcorr,whtscorr)
c 
cccc        call prin2('ts as created*',ts,n)
cccc        call prin2('tscorr as created*',tscorr,ncorr)
cccc        call prin2('u as constructed*',u,n**2)
cccc        call prin2('vcorr as constructed=*',vcorr,ncorr**2)
c 
c       embed the n*n- matrix converting n values into n
c       Legendre coefficients into an ncorr*n matrix converting
c       n values into ncorr coefficients (the last ncorr-n of
c       them being equal to zero).
c 
        do 1400 i=1,n
        do 1200 j=1,ncorr
c 
        uembed(j,i)=0
 1200 continue
 1400 continue
c 
        do 1800 i=1,n
        do 1600 j=1,n
c 
        uembed(j,i)=u(j,i)
 1600 continue
 1800 continue
c 
cccc        call prin2('uembed as constructed*',u,n*ncorr)
cccc        call prin2('while vcorr=*',vcorr,ncorr**2)
c 
c       multiply the embedded matrix by the matrix evaluating an
c       expansion of order ncorr at ncor Gaussian nodes, obtaining
c       the interpolating matrix from n to ncorr nodes
c 
        call unimatmu(vcorr,uembed,ncorr,ncorr,n,ainterp)
  
cccc        call prin2('and ainterp=*',ainterp,ncorr*n)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine unimatmu(a,b,n,m,k,c)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),b(m,k),c(n,k)
c 
        do 1600 i=1,n
        do 1400 j=1,k
        d=0
        do 1200 ii=1,m
c 
        d=d+a(i,ii)*b(ii,j)
 1200 continue
c 
        c(i,j)=d
 1400 continue
 1600 continue
        return
        end
  
c 
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the testing code and the beginning of the
c        prolate bell/taper code proper
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c        This file contains the following user-callable subroutines:
c 
c 
c   protaper2 - evaluates at the user-supplied point x \in [-1,1]
c        the taper function, together with its derivative; the geometry
c        of the taper is compatible with that assumed by most FFT
c        routines
c 
c   protaper -  evaluates at the user-supplied point x \in [-1,1]
c        the taper function, together with its derivative; the
c        geometry of the taper is the intuitive one, and is NOT
c        compatible with that assumed by most FFT routines
c 
c   probelev - evaluates the prolate bell at the user-supplied
c        point x \in [-1,1]. Please note that the prolate bell is NOT
c        the prolate function
c 
c   prolbret - returns to the user the coefficients of the Legendre
c        series of a prolate bell function.
c 
c   proltrap - constructs the "universal" quadrature rule based on
c        the prolate bell.
c 
c   proresf - constructs the image of the point t \in [-1,1] under
c        the "universal reparametrization" based on the the prolate
c        bell.
c 
c 
        subroutine protaper2(ndigits,wing,x,bell,der)
        implicit real *8 (a-h,o-z)
c 
c       This subroutine evaluates at the user-supplied point
c       x \in [-1,1] the taper function, together with its derivative;
c       the geometry of the taper is compatible with that assumed by
c       most FFT routines, so that the graph of the taper looks as follows:
c 
c 
c      *
c     ***
c      *
c      *
c      *
c   1> *               **********                 **********
c      *                         *               *
c      *                          *             *
c      *                           *           *
c   0> *                            ***********
c      *
c      *
c      *
c      *
c      *                                                                 *
c      ********************************************************************
c                      ^        ^        ^        ^         ^            *
c                     -1      -wing      0       wing       1
c 
c 
c 
c                         Input parameters:
c 
c  ndigits - the number of digits with which the bell is a bell, etc.;
c            must be between 1 and 16 (inclusive); however, since the
c            calculations are conducted in double precision, there is
c            no actual accuracy gain after ndigits = 14 or even ndigits=13
c  wing - the size of the wing (see cartoon above)
c  x - the point on the interval [-1,1] where the taper function and its
c            derivative are to be evaluated
c 
c                         Output parameters:
c  bell - the value of the taper function at the point x \in [-1,1]
c  der - the derivative of the taper function at the point x \in [-1,1]
c 
c 
c       . . . if the user-supplied point is not in either of the wings
c             - set the bell to 1 and get out
c 
        save
        if(abs(x) .lt. wing) goto 1400
c 
        bell=1
        der=0
        return
 1400 continue
c 
c       if the user-supplied point is on the left wing - act accordingly
c 
        if(x .gt. 0) goto 1800
c 
        beta=1
        alpha=beta/wing
c 
        t=alpha*x+beta
c 
        call probelev(ndigits,t,bell,der)
  
        der=der/wing*2
c 
        return
 1800 continue
c 
c       the user-supplied point is on the right wing - act accordingly
c 
        beta=1
        alpha=-beta/wing
  
c 
        t=alpha*x+beta
c 
        call probelev(ndigits,t,bell,der)
  
        der=der/wing*2
        return
        end
c 
c 
c 
c 
c 
        subroutine protaper(ndigits,wing,x,bell,der)
        implicit real *8 (a-h,o-z)
c 
c       This subroutine evaluates at the user-supplied point
c       x \in [-1,1] the taper function, together with its derivative;
c       the geometry of the taper is the intuitive one, and is NOT
c       compatible with that assumed by most FFT routines, so that
c       the graph of the taper looks as follows:
c 
c 
c      *
c     ***
c      *
c   1> *                       *******************
c      *                      *                   *
c      *                     *                     *
c      *                    *                       *
c   0> *              ******                         ******
c      *
c      *
c      *
c      *
c      *                                                                *
c      ********************************************************************
c                     ^        ^        ^        ^        ^             *
c                    -1     -1+wing     0      1-wing     1
c 
c 
c 
c                         Input parameters:
c 
c  ndigits - the number of digits with which the bell is a bell, etc.;
c            must be between 1 and 16 (inclusive); however, since the
c            calculations are conducted in double precision, there is
c            no actual accuracy gain after ndigits = 14 or even ndigits=13
c  wing - the size of the wing (see cartoon above)
c  x - the point on the interval [-1,1] where the taper function and its
c            derivative are to be evaluated
c 
c                         Output parameters:
c  bell - the value of the taper function at the point x \in [-1,1]
c  der - the derivative of the taper function at the point x \in [-1,1]
c 
c 
c       . . . if the user-supplied point is not in either of the wings
c             - set the bell to 1 and get out
c 
        save
        if(abs(x-1) .lt. wing) goto 1400
        if(abs(x+1) .lt. wing) goto 1400
c 
        bell=1
        der=0
        return
 1400 continue
c 
c       if the user-supplied point is on the left wing - act accordingly
c 
        if(x .gt. 0) goto 1800
c 
        b=1-wing
        alpha=1/(b-1)
        beta=1+alpha
c 
        t=alpha*x+beta
c 
        call probelev(ndigits,t,bell,der)
        der=der/wing*2
c 
        return
 1800 continue
c 
c       the user-supplied point is on the right wing - act accordingly
c 
        b=1-wing
        alpha=1/(1-b)
        beta=1-alpha
c 
        t=alpha*x+beta
c 
        call probelev(ndigits,t,bell,der)
        der=der/wing*2
        return
        end
c 
c 
c 
c 
c 
        subroutine probelev(ndigits,x,bell,der)
        implicit real *8 (a-h,o-z)
        save
        dimension coesbell(1000)
c 
        data ndigold/-1/
c 
c        This subroutine evaluates the prolate bell at the user-supplied
c        point x \in [-1,1]. Please note that the prolate bell is NOT
c        the prolate function, and has the folloeing appearance:
c 
c 
c      *
c     ***
c      *
c   1> *                           ***********
c      *                          *           *
c      *                         *             *
c      *                        *               *
c      *                       *                 *
c      *                      *                   *
c      *                     *                     *
c      *                    *                       *
c   0> *              ******                         ******
c      *
c      *
c      *
c      *
c      *                                                                *
c      ********************************************************************
c                     ^                 ^                 ^             *
c                    -1                 0                 1
c 
c 
c        with all derivatives equal to zero at the points -1, 0, 1.
c 
c 
c                         Input parameters:
c 
c  ndigits - the number of digits with which the bell is a bell, etc.;
c            must be between 1 and 16 (inclusive); however, since the
c            calculations are conducted in double precision, there is
c            no actual accuracy gain after ndigits = 14 or even ndigits=13
c  x - the point on the interval [-1,1] where the taper function and its
c            derivative are to be evaluated
c 
c                         Output parameters:
c  bell - the value of the bell function at the point x \in [-1,1]
c  der - the derivative of the bell function at the point x \in [-1,1]
c 
c 
c       . . . if the user-requested number of digits is different from
c             the one requested during the preceding call - retrieve
c             the appropriate coefficients from the subroutine prolbret
c 
        if(ndigold .eq. ndigits) goto 1100
c 
        call prolbret(ier,ndigits,coesbell,nterms)
        ndigold=ndigits
 1100 continue
c 
c               if the point is outside the interval [-1,1]
c             - set both the bell and its derivative to zero
c 
        if(abs(x) .lt. 1) goto 1200
        bell=0
        der=0
        return
 1200 continue
c 
c       calculate the bell and its derivative if x < 0
c 
        if(x .gt. 0) goto 2200
c 
        t=x*2+1
c 
        call proserev(t,bell,der,coesbell,Nterms)
c 
        return
 2200 continue
c 
c       calculate the bell and its derivative if x < 0
c 
        if(x .le. 0) goto 2400
c 
        t=-x*2+1
c 
        call proserev(t,bell,der,coesbell,Nterms)
c 
        return
 2400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine prolbret(ier,ndigits,coesbell,nterms)
        implicit real *8 (a-h,o-z)
        save
        dimension coesbell(1)
c 
        dimension bell6(21),bell7(22),bell8(24),bell9(25),
     1      bell10(26),bell11(27),bell12(28),bell13(28),
     2      bell14(29),bell15(30),bell16(31),
  
     3      bell1(12),bell2(15),bell3(17),bell4(18),bell5(20)
c 
c       This subroutine returns to the user the coefficients of the
c       Legendre series of a prolate bell function. The left wing
c       of the bell (i .e. for the values of x in the interval [-1,0])
c       is given by the expression
c 
c        b(x) = \int_{-1}^t \psi_0(\tau) d \tau,                       (1)
c 
c       with t=2*x+1; the right wing of the bell is given by the
c       same formula (1), with t=-2*x+1.
c 
c                        Input parameters:
c 
c  ndigits - the number of digits with which the bell to be constructed
c       is indeed a bell. Permitted values: 1-16; please note that the
c       subroutine treats this parameter conservatively, in that the
c       actuall precision tends to be higher than the one requested.
c 
c                        Output parameters:
c 
c  ier - errorr eturn code. Will be equal to zero for any integer
c       ndigits between 1 and 16 (inclusive). Will be equal to 16
c       for any other vlue of ndigits; please note that this is a
c       fatal error.
c  coesbell - the coefficients of the Legendre expansion of the prolate
c       bell. Following is a cartoon of the bell whose Legendre expansion
c       the array coesbell contains.
c 
c 
c 
c          *
c         ***
c          *
c          *
c          *
c          *
c          *
c          *
c          *
c          *
c          *
c          *
c      1   *                                   ************
c          *                                  *
c          *                                 *
c          *                                *
c          *                               *
c          *                              *
c          *                             *
c          *                            *
c          *                           *
c          *                          *
c          *                         *
c          *                        *
c          *                       *
c          *                      *
c          *                     *
c          *                    *
c      0   *        ************
c          *
c          *                                                      *
c          *********************************************************
c                                                                 *
c                   0                                      1
c 
c  nterms - the order (which is one less than the number of terms) of
c          the Legendre expsnsion returned in array coesbell
c 
c        In this cycle, eps=.100000E+00, nterms= 23
c        c= 0.400592968933142E+01
c 
         data bell1/
     1   0.60146380521739678E+00,-.11407161876618731E+00,
     2   0.13517050808556468E-01,-.95099637906870625E-03,
     3   0.43076235352447691E-04,-.13472944152457946E-05,
     4   0.30702429187461044E-07,-.53119133224436952E-09,
     5   0.72056414387414163E-11,-.78630974259126637E-13,
     6   0.70484836699800012E-15,-.52474422774879081E-17/
c 
c        In this cycle, eps=.100000E-01, nterms= 29,
c        c= 0.683119074050400E+01
c 
         data bell2/
     1   0.65282374627936021E+00,-.19141750864209432E+00,
     2   0.45122452021114502E-01,-.72865198821214535E-02,
     3   0.82086450595348221E-03,-.66943839670080946E-04,
     4   0.40964992557243430E-05,-.19402525396135610E-06,
     5   0.73010291375707074E-08,-.22309417793706759E-09,
     6   0.56385323117075611E-11,-.11972968057404566E-12,
     7   0.21647271671191502E-14,-.33712088197411763E-16,
     8   0.45141417098889022E-18/
c 
c        In this cycle, eps=.100000E-02, nterms= 33,
c        c= 0.947202580651210E+01
c 
         data bell3/
     1   0.67728240463901267E+00,-.23813171072185113E+00,
     2   0.76413180567054787E-01,-.18490113683364890E-01,
     3   0.33386092789302805E-02,-.45703189856104223E-03,
     4   0.48470073868898336E-04,-.40694706145776921E-05,
     5   0.27590122830683477E-06,-.15373497945007464E-07,
     6   0.71500879640970183E-09,-.28135054970914944E-10,
     7   0.94783706114451960E-12,-.27624886843876591E-13,
     8   0.70299448190778179E-15,-.15748628500911991E-16,
     9   0.30739584628451690E-18/
c 
c        In this cycle, eps=.100000E-03,nterms= 35,
c        c= 0.120226677411590E+02
c 
         data bell4/
     1   0.69160039102561541E+00,-.26956152938432798E+00,
     2   0.10339789628119517E+00,-.31880180988271425E-01,
     3   0.77165963648888549E-02,-.14722336781315063E-02,
     4   0.22409631057354861E-03,-.27610049913059973E-04,
     5   0.27937293542265681E-05,-.23536839784397071E-06,
     6   0.16719790988592986E-07,-.10129458216838715E-08,
     7   0.52877576297400086E-10,-.24004335615733336E-11,
     8   0.95551421249843187E-13,-.33600919195352383E-14,
     9   0.10508799304097005E-15,-.28685440219481282E-17/
c 
c        In this cycle, eps=.100000E-04, nterms= 39,
c        c= 0.145195500384292E+02
c 
         data bell5/
     1   0.70106421386446655E+00,-.29222790526085048E+00,
     2   0.12606530030258014E+00,-.45721761909315494E-01,
     3   0.13525407164381889E-01,-.32545522934381341E-02,
     4   0.64085583402647623E-03,-.10423306262155803E-03,
     5   0.14149797470105478E-04,-.16201274920117543E-05,
     6   0.15804337182945451E-06,-.13259577328740125E-07,
     7   0.96515742521948897E-09,-.61441160392107142E-10,
     8   0.34458343893981842E-11,-.17139928975007367E-12,
     9   0.76077856393777156E-14,-.30301942325283824E-15,
     *   0.10886045803195425E-16,-.34418073748126712E-18/
c 
c        In this cycle, eps=.100000E-05, nterms= 41,
c        c= 0.339623643706311E+02
c 
         data bell6/
     1   0.70781525509854548E+00,-.30938574978834029E+00,
     2   0.14510626605369960E+00,-.59148369608864034E-01,
     3   0.20253070337337257E-01,-.57860739292519763E-02,
     4   0.13820553515429167E-02,-.27759730353974528E-03,
     5   0.47232599840761851E-04,-.68618709240993842E-05,
     6   0.85797091383960354E-06,-.93042640854041027E-07,
     7   0.88156222194551613E-08,-.73480322384888451E-09,
     8   0.54227218191984347E-10,-.35642765554735733E-11,
     9   0.20980738472474320E-12,-.11116633714923516E-13,
     *   0.53268535296938006E-15,-.23184759578853445E-16,
     1   0.88792263511040825E-18/
c 
c        In this cycle, eps=.100000E-06, nterms= 43,
c        c= 0.194184111364534E+02
c 
         data bell7/
     1   0.71288918720492832E+00,-.32284814794006126E+00,
     2   0.16122281344246461E+00,-.71781780641734329E-01,
     3   0.27478699876994545E-01,-.89601424754312892E-02,
     4   0.24873539254537715E-02,-.58976533815988292E-03,
     5   0.12004427534253237E-03,-.21100541153305824E-04,
     6   0.32228851493449990E-05,-.43044067318777414E-06,
     7   0.50576695762774310E-07,-.52589440086320866E-08,
     8   0.48659622445645201E-09,-.40274662699573576E-10,
     9   0.29965220223385728E-11,-.20133386485280561E-12,
     *   0.12268448731417556E-13,-.68072476766827661E-15,
     1   0.34520921163033922E-16,-.15394944807872777E-17/
c 
c        In this cycle, eps=.100000E-07, nterms= 47,
c        c= 0.218381888425424E+02
c 
         data bell8/
     1   0.71685031990470713E+00,-.33370690910553134E+00,
     2   0.17499558466632966E+00,-.83491606722962008E-01,
     3   0.34899017915043219E-01,-.12639999135799853E-01,
     4   0.39576220637527893E-02,-.10727983124730570E-02,
     5   0.25262042780578100E-03,-.51902174277792886E-04,
     6   0.93493472772540754E-05,-.14840562131095893E-05,
     7   0.20863679311654354E-06,-.26106931245449034E-07,
     8   0.29216177628327297E-08,-.29375184030075404E-09,
     9   0.26651081599077831E-10,-.21908425186663199E-11,
     *   0.16381473104725301E-12,-.11182139925510939E-13,
     1   0.69922570303905333E-15,-.40181975911906923E-16,
     2   0.21285486368090922E-17,-.99695442497785624E-19/
c 
c        In this cycle, eps=.100000E-08, nterms= 49,
c        c= 0.242452552850318E+02
c 
         data bell9/
     1   0.72003343767234023E+00,-.34265971536139751E+00,
     2   0.18688008071724855E+00,-.94267770930574685E-01,
     3   0.42309529664897900E-01,-.16691420583463494E-01,
     4   0.57672876111191780E-02,-.17456479380086028E-02,
     5   0.46384356380344607E-03,-.10854243962020546E-03,
     6   0.22453146580589192E-04,-.41225880062866465E-05,
     7   0.67467421229515356E-06,-.98824906951249980E-07,
     8   0.13009815218519715E-07,-.15453945035179439E-08,
     9   0.16627862608328004E-09,-.16265027507891156E-10,
     *   0.14514797039221807E-11,-.11856181691385371E-12,
     1   0.88924411081945860E-14,-.61422989467764935E-15,
     2   0.39182869009462291E-16,-.23145779514829928E-17,
     3   0.12074398870557028E-18/
c 
c        In this cycle, eps=.100000E-09, nterms= 51,
c        c= 0.266429920398627E+02
c 
         data bell10/
     1   0.72265024218030161E+00,-.35017379777718885E+00,
     2   0.19722926395934373E+00,-.10415556117192125E+00,
     3   0.49578575336458822E-01,-.20996892345871915E-01,
     4   0.78761430947733715E-02,-.26148099232848047E-02,
     5   0.76920773101401832E-03,-.20095844856769457E-03,
     6   0.46761202481232916E-04,-.97229256089443292E-05,
     7   0.18127477478007131E-05,-.30411392656855757E-06,
     8   0.46070369377057073E-07,-.63240792649081321E-08,
     9   0.78927384619368381E-09,-.89852184786600232E-10,
     *   0.93596837596536784E-11,-.89480561109064473E-12,
     1   0.78736443128404064E-13,-.63942333012066392E-14,
     2   0.48050222879640828E-15,-.33494203692892882E-16,
     3   0.21708556140489336E-17,-.12409046796771768E-18/
c 
c        In this cycle, eps=.100000E-10, nterms= 53,
c        c= 0.290338967451443E+02
c 
         data bell11/
     1   0.72484142119894892E+00,-.35657417068550644E+00,
     2   0.20631721674407225E+00,-.11322296495754801E+00,
     3   0.56625276290538658E-01,-.25459655104575015E-01,
     4   0.10238059736725647E-01,-.36767620026505914E-02,
     5   0.11797438399406762E-02,-.33872898114094399E-03,
     6   0.87219645960996153E-04,-.20193087115546944E-04,
     7   0.42154985646098509E-05,-.79586321304261121E-06,
     8   0.13629420187652191E-06,-.21235846063102567E-07,
     9   0.30192575278250960E-08,-.39284802925673729E-09,
     *   0.46909530580066767E-10,-.51544938436199178E-11,
     1   0.52255387451556703E-12,-.48998405039421854E-13,
     2   0.42596914879735655E-14,-.34412283117310492E-15,
     3   0.25890301679971393E-16,-.18178368133896615E-17,
     4   0.11241276007466111E-18/
c 
c        In this cycle, eps=.100000E-11, nterms= 55,
c        c= 0.314198649407190E+02
c 
         data bell12/
     1   0.72670426557721097E+00,-.36209412808790998E+00,
     2   0.21435841402308278E+00,-.12154444193288677E+00,
     3   0.63403158503401920E-01,-.30002987863130041E-01,
     4   0.12806570318321467E-01,-.49208291720117231E-02,
     5   0.17018686057097149E-02,-.53029331002328495E-03,
     6   0.14911460401568650E-03,-.37916935615348814E-04,
     7   0.87390729170068300E-05,-.18301867367869556E-05,
     8   0.34917145167380390E-06,-.60844917506498784E-07,
     9   0.97090829819096541E-08,-.14223674484750220E-08,
     *   0.19178442827107007E-09,-.23858447724680960E-10,
     1   0.27448898521915346E-11,-.29272038686643836E-12,
     2   0.28998882512687375E-13,-.26744085642087478E-14,
     3   0.23007717818372287E-15,-.18499645976587764E-16,
     4   0.13928576577639397E-17,-.92230030902493926E-19/
c 
c        In this cycle, eps=.100000E-12, nterms= 55,
c        c= 0.338023678685776E+02
c 
         data bell13/
     1   0.72830828725692189E+00,-.36690560304442768E+00,
     2   0.22152234090951156E+00,-.12919309221073432E+00,
     3   0.69888646371097535E-01,-.34567629854762538E-01,
     4   0.15538054952952736E-01,-.63317833395768242E-02,
     5   0.23378101896359016E-02,-.78247173622024081E-03,
     6   0.23769438780614089E-03,-.65639154377989024E-04,
     7   0.16509211200052306E-04,-.37897708850065786E-05,
     8   0.79574590566744040E-06,-.15317486395358602E-06,
     9   0.27091807733202199E-07,-.44127421685386623E-08,
     *   0.66339166893931869E-09,-.92251803074226076E-10,
     1   0.11891940114038251E-10,-.14239906053687908E-11,
     2   0.15871337630553980E-12,-.16497495021326307E-13,
     3   0.16022735538618379E-14,-.14566511275881893E-15,
     4   0.12417406138995284E-16,-.92432572435485986E-18/
c 
c        In this cycle, eps=.100000E-13, nterms= 57,
c        c= 0.361825675685409E+02
c 
         data bell14/
     1   0.72970452194180375E+00,-.37113818321158001E+00,
     2   0.22794432192468493E+00,-.13623711641283905E+00,
     3   0.76073173888658709E-01,-.39108748518235550E-01,
     4   0.18393307094176820E-01,-.78919205446447851E-02,
     5   0.30863134341658496E-02,-.11002595666659075E-02,
     6   0.35785108270332000E-03,-.10631778126111510E-03,
     7   0.28898618622276331E-04,-.71989933421878563E-05,
     8   0.16466487174756885E-05,-.34650482896424909E-06,
     9   0.67213914343019794E-07,-.12042563927484434E-07,
     *   0.19968886912922435E-08,-.30705724877956433E-09,
     1   0.43868840900603996E-10,-.58342705642695605E-11,
     2   0.72362214317549357E-12,-.83851433131330231E-13,
     3   0.90936251920128139E-14,-.92453112843721271E-15,
     4   0.88261067505409799E-16,-.79242918588318514E-17,
     5   0.62042370690424047E-18/
  
  
c        In this cycle, eps=.100000E-14, nterms= 59,
c        c= 0.385613967628619E+02
c 
         data bell15/
     1   0.73093131634804077E+00,-.37489145105707255E+00,
     2   0.23373351933805321E+00,-.14273847584595978E+00,
     3   0.81957848967598258E-01,-.43593079784781377E-01,
     4   0.21338095693388081E-01,-.95825925931449890E-02,
     5   0.39434229425753760E-02,-.14868286534725882E-02,
     6   0.51388789920219987E-03,-.16297032430863914E-03,
     7   0.47481354997545499E-04,-.12727477499146370E-04,
     8   0.31438121487776447E-05,-.71679636845232585E-06,
     9   0.15111730093226813E-06,-.29510643447351294E-07,
     *   0.53476135423298256E-08,-.90079502415992665E-09,
     1   0.14129806087803155E-09,-.20674652859139730E-10,
     2   0.28266013690806569E-11,-.36168504670988642E-12,
     3   0.43384323912721639E-13,-.48859423261401731E-14,
     4   0.51740862850929786E-15,-.51597042277305565E-16,
     5   0.48521817601279962E-17,-.39738166235714609E-18/
c 
c        In this cycle, eps=.100000E-15, nterms= 61,
c        c= 0.409396126634395E+02
c 
         data bell16/
     1   0.73201805993908964E+00,-.37824323060823073E+00,
     2   0.23897885366502557E+00,-.14875264604747033E+00,
     3   0.87549858928615320E-01,-.47996457884380625E-01,
     4   0.24343141972939695E-01,-.11385263416623253E-01,
     5   0.49032248522179760E-02,-.19436649110042717E-02,
     6   0.70937488761531426E-03,-.23853278462416243E-03,
     7   0.73972604360034207E-04,-.21182030464030374E-04,
     8   0.56082167356631754E-05,-.13749056203636324E-05,
     9   0.31258725571379400E-06,-.66008083754632859E-07,
     *   0.12966919487326317E-07,-.23734447194912741E-08,
     1   0.40542452732477674E-09,-.64730445737997258E-10,
     2   0.96748428564570320E-11,-.13557282057622913E-11,
     3   0.17837631028848128E-12,-.22068077322927096E-13,
     4   0.25707813184543927E-14,-.28237916448093274E-15,
     5   0.29284868964257286E-16,-.28711498484446698E-17,
     6   0.24482549070221854E-18/
  
c 
c       process the 1-digit case
c 
        if(ndigits .ne. 1) goto 1100
c 
        nterms=23
        nn=12
        call prolcpy2(bell1,nn,coesbell,nterms)
        return
 1100 continue
c 
c       process the 2-digit case
c 
        if(ndigits .ne. 2) goto 1200
c 
        nterms=29
        nn=15
        call prolcpy2(bell2,nn,coesbell,nterms)
        return
 1200 continue
c 
c       process the 3-digit case
c 
        if(ndigits .ne. 3) goto 1300
c 
        nterms=33
        nn=17
        call prolcpy2(bell3,nn,coesbell,nterms)
        return
 1300 continue
c 
c       process the 4-digit case
c 
        if(ndigits .ne. 4) goto 1400
c 
        nterms=35
        nn=18
        call prolcpy2(bell4,nn,coesbell,nterms)
        return
 1400 continue
c 
c       process the 5-digit case
c 
        if(ndigits .ne. 5) goto 1500
c 
        nterms=39
        nn=20
        call prolcpy2(bell5,nn,coesbell,nterms)
        return
 1500 continue
c 
c       process the 6-digit case
c 
        if(ndigits .ne. 6) goto 1600
c 
        nterms=41
        nn=21
        call prolcpy2(bell6,nn,coesbell,nterms)
        return
 1600 continue
c 
c       process the 7-digit case
c 
        if(ndigits .ne. 7) goto 1700
c 
        nterms=43
        nn=22
        call prolcpy2(bell7,nn,coesbell,nterms)
        return
 1700 continue
c 
c       process the 8-digit case
c 
        if(ndigits .ne. 8) goto 1800
c 
        nterms=47
        nn=24
        call prolcpy2(bell8,nn,coesbell,nterms)
        return
 1800 continue
c 
c       process the 9-digit case
c 
        if(ndigits .ne. 9) goto 1900
c 
        nterms=49
        nn=25
        call prolcpy2(bell9,nn,coesbell,nterms)
        return
 1900 continue
c 
c       process the 10-digit case
c 
        if(ndigits .ne. 10) goto 2000
c 
        nterms=51
        nn=26
        call prolcpy2(bell10,nn,coesbell,nterms)
        return
 2000 continue
c 
c       process the 11-digit case
c 
        if(ndigits .ne. 11) goto 2100
c 
        nterms=53
        nn=27
        call prolcpy2(bell11,nn,coesbell,nterms)
        return
 2100 continue
c 
c       process the 12-digit case
c 
        if(ndigits .ne. 12) goto 2200
c 
        nterms=55
        nn=28
        call prolcpy2(bell12,nn,coesbell,nterms)
        return
 2200 continue
c 
c       process the 13-digit case
c 
        if(ndigits .ne. 13) goto 2300
c 
        nterms=55
        nn=28
        call prolcpy2(bell13,nn,coesbell,nterms)
        return
 2300 continue
c 
c       process the 14-digit case
c 
        if(ndigits .ne. 14) goto 2400
c 
        nterms=57
        nn=29
        call prolcpy2(bell14,nn,coesbell,nterms)
        return
 2400 continue
c 
c       process the 15-digit case
c 
        if(ndigits .ne. 15) goto 2500
c 
        nterms=59
        nn=30
        call prolcpy2(bell15,nn,coesbell,nterms)
        return
 2500 continue
c 
c       process the 16-digit case
c 
        if(ndigits .ne. 16) goto 2600
c 
        nterms=61
        nn=31
        call prolcpy2(bell16,nn,coesbell,nterms)
        return
 2600 continue
c 
        ier=16
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine prolcpy2(bell,nn,coesbell,nterms)
        implicit real *8 (a-h,o-z)
        save
        dimension bell(1),coesbell(1)
c 
        do 1200 i=1,nterms+1
        coesbell(i)=0
 1200 continue
c 
        coesbell(1)=0.5d0
        do 1400 i=1,nn
        coesbell(i*2)=bell(i)
 1400 continue
        return
        end
  
  
  
  
c 
c 
c 
c 
c 
      SUBROUTINE proserev(X,VAL,der,PEXP,N)
      IMPLICIT REAL *8 (A-H,O-Z)
        save
      REAL *8 PEXP(1)
C 
C     This subroutine computes the value and the derivative
c     of a gaussian expansion with coefficients PEXP
C     at point X in interval [-1,1].
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
        done=1
        pjm2=1
        pjm1=x
        derjm2=0
        derjm1=1
c 
        val=pexp(1)*pjm2+pexp(2)*pjm1
        der=pexp(2)
c 
        DO 600 J = 2,N
c 
        pj= ( (2*j-1)*x*pjm1-(j-1)*pjm2 ) / j
        val=val+pexp(j+1)*pj
c 
        derj=(2*j-1)*(pjm1+x*derjm1)-(j-1)*derjm2
c 
        derj=derj/j
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
        subroutine proltrap(ndigits,n,ts,xs,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension coesbell(100),xs(1),ts(1),whts(1)
c 
        data ndigold/-1/
c 
c       This subroutine constructs the "universal" quadrature
c       rule based on the prolate bell. To be more specific,
c       the subroutine re-parametrizes the interval [-1,1],
c       replacing the parameter $x \in [-1,1]$ with the parameter
c       $t \in [-1,1]$, in such a manner that all derivatives
c       of x with respect to t are equal to zero at the points
c       -1,1 (to precision 1.0E-(ndigits)); after that, it
c       constructs the trapezoidal discretization in t on the
c       interval [-1,1], finds the corresponding points in
c       x (obviously, these cluster near the ends), and evaluates
c       the corresponding weights (by multiplying dxtd by the
c       sampling interval in t, what a surprise!)
c 
c       The quadrature produced by this subroutine is "universal"
c       in the sense that it integrates effectively functions with
c       a wide class of singularities at the ends.
c 
c       IMPORTANT NOTE: The n trapezoidal nodes produced by this
c       subroutine include the ends. However, the weights associated
c       with the end-nodes are quite small, so that the end-nodes
c       can be omitted with a fairly minimal penalty in accuracy.
c       When the function to be integrated is infinite at one or
c       both of the points -1, 1, the corresponding end-point HAS
c       TO BE OMITTED (duh!); otherwise, the quadrature will produce
c       an infinite result (of a NaN, or something equally useful).
c       Thus, the produced quadrature can (and often should) be viewed
c       as an (n-2) - point one
c 
c                     Input parameters:
c 
c  ndigits - the number of digits of accuracy to which everything
c       is done; has to be between 1 and 16 (inclusive)
c  n - the number of nodes in the discretization, including the
c       ends (see IMPORTANT NOTE above)
c 
c                     Output parameters:
c 
c  ts - the n trapezoidal nodes on the interval [-1,1], including
c       the ends (see IMPORTANT NOTE above)
c  xs - images of the nodes ts under the parametrization
c  whts - the weights corresponding to the nodes xs
c 
c       . . . if the user-requested number of digits is different from
c             the one requested during the preceding call - retrieve
c             the coefficients of the Legendre expansion of the
c             bell from the subroutine prolbret
c 
        if(ndigold .eq. ndigits) goto 1100
c 
        call prolbret(ier,ndigits,coesbell,nterms)
        ndigold=ndigits
 1100 continue
  
  
        call prin2('in proltrap, coesbell=*',coesbell,nterms)
        call prinf('in proltrap, n=*',n,1)
c 
c       construct the equispaced discretization of the "source"
c       domain
c 
        h=2
        h=h/(n-1)
c 
        do 1200 i=1,n
c 
        ts(i)=(i-1)*h-1
c 
        call proserev(ts(i),xs(i),whts(i),coesbell,Nterms)
c 
        xs(i)=xs(i)*2-1
c 
        whts(i)=whts(i)*h*2
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine proresf(ndigits,t,x,dxdt)
        implicit real *8 (a-h,o-z)
        save
        dimension coesbell(100)
c 
        data ndigold/-1/
c 
c       This subroutine constructs the image of the point
c       t \in [-1,1] under the "universal reparametrization"
c       based on the the prolate bell. To be more specific,
c       the subroutine re-parametrizes the interval [-1,1],
c       replacing the parameter $x \in [-1,1]$ with the parameter
c       $t \in [-1,1]$, in such a manner that all derivatives
c       of x with respect to t are equal to zero at the points
c       -1,1 (to precision 1.0E-(ndigits)); after that, it
c       constructs the image of the point t under the obtained
c       reparametrization, and also the detivative dxdt at the
c       point t.
c 
c       . . . if the user-requested number of digits is different from
c             the one requested during the preceding call - retrieve
c             the coefficients of the Legendre expansion of the
c             bell from the subroutine prolbret
c 
        if(ndigold .eq. ndigits) goto 1100
c 
        call prolbret(ier,ndigits,coesbell,nterms)
        ndigold=ndigits
c 
 1100 continue
c 
c       evaluate the resampling function and its derivative
c 
        call proserev(t,x,dxdt,coesbell,Nterms)
c 
        x=x*2-1
        dxdt=dxdt*2
        return
        end
  
c 
c 
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        this is the end of the debugging code and the beginning
c        of the legendre expansion routines
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c        This file contains a set of subroutines for the handling
c        of Legendre expansions. It contains 19 subroutines that are
c        user-callable. Following is a brief description of these
c        subroutines.
c 
c   legeexps - constructs Legendre nodes, and  corresponding Gaussian
c        weights. Also constructs the matrix v converting the
c         coefficients of a legendre expansion into its values at
c         the n Gaussian nodes, and its inverse u, converting the
c         values of a function at n Gaussian nodes into the
c         coefficients of the corresponding Legendre series.
c 
c   legepol - evaluates a single Legendre polynomial (together
c         with its derivative) at the user-provided point
c 
c   legepols - evaluates a bunch of Legendre polynomials
c         at the user-provided point
c 
c   legeinmt - for the user-specified n, constructs the matrices of
c        spectral indefinite integration differentiation on the n
c        Gaussian nodes on the interval [-1,1].
c 
c   legeinte - computes the indefinite integral of the legendre
c        expansion polin getting the expansion polout
c 
c   legediff -  differentiates the legendre expansion polin getting
c        the expansion polout
c 
c   legefder - computes the value and the derivative of a Legendre
c        expansion at point X in interval [-1,1]; this subroutine
c        is not designed to be very efficient, but it does not
c        use any exdternally supplied arrays
c 
c   legefde2 - the same as legefder, except it is desigmed to be
c        fairly efficient; it uses externally supplied arrays
c        that are precomputed
c 
c   legeexev - computes the value of a Legendre expansion with
c        at point X in interval [-1,1]; same as legefder, but does
c        not compute the derivative of the expansion
c 
c   legeexe2 - the same as legeexev, except it is desigmed to be
c        fairly efficient; it uses externally supplied arrays
c        that are precomputed
c 
c   lematrin - constructs the matrix interpolating functions from
c        the n-point Gaussian grid on the interval [-1,1] to an
c        arbitrary m-point grid (the nodes of the latter are
c        user-provided)
c 
c   levecin - constructs the coefficients of the standard
c        interpolation formula connecting the values of a
c        function at n Gaussian nodes on the interval [a,b] with
c        its value at the point x \in R^1
c 
c   legeodev - evaluates at the point x a Legendre expansion
c        having only odd-numbered elements; this is a fairly
c        efficient code, using external arrays that are
c        precomputed
c 
c   legeevev - evaluates at the point x a Legendre expansion
c        having only even-numbered elements; this is a fairly
c        efficient code, using external arrays that are
c        precomputed
c 
c   legepeven - evaluates even-numbered Legendre polynomials
c        of the argument x; this is a fairly efficient code,
c        using external arrays that are precomputed
c 
c   legepodd - evaluates odd-numbered Legendre polynomials
c        of the argument x; this is a fairly efficient code,
c        using external arrays that are precomputed
c 
C   legefdeq - computes the value and the derivative of a
c        Legendre Q-expansion with coefficients coefs
C     at point X in interval (-1,1); please note that this is
c     the evil twin of the subroutine legefder, evaluating the
c     proper (P-function) Legendre expansion; this subroutine
c        is not designed to be very efficient, but it does not
c        use any exdternally supplied arrays
c 
c   legeq - calculates the values and derivatives of a bunch
c        of Legendre Q-functions at the user-specified point
c        x on the interval (-1,1)
c 
c   legeqs - calculates the value and the derivative of a single
c        Legendre Q-function at the user-specified point
c        x on the interval (-1,1)
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
        subroutine legeinmt(n,ainte,adiff,x,whts,endinter,
     1      itype,w)
        implicit real *8 (a-h,o-z)
        save
        dimension ainte(1),w(1),x(1),whts(1),adiff(1),endinter(1)
c 
c 
c        for the user-specified n, this subroutine constructs
c        the matrices of spectral indefinite integration and/or
c        spectral differentiation on the n Gaussian nodes
c        on the interval [-1,1]. Actually, this is omnly a
c        memory management routine. All the actual work is done
c        by the subroutine legeinm0 (see)
c 
c                           input parameters:
c 
c  n - the number of Gaussian nodes on the interval [-1,1]
c  itype - the type of the calculation to be performed
c          EXPLANATION:
c       itype=1 means that only the matrix ainte will
c               be constructed
c       itype=2 means that only the matrix adiff will
c               be constructed
c       itype=3 means that both matrices ainte and adiff
c               will be constructed
c 
c                           output paramaters:
c 
c  ainte - the matrix of spectral indefinite integration on
c          the Gaussian nodes
c  adiff - the matrix of spectral differentiation on
c          the Gaussian nodes
c  x - the n Gaussian nodes on the intervl [-1,1]
c  whts - the n Gaussian weights on the interval [-1,1]
c  endinter - the interpolation coefficients converting the
c          values of a function at n Gaussian nodes into its
c          value at 1 (the right end of the interval)
c 
c                           work arrays:
c 
c  w - must be 3* n**2 + 2*n +50 *8 locations long
c 
c        . . . allocate memory for the construction of the integrating
c              matrix
c 
        ipolin=1
        lpolin=n+5
c 
        ipolout=ipolin+lpolin
        lpolout=n+5
c 
        iu=ipolout+lpolout
        lu=n**2+1
c 
        iv=iu+lu
        lv=n**2+1
c 
        iw=iv+lv
        lw=n**2+1
c 
        ltot=iw+lw
c 
c        construct the integrating matrix
c 
        call legeinm0(n,ainte,adiff,w(ipolin),w(ipolout),
     1      x,whts,w(iu),w(iv),w(iw),itype,endinter)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine legeinm0(n,ainte,adiff,polin,polout,
     1      x,whts,u,v,w,itype,endinter)
        implicit real *8 (a-h,o-z)
        save
        dimension ainte(n,n),u(n,n),v(n,n),w(n,n),
     1      endinter(1),x(n),whts(n),polin(n),polout(n),
     2      adiff(n,n)
c 
c        for the user-specified n, this subroutine constructs
c        the matrices of spectral indefinite integration and/or
c        spectral differentiation on the n Gaussian nodes
c        on the interval [-1,1]
c 
c                           input parameters:
c 
c  n - the number of Gaussian nodes on the interval [-1,1]
c  itype - the type of the calculation to be performed
c          EXPLANATION:
c       itype=1 means that only the matrix ainte will
c               be constructed
c       itype=2 means that only the matrix adiff will
c               be constructed
c       itype=3 means that both matrices ainte and adiff
c               will be constructed
c 
c                           output paramaters:
c 
c  ainte - the matrix of spectral indefinite integration on
c          the Gaussian nodes
c  adiff - the matrix of spectral differentiation on
c          the Gaussian nodes
c  x - the n Gaussian nodes on the intervl [-1,1]
c  whts - the n Gaussian weights on the interval [-1,1]
c 
c                           work arrays:
c 
c  polin, polout - must be n+3 real *8 locations each
c 
c  u, v, w - must be n**2+1 real *8 locations each
c 
c        . . . construct the matrices of the forward and inverse
c              Legendre transforms
c 
        itype2=2
        call legeexps(itype2,n,x,u,v,whts)
c 
cccc         call prin2('after legeexps, u=*',u,n*n)
c 
c        if the user so requested,
c        construct the matrix converting the coefficients of
c        the Legendre series of a function into the coefficients
c        of the indefinite integral of that function
c 
        if(itype. eq. 2) goto 2000
c 
        do 1600 i=1,n
c 
        do 1200 j=1,n+2
        polin(j)=0
 1200 continue
c 
        polin(i)=1
c 
        call legeinte(polin,n,polout)
c 
        do 1400 j=1,n
        ainte(j,i)=polout(j)
 1400 continue
c 
 1600 continue
c 
cccc         call prin2('ainte initially is*',ainte,n*n)
c 
c        multiply the three, obtaining the integrating matrix
c 
        call matmul(ainte,u,w,n)
        call matmul(v,w,ainte,n)
c 
 2000 continue
c 
c        if the user so requested,
c        construct the matrix converting the coefficients of
c        the Legendre series of a function into the coefficients
c        of the derivative of that function
c 
        if(itype. eq. 1) goto 3000
c 
        do 2600 i=1,n
c 
        do 2200 j=1,n
        polin(j)=0
 2200 continue
c 
        polin(i)=1
c 
        call legediff(polin,n,polout)
c 
        do 2400 j=1,n
        adiff(j,i)=polout(j)
cccc        ainte(i,j)=polout(j)
 2400 continue
c 
 2600 continue
c 
cccc         call prin2('adiff initially is*',adiff,n*n)
c 
c        multiply the three, obtaining the integrating matrix
c 
        call matmul(adiff,u,w,n)
        call matmul(v,w,adiff,n)
c 
 3000 continue
c 
c        construct the vector of interpolation coefficients
c        converting the values of a polynomial at the Gaussian
c        nodes into its value at the right end of the interval
c 
        do 3400 i=1,n
c 
        d=0
        do 3200 j=1,n
        d=d+u(j,i)
 3200 continue
        endinter(i)=d
 3400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine legeinte(polin,n,polout)
        implicit real *8 (a-h,o-z)
        save
        dimension polin(1),polout(1)
c 
c       this subroutine computes the indefinite integral of the
c       legendre expansion polin getting the expansion polout
c 
c 
c                       input parameters:
c 
c  polin - the legendre expansion to be integrated
c  n - the order of the expansion polin
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c         also nothe that the order of the integrated expansion is
c         n+1 (who could think!)
c 
c                       output parameters:
c 
c  polout - the legendre expansion of the integral of the function
c         represented by the expansion polin
c 
        do 1200 i=1,n+2
        polout(i)=0
 1200 continue
c 
        do 2000 k=2,n+1
        j=k-1
c 
cccc        polout(k+1)=polin(k)/(2*j+1)+polout(k+1)
        polout(k+1)=polin(k)/(2*j+1)
        polout(k-1)=-polin(k)/(2*j+1)+polout(k-1)
c 
 2000 continue
c 
        polout(2)=polin(1)+polout(2)
c 
        dd=0
        sss=-1
        do 2200 k=2,n+1
c 
        dd=dd+polout(k)*sss
        sss=-sss
 2200 continue
c 
        call prin2('dd=*',dd,1)
        polout(1)=-dd
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine legediff(polin,n,polout)
        implicit real *8 (a-h,o-z)
        save
        dimension polin(1),polout(1)
c 
c       this subroutine differentiates the legendre
c       expansion polin getting the expansion polout
c 
c 
c                       input parameters:
c 
c  polin - the legendre expansion to be differentiated
c  n - the order of the expansion polin
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c         also nothe that the order of the integrated expansion is
c         n+1 (who could think!)
c 
c                       output parameters:
c 
c  polout - the legendre expansion of the derivative of the function
c         represented by the expansion polin
c 
        do 1200 k=1,n+1
        polout(i)=0
 1200 continue
c 
        pk=polin(n+1)
        pkm1=polin(n)
        pkm2=0
        do 2000 k=n+1,2,-1
c 
        j=k-1
c 
        polout(k-1)=pk*(2*j-1)
        if(k .ge. 3) pkm2=polin(k-2)+pk
c 
        pk=pkm1
        pkm1=pkm2
c 
 2000 continue
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
C     at point X in interval [-1,1].
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
  
  
        done=1
        pjm2=1
        pjm1=x
        derjm2=0
        derjm1=1
c 
        val=pexp(1)*pjm2+pexp(2)*pjm1
        der=pexp(2)
c 
        DO 600 J = 2,N
c 
        pj= ( (2*j-1)*x*pjm1-(j-1)*pjm2 ) / j
        val=val+pexp(j+1)*pj
c 
        derj=(2*j-1)*(pjm1+x*derjm1)-(j-1)*derjm2
c 
        derj=derj/j
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
      SUBROUTINE legeFDE2(X,VAL,der,PEXP,N,
     1    pjcoefs1,pjcoefs2,ninit)
      IMPLICIT REAL *8 (A-H,O-Z)
        save
      REAL *8 PEXP(1),pjcoefs1(1),pjcoefs2(1)
c 
C     This subroutine computes the value and the derivative
c     of a gaussian expansion with coefficients PEXP
C     at point X in interval [-1,1].
c 
c                input parameters:
c 
C  X - evaluation point
C  PEXP - expansion coefficients
C  N  - order of expansion
c  pjcoefs1, pjcoefs2 - two arrays precomputed on a previous call
c      on a previous call to this subroutine. Please note that this
c      is only an input parameter if the parameter ninit (see below)
c      has been set to 0; otherwise, these are output parameters
c  ninit - tells the subroutine whether and to what maximum order the
c       arrays coepnm1,coepnp1,coexpnp1 should be initialized.
c     EXPLANATION: The subroutine will initialize the first ninit
c       elements of each of the arrays pjcoefs1, pjcoefs2. On the first
c       call to this subroutine, ninit should be set to the maximum
c       order n for which this subroutine might have to be called;
c       on subsequent calls, ninit should be set to 0. PLEASE NOTE
c       THAT THAT THESE ARRAYS USED BY THIS SUBROUTINE
c       ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES USED BY THE
c       SUBROUTINE LEGEEXE2. If these arrays have been initialized
c       by one of these two subroutines, they do not need to be
c       initialized by the other one.
c 
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c 
c                output parameters:
c 
C  VAL - computed value
C  der - computed value of the derivative
C 
C 
        if(ninit .eq. 0) goto 1400
c 
        done=1
        do 1200 j=2,ninit
c 
        pjcoefs1(j)=(2*j-done)/j
        pjcoefs2(j)=-(j-done)/j
c 
 1200 continue
c 
        ifcalled=1
 1400 continue
c 
        pjm2=1
        pjm1=x
        derjm2=0
        derjm1=1
c 
        val=pexp(1)*pjm2+pexp(2)*pjm1
        der=pexp(2)
c 
        DO 1600 J = 2,N
c 
cccc        pj= ( (2*j-1)*x*pjm1-(j-1)*pjm2 ) / j
c 
        pj= pjcoefs1(j)*x*pjm1+pjcoefs2(j)*pjm2
  
  
        val=val+pexp(j+1)*pj
c 
cccc        derj=(2*j-1)*(pjm1+x*derjm1)-(j-1)*derjm2
        derj=pjcoefs1(j)*(pjm1+x*derjm1)+pjcoefs2(j)*derjm2
  
ccc         call prin2('derj=*',derj,1)
  
  
cccc        derj=derj/j
        der=der+pexp(j+1)*derj
c 
        pjm2=pjm1
        pjm1=pj
        derjm2=derjm1
        derjm1=derj
 1600   CONTINUE
c 
      RETURN
      END
c 
c 
c 
c 
c 
      SUBROUTINE legeexev(X,VAL,PEXP,N)
      IMPLICIT REAL *8 (A-H,O-Z)
        save
      REAL *8 PEXP(1)
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
        der=pexp(2)
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
      SUBROUTINE legeexe2(X,VAL,PEXP,N,
     1      pjcoefs1,pjcoefs2,ninit)
      IMPLICIT REAL *8 (A-H,O-Z)
        save
      REAL *8 PEXP(1),pjcoefs1(1),pjcoefs2(1)
c 
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
        if(ninit .eq. 0) goto 1400
c 
        done=1
        do 1200 j=2,ninit
c 
        pjcoefs1(j)=(2*j-done)/j
        pjcoefs2(j)=-(j-done)/j
c 
 1200 continue
c 
        ifcalled=1
 1400 continue
c 
        pjm2=1
        pjm1=x
c 
        val=pexp(1)*pjm2+pexp(2)*pjm1
        der=pexp(2)
c 
        DO 600 J = 2,N
c 
        pj= pjcoefs1(j)*x*pjm1+pjcoefs2(j)*pjm2
  
cccc        pj= ( (2*j-1)*x*pjm1-(j-1)*pjm2 ) / j
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
        subroutine lematrin(n,m,xs,amatrint,ts,w)
        implicit real *8 (a-h,o-z)
        save
        dimension amatrint(m,n),xs(1),w(1),ts(1)
c 
c 
c        This subroutine constructs the matrix interpolating
c        functions from the n-point Gaussian grid on the interval [-1,1]
c        to an arbitrary m-point grid (the nodes of the latter are
c        user-provided)
c 
c                 Input parameters:
c 
c  n - the number of interpolation nodes
c  m - the number of nodes to which the functions will be interpolated
c  xs - the points at which the function is to be interpolated
c 
c                  Output parameters:
c 
c  amatrint - the m \times n matrix conerting the values of a function
c        at the n Legendre nodes into its values at m user-specified
c        (arbitrary) nodes
c  ts - the n Gaussian nodes on the interval [-1,1]
c 
c                  Work arrays:
c 
c  w - must be at least 2*n**2+n + 100 real *8 locations long
c 
  
        icoefs=1
        lcoefs=n+2
c 
        iu=icoefs+lcoefs
        lu=n**2+10
c 
        iv=iu+lu
c 
        ifinit=1
        do 2000 i=1,m
c 
        call levecin(n,xs(i),ts,w(iu),w(iv),w(icoefs),ifinit)
c 
        do 1400 j=1,n
        amatrint(i,j)=w(j)
 1400 continue
c 
        ifinit=0
 2000 continue
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine levecin(n,x,ts,u,v,coefs,ifinit)
        implicit real *8 (a-h,o-z)
        save
        dimension u(n,n),v(n,n),ts(1),coefs(1)
c 
c        This subroutine constructs the coefficients of the
c        standard interpolation formula connecting the values of a
c        function at n Gaussian nodes on the interval [a,b] with
c        its value at the point x \in R^1
c 
c                 Input parameters:
c 
c  n - the number of interpolation nodes
c  x - the points at which the function is to be interpolated
c  ts - the n Gaussian nodes on the interval [-1,1]; please note that
c        it is an input parameter only if the parameter ifinit (see
c        below) has been set to 1; otherwise, it is an output parameter
c  u - the n*n matrix converting the  values at of a polynomial of order
c         n-1 at n legendre nodes into the coefficients of its
c        legendre expansion; please note that
c        it is an input parameter only if the parameter ifinit (see
c        below) has been set to 1; otherwise, it is an output parameter
c  ifinit - an integer parameter telling the subroutine whether it should
c        initialize the Legendre expander;
c     ifinit=1 will cause the subroutine to perform the initialization
c     ifinit=0 will cause the subroutine to  skip the initialization
c 
c                  Output parameters:
c 
c  coefs - the interpolation coefficients
c 
c                 Work arrays:
c 
c  v - must be at least n*n real *8 locations long
c 
c       . . . construct the n Gausian nodes on the interval [-1,1];
c             also the corresponding Gaussian expansion-evaluation
c             matrices
c 
        itype=2
        if(ifinit .ne.0) call legeexps(itype,n,ts,u,v,coefs)
c 
c       evaluate the n Legendre polynomials at the point where the
c       functions will have to be interpolated
c 
        call legepols(x,n+1,v)
c 
c       apply the interpolation matrix to the ector of values
c       of polynomials from the right
c 
        call lematvec(u,v,coefs,n)
        return
        end
c 
c 
c 
c 
c 
        subroutine lematvec(a,x,y,n)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),x(n),y(n)
c 
        do 1400 i=1,n
        d=0
        do 1200 j=1,n
        d=d+a(j,i)*x(j)
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
        subroutine matmul(a,b,c,n)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),b(n,n),c(n,n)
c 
        do 2000 i=1,n
        do 1800 j=1,n
        d=0
        do 1600 k=1,n
        d=d+a(i,k)*b(k,j)
 1600 continue
        c(i,j)=d
 1800 continue
 2000 continue
        return
c 
c 
c 
c 
        entry matmua(a,b,c,n)
ccc          call prin2('in matmua, a=*',a,n**2)
ccc          call prin2('in matmua, b=*',b,n**2)
        do 3000 i=1,n
        do 2800 j=1,n
        d=0
        do 2600 k=1,n
        d=d+a(i,k)*b(j,k)
 2600 continue
        c(i,j)=d
 2800 continue
 3000 continue
ccc          call prin2('exiting, c=*',c,n**2)
        return
        end
  
c 
c 
c 
c 
c 
        subroutine legeodev(x,nn,coefs,val,ninit,
     1      coepnm1,coepnp1,coexpnp1)
        implicit real *8 (a-h,o-z)
        save
        dimension coepnm1(1),coepnp1(1),
     1            coexpnp1(1),coefs(1)
c 
c 
c       This subroutine evaluates at the point x a Legendre expansion
c       having only odd-numbered elements
c 
c                  Input parameters:
c 
c  x - point on the interval [-1,1] at which the Legendre expansion
c       is to be evaluated
c  nn - order of the expansion to be evaluated
c  coefs - odd-numbered coefficients of the Legendre expansion
c       to be evaluated at the point x (nn/2+2 of them things)
c  ninit - tells the subroutine whether and to what maximum order the
c       arrays coepnm1,coepnp1,coexpnp1 should be initialized.
c     EXPLANATION: The subroutine will initialize the first ninit/2+2
c                  (or so) elements of each of the arrays  coepnm1,
c       coepnp1, coexpnp1. On the first call to this subroutine, ninit
c       should be set to the maximum order nn for which this subroutine
c       might have to be called; on subsequent calls, ninit should be
c       set to 0. PLEASE NOTE THAT THAT THESE ARRAYS USED BY THIS SUBROUTINE
c       ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES USED BY THE
c       SUBROUTINE LEGEPODD. IF these arrays have been initialized
c       by one of these two subroutines, they do not need to be
c       initialized by the other one.
c  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 real *8 elements long
c       each. Please note that these are input arrays only if ninit
c       (see above) has been set to 0; otherwise, these are output arrays.
c 
c                  Output parameters:
c 
c  val - the value at the point x of the Legendre expansion with
c       coefficients coefs (see above)
c    EXPLANATION: On exit from the subroutine, pols(1) = P_0 (x),
c       pols(2) = P_2(x), pols(3) =P_4(x),  etc.
c  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 real *8 elements long
c       each. Please note that these are output parameters only if ninit
c       (see above) has not been set to 0; otherwise, these are input
c       parameters
c 
c 
        if(ninit .eq. 0) goto 1400
        done=1
        n=0
        i=0
c 
        do 1200 nnn=2,ninit,2
c 
        n=n+2
        i=i+1
c 
        coepnm1(i)=-(5*n+7*(n*done)**2+2*(n*done)**3)
        coepnp1(i)=-(9+24*n+18*(n*done)**2+4*(n*done)**3)
        coexpnp1(i)=15+46*n+36*(n*done)**2+8*(n*done)**3
c 
        d=(2+n*done)*(3+n*done)*(1+2*n*done)
        coepnm1(i)=coepnm1(i)/d
        coepnp1(i)=coepnp1(i)/d
        coexpnp1(i)=coexpnp1(i)/d
c 
 1200 continue
c 
 1400 continue
c 
        x22=x**2
c 
        pi=x
        pip1=x*(2.5d0*x22-1.5d0)
c 
        val=coefs(1)*pi+coefs(2)*pip1
  
        do 2000 i=1,nn/2-2
c 
        pip2 = coepnm1(i)*pi +
     1      (coepnp1(i)+coexpnp1(i)*x22)*pip1
c 
        val=val+coefs(i+2)*pip2
c 
        pi=pip1
        pip1=pip2
  
  
 2000 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine legeevev(x,nn,coefs,val,ninit,
     1      coepnm1,coepnp1,coexpnp1)
        implicit real *8 (a-h,o-z)
        save
        dimension coepnm1(1),coepnp1(1),
     1            coexpnp1(1),coefs(1)
c 
c 
c       This subroutine evaluates at the point x a Legendre expansion
c       having only even-numbered elements
c 
c                  Input parameters:
c 
c  x - point on the interval [-1,1] at which the Legendre expansion
c       is to be evaluated
c  nn - order of the expansion to be evaluated
c  coefs - even-numbered coefficients of the Legendre expansion
c       to be evaluated at the point x (nn/2+2 of them things)
c  ninit - tells the subroutine whether and to what maximum order the
c       arrays coepnm1,coepnp1,coexpnp1 should be initialized.
c     EXPLANATION: The subroutine will initialize the first ninit/2+2
c                  (or so) elements of each of the arrays  coepnm1,
c       coepnp1, coexpnp1. On the first call to this subroutine, ninit
c       should be set to the maximum order nn for which this subroutine
c       might have to be called; on subsequent calls, ninit should be
c       set to 0. PLEASE NOTE THAT THAT THESE ARRAYS USED BY THIS SUBROUTINE
c       ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES USED BY THE
c       SUBROUTINE LEGEPEVEN. IF these aqrrays have been initialized
c       by one of these two subroutines, they do not need to be
c       initialized by the other one.
c  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 real *8 elements long
c       each. Please note that these are input arrays only if ninit
c       (see above) has been set to 0; otherwise, these are output arrays.
c 
c                  Output parameters:
c 
c  val - the value at the point x of the Legendre expansion with
c       coefficients coefs (see above)
c    EXPLANATION: On exit from the subroutine, pols(1) = P_0 (x),
c       pols(2) = P_2(x), pols(3) =P_4(x),  etc.
c  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 real *8 elements long
c       each. Please note that these are output parameters only if ninit
c       (see above) has not been set to 0; otherwise, these are input
c       parameters
c 
        if(ninit .eq. 0) goto 1400
c 
        done=1
        n=-1
        i=0
        do 1200 nnn=1,ninit,2
c 
        n=n+2
        i=i+1
c 
        coepnm1(i)=-(5*n+7*(n*done)**2+2*(n*done)**3)
        coepnp1(i)=-(9+24*n+18*(n*done)**2+4*(n*done)**3)
        coexpnp1(i)=15+46*n+36*(n*done)**2+8*(n*done)**3
c 
        d=(2+n*done)*(3+n*done)*(1+2*n*done)
        coepnm1(i)=coepnm1(i)/d
        coepnp1(i)=coepnp1(i)/d
        coexpnp1(i)=coexpnp1(i)/d
c 
 1200 continue
c 
 1400 continue
c 
        x22=x**2
c 
        pi=1
        pip1=1.5d0*x22-0.5d0
c 
        val=coefs(1)+coefs(2)*pip1
c 
c       n is greater than 2. conduct recursion
c 
        do 2000 i=1,nn/2-2
c 
        pip2 = coepnm1(i)*pi +
     1      (coepnp1(i)+coexpnp1(i)*x22) *pip1
        val=val+coefs(i+2)*pip2
c 
        pi=pip1
        pip1=pip2
c 
 2000 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine legepeven(x,nn,pols,ninit,
     1      coepnm1,coepnp1,coexpnp1)
        implicit real *8 (a-h,o-z)
        save
        dimension pols(1),coepnm1(1),coepnp1(1),
     1            coexpnp1(1)
c 
c       This subroutine evaluates even-numbered Legendre polynomials
c       of the argument x, up to order nn+1
c 
c                  Input parameters:
c 
c  x - the argument for which the Legendre polynomials are
c       to be evaluated
c  nn - the maximum order for which the Legendre polynomials are
c       to be evaluated
c  ninit - tells the subroutine whether and to what maximum order the
c       arrays coepnm1,coepnp1,coexpnp1 should be initialized.
c     EXPLANATION: The subroutine ill initialize the first ninit/2+2
c                  (or so) elements of each of the arrays  coepnm1,
c       coepnp1, coexpnp1. On the first call to this subroutine, ninit
c       should be set to the maximum order nn for which this subroutine
c       might have to be called; on subsequent calls, ninit should be
c       set to 0.
c  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 real *8 elements long
c       each. Please note that these are input arrays only if ninit
c       (see above) has been set to 0; otherwise, these are output arrays.
c       PLEASE NOTE THAT THAT THESE ARRAYS USED BY THIS SUBROUTINE
c       ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES USED BY THE
c       SUBROUTINE LEGEEVEV. IF these aqrrays have been initialized
c       by one of these two subroutines, they do not need to be
c       initialized by the other one.
c 
c                  Output parameters:
c 
c  pols - even-numbered Legendre polynomials of the input parameter x
c         (nn/2+2 of them things)
c    EXPLANATION: On exit from the subroutine, pols(1) = P_0 (x),
c       pols(2) = P_2(x), pols(3) =P_4(x),  etc.
c  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 real *8 elements long
c       each. Please note that these are output parameters only if ninit
c       (see above) has not been set to 0; otherwise, these are input
c       parameters. PLEASE NOTE THAT THAT THESE ARRAYS USED BY THIS SUBROUTINE
c       ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES USED BY THE
c       SUBROUTINE LEGEEVEV. If these arrays have been initialized
c       by one of these two subroutines, they do not need to be
c       initialized by the other one.
c 
c 
        if(ninit .eq. 0) goto 1400
c 
        done=1
        n=-1
        i=0
        do 1200 nnn=1,ninit,2
c 
        n=n+2
        i=i+1
c 
        coepnm1(i)=-(5*n+7*(n*done)**2+2*(n*done)**3)
        coepnp1(i)=-(9+24*n+18*(n*done)**2+4*(n*done)**3)
        coexpnp1(i)=15+46*n+36*(n*done)**2+8*(n*done)**3
c 
        d=(2+n*done)*(3+n*done)*(1+2*n*done)
        coepnm1(i)=coepnm1(i)/d
        coepnp1(i)=coepnp1(i)/d
        coexpnp1(i)=coexpnp1(i)/d
c 
 1200 continue
c 
 1400 continue
c 
        x22=x**2
c 
        pols(1)=1
        pols(2)=1.5d0*x22-0.5d0
c 
c       n is greater than 2. conduct recursion
c 
        do 2000 i=1,nn/2
c 
        pols(i+2) = coepnm1(i)*pols(i) +
     1      (coepnp1(i)+coexpnp1(i)*x22) *pols(i+1)
c 
 2000 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine legepodd(x,nn,pols,ninit,
     1      coepnm1,coepnp1,coexpnp1)
        implicit real *8 (a-h,o-z)
        save
        dimension pols(1),coepnm1(1),coepnp1(1),
     1            coexpnp1(1)
c 
c       This subroutine evaluates odd-numbered Legendre polynomials
c       of the argument x, up to order nn+1
c 
c                  Input parameters:
c 
c  x - the argument for which the Legendre polynomials are
c       to be evaluated
c  nn - the maximum order for which the Legendre polynomials are
c       to be evaluated
c  ninit - tells the subroutine whether and to what maximum order the
c       arrays coepnm1,coepnp1,coexpnp1 should be initialized.
c     EXPLANATION: The subroutine will initialize the first ninit/2+2
c                  (or so) elements of each of the arrays  coepnm1,
c       coepnp1, coexpnp1. On the first call to this subroutine, ninit
c       should be set to the maximum order nn for which this subroutine
c       might have to be called; on subsequent calls, ninit should be
c       set to 0.
c  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 real *8 elements long
c       each. Please note that these are input arrays only if ninit
c       (see above) has been set to 0; otherwise, these are output arrays.
c 
c                  Output parameters:
c 
c  pols - the odd-numbered Legendre polynomials of the input parameter x
c         (nn/2+2 of them things)
c    EXPLANATION: On exit from the subroutine, pols(1) = P_1(x),
c       pols(2) = P_3(x), pols(3) = P_5 (x), etc.
c  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 real *8 elements long
c       each. Please note that these are output parameters only if ninit
c       (see above) has not been set to 0; otherwise, these are input
c       parameters. PLEASE NOTE THAT THAT THESE ARRAYS USED BY THIS
c       SUBROUTINE ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES
c       SUSED BY THE UBROUTINE LEGEODEV. IF these arrays have been
c       initialized by one of these two subroutines, they do not need
c       to be initialized by the other one.
c 
        if(ninit .eq. 0) goto 1400
        done=1
        n=0
        i=0
c 
        do 1200 nnn=2,ninit,2
c 
        n=n+2
        i=i+1
c 
        coepnm1(i)=-(5*n+7*(n*done)**2+2*(n*done)**3)
        coepnp1(i)=-(9+24*n+18*(n*done)**2+4*(n*done)**3)
        coexpnp1(i)=15+46*n+36*(n*done)**2+8*(n*done)**3
c 
        d=(2+n*done)*(3+n*done)*(1+2*n*done)
        coepnm1(i)=coepnm1(i)/d
        coepnp1(i)=coepnp1(i)/d
        coexpnp1(i)=coexpnp1(i)/d
c 
 1200 continue
c 
 1400 continue
c 
        x22=x**2
c 
        pols(1)=x
        pols(2)=x*(2.5d0*x22-1.5d0)
c 
        do 2000 i=1,nn/2
c 
        pols(i+2) = coepnm1(i)*pols(i) +
     1      (coepnp1(i)+coexpnp1(i)*x22)*pols(i+1)
c 
 2000 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine legefdeq(x,val,der,coefs,n)
        implicit real *8 (a-h,o-z)
        save
        dimension coefs(1)
C 
C     This subroutine computes the value and the derivative
c     of a Legendre Q-expansion with coefficients coefs
C     at point X in interval (-1,1); please note that this is
c     the evil twin of the subroutine legefder, evaluating the
c     proper (P-function) Legendre expansion
c 
c                input parameters:
c 
C  X = evaluation point
C  coefs = expansion coefficients
C  N  = order of expansion
c 
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c 
c                output parameters:
c 
C     VAL = computed value
C     der = computed value of the derivative
C 
c 
        val=0
        der=0
c 
        d= log( (1+x) /(1-x) ) /2
        pkm1=d
        pk=d*x-1
c 
        pk=d
        pkp1=d*x-1
  
        derk=(1/(1+x)+1/(1-x)) /2
        derkp1=d + derk *x
c 
        val=coefs(1)*pk+coefs(2)*pkp1
        der=coefs(1)*derk+coefs(2)*derkp1
c 
c        if n=0 or n=1 - exit
c 
        if(n .ge. 2) goto 1200
c 
        if(n .eq. 0) return
c 
        return
 1200 continue
c 
c       n is greater than 2. conduct recursion
c 
        do 2000 k=1,n-1
        pkm1=pk
        pk=pkp1
c 
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
c 
        derkm1=derk
        derk=derkp1
c 
        derkp1= ( (2*k+1)*pk+(2*k+1)*x*derk - k*derkm1 )/(k+1)
c 
        val=val+coefs(k+2)*pkp1
        der=der+coefs(k+2)*derkp1
c 
 2000 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine legeqs(x,n,pols,ders)
        implicit real *8 (a-h,o-z)
        save
        dimension pols(1),ders(1)
c 
c       This subroutine calculates the values and derivatives of
c       a bunch of Legendre Q-functions at the user-specified point
c       x on the interval (-1,1)
c 
c                     Input parameters:
c 
c  x - the point on the interval [-1,1] where the Q-functions and
c       their derivatives are to be evaluated
c  n - the highest order for which the functions are to be evaluated
c 
c                     Output parameters:
c 
c  pols - the values of the Q-functions (the evil twins of the
c       Legeendre polynomials) at the point x (n+1 of them things)
c  ders - the derivatives of the Q-functions (the evil twins of the
c       Legeendre polynomials) at the point x (n+1 of them things)
c 
c 
        d= log( (1+x) /(1-x) ) /2
        pkm1=d
        pk=d*x-1
c 
        pk=d
        pkp1=d*x-1
  
        derk=(1/(1+x)+1/(1-x)) /2
        derkp1=d + derk *x
c 
c        if n=0 or n=1 - exit
c 
        if(n .ge. 2) goto 1200
        pols(1)=pk
        ders(1)=derk
        if(n .eq. 0) return
c 
        pols(2)=pkp1
        ders(2)=derkp1
        return
 1200 continue
c 
        pols(1)=pk
        pols(2)=pkp1
c 
c       n is greater than 2. conduct recursion
c 
        ders(1)=derk
        ders(2)=derkp1
c 
        do 2000 k=1,n-1
        pkm1=pk
        pk=pkp1
c 
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
        pols(k+2)=pkp1
c 
        derkm1=derk
        derk=derkp1
c 
        derkp1= ( (2*k+1)*pk+(2*k+1)*x*derk - k*derkm1 )/(k+1)
        ders(k+2)=derkp1
 2000 continue
c 
        return
        end
c 
c 
c 
c 
c 
  
  
        subroutine legeq(x,n,pol,der)
        implicit real *8 (a-h,o-z)
c 
c       This subroutine calculates the value and derivative of
c       a Legendre Q-function at the user-specified point
c       x on the interval (-1,1)
c 
c 
c                     Input parameters:
c 
c  x - the point on the interval [-1,1] where the Q-functions and
c       their derivatives are to be evaluated
c  n - the order for which the function is to be evaluated
c 
c                     Output parameters:
c 
c  pol - the value of the n-th Q-function (the evil twin of the
c       Legeendre polynomial) at the point x
c  ders - the derivatives of the Q-function at the point x
c 
c 
        save
        d= log( (1+x) /(1-x) ) /2
        pk=d
        pkp1=d*x-1
c 
c        if n=0 or n=1 - exit
c 
        if(n .ge. 2) goto 1200
        pol=d
  
        der=(1/(1+x)+1/(1-x)) /2
  
        if(n .eq. 0) return
c 
        pol=pkp1
        der=d + der *x
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
      SUBROUTINE legecFDE(X,VAL,der,PEXP,N)
      IMPLICIT REAL *8 (A-H,O-Z)
        save
      complex *16 PEXP(1),val,der
C 
C     This subroutine computes the value and the derivative
c     of a gaussian expansion with complex coefficients PEXP
C     at point X in interval [-1,1].
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
        done=1
        pjm2=1
        pjm1=x
        derjm2=0
        derjm1=1
c 
        val=pexp(1)*pjm2+pexp(2)*pjm1
        der=pexp(2)
c 
        DO 600 J = 2,N
c 
        pj= ( (2*j-1)*x*pjm1-(j-1)*pjm2 ) / j
        val=val+pexp(j+1)*pj
c 
        derj=(2*j-1)*(pjm1+x*derjm1)-(j-1)*derjm2
c 
        derj=derj/j
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
c 
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
  
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c       this is the end of the debugging code and the beginning of the
c       least squares subroutines proper
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine leastsq(a,u,w,t,n,m,ncols,rnorms,eps,v)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),u(n,1),v(m,1),w(m,1),t(1),rnorms(1)
c 
c        this subroutine decomposes the input matrix a in the
c        form
c               A=U * T * W^*                                       (1)
c 
c        with the columns of each of the matrices U,W orthonormal,
c        and T a lower triangular matrix of dimensionality (ncols,ncols)
c 
c   NOTE: this subroutine uses the subroutine leastsq1 (see) to perform
c        almost all work. However, when m < n, it tends to be much
c        more efficient than leastsq1, since it performs the two
c        Gram-Schmidt processes in the optimal order (starting with
c        the longer gram-schmidt involving shorter vectors)
c 
c 
c                     input parameters:
c 
c  a - the matrix to be decomposed (note that a is NOT damaged
c       by this subroutine in any way)
c  n,m - the dimensionalities of a
c  eps - the accuracy to which the decomposition (1) will be computed.
c 
c                     output parameters:
c 
c  u - the matrix in (1). note that the matrix u has to be dimensioned
c        u(n,m), but on exit from this subroutine only the first
c        ncols columns of  u  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit u is dimensioned
c        u(n,ncols)
c  w - the third matrix in (1). it is dimensioned w(m,ncols) on
c        exit. on entry the parameter ncols is unknown,
c        so w should be dimensioned w(m,n) by the user. however,
c        on exit only the first ncols columns of w are meaningful
c  t - the second matrix in (1)
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices u, w on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c 
c 
c                      work arrays:
c 
c  v - must be at least n*m real *8 locations long
c 
c       . . . if m > n , construct the decomposition of the matrix as
c             specified by the user
c 
cccc        if(2. ne. 3) goto 2000
        if(m .lt. n-2) goto 2000
c 
        call leastsq1(a,u,w,t,n,m,ncols,rnorms,eps,v)
cccc        call prin2('after first leastsq1, u=*',u,n*ncols)
c 
        return
c 
 2000 continue
c 
c       n is greater than m. transpose the matrix, and decompose
c       the transpose
c 
        call leastran(a,n,m,v)
        call leascopy(v,a,n*m)
c 
        call leastsq1(a,w,u,t,m,n,ncols,rnorms,eps,v)
c 
cccc        call prin2('after leastsq1, u=*',u,n*ncols)
c 
c        transpose back everything that needs to be transposed back
c 
cccc        call leastran(a,n,m,v)
        call leastran(a,m,n,v)
        call leascopy(v,a,n*m)
c 
        call leastran(t,ncols,ncols,v)
        call leascopy(v,t,ncols**2)
c 
        call leasrevr(t,ncols,ncols,v)
        call leascopy(v,t,ncols**2)
c 
        call leasrevc(u,n,ncols,v)
        call leascopy(v,u,n*ncols)
c 
c 
c 
        call leasrevc(w,m,ncols,v)
        call leascopy(v,w,m*ncols)
c 
        call leasrevc(t,ncols,ncols,v)
        call leascopy(v,t,ncols**2)
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine leastran(a,n,m,b)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),b(m,n),x(1),y(1)
c 
c       transpose a
c 
        do 1400 i=1,n
        do 1200 j=1,m
        b(j,i)=a(i,j)
 1200 continue
 1400 continue
        return
c 
c 
c 
c 
        entry leascopy(x,y,n)
c 
        do 2400 i=1,n
        y(i)=x(i)
 2400 continue
        end
c 
c 
c 
c 
c 
        subroutine leastsq2(u,w,t,n,m,ncols,y,x,work)
        implicit real *8 (a-h,o-z)
        save
        dimension u(n,ncols),w(m,ncols),t(ncols,ncols),
     1      x(m),y(n),work(ncols)
c 
c         this subroutine uses the QR-type (not exactly QR) decomposition
c         of a matrix to solve in the least squares sense the linear
c         system
c 
c                    A X = Y                                             (1)
c 
c         The expansion used must have been produced by a prior call
c         to the subroutine leastsq1 (see), and is of the form
c 
c               A=U * T * W^*                                            (2)
c 
c        with the columns of each of the matrices U,W orthonormal,
c        and T a lower triangular matrix of dimensionality (ncols,ncols)
c 
c                     input parameters:
c 
c  u - the matrix in (1). note that the matrix u has to be dimensioned
c        u(n,m), but on exit from this subroutine only the first
c        ncols columns of  u  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit u is dimensioned
c        u(n,ncols)
c  w - the third matrix in (1). it is dimensioned w(m,ncols) on
c        exit. on entry the parameter ncols is unknown,
c        so w should be dimensioned w(m,n) by the user. however,
c        on exit only the first ncols columns of w are meaningful
c  t - the second matrix in (1)
c  n,m - the dimensionalities of a
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices u, w on exit
c  y - the right-hand side in (1)
c 
c                     output parameters:
c 
c  x - the solution of the system (1) in the least squares sense
c 
c                      work arrays:
c 
c  work - must be at least ncols real *8 locations long
c 
c        . . . apply to the right-hand side the matrux  U^*
c 
        do 1400 i=1,ncols
        d=0
        do 1200 j=1,n
        d=d+u(j,i)*y(j)
 1200 continue
        work(i)=d
 1400 continue
c 
c       apply to the vector work the inverse of the matrix t
c 
        x(1)=work(1)/t(1,1)
        do 2000 i=2,ncols
c 
        d=0
        do 1600 j=1,i-1
        d=d+t(i,j)*x(j)
 1600 continue
c 
        x(i)=(work(i)-d)/t(i,i)
 2000 continue
c 
        do 2200 i=1,ncols
        work(i)=x(i)
 2200 continue
c 
c       multiply work by the matrix w
c 
        do 2600 i=1,m
        d=0
        do 2400 j=1,ncols
        d=d+w(i,j)*work(j)
 2400 continue
        x(i)=d
 2600 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine cleastsq2(u,w,t,n,m,ncols,y,x,work,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension u(n,ncols),w(m,ncols),t(ncols,ncols),whts(1)
        complex *16 x(m),y(n),work(ncols),d
c 
c         this subroutine uses the QR-type (not exactly QR) decomposition
c         of a matrix to solve in the least squares sense the linear
c         system
c 
c                    A X = Y                                             (1)
c 
c         The expansion used must have been produced by a prior call
c         to the subroutine leastsq1 (see), and is of the form
c 
c               A=U * T * W^*                                            (2)
c 
c        with the columns of each of the matrices U,W orthonormal,
c        and T a lower triangular matrix of dimensionality (ncols,ncols)
c 
c                     input parameters:
c 
c  u - the matrix in (1). note that the matrix u has to be dimensioned
c        u(n,m), but on exit from this subroutine only the first
c        ncols columns of  u  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit u is dimensioned
c        u(n,ncols)
c  w - the third matrix in (1). it is dimensioned w(m,ncols) on
c        exit. on entry the parameter ncols is unknown,
c        so w should be dimensioned w(m,n) by the user. however,
c        on exit only the first ncols columns of w are meaningful
c  t - the second matrix in (1)
c  n,m - the dimensionalities of a
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices u, w on exit
c  y - the right-hand side in (1)
c 
c                     output parameters:
c 
c  x - the solution of the system (1) in the least squares sense
c 
c                      work arrays:
c 
c  work - must be at least ncols real *8 locations long
c 
c        . . . apply to the right-hand side the matrux  U^*
c 
        do 1400 i=1,ncols
        d=0
        do 1200 j=1,n
        d=d+u(j,i)*y(j)*whts(j)
 1200 continue
        work(i)=d
 1400 continue
c 
c       apply to the vector work the inverse of the matrix t
c 
        x(1)=work(1)/t(1,1)
        do 2000 i=2,ncols
c 
        d=0
        do 1600 j=1,i-1
        d=d+t(i,j)*x(j)
 1600 continue
c 
        x(i)=(work(i)-d)/t(i,i)
 2000 continue
c 
        do 2200 i=1,ncols
        work(i)=x(i)
 2200 continue
c 
c       multiply work by the matrix w
c 
        do 2600 i=1,m
        d=0
        do 2400 j=1,ncols
        d=d+w(i,j)*work(j)
 2400 continue
        x(i)=d
 2600 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine leastsq1(a,u,w,t,n,m,ncols,rnorms,
     1    eps,v)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),u(n,1),v(m,1),w(m,1),t(1),rnorms(1)
c 
c        this subroutine decomposes the input matrix a in the
c        form
c               A=U * T * W^*                                       (1)
c 
c        with the columns of each of the matrices U,W orthonormal,
c        and T a lower triangular matrix of dimensionality (ncols,ncols)
c 
c 
c                     input parameters:
c 
c  a - the matrix to be decomposed (note that a is NOT damaged
c       by this subroutine in any way)
c  n,m - the dimensionalities of a
c  eps - the accuracy to which the decomposition (1) will be computed.
c 
c                     output parameters:
c 
c  u - the matrix in (1). note that the matrix u has to be dimensioned
c        u(n,m), but on exit from this subroutine only the first
c        ncols columns of  u  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit u is dimensioned
c        u(n,ncols)
c  w - the third matrix in (1). it is dimensioned w(m,ncols) on
c        exit. on entry the parameter ncols is unknown,
c        so w should be dimensioned w(m,n) by the user. however,
c        on exit only the first ncols columns of w are meaningful
c  t - the second matrix in (1)
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices u, w on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c 
c 
c                      work arrays:
c 
c  v - must be at least n*m real *8 locations long
c 
c        . . . using gram-schmidt process with pivoting, decompose
c              the matrix a in the form
c 
c          a=U  V^*,
c 
c        with  u an orthogonal matrix of minimum rank,
c        and v whatever it wants to be
c 
        ifpivot=1
        call leaspiv0(a,n,m,u,v,ncols,rnorms,eps,ifpivot)
cccc         call prin2('after first leaspiv0, rnorms=*',rnorms,ncols)
c 
c        using gram-schmidt process without pivoting, decompose
c        the matrix v in the form
c 
c          v=w t^*,
c 
c        with  w an orthogonal matrix of minimum rank,
c        and t a triangular matrix of dimensionality ncols * ncols
c 
        ifpivot=0
        call leaspiv0(v,m,ncols,w,t,ncols2,rnorms,eps,ifpivot)
cccc         call prin2('after second leaspiv0, rnorms=*',rnorms,ncols2)
c 
c       if the need be - restructure the matrix t
c 
        if(ncols2 .eq. ncols) return
c 
        call leastres(t,v,t,ncols,ncols2)
        ncols=ncols2
  
        return
        end
c 
c 
c 
c 
c 
        subroutine leastres(a,b,c,n1,n2)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n1,n1),b(n2,n2),c(n2,n2)
c 
        do 1400 i=1,n2
        do 1200 j=1,n2
        b(j,i)=a(j,i)
 1200 continue
 1400 continue
c 
        do 2400 i=1,n1
        do 2200 j=1,n1
        c(j,i)=b(j,i)
 2200 continue
 2400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine leaspiv0(a,n,m,b,v,ncols,rnorms,eps,ifpivot)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),b(n,m),v(m,1),rnorms(1)
c 
c       this matrix applies the compressing gram-schmidt
c       process to the matrix a, obtaining its decomposition
c       in the form
c 
c             a=b v^T,                                        (1)
c 
c       with the matrices b, v having dimensionalities
c       b(n,ncols), v(m,ncols), respectively. the reason for the
c       existence of this subroutine is the hope that the
c       dimensionality ncols, determined by this subroutine,
c       is comsiderably lower than either n or m
c 
c                     input parameters:
c 
c  a - the matrix to be decomposed (note that a is NOT damaged
c       by this subroutine in any way)
c  n,m - the dimensionalities of a
c  eps - the accuracy to which the decomposition (1) will be computed.
c 
c                     output parameters:
c 
c  b - the matrix in (1). note that the matrix b has to be dimensioned
c        b(n,m), but on exit from this subroutine only the first
c        ncols columns of  b  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit b is dimensioned
c        b(n,ncols)
c  v - the second matrix in (1). it is dimensioned v(m,ncols) on
c        exit. on entry the parameter ncols is unknown,
c        so v should be dimensioned v(m,n) by the user. however,
c        on exit only the first ncols columns of v are meaningful
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices b, v on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c 
c        . . . copy the user-supplied matrix a into b
c 
        do 1400 i=1,m
        do 1200 j=1,n
        b(j,i)=a(j,i)
 1200 continue
 1400 continue
c 
c        apply the gram-schmidt proces (with pivoting) to b
c 
         call leasgrm(b,n,m,rnorms,eps,ncols,ifpivot)
c 
c        project the original matrix on the obtained orthogonal
c        columns
c 
        do 2400 i=1,ncols
        do 2200 j=1,m
c 
        call leascapr(a(1,j),b(1,i),n,prod)
        v(j,i)=prod
 2200 continue
 2400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine leasgrm(b,n,m,rnorms,eps,ncols,ifpivot)
        implicit real *8 (a-h,o-z)
        save
        dimension b(n,m),rnorms(1)
c 
c       this subroutine applies a pivoted double gram-schmidt
c       procedure to the matrix b, resulting in a collection of
c       orthogonal vectors. the number of these vectors is
c       equal to the numerical rank of b, given the precision
c       of computations eps, which is supplied by the user.
c 
c                    input paramneters:
c 
c  b - the matrix to be gram-schmidt'ed. it is destroyed by this
c       subroutine
c  n,m - dimensionalities of the matrix b
c  eps - the machine accuracy
c 
c                     output parameters:
c 
c  b - the matrix of gram-schmidt vectors of the matrix a. note
c        that on exit from this subroutine only the first
c        ncols columns of  b  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit b is dimensioned
c        b(n,ncols)
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrix b on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c 
c       . . . prepare the array of values of norms
c             of columns
c 
        done=1
        dtot=0
        do 1400 i=1,m
        d=0
        do 1200 j=1,n
        d=d+b(j,i)**2
 1200 continue
        rnorms(i)=d
        dtot=dtot+d
 1400 continue
c 
c       . . . conduct gram-schmidt iterations
c 
        thresh=dtot*eps**2
        do 4000 i=1,m
c 
c       find the pivot
c 
         if(ifpivot .eq. 0) goto 2700
        ipivot=i
        rn=rnorms(i)
c 
        do 2200 j=i+1,m
        if(rnorms(j) .le. rn) goto 2200
        rn=rnorms(j)
        ipivot=j
 2200 continue
 2400 continue
c 
c       put the column number ipivot in the i-th place
c 
        do 2600 j=1,n
        d=b(j,i)
        b(j,i)=b(j,ipivot)
        b(j,ipivot)=d
 2600 continue
c 
        d=rnorms(ipivot)
        rnorms(ipivot)=rnorms(i)
        rnorms(i)=d
 2700 continue
c 
c       orthogonalize the i-th column to all preceeding ones
c 
        if(i .eq. 1) goto 2790
        do 2780 j=1,i-1
c 
        call leascapr(b(1,i),b(1,j),n,d)
c 
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*d
 2770 continue
 2780 continue
 2790 continue
c 
c       normalize the i-th column
c 
        call leascapr(b(1,i),b(1,i),n,d)
c 
c       if the norm of the remaining part of the matrix
c       is sufficiently small - terminate the process
c 
        if(d .lt. thresh) return
        ncols=i
c 
        d=done/dsqrt(d)
        do 2800 j=1,n
        b(j,i)=b(j,i)*d
 2800 continue
c 
c        orthogonalize everything else to it
c 
        do 3200 j=i+1,m
c 
        if(rnorms(j) .lt. thresh) goto 3200
c 
        call leascapr(b(1,i),b(1,j),n,d)
c 
         rrn=0
        do 3000 l=1,n
        b(l,j)=b(l,j)-b(l,i)*d
        rrn=rrn+b(l,j)**2
 3000 continue
        rnorms(j)=dsqrt(rrn)
 3200 continue
 3400 continue
 4000 continue
  
 4200 continue
  
  
        return
        end
c 
c 
c 
c 
c 
        subroutine leascapr(x,y,n,prod)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1)
c 
        prod=0
        do 1200 i=1,n
        prod=prod+x(i)*y(i)
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine leasrevc(a,n,m,b)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),b(n,m)
c 
c       reverse the columns of a
c 
        do 1400 i=1,m
        do 1200 j=1,n
        b(j,i)=a(j,m-i+1)
 1200 continue
 1400 continue
        return
c 
c 
c 
c 
        entry leasrevr(a,n,m,b)
c 
c       reverse the rows of a
c 
        do 2400 i=1,m
        do 2200 j=1,n
        b(j,i)=a(n-j+1,i)
 2200 continue
 2400 continue
        return
        end
