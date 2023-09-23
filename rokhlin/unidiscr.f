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
