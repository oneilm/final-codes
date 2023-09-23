        implicit real *8 (a-h,o-z)
        complex *16 ima,rk
        dimension thetas(2000),phis(2000),whts(2000)
         real *8 vectors(3,80 000),
     1       sources(3,1 000),receiv(3,1 000)
cccc         complex *16 hexp(80 000),trhj(80 000),jexp(80 000),
cccc     1       fldcoes(1600 000)
c 
         complex *16 hexp(40 000),trhj(40 000),jexp(40 000),
     1       fldcoes(800 000)
cccc     1       fldcoes(200 000)
c 
        data ima/(0.0d0,1.0d0)/
c 
        CALL PRINI(6,13)
C 
c       . . . construct the testing parameters
c 
c 
         PRINT *, 'ENTER r'
         READ *,r
         CALL PRIn2('r=*',r,1)
c 
         PRINT *, 'ENTER rho'
         READ *,rho
         CALL PRIn2('rho=*',rho,1)
c 
        PRINT *, 'ENTER -log(eps) base 10'
        READ *,db
        CALL PRIN2('db=*',db,1 )
c 
        eps=0.1**db
        call prin2('and eps=*',eps,1)
c 
c        find the number of bessel functions of rk*r that are
c        greater than eps
c 
c 
        rk=1
c 
        coef=1.44
        coef=1.5
cccc        coef=1.732
cccc        coef=2
c 
c        construct the geometry of the test
c 
c        . . . the center of the h-expansion
c 
        xh=0
        yh=0
        zh=0
c 
c        . . . the position of the center of the j-expansion
c 
       xj=xh-rho
ccc       xj=xh
       yj=yh
cccc       zj=zh+rho
       zj=zh
c 
 1200 format(2x,'xh= ',e11.5,' yh= ',e11.5,' zh= ',e11.5)
 1400 format(10x)
         write(6,1400)
         write(13,1400)
         write(6,1200) xh,yh,zh
         write(13,1200) xh,yh,zh
c 
 2400 format(2x,'xj= ',e11.5,' yj= ',e11.5,' zj= ',e11.5)
         write(6,1400)
         write(13,1400)
         write(6,2400) xj,yj,zj
         write(13,2400) xj,yj,zj
c 
c        conduct massive testing
c 
        call getpara(x0,y0,z0,x,y,z,r,rk,eps,coef,
     1    vectors,phis,thetas,whts,nphi,ntheta,nterms,hexp)
  
  
  
        ntest=4
        call masstest2(nphi,ntheta,nterms,sources,receiv,
     1      xh,yh,zh,xj,yj,zj,r,rho,ntest,vectors,rk,
     2      hexp,jexp,trhj,eps,coef,fldcoes,whts)
  
         stop
         end
c 
c 
c 
c 
c 
        subroutine masstest2(nphi,ntheta,nterms,sources,receiv,
     1      xh,yh,zh,xj,yj,zj,r,rho,ntest,vectors,rk,
     2      hexp,jexp,trhj,eps,coef,fldcoes,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension sources(3,1),receiv(3,1),vectors(3,1),
     1      thtest(100),phitest(100),whts(1),errs(20 000)
        complex *16 trhj(1),rk,hexp(1),f1,f2,coetot(1000)
c 
        complex *16 jexp(1),fldcoes(ntheta*nphi,1)
c 
c        construct the sphere and the whole structure on it
c 
cccc             call prinf('in masstest2, nphi=*',nphi,1)
cccc             call prinf('in masstest2, ntheta=*',ntheta,1)
c 
c       construct the arrays of sources and receivers
c 
        call prin2('in masstest2, xh=*',xh,1)
        call prin2('in masstest2, yh=*',yh,1)
        call prin2('in masstest2, zh=*',zh,1)
  
        call prin2('in masstest2, xj=*',xj,1)
        call prin2('in masstest2, yj=*',yj,1)
        call prin2('in masstest2, zj=*',zj,1)
c 
        call prin2('in masstest2, rk=*',rk,1)
c 
        done=1
        pi=datan(done)*4
        htheta=pi/ntest
        hphi=2*pi/ntest
c 
        do 1200 i=1,ntest
        thtest(i)=(i-1)*htheta + htheta/3
        phitest(i)=(i-1)*hphi
 1200 continue
        call prin2('phitest=*',phitest,ntest)
        call prin2('thtest=*',thtest,ntest)
c 
        ii=0
        do 1600 i=1,ntest
        do 1400 j=1,ntest
c 
        ii=ii+1
        sources(1,ii)=dsin(thtest(i))*dcos(phitest(j)) *r +xh
        sources(2,ii)=dsin(thtest(i))*dsin(phitest(j)) *r   +yh
        sources(3,ii)=dcos(thtest(i))               *r   +zh
c 
        receiv(1,ii)=dsin(thtest(i))*dcos(phitest(j)) *r +xj
        receiv(2,ii)=dsin(thtest(i))*dsin(phitest(j)) *r   +yj
        receiv(3,ii)=dcos(thtest(i))               *r   +zj
c 
 1400 continue
 1600 continue
c 
         call prin2('sources as created*',sources,ntest**2*3)
         call prin2('receiv as created*',receiv,ntest**2*3)
c 
c       construct the translation operator connecting
c       the h-expansion with the center at (xh,yh,zh) with
c       the j-expansion with the center at (xj,yj,zj)
c 
        call hjtran30(ier,xh,yh,zh,xj,yj,zj,
     1    vectors,nphi,ntheta,rk,trhj,r,eps,
     2    b,prop,nterms,fldcoes,coetot)
c 
        call prinf('after hjtran30, ier=*',ier,1)
        call prin2('after hjtran30, coetot=*',
     1      coetot,ntheta*2)
c 
        call prin2('after hjtran30, b=*',b,1)
        call prin2('after hjtran30, prop=*',prop,1)
  
c 
c        one receiver after another, construct the coefficients of
c        the linear forms converting the j-expansion into the values
c        of the potential at the receivers. store the latter in the
c        array fldcoes
c 
        do 2600 i=1,ntest**2
c 
c       . . . construct the set of coefficients
c 
        call jexeva0 ( xj,yj,zj,receiv(1,i),receiv(2,i),receiv(3,i),
     1      vectors,rk,nphi,ntheta,whts,fldcoes(1,i)  )
c 
 2600 continue
c 
c       for all pairs of sources and receivers, construct their
c       interactions, both via the expansions and directly;
c       calculate the differences
c 
        ii=0
        errmax=0
c 
        do 3000 i=1,ntest**2
c 
        call prinf('in masstest in the second loop, i=*',i,1)
c 
c       . . . expand the source into an h-expansion
c 
        call hexpan(xh,yh,zh,sources(1,i),sources(2,i),
     1      sources(3,i),vectors,rk,nphi,ntheta,hexp)
c 
cccc         call prin2('hexp as created*',hexp,ntheta*nphi*2)
c 
c        convert the h-expansion into a j-expansion
c 
        call hjtrans(trhj,hexp,nphi,ntheta,jexp )
c 
c        one receiver after another, evaluate the j-expansion at them
c        and compare the result with that for the direct calculation
c 
        do 2800 j=1,ntest**2
c 
         ii=ii+1
c 
c       . . . evaluate the expansion
c 
cccc        call jexeva1(fldcoes(1,j),jexp,nphi,ntheta,f1)
        call scawrong(fldcoes(1,j),jexp,nphi,ntheta,f1)
c 
c       . . . evaluate the potential directly
c 
        call fldchg(sources(1,i),sources(2,i),sources(3,i),
     1     receiv(1,j),receiv(2,j),receiv(3,j),rk,f2)
c 
        errs(ii)=cdabs( (f2-f1) /f2)
  
        if(errmax .lt. errs(ii)) errmax=errs(ii)
c 
 2800 continue
         call prin2('and errsmax=*',errmax,1)
 3000 continue
  
         call prin2('and relative errors are*',errs,ii)
         call prin2('and maximum relative error is*',errmax,1)
  
         call prinf('remember, nphi=*',nphi,1)
         call prinf('remember, ntheta=*',ntheta,1)
         call prinf('and nphi times ntheta=*',nphi*ntheta,1)
         call prin2(
     1     'again, proportion of sphere where operator is non-zero*',
     2     prop,1)
  
  
  
        return
        end
c 
c 
c 
c 
c 
        subroutine getpara(x0,y0,z0,x,y,z,r,rk,eps,coef,
     1    vectors,phis,thetas,whts,nphi,ntheta,nterms,w)
        implicit real *8 (a-h,o-z)
        save
        complex *16 rk
        real *8 vectors(3,nphi,ntheta),w(1),thetas(1),phis(1)
c 
c       this subroutine evaluates the parameters to be used by the
c       subroutine hjtran30 in the construction of the beam-like
c       translation operator. It also constructs the array vectors
c       containing the exterior normals to the spehere at the nodes
c       of the standard (nphi x theta) - discretization (incidentally,
c       these are also the nodes of the standard discretization of
c       the unit sphere itself), and several other parameters.
  
c               input parameters:
c 
c  x0,y0,z0 - coordinates of the center of the h-expansion
c  x,y,z  - coordinates of the center of the j-expansion
c  r - the radius of the group whose potential is to be translated
c  rk - the helmholtz coefficient
c  coef - the oversampling coefficient to be used for the construction
c        of the beam-like translation operator; so far, values between
c        1.44 and 2.2 have been decent
c 
c               output parameters:
c 
c  vectors - the normals to the unity sphere at the standard nodes,
c       as created by subroutine sphecre; also nodes of the standard
c       (nphi x ntheta) discretization of the unit sphere
c  phis - the value of phi at the nodes on each parallel
c        (nphi of them)
c  thetas - the value of theta at the nodes on each meridian
c        (ntheta of them)
c  whts - the weights to be used by the subroutine spheint,
c         (see) to integrate functions on the sphere.
c  nphi - the number of elements in the discretization
c        of each parallel
c  ntheta - the number of elements in the discretization
c        of each meridian; this is also the length of the Legendre
c        expansion of the translation operator as a function of cos(theta)
c  nterms - the order of the spherical expansion to be used for the
c        construction of the transfer function {\bf in the absence of
c        oversampling and other measures for the construction of the
c        beam-like translation operator}
c 
c                work arrays:
c 
c  w - must be at least ????
c 
c        . . . find the numbers of nodes required for
c              discretization of the sphere
c 
        call jfuns3d(rk*r,eps/10,w,njs)
c 
        call prinf(
     1  'in getpara, number of Bessel functions greater than eps*',
     2      njs,1)
c 
        nphi=njs*4
        ntheta=njs*2
        nterms=njs*2
c 
        nphi=nphi*coef
        ntheta=ntheta*coef
c 
        m=ntheta-nterms
c 
        call prinf('nphi as created in getpara*',nphi,1)
        call prinf('ntheta as created in getpara*',ntheta,1)
        call prinf('nterms as created in getpara*',nterms,1)
c 
        call sphecre(nphi,ntheta,phis,thetas,whts,
     1         vectors)
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine hjtran30(ier,x0,y0,z0,x,y,z,
     1    vectors,nphi,ntheta,rk,trans,r,
     2      eps,b,prop,nterms,ww,coetot)
        implicit real *8 (a-h,o-z)
        save
        complex *16 rk,trans(nphi,ntheta),coetot(1)
        real *8 vectors(3,nphi,ntheta),ww(1)
c 
c       this subroutine constructs the transfer coefficients
c       for the conversion of a far-field representation of
c       an h-expansion onto a far-field representation of
c       a j-expansion around a sufficiently shifted center
c 
c               input parameters:
c 
c  x0,y0,z0 - coordinates of the center of the h-expansion
c  x,y,z  - coordinates of the center of the j-expansion
c  vectors - the normals to the unity sphere at the
c           standard nodes, as created by subroutine sphecre
c  nphi - the number of elements in the discretization
c        of each parallel
c  ntheta - the number of nodes in the discretization
c        of each meridian; this is also the length of the Legendre
c        expansion of the translation operator as a function of cos(theta)
c  rk - the helmholtz coefficient
c  r - the radius of the group whose potential is to be translated
c  nterms - the order of the spherical expansion to be used for the
c        construction of the transfer function {\bf in the absence of
c        oversampling and other measures for the construction of the
c        beam-like translation operator}
c 
c   IMOPORTANT NOTE: It is expected that the parameters nphi, ntheta,
c                    nterms, vectors have been obtained as an output of
c                    the subroutine getpara (see).
c 
c               output parameters:
c 
c  ier - error return code; ier=0 means successful conclusion
c                           ier=8 means that bisection couls not be
c          initialized because no value of the parameter b could be
c          found for which the operator was supported by less than
c          0.9 of the area of the sphere. This means trouble; a possible
c          excplanation is that eps is too small for the geometry; this
c          also has been observed for eps > 0.001, for which this
c          subroutine has not been tuned.
c  trans - the beam-like transfer function tabulated at nphi* ntheta
c        points on the sphere
c  b - the parameter indicating the size of the part of the sphere
c        where the translation operator is large. Specifically, the
c        operator should be considered large (or at least, suspected
c        of being large) for all theta such that cos(theta) > b, with
c        theta the angle between the direction of translation and the
c        point on the sphere
c  prop - the proportion of the sphere where the translation operator
c        is large (obvoiously, it has to be between 0 and 1)
c  coetot - the coefficients of the Legendre expansion of the
c        translation operator (as a function of cos(theta) ).
c 
c                work arrays:
c 
c  ww - must be sufficiently large
c 
c        . . . find the distance between the centers
c              of the expansions, the shift vector
c 
        dx=x-x0
        dy=y-y0
        dz=z-z0
        d=dx**2+dy**2+dz**2
        d=dsqrt(d)
        rho=d
        dx=dx/d
        dy=dy/d
        dz=dz/d
c 
c       . . . construct the Legendre expansion of the translation
c             operator
c 
        m=ntheta-nterms
        npts=nterms+m+4
c 
        call prinf('in hjtran30, nphi=*',nphi,1)
        call prinf('in hjtran30, ntheta=*',ntheta,1)
        call prinf('in hjtran30, nterms=*',nterms,1)
        call prinf('in hjtran30, m=*',m,1)
        call prinf('in hjtran30, npts=*',npts,1)
c 
cccc        call prinf('before projcon2, nterms+m=*',nterms+m,1)
cccc        call prinf('and npts x (nterms+m)=*',npts*(nterms+m),1)
  
        call projcon2(ier,nterms,m,npts,
     1      eps,rho,rk,r,coetot,ww,b,prop)
c 
        call prinf('after projcon2, ier=*',ier,1)
        call prin2('after projcon2, b=*',b,1)
        call prin2('after projcon2, prop=*',prop,1)
c 
        if(ier .ne. 0) return
c 
c        one node after another, construct the elements
c        of the diagonal operator converting the user-supplied
c        h-expansion into a j-expansion (both in far-field
c        form)
c 
        do 3400 i=1,ntheta
        do 3200 j=1,nphi
c 
        call onetran3(vectors(1,j,i),rk,ntheta,trans(j,i),
     1        ww(2),ww(n+20),dx,dy,dz,coetot)
c 
 3200 continue
 3400 continue
c 
c        count the elements of the translation operator that are
c        greater than eps
c 
        eps2=eps**2
        eps2=eps2 *1000
        d=cdabs(rho*rk)
        eps2=eps2*25
        eps2=eps2 / d**2
c 
        nonzero=0
        do 3800 i=1,ntheta
        do 3600 j=1,nphi
        dd=trans(j,i)*dconjg(trans(j,i))
        if(dd .gt. eps2) nonzero=nonzero+1
        if(dd .lt. eps2) trans(j,i)=0
 3600 continue
 3800 continue
         call prinf('exiting hjtran30, nonzero=*',nonzero,1)
         call prinf('out of *',nphi*ntheta,1)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine projcon2(ier,nterms,m,npts,
     1      eps,rho,rk,r,coetot,ww,b,prop)
        implicit real *8 (a-h,o-z)
        save
        dimension coetot(1),prot(6,100),ww(1)
c 
c        This subroutine uses an ugly version of the bisection
c        procedure to find the (somewhat) optimal value of the
c        parameter b resulting in the sparsest possible translation
c        operator that still provides accuracy eps (or so). The
c        operator converts an h-expansion of a group of charges of
c        radius r into a j-expansion around a center removed
c        by rho from the center of the h-expansion.
c 
c        For each value of b, the subroutine uses the least squares
c        procedure (implemented by the subroutine projcon (see))
c        to design the translation operator on one meridian; the
c        operator is a function of cos(theta), with
c        cos(theta) \in [-1,1]. The operator is chosen in such a
c        way that it performes the translation with accuracy eps,
c        and is minimized on the interva cos(theta) \in [-1,b].
c 
c        The subroutine also produces the coefficients coetot
c        of the Legendre expansion of the sparse translation operator,
c        and several other parameters (see below)
c 
c                            input parameters:
c 
c  nterms - the order of the spherical expansion to be used for the
c        construction of the transfer function {\bf in the absence of
c        oversampling and other measures for the construction of the
c        beam-like translation operator}
c  m - the number of oversampling nodes along each meridian
c  npts - the number of test points on the interval [-1,b] at
c        which the attempt is to be made to make the operator small
c  eps - the accuracy to which the subroutine will attempt to construct
c        the translation operator
c  rho - the distance between the centers of h- and j-expansions
c  rk - the helmholtz coefficient
c  r - the radius of the group whose potential is to be translated
c 
c               output parameters:
c 
c  ier - error return code; ier=0 means successful conclusion
c                           ier=8 means that bisection couls not be
c          initialized because no value of the parameter b could be
c          found for which the operator was supported by less than
c          0.9 of the area of the sphere. This means trouble; a possible
c          excplanation is that eps is too small for the geometry; this
c          also has been observed for eps > 0.001, for which this
c          subroutine has not been tuned.
c  trans - the beam-like transfer function
c  b - the parameter indicating the size of the part of the sphere
c        where the translation operator is large. Specifically, the
c        operator should be considered large (or at least, suspected
c        of being large) for all theta such that cos(theta) > b, with
c        theta the angle between the direction of translation and the
c        point on the sphere
c  prop - the proportion of the sphere where the translation operator
c        is large (obvoiously, it has to be between 0 and 1)
c  coetot - the coefficients of the Legendre expansion of the
c        translation operator (as a function of cos(theta).
c 
c                      work arrays:
c 
c  ww - must be sufficiently large
c 
c 
c 
c       . . . initialize the bisection process for the finding
c             {\bf some} b for which the operator is supported on less
c             than 100% of the sphere
c 
        ier=0
        b1=-1
        f1=1
        b3=1
        f3=1
c 
c       conduct the bisection
c 
        numit=14
        thresh=0.9
        ifplot=0
        ifprint=0
        do 1200 i=1,numit
c 
cccc        call prinf('in the first loop in projcon2, i=*',i,1)
c 
c       construct the point b2 and the value of the function at it
c 
        b2=(b1+b3)/2
c 
        call projcon(ww,nterms,m,npts,b2,eps,
     1      rho,rk,r,coetot,f2,ifplot,ifprint)
c 
        if(f2. lt. thresh) goto 1400
c 
        b3=b2
        f3=f2
 1200 continue
c 
        ier=8
        return
 1400 continue
c 
c       conduct the bisection to find the point b at which
c       the proportion of the sphere supporting the operator
c       is minimized
c 
        do 2000 i=1,numit
c 
        call prinf('in the second  loop in projcon2, i=*',i,1)
c 
        if(i .eq. numit) ifplot=1
        if(i .eq. numit) ifprint=1
c 
        prot(1,i)=b1
        prot(2,i)=b3
        prot(3,i)=f1
        prot(4,i)=f3
        prot(5,i)=b3-b1
        prot(6,i)=f3-f1
c 
       call bisect(b1,b2,b3,f1,f2,f3,
     1      nterms,m,npts,eps,
     2      rho,rk,r,coetot,ifplot,ifprint,ww)
c 
 2000 continue
c 
        b=b2
        prop=f2
c 
        call prin2('and protocol of bisection is*',prot,0)
c 
 2100 format(7x,'b1',10x,'b3',11x,'f1',10x,'f3',10x,'b3-b1',
     1    5x,'f3-f1')
        write(6,2100)
        write(13,2100)
  
 2200 format(2x,i3,2x,e10.4,2x,E10.4,2x,e9.4,2x,e9.4,2x,e9.4,2x,e10.4)
        write(6,2200) (i,(prot(l,i),l=1,6),i=1,numit)
        write(13,2200) (i,(prot(l,i),l=1,6),i=1,numit)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine bisect(b1,b2,b3,f1,f2,f3,
     1      nterms,m,npts,eps,rho,rk,r,coetot,
     2      ifplot,ifprint,ww)
        implicit real *8 (a-h,o-z)
        save
        dimension coetot(1),ww(1)
c 
c        this subroutine conducts a single step of the bisection
c        procedure to find the (somewhat) optimal value of the
c        parameter b resulting in the sparsest possible translation
c        operator that still provides accuracy eps (or so).
c     IMPORTANT NOTE. This procedure attempts to find the minimum of
c        a convex function, not a zero of a function. Thus, it is a
c        somewhat weird scheme; its single iteration reduces the size
c        of the onterval on which the minimum is localized by a factor
c        of sqrt(2) (!!) - some bisection!. However, the function
c        being minimized is convex, but otherwise fairly ill-behaved
c        - it is somewhat difficult to construct a robust Newton-style
c        scheme for it.
c 
c                            input parameters:
c 
c  b1,b2,b3 - the current three points defining the interval where the
c        minimum is located. Specifically, the minimum must be located
c        on the interval (b1, b3); b2 is some point on the interval
c        (b1, b3) such thta f(b2) < f(b1) and f(b2) < f(b3).
c  f1,f2,f3 - the values of the function being minimized at the points
c        b1,b2,b3 respectively.
c  nterms - the order of the spherical expansion to be used for the
c        construction of the transfer function {\bf in the absence of
c        oversampling and other measures for the construction of the
c        beam-like translation operator}
c  m - the number of oversampling nodes along each meridian
c  npts - the number of test points on the interval [-1,b] at
c        which the attempt is to be made to make the operator small
c  eps - the accuracy to which the subroutine will attempt to construct
c        the translation operator
c  rho - the distance between the centers of h- and j-expansions
c  rk - the helmholtz coefficient
c  r - the radius of the group whose potential is to be translated
c  ifplot - the parameter telling the subroutine whether it should
c        plot the  real and imaginary parts of the obtained translation
c        operator on the FORTRAN files number 21, 22 respectively.
c        ifplot=0 causes the subroutine to suppress the plotting.
c  ifprint - the parameter telling the subroutine whether it should
c        print various types of internal information on the FORTRAN
c        files with numbers 6, 13. ifprint=0 causes the subroutine to
c        suppress the printing
c 
c               output parameters:
c 
c  b1,b2,b3 - the improved three points defining the interval where the
c        minimum is located. Specifically, the minimum must be located
c        on the interval (b1, b3); b2 is some point on the interval
c        (b1, b3) such thta f(b2) < f(b1) and f(b2) < f(b3).
c  f1,f2,f3 - the values of the function being minimized at the points
c        b1,b2,b3 respectively.
c  coetot - the coefficients of the Legendre expansion of the obtained
c        translation operator (as a function of cos(theta).
c 
c                      work arrays:
c 
c  ww - must be sufficiently large
c 
c 
c       . . . find out which of the two intervals is the longer
c 
        d12=b2-b1
        d23=b3-b2
c 
c       if the first interval is the longer one - act accordingly
c 
        if(d23 .gt. d12) goto 2000
c 
        c=(b2+b1)/2
c 
        call projcon(ww,nterms,m,npts,c,eps,
     1      rho,rk,r,coetot,g,ifplot,ifprint)
c 
        if(g .gt. f2) goto 1400
c 
        f3=f2
        b3=b2
        b2=c
        f2=g
        return
 1400 continue
c 
        b1=c
        f1=g
        return
c 
 2000 continue
c 
c        the second interval is the longer one - act accordingly
c 
        c=(b3+b2)/2
c 
        tt1=second(IT1)
        call projcon(ww,nterms,m,npts,c,eps,
     1      rho,rk,r,coetot,g,ifplot,ifprint)
c 
        tt2=second(IT1)
c 
        call prin2('cpu time for projcon0 is*',tt2-tt1,1)
c 
        If(g .gt. f2) goto 2400
c 
        f1=f2
        b1=b2
        b2=c
        f2=g
        return
c 
 2400 continue
c 
        b3=c
        f3=g
        return
        end
c 
c 
c 
c 
c 
        subroutine projcon(w,nterms,m,npts,b,eps,
     1      rho,rk,r,coetot,prop,ifplot,ifprint)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1)
        complex *16 rk,coetot(1)
c 
c        This subroutine constructs the diagonal form of the
c        operator converting an h-expansion of a group of charges of
c        radius r into a j-expansion around a center removed
c        by rho from the center of the h-expansion.
c 
c        Depending on the value of b, the subroutine designes
c        the translation operator on one meridian; the operator
c        is a function of cos(theta), with cos(theta) \in [-1,1].
c        The operator is chosen in such a way that it performes
c        the translation with accuracy eps, and is minimized on
c        the interval cos(theta) \in [-1,b].
c 
c        The subroutine also produces the coefficients coetot
c        of the Legendre expansion of the sparse translation operator,
c        and several other parameters (see below). Actually, this is
c        simply a prinitive memory manager for the subroutine projcon0
c        (see) which does all the work.
c 
c    IMPORTANT NOTE: This subroutine does not optimize the value of
c        the parameter b, which should be chosen depending on r,
c        rho, rk, and eps. The optimization of be is performed by
c        the subroutine projcon2 (see).
c 
c                            input parameters:
c 
c  n - the order of the spherical expansion to be used for the
c        construction of the transfer function {\bf in the absence of
c        oversampling and other measures for the construction of the
c        beam-like translation operator}; IN MOST OTHER SUBROUTINES
C        THIS PARAMETER IS DENOTED BY NTERMS.
c  m - the number of oversampling nodes along each meridian
c  npts - the number of test points on the interval [-1,b] at
c        which the attempt is to be made to make the operator small;
c        should be equal to n+m+4
c  b - the parameter indicating the size of the part of the sphere
c        where the translation operator is large. Specifically, the
c        operator should be considered large (or at least, suspected
c        of being large) for all theta such that cos(theta) > b, with
c        theta the angle between the direction of translation and the
c        point on the sphere
c  eps - the accuracy to which the subroutine will attempt to construct
c        the translation operator
c  rho - the distance between the centers of h- and j-expansions
c  rk - the helmholtz coefficient
c  r - the radius of the group whose potential is to be translated
c  ifplot - the parameter telling the subroutine whether it should
c        plot the  real and imaginary parts of the obtained translation
c        operator on the FORTRAN files number 21, 22 respectively.
c        ifplot=0 causes the subroutine to suppress the plotting.
c  ifprint - the parameter telling the subroutine whether it should
c        print various types of internal information on the FORTRAN
c        files with numbers 6, 13. ifprint=0 causes the subroutine to
c        suppress the printing
c 
c               output parameters:
c 
c  prop - the proportion of the sphere where the translation operator
c        is large (obvoiously, it has to be between 0 and 1)
c  coetot - the coefficients of the Legendre expansion of the
c        translation operator (as a function of cos(theta).
c 
c                      work arrays:
c 
c     w - must be at least 4 * (npts+nterms+m) * (nterms+m) +
c             + (nterms+m)**2 + 50 real *8 locations long
c 
c        . . . allocate memory for the construction of the
c              translation operator
c 
        ivals=1
        lvals=(npts+nterms+m)*(nterms+m) + 10
c 
        iu=ivals+lvals
        lu=(npts+nterms+m)*(nterms+m) + 10
c 
        iv=iu+lu
        lv=(npts+nterms+m)*(nterms+m) + 10
c 
        iw=iv+lv
        lw=(npts+nterms+m)*(nterms+m) + 10
c 
        itria=iw+lw
        ltria=(nterms+m)**2 + 10
c 
c        construct the translation operator
c 
        call projcon0(w(ivals),nterms,m,npts,w(iu),w(iv),b,
     1      eps,rho,rk,r,coetot,prop,ifplot,ifprint,
     2      w(iw),w(itria) )
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine projcon0(vals,n,m,npts,u,v,b,eps,
     1      rho,rk,r,coetot,prop,ifplot,ifprint,w,tria)
        implicit real *8 (a-h,o-z)
        save
        dimension vals(npts+n+m,n+m),t(2000),whts(2000),
     1     u(npts+n+m,1),v(n+m,1),w(1),p(2000),t00(2000),
     2     ttest(8000),
     3     frea(8000),fima(8000),rea(2),tria(1)
c 
        complex *16 coefs(2000),rhs(2000),work(2 000),
     1      coetot(1),sol(2000),com,f(8000),
     2      f2(8000),hs(2000),ima,rk,cd
c 
        equivalence (rea(1),com),(hs(1),rhs(1)),(hs(1),p(1))
        data ima/(0.0d0,1.0d0)/,ifinit/1/
c 
c        This subroutine constructs the diagonal form of the
c        operator converting an h-expansion of a group of charges of
c        radius r into a j-expansion around a center removed
c        by rho from the center of the h-expansion.
c 
c        Depending on the value of b, the subroutine designes
c        the translation operator on one meridian; the operator
c        is a function of cos(theta), with cos(theta) \in [-1,1].
c        The operator is chosen in such a way that it performes
c        the translation with accuracy eps, and is minimized on
c        the interval cos(theta) \in [-1,b].
c 
c        The subroutine also produces the coefficients coetot
c        of the Legendre expansion of the sparse translation operator,
c        and several other parameters (see below).
c 
c    IMPORTANT NOTE: This subroutine does not optimize the value of
c        the parameter b, which should be chosen depending on r,
c        rho, rk, and eps. The optimization of be is performed by
c        the subroutine projcon2 (see)
c 
c                            input parameters:
c 
c  n - the order of the spherical expansion to be used for the
c        construction of the transfer function {\bf in the absence of
c        oversampling and other measures for the construction of the
c        beam-like translation operator}; IN MOST OTHER SUBROUTINES
C        THIS PARAMETER IS DENOTED BY NTERMS.
c  m - the number of oversampling nodes along each meridian
c  npts - the number of test points on the interval [-1,b] at
c        which the attempt is to be made to make the operator small;
c        should be equal to n+m+4
c  b - the parameter indicating the size of the part of the sphere
c        where the translation operator is large. Specifically, the
c        operator should be considered large (or at least, suspected
c        of being large) for all theta such that cos(theta) > b, with
c        theta the angle between the direction of translation and the
c        point on the sphere
c  eps - the accuracy to which the subroutine will attempt to construct
c        the translation operator
c  rho - the distance between the centers of h- and j-expansions
c  rk - the helmholtz coefficient
c  r - the radius of the group whose potential is to be translated
c  ifplot - the parameter telling the subroutine whether it should
c        plot the  real and imaginary parts of the obtained translation
c        operator on the FORTRAN files number 21, 22 respectively.
c        ifplot=0 causes the subroutine to suppress the plotting.
c  ifprint - the parameter telling the subroutine whether it should
c        print various types of internal information on the FORTRAN
c        files with numbers 6, 13. ifprint=0 causes the subroutine to
c        suppress the printing
c 
c               output parameters:
c 
c  prop - the proportion of the sphere where the translation operator
c        is large (obvoiously, it has to be between 0 and 1)
c  coetot - the coefficients of the Legendre expansion of the
c        translation operator (as a function of cos(theta).
c 
c                      work arrays:
c 
c  vals,u,v,w - must be at least (npts+nterms+m) * (nterms+m) +10
c        real *8 locations each
c  ltria - must be at least (nterms+m)**2 + 10 real *8 locations each
c 
        if(ifprint .ne. 0) call prinf('in projcon0, n=*',n,1)
        if(ifprint .ne. 0) call prinf('in projcon0, m=*',m,1)
c 
c        initialize the matrix val
c 
        do 1100 i=1,n+m
        do 1050 j=1,npts+n+m
        vals(j,i)=0
 1050 continue
 1100 continue
c 
c       construct the Gaussian nodes on the interval [-1,b]
c 
cccc         call prinf('before first gauswhts, ifinit=*',ifinit,1)
        ifwhts=1
        if(ifinit .ne. 0) call gauswhts(npts,t00,whts,ifwhts)
c 
        alpha=(b+1)/2
        beta=(b-1)/2
        do 1200 i=1,npts
        t(i)=t00(i)*alpha+beta
 1200 continue
c 
c       fill in the penalty part of the weights
c 
        do 1300 i=npts+1,npts+n+m
        whts(i)=1
 1300 continue
c 
c        evaluate the appropriate Legendre polynomials
c        at these nodes
c 
        do 1600 i=1,npts
c 
        call legpols(t(i),n+m+1,p(2))
c 
        d=dsqrt(whts(i))
        do 1400 j=1,n+m
        vals(i,j)=p(j)*d
 1400 continue
 1600 continue
c 
c        construct the penalty part of the matrix
c 
        eps3=1.0d-16
        call jfuns3d(rk*r*2,eps3,hs(2),npens)
c 
        if(ifprint .ne. 0)
     1      call prinf('constructing penalty, npens=*',npens,1)
c 
        call hfuns3d(rk*r*2,npens,hs(2) )
c 
        pen=1
c 
        do 1650 i=1,n+m
c 
        if(i .le. npens) vals(i+npts,i)=pen/cdabs(hs(i))
 1650 continue
c 
        eps2=1.0d-12
        lw=100 000
c 
cccc        call leastsq1(vals,u,v,tria,npts+n+m,n+m,ncols,rnorms,
        call leastsq1(vals,u,v,tria,npts+n+m,n+m,ncols,sol,
     1    eps2,w)
c 
         if(ifprint .eq. 0) goto 1670
         call prinf(' after leatssq1, ncols=*',ncols,1)
         call prin2(' after leastsq1, rnorms=*',sol,ncols)
 1670 continue
c 
c        scale the elements of the matrix u to compensate
c        for the weights
c 
        do 1800 i=1,npts
        d=1/dsqrt(whts(i))
        do 1700 j=1,m+n
        u(i,j)=u(i,j)*d
 1700 continue
 1800 continue
c 
c        construct the Legendre expansion of the right-hand side
c 
        do 2100 i=1,n+6
        coefs(i)=0
 2100 continue
c 
c       . . . evaluate the hankel functions of the distance
c 
        call hfuns3d(rk*rho,n+4,hs(2) )
c 
cccc        call prin2('after hfuns3d, hs=*',hs,n*2)
c 
        do 2200 i=0,n
c 
        coefs(i+1)=2*i+1
c 
        coefs(i+1)=(2*i+1)*hs(i+1)*ima**i
 2200 continue
c 
c       evaluate the expansion at the Gaussian nodes on
c       the interval [-1,b]
c 
        do 2400 i=1,npts
c 
        call clegeexe(t(i),rhs(i),coefs,n-1)
 2400 continue
c 
c        fill in the penalty part of the right-hand side
c 
        do 2500 i=npts+1,npts+n+m
        rhs(i)=0
 2500 continue
c 
cccc        call prin2('rhs as created*',rhs,(npts+n+m)*2)
c 
c       l2-approximate the rhs by a linear combination of
c       Legendre polynomials of orders n+1, . . . n+m
c 
        call cleastsq2(u,v,tria,npts+n+m,n+m,ncols,rhs,sol,work,whts)
c 
c       now, construct the Legendre expansion of the total transfer function
c 
        do 3200 i=1,n
        coetot(i)=coefs(i)-sol(i)
 3200 continue
c 
        do 3400 i=1,m
        coetot(n+i)=-sol(n+i)
 3400 continue
c 
c       for verification, evaluate the total expansion at the
c       nodes on the interval [-1,b]
c 
        totmax=0
        do 3600 i=1,npts
c 
        call clegeexe(t(i),cd,coetot,N+m)
        if(totmax .lt. cdabs(cd)) totmax=cdabs(cd)
 3600 continue
c 
        if(ifprint .ne. 0) call prin2(
     1      'maximum value of total transfer function on [-1,b]*',
     2      totmax,1)
c 
c        evaluate both transfer functions on the interval [-1,1]
c 
        ntest=npts*6
c 
cccc         call prinf('before second gauswhts, ifinit=*',ifinit,1)
        ifwhts=0
        if(ifinit .ne. 0) call gauswhts(ntest,ttest,dummy,ifwhts)
  
  
        do 3800 i=1,ntest
c 
        call clegeexe(ttest(i),f2(i),coetot,N+m)
c 
        if(ifplot .ne. 0)
     1      call clegeexe(ttest(i),f(i),coefs,N-1)
c 
 3800 continue
c 
c       separate the real and imaginary parts of the transfer functions
c       and plot them
c 
        if(ifplot .eq. 0) goto 4300
c 
cccc          call prinf('ifplot=*',ifplot,1)
        do 4000 i=1,ntest
c 
        com=f(i)
        frea(i)=rea(1)
        fima(i)=rea(2)
 4000 continue
c 
        iw=21
        call RSgraf(ttest,frea,ntest,IW,
     1      'real part of original transfer function*')
c 
        iw=22
        call RSgraf(ttest,fima,ntest,IW,
     1      'imaginary part of original transfer function*')
c 
        do 4200 i=1,ntest
c 
        com=f2(i)
        frea(i)=rea(1)
        fima(i)=rea(2)
 4200 continue
c 
        iw=23
        call RSgraf(ttest,frea,ntest,IW,
     1      'real part of sparse transfer function*')
c 
        iw=24
        call RSgraf(ttest,fima,ntest,IW,
     1      'imaginary part of sparse transfer function*')
c 
 4300 continue
c 
c        calculate the number of elements in the array f2 where
c        the transfer function is greater than eps
c 
        eps3=eps*30
c 
        eps3=eps3*5
  
        d=cdabs(rho*rk)
        eps3=eps3/d
c 
        nzero=0
        istart=20
        do 4400 i=istart,ntest
        d=cdabs(f2(i))
        if(d .gt. eps3) goto 4600
        nzero=i
 4400 continue
 4600 continue
c 
        if(ifprint .ne. 0) call prinf('in projcon0, nzero=*',nzero,1)
        if(ifprint .ne. 0) call prinf('out of *',ntest,1)
c 
c        calculate the proportion of the area of the sphere
c        supporting all of the translation operator that is greater
c        than eps
c 
        thresh=ttest(nzero+1)
        if(ifprint .ne. 0) call prin2('thresh=*',thresh,1)
c 
        if(ifprint .ne. 0) call prinf('with ntest=*',ntest,1)
c 
c       finally, calculate the proportion of the area of the
c       sphere supporting the translation operator to the
c       precision eps
c 
        rnum=1-thresh
        prop=rnum/2
c 
        if(ifprint .ne. 0) call prin2(
     1   'proportion of sphere supporting translation operator*',
     2    prop,1)
c 
        ifinit=0
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine onetran3(vect,rk,n,trans,p,hs,
     1    dx,dy,dz,coetot)
        implicit real *8 (a-h,o-z)
        save
        complex *16 rk,trans,hs(1),ima,coetot(1)
        real *8 vect(1),p(1)
         data ima/(0.0d0,1.0d0)/
c 
c       this subroutine constructs a single transfer coefficient
c       (for a single direction) for the conversion of a far-field
c       representation of an h-expansion onto a far-field
c       representation of a j-expansion
c 
c               input parameters:
c 
c  vect - the direction in which the transfer coefficient
c           is being constructed.
c  rk - the helmholtz coefficient
c  n - the total number of Legendre coefficients in the construction
c        of the translation operator
c        construction of the transfer function
c  hs - the array of hankel functions of orders 0 to n
c        of rk*d, where d is the distance between the centers
c        of expansions. note that h0 must be in hs(0),
c        not hs(1), i.e. the array is assumed to have
c        an element number zero.
c  (dx,dy,dz) - the shift vector from the center of the
c        h-expansion  to the center of the j-expansion
c 
c               output parameters:
c 
c  trans - the transfer function in the direction vect
c 
c                work arrays:
c 
c  p - must be at least n+4 real *8 locations long. it is
c        assumed that it has the element number 0.
c 
c        . . . determine the cosine of the angle between the
c              direction of the shift and the direction of the
c              vector along which we are taking the far field
c 
        cosfi=dx*vect(1)+dy*vect(2)+dz*vect(3)
c 
c       find the legendre polynomials
c 
        call legpols(cosfi,n,p)
c 
c        construct the transformation in this direction
c 
        trans=0
        do 2200 i=0,n
c 
        trans=trans+coetot(i+1)*p(i)
 2200 continue
        return
        end
  
c 
c 
c 
c 
c 
        subroutine fldchg(xs,ys,zs,
     1     xtarg,ytarg,ztarg,rk,f)
        implicit real *8 (a-h,o-z)
        save
        complex *16 rk,f,h1
c 
c       this subroutine evaluates the potential of a
c       charge at a point.
c 
c                input parameters:
c 
c  (xs,ys,zs) - the coordinates of the charge
c  (xtarg,ytarg,ztarg) - the coordinates of the point
c       where the potential is to be evaluated
c  rk - the helmholtz coefficients
c 
        d=(xtarg-xs)**2+(ytarg-ys)**2+(ztarg-zs)**2
        d=dsqrt(d)
        call han3d(d*rk,f,h1)
        return
        end
c 
c 
c 
c 
c 
        subroutine hhtran(x0,y0,z0,x,y,z,
     1    vectors,nphi,ntheta,rk,trans)
        implicit real *8 (a-h,o-z)
        save
        complex *16 rk,trans(nphi,ntheta),
     1    ima
        real *8 vectors(3,nphi,ntheta)
         data ima/(0.0d0,1.0d0)/
c 
        entry jjtran(x0,y0,z0,x,y,z,
     1    vectors,nphi,ntheta,rk,trans)
c 
c       this subroutine constructs the transfer coefficients
c       for the conversion of a far-field representation of
c       an h-expansion onto a far-field representation of
c       a h-expansion around a new center. the same code
c       produces transfer coefficients for converting
c       j-expansions into j-expansions (the formulae are
c       identical). for consistency in notation, the entry
c       jjtran is introduced (see immediately above).
c 
c               input parameters:
c 
c  x0,y0,z0 - coordinates of the center of the original h-expansion
c  x,y,z  - coordinates of the center of the shifted h-expansion
c  vectors - the normals to the unity sphere at the
c           standard nodes, as created by subroutine sphecre
c  nphi - the number of elements in the discretization
c        of each parallel
c  ntheta - the number of elements in the discretization
c        of each meridian
c  rk - the helmholtz coefficient
c 
c               output parameters:
c 
c  trans - the transfer function
c 
c        . . . find the shift vecror between the centers
c              of the expansions
c 
        dx=x-x0
        dy=y-y0
        dz=z-z0
ccc     d=dx**2+dy**2+dz**2
ccc     d=dsqrt(d)
c 
c        one node after another, construct the elements
c        of the diagonal operator converting the user-supplied
c        h-expansion into a new h-expansion (both in far-field
c        form)
c 
        do 2000 i=1,ntheta
        do 1800 j=1,nphi
c 
        d=dx*vectors(1,j,i)+dy*vectors(2,j,i)+
     1        dz*vectors(3,j,i)
        trans(j,i)=cdexp(ima*rk*d)
 1800 continue
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine jexeva1(fldcoe,g,nphi,ntheta,fld)
        implicit real *8 (a-h,o-z)
        save
        complex *16 fldcoe(1),g(1),fld
c 
        fld=0
        do 1200 i=1,nphi*ntheta
        fld=fld+g(i)*fldcoe(i)
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine jexeva0(x0,y0,z0,x,y,z,vectors,
     1    rk,nphi,ntheta,whts,fldcoe)
        implicit real *8 (a-h,o-z)
        save
        complex *16 rk,fldcoe(nphi,ntheta),ima,fint,cd,
     1    cdd
        dimension vectors(3,nphi,ntheta),whts(1)
         data ima/(0.0d0,1.0d0)/
c 
c       this subroutine constructs the set of coefficients fldcoe
c       to be used for the evaluation at a point in R^3 of a
c       j-expansion given by its far-field representation. The
c       actual evaluation is performed by the subroutine jexeva1
c       (see) or its exact clone scawrong (see)
c 
c               input parameters:
c 
c  x0,y0,z0 - coordinates of the center of the expansion
c  x,y,z - coordinates of the point where the expansion
c           is to be evaluated
c  vectors - the normals to the unity sphere at the
c           standard nodes, as created by subroutine sphecre
c  rk - the helmholtz coefficient
c  nphi - the number of elements in the discretization
c        of each parallel
c  ntheta - the number of elements in the discretization
c        of each meridian
c  whts - the weights to be used to integrate along the meridians
c 
c               output parameters:
c 
c  fldcoe - the nphi*ntheta coefficients
c 
        done=1
        pi=datan(done)*4
c 
        cd=1
        dx=x-x0
        dy=y-y0
        dz=z-z0
c 
        fint=0
        do 2000 i=1,ntheta
        cdd=0
        do 1800 j=1,nphi
c 
c       calculate the projection of the location of the
c       charge on the direction in question
c 
        proj=vectors(1,j,i)*dx+vectors(2,j,i)*dy+
     1         vectors(3,j,i)*dz
c 
c       calculate the contribution of this component
c       of the far-field
c 
        cd    =cdexp( ima*rk*proj)
c 
        fldcoe(j,i)=cd*whts(i)/nphi/pi/4*ima
 1800 continue
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine jexeval(x0,y0,z0,x,y,z,vectors,
     1    rk,nphi,ntheta,g,whts,fint)
        implicit real *8 (a-h,o-z)
        save
        complex *16 rk,g(nphi,ntheta),ima,fint,cd,
     1    cdd
        dimension vectors(3,nphi,ntheta),whts(1)
         data ima/(0.0d0,1.0d0)/
c 
c       this subroutine evaluates  at a point in R^3 a
c       j-expansion given by its far-field representation
c 
c               input parameters:
c 
c  x0,y0,z0 - coordinates of the center of the expansion
c  x,y,z - coordinates of the point where the expansion
c           is to be evaluated
c  vectors - the normals to the unity sphere at the
c           standard nodes, as created by subroutine sphecre
c  rk - the helmholtz coefficient
c  nphi - the number of elements in the discretization
c        of each parallel
c  ntheta - the number of elements in the discretization
c        of each meridian
c  g - the far-field representation of the potential
c       to be evaluated
c  whts - the weights to be used to integrate along the meridians
c 
c               output parameters:
c 
c  fint - the potential evaluated at the point (x,y,z)
c 
        cd=1
        dx=x-x0
        dy=y-y0
        dz=z-z0
c 
        fint=0
        do 2000 i=1,ntheta
        cdd=0
        do 1800 j=1,nphi
c 
c       calculate the projection of the location of the
c       charge on the direction in question
c 
        proj=vectors(1,j,i)*dx+vectors(2,j,i)*dy+
     1         vectors(3,j,i)*dz
c 
c       calculate the contribution of this component
c       of the far-field
c 
        cd    =cdexp( ima*rk*proj)
        cdd=  cdd+cd * g(j,i)
 1800 continue
        cdd=cdd/nphi*whts(i)
        fint=fint+cdd
 2000 continue
        done=1
        pi=datan(done)*4
        fint=fint/pi/4*ima
        return
        end
c 
c 
c 
c 
c 
        subroutine hjtrans(trans,exp1,nphi,ntheta,
     1     exp2)
        implicit real *8 (a-h,o-z)
        save
        complex *16 trans(nphi,ntheta),
     1     exp1(nphi,ntheta),exp2(nphi,ntheta)
c 
        entry hhtrans(trans,exp1,nphi,ntheta,
     1     exp2)
c 
        entry jjtrans(trans,exp1,nphi,ntheta,
     1     exp2)
c 
c       this subroutine applies a translation operator
c       to an h- or j-expansion given by its far-field.
c       for consistency in notation, the additional
c       entries hhtrans and jjtrans are provided. the
c       actual calculations performed are identical.
c 
c                  input parameters:
c 
c  trans - the transfer function as created by hhtran,
c       jjtran, or hjtran (see)
c  exp1 - the expansion to be shifted
c  nphi - the number of elements in the discretization
c        of each parallel
c  ntheta - the number of elements in the discretization
c        of each meridian
c 
c                output parameters:
c 
c  exp2 - the shifted expansion
c 
        do 2000 i=1,ntheta
        do 1800 j=1,nphi
        exp2(j,i)=exp1(j,i)*trans(j,i)
 1800 continue
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine han3d(z,h0,h1)
        implicit real *8 (a-h,o-z)
        save
        complex *16 z,h0,h1,ima,cd
        data ima/(0.0d0,1.0d0)/,done/1.0d0/
c 
c        this subroutine calculates h0 and h1 - the spherical
c        hankel functions of orders 0, 1. note that this is
c        not a very efficient subroutine. furthermore, it is
c        unstable near z=0.
c 
cccc    h0=-ima/z * cdexp(ima*z)
cccc    h1=-ima/z**2*cdexp(ima*z)-cdexp(ima*z)/z
c 
        cd=cdexp(ima*z)
        h0=-ima*cd/z
        h1=-(ima/z**2+done/z)*cd
c 
        return
        end
c 
c 
c 
c 
c 
c 
        subroutine hfuns3d(z,n,hs)
        implicit real *8 (a-h,o-z)
        save
        complex *16 hs(1),z
c 
c       this subroutine calculates spherical hankel
c       functions of orders 0 through n of the complex
c       argument z.
c 
c                  input parameters:
c 
c  z - the argument
c  n - the highest order of the hankel function
c        to be returned
c 
c                  output parameters:
c 
c  hs - the spherical hankel functions of z.
c        note that this array is assumed to have element
c        number zero.
c 
c       . . . calculate h0,h1
c 
        i0=0
        call han3d(z,hs(i0),hs(1))
c 
c        calculate the rest of the hankel functions
c 
        do 1200 i=1,n-1
        hs(i+1)=(2*i+1)/z*hs(i)-hs(i-1)
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine legpols(x,n,p)
        implicit real *8 (a-h,o-z)
        save
        dimension p(1)
c 
c       this subroutine calculates legendre polynomials
c       functions of orders 0 through n of the real
c       argument x.
c 
c                  input parameters:
c 
c  x - the argument
c  n - the highest order of the legendre polynomial
c        to be returned
c 
c                  output parameters:
c 
c  p - the spherical hankel functions of z.
c        note that this array is assumed to have element
c        number zero.
c 
c 
c        . . . construct p0(x), p1(1)
c 
        i0=0
        p(i0)=1
        p(1) = x
c 
c       conduct the recursion
c 
        do 1200 i=1,n-1
        p(i+1)=( (2*i+1)*x*p(i)-i*p(i-1) )/(i+1)
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine sphecre(nphi,ntheta,phis,thetas,
     1    whts,vectors)
        implicit real *8 (a-h,o-z)
        save
        dimension phis(1),thetas(1),whts(1),
     1    vectors(3,nphi,ntheta)
c 
c       this subroutine discretizes the two-dimensional
c       sphere, generating a standard grid, equispaced
c       in both phi and theta. note that theta varies from
c       0 to pi, i.e. we count from one pole to the other.
c       in addition, coefficients of a quadrature formula
c       are created, to be used by the subroutine spheint
c       (see) to integrate functions on the sphere.
c 
c                input parameters:
c 
c  nphi - the number of elements in the discretization
c        of each parallel
c  ntheta - the number of elements in the discretization
c        of each meridian
c 
c                output parameters:
c 
c  phis - the value of phi at the nodes on each parallel
c        (nphi of them)
c  thetas - the value of theta at the nodes on each meridian
c        (ntheta of them)
c  whts - the weights to be used by the subroutine spheint,
c         (see) to integrate functions on the sphere.
c  vectors - the array of coordinates of normal (exterior)
c         vectors at the nodes on the surface of the sphere
c 
c          . . . create the discretization in phi
c 
        done =1
        pi=datan(done)*4
        h=2*pi/nphi
        do 1200 i=1,nphi
        phis(i)=(i-1)*h
 1200 continue
c 
c       calculate the gaussian weights and nodes to
c       be used in the construction of the quadrature
c       rule on the sphere
c 
cccc        itype=1
cccc        call legeexps(itype,ntheta,thetas,u,v,whts)
c 
        ifwhts=1
        call gauswhts(ntheta,thetas,whts,ifwhts)
  
c 
c        create the discretization in theta, and
c        the appropriate quadrature formula
c 
        do 1400 i=1,ntheta
        thetas(i)=dacos(thetas(i))
        whts(i)=whts(i)*pi*2
 1400 continue
c 
c        create the outer normals at all nodes
c 
        do 2000 i=1,ntheta
        do 1800 j=1,nphi
c 
        vectors(1,j,i)=dsin(thetas(i))*dcos(phis(j))
        vectors(2,j,i)=dsin(thetas(i))*dsin(phis(j))
        vectors(3,j,i)=dcos(thetas(i))
 1800 continue
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine spheint(nphi,ntheta,phis,thetas,
     1     f,fint,whts)
        implicit real *8 (a-h,o-z)
        save
        complex *16 f(nphi,ntheta),fint,cd,cdd
        real *8 phis(1),thetas(1),whts(1)
c 
c       this subroutine integrates a user-specified
c       function on the two-dimensional sphere. it is
c       assumed that the function is tabulated on a standard
c       grid, constructed by the subroutine sphecre (see).
c       in addition, coefficients of a quadrature formula
c       are used, also created by the subroutine sphecre.
c 
c                input parameters:
c 
c  nphi - the number of elements in the discretization
c        of each parallel
c  ntheta - the number of elements in the discretization
c        of each meridian
c  phis - the value of phi at the nodes on each parallel
c        (nphi of them)
c  thetas - the value of theta at the nodes on each meridian
c        (ntheta of them)
c  f - the values of the function to be integrated, it is
c        assumed to be tabulated on the nphi \times ntheta
c        grid, as constructed by the subroutine sphecre.
c  whts - the weights to be used to integrate along the meridians
c 
c                output parameters:
c 
c  fint - the integral of f (complex)
c 
c         . . . integrate the function f on the sphere
c 
        do 2000 i=1,ntheta
c 
c       integrate along i-th parallel
c 
        cd=0
        do 1200 j=1,nphi
        cd=cd+f(j,i)
 1200 continue
cccc    cdd=  cd*whts(i) /nphi * 8
        cdd=  cd*whts(i) /nphi
        fint=fint+cdd
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine hexpan(x0,y0,z0,x,y,z,vectors,
     1    rk,nphi,ntheta,f)
        implicit real *8 (a-h,o-z)
        save
        complex *16 rk,f(nphi,ntheta),ima
        dimension vectors(3,nphi,ntheta)
         data ima/(0.0d0,1.0d0)/
c 
c       this subroutine constructs the far-field
c       representation of the potential created at a
c       unity charge located at the point (x,y,z).
c 
c               input parameters:
c 
c  x0,y0,z0 - coordinates of the center of the expansion
c  x,y,z - coordinates of the charge whose potential
c           is being expanded
c  vectors - the normals to the unity sphere at the
c           standard nodes, as created by subroutine sphecre
c  rk - the helmholtz coefficient
c  nphi - the number of elements in the discretization
c        of each parallel
c  ntheta - the number of elements in the discretization
c        of each meridian
c 
c               output parameters:
c 
c  f - the far-field representation of the potential of
c       the charge being expanded
c 
        dx=x-x0
        dy=y-y0
        dz=z-z0
c 
        do 2000 i=1,ntheta
        do 1800 j=1,nphi
c 
c       calculate the projection of the location of the
c       charge on the direction in question
c 
        proj=vectors(1,j,i)*dx+vectors(2,j,i)*dy+
     1         vectors(3,j,i)*dz
c 
c       calculate the far-field expansion
c 
cccc    f(j,i)=cdexp(-ima*rk*proj) /ima
        f(j,i)=cdexp(-ima*rk*proj) /ima
 1800 continue
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine jfuns3d(z,eps,fjs,n)
        implicit real *8 (a-h,o-z)
        save
        complex *16 z,fjs(1),zinv,com,fj0,scale,fj1
        dimension rea(2)
        equivalence (rea(1),com)
c 
c        this subroutine evaluates all spherical Bessel functions
c        of the complex argument z that are greater (in absolute
c        value) than eps.
c 
c                        input parameters:
c 
c  z - the argument
c  eps - machine zero - all Besel functions that are smaller
c        than eps are declared to be zero.
c 
c                        output parameters:
c 
c  fjs - the array of Bessel functions.
c 
c      Attention: the subroutie assumes that this array has
c                 element number zero!!!!!!
c 
c  n - the number of elements in fjs that are greater (in
c      absolute value) than eps.
c 
c       . . . conduct the recursion up to determine how high
c             we should go
c 
        done=1
        i0=0
        fjs(i0)=0
        fjs(1)=1
        rlarge=done/eps*1.0d17
        zinv=done/z
c 
        do 1200 i=1,10000000
        nn=i
        fjs(i+1)=(2*i+done)*zinv*fjs(i)-fjs(i-1)
c 
        com=fjs(i+1)
        dd=dabs(rea(1))+dabs(rea(2))
        if(dd .gt. rlarge) goto 1400
 1200 continue
c 
 1400 continue
cccc        call prinf('nn=*',nn,1)
c 
c       conduct recursion down
c 
        fjs(nn+1)=0
        fjs(nn)=1
        do 2200 i=nn,1,-1
c 
        fjs(i-1)=(2*i+done)*zinv*fjs(i)-fjs(i+1)
 2200 continue
cccc         call prin2('fjs before scaling*',fjs(i0),nn*2+2)
c 
c       determine the scaling parameter
c 
        fj0=cdsin(z)*zinv
        fj1=fj0*zinv-cdcos(z)*zinv
c 
        com=fj0
        d0=dabs(rea(1))+dabs(rea(2))
        com=fj1
        d1=dabs(rea(1))+dabs(rea(2))
c 
        if(d1 .gt. d0) goto 2250
        scale=fj0/fjs(i0)
        goto 2300
 2250 continue
c 
        scale=fj1/fjs(1)
 2300 continue
cccc        call prin2('scale=*',scale,2)
c 
c       . . . scale
c 
        do 2400 i=i0,nn
        fjs(i)=fjs(i)*scale
 2400 continue
c 
c       find the number of Bessel functions that are greater than eps
c 
        do 2600 i=nn,1,-1
        n=i
        com=fjs(i)
        dd=dabs(rea(1))+dabs(rea(2))
        if(dd .gt. eps) goto 2800
 2600 continue
 2800 continue
        return
        end
