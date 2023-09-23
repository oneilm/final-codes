        implicit real *8 (a-h,o-z)
        dimension ts(1000),rhos(1000),angles(1000),
     1      apert(1000),errs(1000),nnnodes(1000)
c 
        CALL PRINI(6,13)
C 
c       . . . construct the testing parameters
c 
         PRINT *, 'ENTER r'
         READ *,r
         CALL PRIn2('r=*',r,1)
c 
         PRINT *, 'ENTER rho'
         READ *,rho
         CALL PRIn2('rho=*',rho,1)
c 
         ifinit=1
         ifplot=0
         call alltrtst(r,rho,angtotal,ifinit,errmax,nnodes,
     1       ifplot)
         call prin2('after alltrtst, angtotal=*',angtotal,1)
c 
c        construct a graph of angtotal as a function of rho
c 
         nrhos=60
cccc         nrhos=6
         rho1=r*4
         rho2=r*100
         hrho=(rho2-rho1)/nrhos
         done=1
         pi=datan(done)*4
c 
         ifinit=0
         ifplot=0
c 
         do 1400 i=1,nrhos
  
         if(i .eq. 10) ifplot=1
c 
         call prinf('in the main program, i=*',i,1)
c 
         rhos(i)=(i-1)*hrho+rho1
         ts(i)=rhos(i)/r
         call alltrtst(r,rhos(i),angles(i),ifinit,errs(i),
     1       nnnodes(i),ifplot)
         angles(i)=angles(i)/(2*pi)
c 
         apert(i)=datan(r/rhos(i))/(2*pi)
c 
         ifplot=0
 1400 continue
         call prin2('rhos=*',rhos,nrhos)
         call prin2('ts=*',ts,nrhos)
         call prin2('angles=*',angles,nrhos)
         call prin2('errs=*',errs,nrhos)
         call prinf('nnnodes=*',nnnodes,nrhos)
c 
c        now, plot the part of the circle that is non-zero as
c        a function of rho
c 
        iw=41
        call RSGRAF(ts,angles,nrhos,IW,
     1      'angles as a function of lambda r*')
  
c 
c        for comparison, plot the aperture  as
c        a function of rho
c 
        iw=42
        call RSGRAF(ts,apert,nrhos,IW,
     1      'aperture as a function of lambda r*')
  
  
         stop
         end
  
  
  
  
  
  
        subroutine alltrtst(r,rho,angtotal,ifinit,errmax,
     1      nnodes,ifplot)
        implicit real *8 (a-h,o-z)
        save
        complex *16 ima,rk
         real *8
     1   xsource(10 000),ysource(10 000),xreceiv(10 000),
     2       yreceiv(10 000),relerrs(10 000)
         complex *16 hconv(100 000),
     1       jconv(100 000),
     2       errs(10 000),vals(10 000),work(100 000)
         data ima/(0.0d0,1.0d0)/
c 
        lwork=200 000
        ntest=3
        rk=1
        eps=1.0d-3
        coef=6
ccccc        ifplot=0
        ifprint=0
c 
c        construct the geometry of the test
c 
c        . . . the center of the h-expansion
c 
        xh1=0
        yh1=0
c 
c        . . . the  center of the j-expansion
c 
        xj1=xh1+rho
        yj1=yh1
c 
cccc        yj1=yj1-1000
c 
         if(2 .ne. 3) goto 2500
 1200 format(2x,'xh1= ',e11.5,' yh1= ',e11.5)
 1400 format(10x)
         write(6,1400)
         write(13,1400)
         write(6,1200) xh1,yh1
         write(13,1200) xh1,yh1
c 
 2400 format(2x,'xj1= ',e11.5,' yj1= ',e11.5)
         write(6,1400)
         write(13,1400)
         write(6,2400) xj1,yj1
         write(13,2400) xj1,yj1
 2500 continue
c 
c        check the accuracy for all sorts of pairs of sources and
c        receivers living on two circles
c 
c        . . . create the test circles
c 
         done=1
         pi=datan(done)*4
         h=2*pi/ntest
         do 3200 i=1,ntest
         tt=(i-1)*h
         xsource(i)=r*dcos(tt)+xh1
         ysource(i)=r*dsin(tt)+yh1
c 
         xreceiv(i)=r*dcos(tt)+xj1
         yreceiv(i)=r*dsin(tt)+yj1
 3200 continue
  
          if(ifinit .eq. 0) goto 3300
  
         call prin2('xsource as created*',xsource,ntest)
         call prin2('ysource as created*',ysource,ntest)
         call prin2('xreceiv as created*',xreceiv,ntest)
         call prin2('yreceiv as created*',yreceiv,ntest)
c 
 3300 continue
c 
        rho=(xj1-xh1)**2+(yj1-yh1)**2
        rho=dsqrt(rho)
          if(ifinit .eq. 1)
cccc     1      call getpara(m,r,rho,rk,eps/2,njs,a,b,n,nhs,coef,
     1      call getpara(m,r,rk,eps/2,njs,a,n,nhs,coef,
     2       work)
         nphi=m
         call hjtrtest(rk,nphi,r,eps,
     1    xh1,yh1,xj1,yj1,
     2      errs,vals,
     3      ifinit,hconv,jconv,ntest,xsource,xreceiv,
     4      ysource,yreceiv,ifplot,ifprint,work,lwork,angtotal,
     5       ifinit,nnodes)
c 
cccc          call prin2('vals are*',vals,ntest**2 *2)
cccc          call prin2('end errors are*',errs,ntest**2*2)
c 
          errmax=0
          do 3400 i=1,ntest**2
          relerrs(i)=cdabs(errs(i)/vals(i))
          if(errmax .lt. relerrs(i) ) errmax=relerrs(i)
 3400 continue
c 
ccccc          call prin2('and relative errors are*',relerrs,ntest**2)
          call prin2('and errmax=*',errmax,1)
         return
         end
c 
c 
c 
c 
c 
        subroutine hjtrtest(rk,nphi,r,eps,
     1    xh1,yh1,xj1,yj1,errs,f2s,
     2      ifinit,hconv,jconv,ntest,xsource,xreceiv,
     3      ysource,yreceiv,ifplot,ifprint,work,lwork,angtotal,
     4      ifinit,nnodes)
        implicit real *8 (a-h,o-z)
        save
        complex *16 ima,rk,hconv(nphi,1),
     1      jconv(nphi,1)
         real *8 phis(10 000),
     1    vectors(30 000),xsource(1),xreceiv(1),
     2       ysource(1),yreceiv(1),work(1)
         complex *16 jexp1(10 000),
     1    f1,f2,trhj21(10 000),errs(1),f2s(1)
         data ima/(0.0d0,1.0d0)/
c 
c            initialize the hankel function computer
c 
         call HBFUNI(RK)
c 
c 
c       construct the sparse h->j transformation vector
c 
        coef=6
        rho=(xj1-xh1)**2+(yj1-yh1)**2
        rho=dsqrt(rho)
        xshift=xj1-xh1
        yshift=yj1-yh1
c 
        call hjtrsprs(ier,r,xshift,yshift,rk,nj,nphi,eps,
     1    trhj21,ifprint,nnodes,angtotal,istart,
     2    ifinit,work,lwork,ltot)
c 
cccc            call prinf('after hjtrsprs, istart=*',istart,1)
         if(ifinit .eq. 0) goto 1100
c 
         call prinf('after hjtrsprs, lwork=*',lwork,1)
         call prinf('after hjtrsprs, ltot=*',ltot,1)
         call prin2('for comparison, 120 |r||rk|+300    *',
     1       120*r*cdabs(rk)+300,1)
c 
        call prinf('after hjtrsprs, nphi=*',nphi,1)
        call prinf('after hjtrsprs, nnodes=*',nnodes,1)
        call prin2('after hjtrsprs, angtotal=*',angtotal,1)
 1100 continue
c 
c        set to zero the part of the translation
c        operator that is smaller than eps, in the opinion
c        of trhj21
c 
        do 1200 i=1,nphi
        jexp1(i)=0
 1200 continue
c 
        do 1250 i=1,nnodes
        j=iwrap(istart+i-1,nphi)
        jexp1(j)=trhj21(j)
 1250 continue
c 
        do 1260 i=1,nphi
        trhj21(i)=jexp1(i)
 1260 continue
c 
c        construct the circle and the whole structure on it
c 
        call circre(nphi,phis,vectors)
c 
c        if the user so requested, plot the real and imaginary
c        parts of the translation operator. also, its absolute value
c 
        if(ifplot .eq. 0) goto 1300
  
        iwrea=21
        iwima=22
        iwabs=23
c 
c        if the user so requested, plot the real and imaginary
c        parts of the translation operator. also, its absolute value
c 
        call cmplxgrf(phis,trhj21,nphi,iwrea,iwima,
     1  'real part of the translation operator*',
     2  'imaginary part of the translation operator*',
     3      iwabs,
     4  'absolute value  of the translation operator*')
c 
 1300 continue
c 
c        for each charge, construct the expansion of the
c        potential of the charge around the center of the
c        h-expansion. also construct the vectors to be
c        used to convert the local expansions around the
c        center of the j-expansion into the values at all
c        receivers
c 
        do 1400 i=1,ntest
c 
        call hexp2d(xh1,yh1,xsource(i),ysource(i),vectors,
     1    rk,nphi,hconv(1,i) )
c 
        call jextoval(xj1,yj1,xreceiv(i),yreceiv(i),vectors,
     1    rk,nphi,jconv(1,i) )
c 
 1400 continue
c 
c        now, for each pair (source-receiver), evaluate the
c        interaction vea expansions, and directly for
c        comparison
c 
        istore=0
        do 3000 i=1,ntest
  
  
        call hjtrans(trhj21,hconv(1,i),nphi,jexp1)
c 
        do 2800 j=1,ntest
c 
        istore=istore+1
c 
c       convert the i-th h-expansion into a j-expansion
c       around the j-center
c 
c 
c        evaluate the j-expansion at the target point
c 
        f1=0
        do 1600 k=1,nphi
        f1=f1+jexp1(k)*jconv(k,j)
 1600 continue
cccc         call prin2('f1 via expansions is*',f1,2)
c 
c        evaluate the potential at the target point
c        directly
c 
        call fldch2d(xsource(i),ysource(i),
     1      xreceiv(j),yreceiv(j),rk,f2)
cccc        call prin2('and f2 evaluated directly is*',f2,2)
cccc        call prin2('and the difference  is*',f1-f2,2)
c 
        errs(istore)=f1-f2
        f2s(istore)=f2
 2800 continue
 3000 continue
c 
         return
         end
c 
c 
c 
c 
c 
        subroutine cmplxgrf(t,z,n,iwrea,iwima,mesrea,mesima,
     1      iwabs,mesabs)
        implicit real *8 (a-h,o-z)
        save
        complex *16 z(1),com
        dimension t(1),rea(2),w(10 000)
        character *1 mesrea(1),mesima(1),mesabs(1)
        equivalence (rea(1),com)
c 
c       plot the real part
c 
        do 1200 i=1,n
        com=z(i)
        w(i)=rea(1)
 1200 continue
        call RSGRAF(t,w,n,IWrea,mesrea)
c 
c       plot the imaginary part
c 
        do 1400 i=1,n
        com=z(i)
        w(i)=rea(2)
 1400 continue
        call RSGRAF(t,w,n,IWima,mesima)
c 
c       plot the absolute value
c 
        do 1600 i=1,n
        w(i)=cdabs(z(i))
 1600 continue
        call RSGRAF(t,w,n,IWabs,mesabs)
        return
        end
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c       This is the end of the debugging code and the beginning
c       of the translation operator code proper.
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        subroutine hjtrsprs(ier,r,xshift,yshift,rk,nj,m,eps,
     1    trans,ifprint,nnodes,angtotal,istart,
     2    ifinit,w,lwork,ltot)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1)
        complex *16 rk,trans(1)
c 
c       this subroutine constructs the sparse diagonal form of the
c       h->j translation operator, with parameters specified by
c       the user. The subroutine makes several decisions that
c       are probably less than optimal. However, it is believed
c       that the bell produced by it is within a factor of 2
c       of the optimal one. if you do not like it, please feel
c       free to write your own.
c 
c                           input parameters:
c 
c  r - the radius of the collection of charges whose potential
c       is to be evaluated
c  (xshift,yshift) - the vector by which h->j translation is
c       to be performed.
c  rk - the helmholtz coefficient
c  eps - the precision to which the calculations are to be performed
c  ifprint - tells the subroutine whether the user wants it to describe
c       in print ojn FORTRAN files 6 and 13 the state of certain of
c       its internal variables.
c     ifprint=0 suppresses the output
c     ifprint=1 activates the output
c       This is a debugging/testing feature; normally, ifprint=0.
c  ifinit - tells the subroutine if it should perform certain
c       initializations.
c     ifinit=1 should be used either when the subroutine is being
c       called for the first time, or when one of the following
c       parameters has changed since the preceding call:
c 
c       r, rk, eps.
c 
c     ifinit=0 should be set when repeated calls are made to
c       the subroutine with the parameters r, rk, eps unchanged
c       (in other words, only (xshift,yshift) have changed). Note
c       also that the first  ltot  elements of the work array
c       w must not be altered between to calls to this subroutine.
c  lwork - the length of the work array w. A reasonable estimate for
c       the minimum length is 120 |r||rk|+300. if lwork is too small,
c       ier is set to 16, and the execution is terminated.
c 
c 
c                           output parameters:
c 
c  ier - error return code. ier=0 means successful conclusion.
c                           ier=16 means that too little memory
c                           has been provided in array w (see below).
c  nj - the number of bessel functions of  rk*r that are bigger
c       than eps. provided solely for the amusement of the user.
c  m - the number of nodes in the discretization of the horizon
c  trans - the translation operator
c  nnodes - the number of values in trans that are greater than eps
c       (in relative terms)
c  angtotal - the angle covered by non-zero part of the translation
c       operator (in radians)
c  istart - the first point in array trans whose absolute value is
c       grater than eps (relatively speaking)
c  ltot - the length of the part of the array w that has been used
c       by the subroutine. this part should not be changed between
c       calls to this subroutine. should it be changed, the following
c       call to this subroutine should be with ifinit=1 (see above).
c 
c                          work array:
c 
c  w - must be at least 120 |r||rk|+300 elements long. this is a
c       fairly crude estimate; it tends to be excessively large
c       for large-scale problems. check ltot (see above) on exit.
c 
c 
c       . . . get the parameters for the construction of the spar!se
c             h-to-j translation operator
c 
c 
        iw1=6
        iw2=13
        coef=6
        rho=dsqrt(xshift**2+yshift**2)
        b=rho
        if(ifinit .ne. 0)
cccc     1      call getpara(m,r,rho,rk,eps/2,njs,a,bbb,n,nhs,coef,
     1      call getpara(m,r,rk,eps/2,njs,a,n,nhs,coef,
     2                   w)
c 
        if(ifprint .eq. 0) goto 3000
c 
 1200 format('    ',/,'    ',/,'    ')
 1400 format(' after getpara, njs (number of Bessel functions ',/,
     1       'of rk*r that are greater than eps) is ',i6)
      write(iw1,1200)
      write(iw2,1200)
c 
      write(iw1,1400) njs
      write(iw2,1400) njs
c 
 1600 format(' after getpara, m (total number of rays discretizing ',/,
     1       'the horizon) is',i6)
      write(iw1,1200)
      write(iw2,1200)
c 
      write(iw1,1600) m
      write(iw2,1600) m
c 
c 
 1800 format(' after getpara, a (the parameter in the exponent ',/,
     1       ' of the Gaussian bell) is ',e11.5)
      write(iw1,1200)
      write(iw2,1200)
c 
      write(iw1,1800) a
      write(iw2,1800) a
c 
 2000 format(' after getpara, b (distance between the expansions  ',/,
     1       'to be shifted) is ',e11.5)
      write(iw1,1200)
      write(iw2,1200)
c 
      write(iw1,2000) b
      write(iw2,2000) b
c 
 2200 format(' after getpara, n (number of positive, and also negative,',
     1   /,'exponentials in the truncated approximation to the ',/,
     2   'delta function) is',i6)
      write(iw1,1200)
      write(iw2,1200)
c 
      write(iw1,2200) n
      write(iw2,2200) n
c 
 2400 format(' after getpara, nhs (number of hankel functions of ',
     1   /,'b*k to be constructed) is ',i6)
      write(iw1,1200)
      write(iw2,1200)
c 
      write(iw1,2400) nhs
      write(iw2,2400) nhs
c 
      write(iw1,1200)
      write(iw2,1200)
      write(iw1,1200)
      write(iw2,1200)
c 
 3000 continue
c 
c       allocate memory for the subroutine hjtrcr22
c       finding the sparse form of the diagonal operator
c 
c 
        if(ifinit .eq. 0) goto 3400
c 
        ithetas=1
        lthetas=m+1
c 
        itrans0=ithetas+lthetas
        ltrans0=m*2+1
c 
        iwsave=itrans0+ltrans0
        lwsave=5*m+20
c 
        ihs=iwsave+lwsave
        lhs=m*2+10
c 
        ltot=ihs+lhs
c 
 3400 continue
c 
c       if the user-supplied memory is insufficient - bomb
c 
        if(ltot .le. lwork) goto 3600
        ier=16
        return
 3600 continue
c 
c       construct the translation operator
c 
        call hjtrcr22(n,m,a,rk,nhs,eps,trans,njs,
     1      nnodes,angtotal,ifinit,
     2      w(ithetas),w(itrans0),w(iwsave),w(ihs),
     3      xshift,yshift,istart)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine getpara(m,r,rk,eps,nj,a,n,nhs,coef,
     1      work)
        implicit real *8 (a-h,o-z)
        save
        complex *16 rk,work(1)
c 
c        this subroutine determines the parameters to be used for
c        the design of the ray-type translation operator for the
c        FMM for the Helmholtz equation.
c 
c 
c                       input parameters:
c 
c  r - the radius of the collection of charges whose potential
c       is to be evaluated
c  (xshift,yshift) - the vector by which h->j translation is
c       to be performed.
c  rk - the helmholtz coefficient
c  eps - the precision to which the calculations are to be performed
c  coef - a tuning parameter. should be close to 6.
c 
c                           output parameters:
c 
c  m - the number of nodes in the discretization of the horizon
c       (in relative terms)
c  nj - the number of bessel functions of  rk*r that are bigger
c       than eps. provided solely for the amusement of the user.
c  a - the parameter inside the exponent in (1) above - specifies
c          the width of the Gaussian bell used.
c  n - the number of exponentials in the expansion of the
c          delta-function to be used by the algorithm (see (1) above )
c  nhs - the highest order of Hankel function of the argument rk*b
c          that will be created.
c 
c                           work array:
c 
c  w - must be fairly long (something like |rk*r| *2 + 300
c 
c        . . . determine the number nj of Bessel functions J_m(rk*r)
c              whose absolute values are greater than eps
c 
cccc         call prin2('in getpara, rk=*',rk,2)
cccc         call prin2('in getpara, r=*',r,1)
         call jBFUN(rk*r,work(2),Nj,EPS)
c 
cccc          call prinf('in getpara after jbfun, nj=*',nj,1)
c 
c       set the total number of nodes in the discretizatin
c       of the horizon to nj*6
c 
         m=nj*coef
c 
         call hjmfind(m,nmin,ktwos,kthrees,kfives)
         m=nmin
c 
c        now, the total width of the Fourier transform of the
c        bell is (coef*nj-4*nj)/2=(coef-4)*nj/2. one wing of
c        the bell has the width of (coef-4)*nj/2 /2.
c 
c         . . . determine the parameter of the bell
c 
         lwing=(coef-4)*nj/2 /2
cccc          call prinf('old lwing=*',lwing,1)
         lwing=(m-4*nj)/2 /2
cccc          call prinf('new lwing=*',lwing,1)
c 
         x01=0.001
         x02=lwing**2*2
         call fndxofm(lwing,eps/2,a,work(2),x01,x02)
c 
cccc          call prin2('in getpara after fndxofm a=*',a,1)
c 
c       determine the number of hankel functions to be computed
c 
        nhs=m/2
c 
c       finally, determine the number of ones in the truncated
c       fourier series of the delte-function
c 
        n=2*nj+lwing
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine hjtrcr22(n,m,a,rk,nhs,eps,trans,nj,
     1      nnodes,angtotal,ifinit,
     2      thetas,trans0,wsave,hs,
     3      xshift,yshift,istart)
c 
        implicit real *8 (a-h,o-z)
        save
        complex *16 ima,trans(1),com,trans0(1),
     1      rk,hs(1),cd,cdd,cd1,cdd1
        dimension thetas(1),rea(2),wsave(1)
c 
        equivalence(rea(1),com)
c 
        data ima/(0.0d0,1.0d0)/
c 
c        this subroutine constructs the beam-form translation
c        operator for the helmholtz equation. The operator
c        corresponds to the shift by the vector (xshift,yshift).
c 
c        the transfer function is constructed by specifying its
c        far-field representation, and then bringing it in as
c        a solution of the helmholtz equation. the far-field
c        representation is given by the formula
c 
c        f(theta)= dexp(a*(cos(theta)-1) *
c                                                                 (1)
c        \sum_{j=-n}^n cdexp(ima*j*theta).
c 
c        in other words, the far-field representation of the
c        transfer function is a product of the gaussian
c        bell dexp(a*(cos(theta)-1) with a finite-frequency
c        approximation to the delta-function ("truncated
c        delta-function"). the actual transfer function is
c        calculated via the formula
c 
c        trans(theta)= \sum_{j=-n}^n ( f_j *
c                                                                 (2)
c                      cdexp(ima*j*theta) * *H_j(rk * b) )
c 
c        with f_j the fourier coefficients of the far-field of the
c        function f (see (1) above).
c 
c 
c                          input parameters:
c 
c  n - the number of exponentials in the expansion of the
c          delta-function to be used by the algorithm (see (1) above )
c  m - the number of nodes in the discretization of the horizon
c  a - the parameter inside the exponent in (1) above - specifies
c          the width of the Gaussian bell used.
c  rk - the helmholtz coefficient
c  nhs - the highest order of Hankel function of the argument rk*b
c          that will be created.
c  eps - the relative precision to which the computations will be
c          performed. note that by increasing eps, the user improves
c          sparsity
c 
c  nj - the number of bessel functions of  rk*r that are bigger
c       than eps. provided solely for the amusement of the user.
c  ifinit - tells the subroutine if it should perform certain
c       initializations.
c     ifinit=1 means that the subroutine recomputes all parameters it
c          needs.
c     ifinit=0 means that the subroutine recomputes only parameters
c          having to do with the vector (xshift,yshift). The
c          parameters having to do with the patches themselves
c          (rk, r,m, etc.) are not recomputed.
c  thetas - the m-point discretization of the circle. this is an input
c          parameter only when ifinit=0; otherwise, it is an output
c          parameter.
c  trans0 - fourier series of the function f (see (1) above).
c          this is an input parameter only when ifinit=0; otherwise,
c          it is an output parameter.
c  wsave - the auxiliary arrayused by the fourier transform routines.
c          this is an input parameter only when ifinit=0; otherwise,
c          it is an output parameter.
c  (xshift,yshift) - the vector by which h->j translation is
c       to be performed.
c 
c 
c                           output parameters:
c 
c  trans - the translation operator
c  nnodes - the number of values in trans that are greater than eps
c       (in relative terms)
c  angtotal - the angle covered by non-zero part of the translation
c       operator (in radians)
c  istart - the first point in array trans whose absolute value is
c       grater than eps (relatively speaking)
c 
c                          work array:
c 
c  hs  - must be at least nhs+1 complex *16 elements long
c 
c         . . . if the user so ordered - bypass the initialization
c 
        if(ifinit .eq. 0) goto 2635
c 
c        . . . onstruct the array of angles
c 
        done=1
        pi=datan(done)*4
        h=2*pi/m
        do 1200 i=1,m
        thetas(i)=(i-1)*h
 1200 continue
c 
c        construct the truncated delta-function
c 
        trans0(1)=2*n+1
        do 1800 i=2,m
        trans0(i)=(done-cdexp(ima*(2*n+1)*thetas(i)) )
     1      /(done-cdexp(ima*thetas(i)))
        trans0(i)=trans0(i)*cdexp(-ima*n*thetas(i))
 1800 continue
c 
c       multiply the truncated delta-function by a tapering bell
c 
        do 2200 i=1,m
        bell=dexp(  a*(dcos(thetas(i))-1)    )
        trans0(i)=trans0(i)*bell
 2200 continue
c 
c        now, calculate the fourier series of the tapered delta-function
c 
        call DCFFTI(m,WSAVE)
        CALL DCFFTF(m,trans0,WSAVE)
c 
c         construct the hankel functions of the user-specified argument
c 
        zero=0
        do 2565 i=1,m/2+1
        hs(i)=zero
 2565 continue
c 
        dsig=1
        do 2600 i=1,m/2
        dsig=-dsig
        trans0(i)=trans0(i)*cdexp(ima*(i-1)*pi/2) /m
c 
        trans0(m-i+1)=trans0(m-i+1)*cdexp(-ima*i*pi/2) /m
     1      *dsig
 2600 continue
c 
        call HBFUNlI(RK)
c 
 2635 continue
c 
        b=dsqrt(xshift**2+yshift**2)
        call HBFUNl(b,Hs(2),Nhs)
c 
c       multiply the fourier series of the tapered delta-function
c       by the appropriate hankel functions. at the same time,
c       them by the phase factor to rotate the beam.
c 
        dd=xshift**2+yshift**2
        dd=dsqrt(dd)
        cd=xshift/dd-ima*yshift/dd
        done=1
        cd1=done/cd
        cdd=done
        cdd1=cd1
c 
        do 2650 i=1,m/2
c 
        trans(i)=trans0(i)* hs(i) * cdd
        trans(m-i+1)=trans0(m-i+1)*hs(i+1)*cdd1
c 
        cdd=cdd*cd
        cdd1=cdd1*cd1
 2650 continue
c 
c        fourier-transform the junk obtaining the final transfer
c        function
c 
        CALL DCFFTb(m,trans,WSAVE)
c 
c        calculate the angular width of the beam and the number
c        of nodes in it. also the location of the beam on
c        the circle.
c 
        call hjnnfind(xshift,yshift,trans,m,eps,
     1      istart,nnodes)
c 
        angtotal=h*nnodes
c 
         return
         end
c 
c 
c 
c 
c 
        subroutine hjnnfind(xshift,yshift,f,n,eps,
     1      istart,nnodes)
        implicit real *8 (a-h,o-z)
        save
        complex *16 com,f(1)
        real *8 rea(2)
c 
        data twopi/0.6283185307179586d+01/,
     1      pi/0.3141592653589793d+01/,
     2      halfpi/0.1570796326794897d+01/
c 
        equivalence (rea(1),com)
c 
c        this subroutine finds the chunk in the complex array f
c        where the elements of f are greater than eps.
c 
c 
c        . . . find the point on the circle at which the translation
c              operator is minimized
c 
        if(dabs(yshift) .gt. dabs(xshift)) goto 1200
c 
        phi=datan(yshift/xshift)
        if(xshift .lt. 0) phi=phi-pi
        if(phi .lt. 0) phi=phi+twopi
        goto 1400
 1200 continue
c 
        phi=datan(xshift/yshift)
        phi=halfpi-phi
        if(yshift .lt. 0) phi=phi+pi
 1400 continue
c 
cccc        call prin2('phi as reconstructed =*',phi,1)
c 
        phimax=phi
        phimin=phi+pi
        if(phimin .gt. twopi) phimin=phimin-twopi
cccc        call prin2('phimax=*',phimax,1)
cccc        call prin2('phimin=*',phimin,1)
c 
c 
c       find the number of values of the translation operator that
c       are greater than eps (scaled to the maximum)
c 
        h=twopi/n
        imax=phimax/h
        imin=phimin/h
ccccc         call prinf('imin=*',imin,1)
cccc         call prinf('imax=*',imax,1)
         thresh=cdabs(f(iwrap(imax,n)))*eps
cccc         call prin2('thresh=*',thresh,1)
c 
c        . . . count zeroes forward
c 
        nforward=0
        do 1600 i=imin,n  +imin
        com=f(iwrap(i,n))
        d=dabs(rea(1))+dabs(rea(2))
        if(d. gt. thresh) goto 1800
        nforward=nforward+1
 1600 continue
 1800 continue
ccccc         call prinf('nforward as calculated*',nforward,1)
c 
c       count zeroes backward
c 
 2000 continue
c 
        nback=0
        do 2600 i=imin-1,-(n-imin),-1
        com=f(iwrap(i,n))
        d=dabs(rea(1))+dabs(rea(2))
        if(d. gt. thresh) goto 2800
        nback=nback+1
 2600 continue
 2800 continue
ccccc         call prinf('nback as calculated*',nback,1)
c 
c       calculate the number of non-zero nodes and the start of
c       the non-zero region
c 
        nnodes=n-(nforward+nback)
        istart=imin+(nforward+nback)/2-1
        if(cdabs(f(iwrap(istart,n))) .lt. thresh) istart=istart+1
        if(cdabs(f(iwrap(istart,n))) .lt. thresh) istart=istart+1
        if(cdabs(f(iwrap(istart,n))) .lt. thresh) istart=istart+1
c 
        if(istart .gt. n) istart=istart-n
        return
        end
c 
c 
c 
c 
c 
        function iwrap(i,n)
        save
        if( (i .gt. n) .or. (i .lt. 1) ) goto 1200
        iwrap=i
        return
 1200 continue
c 
        if(i. gt. n) goto 1400
        iwrap=n+i
        return
 1400 continue
c 
        iwrap=i-n
        return
        end
  
c 
c 
c 
c 
c 
        subroutine jextoval(x0,y0,x,y,vectors,
     1    rk,nphi,coefs)
        implicit real *8 (a-h,o-z)
        save
        complex *16 rk,coefs(nphi),ima,cd,
     1    cdd
        dimension vectors(2,nphi)
         data ima/(0.0d0,1.0d0)/
c 
c       this subroutine constructs the coefficients of the
c       linear form connecting the far-field representation of a
c       j-expansion with its values at the point (x,y)
c 
c               input parameters:
c 
c  x0,y0 - coordinates of the center of the expansion
c  x,y - coordinates of the point where the expansion
c           is to be evaluated
c  vectors - the normals to the unity circle at the
c          equispaced nodes
c  rk - the helmholtz coefficient
c  nphi - the number of elements in the discretization
c        of each parallel
c 
c               output parameters:
c 
c  coefs - the coefficients of the linear form connecting the
c        far-field representation of a j-expansion with its values
c        at the point (x,y)
c 
        dx=x-x0
        dy=y-y0
c 
        done=1
        pi=datan(done)*4
        h=2*pi/nphi
c 
        cdd=0
        do 1800 j=1,nphi
c 
c       calculate the projection of the location of the
c       charge on the direction in question
c 
        proj=vectors(1,j)*dx+vectors(2,j)*dy
c 
c       calculate the contribution of this component
c       of the far-field
c 
        cd    =cdexp( ima*rk*proj)
        coefs(j)=cd *h*ima/pi/2
cccc        coefs(j)=-coefs(j)
 1800 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine circre(nphi,phis,vectors)
        implicit real *8 (a-h,o-z)
        save
        dimension phis(1),vectors(2,nphi)
c 
c       this subroutine discretizes the circle in the
c       equispaced manner.
c 
c                input parameters:
c 
c  nphi - the number of elements in the discretization
c        of each parallel
c 
c                output parameters:
c 
c  phis - the value of phi at the nodes on each parallel
c        (nphi of them)
c  vectors - the array of coordinates of normal (exterior)
c         vectors at the nodes on the circle
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
c        create the outer normals at all nodes
c 
        do 1800 j=1,nphi
        vectors(1,j)=dcos(phis(j))
        vectors(2,j)=dsin(phis(j))
 1800 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine hexp2d(x0,y0,x,y,vectors,
     1    rk,nphi,f)
        implicit real *8 (a-h,o-z)
        save
        complex *16 rk,f(1),ima
        dimension vectors(2,nphi)
         data ima/(0.0d0,1.0d0)/
c 
c       this subroutine constructs the far-field
c       representation of the potential created at a
c       unity charge located at the point (x,y,z).
c 
c               input parameters:
c 
c  x0,y0 - coordinates of the center of the expansion
c  x,y - coordinates of the charge whose potential
c           is being expanded
c  vectors - the normals to the unity circle
c  rk - the helmholtz coefficient
c  nphi - the number of elements in the discretization
c        of the circle
c 
c               output parameters:
c 
c  f - the far-field representation of the potential of
c       the charge being expanded
c 
        dx=x-x0
        dy=y-y0
c 
        do 1800 j=1,nphi
c 
c       calculate the projection of the location of the
c       charge on the direction in question
c 
        proj=vectors(1,j)*dx+vectors(2,j)*dy
c 
c       calculate the far-field expansion
c 
        f(j)=cdexp(-ima*rk*proj) /ima
 1800 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine fldch2d(xs,ys,
     1     xtarg,ytarg,rk,f)
        implicit real *8 (a-h,o-z)
        save
        complex *16 rk,f,h17
c 
c       this subroutine evaluates the potential of a
c       charge at a point.
c 
c                input parameters:
c 
c  (xs,ys) - the coordinates of the charge
c  (xtarg,ytarg) - the coordinates of the point
c       where the potential is to be evaluated
c  rk - the helmholtz coefficients
c 
        d=(xtarg-xs)**2+(ytarg-ys)**2
        d=dsqrt(d)
cccc    call han3d(d*rk,f,h1)
cccc    call prin2('in fldch2d, d=*',d,1)
        call hfun(d,f,h17)
cccc       call prin2('after hfun,f=*',f,2)
        return
        end
c 
c 
c 
c 
c 
        subroutine hhtran2d(x0,y0,x,y,
     1    vectors,nphi,rk,trans)
        implicit real *8 (a-h,o-z)
        save
        complex *16 rk,trans(nphi),ima
        real *8 vectors(2,nphi)
         data ima/(0.0d0,1.0d0)/
c 
        entry jjtran(x0,y0,x,y,
     1    vectors,nphi,rk,trans)
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
c  x0,y0 - coordinates of the center of the original h-expansion
c  x,y,z  - coordinates of the center of the shifted h-expansion
c  vectors - the normals to the unity sphere at the
c           standard nodes, as created by subroutine sphecre
c  nphi - the number of elements in the discretization
c        of each parallel
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
c 
c        one node after another, construct the elements
c        of the diagonal operator converting the user-supplied
c        h-expansion into a new h-expansion (both in far-field
c        form)
c 
cccc     call prin2('in hhtran2d, d=*',d,1)
        do 1800 j=1,nphi
cccc     call prinf('in hhtran2d, j=*',j,1)
c 
        d=dx*vectors(1,j)+dy*vectors(2,j)
cccc     call prin2('d=*',d,1)
        trans(j)=cdexp(ima*rk*d)
 1800 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine jexev2d(x0,y0,x,y,vectors,
     1    rk,nphi,g,fint)
        implicit real *8 (a-h,o-z)
        save
        complex *16 rk,g(nphi),ima,fint,cd,
     1    cdd
        dimension vectors(2,nphi)
         data ima/(0.0d0,1.0d0)/
c 
c       this subroutine evaluates  at a point in R^3 a
c       j-expansion given by its far-field representation
c 
c               input parameters:
c 
c  x0,y0 - coordinates of the center of the expansion
c  x,y - coordinates of the point where the expansion
c           is to be evaluated
c  vectors - the normals to the unity circle at the
c          equispaced nodes
c  rk - the helmholtz coefficient
c  nphi - the number of elements in the discretization
c        of each parallel
c  g - the far-field representation of the potential
c       to be evaluated
c 
c               output parameters:
c 
c  fint - the potential evaluated at the point (x,y,z)
c 
cccc           call prin2('in jexev2d, g=*',g,nphi*2)
        dx=x-x0
        dy=y-y0
c 
        done=1
        pi=datan(done)*4
        h=2*pi/nphi
c 
        cdd=0
        do 1800 j=1,nphi
c 
c       calculate the projection of the location of the
c       charge on the direction in question
c 
        proj=vectors(1,j)*dx+vectors(2,j)*dy
cccc      call prin2('proj=*',proj,1)
c 
c       calculate the contribution of this component
c       of the far-field
c 
        cd    =cdexp( ima*rk*proj)
        cdd=  cdd+cd * g(j)
 1800 continue
        fint=cdd *h
        fint=fint*ima/pi/2
        return
        end
c 
c 
c 
c 
c 
        subroutine hjtrans(trans,exp1,nphi,
     1     exp2)
        implicit real *8 (a-h,o-z)
        save
        complex *16 trans(nphi),
     1     exp1(nphi),exp2(nphi)
c 
        entry hhtrans(trans,exp1,nphi,
     1     exp2)
c 
        entry jjtrans(trans,exp1,nphi,
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
c        of the circle
c 
c                output parameters:
c 
c  exp2 - the shifted expansion
c 
cccc         call prin2('in hjtrans, exp1=*',exp1,nphi*2)
        do 1800 j=1,nphi
        exp2(j)=exp1(j)*trans(j)
 1800 continue
        return
        end
C 
C 
C 
C 
C 
        SUBROUTINE HBFUN(X,H,N)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        COMPLEX *16 H(1),H0,H1,Z,RK,RK7,work(1300),
     1    h07,h17
C 
C        CONSTRUCT H0,H1
C 
        CALL HANK71(X,work,H0,H1)
C 
C     . . .  AND THE REST
C 
        Z=RK*X
        I=0
        H(I) =H0
        H(1)=H1
        DO 1200 I=1,N
        H(I+1)=I*2/Z*H(I)-H(I-1)
 1200 CONTINUE
        RETURN
C 
C 
C 
        ENTRY HBFUNI(RK7)
        RK=RK7
        CALL HANI71(RK,work)
        RETURN
c 
c 
c 
c 
        entry hfun(x,h07,h17)
cccc     call prin2('inside hfun, x=*',x,1)
        CALL HANK71(X,work,H07,H17)
cccc     call prin2('inside hfun, h07=*',h07,2)
        return
        end
c 
c 
c 
c 
c 
        SUBROUTINE HBFUNl(X,H,N)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        COMPLEX *16 H(1),H0,H1,Z,RK,WORK(1300),RK7
C 
C        CONSTRUCT H0,H1
C 
        CALL HANK71(X,work,H0,H1)
C 
C     . . .  AND THE REST
C 
        Z=RK*X
        I=0
        H(I) =H0
        H(1)=H1
c 
         rmax=1.0d50
        do 1100 i=2,n
        h(i)=rmax
 1100 continue
c 
        DO 1200 I=1,N
        H(I+1)=I*2/Z*H(I)-H(I-1)
         if(cdabs(h(i+1)) .gt. rmax) goto 1400
 1200 CONTINUE
 1400 CONTINUE
        RETURN
C 
C 
C 
        ENTRY HBFUNlI(RK7)
        RK=RK7
        CALL HANI71(RK,work)
cccc    call hani41(RK)
        RETURN
        END
c 
c 
c 
c 
c 
        subroutine fndxofm(m,eps,x,funs,x01,x02)
        implicit real *8 (a-h,o-z)
        save
        dimension funs(1)
c 
c        this subroutine finds the positive real number x such that
c        for the user-specified m, eps,
c 
c           I_m (x)=eps                                                (1)
c 
c 
c          use bisection to find x solving (1)
c 
c              . . . initialize the bisection
c 
        x1=x01
        x2=x02
        eps2=eps*1.0d-17
        call ifunscom(x1,eps2,n,funs(2))
        f1=funs(m+1)
        call ifunscom(x2,eps2,n,funs(2))
        f2=funs(m+1)
c 
c       conduct the bisection
c 
        imax=60
        do 2000 i=1,imax
c 
        if(i .ne. imax) goto 1200
cccc          call prinf('in fndxofm, i=*',i,1)
cccc          call prin2('in fndxofm, x1=*',x1,1)
cccc          call prin2('in fndxofm, x2=*',x2,1)
cccc          call prin2('in fndxofm, f1=*',f1,1)
cccccc          call prin2('in fndxofm, f2=*',f2,1)
cccccc          call prin2('and f2-eps=*',f2-eps,1)
 1200 continue
c 
        x3=(x1+x2)/2
        call ifunscom(x3,eps2,n,funs(2))
        f3=funs(m+1)
c 
        if(f3 .lt. eps) goto 1400
        x2=x3
        f2=f3
        goto 2000
 1400 continue
        x1=x3
        f1=f3
 2000 continue
        x=(x2+x1)/2
        return
        end
c 
c 
c 
c 
c 
  
        subroutine ifunscom(x,eps,n,funs)
        implicit real *8 (a-h,o-z)
        save
        dimension funs(1)
c 
c       this subroutine evaluates all scaled I-type (modified)
c       Bessel functions of the argument x that are greater
c       than eps. More specifically, on exit, the i-the element
c       of the array funs is
c 
c              e^{-x} \cdot I_{i-1} (x).                            (1)
c 
c                       input parameters:
c 
c  x - the argument of which the modified Bessel functions are
c       to be evaluated.
c  eps - the threshold after which the functions (1) will be viewed
c        as zero.
c 
c                        output parameters:
c 
c  n - the number of scaled Bessel functions (1) returned. In other
c         words, e^{-x} \cdot I_{i-1} (x) < eps for all i>n.
c  funs - the array of functions (1). Note that it must have
c         element number 0 !!)
c 
c       . . . recurse up till the thing becomes fairly large
c 
        funs(1)=1
        funs(2)=0
        rlarge=1/eps*1.0d17
c 
        do 1200 i=1,1000 000
        ni=i
        funs(i+1)=-2*i/x*funs(i)+funs(i-1)
        if(funs(i+1) .gt. rlarge) goto 1400
 1200 continue
 1400 continue
cccc         call prin2('funs after firt recursion*',funs,ni)
c 
c       now, starting with the position ni, recurse down
c 
        funs(ni+1)=0
        funs(ni)=1
c 
        do 1600 i=ni,1,-1
        funs(i-1)=2*i/x*funs(i)+funs(i+1)
 1600 continue
c 
c       sum up the squares of these recursed values, to get
c       the normalizing coefficient
c 
        i0=0
        sum=funs(i0)/2
        do 1800 i=1,ni
        sum=sum+funs(i)
 1800 continue
        sum=1/(sum*2)
c 
c       . . . scale
c 
        do 2000 i=i0,ni
        funs(i)=funs(i)*sum
 2000 continue
        n=ni
        return
        end
c 
c 
c 
c 
c 
        subroutine hjmfind(n,nmin,ktwos,kthrees,kfives)
        implicit real *8 (a-h,o-z)
        save
        dimension twos(0:100),threes(0:100),
     1      fives(0:100)
c 
c       this subroutine finds the smallest integer nmin that is
c 
c    1. Greater than the user-specified number n
c    2. is a product of integer powers of three primes: 2, 3, 5.
c 
c       It also produces the three powers ktwos,kthrees,kfives
c       with which the numbers 2,3,5 enter in nmin.
c 
c        . . . determine the maximum number of primes
c              to be considered
c 
        m=0
        ii=1
        do 1200 i=1,1000
        m=m+1
        ii=ii*2
        if(ii .ge. n) goto 1400
 1200 continue
 1400 continue
c 
cccc        call prinf('m as computed is*',m,1)
c 
c       construct the tables of powers of 2, 3, 5
c 
        twos(0)=1
        threes(0)=1
        fives(0)=1
c 
        do 1500 i=1,m
        twos(i)=twos(i-1)*2
        threes(i)=threes(i-1)*3
        fives(i)=fives(i-1)*5
 1500 continue
c 
cccc        call prin2('ones=*',ones,m)
cccc        call prin2('twos=*',twos,m)
cccc        call prin2('threes=*',threes,m)
cccc        call prin2('fives=*',fives,m)
c 
c       now, test all the permitted combinations, looking for
c       the best approximation to n
c 
        ddmin=1.0d50
c 
        maxtwos=m
c 
        do 1800 ntwos=0,maxtwos
c 
        maxthrees=m-ntwos
c 
        do 1600 nthrees=0,maxthrees
c 
        maxfives=m-(nthwos+nthrees)
c 
        do 1550 nfives=0,maxfives
c 
        d=twos(ntwos)*threes(nthrees)*fives(nfives)
c 
cccc         call prin2('d=*',d,1)
        dd=d-n
        if(dd .lt. 0) goto 1550
        if(dd .ge. ddmin) goto 1550
  
        ktwos=ntwos
        kthrees=nthrees
        kfives=nfives
c 
        ddmin=dd
        dmin=d
 1550 continue
 1600 continue
 1800 continue
 2000 continue
c 
        nmin=dmin+0.1
c 
        return
        end
  
c 
