        implicit real *8 (a-h,o-z)
        real *8 ts(100 000),xs(100 000),ys(100 000),
     1      frea(100 000),fima(100 000),ts2(100 000)
c 
        complex *16 coefs(100 000),rk,ima,vals(100 000),val,
     1      ders(100 000),der,w(5000 000),hanks(100 000)
c 
        data ima/(0.0d0,1.0d0)/
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER ndigits'
        READ *,ndigits
        CALL PRINf('ndigits=*',ndigits,1 )
  
c 
c        evaluate the appropriate value of c
c 
        done=1
        pi=atan(done)*4
  
        phimax=pi/100 *2 /10
        phimax=pi/220
        phimax=pi/170
cccc        phimax=pi/85
  
cccc        phimax=pi/1000
cccc        phimax=pi/1200
  
        rk=1
        r=350
        r=700
cccc        r=5000
  
  
cccc        r=r/5 /2
cccc        phimax=phimax/3 *5 *2 /1.1
  
  
        lenw=1000 000
        call lensemake(ier,r,rk,phimax,ndigits,
     1    coefs,nterms,w,lenw,lused)
  
        call prin2('after lensemake, coefs=*',coefs,nterms*2)
        call prinf('after lensemake, lused=*',lused,1)
        call prinf('after lensemake, nterms=*',nterms,1)
  
        call prinf('after lensemake, ier=*',ier,1)
  
        if(ier .ne. 0) stop
c 
c        evaluate the beam at a bunch of nodes on the circle
c 
  
c 
  
cccc        r=r+24
  
  
  
        ntest=200
  
        done=1
        h=2*pi/ntest
c 
        do 2400 i=1,ntest
c 
        ts(i)=(i-1)*h
        xs(i)=r*cos(ts(i ))
        ys(i)=r*sin(ts(i))
c 
        call lenseeval(rk,coefs,nterms,xs(i),ys(i),vals(i),
     1      ders(i),hanks)
c 
        frea(i)=abs(vals(i))
        fima(i)=abs(ders(i))
  
  
cccc        call prinf('i=*',i,1)
cccc        call prin2('frea(i)=*',frea(i),1)
c 
 2400 continue
c 
        iw=28
        call quagraph2(iw,ts,frea,ntest,3,ts,fima,ntest,3,
     1      'bell on the circle *')
c 
c        construct and plot the beam as a function of physical width
c 
        phi0=pi/2
        aa=phi0-phimax*5
        bb=phi0+phimax*5
  
cccc        aa=phi0-phimax
cccc        bb=phi0+phimax
c 
        h=(bb-aa)/(ntest-1)
  
        amax=0
        eps=0.1**ndigits
c 
        do 2600 i=1,ntest
c 
        ts2(i)=(i-1)*h+aa
  
        xs(i)=r*cos(ts2(i))
        ys(i)=r*sin(ts2(i))
        ts2(i)=(ts2(i)-pi/2)*r
c 
        call lenseeval(rk,coefs,nterms,xs(i),ys(i),val,der,
     1      hanks)
c 
        if(abs(val)+abs(der) .gt. eps*1000) amax=ts2(i)
  
cccc        frea(i)=abs(val)
        frea(i)=log10(abs(val))
        fima(i)=log10(abs(der))
  
  
  
 2600 continue
c 
        iw=29
        call quagraph(iw,ts2,frea,ntest,3,
     1      'absolute value on the interval *')
c 
        iw=31
        call quagraph2(iw,ts2,frea,ntest,3,ts2,fima,ntest,3,
     1      'absolute values of f and df/dr *')
  
        call prin2('and amax=*',amax,1)
  
        stop
        end
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code, and the beginning
c        of the lensemaking code proper.
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c        There are two user-callable subroutines in this file:
c        lensemake and lenseeval. Following is a brief description
c        of these subroutines.
c 
c   lensemake - constructs a solution of the Helmholtz
c        equation resembling a potential focused by a lense.
c        More specifically, it constructs the coefficients of an
c        H-expansion such that the potential represented by the
c        said expansion (together with its gradient) is small on
c        the circle of radius r, except for a small spot around
c        theta=pi/2. The obtaioned expansion and its derivative
c        with respect to r can be evaluated (inter alia) by calling
c        the subroutine lenseeval (see)
c 
c   lenseeval - evaluates at the user-supplied point
c         (x,y) in R^2 the potential represented by an H-expansion
c         with the (user-supplied) coefficients coefs; the
c         derivative of the said potential with respect to the
c         radius sqrt(x**2+y**2) is also retrned.
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        subroutine lensemake(ier,r,rk,phimax,ndigits,
     1    coefs,nterms,w,lenw,lused)
c 
        implicit real *8 (a-h,o-z)
        save
        real *8 w(1)
c 
        complex *16 rk,coefs(1)
c 
c        This subroutine constructs a solution of the Helmholtz
c        equation resembling a potential focused by a lense.
c        More specifically, it constructs the coefficients of an
c        h-expansion such that the potential represented by the
c        said expansion (together with its gradient) is small on
c        the circle of radius r, except for a small spot around
c        theta=pi/2. The obtaioned expansion and its derivative
c        with respect to r can be evaluated (inter alia) by calling
c        the subroutine lenseeval (see)
c 
c        PLEASE NOTE THAT THE PERFORMANCE OF THIS SUBROUTINE DEPENDS
C        ON A CORRECT CHOICE OF THE PARAMETERS R, RK (FOR REAL RK,
C        THE ACTUAL PARAMETER BEING CHOSEN IS R*RK), phimax, ndigits;
c        THE CHOICE OF THE SAID PARAMETERS IS THE USER'S RESPONSIBILITY
c 
c                        Input parameters:
c 
c  r - the radius of the circle at the point theta=pi/2 of which
c        the lsnse is focused.
c  rk - the Helmholtz coefficient
c  phimax - the radius (in radians) pf the illuminated patch
c  ndigits - the accuracy (in numbers of digits) with which the shadow
c        is a shadow. PLEASE NOTE THAT NDIGITS IS THE ACCURACY WITH
C        WHICH THE AMPLITUDE IS ZERO. THE ENERGY DENSITY IS THE
C        SQUARE OF THAT
c  lenw - the amount of storage (in real *8 words) in the user-provided
c        work array w - make it reasonably large (1.0d10 or so).
c 
c                        Output parameters:
c 
c  ier - error return code.
c         ier=0 means (somewhat) successful execution
c         ier=8 means that the user-provided r is too small for the
c                 user-provided ndigits, phimax. Please note that
c                 SMALLER phimax require LARGER r, and LARGER ndigits
c                 require LARGER r.
c  ier=4096 means that the memory is insufficient, but in a somewhat
c                 strange manner. The subroutine is smelling some sort
c                 of rat.
c  ier=8192 means that some sort of serious trouble (memory?)
c 
c           All other ier mean insufficient memory
c  coefs - the coefficients of the H-expansion describing the beam;
c 
c  nterms - the order of expansion (1) above. PLEASE NOTE THAT THE
C         LENGTH OF ARRAY COEFS IS ACTUALLY 2*NTERMS COMPLEX *16
C         ELEMENTS!!!
c  lused - the number of elements of the work array w actually
c         used by the subroutine
c 
c 
c        determine the number of nodes to be used in the
c        discretization of the circle (in order to guarantee
c        the necessary resolution
c 
        done=1
        pi=atan(done)*4
c 
        eps=0.1**ndigits
        call provcget(ier,eps,c)
c 
        call prin2('c as chosen*',c,1)
c 
        ntest_estim=2*c/phimax
  
        call prinf('ntest_estim=*',ntest_estim,1)
c 
        ntest=4
  
        call prinf('ntest=*',ntest,1)
  
        do 1200 i=1,30
c 
        call prinf('i=*',i,1)
  
        ntest=ntest*2
  
        call prinf('ntest=*',ntest,1)
        call prinf('ntest_estim=*',ntest_estim,1)
  
  
        if(ntest .gt. ntest_estim*1.5) goto 1400
 1200 continue
 1400 continue
  
  
        call prinf('ntest=*',ntest,1)
c 
c        allocate memory
c 
        icfs2=1
        lcfs2=ntest*2+10
c 
        iwsave=icfs2+lcfs2
        lwsave=4*ntest+50
c 
        ihanks=iwsave+lwsave
        lhanks=lenw-ihanks
c 
        lw=lhanks
c 
        if(lhanks .lt. lcfs2) then
            ier=256
            return
        endif
c 
        call lensemake0(ier,r,rk,phimax,ndigits,ntest,
     1        coefs,nterms,
     2        w(iwsave),w(icfs2),w(ihanks),lw,ltot)
c 
        lused=ihanks+ltot
        if(lcfs2 .gt. ltot) lused=ihanks+lcfs2
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine lensemake0(ier,r,rk,phimax,ndigits,ntest,
     1        coefs,nterms,
     2        wsave,cfs2,hanks,lenhanks,ltot)
c 
        implicit real *8 (a-h,o-z)
        save
        complex *16 wsave(1),cfs2(1),rk,z,hanks(1),
     1      coefs(1)
c 
c        evaluate the appropriate value of c
c 
        ier=0
        lenw=lenhanks
c 
        done=1
        pi=atan(done)*4
c 
        eps=0.1**ndigits
        call provcget(ier,eps,c)
  
        call prin2('c as chosen*',c,1)
c 
c       initialize the evaluation of the prolate function
c 
        call prol0ini(jer,c,hanks,rlam20,rkhi,lenw,keep,ltot)
c 
        if(jer .ne. 0) then
            call prin2('bombing on jer from prol0ini, jer=*',jer,1)
            ier=64
            if(jer .eq. 1024) ier=4096
            if(jer .eq. 2048) ier=8192
            return
        endif
c 
c       calculate the scaled bell on the interval [-1,1]
c 
        call DCFFTI(Ntest,WSAVE)
c 
c        define the bell as the function on the interval
c        [0, 2*pi]
c 
        u=1/phimax
        v=-u*pi/2
c 
        h=2*pi/ntest
        do 3200 i=1,ntest
c 
        phi=(i-1)*h
c 
        d=u*phi+v
        if(abs(d) .gt. 1) then
            cfs2(i)=0
            goto 3200
        endif
c 
        call prol0eva(d,hanks,dd,derpsi0)
        cfs2(i)=dd
c 
 3200 continue
c 
c        Forier transform the obtained beam
c 
        call DCFFTF(Ntest,cfs2,WSAVE)
c 
c        determine the number of actual Fourier components
c        in the bell
c 
        do 3300 i=1,ntest/2
c 
        d=abs(cfs2(i))
        if(d .gt. eps*100) nterms=i
 3300 continue
c 
c        multiply the coefficients of the Fourier Transform
c        by the appropriate Hankel functions, FFT back,
c        and plot
c 
        z=r*rk
        ifexpon=1
c 
        aa=1.0d40
        call lensehanks(z,aa,hanks,nterms2,ifexpon,nterms+2)
c 
        call prinf('nterms=*',nterms,1)
        call prinf('nterms2=*',nterms2,1)
c 
        if(nterms .gt. nterms2) then
            ier=8
            return
        endif
c 
c        return to the user the coefficients of the
c        Fourier series scaled by the Hankel functionc
c 
        nn=nterms*2
        do 3800 i=1,nterms
c 
        coefs(i)=cfs2(i)/hanks(i)/ntest
 3800 continue
c 
        do 4000 i=1,nterms
c 
        coefs(nn-i+1)=cfs2(ntest-i+1) /hanks(i+1)/ntest *
     1       (-1)**i
 4000 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine lenseeval(rk,coefs,nterms,x,y,val,der,hanks)
        implicit real *8 (a-h,o-z)
        save
        complex *16 z,ima,coefs(1),val,rk,hanks(1),zz,
     1       der
c 
        data ima/(0.0d0,1.0d0)/
c 
c         This subroutine evaluates at the user-supplied point
c         (x,y) in R^2 the potential represented by an H-expansion
c         with the (user-supplied) coefficients coefs; the
c         derivative of the said potential with respect to the
c         radius sqrt(x**2+y**2) is also retrned.
c 
c 
c 
c                        Input parameters:
c 
c  rk - the Helmholtz coefficient
c  coefs - the coefficients of the H-expansion describing the beam;
c 
c  nterms - the order of expansion (1) above. PLEASE NOTE THAT THE
C         LENGTH OF ARRAY COEFS IS ACTUALLY 2*NTERMS COMPLEX *16
C         ELEMENTS!!!
c  x,y - the coordinates of the point in R^2 where the expansion
c         is to be evaluated
c 
c                        Output parameters:
c 
c  val - the value of the potential at the point (x,y)
c  der - the derivative of the potential with respect to \
c         r=sqrt(x**2+y**2)
c 
c                        Work arrays:
c 
c  hanks - must be at least nterms*8+100 real *8 locations long
  
  
c        . . . evaluate the H-expansion at the point (x,y)
c 
        r=sqrt(x**2+y**2)
        z=x/r+ima*y/r
c 
        iders=nterms+100
        zz=r*rk
c 
        ifexpon=1
        call hanks103(zz,hanks,nterms+5,ifexpon)
c 
        do 1100 i=2,nterms+1
c 
        hanks(iders+i-1)=(hanks(i-1)-hanks(i+1))/2
c 
 1100 continue
c 
        hanks(iders)=-hanks(2)
c 
        val=0
        der=0
        n=nterms*2
        do 1200 i=1,nterms-1
c 
        val=val+coefs(i)*z**(i-1) *hanks(i)
        val=val+coefs(n-i+1)/z**i *hanks(i+1) *(-1)**i
c 
        der=der+coefs(i)*z**(i-1) *hanks(iders+i-1)
        der=der+coefs(n-i+1)/z**i *hanks(iders+i) *(-1)**i
 1200 continue
c 
        der=der*rk
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine lensehanks(z,aa,hanks,nout,ifexpon,nlarge)
        implicit real *8 (a-h,o-z)
        save
        complex *16 z,hanks(1),cd,cdd
c 
        call hank103(z,hanks(1),hanks(2),ifexpon)
c 
c       conduct recursion
c 
        cd=2/z
        cdd=cd
        do 1200 i1=2,nlarge
c 
        i=i1-1
c 
        hanks(i1+1)=cdd*hanks(i1)-hanks(i1-1)
c 
        cdd=cdd+cd
        nout=i1
        if(abs(hanks(i1+1)) .gt. aa) goto 1400
c 
 1200 continue
c 
 1400 continue
c 
        return
        end
  
