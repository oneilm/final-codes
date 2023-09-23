c     This set of subroutines enabes one to do two things:
c     1) Compute the Hankel functions of orders 0 and 1
c     on the ray arg/cdabs(arg)=const in the upper half
c     of the complex arg-plane which is done by the subroutine
c     hank71; before any computations in this regime
c     are performed, array wwork must be initialized
c     by the subroutine hani71.
ccc
c     2) Compute the Hankel functions of orders 0 and 1
c     of arbitrary complex argument arg in the upper half
c     of the complex arg-plane which is done by the subroutine
c     hanvo7; before any computations in this regime
c     are performed, array wsave must be initialized
c     by the subroutine haaus7.
ccc
c     NOTE. The h-functions are complex *16 numbers.
ccc
ccc
      subroutine hank71(rxx,wwork,han0,han1)
c 
c     subroutine hank71 evaluates the Hankel functions
c     of orders 0 and 1 of arbitrary complex
c     argument arg  lying on a ray arg/cdabs(arg)=const
c     in the upper half of the complex arg-plain
c     Im(arg).ge.0; before this subroutine
c     is called the array wwork must be initialized
c     by the subroutine hani71
c 
c       note:
c 
c     this code has been designed by G. Matviyenko at Yale
c     University, Department of Computer Science. Mathematical
c     details can be found in
c 
c     Research Report YALEU/DCS/RR-903
c      May 8, 1992.
c 
c          INPUT PARAMETERS:
c 
c     rxx- real*8 value of the 'dimensional modulus'
c     of the argument; the actual argument of
c     the functions is arg=rxx*gk, where
c     gk is (generally) complex *16 dimensional
c     'Helmholtz parameter';
cc
cc
c     wwork(*)- real*8 array; its dimension must be
c     at least 700;
c 
c          OUTPUT PARAMETERS:
c 
c     han0-complex*16 value of the h-function of order 0;
cc
c     han1-complex*16 value of the h-function of order 1.
c 
      implicit real*8 (a-h,o-z)
        save
      complex*16 rk,ima,arg
      real*8 wwork(1),tn(15),han0(1),han1(1),ark(2),aarg(2),gk(2)
      equivalence (ark(1),rk), (aarg(1),arg)
      data ima /(0.0d0,1.0d0)/
      data n /10/
cc
      rgk=wwork(680)
      gk(1)=wwork(681)
      gk(2)=wwork(682)
      ro=rxx*rgk
c 
c     evaluate han0 and han1 for small arguments by means
c     of the Taylor's expansion
c 
  
      if (ro.gt.2.d0) go to 10
      dlfi=wwork(499)
      cc=wwork(498)
      ss=wwork(497)
      call hatsu7(ro,dlfi,cc,ss,han0,han1)
      return
 10   continue
c 
c     evaluate han0 and han1 for large arguments by means
c     of the asymptotic expansion
c 
      if (ro.lt.50.d0) go to 20
      x=1.d0/ro
      call haass7(x,yre0,yim0,yre1,yim1)
      go to 1000
 20   continue
c 
c     start evaluation of han0 and han1 for intermediate
c     arguments by means of the Chebychev expansion
c 
cc
         if (ro.gt.3.d0) go to 200
c 
c     a) first segment
c 
cc
         x=1.d0/ro
         xm=wwork(247)
         xp=wwork(248)
         xd=wwork(249)
         x=xm*x-xd
cc
c 
c     evaluate real and imaginary parts of han0 and han1
c     yre0,yim0,yre1,yim1
c 
         call hits17(n,x,tn,yre0,yim0,yre1,yim1)
         go to 1000
 200     continue
         if (ro.gt.5.d0) go to 300
c 
c     b) second segment
c 
         x=1.d0/ro
         xm=wwork(327)
         xp=wwork(328)
         xd=wwork(329)
         x=xm*x-xd
cc
c 
c     evaluate real and imaginary parts of han0 and han1
c     yre0,yim0,yre1,yim1
c 
         call hits27(n,x,tn,yre0,yim0,yre1,yim1)
         go to 1000
 300     continue
         if (ro.gt.9.d0) go to 400
cc
c 
c     c) third segment
c 
         x=1.d0/ro
         xm=wwork(407)
         xp=wwork(408)
         xd=wwork(409)
         x=xm*x-xd
cc
c 
c     evaluate real and imaginary parts of han0 and han1
c     yre0,yim0,yre1,yim1
c 
         call hits37(n,x,tn,yre0,yim0,yre1,yim1)
         go to 1000
 400     continue
cc
c 
c     d) fourth segment
c 
         x=1.d0/ro
         xm=wwork(487)
         xp=wwork(488)
         xd=wwork(489)
         x=xm*x-xd
cc
c 
c     evaluate real and imaginary parts of han0 and han1
c     yre0,yim0,yre1,yim1
c 
         call hits47(n,x,tn,yre0,yim0,yre1,yim1)
c 
c     finish computing of han0 and han1 of the intermediate
c     arguments
c 
 1000    continue
         sqr=1.d0/dsqrt(ro)
         aarg(1)=gk(2)*(-rxx)
         aarg(2)=gk(1)*(rxx)
cc
         rk=cdexp(arg)*sqr
cc
         han0(1)=ark(1)*yre0-ark(2)*yim0
         han0(2)=ark(2)*yre0+ark(1)*yim0
cc
         han1(1)=ark(1)*yre1-ark(2)*yim1
         han1(2)=ark(2)*yre1+ark(1)*yim1
         return
         end
ccc
ccc
      subroutine hani71(gk,wwork)
c 
c     subroutine hani71 initializes array wwork to be used
c     in the subroutine hank71
c 
c            INPUT PARAMETER:
c 
c     gk- complex*16 value of the dimensional 'Helmholtz
c     parameter';
c 
c            OUTPUT PARAMETERS:
c 
c     wwork(*)- real*8 array to be used by the subroutine hank71;
c     its dimension mast be at least 700.
c 
      implicit real*8 (a-h,o-z)
        save
      real*8 wwork(1),tn(15),fre0(15),fim0(15),fre1(15),fim1(15)
      complex*16 cd,ima,gk
      data ima /(0.0d0,1.0d0)/
cc
      rgk=cdabs(gk)
      wwork(680)=rgk
      wwork(681)=gk
      wwork(682)=(-ima)*gk
      cd=gk/rgk
c 
c     allocate the memory
c 
      iwsave=1
      iroot=151
cc
      i1pre0=171
      i1pim0=191
      i1pre1=211
      i1pim1=231
cc
      i2pre0=251
      i2pim0=271
      i2pre1=291
      i2pim1=311
cc
      i3pre0=331
      i3pim0=351
      i3pre1=371
      i3pim1=391
cc
      i4pre0=411
      i4pim0=431
      i4pre1=451
      i4pim1=471
cc
c*********************************************
      iajr0=501
      iaji0=521
      iayr0=541
      iayi0=561
cc
      iajr1=581
      iaji1=601
      iayr1=621
      iayi1=641
      n=10
      n1=n+1
c 
c     compute 'geometries' of the segments
c 
c 
c     a) 'geometry' of the first segment
c 
      ro1=3.d0
      ro2=2.d0
      x1=1.d0/ro1
      x2=1.d0/ro2
cc      wwork(245)=x1
cc      wwork(246)=x2
      wwork(247)=2.d0/(x2-x1)
      wwork(248)=x1+x2
      wwork(249)=wwork(248)/(x2-x1)
c 
c     b) 'geometry' of the second segment
c 
      ro1=5.d0
      ro2=3.d0
      x1=1.d0/ro1
      x2=1.d0/ro2
ccc      wwork(325)=x1
ccc      wwork(326)=x2
      wwork(327)=2.d0/(x2-x1)
      wwork(328)=x1+x2
      wwork(329)=wwork(328)/(x2-x1)
c 
c     c) 'geometry' of the third segment
c 
      ro1=9.d0
      ro2=5.d0
      x1=1.d0/ro1
      x2=1.d0/ro2
ccc      wwork(405)=x1
ccc      wwork(406)=x2
      wwork(407)=2.d0/(x2-x1)
      wwork(408)=x1+x2
      wwork(409)=wwork(408)/(x2-x1)
c 
c     d) 'geometry' of the fourth segment
c 
      ro1=50.d0
      ro2=9.d0
      x1=1.d0/ro1
      x2=1.d0/ro2
ccc      wwork(485)=x1
ccc      wwork(486)=x2
      wwork(487)=2.d0/(x2-x1)
      wwork(488)=x1+x2
      wwork(489)=wwork(488)/(x2-x1)
c 
c     initialize arrays root (see subroutine hachr7) and
c     wsave (see subroutine haaus7)
c 
      call hachr7(n1,wwork(iroot))
      call haaus7(wwork(iwsave))
c 
c     initialize arrays of the coefficients of the Chebychev
c     expansions for each of the four segments
c 
c     a) first segment
c 
      ro1=3.d0
      ro2=2.d0
      call haexp7(n,cd,ro1,ro2,wwork(iwsave),tn,wwork(iroot),
     *     fre0,fim0,fre1,fim1,wwork(i1pre0),
     *     wwork(i1pim0),wwork(i1pre1),wwork(i1pim1))
cc
c 
c     b) second segment
c 
      ro1=5.d0
      ro2=3.d0
      call haexp7(n,cd,ro1,ro2,wwork(iwsave),tn,wwork(iroot),
     *     fre0,fim0,fre1,fim1,wwork(i2pre0),
     *     wwork(i2pim0),wwork(i2pre1),wwork(i2pim1))
cc
c 
c     c) third segment
c 
      ro1=9.d0
      ro2=5.d0
      call haexp7(n,cd,ro1,ro2,wwork(iwsave),tn,wwork(iroot),
     *     fre0,fim0,fre1,fim1,wwork(i3pre0),
     *     wwork(i3pim0),wwork(i3pre1),wwork(i3pim1))
cc
c 
c     d) fourth segment
c 
      ro1=50.d0
      ro2=9.d0
      call haexp7(n,cd,ro1,ro2,wwork(iwsave),tn,wwork(iroot),
     *     fre0,fim0,fre1,fim1,wwork(i4pre0),
     *     wwork(i4pim0),wwork(i4pre1),wwork(i4pim1))
cc
c 
c     compute dlfi,cc and ss (see subroutine hatsu7)
c 
      wwork(499)=(-ima)*cdlog(cd)
      wwork(498)=cd
      wwork(497)=(-ima)*cd
c 
c     compute arrays to be used in the subroutine hatsu7
c 
      m=10
      call hatac7(m,wwork(1),wwork(21),cd,wwork(iajr0),wwork(iaji0),
     *        wwork(iayr0),wwork(iayi0),wwork(iajr1),
     *     wwork(iaji1),wwork(iayr1),wwork(iayi1))
ccc
c 
c     inicialize the entrees
c 
      call halex7(wwork)
      call hite17(n,wwork,tn)
      call hite27(n,wwork,tn)
      call hite37(n,wwork,tn)
      call hite47(n,wwork,tn)
      call haasy7(wwork,cd)
cc
      return
      end
ccc
ccc
      subroutine haexp7(n,cd,ro1,ro2,wsave,tn,root,
     *     fre0,fim0,fre1,fim1,pre0,pim0,pre1,pim1)
c 
c     subroutine haexp7 evaluates coefficients in Chebyshev
c     expansion of the Hankel functions of orders 0 and 1;
c     these expansions are valid on a given segment of
c     a given ray (arg/cdabs(arg)=const) lying in the upper
c     half of the complex arg-plane (Im(arg).ge.0)
c 
c          INPUT PARAMETERS:
c 
c     n- integer order of the Chebychev approximation;
cc
c     cd- complex*16 constant; cd=arg/cdabs(arg);
cc
c     ro1,ro2- real*8 coordinates of the ends of the segment
c     of the ray;
cc
c     wsave(*),tn(*),fre0(*),fim0(*),fre1(*),fim1(*)- real*8
c     work arrays;
cc
c     root(*)- real*8 array of the roots of the Chebychev
c     polinomial of order n+1 (see subroutine hachr7);
c 
c 
c           OUTPUT PARAMETERS:
c 
c     pre0(*)- real*8 array of the Chebychev coefficients of
c     the real part of the smooth factor of the Hankel function
c     of order 0; its dimension must be at least n+1;
cc
c     pim0(*)- real*8 array of the Chebychev coefficients of
c     the imaginary part of the smooth factor
c     of the Hankel function of order 0; its
c     dimension must be at least n+1;
cc
c     pre1(*)- real*8 array of the Chebychev coefficients of
c     the real part of the smooth factor
c     of the hankel function of order 1; its
c     dimension must be at least n+1;
cc
c     pim1(*)- real*8 array of the Chebychev coefficients of
c     the imaginary part of the smooth factor
c     of the Hankel function of order 1; its
c     dimension must be at least n+1;
cc
c     NOTES. 1. The expansion is sought on the interval [1/ro1,1/ro2].
c            2. By definition the smooth factor=
c               (Hankel function)*cdexp((-ima)*arg)*cdsqrt(arg)
c 
      implicit real*8 (a-h,o-z)
        save
      real*8 tn(1),root(1),fre0(1),fim0(1),fre1(1),fim1(1),
     *     pre0(1),pim0(1),pre1(1),pim1(1),wsave(1)
      complex*16 ima,arg,rk,cd,han0,han1
      data ima /(0.0d0,1.0d0)/
ccc      call prin2('cd=*',cd,2)
      n1=n+1
cc
c 
c     compute ends of the interval
c 
      x1=1.d0/ro1
      x2=1.d0/ro2
c 
c     evaluate real and imaginary parts of the functions
c     to be approximated
c 
      do 150 i=1,n1
c 
c     map the interval [-1,1] on the interval [1/ro1,1/ro2]
c 
         x=0.5d0*((x2-x1)*root(i)+x1+x2)
c 
c     compute actual value of the argument
c 
         ro=1.d0/x
         arg=ro*cd
c 
c     compute (almost) inverse of the dominant term
c     of the asymptotic expansion
c 
         rk=cdexp((-ima)*arg)*cdsqrt(arg)/cdsqrt(cd)
c 
c     compute real and imaginary parts of the
c     smooth factors of han0 and han1
c 
         call hanvo7(wsave,arg,han0,han1)
         fre0(i)=rk*han0
         fim0(i)=(-ima)*rk*han0
         fre1(i)=rk*han1
         fim1(i)=(-ima)*rk*han1
 150  continue
c 
c     evaluate the coefficients
c 
      call hachc7(n,fre0,tn,root,pre0)
      call hachc7(n,fim0,tn,root,pim0)
      call hachc7(n,fre1,tn,root,pre1)
      call hachc7(n,fim1,tn,root,pim1)
      return
      end
ccc
ccc
      subroutine hachp7(x,n,tn)
c 
c     subroutine hachp7 evaluates chebychev polinomials
c     of the argument x
c 
c          INPUT PARAMETERS:
c 
c     x- real*8 argument
c     n- integer index of the last computed polinomial
c 
c          OUTPUT PARAMETER:
c 
c     tn(*)- real*8 array of the values of the polinomials;
c            its dimension must be at least n+1
c 
c      NOTE. tn(i+1)=Ti
c 
      implicit real*8 (a-h,o-z)
        save
      real*8 tn(1)
      tn(1)=1.d0
      tn(2)=x
      y=2.d0*x
      do 100 i=2,n
         tn(i+1)=y*tn(i)-tn(i-1)
 100  continue
      return
      end
ccc
ccc
      subroutine hachr7(n1,root)
c 
c     subroutine hachr7 evaluates roots of Chebychev polinomials
c 
c 
c             INPUT PARAMETERS:
c 
c     n1- integer order of the polinomial
c 
c             OUTPUT PARAMETERS:
c 
c     root(*)- real*8 array of the roots
c 
      implicit real*8 (a-h,o-z)
        save
      real*8 root(1)
      pi=3.14159265358979323846d0
      x=pi/(2.d0*n1)
      do 100 i=1,n1
         root(i)=dcos((2*i-1)*x)
 100  continue
      return
      end
ccc
ccc
      subroutine hachc7(n,fun,tn,root,ach)
c 
c     subroutine hachc7 evaluates coefficients of the
c     Chebychev expansion of order n of a given function
c     fun that was computed at the roots of the Chebychev
c     polinomial of order n+1
c 
c 
c          INPUT PARAMETERS:
c 
c     n- integer order of approximation
c     fun(*)- real*8 array of the values of the function;
c     tn(*),root(*)- real*8 work arrays; their dimension
c     must be at least n+1;
c 
c          OUTPUT PARAMETERS:
c 
c     ach(*)- real*8 array of the expansion coefficients;
c     its dimension must be at least n+1
c 
c     NOTE. ach(i+1)='weight' of the polinomial Ti in the
c     expansion
c 
      implicit real*8 (a-h,o-z)
        save
      real*8 fun(1),tn(1),root(1),ach(1)
      n1=n+1
c 
c     evaluate nonnormalized coefficients
c 
      do 100 i=1,n1
         ach(i)=0.d0
 100  continue
      do 300 i=1,n1
         x=root(i)
         f=fun(i)
         call hachp7(x,n,tn)
         do 200 j=1,n1
            ach(j)=ach(j)+f*tn(j)
 200     continue
 300  continue
c 
c     normalize the result
c 
      x=2.d0/n1
      ach(1)=ach(1)*x/2
      do 400 i=2,n1
         ach(i)=ach(i)*x
 400  continue
      return
      end
ccc
ccc
ccc
ccc
         subroutine hatac7(m,gam,psi,cd,ajr0,aji0,
     *        ayr0,ayi0,ajr1,aji1,ayr1,ayi1)
c 
c     subroutine hatac7 evaluates  coefficients of Taylor
c     expansions of real and imaginary parts of analitic
c     'portions' of  h-functions of orders 0 and 1
c 
c              INPUT PARAMETERS:
c 
c     m- integer number of terms in the expansion;
c 
c     gam(*),psi(*)- real*8 arrays (see subroutine haaux7);
c 
c     cd- complex*16 number; cd=arg/cdabs(arg).
c 
c          OUTPUT PARAMETERS:
c 
c     ajr0(*)- real*8 array of the Taylor coefficients of
c     the real part of the j-function of order 0;
cc
c     aji0(*)- real*8 array of the Taylor coefficients of
c     the imaginary part of the j-function of order 0;
cc
c     ayr0(*)- real*8 array of the Taylor coefficients of
c     the real part of the analitic 'portion' of the
c     y-function of order 0;
cc
c     ayi0(*)- real*8 array of the Taylor coefficients of
c     the imaginary part of the analitic 'portion' of the
c     y-function of order 0;
cccc
c     ajr1(*)- real*8 array of the Taylor coefficients of
c     the real part of the j-function of order 1;
cc
c     aji1(*)- real*8 array of the Taylor coefficients of
c     the imaginary part of the j-function of order 1;
cc
c     ayr1(*)- real*8 array of the Taylor coefficients of
c     the real part of the analitic 'portion' of the
c     y-function of order 1;
cc
c     ayi1(*)- real*8 array of the Taylor coefficients of
c     the imaginary part of the analitic 'portion' of the
c     y-function of order 1;
cc
         implicit real*8 (a-h,o-z)
        save
         complex*16 cd,ima,z,z1
         real*8 gam(1),psi(1),ajr0(1),aji0(1),
     *        ayr0(1),ayi0(1),ajr1(1),aji1(1),ayr1(1),ayi1(1)
         data ima /(0.0d0,1.0d0)/
         pi=3.14159265358979323846d0
cc
         m1=m+1
         z=1.d0
         z1=(-0.25d0)*cd**2
c 
c     evaluate the arrays
c 
         do 100 i=1,m1
            ajr0(i)=z/(gam(i)*gam(i))
            aji0(i)=(-ima)*z/(gam(i)*gam(i))
            ayr0(i)=2.d0/pi*z*psi(i)/(gam(i)*gam(i))
            ayi0(i)=(-ima)*2.d0/pi*z*psi(i)/(gam(i)*gam(i))
cc
            ajr1(i)=z/(gam(i)*gam(i+1))*0.5d0*cd
            aji1(i)=(-ima)*z/(gam(i)*gam(i+1))*0.5d0*cd
            ayr1(i)=z*(psi(i)+psi(i+1))/(pi*gam(i)
     *           *gam(i+1))*0.5d0*cd
            ayi1(i)=(-ima)*z*(psi(i)+psi(i+1))/(pi*gam(i)
     *           *gam(i+1))*0.5d0*cd
cc
            z=z*z1
 100     continue
         return
         end
ccc
ccc
         subroutine hatsu7(ro,dlfi,cc,ss,han0,han1)
c 
c     subroutine hatsu7 evaluates han0 and han1 for
c     cdabs(arg).lt.2.d0 summing Taylor series
c     on the ray
c 
c 
c               INPUT PARAMETERS:
c 
c     m- integer number of terms in the expansion;
cc
c     arg- complex*16 argument of the h-functions;
cc
c     ro- real*8 value of the modulus of arg; ro=cdabs(arg);
cc
c     dlfi- real*8 value of the argument of arg;
c     dlfi=(-ima)*cdlog(arg/ro)
cc
c     ajr0(*),aji0(*),ayr0(*),ayi0(*),ajr1(*),aji1(*),
c     aayr1(*),ayi1(*)- real*8 arrays (see subroutine hatac7);
c     their dimension must be at least m+1;
c 
c              OUTPUT PARAMETERTS:
c 
cc
c     han0- complex*16 value of the h-function of order 0;
cc
c     han1- complex*16 value of the h-function of order 1.
c 
c     NOTE. Before the subroutine is called the arrays in
c     its calling sequence must be initialized by the
c     subroutine hatac7
c 
         implicit real*8 (a-h,o-z)
        save
         real*8 han0(2),han1(2)
         data ima /(0.0d0,1.0d0)/
         data pi /3.14159265358979323846d0/
cc
c 
c     evaluate the sums related to han0 and han1
c 
         ro2=ro**2
         call hatal7(ro2,xr0,xi0,yr0,yi0,xr1,xi1,yr1,yi1)
ccc
c 
c     start computing of  han0 and han1
c 
         con=2.d0/pi
         dl=dlog(0.5d0*ro)
         ro22=con/ro2
cc
c 
c     compute han0
c 
cc
         yr=con*(xr0*dl-xi0*dlfi)-yr0
         yi=con*(xi0*dl+xr0*dlfi)-yi0
         han0(1)=xr0-yi
         han0(2)=(xi0+yr)
c 
c     compute han1
c 
         yr=con*(xr1*dl-xi1*dlfi)-yr1
         yi=con*(xi1*dl+xr1*dlfi)-yi1
         han1(1)=ro*(xr1-yi-ro22*ss)
         han1(2)=ro*(xi1+yr-ro22*cc)
cc
         return
         end
ccc
ccc
         subroutine halex7(wwsave)
c 
c            INPUT PARAMETERS:
c     wwsave(*)- real*8 array of dimension at least 700;
c     it must be initialized by the subroutine hani71 where
c     it has the name wwork
c 
c     subroutine halex7 is used to obtain Taylor's coefficients
c     of the expansions of the h-functions on the ray; its main
c     part is entry hatal7 that sums up those series explicitly
c     in assembler-like manner
c 
         implicit real*8 (a-h,o-z)
        save
         real*8 wwsave(1)
ccc
         a01=wwsave(501)
         a02=wwsave(502)
         a03=wwsave(503)
         a04=wwsave(504)
         a05=wwsave(505)
         a06=wwsave(506)
         a07=wwsave(507)
         a08=wwsave(508)
         a09=wwsave(509)
         a010=wwsave(510)
         a011=wwsave(511)
ccc
         b01=wwsave(521)
         b02=wwsave(522)
         b03=wwsave(523)
         b04=wwsave(524)
         b05=wwsave(525)
         b06=wwsave(526)
         b07=wwsave(527)
         b08=wwsave(528)
         b09=wwsave(529)
         b010=wwsave(530)
         b011=wwsave(531)
ccc
         c01=wwsave(541)
         c02=wwsave(542)
         c03=wwsave(543)
         c04=wwsave(544)
         c05=wwsave(545)
         c06=wwsave(546)
         c07=wwsave(547)
         c08=wwsave(548)
         c09=wwsave(549)
         c010=wwsave(550)
         c011=wwsave(551)
ccc
         d01=wwsave(561)
         d02=wwsave(562)
         d03=wwsave(563)
         d04=wwsave(564)
         d05=wwsave(565)
         d06=wwsave(566)
         d07=wwsave(567)
         d08=wwsave(568)
         d09=wwsave(569)
         d010=wwsave(570)
         d011=wwsave(571)
ccc****************************
ccc****************************
         a11=wwsave(581)
         a12=wwsave(582)
         a13=wwsave(583)
         a14=wwsave(584)
         a15=wwsave(585)
         a16=wwsave(586)
         a17=wwsave(587)
         a18=wwsave(588)
         a19=wwsave(589)
         a110=wwsave(590)
         a111=wwsave(591)
ccc
         b11=wwsave(601)
         b12=wwsave(602)
         b13=wwsave(603)
         b14=wwsave(604)
         b15=wwsave(605)
         b16=wwsave(606)
         b17=wwsave(607)
         b18=wwsave(608)
         b19=wwsave(609)
         b110=wwsave(610)
         b111=wwsave(611)
ccc
         c11=wwsave(621)
         c12=wwsave(622)
         c13=wwsave(623)
         c14=wwsave(624)
         c15=wwsave(625)
         c16=wwsave(626)
         c17=wwsave(627)
         c18=wwsave(628)
         c19=wwsave(629)
         c110=wwsave(630)
         c111=wwsave(631)
ccc
         d11=wwsave(641)
         d12=wwsave(642)
         d13=wwsave(643)
         d14=wwsave(644)
         d15=wwsave(645)
         d16=wwsave(646)
         d17=wwsave(647)
         d18=wwsave(648)
         d19=wwsave(649)
         d110=wwsave(650)
         d111=wwsave(651)
         return
ccc
ccc
         entry hatal7(ro2,xr0,xi0,yr0,yi0,xr1,xi1,yr1,yi1)
         xr0=((((((((((a011*ro2+a010)*ro2+a09)*ro2
     *        +a08)*ro2+a07)*ro2+a06)*ro2+a05)*ro2
     *        +a04)*ro2+a03)*ro2+a02)*ro2+a01)
ccc
         xi0=((((((((((b011*ro2+b010)*ro2+b09)*ro2
     *        +b08)*ro2+b07)*ro2+b06)*ro2+b05)*ro2
     *        +b04)*ro2+b03)*ro2+b02)*ro2+b01)
ccc
         yr0=((((((((((c011*ro2+c010)*ro2+c09)*ro2
     *        +c08)*ro2+c07)*ro2+c06)*ro2+c05)*ro2
     *        +c04)*ro2+c03)*ro2+c02)*ro2+c01)
ccc
         yi0=((((((((((d011*ro2+d010)*ro2+d09)*ro2
     *        +d08)*ro2+d07)*ro2+d06)*ro2+d05)*ro2
     *        +d04)*ro2+d03)*ro2+d02)*ro2+d01)
ccc**********************************************
ccc
         xr1=((((((((((a111*ro2+a110)*ro2+a19)*ro2
     *        +a18)*ro2+a17)*ro2+a16)*ro2+a15)*ro2
     *        +a14)*ro2+a13)*ro2+a12)*ro2+a11)
ccc
         xi1=((((((((((b111*ro2+b110)*ro2+b19)*ro2
     *        +b18)*ro2+b17)*ro2+b16)*ro2+b15)*ro2
     *        +b14)*ro2+b13)*ro2+b12)*ro2+b11)
ccc
         yr1=((((((((((c111*ro2+c110)*ro2+c19)*ro2
     *        +c18)*ro2+c17)*ro2+c16)*ro2+c15)*ro2
     *        +c14)*ro2+c13)*ro2+c12)*ro2+c11)
ccc
         yi1=((((((((((d111*ro2+d110)*ro2+d19)*ro2
     *        +d18)*ro2+d17)*ro2+d16)*ro2+d15)*ro2
     *        +d14)*ro2+d13)*ro2+d12)*ro2+d11)
ccc**********************************************
         return
         end
ccc
ccc
         subroutine hite17(n,wwsave,tn)
c 
c     subroutine hite17 is used to obtain Chebychev's
c     coefficients of the expansions of the h-functions
c     on the first segment of the ray; its main
c     part is entry hits17 that sums up those series explicitly
c     in assembler-like manner
c 
c             INPUT PARAMETERS:
c 
c     n-integer order of the Chebychev approximation
c 
c     wwsave(*)- real *8 array equivalent to the array wwork in
c     the subroutine hani71
c 
c     th(*)- real *8 work array of dimension at least 15
c 
         implicit real*8 (a-h,o-z)
        save
         real*8 tn(15),wwsave(1)
ccc
         a01=wwsave(171)
         a02=wwsave(172)
         a03=wwsave(173)
         a04=wwsave(174)
         a05=wwsave(175)
         a06=wwsave(176)
         a07=wwsave(177)
         a08=wwsave(178)
         a09=wwsave(179)
         a010=wwsave(180)
         a011=wwsave(181)
ccc
         b01=wwsave(191)
         b02=wwsave(192)
         b03=wwsave(193)
         b04=wwsave(194)
         b05=wwsave(195)
         b06=wwsave(196)
         b07=wwsave(197)
         b08=wwsave(198)
         b09=wwsave(199)
         b010=wwsave(200)
         b011=wwsave(201)
ccc
         c01=wwsave(211)
         c02=wwsave(212)
         c03=wwsave(213)
         c04=wwsave(214)
         c05=wwsave(215)
         c06=wwsave(216)
         c07=wwsave(217)
         c08=wwsave(218)
         c09=wwsave(219)
         c010=wwsave(220)
         c011=wwsave(221)
ccc
         d01=wwsave(231)
         d02=wwsave(232)
         d03=wwsave(233)
         d04=wwsave(234)
         d05=wwsave(235)
         d06=wwsave(236)
         d07=wwsave(237)
         d08=wwsave(238)
         d09=wwsave(239)
         d010=wwsave(240)
         d011=wwsave(241)
ccc****************************
ccc****************************
         return
cc
         entry hits17(n,x,tn,yre0,yim0,yre1,yim1)
         call hachp7(x,n,tn)
         yre0=a01*tn(1)+a02*tn(2)+a03*tn(3)+a04*tn(4)+a05*tn(5)+
     *        a06*tn(6)+a07*tn(7)+a08*tn(8)+a09*tn(9)+a010*tn(10)+
     *        a011*tn(11)
ccc
         yim0=b01*tn(1)+b02*tn(2)+b03*tn(3)+b04*tn(4)+b05*tn(5)+
     *        b06*tn(6)+b07*tn(7)+b08*tn(8)+b09*tn(9)+b010*tn(10)+
     *        b011*tn(11)
ccc
         yre1=c01*tn(1)+c02*tn(2)+c03*tn(3)+c04*tn(4)+c05*tn(5)+
     *        c06*tn(6)+c07*tn(7)+c08*tn(8)+c09*tn(9)+c010*tn(10)+
     *        c011*tn(11)
ccc
         yim1=d01*tn(1)+d02*tn(2)+d03*tn(3)+d04*tn(4)+d05*tn(5)+
     *        d06*tn(6)+d07*tn(7)+d08*tn(8)+d09*tn(9)+d010*tn(10)+
     *        d011*tn(11)
         return
         end
ccc
ccc
         subroutine hite27(n,wwsave,tn)
c 
c     subroutine hite27 is used to obtain Chebychev's
c     coefficients of the expansions of the h-functions
c     on the second segment of the ray; its main
c     part is entry hits27 that sums up those series explicitly
c     in assembler-like manner
c 
c 
c             INPUT PARAMETERS:
c 
c     n-integer order of the Chebychev approximation
c 
c     wwsave(*)- real *8 array equivalent to the array wwork in
c     the subroutine hani71
c 
c     th(*)- real *8 work array of dimension at least 15
c 
         implicit real*8 (a-h,o-z)
        save
         real*8 tn(15),wwsave(1)
ccc
         a01=wwsave(251)
         a02=wwsave(252)
         a03=wwsave(253)
         a04=wwsave(254)
         a05=wwsave(255)
         a06=wwsave(256)
         a07=wwsave(257)
         a08=wwsave(258)
         a09=wwsave(259)
         a010=wwsave(260)
         a011=wwsave(261)
ccc
         b01=wwsave(271)
         b02=wwsave(272)
         b03=wwsave(273)
         b04=wwsave(274)
         b05=wwsave(275)
         b06=wwsave(276)
         b07=wwsave(277)
         b08=wwsave(278)
         b09=wwsave(279)
         b010=wwsave(280)
         b011=wwsave(281)
ccc
         c01=wwsave(291)
         c02=wwsave(292)
         c03=wwsave(293)
         c04=wwsave(294)
         c05=wwsave(295)
         c06=wwsave(296)
         c07=wwsave(297)
         c08=wwsave(298)
         c09=wwsave(299)
         c010=wwsave(300)
         c011=wwsave(301)
ccc
         d01=wwsave(311)
         d02=wwsave(312)
         d03=wwsave(313)
         d04=wwsave(314)
         d05=wwsave(315)
         d06=wwsave(316)
         d07=wwsave(317)
         d08=wwsave(318)
         d09=wwsave(319)
         d010=wwsave(320)
         d011=wwsave(321)
ccc****************************
ccc****************************
         return
cc
         entry hits27(n,x,tn,yre0,yim0,yre1,yim1)
         call hachp7(x,n,tn)
         yre0=a01*tn(1)+a02*tn(2)+a03*tn(3)+a04*tn(4)+a05*tn(5)+
     *        a06*tn(6)+a07*tn(7)+a08*tn(8)+a09*tn(9)+a010*tn(10)+
     *        a011*tn(11)
ccc
         yim0=b01*tn(1)+b02*tn(2)+b03*tn(3)+b04*tn(4)+b05*tn(5)+
     *        b06*tn(6)+b07*tn(7)+b08*tn(8)+b09*tn(9)+b010*tn(10)+
     *        b011*tn(11)
ccc
         yre1=c01*tn(1)+c02*tn(2)+c03*tn(3)+c04*tn(4)+c05*tn(5)+
     *        c06*tn(6)+c07*tn(7)+c08*tn(8)+c09*tn(9)+c010*tn(10)+
     *        c011*tn(11)
ccc
         yim1=d01*tn(1)+d02*tn(2)+d03*tn(3)+d04*tn(4)+d05*tn(5)+
     *        d06*tn(6)+d07*tn(7)+d08*tn(8)+d09*tn(9)+d010*tn(10)+
     *        d011*tn(11)
         return
         end
ccc
ccc
         subroutine hite37(n,wwsave,tn)
c 
c             INPUT PARAMETERS:
c 
c     n-integer order of the Chebychev approximation
c 
c     wwsave(*)- real *8 array equivalent to the array wwork in
c     the subroutine hani71
c 
c     th(*)- real *8 work array of dimension at least 15
c 
c 
c     subroutine hite37 is used to obtain Chebychev's
c     coefficients of the expansions of the h-functions
c     on the third segment of the ray; its main
c     part is entry hits37 that sums up those series explicitly
c     in assembler-like manner
c 
         implicit real*8 (a-h,o-z)
        save
         real*8 tn(15),wwsave(1)
ccc
         a01=wwsave(331)
         a02=wwsave(332)
         a03=wwsave(333)
         a04=wwsave(334)
         a05=wwsave(335)
         a06=wwsave(336)
         a07=wwsave(337)
         a08=wwsave(338)
         a09=wwsave(339)
         a010=wwsave(340)
         a011=wwsave(341)
ccc
         b01=wwsave(351)
         b02=wwsave(352)
         b03=wwsave(353)
         b04=wwsave(354)
         b05=wwsave(355)
         b06=wwsave(356)
         b07=wwsave(357)
         b08=wwsave(358)
         b09=wwsave(359)
         b010=wwsave(360)
         b011=wwsave(361)
ccc
         c01=wwsave(371)
         c02=wwsave(372)
         c03=wwsave(373)
         c04=wwsave(374)
         c05=wwsave(375)
         c06=wwsave(376)
         c07=wwsave(377)
         c08=wwsave(378)
         c09=wwsave(379)
         c010=wwsave(380)
         c011=wwsave(381)
ccc
         d01=wwsave(391)
         d02=wwsave(392)
         d03=wwsave(393)
         d04=wwsave(394)
         d05=wwsave(395)
         d06=wwsave(396)
         d07=wwsave(397)
         d08=wwsave(398)
         d09=wwsave(399)
         d010=wwsave(400)
         d011=wwsave(401)
ccc****************************
ccc****************************
         return
cc
         entry hits37(n,x,tn,yre0,yim0,yre1,yim1)
         call hachp7(x,n,tn)
         yre0=a01*tn(1)+a02*tn(2)+a03*tn(3)+a04*tn(4)+a05*tn(5)+
     *        a06*tn(6)+a07*tn(7)+a08*tn(8)+a09*tn(9)+a010*tn(10)+
     *        a011*tn(11)
ccc
         yim0=b01*tn(1)+b02*tn(2)+b03*tn(3)+b04*tn(4)+b05*tn(5)+
     *        b06*tn(6)+b07*tn(7)+b08*tn(8)+b09*tn(9)+b010*tn(10)+
     *        b011*tn(11)
ccc
         yre1=c01*tn(1)+c02*tn(2)+c03*tn(3)+c04*tn(4)+c05*tn(5)+
     *        c06*tn(6)+c07*tn(7)+c08*tn(8)+c09*tn(9)+c010*tn(10)+
     *        c011*tn(11)
ccc
         yim1=d01*tn(1)+d02*tn(2)+d03*tn(3)+d04*tn(4)+d05*tn(5)+
     *        d06*tn(6)+d07*tn(7)+d08*tn(8)+d09*tn(9)+d010*tn(10)+
     *        d011*tn(11)
         return
         end
ccc
ccc
ccc
ccc
         subroutine hite47(n,wwsave,tn)
c 
c     subroutine hite47 is used to obtain Chebychev's
c     coefficients of the expansions of the h-functions
c     on the fourth segment of the ray; its main
c     part is entry hits47 that sums up those series explicitly
c     in assembler-like manner
c 
c 
c             INPUT PARAMETERS:
c 
c     n-integer order of the Chebychev approximation
c 
c     wwsave(*)- real *8 array equivalent to the array wwork in
c     the subroutine hani71
c 
c     th(*)- real *8 work array of dimension at least 15
c 
         implicit real*8 (a-h,o-z)
        save
         real*8 tn(15),wwsave(1)
ccc
         a01=wwsave(411)
         a02=wwsave(412)
         a03=wwsave(413)
         a04=wwsave(414)
         a05=wwsave(415)
         a06=wwsave(416)
         a07=wwsave(417)
         a08=wwsave(418)
         a09=wwsave(419)
         a010=wwsave(420)
         a011=wwsave(421)
ccc
         b01=wwsave(431)
         b02=wwsave(432)
         b03=wwsave(433)
         b04=wwsave(434)
         b05=wwsave(435)
         b06=wwsave(436)
         b07=wwsave(437)
         b08=wwsave(438)
         b09=wwsave(439)
         b010=wwsave(440)
         b011=wwsave(441)
ccc
         c01=wwsave(451)
         c02=wwsave(452)
         c03=wwsave(453)
         c04=wwsave(454)
         c05=wwsave(455)
         c06=wwsave(456)
         c07=wwsave(457)
         c08=wwsave(458)
         c09=wwsave(459)
         c010=wwsave(460)
         c011=wwsave(461)
ccc
         d01=wwsave(471)
         d02=wwsave(472)
         d03=wwsave(473)
         d04=wwsave(474)
         d05=wwsave(475)
         d06=wwsave(476)
         d07=wwsave(477)
         d08=wwsave(478)
         d09=wwsave(479)
         d010=wwsave(480)
         d011=wwsave(481)
ccc****************************
ccc****************************
         return
cc
         entry hits47(n,x,tn,yre0,yim0,yre1,yim1)
         call hachp7(x,n,tn)
         yre0=a01*tn(1)+a02*tn(2)+a03*tn(3)+a04*tn(4)+a05*tn(5)+
     *        a06*tn(6)+a07*tn(7)+a08*tn(8)+a09*tn(9)+a010*tn(10)+
     *        a011*tn(11)
ccc
         yim0=b01*tn(1)+b02*tn(2)+b03*tn(3)+b04*tn(4)+b05*tn(5)+
     *        b06*tn(6)+b07*tn(7)+b08*tn(8)+b09*tn(9)+b010*tn(10)+
     *        b011*tn(11)
ccc
         yre1=c01*tn(1)+c02*tn(2)+c03*tn(3)+c04*tn(4)+c05*tn(5)+
     *        c06*tn(6)+c07*tn(7)+c08*tn(8)+c09*tn(9)+c010*tn(10)+
     *        c011*tn(11)
ccc
         yim1=d01*tn(1)+d02*tn(2)+d03*tn(3)+d04*tn(4)+d05*tn(5)+
     *        d06*tn(6)+d07*tn(7)+d08*tn(8)+d09*tn(9)+d010*tn(10)+
     *        d011*tn(11)
         return
         end
ccc
ccc
         subroutine haasy7(wwsave,cd)
c 
c     subroutine haasy7 is used to obtain coefficients
c     of the asymptotic expansions of the h-functions
c     on the ray; its main part is entry haass7 that sums up
c     those series explicitly in assembler-like manner
c 
c 
c 
c             INPUT PARAMETERS:
c 
c     wwsave(*)- real *8 array equivalent to the array wwork in
c     the subroutine hani71
c 
c     cd- complex *16 variable cd=arg/cdabs(arg), where arg is
c     the actual argument of the h-functions
c 
         implicit real *8 (a-h,o-z)
        save
         complex *16 cd,ccd,ima,d00,d11
         complex *16  pwi(20)
         real *8 wwsave(1)
         data ima /(0.d0,1.d0)/
         pi=3.14159265358979323846d0
cc
         d00=(1.d0-ima)/cdsqrt(pi*cd)
         d11=(-ima)*d00
cc
         ccd=ima/cd
         pwi(1)=ccd
         do 100 i=1,14
            pwi(i+1)=ccd*pwi(i)
 100     continue
cc
         a01=wwsave(41)*d00*pwi(1)
         a02=wwsave(42)*d00*pwi(2)
         a03=wwsave(43)*d00*pwi(3)
         a04=wwsave(44)*d00*pwi(4)
         a05=wwsave(45)*d00*pwi(5)
         a06=wwsave(46)*d00*pwi(6)
         a07=wwsave(47)*d00*pwi(7)
         a08=wwsave(48)*d00*pwi(8)
         a09=wwsave(49)*d00*pwi(9)
         a010=wwsave(50)*d00*pwi(10)
         ac=d00
cc
cc
         b01=wwsave(41)*d00*pwi(1)*(-ima)
         b02=wwsave(42)*d00*pwi(2)*(-ima)
         b03=wwsave(43)*d00*pwi(3)*(-ima)
         b04=wwsave(44)*d00*pwi(4)*(-ima)
         b05=wwsave(45)*d00*pwi(5)*(-ima)
         b06=wwsave(46)*d00*pwi(6)*(-ima)
         b07=wwsave(47)*d00*pwi(7)*(-ima)
         b08=wwsave(48)*d00*pwi(8)*(-ima)
         b09=wwsave(49)*d00*pwi(9)*(-ima)
         b010=wwsave(50)*d00*pwi(10)*(-ima)
         bc=(-ima)*d00
ccc
         c01=wwsave(91)*d11*pwi(1)
         c02=wwsave(92)*d11*pwi(2)
         c03=wwsave(93)*d11*pwi(3)
         c04=wwsave(94)*d11*pwi(4)
         c05=wwsave(95)*d11*pwi(5)
         c06=wwsave(96)*d11*pwi(6)
         c07=wwsave(97)*d11*pwi(7)
         c08=wwsave(98)*d11*pwi(8)
         c09=wwsave(99)*d11*pwi(9)
         c010=wwsave(100)*d11*pwi(10)
         cc=d11
ccc
         d01=wwsave(91)*d11*pwi(1)*(-ima)
         d02=wwsave(92)*d11*pwi(2)*(-ima)
         d03=wwsave(93)*d11*pwi(3)*(-ima)
         d04=wwsave(94)*d11*pwi(4)*(-ima)
         d05=wwsave(95)*d11*pwi(5)*(-ima)
         d06=wwsave(96)*d11*pwi(6)*(-ima)
         d07=wwsave(97)*d11*pwi(7)*(-ima)
         d08=wwsave(98)*d11*pwi(8)*(-ima)
         d09=wwsave(99)*d11*pwi(9)*(-ima)
         d010=wwsave(100)*d11*pwi(10)*(-ima)
         dc=d11*(-ima)
ccc****************************
         return
ccc
         entry haass7(ro2,yre0,yim0,yre1,yim1)
         yre0=((((((((((a010*ro2+a09)*ro2+a08)*ro2
     *        +a07)*ro2+a06)*ro2+a05)*ro2+a04)*ro2
     *        +a03)*ro2+a02)*ro2+a01)*ro2+ac)
ccc
         yim0=((((((((((b010*ro2+b09)*ro2+b08)*ro2
     *        +b07)*ro2+b06)*ro2+b05)*ro2+b04)*ro2
     *        +b03)*ro2+b02)*ro2+b01)*ro2+bc)
ccc
         yre1=((((((((((c010*ro2+c09)*ro2+c08)*ro2
     *        +c07)*ro2+c06)*ro2+c05)*ro2+c04)*ro2
     *        +c03)*ro2+c02)*ro2+c01)*ro2+cc)
ccc
         yim1=((((((((((d010*ro2+d09)*ro2+d08)*ro2
     *        +d07)*ro2+d06)*ro2+d05)*ro2+d04)*ro2
     *        +d03)*ro2+d02)*ro2+d01)*ro2+dc)
ccc**********************************************
ccc**********************************************
         return
         end
ccc
ccc
ccc
      subroutine hanvo7(wsave,arg,ha0,ha1)
c 
c     subroutine hanvo7 evaluates Hankel functions of
c     orders 0 and 1 of complex*16 argument arg
c     in the upper half-plane Im(arg).ge.0;
c     before calling of the subroutine
c     array wsave must be initialized by the
c     subroutine haaus7
c 
c             INPUT PARAMETERS:
c 
c     wsave(*)- real*8 array; before calling of the subroutine
c               this array must be initialized by the
c               subroutine haaus7; the dimension of the array
c               must be at least 150;
c 
c     arg- complex*16 argument of the functions.
c 
c          OUTPUT PARAMETERS:
c 
c     ha0- (complex*16) value of the h-function of order 0;
c     ha1- (complex*16) value of the h-function of order 1;
c 
         implicit real*8 (a-h,o-z)
        save
         real*8 wsave(1)
         complex*16 arg,ha0,ha1
c 
c     compute indexes
c 
         igam=1
         ipsi=21
         ias0=41
         ias1=91
         x=cdabs(arg)
cc
         if (x.gt.3.d0) go to 10
c 
c     compute h-functions of small arguments by means
c     of Taylor's expansion
c 
         call  hanta7(15,wsave(igam),wsave(ipsi),arg,ha0,ha1)
         return
cc
 10      continue
         if (x.gt.17.d0) go to 20
c 
c     compute h-functions of intermediate arguments by means
c     of  evaluation of Poisson's integral
c 
         call  hahnk7(24,arg,ha0,ha1)
         return
cc
 20      continue
c 
c     compute h-functions of large arguments by means
c     of Hankel's asymptotic expansion
c 
         call hanas7(wsave(ias0),wsave(ias1),arg,ha0,ha1)
         return
         end
ccc
ccc
         subroutine haaus7(wsave)
c 
c     subroutine haaus7 initializes array wsave to be used
c     in subroutine hanvo7; it consists of arrays gam,
c     psi, aso and as1 (see subroutine haaux7)
c 
c             INPUT PARAMETERS: NONE
c 
c 
c             OUTPUT PARAMETER:
c 
c     wsave(*)- real*8 array; its dimension must be at least 150
c 
         implicit real*8 (a-h,o-z)
        save
         real*8 wsave(1)
c 
c     allocate the memory
c 
         igam=1
         ipsi=21
         ias0=41
         ias1=91
c 
c     initialize wsave
c 
         call haaux7(wsave(igam),wsave(ipsi),wsave(ias0),wsave(ias1))
         return
         end
ccc
ccc
         subroutine hahnk7(n,arg,ha0,ha1)
c 
c     subroutine hahnk7 evaluates hankel functions of
c     orders 0 and 1 of arbirary complex argument
c     arg under condition Im(arg).ge.0; the method
c     of evaluation is based on computation of the
c     Poisson's integrals for the Hankel functions.
c 
c          INPUT PARAMETERS:
c 
c     n- (integer) number of nodes of the quadrature
c     formula;
c 
c     pi- (real*8) constant pi=3.14525....;
c 
c     arg- (complex*16) argument of the functions.
c 
c          OUTPUT PARAMETERS:
c 
c     ha0- (complex*16) value of the h-function of order 0;
c     ha1- (complex*16) value of the h-function of order 1;
c 
      implicit complex *16 (a-h,o-z)
        save
      COMPLEX *16 ima
      real*8  h,pi,t,w,ex
      data ima /(0.0d0,1.0d0)/
      pi=3.14159265358979323846d0
c 
c     difine the step and auxiliary constants
c 
      h=6.5d0/n
      n1=n+1
      zr=ima/(2.d0*arg)
      con=2.d0*cdexp(ima*arg)/(pi*cdsqrt(arg))
c 
c     start the integration
c 
      ha0=0.d0
      ha1=0.d0
      do 15 i=1,n1
         t=(i-1)*h
         t=t**2
cc
         ex=dexp((-t))
         gen=cdsqrt(1.d0+t*zr)
         w=h
         if (i.eq.1) w=h/2.d0
cc
         ha0=ha0+w*ex/gen
         ha1=ha1+w*t*ex*gen
 15      continue
c 
c     multiply the results of integration by
c     proper constants
c 
         ha0=con*(1.d0-ima)*ha0
         ha1=(-2.d0)*con*(1.d0+ima)*ha1
         return
         end
ccc
ccc
         subroutine haaux7(gam,psi,as0,as1)
        save
  
c 
c     subroutine haaux7 evaluates array of factorials
c     and psi-functions of integer argument.le.16;
c     it also evaluates array of coefficients of
c     the hankel asymptotic expansion of the
c     h-functions of orders 0 and 1
c 
c 
c            INPUT PARAMETERS: NONE
c 
c            OUTPUT PARAMETERS:
c 
c     gam(*)- real*8 array of factorials; gam(i)=(i-1)!
c 
c     psi(*)- real*8 array of psi-functions; psi(i)=
c     psi-function(i);
c 
c     as0(*)- real*8 array of coefficients
c     of the asymptotic expansion of the h-function of order 0
c 
c     as1(*)- real*8 array of coefficients
c     of the asymptotic expansion of the h-function of order 1
c 
c     NOTE. Arrays aso and as1 do not contain the coefficients
c           of zero order that are equal to 1 for all orders
c 
c 
            implicit real*8 (a-h,o-z)
            real*8 gam(1),psi(1),as0(1),as1(1)
c 
c     define the starting values for evaluation of gam and psi
c 
            gam(1)=1.d0
            psi(1)=-0.57721566490153286d0
c 
c     recurce
c 
            do 100 i=1,16
               gam(i+1)=gam(i)*i
               psi(i+1)=psi(i)+1.d0/i
ccc   call prin2('gam(i)=*',gam(i),1)
 100        continue
c 
c     define the starting values for evaluation of as0 and as1
c 
            as0(1)=(-1.d0)/8.d0
            as1(1)=(-3.d0)*as0(1)
c 
c     recurce
c 
            do 200 i=1,41
               x=(2.d0*i+1.d0)**2
               y=8.d0*(i+1.d0)
               as0(i+1)=as0(i)*(-x)/y
               as1(i+1)=as1(i)*(4.d0-x)/y
 200        continue
            return
            end
ccc
ccc
         subroutine hanta7(m,gam,psi,arg,ha0,ha1)
c 
c     subroutine hanta7 evaluates h-functions of orders 0 and 1
c     by means of their Taylor's expansions
c 
c              INPUT PARAMETERS;
c 
c     m- integer number of terms in the expansion;
c 
c     gam(*),psi(*)- real*8 arrays (see subroutine haaux7);
c 
c     arg- (complex*16) argument of the functions.
c 
c          OUTPUT PARAMETERS:
c 
c     ha0- (complex*16) value of the h-function of order 0;
c 
c     ha1- (complex*16) value of the h-function of order 1;
c 
        save
  
  
         implicit real*8 (a-h,o-z)
         complex*16 arg,z,bj0,bj1,by0,by1,cdl,ima,ha0,ha1
 207     real*8 gam(1),psi(1)
         data ima /(0.0d0,1.0d0)/
ccc         pi=4.d0*datan(1.d0)
         pi=3.14159265358979323846d0
cc
         bj0=0.d0
         bj1=0.d0
         by0=0.d0
         by1=0.d0
         z=(-0.25d0)*arg**2
c 
c     start evaluation of the sums
c 
         do 100 i=1,m
c 
c     evaluate coefficients of the expansions
c 
            x0=1.d0/(gam(m+2-i)*gam(m+2-i))
            x1=1.d0/(gam(m+2-i)*gam(m+3-i))
            y0=psi(m+2-i)
            y1=psi(m+2-i)+psi(m+3-i)
c 
c     form the sums for computing j- and y-functions
c 
            bj0=z*(bj0+x0)
            bj1=z*(bj1+x1)
            by0=z*(by0+x0*y0)
            by1=z*(by1+x1*y1)
 100        continue
c 
c     correct the sums
c 
            bj0=1.d0+bj0
            bj1=0.5d0*arg*(1.d0+bj1)
cc
            cdl=2.d0*cdlog(0.5d0*arg)/pi
            by0=2.d0*(psi(1)+by0)/pi
            by1=arg*(psi(1)+psi(2)+by1)/(2*pi)
cc
cc
            by0=cdl*bj0-by0
            by1=(-2.d0)/(arg*pi)+cdl*bj1-by1
c 
c     evaluate h-functions
c 
            ha0=bj0+ima*by0
            ha1=bj1+ima*by1
            return
            end
ccc
ccc
            subroutine hanan7(m,as0,as1,arg,ha0,ha1)
c 
c     subroutine hanan7 evaluates h-functions of orders 0 and 1
c     by means of their asymptotic expansions with the given
c     number of terms
c 
c              INPUT PARAMETERS;
c 
c     m- integer number of terms in the expansion;
c 
c     as0(*),as1(*)- real*8 arrays (see subroutine haaux7);
c 
c     arg- (complex*16) argument of the functions.
c 
c          OUTPUT PARAMETERS:
c 
c     ha0- (complex*16) value of the h-function of order 0;
c 
c     ha1- (complex*16) value of the h-function of order 1;
c 
        save
  
  
            implicit real*8 (a-h,o-z)
            complex*16 arg,z,dom,ima,ha0,ha1
            real*8 as1(1),as0(1)
            data ima /(0.0d0,1.0d0)/
ccc            pi=4.d0*datan(1.d0)
            pi=3.14159265358979323846d0
  
cc
            ha0=0.d0
            ha1=0.d0
            z=(ima)/arg
c 
c     start evaluation of the sums
c 
         do 100 i=1,m
            ha0=z*(ha0+as0(m-i+1))
            ha1=z*(ha1+as1(m-i+1))
 100        continue
c 
c     correct the sums
c 
            ha0=1.d0+ha0
            ha1=1.d0+ha1
c 
c     evaluate the dominant term of the expansions
c 
            dom=(1.d0-ima)*cdexp(ima*arg)/cdsqrt(pi*arg)
c 
c     evaluate ha0 and ha1
c 
            ha0=dom*ha0
            ha1=(-ima)*dom*ha1
            return
            end
ccc
ccc
            subroutine hanas7(as0,as1,arg,ha0,ha1)
c 
c     subroutine hanas7 evaluates h-functions of orders 0 and 1
c     by means of their asymptotic expansions with the chioce
c     of the 'optimal' number of terms to be retained
c 
c 
c 
c              INPUT PARAMETERS;
c 
c 
c     as0(*),as1(*)- real*8 arrays (see subroutine haaux7);
c 
c     arg- (complex*16) argument of the functions.
c 
c          OUTPUT PARAMETERS:
c 
c     ha0- (complex*16) value of the h-function of order 0;
c 
c     ha1- (complex*16) value of the h-function of order 1;
c 
c     NOTE. The subroutine gives double precision accuracy
c           in the whole complex plane outside
c           the circle cdabs(arg).eq.17.d0
        save
  
            implicit real*8 (a-h,o-z)
            complex*16 arg,ha0,ha1
            real*8 as1(1),as0(1)
c 
c     evaluate the 'scaling' parameter
c 
            x=cdabs(arg)
c 
c     choose the 'optimal' number of terms m
c 
            m=25
            if (x.gt.25.d0) go to 10
            call hanan7(m,as0,as1,arg,ha0,ha1)
            return
cc
 10         continue
            m=15
            if (x.gt.50.d0) go to 20
            call hanan7(m,as0,as1,arg,ha0,ha1)
            return
cc
 20         continue
            m=10
            if (x.gt.300.d0) go to 30
            call hanan7(m,as0,as1,arg,ha0,ha1)
            return
cc
 30         continue
            m=5
            call hanan7(m,as0,as1,arg,ha0,ha1)
            return
            end
  
