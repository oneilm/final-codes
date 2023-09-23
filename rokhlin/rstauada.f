        implicit real *8 (a-h,o-z)
        dimension w(500 000),phiout(500 000),errs(200 000)
c 
c        test rstauexp for a linear system of equations
c 
cccc        call test1(w,phiout,errs)
c 
c        test rstauada for a linear system of equations
c 
cccc        call test2(w,phiout,errs)
c 
c        test rstauada for a non-linear system of equations
c 
cccc         call test3(w,phiout,errs)
  
c 
c        test rstauexp for a non-linear system of equations
c 
cccc        call test4(w,phiout,errs)
  
  
c 
c        test rstauexp for a non-linear system of equations
c        that is not an autonomous system
c 
cccc        call test5(w,phiout,errs)
  
  
c 
c        test rstauada for a non-linear system of equations
c        that is not an autonomous system
c 
cccc        call test6(w,phiout,errs)
c 
c        test rstauada for a Jacobi elliptic functions sn, cn, cn
c 
        call testjac(w,phiout,errs)
  
  
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine testjac(w,phiout,errs)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),errs(10 000),
     1      tout(10 000),aa(1000),bb(1000),errs(1)
        real *8 phiout(1),phi00(1000),phi11(1000)
        external funeva4,funder4
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
         PRINT *, 'ENTER n'
         READ *,n
         CALL PRINF('n=*',n,1 )
C 
         PRINT *, 'ENTER ncorrect'
         READ *,ncorrect
         CALL PRINF('ncorrect=*',ncorrect,1 )
c 
         PRINT *, 'ENTER rm'
         READ *,rm
         CALL PRIN2('rm*',rm,1)
  
  
         call prin2('and 1-rm=*',1-rm,1)
c 
         PRINT *, 'ENTER bbb - right end of the interval'
         READ *,bbb
         CALL PRIN2('bbb*',bbb,1)
c 
c        test the adaptive solver of the ode
c 
        m=3
        ifsave=1
        aaa=0
c 
        call ellfunev(rm,aaa,sn,cn,dn)
        phi00(1)=sn
        phi00(2)=cn
        phi00(3)=dn
 1200 continue
  
        call prin2('phi00 as created*',phi00,m)
  
        stepinit=1.0d-10
c 
        eps=1.0d-12
c 
        ifinit=1
c 
        maxinte=1000
c 
        do 3300 i=1,20 000
        w(i)=-77
 3300 continue
c 
        lw=100 000
cccc        tt1=second()
c 
  
        par1=rm
        call rstauada(ier,n,ncorrect,aaa,bbb,phi00,funeva4,
     1    funder4,m,par1,par2,par3,eps,stepinit,maxinte,ifsave,
     2    aa,bb,ninter,tout,phiout,phi11,nfcalls,nfinv,ifinit,w,
     3    lw,ltot)
c 
cccc        tt2=second()
cccc        call prin2('time for rstauada is*',tt2-tt1,1)
        call prinf('after rstauada, ier=*',ier,1)
        call prinf('after rstauada, nfcalls=*',nfcalls,1)
        call prinf('after rstauada, nfinv=*',nfinv,1)
        call prinf('after rstauada, ltot =*',ltot,1)
  
             if(ier .ne. 0) stop
  
        call prin2('after ccauada0, phi11=*',phi11,1)
        call prin2('after ccauada0, aa=*',aa,ninter)
        call prin2('after ccauada0, bb=*',bb,ninter)
c 
c       calculate and print the error
c 
        call errprjac(tout,phiout,ninter,n,m,errs,phi11,
     1      aaa,bbb,ifsave,aa,bb,rm)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine ellfunev(rm,u,sn,cn,dn)
        implicit real *8 (a-h,o-z)
        save
        dimension a(0:12),b(0:12),c(0:12),phis(0:12)
c 
c        this subroutine evaluates three (of the twelve) Jacobian
c        elliptic functions sn(u | rm), cn(u | rm), dn(u | rm)
c        via the standard Arithmetic-geometric mean process.
c 
c       conduct the arithmetic-geomatrical averaging process
c 
        a(0)=1
        b(0)=dsqrt(1-rm)
        c(0)=sqrt(rm)
c 
        do 1200 i=0,11
c 
        a(i+1)=(a(i)+b(i))/2
        b(i+1)=sqrt(a(i)*b(i))
        c(i+1)=(a(i)-b(i))/2
 1200 continue
c 
ccccc        call prin2('c as computed*',c,13)
c 
c        construct the auxiliary angles phi
c 
        done=1
        pi=datan(done)*4
c 
        phis(12)=2**12*a(12)*u
        do 1400 i=12,1,-1
c 
        rhs=c(i)/a(i)*sin(phis(i))
c 
        phis(i-1)=(asin(rhs)+phis(i))/2
 1400 continue
c 
c       constuct the elliptic functions
c 
        sn=sin(phis(0))
        cn=cos(phis(0))
        dn=cn/cos(phis(1)-phis(0))
        return
        end
c 
c 
c 
c 
c 
        subroutine funeva4(t,phi,f,m,rm,pareva2,pareva3)
        implicit real *8 (a-h,o-z)
        save
        real *8 phi(m),f(m)
c 
c       evaluate the rhs of the differential equation
c 
        f(1)=phi(2)*phi(3)
        f(2)=-phi(1)*phi(3)
        f(3)=-rm*phi(1)*phi(2)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine funder4(t,phi,dermat,m,a,pareva2,pareva3)
        implicit real *8 (a-h,o-z)
        save
        real *8 phi(m),dermat(m,m)
c 
c       evalate the differential of the rhs of the differential equation
c 
        done=1
        do 1400 i=1,m
        do 1200 j=1,m
        dermat(j,i)=0
 1200 continue
c 
 1400 continue
c 
        dermat(1,2)=phi(3)
        dermat(1,3)=phi(2)
c 
        dermat(2,1)=-phi(3)
        dermat(2,3)=-phi(1)
c 
        dermat(3,1)=-rm*phi(2)
        dermat(3,2)=-rm*phi(1)
c 
        return
        end
c 
c 
c 
c 
c 
  
        subroutine errprjac(tout,phiout,ninter,n,m,errs,
     1      phi11,a,b,ifsave,aa,bb,rm)
        implicit real *8 (a-h,o-z)
        save
        dimension tout(1),errs(m,1),errsb(100),aa(1),bb(1),
     1      rll(10 000),ccd(3),w1(20 000),w2(20 000),w3(20 000)
        real *8 phiout(m,1),cd,phi11(1)
c 
c       print the values of the solution
c 
        if(ifsave .eq. 0) goto 1500
cccc        call prin2('solution is*',phiout,n*ninter*m)
c 
c        calculate the solution exactly and calculate the error
c 
        errmax=0
        done=1
        do 1400 i=1,n*ninter
c 
        call ellfunev(rm,tout(i),ccd(1),ccd(2),ccd(3))
c 
        do 1200 j=1,m
c 
        cd=ccd(j)
        d=abs(phiout(j,i)-cd)
        errs(j,i)=d
        if(errmax .lt. d) errmax=d
 1200 continue
 1400 continue
c 
cccc        call prin2('and errs are*',errs,n*ninter*m)
c 
 1500 continue
c 
c        check the error at the right end
c 
        done=1
  
        call ellfunev(rm,b,ccd(1),ccd(2),ccd(3))
  
        do 1600 i=1,m
c 
        cd=ccd(i)
        errsb(i)=abs(phi11(i)-cd)
c 
cccc        if(errmax .lt. errsb(i)) errmax=errsb(i)
 1600 continue
  
        call prin2('and errors at right end are*',errsb,m)
  
        if(ifsave .eq. 0) goto 1700
c 
        call prin2('and errmax =*',errmax,1)
c 
c       plot the three components of the solution
c 
        do 1650 i=1,n*ninter
c 
        w1(i)=phiout(1,i)
        w2(i)=phiout(2,i)
        w3(i)=phiout(3,i)
 1650 continue
c 
        iw=21
        call lotagraph(iw,tout,w1,n*ninter,'function sn*')
c 
        iw=22
        call lotagraph(iw,tout,w2,n*ninter,'function cn*')
c 
        iw=23
        call lotagraph(iw,tout,w3,n*ninter,'function dn*')
c 
 1700 continue
        call prin2('and values at the right end =*',phi11,m)
c 
c       print the lenghts of the subintervals
c 
        do 1800 i=1,ninter
        rll(i)=bb(i)-aa(i)
 1800 continue
c 
        call prin2('and lengths of subintervals are*',rll,ninter)
        call prinf('while ninter=*',ninter,1)
        return
        end
c 
c 
c 
c 
c 
        subroutine test6(w,phiout,errs)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),errs(10 000),
     1      tout(10 000),aa(1000),bb(1000),errs(1)
        real *8 phiout(1),phi00(1000),phi11(1000)
        external funeva3,funder3
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
         PRINT *, 'ENTER n'
         READ *,n
         CALL PRINF('n=*',n,1 )
C 
         PRINT *, 'ENTER m'
         READ *,m
         CALL PRINF('m=*',m,1 )
c 
         PRINT *, 'ENTER ncorrect'
         READ *,ncorrect
         CALL PRINF('ncorrect=*',ncorrect,1 )
c 
c        test the adaptive solver of the ode
c 
        apar=-1
ccc        apar=-0.1d-6
        ifsave=0
        aaa=0
        bbb=10
c 
        ifsave=0
c 
        do 1200 i=1,m
        phi00(i)=(-1/apar)**i
 1200 continue
  
        call prin2('phi00 as created*',phi00,m)
  
        stepinit=0.0001
c 
        eps=1.0d-6
c 
        ifinit=1
c 
        maxinte=1000
c 
        do 3300 i=1,20 000
        w(i)=-77
 3300 continue
c 
        lw=100 000
cccc        tt1=second()
c 
  
        par1=apar
        call rstauada(ier,n,ncorrect,aaa,bbb,phi00,funeva3,
     1    funder3,m,par1,m,clams,eps,stepinit,maxinte,ifsave,
     2    aa,bb,ninter,tout,phiout,phi11,nfcalls,nfinv,ifinit,w,
     3    lw,ltot)
c 
cccc        tt2=second()
cccc        call prin2('time for rstauada is*',tt2-tt1,1)
        call prinf('after rstauada, ier=*',ier,1)
        call prinf('after rstauada, nfcalls=*',nfcalls,1)
        call prinf('after rstauada, nfinv=*',nfinv,1)
        call prinf('after rstauada, ltot =*',ltot,1)
  
             if(ier .ne. 0) stop
  
        call prin2('after ccauada0, phi11=*',phi11,1)
        call prin2('after ccauada0, aa=*',aa,ninter)
        call prin2('after ccauada0, bb=*',bb,ninter)
c 
c       calculate and print the error
c 
        call errprnon(tout,phiout,ninter,n,m,errs,phi11,
     1      aaa,bbb,ifsave,aa,bb,apar)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine test5(w,phiout,errs)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),errs(10 000),
     1      tout(10 000),aa(1000),bb(1000),errs(1)
        real *8 phiout(1),phi00(1000),phi11(1000)
        external funeva3,funder3
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
         PRINT *, 'ENTER n'
         READ *,n
         CALL PRINF('n=*',n,1 )
C 
         PRINT *, 'ENTER m'
         READ *,m
         CALL PRINF('m=*',m,1 )
c 
         PRINT *, 'ENTER ncorr'
         READ *,ncorr
         CALL PRINF('ncorr=*',ncorr,1 )
c 
c        test the adaptive solver of the ode
c 
        apar=-1
        apar=-0.1
c 
        ifsave=0
        aaa=0
        bbb=10
c 
        ifsave=0
c 
        do 1200 i=1,m
        phi00(i)=(-1/apar)**i
 1200 continue
  
        call prin2('phi00 as created*',phi00,m)
c 
        eps=1.0d-9
        stepinit=0.0001
c 
        ifinit=1
c 
        maxinte=1000
c 
        do 3300 i=1,20 000
        w(i)=-77
 3300 continue
c 
        lw=100 000
        ncorr2=ncorr
        n2=n-1
c 
        par1=apar
c 
cccc        tt1=second()
        call rrstauexp(ier,n,n2,ncorr,ncorr2,aaa,bbb,phi00,
     1    funeva3,funder3,m,par1,par2,par3,eps,stepinit,maxinte,ifsave,
     2    aa,bb,ninter,tout,phiout,phi11,nfcalls,nfinv,ifinit,w,
     3    lw,ltot)
c 
cccc        tt2=second()
cccc        call prin2('time for rrstauexp is*',tt2-tt1,1)
        call prinf('after rrstauexp, ier=*',ier,1)
        call prinf('after rrstauexp, nfcalls=*',nfcalls,1)
        call prinf('after rrstauexp, nfinv=*',nfinv,1)
        call prinf('after rrstauexp, ltot =*',ltot,1)
        call prinf('after rrstauexp, ninter =*',ninter,1)
  
  
             if(ier .ne. 0) stop
  
        call prin2('after ccauada0, phi11=*',phi11,1)
        call prin2('after ccauada0, aa=*',aa,ninter)
        call prin2('after ccauada0, bb=*',bb,ninter)
c 
c       calculate and print the error
c 
        call errprnon(tout,phiout,ninter,n,m,errs,phi11,
     1      aaa,bbb,ifsave,aa,bb,apar)
  
  
  
        return
        end
c 
c 
c 
c 
c 
        subroutine funeva3(t,phi,f,m,a,pareva2,pareva3)
        implicit real *8 (a-h,o-z)
        save
        real *8 phi(m),f(m)
  
  
ccccc        call prin2('in funeva3, a=*',a,1)
cccc        call prinf('in funeva2, m=*',m,1)
c 
c       evaluate the rhs of the differential equation
c 
        done=1
        do 1200 i=1,m
c 
        f(i)=-i*phi(i)**((i+1)*done/i)
  
cccc        f(i)=-i*phi(i)/(t-a)
 1200 continue
  
  
cccc        call prin2('exiting funeva2, f=*',f,m)
        return
        end
c 
c 
c 
c 
c 
        subroutine funder3(t,phi,dermat,m,a,pareva2,pareva3)
        implicit real *8 (a-h,o-z)
        save
        real *8 phi(m),dermat(m,m)
c 
c       evalate the differential of the rhs of the differential equation
c 
  
cccc        call prin2('in funder3, a=*',a,1)
cccc        call prin2('in funder2, t=*',t,1)
cccc        call prinf('in funder2, m=*',m,1)
        done=1
        do 1400 i=1,m
        do 1200 j=1,m
        dermat(j,i)=0
 1200 continue
c 
cccc        dermat(i,i)=-(i+1)*phi(i)**(done/i)
c 
        dermat(i,i)=-(i+1)/(t-a)
c 
 1400 continue
c 
  
cccc        call prin2('exiting funder2, dermat*',dermat,m*m)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine test4(w,phiout,errs)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),errs(10 000),
     1      tout(10 000),aa(1000),bb(1000),errs(1)
        real *8 phiout(1),phi00(1000),phi11(1000)
        external funeva2,funder2
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
         PRINT *, 'ENTER n'
         READ *,n
         CALL PRINF('n=*',n,1 )
C 
         PRINT *, 'ENTER m'
         READ *,m
         CALL PRINF('m=*',m,1 )
c 
         PRINT *, 'ENTER ncorr'
         READ *,ncorr
         CALL PRINF('ncorr=*',ncorr,1 )
c 
c        test the adaptive solver of the ode
c 
        apar=-1
        apar=-0.1
c 
        ifsave=0
        aaa=0
        bbb=10
c 
        ifsave=0
c 
        do 1200 i=1,m
        phi00(i)=(-1/apar)**i
 1200 continue
  
        call prin2('phi00 as created*',phi00,m)
c 
        eps=1.0d-9
        stepinit=0.0001
c 
        ifinit=1
c 
        maxinte=1000
c 
        do 3300 i=1,20 000
        w(i)=-77
 3300 continue
c 
        lw=100 000
        ncorr2=ncorr
        n2=n-1
c 
cccc        tt1=second()
        call rrstauexp(ier,n,n2,ncorr,ncorr2,aaa,bbb,phi00,
cccc     1    funeva2,funder2,m,par1,m,clams,eps,stepinit,maxinte,ifsave,
     1    funeva2,funder2,m,par1,par2,par3,eps,stepinit,maxinte,ifsave,
     2    aa,bb,ninter,tout,phiout,phi11,nfcalls,nfinv,ifinit,w,
     3    lw,ltot)
c 
cccc        tt2=second()
cccc        call prin2('time for rrstauexp is*',tt2-tt1,1)
        call prinf('after rrstauexp, ier=*',ier,1)
        call prinf('after rrstauexp, nfcalls=*',nfcalls,1)
        call prinf('after rrstauexp, nfinv=*',nfinv,1)
        call prinf('after rrstauexp, ltot =*',ltot,1)
        call prinf('after rrstauexp, ninter =*',ninter,1)
  
  
             if(ier .ne. 0) stop
  
        call prin2('after ccauada0, phi11=*',phi11,1)
        call prin2('after ccauada0, aa=*',aa,ninter)
        call prin2('after ccauada0, bb=*',bb,ninter)
c 
c       calculate and print the error
c 
        call errprnon(tout,phiout,ninter,n,m,errs,phi11,
     1      aaa,bbb,ifsave,aa,bb,apar)
  
  
  
        return
        end
  
  
c 
c 
c 
c 
c 
        subroutine test3(w,phiout,errs)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),errs(10 000),
     1      tout(10 000),aa(1000),bb(1000),errs(1)
        real *8 phiout(1),phi00(1000),phi11(1000)
        external funeva2,funder2
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
         PRINT *, 'ENTER n'
         READ *,n
         CALL PRINF('n=*',n,1 )
C 
         PRINT *, 'ENTER m'
         READ *,m
         CALL PRINF('m=*',m,1 )
c 
         PRINT *, 'ENTER ncorrect'
         READ *,ncorrect
         CALL PRINF('ncorrect=*',ncorrect,1 )
c 
c        test the adaptive solver of the ode
c 
        apar=-1
ccc        apar=-0.1d-6
        ifsave=0
        aaa=0
        bbb=10
c 
        ifsave=0
c 
        do 1200 i=1,m
        phi00(i)=(-1/apar)**i
 1200 continue
  
        call prin2('phi00 as created*',phi00,m*2)
  
        stepinit=0.0001
c 
        eps=1.0d-6
c 
        ifinit=1
c 
        maxinte=1000
c 
        do 3300 i=1,20 000
        w(i)=-77
 3300 continue
c 
        lw=100 000
cccc        tt1=second()
c 
  
        par1=apar
        call rstauada(ier,n,ncorrect,aaa,bbb,phi00,funeva2,
cccc     1    funder2,m,par1,m,clams,eps,stepinit,maxinte,ifsave,
     1    funder2,m,par1,par2,par3,eps,stepinit,maxinte,ifsave,
     2    aa,bb,ninter,tout,phiout,phi11,nfcalls,nfinv,ifinit,w,
     3    lw,ltot)
c 
cccc        tt2=second()
cccc        call prin2('time for rstauada is*',tt2-tt1,1)
        call prinf('after rstauada, ier=*',ier,1)
        call prinf('after rstauada, nfcalls=*',nfcalls,1)
        call prinf('after rstauada, nfinv=*',nfinv,1)
        call prinf('after rstauada, ltot =*',ltot,1)
  
             if(ier .ne. 0) stop
  
        call prin2('after ccauada0, phi11=*',phi11,1)
        call prin2('after ccauada0, aa=*',aa,ninter)
        call prin2('after ccauada0, bb=*',bb,ninter)
c 
c       calculate and print the error
c 
        call errprnon(tout,phiout,ninter,n,m,errs,phi11,
     1      aaa,bbb,ifsave,aa,bb,apar)
  
        return
        end
c 
c 
c 
c 
c 
  
        subroutine errprnon(tout,phiout,ninter,n,m,errs,
     1      phi11,a,b,ifsave,aa,bb,apar)
        implicit real *8 (a-h,o-z)
        save
        dimension tout(1),errs(m,1),errsb(100),aa(1),bb(1),
     1      rll(10 000)
        real *8 phiout(m,1),cd,clams(1),phi11(1)
c 
c       print the values of the solution
c 
        if(ifsave .eq. 0) goto 1500
cccc        call prin2('solution is*',phiout,n*ninter*m*2)
        call prin2('clams=*',clams,2*m)
c 
c        calculate the solution exactly and calculate the error
c 
        errmax=0
        done=1
        do 1400 i=1,n*ninter
c 
        do 1200 j=1,m
        cd=done/(tout(i)-apar)**j
        d=abs(phiout(j,i)-cd)
        errs(j,i)=d
        if(errmax .lt. d) errmax=d
 1200 continue
 1400 continue
c 
cccc        call prin2('and errs are*',errs,n*ninter*m)
c 
 1500 continue
c 
c        check the error at the right end
c 
        done=1
        do 1600 i=1,m
c 
        cd=exp(clams(i)*b)
  
        cd=done/(b-apar)**i
        errsb(i)=abs(phi11(i)-cd)
c 
cccc        if(errmax .lt. errsb(i)) errmax=errsb(i)
 1600 continue
  
        call prin2('and errors at right end are*',errsb,m)
        call prin2('and errmax =*',errmax,1)
        call prin2('and values at the right end =*',phi11,m)
c 
c       print the lenghts of the subintervals
c 
        do 1800 i=1,ninter
        rll(i)=bb(i)-aa(i)
 1800 continue
c 
        call prin2('and lengths of subintervals are*',rll,ninter)
        call prinf('while ninter=*',ninter,1)
        return
        end
c 
c 
c 
c 
c 
        subroutine funeva2(t,phi,f,m,a,pareva2,pareva3)
        implicit real *8 (a-h,o-z)
        save
        real *8 phi(m),f(m)
  
  
cccc        call prin2('in funeva2, phi=*',phi,m*2)
cccc        call prinf('in funeva2, m=*',m,1)
c 
c       evaluate the rhs of the differential equation
c 
        done=1
        do 1200 i=1,m
c 
        f(i)=-i*phi(i)**((i+1)*done/i)
 1200 continue
  
  
cccc        call prin2('exiting funeva2, f=*',f,m*2)
        return
        end
c 
c 
c 
c 
c 
        subroutine funder2(t,phi,dermat,m,a,pareva2,pareva3)
        implicit real *8 (a-h,o-z)
        save
        real *8 phi(m),a(m,m),dermat(m,m)
c 
c       evalate the differential of the rhs of the differential equation
c 
  
cccc        call prin2('in funder2, phi=*',phi,m*2)
cccc        call prin2('in funder2, t=*',t,1)
cccc        call prinf('in funder2, m=*',m,1)
        done=1
        do 1400 i=1,m
        do 1200 j=1,m
        dermat(j,i)=0
 1200 continue
c 
cccc        dermat(i,i)=-(i+1)*phi(i)**(done/i)
        dermat(i,i)=-(i+1)*phi(i)**(done/i)
c 
 1400 continue
c 
  
cccc        call prin2('exiting funder2, dermat*',dermat,m*2*m)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine test2(w,phiout,errs)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),errs(10 000),
     1      tout(10 000),aa(1000),bb(1000),errs(1)
        real *8 clams(1000),phiout(1),phi00(1000),phi11(1000)
        external funeva5,funder5
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
         PRINT *, 'ENTER n'
         READ *,n
         CALL PRINF('n=*',n,1 )
C 
         PRINT *, 'ENTER m'
         READ *,m
         CALL PRINF('m=*',m,1 )
c 
         PRINT *, 'ENTER ncorrect'
         READ *,ncorrect
         CALL PRINF('ncorrect=*',ncorrect,1 )
c 
         PRINT *, 'ENTER clam'
         READ *,clam
c 
c        test the adaptive solver of the ode
c 
        ifsave=0
        aaa=0
        bbb=10
c 
        ifsave=0
c 
        do 1200 i=1,m
        phi00(i)=1
        clams(i)=-clam*i
cccc        clams(i)=10*clams(i)/(1+i**2)
c 
cccc        clams(1)=-10000000
 1200 continue
c 
        stepinit=0.0001
c 
        eps=1.0d-6
c 
        call prin2('clams as created*',clams,m)
        ifinit=1
c 
        maxinte=1000
c 
        do 3300 i=1,20 000
        w(i)=-77
 3300 continue
c 
        lw=100 000
cccc        tt1=second()
c 
        call rstauada(ier,n,ncorrect,aaa,bbb,phi00,funeva5,
     1    funder5,m,par1,m,clams,eps,stepinit,maxinte,ifsave,
     2    aa,bb,ninter,tout,phiout,phi11,nfcalls,nfinv,ifinit,w,
     3    lw,ltot)
c 
cccc        tt2=second()
cccc        call prin2('time for rstauada is*',tt2-tt1,1)
        call prinf('after rstauada, ier=*',ier,1)
        call prinf('after rstauada, nfcalls=*',nfcalls,1)
        call prinf('after rstauada, nfinv=*',nfinv,1)
        call prinf('after rstauada, ltot =*',ltot,1)
  
             if(ier .ne. 0) stop
  
        call prin2('after ccauada0, phi11=*',phi11,1)
        call prin2('after ccauada0, aa=*',aa,ninter)
        call prin2('after ccauada0, bb=*',bb,ninter)
c 
c       calculate and print the error
c 
        call errprin(tout,phiout,clams,ninter,n,m,errs,phi11,
     1      aaa,bbb,ifsave,aa,bb)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine test1(w,phiout,errs)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),errs(10 000),
     1      tout(10 000),aa(1000),bb(1000),errs(1)
        real *8 clams(1000),phiout(1),phi00(1000),phi11(1000)
        external funeva5,funder5
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
         PRINT *, 'ENTER n'
         READ *,n
         CALL PRINF('n=*',n,1 )
C 
         PRINT *, 'ENTER m'
         READ *,m
         CALL PRINF('m=*',m,1 )
c 
         PRINT *, 'ENTER ncorr'
         READ *,ncorr
         CALL PRINF('ncorr=*',ncorr,1 )
c 
         PRINT *, 'ENTER clam'
         READ *,clam
         CALL PRIN2('clam=*',clam,2 )
c 
c        test the adaptive solver of the ode
c 
        aaa=0
        bbb=10
c 
        ifsave=0
c 
        do 1200 i=1,m
        phi00(i)=1
        clams(i)=-clam*i
 1200 continue
c 
        eps=1.0d-6
        stepinit=0.0001
c 
        ifinit=1
c 
        maxinte=1000
c 
        do 3300 i=1,20 000
        w(i)=-77
 3300 continue
c 
        lw=100 000
        ncorr2=ncorr
        n2=n-1
c 
cccc        tt1=second()
        call rrstauexp(ier,n,n2,ncorr,ncorr2,aaa,bbb,phi00,
     1    funeva5,funder5,m,par1,m,clams,eps,stepinit,maxinte,ifsave,
     2    aa,bb,ninter,tout,phiout,phi11,nfcalls,nfinv,ifinit,w,
     3    lw,ltot)
c 
cccc        tt2=second()
cccc        call prin2('time for rrstauexp is*',tt2-tt1,1)
        call prinf('after rrstauexp, ier=*',ier,1)
        call prinf('after rrstauexp, nfcalls=*',nfcalls,1)
        call prinf('after rrstauexp, nfinv=*',nfinv,1)
        call prinf('after rrstauexp, ltot =*',ltot,1)
        call prinf('after rrstauexp, ninter =*',ninter,1)
  
  
             if(ier .ne. 0) stop
  
        call prin2('after ccauada0, phi11=*',phi11,1)
        call prin2('after ccauada0, aa=*',aa,ninter)
        call prin2('after ccauada0, bb=*',bb,ninter)
c 
c       calculate and print the error
c 
        call errprin(tout,phiout,clams,ninter,n,m,errs,phi11,
     1      aaa,bbb,ifsave,aa,bb)
  
        return
        end
c 
c 
c 
c 
c 
  
        subroutine errprin(tout,phiout,clams,ninter,n,m,errs,
     1      phi11,a,b,ifsave,aa,bb)
        implicit real *8 (a-h,o-z)
        save
        dimension tout(1),errs(m,1),errsb(100),aa(1),bb(1),
     1      rll(10 000)
        real *8 phiout(m,1),cd,clams(1),phi11(1)
c 
c       print the values of the solution
c 
        if(ifsave .eq. 0) goto 1500
cccc        call prin2('solution is*',phiout,n*ninter*m*2)
        call prin2('clams=*',clams,2*m)
c 
c        calculate the solution exactly and calculate the error
c 
        errmax=0
        do 1400 i=1,n*ninter
c 
        do 1200 j=1,m
        cd=exp(clams(j)*tout(i))
        d=abs(phiout(j,i)-cd)
        errs(j,i)=d
        if(errmax .lt. d) errmax=d
 1200 continue
 1400 continue
c 
cccc        call prin2('and errs are*',errs,n*ninter*m)
c 
 1500 continue
c 
c        check the error at the right end
c 
        do 1600 i=1,m
c 
        cd=exp(clams(i)*b)
        errsb(i)=abs(phi11(i)-cd)
c 
cccc        if(errmax .lt. errsb(i)) errmax=errsb(i)
 1600 continue
  
        call prin2('and errors at right end are*',errsb,m)
        call prin2('and errmax =*',errmax,1)
        call prin2('and values at the right end =*',phi11,m)
c 
c       print the lenghts of the subintervals
c 
        do 1800 i=1,ninter
        rll(i)=bb(i)-aa(i)
 1800 continue
c 
        call prin2('and lengths of subintervals are*',rll,ninter)
        call prinf('while ninter=*',ninter,1)
        return
        end
c 
c 
c 
c 
c 
        subroutine funeva5(t,phi,f,m,pareva1,pareva2,clams)
        implicit real *8 (a-h,o-z)
        save
        real *8 phi(m),f(m),clams(1)
  
cccc        call prin2('in funeva5, clams=*',clams,m*2)
c 
        do 1200 i=1,m
c 
        f(i)=phi(i)*clams(i)
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine funder5(t,phi,dermat,m,parder1,pareva2,clams)
        implicit real *8 (a-h,o-z)
        save
        real *8 phi(m),dermat(m,m),clams(1)
c 
c       evalate the differential of the rhs of the differential equation
c 
        do 1400 i=1,m
        do 1200 j=1,m
        dermat(j,i)=0
 1200 continue
        dermat(i,i)=clams(i)
 1400 continue
  
cccc        call prin2('in funder5, clams=*',clams,m*2)
c 
        return
        end
c 
c 
c 
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        this is the end of the debugging code and the beginning of
c        the adaptive ode solver
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine rstauada(ier,n,ncorrect,a,b,phi0,funeva,
     1    funder,m,par1,par2,par3,eps,stepinit,maxinte,ifsave,
     2    aa,bb,ninter,tout,phiout,phi1,nfcalls,nfinv,ifinit,
     3      w,lenw,ltot)
        implicit real *8 (a-h,o-z)
        save
        dimension aa(1),bb(1),tout(1),w(1),par1(1),par2(1),
     1      par3(1)
        real *8 phi0(m),phi1(m),phiout(m,1)
c 
c       this subroutine solves adaptively the real differential
c       equation
c 
c       phi'(t) = f(t,phi(t)).                                   (1)
c 
c        the adaptive scheme is based on the spectral defferred
c        corrections, and is of arbitrary (user-specified) order.
c 
c          IMPORTANT NOTE:
c 
c          IN THIS IMPLEMENTATION, THE SUBROUTINE USES A DIRECT
C          LINEAR SOLVER AS A PART OF THE MARCHING SCHEME. THUS,
C          THIS IS A LOW-DIMENSIONAL CODE. DEEP INSIDE, IT HAS A
C          MEMORY ALLOCATION THAT LIMITS THE MAXIMUM DIMENSIONALITY
C          OF THE PROBLEM TO BE SOLVED BY 300 (AND THAT IS WAY TOO
C          HIGH; THE PROGRAM WILL TAKE FOREVER TO RUN). MANY PROBLEMS
C          ARE SPARSE; CONVERTING THIS CODE TO A SPARSE ONE WILL
C          INVOLVE CHANGING THE SUBROUTINES RSTSOL1, RSTSOL2.
c 
c                       input parameters:
c 
c  n - the number of nodes of the spectral discretization of each
c        subinterval
c  ncorrect - the number of defferred corrections.
c        IMPORTANT note. The order of the method is equal to
c        min(n,ncorrect+1); the actual optimum choice of the parameters
c        n, m appropriate to a particular problem and accuracy are
c        non-trivial, and have not been investigated in sufficient detail
c  a,b - the ends of the interval on which the equation is to be solved
c  phi0 - the initial value for the cauchy problem (specified at the
c        point a)
c  funeva - the user-supplied function that evaluates the right-hand
c        side in the differential equation. It has the calling sequence
c 
c            funeva(t,phi,ff,par1,par2,par3).
c 
c        The input parameters for funeva are
c 
c     t - the independent variable in (1),
c     phi - the value of the function in (1)
c     par1, par2, par3 - any parameters (real, integer, variables,
c        arrays, whatever) needed by the subroutine funeva.
c 
c        The output parameter for funeva is
c 
c     ff - the value of the derivative of phi
c 
c 
c*********************************************************************
c 
c  funder - the user-supplied function providing the solution
c     of the equation
c 
c     phi_{k+1}-h*f(phi_{k+1},t_{k+1})=phi_k                          (2)
c 
c      that has to be solved as a part of the backward Euler (which,
c      after all, is what drives this stupid scheme). The calling
c      sequence of the subroutine funder is
c 
c            funder(t(k+1),phi(k-1),phi(k),m,h,par1,par2,par3)
c 
c       The input parameters of funder are:
c 
c  t(k+1) - the values of the time at which the new value of phi
c       is to be determined
c  phi(k-1) - the vlaue of the solution at the point preceding the one where
c       the new value is to be determined
c  m - the dimensionality of the problem
c  h - the step length
c  par1, par2, par3  - whatever parameters the subroutine might need
c 
c       The output parameter of the subroutine funder
c 
c  phi(k) - the solution of the equation  (2)
c 
c********************************************************************
c 
c  m - the dimensionality of the problem
c  par1,par2,par3 - the parameters to be used by the subroutine funeva
c  eps - the precision to which the solution is to be evaluated
c  stepinit - the initial step length to be used by the marching
c        algorithm. Specifying a negative stepinit will cause the
c        subroutine to start with attempting to use the whole
c        interval [a,b]. This is generally not a good idea.
c        Starting with a very small stepinit is safe, though it
c        will increase the execution time somewhat.
c  maxinte - the maximum number of refined subintervals that the
c     subroutine is permitted to create
c  ifsave - the parameter telling the subroutine whether it is to
c         return to the user the parameters tout, phiout (see below).
c    ifsave=0 will cause the subroutine to not return to
c         the user the parameters tout, phiout; the only relevant output
c         in this case is the value phi1 of the solution at the point b.
c    ifsave=1 will cause the subroutine to return to the user the
c         points tout that have been generated in the process of the
c         adaptive solution, and the values  phiout of the solution
c         at these points.
c         in this case is teedless to say, this might means that the
c         subroutine will require a fairly large amount of storage
c         in arrays tout and (especially) phiout.
c  ifinit - the parameter telling the subroutine whether it should
c     construct the spectral integration matrix and the associated
c     parameters. ifinit=1 will cause the subroutine to construct
c     the parameters (it is somewhat expensive). ifinit=0 will cause
c     the subroutine to skip the creation of the parameters. In this
c     case, it the first n**2+10 elements of the array w must not have
c     been changed since the preceding call.
c  w - this is the input parameter only is ifinit has been set to 0.
c     it contains the spectral integration matrix on the n Gaussian
c     nodes.
c  lenw - the amount of storage (in real *8 locations) provided in the
c     array w
c 
c                         output parameters:
c 
c  ier - error return code.
c         ier=0 means successful execution
c         ier=4 means that the subroutine failed to construct the
c                solution to the user-specified precision within the
c                user-specified maximum number of steps maxinte
c         ier=8 means that the subroutinefailed to construct the
c                solution to the user-specified precision within
c                100 000 steps (the maximum permitted)
c         ier=256 means that at some point, the subroutine failed
c                to get the desired accuracy on {\bf one} interval,
c                even though the interval was smaller than
c                (bbb-aaa)/10**8. Obviously, something is very seriously
c                wrong
c         ier=512 means that at some point, the scheme kept exploding
c                numerically at even though the interval was
c                smaller than (bbb-aaa)/10**8. Obviously, something is
c                very seriously wrong with the subroutine funeva.
c         ier=1024 means that the amount of storage provided by the
c                user in the array w (specified by the parameter lenw)
c                is insufficient
c  aa - the left ends of solution subintervals
c  bb - the right ends of solution subintervals
c  ninter - the number of refuinements actually constructed
c  tout - the points at which the solution of (1) has been calculated
c  phiout - the values of the solution of (1) at the points tout
c  phi1 - the value of the solution of (1) at the point b
c  nfcalls - the number of times the user-supplied subroutine funeva
c        (see above) has been called
c  nfcalls - the number of times the user-supplied subroutine funder
c        (see above) has been called
c  w - this is an output parameter only if ifinit=1 has been specified
c 
c                         work arrays:
c 
c  w - must be at least 4*n**2+3*n+8*n*m+10*m + 4*m + 200.
c 
c       . . . allocate memory for the solution of the ode
c 
        iainte=1
        lainte=n*n+2
c 
        iw=iainte+lainte
        lw=3*n**2+3*n+50
c 
        nc=2
c 
        iphik1=iw+lw
        lphik1=m*n*nc+10
c 
        iff=iphik1+lphik1
        lff=m*n*nc+10
c 
        idelta=iff+lff
        ldelta=m*n*nc+10
c 
        irhs=idelta+ldelta
        lrhs=m*n*nc+10
c 
        icdnon=irhs+lrhs
        lcdnon=m*nc+10
c 
        iphi0=icdnon+lcdnon
        lphi0=m*nc+10
c 
        iuuu=iphi0+lphi0
        luuu=n**2+2
c 
        icw1=iuuu+luuu
        lcw1=n*nc+2
        lcw1=lcw1+m*nc
c 
        icw2=icw1+lcw1
        lcw2=n*nc+2
        lcw2=lcw2+m*nc
c 
        ltot=icw2+lcw2
        if(ltot .le. lenw) goto 1400
        ier=1024
        return
 1400 continue
c 
c        solve the ode
c 
        nfcalls=0
        nfinv=0
        maxiter=100 000
c 
        call rstauada1(ier,n,a,b,phi0,funeva,funder,
     1      par1,par2,par3,tout,phiout,ncorrect,
     2      aa,bb,maxiter,eps,stepinit,ninter,ifinit,phi1,
     3      w(iainte),w(iw),maxinte,m,ifsave,
     4      w(iphik1),w(iff),w(idelta),w(irhs),w(icdnon),
     5      w(iphi0),w(iuuu),w(icw1),w(icw2),nfcalls,nfinv)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine rstauada1(ier,n,aaa,bbb,phi00,funeva,funder,
     1    par1,par2,par3,tout,phiout,ncorrect,aa,bb,
     2      maxiter,eps,stepinit,ninter,ifinit,phi11,ainte,w,
     3      maxinte,m,ifsave,phik1,ff,delta,rhs,cdnon,
     4      phi0,uuu,cw1,cw2,nfcalls,nfinv)
        implicit real *8 (a-h,o-z)
        save
        external funeva
        dimension par1(1),par2(1),par3(1)
        dimension ainte(1),ts(100),whts(100),endinter(100),
     1    w(1),t(100),tout(1),aa(1),bb(1),uuu(1)
        real *8 phi00(m),phi0(1),phik1(1),phiout(m,1),
     1    phi11(m),ff(1),delta(1),rhs(1),cdnon(1),
     2      cw1(1),cw2(1)
c 
c       this subroutine solves adaptively the differential
c       equation
c 
c       phi'(t) = f(t,phi(t)).                                   (1)
c 
c        the adaptive scheme is based on the spectral defferred
c        corrections, and is of arbitrary (user-specified) order.
c 
c                       input parameters:
c 
c  n - the number of nodes of the spectral discretization of each
c        subinterval
c  aaa,bbb - the ends of the interval on which the equation is to be solved
c  phi00 - the initial value for the cauchy problem (specified at the
c        point a)
c  funeva - the user-supplied function that evaluates the right-hand
c        side in the differential equation. It has the calling sequence
c 
c            funeva(t,phi,ff,par1,par2,par3).
c 
c        The input parameters for funeva are
c 
c     t - the independent variable in (1),
c     phi - the value of the function in (1)
c     par1, par2, par3 - any parameters (real, integer, variables,
c        arrays, whatever) needed by the subroutine funeva.
c 
c        The output parameter for funeva is
c 
c     ff - the value of the derivative of phi
c 
c  par1,par2,par3 - the parameters to be used by the subroutine funeva
c  ncorrect - the number of defferred corrections.
c        IMPORTANT note. The order of the method is equal to
c        min(n,ncorrect+1); the actual optimum choice of the parameters
c        n, m appropriate to a particular problem and accuracy are
c        non-trivial, and have not been investigated in sufficient detail
c  maxiter - the maximum total number of intervals on which the
c        subroutine is permitted to {\bf attempt} to construct
c        defferred corrections. Should be LARGE
c  eps - the precision to which the solution is to be evaluated
c  stepinit - the initial step length to be used by the marching
c        algorithm. Specifying a negative stepinit will cause the
c        subroutine to start with attempting to use the whole
c        interval [a,b]. This is generally not a good idea.
c        Starting with a very small stepinit is safe, though it
c        will increase the execution time somewhat.
c  ifinit - the parameter telling the subroutine whether it should
c     construct the spectral integration matrix and the associated
c     parameters. ifinit=1 will cause the subroutine to construct
c     the the integration matrix (it is somewhat expensive).
c     ifinit=0 will cause the subroutine to skip the creation of the
c     integration matrix. In this it the parameters ainte, uuu
c     must not have been chenged since the preceding call.
c  ainte - the matrix of indefinite integration on the n Gaussian
c     nodes on the interval [-1,1]. This is an input parameter only
c     if the user has specified ifinit=0
c  maxinte - the maximum number of refined subintervals that the
c     subroutine is permitted to create
c  m - the dimensionality of the problem
c  ifsave - the parameter telling the subroutine whether it is to
c         return to the user the parameters tout, phiout (see below).
c    ifsave=0 will cause the subroutine to not return to
c         the user the parameters tout, phiout; the only relevant output
c         in this case is the value phi1 of the solution at the point b.
c    ifsave=1 will cause the subroutine to return to the user the
c         points tout that have been generated in the process of the
c         adaptive solution, and the values  phiout of the solution
c         at these points.
c         in this case is teedless to say, this might means that the
c         subroutine will require a fairly large amount of storage
c         in arrays tout and (especially) phiout.
c  uuu - the matrix converting the values of a function at n Gaussian
c         nodes into the first n Legendre coefficients of that function
c 
c                         output parameters:
c 
c  ier - error return code.
c         ier=0 means successful execution
c         ier=4 means that the subroutine failed to construct the
c                solution to the user-specified precision within the
c                user-specified maximum number of steps maxinte
c         ier=8 means that the subroutine failed to construct the
c                solution to the user-specified precision within
c                100 000 steps (the maximum permitted)
c         ier=256 means that at some point, the subroutine failed
c                to get the desired accuracy on {\bf one} interval,
c                even though the interval was smaller than
c                (bbb-aaa)/10**8. Obviously, something is very seriously
c                wrong
c         ier=512 means that at some point, the scheme kept exploding
c                numerically at some point even though the interval was
c                smaller than (bbb-aaa)/10**8. Obviously, something is
c                very seriously wrong with the subroutine funeva.
c         ier=1024 means that the amount of storage provided by the
c                user in the array w (specified by the parameter lenw)
c                is insufficient
c  aa - the left ends of solution subintervals
c  bb - the right ends of solution subintervals
c  ninter - the number of refuinements actually constructed
c  tout - the points at which the solution of (1) has been calculated
c  phiout - the values of the solution of (1) at the points tout
c  phi11 - the value of the solution of (1) at the point b
c  nfcalls - the number of times the user-supplied subroutine funeva
c        (see above) has been called
c  w - this is an output parameter only if ifinit=1 has been specified
c 
c                       work arrays:
c 
c  phik1,ff,delta,rhs - must be at least n*m+2 real *8 locations long
c  cd,cw1,cw2 - must be at least m+2 real *8 locations long
c 
c 
c        . . . if needed, construct the matrix of indefinite
c              integration at the Gaussian nodes
c 
        if(ifinit .eq. 0) goto 1200
        itype=1
        call rstauinm(n,ainte,adiff,ts,whts,endinter,
     1      itype,w,uuu)
c 
 1200 continue
c 
c       start the process of adaptive solution of the ode
c 
        aa(1)=aaa
        bb(1)=aaa
        step=bbb-aaa
  
        if(stepinit .gt. 0) step=stepinit
c 
        smallcoe=1.0d-12
        small=step*smallcoe
c 
        call rstarrmo(phi00,phi0,m*2)
c 
        ninter=1
        numpts=0
        ifout=0
c 
c       solve the ode adaptively
c 
        ndouble=2
        nmore=0
        do 3000 i=1,maxiter
c 
c        construct the solution on the new interval
c 
        a=bb(ninter)
        b=a+step
c 
        bbbb=bbb-step/10
        if(b .ge. bbbb) ifout=1
        if(b .ge. bbbb) b=bbb
        if(b .ge. bbbb) step=b-a
c 
        call rststd(jer,ts,n,phi0,funeva,funder,par1,par2,par3,
     1      a,b,t,phik1,ncorrect,ainte,endinter,phi11,m,
     3      ff,delta,rhs,cdnon,err,eps,cw1,cw2,nfcalls,nfinv)
c 
        if(jer .ne. 0) call prinf(
     1      'in ccauada0 after rststd, jer=*',jer,1)
c 
c        if the step is so big that the corrected euler has
c        exploded, reduce the step and try again
c 
        if(jer .eq. 0) goto 1400
        step=step/2
        if(step .gt. small) goto 3000
        ier=512
        return
 1400 continue
c 
c       if the obtained solution is sufficiently accurate -
c       act accordingly
c 
        call rstdisrep(uuu,n,m,phik1,ff,discr,val,cw1,cw2)
c 
        if( (err .gt. eps*val) .or.
     1      (discr .gt. eps*val) ) goto 2000
c 
        aa(ninter+1)=bb(ninter)
        bb(ninter+1)=aa(ninter+1)+step
        ninter=ninter+1
c 
        if(ninter .lt. maxinte) goto 1800
c 
        ier=4
        return
 1800 continue
c 
        call rstarrmo(phi11,phi0,m*2)
c 
        if(ifsave .eq. 0) goto 1900
        call rstarrmo(phik1,phiout(1,numpts+1),n*m*2)
        call rstarrmo(t,tout(numpts+1),n)
 1900 continue
c 
        numpts=numpts+n
        nmore=nmore+1
        if(nmore .eq. ndouble) step=step*2
        if(nmore .eq. ndouble) nmore=0
c 
        if(ifout .eq. 1) goto 4000
c 
        goto 3000
c 
 2000 continue
c 
c       the obtained solution is not sufficiently accurate.
c       keep subdividing
c 
        ifout=0
        step=step/2
c 
        if(step .gt. small) goto 3000
        ier=256
        return
 3000 continue
        ier=8
 4000 continue
c 
c       eliminate the blank interval at the beginning of
c       arrays aa, bb
c 
        do 4200 i=1,ninter-1
        aa(i)=aa(i+1)
        bb(i)=bb(i+1)
 4200 continue
        ninter =ninter-1
c 
        return
        end
c 
c 
c 
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        this is the end of the debugging code and the beginning of
c        the adaptive ode solver
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine rrstauexp(ier,n,n2,ncorr,ncorr2,a,b,phi0,
     1    funeva,funder,m,par1,par2,par3,eps,stepinit,maxinte,
     2    ifsave,aa,bb,ninter,tout,phiout,phi1,nfcalls,nfinv,
     3    ifinit,w,lenw,ltot)
        implicit real *8 (a-h,o-z)
        save
        dimension aa(1),bb(1),tout(1),w(1),par1(1),par2(1),
     1      par3(1)
        real *8 phi0(m),phi1(m),phiout(m,1)
c 
c       this subroutine solves adaptively the differential
c       equation (generally, stiff)
c 
c       phi'(t) = f(t,phi(t)).                                   (1)
c 
c        the adaptive scheme is based on the spectral defferred
c        corrections, and is of arbitrary (user-specified) order.
c        The subroutine drives the defferred corrections by the
c        backward euler techniques, with subsequent extrapolation
c        to enforce the appropriate stiff behavior at infinity.
c 
c                       input parameters:
c 
c  n - the number of nodes of the spectral discretization of each
c        subinterval in the first of the defferred correction
c        schemes to be extrapolated
c  n2 - the number of nodes of the spectral discretization of each
c        subinterval in the second of the defferred correction
c        schemes to be extrapolated
c  ncorr - the number of defferred corrections the first of the
c        defferred correction schemes to be extrapolated
c  ncorr2 - the number of defferred corrections the second of the
c        defferred correction schemes to be extrapolated
c 
c        IMPORTANT note. The order of the method is equal to
c        min(n,ncorr+1,n2,ncorr2); the actual optimum choice of the
c        parameters n, n2, ncorrs, ncorrs2  appropriate to a particular
c        problem and accuracy are non-trivial, and have not been
c        investigated in sufficient detail
c  a,b - the ends of the interval on which the equation is to be solved
c  phi0 - the initial value for the cauchy problem (specified at the
c        point a)
c  funeva - the user-supplied function that evaluates the right-hand
c        side in the differential equation. It has the calling sequence
c 
c            funeva(t,phi,ff,par1,par2,par3).
c 
c        The input parameters for funeva are
c 
c     t - the independent variable in (1),
c     phi - the value of the function in (1)
c     par1, par2, par3 - any parameters (real, integer, variables,
c        arrays, whatever) needed by the subroutine funeva.
c 
c        The output parameter for funeva is
c 
c     ff - the value of the derivative of phi
c 
c 
c*********************************************************************
c 
c  funder - the user-supplied function providing the solution
c     of the equation
c 
c     phi_{k+1}-h*f(phi_{k+1},t_{k+1})=phi_k                          (2)
c 
c      that has to be solved as a part of the backward Euler (which,
c      after all, is what drives this stupid scheme). The calling
c      sequence of the subroutine funder is
c 
c            funder(t(k+1),phi(k-1),phi(k),m,h,par1,par2,par3)
c 
c       The input parameters of funder are:
c 
c  t(k+1) - the values of the time at which the new value of phi
c       is to be determined
c  phi(k-1) - the vlaue of the solution at the point preceding the one where
c       the new value is to be determined
c  m - the dimensionality of the problem
c  h - the step length
c  par1, par2, par3  - whatever parameters the subroutine might need
c 
c       The output parameter of the subroutine funder
c 
c  phi(k) - the solution of the equation  (2)
c 
c********************************************************************
c 
c  m - the dimensionality of the problem
c  par1,par2,par3 - the parameters to be used by the subroutine funeva
c  eps - the precision to which the solution is to be evaluated
c  stepinit - the initial step length to be used by the marching
c        algorithm. Specifying a negative stepinit will cause the
c        subroutine to start with attempting to use the whole
c        interval [a,b]. This is generally not a good idea.
c        Starting with a very small stepinit is safe, though it
c        will increase the execution time somewhat.
c  maxinte - the maximum number of refined subintervals that the
c     subroutine is permitted to create
c  ifsave - the parameter telling the subroutine whether it is to
c         return to the user the parameters tout, phiout (see below).
c    ifsave=0 will cause the subroutine to not return to
c         the user the parameters tout, phiout; the only relevant output
c         in this case is the value phi1 of the solution at the point b.
c    ifsave=1 will cause the subroutine to return to the user the
c         points tout that have been generated in the process of the
c         adaptive solution, and the values  phiout of the solution
c         at these points.
c         in this case is teedless to say, this might means that the
c         subroutine will require a fairly large amount of storage
c         in arrays tout and (especially) phiout.
c  ifinit - the parameter telling the subroutine whether it should
c     construct the spectral integration matrix and the associated
c     parameters. ifinit=1 will cause the subroutine to construct
c     the parameters (it is somewhat expensive). ifinit=0 will cause
c     the subroutine to skip the creation of the parameters. In this
c     case, it the first n**2+10 elements of the array w must not have
c     been changed since the preceding call.
c  w - this is the input parameter only is ifinit has been set to 0.
c     it contains the spectral integration matrix on the n Gaussian
c     nodes.
c  lenw - the amount of storage (in real *8 locations) provided in the
c     array w
c 
c                         output parameters:
c 
c  ier - error return code.
c         ier=0 means successful execution
c         ier=4 means that the subroutine failed to construct the
c                solution to the user-specified precision within the
c                user-specified maximum number of steps maxinte
c         ier=8 means that the subroutinefailed to construct the
c                solution to the user-specified precision within
c                100 000 steps (the maximum permitted)
c         ier=256 means that at some point, the subroutine failed
c                to get the desired accuracy on {\bf one} interval,
c                even though the interval was smaller than
c                (bbb-aaa)/10**8. Obviously, something is very seriously
c                wrong
c         ier=512 means that at some point, the scheme kept exploding
c                numerically at some point even though the interval was
c                smaller than (bbb-aaa)/10**8. Obviously, something is
c                very seriously wrong with the subroutine funeva.
c         ier=1024 means that the amount of storage provided by the
c                user in the array w (specified by the parameter lenw)
c                is insufficient
c  aa - the left ends of solution subintervals
c  bb - the right ends of solution subintervals
c  ninter - the number of refuinements actually constructed
c  tout - the points at which the solution of (1) has been calculated
c  phiout - the values of the solution of (1) at the points tout
c  phi1 - the value of the solution of (1) at the point b
c  nfcalls - the number of times the user-supplied subroutine funeva
c        (see above) has been called
c  nfcalls - the number of times the user-supplied subroutine funder
c        (see above) has been called
c  w - this is an output parameter only if ifinit=1 has been specified
c 
c                         work arrays:
c 
c  w - must be at least 4*n**2+3*n+8*n*m+10*m + 4*m + 200.
c 
c       . . . allocate memory for the solution of the ode
c 
        k=n
c 
        if(k .lt. n2) k=n2
c 
        iainte=1
        lainte=k*k+2
c 
        iw=iainte+lainte
        lw=3*k**2+3*k+50
c 
        nc=2
c 
        iphik1=iw+lw
        lphik1=m*k*nc+10
c 
        iff=iphik1+lphik1
        lff=m*k*nc+10
c 
        idelta=iff+lff
        ldelta=m*k*nc+10
c 
        irhs=idelta+ldelta
        lrhs=m*k*nc+10
c 
        icdnon=irhs+lrhs
        lcdnon=m*nc+10
c 
        iphi0=icdnon+lcdnon
        lphi0=m*nc+10
c 
        iuuu=iphi0+lphi0
        luuu=k**2+2
c 
        icw1=iuuu+luuu
        lcw1=k*nc+2
c 
cccc        lcw1=lcw1+m*nc
c 
        icw2=icw1+lcw1
        lcw2=k*nc+2
c 
cccc        lcw2=lcw2+m*nc
c 
        iainte2=icw2+lcw2
        lainte2=k*k+2
c 
        iuuu2=iainte2+lainte2
        luuu2=k*k+2
c 
        its2=iuuu2+luuu2
        lts2=k+2
c 
        iphi112=its2+lts2
        lphi112=m*nc+2
c 
        ltot=iphi112+lphi112
        if(ltot .le. lenw) goto 1400
        ier=1024
        return
 1400 continue
c 
c        solve the ode
c 
        nfcalls=0
        nfinv=0
        maxiter=100 000
c 
        call stauex1(ier,n,a,b,phi0,funeva,funder,
     1      par1,par2,par3,tout,phiout,ncorr,
     2      aa,bb,maxiter,eps,stepinit,ninter,ifinit,phi1,
     3      w(iainte),w(iw),maxinte,m,ifsave,
     4      w(iphik1),w(iff),w(idelta),w(irhs),w(icdnon),
     5      w(iphi0),w(iuuu),w(icw1),w(icw2),nfcalls,nfinv,
     6      n2,ncorr2,w(iainte2),w(iuuu2),w(its2),w(iphi112) )
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine stauex1(ier,n,aaa,bbb,phi00,funeva,funder,
     1    par1,par2,par3,tout,phiout,ncorr,aa,bb,
     2      maxiter,eps,stepinit,ninter,ifinit,phi11,ainte,w,
     3      maxinte,m,ifsave,phik1,ff,delta,rhs,cdnon,
     4      phi0,uuu,cw1,cw2,nfcalls,nfinv,
     5      n2,ncorr2,ainte2,uuu2,ts2,phi112)
        implicit real *8 (a-h,o-z)
        save
        external funeva
        dimension par1(1),par2(1),par3(1)
        dimension ainte(1),ts(100),whts(100),endinter(100),
     1    w(1),t(100),tout(1),aa(1),bb(1),uuu(1)
        real *8 phi00(m),phi0(1),phik1(1),phiout(m,1),
     1    phi11(m),ff(1),delta(1),rhs(1),cdnon(1),
     2      cw1(1),cw2(1),phi112(1)
        dimension ainte2(1),endinte2(100),ts2(1),uuu2(1)
c 
c       this subroutine solves adaptively the differential
c       equation
c 
c       phi'(t) = f(t,phi(t)).                                   (1)
c 
c        the adaptive scheme is based on the spectral defferred
c        corrections, and is of arbitrary (user-specified) order.
c 
c                       input parameters:
c 
c  n - the number of nodes of the spectral discretization of each
c        subinterval
c  aaa,bbb - the ends of the interval on which the equation is to be solved
c  phi00 - the initial value for the cauchy problem (specified at the
c        point a)
c  funeva - the user-supplied function that evaluates the right-hand
c        side in the differential equation. It has the calling sequence
c 
c            funeva(t,phi,ff,par1,par2,par3).
c 
c        The input parameters for funeva are
c 
c     t - the independent variable in (1),
c     phi - the value of the function in (1)
c     par1, par2, par3 - any parameters (real, integer, variables,
c        arrays, whatever) needed by the subroutine funeva.
c 
c        The output parameter for funeva is
c 
c     ff - the value of the derivative of phi
c 
c  par1,par2,par3 - the parameters to be used by the subroutine funeva
c  ncorr - the number of defferred corrections.
c        IMPORTANT note. The order of the method is equal to
c        min(n,ncorr+1); the actual optimum choice of the parameters
c        n, m appropriate to a particular problem and accuracy are
c        non-trivial, and have not been investigated in sufficient detail
c  maxiter - the maximum total number of intervals on which the
c        subroutine is permitted to {\bf attempt} to construct
c        defferred corrections. Should be LARGE
c  eps - the precision to which the solution is to be evaluated
c  stepinit - the initial step length to be used by the marching
c        algorithm. Specifying a negative stepinit will cause the
c        subroutine to start with attempting to use the whole
c        interval [a,b]. This is generally not a good idea.
c        Starting with a very small stepinit is safe, though it
c        will increase the execution time somewhat.
c  ifinit - the parameter telling the subroutine whether it should
c     construct the spectral integration matrix and the associated
c     parameters. ifinit=1 will cause the subroutine to construct
c     the the integration matrix (it is somewhat expensive).
c     ifinit=0 will cause the subroutine to skip the creation of the
c     integration matrix. In this it the parameters ainte, uuu
c     must not have been chenged since the preceding call.
c  ainte - the matrix of indefinite integration on the n Gaussian
c     nodes on the interval [-1,1]. This is an input parameter only
c     if the user has specified ifinit=0
c  maxinte - the maximum number of refined subintervals that the
c     subroutine is permitted to create
c  m - the dimensionality of the problem
c  ifsave - the parameter telling the subroutine whether it is to
c         return to the user the parameters tout, phiout (see below).
c    ifsave=0 will cause the subroutine to not return to
c         the user the parameters tout, phiout; the only relevant output
c         in this case is the value phi1 of the solution at the point b.
c    ifsave=1 will cause the subroutine to return to the user the
c         points tout that have been generated in the process of the
c         adaptive solution, and the values  phiout of the solution
c         at these points.
c         in this case is teedless to say, this might means that the
c         subroutine will require a fairly large amount of storage
c         in arrays tout and (especially) phiout.
c  uuu - the matrix converting the values of a function at n Gaussian
c         nodes into the first n Legendre coefficients of that function
c 
c                         output parameters:
c 
c  ier - error return code.
c         ier=0 means successful execution
c         ier=4 means that the subroutine failed to construct the
c                solution to the user-specified precision within the
c                user-specified maximum number of steps maxinte
c         ier=8 means that the subroutine failed to construct the
c                solution to the user-specified precision within
c                100 000 steps (the maximum permitted)
c         ier=256 means that at some point, the subroutine failed
c                to get the desired accuracy on {\bf one} interval,
c                even though the interval was smaller than
c                (bbb-aaa)/10**8. Obviously, something is very seriously
c                wrong
c         ier=512 means that at some point, the scheme kept exploding
c                numerically at some point even though the interval was
c                smaller than (bbb-aaa)/10**8. Obviously, something is
c                very seriously wrong with the subroutine funeva.
c         ier=1024 means that the amount of storage provided by the
c                user in the array w (specified by the parameter lenw)
c                is insufficient
c  aa - the left ends of solution subintervals
c  bb - the right ends of solution subintervals
c  ninter - the number of refuinements actually constructed
c  tout - the points at which the solution of (1) has been calculated
c  phiout - the values of the solution of (1) at the points tout
c  phi11 - the value of the solution of (1) at the point b
c  nfcalls - the number of times the user-supplied subroutine funeva
c        (see above) has been called
c  w - this is an output parameter only if ifinit=1 has been specified
c 
c                       work arrays:
c 
c  phik1,ff,delta,rhs - must be at least n*m+2 real *8 locations long
c  cd,cw1,cw2 - must be at least m+2 real *8 locations long
c 
c 
c        . . . if needed, construct the matrix of indefinite
c              integration at the Gaussian nodes
c 
        if(ifinit .eq. 0) goto 1200
        itype=1
        call rstauinm(n,ainte,adiff,ts,whts,endinter,
     1      itype,w,uuu)
c 
        call rstauinm(n2,ainte2,adiff,ts2,whts,endinte2,
     1      itype,w,uuu2)
c 
 1200 continue
c 
        if(2 .ne. 3) goto 1300
c 
        call prin2('ainte=*',ainte,n**2)
        call prin2('ainte2=*',ainte2,n2**2)
c 
        call prin2('ts=*',ts,n)
        call prin2('ts2=*',ts2,n2)
c 
        call prin2('endinter=*',endinter,n)
        call prin2('endinte2=*',endinte2,n2)
c 
        call prin2('uuu=*',uuu,n**2)
        call prin2('uuu2=*',uuu2,n2**2)
c 
        if(2 .ne. 3) stop
 1300 continue
c 
c       start the process of adaptive solution of the ode
c 
        aa(1)=aaa
        bb(1)=aaa
        step=bbb-aaa
  
        if(stepinit .gt. 0) step=stepinit
c 
        smallcoe=1.0d-12
        small=step*smallcoe
c 
        call rstarrmo(phi00,phi0,m*2)
c 
        ninter=1
        numpts=0
        ifout=0
c 
c       solve the ode adaptively
c 
        ndouble=2
        nmore=0
c 
        call stgetlim(n,ncorr,rlim)
        call stgetlim(n2,ncorr2,rlim2)
c 
        call prin2('after stgetlim, rlim=*',rlim,1)
        call prin2('after stgetlim, rlim2=*',rlim2,1)
c 
        do 3000 i=1,maxiter
c 
c        construct the solution on the new interval
c 
        a=bb(ninter)
        b=a+step
c 
        bbbb=bbb-step/10
        if(b .ge. bbbb) ifout=1
        if(b .ge. bbbb) b=bbb
        if(b .ge. bbbb) step=b-a
c 
        call rststd(jer,ts2,n2,phi0,funeva,funder,par1,par2,par3,
     1      a,b,t,phik1,ncorr2,ainte2,endinte2,phi112,m,
     2      ff,delta,rhs,cdnon,err,eps,cw1,cw2,nfcalls,nfinv)
c 
        call rststd(jer,ts,n,phi0,funeva,funder,par1,par2,par3,
     1      a,b,t,phik1,ncorr,ainte,endinter,phi11,m,
     2      ff,delta,rhs,cdnon,err,eps,cw1,cw2,nfcalls,nfinv)
c 
        if(jer .ne. 0) call prinf(
     1      'in ccauada0 after rststd, jer=*',jer,1)
c 
c        if the step is so big that the corrected euler has
c        exploded, reduce the step and try again
c 
        if(jer .eq. 0) goto 1400
        step=step/2
        if(step .gt. small) goto 3000
        ier=512
        return
 1400 continue
c 
c       if the obtained solution is sufficiently accurate -
c       act accordingly
c 
        call rstdisrep(uuu,n,m,phik1,ff,discr,val,cw1,cw2)
c 
        if( (err .gt. eps*val) .or.
     1      (discr .gt. eps*val) ) goto 2000
c 
c        . . . extrapolate the two solutions at the right end
c              to obtain appropriate stiff behavior at infinity
c 
cccc        call prin2('phi11 before extrapolation*',phi11,m*2)
        do 1600 j=1,m
c 
        phi11(j)=(rlim*phi112(j)-rlim2*phi11(j))/(rlim-rlim2)
 1600 continue
c 
cccc        call prin2('phi11 after extrapolation*',phi11,m*2)
  
        aa(ninter+1)=bb(ninter)
        bb(ninter+1)=aa(ninter+1)+step
        ninter=ninter+1
c 
        if(ninter .lt. maxinte) goto 1800
c 
        ier=4
        return
 1800 continue
c 
        call rstarrmo(phi11,phi0,m*2)
c 
        if(ifsave .eq. 0) goto 1900
        call rstarrmo(phik1,phiout(1,numpts+1),n*m*2)
        call rstarrmo(t,tout(numpts+1),n)
 1900 continue
c 
        numpts=numpts+n
        nmore=nmore+1
        if(nmore .eq. ndouble) step=step*2
        if(nmore .eq. ndouble) nmore=0
c 
        if(ifout .eq. 1) goto 4000
c 
        goto 3000
c 
 2000 continue
c 
c       the obtained solution is not sufficiently accurate.
c       keep subdividing
c 
        ifout=0
        step=step/2
c 
        if(step .gt. small) goto 3000
        ier=256
        return
 3000 continue
        ier=8
 4000 continue
c 
c       eliminate the blank interval at the beginning of
c       arrays aa, bb
c 
        do 4200 i=1,ninter-1
        aa(i)=aa(i+1)
        bb(i)=bb(i+1)
 4200 continue
        ninter =ninter-1
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine stgetlim(n,m,rlim)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension a1(2),a2(3),a3(4),a4(5),a5(6),a6(7),a7(8),
     1      a8(9),a9(10),a10(11),a11(12),a12(13),a13(14),a14(15),
     2      a15(16),a16(17),a17(18),a18(19),a19(20),a20(21)
c 
        data a1/
     1  0.00000000000000000E+00, 0.40000000000000000E-39/
c 
        data a2/
     1  -.61602540378443865E+00, -.19527222831138382E+00,
     2  -.61898816047911231E-01/
c 
        data a3/
     1  0.35736828507016316E+00, 0.46125571845918720E+00,
     2  0.25229689308007762E+00, 0.90890532918585648E-01/
c 
        data a4/
     1  -.23718356065183069E+00, -.35431317434978541E+00,
     2  -.39134896500399613E+00, -.25309980095027557E+00,
     3  -.95162652465009118E-01/
c 
        data a5/
     1  0.17035920981006248E+00, 0.27663436323939936E+00,
     2  0.35353652325674961E+00, 0.34897703876750667E+00,
     3  0.22838605460876855E+00, 0.76595971206798955E-01/
c 
        data a6/
     1  -.12900199978175976E+00, -.22217852433120531E+00,
     2  -.30792436215679936E+00, -.34856766441407516E+00,
     3  -.31008214584564913E+00, -.18671031258937550E+00,
     4  -.39263953593127135E-01/
c 
        data a7/
     1  0.10146978528417938E+00, 0.18288360281406904E+00,
     2  0.26803206728081261E+00, 0.32846478127954243E+00,
     3  0.33304789900508925E+00, 0.26321820732624092E+00,
     4  0.13044641605863905E+00, -.11870633832230317E-01/
c 
        data a8/
     1  -.82139366378411698E-01, -.15360699344472441E+00,
     2  -.23489539825188192E+00, -.30413508549826241E+00,
     3  -.33393548006303271E+00, -.30171478119217472E+00,
     4  -.20300930571777586E+00, -.61748192556765660E-01,
     5  0.70666133309460795E-01/
c 
        data a9/
     1  0.68004196351002640E-01, 0.13117188256314248E+00,
     2  0.20755864640051865E+00, 0.28019870702050765E+00,
     3  0.32536010209391840E+00, 0.31976525999808093E+00,
     4  0.25114809680987169E+00, 0.12830396151893976E+00,
     5  -.15206242928038287E-01, -.12934988019133071E+00/
c 
        data a10/
     1  -.57330426110915872E-01, -.11356431779013012E+00,
     2  -.18488988731654148E+00, -.25809101866505999E+00,
     3  -.31273415168196273E+00, -.32624875263415334E+00,
     4  -.28239153174236790E+00, -.18058384105566541E+00,
     5  -.41831261779652052E-01, 0.93484351284250088E-01,
     6  0.17860624342474751E+00/
c 
        data a11/
     1  0.49057586458750191E-01, 0.99463756264069204E-01,
     2  0.16592757306190927E+00, 0.23815446037872477E+00,
     3  0.29859681252067179E+00, 0.32599186101612335E+00,
     4  0.30216723477044950E+00, 0.22048472412015109E+00,
     5  0.92683359246769037E-01, -.49875014707220701E-01,
     6  -.16357558123445428E+00, -.20838962917407044E+00/
c 
        data a12/
     1  -.42505743713499398E-01, -.87975946389020301E-01,
     2  -.14991454427342601E+00, -.22033734477882048E+00,
     3  -.28419278498722200E+00, -.32175582334151145E+00,
     4  -.31409474051776766E+00, -.25059283342204641E+00,
     5  -.13605832884762899E+00, 0.59797271388600461E-02,
     6  0.13687019583327276E+00, 0.21435575923404493E+00,
     7  0.20943383667098227E+00/
c 
        data a13/
     1  0.37221737783838624E-01, 0.78477562423134641E-01,
     2  0.13626629398255348E+00, 0.20445334720848325E+00,
     3  0.27013805068541957E+00, 0.31518278981823288E+00,
     4  0.32058719038598061E+00, 0.27312695911734523E+00,
     5  0.17238093463406678E+00, 0.35255572015599991E-01,
     6  -.10498035543509334E+00, -.20682634857593239E+00,
     7  -.23479129697856911E+00, -.17527570428821811E+00/
c 
        data a14/
     1  -.32893603285603789E-01, -.70522923986373650E-01,
     2  -.12453216875129929E+00, -.19028285176486869E+00,
     3  -.25672784842996659E+00, -.30727795789147056E+00,
     4  -.32327704304845652E+00, -.28982478727574377E+00,
     5  -.20253259058698435E+00, -.72764362306582651E-01,
     6  0.71452680745644072E-01, 0.19087806978899306E+00,
     7  0.24682094043959876E+00, 0.21623317873344701E+00,
     8  0.10445094780952749E+00/
c 
        data a15/
     1  0.29300582941355596E-01, 0.63786161156072696E-01,
     2  0.11436266763745713E+00, 0.17761247400553120E+00,
     3  0.24408874763743976E+00, 0.29866939282277925E+00,
     4  0.32328752101117216E+00, 0.30201350000051057E+00,
     5  0.22744128428837902E+00, 0.10631398703762268E+00,
     6  -.38217598507047410E-01, -.16990950968341482E+00,
     7  -.24894819216003823E+00, -.24581405974107571E+00,
     8  -.15497638078893648E+00, -.23690652214780229E-02/
c 
        data a16/
     1  -.26282716833688835E-01, -.58024264517034762E-01,
     2  -.10548437388444487E+00, -.16624921751976849E+00,
     3  -.23225741163902730E+00, -.28975466061634483E+00,
     4  -.32140397452221413E+00, -.31070206191394105E+00,
     5  -.24794956445035794E+00, -.13604199687182948E+00,
     6  0.63036860644389555E-02, 0.14615538573283915E+00,
     7  0.24390923850065348E+00, 0.26575112798948589E+00,
     8  0.19751757559637274E+00, 0.54577196026173109E-01,
     9  -.11772163845421429E+00/
c 
        data a17/
     1  0.23721685057861682E-01, 0.53052826508608382E-01,
     2  0.97681320922797169E-01, 0.15602404315986795E+00,
     3  0.22122283999540941E+00, 0.28078680355461726E+00,
     4  0.31818339077482711E+00, 0.31665983141685879E+00,
     5  0.26478127279118554E+00, 0.16224014049543557E+00,
     6  0.23772317048577394E-01, -.12107936486664891E+00,
     7  -.23377681348117470E+00, -.27781066284897200E+00,
     8  -.23225317856565383E+00, -.10313770618980209E+00,
     9  0.72710346146583741E-01, 0.23525474474385019E+00/
c 
        data a18/
     1  -.21528369682489430E-01, -.48729719547395194E-01,
     2  -.90781250564012516E-01, -.14679112020899263E+00,
     3  -.21094958135436023E+00, -.27192638069656244E+00,
     4  -.31402565110494486E+00, -.32047699686721209E+00,
     5  -.27854361026634572E+00, -.18525134634252251E+00,
     6  -.51787304099742500E-01, 0.95637118698669397E-01,
     7  0.22009280592115136E+00, 0.28355530612639424E+00,
     8  0.25986466306879451E+00, 0.14696791857943756E+00,
     9  -.26695925953365834E-01, -.20603580163158398E+00,
     *  -.32476153808992430E+00/
c 
        data a19/
     1  0.19634563151303611E-01, 0.44943869414826077E-01,
     2  0.84645460759277476E-01, 0.13842548488563046E+00,
     3  0.20139067503174558E+00, 0.26327382022283019E+00,
     4  0.30922078651037173E+00, 0.32260951131110251E+00,
     5  0.28974070075858653E+00, 0.20542073304584601E+00,
     6  0.77686766873442531E-01, -.70447949204812173E-01,
     7  -.20399710331456465E+00, -.28429450094465926E+00,
     8  -.28119411490450885E+00, -.18572565707286986E+00,
     9  -.18386408128621351E-01, 0.17155364837110257E+00,
     *  0.31810120287605268E+00, 0.36025341853105954E+00/
c 
        data a20/
     1  -.17987310459111454E-01, -.41607378465347184E-01,
     2  -.79161258353477397E-01, -.13082029179769414E+00,
     3  -.19249489426078519E+00, -.25488998963556215E+00,
     4  -.30398089231760648E+00, -.32341218296571286E+00,
     5  -.29878953491859310E+00, -.22307325306901922E+00,
     6  -.10151173355899419E+00, 0.45906806496961750E-01,
     7  0.18633082183831482E+00, 0.28109896213009062E+00,
     8  0.29708964031172990E+00, 0.21946869084044407E+00,
     9  0.61381078267458756E-01, -.13412255398794546E+00,
     *  -.30299768174110233E+00, -.37952961593838959E+00,
     *  -.32115670502377906E+00/
c 
c        return to the user the requested coefficient
c 
        if(n .ne. 1) goto 1100
        rlim=a1(m)
 1100 continue
c 
        if(n .ne. 2) goto 1200
        rlim=a2(m)
 1200 continue
c 
        if(n .ne. 3) goto 1300
        rlim=a3(m)
 1300 continue
c 
        if(n .ne. 4) goto 1400
        rlim=a4(m)
 1400 continue
c 
        if(n .ne. 5) goto 1500
        rlim=a5(m)
 1500 continue
c 
        if(n .ne. 6) goto 1600
        rlim=a6(m)
 1600 continue
c 
        if(n .ne. 7) goto 1700
        rlim=a7(m)
 1700 continue
c 
        if(n .ne. 8) goto 1800
        rlim=a8(m)
 1800 continue
c 
        if(n .ne. 9) goto 1900
        rlim=a9(m)
 1900 continue
c 
        if(n .ne. 10) goto 2000
        rlim=a10(m)
 2000 continue
c 
        if(n .ne. 11) goto 2100
        rlim=a11(m)
 2100 continue
c 
        if(n .ne. 12) goto 2200
        rlim=a12(m)
 2200 continue
c 
        if(n .ne. 13) goto 2300
        rlim=a13(m)
 2300 continue
c 
        if(n .ne. 14) goto 2400
        rlim=a14(m)
 2400 continue
c 
        if(n .ne. 15) goto 2500
        rlim=a15(m)
 2500 continue
c 
        if(n .ne. 16) goto 2600
        rlim=a16(m)
 2600 continue
c 
        if(n .ne. 17) goto 2700
        rlim=a17(m)
 2700 continue
c 
        if(n .ne. 18) goto 2800
        rlim=a18(m)
 2800 continue
c 
        if(n .ne. 19) goto 2900
        rlim=a19(m)
 2900 continue
c 
        if(n .ne. 20) goto 3000
        rlim=a20(m)
 3000 continue
  
  
        return
        end
c 
c 
c 
c 
c 
        subroutine rstdisrep(uuu,n,m,phik,w,discr,val,cw1,cw2)
        implicit real *8 (a-h,o-z)
        save
        dimension uuu(1)
        real *8 w(m,n),phik(m,n),cw1(1),cw2(1)
c 
c       construct the Legendre transformation of the solution
c       tabulated at the Gaussian nodes
c 
        val=0
        do 2000 i=1,m
        do 1600 j=1,n
        val=val+phik(i,j)* (phik(i,j))
        cw1(j)=phik(i,j)
 1600 continue
c 
        call stmatve(uuu,cw1,cw2,n)
c 
        do  1800 j=1,n
        w(i,j)=cw2(j)
 1800 continue
cccc         call prin2('in rstdisrep, cw2=*',cw2,n*2)
  
 2000 continue
c 
c       calculate the sum of squares of the last ncheck
c       Legendre coefficients
c 
        discr=0
        ncheck=2
c 
        do 2400 i=1,m
        do 2200 j=1,ncheck
        discr=discr+w(i,n-j+1)* (w(i,n-j+1))
 2200 continue
 2400 continue
c 
        discr=dsqrt(discr)
        val=dsqrt(val)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine rststd(ier,ts,n,phi0,funeva,funder,
     1      par1,par2,par3,a,b,t,
     2      phik,ncorr,ainte,endinter,phi1,m,
     3      ff,delta,rhs,cd,err,eps,cw1,cw2,nfcalls,nfinv)
        implicit real *8 (a-h,o-z)
        save
        dimension par1(1),par2(1),par3(1)
        dimension t(1),ainte(1),ts(1),endinter(1)
        real *8 phi0(m),phik(m,n),cd(1),cdd,
     1      phi1(m),ff(1),delta(1),rhs(1),cw1(1),cw2(1)
  
        real *8 dermat(100 000),work(300),rhs2(300)
  
c 
c        this subroutine uses the implicit Euler-driven defferred
c        correction technique to solve a (generally) stiff system
c        of ODEs
c 
c        phi'(t)=fun(t,phi)                                          (1)
c 
c        on the interval discretized by the nodes t(i), i=1,2,...,n.
c 
c                       input parameters:
c 
c  ts - the n Gaussian nodes on the interval [-1,1]
c  n - the number of Gaussian nodes into which the interval [a,b]
c        is to be discretized
c  phi0 - the initial value at the point a for the cauchy problem (1)
c        to be solved
c  funeva - the user-supplied subroutine evaluating the
c        right-hand side in the equation (1) to be solved
c  par1,par2,par3 - the parameters to be used by the subroutine funeva
c  a,b - the ends of the interval on which
c  ncorr - the number of defferred corrections to be performed
c  ainte - the matrix of indefinite ingtegration at n Gaussian
c        nodes on the interval [-1,1]
c  endinter - the interpolation coefficients converting the
c          values of a function at n Gaussian nodes into its
c          value at 1 (the right end of the interval)
c  m - the dimensionality of the problem
c  eps - the precision to which the solution is to be evaluated
c  nfcalls - the counter of the number of times the user-supplied
c        subroutine funeva has been called. It will be incremented by
c        this subroutine.
c 
c                       output parameters:
c 
c  ier -  error return code;
c             ier=0 means successful conclusion
c             ier=4 means that the process exploded numerically
c                   during the correction process; probably, the
c                   stability condition has been exceeded
c             ier=16 means that the process exploded numerically
c                   during the initial Euler stage; probably, the
c                   stability condition has been dramatically
c                   exceeded
c  t - the Gaussian discretization of the interval [a,b]
c  phik - the solution of the equation (1) with the initial
c        value phi0 at the Gaussian nodes on the interval [a,b]
c  phi1 - the solution of the equation (1) at the point b
c  err - the precision to which the defferred corrections have
c        converged after ncorr iterations. Note that if err
c        became less than the user-specified eps after fewer than
c        ncorr iterations, the iterations are terminated
c  nfcalls - the counter of the number of times the user-supplied
c        subroutine funeva has been called. It is incremented by
c        this subroutine.
c 
c                       work arrays:
c 
c  ff,delta,rhs - must be at least n*m+2 real *8 locations long
c  cd,cw1,cw2 - must be at least m+2 real *8 locations long
c 
c        . . . construct the Gaussian discretization of the
c              interval [a,b]
c 
        ier=0
        eps22=eps**2
        dlarge=1.0d35
c 
        u=(b-a)/2
        v=(b+a)/2
        do 1100 i=1,n
        t(i)=u*ts(i)+v
 1100 continue
c 
c       . . . use euler to construct the initial approximation on
c             the short interval
c 
        call rstarrmo(phi0,phik(1,1),m*2)
c 
        do 1200 i=2,n
c 
        nfinv=nfinv+1
        nrefine=0
c 
        call rstsol1(funeva,funder,m,t(i)-t(i-1),phik(1,i-1),
     1      t(i),phik(1,i),
     1      dermat,work,rhs2,par1,par2,par3,
     2      par1,par2,par3,nrefine)
c 
        dd=0
        do 1150 j=1,m
c 
        dd=dd+phik(j,i)* (phik(j,i))
 1150 continue
c 
        if(dd .lt. dlarge) goto 1200
        ier=16
        return
 1200 continue
c 
c        use euler iterations to conduct defferred corrections
c 
        h0=t(1)-a
        do 2000 ijk=1,ncorr
c 
        call rstoneeu(ier,t,n,phik,ainte,phi0,
     1      funeva,u,h0,par1,par2,par3,m,ff,delta,rhs,a,
     2      nfcalls,nfinv,erreuler,funder,cd,cw1,cw2,
     3      dermat)
c 
        if(ier .ne. 0) return
c 
c       if this is the next to last correction - store its
c       value at the right end for subsequent extrapolation
c       and error control
c 
        if(erreuler .lt. eps22) goto 2050
 2000 continue
c 
 2050 continue
c 
c        find the solution at the right end of the interval
c 
        do 2200 j=1,m
        cdd=0
        do 2100 i=1,n
c 
        cdd=cdd+endinter(i)*phik(j,i)
 2100 continue
c 
        phi1(j)=cdd
 2200 continue
c 
        err=dsqrt(erreuler)
        return
        end
c 
c 
c 
c 
c 
        subroutine rstoneeu(ier,t,n,phik,ainte,phi0,
     1      funeva,u,h0,par1,par2,par3,m,ff,delta,rhs,a,
     2      nfcalls,nfinv,erreuler,funder,cd,cw1,cw2,
     3      dermat)
        implicit real *8 (a-h,o-z)
        save
        dimension par1(1),par2(1),par3(1)
        real *8 phik(m,n),ff(m,n),
     1      delta(m,n),phi0(m),rhs(m,n),cd(1),
     2      cw1(1),cw2(1)
        dimension t(1),hs(100),ainte(1)
c 
        real *8 dermat(1)
c 
c        one step after another, use implicit euler to solve the
c        ODE for the correction
c 
c        . . . get the integrand
c 
        ier=0
        dlarge=1.0d35
c 
        do 1100 i=1,n-1
        hs(i)=t(i+1)-t(i)
 1100 continue
c 
        do 1200 i=1,n
        nfcalls=nfcalls+1
c 
        call funeva(t(i),phik(1,i),ff(1,i),m,par1,par2,par3)
 1200 continue
c 
c       calculate the discrepancy rhs
c 
        do 1300 i=1,m
        do 1250 j=1,n
        cw1(j)=ff(i,j)
 1250 continue
c 
        call stmatve(ainte,cw1,cw2,n)
c 
        do  1270 j=1,n
        rhs(i,j)=cw2(j)
 1270 continue
 1300 continue
c 
        do 1400 i=1,n
c 
        do 1350 j=1,m
        rhs(j,i)=rhs(j,i)*u +phi0(j)-phik(j,i)
 1350 continue
 1400 continue
c 
c       calculate the correction via euler
c 
        nfinv=nfinv+1
        nrefine=0
c 
            call rstsol2(funeva,funder,m,h0,t(1),phik(1,1),
     1      dermat,delta,rhs,par1,par2,par3,
     2      par1,par2,par3,nrefine)
c 
        erreuler=0
        do 1600 i=2,n
c 
        do 1500 j=1,m
        cw1(j)=delta(j,i-1)+rhs(j,i)-rhs(j,i-1)
 1500 continue
c 
        dd=0
c 
        nfinv=nfinv+1
        nrefine=0
c 
        call rstsol2(funeva,funder,m,hs(i-1),t(i),phik(1,i),
     1      dermat,delta(1,i),cw1,par1,par2,par3,
     2      par1,par2,par3,nrefine)
c 
        do 1550 j=1,m
c 
        dd=dd+delta(j,i)* (delta(j,i))
 1550 continue
        if(erreuler .lt. dd) erreuler=dd
        if(dd .lt. dlarge) goto 1600
        ier=4
        return
 1600 continue
c 
        do 1800 i=1,n
        do 1700 j=1,m
        phik(j,i)=phik(j,i)+delta(j,i)
 1700 continue
 1800 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine stmatve(a,x,y,n)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n)
        real *8 x(1),y(1),d
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code and the beginning of
c        the non-linear codes proper.
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine rstsol1(funeva,funder,m,h,phik,tkp1,phikp1,
     1      dermat,work,rhs,pareva1,pareva2,pareva3,
     2      parder1,parder2,parder3,nrefine)
        implicit real *8 (a-h,o-z)
        save
        real *8 phik(m),phikp1(m),dermat(m,m),work(1),
     1      rhs(1),pareva1(1),pareva2(1),pareva3(1),
     2      parder1(1),parder2(1),parder3(1),rhs2(100)
c 
c        This subroutine solves the (generally, non-linear)
c        equation
c 
c        phi_{k+1} - h * F(t_{k+1}, \phi_{k+1}) = \phi_k       (1)
c 
c        via Newton method.
c 
c                   Input parameters:
c 
c  funeva - the user-specified subroutine for the evaluation of
c       the function F in (1). The calling sequence of funeva is
c 
c        funeva(t,phi,f,m,pareva1,pareva2,pareva3),
c 
c        with t, phi the arguments, m the dimensionality of the
c        problem, f the function F evaluated at the point (t,phi),
c        and the parameters pareva1, pareva2, pareva3
c        any user-specified parameters the user might want to
c        specify.
c 
c  funder - the user-specified subroutine for the evaluation
c        of the derivative of the function  F in (1) with respect
c        to phi at the point (t,phi). Please note that the
c        output of the subroutine funder is an m \times m matrix.
c        The calling sequence of funder is
c 
c      funder(t,phi,dermat,m,parder1,parder2,parder3),
c 
c        with t, phi the arguments, m the dimensionality of the
c        problem, dermat the derivative of the function F evaluated
c        at the point (t,phi), and the parameters pareva1, pareva2,
c        pareva3 any user-specified parameters the user might want
c        to specify.
c  m - the dimensionality of the problem
c  h - the parameter in (1) above (the "time step")
c  phik - the right-hand side of (1)
c  tkp1 - the parameter in (1) above (the "next time point")
c  pareva1, pareva2, pareva3 - parameters to be used by the subroutine
c        funeva
c  parder1, parder2, parder3 - parameters to be used by the subroutine
c        funder
c  nrefine - the number of newton iterations after the fisrt one.
c        Normally, nrefine should be set to 0 (a single Newton
c        step is sufficient); the possibility of more Newton
c        iterations is provided as a debugging aid, and for
c        experimentation.
c 
c                   Output parameters:
c 
c  phikp1 - the solution of the equation (1) above
c 
c                   work arrays:
c 
c  dermat - must be at least m**2 real *8 locations
c  work, rhs - must be at least m real *8 locations each
c 
c       . . . starting with phik, perform a single iteration of
c             Newton to get one step of backward euler
c 
  
  
cccc        call prin2('in rstsol1, parder3=*',parder3,m*2)
cccc        call prin2('in rstsol1, pareva3=*',pareva3,m*2)
  
        do 1100 i=1,m
        phikp1(i)=phik(i)
 1100 continue
c 
      do 2000 ijk=1,nrefine+1
c 
        call funeva(tkp1,phikp1,rhs,m,pareva1,pareva2,pareva3)
c 
        do 1200 i=1,m
c 
        rhs(i)=rhs(i)*h+phik(i)-phikp1(i)
 1200 continue
c 
cccc        call prin2('rhs=*',rhs,m*2)
c 
c       . . . get the derivative of the function funder
c             at the point (t(k+1),phik) and multiply
c             it by h
c 
cccc        call prin2('calling funder from rstsol1, phikp1=*',
cccc     1      phikp1,m*2)
  
        call funder(tkp1,phikp1,dermat,m,parder1,parder2,parder3)
  
cccc        call prin2('in rstsol1, dermat=*',dermat,m*m*2)
  
c 
c       subtract the derivative times h from the identity matrix,
c 
        do 1600 i=1,m
        do 1400 j=1,m
c 
        dermat(j,i)=-h*dermat(j,i)
 1400 continue
c 
        dermat(i,i)=dermat(i,i)+1
 1600 continue
c 
c       invert the matrix and apply the inverse to the right-hand side
c 
c 
        call rstarrmo(rhs,rhs2,m*2)
        call rcgrmsol(dermat,m,rhs2,work,rcond)
  
cccc        call prin2('rcond=*',rcond,1)
  
c 
        do 1800 i=1,m
        phikp1(i)=phikp1(i)+work(i)
 1800 continue
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
        subroutine rstsol2(funeva,funder,m,h,tkp1,phikp1,
     1      dermat,delta,rhs,pareva1,pareva2,pareva3,
     2      parder1,parder2,parder3,nrefine)
        implicit real *8 (a-h,o-z)
        save
        real *8 phikp1(m),dermat(m,m),delta(1),
     1      rhs(1),parder1(1),parder2(1),parder3(1),
     2      work(100),work2(100),work3(100),rhs2(100),delta2(100)
c 
c        This subroutine solves the (generally, non-linear)
c        "correction" equation
c 
c        \delta_{k+1} - h * [ F(t_{k+1}, \phi_{k+1}+ \delta_{k+1})
c         - F(t_{k+1}, \phi_{k+1} ) ] =   RHS            (1)
c 
c        via Newton method.
c 
c                   Input parameters:
c 
c  funeva - the user-specified subroutine for the evaluation of
c       the function F in (1). The calling sequence of funeva is
c 
c        funeva(t,phi,f,m,pareva1,pareva2,pareva3),
c 
c        with t, phi the arguments, m the dimensionality of the
c        problem, f the function F evaluated at the point (t,phi),
c        and the parameters pareva1, pareva2, pareva3
c        any user-specified parameters the user might want to
c        specify.
c 
c  funder - the user-specified subroutine for the evaluation
c        of the derivative of the function  F in (1) with respect
c        to phi at the point (t,phi). Please note that the
c        output of the subroutine funder is an m \times m matrix.
c        The calling sequence of funder is
c 
c      funder(t,phi,dermat,m,parder1,parder2,parder3),
c 
c        with t, phi the arguments, m the dimensionality of the
c        problem, dermat the derivative of the function F evaluated
c        at the point (t,phi), and the parameters pareva1, pareva2,
c        pareva3 any user-specified parameters the user might want
c        to specify.
c  m - the dimensionality of the problem
c  h - the parameter in (1) above (the "time step")
c  tkp1 - the parameter in (1) above (the "next time point")
c  phik - the right-hand side of (1)
c  pareva1, pareva2, pareva3 - parameters to be used by the subroutine
c        funeva
c  parder1, parder2, parder3 - parameters to be used by the subroutine
c        funder
c  nrefine - the number of newton iterations after the fisrt one.
c        Normally, nrefine should be set to 0 (a single Newton
c        step is sufficient); the possibility of more Newton
c        iterations is provided as a debugging aid, and for
c        experimentation.
c 
c                   Output parameters:
c 
c  delta - the solution of the equation (1) above
c 
c                   work arrays:
c 
c  dermat - must be at least m**2 real *8 locations
c 
c       . . . starting with 0, perform a single iteration of
c             Newton to get the correction delta at the point t_{k+1}
c 
        call funder(tkp1,phikp1,dermat,m,parder1,parder2,parder3)
c 
c       subtract the derivative times h from the identity matrix,
c 
        do 1600 i=1,m
        do 1400 j=1,m
c 
        dermat(j,i)=-h*dermat(j,i)
 1400 continue
c 
        dermat(i,i)=dermat(i,i)+1
 1600 continue
c 
c       invert the matrix and apply the inverse to the right-hand side
c 
        if(nrefine .gt. 0) goto 1700
c 
        call rstarrmo(rhs,rhs2,m*2)
        call rcgrmsol(dermat,m,rhs2,delta,rcond)
c 
        return
 1700 continue
c 
c       the user has requested that the solution be refined by
c       additional iterations of Newton. Do so
c 
        call rstarrmo(rhs,rhs2,m*2)
        call rcgrmsol(dermat,m,rhs2,delta,rcond)
c 
c        conduct additional Newton iterations to test the convergence
c 
        call funeva(tkp1,phikp1,work,m,pareva1,pareva2,pareva3)
c 
        do 3000 ijk=1,nrefine
  
cccc        call prinf('in newton, ijk=*',ijk,1)
c 
        do 2200 i=1,m
c 
        work2(i)=phikp1(i)+delta(i)
 2200 continue
c 
        call funeva(tkp1,work2,work3,m,pareva1,pareva2,pareva3)
c 
        call funder(tkp1,work2,dermat,m,parder1,parder2,parder3)
c 
c       subtract the derivative times h from the identity matrix,
c 
        do 2600 i=1,m
        do 2400 j=1,m
c 
        dermat(j,i)=-h*dermat(j,i)
 2400 continue
c 
        dermat(i,i)=dermat(i,i)+1
 2600 continue
c 
c       invert the matrix and apply the inverse to the right-hand side
c 
        do 2800 i=1,m
c 
        rhs2(i)=rhs(i)-(work(i)-work3(i))*h -delta(i)
 2800 continue
  
        call rcgrmsol(dermat,m,rhs2,delta2,rcond)
  
cccc        call prin2('after rcgrmsol, delta2=*',delta2,m)
  
        do 2900 i=1,m
        delta(i)=delta(i)+delta2(i)
 2900 continue
  
  
 3000 continue
  
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine rcgrmsol(b,n,x,y,rcond)
        implicit real *8 (a-h,o-z)
        save
        real *8 b(n,n),cd,x(1),y(1)
c 
c        This subroutine uses the pivoted Gram-schmidt procedure
c        to solve a linear system with complex coefficients.
c        This is a primitive routine in that it solves the system
c        with a single right-hand side, and will perform the
c        Gram-Schmidt again and again for each new right-hand
c        side. It should not be used in this regime!
c 
c                 Input parameters:
c 
c  b - the matrix of the system; destroyed by the subroutine
c  n - dimensionality of the system
c  x - the right-hand side; destroyed by the subroutine
c 
c 
c                 Output parameters:
c 
c  y - the solution of the system
c  rcond - estimate of the condition number of the matrix b
c 
c        . . . transpose the matrix b
c 
        do 1400 i=1,n
        do 1200 j=1,i
        cd=b(i,j)
        b(i,j)=b(j,i)
        b(j,i)=cd
 1200 continue
 1400 continue
c 
c        apply the Gram-Schmidt process to the matrix b
c 
        call rrcgrmso0(b,n,x,y,rcond)
c 
c       apply the adjoint of the Gram-Schmidt'ed matrix to the vector
c 
        do 2400 i=1,n
        cd=0
        do 2200 j=1,n
        cd=cd+ (b(i,j))*x(j)
 2200 continue
        y(i)=cd
 2400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine rrcgrmso0(b,n,x,rnorms,rcond)
        implicit real *8 (a-h,o-z)
        save
        real *8 b(n,n),cd,x(1)
        dimension rnorms(1)
c 
c       . . . prepare the array of values of norms
c             of columns
c 
        done=1
        dtot=0
        do 1400 i=1,n
        d=0
        do 1200 j=1,n
        d=d+b(j,i)* (b(j,i))
 1200 continue
        rnorms(i)=d
        dtot=dtot+d
 1400 continue
c 
c       . . . conduct gram-schmidt iterations
c 
        do 4000 i=1,n
c 
c       find the pivot
c 
        ipivot=i
        rn=rnorms(i)
c 
        do 2200 j=i+1,n
        if(rnorms(j) .le. rn) goto 2200
        rn=rnorms(j)
        ipivot=j
 2200 continue
 2400 continue
c 
c       put the column number ipivot in the i-th place
c 
        do 2600 j=1,n
        cd=b(j,i)
        b(j,i)=b(j,ipivot)
        b(j,ipivot)=cd
 2600 continue
c 
        cd=x(i)
        x(i)=x(ipivot)
        x(ipivot)=cd
c 
        d=rnorms(ipivot)
        rnorms(ipivot)=rnorms(i)
        rnorms(i)=d
 2700 continue
c 
c       orthogonalize the i-th column to all preceeding ones
c 
        if(i .eq. 1) goto 2790
c 
        do 2780 j=1,i-1
c 
        call rcgrmsca(b(1,i),b(1,j),n,cd)
c 
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*cd
 2770 continue
c 
        x(i)=x(i)-x(j)*cd
c 
 2780 continue
 2790 continue
c 
c       normalize the i-th column
c 
        call rcgrmsca(b(1,i),b(1,i),n,cd)
c 
c       if the norm of the remaining part of the matrix
c       is sufficiently small - terminate the process
c 
        d=cd
c 
        d=done/dsqrt(d)
        do 2800 j=1,n
        b(j,i)=b(j,i)*d
 2800 continue
c 
        x(i)=x(i)*d
c 
c        orthogonalize everything else to it
c 
        do 3200 j=i+1,n
c 
        call rcgrmsca(b(1,i),b(1,j),n,cd)
        cd= (cd)
c 
         rrn=0
        do 3000 l=1,n
        b(l,j)=b(l,j)-b(l,i)*cd
        rrn=rrn+b(l,j)* (b(l,j))
 3000 continue
        x(j)=x(j)-x(i)*cd
        rnorms(j)=dsqrt(rrn)
 3200 continue
 4000 continue
  
 4200 continue
c 
        rcond=rnorms(1)/rnorms(n)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine rcgrmsca(x,y,n,prod)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1)
c 
        prod=0
        do 1200 i=1,n
        prod=prod+x(i)* (y(i))
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine rstauinm(n,ainte,adiff,x,whts,endinter,
     1      itype,w,u)
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
c  u - the n \times n matrix converting the values of a function
c          at n Gaussian nodes into the coefficients of its
c          Legendre series
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
     1      x,whts,u,w(iv),w(iw),itype,endinter)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine rstarrmo(x,y,n)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1)
c 
        do 1200 i=1,n
        y(i)=x(i)
 1200 continue
        return
        end
