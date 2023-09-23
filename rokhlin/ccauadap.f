        implicit real *8 (a-h,o-z)
        dimension w(100 000),errs(10 000),
     1      tout(10 000),errs(10 000),
     2      ab(2,100 000)
        complex *16 clam,ima,clams(100),
     1      phiout(2,10 000),phi00(100),clam22,phi11(100)
c 
        complex *16 w2(100 000),coefs(100 000),diff(10),
     1      vals(10),ders(10)
c 
        external funeva2,funeva3,funeva4
        data ima/(0.0d0,1.0d0)/
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
         PRINT *, 'ENTER n'
         READ *,n
         CALL PRINF('n=*',n,1 )
c 
         PRINT *, 'ENTER ncorrect'
         READ *,ncorrect
         CALL PRINF('ncorrect=*',ncorrect,1 )
c 
         PRINT *, 'ENTER clam'
         READ *,clam2
         clam=clam2 *ima
         CALL PRIN2('clam=*',clam,2 )
c 
c        test the adaptive solver of the ode
c 
        aaa=-1
        bbb=-0.0001
cccc        bbb=-0.01
c 
        ifsave=1
c 
        m=2
  
        call jy01rea(-aaa,fj,fy,fj1,fy1)
  
        phi00(1)=fj+ima*fy
        phi00(2)=-(-fj1-ima*fy1)
c 
cccc        phi00(1)=fj
cccc        phi00(2)=-fj1
c 
        eps=1.0d-12
c 
        clam22=ima
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
        call ccauadap(ier,n,ncorrect,aaa,bbb,phi00,funeva4,
     1    m,par1,m,clams,eps,maxinte,ifsave,
     2    ab,ninter,tout,phiout,phi11,nfcalls,ifinit,w,
     3    lw,ltot)
c 
cccc        tt2=second()
        call prin2('time for ccauadap is*',tt2-tt1,1)
        call prinf('after first ccauadap, ier=*',ier,1)
        call prinf('after first ccauadap, nfcalls=*',nfcalls,1)
        call prinf('after first ccauadap, ltot =*',ltot,1)
  
  
        call prin2('after ccauada0, phi11=*',phi11,2)
        call prin2('after ccauada0, ab=*',ab,ninter*2)
c 
c       calculate and print the error
c 
        call errprin(tout,phiout,clams,ninter,n,m,errs,phi11,
     1      aaa,bbb,ifsave,ab)
  
        call prin2('and at the end, ab=*',ab,ninter*2)
  
  
c 
c       test the evaluation of solutions via the interpolation
c 
c 
c       . . . construct the nested Legendre decomposition of the solution
c 
        call ccexpbld(phiout,m,n,ninter,coefs,w2)
  
  
cccc        call prin2('after ccexpbld, coefs=*',coefs,n*m*ninter*2)
c 
c       for all points in tsall, evaluate the solution via
c       the Legendre expansion, and compare it to the one
c       obtained by the ODE solver
c 
        erreval=0
c 
        do 4200 i=1,n*ninter
c 
cccc        call ccexpeva(ier,coefs,ab,ninter,n,m,tsall(i),
        call ccexpeva(ier,coefs,ab,ninter,n,m,tout(i),
     1      vals,ders)
  
cccc        call prinf('i=*',i,1)
cccc        call prin2('and vals=*',vals,m*2)
cccc        call prin2('and vals via ode=*',phisall(1,i),m*2)
  
        diff(1)=vals(1)-phiout(1,i)
        diff(2)=vals(2)-phiout(2,i)
  
cccc        call prin2('and the difference = *',diff,m*2)
  
        d=abs(diff(1))+abs(diff(2))
        if(d .gt. erreval) erreval=d
 4200 continue
  
  
        call prin2('and erreval=*',erreval,1)
  
        call prinf('and again, nfcalls=*',nfcalls,1)
  
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine errprin(tout,phiout,clams,ninter,n,m,errs,
     1      phi11,a,b,ifsave,ab)
        implicit real *8 (a-h,o-z)
        save
        dimension tout(1),errs(m,1),errsb(100),
     1      rll(1000),ab(2,1)
        complex *16 phiout(m,1),cd,clams(1),phi11(1),ima
        data ima/(0.0d0,1.0d0)/
c 
c       print the values of the solution
c 
        if(ifsave .eq. 0) goto 1500
        call prin2('solution is*',phiout,n*ninter*m*2)
        call prin2('clams=*',clams,2*m)
c 
c        calculate the solution exactly and calculate the error
c 
        do 1400 i=1,n*ninter
c 
        call jy01rea(-tout(i),fj,fy,fj1,fy1)
c 
        cd=fj+ima*fy
        errs(1,i)=cdabs(phiout(1,i)-cd)
  
  
        cd=fj1+ima*fy1
        errs(2,i)=cdabs(phiout(2,i)-cd)
 1400 continue
c 
        call prin2('and errs are*',errs,n*ninter*m)
c 
 1500 continue
c 
c        check the error at the right end
c 
        call jy01rea(-b,fj,fy,fj1,fy1)
  
        cd=fj+ima*fy
  
        errsb(1)=cdabs(phi11(1)-cd)
  
        cd=(fj1+ima*fy1)
  
        errsb(2)=cdabs(phi11(2)-cd)
  
  
        call prin2('values at right end are*',phi11,4)
        call prin2('and errors at right end are*',errsb,m)
        call prinf('and ninter=*',ninter,1)
c 
c       print the lengths of the intervals
c 
        do 1800 i=1,ninter
c 
        rll(i)=ab(2,i)-ab(1,i)
 1800 continue
c 
        call prin2('and lengths of intervals are*',rll,ninter)
        return
        end
c 
c 
c 
c 
c 
        subroutine funeva2(t,phi,ff,par1,par2,clam)
        implicit real *8 (a-h,o-z)
        save
        complex *16 phi,ff,clam
c 
        ff=phi*clam
        return
        end
c 
c 
c 
c 
c 
        subroutine funeva3(t,phi,ff,par1,m,clams)
        implicit real *8 (a-h,o-z)
        save
        complex *16 phi(m),ff(m),clams(1)
c 
        do 1200 i=1,m
c 
        ff(i)=phi(i)*clams(i)
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine funeva4(t,phi,ff,par1,m,clams)
        implicit real *8 (a-h,o-z)
        save
        complex *16 phi(m),ff(m),clams(1)
c 
        ff(1)=phi(2)
        ff(2)=-( (t**2*phi(1)+t*phi(2))/t**2 )
        return
        end
c 
c 
c 
c 
c 
c 
        subroutine funeval(t,phi,ff)
        implicit real *8 (a-h,o-z)
        save
        complex *16 phi,ff,ima
        data ima/(0.0d0,1.0d0)/
c 
        call jy0rea(t,fj,fy)
  
        ff=fj+ima*fy
c 
        return
        end
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
        subroutine ccexpbld(phis,m,n,ninter,coefs,w)
        implicit real *8 (a-h,o-z)
        save
        complex *16 phis(m,n,ninter),coefs(n,m,ninter),
     1      w(100000)
c 
c       This subroutine constructs a nested Legendre expansion
c       of the function [a,b] \to \C^m, given to it in the
c       form of a table of values at a collection of nested
c       Legendre nodes. Normally, the nested Legendre tabulation
c       has been produced by the subroutine ccaunon (see) or
c       the subroutine ccauadap (see), or something similar.
c 
c                        Input parameters:
c 
c  phis - the values of the solution of (1) at the points tout
c  m - the dimensionality of the problem
c  n - the number of Gaussian nodes on each of nn subintervals in the
c      nested subdivision of the interval [a,b].
c  ninter - the number of subintervals created on the interval [a,b];
c     also the number of entries in the array ab (see above)
c 
c                    Output parameters:
c 
c  coefs - the nested Legendre expansions of the components of the
c        solution of the ODE.
c 
c                    Work arrays:
c 
c  w - must be at least max ( 4*n**2+4*n+20, 2*n**2+4*n+2*n*m+20)
c      real *8 locations long
c 
c       . . . construct the matrix converting a table of values of
c             a function at Gaussian nodes into the Legendre expansion
c             of the said function
c 
        iu=1
        lu=n*n+2
c 
        iv=iu+lu
        lv=n*n+2
c 
        ix=iv+lv
        lx=n+2
c 
        iwhts=ix+lx
        lwhts=n+2
c 
        itype=2
        call legeexps(itype,n,w(ix),w(iu),w(iv),w(iwhts) )
c 
c       allocate memory for the procedure proper
c 
        icw=iu+lu
        lcw=n*2+2
c 
        icw2=icw+lcw
        lcw2=n*2+2
c 
        icoefs0=icw2+lcw2
        lcoefs0=n*m*2+2
c 
        call ccexpbl0(phis,m,n,ninter,coefs,w(iu),
     1        w(icw),w(icw2),w(icoefs0) )
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine ccexpbl0(phis,m,n,ninter,
     1      coefs,u,cw,cw2,coefs0)
        implicit real *8 (a-h,o-z)
        save
        complex *16 phis(m,n,ninter),coefs(n,m,ninter),
     1      u(1),cw(1),cw2(1),coefs0(n,m)
c 
c       one subinterval after another, convert tables of values
c       of the solution into their Legendre expansions on the
c       subintervals
c 
        do 2000 i=1,ninter
c 
        do 1600 j=1,m
c 
        do 1400 k=1,n
c 
        cw(k)=phis(j,k,i)
 1400 continue
c 
        call ccmatve(u,cw,cw2,n)
c 
        do 1500 k=1,n
c 
        coefs0(k,j)=cw2(k)
 1500 continue
 1600 continue
c 
        call ccarrmo(coefs0,coefs(1,1,i),n*m*2)
 2000 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine ccexpeva(ier,coefs,ab,ninter,n,m,x,vals,ders)
        implicit real *8 (a-h,o-z)
        save
        complex *16 coefs(n,m,1),val,der,vals(1),ders(1)
        dimension ab(2,1)
        data intold/-10/,ninteold/-10/
c 
c        this subroutine evaluates at the point x \in [a,b] the
c        solution of the ODE produced by a prior call to the subroutine
c        ccaunon (see) or ccauadap (see), and post-processed (converted
c        into Legendre expansions) by a prior call to the subroutine
c        ccexpbld (see). More specifically, the solution is specified
c        by nested Legendre expansions of its components, contained in
c        the input array coefs.
c 
c                   Input parameters:
c 
c  coefs - the nested Legendre expansions of the components of the
c        solution of the ODE.
c  ab - the array of ends of subintervals created. For each i=1,2,...,nn,
c     ab(1,i) is the left end of the i-th subinterval; ab(2,i) is
c     the right end of the i-th subinterval
c  ninter - the number of subintervals created on the interval [a,b];
c     also the number of entries in the array ab (see above)
c  n - the number of Gaussian nodes on each of nn subintervals in the
c      nested subdivision of the interval [a,b].
c  m - the dimensionality of the problem
c  x - the point where the solution is to be evaluated; must be on
c      the interval [a,b]
c 
c                    Output parameters:
c 
c  vals - the value of the solution at the point x
c  ders - the derivative of the solution at the point x
c 
c       . . . find the subinterval in which the point x lives
c 
        ier=0
c 
        if(intold .lt. 0) goto 1100
        if(ninteold .ne. ninter) goto 1100
c 
        if( (x .lt. ab(1,intold)) .or. (x .gt. ab(2,intold)) ) goto 1100
c 
        int=intold
        goto 1600
c 
 1100 continue
c 
c       the point is not on the same subinterval where the preceding
c       one was. Use bisection to find the proper subinterval
c 
        i1=1
        i2=ninter
c 
        if( (x .ge. ab(1,1)) .and. (x .le. ab(2,ninter)) ) goto 1120
c 
        ier=4
        return
 1120 continue
c 
        ifout=0
c 
        do 1150 i=1,ninter
c 
        i3=(i1+i2)/2
c 
        if(x .gt. ab(2,i3)) goto 1130
c 
        i2=i3
        goto 1140
c 
 1130 continue
c 
        i1=i3
 1140 continue
c 
        if(ifout .eq. 1) goto 1160
        if(i2-i1 .le. 1) ifout=1
 1150 continue
c 
 1160 continue
c 
        intnum=i2
        intold=intnum
        ninteold=ninter
c 
 1600 continue
c 
c       the point x lives in the interval number intnum.
c       evaluate the expansions at this point
c 
        u=2/(ab(2,intnum)-ab(1,intnum))
        v=1-ab(2,intnum)*u
c 
        t=u*x+v
c 
        do 2200 i=1,m
c 
        call legecFDE(t,VAL,der,coefs(1,i,intnum),n-1)
c 
        der=der*u
c 
        vals(i)=val
        ders(i)=der
 2200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine ccauadap(ier,n,ncorrect,a,b,phi0,funeva,m,
     1    par1,par2,par3,eps,maxinte,ifsave,
     2    ab,ninter,tout,phiout,phi1,nfcalls,ifinit,
     3      w,lenw,ltot)
        implicit real *8 (a-h,o-z)
        save
        dimension tout(1),w(1),par1(1),par2(1),
     1      par3(1),ab(2,1)
        complex *16 phi0(m),phi1(m),phiout(m,1)
c 
c       this subroutine solves adaptively the differential
c       equation
c 
c       phi'(t) = f(t,phi(t)).                                   (1)
c 
c        the adaptive scheme is based on the spectral deferred
c        corrections, and is of arbitrary (user-specified) order.
c 
c                       input parameters:
c 
c  n - the number of nodes of the spectral discretization of each
c        subinterval
c  ncorrect - the number of deferred corrections.
c        IMPORTANT note. The order of the method is equal to
c        min(n,ncorrect+1); the actual optimum choice of the parameters
c        n, m appropriate to a particular problem and accuracy are
c        non-trivial, and have not been investigated in sufficient detail
c  a,b - the ends of the interval on which the equation is to be solved
c  phi0 - the initial value for the Cauchy problem (specified at the
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
c  m - the dimensionality of the problem
c  par1,par2,par3 - the parameters to be used by the subroutine funeva
c  eps - the precision to which the solution is to be evaluated
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
c         in this case is needless to say, this might means that the
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
c  ninter - the number of refinements actually constructed
c  tout - the points at which the solution of (1) has been calculated
c  phiout - the values of the solution of (1) at the points tout
c  phi1 - the value of the solution of (1) at the point b
c  nfcalls - the number of times the user-supplied subroutine funeva
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
c 
        icw2=icw1+lcw1
        lcw2=n*nc+2
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
        maxiter=100 000
c 
        call ccauada1(ier,n,a,b,phi0,funeva,
     1      par1,par2,par3,tout,phiout,ncorrect,
     2      ab,maxiter,eps,ninter,ifinit,phi1,
     3      w(iainte),w(iw),maxinte,m,ifsave,
     4      w(iphik1),w(iff),w(idelta),w(irhs),w(icdnon),
     5      w(iphi0),w(iuuu),w(icw1),w(icw2),nfcalls)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine ccauada1(ier,n,aaa,bbb,phi00,funeva,
     1    par1,par2,par3,tout,phiout,ncorrect,ab,
     2      maxiter,eps,ninter,ifinit,phi11,ainte,w,
     3      maxinte,m,ifsave,phik1,ff,delta,rhs,cdnon,
     4      phi0,uuu,cw1,cw2,nfcalls)
        implicit real *8 (a-h,o-z)
        save
        external funeva
        dimension par1(1),par2(1),par3(1)
        dimension ainte(1),ts(100),whts(100),endinter(100),
     1    w(1),t(100),tout(1),uuu(1),ab(2,1)
        complex *16 phi00(m),phi0(1),phik1(1),phiout(m,1),
     1    phi11(m),ff(1),delta(1),rhs(1),cdnon(1),
     2      cw1(1),cw2(1)
c 
c       this subroutine solves adaptively the differential
c       equation
c 
c       phi'(t) = f(t,phi(t)).                                   (1)
c 
c        the adaptive scheme is based on the spectral deferred
c        corrections, and is of arbitrary (user-specified) order.
c 
c                       input parameters:
c 
c  n - the number of nodes of the spectral discretization of each
c        subinterval
c  aaa,bbb - the ends of the interval on which the equation is to be solved
c  phi00 - the initial value for the Cauchy problem (specified at the
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
c  ncorrect - the number of deferred corrections.
c        IMPORTANT note. The order of the method is equal to
c        min(n,ncorrect+1); the actual optimum choice of the parameters
c        n, m appropriate to a particular problem and accuracy are
c        non-trivial, and have not been investigated in sufficient detail
c  maxiter - the maximum total number of intervals on which the
c        subroutine is permitted to {\bf attempt} to construct
c        deferred corrections. Should be LARGE
c  eps - the precision to which the solution is to be evaluated
c  ifinit - the parameter telling the subroutine whether it should
c     construct the spectral integration matrix and the associated
c     parameters. ifinit=1 will cause the subroutine to construct
c     the the integration matrix (it is somewhat expensive).
c     ifinit=0 will cause the subroutine to skip the creation of the
c     integration matrix. In this it the parameters ainte, uuu
c     must not have been changed since the preceding call.
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
c         in this case is needless to say, this might means that the
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
c  phik1,ff,delta,rhs - must be at least n*m+2 complex *16 locations long
c  cd,cw1,cw2 - must be at least m+2 complex *16 locations long
c 
c 
c        . . . if needed, construct the matrix of indefinite
c              integration at the Gaussian nodes
c 
        if(ifinit .eq. 0) goto 1200
        itype=1
        call ccauinmt(n,ainte,adiff,ts,whts,endinter,
     1      itype,w,uuu)
c 
 1200 continue
c 
c       start the process of adaptive solution of the ode
c 
        ab(1,1)=aaa
        ab(2,1)=aaa
  
  
        step=bbb-aaa
c 
        small=step*1.0d-8
c 
        call ccarrmo(phi00,phi0,m*2)
  
        ifnew=1
        ninter=1
        numpts=0
        ifout=0
c 
c       solve the ode adaptively
c 
        ndouble=3
        nmore=0
        do 3000 i=1,maxiter
c 
c        construct the solution on the new interval
c 
cccc        a=bb(ninter)
        a=ab(2,ninter)
        b=a+step
c 
        bbbb=bbb-step/10
        if(b .ge. bbbb) ifout=1
        if(b .ge. bbbb) b=bbb
        if(b .ge. bbbb) step=b-a
c 
        call ccnonstd(jer,ts,n,phi0,funeva,par1,par2,par3,
     1      a,b,t,phik1,ncorrect,ainte,whts,phi11,m,
     3      ff,delta,rhs,cdnon,err,eps,cw1,cw2,nfcalls)
c 
        if(jer .ne. 0) call prinf(
     1      'in ccauada0 after ccnonstd, jer=*',jer,1)
c 
c        if the step is so big that the corrected Euler has
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
        call ccdisrep(uuu,n,m,phik1,ff,discr,val,cw1,cw2)
c 
        if( (err .gt. eps*val) .or.
     1      (discr .gt. eps*val) ) goto 2000
c 
        ab(1,ninter+1)=ab(2,ninter)
        ab(2,ninter+1)=ab(1,ninter+1)+step
  
  
        ninter=ninter+1
c 
        if(ninter .lt. maxinte) goto 1800
c 
        ier=4
        return
 1800 continue
c 
        call ccarrmo(phi11,phi0,m*2)
c 
        if(ifsave .eq. 0) goto 1900
        call ccarrmo(phik1,phiout(1,numpts+1),n*m*2)
        call ccarrmo(t,tout(numpts+1),n)
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
c       array ab
c 
        do 4200 i=1,ninter-1
c 
        ab(1,i)=ab(1,i+1)
        ab(2,i)=ab(2,i+1)
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
        subroutine ccdisrep(uuu,n,m,phik,w,discr,val,cw1,cw2)
        implicit real *8 (a-h,o-z)
        save
        dimension uuu(1)
        complex *16 w(m,n),phik(m,n),cw1(1),cw2(1)
c 
c       construct the Legendre transformation of the solution
c       tabulated at the Gaussian nodes
c 
        val=0
        do 2000 i=1,m
        do 1600 j=1,n
        val=val+phik(i,j)*dconjg(phik(i,j))
        cw1(j)=phik(i,j)
 1600 continue
c 
        call ccmatve(uuu,cw1,cw2,n)
c 
        do  1800 j=1,n
        w(i,j)=cw2(j)
 1800 continue
cccc         call prin2('in ccdisrep, cw2=*',cw2,n*2)
  
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
        discr=discr+w(i,n-j+1)*dconjg(w(i,n-j+1))
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
        subroutine ccnonstd(ier,ts,n,phi0,funeva,
     1      par1,par2,par3,a,b,t,
     2      phik,ncorrect,ainte,whts,phi1,m,
     3      ff,delta,rhs,cd,err,eps,cw1,cw2,nfcalls)
        implicit real *8 (a-h,o-z)
        save
        dimension par1(1),par2(1),par3(1)
        dimension t(1),ainte(n,n),whts(n),ts(n)
        complex *16 phi0(m),phik(m,n),cd(1),
     1      phi1(m),ff(m,n),delta(m,n),rhs(m,n),cw1(1),cw2(1)
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
c       . . . use euler (with trapezoidal corrections) to construct
c             the initial approximation on the short interval
c 
        d=t(1)-a
        call funeva(a,phi0,cd,par1,par2,par3)
        nfcalls=nfcalls+1
c 
        do 1120 j=1,m
c 
        phik(j,1)=phi0(j)+cd(j)*d
 1120 continue
c 
        call funeva(a,phik(1,1),cw1,par1,par2,par3)
        nfcalls=nfcalls+1
        do 1140 j=1,m
c 
        phik(j,1)=phi0(j)+(cd(j)+cw1(j))*d/2
 1140 continue
c 
        do 1200 i=2,n
c 
        call ccarrmo(cw1,cd,m*2)
c 
        dd=0
        do 1150 j=1,m
c 
        phik(j,i)=phik(j,i-1)+(t(i)-t(i-1))*cd(j)
c 
        dd=dd+phik(j,i)*dconjg(phik(j,i))
c 
 1150 continue
c 
        nfcalls=nfcalls+1
        call funeva(t(i),phik(1,i),cw1,par1,par2,par3)
c 
        dd=0
        do 1170 j=1,m
c 
        phik(j,i)=phik(j,i-1)+(t(i)-t(i-1))*(cd(j)+cw1(j))/2
c 
        dd=dd+phik(j,i)*dconjg(phik(j,i))
c 
 1170 continue
c 
        if(dd .lt. dlarge) goto 1200
        ier=16
        return
 1200 continue
c 
c        use euler iterations to conduct deferred corrections
c 
        do 2000 ijk=1,ncorrect
c 
        call cconeeul(ier,t,n,phik,ainte,phi0,funeva,u,
     1      par1,par2,par3,m,ff,delta,rhs,erreuler,cw1,cw2,
     2      nfcalls,cd)
c 
        if(ier .ne. 0) call prinf('after coneeeuler, ier=*',ier,1)
cccc        call prin2('after coneeeuler, erreuler=*',erreuler,1)
        if(ier .ne. 0) return
c 
        if(erreuler .lt. eps22) goto 2050
 2000 continue
c 
 2050 continue
c 
c        find the solution at the right end of the interval
c 
        call ccarrmo(phi0,phi1,m*2)
c 
        do 2200 i=1,n
c 
        nfcalls=nfcalls+1
        call funeva(t(i),phik(1,i),cd,par1,par2,par3)
c 
        call ccarrmo(cd,ff(1,i),m*2)
c 
        do 2100 j=1,m
c 
        phi1(j)=phi1(j)+whts(i)*cd(j)*u
 2100 continue
c 
 2200 continue
c 
        err=dsqrt(erreuler)
c 
c 
c       calculate the discrepancies (of the Picard equation)
c       at all nodes on the interval
c 
        do 3000 i=1,m
c 
        do 2600 j=1,n
c 
        cw1(j)=ff(i,j)
 2600 continue
c 
        call ccmatve(ainte,cw1,cw2,n)
c 
        do 2800 j=1,n
c 
        ff(i,j)=cw2(j)*u
 2800 continue
 3000 continue
c 
        errpic=0
  
        do 3400 i=1,n
        do 3200 j=1,m
c 
        ff(j,i)=ff(j,i)+phi0(j)-phik(j,i)
c 
        if(abs(ff(j,i)) .gt. errpic) errpic=abs(ff(j,i))
 3200 continue
 3400 continue
c 
        err=errpic
        return
        end
c 
c 
c 
c 
c 
        subroutine cconeeul(ier,t,n,phik,ainte,phi0,
     1      funeva,u,par1,par2,par3,m,ff,delta,rhs,
     2      erreuler,cw1,cw2,nfcalls,cd)
        implicit real *8 (a-h,o-z)
        save
        dimension par1(1),par2(1),par3(1)
        complex *16 phik(m,n),ff(m,n),delta(m,n),phi0(m),
     2      cw1(1),cw2(1),rhs(m,n),cd(1)
        dimension t(1),hs(200),ainte(1)
c 
c        one step after another, use euler to solve the
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
        call funeva(t(i),phik(1,i),ff(1,i),par1,par2,par3)
 1200 continue
c 
c       calculate the discrepancy rhs
c 
        do 1300 i=1,m
        do 1250 j=1,n
        cw1(j)=ff(i,j)
 1250 continue
c 
        call ccmatve(ainte,cw1,cw2,n)
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
        call ccarrmo(rhs,delta,m*2)
c 
        erreuler=0
c 
        do 1600 i=2,n
  
        if(i .ne. 2) goto 1520
c 
        do 1500 j=1,m
        cw1(j)=phik(j,i-1)+delta(j,i-1)
 1500 continue
c 
        nfcalls=nfcalls+1
        call funeva(t(i-1),cw1,cd,par1,par2,par3)
 1520 continue
c 
        if(i .ne. 2) call ccarrmo(cw2,cd,m*2)
c 
        dd=0
        do 1550 j=1,m
        delta(j,i)=delta(j,i-1)+hs(i-1)*(cd(j)-ff(j,i-1))
c 
        delta(j,i)=delta(j,i)+(rhs(j,i)-rhs(j,i-1))
c 
        dd=dd+delta(j,i)*dconjg(delta(j,i))
 1550 continue
c 
        do 1560 j=1,m
        cw1(j)=phik(j,i)+delta(j,i)
 1560 continue
c 
        nfcalls=nfcalls+1
        call funeva(t(i),cw1,cw2,par1,par2,par3)
c 
        dd=0
        do 1580 j=1,m
c 
        delta(j,i)=delta(j,i-1)+hs(i-1)*
     1      ( (cd(j)-ff(j,i-1)) + (cw2(j)-ff(j,i)) )/2
c 
        delta(j,i)=delta(j,i)+(rhs(j,i)-rhs(j,i-1))
c 
        dd=dd+delta(j,i)*dconjg(delta(j,i))
 1580 continue
c 
        if(erreuler .lt. dd)  erreuler=dd
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
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine ccauinmt(n,ainte,adiff,x,whts,endinter,
     1      itype,w,u)
        implicit real *8 (a-h,o-z)
        save
        dimension ainte(1),w(1),x(1),whts(1),adiff(1),endinter(1)
c 
c        for the user-specified n, this subroutine constructs
c        the matrices of spectral indefinite integration and/or
c        spectral differentiation on the n Gaussian nodes
c        on the interval [-1,1]. Actually, this is only a
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
c                           Output parameters:
c 
c  ainte - the matrix of spectral indefinite integration on
c          the Gaussian nodes
c  adiff - the matrix of spectral differentiation on
c          the Gaussian nodes
c  x - the n Gaussian nodes on the interval [-1,1]
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
        subroutine ccarrmo(x,y,n)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1)
c 
        do 1200 i=1,n
        y(i)=x(i)
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine ccmatve(a,x,y,n)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n)
        complex *16 x(1),y(1),d
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
  
