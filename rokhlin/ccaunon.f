        implicit real *8 (a-h,o-z)
        dimension tsall(10 000),errend(100),w(500 000),
     1      errsall(2,100 000),ab(2,10000)
        complex *16 phisall(2,100 000),phi00(2),phi1(2),
     1      clams(2),ima,valrend(100),phisall2(2,100 000),
     2      coefs(100 000),vals(100),ders(100),diff(10),
     3      w2(100 000)
c 
        external funeva3,funeva4
c 
        data ima/(0.0d0,1.0d0)/
c 
  
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
c 
         PRINT *, 'ENTER ninter'
         READ *,ninter
         CALL PRINF('ninter=*',ninter,1 )
c 
         a=0
         b=1
c 
         ifinit=1
  
        m=2
  
c 
        clams(1)=1
        clams(2)=ima*2
        eps=1.0d-12
        eps=1.0d-18
  
        phi00(1)=1
        phi00(2)=1
c 
  
        call prin2('before ccaunon, clams=*',clams,4)
c 
        lenw=500 000
  
  
         call ccaunon(ier,funeva4,m,a,b,phi00,
     1       ninter,n,ncorrect,eps,ifinit,
     2       ab,tsall,phisall,phi1,errpic,
     3       par1,m,clams,nfcalls,w,lenw,keep,lused)
c 
        call prinf('after ccaunon, ier=*',ier,1)
  
  
        call funinfo4(ncalled)
  
        call prinf('and ncalled obtained via funinfo is*',ncalled,1)
        call prinf('and nfcalls from ccaunon is*',nfcalls,1)
        call prin2('and errpic from ccaunon is*',errpic,1)
c 
c       evaluate the error
c 
        valrend(1)=exp(clams(1)/3)
        valrend(2)=exp(clams(2)/3)
c 
        errend(1)=abs(phi1(1)-valrend(1))
        errend(2)=abs(phi1(2)-valrend(2))
c 
        call prin2('and errend=*',errend,2)
        call prin2('and phi1=*',phi1,4)
        call prin2('and valrend=*',valrend,4)
c 
c       now, evaluate the solution analytically at all nodes in the
c       discretization, and calculate the error
c 
        nall=n*ninter
        errstot=0
c 
        do 3200 i=1,nall
c 
        phisall2(1,i)=exp(clams(1)*tsall(i))
        phisall2(2,i)=exp(clams(2)*tsall(i))
c 
        phisall2(1,i)=exp(clams(1)*tsall(i)**3/3)
        phisall2(2,i)=exp(clams(2)*tsall(i)**3/3)
c 
        errsall(1,i)=abs(phisall2(1,i)-phisall(1,i))
        errsall(2,i)=abs(phisall2(2,i)-phisall(2,i))
c 
        errstot=errstot+errsall(1,i)**2+errsall(2,i)**2
c 
 3200 continue
c 
        errstot=sqrt(errstot)
  
  
        call ccarrmo(phisall,phisall2,n*m*2*ninter+10)
  
  
        call prin2('and errtot =*',errstot,1)
  
c 
c       construct the nested Legendre decomposition of the solution
c 
        call ccexpbld(phisall2,m,n,ninter,
     1      phisall2,w2)
  
  
        call prin2('after ccexpbld, phisall=*',phisall,n*m*ninter*2)
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
        call ccexpeva(ier,phisall2,ab,ninter,n,m,tsall(i),
     1      vals,ders)
  
cccc        call prinf('i=*',i,1)
cccc        call prin2('and vals=*',vals,m*2)
cccc        call prin2('and vals via ode=*',phisall(1,i),m*2)
  
        diff(1)=vals(1)-phisall(1,i)
        diff(2)=vals(2)-phisall(2,i)
  
cccc        call prin2('and the difference = *',diff,m*2)
  
        d=abs(diff(1))+abs(diff(2))
        if(d .gt. erreval) erreval=d
 4200 continue
  
  
        call prin2('and erreval=*',erreval,1)
        stop
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
        data ncalled/0/
c 
        ncalled=ncalled+1
c 
        do 1200 i=1,m
c 
        ff(i)=phi(i)*clams(i)
 1200 continue
        return
c 
c 
c 
c 
        entry funinfo3(ncalled7)
c 
        ncalled7=ncalled
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
        data ncalled/0/
c 
        ncalled=ncalled+1
c 
        do 1200 i=1,m
c 
        ff(i)=clams(i)*t**2 * phi(i)
 1200 continue
        return
c 
c 
c 
c 
        entry funinfo4(ncalled7)
c 
        ncalled7=ncalled
        return
        end
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c       This is the end of the debugging code and the beginning of the
c       ODE solver code proper
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
c       Legendre nodes. Normally, the nested Legnendre tabulation
c       has been produced by the subroutine ccaunon (see) or
c       something similar.
c 
c                        Input parameters:
c 
c 
c  phisall - the values of the solution of (1) at the points tout
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
c       . . . construct the matrix of Legendre expansion
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
c        this subroutine evaluates at the point x \in [a,b]
c        the solution of the ODE produced by a prior call to
c        the subroutine ccaunon (see) and post-processed (converted
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
  
cccc        call prinf('doinng bisection, ninteold=*',ninteold,1)
  
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
        subroutine ccaunon(ier,funeva,m,a,b,phi00,
     1     ninter,n,ncorrect,eps,ifinit,
     2      ab,tsall,phisall,phi1,errpic,
     3      par1,par2,par3,nfcalls,w,lenw,keep,lused)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension w(1),ab(2,1),tsall(n,1),
     1      par1(1),par2(1),par3(1)
c 
        complex *16 phi0(2),phi1(1),phi00(2),phisall(m,1)
c 
c       this subroutine solves non-adaptively the differential
c       equation
c 
c       phi'(t) = f(t,phi(t)).                                   (1)
c 
c        the scheme is based on the spectral defferred
c        corrections, and is of arbitrary (user-specified) order.
c 
c                       input parameters:
c 
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
c         END OF DESCRIPTION OF PARAMETERS OF FUNEVA
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
  
c  m - the dimensionality of the problem
c  a,b - the ends of the interval on which the equation is to be solved
c  phi00 - the initial value for the cauchy problem (specified at the
c        point a)
c  ninter - the number of subintervals into which the interval
c        [a,b] will be subdivided; on each subinterval, deferred
c        corrections are used to solve the Picard equation associated
c        with the equation (1) above
c  n - the number of nodes of the spectral discretization of each
c        subinterval
c  ncorrect - the number of defferred corrections.
c        IMPORTANT note. The order of the method is equal to
c        min(n,ncorrect+1); the actual optimum choice of the parameters
c        n, m appropriate to a particular problem and accuracy are
c        non-trivial, and have not been investigated in sufficient detail
c  eps - the precision to which the solution to the discretized Picard
c        equation is to be evaluated on each subinterval. PLEASE NOTE
c        THAT THE PRECISION IS ONLY DETERMINED IN TERMS OF THE DISCRETIZED
c        EQUATION; THE FINENESS (AND THEREFORE, ACCURACY) OF DISCRETIZATION
C        DETERMINED BY THE USER-PROVIDED PARAMETERS NINTER, N, AND IS NOT
C        UNDER THE CONTROL OF THE SUBROUTINE.
c  ifinit - the parameter telling the subroutine whether it should
c        construct the spectral integration matrix and the associated
c        parameters. ifinit=1 will cause the subroutine to construct
c        the parameters (it is somewhat expensive). ifinit=0 will cause
c        the subroutine to skip the creation of the parameters. In this
c        case, the first keep elements of the array w must not have
c        been chenged since the preceding call.
c  par1,par2,par3 - the parameters (variables, arrays, integer, real,
c        whatever) to be used by the subroutine funeva
c  lenw - the amount of storage (in real *8 locations) provided in the
c        array w; it is used by the subroutine to bomb if the storage is
c        insufficient (also, see the parameter ier below)
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
c  w - this is the input parameter only is ifinit has been set to 0.
c         It contains the spectral integration matrix on the n Gaussian
c         nodes, and several related parameters
c 
c                         output parameters:
c 
c  ier - error return code.
c     ier=0 means successful execution
c     ier=4 means that on at least one subinterval, the requested
c         precision was not achieved after tncorrect corrections.
c         This error is not viewed as a fatal one
c     ier=16 means that at some point, the process exploded for
c         stability reasons. This is a fatal error
c     ier=1024 and ier=2048 mean that the amount of storage provided by
c                the user in the array w (specified by the parameter
c                lenw) is insufficient. Obviously, this fatal
c  ab - a 2 \times ninter real *8 array containing the ends of the
c       intervals into which the interval [a,b] has been subdivided
c  tsall - the points at which the solution of (1) has been calculated
c  phisall - the values of the solution of (1) at the points tout
c  phi1 - the value of the solution of (1) at the point b
c  errpic - the maximum error of the solution of the DISCRETIZED Picard
c       equation at any node on any of the subintervals
c  nfcalls - the number of times the user-supplied subroutine funeva
c        (see above) has been called
c  w - this is an output parameter only if ifinit=1 has been specified
c  keep - the number of elements in array w that should be kept
c       between a call to this subroutine with ifinit=1 and subsequent
c       calls with ifinit=0
c  lused - the maximum number of elements in array w that was used by
c       this subroutine at any moment during execution
c 
c                         work arrays:
c 
c  w - must be sufficiently long
c 
c 
        ier=0
        nfcalls=0
c 
c        . . . if the user so requested, construct the matrix
c              of indefinite integration at the Gaussian nodes
c 
        if(ifinit .eq. 0) goto 1200
c 
c           . . . allocate memory
c 
        iainte=10
        lainte=n**2+2
c 
        its=iainte+lainte
        lts=n+2
c 
        iwhts=its+lts
        lwhts=n+2
c 
        keep=iwhts+lwhts
c 
        iendinte=iwhts+lwhts
        lendinte=n+2
c 
        iw=iendinte+lendinte
        lw=3* n**2 + 2*n +50
c 
        lused=iw+lw
c 
        if(lused .le. lenw) goto 1100
c 
        ier=2048
        return
 1100 continue
c 
        itype=1
c 
        call legeinmt(n,w(iainte),adiff,w(its),w(iwhts),w(iendinte),
     1      itype,w(iw))
c 
 1200 continue
c 
c       construct the subdivision of the big interval into elementary
c       subintervals
c 
        h=(b-a)/ninter
c 
        ab(1,1)=a
        ab(2,1)=+h
  
        do 1400 i=2,ninter
c 
        ab(1,i)=ab(2,i-1)
        ab(2,i)=ab(1,i)+h
 1400 continue
c 
cccc        call prin2('ab as created*',ab,ninter*2)
c 
c       construct Gaussian discretizations of the subdivided
c       subintervals
c 
        do 1800 i=1,ninter
c 
        alpha=(ab(2,i)-ab(1,i))/2
        beta=(ab(2,i)+ab(1,i))/2
c 
        do 1600 j=1,n
c 
        tsall(j,i)=alpha*w(its+j-1)+beta
  
 1600 continue
 1800 continue
c 
c       one interval after another, use the corrected euler
c       combined with deferred correction to solve the cauchy
c       problem
c 
        call ccarrmo(phi00,phi1,m*2)
c 
c       . . . allocate memory
c 
        it=iwhts+lwhts
        lt=n+2
c 
        iff=it+lt
        lff=m*n*2+2
c 
        idelta=iff+lff
        ldelta=m*n*2+2
c 
        irhs=idelta+ldelta
        lrhs=m*n*2+2
c 
        icd=irhs+lrhs
        lcd=m*2+2
c 
        icw1=icd+lcd
        lcw1=(m+n)*2+2
c 
        icw2=icw1+lcw1
        lcw2=(m+n)*2+2
c 
        lused2=icw2+lcw2
c 
        if(lused2 .gt. lused) lused=lused2
c 
        if(lused .le. lenw) goto 2200
c 
        ier=1024
        return
 2200 continue
c 
        errpic=0
        do 3000 k=1,ninter
c 
        call ccarrmo(phi1,phi0,m*2)
c 
        iinter=(k-1)*n+1
        call ccaunon0(jer,w(its),n,phi0,funeva,
     1      par1,par2,par3,ab(1,k),ab(2,k),w(it),
     2      phisall(1,iinter),ncorrect,w(iainte),w(iwhts),
     3      phi1,m,w(iff),w(idelta),w(irhs),w(icd),err,eps,
     4      w(icw1),w(icw2),nfcalls)
c 
cccc        call prin2('after ccaunon0, err=*',err,1)
c 
        if(errpic .lt. err) errpic=err
c 
        if(jer .eq. 4) ier=4
c 
        if(jer .lt. 5) goto 3000
c 
        ier=128
        return
c 
 3000 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine ccaunon0(ier,ts,n,phi0,funeva,
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
cccc        call prinf('in ccaunon0, ijk=*',ijk,1)
  
        call cconeeul(jer,t,n,phik,ainte,phi0,funeva,u,
     1      par1,par2,par3,m,ff,delta,rhs,erreuler,cw1,cw2,
     2      nfcalls,cd)
c 
        if(jer .eq. 0) goto 1400
  
c 
        call prinf('after coneeeuler, jer=*',jer,1)
        ier=16
        return
 1400 continue
c 
cccc       call prinf(' ijk=*',ijk,1)
cccc       call prin2(' erreuler=*',erreuler,1)
c 
        if(erreuler .lt. eps22) goto 2050
 2000 continue
c 
        call prin2('corrections failed to converge, erreuler=*',
     1      erreuler,1)
        call prin2('while eps22=*',eps22,1)
c 
        ier=4
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
c 
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
c 
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
        dd=dd+delta(j,i)*dconjg(delta(j,i))
 1580 continue
c 
cccc          call prin2('dd=*',dd,1)
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
