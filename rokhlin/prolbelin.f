        implicit real *8 (a-h,o-z)
        real *8 w(200 000),
     1      ts(10 000),fs(10 000),fsm1(10 000),fs0(10 000),
     2      coesbell(10 000),coesder(10 000),der2(10 000),
     3      diffs(10 000),tstest(10 000),bells(10 000),ders(10 000),
     4      rnums(10 000),coefs(10 000),abspro(10 000)
c 
        complex *16 ima,cd,project(10 000)
c 
        data ima/(0.0d0,1.0d0)/
  
  
cccc        data coesbell/10000*77.0/,coesder/10000*77.0/
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
cccc        PRINT *, 'ENTER c'
cccc        READ *,c
cccc        CALL PRIN2('c=*',c,1 )
c 
        PRINT *, 'ENTER epslog'
        READ *,epslog
        CALL PRIN2('epslog=*',epslog,1 )
  
        ntest=100
        lenw=200 000
  
        eps=10.0**(-epslog)
  
        do 1200 i=1,10 000
c 
        coesbell(i)=77
        coesder(i)=77
 1200 continue
  
c 
c       construct the Legendre expansions of the bell and its derivative
c 
  
        call prolbelin(ier,eps,coesbell,coesder,nterms,c,w,lenw,lused)
c 
        call prinf('after prolbelin, ier=*',ier,1)
c 
        if(ier .ne. 0) stop
c 
c        evaluate the integrated expasnion on the interval
c        [-1,1], and plot the result
c 
        h=2
        h=h/(ntest-1)
        do 1400 i=1,ntest
c 
        ts(i)=(i-1)*h-1
  
cccc        if(i .eq. ntest) ts(i)=ts(i)-1.0d-12
c 
        call legeexev(ts(i),fs(i),coesbell,nterms)
        call legeexev(ts(i),fs0(i),coesder,nterms)
  
        fsm1(i)=fs(i)-1
  
  
 1400 continue
c 
        iw=21
        call lotagraph(iw,ts,fs,ntest,
     1      'indefinite integral of \psi_0*')
  
  
        call prin2('and the resulting bell is*',fs,ntest)
  
        call prin2('and the derivative of the bell is*',fs0,ntest)
        call prin2('and the resulting bell-1 is*',fsm1,ntest)
        call prin2('and coesbell=*',coesbell,nterms+5)
        call prin2('and coesder=*',coesder,nterms+5)
        call prinf('and nterms=*',nterms,1)
        call prin2('and c=*',c,1)
c 
c        test the derivative of the bell
c 
        do 1600 i=2,ntest-1
c 
        der2(i)=(fs(i+1)-fs(i-1))/h/2
        diffs(i)=fs0(i)-der2(i)
 1600 continue
  
cccc        call prin2('and the derivative obtained numerically is*',
cccc     1       der2,ntest)
  
        call prin2('and the differences are*',diffs,ntest)
        call prin2('while h=*',h,1)
  
        call prin2('and again, derivatives are*',fs0,ntest)
c 
c       calculate the scaled bell on the interval [-1,1]
c 
        h=2
        h=h/ntest
        do 2200 i=1,ntest
c 
        tstest(i)=(i-1)*h-1
 2200 continue
c 
cccc        coesbell(1)=0
  
        call prin2('and tstest as created*',tstest,ntest)
c 
        do 2400 i=1,ntest
c 
        call prolbell(coesbell,coesder,nterms,tstest(i),
     1      bells(i),ders(i) )
c 
        ders(i)=-ders(i)
 2400 continue
  
        call prin2('and bells=*',bells,ntest)
  
        call prin2('and ders=*',ders,ntest)
c 
c       calculate projections of the bell on a bunch of Fourier modes
c 
        done=1
        pi=atan(done)*4
c 
        nfour=30
c 
        do 3400 i=0,nfour
c 
        cd=0
        do 3200 j=1,ntest
c 
        cd=cd+exp(ima*pi*i*tstest(j))*bells(j)
cccc        cd=cd+exp(ima*pi*i*tstest(j))*ders(j)
 3200 continue
  
c 
        project(i+1)=cd*h
 3400 continue
c 
        call prin2('and project=*',project,nfour*2)
c 
c       plot both the bell and the coefficients of the Fourier series
c 
c 
        iw=22
        call lotagraph(iw,tstest,bells,ntest,
     1      'and the bell *')
c 
  
        do 3600 i=1,ntest/2
c 
        rnums(i)=i*2-1
c 
        coefs(i)=project(i*2)
  
 3600 continue
  
        nn=20
        iw=23
        call lotagraph(iw,rnums,coefs,nn,
     1      'and Fourier coefficients of the bell *')
  
c 
        call prin2('and fourier coefficients are=*',coefs,nfour)
  
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine prolbell(coesbell,coesder,nterms,x,
     1      bell,der)
        implicit real *8 (a-h,o-z)
        save
        dimension coesbell(1),coesder(1)
c 
c       This subroutine evaluates a prolate bell and its
c       derivative at the user-supplied point x. Please note that
c       prior to calling this subroutine, the user is expected to
c       call the initialization subroutine prolbelin (see), in order
c       to generate the arrays coesbell, coesder,and their lengths
c       nterms. Also, please note that if the point x is located
c       outside the interval [-1,1], both the bell and its derivative
c       are set to zero
c 
c                       Input parameters:
c 
c  coesbell - the Legendre expansion on the interval [-1,1] of
c       the left half of the filter
c  coesder - the Legendre expansion on the interval [-1,1] of
c       the derivetive of the left half of the filter
c  nterms - the lengths of the arrays coesbell, coesder
c  x - the point at which the bell and its derivative are to be
c       evaluated
c 
c                       Output parameters:
c 
c  bell - the value of the bell at the point x
c  der - the value of the derivative of the bell at the point x
c 
c       . . . if the user-supplied point x is outside the interval [-1,1]
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
        call legeexev(t,bell,coesbell,nterms)
        call legeexev(t,der,coesder,nterms)
c 
        return
 2200 continue
c 
c 
c       calculate the bell and its derivative if x < 0
c 
        if(x .le. 0) goto 2400
c 
        t=-x*2+1
c 
        call legeexev(t,bell,coesbell,nterms)
        call legeexev(t,der,coesder,nterms)
c 
        return
 2400 continue
c 
c 
  
        return
        end
  
  
c 
c 
c 
c 
c 
        subroutine prolbelin(ier,eps,coesbell,coesder,nterms,c,
     1      w,lenw,lused)
        implicit real *8 (a-h,o-z)
        save
        dimension coesbell(1),coesder(1),w(1)
c 
c        This subroutine constructs the Legendre expansions of
c        the Prolate bell and its derivative corresponding
c        to the user-specified precision eps.
c 
c    Explanation:
c 
c        The bell is defined on the interval [-1,1]; it is equal to zero
c        at x=-1, and to 1 at x=1. The bell is equal to the indefinite
c        integral of \psi_^c, with positive real c chosen in such a
c        manner that
c 
c        (\psi_0^c)' (1) < \sim eps                                    (1)
c 
c        (which seems to be as good a definition of "precision eps"
c        as any.
c 
c                      Input parameters:
c 
c  eps - the precision to which the bell is a bell, in the
c        sense defined by (1) above
c  lenw - the length of the user-supplied array w (in real *8 words);
c       please note that lenw = 7*(3/2*c+300) +100 is always more than
c       sufficient (in fact, grossly excessive).
c 
c                      Output parameters:
c 
c  ier - error return code
c  lused - the maximum number of elements of the array w used by
c       the subroutine at any time
c 
c       determine the parameter c needed to obtain the
c       user-specified precision epsilon
c 
  
        epslolog=log(-log(eps))
  
  
        adjust= -.94818E-01 *epslolog+0.63178E+00
c 
        alpha=  0.32856E-01
        beta=  0.69242E+00
        gamma=  0.78742E+00
  
        ccclog=alpha*epslolog**2+beta*epslolog+gamma
        c=exp(ccclog)
c 
c       construct the coefficients of the Legendre expansion
c       of the function \psi_0^c on the interval [-1,1]
c 
        call prolps0i(jer,c,w,lenw,nterms,lused,rkhi)
c 
        if(jer .eq. 0) goto 1150
c 
        ier=4
        return
 1150 continue
c 
        do 1200 i=1,nterms
        coesder(i)=w(i)
 1200 continue
c 
c        construct Legendre coefficients of the indefinite
c        integral of \psi_0^c on the interval [-1,1]
c 
        call legeinte(w,nterms-2,coesbell)
c 
c       scale the coefficients of both the bell and its derivative
c       so that the bell is normalized
c 
        do 1400 i=1,nterms
c 
        coesbell(i)=coesbell(i)/w(1)/2
        coesder(i)=coesder(i)/w(1)/2
 1400 continue
c 
        nterms=nterms-1
c 
        return
        end
  
