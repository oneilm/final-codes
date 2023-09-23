        implicit real *8 (a-h,o-z)
        dimension ab(2,10 000),xs(10 000),whts(10 000),
     1      expons(100 000),xstest(10000)
        complex *16 coesexp(100 000),
     1      fs(100 000),coesexp2(100 000),diffs(10 000),
     3      fstest(10 000),fstest2(10 000),fxx,fxxp,fxxm,
     4      derout2,fout,derout,errs(10 000)
c 
         external fun1,fun2,fun3,fun4
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER k '
        READ *,k
        CALL PRINF('k=*',k,1 )
C 
c        construct the subdivision of the interval [a,b]
c 
  
  
        eps=1.0d-12
cccc        eps=1.0/4
  
        a=100
  
        b=1000000 * 100.0d0  *100
cccc        b=200
        lab=10 000
  
        iferrs=1
  
        time1=clotatim()
  
        do 2200 i=1,100
  
  
         call exponexpand(ier,a,b,fun4,par1,par2,eps,
     1    ab,lab,nn,k,kk,coesexp,xs,whts,fs,
     2       coesexp2,expons,errs,iferrs)
  
  
cccc        call prinf('kk=*',kk,1)
  
  
 2200 continue
  
        time2=clotatim()
  
        call prin2('after expoexpand, time2-time1=*',time2-time1,1)
c 
         call prinf('after exponexpand, kk=*',kk,1)
         call prinf('after exponexpand, nn=*',nn,1)
c 
c        evaluate the obtained expansions at all points
c 
        call prin2('coesexp2=*',coesexp2,nn*kk*2)
        call prin2('expons=*',expons,nn*kk)
        call prin2('errs=*',errs,nn*k*2)
  
c 
c        test the "large interval" expansion evaluator
c 
        call exponeval(ier,xs(50),ab,nn,kk,expons,coesexp2,
     1      fout,derout)
c 
        call prin2('after exponeval, fout=*',fout,2)
        call prin2('and fs(50)=*',fs(50),2)
  
        call prin2('and fout-fs(50)=*',fout-fs(50),2)
  
c 
c       conduct mass testing of the evaluation
c 
        ntest=500
        h=(b-a)/ntest
        dd=0
        ddd=0
        do 4200 i=1,ntest
c 
        xstest(i)=(i-1)*h+a+h/2
c 
        call exponeval(ier,xstest(i),ab,nn,kk,expons,coesexp2,
     1      fstest(i),derout)
  
        call fun4(xstest(i),par1,par2,fstest2(i))
  
        diffs(i)=fstest(i)-fstest2(i)
c 
        dd=dd+abs(fstest(i))**2
        ddd=ddd+abs(diffs(i))**2
 4200 continue
c 
cccc        call prin2('xstest as created*',xstest,ntest)
cccc        call prin2('fstest as created*',fstest,ntest*2)
cccc        call prin2('fstest2 as created*',fstest2,ntest*2)
cccc        call prin2('diffs=*',diffs,ntest*2)
  
        errel=sqrt(ddd/dd)
  
        call prin2('and errel=*',errel,1)
  
cccc        stop
  
c 
c        test the edrivative
c 
        xx=0.4
        h=0.0001/100
        call exponeval(ier,xx,ab,nn,kk,expons,coesexp2,
     1      fxx,derout)
  
        call prin2('testing derivative, derout=*',derout,2)
  
        call exponeval(ier,xx+h,ab,nn,kk,expons,coesexp2,
     1      fxxp,derout)
  
        call exponeval(ier,xx-h,ab,nn,kk,expons,coesexp2,
     1      fxxm,derout)
  
        derout2=(fxxp-fxxm)/h/2
        call prin2('and derout2=*',derout2,2)
  
        d=0
        do 5200 i=1,nn*kk
c 
        d=d+whts(i)
 5200 continue
c 
        call prin2('d=*',d,1)
        call prin2('1-d=*',1-d,1)
  
  
  
  
        stop
        end
  
c 
c 
c 
c 
c 
        subroutine fun2(xx,par1,par2,f)
        implicit real *8 (a-h,o-z)
        save
        dimension par1(1),par2(1)
        real *8 xx
c 
        complex *16 f,ima
c 
        data ima/(0.0d0,1.0d0)/
  
  
  
        x=xx
cccc        f=1/( sqrt(sqrt(x)) )**3 +exp(ima*1500*x)
        f=1/( sqrt(sqrt(x)) )**3+ sin(5*x)*ima
cccc        f=x**5+cos(100*x)
cccc        f=x**5+cos(50*x)
  
  
  
cccc        f=x**5+cos(10*x)  + ima*x**2/100
  
        f=x**5+cos(10*x)
  
cccc        f=cos(2*x) *ima
  
  
cccc            call prin2('f=*',f,2)
  
  
cccc           f=x
  
  
        f=( sqrt(sqrt(x)) )+ sin(5*x)*ima
cccc        f= sin(5*x)*ima
  
  
        f=sqrt(x)
cccc        call prin2('x=*',x,1)
  
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine fun4(x,par1,par2,f)
        implicit real *8 (a-h,o-z)
        save
        dimension par1(1),par2(1)
c 
        complex *16 f,ima,h1,z,rk
c 
        data ima/(0.0d0,1.0d0)/
  
  
        b=1000
  
  
        rk=1+ima*1.0d-13
  
  
  
        z=sqrt(b**2+x**2)
  
        z=z*rk
  
  
        ifexpon=0
  
  
  
        call hank103(z,f,h1,ifexpon)
  
  
  
  
  
        dd=b**2/ (sqrt(b**2+x**2)+x)
  
        f=f*exp(ima*rk*dd)
  
  
  
  
  
c 
        return
        end
  
  
c 
c 
c 
c 
c 
        subroutine fun3(x,par1,par2,f)
        implicit real *8 (a-h,o-z)
        save
        dimension par1(1),par2(1)
c 
        complex *16 f,ima,h1,z,rk
c 
        data ima/(0.0d0,1.0d0)/
  
  
        b=1000
  
cccc        b=0.1d0
  
        rk=1+ima*1.0d-10
  
  
  
        z=sqrt(b**2+(x*1.0d8)**2)
  
        z=z*rk
  
  
        ifexpon=1
  
  
  
        call hank103(z,f,h1,ifexpon)
  
        f=f*exp(-ima*rk*x*1.0d8)
  
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine fun1(x,par1,par2,f)
        implicit real *8 (a-h,o-z)
        save
        dimension par1(1),par2(1)
c 
        complex *16 f,ima,h1,z,rk
c 
        data ima/(0.0d0,1.0d0)/
  
cccc        f=1/( sqrt(sqrt(x)) )**3 +exp(ima*1500*x)
        f=1/( sqrt(sqrt(x)) )**3+ sin(5*x)*ima
cccc        f=x**5+cos(100*x)
cccc        f=x**5+cos(50*x)
  
  
  
cccc        f=x**5+cos(10*x)  + ima*x**2/100
  
cccc        f=x**5+cos(10*x)
  
  
  
  
cccc        f=x
  
  
        rk=1+ima*1.0d-10
  
        ifexpon=1
        z=x*100000000.0d0 * rk
        call hank103(z,f,h1,ifexpon)
  
        f=f*exp(-ima*z)
  
  
  
  
  
  
  
c 
        return
        end
  
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c         This is the end of the debugging code and the beginning of
c         the actual exponential expansion routines
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c        This file contains five user-callable subroutines:
c        exponexpand, exponeval, expoexpand, exporetr, exponeval
c        (actually, exporetr is an entry in the subroutine
c        expoexpand). Following is a brief discussion of these
c        subroutines.
c 
c   exponexpand - constructs a nested exponential expansion
c        of the user-supplied complex-valued function fun; the
c        latter is specified via the subroutine fun. The subroutine
c        subdivides the user-supplied interval [a,b] into a
c        collection of subintervals, and constructs an exponential
c        expansion on each subinterval.
c 
c   exponeval - evaluates at a user-specified point x a complex-
c        valued function given by a nested exponential expansion;
c        the derivative of the said function is also returned. The
c        subroutine uses the results of a prior call to the
c        subroutine exponexpand (see).
c 
c   expoexpand - this subroutine constructs the coefficients
c        coesexp of an exponential expansion of a user-specified
c        (complex-valued) function. The function is supplied to
c        the subroutine via collection fs of its values tabulated
c        at 40 nodes tt on the interval [-1,1]; the user can
c        obtain the latter by calling the entry exporetr of this
c        subroutine. The expansion has the form
c 
c          f(x)=\sum_{j=1}^{72} coefsexp_j * exp(ima*expon_j),        (1)
c 
c        with the coefficients expon_j also obtained by a call to
c        the entry exporetr of this subroutine; all of the
c        coefficients {expon_j} are between -40 and 40.
c 
c   exporetr - returns to the user the parameters {expon_j},
c        j=1,2,...,72 in the expansion (1) above, the 40 nodes tt
c        to be used in the construction of the expansion (1), and
c        the weights ww corresponding to the nodes tt (these are
c        supplied for no apparent reason).
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
  
         subroutine exponexpand(ier,a,b,fun,par1,par2,
     1       eps,ab,lab,nn,k,kk,coesexp,xs,whts,fs,
     2       coesexp2,expons,errs,iferrs)
         implicit real *8 (a-h,o-z)
        save
         dimension ab(2,1),xs(1),whts(1),par1(1),par2(1),
     1       expons(1),expons72(75)
         complex *16 coesexp(1),fs(1),coesexp2(1),
     1       ima,cd
c 
         external fun
c 
         data ima/(0.0d0,1.0d0)/
c 
c        This subroutine constructs a nested exponential expansion
c        of the user-supplied complex-valued function fun; the
c        latter is specified via the subroutine fun. The subroutine
c        subdivides the user-supplied interval [a,b] into a
c        collection of subintervals; on the m-th subinterval,
c        the expansion has the form
c 
c            f(x)=\sum_{j=1}^k coesexp2(j,m) \cdot
c                  e^{i \dot expons(j,m) \cdot x},                       (1)
c 
c        and can be evaluated by calling the subroutine exponeval
c        (see); the complex coefficients coesexp2 and the exponents
c        expons are returned by the subroutine, as well certain
c        other types of information (see below).
c 
c 
c                      Input parameters:
c 
c  a,b - the ends of the interval on which the function is to be
c        approximated
c  fun - the function to be approximated; given by a subroutine
c        with the calling sequence
c 
c        fun(x,par1,par2,f)
c 
c  par1,par2 - whatever the user-supplied subroutine needs
c  lab - the length of the user-provided array ab (see below)
c  eps - accuracy to which the subroutine will attempt to approximate
c        the function fun
c  iferrs - integer parameter telling the subroutine whether the
c        errors of the approximation at the support nodes (nn*kk of
c        them things) are to be evaluated and returned to the user:
c    iferrs=1 tells the subroutine to evaluate the errors
c    iferrs=0 tells the subroutine to skip the error evaluation
c 
c                     Output parameters:
c 
c 
c  ier - error return code.
c         ier=0 means successful execution
c         ier=1024 means that the space allocated in the array ab is
c         insufficient
c  ab - the array of ends of subintervals created. For each i=1,2,...,nn,
c     ab(1,i) is the left end of the i-th subinterval; ab(2,i) is
c     the right end of the i-th subinterval
c  nn - the number of subintervals created on the interval [a,b];
c         also the number of entries in the array ab (see above)
c  k - the number of nodes at which the function has been evaluated
c         on each subinterval. ALWAYS RETURNED EQUAL TO 40
c  kk - the number of exponentials into which the function will
c         be decomposed on each of nn subintervals. ALWAYS RETURNED
c         EQUAL TO 72
c  coesexp - the coefficients of all nn expansions, each expansion
c         containing 72 coefficients. THESE COEFFICIENTS ARE FOR
C         EXPONENTS DEFINED ON THE ASSUMPTION THAT EACH OF THE
C         SUBINTERVALS IS THE INTERVAL [-1,1]. IN MOST CASES, IT
C         IS MUCH SIMPLER TO USE THE ARRAY COESEXP2   !!!
C  xs - the nodes at which the function has been evaluated; obviously,
c         there are k*nn=40*nn of them things
c  whts - the weights corresponding to the nodes in array xs
c  fs - the values of the function fun at the nodes xs
c  coesexp - the coefficients of all nn expansions, each expansion
c         containing 72 coefficients. THESE ARE THE COEFFICIENTS
C         IN THE EXPANSION (1) ABOVE.
c  expons - the array of coefficients in (1)
c  errs - the array of approximation errors at all k*nn=40*nn points
c 
c 
         ier=0
C 
c        construct the subdivision of the interval [a,b]
c 
         call expo_callinte(jer,a,b,fun,par1,par2,
     1      ab,lab,nn,lused,eps)
c 
         if(jer .ne. 0) then
c 
             ier=1024
             return
         endif
c 
c        discretize all of the obtained subintervals
c 
        call expo_cnodewht(ab,nn,xs,whts,k,expons72)
c 
c       on each of the obtained subintervals, construct the exponential
c       expansion of the user-supplied function fun
c 
cccc        call prinf('before expandem, k=*',k,1)
c 
        call expandem(ab,nn,xs,whts,k,fun,par1,par2,coesexp,fs,
     1      errs,iferrs)
c 
c        rescale the exponential coefficients so that they
c        correspond to the exponentials on the interval [a,b]
c 
        kk=72
        do 3400 i=1,nn
c 
        alpha=2/(ab(2,i)-ab(1,i))
        beta=1-alpha*ab(2,i)
c 
        j0=(i-1)*kk
c 
        do 3200 j=1,kk
c 
        cd=ima*expons72(j)*beta
        coesexp2(j0+j)=coesexp(j0+j)*exp(cd)
c 
        expons(j0+j)=expons72(j)*alpha
 3200 continue
c 
 3400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine exponeval(ier,x,ab,nn,kk,expons,coefs,
     1      fout,derout)
        implicit real *8 (a-h,o-z)
        save
        dimension expons(kk,1),ab(2,1)
        complex *16 coefs(kk,1),fout,derout
c 
c        This subroutine evaluates at a user-specified point x a
c        complex-valued function given by a nested exponential
c        expansion; the derivative of the said function is also
c        returned. he subroutine uses the results of a prior call
c        to the subroutine exponexpand (see).
c 
c                      Input parameters:
c 
c  x - the point where the function is to be evaluated
c  ab,nn,kk,expons,coefs - parameters produced by a prior call
c        to the subroutine exponexpand.
c 
c                      Output parameters:
c 
c  fout - the value of the function at the point x
c  derout - the value of the derivative of the function at the
c        point x
c 
c        find the subinterval on which the point x lives
c 
        call expo_findinte(ier,x,ab,nn,intnum)
c 
c       evaluate the function at the point x
c 
        call expoeval(expons(1,intnum),coefs(1,intnum),
     1      kk,x,fout,derout)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine expoeval(expons,coefs,nn,x,fout,derout)
        implicit real *8 (a-h,o-z)
        save
        dimension expons(1)
        complex *16 coefs(1),cd,ima,fout,derout
c 
        data ima/(0.0d0,1.0d0)/
c 
        fout=0
        derout=0
c 
        do 1200 i=1,nn
c 
        cd=exp(x*expons(i)*ima)
        fout=fout+cd*coefs(i)
c 
        derout=derout+cd*coefs(i)*ima*expons(i)
  
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine expoexpand(fs,coefsexp,diffs,iferror)
        implicit real *8 (a-h,o-z)
        save
        complex *16 coefs(100),fs(1),diffs(1),
     1      cd,coefsexp(1),coefs3(72,37),
     2      cfs3(50),ima,evalcoefs(72,40)
        real *8 tt40(100),ww40(100),arrout(36,40),
     1      expons(73),expon(1),tt(1),ww(1),rcoefs3(72,37)
c 
        data ima/(0.0d0,1.0d0)/,ifcalled/0/
c 
c       this subroutine constructs the coefficients coesexp of
c       an exponential expansion of a user-specified (complex-
c       valued) function. The function is supplied to the
c       subroutine via collection fs of its values tabulated at
c       40 nodes tt on the interval [-1,1]; the the user can
c       obtain the latter by calling the entry exporetr of this
c       subroutine. The expansion has the form
c 
c          f(x)=\sum_{j=1}^{72} coefsexp_j * exp(ima*expon_j),        (1)
c 
c       with the coefficients expon_j also obtained by a call to
c       the entry exporetr of this subroutine; all of the
c       coefficients {expon_j} are between -40 and 40.
c 
c                   Input parameters:
c 
c  fs - 40 complex values of the function to be expanded, tabulated
c       at the nodes tt on the interval [-1,1]
c  iferror - tells the subroutine whether to return the errors of the
c       calculation. This is a CPU time-saving option, since the
c       calculation of the errors is about 3 times more expensive
c       than the calculation of the band-limited approximation
c 
c                   Output parameters:
c 
c  coefsexp - the 72 complex coefficients in (1)
c  diffs - the errors (40 of them) of the approximation at the nodes
c       tt (where the input points are tabulated)
c 
c 
c       . . . if this is the first call to this subroutine,
c             retrieve the coefficients of the exponential
c             expansions and the prolate expansion matrix
c 
        if(ifcalled .eq. 1) goto 1600
c 
        call expons_ret(expons,coefs3,tt40,ww40)
c 
        do 1160 i=1,37
c 
        ii=i/2
        jj=i-ii*2
        cd=-ima
        if(jj .ne. 0) cd=1
        do 1150 j=1,72
c 
        rcoefs3(j,i)=coefs3(j,i)*cd
 1150 continue
 1160 continue
c 
        call expomatr_ret(arrout)
c 
        do 1400 i=1,72
        do 1200 j=1,40
c 
        cd=exp(tt40(j)*expons(i)*ima)
c 
        evalcoefs(i,j)=cd
 1200 continue
 1400 continue
c 
        ifcalled=1
 1600 continue
c 
c        construct the coefficients of the prolate expansion
c 
        call expomatvec2(arrout,36,40,fs,coefs)
c 
c       convert the obtained prolate coefficients into
c       coefficients of the exponential expansion
c 
        call expoappl(coefs,coefsexp,rcoefs3)
c 
c       evaluate the obtained expansion at the nodes
c       tt40, and compate the obtained result with
c       the original function
c 
        if(iferror .eq. 0) return
c 
        do 2800  i=1,40
c 
        cd=0
        do 2700 j=1,72
c 
        cd=cd+evalcoefs(j,i)*coefsexp(j)
 2700 continue
c 
        cfs3(i)=cd
c 
        diffs(i)=cfs3(i)-fs(i)
 2800 continue
c 
        return
c 
c 
c 
c 
        entry exporetr(expon,tt,ww)
c 
c        this entry returns to the user the 72 exponentials to be
c        used by the subroutine ; also the 40 standard nodes on the
c        interval [-1,1], and the corresponding weights
c 
        call points_ret(expon,tt,ww)
        return
        end
c 
c 
c 
c 
c 
        subroutine expandem(ab,nn,xs,whts,k,fun,par1,par2,
     1       coesexp,fs,errs,iferrs)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),whts(1),ab(2,1),par1(1),par2(1)
        complex *16 fs(1),coesexp(72,1),diffs(100),errs(k,1)
c 
c        evaluate the user-supplied function at all nodes on
c        all subintervals
c 
        n=nn*k
        do 2200 i=1,n
c 
        call fun(xs(i),par1,par2,fs(i))
 2200 continue
c 
c        one subinterval after another, decompose the function
c        fun into exponentials
c 
        do 2400 i=1,nn
c 
        iff=(i-1)*k+1
c 
        call expoexpand(fs(iff),coesexp(1,i),diffs,iferrs)
c 
        if(iferrs .eq. 0) goto 2400
c 
        do 2300 j=1,k
c 
        errs(j,i)=diffs(j)
 2300 continue
c 
 2400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine expo_cnodewht(ab,nn,x,whts,k7,expons72)
        implicit real *8 (a-h,o-z)
        save
        dimension ab(2,1),x(1),whts(1),tt40(50),ww40(50),
     1      expons72(1)
c 
c 
        data ifcalled/0/
c 
c       construct the nodes and weights of the prolate
c       quadrature on the interval [-1,1]
c 
        k=40
        if(ifcalled .eq. 0) call points_ret(expons72,tt40,ww40)
        ifcalled=1
        k7=k
c 
c       one subinterval after another, construct nodes and weights
c       of the quadrature on the big interval
c 
        do 2000 i=1,nn
        u=(ab(2,i)-ab(1,i))/2
        v=(ab(2,i)+ab(1,i))/2
        ik=(i-1)*k+1
c 
        call expo_discr(ab(1,i),ab(2,i),tt40,ww40,k,
     1      x(ik),whts(ik) )
 2000 continue
c 
        return
        end
c 
c 
c 
c 
c 
         subroutine expo_callinte(ier,a,b,fun,par1,par2,
     1      ab,lab,nn,lused,eps)
        implicit real *8 (a-h,o-z)
        save
        dimension par1(1),par2(1),ab(2,1)
        equivalence (rea,intint)
c 
c       this subroutine finds a subdivision of the interval [a,b],
c       such that on every subinterval, the user-provided function
c       fun is integrated to precision eps by a k-point Gaussian
c       quadrature.
c 
c                          input parameters:
c 
c  a,b - the ends of the interval to be subdivided
c  fun - the user-supplied function to be integrated
c      The calling sequence of fun must be:
c 
c          fun(x,par1,par2,f),
c 
c       where x is the input argument where the function is to be
c       evaluated, f  (output) is the function, and par1,par2 are
c       both input parameters, carrying whatever information funget
c       might need.
c 
c  par1,par2 - whatever the user-supplied subroutine needs
c  lab - the length of the user-provided array ab (see below)
c 
c                     output parameters:
c 
c  ier - error return code.
c         ier=0 means successful execution
c         ier=4 means that the space allocated in the array ab is
c         insufficient
c  ab - the array of ends of subintervals created. For each i=1,2,...,nn,
c     ab(1,i) is the left end of the i-th subinterval; ab(2,i) is
c     the right end of the i-th subinterval
c  nn - the number of subintervals created on the interval [a,b];
c     also the number of entries in the array ab (see above)
c  lused - the total amount of space used by the subroutine in the
c     array ab; note that normally, it  is somewhat greater than
c     the expected nn*2; The difference is not large, but should
c     be taken into account when allocating memory.
c 
c       . . . start the recursive subdivision process
c 
        ab(1,1)=a
        ab(2,1)=b
        nn=0
        istore=1
        ier=0
c 
c       recursively subdivide the intervals till none need
c       to be subdivided
c 
        lused=0
        do 2000 i=1,1 000 000
c 
        nn=nn+1
c 
c       check if any more intervals need to be subdivided
c 
        if(nn .gt. istore) goto 2200
c 
c       check if this interval needs to be subdivided
c 
        d=ab(2,nn)-ab(1,nn)
c 
        call expo_ifsplit2(fun,par1,par2,ab(1,nn),ab(2,nn),
     1      eps,ier)
c 
        if(ier .eq. 0) goto 2000
c 
c       this interval needs to be subdivided. act accordingly
c 
        if(lused .lt. istore*2+4) lused=istore*2+4
        if(istore*2 +4 .lt. lab) goto 1800
        ier=4
        return
 1800 continue
c 
        istore=istore+1
        ab(1,istore)=ab(1,nn)
        ab(2,istore)=(ab(1,nn)+ab(2,nn))/2
c 
        istore=istore+1
        ab(1,istore)=ab(2,istore-1)
        ab(2,istore)=ab(2,nn)
c 
        intint=9654236
        ab(1,nn)=rea
        ab(2,nn)=1.0d35
 2000 continue
 2200 continue
c 
        nn=nn-1
c 
c       scan the intervals we have created, discarding those that have kids
c 
        ii=0
        do 2400 i=1,nn
        rea=ab(1,i)
c 
        if( (ab(2,i) .gt. 1.0d34 ) .and.
     1      (intint .eq. 9654236) ) goto 2400
        ii=ii+1
        ab(1,ii)=ab(1,i)
        ab(2,ii)=ab(2,i)
 2400 continue
        nn=ii
c 
c        bubble-sort the array ab
c 
        call expo_cbubble(ab,nn)
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine expo_discr(a,b,tt,ww,n,xx,wwxx)
        implicit real *8 (a-h,o-z)
        save
        dimension tt(1),ww(1),xx(1),wwxx(1)
c 
        u=(b-a)/2
        v=(b+a)/2
c 
        do 1200 i=1,n
c 
        xx(i)=u*tt(i)+v
        wwxx(i)=ww(i)*u
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine expo_findinte(ier,x,ab,nn,intnum)
        implicit real *8 (a-h,o-z)
        save
        dimension ab(2,nn)
c 
        data intold/-10/,ithresh/10/
c 
c       check if the point is on the subinterval as the preceding one
c 
        ier=0
cccc        if(intold .lt. 0) goto 2000
        if( (intold .lt. 0) .or. (intold .gt. nn) ) goto 2000
c 
        intnum=intold
        if( (x .ge. ab(1,intnum) ) .and. (x .le. ab(2,intnum) ) ) return
 2000 continue
c 
c      the point is not on the same subinterval as the preceding one.
c      if nn is less than ithresh, use direct scan to find the proper
c      interval
c 
       if(nn .gt. ithresh) goto 3000
c 
       if(x .lt. ab(1,1)) then
           ier=4
           return
       endif
c 
        do 2200 j=1,nn
c 
           intnum=j
c 
        if(ab(2,j) .ge. x) goto 2400
 2200 continue
c 
        ier=4
        return
 2400 continue
c 
        intold=intnum
        return
c 
 3000 continue
c 
c      The point is not on the same subinterval as the preceding one,
c      and nn is greater than ithresh; use bisection to find the proper
c      interval
c 
       i1=1
       i2=nn
       i3=(i1+i2)/2
c 
       do 3400 i=1,100
c 
       if(x .ge. ab(1,i3)) i1=i3
       if(x .le. ab(2,i3)) i2=i3
c 
       if(i2 .eq. i1) goto 3600
c 
       i3=(i1+i2)/2
 3400 continue
c 
 3600 continue
  
       if(x .lt. ab(1,i3)) i3=i3-1
       if(x .gt. ab(2,i3)) i3=i3+1
  
       intnum=i3
       intold=intnum
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine expo_cbubble(ab,nn)
        implicit real *8 (a-h,o-z)
        save
        dimension ab(2,1)
c 
c       sort the array ab with respect to the first coordinate
c 
        do 2000 i=1,nn
        do 1200 j=1,i-1
        if(ab(1,i) .ge. ab(1,j) ) goto 1200
        d=ab(1,i)
        ab(1,i)=ab(1,j)
        ab(1,j)=d
c 
        d=ab(2,i)
        ab(2,i)=ab(2,j)
        ab(2,j)=d
 1200 continue
 2000 continue
        return
        end
  
c 
c 
c 
c 
c 
        subroutine expo_ifsplit2(fun,par1,par2,a,b,eps,ifsplit)
        implicit real *8 (a-h,o-z)
        save
        dimension par1(1),par2(1),xx1(100),ww0(100),xx0(100)
        complex *16 fs0(100),fs1(100),fs2(100),diffs(100)
        dimension pnts(20),expons(100),tt40(50),
     1      ww40(50),ww1(100)
c 
        data ifcalled/0/
        data pnts/
     1  -.99312859918509492479D+00,-.96397192727791379127D+00,
     2  -.91223442825132590587D+00,-.83911697182221882339D+00,
     3  -.74633190646015079261D+00,-.63605368072651502545D+00,
     4  -.51086700195082709800D+00,-.37370608871541956067D+00,
     5  -.22778585114164507808D+00,-.76526521133497333755D-01,
     6  0.76526521133497333755D-01,0.22778585114164507808D+00,
     7  0.37370608871541956067D+00,0.51086700195082709800D+00,
     8  0.63605368072651502545D+00,0.74633190646015079261D+00,
     9  0.83911697182221882339D+00,0.91223442825132590587D+00,
     *  0.96397192727791379127D+00,0.99312859918509492479D+00/
c 
c        if needed, retrieve from the subroutine interp40to20 the
c        interpolation matrix to be used to test the smoothness of
c        the function
c 
        if(ifcalled .eq. 0) call exporetr(expons,tt40,ww40)
        ifcalled=1
c 
c        discretize the interval [a,b] using the 40 prolate nodes
c 
        call expo_discr(a,b,tt40,ww40,40,xx0,ww0)
c 
c        discretize the interval [a,b] using the 20 Gaussian nodes
c 
        call expo_discr(a,b,pnts,ww40,20,xx1,ww1)
c 
c        evaluate the function to be analyzed at both grids
c 
        do 1400 i=1,40
c 
        call fun(xx0(i),par1,par2,fs0(i))
 1400 continue
c 
        do 1600 i=1,20
c 
        call fun(xx1(i),par1,par2,fs1(i))
 1600 continue
c 
c        interpolate the function from the prolate grid to the
c        Gaussian grid
c 
        call expomatvec4(fs0,fs2)
c 
        d=0
        dd=0
        do 1800 i=1,20
c 
        diffs(i)=fs2(i)-fs1(i)
c 
        d=d+fs1(i)*conjg(fs1(i))
        dd=dd+diffs(i)*conjg(diffs(i))
 1800 continue
c 
        errel=sqrt(dd/d)
        ifsplit=0
        if(errel .gt. eps) ifsplit=1
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine expomatvec4(x,y)
        implicit real *8 (a-h,o-z)
        save
        dimension asymm(10,20),aasymm(10,20)
        complex *16 x(40),y(20),xsymm(100),y2(100),y1(100),
     1      xasymm(100),cd1,cd2
c 
        data ifcalled/0/
c 
c       symmetrize and antisymmetrize x
c 
        n=20
        m=40
        do 1200 i=1,20
c 
        xsymm(i)=(x(i)+x(41-i))/2
        xasymm(i)=(x(i)-x(41-i))/2
 1200 continue
c 
        if(ifcalled .eq. 0) call exposymm(asymm,aasymm)
        ifcalled=1
c 
c       apply the operator to the symmetrized part of x
c 
        do 2400 i=1,10
        cd1=0
        cd2=0
        do 2200 j=1,20
c 
        cd1=cd1+asymm(i,j)*xsymm(j)
        cd2=cd2+aasymm(i,j)*xasymm(j)
 2200 continue
c 
        y1(i)=cd1
        y2(i)=cd2
 2400 continue
c 
        do 4600 i=1,10
c 
        y(i)=y1(i)+y2(i)
        y(n-i+1)=y1(i)-y2(i)
 4600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine exposymm(symm,asymm)
        implicit real *8 (a-h,o-z)
        save
        dimension symm1(180),symm2(20),sbuff(2),
     1       asymm1(180),asymm2(20),abuff(2),
     2       symm(200),asymm(200)
c 
        equivalence (symm1(180),sbuff(1)),
     1      (symm2(1),sbuff(2)),(asymm1(180),abuff(1)),
     1      (asymm2(1),abuff(2))
  
        data symm1/
     1  0.2580509808887957D+00,0.3577470693252816D-01,
     2  0.2735139045508189D-02,-.7982852151286928D-02,
     3  -.7922513751392023D-02,-.4030904684909979D-02,
     4  -.2325330821794799D-03,0.2276488259606791D-02,
     5  0.3610306863122156D-02,0.4155632600919709D-02,
     6  0.8997492559705167D+00,-.1518001254251797D+00,
     7  -.1307058206205190D-01,0.2572460861874014D-01,
     8  0.2631107001100516D-01,0.1362220033353475D-01,
     9  0.1064842400028114D-02,-.7268200347919160D-02,
     *  -.1171077458195206D-01,-.1353077496145818D-01,
     1  -.1990402123928380D+00,0.6292108015140060D+00,
     2  0.4117676677324924D-01,-.4273948769057934D-01,
     3  -.4683889896169449D-01,-.2508099543785707D-01,
     4  -.2899122570997356D-02,0.1195857436422535D-01,
     5  0.1992421969344967D-01,0.2320034011492263D-01,
     6  0.1457508727373716D-01,0.6050997447741133D+00,
     7  -.1261790084043618D+00,0.4987628955240654D-01,
     8  0.6400247026959185D-01,0.3622665059868504D-01,
     9  0.6202752452214563D-02,-.1424120790956031D-01,
     *  -.2530557941666655D-01,-.2988540452498821D-01,
     1  0.9745417843457592D-01,-.1611105952486528D+00,
     2  0.7214973752153169D+00,-.3401982126078124D-01,
     3  -.7288972406217402D-01,-.4519070473336552D-01,
     4  -.1119864198786416D-01,0.1263388008591917D-01,
     5  0.2573264826284816D-01,0.3120939150520425D-01,
     6  -.1660283209280554D+00,0.5118306703063724D-01,
     7  0.5190800343122839D+00,-.4289144435351734D-01,
     8  0.6844047916966269D-01,0.5026244157201062D-01,
     9  0.1768776403745589D-01,-.6495174965025064D-02,
     *  -.2014205779810308D-01,-.2594069239227994D-01,
     1  0.1971270406482897D+00,0.1603814513071780D-02,
     2  -.2349698176201908D+00,0.7634516100353986D+00,
     3  -.4232163648130331D-01,-.4962957874227269D-01,
     4  -.2494014567647314D-01,-.3862638853342545D-02,
     5  0.8651803124011546D-02,0.1411942024165957D-01,
     6  -.1951472203690213D+00,-.2868451070518940D-01,
     7  0.1484512167573721D+00,0.4208980519758579D+00,
     8  -.3666447508992361D-01,0.4054757711037240D-01,
     9  0.3160299943108140D-01,0.1720182732147611D-01,
     *  0.7503538306386868D-02,0.3019574832104095D-02,
     1  0.1661282962341295D+00,0.3969296934856287D-01,
     2  -.9362374991163006D-01,-.2145458005727873D+00,
     3  0.8352411691268791D+00,-.1621290922048735D-01,
     4  -.3548917212766903D-01,-.3147727266553179D-01,
     5  -.2618289991538294D-01,-.2328667074648362D-01,
     6  -.1177873136199568D+00,-.3975837887017520D-01,
     7  0.4947182255986316D-01,0.1314845265543138D+00,
     8  0.3109165899166928D+00,-.5301478398637403D-01,
     9  0.3282138381636216D-01,0.4394404587795261D-01,
     *  0.4461691543423091D-01,0.4386439457926106D-01,
     1  0.5873338233034877D-01,0.3269246272265095D-01,
     2  -.1168748080441387D-01,-.7192385597338881D-01,
     3  -.1584087803730584D+00,0.8617535561248862D+00,
     4  -.1507390698956130D-01,-.5107136727086996D-01,
     5  -.5966728884488158D-01,-.6164411156927413D-01,
     6  0.2471082519865202D-02,-.2169956514645479D-01,
     7  -.2004140627208748D-01,0.2339559414997943D-01,
     8  0.9686893660500131D-01,0.2858996668318588D+00,
     9  -.5062892504568567D-01,0.4777436854039228D-01,
     *  0.6796644821119872D-01,0.7350718092452546D-01,
     1  -.5814661310333544D-01,0.9483238660469414D-02,
     2  0.4494722419404636D-01,0.1569259932893883D-01,
     3  -.5651444626132759D-01,-.1621684687763642D+00,
     4  0.8260852748766700D+00,-.2390255286993477D-01,
     5  -.6576337012125371D-01,-.7649332717559284D-01,
     6  0.1021289825729208D+00,0.1795795421731330D-02,
     7  -.6217646783187350D-01,-.4473600411735442D-01,
     8  0.2714659264376961D-01,0.1183598971444077D+00,
     9  0.3414170373002988D+00,-.5699936602105600D-01,
     *  0.4797550740352282D-01,0.6778275026259538D-01,
     1  -.1301700978063432D+00,-.1057400734914689D-01,
     2  0.7120513778589381D-01,0.6300353692377701D-01,
     3  -.6286915709444757D-02,-.9127714820573703D-01,
     4  -.2011865576976419D+00,0.7729255173248614D+00,
     5  -.4300948909165950D-02,-.4432942407741838D-01,
     6  0.1401175216245989D+00,0.1590898781035395D-01,
     7  -.7201556473440635D-01,-.7030922623155100D-01,
     8  -.6985518810447151D-02,0.7042107886331046D-01,
     9  0.1501142388519320D+00,0.4177099542843928D+00,
     *  -.1024955810348785D+00,0.1632447192304017D-02,
     1  -.1319000924766501D+00,-.1746263683740307D-01,
     2  0.6514831665397278D-01,0.6722940491229996D-01,
     3  0.1344287385426190D-01,-.5261573533129723D-01,
     4  -.1141404548465400D+00,-.2149671190087338D+00,
     5  0.7368906357835614D+00,0.7049412451262988D-01,
     6  0.1073533435425348D+00,0.1544541291916771D-01,
     7  -.5167771772271755D-01,-.5511951733550453D-01,
     8  -.1410633655180233D-01,0.3655594064370934D-01,
     9  0.8149761353745759D-01,0.1367394988413669D+00,
     *  0.4712261265257212D+00,-.2087682124847050D+00/
c 
        data symm2/
     1  -.6991885372025381D-01,-.1052782657799973D-01,
     2  0.3313678855681406D-01,0.3602231297963686D-01,
     3  0.1033027131881739D-01,-.2154828776474188D-01,
     4  -.4915225207797027D-01,-.7880633418781158D-01,
     5  -.1671287551104641D+00,0.7657700721766310D+00,
     6  0.2424957237614100D-01,0.3726644512908990D-02,
     7  -.1140802649058719D-01,-.1251052534459817D-01,
     8  -.3761206863114190D-02,0.7120507660631717D-02,
     9  0.1644780539908187D-01,0.2592707919959160D-01,
     *  0.4859910612469506D-01,0.3851232889894433D+00/
c 
        data asymm1/
     1  0.2461869775064312D+00,0.3319082015632594D-01,
     2  0.7179941862860256D-02,-.2417846355676279D-02,
     3  -.5299275005951385D-02,-.5032275139179321D-02,
     4  -.3732135070742633D-02,-.2389875142794903D-02,
     5  -.1290964328226477D-02,-.4044288387018381D-03,
     6  0.9387025933601706D+00,-.1433679092748919D+00,
     7  -.2771880350074941D-01,0.7461081173040961D-02,
     8  0.1777948184193434D-01,0.1702714072103716D-01,
     9  0.1266299656279136D-01,0.8120400453003160D-02,
     *  0.4390382844314759D-02,0.1376035929585597D-02,
     1  -.2669917717139593D+00,0.6146747243536112D+00,
     2  0.6691214531677488D-01,-.1090949749790273D-01,
     3  -.3223173531557257D-01,-.3141477979382933D-01,
     4  -.2349095548585101D-01,-.1510547053766874D-01,
     5  -.8180570230445489D-02,-.2566147325634103D-02,
     6  0.1048819146644494D+00,0.6240421421817753D+00,
     7  -.1607735336384730D+00,0.7645582905667775D-02,
     8  0.4519356949711489D-01,0.4549807396397893D-01,
     9  0.3433699335685958D-01,0.2217728119099177D-01,
     *  0.1204154705796853D-01,0.3782212322463146D-02,
     1  -.2421425583913771D-02,-.1813998387910827D+00,
     2  0.7604393318590011D+00,0.1254438415773607D-01,
     3  -.5316323957421262D-01,-.5694252783899279D-01,
     4  -.4363365989391616D-01,-.2837250163736807D-01,
     5  -.1546374954769501D-01,-.4866133165788830D-02,
     6  -.7259893321253936D-01,0.6912334696887475D-01,
     7  0.4815930026830864D+00,-.8619584929923401D-01,
     8  0.5170796551444335D-01,0.6360838536115853D-01,
     9  0.5004676766178230D-01,0.3288580138034292D-01,
     *  0.1802260756193673D-01,0.5686148669124206D-02,
     1  0.1262433625039162D+00,-.1042730307263558D-01,
     2  -.2049456829179698D+00,0.7958646896027384D+00,
     3  -.3228822779818337D-01,-.6331701764016337D-01,
     4  -.5244736449646421D-01,-.3506768208111594D-01,
     5  -.1937935655228358D-01,-.6137224398154217D-02,
     6  -.1601466361880900D+00,-.2530706673711001D-01,
     7  0.1310976352499634D+00,0.4056937951478451D+00,
     8  -.3723341346927654D-01,0.5305557086757972D-01,
     9  0.4979986954070191D-01,0.3439778286480216D-01,
     *  0.1927303106373762D-01,0.6139164538720217D-02,
     1  0.1754552812321647D+00,0.4641528178341560D-01,
     2  -.9256259805406531D-01,-.2205353817558143D+00,
     3  0.8250376837582272D+00,-.2588399791644544D-01,
     4  -.4083980307942125D-01,-.3040687075298767D-01,
     5  -.1749457513126867D-01,-.5630040536260704D-02,
     6  -.1737983573778172D+00,-.5649432330941651D-01,
     7  0.6625339420479326D-01,0.1597414191614889D+00,
     8  0.3314979259741657D+00,-.4781742625223999D-01,
     9  0.2315045593240875D-01,0.2250511761312894D-01,
     *  0.1383563809241951D-01,0.4553690485268342D-02,
     1  0.1574660304651850D+00,0.5783987081308820D-01,
     2  -.4562194618428838D-01,-.1204862497465295D+00,
     3  -.1873002096733972D+00,0.8624759065374686D+00,
     4  0.1018406451493059D-01,-.9592436957380454D-02,
     5  -.7995061556894469D-02,-.2838869365016828D-02,
     6  -.1293219211842793D+00,-.5237067018277216D-01,
     7  0.2823587589916924D-01,0.8754740871410074D-01,
     8  0.1305954631475690D+00,0.2781720978889769D+00,
     9  -.9036174034764139D-01,-.1099227531734969D-01,
     *  -.6039901551268566D-03,0.3638057107395178D-03,
     1  0.9262299253232555D-01,0.4189320019391354D-01,
     2  -.1309417867773709D-01,-.5717466806119242D-01,
     3  -.9063612492071654D-01,-.1468665922118074D+00,
     4  0.8775912388971752D+00,0.4721573783887807D-01,
     5  0.1322628245777167D-01,0.3109722583988274D-02,
     6  -.5081676847198828D-01,-.2815437733913137D-01,
     7  -.2033462355677994D-03,0.2862188436579152D-01,
     8  0.5680899834604864D-01,0.9556211791132045D-01,
     9  0.2822194075552656D+00,-.1337594128086685D+00,
     *  -.3294836174536327D-01,-.8074861957693399D-02,
     1  0.7346433213811150D-02,0.1282061803242843D-01,
     2  0.1177723288599354D-01,-.2183163171032802D-02,
     3  -.2681364512531532D-01,-.6180684184724647D-01,
     4  -.1394645415983337D+00,0.8473568928924960D+00,
     5  0.6863133736341753D-01,0.1561428514152476D-01,
     6  0.3452303599440878D-01,0.2566795456944246D-02,
     7  -.2161387583406542D-01,-.2146512392335917D-01,
     8  0.4832197501921607D-03,0.3590450797035167D-01,
     9  0.9177882677085276D-01,0.3542866014799127D+00,
     *  -.1566207602825290D+00,-.2835859584787611D-01,
     1  -.7185282267869619D-01,-.1664026842279134D-01,
     2  0.2963989023834865D-01,0.4154662418376158D-01,
     3  0.2168895782336859D-01,-.1549724442981292D-01,
     4  -.6548978217288967D-01,-.1714129234043166D+00,
     5  0.7609198334363377D+00,0.5396484052934684D-01,
     6  0.1021513218852031D+00,0.2825021657810611D-01,
     7  -.3576367542343244D-01,-.5732950195269509D-01,
     8  -.3904923116557260D-01,0.7012174107203829D-04,
     9  0.4887719451615142D-01,0.1216597971242549D+00,
     *  0.4889582515459839D+00,-.1239176548761912D+00/
c 
        data asymm2/
     1  -.1234689506289506D+00,-.3650182897739075D-01,
     2  0.3990005903661213D-01,0.6820939078571054D-01,
     3  0.5099655683669333D-01,0.1034809293805222D-01,
     4  -.3866295140606995D-01,-.1002830921560281D+00,
     5  -.2392643731435241D+00,0.5761253074703486D+00,
     6  0.1344678909596300D+00,0.4078299747321511D-01,
     7  -.4198520151521894D-01,-.7375943481137699D-01,
     8  -.5708697824537674D-01,-.1560921550790164D-01,
     9  0.3373678077679900D-01,0.9155097755884822D-01,
     *  0.1910713725773951D+00,0.7724969545978549D+00/
c 
        do 1200 i=1,200
c 
        symm(i)=symm1(i)
        asymm(i)=asymm1(i)
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine expoappl(coefs,coefsexp,rcoefs3)
        implicit real *8 (a-h,o-z)
        save
        complex *16 coefs(1),cd,coefsexp(1),ima,
     1      coefsexp2(100),coefsexp3(100),cd1
c 
        real *8 rcoefs3(72,37)
c 
        data ima/(0.0d0,1.0d0)/
c 
c       convert the obtained prolate coefficients into
c       coefficients of the exponential expansion
c 
        do 3600 i=1,36
c 
        cd=0
        cd1=0
        do 3400 j=1,36,2
c 
        cd=cd+rcoefs3(i,j)*coefs(j)
c 
        j1=j+1
        cd1=cd1+rcoefs3(i,j1)*coefs(j1)
 3400 continue
c 
        coefsexp2(i)=cd
        coefsexp3(i)=cd1
c 
 3600 continue
c 
        do 4700 i=1,36
c 
        coefsexp3(73-i)=-coefsexp3(i)
        coefsexp2(73-i)=coefsexp2(i)
 4700 continue
c 
        do 4800 i=1,72
c 
        coefsexp(i)=ima*coefsexp3(i)+coefsexp2(i)
 4800 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine expomatvec2(a,n,m,x,y)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m)
        complex *16 x(m),y(n),d,symmetr(100),asymmetr(100),d1
c 
        do 1200 i=1,m/2
c 
        symmetr(i)=x(i)+x(m-i+1)
        asymmetr(i)=x(i)-x(m-i+1)
 1200 continue
c 
        do 3400 i=1,n,2
        d=0
        d1=0
        do 3200 j=1,m/2
c 
        d=d+a(i,j)*symmetr(j)
 3200 continue
c 
        y(i)=d
c 
        if(i .eq. n) goto 3600
        i1=i+1
        do 3300 j=1,m/2
c 
        d1=d1+a(i1,j)*asymmetr(j)
 3300 continue
c 
        y(i1)=d1
 3400 continue
 3600 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine expons_cpy(cs,coefs,cd)
        implicit real *8 (a-h,o-z)
        save
        dimension cs(36)
        complex *16 coefs(1),cd
c 
        n=36
        do 1200 i=1,n
c 
        coefs(i)=cd*cs(i)
c 
        coefs(2*n-i+1)=conjg(coefs(i))
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine expomatr_ret(arrout)
        implicit real *8 (a-h,o-z)
        save
        dimension arr1(180),arr2(180),arr3(180),arr4(180),
     1      arr(36,20),arrout(36,40)
c 
        equivalence (arr1(1),arr(1,1)),(arr2(1),arr(1,6)),
     1      (arr3(1),arr(1,11)),(arr4(1),arr(1,16))
c 
        data arr1/
     1  0.3796510055179544D-09,-.4489342637615922D-08,
     2  0.3674108858643277D-07,-.2398829945631423D-06,
     3  0.1322444393439877D-05,-.6340904627475010D-05,
     4  0.2689708614247764D-04,-.1018934440895213D-03,
     5  0.3459006958733788D-03,-.1047673159539303D-02,
     6  0.2774393140782791D-02,-.6112560909540013D-02,
     7  0.1059518814115923D-01,-.1456882549192189D-01,
     8  0.1694397221936550D-01,-.1799724298223695D-01,
     9  0.1841947897103129D-01,-.1855627192215197D-01,
     *  0.1852036703847816D-01,-.1835201779274609D-01,
     1  0.1807102597084035D-01,-.1768995416712813D-01,
     2  0.1721801843684435D-01,-.1666262731219547D-01,
     3  0.1603014105374354D-01,-.1532629416376204D-01,
     4  0.1455644673493097D-01,-.1372574004358692D-01,
     5  0.1283919497774648D-01,-.1190177460801978D-01,
     6  0.1091842335039566D-01,-.9894090331081892D-02,
     7  0.8833741785966724D-02,-.7742365666610082D-02,
     8  0.6624970597285332D-02,-.5486580673106499D-02,
     9  0.2658959770175387D-08,-.2887008931537122D-07,
     *  0.2165766676379078D-06,-.1293865298096146D-05,
     1  0.6515143607803900D-05,-.2848347071618180D-04,
     2  0.1099826917531801D-03,-.3787261719430924D-03,
     3  0.1167543938874723D-02,-.3211193968988061D-02,
     4  0.7735174067855262D-02,-.1556993371293152D-01,
     5  0.2477143878004770D-01,-.3115256186199254D-01,
     6  0.3256411728300928D-01,-.3023360933552281D-01,
     7  0.2610367444193899D-01,-.2115138054005957D-01,
     8  0.1578724852105240D-01,-.1024669273558446D-01,
     9  0.4707463707086039D-02,0.6810972771880658D-03,
     *  -.5789580736801575D-02,0.1050563087708354D-01,
     1  -.1473304016637192D-01,0.1839164315557690D-01,
     2  -.2141750657706747D-01,0.2376316109569012D-01,
     3  -.2539773348096567D-01,0.2630689749627271D-01,
     4  -.2649259674342963D-01,0.2597251533994158D-01,
     5  -.2477928828399948D-01,0.2295945543499140D-01,
     6  -.2057217260584778D-01,0.1768770105595310D-01,
     7  0.1582958228793282D-07,-.1574184948310640D-06,
     8  0.1078505903161532D-05,-.5865163444661474D-05,
     9  0.2678298633953245D-04,-.1057250988782757D-03,
     *  0.3667397053134318D-03,-.1127841167422737D-02,
     1  0.3084217273697037D-02,-.7468826183852533D-02,
     2  0.1572653731614736D-01,-.2752044659590968D-01,
     3  0.3779930547159007D-01,-.3996750043342569D-01,
     4  0.3261525187359222D-01,-.1983450833171404D-01,
     5  0.5875223429899551D-02,0.7109216175507195D-02,
     6  -.1807129464761980D-01,0.2642446007051132D-01,
     7  -.3186356358348987D-01,0.3430761282757707D-01,
     8  -.3386656210619405D-01,0.3081185173527758D-01,
     9  -.2554596789317805D-01,0.1857004324633509D-01,
     *  -.1044987717094961D-01,0.1781374160312755D-02,
     1  0.6843316018527925D-02,-.1486758234187469D-01,
     2  0.2179872927724099D-01,-.2723175715828734D-01,
     3  0.3086757577638457D-01,-.3252495585333008D-01,
     4  0.3214583545004752D-01,-.2979392922456717D-01,
     5  0.8572652453400550D-07,-.7809786584263233D-06,
     6  0.4887335531037430D-05,-.2419459171631098D-04,
     7  0.1001707188762361D-03,-.3567830124413700D-03,
     8  0.1110062683798275D-02,-.3039142423926945D-02,
     9  0.7327710465879069D-02,-.1544875777449472D-01,
     *  0.2785564253446697D-01,-.4085192851795568D-01,
     1  0.4533377963283673D-01,-.3487472867338052D-01,
     2  0.1300974009521368D-01,0.1047957242616746D-01,
     3  -.2853474958121006D-01,0.3859210363495981D-01,
     4  -.4042994371844119D-01,0.3506252093645826D-01,
     5  -.2425949909447994D-01,0.1019979518016608D-01,
     6  0.4830511884749106D-02,-.1869947514101199D-01,
     7  0.2964032880842814D-01,-.3640261288096374D-01,
     8  0.3833708924498554D-01,-.3541311243501086D-01,
     9  0.2817252433618403D-01,-.1762924218677110D-01,
     *  0.5127927540888469D-02,0.7821956365292568D-02,
     1  -.1971902726438236D-01,0.2923241913232593D-01,
     2  -.3533824599404936D-01,0.3741779535811018D-01,
     3  0.4266463797509277D-06,-.3556810839655804D-05,
     4  0.2030495375018116D-04,-.9136195676690436D-04,
     5  0.3423025983010127D-03,-.1097455676414509D-02,
     6  0.3053152748947737D-02,-.7409771741328076D-02,
     7  0.1565115186068663D-01,-.2841897602151205D-01,
     8  0.4299256270093329D-01,-.5056204661088708D-01,
     9  0.4036797599917526D-01,-.1255490631719295D-01,
     *  -.1940605334223831D-01,0.4014513785785954D-01,
     1  -.4393021646099590D-01,0.3306303607173676D-01,
     2  -.1316345798740762D-01,-.9305862862075582D-02,
     3  0.2852602007787070D-01,-.4030250263665865D-01,
     4  0.4257550023934506D-01,-.3548380094771323D-01,
     5  0.2104817432960468D-01,-.2578044654019467D-02,
     6  -.1607042073278012D-01,0.3125018711285889D-01,
     7  -.4017062658631888D-01,0.4134182791983086D-01,
     8  -.3476112119104525D-01,0.2183158148466883D-01,
     9  -.5045149673756690D-02,-.1250201680936372D-01,
     *  0.2767200525909038D-01,-.3782838565566632D-01/
c 
        data arr2/
     1  0.1952011240425843D-05,-.1486317758203030D-04,
     2  0.7722454476763477D-04,-.3149236073392990D-03,
     3  0.1064000533238828D-02,-.3056945068785071D-02,
     4  0.7560021032237549D-02,-.1613458192101655D-01,
     5  0.2950954363398796D-01,-.4528906952641869D-01,
     6  0.5546005765370116D-01,-.4779662688848887D-01,
     7  0.1781179279543278D-01,0.2025797772854646D-01,
     8  -.4367937704380723D-01,0.4048558533658485D-01,
     9  -.1683236124154897D-01,-.1330533783729819D-01,
     *  0.3696632018353829D-01,-.4611328076279281D-01,
     1  0.3888324261485107D-01,-.1891874914390619D-01,
     2  -.6609720619224401D-02,0.2955629475388895D-01,
     3  -.4321169604300775D-01,0.4401183790383351D-01,
     4  -.3222351619697626D-01,0.1156522372139805D-01,
     5  0.1202656131318506D-01,-.3209800433647443D-01,
     6  0.4338601106696698D-01,-.4311367456002362D-01,
     7  0.3158508025545154D-01,-.1197940378328220D-01,
     8  -.1055847877633008D-01,0.3027227209902396D-01,
     9  0.8185361610219380D-05,-.5677262101666198D-04,
     *  0.2675831209429892D-03,-.9849634278297064D-03,
     1  0.2985335565986729D-02,-.7634168377594799D-02,
     2  0.1662934203201352D-01,-.3080152754494838D-01,
     3  0.4779483484458330D-01,-.5980309399842045D-01,
     4  0.5466293247655700D-01,-.2501457136094244D-01,
     5  -.1684767500073817D-01,0.4314992241264827D-01,
     6  -.3565116180290491D-01,0.2391909610327805D-02,
     7  0.3197515479729560D-01,-.4765562923422252D-01,
     8  0.3846441854598900D-01,-.1097815615362218D-01,
     9  -.2093757468691530D-01,0.4309145698756375D-01,
     *  -.4670556656821714D-01,0.3118722236158591D-01,
     1  -.3487711903149574D-02,-.2521070231586074D-01,
     2  0.4397702214097388D-01,-.4610045823126692D-01,
     3  0.3123707949307739D-01,-.5169206855437622D-02,
     4  -.2254374278937073D-01,0.4207264899082132D-01,
     5  -.4672732532263060D-01,0.3513169481882301D-01,
     6  -.1150903131855458D-01,-.1592212149721992D-01,
     7  0.3133600482924273D-04,-.1972792444679353D-03,
     8  0.8397862563126347D-03,-.2774718187477264D-02,
     9  0.7489749041180374D-02,-.1688027418638606D-01,
     *  0.3193542392430563D-01,-.5024136179383464D-01,
     1  0.6371614352688861D-01,-.6000004475940435D-01,
     2  0.3097444698309431D-01,0.1336830153651855D-01,
     3  -.4257569402578658D-01,0.3337053758118367D-01,
     4  0.5832229007608353D-02,-.4103701301007031D-01,
     5  0.4547084760041353D-01,-.1882916194940900D-01,
     6  -.1944800000142068D-01,0.4566543242590305D-01,
     7  -.4556939228395854D-01,0.2052183333949101D-01,
     8  0.1511867152719688D-01,-.4258624044340990D-01,
     9  0.4834859801322455D-01,-.3022654400209297D-01,
     *  -.2373365126402273D-02,0.3359539794362465D-01,
     1  -.4883359226323301D-01,0.4135176652532305D-01,
     2  -.1497903112523978D-01,-.1803824197609937D-01,
     3  0.4279466810008464D-01,-.4839453373821288D-01,
     4  0.3260488540698514D-01,-.2585768032265058D-02,
     5  0.1091149071070879D-03,-.6207015586696308D-03,
     6  0.2372635965366946D-02,-.6983955834391582D-02,
     7  0.1661897476706763D-01,-.3253626754632150D-01,
     8  0.5228713847707459D-01,-.6724659061647585D-01,
     9  0.6430180842120281D-01,-.3491475474863625D-01,
     *  -.1156909887938098D-01,0.4414388276684511D-01,
     1  -.3464989029786402D-01,-.7860269475681752D-02,
     2  0.4319158185042301D-01,-.3811969882660149D-01,
     3  -.6231439324998313D-03,0.3922714079370673D-01,
     4  -.4837673234396909D-01,0.2294604926101595D-01,
     5  0.1809356772750290D-01,-.4691526235792252D-01,
     6  0.4536249238002135D-01,-.1532401266462646D-01,
     7  -.2396653885780871D-01,0.4858080628401451D-01,
     8  -.4422957825993125D-01,0.1402290058944424D-01,
     9  0.2414631499181284D-01,-.4837300671896481D-01,
     *  0.4518146984139197D-01,-.1672068756238698D-01,
     1  -.2092011574020259D-01,0.4692722346362039D-01,
     2  -.4721015833680703D-01,0.2185687050517833D-01,
     3  0.3444982152080801D-03,-.1760461459490887D-02,
     4  0.5997220621564376D-02,-.1556670166063313D-01,
     5  0.3218315725455715D-01,-.5352408702922477D-01,
     6  0.7031337934117960D-01,-.6819828176228938D-01,
     7  0.3783158459766552D-01,0.1115824221801222D-01,
     8  -.4698077713394544D-01,0.3841952077730244D-01,
     9  0.6606703978579123D-02,-.4199263906826066D-01,
     *  0.3111246614058414D-01,0.1428397196107758D-01,
     1  -.4722553176505658D-01,0.3721933412110873D-01,
     2  0.5321567812944354D-02,-.4363503610471241D-01,
     3  0.4693880584310554D-01,-.1361241471448690D-01,
     4  -.3004175272987504D-01,0.5126011408506556D-01,
     5  -.3498357072899279D-01,-.6401807979591749D-02,
     6  0.4316286161780340D-01,-.4966436550869096D-01,
     7  0.2182561545831302D-01,0.2088459941991257D-01,
     8  -.4935554221697711D-01,0.4460465876915506D-01,
     9  -.1011642828416458D-01,-.3101157203127119D-01,
     *  0.5166453957389062D-01,-.3848438690157288D-01/
c 
        data arr3/
     1  0.9837077354227107D-03,-.4482628758323128D-02,
     2  0.1347364864839326D-01,-.3040339646402474D-01,
     3  0.5343212271969449D-01,-.7270719874170155D-01,
     4  0.7216424042092520D-01,-.4094270503987947D-01,
     5  -.1096813088737104D-01,0.4982073711188005D-01,
     6  -.4234385288824020D-01,-.5126114103930479D-02,
     7  0.4138507841824411D-01,-.2707559048110186D-01,
     8  -.2123300873465985D-01,0.4750340574742344D-01,
     9  -.2213447596044483D-01,-.2705481982057326D-01,
     *  0.5034268813131099D-01,-.2594422125991633D-01,
     1  -.2259556918000686D-01,0.5105390237084441D-01,
     2  -.3467542264560530D-01,-.1168709047387515D-01,
     3  0.4819006320771092D-01,-.4437645211127311D-01,
     4  0.3876181719917722D-02,0.3980050530814352D-01,
     5  -.5127499811578529D-01,0.2163813279160073D-01,
     6  0.2519734875346480D-01,-.5212853280295367D-01,
     7  0.3818038451746877D-01,0.5512071544270471D-02,
     8  -.4491426776444923D-01,0.4967566784100584D-01,
     9  0.2535756373093010D-02,-.1020649027377889D-01,
     *  0.2670506074151683D-01,-.5131525616736328D-01,
     1  0.7399142009798100D-01,-.7650162698490887D-01,
     2  0.4535165767948977D-01,0.9736351326357898D-02,
     3  -.5214657809616185D-01,0.4553517844140779D-01,
     4  0.4466458646514545D-02,-.4332027607553873D-01,
     5  0.2662528744417832D-01,0.2304500380756932D-01,
     6  -.4434177859731011D-01,0.9306095917740648D-02,
     7  0.3928427155201588D-01,-.4423465349802429D-01,
     8  0.3086879710890560D-03,0.4469996469301942D-01,
     9  -.4365688867292943D-01,-.1915655548520447D-02,
     *  0.4590049768445818D-01,-.4495274652083747D-01,
     1  0.3797088311831092D-03,0.4471625043322712D-01,
     2  -.4756270673553850D-01,0.5780561338958014D-02,
     3  0.4145932871857779D-01,-.5051206654607384D-01,
     4  0.1329117335235249D-01,0.3607159070067410D-01,
     5  -.5283185846599396D-01,0.2204213493548696D-01,
     6  0.2851532257410008D-01,-.5366911242404721D-01,
     7  0.5892940032421729D-02,-.2069368405197410D-01,
     8  0.4625361906403562D-01,-.7329789201245010D-01,
     9  0.8116732853721487D-01,-.5205640937780353D-01,
     *  -.6148713879131136D-02,0.5355449078707404D-01,
     1  -.4872570620068247D-01,-.3909612950252213D-02,
     2  0.4635170949575506D-01,-.2873985359254378D-01,
     3  -.2307505310915415D-01,0.4098506286635132D-01,
     4  -.1175035616627304D-02,-.4388158892815317D-01,
     5  0.3327918780526888D-01,0.2008520275596948D-01,
     6  -.5047135339682853D-01,0.2151168512341432D-01,
     7  0.3292408360522892D-01,-.5053930380700611D-01,
     8  0.1165645438780848D-01,0.4044946681526567D-01,
     9  -.4872920868386364D-01,0.4385558116729210D-02,
     *  0.4487440912547175D-01,-.4669860552674663D-01,
     1  -.6794091119888206D-03,0.4749553359677098D-01,
     2  -.4503243587774336D-01,-.4023260452874789D-02,
     3  0.4904148533086739D-01,-.4389544241191863D-01,
     4  -.6046882776546230D-02,0.4991734298226533D-01,
     5  0.1233530026658092D-01,-.3717632016637554D-01,
     6  0.6905324318463481D-01,-.8536763593479568D-01,
     7  0.6178892478909666D-01,-.1328203295540481D-02,
     8  -.5330992161725214D-01,0.5280838524115411D-01,
     9  0.2041298249018407D-02,-.4889155901627522D-01,
     *  0.3164476463706579D-01,0.2372721801717911D-01,
     1  -.4015685273075909D-01,-.2283073173364448D-02,
     2  0.4375334984940402D-01,-.2293888005965107D-01,
     3  -.3257133211094024D-01,0.4579420593373669D-01,
     4  0.1792631277156949D-02,-.4826157331925122D-01,
     5  0.3387227012458557D-01,0.2342074885993878D-01,
     6  -.5236617670352083D-01,0.1732082958283396D-01,
     7  0.3910230517980972D-01,-.4886852453997099D-01,
     8  0.3030841503068422D-03,0.4890729814415658D-01,
     9  -.4070704591806260D-01,-.1525164182114987D-01,
     *  0.5364523330816774D-01,-.2989969641072996D-01,
     1  -.2848352637763642D-01,0.5423321516863871D-01,
     2  -.1781916369690957D-01,-.3905903628189996D-01,
     3  0.2324426142847839D-01,-.5878740754468521D-01,
     4  0.8685934354197085D-01,-.7447788685546907D-01,
     5  0.1452888595580472D-01,0.4991959733174101D-01,
     6  -.5823953514444603D-01,0.2643323231191307D-02,
     7  0.5033528155306945D-01,-.3506034976779021D-01,
     8  -.2443418982417763D-01,0.4253635659368703D-01,
     9  0.2767030049810511D-02,-.4175084870540454D-01,
     *  0.1577435215205348D-01,0.3829054795656895D-01,
     1  -.3729564777643297D-01,-.1936448959376451D-01,
     2  0.5008994083944752D-01,-.1053217845960076D-01,
     3  -.4465304244019980D-01,0.3928253598473791D-01,
     4  0.1973446921455460D-01,-.5287363736265731D-01,
     5  0.1559274964131275D-01,0.4276505679639985D-01,
     6  -.4515244233975121D-01,-.1170055352161893D-01,
     7  0.5364158792635596D-01,-.2611866633432398D-01,
     8  -.3533742906285521D-01,0.5152363166697279D-01,
     9  -.1597428946164351D-02,-.5055641755035650D-01,
     *  0.3834159230606467D-01,0.2272711764414527D-01/
c 
        data arr4/
     1  0.3941761473755601D-01,-.8100979154245714D-01,
     2  0.8798076028821071D-01,-.3524526301211765D-01,
     3  -.4053399136510526D-01,0.6449564758142182D-01,
     4  -.1172462348100541D-01,-.4990041555946889D-01,
     5  0.3996870425085347D-01,0.2355585250978759D-01,
     6  -.4624990268929042D-01,-.1816781215201705D-02,
     7  0.4110566993377548D-01,-.1240171061017832D-01,
     8  -.3930163835435092D-01,0.2935233445879836D-01,
     9  0.2989808567069546D-01,-.4471865809510411D-01,
     *  -.9678759688545088D-02,0.5092751207344006D-01,
     1  -.1608422865812041D-01,-.4357194897169763D-01,
     2  0.3959589839273470D-01,0.2273737859335781D-01,
     3  -.5260948351330229D-01,0.6343195603888621D-02,
     4  0.4953549753186279D-01,-.3460621201468588D-01,
     5  -.3001989029197717D-01,0.5222185716921961D-01,
     6  -.2651136243973108D-03,-.5234378953761934D-01,
     7  0.3117825814853024D-01,0.3406037364780765D-01,
     8  -.5162111907622234D-01,-.3249592407143806D-02,
     9  0.6014578549539706D-01,-.9561991828848516D-01,
     *  0.6360132378623571D-01,0.2045685692878004D-01,
     1  -.6890288804880434D-01,0.2665751757462844D-01,
     2  0.4575475814228475D-01,-.4711487680206381D-01,
     3  -.1951358967566157D-01,0.5004856241601863D-01,
     4  -.7321029921604985D-03,-.4309122252676673D-01,
     5  0.1244667538306833D-01,0.3786197277972862D-01,
     6  -.2396402125750501D-01,-.3460934327934422D-01,
     7  0.3759108927257137D-01,0.2329251128232017D-01,
     8  -.4787458929429280D-01,-.5615817822081633D-02,
     9  0.5139018555116722D-01,-.1513716829075257D-01,
     *  -.4607796435967876D-01,0.3481585242473407D-01,
     1  0.3180523357564476D-01,-.4903380293595016D-01,
     2  -.1061759268135548D-01,0.5419854952835847D-01,
     3  -.1355077856487249D-01,-.4848393606817652D-01,
     4  0.3567920740510192D-01,0.3244634465021293D-01,
     5  -.5080840590894331D-01,-.9086344043011617D-02,
     6  0.5526161036173820D-01,-.1670276052986138D-01,
     7  0.8257267427676158D-01,-.9338819484434835D-01,
     8  0.1594969662694110D-01,0.6456488704953399D-01,
     9  -.4753686482170637D-01,-.3434788933595591D-01,
     *  0.5598833014480003D-01,0.1038897583933163D-01,
     1  -.5353026368776176D-01,0.6214239472984649D-02,
     2  0.4563477644466296D-01,-.1543594666403760D-01,
     3  -.3687281940033222D-01,0.2171642473278938D-01,
     4  0.3509831872196665D-01,-.3183361647255026D-01,
     5  -.3060450046488172D-01,0.4152860040358961D-01,
     6  0.2095398415917229D-01,-.4873429144101229D-01,
     7  -.7750150116350034D-02,0.5212315992783278D-01,
     8  -.7525271197234455D-02,-.5073981867898071D-01,
     9  0.2315794849082880D-01,0.4421000432657884D-01,
     *  -.3728857687139481D-01,-.3282965903944610D-01,
     1  0.4812496802454008D-01,0.1756139516489869D-01,
     2  -.5417799031030128D-01,0.5987155505130625D-04,
     3  0.5447427922009996D-01,-.1810263590204717D-01,
     4  -.4871099304736002D-01,0.3447505556794606D-01,
     5  0.1019960169348262D+00,-.6902746447807947D-01,
     6  -.3827627260063420D-01,0.6939213684799525D-01,
     7  0.1041673867709182D-01,-.6320931138121996D-01,
     8  0.6131614541410737D-02,0.5509578763016761D-01,
     9  -.1618002988999859D-01,-.4639564591824219D-01,
     *  0.2133913308464419D-01,0.3729394556021993D-01,
     1  -.2297991500384595D-01,-.3306156782280880D-01,
     2  0.2857559459355736D-01,0.3303999703428855D-01,
     3  -.3584142220550439D-01,-.2965349438125879D-01,
     4  0.4219294535254349D-01,0.2375871336617241D-01,
     5  -.4745014990501609D-01,-.1608116764944125D-01,
     6  0.5127782764536160D-01,0.7095709291585869D-02,
     7  -.5336253846331744D-01,0.2738723527763038D-02,
     8  0.5347399835945724D-01,-.1294394455711389D-01,
     9  -.5148383252961454D-01,0.2302679176267023D-01,
     *  0.4737345303758795D-01,-.3249549210593931D-01,
     1  -.4123550344223615D-01,0.4088042520515507D-01,
     2  0.3326976433506023D-01,-.4775588050628167D-01,
     3  0.1133583495309738D+00,-.2553808507604230D-01,
     4  -.7449791396225798D-01,0.2964855585413100D-01,
     5  0.5979036357471293D-01,-.3124640177447832D-01,
     6  -.5035091970232152D-01,0.3152025177929101D-01,
     7  0.4305882651538799D-01,-.3062140848678901D-01,
     8  -.3640587293424975D-01,0.2806075616367769D-01,
     9  0.3071099652489073D-01,-.2798640569072092D-01,
     *  -.3178945137393960D-01,0.3255486082937727D-01,
     1  0.3310526787977260D-01,-.3627161801465869D-01,
     2  -.3271386995325684D-01,0.3938029313586181D-01,
     3  0.3145033775435960D-01,-.4217036555770997D-01,
     4  -.2961019939486667D-01,0.4470926234278143D-01,
     5  0.2733460559094777D-01,-.4701304233019133D-01,
     6  -.2471100394483223D-01,0.4908123464726923D-01,
     7  0.2180141573960736D-01,-.5090725599671978D-01,
     8  -.1865397398824384D-01,0.5248246835447075D-01,
     9  0.1530858457535782D-01,-.5379799459303460D-01,
     *  -.1180001247118369D-01,0.5484560620297043D-01/
c 
c        return the expansion matrix
c 
        do 1400 i=1,36
        do 1200 j=1,20
c 
        arrout(i,j)=arr(i,j)
c 
 1200 continue
 1400 continue
c 
        d=-1
        do 2400 i=1,36
c 
        d=-d
        do 2200 j=1,20
c 
        arrout(i,41-j)=arr(i,j)*d
c 
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
        subroutine expons_ret(expons,coefs,tt40,ww40)
        implicit real *8 (a-h,o-z)
        save
        dimension cs1(36),cs2(36),cs3(36),cs4(36),
     1       cs5(36),cs6(36),cs7(36),cs8(36),cs9(36),
     2       cs10(36),cs11(36),cs12(36),cs13(36),cs14(36),
     3       cs15(36),cs16(36),cs17(36),cs18(36),cs19(36),
     4       cs20(36),cs21(36),cs22(36),cs23(36),cs24(36),
     5       cs25(36),cs26(36),cs27(36),cs28(36),cs29(36),
     6       cs30(36),expons0(36),cs0(36),tt0(20),ww0(20),
     7       cs31(36),cs32(36),cs33(36),cs34(36),cs35(36),
     8       cs36(36)
c 
        dimension expons(72),tt40(40),ww40(40),expons72(72)
        complex *16 coefs(72,37),cd0,cd1,cd2,cd3,ima
c 
        data ima/(0.0d0,1.0d0)/
c 
        data tt0/
     1  -.99796430969563536201D+00,-.98930906979315713308D+00,
     2  -.97387786032577240943D+00,-.95189760341409325544D+00,
     3  -.92368031722547892847D+00,-.88960200962536979564D+00,
     4  -.85008309626971718513D+00,-.80557018095505792891D+00,
     5  -.75652045315962213320D+00,-.70338938289267660626D+00,
     6  -.64662181531460088732D+00,-.58664617113998871016D+00,
     7  -.52387125267017717418D+00,-.45868509646901049738D+00,
     8  -.39145534786881832704D+00,-.32253071322792858786D+00,
     9  -.25224314089450271412D+00,-.18091047228154398668D+00,
     *  -.10883938130975203874D+00,-.36328480859055730185D-01/
c 
        data ww0/
     *  0.52202210263669779144D-02,0.12072664924297736051D-01,
     *  0.18751412133019550035D-01,0.25156768524833940654D-01,
     *  0.31214806062386796015D-01,0.36871432999101159839D-01,
     *  0.42091698120202544379D-01,0.46857729013469706818D-01,
     *  0.51165667264918386616D-01,0.55022236998821941073D-01,
     *  0.58441442839184497621D-01,0.61441693937994097497D-01,
     *  0.64043477544761927523D-01,0.66267585668338907402D-01,
     *  0.68133830764931244004D-01,0.69660158402707219656D-01,
     *  0.70862061968500511279D-01,0.71752214889286323233D-01,
     *  0.72340251621454719922D-01,0.72632645295421812467D-01/
c 
        data expons0/
     *  -.39965454976041159016D+02,-.39818629074803446199D+02,
     *  -.39557049020318090863D+02,-.39184834122463674361D+02,
     *  -.38707552668848131108D+02,-.38131775753315235134D+02,
     *  -.37464669874607238146D+02,-.36713636447527061786D+02,
     *  -.35886021783664774005D+02,-.34988905024461485255D+02,
     *  -.34028958724322872238D+02,-.33012369826673918062D+02,
     *  -.31944806618911777646D+02,-.30831418103985343988D+02,
     *  -.29676854510115180153D+02,-.28485300313340347219D+02,
     *  -.27260513592184186948D+02,-.26005867533080991185D+02,
     *  -.24724391420427078905D+02,-.23418809531347891300D+02,
     *  -.22091577099345482900D+02,-.20744912999652091310D+02,
     *  -.19380829115885856237D+02,-.18001156529565628750D+02,
     *  -.16607568773120940047D+02,-.15201602432849922464D+02,
     *  -.13784675401138727038D+02,-.12358103070735553137D+02,
     *  -.10923112746834613350D+02,-.94808565308181889580D+01,
     *  -.80324229062267534860D+01,-.65788472349674645540D+01,
     *  -.51211213511154765680D+01,-.36602024215301279528D+01,
     *  -.21970212271769737332D+01,-.73249000659931290264D+00/
c 
        data cs0/
     8  0.16598133393340511474D-11,-.54912924064271864264D-11,
     9  0.96418426970218691493D-11,-.12682004541047750731D-10,
     *  0.13153810778713020803D-10,-.97180491334235549588D-11,
     1  0.15277007138415216008D-11,0.11229099474707802179D-10,
     2  -.26744019062154172536D-10,0.41205170156186770995D-10,
     3  -.48948962892386682928D-10,0.43466445158986411916D-10,
     4  -.19396152811312900121D-10,-.24738507566498950409D-10,
     5  0.83760671915550742529D-10,-.14411384730725686453D-09,
     6  0.18421995355577919964D-09,-.17731539836539055796D-09,
     7  0.96113873959994667402D-10,0.79113775099066823004D-10,
     8  -.29525635742936850797D-09,-.66464719597614951900D-09,
     9  0.14576872258149030410D-06,0.44450515847403288503D-05,
     *  0.43031929739024316500D-04,0.24833140364000258206D-03,
     1  0.10309254938134704719D-02,0.33588210804566049593D-02,
     2  0.90164787818910623953D-02,0.20554269146662610561D-01,
     3  0.40596055196286689345D-01,0.70434251335865259616D-01,
     4  0.10839812808891954803D+00,0.14899070280046933905D+00,
     5  0.18374338541866243022D+00,0.20390880589726110440D+00/
c 
        data cs1/
     4  -.74092634433400728197D-10,0.24881840387591262933D-09,
     5  -.45082770047533665727D-09,0.62895174636430194813D-09,
     6  -.73217540442804008428D-09,0.71245092776081305698D-09,
     7  -.53462904785463087194D-09,0.19092922886089786368D-09,
     8  0.28388614049002325084D-09,-.80296102385551594252D-09,
     9  0.12311105274727797646D-08,-.14080591246054309270D-08,
     *  0.11918526110980556496D-08,-.51442543945593616997D-09,
     1  -.56776782207100196393D-09,0.18499719234616667502D-08,
     2  -.30003669261061244405D-08,0.36314516679153067421D-08,
     3  -.34142469519922814577D-08,0.22779075970981808466D-08,
     4  -.14133115398946349369D-08,0.15105906583396598120D-07,
     5  -.14201012899240717112D-05,-.34754132519981302755D-04,
     6  -.28571133057469950224D-03,-.14229108651106229268D-02,
     7  -.51204102608572521234D-02,-.14430870859118978291D-01,
     8  -.33266712433425961102D-01,-.64300802464786902238D-01,
     9  -.10558349471757681014D+00,-.14780197743838334340D+00,
     *  -.17502701912937580581D+00,-.17049534008566768821D+00,
     1  -.12551429849197205795D+00,-.46312438664976915480D-01/
c 
        data cs2/
     *  0.52938213377820932580D-09,-.18010877698574248990D-08,
     1  0.33448475940421365375D-08,-.48523804250921932726D-08,
     2  0.59960956962809078795D-08,-.64212061592198192692D-08,
     3  0.57889876995258604784D-08,-.38641994398582319632D-08,
     4  0.63577371789067570443D-09,0.35588374329156940038D-08,
     5  -.79560625676270883701D-08,0.11396823392864396250D-07,
     6  -.12515212351448307333D-07,0.10092505176184554259D-07,
     7  -.35141852508384782155D-08,-.68083022972742191349D-08,
     8  0.19239082565541615717D-07,-.31131328382308408486D-07,
     9  0.39649026924198731746D-07,-.43913172750694549453D-07,
     *  0.54523316159571709852D-07,-.17854632277262712178D-06,
     1  0.95277329540681743842D-05,0.18561130480905037320D-03,
     2  0.12863187739779792560D-02,0.54812989744750289047D-02,
     3  0.16912946780055230383D-01,0.40632970939423448333D-01,
     4  0.78778304354471557795D-01,0.12508139038168978837D+00,
     5  0.16204505895769234759D+00,0.16597725088644784770D+00,
     6  0.12060634663958623509D+00,0.31912955088473798991D-01,
     7  -.67859019427010405272D-01,-.13388941509307075280D+00/
c 
        data cs3/
     *  -.96689949158090696718D-09,0.34027553878356584621D-08,
     *  -.67097056660651853630D-08,0.10601553866664881890D-07,
     *  -.14651440086519924169D-07,0.18129890027625994764D-07,
     *  -.19931004679211422811D-07,0.18635658686532446440D-07,
     *  -.12788168074573986049D-07,0.14347928538841069527D-08,
     *  0.15116136685589584712D-07,-.34495661407839902206D-07,
     *  0.51775335785511916890D-07,-.59690528652116642403D-07,
     *  0.49921492526074927887D-07,-.15438553895670196844D-07,
     *  -.46648616684803894311D-07,0.13304668448704433750D-06,
     *  -.23553851977702506526D-06,0.35252712257453790663D-06,
     *  -.54973960200720166552D-06,0.14204861582316877873D-05,
     *  -.50715669513947918033D-04,-.77895913733360480796D-03,
     *  -.45110028087796618874D-02,-.16269828754680947894D-01,
     *  -.42425316505771289763D-01,-.85079010069275288649D-01,
     *  -.13419169584849875583D+00,-.16499324032105547438D+00,
     *  -.14864165468873245893D+00,-.74472002360065902841D-01,
     *  0.32110271253478565550D-01,0.11589295552997740766D+00,
     *  0.12576446215074583014D+00,0.53751851961329839810D-01/
c 
        data cs4/
     *  0.13711065247621895776D-07,-.47230474569116775496D-07,
     *  0.89788762928124961876D-07,-.13512151311147546196D-06,
     *  0.17633184040492636788D-06,-.20504087322678692701D-06,
     *  0.21162613471160500918D-06,-.18650587936815527107D-06,
     *  0.12254589647352787041D-06,-.18390685005933277305D-07,
     *  -.11794036496328529871D-06,0.26634650780933296055D-06,
     *  -.39428047778802821500D-06,0.46000369050309297851D-06,
     *  -.41989879686823636849D-06,0.23875077862989384591D-06,
     *  0.10012825532926618519D-06,-.59063875238510733595D-06,
     *  0.12165323381534735954D-05,-.20151747147982443603D-05,
     *  0.33774944519104237730D-05,-.80023377309793202203D-05,
     *  0.22573060966153806042D-03,0.27126623887652965160D-02,
     *  0.12983030456344626372D-01,0.39073155032659181874D-01,
     *  0.84358107970521108287D-01,0.13661224062450022440D+00,
     *  0.16478892368213881954D+00,0.13535282053607959639D+00,
     *  0.43534573520437312692D-01,-.67309984953854251122D-01,
     *  -.12607693140378493582D+00,-.88823824744649530985D-01,
     *  0.17231299952005125951D-01,0.10736146382141686537D+00/
c 
        data cs5/
     *  0.92956798927526196345D-07,-.30768032469216100814D-06,
     *  0.54183038913886176422D-06,-.72039690947736662273D-06,
     *  0.77293435205288576819D-06,-.64290003307418478879D-06,
     *  0.30593124120958919256D-06,0.20867962916273875925D-06,
     *  -.80189576975307681747D-06,0.13059857386458986524D-05,
     *  -.15139394724920339494D-05,0.12415586927783101917D-05,
     *  -.41286352163697975587D-06,-.85732332705453208542D-06,
     *  0.22218051821002643677D-05,-.31362540811453162959D-05,
     *  0.29832309597437316357D-05,-.12512639221734453121D-05,
     *  -.23307500734669258641D-05,0.79306393657879751823D-05,
     *  -.16819684722223876167D-04,0.38205432671765961423D-04,
     *  -.86781570260922582214D-03,-.80474400986228592747D-02,
     *  -.31390460028550003224D-01,-.77297751893895738796D-01,
     *  -.13388911782794141823D+00,-.16487598292106393920D+00,
     *  -.13040551190892009072D+00,-.28394688090185583369D-01,
     *  0.83428717990419931549D-01,0.12216603310914471097D+00,
     *  0.54678555685193477566D-01,-.59568192055305999470D-01,
     *  -.11415015222649200773D+00,-.56635856992618942775D-01/
c 
        data cs6/
     *  -.12111952784073076664D-06,0.39874537448245023285D-06,
     *  -.69307449575648154365D-06,0.89464175131379202768D-06,
     *  -.89332113626546796751D-06,0.59167032646931314212D-06,
     *  0.64789231450992392602D-07,-.10455025631131355796D-05,
     *  0.21922291463777036166D-05,-.31961136412507771373D-05,
     *  0.36210126334530127212D-05,-.29966612932327976005D-05,
     *  0.98787414725562920195D-06,0.23861315017200396440D-05,
     *  -.65464784022729979115D-05,0.10253649249603382119D-04,
     *  -.11689842145618013033D-04,0.87288299884542318032D-05,
     *  0.78159533020614994522D-06,-.19359422766896186612D-04,
     *  0.53378792215166969138D-04,-.13068076887244887714D-03,
     *  0.29045985822808107941D-02,0.20607581351121801412D-01,
     *  0.64270899002363460141D-01,0.12575034885968875098D+00,
     *  0.16543418161048340138D+00,0.13428348041701800025D+00,
     *  0.27151063250793306970D-01,-.87980529066140190735D-01,
     *  -.11655601688358545277D+00,-.32488148533783822351D-01,
     *  0.79872257486689383667D-01,0.10301959949288450969D+00,
     *  0.12740251666528436369D-01,-.90335853646643192133D-01/
c 
        data cs7/
     *  0.22050008034035438315D-05,-.73907530106590841601D-05,
     *  0.13342811335009668956D-04,-.18508207553511885397D-04,
     *  0.21357833637856440494D-04,-.20490656125077453401D-04,
     *  0.14947657339613329317D-04,-.46482336486193748389D-05,
     *  -.91800282591609650576D-05,0.23780620751480157743D-04,
     *  -.35089322929022567312D-04,0.38563540307388819401D-04,
     *  -.30595810296466511716D-04,0.10205421902388300813D-04,
     *  0.19567711756452411434D-04,-.51245296530005358951D-04,
     *  0.73938766703212555997D-04,-.75712289982429838334D-04,
     *  0.46232754792089455889D-04,0.22889551903677446185D-04,
     *  -.14731218307500814924D-03,0.38330681393671185218D-03,
     *  -.85812326358247379604D-02,-.45721582810451900276D-01,
     *  -.11089732494395502205D+00,-.16429389014411892821D+00,
     *  -.14552543797854048476D+00,-.38246613910983555329D-01,
     *  0.84551725353248459486D-01,0.11385478053929716242D+00,
     *  0.22717557989888327149D-01,-.87302608105132681809D-01,
     *  -.90160006450337903632D-01,0.15434458202317047247D-01,
     *  0.99110163933073492648D-01,0.57119935910139331164D-01/
c 
        data cs8/
     *  -.59306456331021056743D-05,0.20126842105447333431D-04,
     *  -.37199922079373488591D-04,0.53558136894026259156D-04,
     *  -.65423315651560691368D-04,0.68808853020366501538D-04,
     *  -.60087691589515243503D-04,0.37051286839581737699D-04,
     *  -.29376039626440871876D-06,-.45473673960454742008D-04,
     *  0.90807171863808454708D-04,-.12233387236358129103D-03,
     *  0.12548959709457787903D-03,-.89123211254196286743D-04,
     *  0.11028594331572446181D-04,0.97386289452129098178D-04,
     *  -.20942912759052100763D-03,0.28591247210689024326D-03,
     *  -.28174018818125405947D-03,0.15128601587169965685D-03,
     *  0.16266069807851864440D-03,-.73412771711514636298D-03,
     *  0.22290907498638338071D-01,0.87488834745439734197D-01,
     *  0.15801498142229657428D+00,0.16081589175968006661D+00,
     *  0.60804096784808511079D-01,-.73449981110184660999D-01,
     *  -.11452313661879593127D+00,-.23484630096802347460D-01,
     *  0.87511824371285084503D-01,0.81811896809516033694D-01,
     *  -.29781954913385945072D-01,-.97461858887894542908D-01,
     *  -.30733001258203407246D-01,0.77159069316987931715D-01/
c 
        data cs9/
     *  0.11883455710686210677D-04,-.40492063852125449459D-04,
     *  0.75412983332076096034D-04,-.10988054726785931127D-03,
     *  0.13663365612141838846D-03,-.14763738767021744130D-03,
     *  0.13493513742222018763D-03,-.92521169453226134464D-04,
     *  0.19075423149556117416D-04,0.78996866411082897848D-04,
     *  -.18527869837004788851D-03,0.27323774317117001407D-03,
     *  -.30969268075136687946D-03,0.26253651024456079249D-03,
     *  -.11200966377672958226D-03,-.13695536877260545811D-03,
     *  0.44522701610875180786D-03,-.73679426908328864700D-03,
     *  0.90446963262253842652D-03,-.81898848372031164528D-03,
     *  0.32531780336425922671D-03,0.52969185624445714792D-03,
     *  -.50595528335901788023D-01,-.14250891883514424176D+00,
     *  -.17619042858709975303D+00,-.92530651392561777930D-01,
     *  0.52179706892570063823D-01,0.11681969971878731466D+00,
     *  0.32890493083598926424D-01,-.83224637037862154226D-01,
     *  -.78588775133349173860D-01,0.33645296837996448078D-01,
     *  0.93363923732458834404D-01,0.14669057958038840590D-01,
     *  -.83430354089774781794D-01,-.55156172340294559770D-01/
c 
        data cs10/
     *  -.51310501646367246379D-04,0.17540753286322306379D-03,
     *  -.32875320685874649230D-03,0.48401280035308422450D-03,
     *  -.61190995737815356031D-03,0.67961201408210559217D-03,
     *  -.65372438457956629070D-03,0.50736695419411751749D-03,
     *  -.23071092620459645932D-03,-.15718080274864777903D-03,
     *  0.59926463923480256317D-03,-.99955074564125419943D-03,
     *  0.12352579709587062381D-02,-.11841689549987761730D-02,
     *  0.76385836695187678437D-03,0.26982264195352587914D-04,
     *  -.10796139086996610450D-02,0.21729172481366849972D-02,
     *  -.30042750325398698120D-02,0.32462881658844937517D-02,
     *  -.26272129977410357834D-02,0.20114139747449919664D-02,
     *  0.10021438075172107661D+00,0.19084198163461843265D+00,
     *  0.13225724427914633327D+00,-.21294375496387329122D-01,
     *  -.11477104155851989691D+00,-.50626252480628151503D-01,
     *  0.75399636300863800795D-01,0.80383907161813452457D-01,
     *  -.31453564078078910787D-01,-.88937510304112521863D-01,
     *  -.80959286406293930417D-02,0.83365740343035567236D-01,
     *  0.40123425163454311945D-01,-.65810719871614436014D-01/
c 
        data cs11/
     *  -.42207611890087579436D-04,0.13585295304807089960D-03,
     *  -.22564485502810727010D-03,0.26836559296338366827D-03,
     *  -.22675586790948216176D-03,0.78696811314301571631D-04,
     *  0.16971171350746075705D-03,-.47200547772982128494D-03,
     *  0.73725796108206075540D-03,-.84241613633752501422D-03,
     *  0.66698559874454790237D-03,-.14822685221870116671D-03,
     *  -.65719209975405766740D-03,0.15315560757179751620D-02,
     *  -.21063361637904949790D-02,0.19471777846441045075D-02,
     *  -.71034744139741968036D-03,-.16701234509616442914D-02,
     *  0.48393043866548811688D-02,-.80146939423323861464D-02,
     *  0.10287667480112147785D-01,-.14070681280708892643D-01,
     *  -.16716472186217581430D+00,-.20225980203839024476D+00,
     *  -.19613206375086202799D-01,0.11348473593416150593D+00,
     *  0.68099650729461062181D-01,-.58328935565350546811D-01,
     *  -.92000662916995313698D-01,0.25336566539943563060D-01,
     *  0.91989214131239369860D-01,0.34715620534019523142D-02,
     *  -.79899403369304093152D-01,-.34273738684236144932D-01,
     *  0.69449497598363331481D-01,0.57464958572242988069D-01/
c 
        data cs12/
     *  -.13880592897751296975D-03,0.47996386840310922638D-03,
     *  -.91884641100609615286D-03,0.13973748303710358314D-02,
     *  -.18506404656621930478D-02,0.21963805710821017044D-02,
     *  -.23348771990696020052D-02,0.21593577325119617702D-02,
     *  -.15788045543085013707D-02,0.55280268012550701010D-03,
     *  0.86666127740371145169D-03,-.24976082708259036294D-02,
     *  0.40111567259752571788D-02,-.49532653008780755508D-02,
     *  0.48177549554199065156D-02,-.31695918109151768658D-02,
     *  -.20292146529610883287D-03,0.51505503064633075509D-02,
     *  -.11103297658828224029D-01,0.17187552251136629581D-01,
     *  -.22925330973470871420D-01,0.36260934075529592327D-01,
     *  0.25695461207229648754D+00,0.14880628305210794306D+00,
     *  -.13311936815230085997D+00,-.10972487785056411509D+00,
     *  0.54302442581896793032D-01,0.10862981599595375870D+00,
     *  -.64394904417490413337D-02,-.11494168764982161517D+00,
     *  -.10368938518646568545D-01,0.98491045670363124552D-01,
     *  0.32210781819588841708D-01,-.79236728416626556839D-01,
     *  -.60881303471304901343D-01,0.73642285595279412301D-01/
c 
        data cs13/
     *  -.63756131833678222696D-03,0.21153706748016868328D-02,
     *  -.37437009195123801268D-02,0.50225071052896802699D-02,
     *  -.54810409663607684909D-02,0.47359776512352925456D-02,
     *  -.26090744283639938752D-02,-.73362102746902432616D-03,
     *  0.46767265462330515918D-02,-.81605427559353699923D-02,
     *  0.98680075594975828769D-02,-.86197278750134946240D-02,
     *  0.39137859025895071909D-02,0.35620626042713407306D-02,
     *  -.11695817316884917551D-01,0.17228778976484754204D-01,
     *  -.16631094620813709890D-01,0.73954684193543872833D-02,
     *  0.10734333364068069305D-01,-.35360787644008408585D-01,
     *  0.63406978110105166941D-01,-.11154044643673834113D+00,
     *  -.36903032911166459608D+00,0.57352829965781680550D-02,
     *  0.29940883028267370117D+00,-.38382892613290485471D-01,
     *  -.21518788554500445999D+00,0.49061508102214314122D-02,
     *  0.15986737404708747015D+00,0.53708129843930733503D-01,
     *  -.15360354253888452500D+00,-.81821442849998133249D-01,
     *  0.15346662696736635481D+00,0.82200194246388292991D-01,
     *  -.12741236118613421679D+00,-.93784861242995310125D-01/
c 
        data cs14/
     *  0.27152848558123503225D-03,-.87882862276024373005D-03,
     *  0.14724248328424880008D-02,-.17645758336235842681D-02,
     *  0.14727763340960278195D-02,-.36763339642495631930D-03,
     *  -.16366550116868740162D-02,0.43662509797730148862D-02,
     *  -.72826909577713819570D-02,0.94506502959359316343D-02,
     *  -.96533825231461938359D-02,0.67136540433607616259D-02,
     *  -.16693897855867975705D-04,-.98709390258621231106D-02,
     *  0.20706563933991944331D-01,-.28488451531024092009D-01,
     *  0.28025749568015947117D-01,-.14189628778430358733D-01,
     *  -.16527018708758074622D-01,0.65236275154417612865D-01,
     *  -.13448541120101400064D+00,0.27858100810116806035D+00,
     *  0.52461624139057266933D+00,-.44866394089156159975D+00,
     *  -.37827297983476303366D+00,0.55067338618579327956D+00,
     *  0.16369174818674125305D+00,-.40421295864480481337D+00,
     *  -.14451379868367447280D+00,0.27876022773640457856D+00,
     *  0.24495750736108237319D+00,-.30570213536611016066D+00,
     *  -.24921133958825635403D+00,0.32122273270510464543D+00,
     *  0.19487789699488119574D+00,-.23746608327441316457D+00/
c 
        data cs15/
     *  -.59855098484896838597D-02,0.20027597162422587224D-01,
     *  -.36038016838475063522D-01,0.49729478039412616919D-01,
     *  -.56929930177323388195D-01,0.53915137732256366422D-01,
     *  -.38301888258816874819D-01,0.10237438254455097311D-01,
     *  0.26501596415588658961D-01,-.64080332034495922067D-01,
     *  0.91476913553641741611D-01,-.97043609446552988194D-01,
     *  0.72557755370717859234D-01,-.17742756455720902847D-01,
     *  -.56484586103899859383D-01,0.12796169538002064180D+00,
     *  -.16760522285075357541D+00,0.14759266726060659128D+00,
     *  -.50565999458155349954D-01,-.12570622091240406808D+00,
     *  0.38132496758257580668D+00,-.86175995680607671994D+00,
     *  -.40836158394140640857D+00,0.14595153195523633775D+01,
     *  -.33818071168907580590D+00,-.13780189768011925231D+01,
     *  0.82791014078656586180D+00,0.93412406010175130107D+00,
     *  -.79478691393517169648D+00,-.72445306632343396237D+00,
     *  0.65409344445904437169D+00,0.73973335254640331419D+00,
     *  -.65708131675042207761D+00,-.72841714305024808780D+00,
     *  0.68637972165125240202D+00,0.69016829308143045602D+00/
c 
        data cs16/
     *  0.12969108095960423206D-01,-.43947837966615436553D-01,
     *  0.80999983165099621565D-01,-.11610645691500397495D+00,
     *  0.14089940266921432735D+00,-.14670793655415964759D+00,
     *  0.12590757548997658468D+00,-.74318812668863473491D-01,
     *  -.57774961426047057732D-02,0.10296155155584737955D+00,
     *  -.19578377663849563957D+00,0.25523805246834344239D+00,
     *  -.25158574512520509498D+00,0.16506297291363039914D+00,
     *  0.20536314350316635201D-02,-.21643285779720133749D+00,
     *  0.41271201126053089713D+00,-.50371582861342592118D+00,
     *  0.40013014582609810717D+00,-.28690559688470733346D-01,
     *  -.68945336889735242557D+00,0.22462990504547219805D+01,
     *  -.91034484355833242852D+00,-.28800266559077344225D+01,
     *  0.34763980607424565173D+01,0.11973377614248634449D+01,
     *  -.43441312538036923161D+01,0.88554115973340470871D+00,
     *  0.34229852337339691630D+01,-.16085497865699096803D+01,
     *  -.24765132014630010722D+01,0.14061541682606051097D+01,
     *  0.23936099119246397179D+01,-.15144845207224339996D+01,
     *  -.25189012010289238619D+01,0.21273105690380762538D+01/
c 
        data cs17/
     *  -.30937428911792924522D-01,0.10506308129451802564D+00,
     *  -.19445721079209546930D+00,0.28066799537251250307D+00,
     *  -.34432802627823171046D+00,0.36491498918200932370D+00,
     *  -.32346333592250063972D+00,0.20772770567620238275D+00,
     *  -.19054731228223738317D-01,-.22084569404477552441D+00,
     *  0.46617936519354590281D+00,-.64990555460878720355D+00,
     *  0.69592369361374006780D+00,-.54097512836249336511D+00,
     *  0.16273738510228662138D+00,0.39359396833150985613D+00,
     *  -.10004044516956722437D+01,0.14529590946396812521D+01,
     *  -.14950822415357409714D+01,0.84415762180601957733D+00,
     *  0.88555485128679657442D+00,-.53860314348894980614D+01,
     *  0.63350990112444944918D+01,0.22690774964943666836D+01,
     *  -.12068708986034385483D+02,0.74221273831011530252D+01,
     *  0.83588857447948789455D+01,-.13017300154941300588D+02,
     *  -.50732357605973304414D+00,0.10998656141635743186D+02,
     *  -.28726323391741215646D+01,-.82633730648991592835D+01,
     *  0.28951856242124619939D+01,0.87069893497701935936D+01,
     *  -.50487396700990723946D+01,-.83217051373327542686D+01/
c 
        data cs18/
     *  0.16427934139640404724D+00,-.56149614503388193654D+00,
     *  0.10520555592272768760D+01,-.15483474953732957630D+01,
     *  0.19569024275514692524D+01,-.21735404216083584263D+01,
     *  0.20932823468993994461D+01,-.16331659968203357798D+01,
     *  0.76568583964960599074D+00,0.44442391041022446698D+00,
     *  -.18133788663773257728D+01,0.30387833946128804799D+01,
     *  -.37428822502165448000D+01,0.35641593422003102551D+01,
     *  -.22848957442787493076D+01,-.38460540385598303439D-01,
     *  0.29894080390609989208D+01,-.57935099010194598727D+01,
     *  0.74199704502442046938D+01,-.67161717549016004253D+01,
     *  0.22532006924096773987D+01,0.11254973988896362735D+02,
     *  -.25457930809032400896D+02,0.16319148526427667181D+02,
     *  0.20191221503370064120D+02,-.42837367620017525916D+02,
     *  0.13983888030388161777D+02,0.35045786188836915831D+02,
     *  -.37012438934017103747D+02,-.10617280192618925551D+02,
     *  0.36138871418361937808D+02,-.31756448817144187117D+01,
     *  -.31685940988502414442D+02,0.10766097241079246541D+02,
     *  0.30419735711442073329D+02,-.22425097204772994459D+02/
c 
        data cs19/
     *  0.55065647552885792054D+00,-.18211368505428812884D+01,
     *  0.32025820909518982107D+01,-.42506032183179496520D+01,
     *  0.45540259315679373963D+01,-.37929027079047741260D+01,
     *  0.18462998750107555976D+01,0.10872648964443293790D+01,
     *  -.44102175933621380614D+01,0.71564047126936902740D+01,
     *  -.81916186441807034543D+01,0.65948056766576097443D+01,
     *  -.21424531054133697852D+01,-.42814731046165743797D+01,
     *  0.10588480942151731348D+02,-.13904025602971999530D+02,
     *  0.11542186918251635961D+02,-.23228042032296472413D+01,
     *  -.12214690773324811428D+02,0.27152180530070124279D+02,
     *  -.33778306655603910382D+02,0.10180016025161264395D+02,
     *  0.48682852452385471176D+02,-.86549351793578535502D+02,
     *  0.27666185264008315707D+02,0.95893410342395713251D+02,
     *  -.13066977888967494785D+03,-.88383467344189454040D+00,
     *  0.14980957792038912969D+03,-.11292152912788085036D+03,
     *  -.70396391085985277377D+02,0.15025008864933896663D+03,
     *  -.22430952043682099165D+02,-.11610344650981165992D+03,
     *  0.64265679346217918531D+02,0.79848784701011516636D+02/
c 
        data cs20/
     *  -.47844575600005078222D+00,0.15614152912316813132D+01,
     *  -.26647981178412117490D+01,0.33217331230963855831D+01,
     *  -.30746291659884561138D+01,0.15461548010341329848D+01,
     *  0.14215715500898675491D+01,-.55882388764203255454D+01,
     *  0.10155368496655531048D+02,-.13721605181507102718D+02,
     *  0.14468047434804145331D+02,-.10657359554546478629D+02,
     *  0.14290496004467046519D+01,0.12291826914668723368D+02,
     *  -.27134254840593623473D+02,0.37332534401822589147D+02,
     *  -.35875856202744106513D+02,0.16836216816660311579D+02,
     *  0.21506478548495941940D+02,-.73412780557075948768D+02,
     *  0.12166918548554666192D+03,-.10631880572099665243D+03,
     *  -.47248482092640187633D+02,0.27974272014410715637D+03,
     *  -.32567814906399707472D+03,-.11806502967668104449D+02,
     *  0.47466225623730665493D+03,-.49712348347358222769D+03,
     *  -.71593248359459842023D+02,0.61534498846903418283D+03,
     *  -.45024214847502669460D+03,-.23540114669170120207D+03,
     *  0.56351453677008030863D+03,-.14499254774045190420D+03,
     *  -.38191499322167841202D+03,0.26847787483853688583D+03/
c 
        data cs21/
     *  0.14179249339895559064D+02,-.47560408097179070536D+02,
     *  0.85994426441294944984D+02,-.11963201154155296908D+03,
     *  0.13882188771191024758D+03,-.13476263155383463622D+03,
     *  0.10146683157648713824D+03,-.38458019860897017423D+02,
     *  -.46637766014046313650D+02,0.13688975694520963280D+03,
     *  -.20768616963382692775D+03,0.23213779135671138343D+03,
     *  -.18994696799003741804D+03,0.77603304152012615053D+02,
     *  0.83920899888845537287D+02,-.24836894431763826570D+03,
     *  0.35329793544355163222D+03,-.33864248351376998818D+03,
     *  0.16986037842205636814D+03,0.14243000798537819991D+03,
     *  -.53147793699236109582D+03,0.79825350662854550104D+03,
     *  -.57696952668789858900D+03,-.33909294995686552414D+03,
     *  0.14151012883188285246D+04,-.14597862533822350638D+04,
     *  -.80209197279818740029D+02,0.19396055973867830407D+04,
     *  -.19829659458391923694D+04,-.18408181385747171161D+03,
     *  0.22453351658658795944D+04,-.17300493554021982213D+04,
     *  -.77324810790082312322D+03,0.21395507190689150740D+04,
     *  -.68754501617795614924D+03,-.15676930082313495291D+04/
c 
        data cs22/
     *  -.57750796788614810294D+02,0.19657607396599510938D+03,
     *  -.36546398085436494969D+03,0.53135141830987665096D+03,
     *  -.65950420387040696506D+03,0.71280066311624120091D+03,
     *  -.65645929262275146725D+03,0.46725477795641749057D+03,
     *  -.14578721779096665686D+03,-.27164482338928623879D+03,
     *  0.70710154848438646579D+03,-.10473451205465403106D+04,
     *  0.11670056902246476621D+04,-.96791425971427142998D+03,
     *  0.42613559524539452849D+03,0.36949952193768585894D+03,
     *  -.12070867380425682062D+04,0.17846815365654250083D+04,
     *  -.17876836022353842379D+04,0.99120247215049272296D+03,
     *  0.64643759359029203067D+03,-.27297235233397369290D+04,
     *  0.40116896807742537972D+04,-.26764713179228620421D+04,
     *  -.17707864520633416127D+04,0.64262055843483308818D+04,
     *  -.62907484589025349280D+04,-.32120362376112922877D+03,
     *  0.79596631221544515814D+04,-.83024045905710100456D+04,
     *  -.20228497299184287841D+03,0.88736845901598615888D+04,
     *  -.77995618616700552887D+04,-.21167851340230384750D+04,
     *  0.93331694956839052578D+04,-.52274875089858484525D+04/
c 
        data cs23/
     *  -.15234107159632294784D+03,0.50211389695641321201D+03,
     *  -.87707061676787851202D+03,0.11506832316516718783D+04,
     *  -.12079421676912061304D+04,0.96362093842478433974D+03,
     *  -.39375326068245809396D+03,-.43119737623161679751D+03,
     *  0.13277423831743566974D+04,-.20153629963005600669D+04,
     *  0.21827291566383206112D+04,-.16032624215036328262D+04,
     *  0.27189460395887301756D+03,0.14938124228752306118D+04,
     *  -.30587785974653392161D+04,0.36342324685174129697D+04,
     *  -.25960074716370022397D+04,-.12862007902416296214D+03,
     *  0.37614160868398151450D+04,-.66250321579857786636D+04,
     *  0.64184857921430445032D+04,-.12644015613536503474D+04,
     *  -.78728032831082983068D+04,0.14894169091475538615D+05,
     *  -.10991419472447870017D+05,-.61235275727288084238D+04,
     *  0.24222367079716810557D+05,-.23261134072763513477D+05,
     *  -.31024689008096525513D+04,0.33464176192919793639D+05,
     *  -.34861341013849122519D+05,-.80608275052853156985D+03,
     *  0.40996514837262809759D+05,-.42609313245126908009D+05,
     *  -.23232184059632199529D+02,0.44620039500578638746D+05/
c 
        data cs24/
     *  0.12425242469463896796D+03,-.39731361424467956354D+03,
     *  0.64801363680521704856D+03,-.73251858311704129861D+03,
     *  0.51352317101844478183D+03,0.11157995412290074448D+03,
     *  -.11627647139447969222D+04,0.25203794149970996142D+04,
     *  -.38805819702416151118D+04,0.47565452181176943070D+04,
     *  -.45620678667906980961D+04,0.27991179249290245634D+04,
     *  0.66862542129515033916D+03,-.53332453520795868565D+04,
     *  0.98864793739974007855D+04,-.12349387563198029036D+05,
     *  0.10602249580928585301D+05,-.33023884045021229493D+04,
     *  -.89980295480644807370D+04,0.22846983653662998838D+05,
     *  -.31193415714077367991D+05,0.24159537008408988049D+05,
     *  0.52680907377438577327D+04,-.49227530578204305172D+05,
     *  0.76823482270230059627D+05,-.49353104062575775699D+05,
     *  -.38277401789913397421D+05,0.12721459461537315898D+06,
     *  -.12652890137367038695D+06,0.32522032836237723850D+04,
     *  0.15585533601738727618D+06,-.20009425427952309451D+06,
     *  0.61843800015949797672D+05,0.15245187437147725852D+06,
     *  -.24452577290223288207D+06,0.11757237774319487690D+06/
c 
        data cs25/
     *  -.70988030354808558350D+04,0.23828374385260426699D+05,
     *  -.43147498068273895141D+05,0.60178161027763194703D+05,
     *  -.70140369271644502973D+05,0.68661198283129192672D+05,
     *  -.52733453369439770607D+05,0.22036091823743492533D+05,
     *  0.19795542907511897439D+05,-.64578980456965650377D+05,
     *  0.10037498698525173531D+06,-.11411953886606863207D+06,
     *  0.95930515271051020696D+05,-.44023817443868634471D+05,
     *  -.31652844563053565490D+05,0.10922860020460870968D+06,
     *  -.15962473011965629340D+06,0.15594336818424145221D+06,
     *  -.85134411089354183230D+05,-.41967593721827858129D+05,
     *  0.18537248058208358757D+06,-.27793015660992509681D+06,
     *  0.24255470944757012194D+06,-.37892546889834462849D+05,
     *  -.27494445483317182322D+06,0.50428197331419509877D+06,
     *  -.42090282054267101918D+06,-.33829098336225825983D+05,
     *  0.59503927098210915297D+06,-.79458209754144932316D+06,
     *  0.35542043059241817115D+06,0.46997852538949050299D+06,
     *  -.10105422865581880704D+07,0.74384066485439166978D+06,
     *  0.17735234895595036575D+06,-.99103142128246451319D+06/
c 
        data cs26/
     *  0.37422966498527103958D+05,-.12757958095711320276D+06,
     *  0.23789790811697789974D+06,-.34756830269090851239D+06,
     *  0.43470825762387133779D+06,-.47571237825384107407D+06,
     *  0.44804212391243652463D+06,-.33585920196827819444D+06,
     *  0.13765962792018630309D+06,0.12616954339548852264D+06,
     *  -.40925395985122952934D+06,0.64253331273856475850D+06,
     *  -.74784404709398693969D+06,0.66153224548647597960D+06,
     *  -.36305845560435910813D+06,-.10126573356664479966D+06,
     *  0.61162033758202703832D+06,-.99481856899793543022D+06,
     *  0.10729092969432844154D+07,-.73068463288659912357D+06,
     *  -.14519363095905828613D+05,0.95806685541943423457D+06,
     *  -.16974781204822868708D+07,0.17303981495159368832D+07,
     *  -.72839097070529701432D+06,-.10792308781318713904D+07,
     *  0.27232579536770466657D+07,-.29038116426922567684D+07,
     *  0.98566568107119499520D+06,0.21298270565002439752D+07,
     *  -.42421510231109624781D+07,0.34295183970049656718D+07,
     *  0.15633548982466952077D+06,-.39954701222990010114D+07,
     *  0.50325876115841825753D+07,-.22614554604863803885D+07/
c 
        data cs27/
     *  0.28443828812355817832D+06,-.94884505492885811580D+06,
     *  0.16976747216594324101D+07,-.23220267126820485875D+07,
     *  0.26236633525890454143D+07,-.24349190954009345214D+07,
     *  0.16636711058891848955D+07,-.34811963019297140155D+06,
     *  -.12999715006046823797D+07,0.28925140979970280512D+07,
     *  -.39300253916251658519D+07,0.39437166511725559660D+07,
     *  -.26900289892298894753D+07,0.33207366318126988417D+06,
     *  0.24850313984040459367D+07,-.47486734711738246321D+07,
     *  0.54072818837034355075D+07,-.38371485748719873591D+07,
     *  0.27012257567847982841D+06,0.40433493065063306946D+07,
     *  -.71183662330545288192D+07,0.70164235438791031688D+07,
     *  -.29548146603919867729D+07,-.37338555446703981659D+07,
     *  0.95614891456366204298D+07,-.10416412913317826657D+08,
     *  0.43816993268797079597D+07,0.59141230742894415329D+07,
     *  -.13812577246265776838D+08,0.12839713796569454474D+08,
     *  -.21508706843673227835D+07,-.11443299803604001637D+08,
     *  0.17701382593711001766D+08,-.10998074997873411160D+08,
     *  -.44016550059605758603D+07,0.16974849412493908053D+08/
c 
        data cs28/
     *  -.12341981538446633070D+07,0.41910187280742471109D+07,
     *  -.77569955210004685695D+07,0.11200466311745457170D+08,
     *  -.13763078365356788917D+08,0.14658284717513912125D+08,
     *  -.13187484687662283628D+08,0.89474513284330197556D+07,
     *  -.20882653538685327570D+07,-.64637686527792270330D+07,
     *  0.14943930692510987515D+08,-.20977416074076840798D+08,
     *  0.22166804097845624487D+08,-.16971427977189725744D+08,
     *  0.56301200113636095026D+07,0.92738766643584866522D+07,
     *  -.23031324589171666783D+08,0.30051118488418356436D+08,
     *  -.25982121691951940446D+08,0.10082693324952662273D+08,
     *  0.13175206463092358582D+08,-.34394867616674404836D+08,
     *  0.42257160947209655395D+08,-.29074586957601316866D+08,
     *  -.30985580851043658161D+07,0.39805221035703798560D+08,
     *  -.59109463968636791807D+08,0.44376171575694691286D+08,
     *  0.17801446490265387312D+07,-.53664647032342176167D+08,
     *  0.75941423700959881750D+08,-.48544507989679404219D+08,
     *  -.14961960611969908723D+08,0.72602495268849597512D+08,
     *  -.82760611859088792880D+08,0.35981697557726329468D+08/
c 
        data cs29/
     *  -.59975552097288225608D+07,0.19949113843016380681D+08,
     *  -.35493446245331931778D+08,0.48101177214294770058D+08,
     *  -.53541958388317334156D+08,0.48377760421161010549D+08,
     *  -.30968286379233814869D+08,0.26275495577744430273D+07,
     *  0.31565200214962942361D+08,-.62938302398582672992D+08,
     *  0.80930235431646733490D+08,-.76407147086321421180D+08,
     *  0.45876132755310160653D+08,0.50066135989706695435D+07,
     *  -.60706039192114381198D+08,0.99474951408639761355D+08,
     *  -.10142490549941511975D+09,0.58612638061679170136D+08,
     *  0.17396377200769183980D+08,-.95573793430825961750D+08,
     *  0.13608259227211315436D+09,-.10972348262672930951D+09,
     *  0.18137536663167586486D+08,0.97773305375896259266D+08,
     *  -.17200162371394943149D+09,0.14978642753715077144D+09,
     *  -.28336557575968515641D+08,-.12761336200122602163D+09,
     *  0.21649286671511387043D+09,-.16690255990108351568D+09,
     *  -.31326360277293653129D+07,0.18483810890643581306D+09,
     *  -.24848097168658677997D+09,0.13975374562394469141D+09,
     *  0.71900909299803185999D+08,-.23752885851887535325D+09/
c 
        data cs30/
     8  0.24750736543183966078D+08,-.83858331707255957821D+08,
     9  0.15454908262942091641D+09,-.22165159066929613836D+09,
     *  0.26958653055725691638D+09,-.28257224635558834026D+09,
     1  0.24722220060985434986D+09,-.15691482268837697167D+09,
     2  0.16957997552036603484D+08,0.15135725837859098984D+09,
     3  -.31046029175766960106D+09,0.41248352418414845323D+09,
     4  -.41248079038484413168D+09,0.28675027379882711274D+09,
     5  -.50324372982065890340D+08,-.23508671494190890866D+09,
     6  0.47097617562048523157D+09,-.55287493764508956784D+09,
     7  0.41584467214255201556D+09,-.77706206080056888305D+08,
     8  -.34356420656233830995D+09,0.65874691759708951234D+09,
     9  -.68631452198016708828D+09,0.35627011315450790122D+09,
     *  0.21541903392802501547D+09,-.74333148128121944363D+09,
     1  0.90568409874974814228D+09,-.54391697520385380414D+09,
     2  -.19316976166592788390D+09,0.88480925083831376609D+09,
     3  -.10766485461001150849D+10,0.59155031066112552627D+09,
     4  0.30742543053254633671D+09,-.10511529330911003367D+10,
     5  0.11368003623161922944D+10,-.48515523831866775590D+09/
c 
        data cs31/
     8  -.26072706535612945326D+08,0.85688703576106805063D+08,
     9  -.14796324022895728072D+09,0.18772570462435088013D+09,
     *  -.18051769429595397265D+09,0.10960386174595955152D+09,
     1  0.22773692004621239875D+08,-.18575514220441113371D+09,
     2  0.32013926022018119582D+09,-.35555548482083362399D+09,
     3  0.24400662826919805682D+09,0.28923159292830572923D+07,
     4  -.29587811688322318657D+09,0.49107595211574729528D+09,
     5  -.45784212584926284053D+09,0.16463284489990302733D+09,
     6  0.26930636405198263887D+09,-.60637250825441603097D+09,
     7  0.61552160977352257259D+09,-.23185741169895478268D+09,
     8  -.35550131560875129747D+09,0.77697073150188625627D+09,
     9  -.71092575465280331267D+09,0.13355149916572023833D+09,
     *  0.60289564887226594685D+09,-.97079968867526152754D+09,
     1  0.64900430305557682391D+09,0.19054300247051185129D+09,
     2  -.94930507418856576529D+09,0.10289760137617743120D+10,
     3  -.31578506177237627732D+09,-.67931838323509832055D+09,
     4  0.11846661198254273874D+10,-.77914011533891846807D+09,
     5  -.24534295419374019808D+09,0.10921556783712240178D+10/
c 
        data cs32/
     8  0.21655350683197056369D+08,-.50842529392481330025D+08,
     9  0.29470776185549082826D+08,0.59986042711654682593D+08,
     *  -.17714100304130979870D+09,0.23713460485737805116D+09,
     1  -.16092553570476201748D+09,-.60691249809454595217D+08,
     2  0.33151688941078393641D+09,-.48082701420067410399D+09,
     3  0.36325400506110231969D+09,0.22965948593216253847D+08,
     4  -.48525774957930568617D+09,0.71879796632319572519D+09,
     5  -.50341978169769102888D+09,-.10371976218534444048D+09,
     6  0.73734064612104836403D+09,-.93082154971434676424D+09,
     7  0.46038670186150332257D+09,0.42156729913388710580D+09,
     8  -.10902228344138398246D+10,0.98601931146092117219D+09,
     9  -.89541112419377993599D+08,-.96975887236613825684D+09,
     *  0.13417228227856847779D+10,-.64644291446802520247D+09,
     1  -.63750651178105576168D+09,0.14815474091751195812D+10,
     2  -.11380593394467611767D+10,-.18220091717514697736D+09,
     3  0.14118304081927909663D+10,-.14956337312745078542D+10,
     4  0.31388485532597174797D+09,0.11643097358965653185D+10,
     5  -.16822212321871229922D+10,0.78184286745852119727D+09/
c 
        data cs33/
     8  0.14123965134047758978D+09,-.38422673004387102795D+09,
     9  0.43529337955630440317D+09,-.17438894383220017973D+09,
     *  -.29197810746138061431D+09,0.67625787801907259594D+09,
     1  -.68475076498498224465D+09,0.23291152441027559336D+09,
     2  0.44473705818058912756D+09,-.89986139677603094227D+09,
     3  0.76558306340675658508D+09,-.56905548290106187625D+08,
     4  -.76568684529276338001D+09,0.10801531810070877416D+10,
     5  -.57857273771266527515D+09,-.42654317292235460586D+09,
     6  0.11655179452398997431D+10,-.99692135335832834773D+09,
     7  -.21813735937321305180D+08,0.10919930173828441075D+10,
     8  -.12874654514554701909D+10,0.36767499314014215268D+09,
     9  0.93077523330424285924D+09,-.14671224366916640334D+10,
     *  0.70707836060336520406D+09,0.73079441659714720443D+09,
     1  -.15608905217565154453D+10,0.98660318302814878841D+09,
     2  0.51947854766378493080D+09,-.15890460228856870763D+10,
     3  0.12074749309677344439D+10,0.30905812956271812196D+09,
     4  -.15642915794015183807D+10,0.13745195544845011664D+10,
     5  0.10242893214365992778D+09,-.14925308634925610424D+10/
c 
        data cs34/
     8  0.74025932556557055203D+08,-.15153770833293090448D+09,
     9  0.24978865380489894099D+08,0.28760854656487578308D+09,
     *  -.52518022949003128748D+09,0.38148272311023562074D+09,
     1  0.19532609611387893179D+09,-.85454717906474855479D+09,
     2  0.10325697024487016255D+10,-.40799601945381117818D+09,
     3  -.72242851072945225854D+09,0.15258895313517372855D+10,
     4  -.12361565088066826866D+10,-.13653240576979118331D+09,
     5  0.16247147625202412567D+10,-.19613385462158360283D+10,
     6  0.66539068452947580700D+09,0.13682141900690711331D+10,
     7  -.24402412420185932695D+10,0.14745595859252298409D+10,
     8  0.89421746825490661750D+09,-.26697959379679397018D+10,
     9  0.21785319767444245813D+10,0.32747457749385575692D+09,
     *  -.27006623889934103216D+10,0.27362512775935597899D+10,
     1  -.24996215886896074077D+09,-.25873757362668195563D+10,
     2  0.31426926745359344699D+10,-.79334849140990613001D+09,
     3  -.23704491136440518448D+10,0.34057509814346112814D+10,
     4  -.12824170812366352638D+10,-.20741620744549976277D+10,
     5  0.35345699872289499099D+10,-.17101182307493854990D+10/
c 
        data cs35/
     8  -.87333811581545478468D+09,0.23682418905354603957D+10,
     9  -.26607534712031674081D+10,0.10214757904476817387D+10,
     *  0.18580871639593807678D+10,-.41749156789380624939D+10,
     1  0.41318131318562650498D+10,-.12766414065908039331D+10,
     2  -.28624585579452291509D+10,0.55052419403118163945D+10,
     3  -.44894662785425857060D+10,0.75109829649409165086D+08,
     4  0.48046695445988637018D+10,-.64193995444989147387D+10,
     5  0.31418410072112795846D+10,0.28831083502783908757D+10,
     6  -.69865689683617010308D+10,0.55940948508547168406D+10,
     7  0.56722616345546733205D+09,-.65982772083516614671D+10,
     8  0.72789007833394023137D+10,-.16654813384548154902D+10,
     9  -.56697383361298412097D+10,0.82940897274076035870D+10,
     *  -.36127619205858139454D+10,-.44870172629153190366D+10,
     1  0.87903688234629211305D+10,-.52221942536627998108D+10,
     2  -.32122756483527255507D+10,0.88942562352761849636D+10,
     3  -.65061042810505055780D+10,-.19214653745822233509D+10,
     4  0.86894582323772025804D+10,-.74959088268750497699D+10,
     5  -.63868813458970346992D+09,0.82195344482858311578D+10/
c 
        data cs36/
     8  -.93839655549269080473D+09,0.24310192656886065572D+10,
     9  -.23967578597043429794D+10,0.24431625157139749976D+09,
     *  0.27966960876371541378D+10,-.44250430004894142600D+10,
     1  0.29878470303501310189D+10,0.98693795623979810672D+09,
     2  -.47617147139069686060D+10,0.52460871776736632030D+10,
     3  -.15672956112892534463D+10,-.37620139975200601955D+10,
     4  0.64440926349851023288D+10,-.38953418957454068858D+10,
     5  -.21946433107574384748D+10,0.67937648794769112383D+10,
     6  -.56797785852989643082D+10,-.58586252075006668421D+09,
     7  0.66621251041676375272D+10,-.69464400183604929413D+10,
     8  0.82669713539939177416D+09,0.63235315609456886020D+10,
     9  -.78169645641638186321D+10,0.19767285574751295679D+10,
     *  0.59322573636319880556D+10,-.84053821329056332497D+10,
     1  0.28759648477206180986D+10,0.55554890385333943484D+10,
     2  -.87925459141333125398D+10,0.35655400152246269087D+10,
     3  0.52072150950534142109D+10,-.90284207462875158030D+10,
     4  0.40950374253363214891D+10,0.48711869667125008941D+10,
     5  -.91404239465312649069D+10,0.45144909798977614494D+10/
c 
c        return the exponents
c 
        n=36
c 
        do 1200 i=1,n
c 
        expons(i)=expons0(i)
        expons(2*n-i+1)=-expons(i)
 1200 continue
c 
c       return the nodes and weights
c 
c 
        do 1400 i=1,20
c 
        tt40(i)=tt0(i)
        ww40(i)=ww0(i)
c 
        tt40(41-i)=-tt40(i)
        ww40(41-i)=ww40(i)
 1400 continue
c 
c       return the coefficients
c 
        cd0=1
        cd1=-cd0*ima
        cd2=-cd1*ima
        cd3=-cd2*ima
c 
        call expons_cpy(cs0,coefs(1,1),cd0)
        call expons_cpy(cs1,coefs(1,2),cd1)
        call expons_cpy(cs2,coefs(1,3),cd2)
        call expons_cpy(cs3,coefs(1,4),cd3)
c 
        call expons_cpy(cs4,coefs(1,5),cd0)
        call expons_cpy(cs5,coefs(1,6),cd1)
        call expons_cpy(cs6,coefs(1,7),cd2)
        call expons_cpy(cs7,coefs(1,8),cd3)
c 
        call expons_cpy(cs8,coefs(1,9),cd0)
        call expons_cpy(cs9,coefs(1,10),cd1)
        call expons_cpy(cs10,coefs(1,11),cd2)
        call expons_cpy(cs11,coefs(1,12),cd3)
c 
        call expons_cpy(cs12,coefs(1,13),cd0)
        call expons_cpy(cs13,coefs(1,14),cd1)
        call expons_cpy(cs14,coefs(1,15),cd2)
        call expons_cpy(cs15,coefs(1,16),cd3)
c 
        call expons_cpy(cs16,coefs(1,17),cd0)
        call expons_cpy(cs17,coefs(1,18),cd1)
        call expons_cpy(cs18,coefs(1,19),cd2)
        call expons_cpy(cs19,coefs(1,20),cd3)
c 
        call expons_cpy(cs20,coefs(1,21),cd0)
        call expons_cpy(cs21,coefs(1,22),cd1)
        call expons_cpy(cs22,coefs(1,23),cd2)
        call expons_cpy(cs23,coefs(1,24),cd3)
c 
        call expons_cpy(cs24,coefs(1,25),cd0)
        call expons_cpy(cs25,coefs(1,26),cd1)
        call expons_cpy(cs26,coefs(1,27),cd2)
        call expons_cpy(cs27,coefs(1,28),cd3)
c 
        call expons_cpy(cs28,coefs(1,29),cd0)
        call expons_cpy(cs29,coefs(1,30),cd1)
        call expons_cpy(cs30,coefs(1,31),cd2)
        call expons_cpy(cs31,coefs(1,32),cd3)
c 
        call expons_cpy(cs32,coefs(1,33),cd0)
        call expons_cpy(cs33,coefs(1,34),cd1)
        call expons_cpy(cs34,coefs(1,35),cd2)
        call expons_cpy(cs35,coefs(1,36),cd3)
c 
        call expons_cpy(cs36,coefs(1,37),cd0)
c 
        return
c 
c 
c 
c 
        entry points_ret(expons72,tt40,ww40)
c 
c       This entry returns to the user the 40 prolate nodes tt40
c       and their corresponding weights ww40, discretizing to
c       16-digit (or so) accuracy band-limited functions with
c       c=20 on the interval [-1,1]. It also returns the 72
c       exponents expons72 (all between -40 and 40) to be used
c       for the approximation of functions (for example, via the
c       subroutine expoeval).
c 
c        return the exponents
c 
        n=36
c 
        do 2200 i=1,n
c 
        expons72(i)=expons0(i)
        expons72(2*n-i+1)=-expons72(i)
 2200 continue
c 
c       return the nodes and weights
c 
c 
        do 2400 i=1,20
c 
        tt40(i)=tt0(i)
        ww40(i)=ww0(i)
c 
        tt40(41-i)=-tt40(i)
        ww40(41-i)=ww40(i)
 2400 continue
c 
        return
        end
