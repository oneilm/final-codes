        implicit real *8 (a-h,o-z)
        dimension ab(2,10 000),xs(100 000),
     1      whts(100 000),
     2      adiff(100 000),endinter(1000),ts(1000),ws(1000),
     3      ainte(1000000)
c 
        dimension w(1000 000),w2(1000 000)
  
        complex *16 sigma(100000),ima,rlam,
     1      fsout(100 000),errs(100 000),diffs(100 000),
     2      fs2(100000),fs3(100 000),fs4(100 000),oper(100 000)
  
  
c 
        external fun3
c 
        data ima/(0.0d0,1.0d0)/
c 
        call prini(6,13)
c 
C 
        PRINT *, 'ENTER k'
        READ *,k
        CALL PRINf('k=*',k,1 )
c 
        rlam=10000+ima
        rlam=3+ima
  
cccc        rlam=1.0d-12
cccc        rlam=1
c 
c       construct the composite chebychev discretization
c       of the interval [-b,-a]
c 
        n=10
cccc        n=1000
        a=1
        b=2
c 
        a2=-a
        b2=-b
c 
        call constcoef_crea(k,n,a,b,
     1      ab,xs,whts,nn,w,keep,lused)
  
c 
        call prinf('after constcoef_crea, lused=*',lused,1)
        call prinf('after constcoef_crea, keep=*',keep,1)
        call prin2('after constcoef_crea, whts=*',whts,nn)
        call prin2('after constcoef_crea, xs=*',xs,nn)
c 
c       create the array of values of sigma, and integrate
c       the product sigma*exp(rlam*x)
c 
        do 2200 i=1,nn
c 
        sigma(i)=sigma3(xs(i))
cccc        sigma(i)=1
c 
 2200 continue
c 
c       apply the Lippman-Schwinger operator to sigma
c 
        call constcoef_eval(rlam,sigma,w,fsout,w2)
c 
        call prin2('after constcoef_eval, fsout=*',
     1      fsout,nn*2)
c 
c        substitute the result into the ODE
c 
        itype=2
        call legeinmt(k,ainte,adiff,ts,ws,endinter,
     1      itype,w2)
  
        call prin2('after legeinmt, ts=*',ts,k)
c 
        do 3200 i=1,nn
c 
        fs2(i)=fsout(i)
 3200 continue
c 
        call prin2('fs2 as created*',fs2,nn*2)
c 
c       differentiate it twice
c 
        do 3300 i=1,n
cccc        call constcoef_matvec(adiff,fs2,fs3,k)
cccc        call constcoef_matvec(adiff,fs3,fs4,k)
c 
        ii=(i-1)*k+1
        call constcoef_matvec(adiff,fs2(ii),fs3(ii),k)
        call constcoef_matvec(adiff,fs3(ii),fs4(ii),k)
c 
 3300 continue
c 
        do 3400 i=1,nn
c 
        fs4(i)=fs4(i)/(ab(2,1)-ab(1,1))**2 *4
 3400 continue
  
        call prin2('second derivative*',fs4,nn*2)
c 
c        evaluate the operator
c 
        do 3600 i=1,nn
c 
        oper(i)=fs4(i)-rlam**2*fsout(i)
c 
        diffs(i)=oper(i)-sigma(i)
 3600 continue
  
        call prin2('and oper=*',oper,nn*2)
        call prin2('while sigma=*',sigma,nn*2)
        call prin2('and diffs=*',diffs,nn*2)
  
c 
        call convert_test(k)
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine convert_test(k)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1000 000),
     1      xsin(100 000),xsout(100 000),
     3      ts(1000),ab(2,1000)
c 
        complex *16 fsin(100 000),fsout(100 000),
     1      fsout2(100 000),errs(100 000)
  
c 
        a=1
        b=2
        kin=6
        nin=30
c 
        kout=12
        nout=2
c 
cccc        kout=7
cccc        nout=10
c 
        call constcoef_convert_init(a,b,kin,nin,kout,nout,
     1      w,keep,lused)
c 
c        construct the test function at the nodes of both
c        discretizations
c 
        call constcoef_discr1(a,b,kin,nin,
     1      ab,xsin,nnin)
c 
        call constcoef_discr1(a,b,kout,nout,
     1      ab,xsout,nnout)
c 
        do 2200 i=1,kin*nin
c 
        call funtest(xsin(i),fsin(i))
 2200 continue
c 
        do 2400 i=1,kout*nout
c 
        call funtest(xsout(i),fsout(i))
 2400 continue
c 
c        . . . interpolate
c 
        call constcoef_convert_evalc(fsin,fsout2,w)
c 
        call prin2('and fsout2=*',fsout2,nout*kout)
        call prin2('and fsout=*',fsout,nout*kout)
c 
        do 2600 i=1,kout*nout
c 
        errs(i)=fsout2(i)-fsout(i)
 2600 continue
  
        call prin2('and errs=*',errs,kout*nout)
        call prinf('while keep=*',keep,1)
        call prinf('and lused=*',lused,1)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine funtest(x,f)
        implicit real *8 (a-h,o-z)
        save
        complex *16 f,ima
        data ima/(0.0d0,1.0d0)/
c 
        f=x**3+ima*sin(x)
        return
        end
c 
c 
c 
c 
c 
        function fun3(x,rlam,par2)
        implicit real *8 (a-h,o-z)
        save
        complex *16 fun3,rlam
c 
        fun3=sigma3(x)*exp(-rlam*x)
  
cccc        f=1
cccc        call prin2('f=*',f,1)
  
cccc        fun3=1
  
cccc        fun3=exp(-rlam*x)
  
  
  
  
        return
        end
  
c 
c 
c 
c 
c 
        function sigma3(x)
        implicit real *8 (a-h,o-z)
c 
cccc        sigma3=x**2
        save
        sigma3=x-x**2
cccc        sigma3=1
        return
        end
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c       This is the end of the debugging code, and the beginning of
c       the ODE solver code proper
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c       This file contains 6 user-callable subroutines (or,
c       rather, entries): constcoef_crea, constcoef_eval,
c       constcoef_convert_init, constcoef_convert_eval,
c       constcoef_convert_evalc, constcoef_discr0. Following
c       is a brief description of these subroutines.
c 
c   constcoef_crea - constructs the various types of data
c        to be used in the (numerical) evaluation of integrals
c        of the form
c 
c        -1/(2*rlam)* \int_a^b e^{-rlam*|x-t|} * sigma(t) dt,          (1)
c 
c        with rlam a complex number, [a,b] a user specified
c        interval, and sigma a user supplied function to be
c        tabulated at the nodes xs provided by this subroutine.
c 
c   constcoef_eval - uses various types of data constructed
c        via a preceding call to the entry constcoef_crea
c        to evaluate the expression (1) with x=x(1), x=x(2),
c        ... x=x(n*k) constructed by the entry constcoef_crea
c 
c   constcoef_convert_init - constructs various types of data
c        to be used (by the entries constcoef_convert_eval,
c        cconstcoef_convert_eval) to interpolate functions
c        from one composite Legendre discretization of the
c        interval [a,b] to another. Please note that all
c        subintervals within each of the discretizations are
c        of the same length.
c 
c   constcoef_convert_eval - interpolates a real function from
c        one composite Gaussian discretization of the interval
c        [a,b] to another. It uses data generated via a
c        preceding call to the entry constcoef_convert_init.
c 
c   constcoef_convert_evalc - interpolates a complex function
c        from one composite Gaussian discretization of the
c        interval [a,b] to another. It uses data generated via
c        a preceding call to the entry constcoef_convert_init.
c 
c   constcoef_discr0 - constructs a composite Gaussian discre-
c        tization of the interval [-1,1]
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        subroutine constcoef_crea(k,n,a,b,
     1      ab,xs,whts,nn,w,keep,lused)
        implicit real *8 (a-h,o-z)
        save
        dimension ab(2,1),xs(1),whts(1),w(1),w2(1)
        complex *16 rlam,fsout(1),sigma(1)
c 
c        This subroutine constructs the various types of data
c        to be used in the (numerical) evaluation of integrals
c        of the form
c 
c        -1/(2*rlam)* \int_a^b e^{-rlam*|x-t|} * sigma(t) dt,          (1)
c 
c        with rlam a complex number, [a,b] a user specified
c        interval, and sigma a user supplied function to be
c        tabulated at the nodes xs provided by this subroutine.
c 
c        Obviously, the evaluation of the integral (1) above
c        produces the solution of the differential equation
c 
c             phi''(x)-rlam**2 * phi(x) = sigma (x),                  (1a)
c 
c        satisfying the radiation boundary condtions
c 
c             phi'(a)=-rlam*phi(a),
c                                                                     (1b)
c             phi'(b)=rlam*phi(b).
c 
c 
c        Please note that this is the initialization entry
c        point of this subroutine; the actual evaluation of the
c        integrals of the form (1) is to be performed by the
c        entry constcoef_eval (see below). Thus, the
c        coefficient lambda is not supplied to this entry;
C        IN OTHER WORDS, THIS ENTRY IS NOT LAMBDA-SPECIFIC.
c 
c                      Input parameters:
c 
c  k - the number of Chebychev nodes on each of the intervals to
c        be constructed
c  n - the number of equispaced intervals into which the interval
c        [a,b] is to be subdivided
c  a,b - the ends of the interval on which the expression (1) is
c        to be evaluated; please note that a < b is a must, unless
c        the real part of rlam (see the entry constcoef_evalall
c        of this subroutine) is small.
c 
c                      Output parameters:
c 
c  ab - the ends of the n subintervals into which the interval
c        [a,b] has been subdivided: ab(1,i) is the left end of
c        the i-th subinterval, and a(2,i) is the right end of
c        the said subinterval.
c  xs - the nodes (n*k of them things) into which the interval
c        [a,b] has been discretized
c  whts - the quadrature weights corresponding to the nodes xs
c  nn - the total number of nodes in the array xs (always equal
c        to n*k)
c  w - array containing all of the data created by this entry
c        as the fodder for the entry constcoef_evalall (see)
c  keep - the number of real *8 elements in array w that should
c        not be changed between the call to this entry and the
c        subsequent call(s) to the entry constcoef_evalall
c  lused - the total amount of space in array w (in real *8
c        elements) actually used by this entry
c 
c                      Work array:
c 
c  w - must be ????? real *8 elements long
c 
c       . . . construct the data to be used for the evaluation
c             of the expressions of the form
c 
c 
c           e^{-rlam*x} * \int_a^x e^{rlam*t) * sigma(t) dt,          (2)
c 
        iw1=11
c 
        call constcoef_crea0(k,n,a,b,
     1      ab,xs,whts,nn,w(iw1),keep1,lused)
c 
c       construct the data to be used for the evaluation
c       of the expressions of the form
c 
c 
c           e^{rlam*x} * \int_a^x e^{-rlam*t) * sigma(t) dt,          (3)
c 
        iw2=iw1+keep1+2
c 
        call constcoef_rev(k,n,a,b,
     1      ab,xs,whts,nn,w(iw2),keep2,lused)
c 
        w(1)=iw1+0.1
        w(2)=iw2+0.1
        keep=iw2+keep2
c 
        return
c 
c 
c 
c 
        entry constcoef_eval(rlam,sigma,w,fsout,w2)
c 
c        This ENTRY uses various types of data constructed
c        via a preceding call to the entry constcoef_crea
c        of this subroutine to evaluate the expression
c 
c 
c        -1/(2*rlam)* \int_a^b e^{-rlam*|x-t|} * sigma(t) dt,          (4)
c 
c        with x=x(1), x=x(2), ... x=x(n*k) constructed by the
c        entry constcoef_crea of this subroutine, rlam
c        a complex number with positive real part, and sigma a
c        user-supplied array of complex values at the nodes x(1),
c        x(1),...,x(n*k) of whatever function he likes most.
c        The claim to fame of this subroutine the fact that the
c        scheme it uses is stable as long as the real part of
c        rlam is positive. It does not overflow for large real
c        part of rlam.
c 
c                          Input parameters:
c 
c  rlam - the complex parameter in (4) above
c  sigma - the function in (4), tabulated at nn nodes {xs(i)}
c        returned by the preceding call to the entry
c        constcoef_crea of this subroutine
c  w - the array containing the various types of data created
c        during a preceding call to constcoef_crea
c 
c                          Output parameters:
c 
c  fsout - the values of the integral (4) tabulated at the node xs
c 
c                          Work arrays:
c 
c  w2 - must be ????? real *8 eloements long
c 
c 
c       . . . evaluate the expression
c 
c           e^{-rlam*x} * \int_a^x e^{rlam*t) * sigma(t) dt
c 
        kk=w(15)
        nnn=w(16)
        ifsout2=1
        lfsout2=kk*nnn*2+4
        iww2=ifsout2+lfsout2
        iw1=w(1)
c 
        call constcoef_evalall(rlam,sigma,w(iw1),
     1      w2(ifsout2),w2(iww2))
c 
c       evaluate the expression
c 
c           e^{rlam*x} * \int_a^x e^{-rlam*t) * sigma(t) dt
c 
        iw2=w(2)
        call constcoef_rev0(rlam,sigma,w(iw2),fsout,w2(iww2))
c 
c        add them things up and scale them properly
c 
        call constcoef_add(fsout,w2(ifsout2),nnn*kk*2)
c 
        do 3200 i=1,nnn*kk
c 
        fsout(i)=-fsout(i)/(2*rlam)
 3200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine constcoef_add(x,y,n)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1)
c 
        do 1200 i=1,n
        x(i)=x(i)+y(i)
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine constcoef_rev(k,n,a,b,
     1      ab,xs,whts,nn,w,keep,lused)
        implicit real *8 (a-h,o-z)
        save
        dimension ab(2,1),xs(1),whts(1),w(1),w2(1)
        complex *16 rlam,fsout(1),sigma(1),cd
c 
c        This subroutine constructs various types of data to
c        be used for the evaluation of expressions of the form
c 
c           e^{rlam*x} * \int_a^x e^{-rlam*t) * sigma(t) dt,          (1)
c 
c        with x \in [a,b], rlam a complex number, and sigma
c        the user-supplied function. The claim to fame of
c        this subroutine is the fact that the scheme it uses
c        is stable as long as the real part of rlam is positive.
c        THIS SUBROUTINE DOES NOT BOMB FOR VERY LARGE RLAM.
c        The output from this subroutine is used as a fodder
c        for the entry constcoef_evalall of this subroutine.
c        This entry has no known uses as a stand-alone device.
c 
c        PLEASE NOTE THAT THIS SUBROUTINE IS MERELY AN INTEFACE.
C        ALL OF THE HONEST WORK IS DONE BY THE SUBROUTINE
C        CONSTCOEF_CREA0 (SEE).
c 
c       reverse what needs to be reversed
c 
        a2=-a
        b2=-b
c 
        call constcoef_crea0(k,n,b2,a2,
     1      ab,xs,whts,nn,w,keep,lused)
c 
c        reverse the arrays xs, whts
c 
        do 1400 i=1,nn/2
c 
        d=xs(i)
        xs(i)=-xs(nn-i+1)
        xs(nn-i+1)=-d
c 
        d=whts(i)
        whts(i)=whts(nn-i+1)
        whts(nn-i+1)=d
 1400 continue
c 
        return
c 
c 
c 
c 
        entry constcoef_rev0(rlam,sigma,w,fsout,w2)
c 
c        This ENTRY uses various types of data constructed
c        via a preceding call to the entry constcoef_crea0 of
c        this subroutine to evaluate the expression
c 
c           e^{rlam*x} * \int_a^x e^{-rlam*t) * sigma(t) dt,          (1)
c 
c        with x=x(1), x=x(2), ... x=x(n*k) constructed by
c        constcoef_rev, rlam a complex number with positive real
c        part, and sigma a user-supplied array of complex values
c        at the nodes x(1),x(1),...,x(n*k) of whatever function
c        he likes most. The claim to fame of this subroutine the
c        fact that the scheme it uses is stable as long as the
c        real part of rlam is positive.
c 
c        This entry uses data constructed by a preceding call to
c        the entry constcoef_crea0 of this subroutine. This entry
c        has no known uses as a stand-alone device.
c 
c        PLEASE NOTE THAT THIS SUBROUTINE IS MERELY AN INTEFACE.
C        ALL OF THE HONEST WORK IS DONE BY THE SUBROUTINE
C        CONSTCOEF_EVALALL (SEE).
c 
c       construct the indefinite integral of fs
c 
        kk=w(5)
        nnn=w(6)
        nnnn=kk*nnn
c 
        do 2200 i=1,nnnn/2
c 
        cd=sigma(i)
        sigma(i)=sigma(nnnn-i+1)
        sigma(nnnn-i+1)=cd
 2200 continue
c 
        call constcoef_evalall(rlam,sigma,w,fsout,w2)
c 
c       reverse the obtained values
c 
        kk=w(5)
        nnn=w(6)
        nnnn=kk*nnn
c 
        do 2400 i=1,nnnn/2
c 
        cd=fsout(i)
        fsout(i)=fsout(nnnn-i+1)
        fsout(nnnn-i+1)=cd
c 
        cd=sigma(i)
        sigma(i)=sigma(nnnn-i+1)
        sigma(nnnn-i+1)=cd
 2400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine constcoef_crea0(k,n,a,b,
     1      ab,xs,whts,nn,w,keep,lused)
        implicit real *8 (a-h,o-z)
        save
        dimension ab(2,1),xs(1),whts(1),w(1),w2(1)
        complex *16 rlam,fsout(1),sigma(1)
c 
c        This subroutine constructs various types of data to
c        be used for the evaluation of expressions of the form
c 
c           e^{-rlam*x} * \int_a^x e^{rlam*t) * sigma(t) dt,          (1)
c 
c        with x \in [a,b], rlam a complex number, and sigma
c        the user-supplied function. The claim to fame of
c        this subroutine is the fact that the scheme it uses
c        is stable as long as the real part of rlam is positive.
c        THIS SUBROUTINE DOES NOT BOMB FOR VERY LARGE RLAM.
c        The output from this subroutine is used as a fodder
c        for the entry constcoef_evalall of this subroutine.
c        This entry has no known uses as a stand-alone device.
c 
c        This subroutine starts with constructing a composite
c        Chebychev discretization of the used-supplied interval
c        [a,b], consisting of n equal subintervals, each of the
c        latter discretized into k Chebychev nodes. Then, the
c        integration machinery for the evaluation of (1) at all
c        n*k of the said nodes is built, and returned to the
c        user in the array w.
c 
c                      Input parameters:
c 
c  k - the number of Chebychev nodes on each of the intervals to
c        be constructed
c  n - the number of equispaced intervals into which the interval
c        [a,b] is to be subdivided
c  a,b - the ends of the interval on which the expression (1) is
c        to be evaluated; please note that a < b is a must, unless
c        the real part of rlam (see the entry constcoef_evalall of
c        this subroutine) is small.
c 
c                      Output parameters:
c 
c  ab - the ends of the n subintervals into which the interval
c        [a,b] has been subdivided: ab(1,i) is the left end of
c        the i-th subinterval, and a(2,i) is the right end of
c        the said subinterval.
c  xs - the nodes (n*k of them things) into which the interval
c        [a,b] has been discretized
c  whts - the quadrature weights corresponding to the nodes xs
c  nn - the total number of nodes in the array xs (always equal
c        to n*k)
c  w - array containing all of the data created by this entry
c        as the fodder for the entry constcoef_evalall (see)
c  keep - the number of real *8 elements in array w that should
c        not be changed between the call to this entry and the
c        subsequent call(s) to the entry constcoef_evalall
c  lused - the total amount of space in array w (in real *8
c        elements) actually used by this entry
c 
c        allocate memory
c 
        iab=11
        lab=n*2+4
c 
        iainte=iab+lab
        lainte=k**2+10
c 
        ixs=iainte+lainte
        lxs=n*k+10
c 
        iendinter=ixs+lxs
        lendinter=k+2
c 
        keep=iendinter+lendinter
c 
        iw2=iendinter+lendinter
        lw2=3*k**2+5*k+50
c 
        lused=iw2+lw2
c 
        call constcoef_discr(k,n,a,b,w(iab),w(ixs),whts,nn,
     1      w(iainte),w(iendinter),w(iw2) )
c 
c        return to the user the parameters that are believed
c        to be of use to him (or her, or it - we do not
c        discriminate)
c 
        call constcoef_move(w(iab),ab,k*2)
        call constcoef_move(w(ixs),xs,n*k)
c 
c        store the memory map in the beginning of array w
c 
        w(1)=iab+0.1
        w(2)=iainte+0.1
        w(3)=ixs+0.1
        w(4)=iendinter+0.1
        w(5)=k+0.1
        w(6)=n+0.1
c 
        return
c 
c 
c 
c 
        entry constcoef_evalall(rlam,sigma,w,fsout,w2)
c 
c        This ENTRY uses various types of data constructed
c        via a preceding call to the entry constcoef_crea0 of
c        this subroutine to evaluate the expression
c 
c           e^{-rlam*x} * \int_a^x e^{rlam*t) * sigma(t) dt,          (1)
c 
c        with x=x(1), x=x(2), ... x=x(n*k) constructed by
c        constcoef_crea0, rlam a complex number with positive
c        real part, and sigma a user-supplied array of complex
c        values at the nodes x(1),x(1),...,x(n*k) of whatever
c        function he likes most. The claim to fame of this
c        subroutine the fact that the scheme it uses is stable
c        as long as the real part of rlam is positive.
c 
c        This entry uses data constructed by a preceding call to
c        the entry constcoef_crea0 of this subroutine. This entry
c        has no known uses as a stand-alone device.
c 
c        . . . unpack the memory map
c 
        iab=w(1)
        iainte=w(2)
        ixs=w(3)
        iendinter=w(4)
        kk=w(5)
        nnn=w(6)
c 
        iadds=1
        ladds=nnn*2+4
c 
        ifs=iadds+ladds
        lfs=nnn*kk*2+10
c 
        iw2=ifs+lfs
c 
        call constcoef_integr(kk,nnn,w2(ifs),fsout,
     1      w(ixs),w(iainte),w(iendinter),w2(iadds),
     2      w2(iw2),sigma,rlam,w(iab))
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine constcoef_integr(k,m,fs,fsout,
     1      xs,ainte,endinter,adds,w,sigma,rlam,ab)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),ainte(1),endinter(1),
     1      w(1),ab(2,1)
c 
        complex *16 sigma(1),adds(1),fsout(1),fs(1),
     1      rlam,coef
c 
c        evaluate the product sigma * exp(rlam * t) at
c        the nodes of discretization
c 
        nn=m*k
        do 1400 i=1,m
        ii=(i-1)*k
        shift=ab(1,i)
c 
        do 1200 j=1,k
c 
        fs(ii+j)=exp(rlam*(xs(ii+j)-shift)) * sigma(ii+j)
 1200 continue
 1400 continue
  
cccc        call prin2('fs as created*',fs,nn)
  
cccc        stop
  
c 
c        one subinterval after another, construct the indefinite
c        integral
c 
        adds(1)=0
        do 2000 i=1,m
c 
        ii=(i-1)*k+1
        call constcoef_matvec(ainte,fs(ii),fsout(ii),k)
        call constcoef_scap(fsout(ii),endinter,k,adds(i+1))
 2000 continue
c 
c       splice together the subintervals, obtaining the indefinite
c       integral on the whole interval
c 
        do 2200 i=2,m
c 
        d=ab(1,i)-ab(1,i-1)
        coef=exp(-rlam*d)
        adds(i)=adds(i-1)*coef+adds(i)
 2200 continue
c 
  
cccc        call prin2('fsout before scaling*',fsout,nn)
  
cccc        call prin2('adds before scaling*',adds,m)
  
cccc        stop
  
  
        do 2600 i=1,m
c 
        ii=(i-1)*k
        shift=ab(1,i)
c 
        d=0
        if(i .ne. 1) d=ab(1,i)-ab(1,i-1)
c 
        coef=exp(-rlam*d)
  
cccc        call prinf('i=*',i,1)
cccc        call prin2('and coef=*',coef,2)
c 
        alpha=(ab(2,i)-ab(1,i))/2
  
  
        do 2400 j=1,k
c 
        fsout(ii+j)=fsout(ii+j)+adds(i)*coef
        fsout(ii+j)=fsout(ii+j)*alpha
        fsout(ii+j)=fsout(ii+j)*exp(-rlam*(xs(ii+j)-shift))
c 
 2400 continue
c 
 2600 continue
c 
  
cccc        call prin2('fsout after scaling*',fsout,nn)
  
cccc        stop
  
  
  
        return
        end
c 
c 
c 
c 
c 
        subroutine constcoef_discr(k,n,a,b,
     1      ab,xs,whts,nn,ainte,endinter,w)
        implicit real *8 (a-h,o-z)
        save
        dimension ab(2,1),xs(1),whts(1),ainte(1),
     1      endinter(1),w(1)
c 
        its=1
        lts=k+2
c 
        iws=its+lts
        lws=k+2
c 
        iw=iws+lws
c 
c       construct the composite chebychev discretization
c       of the interval [a,b]
c 
        call constcoef_discr0(a,b,k,n,
     1      ab,xs,whts,nn)
c 
c       construct the matrix of the indefinite integral
c 
        itype=1
        call legeinmt(k,ainte,adiff,w(its),w(iws),endinter,
     1      itype,w(iw))
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine constcoef_matvec(a,x,y,n)
        implicit real *8 (a-h,o-z)
        save
        complex *16 x(1),y(1),cd
        dimension a(n,n)
c 
        do 2000 i=1,n
        cd=0
        do 1400 j=1,n
        cd=cd+a(i,j)*x(j)
 1400 continue
        y(i)=cd
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine constcoef_scap(x,y,n,prod)
        implicit real *8 (a-h,o-z)
        save
        complex *16 x(1),prod
        dimension y(1)
c 
        prod=0
        do 1200 i=1,n
        prod=prod+x(i)*y(i)
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine constcoef_discr0(a,b,k,n,
     1      ab,xs,whts,nn)
        implicit real *8 (a-h,o-z)
        save
        dimension ab(2,1),xs(1),whts(1)
c 
c        This subroutine constructs a composite Legendre
c        discretization of the user-specified interval [a,b]
c 
c                  Input parameters:
c 
c  a,b - the ends of the interval on which the expression (1) is
c        to be evaluated; please note that a < b is a must, unless
c        the real part of rlam (see the entry constcoef_evalall
c        of this subroutine) is small.
c  k - the number of Chebychev nodes on each of the intervals to
c        be constructed
c  n - the number of equispaced intervals into which the interval
c        [a,b] is to be subdivided
c 
c                      Output parameters:
c 
c  ab - the ends of the n subintervals into which the interval
c        [a,b] has been subdivided: ab(1,i) is the left end of
c        the i-th subinterval, and a(2,i) is the right end of
c        the said subinterval.
c  xs - the nodes (n*k of them things) into which the interval
c        [a,b] has been discretized
c  whts - the quadrature weights corresponding to the nodes xs
c  nn - the total number of nodes in the array xs (always equal
c        to n*k)
c 
c        construct Chebychev nodes and weights on the
c        interval [-1,1]
c 
        itype=1
        call legeexps(itype,k,xs,u,v,whts)
c 
c        construct the nodes and weights of the composite
c        quadrature on the interval [a,b]
c 
        hab=(b-a)/n
c 
        do 1400 i=n,1,-1
c 
        ab(1,i)=(i-1)*hab+a
        ab(2,i)=ab(1,i)+hab
c 
        alpha=(ab(2,i)-ab(1,i))/2
        beta=(ab(2,i)+ab(1,i))/2
c 
        mmm=(i-1)*k
c 
        do 1200 j=1,k
c 
        xs(mmm+j)=alpha*xs(j)+beta
        whts(mmm+j)=whts(j)*alpha
 1200 continue
c 
 1400 continue
c 
        nn=n*k
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine constcoef_move(x,y,n)
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
        subroutine constcoef_convert_init(a,b,kin,nin,kout,nout,
     1      w,keep,lused)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1)
c 
c        This subroutine constructs various types of data to be used
c        (by the entries of this subroutine constcoef_convert_eval,
c        cconstcoef_convert_eval) to interpolate functions from one
c        composite Legendre discretization of the interval [a,b] to
c        another. Please note that all subintervals within each of
c        the discretizations are of the same length.
c 
c                   Input parameters:
c 
c  a,b - the ends of the interval on which the discretizations are
c        defined
c  kin - the number of Gaussian nodes on each of the nin
c        subintervals in the discretization FROM which we will be
c        interpolating,
c  nin - the number of subintervals in the discretization FROM which
c        we will be interpolating
c  kout - the number of Gaussian nodes on each of the nin
c        subintervals in the discretization TO which we will be
c        interpolating,
c  nout - the number of subintervals in the discretization TO which
c        we will be interpolating
c 
c                   Output parameters:
c 
c  w - the array in which the subroutine sored the various data
c        to be used by the entries constcoef_convert_eval,
c        cconstcoef_convert_eval. Must be at least ???? real *8
c        locations long.
c  keep - the number of real *8 elements in array w that have
c        to be unchanged between the call to this subroutine
c        and the subsequent calls to the entries constcoef_convert_eval,
c        cconstcoef_convert_eval.
c  lused - the amount of space (in real *8 locations) in array
c        w that was actually used by this subroutine
c 
c       . . . allocate memory
c 
        nn=kout*nout
        iinwhich=11
        linwhich=nn+1
c 
        icoefs=iinwhich+linwhich
        lcoefs=nn*kin+10
c 
        keep=icoefs+lcoefs
c 
        iabin=icoefs+lcoefs
        labin=nin*2+4
c 
        iabout=iabin+labin
        labout=nout*2+4
c 
        itsin=iabout+labout
        ltsin=nin+2
c 
        ixsin=itsin+ltsin
        lxsin=kin*nin+10
c 
        ixsout=ixsin+lxsin
        lxsout=kout*nout+10
c 
        iuin=ixsout+lxsout
        luin=kin**2+4
c 
        ivin=iuin+luin
        lvin=kin**2+4
c 
        lused=ivin+lvin
c 
c       construct the data to be used for the interpolation
c 
        call constcoef_convert_init0(a,b,kin,nin,kout,nout,
     1      w(iinwhich),w(icoefs),w(ixsin),w(ixsout),w(iabin),
     2      w(iuin),w(ivin),w(iabout),w(itsin))
c 
c       store various integer data in the beginning
c       of array w
c 
        w(1)=kin+0.1
        w(2)=nin+0.1
        w(3)=kout+0.1
        w(4)=nout+0.1
c 
        return
c 
c 
c 
c 
        entry constcoef_convert_eval(fsin,fsout,w)
c 
c        This entry interpolates a real function from one
c        composite Gaussian discretization of the interval
c        [a,b] to another. It uses data generated via a
c        preceding call to the entry constcoef_convert_init
c        of this subroutine.
c 
c                   Input parameters:
c 
c  fsin - the function to be interpolated (tabulated at the
c        nodes of the first discretization)
c  w - array containing data generated by the preceding
c        call to the entry constcoef_convert_init of this subroutine
c 
c        output parameters:
c 
c  fsout - the interpolated function
c 
c 
        kin2=w(1)
        nin2=w(2)
        kout2=w(3)
        nout2=w(4)
c 
        call constcoef_convert_eval0(kin2,nin2,fsin,
     1      w(iinwhich),w(icoefs),fsout,kout2,nout2)
c 
        return
c 
c 
c 
c 
        entry constcoef_convert_evalc(fsin,fsout,w)
c 
c        This entry interpolates a complex function from one
c        composite Gaussian discretization of the interval
c        [a,b] to another. It uses data generated via a
c        preceding call to the entry constcoef_convert_init
c        of this subroutine.
c 
c                   Input parameters:
c 
c  fsin - the function to be interpolated (tabulated at the
c        nodes of the first discretization)
c  w - array containing data generated by the preceding
c        call to the entry constcoef_convert_init of this subroutine
c 
c        output parameters:
c 
c  fsout - the interpolated function
c 
        kin2=w(1)
        nin2=w(2)
        kout2=w(3)
        nout2=w(4)
c 
        call cconstcoef_convert_eval0(kin2,nin2,fsin,
     1      w(iinwhich),w(icoefs),fsout,kout2,nout2)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine constcoef_convert_init0(a,b,kin,nin,
     1      kout,nout,inwhich,coefs,xsin,xsout,abin,
     2      uin,vin,about,tsin)
        implicit real *8 (a-h,o-z)
        save
        dimension abin(2,1),xsin(1),tsin(1),uin(1),vin(1),
     1      inwhich(1),about(2,1),xsout(1),
     2      coefs(kin,1)
c 
c       construct the input discretization of the interval [a,b]
c 
        call constcoef_discr1(a,b,kin,nin,
     1      abin,xsin,nnin)
c 
c       construct the output discretization of the interval [a,b]
c 
        call constcoef_discr1(a,b,kout,nout,
     1      about,xsout,nnout)
c 
c        for each output node, determine the subinterval of
c        the input discretization on which it lives
c 
        h=abin(2,1)-abin(1,1)
        do 1400 i=1,nnout
c 
        inode=(xsout(i)-a)/h
        inode=inode+1
c 
        inwhich(i)=inode
 1400 continue
c 
c        for each output node, construct the interpolation
c        formula on the appropriate subinterval of the
c        input discretization
c 
        h=abin(2,1)-abin(1,1)
        ifinit=1
        do 1800 i=1,nnout
c 
        inn=inwhich(i)
        alpha=2/h
        beta=1-abin(2,inn)*alpha
c 
        x=alpha*xsout(i)+beta
c 
        call levecin(kin,x,tsin,uin,vin,coefs(1,i),ifinit)
c 
        ifinit=0
c 
 1800 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine constcoef_convert_eval0(kin,nin,fsin,
     1      inwhich,coefs,fsout,kout,nout)
        implicit real *8 (a-h,o-z)
        save
        dimension inwhich(1),fsout(1),
     1      fsin(1),coefs(kin,1)
c 
c       interpolate the values of the function fsin from the
c       input grid to the output
c 
        do 2000 ipoint=1,kout*nout
c 
        ii=inwhich(ipoint)
c 
        i1=(ii-1)*kin
        fout=0
        do 1200 i=1,kin
c 
        fout=fout+coefs(i,ipoint)*fsin(i1+i)
 1200 continue
c 
        fsout(ipoint)=fout
 2000 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cconstcoef_convert_eval0(kin,nin,fsin,
     1      inwhich,coefs,fsout,kout,nout)
        implicit real *8 (a-h,o-z)
        save
        dimension inwhich(1),coefs(kin,1)
c 
        complex *16 fsout(1),fsin(1),fout
c 
c       interpolate the values of the function fsin from the
c       input grid to the output
c 
        do 2000 ipoint=1,kout*nout
c 
        ii=inwhich(ipoint)
c 
        i1=(ii-1)*kin
        fout=0
        do 1200 i=1,kin
c 
        fout=fout+coefs(i,ipoint)*fsin(i1+i)
 1200 continue
c 
        fsout(ipoint)=fout
 2000 continue
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine constcoef_discr1(a,b,k,n,
     1      ab,xs,nn)
        implicit real *8 (a-h,o-z)
        save
        dimension ab(2,1),xs(1)
c 
c        construct Legendre nodes and weights on the interval [-1,1]
c 
        itype=0
        call legeexps(itype,k,xs,u,v,ws)
c 
c        construct the nodes and weights of the composite
c        quadrature on the interval [a,b]
c 
        hab=(b-a)/n
c 
        do 1400 i=n,1,-1
c 
        ab(1,i)=(i-1)*hab+a
        ab(2,i)=ab(1,i)+hab
c 
        alpha=(ab(2,i)-ab(1,i))/2
        beta=(ab(2,i)+ab(1,i))/2
c 
        mmm=(i-1)*k
c 
        do 1200 j=1,k
c 
        xs(mmm+j)=alpha*xs(j)+beta
 1200 continue
c 
 1400 continue
c 
        nn=n*k
c 
        return
        end
