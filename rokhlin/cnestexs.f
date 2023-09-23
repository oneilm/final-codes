         implicit real *8 (a-h,o-z)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER k '
        READ *,k
        CALL PRINF('k=*',k,1 )
C 
        PRINT *, 'ENTER nsings '
        READ *,nsings
        CALL PRINF('nsings=*',nsings,1 )
  
  
cccc        call testex(k,nsings)
  
        call testexs(k,nsings)
  
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine testexs(k,nsings)
        implicit real *8 (a-h,o-z)
        save
        dimension par1(100),ab(20 000),
     1      w(100 000),xstest(10000),errsmax(10000)
c 
        complex *16 fs(4 000 000),
     1  vals0(1000),vals1(1000)
  
        external fun1
c 
c      construct all parameters
c 
        call fun1ini(nsings)
        a=0
        b=1
        eps=1.0d-12
        lab=10 000
c 
  
        lenw=100 000
  
        nfuns=nsings
  
        lenfs=4 000 000
        ifinit=1
        call cnestexs(ier,a,b,k,fun1,nfuns,eps,ab,lab,nn,
     1      w,lenw,ipar1,par2,fs,lenfs,ltotfs,ifinit,ltot,keep)
  
        call prinf('after cnestecxs, nn=*',nn,1)
c 
c 
c       test the obtained expansions
c 
        ipar1=1
        ntest=100
        h=(b-a)/ntest
        do 2200 i=1,ntest
c 
        xstest(i)=a+h/2+(i-1)*h
 2200 continue
c 
        call prin2('and xstest as created*',xstest,ntest)
  
  
        errmax=0
  
  
  
        do 2400 i=1,ntest
c 
        errsmax(i)=0
cccc        call prinf('in main program, i=*',i,1)
  
  
        call cnestem(ier,fs,ab,nn,k,nfuns,xstest(i),vals0)
  
  
cccc        call prin2('and vals0=*',vals0,nfuns*2)
  
  
        do 2300 j=1,nfuns
c 
        call fun1(xstest(i),j,par1,par2,vals1(j))
c 
        d=abs(vals1(j)-vals0(j))
        if(d .gt. errmax) errmax=d
        if(d .gt. errsmax(i)) errsmax(i)=d
  
  
 2300 continue
c 
cccc        call prin2('and vals1=*',vals1,nfuns*2)
cccc        call prin2('and errmax=*',errmax,1)
  
cccc         stop
  
  
 2400 continue
  
        call prin2('and errsmax=*',errsmax,ntest)
  
c 
c       now, time the evaluation
c 
  
  
        tt1=clotatim()
  
        do 3600 ijk=1,100
        do 3400 i=1,ntest
c 
        call cnestem(ier,fs,ab,nn,k,nfuns,xstest(i),vals0)
  
 3400 continue
  
 3600 continue
  
        tt2=clotatim()
  
        call prin2('and CPU time for all functions*',tt2-tt1,1)
  
        return
        end
  
c 
c 
c 
c 
c 
        subroutine testex(k,nsings)
        implicit real *8 (a-h,o-z)
        save
        dimension par1(100),ab(20 000),
     1      w(100 000),xstest(10000)
c 
        complex *16 f(100 000),val,der,val2,vals(10 000),
     1   vals2(10 000),errs(10 000),valp,valm,derp,derm,der2
  
        external fun1
c 
c      construct all parameters
c 
        call fun1ini(nsings)
        a=0
        b=1
        eps=1.0d-12
        lab=10 000
        ifinit=1
c 
  
        lenw=100 000
        lenf=100 000
  
         ipar1=11
  
        call cnestex(ier,a,b,k,fun1,eps,ab,lab,nn,w,lenw,
     1      ipar1,par2,par3,f,lenf,ifinit,ltot,keep)
c 
        call prinf('after cnestex, ltot=*',ltot,1)
        call prinf('after cnestex, keep=*',keep,1)
  
        ifinit=0
  
        call cnestex(ier,a,b,k,fun1,eps,ab,lab,nn,w,lenw,
     1      ipar1,par2,par3,f,lenf,ifinit,ltot,keep)
  
  
         call prin2('after cnestex, ab=*',ab,nn*2)
c 
c       test the evaluation of the function
c 
        ntest=100
        h=(b-a)/ntest
        do 2200 i=1,ntest
c 
        xstest(i)=a+h/2+(i-1)*h
 2200 continue
c 
        call prin2('and xstest as created*',xstest,ntest)
  
  
cccc          stop
  
  
  
  
  
        xtest=0.7
cccc        xtest=0.2
        call cneste(ier,f,ab,nn,k,xtest,val,der)
  
        call prin2('xtest=*',xtest,1)
        call prin2('val=*',val,2)
  
        call fun1(xtest,ipar1,par1,par2,val2)
  
        call prin2('and val2=*',val2,2)
        call prin2('and val-val2=*',val-val2,2)
  
        do 2400 i=1,ntest
c 
cccc        call cneste(ier,f,ab,nn,k,xstest(i),vals(i),der)
        call cneste2(ier,f,ab,nn,k,xstest(i),vals(i))
c 
        call fun1(xstest(i),ipar1,par1,par2,vals2(i))
  
        errs(i)=vals(i)-vals2(i)
 2400 continue
  
        call prin2('and xstest=*',xstest,ntest)
        call prin2('and vals=*',vals,ntest*2)
        call prin2('and errs=*',errs,ntest*2)
c 
c       test the evaluation of the derivative
c 
        h=0.00001
        call cneste(ier,f,ab,nn,k,xtest,val,der)
  
        call cneste(ier,f,ab,nn,k,xtest+h,valp,derp)
        call cneste(ier,f,ab,nn,k,xtest-h,valm,derm)
  
        der2=(valp-valm)/h/2
  
        call prin2('testing derivative, der=*',der,2)
        call prin2('testing derivative, der2=*',der2,2)
        call prin2('testing derivative, der-der2=*',der-der2,2)
  
  
        return
        end
c 
c 
c 
c 
c 
        subroutine fun1(x,i,par1,par2,f)
        implicit real *8 (a-h,o-z)
        save
        dimension par1(1),par2(1),sings(10000)
c 
        complex *16 f,ima
c 
        data ima/(0.0d0,1.0d0)/
c 
        f=x**sings(i)+ima*(1-x)**sings(i)
cccc        f=x**sings(i)+(1-x)**sings(i)
  
cccc        call prin2('x=*',x,1)
cccc        call prin2('f=*',f,2)
c 
cccc        f=sin(x)
  
        return
c 
c 
c 
c 
        entry fun1ini(nsings7)
c 
        nsings=nsings7
        done=1
        h=done/nsings
        do 2200 j=1,nsings
        sings(j)=j*h
 2200 continue
c 
        call prin2('in fun1ini, sings=*',sings,nsings)
c 
        return
        end
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c       This is the end of the debugging code and the beginning of the
c       interpolation code proper
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c       This file contains numerical machinery for the construction and
c       evaluation of nested Legendre expansions of complex functions
c       of one real variable. There are 6 user-callable subroutines in
c       this file: cnestex, cneste, cneste2, cnestexs, cnesteva, cnestev2.
c       Following is a brief description of these subroutines.
c 
c   cnestex - finds a subdivision of the interval [a,b], such that on
c       every subinterval, the user-provided function given by the
c       subroutine fun (see description below) is interpolated to
c       precision eps by a k-point Legendre interpolation formula; then
c       it evaluates on each subinterval the coefficients of the Legendre
c       series approximation of f. The resulting is gauaranteed to accuracy
c       eps, and the output of this subroutine can be used by the
c       subroutines cneste, cneste2 (see) to evaluate the interpolated
c       function at arbitrary points on the interval [a,b].
c 
c   cneste - evaluates at the point x \in [a,b] the complex
c        function (together with its derivative) specified by its nested
c        Legendre expansion, whose coefficients are stored in the input
c        array coefs. Note that this subroutine uses output of the
c        subroutine cnestex (see) as its input; it has no known uses as
c        a stand-alone device
c 
c   cneste2 - evaluates at the point x \in [a,b] the complex function
c        specified by its nested Legendre expansion, whose coefficients
c        are provideds in the input array coefs. Please note that this
c        subroutine uses output of the subroutine cnestex (see) as its
c        input; it has no known uses as a stand-alone device
c 
c   cnestexs - finds a subdivision of the interval [a,b], such that on
c       every subinterval, the user-provided functions (nfuns of them
c       things) given by the subroutine fun is interpolated to precision
c       eps by a k-point Legendre interpolation formula; then it evaluates
c       on each subinterval the coefficients of the Legendre series
c       approximation of f. The resulting approximation is gauaranteed
c       to accuracy eps, and the output of this subroutine can be used
c       by the subroutines cnesteva, cnestev2, cnestem (see) to evaluate
c       the interpolated functions at arbitrary points on the interval
c       [a,b].
c 
c   cnesteva - evaluates at the point x \in [a,b] the function specified
c       by its nested Legendre expansion, together with the derivativce
c       of the said function. The coefficients of the Legendre expansion
c       are supplied in the input array coefs. Normally, the latter has
c       been produced by a prior call to the subroutine cnestex (see)
c 
c   cnestev2 - evaluates at the point x \in [a,b] the function specified
c       by its nested Legendre expansion, whose coefficients are supplied
c       in the input array coefs. Normally, the latter has been produced
c       by a prior call to the subroutine cnestex (see)
c 
c   cnestem - evaluates at the point x \in [a,b] the nfuns functions
c        specified by their nested Legendre expansions, whose coefficients
c        are supplied in the input array coefs. Normally, the latter has
c        been produced by a prior call to the subroutine cnestexs (see)
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine cnestem(ier,coefs,ab,nn,k,nfuns,x,vals)
        implicit real *8 (a-h,o-z)
        save
        dimension ab(2,1),pols(300)
        complex *16  coefs(nn*k,nfuns),vals(1)
c 
c        this subroutine evaluates at the point x \in [a,b] the
c        nfuns functions specified by their nested Legendre
c        expansions, whose coefficients are supplied in the input
c        array coefs. Normally, the latter has been produced by a
c        prior call to the subroutine cnestexs (see)
c 
c                   Input parameters:
c 
c  coefs - the nested Legendre expansions of the functions stored
c      in the array fs.
c  ab - the array of ends of subintervals created. For each i=1,2,
c      ...,nn, ab(1,i) is the left end of the i-th subinterval;
c      ab(2,i) is the right end of the i-th subinterval
c  nn - the number of subintervals created on the interval [a,b];
c     also the number of entries in the array ab (see above)
c  k - the number of Gaussian nodes on each of nn subintervals in
c      the nested subdivision of the interval [a,b].
c  nfuns - the number of functions whose expansions are stored in
c      the array coefs; also the number of functions whose values
c      will be returned in array vals
c  x - the point where the functions are to be evaluated; must be
c      on the interval [a,b]
c 
c                    Output parameters:
c 
c  ier - error return code. ier = 0 means successful execution
c                           ier > 0 means that the point x is
c                           outside the interval [a,b]
c  vals - values of the nfuns functions evaluated at the point x
c 
c       . . . find the subinterval in which the point x lives
c 
        ier=0
c 
        call findinte(ier,x,ab,nn,intnum)
        if(ier .ne. 0) return
c 
c       the point x lives in the interval number intnum.
c       evaluate the expansion at this point
c 
        u=2/(ab(2,intnum)-ab(1,intnum))
        v=1-ab(2,intnum)*u
c 
        t=u*x+v
c 
        jj=(intnum-1)*k+1
c 
c        evaluate all Legendre polynomials of t
c 
        call legepls2(t,k,pols)
c 
c       evaluate all expansions
c 
        do 1200 i=1,nfuns
c 
        vals(i)=0
 1200 continue
c 
        jj0=(intnum-1)*k
        do 1600 j=1,k
c 
        jj=jj0+j
        do 1400 i=1,nfuns
c 
        vals(i)=vals(i)+coefs(jj,i)*pols(j)
 1400 continue
 1600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cnestev2(ier,coefs,ab,nn,k,x,i,val)
        implicit real *8 (a-h,o-z)
        save
        dimension ab(2,1),pjcoefs1(300),pjcoefs2(300)
        complex *16  coefs(nn*k,1),val,der
c 
c        This subroutine evaluates at the point x \in [a,b] one of the
c        functions (specifically, the function number i) specified by
c        their nested Legendre expansions, whose coefficients are
c        supplied in the input array coefs. Normally, the latter has
c        been produced by a prior call to the subroutine cnestexs (see)
c 
c                   Input parameters:
c 
c  coefs - the nested Legendre expansions of the functions stored in the
c      array fs.
c  ab - the array of ends of subintervals created. For each i=1,2,...,nn,
c     ab(1,i) is the left end of the i-th subinterval; ab(2,i) is
c     the right end of the i-th subinterval
c  nn - the number of subintervals created on the interval [a,b];
c     also the number of entries in the array ab (see above)
c  k - the number of Gaussian nodes on each of nn subintervals in the
c      nested subdivision of the interval [a,b].
c  x - the point where the function is to be evaluated; must be on
c      the interval [a,b]
c  i - the sequence number of the function to be evaluated; must be
c      between 1 and ncols.
c 
c                    Output parameters:
c 
c  ier - error return code. ier = 0 means successful execution
c                           ier > 0 means that the point x is
c                           outside the interval [a,b]
c  val - the value of the i-th function at the point x
c 
c       . . . find the subinterval in which the point x lives
c 
        ier=0
c 
        call findinte(ier,x,ab,nn,intnum)
        if(ier .ne. 0) return
c 
c       the point x lives in the interval number intnum.
c       evaluate the expansion at this point
c 
        u=2/(ab(2,intnum)-ab(1,intnum))
        v=1-ab(2,intnum)*u
c 
        t=u*x+v
c 
        jj=(intnum-1)*k+1
c 
        ninit=290
        if(ifcalled .eq. 1) ninit=0
        ifcalled=1
c 
        call legecva2(t,VAL,coefs(jj,i),k-1,
     1    pjcoefs1,pjcoefs2,ninit)
c 
        der=der*u
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cnesteva(ier,coefs,ab,nn,k,x,i,val,der)
        implicit real *8 (a-h,o-z)
        save
        dimension ab(2,1),pjcoefs1(300),pjcoefs2(300)
        complex *16  coefs(nn*k,1),val,der
c 
c        This subroutine calculates the value and the derivative at the
c        point x \in [a,b] of one of the functions (specifically, the
c        function number i) specified by their nested Legendre expansions,
c        The coefficients of the latter are supplied in the input array
c        coefs. Normally, the latter has been produced by a prior call to
c        the subroutine cnestexs (see)
c 
c                   Input parameters:
c 
c  coefs - the nested Legendre expansions of the functions stored in the
c      array fs.
c  ab - the array of ends of subintervals created. For each i=1,2,...,nn,
c     ab(1,i) is the left end of the i-th subinterval; ab(2,i) is
c     the right end of the i-th subinterval
c  nn - the number of subintervals created on the interval [a,b];
c     also the number of entries in the array ab (see above)
c  k - the number of Gaussian nodes on each of nn subintervals in the
c      nested subdivision of the interval [a,b].
c  x - the point where the function is to be evaluated; must be on
c      the interval [a,b]
c  i - the sequence number of the function to be evaluated; must be
c      between 1 and ncols.
c 
c                    Output parameters:
c 
c  ier - error return code. ier = 0 means successful execution
c                           ier > 0 means that the point x is
c                           outside the interval [a,b]
c  val - the value of the i-th function at the point x
c  der - the derivative of the i-th function at the point x
c 
c       . . . find the subinterval in which the point x lives
c 
        call findinte(ier,x,ab,nn,intnum)
        if(ier .ne. 0) return
c 
c       the point x lives in the interval number intnum.
c       evaluate the expansion at this point
c 
        u=2/(ab(2,intnum)-ab(1,intnum))
        v=1-ab(2,intnum)*u
c 
        t=u*x+v
c 
        jj=(intnum-1)*k+1
c 
        ninit=290
        if(ifcalled .eq. 1) ninit=0
        ifcalled=1
c 
        call legecFD2(t,VAL,der,coefs(jj,i),k-1,
     1    pjcoefs1,pjcoefs2,ninit)
c 
        der=der*u
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cnestexs(ier,a,b,k,fun,nfuns,eps,ab,lab,nn,
     1      w,lenw,par1,par2,fs,lenfs,ltotfs,ifinit,ltot,keep)
        implicit real *8 (a-h,o-z)
        save
        dimension ab(2,1),w(1),par1(1),par2(1)
c 
        real *8 fs(1)
c 
        external fun
c 
c       this subroutine finds a subdivision of the interval [a,b],
c       such that on every subinterval, the user-provided functions
c       (nfuns of them things) given by the subroutine fun (see
c       description below) is interpolated to precision eps by a
c       k-point Legendre interpolation formula; then it evaluates
c       on each subinterval the coefficients of the Legendre series
c       approximation of f. The resulting approximation is gauaranteed
c       to accuracy eps. The output of this subroutine can be used by
c       the subroutines cnesteva, cnestev2, cnestem (see) to evaluate
c       the interpolated functions at arbitrary points on the interval
c       [a,b].
c 
c                          Input parameters:
c 
c  a,b - the ends of the interval to be subdivided
c  k - the number of Gaussian nodes on each subinterval
c  fun - the user-supplied subroutine defining the function to be
c      interpolated. The calling sequence of fun must be:
c 
c      fun(x,par1,par2,f),
c 
c      where x is the input argument where the function is to be
c      evaluated, apar1, par2, par3 are whatever parameters the
c      subroutine fun needs, and f (output) is the value of the
c      function to be interpolated at the point x.
c  nfuns - the number of functions to be interpolated
c  eps - the accuracy to which the legendre expansion must represent
c     the user-supplied function for the subroutine to decide that
c     the representation is accurate
c  lab - the length of the user-provided array ab (see below)
c  w - the array containing the data to be used by the subroutine,
c     created by a preceding call to it. It must be at
c     least 16*k+8*k**2+30 + 20 000 real *8 locations long. Please
c     note that this is an input parameter only if ifinit (see below)
c     has been set to 0; otherwise, it is an output parameter
c  lenw - the length of the array w (in real *8 location; should not
c     be less than 16*k+8*k**2+30 real *8 locations
c  par1,par2 - whatever the user-supplied subroutine fun needs
c  lenfs - the length of the user-supplied array fs; must be large
c  ifinit - the index telling the subroutine whether it should construct
c     the parameter w (see above), or view it as an input parameter.
c     In the latter case, this  must have been created by a prior call
c     to this subroutine with the same parameter k.
c     ifinit=1 means that the parameter w will be created.
c     ifinit=0 means that the parameter w will be viewed
c     as an input parameter.
c 
c                     Output parameters:
c 
c  ier - error return code.
c         ier=0 means successful execution
c         ier=4 means that the space allocated in the array ab is
c               insufficient
c         ier=16  means that the space allocated in the array fs is
c               insufficient
c  ab - the array of ends of subintervals created. For each i=1,2,...,nn,
c     ab(1,i) is the left end of the i-th subinterval; ab(2,i) is
c     the right end of the i-th subinterval
c  nn - the number of subintervals created on the interval [a,b];
c     also the number of entries in the array ab (see above)
c  w - the array containing the data used by the subroutine,
c     It is only an output parameter if the user has set the
c     parameter ifinit (see above) to 1. Otherwise, it is an input
c     parameter, and must have been created  by a prior call to
c     this subroutine. Note that it should not be altered between
c     consecutive calls to this subroutine; the piece used by this
c     subroutine is 16*k+8*k**2+30 real *8 locations long.
c  fs - coefficients of all nn nested Legendre expansions; can be used
c     by the subroutines cnesteva, cnestev2 (see) for the evaluation of
c     the functions fun at arbitrary points on the interval [a,b]. The
c     number of elements of fs (viewed as real *8 elements) used by
c     this subroutine is nn*k*nfuns*2+nn*k+100; these elements should
c     not be changed between the call to this subroutine and subsequent
c     calls to the subroutines cnesteva, cnestev, cnestem (see).
c  ltot - the length of the part of the array w (in real *8 words)
c     actually used by the subroutine
c  keep - the length of the part of the array w (in real *8 words)
c     that should not be changed between the call to this subroutine and
c     the subsequent calls to the subroutines cneste, cneste2 (see)
c 
c 
c        . . . construct a subdivision of the interval [a,b], such that
c              on every subinterval, the user-provided function fun is
c              interpolated to to precision eps by a k-point Legendre
c              expansion
c 
         ier=0
c 
         iw=1
         lw=16*k+8*k**2+30
c 
         it=iw+lw
         lt=k+2
c 
         iu=it+lt
         lu=k**2+10
c 
         iv=iu+lu
         lv=k**2+10
c 
         iww=iv+lv
         lww=k*2+4
c 
         ltot=iww+lww
         keep=iv
c 
         if(ltot .le. lenw) goto 1200
c 
         ier=4
         return
 1200 continue
c 
         call cmaninco(jer,a,b,k,fun,nfuns,par1,par2,
     1      eps,ab,lab,nn,fs,lenfs)
c 
         if(jer .ne. 0) then
            ier=4
            return
         endif
c 
         n=k*nn
c 
         ifs=1
         lfs=n*nfuns*2+10
c 
         ix=ifs+lfs
         lx=n+2
c 
         ltotfs=n*nfuns*2+n+100
c 
         call prinf('lenfs=*',lenfs,1)
         call prinf('and ltotfs=*',ltotfs,1)
c 
         if(ltotfs .gt. lenfs) then
             ier=16
             return
         endif
c 
         iwhts=ifs
c 
c       construct all nodes and weights on all subintervals
c 
        call cnodewht(ab,nn,k,fs(ix),fs(iwhts) )
c 
c        construct the Legendre expansions of the function on all
c        subintervals
c 
        itype=2
        if(ifinit .eq. 1)
     1      call legeexps(itype,k,w(it),w(iu),w(iv),w(iww))
c 
        call callexpa(fs(ix),fun,par1,par2,nfuns,k,nn,fs(ifs),
     1      w(iu),w(iww) )
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cnestex(ier,a,b,k,fun,eps,ab,lab,nn,
     1      w,lenw,par1,par2,par3,f,lenf,ifinit,ltot,keep)
        implicit real *8 (a-h,o-z)
        save
        dimension ab(2,1),w(1),par1(1),
     1      par2(1),par3(1),f(1)
c 
        external fun
c 
c       This subroutine finds a subdivision of the interval [a,b],
c       such that on every subinterval, the user-provided function
c       given by the subroutine fun (see description below) is
c       interpolated to precision eps by a k-point Legendre interpolation
c       formula; then it evaluates on each subinterval the coefficients
c       of the Legendre series approximation of f. The resulting
c       is gauaranteed to accuracy eps, and the output of this subroutine
c       can be used by the subroutines cneste, cneste2 (see) to evaluate
c       the function fun at arbitrary points on the interval [a,b].
c 
c                          input parameters:
c 
c  a,b - the ends of the interval to be subdivided
c  k - the number of Gaussian nodes on each subinterval
c  fun - the user-supplied subroutine defining the function to be
c      interpolated. The calling sequence of fun must be:
c 
c      fun(x,par1,par2,par3,f),
c 
c      where x is the input argument where the function is to be
c      evaluated, apar1, par2, par3 are whatever parameters the
c      subroutine fun needs, and f (output) is the value of the
c      function to be interpolated at the point x.
c  eps - the accuracy to which the legendre expansion must represent
c     the user-supplied function for the subroutine to decide that
c     the representation is accurate
c  lab - the length of the user-provided array ab (see below)
c  w - the array containing the data to be used by the subroutine,
c     created by a preceding call to it. It must be at
c     least 10*k+8*k**2+30 + 20 000 real *8 locations long.
c  par1,par2,par3 - whatever the user-supplied subroutine fun needs
c  ifinit - the index telling the subroutine whether it should construct
c     the parameter w (see above), or view it as an input parameter.
c     In the latter case, this  must have been created by a prior call
c     to this subroutine with the same parameter k.
c     ifinit=1 means that the parameter w will be created.
c     ifinit=0 means that the parameter w will be viewed
c     as an input parameter.
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
c  x - the nodes of the nested discretization (k*nn of them things)
c  whts - the weights corresponding to the nodes x
c  w - the array containing the data used by the subroutine,
c     It is only an output parameter if the user has set the
c     parameter ifinit (see above) to 1. Otherwise, it is an input
c     parameter, and must have been created  by a prior call to
c     this subroutine. Note that it should not be altered between
c     consecutive calls to this subroutine; the piece
c     used by this subroutine is 10*k+8*k**2+30 real *8 locations long.
c  f - coefficients of all nn nested Legendre expansions; can be used
c     by the subroutines cneste, cneste2 (see) for the evaluation of the
c     function fun at arbitrary points on the interval [a,b]
c  ltot - the length of the part of the array w (in real *8 words)
c     actually used by the subroutine
c  keep - the length of the part of the array w (in real *8 words)
c     that should not be changed between the call to this subroutine and
c     the subsequent calls to the subroutines cneste, cneste2 (see)
c 
c        . . . construct a subdivision of the interval [a,b], such that
c              on every subinterval, the user-provided function fun is
c              interpolated to to precision eps by a k-point Legendre
c              expansion
c 
         ier=0
c 
         iw=1
         lw=16*k+8*k**2+30
c 
         it=iw+lw
         lt=k+2
c 
         iu=it+lt
         lu=k**2+10
c 
         iv=iu+lu
         lv=k**2+10
c 
         iww=iv+lv
         lww=k*2+4
c 
         ltot=iww+lww
         keep=iv
c 
         if(ltot .le. lenw) goto 1200
c 
         ier=4
         return
 1200 continue
c 
         call callinco(ier,a,b,k,fun,par1,par2,par3,eps,ab,
     1       lab,nn,labused,w(iw),ifinit)
c 
         n=k*nn
c 
         ifs=1
         lfs=n*2+10
c 
         ix=ifs+lfs
         lx=n+2
c 
         ltotfs=n*2+n+100
c 
         call prinf('lenf=*',lenf,1)
         call prinf('and ltotfs=*',ltotfs,1)
c 
         if(ltotfs .gt. lenf) then
             ier=16
             return
         endif
c 
         iwhts=ifs
c 
c       construct all nodes and weights on all subintervals
c 
        call cnodewht(ab,nn,k,f(ix),f(iwhts))
c 
c        construct the Legendre expansions of the function on all
c        subintervals
c 
        itype=2
        if(ifinit .eq. 1)
     1      call legeexps(itype,k,w(it),w(iu),w(iv),w(iww))
c 
        do 1400 i=1,n
        call fun(f(ix+i-1),par1,par2,par3,f(2*i-1) )
 1400 continue
c 
        do 2200 i=1,nn
c 
        ii=(i-1)*k+1
        call callqrma(w(iu),f(2*ii-1),w(iww),k)
c 
        do 1800 j=1,k*2
        f(2*ii+j-2)=w(iww+j-1)
 1800 continue
c 
 2200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cneste(ier,coefs,ab,nn,k,x,val,der)
        implicit real *8 (a-h,o-z)
        save
        dimension ab(2,1),pjcoefs1(300),pjcoefs2(300)
c 
        complex *16 coefs(1),val,der
c 
        data ifcalled/0/
c 
c        This subroutine evaluates at the point x \in [a,b] the complex
c        function (together with its derivative) specified by its nested
c        Legendre expansion, whose coefficients are supplied in the input
c        array coefs. Please note that this subroutine uses output of the
c        subroutine cnestex (see) as its input; it has no known uses as
c        a stand-alone device.
c 
c                   Input parameters:
c 
c  coefs - the nested Legendre expansions of the functions stored in the
c      array fs.
c  ab - the array of ends of subintervals created. For each i=1,2,...,nn,
c     ab(1,i) is the left end of the i-th subinterval; ab(2,i) is
c     the right end of the i-th subinterval
c  nn - the number of subintervals created on the interval [a,b];
c     also the number of entries in the array ab (see above)
c  k - the number of Gaussian nodes on each of nn subintervals in the
c      nested subdivision of the interval [a,b].
c  x - the point where the function is to be evaluated; must be on
c      the interval [a,b]
c 
c                    Output parameters:
c 
c  ier - error return code; ier=0 means successful conclusion
c                           ier=4 means that the user-supplied foint x
c                           is outside the interval [a,b]
c  val - the value of the function fun at the point x
c  der - the derivative of the function fun at the point x
c 
c       . . . find the subinterval in which the point x lives
c 
        ier=0
        call findinte(ier,x,ab,nn,intnum)
        if(ier .ne. 0) return
c 
c       the point x lives in the interval number intnum.
c       evaluate the expansion at this point
c 
        u=2/(ab(2,intnum)-ab(1,intnum))
        v=1-ab(2,intnum)*u
c 
        t=u*x+v
c 
        jj=(intnum-1)*k+1
c 
        ninit=290
        if(ifcalled .eq. 1) ninit=0
        ifcalled=1
c 
        call legecFD2(t,VAL,der,coefs(jj),k-1,
     1    pjcoefs1,pjcoefs2,ninit)
c 
        der=der*u
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cneste2(ier,coefs,ab,nn,k,x,val)
        implicit real *8 (a-h,o-z)
        save
        dimension ab(2,1),pjcoefs1(300),pjcoefs2(300)
c 
        complex *16 coefs(1),val
c 
        data ifcalled/0/
c 
c        This subroutine evaluates at the point x \in [a,b] the complex
c        function specified by its nested Legendre expansion, whose
c        coefficients are supplied in the input array coefs. Please
c        note that this subroutine uses output of the subroutine cnestex
c        (see) as its input; it has no known uses as a stand-alone device.
c 
c                   Input parameters:
c 
c  coefs - the nested Legendre expansions of the functions stored in the
c      array fs.
c  ab - the array of ends of subintervals created. For each i=1,2,...,nn,
c     ab(1,i) is the left end of the i-th subinterval; ab(2,i) is
c     the right end of the i-th subinterval
c  nn - the number of subintervals created on the interval [a,b];
c     also the number of entries in the array ab (see above)
c  k - the number of Gaussian nodes on each of nn subintervals in the
c      nested subdivision of the interval [a,b].
c  x - the point where the function is to be evaluated; must be on
c      the interval [a,b]
c 
c                    Output parameters:
c 
c  ier - error return code; ier=0 means successful conclusion
c                           ier=4 means that the user-supplied foint x
c                           is outside the interval [a,b]
c  val - the value of the function fun at the point x
c 
c       . . . find the subinterval in which the point x lives
c 
        ier=0
c 
        call findinte(ier,x,ab,nn,intnum)
        if(ier .ne. 0) return
c 
c       the point x lives in the interval number intnum.
c       evaluate the expansion at this point
c 
        u=2/(ab(2,intnum)-ab(1,intnum))
        v=1-ab(2,intnum)*u
c 
        t=u*x+v
c 
        jj=(intnum-1)*k+1
c 
        ninit=290
        if(ifcalled .eq. 1) ninit=0
        ifcalled=1
c 
        call legecva2(t,VAL,coefs(jj),k-1,
     1    pjcoefs1,pjcoefs2,ninit)
c 
        der=der*u
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine callexpa(xs,fun,par1,par2,nfuns,k,nn,fs,u,ww)
        implicit real *8 (a-h,o-z)
        save
        dimension par1(1),par2(1),xs(1),ww(1),u(k,k)
        complex *16 fs(nn*k,nfuns)
c 
        external fun
c 
c        construct the Legendre expansions of all functions on all
c        subintervals
c 
        n=nn*k
        do 1600 i=1,nfuns
  
        call prinf('in callexpa, i=*',i,1)
  
        do 1400 j=1,n
        call fun(xs(j),i,par1,par1,fs(j,i) )
 1400 continue
 1600 continue
c 
        do 2200 i=1,nfuns
c 
        do 1800 j=1,nn
c 
        jj=(j-1)*k+1
        call callqrma(u,fs(jj,i),ww,k)
c 
        call arrmove(ww,fs(jj,i),k*2)
c 
 1800 continue
 2200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cmaninco(ier,a,b,k,fun,nfuns,par1,par2,
     1    eps,ab,lab,nn,w,lenw)
        save
        dimension par1(1),ab(2,1),w(1),par2(1)
        implicit real *8 (a-h,o-z)
c 
        external fun
c 
c       this subroutine finds a subdivision of the interval [a,b],
c       such that on every subinterval, every one of the
c       user-provided collection of functions given by the function
c       fun (see description below) is integrated to precision eps
c       by a k-point Gaussian quadrature.
c 
c                          input parameters:
c 
c  a,b - the ends of the interval to be subdivided
c  k - the number of Gaussian nodes on each subinterval
c  fun - the user-supplied function to be integrated
c      The calling sequence of fun must be:
c 
c      fun(x,i,par1,par2,f),
c 
c      where x is the input argument where the function is to be
c      evaluated, i the sequence number of function, f (output)
c      is the value of the i-th function at the point x, and par
c      is an input parameter, carrying whatever information fun
c      might need.
c 
c  par1,par2 - whatever the user-supplied subroutine fun needs
c  eps - the accuracy to which the legendre expansion must represent
c     the user-supplied function for the subroutine to decide that
c     the representation is accurate
c  lab - the length of the user-provided array ab (see below)
c  w - the array containing the data to be used by the subroutine,
c     created by a preceding call to it. It must be at
c     least 10*k+8*k**2+30 + 20 000 real *8 locations long.
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
c  labused - the total amount of space used by the subroutine in the
c     array ab; note that normally, it  is somewhat greater than
c     the expected nn*2; The difference is not large, but should
c     be taken into account when allocating memory.
c  w - the array containing the data used by the subroutine,
c     It is only an output parameter if the user has set the
c     parameter ifinit (see above) to 1. Otherwise, it is an input
c     parameter, and must have been created  by a prior call to
c     this subroutine. Note that it should not be altered between
c     consecutive calls to this subroutine; the piece
c     used by this subroutine is 10*k+8*k**2+30 real *8 locations long.
c 
c 
c        one after another, construct the subdivisions
c        of the interval [a,b] corresponding to the
c        user-specified functions to be integrated, and merge
c        all these intervals
c 
c        . . . initialize the merging process
c 
        ier=0
        ipar1=1
        ifinit=1
c 
        iw=1
        lw=16*k+8*k**2+30
c 
        iab1=iw+lw
        lab1=lenw-iab1
c 
         call callinco(ier,a,b,k,fun,ipar1,par1,par2,eps,
     1    w(iab1),lab1,nn1,labused,w(iw),ifinit)
c 
        if(ier .ne. 0) return
c 
c        . . . merge
c 
        do 2000 i=2,nfuns
c 
        ipar1=i
        ifinit=0
c 
        iab2=iab1+nn1*2+2
        lab2=lenw-iab2
c 
         call callinco(ier,a,b,k,fun,ipar1,par1,par2,eps,
     1    w(iab2),lab2,nn2,labused,w(iw),ifinit)
c 
        if(ier .ne. 0) return
c 
        iab3=iab2+nn2*2+2
        lab3=(nn1+nn2)*2+4
c 
        iww=iab3+lab3
        lleft=lenw-iww
        if(lleft .gt. (nn1+nn2)*2+4) goto 1600
        ier=4
        return
 1600 continue
c 
        call  cmgrstru(w(iab1),nn1,w(iab2),nn2,
     1      w(iab3),nn3,eps,w(iww))
c 
        call arrmove(w(iab3),w(iab1),nn3*2)
        nn1=nn3
c 
 2000 continue
c 
        call arrmove(w(iab3),ab,nn3*2)
        nn=nn3
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine arrmove(x,y,n)
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
        subroutine cmgrstru(ab1,nn1,ab2,nn2,ab3,nn3,eps,w)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),ab1(1),ab2(1),ab3(1)
c 
c        this subroutine merges the user-provided subdivisions to
c        get one single subdivision that is as dense as either of
c        the original ones at any point on the (common) interval
c        of their definition
c 
c                        input parameters:
c 
c  ab1 - the left and right end of subintervals in the
c        first subdivision
c  nn1 - the number of subintervals in the array ab1
c  ab2 - the left and right end of subintervals in the
c        second subdivision
c  nn2 - the number of subintervals in the array ab2
c  eps - tells the subroutine how short a subinterval must
c        be to be considered of length zero
c 
c                        input parameters:
c 
c  ab3 - the left and right end of subintervals in the
c        merged subdivision
c  nn3 - the number of subintervals in the array ab3
c 
c                        work arrays:
c 
c  w - must be at least (nn1+nn2)*2+2 real *8 elements long
c 
c 
c        . . . merge the arrays ab1,ab2 (both left and right ends)
c              into one long array with non-decreasing elements
c 
        call ctwomrge(ab1,nn1*2,ab2,nn2*2,w)
c 
c        now, scan the obtained merged array and discard repeating
c        elements
c 
        nnn=nn1*2+nn2*2
        ii=1
        ab3(ii)=w(1)
        do 1400 i=2,nnn
        d=w(i)-w(i-1)
        if(d .lt. eps) goto 1400
        ii=ii+1
        ab3(ii)=w(i)
 1400 continue
c 
c       now, duplicate all interior elements in array ab3
c 
        w(1)=ab3(1)
        iii=1
        do 1600 i=2,ii-1
        iii=iii+1
        w(iii)=ab3(i)
        iii=iii+1
        w(iii)=ab3(i)
 1600 continue
        w(iii+1)=ab3(ii)
        iii=iii+1
c 
c       copy the final array to where it belongs
c 
        nn3=iii/2
        do 1800 i=1,iii
        ab3(i)=w(i)
 1800 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine ctwomrge(a,na,b,nb,c)
        implicit real *8 (a-h,o-z)
        save
        real *8  a(1),b(1),c(1)
c 
c       merge the sorted arrays ia, ib obtaining the sorted array ic
c 
        ia=1
        ib=1
        ic=1
c 
        do 1600 i=1,na+nb
        ic=i
c 
        if(a(ia) .gt. b(ib) ) goto 1200
        c(i)=a(ia)
        ia=ia+1
        goto 1400
c 
 1200 continue
        c(i)=b(ib)
        ib=ib+1
 1400 continue
c 
        if(ib .gt. nb) goto 1800
        if(ia .gt. na) goto 2200
 1600 continue
c 
 1800 continue
c 
c       we have run out of elemnets in b. finish off the array a
c 
        do 2000 i=1,na-ia+1
        c(ic+i)=a(ia+i-1)
 2000 continue
        return
c 
 2200 continue
c 
c       we have run out of elements in a. finish off the array b
c 
        do 2400 i=1,nb-ib+1
        c(ic+i)=b(ib+i-1)
 2400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine cnodewht(ab,nn,k,x,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension ab(2,1),x(1),whts(1),t(2000),w(2000)
c 
c 
        data kold/-10/
c 
c       construct the nodes and weights of the k-point Gaussian
c       quadrature on the interval [-1,1]
c 
        itype=1
c 
        if(k .ne. kold) call legeexps(itype,k,t,u,v,w)
c 
        kold=k
c 
c       one subinterval after another, construct nodes and weights
c       of the quadrature on the big interval
c 
        do 2000 i=1,nn
        u=(ab(2,i)-ab(1,i))/2
        v=(ab(2,i)+ab(1,i))/2
        ik=(i-1)*k
c 
        do 1200 j=1,k
        x(ik+j)=t(j)*u+v
        whts(ik+j)=w(j)*u
 1200 continue
 2000 continue
c 
        return
        end
c 
c 
c 
c 
c 
         subroutine callinco(ier,a,b,k,fun,ipar1,par1,
     1    par2,eps,ab,lab,nn,labused,w,ifinit)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),par1(1),par2(1),ab(2,1)
c 
c       this subroutine finds a subdivision of the interval [a,b],
c       such that on every subinterval, the user-provided function
c       fun is integrated to precision eps by a k-point Gaussian
c       quadrature. Actually, this is a memory management routine;
c       all the work is done by the routine callinte (see).
c 
c                          input parameters:
c 
c  a,b - the ends of the interval to be subdivided
c  k - the number of Gaussian nodes on each subinterval
c  fun - the user-supplied function to be integrated
c      The calling sequence of funget must be:
c 
c      fun(x,i,par1,par2,f),
c 
c      where x is the input argument where the function is to be
c      evaluated, i is the sequence  number of the function to be
c      evaluated, f  (output) is the function, and par1,par2,apr3
c      are input parameters, carrying whatever information fun
c      might need.
c 
c  par1,par2 - whatever the user-supplied subroutine fun needs
c  eps - the accuracy to which the legendre expansion must represent
c     the user-supplied function for the subroutine to decide that
c     the representation is accurate
c  lab - the length of the user-provided array ab (see below)
c  w - the array containing the data to be used by the subroutine,
c     created by a preceding call to it. It is only an input parameter
c     if the user has set the parameter ifinit (see below) to 0.
c     Otherwise, it is an output parameter, and must be at
c     least 10*k+8*k**2+30 real *8 locations long.
c  ifinit - the index telling the subroutine whether it should construct
c     the parameter w (see above), or view it as an input parameter.
c     In the latter case, this  must have been created by a prior call
c     to this subroutine with the same parameter k.
c     ifinit=1 means that the parameter w will be created.
c     ifinit=0 means that the parameter w will be viewed
c     as an input parameter.
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
c  labused - the total amount of space used by the subroutine in the
c     array ab; note that normally, it  is somewhat greater than
c     the expected nn*2; The difference is not large, but should
c     be taken into account when allocating memory.
c  w - the array containing the data used by the subroutine,
c     It is only an output parameter if the user has set the
c     parameter ifinit (see above) to 1. Otherwise, it is an input
c     parameter, and must have been created  by a prior call to
c     this subroutine. Note that it should not be altered between
c     consecutive calls to this subroutine; the piece
c     used by this subroutine is 16*k+8*k**2+30 real *8 locations long.
c 
c       . . . allocate memory for the subroutine callinte
c 
        ier=0
c 
        it=1
        lt=k*2+2
c 
        ix=it+lt
        lx=k*2+2
c 
        if=ix+lx
        lf=k*2+2
        lf=lf*2
c 
        icoefs=if+lf
        lcoefs=k*2+2
        lcoefs=lcoefs*2
c 
        iwhts=icoefs+lcoefs
        lwhts=k*2+1
c 
        iu=iwhts+lwhts
        lu=k**2*4 +10
c 
        iv=iu+lu
        lv=k**2*4 +10
c 
c       recursively subdivide the interval [a,b] into subintervals
c       such that on each of them, the k-point Gaussian rule
c       will obtain the accuracy epc
c 
         call callinte(ier,a,b,fun,ipar1,par1,par2,k,ab,lab,
     1    eps,w(it),w(iu),w(iv),w(iwhts),w(ix),w(if),
     2       w(icoefs),nn,labused,ifinit)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cbubble(ab,nn)
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
         subroutine callinte(ier,a,b,fun,ipar1,par1,par2,
     1    k,ab,lab,eps,t,u,v,whts,x,f,coefs,nn,lused,ifinit7)
        implicit real *8 (a-h,o-z)
        save
        dimension t(1),u(2*k,2*k),v(2*k,2*k),whts(1),x(1),
     1      par1(1),par2(1),ab(2,1)
        equivalence (rea,intint)
c 
        complex *16 f(1),coefs(1)
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
c      The calling sequence of funget must be:
c 
c      funget(x,par1,par2,f),
c 
c      where x is the input argument where the function is to be
c      evaluated, f  (output) is the function, and par1,par2 are
c      both input parameters, carrying whatever information funget
c      might need.
c 
c  par1,par2 - whatever the user-supplied subroutine funget needs
c  k - the number of Gaussian nodes on each subinterval
c  lab - the length of the user-provided array ab (see below)
c  eps - the accuracy to which the legendre expansion must represent
c     the user-supplied function for the subroutine to decide that
c     the representation is accurate
c  t - the Gaussian nodes on the interval [-1,1]. This is only an
c     input parameter is ifinit7 (see below) has been set to 0;
c     otherwise, it is an output parameter.
c  u - the 2*k * 2*k  matrix converting the  values at of a polynomial
c     of order 2*k-1 at 2*k legendre nodes into the coefficients of its
c     legendre expansion. This is only an input parameter is ifinit7
c     (see below) has been set to 0; otherwise, it is an output parameter.
c  v - the 2*k * 2*k matrix converting the coefficients
c     of an 2*k-term legendre expansion into its values at 2*k legendre
c     nodes (note that v is the inverse of u). This is only an input
c     parameter is ifinit7 (see below) has been set to 0; otherwise, it
c     is an output parameter.
c  whts - the Gaussian quadrature weights (2*k of them).  This is only
c     an input parameter is ifinit7 (see below) has been set to 0;
c     otherwise, its an output parameter.
c  ifinit7 - the index telling the subroutine whether it should construct
c     the parameters t,u,v,whts (by calling the subroutine legeexps),
c     or view them as input parameters. In the latter case, these
c     must have been created by a prior call to this subroutine (or to
c     subroutine legeexps, for an adventurous and creative user).
c     ifinit7=1 means that the parameters t,u,v,whts will be created.
c     ifinit7=0 means that the parameters t,u,v,whts will be viewed
c     as input parameters.
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
c  t,u,v,whts - are only output is the parameter ifinit (see above)
c     has been set by the user to 1.
c  x - the Gaussian nodes on the interval [a,b]
c  f - the values of the user-specified function fun at
c     the nodes x
c  coefs - the coefficients of the Legendre expansion of the
c     function fun on the interval [a,b]
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
        ifinit=ifinit7
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
        call cifsplit(ier,ab(1,nn),ab(2,nn),k,fun,ipar1,par1,
     1    par2,eps/d,t,u,v,whts,x,ifinit,f,coefs)
c 
        ifinit=0
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
        call cbubble(ab,nn)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine cifsplit(ier,a,b,k,fun,ipar1,par1,
     1    par2,eps,t,u,v,whts,x,ifinit,f,coefs)
        implicit real *8 (a-h,o-z)
        save
        dimension t(1),u(2*k,2*k),v(2*k,2*k),whts(1),x(1),
     1      par1(1),par2(1)
c 
        complex *16 f(1),coefs(1)
c 
c       this subroutine determines whether a k-point legendre
c       expansion of the user-supplied function f on the interval
c       [a,b] represents it with accuracy eps.
c 
c                          input parameters:
c 
c  a,b - the ends of the interval on which the approximation is checked
c  k - the order of the approximation for which the accuracy is checked
c  funget - the user-supplied function for which the accuracy is checked
c      The calling sequence of funget must be:
c 
c      fun(x,i,par1,par2,f),
c 
c      where x is the input argument where the function is to be
c      evaluated, f  (output) is the function, and par1,par2 are
c      both input parameters, carrying whatever information funget
c      might need.
c 
c  par1,par2 - whatever the user-supplied subroutine funget needs
c  eps - the accuracy to which the legendre expansion must represent
c     the user-supplied function for the subroutine to decide that
c     the representation is accurate
c  t - the Gaussian nodes on the interval [-1,1]. This is only an
c     input parameter is ifinit (see below) has been set to 0;
c     otherwise, it is an output parameter.
c  u - the 2*k * 2*k  matrix converting the  values at of a polynomial
c     of order 2*k-1 at 2*k legendre nodes into the coefficients of its
c     legendre expansion. This is only an input parameter is ifinit
c     (see below) has been set to 0; otherwise, it is an output parameter.
c  v - the 2*k * 2*k matrix converting the coefficients
c     of an 2*k-term legendre expansion into its values at 2*k legendre
c     nodes (note that v is the inverse of u). This is only an input
c     parameter is ifinit (see below) has been set to 0; otherwise, it
c     is an output parameter.
c  whts - the Gaussian quadrature weights (2*k of them).  This is only
c     an input parameter is ifinit (see below) has been set to 0;
c     otherwise, its an output parameter.
c  ifinit - the index telling the subroutine whether it sholud construct
c     the parameters t,u,v,whts (by calling the subroutine legeexps),
c     or view them as input parameters. In the latter case, these
c     must have been created by a prior call to this subroutine (or to
c     subroutine legeexps, for an adventurous and creative user).
c     ifinit=1 means that the parameters t,u,v,whts will be created.
c     ifinit=0 means that the parameters t,u,v,whts will be viewed
c     as input parameters.
c 
c                     output parameters:
c 
c  t,u,v,whts - are only output is the parameter ifinit (see above)
c     has been set by the user to 1.
c  x - the Gaussian nodes on the interval [a,b]
c  f - the values of the user-specified function funget at
c     the nodes x
c  coefs - the coefficients of the Legendre expansion of the
c     function f on the interval [a,b]
c 
c       . . . if the need be, construct the gaussian nodes on the
c             interval [-1,1], and the matrix conevrting the values
c             of the function at these nodes into coefficients of
c             its legendre expansion
c 
        if(ifinit .eq. 0) goto 1200
c 
        itype=2
        call legeexps(itype,k*2,t,u,v,whts)
c 
 1200 continue
c 
c        discretize the interval and evaluate the
c        user-supplied function at its nodes
c 
        alpha=(b-a)/2
        beta=(b+a)/2
c 
        do 1400 i=1,k*2
        x(i)=alpha*t(i)+beta
        call fun(x(i),ipar1,par1,par2,f(i) )
 1400 continue
c 
c       construct the legendre expansion of the function on
c       the interval
c 
        call calqrma2(u,f,coefs,k*2,k+1)
c 
c       scan the coefficients and see if the last k
c       of them are small
c 
        eps22=eps**2
c 
        ier=0
        do 1800 i=1,k
        d=coefs(k+i)*conjg(coefs(k+i))
        if(d .gt. eps22) ier=4
 1800 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine calqrma2(a,x,y,n,n0)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n)
        complex *16 x(1),y(1),d
c 
        do 2000 i=n0,n
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
c 
c 
c 
        subroutine callqrma(a,x,y,n)
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
  
c 
c 
c 
c 
c 
        subroutine findinte(ier,x,ab,nn,intnum)
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
