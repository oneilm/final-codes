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
        call testex(k)
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine testex(k)
        implicit real *8 (a-h,o-z)
        save
        dimension par1(100),
     1      w(100 000),xstest(10000),w2(100000)
c 
        real *8 val,der,val2,vals(10 000),
     1   vals2(10 000),errs(10 000),valp,valm,derp,derm,der2
  
        external fun1
c 
c      construct all parameters
c 
        a=0
        b=1
        eps=1.0d-12
        lab=10 000
        ifinit=1
c 
  
        lenw2=100 000
        lenf=100 000
  
         ipar1=11
  
c 
        call rnestex_cheb(ier,a,b,k,fun1,
     1      ipar1,par2,par3,eps,nn,w,ifinit,w2,lenw2,keepw2)
c 
        ifinit=0
  
        call rnestex_cheb(ier,a,b,k,fun1,
     1      ipar1,par2,par3,eps,nn,w,ifinit,w2,lenw2,keepw2)
c 
c       test the evaluation of the function
c 
        iab=w2(1)
        ifs=w2(4)
  
        ntest=100
        h=(b-a)/ntest
        do 2200 i=1,ntest
c 
        xstest(i)=a+h/2+(i-1)*h
 2200 continue
c 
        call prin2('and xstest as created*',xstest,ntest)
c 
        xtest=0.7
        call rneste_cheb(ier,w2,xtest,val,der)
c 
        call prin2('xtest=*',xtest,1)
        call prin2('val=*',val,1)
  
        call fun1(xtest,ipar1,par1,par2,val2)
  
        call prin2('and val2=*',val2,1)
        call prin2('and val-val2=*',val-val2,1)
c
        do 2400 i=1,ntest
c 
        call rneste2_cheb(ier,w2,xstest(i),vals(i))
c 
        call fun1(xstest(i),ipar1,par1,par2,vals2(i))
  
        errs(i)=vals(i)-vals2(i)
 2400 continue
  
        call prin2('and xstest=*',xstest,ntest)
        call prin2('and vals=*',vals,ntest)
        call prin2('and errs=*',errs,ntest)
c 
c       test the evaluation of the derivative
c 
        h=0.00001
  
        call rneste_cheb(ier,w2,xtest,val,der)
        call rneste_cheb(ier,w2,xtest+h,valp,derp)
        call rneste_cheb(ier,w2,xtest-h,valm,derm)
c 
        der2=(valp-valm)/h/2
  
        call prin2('testing derivative, der=*',der,1)
        call prin2('testing derivative, der2=*',der2,1)
        call prin2('testing derivative, der-der2=*',der-der2,1)
  
        call prinf('and nn=*',nn,1)
  
  
        return
        end
c 
c 
c 
c 
c 
c 
        subroutine arrcpy(x,y,n)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1)
c 
        do 1200 i=1,n
        y(i)=x(i)
 1200 continue
c 
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
        dimension par1(1),par2(1)
c 
        real *8 f,ima
c 
        data ima/(0.0d0,1.0d0)/
c 
        f=log(x)
cccc        f=sin(10*x)
c 
        return
        end
c 
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c       This is the end of the debugging code, and the beginning
c       of the nested Chebychev interpolation code proper
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c       This file contains numerical machinery for the construction and
c       evaluation of nested Chebychev expansions of real functions
c       of one real variable. There are 3 user-callable subroutines in
c       this file: rnestex_cheb, rneste_cheb, rneste2_cheb (actually,
c       rneste2_cheb is an entry in the subroutine rneste_cheb. Following
c       is a brief description of these subroutines.
c 
c   rnestex_cheb - finds a subdivision of the interval [a,b], such that on
c       every subinterval, the user-provided function given by the
c       subroutine fun (see description below) is interpolated to
c       precision eps by a k-point Chebychev interpolation formula; then
c       it evaluates on each subinterval the coefficients of the chebychev
c       series approximation of f. The resulting is gauaranteed to accuracy
c       eps, and the output of this subroutine can be used by the
c       subroutines rneste_cheb, rneste2_cheb (see) to evaluate the
c       interpolated function at arbitrary points on the interval [a,b].
c 
c   rneste_cheb - evaluates at the point x \in [a,b] the real
c        function (together with its derivative) specified by its nested
c        Chebychev expansion, whose coefficients are stored in the input
c        array coefs. Note that this subroutine uses output of the
c        subroutine rnestex_cheb (see) as its input; it has no known uses
c        as a stand-alone device
c 
c   rneste2_cheb - evaluates at the point x \in [a,b] the real function
c        specified by its nested Chebychev expansion, whose coefficients
c        are provideds in the input array coefs. Please note that this
c        subroutine uses output of the subroutine rnestex_cheb (see) as its
c        input; it has no known uses as a stand-alone device
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        subroutine rnestex_cheb(ier,a,b,k,fun,par1,par2,par3,
     1      eps,nn,w,ifinit,w2,lenw2,keepw2)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),par1(1),par2(1),par3(1),w2(1)
c 
        external fun
c 
c       This subroutine finds a subdivision of the interval [a,b],
c       such that on every subinterval, the user-provided function
c       given by the subroutine fun (see description below) is
c       interpolated to precision eps by a k-point Chebychev interpolation
c       formula; then it evaluates on each subinterval the coefficients
c       of the Chebychev series approximation of f. The resulting
c       interpolation scheme is gauaranteed to accuracy eps, and the
c       output of this subroutine can be used by the subroutines rneste,
c       rneste2_cheb (see) to evaluate the function fun at arbitrary points
c       on the interval [a,b].
c 
c                          input parameters:
c 
c  a,b - the ends of the interval to be subdivided
c  k - the number of Chebychev nodes on each subinterval
c  fun - the user-supplied subroutine defining the function to be
c      interpolated. The calling sequence of fun must be:
c 
c      fun(x,par1,par2,par3,f),
c 
c      where x is the input argument where the function is to be
c      evaluated, apar1, par2, par3 are whatever parameters the
c      subroutine fun needs, and f (output) is the value of the
c      function to be interpolated at the point x.
c  par1,par2,par3 - whatever the user-supplied subroutine fun needs
c  eps - the accuracy to which the Chebychev expansion must represent
c     the user-supplied function for the subroutine to decide that
c     the representation is accurate
c  w - the array containing the data to be used by the subroutine,
c     created by a preceding call to it. It must be at
c     least 10*k+8*k**2+30 + 20 000 real *8 locations long.
c  ifinit - the index telling the subroutine whether it should construct
c     the parameter w (see above), or view it as an input parameter.
c     In the latter case, this  must have been created by a prior call
c     to this subroutine with the same parameter k.
c     ifinit=1 means that the parameter w will be created.
c     ifinit=0 means that the parameter w will be viewed
c     as an input parameter.
c  lenw2 - the length of the user-provided array w2 (see below)
c 
c                     output parameters:
c 
c  ier - error return code.
c         ier=0 means successful execution
c         ier=4 means that the space allocated in the array ab is
c         insufficient
c  nn - the number of subintervals created by the subroutine on
c     the interval [a,b]
c  w2 - coefficients of all nn nested Chebychev expansions (also,
c     certain other types of information); can be used by the
c     subroutines rneste, rneste2_cheb (see) for the evaluation
c     of the function fun at arbitrary points on the interval [a,b]
c  keep - the length of the part of the array w2 (in real *8 words)
c     that should not be changed between the call to this subroutine
c     and the subsequent calls to the subroutines rneste, rneste2_cheb
c     (see)
c 
c        . . . construct a subdivision of the interval [a,b], such that
c              on every subinterval, the user-provided function fun is
c              interpolated to to precision eps by a k-point Chebychev
c              expansion
c 
         ier=0
c 
         iw=1
         lw=16*k+8*k**2+30
  
         call prinf('lw=*',lw,1)
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
c 
         call prinf('ltot=*',ltot,1)
  
         ltotana=10*k**2+20*k+100
  
         call prinf('ltotana=*',ltotana,1)
c 
         iab=21
         lenab=lenw2-iab
c 
         call rallinco_cheb(ier,a,b,k,fun,par1,par2,par3,eps,
     1       w2(iab),lenab,nn,labused,w(iw),ifinit)
c
         lab=nn*2+2
         n=k*nn
c 
         ifs=iab+lab
cccc         lfs=n*2+10
         lfs=n+10
c 
         ix=ifs+lfs
         lx=n+2
c 
         ltotw2=ix+lx+100
c 
         if(ltotw2 .lt. lenw2) goto 1300
c 
         ier=16
         return
 1300  continue
c 
c       construct all nodes and weights on all subintervals
c 
        call rnodewht_cheb(w2(iab),nn,k,w2(ix))
c 
c        construct the Chebychev expansions of the function on all
c        subintervals
c 
        itype=2
        if(ifinit .eq. 1)
     1      call chebexps(itype,k,w(it),w(iu),w(iv),w(iww))
c 
        do 1400 i=1,n
c 
        call fun(w2(ix+i-1),par1,par2,par3,w2(ifs+i-1) )
 1400 continue
c 
        do 2200 i=1,nn
c 
        ii=(i-1)*k+1
        call ralqrma2_cheb(w(iu),w2(ifs+ii-1),w(iww),k,1)
c 
        do 1800 j=1,k
c 
        w2(ifs+ii+j-2)=w(iww+j-1)
 1800 continue
c 
 2200 continue
c 
c       store in the beginning of the array w2 the
c       map of the said array and other integer information
c 
        w2(1)=iab+0.1
        w2(2)=nn+0.1
        w2(3)=k+0.1
        w2(4)=ifs+0.1
c 
        keepw2=ltotw2
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine rneste_cheb(ier,w2,x,val,der)
        implicit real *8 (a-h,o-z)
        save
        real *8 val,der
        dimension w2(1)
  
c 
c        This subroutine evaluates at the point x \in [a,b] the
c        real function (together with its derivative) specified
c        by its nested Chebychev expansion, whose coefficients are
c        supplied in the input array w2 (together with certain other
c        types of data). Please note that this subroutine uses output
c        of the subroutine rnestex_cheb (see) as its input; it has
c        no known uses as a stand-alone device.
c 
c                   Input parameters:
c 
c  w2 - array produced via a preceding call to the subroutine
c        rnestex_cheb (see)
c      array fs.
c  x - the point where the function is to be evaluated; must be on
c      the interval [a,b]
c 
c                    Output parameters:
c 
c  ier - error return code; ier=0 means successful conclusion
c                           ier=4 means that the user-supplied foint x
c                           is outside for which the array w2 has been
c                           constructed (see subroutine rnestex_cheb)
c  val - the value of the function fun at the point x
c  der - the derivative of the function fun at the point x
c 
c 
c        . . . decode the beginning of the array w2
c 
        iab=w2(1)
        nn=w2(2)
        k=w2(3)
        icoefs=w2(4)
c 
c       evaluate the function and its derivative
c 
        call rneste_cheb0(ier,w2(icoefs),w2(iab),nn,k,x,val,der)
c 
        return
c 
c 
c 
c 
        entry rneste2_cheb(ier,w2,x,val)
c 
c        This entry evaluates at the point x \in [a,b] the real
c        function specified by its nested Chebychev expansion, whose
c        coefficients are supplied in the input array w2 (together
c        with certain other types of data). Please note that this
c        subroutine uses output of the subroutine rnestex_cheb (see)
c        as its input; it has no known uses as a stand-alone device.
c 
c                   Input parameters:
c 
c  w2 - array produced via a preceding call to the subroutine
c        rnestex_cheb (see)
c      array fs.
c  x - the point where the function is to be evaluated; must be on
c      the interval [a,b]
c 
c                    Output parameters:
c 
c  ier - error return code; ier=0 means successful conclusion
c                           ier=4 means that the user-supplied foint x
c                           is outside for which the array w2 has been
c                           constructed (see subroutine rnestex_cheb)
c  val - the value of the function fun at the point x
c 
c 
c        . . . decode the beginning of the array w2
c 
        iab=w2(1)
        nn=w2(2)
        k=w2(3)
        icoefs=w2(4)
c 
c       evaluate the function and its derivative
c 
        call rneste2_cheb0(ier,w2(icoefs),w2(iab),nn,k,x,val)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine rneste_cheb0(ier,coefs,ab,nn,k,x,val,der)
        implicit real *8 (a-h,o-z)
        save
        dimension ab(2,1)
c 
        real *8 coefs(1),val,der
c 
c       . . . find the subinterval in which the point x lives
c 
        ier=0
        call rfindinte_cheb(ier,x,ab,nn,intnum)
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
cccc        call chfuncder(t,VAL,der,coefs(jj),k-1)
        call chfunder(t,VAL,der,coefs(jj),k-1)
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
        subroutine rneste2_cheb0(ier,coefs,ab,nn,k,x,val)
        implicit real *8 (a-h,o-z)
        save
        dimension ab(2,1)
c 
        real *8 coefs(1),val
c 
c        This subroutine evaluates at the point x \in [a,b] the real
c        function specified by its nested Chebychev expansion, whose
c        coefficients are supplied in the input array coefs. Please
c        note that this subroutine uses output of the subroutine rnestex_cheb
c        (see) as its input; it has no known uses as a stand-alone device.
c 
c                   Input parameters:
c 
c  coefs - the nested Chebychev expansions of the functions stored in the
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
        call rfindinte_cheb(ier,x,ab,nn,intnum)
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
        call chebval(t,VAL,coefs(jj),k-1)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine rnodewht_cheb(ab,nn,k,x)
        implicit real *8 (a-h,o-z)
        save
        dimension ab(2,1),x(1),t(2000)
c 
c 
        data kold/-10/
c 
c       construct the nodes and weights of the k-point Gaussian
c       quadrature on the interval [-1,1]
c 
        itype=0
c 
        if(k .ne. kold) call chebexps(itype,k,t,u,v,w)
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
         subroutine rallinco_cheb(ier,a,b,k,fun,ipar1,par1,
     1    par2,eps,ab,lab,nn,labused,w,ifinit)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),par1(1),par2(1),ab(2,1)
c 
c       this subroutine finds a subdivision of the interval [a,b],
c       such that on every subinterval, the user-provided function
c       fun is integrated to precision eps by a k-point Gaussian
c       quadrature. Actually, this is a memory management routine;
c       all the work is done by the routine rallinte_cheb (see).
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
c  eps - the accuracy to which the Chebychev expansion must represent
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
c       . . . allocate memory for the subroutine rallinte_cheb
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
         call rallinte_cheb(ier,a,b,fun,ipar1,par1,par2,k,ab,lab,
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
        subroutine rbubble_cheb(ab,nn)
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
         subroutine rallinte_cheb(ier,a,b,fun,ipar1,par1,par2,
     1    k,ab,lab,eps,t,u,v,whts,x,f,coefs,nn,lused,ifinit7)
        implicit real *8 (a-h,o-z)
        save
        dimension t(1),u(2*k,2*k),v(2*k,2*k),whts(1),x(1),
     1      par1(1),par2(1),ab(2,1)
        equivalence (rea,intint)
c 
        real *8 f(1),coefs(1)
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
c  eps - the accuracy to which the Chebychev expansion must represent
c     the user-supplied function for the subroutine to decide that
c     the representation is accurate
c  t - the Gaussian nodes on the interval [-1,1]. This is only an
c     input parameter is ifinit7 (see below) has been set to 0;
c     otherwise, it is an output parameter.
c  u - the 2*k * 2*k  matrix converting the  values at of a polynomial
c     of order 2*k-1 at 2*k Chebychev nodes into the coefficients of its
c     Chebychev expansion. This is only an input parameter is ifinit7
c     (see below) has been set to 0; otherwise, it is an output parameter.
c  v - the 2*k * 2*k matrix converting the coefficients
c     of an 2*k-term Chebychev expansion into its values at 2*k Chebychev
c     nodes (note that v is the inverse of u). This is only an input
c     parameter is ifinit7 (see below) has been set to 0; otherwise, it
c     is an output parameter.
c  whts - the Gaussian quadrature weights (2*k of them).  This is only
c     an input parameter is ifinit7 (see below) has been set to 0;
c     otherwise, its an output parameter.
c  ifinit7 - the index telling the subroutine whether it should construct
c     the parameters t,u,v,whts (by calling the subroutine chebexps),
c     or view them as input parameters. In the latter case, these
c     must have been created by a prior call to this subroutine (or to
c     subroutine chebexps, for an adventurous and creative user).
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
c  coefs - the coefficients of the Chebychev expansion of the
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
        call rifsplit_cheb(ier,ab(1,nn),ab(2,nn),k,fun,ipar1,par1,
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
        call rbubble_cheb(ab,nn)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine rifsplit_cheb(ier,a,b,k,fun,ipar1,par1,
     1    par2,eps,t,u,v,whts,x,ifinit,f,coefs)
        implicit real *8 (a-h,o-z)
        save
        dimension t(1),u(2*k,2*k),v(2*k,2*k),whts(1),x(1),
     1      par1(1),par2(1)
c 
        real *8 f(1),coefs(1)
c 
c       this subroutine determines whether a k-point Chebychev
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
c  eps - the accuracy to which the Chebychev expansion must represent
c     the user-supplied function for the subroutine to decide that
c     the representation is accurate
c  t - the Gaussian nodes on the interval [-1,1]. This is only an
c     input parameter is ifinit (see below) has been set to 0;
c     otherwise, it is an output parameter.
c  u - the 2*k * 2*k  matrix converting the  values at of a polynomial
c     of order 2*k-1 at 2*k Chebychev nodes into the coefficients of its
c     Chebychev expansion. This is only an input parameter is ifinit
c     (see below) has been set to 0; otherwise, it is an output parameter.
c  v - the 2*k * 2*k matrix converting the coefficients
c     of an 2*k-term Chebychev expansion into its values at 2*k Chebychev
c     nodes (note that v is the inverse of u). This is only an input
c     parameter is ifinit (see below) has been set to 0; otherwise, it
c     is an output parameter.
c  whts - the Gaussian quadrature weights (2*k of them).  This is only
c     an input parameter is ifinit (see below) has been set to 0;
c     otherwise, its an output parameter.
c  ifinit - the index telling the subroutine whether it sholud construct
c     the parameters t,u,v,whts (by calling the subroutine chebexps),
c     or view them as input parameters. In the latter case, these
c     must have been created by a prior call to this subroutine (or to
c     subroutine chebexps, for an adventurous and creative user).
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
c  coefs - the coefficients of the Chebychev expansion of the
c     function f on the interval [a,b]
c 
c       . . . if the need be, construct the gaussian nodes on the
c             interval [-1,1], and the matrix conevrting the values
c             of the function at these nodes into coefficients of
c             its Chebychev expansion
c 
        if(ifinit .eq. 0) goto 1200
c 
        itype=2
        call chebexps(itype,k*2,t,u,v,whts)
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
c       construct the Chebychev expansion of the function on
c       the interval
c 
        call ralqrma2_cheb(u,f,coefs,k*2,k+1)
c 
c       scan the coefficients and see if the last k
c       of them are small
c 
        eps22=eps**2
c 
        ier=0
        do 1800 i=1,k
cccc        d=coefs(k+i)*conjg(coefs(k+i))
        d=coefs(k+i)**2
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
        subroutine ralqrma2_cheb(a,x,y,n,n0)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n)
        real *8 x(1),y(1),d
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
        subroutine rfindinte_cheb(ier,x,ab,nn,intnum)
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
