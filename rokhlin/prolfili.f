        implicit real *8 (a-h,o-z)
        real *8 proj2(300 000),rinterp2(300 000),
     1      xs(10 000),ws(10 000),xs2(10 000),ws2(10 000)
c 
        real *8 w(10 000 000),u(500 000),v(500 000),
     1      riprfast(300 000),errs1(10 000),errs2(10 000),
     2      coefstest(10 000),errs3(10 000),errs4(10 000)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER wavesnum'
        READ *,wavesnum
        CALL PRIN2('wavesnum=*',wavesnum,1 )
  
c 
c       initialize the filtering/interpolation subroutine
c 
        eps=1.0d-12
        eps=1.0d-13
        eps=1.0d-7
  
  
        ifrefine=1
        lenw=10 000 000
c 
  
        call fileflush(13)
  
  
        call prolfili(ier,wavesnum,eps,ifrefine,
     1      xs,ws,xs2,ws2,riprfast,n,rinterp2,proj2,
     2      u,v,w,lenw,keep,lused)
  
        call fileflush(13)
  
  
  
        call prinf('and after prolfili, n=*',n,1)
        call prin2('and after prolfili, xs=*',xs,n)
        call prin2('and after prolfili, xs2=*',xs2,n*2)
        call prinf('and after prolfili, keep=*',keep,1)
        call prinf('and after prolfili, lused=*',lused,1)
  
        call fileflush(13)
  
  
c 
c        conduct mass testing of the interpolation subroutine
c 
        done=1
        pi=atan(done)*4
  
        coefmax=wavesnum*pi
c 
         ntest=32
         hcoef=coefmax/ntest
c 
         do 2400 i=1,ntest
c 
         coefstest(i)=i*hcoef
         call funtesin(coefstest(i))
  
        call cintetest(xs,ws,xs2,ws2,n,rinterp2,riprfast,
     1      err1max,err2max)
c 
        errs1(i)=err1max
        errs2(i)=err2max
  
 2400 continue
c 
        call prin2(
     1      'after mass testing of interpolation, errs1=*',
     2      errs1,ntest)
c 
        call prin2(
     1      'after mass testing of interpolation, errs2=*',
     2      errs2,ntest)
c 
        call prin2('while coefstest=*',coefstest,ntest)
  
        call fileflush(13)
  
  
  
c 
c       test the obtained filtering matrix
c 
         coef=wavesnum*pi *3
c 
         coef2=wavesnum*pi
cccc         coef2=wavesnum*pi /4
c 
         ntest=32
         hcoef=coef/ntest
         ifinit=1
  
         do 2800 i=1,ntest
c 
  
cccc        call prinf('i=*',i,1)
  
         coefstest(i)=i*hcoef
         call funtesin(coefstest(i))
c 
        call cproj2te(xs,ws,xsr,wsr,n,
     1    v,w,xs2,ws2,proj2,riprfast,coef2,errmax,err0max,
     2      ifinit)
  
cccc  call prin2('after cproj2te, errmax=*',errmax,1)
cccc        call prin2('after cproj2te, err0max=*',err0max,1)
  
c 
         ifinit=0
c 
        errs3(i)=err0max
        errs4(i)=errmax
 2800 continue
  
  
        call prin2('after mass test of filtering, errs3=*',
     1      errs3,ntest)
  
        call prin2('after mass test of filtering, errs4=*',
     1      errs4,ntest)
  
  
        call fileflush(13)
  
  
  
        stop
        end
  
c 
c 
c 
c 
c 
        subroutine cproj2te(xs,ws,xsr,wsr,n,
     1    v,w,xs2,ws2,proj2,riprfast,coef,errmax,err0max,
     2      ifinit)
        implicit real *8 (a-h,o-z)
        save
        real *8 xs(1),ws(1),
     3      xsr(1),wsr(1),ts(10 000),whts(10 000),
     4      tsr(10 000),whtsr(10 000),www(10 000)
c 
        dimension xs2(1),whts2(10 000),proj2(n*2,n*2),ts2(10 000),
     1      riprfast(1),freqs(10 000),errs(10 000)
c 
        complex *16 fproj(10 000),ftest2(10 000),
     2      expprol(10 000),explege(10 000),vals(10 000),
     5      vals7(10 000),fproj2(10 000),vals8(10 000),
     3      rint1,rint2,der
c 
c 
c       construct Gaussian nodes on the intervals [-1,1], [0,1]
c 
        if(ifinit .eq. 0) goto 1130
  
        ngauss=1000
cccc        ngauss=500
c 
        itype=1
        call legeexps(itype,ngauss,ts,u,v,whts)
c 
        do 1100 i=1,ngauss
c 
        tsr(i)=(ts(i)+1)/2
        whtsr(i)=whts(i)/2
c 
 1100 continue
 1130 continue
  
cccc        call prin2('in cproj2te, tsr=*',tsr,ngauss)
  
c 
        do 1150 i=1,ngauss
c 
        ts2(i)=-tsr(ngauss-i+1)
        ts2(ngauss+i)=tsr(i)
c 
        whts2(i)=whtsr(i)
        whts2(ngauss+i)=whtsr(i)
 1150 continue
  
cccc        call prin2('in cproj2te, whts2=*',whts2,n*2)
  
c 
c       construct the test function
c 
        do 1200 i=1,n*2
c 
        call funtest(xs2(i),ftest2(i))
 1200 continue
  
cccc        call prin2('in cproj2te, ftest2=*',funtest2,n*4)
  
  
c 
c       calculate the projection of the test function on the
c       subintervals on band-limited function on the whole interval
c 
        call cprolma2(proj2,ftest2,n,n*2,fproj)
c 
        call cprolfia(riprfast,n,ftest2,fproj2,www)
        err0max=0
        do 1300 i=1,n
c 
        d=abs(fproj(i)-fproj2(i))
        if(d .gt. err0max) err0max=d
 1300 continue
c 
c       expand the projected function into prolate series
c 
        call cprolmat(v,fproj,n,expprol)
c 
c        convert the prolate series of the interpolated function
c        into a legendre series
c 
        call cprotole(expprol,n,w,explege,nlege)
c 
c       evaluate and plot the projected function at the Gaussian nodes
c       on the interval [-1,1]
c 
        do 1400 i=1,ngauss
c 
        call legecFDE(ts(i),VALs(i),der,explege,Nlege-1)
c 
        call funtest(ts(i),vals7(i))
 1400 continue
c 
cccc        iw=21
cccc        call lotagraph2(iw,ts,vals,ngauss,ts,vals7,ngauss,
cccc     1      'filtered against unfiltered function*')
c 
c       calculate the inner products between the projected function
c       and various complex exponentials, and compare the obtained
c       products with those for the original function
c 
        freqmax=coef
c 
        ntest=32
        hfreq=freqmax/ntest
        errmax=0
c 
        do 2600 i=1,ntest
c 
        freqs(i)=i*hfreq
        call funtes2in(freqs(i))
c 
        rint1=0
        rint2=0
        do 2200 j=1,ngauss
c 
        call funtest2(ts(j),vals8(j))
c 
        rint1=rint1+vals7(j)*vals8(j)*whts(j)
        rint2=rint2+vals(j)*vals8(j)*whts(j)
 2200 continue
c 
        errs(i)=abs(rint2-rint1)
  
        if(abs(errs(i)) .gt. errmax) errmax=abs(errs(i))
c 
 2600 continue
c 
cccc        call prin2('and errs in cproj2te are*',errs,ntest)
cccc        call prin2('and freqs are*',freqs,ntest)
cccc        call prin2('and errmax=*',errmax,1)
        return
        end
c 
c 
c 
c 
c 
        subroutine cprolmat(a,x,n,y)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n)
        complex *16 x(1),y(1),d
c 
        do 1400 i=1,n
        d=0
        do 1200 j=1,n
        d=d+a(i,j)*x(j)
 1200 continue
        y(i)=d
 1400 continue
  
        return
        end
c 
c 
c 
c 
c 
        subroutine prolmatv(a,x,n,y)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),x(1),y(1)
c 
        do 1400 i=1,n
        d=0
        do 1200 j=1,n
        d=d+a(i,j)*x(j)
 1200 continue
        y(i)=d
 1400 continue
  
        return
        end
c 
c 
c 
c 
c 
        subroutine cintetest(xs,ws,xs2,ws2,n,rinterp2,riprfast,
     1      err1max,err2max)
        implicit real *8 (a-h,o-z)
        save
        real *8 xs(1),ws(1),errs(10 000),
     1      rinterp2(2*n,n),xs2(1),ws2(1),riprfast(1),www(10 000)
c 
        complex *16 ftest(10 000),ftest2(10 000),fproj(10 000),
     1      fproj2(10 000)
c 
c       construct the test function
c 
        do 1200 i=1,n
c 
        call funtest(xs(i),ftest(i))
  
        call funtest(xs2(i),ftest2(i))
        call funtest(xs2(n+i),ftest2(n+i))
c 
 1200 continue
c 
c       calculate the interpolation of the test function on the
c       interval on band-limited functions on the subinterval
c 
        call cprolma2(rinterp2,ftest,n*2,n,fproj)
c 
cccc        call prin2('in cintetest, fproj=*',fproj,n*2)
cccc        call prin2('in cintetest, ftest2=*',ftest2,n*2)
  
  
        call cprolina(riprfast,n,ftest,fproj2,www)
  
cccc        call prin2('in cintetest, fproj2 from cprolina=*',fproj2,n*4)
cccc        call prin2('in cintetest, ftest2 from cprolina=*',ftest2,n*4)
c 
        err1max=0
        err2max=0
c 
        do 1400 i=1,n*2
        errs(i)=fproj(i)-ftest2(i)
  
c 
        d1=abs(fproj(i)-ftest2(i))
        d2=abs(fproj2(i)-ftest2(i))
c 
        if(err1max .lt. d1) err1max=d1
        if(err2max .lt. d2) err2max=d2
 1400 continue
  
cccc         call prin2('in cintetest, errs=*',errs,n*2)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine funtest2(x,f)
        implicit real *8 (a-h,o-z)
c 
        save
        complex *16 f,ima
c 
        data ima/(0.0d0,1.0d0)/
c 
        f=sin(x*coef)+ima*cos(x*coef)
c 
        return
c 
c 
c 
c 
        entry funtes2in(coef7)
        coef=coef7
        return
        end
c 
c 
c 
c 
c 
        subroutine funtest(x,f)
        implicit real *8 (a-h,o-z)
c 
        save
        complex *16 f,ima
c 
        data ima/(0.0d0,1.0d0)/
  
c 
c 
        f=sin(x*coef)+ima*cos(x*coef)
c 
        return
c 
c 
c 
c 
        entry funtesin(coef7)
        coef=coef7
        return
        end
  
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c       This is the end of the testing code and the beginning of the
c       prolate expansion code proper.
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c       This file contains a set of subroutines for the interpolation
c       and filtering of band-limited functions on the interval (for
c       convenience, we will be assuming that the interval is [-1,1]),
c       and on the square (again, we will be assuming that the square
c       is [-1,1] \times [-1,1]).
c 
c       The real functions on the interval are handled by the subroutine
c       prolfiap (it has two entries, prolfiap for the filtering, and
c       prolinap for the interpolation). The complex functions on the
c       interval are handled by the subroutine cprolfia (it has two
c       entries, cprolfia for the filtering and cprolina for the
c       interpolation).
c 
c       The complex functions on the square are  handled by the subroutine
c       prosqint (it has two entries, prolsqint for the interpolation,
c       and prosqfil for the filtering). At this time (4.20.00) there is
c       no code for the handling of real functions on the square.
c 
c       Before any interpolation and filtering can be performed, an
c       initialization subroutine has to be called. There are two such
c       subroutines. The initialization subroutine for the treatment of
c       functions on the square is prolfili (see); the initialization
c       subroutine for dealing with functions on the square is prosqbld
c       (see). While the subroutine prosqbld prepares all the data that
c       are needed for dealing with both squares and intervals, the
c       subroutine prolfili only prepares data for dealing with intervals
c       (which is only fair, since the subroutine prosqbld actually calls
c       the subroutine prolfili).
c 
c   PLEASE NOTE THAT THE SUBROUTINE PROLFILI BELOW IS THE MOST EXTENSIVELY
c   DOCUMENTED SUBROUTINE IN THIS FILE. IF YOU ARE LOOKING FOR INFORMATION,
c   PROLFILI IS YOUR BEST BET (not to say that it is a very good bet!)
c 
c    There are 8 subroutines and entries in this file that are supposed
c       to be user-callable, as follows.
c 
c   prosqbld - constructs matrices of filtering and
c         interpolation operators for band-limited functions,
c         both for functions defined on the interval and for
c         functions defined on the square
c 
c   prolfili - constructs matrices of filtering and
c         interpolation operators for band-limited functions
c         defined on the interval
c   prosqint - interpolates a complex-valued function on the
c       square, from an  n \times n prolate tensor product grid
c       onto an 2*n \times 2*n grid.
c 
c   prosqfil - filters a complex-valued function on the square,
c       and evaluates the filtered function at the nodes of a sparser
c       grid on the said square.
c 
c   cprolfia - applies the filtering operator to the (complex) input
c       vector of dimensionality 2*n, obtaining the (complex) output
c       vector of dimensionality n.
c 
c   cprolina - applies the interpolating operator to the (complex) input
c       vector of dimensionality n, obtaining the (complex) output
c       vector of dimensionality 2*n.
c 
c   prolfiap - applies the filtering operator to the (real input
c       vector of dimensionality 2*n, obtaining the (real) output
c       vector of dimensionality n.
c 
c   prolinap - applies the interpolating operator to the (real) input
c       vector of dimensionality n, obtaining the (real) output
c       vector of dimensionality 2*n.
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine prosqint(riprfast,n,fssqin,fssqout,w)
        implicit real *8 (a-h,o-z)
        save
        real *8 riprfast(1),fssqin(2*n,2*n),fssqout(n,n),w(1)
c 
c       This entry interpolates a complex-valued function on the
c       square, from an  n \times n prolate tensor product grid
c       (which has been (hopefully) obtained via a preceding call to
c       the subroutine prosqbld (see) onto an 2*n \times 2*n grid
c       xs2sq, which has been obtained via the same call to the
c       subroutine proqsbld (see). The interpolation operator is
c       supplied by the user in the array riprfast, which (hopefully)
c       has been generated by the same prior call to the subroutine
c       prosqbld (see).
c 
c                     Input parameters:
c 
c  riprfast - the array containing symmetrized versions of the operators
c         rinterp2, proj2, to be used by this entry for the application
c         the operator proj2 to arbitrary vectors
c  n - the number of nodes discretizing the interval [-1,1] with precision
c         eps and band-limit c; hopefully has been generated by the same
c         prior call to the subroutine prosqbld (see)
c  fssqin - the input function, tabulated on the n \times n grid xssq
c 
c                     Output parameters:
c 
c  fssqout - the output function, tabulated on the 2*n \times 2*n grid xs2sq
c 
c                     Work arrays:
c 
c  w - must be at least 20*n+100  real *8 elements long
c 
c 
c       . . . allocate memory for the interpolation
c             of the function on the square
c 
        iw=1
        lw=n*8+8
c 
        iww=iw+lw
        lww=n*4+4
c 
        iwww=iww+lww
        lwww=n*4+4
c 
c       . . . interpolate
c 
        call prosqin0(riprfast,n,fssqin,fssqout,w(iw),w(iww),w(iwww) )
c 
        return
c 
c 
c 
c 
        entry prosqfil(riprfast,n,fssqin,fssqout,w)
c 
c       This entry filters a complex-valued function on the square,
c       and evaluates the filtered function at the nodes of a sparser
c       grid on the said square. The input function fssqin is tabulated
c       on the 2*n \times 2*n grid xs2sq (which has been (hopefully)
c       obtained via a preceding call to the subroutine prosqbld (see);
c       the output function fssquot is tabulated at the nodes of the
c       grid xssq, also obtained via a preceding call to prosqbld.
c       The filtering operator is supplied by the user in the array
c       riprfast, which (hopefully) has been generated by the same
c       prior call to the subroutine prosqbld (see)
c 
c                     Input parameters:
c 
c  riprfast - the array containing symmetrized versions of the operators
c         rinterp2, proj2, to be used by this entry for the application
c         the operator proj2 to arbitrary vectors
c  n - the number of nodes discretizing the interval [-1,1] with precision
c         eps and band-limit c; hopefully has been generated by the same
c         prior call to the subroutine prosqbld (see)
c  fssqin - the input function, tabulated on the 2*n \times 2*n grid xs2sq
c 
c                     Output parameters:
c 
c  fssqout - the output function, tabulated on the n \times n grid xssq
c 
c                     Work arrays:
c 
c  w - must be at least n**2*4+20*n+100  real *8 elements long
c 
c 
c       allocate memory for the filtering
c       of the function on the square
c 
        iw=1
        lw=n*8+8
c 
        iww=iw+lw
        lww=n*4+4
c 
        iwww=iww+lww
        lwww=n*4+4
c 
        iwork=iwww+lwww
        lwork=n**2*4 +10
c 
c       . . . filter
c 
        call prosqfi0(riprfast,n,fssqin,fssqout,
     1      w(iww),w(iwww),w(iw),w(iwork) )
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine prosqin0(riprfast,n,fssqin,fssqout,w,ww,www)
        implicit real *8 (a-h,o-z)
        save
        real *8 riprfast(1)
        complex *16 fssqin(n,n),fssqout(2*n,2*n),
     1          ww(1),www(1),w(1)
c 
c       interpolate the function on the square in the horizontal
c       direction
c 
        do 1600 i=1,n
c 
        do 1200 j=1,n
c 
        ww(j)=fssqin(i,j)
 1200 continue
c 
        call cprolina(riprfast,n,ww,www,w)
c 
        do 1400 j=1,n*2
c 
        fssqout(i,j)=www(j)
 1400 continue
c 
 1600 continue
c 
c       interpolate the function on the square in the vertical
c       direction
c 
        do 2600 i=1,n*2
c 
        do 2200 j=1,n
c 
        ww(j)=fssqout(j,i)
 2200 continue
c 
        call cprolina(riprfast,n,ww,www,w)
c 
        do 2400 j=1,n*2
c 
        fssqout(j,i)=www(j)
 2400 continue
c 
 2600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine prosqfi0(riprfast,n,fssqin,fssqout,
     1      ww,www,w,work)
        implicit real *8 (a-h,o-z)
        save
        real *8 riprfast(1)
c 
        complex *16 fssqin(2*n,2*n),fssqout(n,n),work(n*2,n),
     1          ww(1),www(1),w(1)
c 
c       filter the function on the square in the horizontal
c       direction
c 
        do 1600 i=1,n*2
c 
        do 1200 j=1,n*2
c 
        ww(j)=fssqin(i,j)
 1200 continue
c 
        call cprolfia(riprfast,n,ww,www,w)
c 
        do 1400 j=1,n
c 
        work(i,j)=www(j)
 1400 continue
c 
 1600 continue
c 
c       filter the function on the square in the horizontal
c       direction
c 
        do 2600 i=1,n
c 
        do 2200 j=1,n*2
c 
        ww(j)=work(j,i)
 2200 continue
c 
        call cprolfia(riprfast,n,ww,www,w)
c 
        do 2400 j=1,n
c 
        fssqout(j,i)=www(j)
 2400 continue
c 
 2600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine prosqbld(ier,wavesnum,eps,ifrefine,
     1      xs,ws,xs2,ws2,riprfast,n,rinterp2,proj2,
     2      xssq,wssq,xs2sq,ws2sq,u,v,w,lenw,keep,lused)
        implicit real *8 (a-h,o-z)
c 
        save
        real *8 w(1),u(1),v(1),xs(1),ws(1),xs2(1),ws2(1),
     1          proj2(1),rinterp2(1),riprfast(1),
     2          xssq(1),wssq(1),xs2sq(1),ws2sq(1)
c 
c         This subroutine constructs matrices of filtering and
c         interpolation operators for band-limited functions,
c         both for functions defined on the interval (for
c         simplicity, we will be assuming that the interval is
c         [-1,1]), and on the square (again, we assume that the
c         square is [-1,1] \times [-1,1]).
c 
c         The environment is not assumed to be periodic,
c         and the operators constructed by this subroutine are
c         filtering and interpolation operators in the weak sense
c         only (actually, the interpolation operator is a true
c         interpolation operator, in a strong sense; the filtering
c         operator filters things in a weak sense only). This
c         subroutine produces two operators: interpolation and
c         filtering.
c 
c         The filtering and interpolation operators are returned
c         by the subroutine in two forms: matrices rinterp2, proj2,
c         and the array riprfast. The matrices rinterp2, proj2 are
c         truly matrices representing the operators of interpolation
c         and filtering, respectively (see detailed descriptions
c         below). The user can be apply them to vectors of function
c         values (for example, via the subroutine prolmav2 present
c         in this package; on the other hand, the user might be
c         able to write his own subroutine applying a matrix to a
c         vector).
c 
c         IMPORTANT NOTE: THE SIZES OF OUTPUT ARRAYS PRODUCED BY THIS
c         SUBROUTINE DEPEND ON THE OUTPUT PARAMETER N, UNKNOWN BEFORE
c         THE SUBROUTINE HAS BEEN CALLED. THIS IS A SIGNIFICANT
c         DRAWBACK, SOMEWHAT ALLEVIATED BY THE FACT THAT IN ANTICIPATED
c         APPLICATIONS, EVERYTHING IS FAIRLY SMALL-SCALE. THE FOLLOWING
c         SIZES ARE ALWAYS SUFFICIENT:
c 
c      xs,ws :     c+200 real *8 elements long each
c      xs2,ws2 :   2*c+400 real *8 elements long each
c      riprfast:   4*(c+200)**2+100 real *8 elements long
c      rinterp2:   2*(c+200)**2 real *8 elements long
c      proj2:      2*(c+200)**2 real *8 elements long
c      u,v:        (c+200)**2 real *8 elements long each
c      w  -  make it long, and tell the subroutine honestly how long
c                  it is, via the parameter lenw.
c 
c         IN THE ABOVE, c IS GIVEN BY THE FORMULA
c 
c                    c=\pi * wavesnum
c 
c         The array riprfast contains symmetrized versions
c         of matrices rinterp2, proj2; these versions are used by
c         the entries prolfiap, prolinap of the subroutine prolfiap
c         (see) to apply the filtering and interpolation operators
c         for about 1/2 of the cost of applying the matrices rinterp2,
c         proj2 directly. The reason for the existence of the array
c         riprfast and the subroutine prolfiap is the fact that the
c         matrices rinterp2, proj2 have a symmetry in them; thus,
c         a factor of two can be saved compared to their direct
c         application to vectors.
c 
c         To look at the matters in finer detail, the interpolation
c         operator is represented by the matrix rinterp2, converting
c         a band-limited function on the n-point prolate grid on
c         the interval [-1,1] (nodes of this grid are returned to the
c         user in the array xs) into a function on the 2*n-point grid
c         on the interval [-1,1] (returned to the user in the array
c         xs2). The latter grid is obtained by concatenating two
c         n-point prolate grids on the intervals [-1,0], [0,1].
c 
c         Similarly, the filtering operator is represented by the
c         matrix proj2, converting a band-limited function on the
c         2*n-point grid on the interval [-1,1] (returned to the
c         user in the array xs2), into a filtered function at the
c         nodes of an  n-point grid on the interval [-1,1] (returned
c         to the user in the array xs).
c 
c 
c                      Input parameters:
c 
c  wavesnum - the number of wavelengths on the interval [-1,1]
c         at the band-limit of the functions TO BE INTERPOLATED;
c         at the band-limit of the functions to be filtered, there
c         are 3*wavesnum wavelength on the interval [-1,1]
c  eps - the precision to which filtering and interpolation are to
c         be performed. Generally the accuracy of the interpolation
c         will be eps (a little better, as a matter of fact); the
c         accuracy of filtering will be considerably (1-2 digits)
c         better.
c  ifrefine - the parameter telling the subroutine which type of
c       prolate nodes and weights to generate:
c    ifrefine=0 will cause the nodes based on the division theorem
c       to be generated (i.e. the resulting nodes will be th eroots
c       of the appropriately chosen prolate function)
c    ifrefine=1 will cause the optimal prolate nodes to be generated.
c       PLEASE NOTE THAT THE TWO SETS OF NODES ARE VERY CLOSE TO EACH
c       OTHER. AS LONG AS THEY ARE USED STRICTLY FOR INTERPOLATION,
c       THERE IS NO SERIOUS ADVANTAGE FOR EITHER SET. HOWEVER, WHEN
c       USED AS QUADRATURE NODES, THE OPTIMAL NODES GAIN A COUPLE
c       OF DIGITS ON THE ONES PRODUCED WITH IFREFINE=0. On the other
c       hand, the accuracy of both quadratures is better than eps
c       (see above), and ifrefine=1 causes a MUCH greater expenditure
c       of CPU time than ifrefine=0.
c  lenw - the amount of space (n real *8 words) provided by the user
c         in the work array w
c 
c                      Output parameters:
c 
c  xs - the n prolate nodes on the interval [-1,1]; the functions will
c         be interpolated from these nodes onto nodes xs2, and filtered
c         from the nodes xs2 to the nodes xs.
c  ws - quadrature weights corresponding to the nodes {xs}
c  xs2 - the nodes on the interval [-,1], by concatenating sets on n
c         prolate nodes on the intervals [-1,0], [0,1]; the functions
c         will be filtered from the nodes {xs2} onto nodes xs, and
c         interpolated from the nodes xs to the nodes {xs2}.
c  ws2 - quadrature weights corresponding to the nodes {xs2}
c  riprfast - the array containing symmetrized versions of the operators
c         rinterp2, proj2, to be used by the subroutine prolfiap (see)
c         for the application of these operators to arbitrary vectors.
c         Please note that 2*n**2+100 elements (real *8 elements, that
c         is) of the array riprfast have to remain unchanged between
c         the call to this subroutine and subsequent calls to the
c         subroutines prosqint, prosqfil, cprolfia, cprolina, prolfiap,
c         prolinap, each of which uses this array for its own dark
c         purposes.
c  rinterp2 - the 2*n \times n matrix interpolating band-limited functions
c         from the nodes xs onto nodes xs2
c  proj2 - the n \times 2*n matrix filtering band-limited functions
c         from the nodes xs2 onto nodes xs
c  xssq - the nodes on the square [-1,1] \times [-1,1] (n**2 of them
c         things), obtained as the tensor square of the prolate nodes on
c         the interval [-1,1]
c  wssq - the weights corresponding to the nodes xssq
c 
c  xs2sq - the nodes on the square [-1,1] \times [-1,1] (n**2 * 4  of
c         them), obtained as the tensor square of the nodes xs2 on
c         the interval [-1,1]
c  ws2sq - the weights corresponding to the nodes xs2sq
c 
c  n - the number of nodes discretizing the interval [-1,1] with precision
c         eps and band-limit c; please note that there are n nodes xs,
c         and 2*n nodes xs2
c  u - the matrix converting the n coefficients of a prolate expansion
c       on the interval [-1,1] into the values of the said expansion
c       at prolate nodes; provided for user's entertainment
c  v - the matrix converting the values of a prolate expansion
c         at prolate nodes into  the n coefficients of the said prolate
c         expansion; provided for user's entertainment. Needless to say,
c         u,v are inverses of each other.
c  keep - the number of elements in the array w that have to be
c       unchanged between the call to this subroutine and the subsequent
c       calls to the subroutine proleva and its relatives (see)
c  lused - integer number telling the user how much space in the
c       array w (in real *8 words) the subroutine actually used.
c  w - contains data that can be used by the subroutine and its relatives
c         to evaluate prolate functions at arbitrary points on R^1, and
c         for whatever other needs of the subroutines in the prolcrea
c         package; otherwise, it is a work array; do not skimp on its
c         size!
c 
c                      Work arrays:
c 
c  w - must be fairly large (its size grows quadratically with wavesnum)
c 
c 
c         construct the one-dimensional nodes and transformations
c 
        call prolfili(ier,wavesnum,eps,ifrefine,
     1      xs,ws,xs2,ws2,riprfast,n,rinterp2,proj2,
     2      u,v,w,lenw,keep,lused)
c 
        if(ier .ne. 0) return
c 
c         construct the two-dimensional nodes
c 
        call prosqbl0(xs,ws,xs2,ws2,n,
     1      xssq,wssq,xs2sq,ws2sq)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine prosqbl0(xs,ws,xs2,ws2,n,
     1      xssq,wssq,xs2sq,ws2sq)
        implicit real *8 (a-h,o-z)
c 
        save
        real *8 xs(n),ws(n),xs2(n*2),ws2(n*2),
     1          xssq(2,n,n),wssq(n,n),xs2sq(2,n*2,n*2),
     2      ws2sq(2*n,2*n)
c 
c       construct the small square
c 
        do 1400 i=1,n
        do 1200 j=1,n
c 
        xssq(1,i,j)=xs(i)
        xssq(2,i,j)=xs(j)
        wssq(i,j)=ws(i)*ws(j)
 1200 continue
 1400 continue
c 
c       construct the small square
c 
        do 2400 i=1,n*2
        do 2200 j=1,n*2
c 
        xs2sq(1,i,j)=xs2(i)
        xs2sq(2,i,j)=xs2(j)
        ws2sq(i,j)=ws2(i)*ws2(j)
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
        subroutine cprolfia(riprfast,n,fin,fout,w)
        implicit real *8 (a-h,o-z)
        save
        dimension riprfast(1)
        complex *16 w(1),fin(1),fout(1)
c 
c    ATTENTION: THIS IS THE COMPLEX VERSION OF THE REAL SUBROUTINE
c       PROLFIAP, IDENTICAL TO THE LATTER IN ALL ESSENTIAL DETAILS
c 
c       This entry applies the projection operator proj2
c       to the (complex) input vector fin of dimensionality 2*n,
c       obtaining (complex) the output vector fout of dimensionality n.
c       The projection operator is supplied by the user in the
c       array riprfast, which (hopefully) has been generated by
c       a prior call to the subroutine prolfili (see)
c 
c                     Input parameters:
c 
c  riprfast - the array containing symmetrized versions of the operators
c         rinterp2, proj2, to be used by this entry for the application
c         the operator proj2 to arbitrary vectors
c  n - the number of nodes discretizing the interval [-1,1] with precision
c         eps and band-limit c; please note that there are n nodes xs,
c         and 2*n nodes xs2
c  fin - the (complex) vector of dimensionality n*2 to which the operator
c         proj2 is to be applied
c 
c                     Output parameters:
c 
c  fout - the (complex) vector of dimensionality n which is the result of
c         applying the operator proj2 to the input vector fin
c 
c                     Work arrays:
c 
c  w - must be at least n*4+20 complex *16 elements long
c 
c 
c       . . . apply the symmetrized filtering matrix to the
c             input vector fin
c 
        isproj=riprfast(1)
        iaproj=riprfast(2)
c 
        iwp=1
        lwp=n+2
c 
        iwm=iwp+lwp
        lwm=n+2
c 
        ifp=iwm+lwm
        lfp=n+2
c 
        ifm=ifp+lfp
        lfm=n+2
c 
cccc        call prolproa(riprfast(isproj),riprfast(iaproj),
        call cprolpro(riprfast(isproj),riprfast(iaproj),
     1      fin,fout,n,w(iwp),w(iwm),w(ifp),w(ifm) )
        return
c 
c 
c 
c 
        entry cprolina(riprfast,n,fin,fout,w)
c 
c 
c    ATTENTION: THIS IS THE COMPLEX VERSION OF THE REAL SUBROUTINE
c       PROLINAP, IDENTICAL TO THE LATTER IN ALL ESSENTIAL DETAILS
c 
c       This entry applies the interpolation operator rinterp2
c       to the input vector fin of dimensionality n, obtaining
c       the output vector fout of dimensionality 2*n.
c       The interpolation operator is supplied by the user in the
c       array riprfast, which (hopefully) has been generated by
c       a prior call to the subroutine prolfili (see)
c 
c                     Input parameters:
c 
c  riprfast - the array containing symmetrized versions of the operators
c         rinterp2, proj2, to be used by this entry for the application
c         the operator proj2 to arbitrary vectors
c  n - the number of nodes discretizing the interval [-1,1] with precision
c         eps and band-limit c; please note that there are n nodes xs,
c         and 2*n nodes xs2
c  fin - the (complex) vector of dimensionality n to which the operator
c         rinterp2 is to be applied
c 
c                     Output parameters:
c 
c  fout - the (complex) vector of dimensionality 2*n which is the result
c         of applying the operator rinterp2 to the input vector fin
c 
c                     Work arrays:
c 
c  w - must be at least n*4+20 complex *16 elements long
c 
c       . . . apply the symmetrized interpolation matrix to the
c             input vector fin
c 
        isinterp=riprfast(3)
        iainterp=riprfast(4)
c 
        iwp=1
        lwp=n+2
c 
        iwm=iwp+lwp
        lwm=n+2
c 
        ifp=iwm+lwm
        lfp=n+2
c 
        ifm=ifp+lfp
        lfm=n+2
c 
cccc        call prolinta(riprfast(isinterp),riprfast(iainterp),
        call cprolint(riprfast(isinterp),riprfast(iainterp),
     1      fin,fout,n,w(iwp),w(iwm),w(ifp),w(ifm) )
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine prolfiap(riprfast,n,fin,fout,w)
        implicit real *8 (a-h,o-z)
        save
        dimension fin(1),fout(1),riprfast(1),w(1)
c 
c       This entry applies the projection operator proj2
c       to the input vector fin of dimensionality 2*n, obtaining
c       the output vector fout of dimensionality n. The projection
c       operator is supplied by the user in the array riprfast,
c       which (hopefully) has been generated by a prior call to the
c       subroutine prolfili (see). Please note that this subroutine
c       has the companion entry prolinap (see below) that performs
c       the dual operation of interpolation.
c 
c                     Input parameters:
c 
c  riprfast - the array containing symmetrized versions of the operators
c         rinterp2, proj2, to be used by this entry for the application
c         the operator proj2 to arbitrary vectors
c  n - the number of nodes discretizing the interval [-1,1] with precision
c         eps and band-limit c; please note that there are n nodes xs,
c         and 2*n nodes xs2
c  fin - the vector of dimensionality n*2 to which the operator proj2
c         is to be applied
c 
c                     Output parameters:
c 
c  fout - the vector of dimensionality n which is the result of
c         applying the operator proj2 to the input vector fin
c 
c                     Work arrays:
c 
c  w - must be at least n*4+20 real *8 elements long
c 
c 
c       . . . apply the symmetrized filtering matrix to the
c             input vector fin
c 
        isproj=riprfast(1)
        iaproj=riprfast(2)
c 
        iwp=1
        lwp=n+2
c 
        iwm=iwp+lwp
        lwm=n+2
c 
        ifp=iwm+lwm
        lfp=n+2
c 
        ifm=ifp+lfp
        lfm=n+2
c 
        call prolproa(riprfast(isproj),riprfast(iaproj),
     1      fin,fout,n,w(iwp),w(iwm),w(ifp),w(ifm) )
        return
c 
c 
c 
c 
        entry prolinap(riprfast,n,fin,fout,w)
c 
c       This entry applies the interpolation operator rinterp2
c       to the input vector fin of dimensionality n, obtaining
c       the output vector fout of dimensionality 2*n. The interpolation
c       operator is supplied by the user in the array riprfast, which
c       (hopefully) has been generated by a prior call to the subroutine
c       prolfili (see). Please note that this subroutine has the
c       companion entry prolfiap (see above) that performs the dual
c       operation of filtering.
c 
c                     Input parameters:
c 
c  riprfast - the array containing symmetrized versions of the operators
c         rinterp2, proj2, to be used by this entry for the application
c         the operator proj2 to arbitrary vectors
c  n - the number of nodes discretizing the interval [-1,1] with precision
c         eps and band-limit c; please note that there are n nodes xs,
c         and 2*n nodes xs2
c  fin - the vector of dimensionality n to which the operator rinterp2
c         is to be applied
c 
c                     Output parameters:
c 
c  fout - the vector of dimensionality 2*n which is the result of
c         applying the operator rinterp2 to the input vector fin
c 
c                     Work arrays:
c 
c  w - must be at least n*4+20 real *8 elements long
c 
c       . . . apply the symmetrized interpolation matrix to the
c             input vector fin
c 
        isinterp=riprfast(3)
        iainterp=riprfast(4)
c 
        iwp=1
        lwp=n+2
c 
        iwm=iwp+lwp
        lwm=n+2
c 
        ifp=iwm+lwm
        lfp=n+2
c 
        ifm=ifp+lfp
        lfm=n+2
c 
        call prolinta(riprfast(isinterp),riprfast(iainterp),
     1      fin,fout,n,w(iwp),w(iwm),w(ifp),w(ifm) )
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine prolfili(ier,wavesnum,eps,ifrefine,
     1      xs,ws,xs2,ws2,riprfast,n,rinterp2,proj2,
     2      u,v,w,lenw,keep,lused)
        implicit real *8 (a-h,o-z)
c 
        save
        real *8 w(1),u(1),v(1),xs(1),ws(1),xs2(1),ws2(1),
     1          proj2(1),rinterp2(1),riprfast(1)
c 
c         This subroutine constructs matrices of filtering and
c         interpolation operators for band-limited functions.
c         The environment is not assumed to be periodic,
c         and the operators constructed by this subroutine are
c         filtering and interpolation operators in the weak sense
c         only (actually, the interpolation operator is a true
c         interpolation operator, in the strong sense; the filtering
c         operator filters things in a weak sense only). This
c         subroutine produces two operators: interpolation and
c         filtering.
c 
c    IMPORTANT NOTE ON BAND-LIMITS AND ACCURACY:
c 
c    The frequency behavior of the interpolation operator is simple:
c         it interpolates with respect to t all functions of the form
c 
c                     e^{i \cdot c \cdot x \cdot t},                    (1)
c 
c         with |x|< 1, and
c 
c         c = \pi * wavesnum.                                           (2)
c 
c         The interpolation is from the nodes xs (n of them)  to the
c         nodes xs2 (n \cdot 2 of them); all of the nodes {xs}, {xs2}
c         are on the interval [-1,1], and the precision of interpolation
c         is at least eps.
c 
c    The frequency behavior of the filtering operator is much more
c         complicated, and is as follows: given a function of the
c         form (1) with |x| < 3 (please note that "3" in the preceding
c         formula is not a typo) tabulated at the nodes {xs2}, the
c         filtering operator constructs a function f and evaluates
c         it at the nodes {xs}. The function f possesses the following
c         property:
c 
c 
c     | \int_{-1}^1 f(t) \cdot e^{i \cdot c \cdot y \cdot t} dt
c 
c     - \int_{-1}^1 f(t) \cdot e^{i \cdot c \cdot y \cdot t}            (2)
c 
c     \cdot e^{i \cdot c \cdot x \cdot t}  dt | < eps,
c 
c 
c         for all real y such that  |y|< 1. Furthermore, if the integral
c 
c 
c        \int_{-1}^1 f(t) \cdot e^{i \cdot c \cdot y \cdot t} dt        (3)
c 
c         is approximated via the quadrature formula with the nodes
c         {xs} and weights {ws}, the error of the resulting approximation
c         is also less than eps.
c 
c         IN OTHER WORDS, THE "FILTERING OPERATOR" IS A LOW-PASS FILTER
c         ON A FINITE INTERVAL ([-1,1]), AS IT HAPPENS), THAT IS ONLY A
c         FILTER IN A WEAK SENSE, WHEN TESTED AGAINST TEST FUNCTIONS OF
c         THE FORM (1), WITH |x| < 1; THE SAID LOWPASS FILTER ONLY WORKS
c         FOR BAND-LIMITED FUNCTIONS, WITH THE BAND-LIMIT c \cdot 3, AND
c         C DEFINED BY (2). ALSO, PLEASE NOTE THAT THE NUMBER N OF NODES
c         IS DETERMINED BY THIS SUBROUTINE, AND DEPENDS ON BOTH WAVESNUM
c         AND EPS, AS ONE WOULD EXPECT.
c 
c         The filtering and interpolation operators are returned
c         by this subroutine in two forms: matrices rinterp2, proj2,
c         and the array riprfast. The matrices rinterp2, proj2 are
c         truly matrices representing the operators of interpolation
c         and filtering, respectively (see detailed descriptions
c         below). The user can be apply them to vectors of function
c         values (for example, via the subroutine prolmav2 present
c         in this package; on the other hand, the user might be
c         able to write his own subroutine applying a matrix to a
c         vector).
c 
c         IMPORTANT NOTE: THE SIZES OF OUTPUT ARRAYS PRODUCED BY THIS
c         SUBROUTINE DEPEND ON THE OUTPUT PARAMETER N, UNKNOWN BEFORE
c         THE SUBROUTINE HAS BEEN CALLED. THIS IS A SIGNIFICANT
c         DRAWBACK, SOMEWHAT ALLEVIATED BY THE FACT THAT IN ANTICIPATED
c         APPLICATIONS, EVERYTHING IS FAIRLY SMALL-SCALE. THE FOLLOWING
c         SIZES ARE ALWAYS SUFFICIENT:
c 
c      xs,ws :     c+200 real *8 elements long each
c      xs2,ws2 :   2*c+400 real *8 elements long each
c      riprfast:   4*(c+200)**2+100 real *8 elements long
c      rinterp2:   2*(c+200)**2 real *8 elements long
c      proj2:      2*(c+200)**2 real *8 elements long
c      u,v:        (c+200)**2 real *8 elements long each
c      w  -  make it long, and tell the subroutine honestly how long
c                  it is, via the parameter lenw.
c 
c         IN THE ABOVE, c IS GIVEN BY THE FORMULA
c 
c                    c=\pi * wavesnum
c 
c         The array riprfast contains symmetrized versions
c         of matrices rinterp2, proj2; these versions are used by
c         the entries prolfiap, prolinap of the subroutine prolfiap
c         (see) to apply the filtering and interpolation operators
c         for about 1/2 of the cost of applying the matrices rinterp2,
c         proj2 directly. The reason for the existence of the array
c         riprfast and the subroutine prolfiap is the fact that the
c         matrices rinterp2, proj2 have a symmetry in them; thus,
c         a factor of two can be saved compared to their direct
c         application to vectors.
c 
c         To look at the matters in finer detail, the interpolation
c         operator is represented by the matrix rinterp2, converting
c         a band-limited function on the n-point prolate grid on
c         the interval [-1,1] (nodes of this grid are returned to the
c         user in the array xs) into a function on the 2*n-point grid
c         on the interval [-1,1] (returned to the user in the array
c         xs2). The latter grid is obtained by concatenating two
c         n-point prolate grids on the intervals [-1,0], [0,1].
c 
c         Similarly, the filtering operator is represented by the
c         matrix proj2, converting a band-limited function on the
c         2*n-point grid on the interval [-1,1] (returned to the
c         user in the array xs2), into a filtered function at the
c         nodes of an  n-point grid on the interval [-1,1] (returned
c         to the user in the array xs).
c 
c 
c                      Input parameters:
c 
c  wavesnum - the number of wavelengths on the interval [-1,1]
c         at the band-limit of the functions TO BE INTERPOLATED;
c         at the band-limit of the functions to be filtered, there
c         are 3*wavesnum wavelength on the interval [-1,1]
c  eps - the precision to which filtering and interpolation are to
c         be performed. Generally the accuracy of the interpolation
c         will be eps (a little better, as a matter of fact); the
c         accuracy of filtering will be considerably (1-2 digits)
c         better.
c  ifrefine - the parameter telling the subroutine which type of
c       prolate nodes and weights to generate:
c    ifrefine=0 will cause the nodes based on the division theorem
c       to be generated (i.e. the resulting nodes will be th eroots
c       of the appropriately chosen prolate function)
c    ifrefine=1 will cause the optimal prolate nodes to be generated.
c       PLEASE NOTE THAT THE TWO SETS OF NODES ARE VERY CLOSE TO EACH
c       OTHER. AS LONG AS THEY ARE USED STRICTLY FOR INTERPOLATION,
c       THERE IS NO SERIOUS ADVANTAGE FOR EITHER SET. HOWEVER, WHEN
c       USED AS QUADRATURE NODES, THE OPTIMAL NODES GAIN A COUPLE
c       OF DIGITS ON THE ONES PRODUCED WITH IFREFINE=0. On the other
c       hand, the accuracy of both quadratures is better than eps
c       (see above), and ifrefine=1 causes a MUCH greater expenditure
c       of CPU time than ifrefine=0.
c  lenw - the amount of space (n real *8 words) provided by the user
c         in the work array w
c 
c                      Output parameters:
c 
c  ier - error return code
c 
c    ier=0 means successful execution
c    ier=16 or ier=8 means that the amount of memory provided by the
c         user in the array w (and supplied to the subroutine in the
c         parameter lenw) is insufficient. Obviously, it is a fatal
c         error.
c  xs - the n prolate nodes on the interval [-1,1]; the functions will
c         be interpolated from these nodes onto nodes xs2, and filtered
c         from the nodes xs2 to the nodes xs.
c  ws - quadrature weights corresponding to the nodes {xs}
c  xs2 - the nodes on the interval [-,1], by concatenating sets on n
c         prolate nodes on the intervals [-1,0], [0,1]; the functions
c         will be filtered from the nodes {xs2} onto nodes xs, and
c         interpolated from the nodes xs to the nodes {xs2}.
c  ws2 - quadrature weights corresponding to the nodes {xs2}
c  riprfast - the array containing symmetrized versions of the operators
c         rinterp2, proj2, to be used by the subroutine prolfiap (see)
c         for the application of these operators to arbitrary vectors
c  rinterp2 - the 2*n \times n matrix interpolating band-limited functions
c         from the nodes xs onto nodes xs2
c  proj2 - the n \times 2*n matrix filtering band-limited functions
c         from the nodes xs2 onto nodes xs
c  n - the number of nodes discretizing the interval [-1,1] with precision
c         eps and band-limit c; please note that there are n nodes xs,
c         and 2*n nodes xs2; also note that n is always even.
c  u - the matrix converting the n coefficients of a prolate expansion
c       on the interval [-1,1] into the values of the said expansion
c       at prolate nodes; provided for user's entertainment
c  v - the matrix converting the values of a prolate expansion
c         at prolate nodes into  the n coefficients of the said prolate
c         expansion; provided for user's entertainment. Needless to say,
c         u,v are inverses of each other.
c  keep - the number of elements in the array w that have to be
c       unchanged between the call to this subroutine and the subsequent
c       calls to the subroutine proleva and its relatives (see)
c  lused - integer number telling the user how much space in the
c       array w (in real *8 words) the subroutine actually used.
c  w - contains data that can be used by the subroutine proleva and its
c         relatives to evaluate prolate functions at arbitrary points on
c         the interval [-1,1], and for whatever other needs the subroutines
c         in the prolcrea might have; otherwise, it is a work array; do not
c         skimp on its size!
c  keep - the number of elements in the array w that have to be
c       unchanged between the call to this subroutine and the subsequent
c       calls to the subroutine proleva and its relatives (see)
c  lused - integer number telling the user how much space in the
c       array w (in real *8 words) the subroutine has actually used.
c 
c                      Work arrays:
c 
c  w - must be fairly large (its size grows quadratically with wavesnum)
c 
c 
        ier=0
c 
        done=1
        pi=atan(done)*4
c 
        c=wavesnum*pi
c 
c       construct the matrices of projection and interpolation
c 
        call prolinfi(jer,c,eps,ifrefine,
     1    xs,ws,xs2,ws2,rinterp2,proj2,n,u,v,w,lenw,keep,lused)
c 
        if(jer .eq. 0) goto 1400
        ier=16
        return
 1400 continue
c 
c        allocate memory for the remainder of the procedure, and
c        store the data where they have to be stored
c 
        ixsr=keep+2
        lxsr=n+2
c 
        iwsr=ixsr+lxsr
        lwsr=n+2
c 
        iproj=iwsr+lwsr
        lproj=n**2*2+10
c 
        irinterp=iproj+lproj
        lrinterp=n**2*2+10
c 
        lused2=irinterp+lrinterp
  
        call prinf('lused2=*',lused2,1)
        call prinf('while lused=*',lused,1)
        call prinf('and keep=*',keep,1)
  
        if(lused2 .gt. lused) lused=lused2
c 
        if(lused .lt. lenw) goto 1600
        ier=8
        return
 1600 continue
c 
        call prolarrm(xs2,w(ixsr),n)
        call prolarrm(ws2,w(iwsr),n)
        call prolarrm(rinterp2,w(irinterp),n**2)
        call prolarrm(proj2,w(iproj),n**2)
  
cccc          stop
c 
c       construct the projection and interpolation operators for
c       the pair of subintervals, starting with the projection and
c       interpolation operators for the single right subinterval,
c       generated by the prior call to prolinfi
c 
        call prolinpa(xs,ws,w(ixsr),w(iwsr),n,w(iproj),w(irinterp),
     1      xs2,ws2,proj2,rinterp2)
c 
c       allocate memory for the final restructuring and symmetrization
c       of the projection and interpolation operators
c 
        iww=keep+2
        lww=n*n+10
c 
        lused3=iww+lww
        if(lused3 .gt. lused) lused=lused3
  
        call prinf('and lused3=*',lused3,1)
c 
        if(lused .lt. lenw) goto 1800
        ier=8
        return
 1800 continue
c 
        isproj=11
        lsproj=n/2*n+10
c 
        iaproj=isproj+lsproj
        laproj=n/2*n+10
c 
        isinterp=iaproj+laproj
        lsinterp=n/2*n+10
c 
        iainterp=isinterp+lsinterp
        lainterp=n/2*n+10
c 
        lripr=iainterp+lainterp
        call prinf('lripr=*',lripr,1)
        n22p=2*n**2+100
c 
        call prinf('and n22p=*',n22p,1)
c 
        riprfast(1)=isproj+0.1
        riprfast(2)=iaproj+0.1
        riprfast(3)=isinterp+0.1
        riprfast(4)=iainterp+0.1
c 
c       construct the symmetrized projection and interpolation
c       matrices, in order to save a factor of 2 in applying these
c       matrices to vectors
c 
       call prolsym(proj2,riprfast(isproj),riprfast(iaproj),
     1     w(iww),n,rinterp2,riprfast(isinterp),riprfast(iainterp))
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine prolinta(sinter,ainter,fin,fout,n,
     1      wp,wm,fp,fm)
        implicit real *8 (a-h,o-z)
        save
        real *8 sinter(n,n/2),ainter(n,n/2),fin(1),fout(1),
     1           wp(1),wm(1),fp(1),fm(1)
c 
c       symmetrize and antisymmetrize the input vector fin
c 
  
        do 2400 i=1,n/2
c 
        wp(i)=fin(i)+fin(n-i+1)
        wm(i)=fin(i)-fin(n-i+1)
 2400 continue
c 
c       apply the symmetrized and antisymmetrized projection
c       matrices to the symmetrized and antisymmetrized input
c       vectors, respectively
c 
        call prolmav2(sinter,wp,n,n/2,fp)
        call prolmav2(ainter,wm,n,n/2,fm)
c 
c 
c 
        do 2600 i=1,n
c 
        fout(i)=fp(i)+fm(i)
        fout(2*n-i+1)=fp(i)-fm(i)
 2600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine prolproa(sproj,aproj,fin,fout,n,
     1      wp,wm,fp,fm)
        implicit real *8 (a-h,o-z)
        save
        real *8 sproj(n/2,n),aproj(n/2,n),fin(1),fout(1),
     1           wp(1),wm(1),fp(1),fm(1)
c 
c       symmetrize and antisymmetrize the input vector fin
c 
        do 2400 i=1,n
c 
        wp(i)=fin(i)+fin(2*n-i+1)
        wm(i)=fin(i)-fin(2*n-i+1)
 2400 continue
c 
c       apply the symmetrized and antisymmetrized projection
c       matrices to the symmetrized and antisymmetrized input
c       vectors, respectively
c 
        call prolmav2(sproj,wp,n/2,n,fp)
        call prolmav2(aproj,wm,n/2,n,fm)
c 
c       combine the symmetrized vectors appropriately
c 
        do 2600 i=1,n/2
c 
        fout(i)=fp(i)+fm(i)
        fout(n-i+1)=fp(i)-fm(i)
 2600 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine prolsym(proj2,sproj,aproj,ww,n,
     1      rinterp2,sinterp,ainterp)
        implicit real *8 (a-h,o-z)
        save
        real *8 proj2(n,n*2),ww(n,n),sproj(n/2,n),aproj(n/2,n),
     1      rinterp2(n*2,n),sinterp(n,n/2),ainterp(n,n/2)
c 
c           first, deal with the matrix proj2 as it deserves
c 
c 
c       . . . symmetrize the matrix proj2
c 
        do 2400 i=1,n
        do 2200 j=1,n
c 
        ww(i,j)=proj2(i,j)+proj2(i,2*n-j+1)
 2200 continue
 2400 continue
c 
        do 2800 i=1,n/2
        do 2600 j=1,n
c 
        sproj(i,j)=(ww(i,j)+ww(n-i+1,j))/4
 2600 continue
 2800 continue
c 
c       antisymmetrize the matrix proj2
c 
        do 3400 i=1,n
        do 3200 j=1,n
c 
        ww(i,j)=proj2(i,j)-proj2(i,2*n-j+1)
 3200 continue
 3400 continue
c 
        do 3800 i=1,n/2
        do 3600 j=1,n
c 
        aproj(i,j)=(ww(i,j)-ww(n-i+1,j))/4
 3600 continue
 3800 continue
c 
c           now, deal with the matrix rinterp2
c 
c       . . . symmetrize the matrix rinterp2
c 
        do 4400 i=1,n
        do 4200 j=1,n
c 
        ww(i,j)=rinterp2(i,j)+rinterp2(2*n-i+1,j)
 4200 continue
 4400 continue
c 
        do 5800 i=1,n
        do 5600 j=1,n/2
c 
        sinterp(i,j)=(ww(i,j)+ww(i,n-j+1))/4
 5600 continue
 5800 continue
c 
c       antisymmetrize the matrix rinterp2
c 
        do 6400 i=1,n
        do 6200 j=1,n
c 
        ww(i,j)=rinterp2(i,j)-rinterp2(2*n-i+1,j)
 6200 continue
 6400 continue
c 
        do 7800 i=1,n
        do 7600 j=1,n/2
c 
        ainterp(i,j)=(ww(i,j)-ww(i,n-j+1))/4
 7600 continue
 7800 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine prolinpa(xs,ws,xsr,wsr,n,proj,rinterp,
     1      xs2,ws2,proj2,rinterp2)
        implicit real *8 (a-h,o-z)
        save
        real *8 proj(n,n),proj2(n,n*2),rinterp2(2*n,n),
     1      xs(1),ws(1),xsr(1),wsr(1),xs2(1),ws2(1),rinterp(n,n)
c 
c        This subroutine converts the one-sided interpolation and
c        filtering matrices (produced by the subroutine prolinfi)
c        into full two-sided matrices. It is simply a reordering
c        procedure; it contains logic but no analysis.
c 
c                      Input parameters:
c 
c  xs - the n prolate nodes on the interval [-1,1]
c  ws - the weights corresponding to the nodes xs
c  xsr - the n prolate nodes on the interval [0,1]
c  wsr - the weights corresponding to the nodes xsr
c  n - the numbers of elements in the arrays xs,ws,xsr,wsr
c  proj - the n*n-matrix converting functions on the interval [0,1]
c        tabulated at the n prolate nodes into functions on the
c        interval [-1,1], appropriately filtered and tabulated at the
c        n prolate nodes
c  rinterp - the n*n-matrix converting functions on the interval [-1,1]
c        tabulated at the n prolate nodes into functions on the
c        interval [0,1], tabulated at the n prolate nodes
c 
c                      Output parameters:
c 
c  xs2 - 2*n nodes on the interval [-1,1] obtained by concatenating
c        n prolate nodes on the interval [-1,0] with n prolate nodes
c        on the interval [0,1]
c  ws2 - the weights corresponding to the nodes xs2
c  proj2 - the n \times n*2 -matrix converting functions with band-limit
c        2*c on the interval [-1,1] tabulated at the nodes xs2 into
c        functions with band-limit c on the same interval [-1,1],
c        appropriately filtered and tabulated at the n prolate nodes xs
c  rinterp - the 2*n \time n - matrix interpolating functions with
c        band-limit c on the interval [-1,1] tabulated at the n prolate
c        nodes xs the nodes xs2
c 
c        . . . construct the discretizations of the two subintervals
c              (in one array)
c 
        do 1100 i=1,n
c 
        xs2(i)=-xsr(n-i+1)
        xs2(n+i)=xsr(i)
c 
        ws2(i)=wsr(i)
        ws2(n+i)=wsr(i)
 1100 continue
c 
c       construct the big projection matrix
c 
        do 1400 i=1,n
        do 1200 j=1,n
c 
        proj2(i,n+j)=proj(i,j)
 1200 continue
 1400 continue
c 
        do 2400 i=1,n
        do 2200 j=1,n
c 
        proj2(i,j)=proj2(n-i+1,2*n-j+1)
 2200 continue
 2400 continue
c 
c       construct the big interpolation matrix
c 
        do 3400 i=1,n
        do 3200 j=1,n
c 
        rinterp2(n+i,j)=rinterp(i,j)
 3200 continue
 3400 continue
c 
        do 4400 i=1,n
        do 4200 j=1,n
c 
        rinterp2(i,j)=rinterp2(2*n-i+1,n-j+1)
 4200 continue
 4400 continue
c 
  
        return
        end
c 
c 
c 
c 
c 
        subroutine prolinfi(ier,c,eps,ifrefine,
     1    xs,ws,xsr,wsr,rinterp,proj,n,u,v,w,lenw,keep,lused)
        implicit real *8 (a-h,o-z)
        save
        real *8 u(1),v(1),proj(1),rinterp(1),
     1      w(1),xsr(1),wsr(1),xs(1),ws(1)
c 
c         This subroutine constructs matrices of filtering and
c         interpolation operators for band-limited functions.
c         The environment is not assumed to be periodic,
c         and the operators constructed by this subroutine are
c         filtering and interpolation operators in the weak sense
c         only (actually, the interpolation operator is a true
c         interpolation operator, in a strong sense; the filtering
c         operator filters things in a weak sense only).
c 
c         Please note that this subroutine produces interpolation
c         coefficients that interpolate from an n-point prolate grid on
c         the interval [-1,1] onto a prolate grid (of the same size) on
c         the interval [0,1]; the filtering coefficients it produces
c         filter from n prolate nodes on the interval [0,1] onto n
c         prolate nodes on the interval [-1,1].
c 
c 
c                      Input parameters:
c 
c  c - the band-limit
c  eps - the precision to which filtering and interpolation are to
c         be performed
c  ifrefine - the parameter telling the subroutine which type of
c       prolate nodes and weights to generate:
c    ifrefine=0 will cause the nodes based on the division theorem
c       to be generated (i.e. the resulting nodes will be th eroots
c       of the appropriately chosen prolate function)
c    ifrefine=1 will cause the optimal prolate nodes to be generated.
c       PLEASE NOTE THAT THE TWO SETS OF NODES ARE VERY CLOSE TO EACH
c       OTHER. AS LONG AS THEY ARE USED STRICTLY FOR INTERPOLATION,
c       THERE IS NO SERIOUS ADVANTAGE FOR EITHER SET. HOWEVER, WHEN
c       USED AS QUADRATURE NODES, THE OPTIMAL NODES GAIN A COUPLE
c       OF DIGITS ON THE ONES PRODUCED WITH IFREFINE=0; on the other
c       hand, ifrefine=1 causes a much greater expenditure of CPU time
c       than ifrefine=0.
c  lenw - the amount of space (n real *8 words) provided by the user
c         in the work array w
c 
c                      Output parameters:
c 
c  xs - the nodes on the interval [-1,1]; the functions will be
c         interpolated from these nodes onto nodes xsr, and filtered
c         from the nodes xsr to the nodes xs.
c  ws - quadrature weights corresponding to the nodes {xs}
c  xsr - the nodes on the interval [0,1]; the functions will be
c         filtered from these nodes onto nodes xsr, and interpolated
c         from the nodes xsr to the nodes xs.
c  wsr - quadrature weights corresponding to the nodes {xsr}
c  rinterp - the n \times n matrix interpolating band limited functions
c         from the nodes xsr onto nodes xs
c  proj - the n \times n matrix filtering band-limited functions
c         from the nodes xs onto nodes xsr
c  n the number of nodes discretizing the interval [-1,1] with precision
c         eps and band-limit c
c  u - the matrix converting the n coefficients of a prolate expansion
c       on the interval [-1,1] into the values of the said expansion
c       at prolate nodes.
c  v - the matrix converting the values of a prolate expansion
c       at prolate nodes into  the n coefficients of the said prolate
c       expansion. Needless to say, u,v are inverses of each other.
c  keep - the number of elements in the array w that have to be
c       unchanged between the call to this subroutine and the subsequent
c       calls to the subroutine proleva and its relatives (see)
c  lused - integer number telling the user how much space in the
c       array w (in real *8 words) the subroutine actually used.
c 
c                      Work arrays:
c 
c  w - must be fairly large (its size grows quadratically with c)
c 
        ier=0
c 
        call prolexps(jer,c,eps,ifrefine,n,xs,ws,u,v,w,lenw,keep,lused)
c 
        if(jer .eq. 0) goto 1400
c 
        call prinf('after prolexps, jer=*',jer,1)
        ier=32
        return
 1400 continue
c 
c       allocate memory for the construction of the matrix of inner
c       products
c 
        nlegend=c+100
c 
        iprods=keep+2
        lprods=n**2+2
c 
        ifuns1=iprods+lprods
        lfuns1=nlegend*n+10
c 
        ifuns2=ifuns1+lfuns1
        lfuns2=nlegend*n+10
c 
        its=ifuns2+lfuns2
        lts=nlegend+2
c 
        its2=its+lts
        lts2=nlegend+2
c 
        iwhts=its2+lts2
        lwhts=nlegend+2
c 
        lused=iwhts+lwhts
c 
        if(lused .lt. lenw) goto 1600
        ier=16
        return
 1600 continue
c 
c       construct the matrix of inner products on the interval [0,1]
c       between prolate functions and those same functions scaled
c       to the interval [0,1]
c 
        call prodseva(c,nlegend,w,w(ifuns1),w(ifuns2),w(iprods),
     1      n,w(its),w(its2),w(iwhts) )
c 
c       allocate memory for the remainder of the subroutine
c 
        iwork=iprods+lprods
        lwork=n*n+2
c 
c       construct the matrix of the projection operator
c 
        call promama(u,w(iprods),n,w(iwork) )
        call promama(w(iwork),v,n,proj)
c 
c       construct the nodes and weights on the subintervals
c 
        do 2200 i=1,n
        xsr(i)=(xs(i)+1)/2
        wsr(i)=ws(i)/2
 2200 continue
c 
c        scale the projection matrix by a factor of 2
c 
        do 2400 i=1,n**2
c 
        proj(i)=proj(i)/2
 2400 continue
c 
c       construct the matrix of the interpolation operator
c 
        call prointemat(w,xsr,n,w(iwork),u,v,rinterp)
c 
        return
        end
c 
c 
c 
c 
c 
c 
        subroutine prointemat(w,xsr,n,a,u,v,rinterp)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),xsr(1),w(1),u(1),v(1),rinterp(100000)
c 
c        construct the matrix of values of prolate functions at the
c        prolate nodes on the interval [0,1]
  
  
        do 1600 i=1,n
        do 1400 j=1,n
c 
        call prolevv(i-1,xsr(j),w,val)
c 
        a(j,i)=val
c 
 1400 continue
 1600 continue
c 
c       multiply the matrix a by the decompusition matrix u
c       from the right, obtaining the interpolation matrix
c 
        call promama(a,v,n,rinterp)
  
        call prin2('in prouintemat, rinterp=*',rinterp,100)
  
  
        return
        end
c 
c 
c 
c 
        subroutine prodseva(c,n,w,funs1,funs2,prods,nprods,
     1      ts,ts2,whts)
        implicit real *8 (a-h,o-z)
        save
        real *8 w(1),funs1(n,1),funs2(n,1),prods(nprods,1)
c 
        real *8 ts(1),ts2(1),whts(1)
c 
c       discretize the intervals [-1,1], [0,1]
c 
        itype=1
        call legeexps(itype,n,ts,u,v,whts)
c 
        do 1200 i=1,n
c 
        ts2(i)=(ts(i)+1)/2
 1200 continue
c 
c       evaluate prolate functions at both sets of points
c 
        do 1600 i=1,nprods
c 
        do 1400 j=1,n
c 
        call prolevv(i-1,ts(j),w,val)
        funs1(j,i)=val
c 
        call prolevv(i-1,ts2(j),w,val)
        funs2(j,i)=val
 1400 continue
 1600 continue
c 
c       construct the matrix of inner products between thge
c       obtained functions
c 
        do 2400 i=1,nprods
        do 2200 j=1,nprods
c 
        call prolsca(funs1(1,i),funs2(1,j),n,prod,whts)
c 
        prods(j,i)=prod
 2200 continue
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
        subroutine prolsca(x,y,n,prod,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension y(1),whts(1)
        real *8 x(1),prod
c 
        prod=0
        do 1800 i=1,n
        prod=prod+x(i)*y(i)*whts(i)
 1800 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine prolmav2(a,x,n,m,y)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),x(m),y(n)
c 
cccc        call prin2('in prolmav2, a=*',a,n*m)
cccc        call prin2('in prolmav2, x=*',x,m)
        do 1400 i=1,n
        d=0
        do 1200 j=1,m
        d=d+a(i,j)*x(j)
 1200 continue
  
cccc        call prin2('in prolmav2, d=*',d,1)
        y(i)=d
 1400 continue
  
        return
        end
c 
c 
c 
c 
c 
       subroutine promamaa(a,b,n,c)
       implicit real *8 (a-h,o-z)
        save
       dimension a(n,n),b(n,n),c(n,n)
c 
       do 1600 i=1,n
       do 1400 j=1,n
       d=0
       do 1200 k=1,n
c 
       d=d+a(i,k)*b(j,k)
 1200 continue
c 
       c(i,j)=d
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
       subroutine promama(a,b,n,c)
       implicit real *8 (a-h,o-z)
        save
       dimension a(n,n),b(n,n),c(n,n)
c 
       do 1600 i=1,n
       do 1400 j=1,n
       d=0
       do 1200 k=1,n
c 
       d=d+a(i,k)*b(k,j)
 1200 continue
c 
       c(i,j)=d
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
        subroutine cprolint(sinter,ainter,fin,fout,n,
     1      wp,wm,fp,fm)
        implicit real *8 (a-h,o-z)
        save
        real *8 sinter(n,n/2),ainter(n,n/2)
        complex * 16 fin(1),fout(1),wp(1),wm(1),fp(1),fm(1)
c 
c    ATTENTION: THIS IS THE COMPLEX VERSION OF THE REAL SUBROUTINE
c       PROLINTA, IDENTICAL TO THE LATTER IN ALL ESSENTIAL DETAILS
c 
c       symmetrize and antisymmetrize the input vector fin
c 
        do 2400 i=1,n/2
c 
        wp(i)=fin(i)+fin(n-i+1)
        wm(i)=fin(i)-fin(n-i+1)
 2400 continue
c 
c       apply the symmetrized and antisymmetrized projection
c       matrices to the symmetrized and antisymmetrized input
c       vectors, respectively
c 
        call cprolma2(sinter,wp,n,n/2,fp)
        call cprolma2(ainter,wm,n,n/2,fm)
c 
        do 2600 i=1,n
c 
        fout(i)=fp(i)+fm(i)
        fout(2*n-i+1)=fp(i)-fm(i)
 2600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cprolpro(sproj,aproj,fin,fout,n,
     1      wp,wm,fp,fm)
        implicit real *8 (a-h,o-z)
        save
        real *8 sproj(n/2,n),aproj(n/2,n)
        complex *16 fin(1),fout(1),wp(1),wm(1),fp(1),fm(1)
c 
c    ATTENTION: THIS IS THE COMPLEX VERSION OF THE REAL SUBROUTINE
c       PROLPROA, IDENTICAL TO THE LATTER IN ALL ESSENTIAL DETAILS
c 
c       symmetrize and antisymmetrize the input vector fin
c 
        do 2400 i=1,n
c 
        wp(i)=fin(i)+fin(2*n-i+1)
        wm(i)=fin(i)-fin(2*n-i+1)
 2400 continue
c 
c       apply the symmetrized and antisymmetrized projection
c       matrices to the symmetrized and antisymmetrized input
c       vectors, respectively
c 
        call cprolma2(sproj,wp,n/2,n,fp)
        call cprolma2(aproj,wm,n/2,n,fm)
c 
c       combine the symmetrized vectors appropriately
c 
        do 2600 i=1,n/2
c 
        fout(i)=fp(i)+fm(i)
        fout(n-i+1)=fp(i)-fm(i)
 2600 continue
        return
        end
  
c 
c 
c 
c 
c 
        subroutine cprolma2(a,x,n,m,y)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m)
        complex * 16 x(m),y(n),d
c 
        do 1400 i=1,n
        d=0
        do 1200 j=1,m
           d=d+a(i,j)*x(j)
 1200 continue
c 
        y(i)=d
 1400 continue
  
        return
        end
