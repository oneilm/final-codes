        implicit real *8 (a-h,o-z)
        real *8 xs1(10 000),ws1(10 000),xs2(10 000),
     1      ws2(10 000),rinterp(100 000),filter(100 000),
     2      w1(500 000),
     3      xs3(1000),ws3(1000),u3(100 000),v3(100 000),
     4      w3(1 000 000),riprfast(100 000)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, '    '
        PRINT *, 'ENTER wavesnum - not used in this code'
        READ *,wavesnum
        CALL PRIN2('wavesnum=*',wavesnum,1)
  
c 
c       initialize the filtering/interpolation subroutine
c 
        done=1
        pi=atan(done)*4
        wavesnum1=3
        wavesnum2=5
  
  
        c1=wavesnum1*pi
        c2=wavesnum2*pi
  
ccc        c2=40
c 
        eps1=1.0d-10
        eps2=1.0d-10
  
        call prin2('in main program, c1=*',c1,1)
        call prin2('in main program, c2=*',c2,1)
  
  
  
        call fileflush(13)
c 
        lenw1=500 000
c 
  
        ifrefine=1
  
cccc        call prolfltinte(ier,c1,eps1,xs1,ws1,n1,
        call prolfltinte(ier,wavesnum1,eps1,xs1,ws1,n1,
     1      wavesnum2,eps2,xs2,ws2,n2,ifrefine,rinterp,filter,
     2      w1,lenw1,riprfast)
c 
c       test the filtering and the interpolation
c 
        lenw3=500 000
        call prolexps(jer1,c1,eps1,ifrefine,n3,xs3,ws3,
     1      u3,v3,w3,lenw3,keep1,lused1)
  
c 
cccc        call cfilttst(c1,xs3,ws3,n3,c2,xs2,ws2,n2,filter,v3,
        call filttst(c1,xs3,ws3,n3,c2,xs2,ws2,n2,filter,v3,
     1       w3,rinterp,riprfast)
  
  
        call prinf('and by the way, n1=*',n1,1)
        call prinf('and by the way, n2=*',n2,1)
  
        stop
        end
  
c 
c 
c 
c 
c 
        subroutine cfilttst(c1,xs1,ws1,n1,
     1      c2,xs2,ws2,n2,filter,v1,w1,rinterp,riprfast)
        implicit real *8 (a-h,o-z)
        save
        real *8 xs1(1),ws1(1),ws2(1),xs2(1),
     1      filter(n1,n2),rinterp(n2,n1),
     2      v1(1),t22(10 000),
     3      coefs(10 000),w1(1),coefslege(10000),
     4      ww(100000),fsrea1(10 000),fsima1(10 000),
     5      fsrea2(10 000),fsima2(10 000)
c 
        complex *16 fsin(1000),fsout(1000),diffs(1000),
     1      fs2(1000),f22(10 000),f23(10 000),fs1(1000),
     2      fs3(1000),errs(1000),cd,cdd1,cdd2,ima,der
c 
        data ima/(0.0d0,1.0d0)/
c 
c       construct the function to be interpolated
c 
        call cfuntestin(c2 * 0.9)
c 
        do 1200 i=1,n2
c 
        call cfuntest(xs2(i),fsin(i))
 1200 continue
c 
        call cprolc1c2filter(riprfast,fsin,fsout,ww)
  
  
        call prin2('and fsout=*',fsout,n1*2)
  
  
cccc        stop
c 
        do 1400 i=1,n1
c 
        call cfuntest(xs1(i),fs2(i))
c 
        diffs(i)=fs2(i)-fsout(i)
c 
 1400 continue
c 
        call prin2('and fs2=*',fs2,n1*2)
        call prin2(
     1   'differences between filtered and original functions=*',
     2      diffs,n1*2)
  
  
cccc        stop
c 
c       plot both of them things
c 
        do 1450 i=1,n1
c 
        fsrea2(i)=fs2(i)
        fsrea1(i)=fsout(i)
c 
        fsima2(i)=fs2(i)*ima
        fsima1(i)=fsout(i)*ima
 1450 continue
c 
        iw=21
        call quagraph2(iw,xs1,fsrea2,n1,3,xs1,fsrea1,n1,3,
     1      'filtered real part vs. original real part*')
c 
        iw=22
        call quagraph2(iw,xs1,fsima2,n1,3,xs1,fsima1,n1,3,
     1  'filtered imaginary part vs. original imaginary part*')
  
cccc        stop
c 
c       evaluate inner products of the filtered function
c       with a bunch of complex exponentials, and compare
c       their inner products with those of the unfiltered
c       function
c 
        ntest=100
        h=c1/(ntest-1)
        do 1800 ii=1,ntest
c 
        cc=(ii-1)*h
  
cccc        call prin2('cc=*',cc,1)
  
        cdd1=0
        do 1600 i=1,n1
c 
        cd=sin(cc*xs1(i))
        cdd1=cdd1+ws1(i)*cd*fsout(i)
 1600 continue
  
cccc        call prin2('and dd1=*',dd1,1)
  
  
        cdd2=0
        do 1700 i=1,n2
c 
        cd=sin(cc*xs2(i))
        cdd2=cdd2+ws2(i)*cd*fsin(i)
 1700 continue
  
cccc        call prin2('and dd2=*',dd2,1)
cccc        call prin2('and dd2-dd1=*',dd2-dd1,1)
  
        errs(ii)=cdd2-cdd1
  
 1800 continue
  
  
        call prin2('errors in projection are=*',errs,ntest*2)
  
  
cccc          stop
  
c 
c       construct a much finer discretization of the functiun
c       fsout, and plot it
c 
cccc        call prolmat(v1,fsout,n1,n1,coefs)
        call cprolma2(v1,fsout,n1,n1,coefs)
  
        call prin2('and coefs=*',coefs,n1*2)
c 
        call cprotole(coefs,n1,w1,coefslege,nlege)
  
        call prin2('and coefslege=*',coefslege,nlege*2)
  
cccc        stop
  
  
  
        ntest=1000
  
  
        call prinf('ntest as created=*',ntest,1)
        done=1
        h=2*done/(ntest-1)
        do 2200 i=1,ntest
c 
        t22(i)=-done+(i-1)*h
c 
        call legecFDE(t22(i),f22(i),der,coefslege,Nlege-1)
  
  
cccc      SUBROUTINE legecFDE(X,VAL,der,PEXP,N)
  
c 
        call cfuntest(t22(i),f23(i))
  
  
c 
        fsrea2(i)=f22(i)
        fsrea1(i)=f23(i)
c 
        fsima2(i)=f22(i)*ima
        fsima1(i)=f23(i)*ima
  
  
  
 2200 continue
  
        call prinf('ntest=*',ntest,1)
  
cccc        call prin2('f22=*',f22,ntest*2)
cccc        call prin2('f23=*',f23,ntest*2)
cccc        call prin2('fsrea1=*',fsrea1,ntest)
cccc        call prin2('fsrea2=*',fsrea2,ntest)
c 
        iw=23
        call quagraph2(iw,t22,fsrea2,ntest,3,t22,fsrea1,ntest,3,
     1      'real part of filtered vs. original after resampling*')
c 
c 
        iw=24
        call quagraph2(iw,t22,fsima2,ntest,3,t22,fsima1,ntest,3,
     1  'imaginary part of filtered vs. original after resampling*')
c 
        call prinf('and nlege=*',nlege,1)
  
  
  
cccc        stop
c 
c       test the interpolation operator
c 
  
        call cfuntestin(c1 * 0.9)
  
  
c 
        do 3200 i=1,n1
c 
        call cfuntest(xs1(i),fsin(i))
 3200 continue
c 
  
  
        call cprolc1c2interp(riprfast,fsin,fsout,ww)
  
  
        call prin2('after prolproc1c2, fsout=*',fsout,n2*2)
  
  
cccc        stop
c 
        do 3400 i=1,n2
c 
        call cfuntest(xs2(i),fs2(i))
c 
        diffs(i)=fsout(i)-fs2(i)
 3400 continue
c 
  
        call prin2('and fs2=*',fs2,n2*2)
        call prin2('and interpolation errors are=*',diffs,n2*2)
  
  
  
  
  
        return
  
        end
c 
c 
c 
c 
c 
        subroutine cfuntest(x,f)
        implicit real *8 (a-h,o-z)
        save
        complex *16 f,ima
c 
        data ima/(0.0d0,1.0d0)/
c 
        cc=0.9
c 
        f=sin(x*coef) + cos(x*coef) * ima + sin(x*coef*cc) * ima
cccc        f=sin(x*coef)
        return
c 
c 
c 
c 
        entry cfuntestin(coef7)
        coef=coef7
        return
        end
  
c 
c 
c 
c 
c 
        subroutine filttst(c1,xs1,ws1,n1,
     1      c2,xs2,ws2,n2,filter,v1,w1,rinterp,riprfast)
        implicit real *8 (a-h,o-z)
        save
        real *8 xs1(1),ws1(1),ws2(1),xs2(1),
     1      filter(n1,n2),fsin(1000),fsout(1000),rinterp(n2,n1),
     2      diffs(1000),fs2(1000),v1(1),t22(10 000),f22(10 000),
     3      f23(10 000),coefs(10 000),w1(1),coefslege(10000)
c 
        dimension ww(100000),fs1(1000),fs3(1000),errs(1000)
c 
c       construct the function to be interpolated
c 
        call funtestin(c2 * 0.9)
c 
        do 1200 i=1,n2
c 
        fsin(i)=xs2(i)
        fsin(i)=1
        call funtest(xs2(i),fsin(i))
 1200 continue
c 
cccc        call prolmat(filter,fsin,n1,n2,fsout)
  
        call prolc1c2filter(riprfast,fsin,fsout,ww)
  
  
        call prin2('and fsout=*',fsout,n1)
  
  
c 
        do 1400 i=1,n1
c 
        fs2(i)=xs1(i)
        fs2(i)=1
  
  
        call funtest(xs1(i),fs2(i))
c 
        diffs(i)=fs2(i)-fsout(i)
c 
 1400 continue
c 
        call prin2('and fs2=*',fs2,n1)
        call prin2(
     1   'differences between filtered and original functions=*',
     2      diffs,n1)
  
  
cccc        stop
c 
c       plot both of them things
c 
        iw=21
        call quagraph2(iw,xs1,fs2,n1,3,xs1,fsout,n1,3,
     1      'filterd vs. original*')
c 
c       evaluate inner products of the filtered function
c       with a bunch of complex exponentials, and compare
c       their inner products with those of the unfiltered
c       function
c 
        ntest=100
        h=c1/(ntest-1)
        do 1800 ii=1,ntest
c 
        cc=(ii-1)*h
  
cccc        call prin2('cc=*',cc,1)
  
        dd1=0
        do 1600 i=1,n1
c 
        d=sin(cc*xs1(i))
        dd1=dd1+ws1(i)*d*fsout(i)
 1600 continue
  
cccc        call prin2('and dd1=*',dd1,1)
  
  
        dd2=0
        do 1700 i=1,n2
c 
        d=sin(cc*xs2(i))
        dd2=dd2+ws2(i)*d*fsin(i)
 1700 continue
  
cccc        call prin2('and dd2=*',dd2,1)
cccc        call prin2('and dd2-dd1=*',dd2-dd1,1)
  
        errs(ii)=dd2-dd1
  
 1800 continue
  
  
        call prin2('errors in projection are=*',errs,ntest)
c 
c       construct a much finer discretization of the functiun
c       fsout, and plot it
c 
        call prolmat(v1,fsout,n1,n1,coefs)
  
        call prin2('and coefs=*',coefs,n1)
c 
        call protoleg(coefs,n1,w1,coefslege,nlege)
  
        call prin2('and coefslege=*',coefslege,nlege)
  
  
        ntest=1000
        done=1
        h=2*done/(ntest-1)
        do 2200 i=1,ntest
c 
        t22(i)=-done+(i-1)*h
c 
        call legeFDER(t22(i),f22(i),der,coefslege,Nlege-1)
c 
cccc        f23(i)=cos(t22(i)* c1)
        call funtest(t22(i),f23(i))
  
 2200 continue
c 
        iw=22
        call quagraph2(iw,t22,f22,ntest,3,t22,f23,ntest,3,
     1      'filterd vs. original after resampling*')
c 
        call prinf('and nlege=*',nlege,1)
  
  
  
cccc        stop
c 
c       test the interpolation operator
c 
  
        call funtestin(c1 * 0.9)
  
  
c 
        do 3200 i=1,n1
c 
        call funtest(xs1(i),fsin(i))
 3200 continue
c 
  
  
        call prolc1c2interp(riprfast,fsin,fsout,ww)
  
  
        call prin2('after prolproc1c2, fsout=*',fsout,n2)
  
  
c 
        do 3400 i=1,n2
c 
        call funtest(xs2(i),fs2(i))
c 
        diffs(i)=fsout(i)-fs2(i)
 3400 continue
c 
  
        call prin2('and fs2=*',fs2,n2)
        call prin2('and interpolation errors are=*',diffs,n2)
  
  
        stop
  
  
  
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
        f=sin(x*coef)
        return
c 
c 
c 
c 
        entry funtestin(coef7)
        coef=coef7
        return
        end
c 
c 
c 
c 
c 
        subroutine prolmat(a,x,n,m,y)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m)
        real *8 x(m),y(n),d
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
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c       This is the end of the debugging code, and the beginning
c       of the filtering/interpolation code proper
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c        This file contains 5 user-callable subroutines: prolfltinte,
c        prolc1c2interp, prolc1c2filter, cprolc1c2interp,
c        cprolc1c2filter; actually, prolc1c2filter, cprolc1c2filter
c        are entries rather than subroutines. Following are brief
c        descriptions of the fife user-callable entries.
c 
c   prolfltinte - produces filtering and interpolation matrices for
c        (non-periodic) band-limited functions on
c        the interval [-1,1]. It also generates the requisite
c        discretization nodes and the corresponding weights.
c        The band-limits c1, c2 between which the interpolation/
c        filtering are to be performed are more or less arbitrary.
c 
c   prolc1c2interp - applies the interepolation operator, converting
c        a real vector f1 of length n1 into a real vector f2
c        of length n2. The operator is in the array riprfast, where
c        it has been stored (one fervently hopes) by a prior call to
c        the subroutine prolfltinte (see)
c 
c   prolc1c2filter - applies the filtering operator, converting a
c        real vector f2 of length n2 into a real vector f1 of
c        length n1. The operator is in the array riprfast, where it
c        has been stored (one fervently hopes) by a prior call to
c        the subroutine prolfltinte (see)
c 
c   cprolc1c2interp - applies the interepolation operator, converting
c        a complex  vector f1 of length n1 into a complex vector f2
c        of length n2. The operator is in the array riprfast, where
c        it has been stored (one fervently hopes) by a prior call to
c        the subroutine prolfltinte (see)
c 
c   cprolc1c2filter - applies  the filtering operator, converting a
c        complex vector f2 of length n2 into a complex vector f1 of
c        length n1. The operator is in the array riprfast, where it
c        has been stored (one fervently hopes) by a prior call to
c        the subroutine prolfltinte (see)
c 
        subroutine prolfltinte(ier,wavesnum1,eps1,xs1,ws1,n1,
     1      wavesnum2,eps2,xs2,ws2,n2,ifrefine,rinterp,filter,
     2      w1,lenw1,riprfast)
        implicit real *8 (a-h,o-z)
        save
        real *8 w1(1),rinterp(1),filter(1)
        dimension xs1(1),ws1(1),xs2(1),ws2(1),riprfast(1)
c 
c       This subroutine produces filtering and interpolation
c       matrices for (non-periodic) band-limited functions on
c       the interval [-1,1]. It also generates the requisite
c       discretization nodes and the corresponding weights.
c       The band-limits c1, c2 between which the interpolation/
c       filtering are to be performed are more or less arbitrary.
c 
  
c 
c                      Input parameters:
c 
c  wavesnum1 - the lower band-limit (from which the functions are to be
c       interpolated, and to which they will be filtered), in
c       wavelengths
c  eps1 - the accuracy with which the functions with band-limit c1
c       are to be represented
c  wavesnum2 - the upper band-limit (to which the functions are to be
c       interpolated, and from which they will be filtered), in
c       wavelengths
c  eps2 - the accuracy with which the functions with band-limit c2
c       are to be represented
c 
c     IMPORTANT NOTE: C2 MUST BE GREATER THAN OR EQUAL TO C1  (!!!!)
c 
c  ifrefine - tells the subroutine whether the interpolation/quadrature
c       nodes are to be "refined":
c     ifrefine=1 will produce refined nodes
c     ifrefine=0 will suppress the refining
c  lenw1 - the length of the work array w1 in real *8 words
c 
c                      Output parameters:
c 
c  ier - error return code:
c     ier=0 means successful conclusion
c     ier=4 means that c2 < c1
c     ier>4 means that the memory available in the array w1 is
c       insufficient
c  xs1 - the discretization of the interval [-1,1] corresponding to
c       the band-limit c1
c  ws1 - the quadrature node scorresponding to the nodes xs1
c  n1 - the number of nodes in the array xs1
c  xs2 - the discretization of the interval [-1,1] corresponding to
c       the band-limit c2
c  ws2 - the quadrature node scorresponding to the nodes xs2
c  n2 - the number of nodes in the array xs2
c  rinterp - the matrix of interpolation from the nodes in array xs1
c        to the nodes in array xs2, dimensioned n2 \times n1
c  filter - the matrix of filtering from the nodes in array xs2
c        to the nodes in array xs1, dimensioned n1 \times n2
c  riprfast - the array containing symmetrized versions of the
c        operators rinterp, filter (se above), to be applied to
c        arbitrary vectors by the subroutine prolc1c2interp, and its
c        entry prolc1c2filter (see). Using these entries saves about
c        a factor of 2 in the CPU time over the direct application of
c        the operators rinterp, filter
c 
        ier=0
c 
c        construct the prolate functions and related machinery for
c        the first band-limit
c 
        done=1
        pi=atan(done)*4
        c1=wavesnum1*pi
        c2=wavesnum2*pi
  
  
  
cccc        call prin2('in prolfltinte, c1=*',c1,1)
cccc        call prin2('in prolfltinte, c2=*',c2,1)
  
  
cccc        stop
  
        nn1=c1+100
        iu1=1
        lu1=nn1**2
c 
        iv1=iu1+lu1
        lv1=nn1**2
c 
        iw1=iv1+lv1
        lw1=lenw1-iw1
c 
        call prolexps(jer1,c1,eps1,ifrefine,n1,xs1,ws1,
     1      w1(iu1),w1(iv1),w1(iw1),lw1,keep1,lused1)
c 
        call prinf('after first prolexps, jer1=*',jer1,1)
c 
c        construct the prolate functions and related machinery for
c        the second band-limit
c 
        nn2=c2+100
c 
        iu2=iw1+keep1+10
        lu2=nn2**2
c 
        iv2=iu2+lu2
        lv2=nn2**2
c 
        iw2=iv2+lv2
        lw2=lenw1-iw2
c 
        call prolexps(jer2,c2,eps2,ifrefine,n2,xs2,ws2,
     1      w1(iu2),w1(iv2),w1(iw2),lw2,keep2,lused2)
c 
        call prinf('after second prolexps, jer2=*',jer2,1)
c 
c 
        call prinf('after both calls to prolexps, n1=*',n1,1)
        call prinf('after both calls to prolexps, n2=*',n2,1)
c 
c       construct the interpolation scheme from the n1 nodes
c       to n2 nodes
c 
        ia=iw2+keep2+10
c 
        call prointc1c2(w1(iw1),xs2,n1,n2,w1(ia),w1(iv1),rinterp)
c 
c        construct the filtering scheme from n2 nodes and frequency
c        c2 to n1 nodes and frequency c1
c 
        nlegend=c2+100
c 
        ifuns1=iw2+keep2+10
        lfuns1=n1*nlegend+10
c 
        ifuns2=ifuns1+lfuns1
        lfuns2=n2*nlegend+10
c 
        iprods=ifuns2+lfuns2
        lprods=n1*n2+10
c 
        iwork=iprods+lprods
        lwork=n1*n2+10
c 
        call prodsc1c2(nlegend,w1(iw1),w1(iw2),w1(ifuns1),
     1      w1(ifuns2),n1,n2,w1(iprods) )
c 
        call promamac1c2(w1(iu1),w1(iprods),n1,n1,n2,w1(iwork) )
        call promamac1c2(w1(iwork),w1(iv2),n1,n2,n2,filter)
c 
c 
c       construct the symmetrized projection and interpolation
c       matrices, in order to save a factor of 2 in applying these
c       matrices to vectors
c 
c       . . . allocate memory
c 
  
        isproj=11
        lsproj=n1*n2/4+10
c 
        iaproj=isproj+lsproj
        laproj=n1*n2/4+10
c 
        isrint=iaproj+laproj
        lsrint=n1*n2/4+10
c 
        iarint=isrint+lsrint
        larint=n1*n2/4+10
c 
        riprfast(1)=isproj+0.1
        riprfast(2)=iaproj+0.1
        riprfast(3)=isrint+0.1
        riprfast(4)=iarint+0.1
c 
        riprfast(5)=n1+0.1
        riprfast(6)=n2+0.1
c 
c       . . . symmetrize the matrix filter
c 
        call prolsymc1c2(filter,riprfast(isproj),
     1      riprfast(iaproj),n1,n2,w1)
c 
c        now, symmetrize the matrix rinterp
c 
        call prolsymc1c2(rinterp,riprfast(isrint),
     1      riprfast(iarint),n2,n1,w1)
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine cprolc1c2interp(riprfast,fs1,fs2,w)
        implicit real *8 (a-h,o-z)
        save
        real *8 riprfast(1),fs1(1),fs2(1),w(1)
c 
c        This subroutine applies the interepolation operator,
c        converting a complex  vector f1 of length n1 into a complex
c        vector f2 of length n2. The operator is in the array
c        riprfast, where it has been stored (one fervently hopes)
c        by a prior call to the subroutine prolfltinte (see)
c 
c                    Input patrameters:
c 
c  riprfast - (hopefully) created by a prior call to prolfltinte (see)
c  fs1 - the vector to be interpolated (of length n1)
c 
c                    Output parameters:
c 
c  fs2 - the interpolated version of the vector n2
c 
c                    Work array:
c 
c  w - must be at least n1*n2+100 real *8 elements
c 
c 
c       . . . decode the first few elements of array riprfast
c 
        isrint=riprfast(3)
        iarint=riprfast(4)
c 
        n1=riprfast(5)
        n2=riprfast(6)
c 
        iwp=1
        lwp=n1+2
c 
        iwm=iwp+lwp
        lwm=n1+2
c 
        ifp=iwm+lwm
        lfp=n2+2
c 
        ifm=ifp+lfp
        lfm=n2+2
c 
c       interpolate the user-supplied array fs1 of length n1,
c       obtaining the array fs2 of length n2
c 
        call cprolproc1c2(riprfast(isrint),riprfast(iarint),
     1      fs1,fs2,n2,n1,w(iwp),w(iwm),w(ifp),w(ifm) )
        return
c 
c 
c 
c 
        entry cprolc1c2filter(riprfast,fs2,fs1,w)
c 
c        This entry applies the filtering operator,
c        converting a complex vector f2 of length n2 into a complex
c        vector f1 of length n1. The operator is in the array
c        riprfast, where it has been stored (one fervently hopes)
c        by a prior call to the subroutine prolfltinte (see)
c 
c                    Input patrameters:
c 
c  riprfast - (hopefully) created by a prior call to prolfltinte (see)
c  fs2 - the vector to be filtered (of length n2)
c 
c                    Output parameters:
c 
c  fs1 - the filtered version of the vector n1
c 
c                    Work array:
c 
c  w - must be at least n1*n2+100 real *8 elements
c 
c 
c       . . . decode the first few elements of array riprfast
c 
        isproj=riprfast(1)
        iaproj=riprfast(2)
c 
        n1=riprfast(5)
        n2=riprfast(6)
c 
        iwp=1
        lwp=n2+2
        lwp=lwp*2
c 
        iwm=iwp+lwp
        lwm=n2+2
        lwm=lwm*2
c 
        ifp=iwm+lwm
        lfp=n1+2
        lfp=lfp*2
c 
        ifm=ifp+lfp
        lfm=n1+2
        lfm=lfm*2
c 
c       filter the user-supplied array fs2 of length n2,
c       obtaining the array fs1 of length n1
c 
        call cprolproc1c2(riprfast(isproj),riprfast(iaproj),
     1      fs2,fs1,n1,n2,w(iwp),w(iwm),w(ifp),w(ifm) )
  
        return
        end
  
c 
c 
c 
c 
c 
        subroutine prolc1c2interp(riprfast,fs1,fs2,w)
        implicit real *8 (a-h,o-z)
        save
        real *8 riprfast(1),fs1(1),fs2(1),w(1)
c 
c        This subroutine applies the interepolation operator,
c        converting a real vector f1 of length n1 into a real
c        vector f2 of length n2. The operator is in the array
c        riprfast, where it has been stored (one fervently hopes)
c        by a prior call to the subroutine prolfltinte (see)
c 
c                    Input patrameters:
c 
c  riprfast - (hopefully) created by a prior call to prolfltinte (see)
c  fs1 - the vector to be interpolated (of length n1)
c 
c                    Output parameters:
c 
c  fs2 - the interpolated version of the vector n2
c 
c                    Work array:
c 
c  w - must be at least n1*n2+100 real *8 elements
c 
c 
c       . . . decode the first few elements of array riprfast
c 
        isrint=riprfast(3)
        iarint=riprfast(4)
c 
        n1=riprfast(5)
        n2=riprfast(6)
c 
        iwp=1
        lwp=n1+2
c 
        iwm=iwp+lwp
        lwm=n1+2
c 
        ifp=iwm+lwm
        lfp=n2+2
c 
        ifm=ifp+lfp
        lfm=n2+2
c 
c       interpolate the user-supplied array fs1 of length n1,
c       obtaining the array fs2 of length n2
c 
        call prolproc1c2(riprfast(isrint),riprfast(iarint),
     1      fs1,fs2,n2,n1,w(iwp),w(iwm),w(ifp),w(ifm) )
        return
c 
c 
c 
c 
        entry prolc1c2filter(riprfast,fs2,fs1,w)
c 
c        This entry applies the filtering operator,
c        converting a real vector f2 of length n2 into a real
c        vector f1 of length n1. The operator is in the array
c        riprfast, where it has been stored (one fervently hopes)
c        by a prior call to the subroutine prolfltinte (see)
c 
c                    Input patrameters:
c 
c  riprfast - (hopefully) created by a prior call to prolfltinte (see)
c  fs2 - the vector to be filtered (of length n2)
c 
c                    Output parameters:
c 
c  fs1 - the filtered version of the vector n1
c 
c                    Work array:
c 
c  w - must be at least n1*n2+100 real *8 elements
c 
c 
c       . . . decode the first few elements of array riprfast
c 
        isproj=riprfast(1)
        iaproj=riprfast(2)
c 
        n1=riprfast(5)
        n2=riprfast(6)
c 
        iwp=1
        lwp=n2+2
c 
        iwm=iwp+lwp
        lwm=n2+2
c 
        ifp=iwm+lwm
        lfp=n1+2
c 
        ifm=ifp+lfp
        lfm=n1+2
c 
c       filter the user-supplied array fs2 of length n2,
c       obtaining the array fs1 of length n1
c 
        call prolproc1c2(riprfast(isproj),riprfast(iaproj),
     1      fs2,fs1,n1,n2,w(iwp),w(iwm),w(ifp),w(ifm) )
  
        return
        end
c 
c 
c 
c 
c 
        subroutine prolproc1c2(sproj,aproj,fin,fout,n1,n2,
     1      wp,wm,fp,fm)
        implicit real *8 (a-h,o-z)
        save
        real *8 sproj(n1/2,n2/2),aproj(n1/2,n2/2),fin(1),fout(1),
     1           wp(1),wm(1),fp(1),fm(1)
c 
c       symmetrize and antisymmetrize the input vector fin
c 
        do 2400 i=1,n2/2
c 
        wp(i)=fin(i)+fin(n2-i+1)
        wm(i)=fin(i)-fin(n2-i+1)
 2400 continue
c 
c       apply the symmetrized and antisymmetrized projection
c       matrices to the symmetrized and antisymmetrized input
c       vectors, respectively
c 
        call prolmav2(sproj,wp,n1/2,n2/2,fp)
        call prolmav2(aproj,wm,n1/2,n2/2,fm)
c 
c       combine the symmetrized vectors appropriately
c 
        do 2600 i=1,n1/2
c 
        fout(i)=fp(i)+fm(i)
        fout(n1-i+1)=fp(i)-fm(i)
 2600 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine cprolproc1c2(sproj,aproj,fin,fout,n1,n2,
     1      wp,wm,fp,fm)
        implicit real *8 (a-h,o-z)
        save
        real *8 sproj(n1/2,n2/2),aproj(n1/2,n2/2)
c 
        complex *16 fin(1),fout(1),wp(1),wm(1),fp(1),fm(1)
c 
c       symmetrize and antisymmetrize the input vector fin
c 
        do 2400 i=1,n2/2
c 
        wp(i)=fin(i)+fin(n2-i+1)
        wm(i)=fin(i)-fin(n2-i+1)
 2400 continue
c 
c       apply the symmetrized and antisymmetrized projection
c       matrices to the symmetrized and antisymmetrized input
c       vectors, respectively
c 
        call cprolma2(sproj,wp,n1/2,n2/2,fp)
        call cprolma2(aproj,wm,n1/2,n2/2,fm)
c 
c       combine the symmetrized vectors appropriately
c 
        do 2600 i=1,n1/2
c 
        fout(i)=fp(i)+fm(i)
        fout(n1-i+1)=fp(i)-fm(i)
 2600 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine prolsymc1c2(proj,sproj,aproj,n1,n2,ww)
        implicit real *8 (a-h,o-z)
        save
        real *8 proj(n1,n2),sproj(n1/2,n2/2),aproj(n1/2,n2/2),
     1      ww(n1,n2)
c 
c       symmetrize the matrix proj
c 
        do 2400 i=1,n1
        do 2200 j=1,n2/2
c 
        ww(i,j)=proj(i,j)+proj(i,n2-j+1)
 2200 continue
 2400 continue
c 
        do 2800 i=1,n1/2
        do 2600 j=1,n2/2
c 
        sproj(i,j)=(ww(i,j)+ww(n1-i+1,j))/4
 2600 continue
 2800 continue
c 
c       antisymmetrize the matrix proj
c 
        do 3400 i=1,n1
        do 3200 j=1,n2/2
c 
        ww(i,j)=proj(i,j)-proj(i,n2-j+1)
 3200 continue
 3400 continue
c 
        do 3800 i=1,n1/2
        do 3600 j=1,n2/2
c 
        aproj(i,j)=(ww(i,j)-ww(n1-i+1,j))/4
 3600 continue
 3800 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine prodsc1c2(n,w1,w2,funs1,funs2,n1,n2,prods)
        implicit real *8 (a-h,o-z)
        save
        real *8 w1(1),w2(1),funs1(n,1),funs2(n,1),prods(n1,n2)
c 
        real *8 ts(10 000),whts(10 000)
c 
c       discretize the intervals [-1,1], [0,1]
c 
        itype=1
        call legeexps(itype,n,ts,u,v,whts)
c 
c       evaluate both sets of prolate functions at the quadrature nodes
c 
        do 1600 i=1,n1
c 
        do 1400 j=1,n
c 
        call prolevv(i-1,ts(j),w1,val)
        funs1(j,i)=val
 1400 continue
 1600 continue
c 
c 
        do 2000 i=1,n2
c 
        do 1800 j=1,n
c 
        call prolevv(i-1,ts(j),w2,val)
        funs2(j,i)=val
 1800 continue
 2000 continue
c 
c       construct the matrix of inner products between thge
c       obtained functions
c 
        do 2400 i=1,n1
        do 2200 j=1,n2
c 
        call prolsca(funs1(1,i),funs2(1,j),n,prod,whts)
c 
        prods(i,j)=prod
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
        subroutine prointc1c2(w1,xs2,n1,n2,a,v1,rinterp)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n2,n1),xs2(1),w1(1),v1(n1,n1),rinterp(n2,n1)
c 
c        construct the matrix of values of prolate functions
c        with c=c1 at the prolate nodes on the interval [-1,1]
c        corresponding to c=c2
c 
        do 1600 i=1,n1
        do 1400 j=1,n2
c 
        call prolevv(i-1,xs2(j),w1,val)
c 
        a(j,i)=val
 1400 continue
 1600 continue
c 
c       multiply the matrix a by the decomposition matrix u
c       from the right, obtaining the interpolation matrix
c 
        call promamac1c2(a,v1,n2,n1,n1,rinterp)
c 
        return
        end
c 
c 
c 
c 
c 
       subroutine promamac1c2(a,b,n,m,k,c)
       implicit real *8 (a-h,o-z)
        save
       dimension a(n,m),b(m,k),c(n,k)
c 
       do 1600 i=1,n
       do 1400 j=1,k
       d=0
       do 1200 jj=1,m
c 
       d=d+a(i,jj)*b(jj,j)
 1200 continue
c 
       c(i,j)=d
 1400 continue
 1600 continue
c 
       return
       end
