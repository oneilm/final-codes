        implicit real *8 (a-h,o-z)
        real *8 tt(10 000),ww(10 000),www(100 000)
c 
        complex *16 ima,
     1      ff(10 000),ff2(10 000),
     2      diffs(10 000),derout,ff1(10 000),dersout(10 000),
     3      pattexp(10 000)
c 
        dimension tt1(10 000),
     1      ff1rea(10 000),ff1ima(10 000),ww1(10 000)
c 
c 
c 
c 
        data ima/(0.0d0,1.0d0)/
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER supergain'
        READ *,sgain
        CALL PRIN2('sgain=*',sgain,1 )
C 
  
  
c 
c       retrieve the nodes and weights of the prolate quadrature
c 
        call expo_ttww_ret63(tt,ww,nn,c)
cccc        call expo_ttww_ret20(tt,ww,nn,c)
c 
c       construct the function to be expanded
c 
        do 1200 i=1,nn
c 
cccc        ff(i)=1 +ima*tt(i)**5
        ff(i)=1/sqrt(3+tt(i)) + ima*tt(i)**2
        ff(i)=1/sqrt(3+tt(i))
cccc        ff(i)=1
cccc        ff(i)=1/sqrt(1.2+tt(i))
  
 1200 continue
  
        call prin2('ff as created*',ff,nn*2)
c 
c       expand the function into exponentials
c 
        nexpons=nn
c 
        time1=clotatim()
  
  
        ifinit=1
        do 1400 i=1,1000
cccc        do 1400 i=1,1
        call expoapp63(ier,sgain,ff,sgainout,www,pattexp,ifinit)
cccc        call expoapp20(ier,sgain,ff,sgainout,www,pattexp,ifinit)
c 
        ifinit=0
 1400 continue
  
  
        time2=clotatim()
  
        call prin2('after expoapp63, time2-time1=*',time2-time1,1)
        call prinf('after expoapp63, ier=*',ier,1)
        call prin2('actual supergain of obtained pattern*',sgainout,1)
c 
c       evaluate the obtained exponential expansion at the nodes of
c       the discretization and compare it to the original function ff
c 
        diffmax=0
        do 5200 i=1,nn
c 
        call expoextr(c,tt,pattexp,nexpons,tt(i),ff2(i),derout )
c 
        diffs(i)=ff2(i)-ff(i)
  
        d=abs(diffs(i))
  
        if(d .gt. diffmax) diffmax=d
 5200 continue
  
        call prin2('errors of exponential expansion are*',
     1       diffs,nn*2)
        call prin2('maximum error of exponential expansion is*',
     1       diffmax,1)
c 
c       discretize the interval [-2,2] and plot the real and imaginary parts
c 
        nplot=500
        h=1
        h=2*h/(nplot-1)
c 
        do 5400 i=1,nplot
c 
        tt1(i)=-1+(i-1)*h
        ww1(i)=h
c 
        tt1(i)=tt1(i)*3
        ww1(i)=ww1(i)*3
c 
        call expoextr(c,tt,pattexp,nexpons,tt1(i),ff1(i),dersout(i))
c 
        ff1rea(i)=ff1(i)
        ff1ima(i)=-ff1(i)*ima
 5400 continue
c 
c 
        iw=25
        call quagraph(iw,tt1,ff1rea,nplot,3,
     1      'extended real part*')
c 
c 
        iw=26
        call quagraph(iw,tt1,ff1ima,nplot,3,
     1      'extended imaginary part*')
  
  
  
        call prinf('and nexpons=*',nexpons,1)
        call prin2('and pattexp=*',pattexp,nexpons*2)
  
        stop
        end
c 
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c         This is the end of the debugging code and the beginning of the
c         actual band-limited approximatipn routines
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c      This file contains 5 user-callable subroutines and entries:
c      expoapp20, expo_ttww_ret20, expoapp63, expo_ttww_ret63,
c      expoextr. Following is a brief description of these entries.
c 
c   expoapp20 - constructs an expansion of the user-supplied function
c      into a linear combination of 38 exponentials, of the form
c 
c          f(x)=\sum_{j=1}^{38} coesexp_j * exp(ima*c*tt_j),        (1)
c 
c      with c=20.
c 
c   expo_ttww_ret20 - returns to the user the 38 prolate nodes tt
c      and their corresponding weights ww, discretizing to 16-digit
c      (or so) accuracy band-limited functions with c=20 on the
c      interval [-1,1].
c 
c   expoapp63 - constructs an expansion of the user-supplied
c      function into a linear combination of 74
c      exponentials, of the form
c 
c          f(x)=\sum_{j=1}^{74} coesexp_j * exp(ima*c*tt_j),        (2)
c 
c        with c=63.
c 
c   expo_ttww_ret63 - returns to the user the 74 prolate nodes tt
c      and their corresponding weights ww, discretizing to 16-digit
c      (or so) accuracy band-limited functions with c=63 on the
c      interval [-1,1].
c 
c   expoextr - given a collection of user-provided frequencies ttt
c       and (complex) coefficients vals(j),this subroutine evaluates
c       the expansion of the form
c 
c       fout = sum_{j=1}^nnn
c              vals(i) \cdot e^{i \cdot c \cdot ttt(j) \cdot x}      (3)
c 
c       at the point x \in R^2; the subroutine also returns the
c       derivative of fout with respect to x
c 
c 
        subroutine expoapp20(ier,sgain,ff,
     1      sgainout,www,coesexp,ifinit)
        implicit real *8 (a-h,o-z)
        save
        dimension tt(50),cm(50),ww(50),www(38,1),rcoefs(100)
c 
        complex *16 rlams(50),ff(1),cd,alphas(50),
     1      coesexp(1)
c 
c        This subroutine constructs an expansion of the user
c        -supplied function into a linear combination of 38
c        exponentials, of the form
c 
c             f(x)=\sum_{j=1}^{nn} coesexp_j * exp(ima*c*tt_j),        (1)
c 
c        with c=20. The function is specified by the user in the
c        form of a table of values tabulated at the nodes tt; the
c        nodes tt have been (hopefully) retrieved by the user from
c        the entry expo_ttww_ret20 of the subroutine exporet20 (see).
c 
c                      Input parameters:
c 
c  sgain - the supergain the approximation is permitted to have
c  ff - the values at the nodes tt (obtained from the entry
c       expo_ttww_ret20 of the subroutine exporet20) of the
c       function to be expanded
c  ifinit - the integer parameter telling the subroutine whether
c       the array www has to be created (it is a minor CPU time
c       saving feature).
c     ifinit =1 means that www has to be initialized
c     ifinit =0 means that the initialization is to be skipped
c 
c 
c                      Output parameters:
c 
c  ier - error return code
c  sgainout - the actual supergain of the obtained approximation
c  coesexp - the (complex) coefficients of the expansion (38 of
c       them things)
c 
c                      Work array:
c 
c  www - must be at least 1600 real *8 locations long
c 
        nn=38
        n=nn+1
        c=20
c 
c       find the L^2 - norm of the user supplied right-hand
c       side, in order to evaluate the total (visible+invisible)
c       energy of the approximation to be constructed
c 
        if(ifinit .eq. 1) call exporet20(tt,ww,cm,rlams,www)
c 
        ier=0
        d=0
c 
        do 1600 i=1,nn
c 
        d=d+ff(i)*conjg(ff(i))*ww(i)
 1600 continue
c 
        totpowr=d*sgain
c 
c       construct the projections of the function ff on all
c       relevant prolate functions
c 
        call expoproj(www,ff,ww,alphas,nn,n)
c 
        do 2200 i=1,n
c 
        rcoefs(i)= cm(i)*conjg(alphas(i))*alphas(i)
 2200 continue
c 
c       find the Lagrange multiplier corresponding to the
c       user-supplied total power
c 
        call expobis(jer,cm,alphas,n,totpowr,rlout,rcoefs)
  
cccc        call prinf('after expobis, jer=*',jer,1)
c 
        if(jer .ne. 0) ier=2
c 
c       evaluate all prolate coefficients
c 
        sumtes=0
        denom=0
        do 2400 i=1,n
c 
        alphas(i)=alphas(i)/(1+rlout*cm(i))
  
        if(cm(i) .lt. 0.1) alphas(i)=0
c 
        sumtes=sumtes+cm(i)*(conjg(alphas(i))*alphas(i))
c 
        denom=denom+(conjg(alphas(i))*alphas(i))
 2400 continue
c 
        sgainout=sumtes/denom
c 
c       determine the number of alphas that are not zero,
c       and scale these by 1/rlams(i)
c 
        do 2600 i=n,1,-1
c 
        nbetas=i
        if(abs(alphas(i)) .gt. 1.0d-60) goto 2800
 2600 continue
 2800 continue
c 
        n=nbetas-1
c 
        do 3200 i=1,nn
c 
        alphas(i)=alphas(i)/rlams(i)
 3200 continue
c 
c        construct the exponential expansion
c 
        do 4400 i=1,nn
c 
        cd=0
        do 4200 j=1,n-1
        cd=cd+alphas(j)*(www(i,j)*ww(i))
 4200 continue
        coesexp(i)=cd
 4400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine expoapp63(ier,sgain,ff,
     1      sgainout,www,coesexp,ifinit)
        implicit real *8 (a-h,o-z)
        save
        dimension tt(100),cm(100),ww(100),www(74,1),rcoefs(200)
c 
        complex *16 rlams(100),ff(1),cd,alphas(100),
     1      coesexp(1)
c 
c        This subroutine constructs an expansion of the user
c        -supplied function into a linear combination of 74
c        exponentials, of the form
c 
c             f(x)=\sum_{j=1}^{nn} coesexp_j * exp(ima*c*tt_j),        (1)
c 
c        with c=63. The function is specified by the user in the
c        form of a table of values tabulated at the nodes tt; the
c        nodes tt have been (hopefully) retrieved by the user from
c        the entry expo_ttww_ret63 of the subroutine exporet63 (see).
c 
c                      Input parameters:
c 
c  sgain - the supergain the approximation is permitted to have
c  ff - the values at the nodes tt (obtained from the entry
c       expo_ttww_ret63 of the subroutine exporet63) of the
c       function to be expanded
c  ifinit - the integer parameter telling the subroutine whether
c       the array www has to be created (it is a minor CPU time
c       saving feature).
c     ifinit =1 means that www has to be initialized
c     ifinit =0 means that the initialization is to be skipped
c 
c 
c                      Output parameters:
c 
c  ier - error return code
c  sgainout - the actual supergain of the obtained approximation
c  coesexp - the (complex) coefficients of the expansion (74 of
c       them things)
c 
c                      Work array:
c 
c  www - must be at least 1600 real *8 locations long
c 
        nn=74
        n=nn+1
        c=63
c 
c       find the L^2 - norm of the user supplied right-hand
c       side, in order to evaluate the total (visible+invisible)
c       energy of the approximation to be constructed
c 
        if(ifinit .eq. 1) call exporet63(tt,ww,cm,rlams,www)
c 
        ier=0
        d=0
c 
        do 1600 i=1,nn
c 
        d=d+ff(i)*conjg(ff(i))*ww(i)
 1600 continue
c 
        totpowr=d*sgain
c 
c       construct the projections of the function ff on all
c       relevant prolate functions
c 
        call expoproj(www,ff,ww,alphas,nn,n)
c 
cccc        call prin2('alphas via expoproj=*',alphas,n*2)
c 
        do 2200 i=1,n
c 
        rcoefs(i)= cm(i)*conjg(alphas(i))*alphas(i)
 2200 continue
c 
c       find the Lagrange multiplier corresponding to the
c       user-supplied total power
c 
        call expobis(jer,cm,alphas,n,totpowr,rlout,rcoefs)
  
cccc        call prinf('after expobis, jer=*',jer,1)
c 
        if(jer .ne. 0) ier=2
c 
c       evaluate all prolate coefficients
c 
        sumtes=0
        denom=0
        do 2400 i=1,n
c 
        alphas(i)=alphas(i)/(1+rlout*cm(i))
  
        if(cm(i) .lt. 0.1) alphas(i)=0
c 
        sumtes=sumtes+cm(i)*(conjg(alphas(i))*alphas(i))
c 
        denom=denom+(conjg(alphas(i))*alphas(i))
 2400 continue
c 
        sgainout=sumtes/denom
c 
c       determine the number of alphas that are not zero,
c       and scale these by 1/rlams(i)
c 
        do 2600 i=n,1,-1
c 
        nbetas=i
        if(abs(alphas(i)) .gt. 1.0d-60) goto 2800
 2600 continue
 2800 continue
c 
        n=nbetas-1
c 
        do 3200 i=1,nn
c 
        alphas(i)=alphas(i)/rlams(i)
 3200 continue
c 
c        construct the exponential expansion
c 
        do 4400 i=1,nn
c 
        cd=0
        do 4200 j=1,n-1
        cd=cd+alphas(j)*(www(i,j)*ww(i))
 4200 continue
        coesexp(i)=cd
 4400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine expoproj(www,ff,ww,alphas,nn,n)
        implicit real *8 (a-h,o-z)
        save
        dimension ww(1),www(nn,1)
        complex *16 ff(1),cd,alphas(1),symm(100),antisymm(100)
c 
c        symmetrize and antisymmetrize the function ff
c 
        do 1200 i=1,nn/2
c 
        symm(i)=ff(i)+ff(nn-i+1)
        antisymm(i)=ff(i)-ff(nn-i+1)
c 
 1200 continue
c 
cccc        call prin2('symm=*',symm,nn)
cccc        call prin2('antisymmy=*',antisymm,nn)
c 
c       construct the projections of the function ff on all
c       relevant prolate functions
c 
        do 1800 i=1,n/2
c 
        i1=i*2-1
        i2=i1+1
c 
        call expoinsca(symm,www(1,i1),nn/2,cd,ww)
c 
        alphas(i1)=cd
c 
        call expoinsca(antisymm,www(1,i2),nn/2,cd,ww)
c 
        alphas(i2)=cd
 1800 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine expobis(ier,cm,prods,n,totpowr,rlout,rcoefs)
        implicit real *8 (a-h,o-z)
        save
        dimension cm(1),prods(1),rcoefs(1)
c 
c       initialize the bisection process
c 
        ier=0
        rl1=0
        rl2=1.0d10
        rlout=0
c 
        call expoinsum2(cm,n,rl1,f1,rcoefs)
c 
        if( (f1 .gt. totpowr) .and. (f2 .lt. totpowr) ) goto 1200
c 
        if(f1 .lt. totpowr) ier=4
        if(f2 .gt. totpowr) ier=16
        return
c 
 1200 continue
c 
c        reduce rl2 by decimal scaling as far as it will go
c 
  
        do 1300 i=1,100
c 
cccc        call prinf('i=*',i,1)
  
        rl22=rl2/1000
        call expoinsum2(cm,n,rl22,f22,rcoefs)
c 
        if(f22 .ge. totpowr) goto 1350
c 
        rl2=rl22
        f2=f22
 1300 continue
 1350 continue
c 
        do 2000 i=1,15
c 
        rl3=(rl1+rl2)/2
c 
        call expoinsum2(cm,n,rl3,f3,rcoefs)
c 
        if(f3 .lt. totpowr) goto 1400
c 
        f1=f3
        rl1=rl3
        goto 2000
 1400 continue
c 
        f2=f3
        rl2=rl3
c 
 2000 continue
c 
        rlout=(rl1+rl2)/2
c 
c       use Newton to tighten the obtained value of rlout
c 
        do 1600 i=1,5
c 
        call expoinsum2(cm,n,rlout,sum,rcoefs)
        call expoinder2(cm,n,rlout,der,rcoefs)
c 
        ff=sum-totpowr
c 
        rlout=rlout-ff/der
 1600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine expoinder2(cm,n,rlagr,der,rcoefs)
        implicit real *8 (a-h,o-z)
        save
        dimension cm(1),rcoefs(1)
c 
c       evaluate the sum that is supposed to be equal to the
c       total power
c 
        der=0
        do 1200 i=1,n
c 
        der=der + cm(i)*rcoefs(i)/(1+rlagr*cm(i))**3
 1200 continue
c 
        der=-der*2
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine expoinsum2(cm,n,rlagr,sum,rcoefs)
        implicit real *8 (a-h,o-z)
        save
        dimension cm(1),rcoefs(1)
c 
c       evaluate the sum that is supposed to be equal to the
c       total power
c 
        done=1
        sum=0
        do 1200 i=1,n
c 
        sum=sum + rcoefs(i)/(1+rlagr*cm(i))**2
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine expoinsca(cx,y,n,prod,ww)
        implicit real *8 (a-h,o-z)
        save
        dimension y(1),ww(1)
        complex *16 cx(1),prod
c 
        prod=0
        do 1800 i=1,n
        prod=prod+cx(i)*(y(i)*ww(i))
 1800 continue
  
        return
        end
c 
c 
c 
c 
c 
        subroutine expoextr(c,ttt,vals,nnn,x,fout,derout)
        implicit real *8 (a-h,o-z)
        save
        dimension ttt(1)
        complex *16 vals(1),cd,ima,fout,derout,xcima
c 
        data ima/(0.0d0,1.0d0)/
c 
c       Given a collection of user-provided frequencies ttt and
c       (complex) coefficients vals(j),this subroutine evaluates
c       the expansion of the form
c 
c       fout = sum_{j=1}^nnn
c              vals(i) \cdot e^{i \cdot c \cdot ttt(j) \cdot x}      (1)
c 
c       at the point x \in R^2; the subroutine also returns the
c       derivative derout of fout with respect to x.
c 
c                       Input parameters:
c 
c  c - the parameter in (1) above
c  ttt - the exponents
c  vals - the intensities
c  nnn - the number of terms in (1)
c  x - the point where the expansion (1) is to be evaluated
c 
c                       Output parameters:
c 
c  fout - the value of the expansion (1) at the point x
c  derout  - the derivative of fout with respect to x
c 
c       . . . evaluate the function and the derivative at the user-specified
c             point x
c 
        xcima=ima*c*x
        fout=0
        derout=0
c 
        do 1200 i=1,nnn
c 
        cd=exp(xcima*ttt(i))*vals(i)
c 
        fout=fout+cd
        derout=derout+cd*ttt(i)
 1200 continue
c 
        derout=derout*ima*c
        return
        end
c 
c 
c 
c 
c 
        subroutine exporet63(tt,ww,cm,rlams,www)
        implicit real *8 (a-h,o-z)
        save
        dimension tt0(37),ww0(37),cm0(74),rlams0(78),
     1      tt(1),ww(1),cm(1),www1(180),www2(180),www3(180),
     2      www4(180),www5(180),www6(180),www7(180),www8(180),
     3      www9(180),www10(180),www11(180),www12(180),
     4      www13(180),www14(180),www15(180),www16(112),
     5      www(74,76),www0(3000),www101(37,2)
c 
        equivalence (www0(1),www1(1)),(www0(181),www2(1)),
     1       (www0(361),www3(1)),(www0(541),www4(1)),
     2       (www0(721),www5(1)),(www0(901),www6(1)),
     3       (www0(1081),www7(1)),(www0(1261),www8(1)),
     4       (www0(1441),www9(1)),(www0(1621),www10(1)),
     5       (www0(1801),www11(1)),(www0(1981),www12(1)),
     6       (www0(2161),www13(1)),(www0(2341),www14(1)),
     7       (www0(2521),www15(1)),(www0(2701),www16(1))
c 
        equivalence(www0(1),www101(1,1))
  
        complex *16 rlams(1),ima
c 
        data nn/74/,ima/(0.0d0,1.0d0)/
c 
        data tt0/
     5  -.9991622975568036D+00,-.9956022276930717D+00,
     6  -.9892612773641241D+00,-.9802416939782934D+00,
     7  -.9686813577261812D+00,-.9547425438011088D+00,
     8  -.9386016909497066D+00,-.9204404225409142D+00,
     9  -.9004383983572652D+00,-.8787681572645679D+00,
     *  -.8555917894531944D+00,-.8310591093806521D+00,
     1  -.8053069555497094D+00,-.7784592725411463D+00,
     2  -.7506276934219378D+00,-.7219124101547703D+00,
     3  -.6924031820647918D+00,-.6621803825363628D+00,
     4  -.6313160215002649D+00,-.5998747076843717D+00,
     5  -.5679145324284509D+00,-.5354878683840564D+00,
     6  -.5026420835124057D+00,-.4694201748877228D+00,
     7  -.4358613289455146D+00,-.4020014157056263D+00,
     8  -.3678734246365461D+00,-.3335078495321746D+00,
     9  -.2989330292523560D+00,-.2641754505637434D+00,
     *  -.2292600186868278D+00,-.1942103005538710D+00,
     1  -.1590487452357610D+00,-.1237968855152411D+00,
     2  -.8847552417364960D-01,-.5310490821808625D-01,
     3  -.1770489400349689D-01/
c 
        data ww0/
     4  0.2147971382760330D-02,0.4963923329520173D-02,
     5  0.7700680075925467D-02,0.1031532691768686D-01,
     6  0.1277812702795207D-01,0.1506994328373336D-01,
     7  0.1718137239588103D-01,0.1911113561881574D-01,
     8  0.2086410164066793D-01,0.2244933598090717D-01,
     9  0.2387841972778056D-01,0.2516413604028332D-01,
     *  0.2631952440050601D-01,0.2735725096407706D-01,
     1  0.2828922565416356D-01,0.2912639858186779D-01,
     2  0.2987867915786881D-01,0.3055493432223061D-01,
     3  0.3116303437307834D-01,0.3170992467826816D-01,
     4  0.3220170894679068D-01,0.3264373502890007D-01,
     5  0.3304067785461967D-01,0.3339661654079579D-01,
     6  0.3371510425375067D-01,0.3399923037954402D-01,
     7  0.3425167512378875D-01,0.3447475697566057D-01,
     8  0.3467047362009134D-01,0.3484053693076344D-01,
     9  0.3498640266536309D-01,0.3510929543982443D-01,
     *  0.3521022949639451D-01,0.3529002571119960D-01,
     1  0.3534932521631922D-01,0.3538859994221111D-01,
     2  0.3540816032000867D-01/
c 
        data cm0/
     3  0.1000000000000000D+01,0.1000000000000000D+01,
     4  0.1000000000000000D+01,0.1000000000000000D+01,
     5  0.1000000000000000D+01,0.1000000000000000D+01,
     6  0.1000000000000000D+01,0.1000000000000000D+01,
     7  0.1000000000000000D+01,0.1000000000000000D+01,
     8  0.1000000000000000D+01,0.1000000000000000D+01,
     9  0.1000000000000000D+01,0.1000000000000000D+01,
     *  0.1000000000000000D+01,0.1000000000000000D+01,
     1  0.1000000000000000D+01,0.1000000000000000D+01,
     2  0.1000000000000000D+01,0.1000000000000000D+01,
     3  0.1000000000000000D+01,0.1000000000000000D+01,
     4  0.1000000000000000D+01,0.1000000000000001D+01,
     5  0.1000000000000015D+01,0.1000000000000201D+01,
     6  0.1000000000002487D+01,0.1000000000028861D+01,
     7  0.1000000000313518D+01,0.1000000003186544D+01,
     8  0.1000000030269325D+01,0.1000000268292645D+01,
     9  0.1000002213851672D+01,0.1000016952706787D+01,
     *  0.1000119926128933D+01,0.1000778489966977D+01,
     1  0.1004588453312943D+01,0.1024148876293116D+01,
     2  0.1111176944490417D+01,0.1447684887610384D+01,
     3  0.2681465229046134D+01,0.7605909019388684D+01,
     4  0.3036684495533352D+02,0.1500498322706982D+03,
     5  0.8470753247604217D+03,0.5270020477690781D+04,
     6  0.3554112439779423D+05,0.2574356954798891D+06,
     7  0.1990311395590338D+07,0.1634749588070921D+08,
     8  0.1421125711649668D+09,0.1303524104849669D+10,
     9  0.1258260733258334D+11,0.1275273170612599D+12,
     *  0.1354425805666507D+13,0.1504745819926326D+14,
     1  0.1745997412369720D+15,0.2112903070773511D+16,
     2  0.2663237496276053D+17,0.3492374686544143D+18,
     3  0.4759259750603433D+19,0.6733300316157645D+20,
     4  0.9880519757858019D+21,0.1502502576291876D+23,
     5  0.2365803990919435D+24,0.3854218411418182D+25,
     6  0.6491918565714795D+26,0.1129773685769976D+28,
     7  0.2030063555993667D+29,0.3764089815313410D+30,
     8  0.7197602313308304D+31,0.1418569108978058D+33,
     9  0.2880156668739321D+34,0.6020921140143173D+35/
c 
        data rlams0/
     3  0.3158054782836449D+00,0.3158054782836449D+00,
     4  0.3158054782836449D+00,0.3158054782836449D+00,
     5  0.3158054782836449D+00,0.3158054782836449D+00,
     6  0.3158054782836449D+00,0.3158054782836449D+00,
     7  0.3158054782836449D+00,0.3158054782836449D+00,
     8  0.3158054782836449D+00,0.3158054782836449D+00,
     9  0.3158054782836449D+00,0.3158054782836449D+00,
     *  0.3158054782836449D+00,0.3158054782836449D+00,
     1  0.3158054782836449D+00,0.3158054782836449D+00,
     2  0.3158054782836449D+00,0.3158054782836449D+00,
     3  0.3158054782836449D+00,0.3158054782836449D+00,
     4  0.3158054782836449D+00,0.3158054782836447D+00,
     5  0.3158054782836425D+00,0.3158054782836132D+00,
     6  0.3158054782832521D+00,0.3158054782790877D+00,
     7  0.3158054782341395D+00,0.3158054777804808D+00,
     8  0.3158054735040356D+00,0.3158054359195098D+00,
     9  0.3158051287109822D+00,0.3158028014388421D+00,
     *  0.3157865433224720D+00,0.3156826243112518D+00,
     1  0.3150834327931803D+00,0.3120600063818407D+00,
     2  0.2995905074370301D+00,0.2624717936187662D+00,
     3  0.1928561863700381D+00,0.1145101797457262D+00,
     4  0.5730860320261958D-01,0.2578112724092602D-01,
     5  0.1085072217871872D-01,0.4350242911963569D-02,
     6  0.1675151450874822D-02,0.6224224903979428D-03,
     7  0.2238510545505818D-03,0.7810773446274666D-04,
     8  0.2649130854817351D-04,0.8747020133427219D-05,
     9  0.2815362603624995D-05,0.8843375131067677D-06,
     *  0.2713576375379882D-06,0.8141193677888181D-07,
     1  0.2389999783389688D-07,0.6870365980658344D-08,
     2  0.1935150342236489D-08,0.5343911332440496D-09,
     3  0.1447604602666306D-09,0.3848625622534704D-10,
     4  0.1004684643070356D-10,0.2576392621882187D-11,
     5  0.6492773897864849D-12,0.1608612720283756D-12,
     6  0.3919522192121082D-13,0.9395587293618134D-14,
     7  0.2216485232246792D-14,0.5147420770431821D-15,
     8  0.1177133544657230D-15,0.2651516965646360D-16,
     9  0.5884527588932938D-17,0.1287028576768551D-17,
     *  0.2774791726781680D-18,0.5898480025109839D-19,
     1  0.1236560169495527D-19,0.2557113317194904D-20/
  
  
c 
        data www1/
     1  0.8880073727918875D-25,0.1523622800184043D-23,
     2  0.3220121603046284D-22,0.6979986624424650D-21,
     3  0.1466527036037469D-19,0.2898433716812078D-18,
     4  0.5286592615410806D-17,0.8784857740344886D-16,
     5  0.1318526236450321D-14,0.1777602405895916D-13,
     6  0.2145743690982233D-12,0.2315675399127602D-11,
     7  0.2233822143656476D-10,0.1927547370359204D-09,
     8  0.1489845452021692D-08,0.1033328702393076D-07,
     9  0.6444594869151750D-07,0.3622232623132090D-06,
     *  0.1838953497441633D-05,0.8452150657357162D-05,
     1  0.3524782554508543D-04,0.1336585213782661D-03,
     2  0.4617907370822551D-03,0.1456502110910582D-02,
     3  0.4201182752429197D-02,0.1110058089666580D-01,
     4  0.2690875105011746D-01,0.5992607675863471D-01,
     5  0.1227593365008930D+00,0.2315759203673169D+00,
     6  0.4026790878803939D+00,0.6459844299112727D+00,
     7  0.9567482603913051D+00,0.1309022643776508D+01,
     8  0.1655306242952144D+01,0.1935292080918204D+01,
     9  0.2092448855972035D+01,-.1911820883388790D-23,
     *  -.3106716759015617D-22,-.6222119476038028D-21,
     1  -.1278307606451691D-19,-.2545921090645658D-18,
     2  -.4770248267329179D-17,-.8249089199831642D-16,
     3  -.1299628518876162D-14,-.1849237339214288D-13,
     4  -.2363078138399505D-12,-.2702892880970082D-11,
     5  -.2762771807336767D-10,-.2522755656240026D-09,
     6  -.2059022414785314D-08,-.1503892001720010D-07,
     7  -.9845434250003429D-07,-.5787838214998108D-06,
     8  -.3061366172466910D-05,-.1459815812325327D-04,
     9  -.6287996378076373D-04,-.2451106389538300D-03,
     *  -.8661415452957617D-03,-.2778774847075097D-02,
     1  -.8104379966408058D-02,-.2150985369678320D-01,
     2  -.5198983776043670D-01,-.1144741894200438D+00,
     3  -.2295822130437986D+00,-.4190488912958173D+00,
     4  -.6949468871483042D+00,-.1043960039427650D+01,
     5  -.1413238913335048D+01,-.1708685360990682D+01,
     6  -.1815043878307117D+01,-.1637220392687096D+01,
     7  -.1147462960114806D+01,-.4133634478524352D+00,
     8  0.2892181710846558D-22,0.4449002190728648D-21,
     9  0.8439943833861831D-20,0.1642658180999865D-18,
     *  0.3099703791230851D-17,0.5503246139237322D-16,
     1  0.9017876336636541D-15,0.1346246104865332D-13,
     2  0.1814894480371016D-12,0.2196788547131956D-11,
     3  0.2379204782233568D-10,0.2301536896635340D-09,
     4  0.1987573656447505D-08,0.1532884277117282D-07,
     5  0.1056820690421640D-06,0.6522177616089640D-06,
     6  0.3608878033653985D-05,0.1793359441206821D-04,
     7  0.8016760940830387D-04,0.3228829795606753D-03,
     8  0.1173303369112424D-02,0.3851197583433327D-02,
     9  0.1142798579228617D-01,0.3067123145867186D-01,
     *  0.7444992958331015D-01,0.1633321054535487D+00,
     1  0.3233506361786171D+00,0.5760310502859016D+00,
     2  0.9190368316810619D+00,0.1302902906205909D+01,
     3  0.1619069879996334D+01,0.1718830537890309D+01,
     4  0.1472177810404408D+01,0.8486834738901383D+00,
     5  -.2529109781075723D-01,-.8819324422523910D+00,
     6  -.1413106253674420D+01,-.3549470067671627D-21,
     7  -.5166113598742685D-20,-.9278365219365948D-19,
     8  -.1709926438881821D-17,-.3055585637598904D-16,
     9  -.5137671623001919D-15,-.7973136833084040D-14,
     *  -.1127189808952422D-12,-.1438791880822731D-11,
     1  -.1648475144062559D-10,-.1689217530801441D-09,
     2  -.1545159472044355D-08,-.1260791573310902D-07,
     3  -.9178435898288870D-07,-.5965885287776843D-06,
     4  -.3466115548522703D-05,-.1802327159775486D-04,
     5  -.8398964624001869D-04,-.3512099078028010D-03,
     6  -.1319253643209890D-02,-.4455171382860429D-02,
     7  -.1353215148500085D-01,-.3696677636725438D-01,
     8  -.9075854114766517D-01,-.1999367527259824D+00,
     9  -.3940676946146334D+00,-.6915690114137975D+00,
     *  -.1072165249741431D+01,-.1448939241331172D+01,
     1  -.1665673395378090D+01,-.1546417882837516D+01,
     2  -.9974865447314657D+00,-.1127789975628871D+00,
     3  0.8017024114365698D+00,0.1342654278734629D+01,
     4  0.1229169537962837D+01,0.4945609669830609D+00,
     5  0.3747805252773753D-20,0.5158375106238843D-19,
     6  0.8766630719519895D-18,0.1529029679757545D-16,
     7  0.2586117372577569D-15,0.4115792403860796D-14,
     8  0.6045603963614355D-13,0.8088779203633727D-12,
     9  0.9769303253195242D-11,0.1058703625272422D-09,
     *  0.1025613230201217D-08,0.8862991287233419D-08,
     1  0.6826140602128754D-07,0.4685359835073761D-06,
     2  0.2867457139213452D-05,0.1565993308489240D-04,
     3  0.7638947786411837D-04,0.3331438699742092D-03,
     4  0.1299937732923021D-02,0.4540687497461476D-02,
     5  0.1419919770820721D-01,0.3973121854663923D-01,
     6  0.9934856977602132D-01,0.2214633532896174D+00,
     7  0.4383189682230255D+00,0.7651775395714934D+00,
     8  0.1165470381266563D+01,0.1519891287436058D+01,
     9  0.1636099124304333D+01,0.1331481663666250D+01,
     *  0.5738891132504533D+00,-.4083410530207641D+00/
c 
c        following is block number   2 of length 180 elements
c 
        data www2/
     1  -.1167956027864957D+01,-.1275652134617670D+01,
     2  -.6297789333609238D+00,0.4012961490746016D+00,
     3  0.1167964563499795D+01,-.3515743755294047D-19,
     4  -.4573528568056526D-18,-.7351161013280311D-17,
     5  -.1212790751369735D-15,-.1940415145165884D-14,
     6  -.2921342056178768D-13,-.4059060224046192D-12,
     7  -.5136390241961487D-11,-.5865537260252496D-10,
     8  -.6007653315886471D-09,-.5497231568936410D-08,
     9  -.4483632581012227D-07,-.3255945202679543D-06,
     *  -.2104502981713629D-05,-.1210971307919021D-04,
     1  -.6206331649606334D-04,-.2834589936791624D-03,
     2  -.1154234887698819D-02,-.4191099116559179D-02,
     3  -.1356703870388772D-01,-.3911869621721552D-01,
     4  -.1002890861695480D+00,-.2279074551639218D+00,
     5  -.4567987656231599D+00,-.8010264317417925D+00,
     6  -.1212590599667439D+01,-.1547413715403922D+01,
     7  -.1586135497644613D+01,-.1147495925669038D+01,
     8  -.2619336224692323D+00,0.7254884387264400D+00,
     9  0.1275799086540079D+01,0.1009544859657504D+01,
     *  0.5872842602079024D-01,-.9214468827659209D+00,
     1  -.1195864010823278D+01,-.5400067491377234D+00,
     2  0.2990075178701321D-18,0.3674239852667059D-17,
     3  0.5582379029114391D-16,0.8706736052802324D-15,
     4  0.1317015728931577D-13,0.1874556380913761D-12,
     5  0.2462173009004960D-11,0.2944684822113311D-10,
     6  0.3177090107499908D-09,0.3072924001958064D-08,
     7  0.2653504765510124D-07,0.2040540142620958D-06,
     8  0.1395514635797423D-05,0.8482572211561719D-05,
     9  0.4582114774564331D-04,0.2199798358721650D-03,
     *  0.9386615424387600D-03,0.3559468518892944D-02,
     1  0.1198871783019947D-01,0.3582210548402981D-01,
     2  0.9475087702232646D-01,0.2210647357794457D+00,
     3  0.4523552795156705D+00,0.8044223159311751D+00,
     4  0.1224365189378082D+01,0.1551777191082976D+01,
     5  0.1545937280770670D+01,0.1024190814693798D+01,
     6  0.6493429161895770D-01,-.8966383239498067D+00,
     7  -.1263887874020775D+01,-.7330192934563689D+00,
     8  0.3463484718540720D+00,0.1124455379062621D+01,
     9  0.9312402371667902D+00,-.8360311254495984D-01,
     *  -.1016809348080669D+01,-.2337863238040360D-17,
     1  -.2712052727634616D-16,-.3892683173221800D-15,
     2  -.5736403753981857D-14,-.8198619649203933D-13,
     3  -.1102535131249222D-11,-.1368028213462553D-10,
     4  -.1545200532171133D-09,-.1573869605021194D-08,
     5  -.1436256788292239D-07,-.1169230405614500D-06,
     6  -.8467995358639831D-06,-.5447030301535735D-05,
     7  -.3109111723905320D-04,-.1573925475711905D-03,
     8  -.7063744907915038D-03,-.2809129563412175D-02,
     9  -.9890605975161987D-02,-.3078538160649069D-01,
     *  -.8450279098859780D-01,-.2037513564906484D+00,
     1  -.4288944617201610D+00,-.7804448432988380D+00,
     2  -.1207699111684665D+01,-.1542800059112418D+01,
     3  -.1527342793871598D+01,-.9679602970666029D+00,
     4  0.3218609225565594D-01,0.9745836410455643D+00,
     5  0.1218844414854994D+01,0.5238221219005458D+00,
     6  -.5902159840469066D+00,-.1152980499820211D+01,
     7  -.6205746493986929D+00,0.5124661286308823D+00,
     8  0.1117466848502673D+01,0.5694612644719520D+00,
     9  0.1697587114331162D-16,0.1857958149727147D-15,
     *  0.2517829114518109D-14,0.3503549129875573D-13,
     1  0.4728252575337118D-12,0.6003532550506734D-11,
     2  0.7032047097125326D-10,0.7495569098070264D-09,
     3  0.7201325514868638D-08,0.6194510374403848D-07,
     4  0.4749138249006231D-06,0.3235400240676689D-05,
     5  0.1954774015182072D-04,0.1046061987337046D-03,
     6  0.4953243803655547D-03,0.2073412513000543D-02,
     7  0.7663496457660214D-02,0.2496612142754901D-01,
     8  0.7149722101438237D-01,0.1792416461631799D+00,
     9  0.3908474687964156D+00,0.7338004992515103D+00,
     *  0.1166385307228404D+01,0.1522669987445746D+01,
     1  0.1530464399504841D+01,0.9741549737898157D+00,
     2  -.4563971300262342D-01,-.9930280215568766D+00,
     3  -.1183459486595505D+01,-.4018997897852153D+00,
     4  0.7155180300195878D+00,0.1116647731335210D+01,
     5  0.3733049905335183D+00,-.7539340969732559D+00,
     6  -.1035938980797829D+01,-.1459567105340389D+00,
     7  0.9063714733700704D+00,-.1153620930140625D-15,
     8  -.1190449657643477D-14,-.1522188843000292D-13,
     9  -.1998783622106432D-12,-.2545434452508591D-11,
     *  -.3049423347972374D-10,-.3369281377631038D-09,
     1  -.3386386558710975D-08,-.3066021647445334D-07,
     2  -.2483489984163724D-06,-.1791079239639839D-05,
     3  -.1146289604304087D-04,-.6495253504069382D-04,
     4  -.3252896650810942D-03,-.1437694844560014D-02,
     5  -.5598733580985650D-02,-.1917117339825849D-01,
     6  -.5755490954129034D-01,-.1508428518399073D+00,
     7  -.3428820373287085D+00,-.6691475198772473D+00,
     8  -.1102582585104863D+01,-.1488928797620601D+01,
     9  -.1548251195606072D+01,-.1032844579355623D+01,
     *  -.1229987394949419D-01,0.9675003336610446D+00/
c 
c        following is block number   3 of length 180 elements
c 
        data www3/
     1  0.1172087756901584D+01,0.3642291461146100D+00,
     2  -.7597902635683045D+00,-.1076305234630621D+01,
     3  -.2192905879557470D+00,0.8685066415019972D+00,
     4  0.9109883658920778D+00,-.1559885945716400D+00,
     5  -.1018482865904952D+01,-.5894968345991693D+00,
     6  0.7381213620680036D-15,0.7176722389717785D-14,
     7  0.8652966803961102D-13,0.1071489052827784D-11,
     8  0.1286722422758554D-10,0.1453344963741082D-09,
     9  0.1513513367231833D-08,0.1433121593213751D-07,
     *  0.1221611206016250D-06,0.9307695229202702D-06,
     1  0.6306727539642957D-05,0.3786461204354963D-04,
     2  0.2008851553050573D-03,0.9396789317473140D-03,
     3  0.3867279044919931D-02,0.1396943744651553D-01,
     4  0.4415101572066002D-01,0.1215567513369134D+00,
     5  0.2896506257249598D+00,0.5914260760533314D+00,
     6  0.1018279263086321D+01,0.1437133377513913D+01,
     7  0.1569537865713116D+01,0.1130827367605331D+01,
     8  0.1327114313927507D+00,-.9003351625544266D+00,
     9  -.1181867050370876D+01,-.4009412403961977D+00,
     *  0.7445326205976462D+00,0.1056709353445044D+01,
     1  0.1562458027887061D+00,-.9087178712598413D+00,
     2  -.8127797180402835D+00,0.3342135711327505D+00,
     3  0.1021234063797856D+01,0.3169895893247709D+00,
     4  -.8187227644880437D+00,-.4468086495054536D-14,
     5  -.4090374965530345D-13,-.4647108925040600D-12,
     6  -.5422831902212924D-11,-.6136281681713060D-10,
     7  -.6529461237414793D-09,-.6403615412811975D-08,
     8  -.5707106884051002D-07,-.4575389619577291D-06,
     9  -.3275308572176554D-05,-.2082300162491287D-04,
     *  -.1170977108141752D-03,-.5806030081015672D-03,
     1  -.2531114010756662D-02,-.9673726185808030D-02,
     2  -.3230300497114366D-01,-.9382083001935143D-01,
     3  -.2354935404254855D+00,-.5059540196512244D+00,
     4  -.9163956604376309D+00,-.1363102221650489D+01,
     5  -.1581380562115433D+01,-.1252101005965133D+01,
     6  -.3073007186587574D+00,0.7861239863143732D+00,
     7  0.1199106913507099D+01,0.5000733190849755D+00,
     8  -.6766776839530192D+00,-.1059992238456312D+01,
     9  -.1735712408135104D+00,0.9043360262767195D+00,
     *  0.7668523063640674D+00,-.4132700315288168D+00,
     1  -.9922205881388528D+00,-.1365861418408448D+00,
     2  0.9106463343191834D+00,0.6031797835529426D+00,
     3  0.2568942819512762D-13,0.2212672806811264D-12,
     4  0.2367015500682843D-11,0.2601020150700934D-10,
     5  0.2771201210827631D-09,0.2775673970980703D-08,
     6  0.2561272815445071D-07,0.2146394687558588D-06,
     7  0.1616595460333741D-05,0.1085906484615725D-04,
     8  0.6468155584583275D-04,0.3401102590074369D-03,
     9  0.1572830970189749D-02,0.6374367047919212D-02,
     *  0.2255435406475972D-01,0.6934725306608409D-01,
     1  0.1841128982006138D+00,0.4182202105146390D+00,
     2  0.8013424315697251D+00,0.1264675568128596D+01,
     3  0.1571218247713123D+01,0.1378247880264359D+01,
     4  0.5253529804403189D+00,-.6169328299567909D+00,
     5  -.1202943260184660D+01,-.6468182550294477D+00,
     6  0.5538621298865013D+00,0.1074034583405019D+01,
     7  0.2600792086435986D+00,-.8645063755951172D+00,
     8  -.7750156898949268D+00,0.4147694575040932D+00,
     9  0.9691834884617889D+00,0.5017949881451246D-01,
     *  -.9338162029404410D+00,-.4453120643464170D+00,
     1  0.7457296516794740D+00,-.1407478022622580D-12,
     2  -.1139693017281120D-11,-.1147100973076178D-10,
     3  -.1186053737745664D-09,-.1188823840763259D-08,
     4  -.1119861193116519D-07,-.9713471398927297D-07,
     5  -.7645914819918165D-06,-.5403572474665577D-05,
     6  -.3401292885026767D-04,-.1895119929196742D-03,
     7  -.9300150785166224D-03,-.4002140515033297D-02,
     8  -.1503638642715928D-01,-.4907816628524496D-01,
     9  -.1382879370121109D+00,-.3334375108969801D+00,
     *  -.6789742941895874D+00,-.1142843989792418D+01,
     1  -.1529000535967868D+01,-.1489543603213697D+01,
     2  -.7714478646317131D+00,0.3872722697643476D+00,
     3  0.1168402434592632D+01,0.8209200396531326D+00,
     4  -.3707694418240017D+00,-.1077050824864925D+01,
     5  -.4034488200484534D+00,0.7840020019893693D+00,
     6  0.8261609666558909D+00,-.3500715158037914D+00,
     7  -.9617620568217888D+00,-.4691020025004011D-01,
     8  0.9276992262793366D+00,0.3652653427592031D+00,
     9  -.8004511424358123D+00,-.6122073443912441D+00,
     *  0.7368301916108836D-12,0.5604573649068069D-11,
     1  0.5303132647323197D-10,0.5155091325475600D-09,
     2  0.4856912827523582D-08,0.4298794848244829D-07,
     3  0.3501341943020961D-06,0.2585797249569608D-05,
     4  0.1712549758656655D-04,0.1008624241162873D-03,
     5  0.5247646011786898D-03,0.2398382171245449D-02,
     6  0.9579455928478713D-02,0.3325660471830247D-01,
     7  0.9971019173643705D-01,0.2560024644195797D+00,
     8  0.5560086178322281D+00,0.1002086121966002D+01,
     9  0.1449070835242646D+01,0.1567034085602484D+01,
     *  0.1024705894878996D+01,-.9856964988061317D-01/
  
c 
c        following is block number   4 of length 180 elements
c 
        data www4/
     1  -.1070299788491014D+01,-.9945761293441077D+00,
     2  0.1258564988763547D+00,0.1041094463797768D+01,
     3  0.5863119707139069D+00,-.6501003885181206D+00,
     4  -.8998594790162545D+00,0.2227507600324419D+00,
     5  0.9613154275704662D+00,0.1155671463804413D+00,
     6  -.9055898113391026D+00,-.3608008813410047D+00,
     7  0.8042789273623734D+00,0.5408557602596370D+00,
     8  -.6830266872987908D+00,-.3694289824627343D-11,
     9  -.2637316670201494D-10,-.2343976446193109D-09,
     *  -.2140314372601869D-08,-.1893696778284434D-07,
     1  -.1573267931809212D-06,-.1201961270953944D-05,
     2  -.8318020730321985D-05,-.5155308710033275D-04,
     3  -.2836304451670758D-03,-.1375264107861886D-02,
     4  -.5840034144719909D-02,-.2158641842428872D-01,
     5  -.6898608938850009D-01,-.1890361740929714D+00,
     6  -.4391056921281539D+00,-.8498550392220914D+00,
     7  -.1331392843080391D+01,-.1595307027834172D+01,
     8  -.1260029918081376D+01,-.2375738107075465D+00,
     9  0.8884338899785981D+00,0.1132568082788345D+01,
     *  0.1721493946821226D+00,-.9367529802189891D+00,
     1  -.7820550186463802D+00,0.4500627644638497D+00,
     2  0.9669817485995405D+00,-.3497002433051473D-01,
     3  -.9469497160699611D+00,-.2453824100976222D+00,
     4  0.8620879809529084D+00,0.4216136334531967D+00,
     5  -.7718078913388138D+00,-.5350036704478940D+00,
     6  0.6917736979418472D+00,0.6176185143186129D+00,
     7  0.1777378740639995D-10,0.1189803308543767D-09,
     8  0.9923580913570955D-09,0.8503709175195583D-08,
     9  0.7058652013480918D-07,0.5498642025207562D-06,
     *  0.3935766730402899D-05,0.2548875462974385D-04,
     1  0.1476052892060783D-03,0.7572304453085924D-03,
     2  0.3414424941992403D-02,0.1343574235099590D-01,
     3  0.4580323726663854D-01,0.1341467776146820D+00,
     4  0.3338910181536803D+00,0.6953414341191235D+00,
     5  0.1181724329002706D+01,0.1565267318963269D+01,
     6  0.1451318145719322D+01,0.5983304234494403D+00,
     7  -.6134792010631415D+00,-.1195746332013892D+01,
     8  -.4992886569621275D+00,0.7401911161401029D+00,
     9  0.9528206616740936D+00,-.1795370847428977D+00,
     *  -.9909846122439361D+00,-.2046851053868811D+00,
     1  0.8901511734673179D+00,0.4215371094315929D+00,
     2  -.7801566858619253D+00,-.5323535823357039D+00,
     3  0.7005710851239707D+00,0.5852425275277089D+00,
     4  -.6533725798332196D+00,-.6104964154962510D+00,
     5  0.6279986558212396D+00,-.8219320110144595D-10,
     6  -.5154442296889008D-09,-.4030436694167276D-08,
     7  -.3237994929355050D-07,-.2518906056612360D-06,
     8  -.1837758090513582D-05,-.1230816652419308D-04,
     9  -.7448547544945482D-04,-.4023566114521196D-03,
     *  -.1920858929418320D-02,-.8034930270047273D-02,
     1  -.2920880161760330D-01,-.9147350614617693D-01,
     2  -.2442009244663335D+00,-.5478413984045400D+00,
     3  -.1010574254663911D+01,-.1476038936524016D+01,
     4  -.1575988545000137D+01,-.9515951269254367D+00,
     5  0.2520455571104147D+00,0.1148022300874036D+01,
     6  0.8148803868869218D+00,-.4421388830146149D+00,
     7  -.1051927287836711D+01,-.1489283962777269D+00,
     8  0.9326172258080271D+00,0.4737163076577744D+00,
     9  -.7608890718353171D+00,-.6192781585978773D+00,
     *  0.6383226482037594D+00,0.6702132414860110D+00,
     1  -.5783477386771121D+00,-.6737844024114142D+00,
     2  0.5674312223733713D+00,0.6528115986945231D+00,
     3  -.5870173274452998D+00,-.6200889155949484D+00,
     4  0.3658582769020855D-09,0.2147204321840412D-08,
     5  0.1572412444401086D-07,0.1183077201603988D-06,
     6  0.8615490257824169D-06,0.5879786537207020D-05,
     7  0.3679558875971562D-04,0.2077511259447160D-03,
     8  0.1044866103975250D-02,0.4631601056097609D-02,
     9  0.1792310028537459D-01,0.5997861446921160D-01,
     *  0.1717434712854339D+00,0.4151719931687760D+00,
     1  0.8311401851888702D+00,0.1335262262996261D+01,
     2  0.1619581119435315D+01,0.1261007841103145D+01,
     3  0.1710811662096471D+00,-.9656281683095705D+00,
     4  -.1065358885309242D+01,0.5663779253284116D-01,
     5  0.1032096668696560D+01,0.4997022573371339D+00,
     6  -.7594141975846935D+00,-.7311302265353488D+00,
     7  0.5369448093044025D+00,0.7998717953768159D+00,
     8  -.4196197683056463D+00,-.8016286729623616D+00,
     9  0.3916716023193492D+00,0.7737071448434782D+00,
     *  -.4255205362406647D+00,-.7256579142243688D+00,
     1  0.4953107228987033D+00,0.6593004823176414D+00,
     2  -.5789505007729807D+00,-.1569395289538605D-08,
     3  -.8610874822311459D-08,-.5898961462861121D-07,
     4  -.4151925614700964D-06,-.2826953825598411D-05,
     5  -.1802286642643117D-04,-.1052280135850978D-03,
     6  -.5533379091795049D-03,-.2585751670752947D-02,
     7  -.1061569059756303D-01,-.3788327034635048D-01,
     8  -.1162198835702757D+00,-.3025464775785302D+00,
     9  -.6567586801854982D+00,-.1157522757519471D+01,
     *  -.1578980196986266D+01,-.1493063209809347D+01/
  
c 
c        following is block number 5 of length 180 elements
c 
        data www5/
     1  -.6151344706337315D+00,0.6459865977193222D+00,
     2  0.1194085984588603D+01,0.3740433845871656D+00,
     3  -.8589857817691345D+00,-.8138236823606808D+00,
     4  0.4590785454795069D+00,0.9191082546332291D+00,
     5  -.2160422148556008D+00,-.9117616958616214D+00,
     6  0.1226288378130736D+00,0.8819872506196583D+00,
     7  -.1359701794959343D+00,-.8478426671108373D+00,
     8  0.2202524362782217D+00,0.8003217327891143D+00,
     9  -.3464006176552922D+00,-.7263239272701271D+00,
     *  0.4876901436252612D+00,0.6200707818023282D+00,
     1  0.6494322539993899D-08,0.3327489750435725D-07,
     2  0.2129881714229042D-06,0.1400630922999393D-05,
     3  0.8904788813501369D-05,0.5295662032793951D-04,
     4  0.2879935577141799D-03,0.1407715603260948D-02,
     5  0.6097955395902086D-02,0.2312015503507545D-01,
     6  0.7580664397851095D-01,0.2121487458970634D+00,
     7  0.4985577708232785D+00,0.9612692604182907D+00,
     8  0.1463113256673591D+01,0.1624376478847931D+01,
     9  0.1030353356314536D+01,-.2127918081365393D+00,
     *  -.1155582875980852D+01,-.7820580939617058D+00,
     1  0.5267563713240564D+00,0.1017945924064839D+01,
     2  -.5289230518849239D-01,-.9737629921499054D+00,
     3  -.1725791683074270D+00,0.8998316440924969D+00,
     4  0.2275597586597319D+00,-.8617285976218480D+00,
     5  -.1735507235127739D+00,0.8514121443422837D+00,
     6  0.4424021561477054D-01,-.8385467570536987D+00,
     7  0.1360401190783044D+00,0.7918447863805431D+00,
     8  -.3397620818729323D+00,-.6911567516943148D+00,
     9  0.5347148153794550D+00,-.2594677467858536D-07,
     *  -.1239993556265821D-06,-.7406307646917918D-06,
     1  -.4544517640366184D-05,-.2694004292236002D-04,
     2  -.1492093576453076D-03,-.7544358734310694D-03,
     3  -.3420533950468783D-02,-.1369933990148277D-01,
     4  -.4780994720926139D-01,-.1434060689966377D+00,
     5  -.3639026219980785D+00,-.7650960943112255D+00,
     6  -.1290787923116648D+01,-.1646842040920935D+01,
     7  -.1368335776144818D+01,-.2855133764618137D+00,
     8  0.9302042345180822D+00,0.1086145932809213D+01,
     9  -.6985749952601580D-01,-.1043127392744593D+01,
     *  -.3972513666223539D+00,0.8450843418683414D+00,
     1  0.5631749683844298D+00,-.7233549967077853D+00,
     2  -.5734466246756556D+00,0.7005808487317778D+00,
     3  0.4937140753165585D+00,-.7414672103334288D+00,
     4  -.3421743367308659D+00,0.7969731198405260D+00,
     5  0.1271877713227948D+00,-.8172802774063322D+00,
     6  0.1309981689875487D+00,0.7630589884496917D+00,
     7  -.3947290088725394D+00,-.6178645960581906D+00,
     8  0.1001565148384067D-06,0.4458853947204619D-06,
     9  0.2481613041053396D-05,0.1418770064296700D-04,
     *  0.7829910544377395D-04,0.4031787234823570D-03,
     1  0.1891517436439809D-02,0.7935479510561948D-02,
     2  0.2929670369207961D-01,0.9375187591460384D-01,
     3  0.2559192902522168D+00,0.5844431169730671D+00,
     4  0.1086315924593825D+01,0.1568981389043754D+01,
     5  0.1592809480987497D+01,0.7846431663250402D+00,
     6  -.5339060303440850D+00,-.1211691723586929D+01,
     7  -.4355252648540299D+00,0.8496038635644539D+00,
     8  0.7949093940066053D+00,-.5213825422368942D+00,
     9  -.8617072746965648D+00,0.3787189301165352D+00,
     *  0.8307022837964134D+00,-.3879848315471542D+00,
     1  -.7535871363430242D+00,0.4936607939866470D+00,
     2  0.6205566127485724D+00,-.6385001017386333D+00,
     3  -.4122511232409023D+00,0.7599417237913452D+00,
     4  0.1295783348098046D+00,-.7970129357776876D+00,
     5  0.1931574743002435D+00,0.7091293754297043D+00,
     6  -.4944458507735643D+00,-.3737281908937944D-06,
     7  -.1547854156059518D-05,-.8014916572790200D-05,
     8  -.4262710109606995D-04,-.2186368956245627D-03,
     9  -.1044647062037252D-02,-.4537210176258745D-02,
     *  -.1756529577661705D-01,-.5957230797540552D-01,
     1  -.1740094743756424D+00,-.4295403249449734D+00,
     2  -.8743107772104877D+00,-.1413079652973657D+01,
     3  -.1687398504545303D+01,-.1218709909320143D+01,
     4  0.1854356816987856D-01,0.1113019849772143D+01,
     5  0.8844741429768195D+00,-.4490097665879411D+00,
     6  -.1031505972041777D+01,0.4851586180481978D-01,
     7  0.9705700268312298D+00,0.8229723499145267D-01,
     8  -.9096296336681647D+00,-.3752048503310198D-01,
     9  0.8687606610122161D+00,-.1223222106236938D+00,
     *  -.8041305240644152D+00,0.3501979508872338D+00,
     1  0.6617035997649543D+00,-.5852871864326525D+00,
     2  -.4100098174955356D+00,0.7495137077527293D+00,
     3  0.6425787771574778D-01,-.7700365050877849D+00,
     4  0.3086966887214634D+00,0.6136556228521829D+00,
     5  0.1348609595310324D-05,0.5188928404899194D-05,
     6  0.2495574676984979D-04,0.1232584116107033D-03,
     7  0.5864477175732775D-03,0.2594469696174952D-02,
     8  0.1040568373556854D-01,0.3705793689038143D-01,
     9  0.1149906475448936D+00,0.3049140744377468D+00,
     *  0.6752112673199199D+00,0.1209200100753764D+01/
  
c 
c        following is block number 6 of length 180 elements
c 
        data www6/
     1  0.1657914166196813D+01,0.1535413208928588D+01,
     2  0.5387486804070905D+00,-.7898533646149452D+00,
     3  -.1169460282996070D+01,-.8647723672031326D-01,
     4  0.1021319746732472D+01,0.4678205727438545D+00,
     5  -.8256723736698031D+00,-.5508604058997288D+00,
     6  0.7503525712649021D+00,0.4815816511613217D+00,
     7  -.7710045523704610D+00,-.3055779343781105D+00,
     8  0.8154455017103572D+00,0.3767328930090577D-01,
     9  -.7994303122755136D+00,0.2854741242908678D+00,
     *  0.6516585038426067D+00,-.5846372129978169D+00,
     1  -.3485316470138259D+00,0.7552246175444275D+00,
     2  -.5959880503752474D-01,-.7156682779768237D+00,
     3  0.4575014116104749D+00,-.4707408899970804D-05,
     4  -.1680120165700808D-04,-.7491151293187084D-04,
     5  -.3429460771279318D-03,-.1510446645906995D-02,
     6  -.6172373618335044D-02,-.2279403400289728D-01,
     7  -.7440455566239474D-01,-.2102278968494165D+00,
     8  -.5026393515683927D+00,-.9881670411970339D+00,
     9  -.1528719738779371D+01,-.1707171469976906D+01,
     *  -.1054246805053616D+01,0.2908667849281890D+00,
     1  0.1211855260845011D+01,0.6355932073302089D+00,
     2  -.7373429713247402D+00,-.8844882221012156D+00,
     3  0.4327067206183988D+00,0.8851312318063487D+00,
     4  -.3600421830690267D+00,-.8112021905095620D+00,
     5  0.4458900167966013D+00,0.6749212914216471D+00,
     6  -.6085962643638603D+00,-.4422776036080374D+00,
     7  0.7540271351475380D+00,0.1020127114113426D+00,
     8  -.7777879314716272D+00,0.2931451170835424D+00,
     9  0.6024884237026652D+00,-.6239557460187662D+00,
     *  -.2308831291027487D+00,0.7535925220849591D+00,
     1  -.2299110567064298D+00,-.6075285212007299D+00,
     2  0.1589613129029952D-04,0.5254445823114054D-04,
     3  0.2167477901811808D-03,0.9177875782113657D-03,
     4  0.3733118520633045D-02,0.1405277887074118D-01,
     5  0.4762480088670929D-01,0.1418811206525071D+00,
     6  0.3629091028653561D+00,0.7756969780612405D+00,
     7  0.1334872641503207D+01,0.1734732877191047D+01,
     8  0.1458115934498565D+01,0.2981440717184017D+00,
     9  -.9871544449064021D+00,-.1059622862080848D+01,
     *  0.2321228540945412D+00,0.1065477119947861D+01,
     1  0.1158398439518964D+00,-.9585599449519617D+00,
     2  -.1669849530308642D+00,0.8987001995995444D+00,
     3  0.3870029755556563D-01,-.8539968761181683D+00,
     4  0.2052365863455155D+00,0.7408311445828815D+00,
     5  -.4945280986032723D+00,-.4872567606980925D+00,
     6  0.7183841017481599D+00,0.8803301074826193D-01,
     7  -.7482372628453276D+00,0.3595973829794742D+00,
     8  0.5103841365526674D+00,-.6799071365520021D+00,
     9  -.5848006554747559D-01,0.7127369932207470D+00,
     *  -.4233684173425961D+00,-.5192737286662271D-04,
     1  -.1587001242559871D-03,-.6042625680286599D-03,
     2  -.2360968602227872D-02,-.8845475171369162D-02,
     3  -.3057737499212362D-01,-.9473003668581339D-01,
     4  -.2562542688254920D+00,-.5890937086504014D+00,
     5  -.1113036784301676D+01,-.1642802055863754D+01,
     6  -.1709648013813761D+01,-.8776726007945510D+00,
     7  0.5334626308447211D+00,0.1244164402252787D+01,
     8  0.3688996025393834D+00,-.9336584272434046D+00,
     9  -.6575817859326154D+00,0.7163302073499136D+00,
     *  0.6642690555806890D+00,-.6802478548208014D+00,
     1  -.5340933474765653D+00,0.7492927849438955D+00,
     2  0.2873996270865710D+00,-.8112931152870011D+00,
     3  0.6447536621415138D-01,0.7449057552383721D+00,
     4  -.4472625069761487D+00,-.4681560126134094D+00,
     5  0.7117728185544972D+00,0.8854217333132047D-02,
     6  -.7005221324287463D+00,0.4656421171115588D+00,
     7  0.3648884172325551D+00,-.7193047732538397D+00,
     8  0.1585372389915567D+00,0.5994632126551478D+00,
     9  0.1640606067907842D-03,0.4627600523213530D-03,
     *  0.1622195730549795D-02,0.5832766160745940D-02,
     1  0.2006758571575796D-01,0.6347342346928671D-01,
     2  0.1789360893228807D+00,0.4367719857777806D+00,
     3  0.8940823776263601D+00,0.1470194196451574D+01,
     4  0.1802610538962184D+01,0.1360871902085931D+01,
     5  0.6067244037461318D-01,-.1134881661168101D+01,
     6  -.9029511045121607D+00,0.5071443194846208D+00,
     7  0.1011679674002935D+01,-.2133330577155667D+00,
     8  -.9477663310309129D+00,0.2021760992897097D+00,
     9  0.8547358363527003D+00,-.3624343515247819D+00,
     *  -.6966506358221596D+00,0.5919027537196529D+00,
     1  0.4078377620868296D+00,-.7599126519583254D+00,
     2  0.1604249819252408D-01,0.7202552375961910D+00,
     3  -.4613578294183719D+00,-.3957537165449317D+00,
     4  0.7193261583187537D+00,-.1245179834548382D+00,
     5  -.6137255903095414D+00,0.5830611219769038D+00,
     6  0.1598111649769613D+00,-.7018831268381760D+00,
     7  0.3916088451995342D+00,-.5011159064427980D-03,
     8  -.1302092537529811D-02,-.4190047370295036D-02,
     9  -.1382154670218891D-01,-.4351607524502258D-01,
     *  -.1254050939573858D+00,-.3199070740600934D+00/
c 
c        following is block number 7 of length 180 elements
c 
        data www7/
     1  -.6990934561970858D+00,-.1258594296900062D+01,
     2  -.1760306327988728D+01,-.1693444306723797D+01,
     3  -.6841959287505903D+00,0.7525687174068808D+00,
     4  0.1221958903783947D+01,0.9883673762564373D-01,
     5  -.1046940136913611D+01,-.3905180172991969D+00,
     6  0.8925746375539425D+00,0.3749345976651137D+00,
     7  -.8554424816104453D+00,-.1834582571520208D+00,
     8  0.8480835411910765D+00,-.1322846180235744D+00,
     9  -.7420400521486957D+00,0.4894650330536664D+00,
     *  0.4386686879279083D+00,-.7288777405158940D+00,
     1  0.4400636301005340D-01,0.6731543466300474D+00,
     2  -.5203313479239026D+00,-.2702330802988185D+00,
     3  0.7143092223762893D+00,-.2956839955577534D+00,
     4  -.4654679874088817D+00,0.6719763901602314D+00,
     5  -.9466222130950613D-01,-.5893062171359880D+00,
     6  0.1478833469762111D-02,0.3532758828157223D-02,
     7  0.1040094340104506D-01,0.3136326005784848D-01,
     8  0.8999227408890160D-01,0.2350833063351987D+00,
     9  0.5389389422547142D+00,0.1043647822429611D+01,
     *  0.1623858056380115D+01,0.1861104873322325D+01,
     1  0.1235921435732545D+01,-.1817963654989027D+00,
     2  -.1239769285642340D+01,-.7088297124125551D+00,
     3  0.7363544795342768D+00,0.8837517611548340D+00,
     4  -.4963521952042582D+00,-.8287200917697791D+00,
     5  0.5084257168019726D+00,0.6747126776254917D+00,
     6  -.6514535709593858D+00,-.3991757340029964D+00,
     7  0.7809437491216877D+00,-.9826495384195670D-02,
     8  -.7340355534591797D+00,0.4543730098827148D+00,
     9  0.4050399288817407D+00,-.7176135080906227D+00,
     *  0.1337754039428894D+00,0.5930196675449808D+00,
     1  -.5988732072050938D+00,-.9126754624810520D-01,
     2  0.6631326192615455D+00,-.4749221111022150D+00,
     3  -.2438833105546876D+00,0.6842528311914117D+00,
     4  -.3618088382031753D+00,-.4212513130569039D-02,
     5  -.9232674979613458D-02,-.2477432436955186D-01,
     6  -.6800287839431725D-01,-.1769582662271223D+00,
     7  -.4164060471846654D+00,-.8503885173951816D+00,
     8  -.1439020808372157D+01,-.1884409364033647D+01,
     9  -.1649251033148742D+01,-.4608090287386733D+00,
     *  0.9550011598105018D+00,0.1148166799017688D+01,
     1  -.1709459346115384D+00,-.1085784729425492D+01,
     2  -.1065832541649197D+00,0.9681315792231672D+00,
     3  0.6412898417248127D-01,-.8946496878124789D+00,
     4  0.1643472748176547D+00,0.7728230907616039D+00,
     5  -.4783090199198544D+00,-.4887478905918898D+00,
     6  0.7225531383349792D+00,0.1997004588187136D-01,
     7  -.7009316026694230D+00,0.4772651443867975D+00,
     8  0.3163898654328724D+00,-.7079966758201645D+00,
     9  0.2692266840090674D+00,0.4625127993999689D+00,
     *  -.6618897660933354D+00,0.1305744576847783D+00,
     1  0.5338876971888708D+00,-.6156370051209955D+00,
     2  0.3837005629579650D-01,0.5766942382908483D+00,
     3  0.1156719955491533D-01,0.2321011274343100D-01,
     4  0.5651105159101926D-01,0.1404863829803059D+00,
     5  0.3295549322462606D+00,0.6930499516542948D+00,
     6  0.1246101570790276D+01,0.1805642015325810D+01,
     7  0.1902033177999305D+01,0.1064414449546494D+01,
     8  -.4417489839406095D+00,-.1301353738389627D+01,
     9  -.4752538949760369D+00,0.9201605364212545D+00,
     *  0.6944833385786625D+00,-.7233160652907645D+00,
     1  -.6301182156105312D+00,0.7302852343573940D+00,
     2  0.4116614900326117D+00,-.8051448736434581D+00,
     3  -.5230968984104097D-01,0.7759120470005064D+00,
     4  -.3809323193105565D+00,-.4969598810802904D+00,
     5  0.6934988935898223D+00,-.2457305880035294D-01,
     6  -.6445542527903632D+00,0.5381201904391619D+00,
     7  0.1743313251979230D+00,-.6714242810886409D+00,
     8  0.4256343304223285D+00,0.2689649233212775D+00,
     9  -.6666250731715980D+00,0.3639208999349503D+00,
     *  0.3105512202316869D+00,-.6605128716333911D+00,
     1  0.3335054116743334D+00,-.3056053276089139D-01,
     2  -.5601863478815348D-01,-.1231142770747138D+00,
     3  -.2754757485986950D+00,-.5781464408011263D+00,
     4  -.1075351975752594D+01,-.1674416546474150D+01,
     5  -.2011741119801901D+01,-.1553212838136637D+01,
     6  -.1845949268429354D+00,0.1143179984363706D+01,
     7  0.1012658111954592D+01,-.4425782915476415D+00,
     8  -.1050502597530857D+01,0.1831075498677203D+00,
     9  0.9492495097620265D+00,-.2391644566872959D+00,
     *  -.8157839533119837D+00,0.4623695571691412D+00,
     1  0.5705270683315868D+00,-.6987217568908668D+00,
     2  -.1497495745165436D+00,0.7468708850750366D+00,
     3  -.3597776850914934D+00,-.4476718354992130D+00,
     4  0.6866988560968892D+00,-.1278377780913367D+00,
     5  -.5516231287957181D+00,0.6077066222326141D+00,
     6  -.1620264447570847D-01,-.5761036445470166D+00,
     7  0.5641796341070490D+00,0.1758762047469329D-01,
     8  -.5724618127137065D+00,0.5534956978875837D+00,
     9  0.1015010577206929D-01,-.5608471230303807D+00,
     *  0.7746936404418255D-01,0.1294514700897340D+00/
  
c 
c        following is block number 8 of length 180 elements
c 
        data www8/
     1  0.2552355944093297D+00,0.5100252544346579D+00,
     2  0.9482087563222563D+00,0.1537641902837798D+01,
     3  0.2022077107944065D+01,0.1898808753682047D+01,
     4  0.8092566326490843D+00,-.7328075515720888D+00,
     5  -.1304042752368791D+01,-.1890518167074209D+00,
     6  0.1052553453813897D+01,0.4436868964871803D+00,
     7  -.8860208089669544D+00,-.3686101626122485D+00,
     8  0.8569478419965463D+00,0.1052337716461684D+00,
     9  -.8229288921300610D+00,0.2780570514584103D+00,
     *  0.6143261863662071D+00,-.6308142499834145D+00,
     1  -.1552155201931309D+00,0.7070215086875513D+00,
     2  -.3992710706892087D+00,-.3470813335671321D+00,
     3  0.6791659247125093D+00,-.2716457441331986D+00,
     4  -.4050333475340314D+00,0.6475163662211284D+00,
     5  -.2355308952453938D+00,-.3986187620935694D+00,
     6  0.6335342264172942D+00,-.2568531931090506D+00,
     7  -.3594819379919693D+00,0.6305018783012455D+00,
     8  -.3060071958652730D+00,-.1875820118055840D+00,
     9  -.2851975595830044D+00,-.5008512260053845D+00,
     *  -.8847842084376594D+00,-.1437470756794964D+01,
     1  -.1989433405554128D+01,-.2116664808899325D+01,
     2  -.1350555058335013D+01,0.1795526339307495D+00,
     3  0.1303762150014085D+01,0.7848863254800752D+00,
     4  -.7160913271757107D+00,-.9243491542239299D+00,
     5  0.4718943065212517D+00,0.8307522118103458D+00,
     6  -.5148761309221716D+00,-.6309037865025858D+00,
     7  0.6805380109392950D+00,0.2838497396992169D+00,
     8  -.7736807628068010D+00,0.1944362286422398D+00,
     9  0.5936466742267666D+00,-.6090566836433076D+00,
     *  -.8957925838247913D-01,0.6507380581666180D+00,
     1  -.4757266002557630D+00,-.1945710284913175D+00,
     2  0.6387958426648291D+00,-.4272186844500334D+00,
     3  -.1966880024588833D+00,0.6146128342448148D+00,
     4  -.4415317257569420D+00,-.1400062042618868D+00,
     5  0.5831618288029610D+00,-.4875901601297967D+00,
     6  -.5034775209032997D-01,0.5398422477823570D+00,
     7  0.4302879771860331D+00,0.5944066888437048D+00,
     8  0.9220421041767678D+00,0.1420701962211284D+01,
     9  0.1977275173638194D+01,0.2254691362569372D+01,
     *  0.1774048931778339D+01,0.3995100059986310D+00,
     1  -.1057574042355641D+01,-.1197796887207693D+01,
     2  0.1726692267932930D+00,0.1104697296850059D+01,
     3  0.1177121368247339D+00,-.9619644967024460D+00,
     4  -.4925490559793620D-01,0.8704786730976642D+00,
     5  -.2156758816510566D+00,-.7051760413882513D+00,
     6  0.5440009400808059D+00,0.3387074954869055D+00,
     7  -.7229348125385055D+00,0.1965194905141914D+00,
     8  0.5277017973911690D+00,-.6169111529370106D+00,
     9  0.3240894292088457D-01,0.5583692136954410D+00,
     *  -.5561105323498039D+00,0.4120030735203543D-02,
     1  0.5326881318526853D+00,-.5504774433171196D+00,
     2  0.5838973576996338D-01,0.4754583499950482D+00,
     3  -.5720765522660159D+00,0.1592870876931934D+00,
     4  0.3886390373687218D+00,-.5917085420885882D+00,
     5  0.2777140069445294D+00,-.9191081685009258D+00,
     6  -.1153214351274300D+01,-.1565534564067097D+01,
     7  -.2065486007753757D+01,-.2385346212013025D+01,
     8  -.2091210215327648D+01,-.9117598982843072D+00,
     9  0.6737235477108288D+00,0.1371481967904694D+01,
     *  0.3985614009629245D+00,-.9673439611708624D+00,
     1  -.6572381911479182D+00,0.7398922909030486D+00,
     2  0.5828564496649169D+00,-.7330617711527658D+00,
     3  -.3371282113571682D+00,0.7837472244706101D+00,
     4  -.5314163207689307D-01,-.6889067632664572D+00,
     5  0.4787933950028658D+00,0.3072771252489702D+00,
     6  -.6788899672240445D+00,0.2617718405095195D+00,
     7  0.4154731839418278D+00,-.6236541251464960D+00,
     8  0.1921479807522147D+00,0.4090603004028455D+00,
     9  -.5973065746594473D+00,0.2238890894184803D+00,
     *  0.3416639026099896D+00,-.5872758686634310D+00,
     1  0.3118718153495547D+00,0.2281334058421808D+00,
     2  -.5643245471361249D+00,0.4174195855950420D+00,
     3  0.8040505628316702D-01,-.5080099307913500D+00,
     4  0.1766814614644337D+01,0.2016104317495241D+01,
     5  0.2372061188984965D+01,0.2609540158320411D+01,
     6  0.2356961405060046D+01,0.1310449164108145D+01,
     7  -.2812069840034042D+00,-.1352579769279805D+01,
     8  -.8556608918314116D+00,0.6216165745315250D+00,
     9  0.9880278266183715D+00,-.2990041411795834D+00,
     *  -.8854561354290778D+00,0.3208802189422670D+00,
     1  0.7238915753764700D+00,-.5141498743439094D+00,
     2  -.4388546715049648D+00,0.6924409650386187D+00,
     3  -.1030622743461405D-01,-.6317379845563086D+00,
     4  0.4752465877650908D+00,0.2110710503888363D+00,
     5  -.6255159619429039D+00,0.3622279339814799D+00,
     6  0.2527620572696162D+00,-.5922967961616812D+00,
     7  0.3570594572425016D+00,0.1976918176845541D+00,
     8  -.5558178787510763D+00,0.4147226933387033D+00,
     9  0.8006447960746360D-01,-.4955832001794641D+00,
     *  0.4885244423982859D+00,-.7742375781275533D-01/
  
c 
c        following is block number 9 of length 180 elements
c 
        data www9/
     1  -.3921727629927899D+00,0.5378578814074837D+00,
     2  -.2456075710481404D+00,-.2918338048391883D+01,
     3  -.3034610246953614D+01,-.3057380887498579D+01,
     4  -.2681022503968316D+01,-.1629778399004426D+01,
     5  -.3416838902164924D-01,0.1260621666368625D+01,
     6  0.1171817898036141D+01,-.2161153939147167D+00,
     7  -.1098496797624310D+01,-.1822089159893226D+00,
     8  0.9180390462984134D+00,0.1691537023518981D+00,
     9  -.8303353624701186D+00,0.6096172004738551D-01,
     *  0.7218282601828371D+00,-.3933534210683229D+00,
     1  -.4370400255465089D+00,0.6401650662848070D+00,
     2  -.5265442180188327D-01,-.5533448245900135D+00,
     3  0.5107000340341383D+00,0.6728312112761596D-01,
     4  -.5446687216002560D+00,0.4658382915694571D+00,
     5  0.4905966959878381D-01,-.4970871618106331D+00,
     6  0.4849533728691938D+00,-.5156291382102654D-01,
     7  -.4125003043069974D+00,0.5198031711324057D+00,
     8  -.1976538201447958D+00,-.2784558805807067D+00,
     9  0.5239039825407135D+00,-.3499615573702791D+00,
     *  -.9900235422584095D-01,0.4676617488528993D+00,
     1  0.4049700831592706D+01,0.3830033369983175D+01,
     2  0.3223837411415272D+01,0.2022717230795254D+01,
     3  0.3207551582237040D+00,-.1182547770602693D+01,
     4  -.1438489020139304D+01,-.2007692076354944D+00,
     5  0.1070136134057731D+01,0.6323324722594845D+00,
     6  -.7504827398502840D+00,-.6175431268265433D+00,
     7  0.6876495263223741D+00,0.4149385269921435D+00,
     8  -.7396592254784236D+00,-.6265783388088699D-01,
     9  0.7007224296876174D+00,-.3656458453775710D+00,
     *  -.3932616136204247D+00,0.6327309264655472D+00,
     1  -.1494810449927383D+00,-.4609898201294823D+00,
     2  0.5638882135902659D+00,-.1079662773914290D+00,
     3  -.4266640956029669D+00,0.5459048294005962D+00,
     4  -.1754874483500837D+00,-.3323258989712838D+00,
     5  0.5425765444951590D+00,-.2988108853362352D+00,
     6  -.1803982179428835D+00,0.5049423463161422D+00,
     7  -.4263956909118109D+00,0.1760057311221846D-01,
     8  0.3998028098245213D+00,-.5071666354553527D+00,
     9  0.2265340959457309D+00,-.4831341472939090D+01,
     *  -.4109651143384095D+01,-.2723380350788557D+01,
     1  -.7974348724261039D+00,0.1032799503768320D+01,
     2  0.1715484638016471D+01,0.7105438746300046D+00,
     3  -.8941235598507432D+00,-.1061516058280806D+01,
     4  0.3833825329487314D+00,0.9756010188097783D+00,
     5  -.3072577266448802D+00,-.8151965600846212D+00,
     6  0.4729246415466327D+00,0.5524423106335208D+00,
     7  -.6836753971272618D+00,-.1138473315339578D+00,
     8  0.6963405020624226D+00,-.3992805439526372D+00,
     9  -.3270215807414314D+00,0.6489881642495509D+00,
     *  -.2815253299222028D+00,-.3413781149278530D+00,
     1  0.6056161470030407D+00,-.3021910644529106D+00,
     2  -.2545816594116374D+00,0.5690669759556401D+00,
     3  -.3957585329764919D+00,-.9457708286143845D-01,
     4  0.4945263166520010D+00,-.4962397458901114D+00,
     5  0.1151319597618804D+00,0.3462933113291743D+00,
     6  -.5398141712163532D+00,0.3279661371047520D+00,
     7  0.1255036197648641D+00,-.4837730763458327D+00,
     8  0.5211734801126140D+01,0.3905520999442787D+01,
     9  0.1775516093367901D+01,-.5261038890637303D+00,
     *  -.1871699617959158D+01,-.1362435393988294D+01,
     1  0.4285610449515197D+00,0.1371650885574276D+01,
     2  0.2336076081636103D+00,-.1096142856264671D+01,
     3  -.3071841448642183D+00,0.9628958496411746D+00,
     4  0.8719548020957059D-01,-.8761886418652211D+00,
     5  0.2941473262134436D+00,0.6316439088328793D+00,
     6  -.6545695971334199D+00,-.1197488692712087D+00,
     7  0.6923754146235923D+00,-.4661688592383947D+00,
     8  -.2278780934177084D+00,0.6502880553281526D+00,
     9  -.4288637749511712D+00,-.1723210922168585D+00,
     *  0.5908447879318905D+00,-.4856168391114459D+00,
     1  -.1743336895045146D-01,0.4872709511965409D+00,
     2  -.5575653117058082D+00,0.1961410756256412D+00,
     3  0.3040749357883315D+00,-.5666138199594549D+00,
     4  0.4080215536658739D+00,0.4453217957755151D-01,
     5  -.4595428361996576D+00,0.5417732055349237D+00,
     6  -.2367172686564486D+00,-.5342731454993147D+01,
     7  -.3430262854047045D+01,-.7033946469637327D+00,
     8  0.1560861745045864D+01,0.1997738745656728D+01,
     9  0.4273196470681696D+00,-.1311390761246437D+01,
     *  -.9970530681562964D+00,0.7645193079017206D+00,
     1  0.9666467630534920D+00,-.6430818516619797D+00,
     2  -.7509266612232672D+00,0.7586389440476989D+00,
     3  0.3736023423960686D+00,-.8573274139214244D+00,
     4  0.1651292838909696D+00,0.6648377680323878D+00,
     5  -.6472215051745156D+00,-.8163723476375782D-01,
     6  0.6693675459739860D+00,-.5490059091603961D+00,
     7  -.8583736521721505D-01,0.6040341557170333D+00,
     8  -.5634049426199418D+00,0.4777535035451208D-01,
     9  0.4830984132831847D+00,-.6072899129859108D+00,
     *  0.2563899501913887D+00,0.2764850963518055D+00/
c 
c        following is block number 10 of length 180 elements
c 
        data www10/
     1  -.5886414909920528D+00,0.4647817282783468D+00,
     2  -.1125628104066205D-01,-.4415829222477316D+00,
     3  0.5790222154524389D+00,-.3137514313130358D+00,
     4  -.1653309766187957D+00,0.5284569900925752D+00,
     5  0.5357943685186967D+01,0.2843253764288105D+01,
     6  -.2972137276308305D+00,-.2163706795980860D+01,
     7  -.1536040586445373D+01,0.6199137060566464D+00,
     8  0.1558372720052655D+01,0.8266830430302820D-01,
     9  -.1277956597477451D+01,-.1528124158022829D+00,
     *  0.1104992520569898D+01,-.1061129190915107D+00,
     1  -.9297170203659215D+00,0.5138362531090301D+00,
     2  0.5433326569321447D+00,-.8117988300571108D+00,
     3  0.9860750702421187D-01,0.6610095840401607D+00,
     4  -.6632445182824290D+00,0.3177453982593444D-02,
     5  0.6132990573492814D+00,-.6318632551463060D+00,
     6  0.9824894900738804D-01,0.4883044079935156D+00,
     7  -.6496803517829702D+00,0.2942197341481302D+00,
     8  0.2703277122958035D+00,-.6131470379439777D+00,
     9  0.5005353782068456D+00,-.3652074588924764D-01,
     *  -.4407517245941232D+00,0.6053493411876363D+00,
     1  -.3575884445487741D+00,-.1255314808665815D+00,
     2  0.5195872741543334D+00,-.5675069000338581D+00,
     3  0.2419883983343572D+00,-.5317848079545108D+01,
     4  -.2217752055738074D+01,0.1143821928799593D+01,
     5  0.2339063961315599D+01,0.7256412669251432D+00,
     6  -.1409283037095012D+01,-.1155250005787810D+01,
     7  0.8524310940356937D+00,0.1054222686307659D+01,
     8  -.7583773008533094D+00,-.7699400483999265D+00,
     9  0.8884734851727772D+00,0.3012516959850583D+00,
     *  -.9425972298626577D+00,0.3173414033605738D+00,
     1  0.6219837392690592D+00,-.7759960542154201D+00,
     2  0.9671395237404079D-01,0.6236706571432725D+00,
     3  -.6968141617042254D+00,0.1338187443000640D+00,
     4  0.5084717452793898D+00,-.6898451052058332D+00,
     5  0.3096287024510655D+00,0.2904352764371268D+00,
     6  -.6457967677413750D+00,0.5174979276120817D+00,
     7  -.2812077504956665D-01,-.4642027124468686D+00,
     8  0.6279062937278663D+00,-.3703603353557414D+00,
     9  -.1235408116223296D+00,0.5284014055759818D+00,
     *  -.5896496836270778D+00,0.2766971255598840D+00,
     1  0.2078719368264908D+00,-.5594649297339551D+00,
     2  0.5244946737686191D+01,0.1586467957461671D+01,
     3  -.1804530588786986D+01,-.2151581542538152D+01,
     4  0.1926146609482027D+00,0.1743086924136421D+01,
     5  0.3336674659916305D+00,-.1379020487018884D+01,
     6  -.2746984884197629D+00,0.1203541453532129D+01,
     7  -.8802841753261702D-01,-.1007825852134271D+01,
     8  0.5799189475963334D+00,0.5426501951392433D+00,
     9  -.8913138386398881D+00,0.2003364926655157D+00,
     *  0.6378410103570734D+00,-.7639953816366310D+00,
     1  0.1555267956560215D+00,0.5470605672458669D+00,
     2  -.7304566642660586D+00,0.3016384115586381D+00,
     3  0.3395338578691357D+00,-.6877225060838049D+00,
     4  0.5140865482209004D+00,0.1677067167038701D-01,
     5  -.5140051244341442D+00,0.6470992535311515D+00,
     6  -.3508108623025670D+00,-.1623253057944903D+00,
     7  0.5596020268615715D+00,-.5999508443748551D+00,
     8  0.2711951021300365D+00,0.2173029259120175D+00,
     9  -.5690352590480619D+00,0.5767225453117961D+00,
     *  -.2397722775650518D+00,-.5148453858574672D+01,
     1  -.9675364152223717D+00,0.2269841311793331D+01,
     2  0.1690225787213589D+01,-.1019424677591384D+01,
     3  -.1588364471349712D+01,0.5756957871292173D+00,
     4  0.1310328676601441D+01,-.6287896141712853D+00,
     5  -.9556778482342932D+00,0.8858942761196590D+00,
     6  0.4060431933710348D+00,-.1010890382532522D+01,
     7  0.3224828209627383D+00,0.6565143622518501D+00,
     8  -.8389163766037358D+00,0.1646963145143220D+00,
     9  0.6062757068590467D+00,-.7716596499489095D+00,
     *  0.2674730248441088D+00,0.4188509886973790D+00,
     1  -.7353649630830229D+00,0.4842398967269206D+00,
     2  0.1012649949874701D+00,-.5872056923821412D+00,
     3  0.6559544688105250D+00,-.2930540973469822D+00,
     4  -.2429939498134164D+00,0.6099347352008233D+00,
     5  -.5941315461670513D+00,0.2228931187699098D+00,
     6  0.2715603106726231D+00,-.5986467026469435D+00,
     7  0.5757634055402846D+00,-.2239621029112691D+00,
     8  -.2522281087314475D+00,0.5830451257574513D+00,
     9  0.5033145322114767D+01,0.3732226803493336D+00,
     *  -.2543021748101097D+01,-.1051726675961941D+01,
     1  0.1615201735949579D+01,0.1038902984176996D+01,
     2  -.1268948951501257D+01,-.7219407224606583D+00,
     3  0.1220298619640211D+01,0.1899318860514999D+00,
     4  -.1132836355919105D+01,0.4738301357917219D+00,
     5  0.6793183098748046D+00,-.9275631079057153D+00,
     6  0.1624464848567670D+00,0.6880760107021496D+00,
     7  -.8115663736850078D+00,0.2026015616158504D+00,
     8  0.5280041238389543D+00,-.7804195564143127D+00,
     9  0.4180011651095032D+00,0.2275132549807120D+00,
     *  -.6747824096530065D+00,0.6402198158749146D+00/
c 
c        following is block number 11 of length 180 elements
c 
        data www11/
     1  -.1878631739345199D+00,-.3632992804191695D+00,
     2  0.6675519030652408D+00,-.5589166518134097D+00,
     3  0.1252610608736164D+00,0.3678247552289286D+00,
     4  -.6400437177545594D+00,0.5494256263789092D+00,
     5  -.1580557764901868D+00,-.3135068068838113D+00,
     6  0.6100584171492838D+00,-.5761266001848216D+00,
     7  0.2335515598132133D+00,-.4902067321622479D+01,
     8  0.1873439365138454D+00,0.2635959412533801D+01,
     9  0.3302875298186422D+00,-.1904886329353431D+01,
     *  -.2642468633790321D+00,0.1552306951500251D+01,
     1  -.1253082107583728D+00,-.1252847519166199D+01,
     2  0.6658995510496086D+00,0.7066438466244374D+00,
     3  -.1037215937562952D+01,0.1488302876719776D+00,
     4  0.7952978728584508D+00,-.8465468656462475D+00,
     5  0.1008717574928904D+00,0.6643525523888697D+00,
     6  -.8095452621543996D+00,0.3028787612829415D+00,
     7  0.3939922342044226D+00,-.7596766666815076D+00,
     8  0.5791650768878351D+00,-.2621113380496269D-01,
     9  -.5136103930336851D+00,0.7102683218576795D+00,
     *  -.4737445212892812D+00,-.2827368243807729D-01,
     1  0.4957768005962510D+00,-.6734233240367458D+00,
     2  0.4805140504594467D+00,-.3635767481989899D-01,
     3  -.4171325520382231D+00,0.6459052839875501D+00,
     4  -.5396280037457474D+00,0.1596815529334961D+00,
     5  0.2980195472101586D+00,-.6028864425683032D+00,
     6  0.4757423154577068D+01,-.7071462774250689D+00,
     7  -.2566632387366494D+01,0.3894355896155893D+00,
     8  0.1873523166568159D+01,-.5419943566513719D+00,
     9  -.1375502039848850D+01,0.9082892117089581D+00,
     *  0.7448178563591810D+00,-.1177289567420743D+01,
     1  0.1204093172827958D+00,0.9332572859081395D+00,
     2  -.8708000076437313D+00,-.4580241701813980D-01,
     3  0.8220376090220509D+00,-.8039902526592157D+00,
     4  0.1261121810046919D+00,0.5916124380272208D+00,
     5  -.8148546399266269D+00,0.4482278856733146D+00,
     6  0.1945277701016943D+00,-.6712633477326484D+00,
     7  0.7045262276548390D+00,-.3153821554349103D+00,
     8  -.2367350063023551D+00,0.6315569415098342D+00,
     9  -.6663881422691869D+00,0.3460210224290730D+00,
     *  0.1435431676998054D+00,-.5458149241960187D+00,
     1  0.6647575438212744D+00,-.4533683070293716D+00,
     2  0.2509989783401604D-01,0.4088644629647716D+00,
     3  -.6421909903807955D+00,0.5679788138137878D+00,
     4  -.2246428889236793D+00,-.4600926702874416D+01,
     5  0.1180785409936567D+01,0.2357228813854954D+01,
     6  -.1036937694350650D+01,-.1555734284022177D+01,
     7  0.1205695118806253D+01,0.8207663791410504D+00,
     8  -.1360308564792966D+01,0.6376939372803637D-01,
     9  0.1112161100843211D+01,-.8739788693181969D+00,
     *  -.2488336503815124D+00,0.9905091731935248D+00,
     1  -.7390488850378399D+00,-.1219641831766721D+00,
     2  0.7990090023318421D+00,-.8027540732126996D+00,
     3  0.2250003782341705D+00,0.4615762417247108D+00,
     4  -.7951897469074658D+00,0.6081115260195326D+00,
     5  -.6769489101214745D-01,-.4804556490244717D+00,
     6  0.7317502401090305D+00,-.5764253001014310D+00,
     7  0.1266194285039684D+00,0.3694346433915527D+00,
     8  -.6640944700970472D+00,0.6268129948195552D+00,
     9  -.2922574797435663D+00,-.1706104941228481D+00,
     *  0.5447640081447942D+00,-.6639460843401327D+00,
     1  0.4818948690398120D+00,-.8621994709885659D-01,
     2  -.3446037004351557D+00,0.6202999466574006D+00,
     3  0.4433976017166590D+01,-.1604190914152862D+01,
     4  -.2032627778707010D+01,0.1558511865879763D+01,
     5  0.1021794138344792D+01,-.1602548794927620D+01,
     6  -.6101136336173155D-01,0.1350123304371822D+01,
     7  -.8337952582410369D+00,-.5260522794807554D+00,
     8  0.1151207118538089D+01,-.5827413256382536D+00,
     9  -.4431203021899939D+00,0.9774236065804665D+00,
     *  -.6778175719687742D+00,-.1001206389580730D+00,
     1  0.7360557531372195D+00,-.8244683513866168D+00,
     2  0.3805899762875923D+00,0.2620010107582771D+00,
     3  -.7088595913416190D+00,0.7326833984190006D+00,
     4  -.3621850226496059D+00,-.1761609201166668D+00,
     5  0.6001470469342608D+00,-.7135502747751819D+00,
     6  0.4846930771741500D+00,-.4024686306402485D-01,
     7  -.4080141728433026D+00,0.6625843466063956D+00,
     8  -.6219922635906638D+00,0.3150639212305571D+00,
     9  0.1200344893875997D+00,-.4982581781828246D+00,
     *  0.6642635972837399D+00,-.5531788680641278D+00,
     1  0.2136704004259053D+00,-.4257752005205002D+01,
     2  0.1974428410509214D+01,0.1619127392409469D+01,
     3  -.1918287929115848D+01,-.3625496231082449D+00,
     4  0.1673456438582815D+01,-.6961636982319397D+00,
     5  -.9059010314845774D+00,0.1268289663195542D+01,
     6  -.2927441363790122D+00,-.8235743246310451D+00,
     7  0.1065294469301931D+01,-.3941670633142216D+00,
     8  -.5079922263146237D+00,0.9439047383651727D+00,
     9  -.6864678095393837D+00,0.5041930847990198D-02,
     *  0.6219475378414343D+00,-.8342622812831679D+00/
c 
c        following is block number 12 of length 180 elements
c 
        data www12/
     1  0.5633825429426618D+00,-.8419680168496886D-02,
     2  -.5164073399580941D+00,0.7550999553401466D+00,
     3  -.6196462172969029D+00,0.2033671361908953D+00,
     4  0.2862182424896045D+00,-.6308176082541713D+00,
     5  0.6930370699670428D+00,-.4629232990519493D+00,
     6  0.4906904301952753D-01,0.3743338278015838D+00,
     7  -.6394772409868817D+00,0.6479702510233816D+00,
     8  -.4035402447778481D+00,0.5590859295314609D-02,
     9  0.3912956335795963D+00,-.6358973263814614D+00,
     *  0.4073279754195837D+01,-.2289560072581684D+01,
     1  -.1143365089025919D+01,0.2097961338021169D+01,
     2  -.3250450468285903D+00,-.1425846967018422D+01,
     3  0.1263096428852039D+01,0.1877080600018870D+00,
     4  -.1218351967922733D+01,0.9722510931759791D+00,
     5  0.7838420167246698D-01,-.9299018313782093D+00,
     6  0.9747818963357766D+00,-.3215257470154691D+00,
     7  -.4787437971080099D+00,0.8991522699328565D+00,
     8  -.7446258838185683D+00,0.1810658840665644D+00,
     9  0.4362973352275288D+00,-.7844249854005845D+00,
     *  0.7224316498815129D+00,-.3220111526801283D+00,
     1  -.2026319973653056D+00,0.6105465811262853D+00,
     2  -.7381252652706928D+00,0.5537690911240925D+00,
     3  -.1527887289503569D+00,-.2942483281847913D+00,
     4  0.6132858789802472D+00,-.6906027393728357D+00,
     5  0.5081847529433683D+00,-.1425399830257640D+00,
     6  -.2684533832767632D+00,0.5767012816834038D+00,
     7  -.6750632796151515D+00,0.5322067765360629D+00,
     8  -.2010095309591970D+00,-.3881468853644003D+01,
     9  0.2548536283599024D+01,0.6314010805519841D+00,
     *  -.2095484215594645D+01,0.9504730647937030D+00,
     1  0.9236159159599852D+00,-.1513596516126979D+01,
     2  0.5745287190796243D+00,0.7223193017263534D+00,
     3  -.1223729898832305D+01,0.6995764759432501D+00,
     4  0.2697588840259852D+00,-.9378106297020039D+00,
     5  0.9256072639731221D+00,-.3513279418042344D+00,
     6  -.3695489511947173D+00,0.8247071148140596D+00,
     7  -.8156575728868690D+00,0.4056351818620311D+00,
     8  0.1659058741127279D+00,-.6247772305904297D+00,
     9  0.7847588983397223D+00,-.6077895484668047D+00,
     *  0.1944346694787679D+00,0.2754879066932548D+00,
     1  -.6202387978166873D+00,0.7213643274163920D+00,
     2  -.5574178885647743D+00,0.2007251983957383D+00,
     3  0.2159714537419334D+00,-.5483812708801037D+00,
     4  0.6882107945800941D+00,-.5953004654173582D+00,
     5  0.3063716031196576D+00,0.8019382588854498D-01,
     6  -.4374140518566665D+00,0.6500152786378870D+00,
     7  0.3683140992754788D+01,-.2751106438587011D+01,
     8  -.1079484618185755D+00,0.1922973008055200D+01,
     9  -.1439417995809242D+01,-.2693970269195313D+00,
     *  0.1404805788486647D+01,-.1154677563902556D+01,
     1  0.2517329938570768D-01,0.9555919056837775D+00,
     2  -.1129640629758952D+01,0.5358734890464424D+00,
     3  0.3197093127424656D+00,-.8942678357828589D+00,
     4  0.9177535311267683D+00,-.4616373711063394D+00,
     5  -.1815585571573859D+00,0.6865513605172679D+00,
     6  -.8460060752630993D+00,0.6325259747929235D+00,
     7  -.1737345284302524D+00,-.3257632672375817D+00,
     8  0.6734203291141050D+00,-.7553700869455215D+00,
     9  0.5637625248566405D+00,-.1845751462591582D+00,
     *  -.2423971186728088D+00,0.5741973345622290D+00,
     1  -.7093715345379232D+00,0.6148013089312815D+00,
     2  -.3296906764509375D+00,-.5151244528912973D-01,
     3  0.4102740661436496D+00,-.6394857006220044D+00,
     4  0.6735359574705580D+00,-.5053900393829489D+00,
     5  0.1869170653315365D+00,-.3479049433839138D+01,
     6  0.2897741858469804D+01,-.4042603495046755D+00,
     7  -.1604071559888224D+01,0.1740293621209632D+01,
     8  -.4167966381199335D+00,-.9780350818719702D+00,
     9  0.1393875648837193D+01,-.7539694205970104D+00,
     *  -.2953920332380405D+00,0.1015941653784270D+01,
     1  -.1045340971458526D+01,0.4906970709990586D+00,
     2  0.2599375220181004D+00,-.8056451988943580D+00,
     3  0.9229306296809998D+00,-.6203590023191889D+00,
     4  0.8057621525321797D-01,0.4507084825172242D+00,
     5  -.7707334004537515D+00,0.7844549607211647D+00,
     6  -.5175988548585850D+00,0.8728709309261024D-01,
     7  0.3502622238273396D+00,-.6541346128957557D+00,
     8  0.7391418247806055D+00,-.5937909420336242D+00,
     9  0.2750456507361483D+00,0.1148543960139574D+00,
     *  -.4608161412705561D+00,0.6671113201122851D+00,
     1  -.6814378699878402D+00,0.5058488100971587D+00,
     2  -.1931729011769113D+00,-.1690260358319377D+00,
     3  0.4822914367929913D+00,-.6628558851663823D+00,
     4  0.3269893089441535D+01,-.2989566516668425D+01,
     5  0.8849782719695163D+00,0.1170991107741274D+01,
     6  -.1827451057537405D+01,0.1018784239333951D+01,
     7  0.3416316444441060D+00,-.1238364470494732D+01,
     8  0.1215829553335102D+01,-.4776450561068184D+00,
     9  -.4264005957527706D+00,0.9933373841902246D+00,
     *  -.1003728749465612D+01,0.5455852326656521D+00/
c 
c        following is block number 13 of length 180 elements
c 
        data www13/
     1  0.1064536302008601D+00,-.6528330080335500D+00,
     2  0.8927420770041758D+00,-.7760314573295088D+00,
     3  0.3880340878828887D+00,0.1071386408031718D+00,
     4  -.5348198798350502D+00,0.7657996541265381D+00,
     5  -.7471419036958569D+00,0.5053904028441462D+00,
     6  -.1275466708962384D+00,-.2706394777524458D+00,
     7  0.5787364923482757D+00,-.7194295030974725D+00,
     8  0.6645614585243132D+00,-.4375035616280473D+00,
     9  0.1031777847309019D+00,0.2505625706593139D+00,
     *  -.5352089032057834D+00,0.6823924104360428D+00,
     1  -.6588835037440458D+00,0.4730073388481372D+00,
     2  -.1715834330478730D+00,-.3056326898204173D+01,
     3  0.3028292880757042D+01,-.1316825667347118D+01,
     4  -.6614283406372837D+00,0.1701275088837427D+01,
     5  -.1442191214622543D+01,0.3580860467931627D+00,
     6  0.7452191968374065D+00,-.1264963622850387D+01,
     7  0.1053273342317959D+01,-.3551119845801965D+00,
     8  -.4196621916125094D+00,0.9222275532006212D+00,
     9  -.9917957026231073D+00,0.6686134341433382D+00,
     *  -.1286751695865894D+00,-.4098092256207956D+00,
     1  0.7677937453485452D+00,-.8530301151015835D+00,
     2  0.6700551947853393D+00,-.3004857477196885D+00,
     3  -.1331269156698656D+00,0.5065240599626354D+00,
     4  -.7257608852531037D+00,0.7461211529830819D+00,
     5  -.5763552724398393D+00,0.2701195220639930D+00,
     6  0.9130903204119998D-01,-.4200518795319643D+00,
     7  0.6415760950651518D+00,-.7097503968656178D+00,
     8  0.6146528263269551D+00,-.3823912219882175D+00,
     9  0.6770534029756587D-01,0.2586861541640525D+00,
     *  -.5252787281500638D+00,0.6745468073971779D+00,
     1  0.2838969593055321D+01,-.3016161156477520D+01,
     2  0.1685556155481861D+01,0.1155390369999321D+00,
     3  -.1385614811874917D+01,0.1627047955564231D+01,
     4  -.9704535822636003D+00,-.5851887404152763D-01,
     5  0.8980249810297985D+00,-.1212416012608392D+01,
     6  0.9696969618574054D+00,-.3678206195853837D+00,
     7  -.3045324436715156D+00,0.7959333781496548D+00,
     8  -.9677385641476703D+00,0.8106198470691771D+00,
     9  -.4169946736198138D+00,-.6826344734698281D-01,
     *  0.4985985073442328D+00,-.7644900454065994D+00,
     1  0.8140361079304543D+00,-.6554216079172465D+00,
     2  0.3452367031580245D+00,0.3157006018658195D-01,
     3  -.3841843486685896D+00,0.6359313354263054D+00,
     4  -.7384359584920258D+00,0.6783156349104264D+00,
     5  -.4763721753431980D+00,0.1805031069134837D+00,
     6  0.1456221967658919D+00,-.4356962477087487D+00,
     7  0.6333194020682686D+00,-.7019217205861497D+00,
     8  0.6306241697915991D+00,-.4353368803349220D+00,
     9  0.1551581654677127D+00,-.2618409589704171D+01,
     *  0.2955880850787279D+01,-.1980212854066665D+01,
     1  0.4268817476593159D+00,0.9231546245548432D+00,
     2  -.1553932455548067D+01,0.1371502802282013D+01,
     3  -.6345493781722385D+00,-.2454882180900756D+00,
     4  0.9049340645857451D+00,-.1147000485593757D+01,
     5  0.9633112364692904D+00,-.4854042023879771D+00,
     6  -.9270237932226257D-01,0.5869320756737445D+00,
     7  -.8715813786746013D+00,0.8988065517281133D+00,
     8  -.6935524841822630D+00,0.3326372532743379D+00,
     9  0.8274120143540498D-01,-.4534693171010256D+00,
     *  0.7029179925858820D+00,-.7888571046606584D+00,
     1  0.7066383653045475D+00,-.4849427399034482D+00,
     2  0.1762832291174292D+00,0.1552158315153777D+00,
     3  -.4463492388337561D+00,0.6456670013374583D+00,
     4  -.7211972919966263D+00,0.6642909626630698D+00,
     5  -.4894576685290473D+00,0.2306168735203062D+00,
     6  0.6538531416128287D-01,-.3468651870757191D+00,
     7  0.5657495785822275D+00,-.6851712639478468D+00,
     8  0.2395209493189977D+01,-.2850573986000615D+01,
     9  0.2193182495187298D+01,-.9288906664307295D+00,
     *  -.3694208245645470D+00,0.1243597237785603D+01,
     1  -.1486753728477276D+01,0.1154520199369350D+01,
     2  -.4779294111731727D+00,-.2584768899157736D+00,
     3  0.8187204488874777D+00,-.1071129316092345D+01,
     4  0.9971289968795518D+00,-.6670831747820946D+00,
     5  0.2009129822734307D+00,0.2717539769818922D+00,
     6  -.6426674431784648D+00,0.8433882180549319D+00,
     7  -.8511679385405859D+00,0.6847391024786847D+00,
     8  -.3934859672815645D+00,0.4351540631620356D-01,
     9  0.2963766613575842D+00,-.5666407351459049D+00,
     *  0.7254624283646328D+00,-.7532411090005811D+00,
     1  0.6531421912803267D+00,-.4482862953093341D+00,
     2  0.1764477038517791D+00,0.1167134538747974D+00,
     3  -.3848440624076976D+00,0.5875667102507273D+00,
     4  -.6958479503925024D+00,0.6955162360222639D+00,
     5  -.5886334081309268D+00,0.3926811014870490D+00,
     6  -.1377638213149680D+00,-.2169909570833715D+01,
     7  0.2703719578492528D+01,-.2320155714428504D+01,
     8  0.1358465899295843D+01,-.2138588004486023D+00,
     9  -.7508927262267700D+00,0.1302351923177802D+01,
     *  -.1374242232938581D+01,0.1044778097390988D+01/
c 
c        following is block number 14 of length 180 elements
c 
        data www14/
     1  -.4786010390966556D+00,-.1377694178762489D+00,
     2  0.6468060850925682D+00,-.9487440199760855D+00,
     3  0.1008951799442107D+01,-.8494962582305758D+00,
     4  0.5318377758869662D+00,-.1367891738482826D+00,
     5  -.2538474932524369D+00,0.5713844504869592D+00,
     6  -.7690528156348810D+00,0.8255860957286221D+00,
     7  -.7443622755719464D+00,0.5492443677742561D+00,
     8  -.2783662615200634D+00,-.2295419986380634D-01,
     9  0.3091637313883942D+00,-.5404991731605462D+00,
     *  0.6875520982067541D+00,-.7340758268098667D+00,
     1  0.6779577038898554D+00,-.5304664464758835D+00,
     2  0.3140325897193398D+00,-.5892667222141530D-01,
     3  -.2007355696912026D+00,0.4312075099637693D+00,
     4  -.6031050199108406D+00,0.6947844838216287D+00,
     5  0.1943030440443406D+01,-.2519099189984527D+01,
     6  0.2360002651200895D+01,-.1690093634574746D+01,
     7  0.7653620181141052D+00,0.1544438323607497D+00,
     8  -.8636896339417402D+00,0.1247079723074330D+01,
     9  -.1284900602817450D+01,0.1034113207904722D+01,
     *  -.5975482091663620D+00,0.9228550337143166D-01,
     1  0.3756985501623511D+00,-.7276689095608679D+00,
     2  0.9191241851437379D+00,-.9393886087543535D+00,
     3  0.8061587559677216D+00,-.5572995465479350D+00,
     4  0.2417868773632755D+00,0.8881504386756762D-01,
     5  -.3874480034682742D+00,0.6167525496033913D+00,
     6  -.7522141322649609D+00,0.7833557707447949D+00,
     7  -.7132304176916944D+00,0.5565621004862034D+00,
     8  -.3369302825685766D+00,0.8339849183424478D-01,
     9  0.1730340763782875D+00,-.4028225264920639D+00,
     *  0.5808299349469220D+00,-.6886333520767181D+00,
     1  0.7160415852216838D+00,-.6617555016411382D+00,
     2  0.5331700167716153D+00,-.3453796613174453D+00,
     3  0.1195042515095535D+00,-.1715075154832002D+01,
     4  0.2300743495834933D+01,-.2314573766433793D+01,
     5  0.1905838449784636D+01,-.1230150943997388D+01,
     6  0.4560873794156834D+00,0.2630027321990933D+00,
     7  -.8135314797411937D+00,0.1133418621743543D+01,
     8  -.1211097831234884D+01,0.1075669220678813D+01,
     9  -.7828678382843283D+00,0.4006751886258035D+00,
     *  0.2980985483983102D-02,-.3692778292071218D+00,
     1  0.6540737588908623D+00,-.8303500661116594D+00,
     2  0.8881318725180174D+00,-.8325823998402929D+00,
     3  0.6809314059903732D+00,-.4588095489634584D+00,
     4  0.1964507162200668D+00,0.7488628004244370D-01,
     5  -.3260270479000151D+00,0.5323614307115028D+00,
     6  -.6756085669660320D+00,0.7448441691372308D+00,
     7  -.7368171190546012D+00,0.6556216919365309D+00,
     8  -.5118211084657555D+00,0.3211386205248276D+00,
     9  -.1028447872724186D+00,-.1220252623155234D+00,
     *  0.3324922532256477D+00,-.5093659253710130D+00,
     1  0.6367785479342465D+00,-.7034233944481632D+00,
     2  0.1486530816571395D+01,-.2052879911042778D+01,
     3  0.2188436425164817D+01,-.1995892522898998D+01,
     4  0.1564549246245671D+01,-.9932183992541947D+00,
     5  0.3807660856565230D+00,0.1866150236357769D+00,
     6  -.6445021026135712D+00,0.9541464129813289D+00,
     7  -.1101816801656939D+01,0.1095025924133427D+01,
     8  -.9571356586670679D+00,0.7215399469927417D+00,
     9  -.4262336965362963D+00,0.1092281440651648D+00,
     *  0.1949875627988181D+00,-.4579117161747128D+00,
     1  0.6584985622668275D+00,-.7837603131583302D+00,
     2  0.8286556101499914D+00,-.7954211662000594D+00,
     3  0.6924965360021341D+00,-.5331796954224246D+00,
     4  0.3341362594685924D+00,-.1138692923115681D+00,
     5  -.1087594458212622D+00,0.3158819326643375D+00,
     6  -.4917907482154463D+00,0.6238763070859672D+00,
     7  -.7033050400282626D+00,0.7254230853490204D+00,
     8  -.6898832236744045D+00,0.6005049913631724D+00,
     9  -.4648889619212472D+00,0.2938159117058777D+00,
     *  -.1004698145370018D+00,-.1257869823533538D+01,
     1  0.1779881382479340D+01,-.1988558272583060D+01,
     2  0.1958623743433523D+01,-.1739588789163545D+01,
     3  0.1383555578864706D+01,-.9455281331069459D+00,
     4  0.4788089750764386D+00,-.3036401958868544D-01,
     5  -.3623929263793325D+00,0.6731418150659164D+00,
     6  -.8867645482955178D+00,0.9984213734509454D+00,
     7  -.1011992277087438D+01,0.9382066424559074D+00,
     8  -.7926961776343464D+00,0.5941279625832955D+00,
     9  -.3625142357772393D+00,0.1177527860283742D+00,
     *  0.1215769801491735D+00,-.3391559997138879D+00,
     1  0.5216113105341785D+00,-.6589591540498923D+00,
     2  0.7448241032823253D+00,-.7764564004403107D+00,
     3  0.7545721638585554D+00,-.6830428926160691D+00,
     4  0.5684619110608112D+00,-.4196161098675868D+00,
     5  0.2468915608241012D+00,-.6164129492744682D-01,
     6  -.1244572788913709D+00,0.3000291823600758D+00,
     7  -.4545886448128301D+00,0.5790639694604622D+00,
     8  -.6662414595399823D+00,0.7111126416760064D+00,
     9  0.1029550821751628D+01,-.1486216502883956D+01,
     *  0.1723948730694752D+01,-.1800164662318318D+01/
c 
c        following is block number 15 of length 180 elements
c 
        data www15/
     1  0.1742850697063349D+01,-.1576957593166293D+01,
     2  0.1328553022313204D+01,-.1024445816400476D+01,
     3  0.6907715383386162D+00,-.3515833803078953D+00,
     4  0.2778395018567716D-01,0.2635371536530407D+00,
     5  -.5094027404459992D+00,0.7009611136246944D+00,
     6  -.8332901943610382D+00,0.9050579568148933D+00,
     7  -.9180897924130920D+00,0.8768804657630591D+00,
     8  -.7880783316678217D+00,0.6599623376010509D+00,
     9  -.5019274407629226D+00,0.3239907611536930D+00,
     *  -.1363285105125777D+00,-.5114793822219861D-01,
     1  0.2291698338937030D+00,-.3894369005244720D+00,
     2  0.5248866803572404D+00,-.6299040189292290D+00,
     3  0.7004686295459566D+00,-.7342400078294323D+00,
     4  0.7305803304961963D+00,-.6905173393123865D+00,
     5  0.6166505498591056D+00,-.5130054022506439D+00,
     6  0.3848411546301224D+00,-.2384193794137912D+00,
     7  0.8074082924260898D-01,-.8020194245897756D+00,
     8  0.1176401133618858D+01,-.1405270046466130D+01,
     9  0.1533602263669746D+01,-.1578676997871027D+01,
     *  0.1551998872888544D+01,-.1464039049676816D+01,
     1  0.1325393776286283D+01,-.1146896338312140D+01,
     2  0.9394027489621884D+00,-.7135020579738462D+00,
     3  0.4792429011151923D+00,-.2459057980281931D+00,
     4  0.2182600530847167D-01,0.1857370559668586D+00,
     5  -.3706898164049277D+00,0.5281595103819739D+00,
     6  -.6545183510729532D+00,0.7473801016477688D+00,
     7  -.8055714393712056D+00,0.8290800878531437D+00,
     8  -.8189815598970759D+00,0.7773464242862997D+00,
     9  -.7071302067488548D+00,0.6120482957080132D+00,
     *  -.4964385000522910D+00,0.3651141648118001D+00,
     1  -.2232109672344457D+00,0.7603067397472501D-01,
     2  0.7111477054321963D-01,-.2130554400285739D+00,
     3  0.3449085847180837D+00,-.4622135701077351D+00,
     4  0.5610536274785586D+00,-.6381617303890680D+00,
     5  0.6910082564994172D+00,-.7178684936465191D+00,
     6  0.5757087442188648D+00,-.8549517426919102D+00,
     7  0.1044429226330316D+01,-.1177842374802055D+01,
     8  0.1266866577678294D+01,-.1317327410837189D+01,
     9  0.1333062709277839D+01,-.1317217548805642D+01,
     *  0.1272753191667120D+01,-.1202659206542256D+01,
     1  0.1110037659077424D+01,-.9981264206718233D+00,
     2  0.8702907230154418D+00,-.7299965437088290D+00,
     3  0.5807726989446132D+00,-.4261655577827685D+00,
     4  0.2696889770693858D+00,-.1147714662141038D+00,
     5  -.3529768727260483D-01,0.1774217547112319D+00,
     6  -.3087479025212987D+00,0.4267151740121776D+00,
     7  -.5290978824173337D+00,0.6140433154301013D+00,
     8  -.6801028014750890D+00,0.7262553518400833D+00,
     9  -.7519232681549939D+00,0.7569792886182506D+00,
     *  -.7417450356424329D+00,0.7069807188541908D+00,
     1  -.6538662371771776D+00,0.5839740085626479D+00,
     2  -.4992340323664271D+00,0.4018918540753846D+00,
     3  -.2944602519222616D+00,0.1796655970139540D+00,
     4  -.6038995037742634D-01,-.3510397716565101D+00,
     5  0.5263406746993044D+00,-.6541618981366773D+00,
     6  0.7562319452598147D+00,-.8400843865735621D+00,
     7  0.9090539927848734D+00,-.9647411602326444D+00,
     8  0.1007905446810812D+01,-.1038873906158425D+01,
     9  0.1057759943365566D+01,-.1064594415423373D+01,
     *  0.1059410107239815D+01,-.1042298449701036D+01,
     1  0.1013448088326349D+01,-.9731706437013304D+00,
     2  0.9219168563202508D+00,-.8602851456576862D+00,
     3  0.7890239498228440D+00,-.7090288215458890D+00,
     4  0.6213350241587654D+00,-.5271062352789233D+00,
     5  0.4276198896153743D+00,-.3242496531308033D+00,
     6  0.2184455041277785D+00,-.1117118929460496D+00,
     7  0.5584454213916678D-02,0.9839425061273443D-01,
     8  -.1986994805653064D+00,0.2938490742211013D+00,
     9  -.3824277494973146D+00,0.4631104865580793D+00,
     *  -.5346845563836672D+00,0.5960697846638253D+00,
     1  -.6463366746036141D+00,0.6847220526050230D+00,
     2  -.7106419470468314D+00,0.7237014617979318D+00,
     3  0.1284216342422759D+00,-.1949535719866160D+00,
     4  0.2476186573737516D+00,-.2950273903776195D+00,
     5  0.3402993633893572D+00,-.3847848998290912D+00,
     6  0.4290135544690596D+00,-.4730594236909377D+00,
     7  0.5167215514760596D+00,-.5596278980609735D+00,
     8  0.6013010032274149D+00,-.6412018456876520D+00,
     9  0.6787603856689003D+00,-.7133977104486145D+00,
     *  0.7445428162678904D+00,-.7716459291410001D+00,
     1  0.7941895425945407D+00,-.8116978760559119D+00,
     2  0.8237451483639688D+00,-.8299628630567932D+00,
     3  0.8300461794866040D+00,-.8237593715539745D+00,
     4  0.8109403369757258D+00,-.7915041030706186D+00,
     5  0.7654452724471219D+00,-.7328393586646012D+00,
     6  0.6938429745701135D+00,-.6486928523482110D+00,
     7  0.5977036928278812D+00,-.5412648611842929D+00,
     8  0.4798359660524061D+00,-.4139413786149240D+00,
     9  0.3441637669376922D+00,-.2711367382829591D+00,
     *  0.1955366979782599D+00,-.1180740473506653D+00/
c 
c        following is block number 16 of length 112 elements
c 
        data www16/
     1  0.3948385499392695D-01,0.9174824624319424D-01,
     2  -.1349508371215745D+00,0.1620362079930221D+00,
     3  -.1782021517764676D+00,0.1853793086188237D+00,
     4  -.1847301693877507D+00,0.1771861969258051D+00,
     5  -.1635845923775079D+00,0.1446962501634512D+00,
     6  -.1212256720213820D+00,0.9380829025029251D-01,
     7  -.6301210985459161D-01,0.2934397731005927D-01,
     8  0.6741201814036769D-02,-.4483014019749116D-01,
     9  0.8454114200064233D-01,-.1255152373318130D+00,
     *  0.1674089035720270D+00,-.2098884957382540D+00,
     1  0.2526263236347159D+00,-.2952982156309132D+00,
     2  0.3375823671676245D+00,-.3791592635405441D+00,
     3  0.4197124756221066D+00,-.4589301443196658D+00,
     4  0.4965069892061485D+00,-.5321466960181988D+00,
     5  0.5655645553138636D+00,-.5964902401349509D+00,
     6  0.6246706241123417D+00,-.6498725533781006D+00,
     7  0.6718854962621738D+00,-.6905240043977252D+00,
     8  0.7056299278345789D+00,-.7170743353022603D+00,
     9  0.7247590990507156D+00,-.7286181118522980D+00,
     *  -.3090840759906456D+00,0.4592783756378569D+00,
     1  -.5620045441754276D+00,0.6361511654604086D+00,
     2  -.6890773220047616D+00,0.7251488725410975D+00,
     3  -.7476294728422488D+00,0.7592088636425795D+00,
     4  -.7621541634944084D+00,0.7583561200560757D+00,
     5  -.7493560736522343D+00,0.7363778742969507D+00,
     6  -.7203674783716560D+00,0.7020365493127541D+00,
     7  -.6819057534875589D+00,0.6603446515443267D+00,
     8  -.6376065286187083D+00,0.6138576331808867D+00,
     9  -.5892010101924205D+00,0.5636954915851083D+00,
     *  -.5373705531627361D+00,0.5102377549498618D+00,
     1  -.4822994183343388D+00,0.4535551003545334D+00,
     2  -.4240063274586804D+00,0.3936599604284078D+00,
     3  -.3625304841300833D+00,0.3306414515323948D+00,
     4  -.2980262601689175D+00,0.2647283992257854D+00,
     5  -.2308012747647136D+00,0.1963076974108941D+00,
     6  -.1613190995348386D+00,0.1259145361986813D+00,
     7  -.9017951484418188D-01,0.5420469202750022D-01,
     8  -.1808447080797325D-01,0.5232116767081041D+00,
     9  -.7741342860149014D+00,0.9402046065150807D+00,
     *  -.1053255171435532D+01,0.1126362852468384D+01,
     1  -.1168041789455170D+01,0.1185186400899267D+01,
     2  -.1183746567601229D+01,0.1168829890830902D+01,
     3  -.1144694992752059D+01,0.1114768791977114D+01,
     4  -.1081710432684130D+01,0.1047509074659856D+01,
     5  -.1013596190599687D+01,0.9809567732498248D+00,
     6  -.9502298789867737D+00,0.9217941085366299D+00,
     7  -.8958371457319190D+00,0.8724104869270404D+00,
     8  -.8514714182284320D+00,0.8329145452434759D+00,
     9  -.8165950591430609D+00,0.8023456352113559D+00,
     *  -.7899885241368606D+00,0.7793440741436697D+00,
     1  -.7702366411798358D+00,0.7624986128556170D+00,
     2  -.7559730878718678D+00,0.7505156102357613D+00,
     3  -.7459952493320651D+00,0.7422952358803207D+00,
     4  -.7393133038459493D+00,0.7369618444524820D+00,
     5  -.7351679465609392D+00,0.7338733747298078D+00,
     6  -.7330345198590142D+00,0.7326223456343998D+00/
  
c 
c       return the varous types of data to the user
c 
        do 1200 i=1,nn/2
c 
        tt(i)=tt0(i)
        tt(nn-i+1)=-tt0(i)
c 
        ww(i)=ww0(i)
        ww(nn-i+1)=ww0(i)
 1200 continue
c 
        do 1400 i=1,74
c 
        cm(i)=cm0(i)
c 
        rlams(i)=rlams0(i)*ima**(i-1)
 1400 continue
c 
        do 1500 i=1,74
c 
        rlams(i)=rlams0(i)*ima**(i-1)
 1500 continue
c 
        d=1
        do 1800 i=1,74
c 
        do 1600 j=1,37
c 
        www(j,i)=www101(j,i)
        www(74-j+1,i)=www101(j,i)*d
c 
 1600 continue
        d=-d
 1800 continue
c 
        return
c 
c 
c 
c 
        entry expo_ttww_ret63(tt,ww,nnout,c)
c 
c       This entry returns to the user the 74 prolate nodes tt and
c       their corresponding weights ww, discretizing to 16-digit
c       (or so) accuracy band-limited functions with c=63 on the
c       interval [-1,1]. It also returns nnout=74, and c=63
c 
c       return the nodes and weights to the user
c 
        do 2200 i=1,nn/2
c 
        tt(i)=tt0(i)
        tt(nn-i+1)=-tt0(i)
c 
        ww(i)=ww0(i)
        ww(nn-i+1)=ww0(i)
 2200 continue
c 
        c=63
c 
        nnout=nn
        return
        end
  
  
  
c 
c 
c 
c 
c 
        subroutine exporet20(tt,ww,cm,rlams,www)
        implicit real *8 (a-h,o-z)
        save
        dimension tt0(19),ww0(19),cm0(39),rlams0(42),
     1      tt(1),ww(1),cm(1),www1(180),www2(180),www3(180),
     2      www4(180),www5(40),www(38,40),www0(2270),
     3      www10(19,40)
c 
        equivalence (www0(1),www1(1)),(www0(181),www2(1)),
     1       (www0(361),www3(1)),(www0(541),www4(1)),
     2       (www0(721),www5(1))
c 
        equivalence(www0(1),www10(1,1))
  
        complex *16 rlams(1),ima
c 
        data nn/38/,ima/(0.0d0,1.0d0)/
c 
        data tt0/
     1  -.9977099855266336D+00,-.9879814000447955D+00,
     2  -.9706678464542580D+00,-.9460728575853053D+00,
     3  -.9146071418462681D+00,-.8767569865116678D+00,
     4  -.8330553189029443D+00,-.7840562169927318D+00,
     5  -.7303145588934070D+00,-.6723713303736169D+00,
     6  -.6107441965326103D+00,-.5459224550598787D+00,
     7  -.4783653442533510D+00,-.4085027481873641D+00,
     8  -.3367375104689089D+00,-.2634487611861341D+00,
     9  -.1889958382343864D+00,-.1137225277646114D+00,
     *  -.3796145638706746D-01/
c 
        data ww0/
     1  0.5871473778196759D-02,0.1356088858562982D-01,
     2  0.2101491843545797D-01,0.2810610167119679D-01,
     3  0.3474379339274391D-01,0.4086726705215778D-01,
     4  0.4644338157697391D-01,0.5146209880237784D-01,
     5  0.5593089820335008D-01,0.5986919973500332D-01,
     6  0.6330349244899238D-01,0.6626346441743553D-01,
     7  0.6877914982179766D-01,0.7087895814127994D-01,
     8  0.7258839494384392D-01,0.7392928503230201D-01,
     9  0.7491933659015427D-01,0.7557192088056520D-01,
     *  0.7589597649054098D-01/
c 
        data cm0/
     1  0.9999999999999994D+00,0.1000000000000020D+01,
     2  0.1000000000001491D+01,0.1000000000070044D+01,
     3  0.1000000002349488D+01,0.1000000059748815D+01,
     4  0.1000001192744885D+01,0.1000019072694090D+01,
     5  0.1000246613986270D+01,0.1002574095223061D+01,
     6  0.1021329771732463D+01,0.1136733579423376D+01,
     7  0.1698388652621668D+01,0.4367027579939473D+01,
     8  0.1990208326013668D+02,0.1347484830107972D+03,
     9  0.1163008854989965D+04,0.1194178983667014D+05,
     *  0.1416416507205055D+06,0.1909311808961025D+07,
     1  0.2892305630924936D+08,0.4881259789907936D+09,
     2  0.9113428127625433D+10,0.1871214793973914D+12,
     3  0.4203853159635022D+13,0.1028774513167125D+15,
     4  0.2731666633291523D+16,0.7842134208199461D+17,
     5  0.2426340384738884D+19,0.8067076142706473D+20,
     6  0.2874549152280671D+22,0.1095081091869021D+24,
     7  0.4449994797369014D+25,0.1924834707056343D+27,
     8  0.8844930933467107D+28,0.4309872993636294D+30,
     9  0.2223081416393892D+32,0.1211892195300464D+34,
     *  0.6971527926813384D+35/
c 
        data rlams0/
     1  0.5604991216397930D+00,0.5604991216397874D+00,
     2  0.5604991216393751D+00,0.5604991216201632D+00,
     3  0.5604991209813498D+00,0.5604991048952143D+00,
     4  0.5604987873738617D+00,0.5604937766021087D+00,
     5  0.5604300209590818D+00,0.5597791222963892D+00,
     6  0.5546154204401222D+00,0.5257092173750537D+00,
     7  0.4300870802206489D+00,0.2682143248247087D+00,
     8  0.1256393459693961D+00,0.4828508430715566D-01,
     9  0.1643551910912521D-01,0.5129088875375108D-02,
     *  0.1489290593986078D-02,0.4056360145345845D-03,
     1  0.1042204261898806D-03,0.2536932868557876D-04,
     2  0.5871297001984451D-05,0.1295725826243741D-05,
     3  0.2733701980937441D-06,0.5526050355923076D-07,
     4  0.1072410573798299D-07,0.2001510169839282D-08,
     5  0.3598314119335601D-09,0.6240463656638248D-10,
     6  0.1045418224565700D-10,0.1693759703123382D-11,
     7  0.2657022210276245D-12,0.4039970704669495D-13,
     8  0.5959745485545675D-14,0.8537735579324383D-15,
     9  0.1188768400342671D-15,0.1610063267752776D-16,
     *  0.2122809148651942D-17,0.2726552743763464D-18,
     1  0.3413873914438881D-19,0.4169575250013808D-20/
c 
        data www1/
     1  0.7573563373310281D-07,0.2527297727170619D-06,
     2  0.1063315163136208D-05,0.4668141405173276D-05,
     3  0.2014562499964947D-04,0.8294048041387841D-04,
     4  0.3196015929886324D-03,0.1137399186786685D-02,
     5  0.3702578059192846D-02,0.1094877755146564D-01,
     6  0.2926424948671362D-01,0.7045315235101572D-01,
     7  0.1524091805958191D+00,0.2957830408847592D+00,
     8  0.5144449673179284D+00,0.8013689319107636D+00,
     9  0.1117618558196425D+01,0.1395195923728090D+01,
     *  0.1558886430899670D+01,-.8923689579886855D-06,
     1  -.2718294812431669D-05,-.1042146553506025D-04,
     2  -.4169111317203468D-04,-.1637287630594050D-03,
     3  -.6120226522265641D-03,-.2134446352349163D-02,
     4  -.6845879674346563D-02,-.1997424794764854D-01,
     5  -.5256405071025874D-01,-.1238810581468614D+00,
     6  -.2598073307105063D+00,-.4817446806278943D+00,
     7  -.7837291921148170D+00,-.1106703012978779D+01,
     8  -.1332670202032805D+01,-.1321526589764515D+01,
     9  -.9868591423188869D+00,-.3669934345300187D+00,
     *  0.7277074490010725D-05,0.2019460720754123D-04,
     1  0.7033423974181475D-04,0.2555974311215893D-03,
     2  0.9103417705835859D-03,0.3077560906235868D-02,
     3  0.9668848432966011D-02,0.2778810807629992D-01,
     4  0.7213732998151748D-01,0.1673095583629131D+00,
     5  0.3430768600889518D+00,0.6148804284773143D+00,
     6  0.9490738810132043D+00,0.1233206439098972D+01,
     7  0.1292249285145860D+01,0.9795075384757406D+00,
     8  0.3080034789372892D+00,-.4827257382465756D+00,
     9  -.1019445737585801D+01,-.4734210511682510D-04,
     *  -.1194387042132895D-03,-.3765763049096955D-03,
     1  -.1238620111320059D-02,-.3984919876317772D-02,
     2  -.1212725247281002D-01,-.3412868843878981D-01,
     3  -.8726062274917107D-01,-.1996365708801532D+00,
     4  -.4027162584253929D+00,-.7046741621373096D+00,
     5  -.1046627787995198D+01,-.1273915238842010D+01,
     6  -.1180332209427835D+01,-.6527075985046640D+00,
     7  0.1637450495412428D+00,0.8534855384587676D+00,
     8  0.9761158613910047D+00,0.4254077988879231D+00,
     9  0.2600596205019238D-03,0.5951865711151749D-03,
     *  0.1691901940603707D-02,0.5015190051271380D-02,
     1  0.1450547791627262D-01,0.3951484189499765D-01,
     2  0.9890839715322876D-01,0.2228946901868231D+00,
     3  0.4436237973283728D+00,0.7635090312506666D+00,
     4  0.1105021856066826D+01,0.1283967597660786D+01,
     5  0.1078560552982060D+01,0.4179171476353463D+00,
     6  -.4365852355844033D+00,-.9545426461045412D+00,
     7  -.7309606919218861D+00,0.8558395260284081D-01,
     8  0.8141052083524086D+00,-.1242519321743034D-02,
     9  -.2574060394724514D-02,-.6565940775257672D-02,
     *  -.1744946137078246D-01,-.4510843287302408D-01,
     1  -.1092198251087665D+00,-.2409439494435535D+00,
     2  -.4725665825565007D+00,-.8029865960124514D+00,
     3  -.1143274735268904D+01,-.1290581104081166D+01,
     4  -.1014166653697318D+01,-.2761209611167692D+00,
     5  0.5779234581198878D+00,0.9494557500923620D+00,
     6  0.5017833436873305D+00,-.3917257706943882D+00,
     7  -.8734816828215742D+00,-.4476797569784074D+00,
     8  0.5252050826626848D-02,0.9827729316047955D-02,
     9  0.2237000061095420D-01,0.5295897827467638D-01,
     *  0.1214661372086720D+00,0.2590206726181624D+00,
     1  0.4973997897041962D+00,0.8336736829746914D+00,
     2  0.1173661669602186D+01,0.1304751410602755D+01,
     3  0.9897743619756058D+00,0.2083958766642108D+00,
     4  -.6418472622350806D+00,-.9205830715308203D+00,
     5  -.3471611774839628D+00,0.5450981359068335D+00,
     6  0.8128101823914210D+00,0.1459278605247662D+00,
     7  -.6821552438636614D+00,-.1982737761592436D-01,
     8  -.3344636195421829D-01,-.6748936789659341D-01,
     9  -.1411768354353321D+00,-.2845366740001647D+00,
     *  -.5277569180984867D+00,-.8665055767736538D+00,
     1  -.1205535828562192D+01,-.1329792474609819D+01,
     2  -.9981588424686301D+00,-.1940765962325170D+00,
     3  0.6609295613159634D+00,0.8991490109984123D+00,
     4  0.2662868781081671D+00,-.6089182592853220D+00,
     5  -.7386515048219099D+00,0.4195868080934234D-01,
     6  0.7466835339617809D+00,0.4509867981498777D+00,
     7  0.6708118222473506D-01,0.1018538399068630D+00,
     8  0.1807907136303741D+00,0.3306273231387975D+00,
     9  0.5778465571960991D+00,0.9154257781323152D+00,
     *  0.1248843961339039D+01,0.1366574743102636D+01,
     1  0.1028383194534675D+01,0.2138374958726836D+00,
     2  -.6535529965451988D+00,-.8910122932859696D+00,
     3  -.2422942365154620D+00,0.6232258515565059D+00,
     4  0.6876568814599800D+00,-.1360433727096175D+00,
     5  -.7424638071297302D+00,-.2821012781498028D+00,
     6  0.5805065947472271D+00,-.2025148455109162D+00,
     7  -.2766016887826660D+00,-.4280941711540518D+00,
     8  -.6747456024640888D+00,-.1002877583328698D+01,
     9  -.1317756673968590D+01,-.1416148685436081D+01,
     *  -.1067071351873962D+01,-.2473652970165125D+00/
c 
c        following is block number   2 of length 180 elements
c 
        data www2/
     1  0.6324465824661565D+00,0.8908743491970442D+00,
     2  0.2550217397079575D+00,-.6086016197195703D+00,
     3  -.6626560549070983D+00,0.1651084874118315D+00,
     4  0.7180857132423628D+00,0.1872150310506188D+00,
     5  -.6157047512471094D+00,-.4375821615977410D+00,
     6  0.5346365875046633D+00,0.6576858555026423D+00,
     7  0.8795175644365801D+00,0.1171894117397706D+01,
     8  0.1435842270300742D+01,0.1480332565196199D+01,
     9  0.1096177912702323D+01,0.2685693831718437D+00,
     *  -.6120698784705333D+00,-.8894491111846026D+00,
     1  -.2818378791301778D+00,0.5765556737707919D+00,
     2  0.6517815278854878D+00,-.1534803671913478D+00,
     3  -.6879063498588831D+00,-.1494503427338255D+00,
     4  0.6050601066212892D+00,0.3468453469674721D+00,
     5  -.4885603644329505D+00,-.1174587442881460D+01,
     6  -.1306821663551457D+01,-.1498191668249683D+01,
     7  -.1642865803761555D+01,-.1566464740113204D+01,
     8  -.1096782160709069D+01,-.2456847629452971D+00,
     9  0.6154588168607213D+00,0.8836425292455696D+00,
     *  0.2997870454918602D+00,-.5414160101001283D+00,
     1  -.6427866941398885D+00,0.1234820902269934D+00,
     2  0.6549422798398004D+00,0.1465161743904082D+00,
     3  -.5730804360177194D+00,-.3014205161397220D+00,
     4  0.4844745149942041D+00,0.4005657794853648D+00,
     5  0.2030726309007622D+01,0.2052581194508255D+01,
     6  0.1996191668115497D+01,0.1716132094407600D+01,
     7  0.1089629801471987D+01,0.1715753425877852D+00,
     8  -.6732696990277240D+00,-.9024357061038303D+00,
     9  -.3082671796114999D+00,0.5291822024937003D+00,
     *  0.6502225472476071D+00,-.9217796197341171D-01,
     1  -.6393284279242578D+00,-.1662821410413297D+00,
     2  0.5474380376088241D+00,0.2949645167578165D+00,
     3  -.4674018232770187D+00,-.3619082811843943D+00,
     4  0.4105103736741037D+00,-.2785154583589062D+01,
     5  -.2545250212793360D+01,-.2028549642533136D+01,
     6  -.1179596957061272D+01,-.1118852715354544D+00,
     7  0.7872556217661573D+00,0.9982943278048069D+00,
     8  0.3515539632148817D+00,-.5490888289481233D+00,
     9  -.7142815602177789D+00,0.4844753231799871D-01,
     *  0.6684629263739596D+00,0.2220039038525612D+00,
     1  -.5507429683982805D+00,-.3369972479718129D+00,
     2  0.4654763291657691D+00,0.3803723678430881D+00,
     3  -.4204823566183548D+00,-.3991498435857849D+00,
     4  0.3229748331892971D+01,0.2613390973602822D+01,
     5  0.1547714522780814D+01,0.2367530664408588D+00,
     6  -.8606220351767815D+00,-.1173124909504114D+01,
     7  -.4987896906230660D+00,0.5399972217598014D+00,
     8  0.8390682742634288D+00,0.5798019245500285D-01,
     9  -.7161273967263700D+00,-.3397728018826645D+00,
     *  0.5519534798425375D+00,0.4352676152569210D+00,
     1  -.4547804232947165D+00,-.4521971972778701D+00,
     2  0.4203033405389235D+00,0.4402754764322203D+00,
     3  -.4232064749295588D+00,-.3418671066367719D+01,
     4  -.2367521875008576D+01,-.7982291546323745D+00,
     5  0.6876131261380227D+00,0.1347103714937450D+01,
     6  0.8039902275647222D+00,-.3914511143669234D+00,
     7  -.9730444825039706D+00,-.2792075173205382D+00,
     8  0.7088044329819487D+00,0.5284241425211344D+00,
     9  -.4860003050273246D+00,-.5794420378156649D+00,
     *  0.3833141632555523D+00,0.5617984916965793D+00,
     1  -.3703413287034019D+00,-.5190323745591689D+00,
     2  0.4074757310369858D+00,0.4638156726328138D+00,
     3  0.3484896112448440D+01,0.1973892435414503D+01,
     4  0.1934095918821288D-01,-.1315714707205642D+01,
     5  -.1238108817568239D+01,-.2742146493716400D-02,
     6  0.1011573893628038D+01,0.6226075498293891D+00,
     7  -.5619117935491863D+00,-.7571447273390933D+00,
     8  0.2957226260595604D+00,0.7339171761306533D+00,
     9  -.2101863364334722D+00,-.6816615798941301D+00,
     *  0.2386080004663811D+00,0.6202982864440537D+00,
     1  -.3279735166487067D+00,-.5403241425961000D+00,
     2  0.4382845001178185D+00,-.3494858006834202D+01,
     3  -.1516501847406266D+01,0.6698654162758537D+00,
     4  0.1567934919329674D+01,0.6899378366685549D+00,
     5  -.7728488870546077D+00,-.1006264352093023D+01,
     6  0.1877360708604334D+00,0.9301857065924341D+00,
     7  0.5761626380834304D-01,-.8197417666269410D+00,
     8  -.8814494920989526D-01,0.7539358378557571D+00,
     9  -.2280788241721384D-02,-.7076762070815967D+00,
     *  0.1615281586033398D+00,0.6374503642928676D+00,
     1  -.3470663237707857D+00,-.5160752477930560D+00,
     2  0.3470326526476252D+01,0.1031882707198996D+01,
     3  -.1213443543340853D+01,-.1466162432856861D+01,
     4  0.5214983984517517D-01,0.1184241624919244D+01,
     5  0.4368914462571924D+00,-.8653144856822403D+00,
     6  -.5453739035808656D+00,0.7062296041324561D+00,
     7  0.4853195712428049D+00,-.6801598502142920D+00,
     8  -.3326418604686691D+00,0.7083785131967673D+00,
     9  0.1066603482931279D+00,-.7096348024044857D+00,
     *  0.1658050561925785D+00,0.6226409659683545D+00/
c 
c        following is block number   3 of length 180 elements
c 
        data www3/
     1  -.4299712754935972D+00,-.3419185076636970D+01,
     2  -.5417201206992392D+00,0.1583593284852183D+01,
     3  0.1084482504110075D+01,-.7462124370484869D+00,
     4  -.1100702085575433D+01,0.3519570679068306D+00,
     5  0.9759220406927361D+00,-.2619954209061519D+00,
     6  -.8415261605929214D+00,0.3512649943356449D+00,
     7  0.6830313133438107D+00,-.5217643482027434D+00,
     8  -.4482627639379957D+00,0.6741809443159581D+00,
     9  0.1191071925852134D+00,-.7071855106331635D+00,
     *  0.2522025664430541D+00,0.5594163031930549D+00,
     1  0.3345397574048634D+01,0.6251501073375510D-01,
     2  -.1770654638581352D+01,-.5248568132553621D+00,
     3  0.1207058569141622D+01,0.5959954351785003D+00,
     4  -.9508861718926831D+00,-.4794119586097699D+00,
     5  0.8731625719956651D+00,0.2453630247539127D+00,
     6  -.8470661573816469D+00,0.7953065435255826D-01,
     7  0.7428507991211620D+00,-.4310549243200154D+00,
     8  -.4693056890579941D+00,0.6791572694096249D+00,
     9  0.3876284782134648D-01,-.6854220505296192D+00,
     *  0.4095741649981866D+00,-.3251540030463391D+01,
     1  0.3920238693937616D+00,0.1779361173771212D+01,
     2  -.1012986036782702D+00,-.1335146326153554D+01,
     3  0.1152223436041261D+00,0.1083606989138183D+01,
     4  -.2997028229127054D+00,-.8527354206616420D+00,
     5  0.5602425674505525D+00,0.5232043102720318D+00,
     6  -.7637816686449153D+00,-.5800910450940894D-01,
     7  0.7446350289851109D+00,-.4313408167842108D+00,
     8  -.4101641504842388D+00,0.7119209071218954D+00,
     9  -.1337203947952425D+00,-.5979600936173106D+00,
     *  0.3139548502820897D+01,-.8102923991309618D+00,
     1  -.1626245047423932D+01,0.6885237461380238D+00,
     2  0.1123163195740552D+01,-.7703112458642766D+00,
     3  -.7119986932545513D+00,0.8978969908071542D+00,
     4  0.2334831511153674D+00,-.9039308035084066D+00,
     5  0.3010422372827193D+00,0.6298953262703849D+00,
     6  -.7053329709967272D+00,-.7013294759927101D-01,
     7  0.7177277130303944D+00,-.5140066727750421D+00,
     8  -.2597729747519309D+00,0.7242555017540853D+00,
     9  -.3811308079938524D+00,-.3011015508080565D+01,
     *  0.1182561884019044D+01,0.1337184541415791D+01,
     1  -.1149513784010782D+01,-.6436832284670040D+00,
     2  0.1149588253969849D+01,0.2506549707649848D-01,
     3  -.9831984351725883D+00,0.5437545298520552D+00,
     4  0.5154212402172251D+00,-.8469955183444879D+00,
     5  0.1750388306502878D+00,0.6360557466219570D+00,
     6  -.7053811223123107D+00,0.3061515587547795D-01,
     7  0.6432771037597356D+00,-.6421339699357578D+00,
     8  -.7319913013334102D-03,0.6326576997949838D+00,
     9  0.2867336198674757D+01,-.1500907014424504D+01,
     *  -.9448901360957919D+00,0.1423625310020874D+01,
     1  0.2328664731812178D-01,-.1140633616082532D+01,
     2  0.6643045471290169D+00,0.5242539709216621D+00,
     3  -.9523693000745908D+00,0.2834665138674417D+00,
     4  0.6135511369324849D+00,-.8011504852875482D+00,
     5  0.1835777542758723D+00,0.5666433778086783D+00,
     6  -.7441922940537973D+00,0.2331058170562830D+00,
     7  0.4688608967769730D+00,-.7351412030208522D+00,
     8  0.3465995513505992D+00,-.2709789412652735D+01,
     9  0.1759191490216398D+01,0.4863248176944702D+00,
     *  -.1481384938157860D+01,0.5902604835167086D+00,
     1  0.7623492485265085D+00,-.1060903500309131D+01,
     2  0.2154346474656724D+00,0.7300794373820998D+00,
     3  -.8597956242272939D+00,0.1789261601726226D+00,
     4  0.5986595130658980D+00,-.7948683994264249D+00,
     5  0.3086555035459126D+00,0.4027032521197360D+00,
     6  -.7566132044824051D+00,0.4983952948710378D+00,
     7  0.1436837944783530D+00,-.6637286340882520D+00,
     8  0.2539585525226593D+01,-.1953068112308434D+01,
     9  -.1272013330719693D-03,0.1325002903806510D+01,
     *  -.1060962634376401D+01,-.1476597126740279D+00,
     1  0.1007146569143784D+01,-.8353383614700629D+00,
     2  -.3352753930340484D-01,0.7707189683532868D+00,
     3  -.8120146150006206D+00,0.2200730679362182D+00,
     4  0.4896838865901735D+00,-.7931582870797675D+00,
     5  0.5138190652232347D+00,0.1115553116750277D+00,
     6  -.6375118096294711D+00,0.7154893819235754D+00,
     7  -.3072135990552051D+00,-.2357895681794304D+01,
     8  0.2079971544329906D+01,-.4758641922400914D+00,
     9  -.9852350898329978D+00,0.1291696800402604D+01,
     *  -.5052705189993485D+00,-.5374953200738838D+00,
     1  0.1013310402173437D+01,-.6821327459862665D+00,
     2  -.9273227044556077D-01,0.7232092195455261D+00,
     3  -.8175532264325585D+00,0.3798837966689034D+00,
     4  0.2652018140710691D+00,-.7118090489849438D+00,
     5  0.7118558399417437D+00,-.2921299572394262D+00,
     6  -.2874682484988102D+00,0.6911589663483335D+00,
     7  0.2165869890028855D+01,-.2139091995556370D+01,
     8  0.9068789072236498D+00,0.5152510291387089D+00,
     9  -.1241132612278658D+01,0.9963217575663528D+00,
     *  -.1469835334575772D+00,-.6661985427367746D+00/
c 
c        following is block number   4 of length 180 elements
c 
        data www4/
     1  0.9629739574963783D+00,-.6547251052724243D+00,
     2  0.3529526885999338D-02,0.5905610259230950D+00,
     3  -.8199084773697679D+00,0.6037560686505530D+00,
     4  -.9013000258602494D-01,-.4453942364350668D+00,
     5  0.7408213356259571D+00,-.6644069416563595D+00,
     6  0.2638626370279016D+00,-.1964648108806327D+01,
     7  0.2131323527404154D+01,-.1263173216276063D+01,
     8  0.1758895212774843D-01,0.9283452888004270D+00,
     9  -.1182238554847524D+01,0.7671072414842106D+00,
     *  -.1934080716999307D-01,-.6391031159046838D+00,
     1  0.9112112403562429D+00,-.7278520824394803D+00,
     2  0.2296536160387343D+00,0.3314395492526179D+00,
     3  -.7122042476325641D+00,0.7709402623890110D+00,
     4  -.5070162034491158D+00,0.4463162302117599D-01,
     5  0.4243743609622472D+00,-.7148504803988695D+00,
     6  0.1755366757294460D+01,-.2059184008757356D+01,
     7  0.1521592199183113D+01,-.5409948339008957D+00,
     8  -.4253522894421097D+00,0.1015087847991818D+01,
     9  -.1078532630580577D+01,0.6916956749298072D+00,
     *  -.8062419090345268D-01,-.4910781653803836D+00,
     1  0.8209177023882905D+00,-.8227588461626221D+00,
     2  0.5328327039255889D+00,-.7744213782512729D-01,
     3  -.3792712455355352D+00,0.6891280097371445D+00,
     4  -.7612350497430138D+00,0.5827882701730380D+00,
     5  -.2172488062319981D+00,-.1539162129981437D+01,
     6  0.1926705979383933D+01,-.1666668825871208D+01,
     7  0.9869546146298928D+00,-.1600761564052157D+00,
     8  -.5513262127689966D+00,0.9648383187096152D+00,
     9  -.1015023522536124D+01,0.7487543544406761D+00,
     *  -.2902720890746753D+00,-.2063345744939509D+00,
     1  0.6015025853319685D+00,-.8015314869173176D+00,
     2  0.7731508823540570D+00,-.5422525829052936D+00,
     3  0.1804824792732120D+00,0.2154638851296148D+00,
     4  -.5469562161967131D+00,0.7346788637192610D+00,
     5  0.1317171664667833D+01,-.1739299302355385D+01,
     6  0.1691218951334923D+01,-.1299869680965536D+01,
     7  0.7074308625816795D+00,-.6895657138829051D-01,
     8  -.4771205352278049D+00,0.8339294087295244D+00,
     9  -.9581274108805436D+00,0.8585804555776558D+00,
     *  -.5852854339842515D+00,0.2131911043026928D+00,
     1  0.1749222957286029D+00,-.5035223878518467D+00,
     2  0.7156253830704402D+00,-.7796011382329548D+00,
     3  0.6914885135063064D+00,-.4732875909908397D+00,
     4  0.1679632038644692D+00,-.1090533688013171D+01,
     5  0.1503587726735387D+01,-.1596418939864032D+01,
     6  0.1442772336187699D+01,-.1107635636395163D+01,
     7  0.6650380195296896D+00,-.1915258117199249D+00,
     8  -.2436513769908789D+00,0.5864643037651188D+00,
     9  -.8026016233963401D+00,0.8785476568997846D+00,
     *  -.8199332462210695D+00,0.6480500000363666D+00,
     1  -.3953559093741969D+00,0.1006339869287407D+00,
     2  0.1957146958542908D+00,-.4559156183855197D+00,
     3  0.6483904815552997D+00,-.7505197895016527D+00,
     4  0.8603860635533427D+00,-.1227222456663468D+01,
     5  0.1391381068221115D+01,-.1401032952642906D+01,
     6  0.1283811258258360D+01,-.1067787755050902D+01,
     7  0.7836920778285456D+00,-.4633351714244319D+00,
     8  0.1373085242233652D+00,0.1670757615140190D+00,
     9  -.4272937642611249D+00,0.6265258560448591D+00,
     *  -.7540495522085835D+00,0.8052123702679723D+00,
     1  -.7810854691676136D+00,0.6878824459640406D+00,
     2  -.5362088786542344D+00,0.3401919769597060D+00,
     3  -.1165274768889852D+00,-.6278640397687938D+00,
     4  0.9186765913481395D+00,-.1092269144084210D+01,
     5  0.1183267676864186D+01,-.1204959342628096D+01,
     6  0.1166430208834496D+01,-.1076266422677171D+01,
     7  0.9434476115781312D+00,-.7774121123368449D+00,
     8  0.5878566500535651D+00,-.3844675260836291D+00,
     9  0.1766564431598220D+00,0.2667438321535171D-01,
     *  -.2173326012111058D+00,0.3880033686047275D+00,
     1  -.5323877388782083D+00,0.6453148925146611D+00,
     2  -.7228303837892069D+00,0.7622616868488053D+00,
     3  0.3940975124687553D+00,-.5870248926048114D+00,
     4  0.7210213005577220D+00,-.8194781465770417D+00,
     5  0.8903818999041364D+00,-.9374388741213849D+00,
     6  0.9627680040269463D+00,-.9678375197938190D+00,
     7  0.9538626752714523D+00,-.9219971940581922D+00,
     8  0.8734328488803835D+00,-.8094527992673186D+00,
     9  0.7314588911604609D+00,-.6409826739610458D+00,
     *  0.5396852773933695D+00,-.4293491499311388D+00,
     1  0.3118636299336174D+00,-.1892058095707313D+00,
     2  0.6341789197385484D-01,-.1602078599843602D+00,
     3  0.2417137680119200D+00,-.3037672886055286D+00,
     4  0.3567580614691926D+00,-.4044447387790943D+00,
     5  0.4485515721135924D+00,-.4899134740548339D+00,
     6  0.5288981440838054D+00,-.5656030607348668D+00,
     7  0.5999632614407025D+00,-.6318170397428307D+00,
     8  0.6609489981804520D+00,-.6871196786304904D+00,
     9  0.7100866633538605D+00,-.7296199645067449D+00,
     *  0.7455134123434719D+00,-.7575931103593522D+00/
c 
c        following is block number 5 of length 40 elements
c 
        data www5/
     1  0.7657236317919674D+00,-.7698123812623279D+00,
     2  -.7269552845737379D-01,0.1076733629754942D+00,
     3  -.1309573367185772D+00,0.1468436979298362D+00,
     4  -.1569372765745747D+00,0.1621605664603444D+00,
     5  -.1632038184565741D+00,0.1606537733763540D+00,
     6  -.1550309965027610D+00,0.1467996323634537D+00,
     7  -.1363704985953814D+00,0.1241042159291339D+00,
     8  -.1103158761014064D+00,0.9528096779055763D-01,
     9  -.7924186461600833D-01,0.6241426898636036D-01,
     *  -.4499323147416242D-01,0.2715858181039333D-01,
     1  -.9079770676491171D-02,0.3035189373743966D+00,
     2  -.4516500199158431D+00,0.5540887472691450D+00,
     3  -.6295327153679563D+00,0.6852749134557665D+00,
     4  -.7256190178770777D+00,0.7537641743805975D+00,
     5  -.7723452317311426D+00,0.7835950924065715D+00,
     6  -.7893966683213904D+00,0.7913102944706223D+00,
     7  -.7906023703137268D+00,0.7882795837331825D+00,
     8  -.7851262578294904D+00,0.7817412183498720D+00,
     9  -.7785713599899180D+00,0.7759402274569689D+00,
     *  -.7740708533241166D+00,0.7731026802733358D+00/
c 
c       return the varous types of data to the user
c 
        do 1200 i=1,nn/2
c 
        tt(i)=tt0(i)
        tt(nn-i+1)=-tt0(i)
c 
        ww(i)=ww0(i)
        ww(nn-i+1)=ww0(i)
 1200 continue
c 
        do 1400 i=1,40
c 
        cm(i)=cm0(i)
c 
        rlams(i)=rlams0(i)*ima**(i-1)
 1400 continue
c 
        do 1500 i=1,40
c 
        rlams(i)=rlams0(i)*ima**(i-1)
 1500 continue
c 
        d=1
        do 1800 i=1,40
c 
        do 1600 j=1,19
c 
        www(j,i)=www10(j,i)
        www(38-j+1,i)=www10(j,i)*d
c 
 1600 continue
        d=-d
 1800 continue
c 
        return
c 
c 
c 
c 
        entry expo_ttww_ret20(tt,ww,nnout,c)
c 
c       This entry returns to the user the 38 prolate nodes tt and
c       their corresponding weights ww, discretizing to 16-digit
c       (or so) accuracy band-limited functions with c=20 on the
c       interval [-1,1]. It also returns nnot=38, and c=20
c 
c       return the nodes and weights to the user
c 
        do 2200 i=1,nn/2
c 
        tt(i)=tt0(i)
        tt(nn-i+1)=-tt0(i)
c 
        ww(i)=ww0(i)
        ww(nn-i+1)=ww0(i)
 2200 continue
c 
        c=20
c 
        nnout=nn
        return
        end
  
  
  
  
  
  
