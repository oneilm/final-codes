        implicit real *8 (a-h,o-z)
        real *8 w(400 000),x(10 000),whts(10 000),
     1      expons(10 000)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER c'
        READ *,c
        CALL PRIN2('c=*',c,1 )
C 
        PRINT *, 'ENTER cc'
        READ *,cc
        CALL PRIN2('cc=*',cc,1 )
C 
        PRINT *, 'ENTER -log(eps) base 10'
        READ *,db
        CALL PRIN2('db=*',db,1 )
c 
        eps=0.1**db
        call prin2('and eps=*',eps,1)
c 
        PRINT *, 'ENTER rl'
        READ *,rl
        CALL PRIN2('rl=*',rl,1 )
c 
c       construct the whole prolate junk
c 
        lenw=400 000
c 
        call prolall(ier,c,cc,eps,w,lenw,x,whts,expons,n,
     1      valmax,lused7,keep)
  
         call prinf('after prolall, ier=*',ier,1)
         call prinf('after prolall, lused7=*',lused7,1)
         call prinf('after prolall, keep=*',keep,1)
         call prin2('after prolall, valmax=*',valmax,1)
c 
c        test the obtained exponsntial expander
c 
        nn=200
        call proltst(w,x,n,nvects,expons,rl,nn)
  
  
  
  
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine proltst(w,x,n,nvects,expons,rl,nn)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),expons(1),w(1),ww(10000),rea(2),
     1      wrea(2000),wima(2000),tt(2000)
c 
        complex *16 coefs(1000),f(1000),ima,vals(1000),com
        data ima/(0.0d0,1.0d0)/
        equivalence (rea(1),com)
c 
c        construct the test function to be expanded
c 
        do 1200 i=1,n
        f(i)=cdexp(ima*x(i)*rl)
  
ccccc        f(i)=x(i)**10
  
  
ccc             f(i)=x(i)
 1200 continue
  
cccc        f(n/4)=f(n/4)+1.0d-2
c 
c       construct the exponential representation of the test function
c 
        call prolalev(f,w,coefs)
c 
        call prin2('the exponential coefficients are*',coefs,n*2)
c 
c        plot the exponential coefficients
c 
        do 1250 i=1,n
        com=coefs(i)
        ww(i)=rea(1)
 1250 continue
c 
        iw=21
  
  
  
  
        call lotagraph(iw,expons,ww,n,'exponential coefficients*')
  
cccc        call RSGRAF(expons,ww,n,IW,'exponential coefficients*')
c 
c       evaluate the norm of the exponential coefficients
c 
        d=0
        do 1300 i=1,n
        d=d+cdabs(coefs(i))**2
 1300 continue
          call prin2('norm of the vector of exponential coefficients*',
     1      dsqrt(d),1)
c 
c       evaluate the exponential expansion at the test points
c       to check the accuracy
c 
        do 1400 i=1,n
c 
        call prolevex(x(i),coefs,expons,n,vals(i))
 1400 continue
cccc        call prin2('and f recomputed from exponential expansion*',
cccc     1      vals,n*2)
c 
cccc        call prin2('original f is*',f,n*2)
c 
        do 1600 i=1,n
        vals(i)=vals(i)-f(i)
 1600 continue
cccc        call prin2('and the difference is*',vals,n*2)
c 
c        calculate the relative l^2 error
c 
        d1=0
        d2=0
        do 1800 i=1,n
        d1=d1+cdabs(vals(i))**2
        d2=d2+cdabs(f(i))**2
 1800 continue
        call prin2('and relative l^2 error*',dsqrt(d1/d2),1)
c 
c       plot the extended function on the interval [-2,2]
c 
        nn=200
        h=4
        h=h/(nn-1)
        do 2200 i=1,nn
        tt(i)=(i-1)*h-2
 2200 continue
c 
cccc        call prin2('tt as created in proltst*',tt,nn)
c 
c       construct the extended function
c 
        do 2400 i=1,nn
c 
        call prolevex(tt(i),coefs,expons,n,com)
c 
        wrea(i)=rea(1)
        wima(i)=rea(2)
 2400 continue
c 
        iw=22
  
        call lotagraph(iw,tt,wrea,nn,'real part of extended function*')
  
cccc        call RSGRAF(tt,wrea,nn,IW,'real part of extended function*')
c 
        iw=23
  
        call lotagraph(iw,tt,wrea,nn,
     1      'imaginary part of extended function*')
  
cccc        call RSGRAF(tt,wima,nn,IW,
cccc     1      'imaginary part of extended function*')
c 
  
c 
        return
        end
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c        this is the end of the debugging code and the beginning of
c        the actual routines for the exponantial expansion of
c        non-periodic functions
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
c 
c 
        subroutine prolalev(f,w,coefs)
        implicit real *8 (a-h,o-z)
        save
        dimension f(1),w(1),coefs(1)
c 
c        this subroutine constructs an exponential expansion for the
c        user-supplied function on the interval [-1,1]. Actually,
c        almost all the work is performed by the subroutine prolall
c        (see). This subroutine essentially applies the matrix prepared
c        by prolall to the user-supplied values of the function at the
c        Gaussian nodes on the interval [-1,1].
c 
c        Specifically, the user supplies the values f(i) at the Gaussian
c        nodes of some (hopefully, smooth) function f. This subroutine
c        constructs the set of coefficients coefs(k) such that
c 
c        f(x) \sim  \sum{k=1}^n coefs(k) * cdexp(ima*expons(k) *x),       (1)
c 
c        where both the array expons and the number n are produced by
c        the subroutine prolall. It is assumed that the function f is
c        band-limited, i.e.
c 
c        f(x)= \int _{-1}^1 \sigma(t) * cdexp(c * ima*t*x) dt,             (2)
c 
c        with c the user-supplied parameter. The acutal evaluation of
c        the expansion (1) can be performed by the subroutine prolevex
c        (see).
c 
c                          input parameters:
c 
c  f - the values at the Gaussian nodes on the interval [-1,1] of
c        the function to be decomposed
c  w - the array containing various types of information produced by
c        the subroutine prolal (see); the first keep elements
c        of this array should not be changed between the call to
c        prolalev and the subsequent calls to this subroutine
c 
c                           output parameters:
c 
c  coefs - the ceofficients in the expansion (1).
c 
c       . . . unpack the map of the array w
c 
        ialphas=w(1)
        ipsis=w(2)
        iwhts=w(3)
        nvects=w(4)
        n=w(5)
c 
cccc        call prinf('in prolalev, ialphas=*',ialphas,1)
cccc        call prinf('in prolalev, ipsis=*',ipsis,1)
cccc        call prinf('in prolalev, nvects=*',nvects,1)
cccc        call prinf('in prolalev, n=*',n,1)
c 
cccc        call prin2('in prolalev, alphas=*',w(ialphas),n*nvects*2)
cccc        call prin2('in prolalev, psis=*',w(ipsis),n*nvects)
  
c 
c        construct the exponential expansion of the user-supplied
c        function
c 
        call prolexco(f,w(ipsis),n,nvects,w(iwhts),coefs,w(ialphas))
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine prolevex(t,coefs,expons,n,cout)
        implicit real *8 (a-h,o-z)
        save
        complex *16 coefs(1),ima,cout
        real *8 expons(1)
        data ima/(0.0d0,1.0d0)/
c 
c        this subroutine evaluates at the point t the user-supplied
c        exponential expansion with n terms
c 
c        \sum{k=1}^n coefs(k) * cdexp(ima*expons(k) *t)              (1)
c 
c 
c                          input parameters:
c 
c  t - the poiint (real) where the expansion (1) is to be avaluated
c  coefs - the ceofficients in the expansion (1).
c  expons - the coefficients in (1) above (n of them)
c  n - the number of elements in the expansion (1)
c 
c                           output parameters:
c 
c  cout - the sum of the expansion (1)
c 
        cout=0
        do 1200 i=1,n
        cout=cout+coefs(i)*cdexp(ima*t*expons(i))
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine prolall(ier,c,cc,eps,w,lenw,
     1      x,whts,expons,n,valmax,lused7,keep)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),x(1),whts(1),expons(1)
c 
c        this subroutine is the initialization entry point for the
c        subroutine prolalev (see) that constructs an exponential
c        expansion for a user-supplied function on the interval [-1,1].
c        actually, this subroutine performs all the work; prolalev
c        essentially applies the matrix prepared by this subroutine
c        to the user-supplied values of the function at the Gaussian
c        nodes on the interval [-1,1].
c 
c        Specifically, the user supplies the values f(i) at the Gaussian
c        nodes of some (hopefully, smooth) function f. The subroutine
c        prolalev will construct the set of coefficients coefs(k)
c        such that
c 
c        f(x) \sim  \sum{k=1}^n coefs(k) * cdexp(ima*expons(k) *x),       (1)
c 
c 
c        where both the array expons and the number n are produced by
c        this subroutine. It is assumed that the function f is
c        band-limited, i.e.
c 
c        f(x)= \int _{-1}^1 \sigma(t) * cdexp(c * ima*t*x) dt,             (2)
c 
c        with c the user-supplied parameter.
c 
c    IMPORTANT note: The elements of the array expons (see (1) above)
c        range from -c*cc to +c*cc, with cc the user supplied
c        "stabilization constant". Thbe existence of cc has to do with
c        the funny behaviour of the problem (1) of exponential
c        interpolation. Due to space limitation, we do not discuss
c        the details here. Normally cc=2 is sugficient; under no
c        conditions should one use cc < 1.
c 
c                          input parameters:
c 
c  c - the parameter in formula (2) above: the highest permitted
c        frequency in the user-supplied function (see (2) above)
c  cc - "stabilization constant" to be used by the subroutine
c        in constructing array expons (see (1) above)
c  eps - the accuracy to which the subroutine will attempt to approximate
c        the user-supplied function
c  lenw - the amount of space provided in array w (in real *8 words);
c        should be sufficiently large
c 
c                           output parameters:
c 
c  ier - error return code. ier=0 means successful conclusion. Any
c        non-zero value returned means that the subroutine has run
c        out of memory in array w (the subroutine is informed about
c        this amount via the parameter lenw; the magnitude of ier is
c        guide to the severity of the problem. ier=2048 means that
c        the amount is GROOOSly insufficient. ier=32 means that lenw
c        was almost large enough.
c  w - the array containing various types of information to be used
c        by the subroutine prolalev (see); the first keep elements
c        of this array should not be changed between the call to
c        this subroutine and subsequent calls to the subroutine prolalev
c  x - the Gaussian nodes on the interval [-1,1]
c  whts - the Gaussian weights on the interval [-1,1]
c  expons - the coefficients in (1) above (n of them)
c  n - the number of elements in the expansion (1)
c  lused7 - the total amount of space in w that has been used by the
c        subroutine. Note that this amount is normally very much greater
c        than the amount keep (parameter below) that has to be kept
c        between the call to this subroutine and subsequent calls to the
c        subroutine prolalev (see). In other words, most of the array
c        w is used as the work array; a part of it of length keep is
c        viewed as output, to be used by the subroutine prolalev.
c  keep - the part of array w (in real *8 words) that has to be unchanged
c        between the call to this subroutine and the subsequent calls
c        to the subroutine prolalev (see)
c 
c 
c        . . . allocate memory for the first prolcrea
c 
        ier=0
        lused7=0
        c2=c*cc
c 
        ikhi1=1
        lkhi1=c*2+100
c 
        irlams1=ikhi1+lkhi1
        lrlams1=(c*2+100)*2
c 
        irlams12=irlams1+lrlams1
        lrlams12=c*2+100
c 
        iw1=irlams12+lrlams12
c 
c       construct the prolate functions for the coefficient c
c 
         lenw1=lenw-iw1
         if(lenw1 .gt. 500) goto 1100
         ier=2048
         return
 1100 continue
c 
        call prolcrea(ier1,c,w(iw1),lenw1,nvects1,nhigh1,
     1    w(ikhi1),w(irlams1),w(irlams12),keep1,lused1)
c 
        call prin2('after prolcrea khi1 are*',w(ikhi1),nvects1)
        call prin2('after prolcrea, rlams12=*',w(irlams12),nvects1)
        call prin2('after prolcrea, rlams1=*',w(irlams1),nvects1*2)
c 
        call prinf('ier1=*',ier1,1)
        call prinf('keep1=*',keep1,1)
        call prinf('lused1=*',lused1,1)
c 
        lused7=iw1+lused1
c 
        if(ier1 .eq. 0) goto 1150
        ier=1024
        return
 1150 continue
c 
c       collect the garbage
c 
        irlams121=1
        do 1200 i=1,nvects1
        w(irlams121+i-1)=w(irlams12+i-1)
 1200 continue
        irlams12=irlams121
        lrlams12=nvects1+1
c 
        iw11=irlams12+lrlams12
c 
        do 1400 i=1,keep1
        w(iw11+i-1)=w(iw1+i-1)
 1400 continue
        iw1=iw11
c 
        lw1=keep1+1
        lused=iw1+lw1
c 
c        allocate memory for the second prolcrea
c 
        ikhi2=lused+1
        lkhi2=c2*2+100
c 
        irlams2=ikhi2+lkhi2
        lrlams2=(c2*2+100)*2
c 
        irlams22=irlams2+lrlams2
        lrlams22=c2*2+100
c 
        iw2=irlams22+lrlams22
c 
c       construct the prolate functions for the coefficient c2
c 
        lenw2=lenw-iw2
        if(lenw2 .gt. 500) goto 1600
        ier=256
        return
 1600 continue
c 
        call prolcrea(ier2,c2,w(iw2),lenw2,nvects2,nhigh2,
     1    w(ikhi2),w(irlams2),w(irlams22),keep2,lused2)
c 
        call prin2('after prolcrea khi2 are*',w(ikhi2),nvects2)
        call prin2('after prolcrea, rlams22=*',w(irlams22),nvects2)
        call prin2('after prolcrea, rlams2=*',w(irlams2),nvects2*2)
c 
        call prinf('ier2=*',ier2,1)
        call prinf('keep2=*',keep2,1)
        call prinf('lused2=*',lused2,1)
c 
        lused8=iw2+lused2
        if(lused8 .gt. lused7) lused7=lused8
c 
        if(ier2 .eq. 0) goto 1800
        ier=128
        return
 1800 continue
c 
c       collect the garbage
c 
        irlams221=lused+1
        do 2200 i=1,nvects2
        w(irlams221+i-1)=w(irlams22+i-1)
 2200 continue
        irlams22=irlams221
        lrlams22=nvects2+1
c 
        iw21=irlams22+lrlams22
c 
        do 2400 i=1,keep2
        w(iw21+i-1)=w(iw2+i-1)
 2400 continue
        iw2=iw21
c 
        lw2=keep2+1
        lused=iw2+lw2
c 
c       adjust nvects1, nvects2 to eliminate those eigenvectors
c       coresponding to very small eigenvalues
c 
        eps7=eps/10
        do 2600 i=1,nvects1
        w(irlams12+i-1)=dsqrt(w(irlams12+i-1))
c 
        if(w(irlams12+i-1) .gt. eps7) n1=i
 2600 continue
c 
        do 2800 i=1,nvects2
        w(irlams22+i-1)=dsqrt(w(irlams22+i-1))
        if(w(irlams22+i-1) .gt. eps7) n2=i
 2800 continue
c 
        nvects1=n1
        nvects2=n2
c 
        call prinf('new nvects1 is*',nvects1,1)
        call prinf('new nvects2 is*',nvects2,1)
c 
        call prin2('w(irlams12) are*',w(irlams12),n1)
        call prin2('w(irlams22) are*',w(irlams22),n2)
c 
c       allocate memory for the subroutine projeva
c 
        iprojs=lused+2
        lprojs=nvects1*nvects2+2
c 
        iwork=iprojs+lprojs
        lwork=nhigh2+1
c 
        iwork1=iwork+lwork
        lwork1=nhigh2+1
        lwork1=lwork1*2+100
c 
        iwork2=iwork1+lwork1
        lwork2=nhigh2+1
        lwork2=lwork2*2+100
c 
        inproj=iwork2+lwork2
        lnproj=nvects1+2
c 
        ltot=inproj+lnproj
c 
        if(ltot .gt. lused7) lused7=ltot
c 
        if(ltot .lt. lenw) goto 3000
        ier=64
        return
 3000 continue
c 
c       evaluate the projections of the eigenfunctions of the
c       first operator on the eigenfunctions of the second
c 
          call prinf('before projeva, iw1=*',iw1,1)
c 
        call projeva(w(iw1),nvects1,w(iw2),nvects2,w(irlams12),
     1    w(irlams22),eps/5,w(iprojs),w(iwork),w(iwork1),
     2      w(iwork2),w(inproj),nhigh2,valmax)
  
        call prin2('after projeva, valmax=*',valmax,1)
c 
c        allocate memory for the subroutine prolexex
c 
        n=nhigh2+2
c 
        ialphas=iprojs+lprojs
        lalphas=nvects1*n*2+10
c 
        ipsis=ialphas+lalphas
        lpsis=nvects1*n+5
c 
        icoelege=ipsis+lpsis
        lcoelege=nhigh2+1
c 
        iphi=icoelege+lcoelege
        lphi=n+2
c 
        iwork=iphi+lphi
        lwork=nhigh2+2
c 
        ltot=iwork+lwork
c 
        if(ltot .gt. lused7) lused7=ltot
c 
        if(ltot .lt. lenw) goto 3200
        ier=32
        return
 3200 continue
c 
c       construct the exponential expansions for the prolate
c       functions corresponding to the coefficient c, and
c       evaluate them at the test nodes
c 
        call prolexex(n,w(iprojs),w(iw2),nvects1,nvects2,
     1      w(iphi),c*cc,nhigh2,w(ialphas),w(ipsis),w(iw1),
     2      x,whts,expons,w(icoelege),w(iwork) )
  
c 
c       perform the final garbage collection
c 
        ialphas2=11
        do 3400 i=1,nvects1*n*2+4
        w(ialphas2+i-1)=-w(ialphas+i-1)
 3400 continue
c 
        ialphas=ialphas2
        lalphas=nvects1*n*2+4
c 
        ipsis2=ialphas+lalphas
        do 3600 i=1,lpsis
        w(ipsis2+i-1)=w(ipsis+i-1)
 3600 continue
c 
        ipsis=ipsis2
c 
        iwhts=ipsis+lpsis
        do 3800 i=1,n
        w(iwhts+i-1)=whts(i)
 3800 continue
c 
c       pack into the beginning of the array w the map of w
c 
        w(1)=ialphas+0.1
        w(2)=ipsis+0.1
        w(3)=iwhts+0.1
        w(4)=nvects1+0.1
        w(5)=n+0.1
c 
        keep=iwhts+lwhts+2
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine prolexco(f,psis,n,nvects,
     1    whts,coefs,alphas)
        implicit real *8 (a-h,o-z)
        save
        complex *16 f(1),coefs(1),cd,alphas(n,1)
        dimension psis(n,1),whts(1)
c 
c        this subroutine constructs the exponential expansion of the
c        user-supplied function f. It uses the parameters psis, whts,
c        alphas produced by a prior call to the subroutine prolexex
c        (see), and has no function as a stand-alone device.
c 
c 
c                          input parameters:
c 
c  f - the values at the n Gaussian nodes of the function whose exponential
c        expansion is to be constructed
c  psis - the values of the first nvects1 prolate functions corresponding
c        to the parameter c. Specifically, i-th column of the array psis
c        contains values of the function of order i-1 at the n Gaussian
c        nodes on the interval [-1,1]. This is an output of the subroutine
c        prolexex (see).
c  n - the number of Gaussian nodes to be used in all calculations
c  nvects - the number of prolate eigenvectors corresponding to the
c        parameter c whose eigenvalues are greater than eps. Must have been
c        produced by the subroutine prolcrea (see)
c  whts - the Gaussian weights on the interval [-1,1]
c  alphas - the coefficients of the exponential expansions of the first
c        nvects1 prolate functions corresponding to the coefficient c1.
c         This is an output of the subroutine prolexex (see).
c 
c                           output parameters:
c 
c  coefs - coefficients of the exponential expansion of the function f
c 
c         . . .  initialize the output array to zero
c 
        do 1200 i=1,n
        coefs(i)=0
 1200 continue
c 
        do 2000 i=1,nvects
c 
c         . . . calculate the projection of f on the
c               i-th prolate function
c 
        cd=0
        do 1400 j=1,n
        cd=cd+f(j)*psis(j,i)*whts(j)
 1400 continue
c 
c       augment the output array
c 
        do 1600 j=1,n
        coefs(j)=coefs(j)+cd*alphas(j,i)
 1600 continue
c 
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine prolexex(n,projs,w2,nvects1,nvects2,
     1      phi,c2,nhigh,alphas,psis,w1,x,whts,expons,
     2      coelege,work)
        implicit real *8 (a-h,o-z)
        save
        dimension w2(1),projs(nvects2,1),w2(1),x(1),whts(1),phi(1),
     1      coelege(1000),work(1000),psis(n,1),w1(1),expons(1)
        complex *16 alphas(n,1)
c 
c        this subroutine constructs the matrix converting the
c        coefficients of the prolate expansion of a function into
c        the coefficients of that function's exponentiual expansion.
c        it also constructs the values of the first prolate functions
c        at the Gaussian nodes on the interval [-1,1], and several
c        other parameters. This subroutine uses as input various data
c        produvced by the subroutines prolcrea, projeva (see). It has
c        no function as a stand-alone device.
c 
c                          input parameters:
c 
c  n - the number of Gaussian nodes to be used in all calculations
c  projs - the array of projections of prolate functions corresponding
c        to the parameter c on the prolate functions corresponding
c        to the parameter c2. This parameter must have been produced by
c        a prior call to the subroutine projeva (see)
c  w2 - the principal output of the subroutine prolcrea (see) correspo-
c        nding  to the parameter c. Its name is w as the parameter
c        in prolcrea
c  nvects1 - the number of prolate eigenvectors corresponding  to the
c        parameter c whose eigenvalues are greater than eps
c  nvects2 - the number of prolate eigenvectors corresponding  to the
c        parameter c2 whose eigenvalues are greater than eps
c  c2 - the greter of the two coefficients for which the prolate functions
c        have been created (by prior calls to the subroutine prolcrea)
c  nhigh - the highest order of the Legendre expansion used to express
c       any of the nvects2 prolate functions (corresponding to the
c       coefficient c2)
c  w1 - the principal output of the subroutine prolcrea (see) correspo-
c        nding  to the parameter c2. Its name is w as the parameter
c        in prolcrea.
c 
c                           output parameters:
c 
c  alphas - the coefficients of the exponential expansions of the first
c        nvects1 prolate functions corresponding to the coefficient c1.
c        More specifically, the i-th column of the matrix alphas contains
c        the coefficients of the exponential expansion of the i-1-st
c        prolate function, so that
c 
c   \psi_{i-1} (x) = \sum_{j=1}^n alphas(j,i) * e^{ima * expons(j) *x)},   (1)
c 
c        with the coefficients expons(j) also produced by this subroutine
c        (see below)
c  psis - the values of the first nvects1 prolate functions corresponding
c        to the parameter c. Specifically, i-th column of the array psis
c        contains values of the function of order i-1 at the n Gaussian
c        nodes on the interval [-1,1].
c  x - the Gaussian nodes on the interval [-1,1]
c  whts - the Gaussian weights on the interval [-1,1]
c  expons - the coefficients in (1) above (n of them)
c 
c                           work arrays:
c 
c  phi - must be at least n+1 real *8 locations long
c  coelege, work - must be at least nhigh+1 real *8 locations long
c 
c        . . . construct the Gaussian nodes and weights on
c              the interval [-1,1]
c 
        itype=1
        call legeexps(itype,n,x,u,v,whts)
c 
c        one prolate function psi_{c,j) after another, construct the
c        values at the legendre nodes on the interval [-1,1]
c        of such function \phi_{c,j} that
c 
c        psi_{c,j}(x)= \int_{-1}^1 phi_{c,j)(t) * cdexp(ima*c*x*t) dr
c 
        done=1
        pi=datan(done)*4
        coefix=c2/pi/2
        coefix=dsqrt(coefix)
  
        do 2000 i=1,nvects1
c 
c       construct the Legendre expansion of the function phi_{c,i}
c 
        call prolcons(projs(1,i),w2,nvects2,coelege,nnn,work,nhigh)
c 
c       evaluate the resulting expansion at the Legendre nodes
c       and construct the matrix converting prolate expansions
c       into exponential ones
c 
        do 1400 j=1,n
c 
        call legeexev(X(j),phi(j),coelege,Nnn)
        phi(j)=phi(j)*coefix
 1400 continue
c 
        jj=i-1
        call proltoex(phi,jj,n,c2,t,x,whts,alphas(1,i) )
 2000 continue
c 
c        evaluate each of the prolate functions at all legendre nodes
c 
        do 2600 i=1,nvects1
c 
        do 2400 j=1,n
        call proleva(i-1,x(j),w1,psis(j,i),der)
 2400 continue
 2600 continue
c 
c        construct the array of exponents in the exponential expansion
c        for the prolate spheroidal functions
c 
        do 2800 i=1,n
        expons(i)=c2*x(i)
 2800 continue
        return
        end
  
c 
c 
c 
c 
c 
        subroutine proltoex(phi,j,n,c2,t,x,whts,coefs)
        implicit real *8 (a-h,o-z)
        save
        dimension phi(1),x(1),whts(1)
        complex *16 ima,cd,coefs(1)
        data ima/(0.0d0,1.0d0)/
c 
c        evaluate the integral expression for the prolate function
c 
        jj=(j+1)/4
        jjj=j+1-jj*4
        cd=(-ima)**j
        if(jjj .eq. 0) cd=-cd
        if(jjj .eq. 1) cd=-cd
        do 2000 i=1,n
c 
        coefs(i)=phi(i)*whts(i)*cd
 2000 continue
        val=cd
        return
        end
c 
c 
c 
c 
c 
        subroutine prolcons(coefs,w,nvects,coelege,n,work,nhigh)
        implicit real *8 (a-h,o-z)
        save
        dimension coefs(1),w(1),coelege(1),work(1)
c 
c       This subroutine converts the user-supplied prolate expansion
c       into a Legendre expansion. Note that this subroutine uses
c       parameters produced by a prior call to the subroutine prolcrea
c       (see). It has no stand-alone function
c 
c                         Input parameters:
c 
c  coefs - the coefficients of a Legendre expansion
c  w - the principal output of the subroutine prolcrea, produced by
c       a prior call to it
c  nvects - the number of solutions of the equation (1) whose Legendre
c       expansions have been computed; is an output of subroutine prolcrea
c  nhigh - the highest order of the Legendre expansion used to express
c        any of the prolate function; is an output of subroutine prolcrea
c 
c                        output parameters:
c 
c  coelege - the coefficients of the Legendre expansion that is equal to
c        the prolate expansion with the coefficients coefs
c  n - the order of the Legendre expansion with coefficients coelege
c 
c                        Work arrays:
c 
c  work - must be at least nhigh+1 real *8 locations long
c 
        do 1200 i=1,nhigh
        coelege(i)=0
 1200 continue
c 
        nn=0
        done=1
c 
        isign=1
c 
        do 2000 i=1,nvects
  
        j=i/2
        jj=i-j*2
        if(jj .eq. 0) isign=-isign
  
c 
        call prolunpk(i-1,w,work,n)
c 
        do 1400 j=1,n
        coelege(j)=coelege(j)+work(j)*coefs(i) * isign
 1400 continue
c 
        if(nn .lt. n) nn=n
 2000 continue
        n=nn
        return
        end
c 
c 
c 
c 
c 
        subroutine projeva(w1,nvects1,w2,nvects2,rlams12,
     1    rlams22,eps,projs,whts,coefs1,coefs2,nproj,nhigh,
     2      valmax)
        implicit real *8 (a-h,o-z)
        save
        dimension w1(1),w2(1),rlams22(1),whts(1),rlams12(1),
     1      coefs1(1),coefs2(1),projs(nvects2,1),nproj(1)
c 
c       for each of the prolate functions corresponding to the
c       coefficient c, this subroutine constructs its projections
c       on all prolate functions corresponding to the parameter c2
c 
c 
c                          input parameters:
c 
c  w1 - the principal output of the subroutine prolcrea (see) correspo-
c        nding  to the parameter c2. Its name is w as the parameter
c        in prolcrea.
c  nvects1 - the number of prolate eigenvectors corresponding  to the
c        parameter c whose eigenvalues are greater than eps
c  w2 - the principal output of the subroutine prolcrea (see) correspo-
c        nding  to the parameter c. Its name is w as the parameter
c        in prolcrea
c  nvects2 - the number of prolate eigenvectors corresponding  to the
c        parameter c2 whose eigenvalues are greater than eps
c  rlams12 - the absolute values of the eigenvalues corresponding
c        to the prolate functions corresponding to the parameter c
c  rlams22 - the absolute values of the eigenvalues corresponding
c        to the prolate functions corresponding to the parameter c2
c  eps - the accuracy to which the projections are to be computed.
c        More specifivally, any projection that is smaller than eps
c        is declared to be zero.
c  nhigh - the highest order of the Legendre expansion used to express
c       any of the nvects2 prolate functions (corresponding to the
c       coefficient c2)
c 
c                           output parameters:
c 
c  projs - the array of projections of prolate functions corresponding
c        to the parameter c on the prolate functions corresponding
c        to the parameter c2.
c 
c                           work arrays:
c 
c  whts,coefs1,coefs2 - must be at least nhigh+1 real *8 locations each
c  nproj - must be at least nvects1 integer *4 locations each
c 
        do 1200 i=1,nhigh+1
        whts(i)=2*i-1
        whts(i)=1/whts(i)
        whts(i)=whts(i)*2
 1200 continue
c 
        do 3400 i=1,nvects1
  
cccc         call prinf('in projeva, i=*',i,1)
c 
        call prolunpk(i-1,w1,coefs1,n1)
c 
        do 3200 j=1,nvects2
c 
cccc         call prinf('in projeva, j=*',j,1)
  
        call prolunpk(j-1,w2,coefs2,n2)
c 
        n=n2
        if(n1 .lt. n2) n=n1
        call proscapr(coefs1,coefs2,whts,n,projs(j,i))
c 
        d=dabs(projs(j,i))
        if(d*rlams12(i) .lt. eps) goto 3200
c 
        nproj(i)=j
 3200 continue
 3400 continue
c 
c        now, construct the array of projections scaled by
c        appropriate eigenvalues
c 
        valmax=0
        do 3800 i=1,nvects1
        do 3600 j=1,nvects2
        projs(j,i)=projs(j,i)/rlams22(j)
        if(j .gt. nproj(i)) projs(j,i)=0
c 
        if(valmax .lt. dabs(projs(j,i))) valmax=dabs(projs(j,i))
 3600 continue
 3800 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine prolexev(x,n,coefs,w,fun,der)
        implicit real *8 (a-h,o-z)
        save
        dimension coefs(1),w(1)
c 
c       evaluate at the point x \in [-1,1] the prolate
c       expansion with coefficients coefs
c 
        fun=0
        der=0
c 
        do 1200 i=1,n+1
        call proleva(i-1,x,w,val,valder)
c 
        fun=fun+val*coefs(i)
        der=der+valder*coefs(i)
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine proscapr(x,y,w,n,prod)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1),w(1)
c 
        prod=0
        do 1200 i=1,n
        prod=prod+x(i)*y(i)*w(i)
 1200 continue
        return
        end
