        implicit real *8 (a-h,o-z)
        dimension xs(10 000),ys(10 000),w(1000 000),
     1      wsxs(10000),
     2      w3(1000 000)
c
        complex *16 sigma(10 000),
     1      vals(10000),vals2(10000),diffs(10000)
c
        call prini(6,13)
C
C       SET ALL PARAMETERS
C
        PRINT *, 'ENTER n'
        READ *,n
        CALL PRINf('n=*',n,1 )
c
c       prepare the test parameters
c
        mm=5
        ntest=6
        nlege=100

        ntest=90
        nlege=100
cccc        nlege=70


        m=n/2
cccc        m=n*3
cccc        m=n-1

c
c       construct the test function to be filtered
c
        itype=1
        call legeexps(itype,n,xs,u,v,wsxs)
c
        call pnmini(nlege+ntest,w3)
c
        done=1
        do 1200 i=1,n
c
        call sigma_eval(xs(i),ntest,mm,sigma(i),w3)
 1200 continue

        call prin2('sigma as created sigma=*',sigma,n*2)

cccc        stop
c
        lenw=1000 000
        ndigits=14
        call pnmproj_init(ier,ndigits,n,m,nlege,
     1      xs,wsxs,ys,w,lenw,keep)

        call prinf('after pnmproj_init, keep=*',keep,1)
c
        call pnmproj_evalc(mm,sigma,vals,w)

        call prin2('after onedfmm_proj_pnm, vals=*',vals,m*2)

cccc        stop


        done=1
        do 2200 i=1,m
c
        call sigma_eval(ys(i),ntest,mm,vals2(i),w3)
        diffs(i)=vals2(i)-vals(i)
 2200 continue

        iw=21
        call quagraph2(iw,xs,sigma,n,3,ys,vals,m,3,
     1      'vals vs. sigma*')


        call prin2('and vals2=*',vals2,m)


        call prin2('and diffs=*',diffs,m)


        stop
        end
c
c
c
c
c
        subroutine sigma_eval(t,n,mm,f,w)
        implicit real *8 (a-h,o-z)
        dimension rs(10000),ders(10000),w(1)
        complex *16 f,ima
        save
        data ima/(0.0d0,1.0d0)/
c

        call pnmeva(ier,t,n+4,mm,rs(2),ders(2),w)
c
        f=rs(n+1)*ima

cccc        call prinf('n=*',n,1)
cccc        call prin2('t=*',t,1)
cccc        call prin2('rs=*',rs,n+1)

c
        return
        end
c     
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging code, and the beginning of
c       the Spherical filtering code proper. 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     
C        This file contains three user-callable subroutines: 
c        pnmproj_init and pnmproj_eval, pnmproj_evalc. Following 
c        is a brief description of these three subroutines.
c
c   pnmproj_init - initializes the application of the spherical
c         filter to user-specified functions; the latter
c        procedure is performed by the subroutine pnmproj_eval
c        (see). Specifically, given the function 
c
c          sigma(x)=\sum{k=m}^N \alpha_k * P_k^m (x)            (1)
c
c       tabulated at the nodes xs(1),xs(2),...,xs(nxs),
c       the subroutine pnmproj_eval produces the function
c
c          sigma(x)=\sum{k=m}^{nlege} \alpha_k * P_k^m (x),     (2)
c
c       tabulated at the points ys(1),ys(2),...,ys(nys).
c
c   pnmproj_eval - applies the sphrical filter to the 
c        user-specified function sigma. Specifically, given 
c        the function 
c
c          sigma(x)=\sum{k=m}^N \alpha_k * P_k^m (x)            (1)
c
c       tabulated at the nodes xs(1),xs(2),...,xs(nxs),
c       the subroutine pnmproj_eval produces the function
c
c          sigma(x)=\sum{k=m}^{nlege} \alpha_k * P_k^m (x),     (2)
c
c       tabulated at the points ys(1),ys(2),...,ys(nys).
c
c       PLEASE NOTE THAT THIS SUBROUTINE USES THE ARRAY W THAT
C       (hopefully) has been created by a preceding call to
c       the subroutine pnmproj_init (see). This subroutine
c       has no uses known to man as a stand-alone device.       
c
c   pnmproj_evalc - a complex *16 version of the subroutine 
c       pnmproj_eval
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
c
c
        subroutine pnmproj_evalc(m,sigma,vals,w)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1)
        complex *16 sigma(1),vals(1)
c
c        PLEASE NOTE THAT THIS IS A COMPLEX VERSION OF THE
C        SUBROUTINE PNMPROJ_EVAL. IN FACT, THESE TWO FUNCTIONS
C        ARE VIRTUALLY IDENTICAL (EXCEPT FOR A FEW DECLARATIONS)
C
c        This subroutine applies the sphrical filter to the 
c        user-specified complex function sigma. Specifically, 
c        given the function 
c
c          sigma(x)=\sum{k=m}^N \alpha_k * P_k^m (x)            (1)
c
c       tabulated at the nodes xs(1),xs(2),...,xs(nxs),
c       the subroutine pnmproj_evalc produces the function
c
c          sigma(x)=\sum{k=m}^{nlege} \alpha_k * P_k^m (x),     (2)
c
c       tabulated at the points ys(1),ys(2),...,ys(nys).
c
c       PLEASE NOTE THAT THIS SUBROUTINE USES THE ARRAY W THAT
C       (hopefully) has been created by a preceding call to
c       the subroutine pnmproj_init (see). This subroutine
c       has no uses known to man as a stand-alone device.       
c
c       POTENTIALLY IMPORTANT: The value of the upper limit
c       N in the summation in (1) above depends on the nature
c       of the nodes xs and the corresponding weights wsxs. 
c       If xs are the Gaussian nodes and wsxs are the 
c       corresponding Gaussian weights, then 
c
c                    N=nxs.                                     (3)
c
c       Any other selection of nodes and weights leads to a
c       smaller N.
c
c                        Input parameters:
c
c  m - the parameters in (1) above
c  sigma - the function being filtered, tabulated at the nxs
c        nodes 
c  w - the big array generated by a preceding call to the
c       subroutine pnmproj_init (see). 
c
c                         Output parameters:
c
c  vals - the filtered version of the function sigma, tabulated
c       at the nys nodes ys (see subroutine pnmproj_init )
c
c        . . . allocate memory
c
        iw=w(1)
        iiis=w(3)
        iw2=w(2)
c
        nxs=w(4)
        nys=w(5)
c
        ii=w(iiis+m)
        ii=ii+iw2-1
c
        iw3=w(6)

c
c        apply the projection operator
c
        call pnmproj_eval0c(w(iw),sigma,vals,nxs,nys,
     1      w(ii),m,w(iw3) )
c
        return
        end
c
c
c
c
c
        subroutine pnmproj_eval0c(w,sigma,vals,
     1      nxs,nys,w2,m,w4)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),w2(1)
        complex *16 sigma(1),vals(1),w4(1)
c
c        allocate memory
c
        ipoln=1
        lpoln=nxs+4
c
        ipolnp1=ipoln+lpoln
        lpolnp1=nxs+4
c
        ipolyn=ipolnp1+lpolnp1
        lpolyn=nys+4
c
        ipolynp1=ipolyn+lpolyn
        lpolynp1=nys+4
c
        ivalsnp1=1
        lvalsnp1=nys+4
c
        ietan=ivalsnp1+lvalsnp1
        letan=nxs+4
c
        ietanp1=ietan+letan
        letanp1=nxs+4
c
c        . . . evaluate
c
        call pnmproj_eval00c(w,sigma,vals,nxs,nys,
     1        w4(ietan),w4(ietanp1),w2(ipoln),w2(ipolnp1),
     2        w2(ipolyn),w2(ipolynp1),w4(ivalsnp1))
c
        return
        end
c
c
c
c
c
        subroutine pnmproj_eval00c(w,sigma,vals,n,m,
     1      etan,etanp1,poln,polnp1,polyn,polynp1,valsnp1)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),poln(1),polnp1(1),polyn(1),polynp1(1)
        complex *16 sigma(1),vals(1),valsnp1(1),etan(1),etanp1(1)
     1      
c
c        construct the functions to be FMM'ed
c
        do 2200 i=1,n
c
        etan(i)=poln(i)*sigma(i)
        etanp1(i)=polnp1(i)*sigma(i)
 2200 continue
c
c        FMM them things
c
        call onedfmm_eval2c(w,etan,etanp1,vals,valsnp1,n,m)
c
c        . . . recombine
c
        do 2600 i=1,m
c
        vals(i)=vals(i)*polynp1(i)-valsnp1(i)*polyn(i)
 2600 continue
c
        return
        end
c
c
c
c
c
        subroutine pnmproj_eval(m,sigma,vals,w)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1)
        real *8 sigma(1),vals(1)
c
c        This subroutine applies the sphrical filter to the 
c        user-specified function sigma. Specifically, given the 
c        function 
c
c          sigma(x)=\sum{k=m}^N \alpha_k * P_k^m (x)            (1)
c
c       tabulated at the nodes xs(1),xs(2),...,xs(nxs),
c       the subroutine pnmproj_eval produces the function
c
c          sigma(x)=\sum{k=m}^{nlege} \alpha_k * P_k^m (x),     (2)
c
c       tabulated at the points ys(1),ys(2),...,ys(nys).
c
c       PLEASE NOTE THAT THIS SUBROUTINE USES THE ARRAY W THAT
C       (hopefully) has been created by a preceding call to
c       the subroutine pnmproj_init (see). This subroutine
c       has no uses known to man as a stand-alone device.       
c
c       POTENTIALLY IMPORTANT: The value of the upper limit
c       N in the summation in (1) above depends on the nature
c       of the nodes xs and the corresponding weights wsxs. 
c       If xs are the Gaussian nodes and wsxs are the 
c       corresponding Gaussian weights, then 
c
c                    N=nxs.                                     (3)
c
c       Any other selection of nodes and weights leads to a
c       smaller N.
c
c                        Input parameters:
c
c  m - the parameters in (1) above
c  sigma - the function being filtered, tabulated at the nxs
c        nodes 
c  w - the big array generated by a preceding call to the
c       subroutine pnmproj_init (see). 
c
c                         Output parameters:
c
c  vals - the filtered version of the function sigma, tabulated
c       at the nys nodes ys (see subroutine pnmproj_init )
c
c        . . . allocate memory
c
        iw=w(1)
        iiis=w(3)
        iw2=w(2)
c
        nxs=w(4)
        nys=w(5)
c
        ii=w(iiis+m)
        ii=ii+iw2-1
c
        iw3=w(6)

c
c        apply the projection operator
c
        call pnmproj_eval0(w(iw),sigma,vals,nxs,nys,w(ii),m,w(iw3) )
c
        return
        end
c
c
c
c
c
        subroutine pnmproj_eval0(w,sigma,vals,
     1      nxs,nys,w2,m,w4)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),w2(1)
        real *8 sigma(1),vals(1),w4(1)
c
c        allocate memory
c
        ipoln=1
        lpoln=nxs+4
c
        ipolnp1=ipoln+lpoln
        lpolnp1=nxs+4
c
        ipolyn=ipolnp1+lpolnp1
        lpolyn=nys+4
c
        ipolynp1=ipolyn+lpolyn
        lpolynp1=nys+4
c
        ivalsnp1=ipolynp1+lpolynp1
        ivalsnp1=1
        lvalsnp1=nys+4
c
        ietan=ivalsnp1+lvalsnp1
        letan=nxs+4
c
        ietanp1=ietan+letan
        letanp1=nxs+4
c
c        . . . evaluate
c
        call pnmproj_eval00(w,sigma,vals,nxs,nys,
     1        w4(ietan),w4(ietanp1),w2(ipoln),w2(ipolnp1),
     2        w2(ipolyn),w2(ipolynp1),w4(ivalsnp1))
c
        return
        end
c
c
c
c
c
        subroutine pnmproj_eval00(w,sigma,vals,n,m,
     1      etan,etanp1,poln,polnp1,polyn,polynp1,valsnp1)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1)
        real *8 sigma(1),vals(1),valsnp1(1),etan(1),etanp1(1),
     1      poln(1),polnp1(1),polyn(1),polynp1(1)
c
c        construct the functions to be FMM'ed
c
        do 2200 i=1,n
c
        etan(i)=poln(i)*sigma(i)
        etanp1(i)=polnp1(i)*sigma(i)
 2200 continue
c
c        FMM them things
c
        call onedfmm_eval2(w,etan,etanp1,vals,valsnp1,n,m)
c
c        . . . recombine
c
        do 2600 i=1,m
c
        vals(i)=vals(i)*polynp1(i)-valsnp1(i)*polyn(i)
 2600 continue
c
        return
        end
c
c
c
c
c
        subroutine pnmproj_init(ier,ndigits,nxs,nys,nlege,
     1      xs,wsxs,ys,w,lenw,keep)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),ys(1),w(1),wsxs(1)
c
c        This subroutine initializes the application of the 
c        spherical filter to user-specified functions; the latter
c        procedure is performed by the subroutine pnmproj_eval
c        (see). Specifically, given the function 
c
c          sigma(x)=\sum{k=m}^N \alpha_k * P_k^m (x)            (1)
c
c       tabulated at the nodes xs(1),xs(2),...,xs(nxs),
c       the subroutine pnmproj_eval produces the function
c
c          sigma(x)=\sum{k=m}^{nlege} \alpha_k * P_k^m (x),     (2)
c
c       tabulated at the points ys(1),ys(2),...,ys(nys).
c
c       POTENTIALLY IMPORTANT: The value of the upper limit
c       N in the summation in (1) above depends on the nature
c       of the nodes xs and the corresponding weights wsxs. 
c       If xs are the Gaussian nodes and wsxs are the 
c       corresponding Gaussian weights, then 
c
c                    N=nxs.                                     (3)
c
c       Any other selection of nodes and weights leads to a
c       smaller N.
c
c                        Input parameters:
c
c  ndigits - the number of digits with which the calculations
c       will be performed. Permitted values: 2-14.
c  nxs - the number of nodes in the original discretization 
c       of the interval [-1,1]; alos, the number of elements 
c       in each of the arrays xs,wsxs
c  nys - the number nodes at which the filterded functions 
c       will be evaluated. Also, the number of nodes in array ys
c  nlege - the highest order of P_n^m that is retained by the
c        filter
c  xs - the nodes at which the function to be filteredwill be
c        tabulated
c  wsxs - the quadrature weights at the nodes in xs
c  ys - the nodes at which the filtered function will be 
c        evaluated
c  lenw - the length of the user-supplied array w, in real *8 
c        words
c  
c                         Output parameters:
c
c  ier - error return code. 
c     ier=0 means successful conclusion
c     ier=16 means that the amount of storage in the user-supplied
c        array w is insufficient
c     ier=32 means that the amount of storage in the user-supplied
c        array w is grossly insufficient
c  w - the array containing various sorts of data to be used bt the
c        subroutine pnmproj_eval (see)
c  keep - the piece of array w that should not be changed between 
c        the call to this subroutine and the subsequent calls to 
c        the subroutine pnmproj_eval (see)
c
        ier=0
c
        itype=1
        call legeexps(itype,nxs,xs,u,v,wsxs)
c
        itype=0
        call legeexps(itype,nys,ys,u,v,wsys)
c
        iiis=21
        liis=nlege+3
c
        iw=iiis+liis
c
c        initialize the one-dimensional FMM
c
        lenw2=lenw-iw-1
        call onedfmm_init(jer,xs,nxs,ys,nys,ndigits,
     1      w(iw),lenw2,keep)
c
        if(jer .ne. 0) then
            ier=32
            return
        endif
c
c        . . . initialize the evaluation of P_n^m's
c
        iw2=iw+keep+1
        lw2=(nlege+1)*(nxs+nys+10)*2
c
        iw3=iw2+lw2
        lw3=6*nlege+200
c
        lw33=4*nxs+3*nys+20
        if(lw33 .gt. lw3) lw3=lw33
c
        lw3=lw3*2
c
        lused=iw3+lw3
        if(lused .gt. lenw) then
            ier=16
            return
        endif
c
        call pnmini(nlege+10,w(iw3) )
c
c        construct projection operators for all 
c        values of mm
c
        ii=1
        do 2400 mm7=0,nlege
c        
        w(iiis+mm7)=ii+0.1   
c
        call pnmproj_init0(w(iw),nlege,
     1      xs,wsxs,nxs,ys,nys,w(iw2+ii-1),keep,mm7,w(iw3))
c
        ii=ii+keep
 2400 continue
c
        w(1)=iw+0.1
        w(2)=iw2+0.1
        w(3)=iiis+0.1
        w(4)=nxs+0.1
        w(5)=nys+0.1
        w(6)=iw3+0.1
c
        keep=iw2+lw2
        keep=lused
c
        return
        end
c
c
c
c
c
        subroutine pnmproj_init0(w,nlege,
     1      xs,wsxs,nxs,ys,nys,w2,keepw2,m,w3)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),xs(1),ys(1),wsxs(1),w2(1),w3(1)
c
c        allocate memory
c
        ipoln=1
        lpoln=nxs+4
c
        ipolnp1=ipoln+lpoln
        lpolnp1=nxs+4
c
        ipolyn=ipolnp1+lpolnp1
        lpolyn=nys+4
c
        ipolynp1=ipolyn+lpolyn
        lpolynp1=nys+4
c
        keepw2=ipolynp1+lpolynp1
c
        call pnmproj_init00(w,nlege,xs,wsxs,nxs,ys,nys,
     1        w2(ipoln),w2(ipolnp1),w2(ipolyn),w2(ipolynp1),m,w3)
c
        return
        end
c
c
c
c
c
        subroutine pnmproj_init00(w,nlege,xs,wsxs,n,ys,m,
     1      poln,polnp1,polyn,polynp1,mm,w3)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),xs(1),ys(1),wsxs(1),poln(1),polnp1(1),
     1      polyn(1),polynp1(1),w3(1)
c
c        construct the values of P_n^{mm}(x), P_{n+1}^{mm}(x)
c
        done=1
        d=((nlege+1)**2-mm**2)/( 4*(nlege+done)**2-done)
        d=sqrt(d)
c
        do 1200 i=1,n
c
        call pnmeva_single(ier,xs(i),nlege,mm,
     1      w3,rsn,rsnp1)
c
        poln(i)=rsn*wsxs(i)*d
c
        polnp1(i)=rsnp1*wsxs(i)*d
c
 1200 continue
c
        do 1400 i=1,m
c
        call pnmeva_single(ier,ys(i),nlege,mm,
     1      w3,rsn,rsnp1)

        polyn(i)=rsn
        polynp1(i)=rsnp1
 1400 continue
c
        return
        end
c
c 
c 
c 
c 
c 
        subroutine pnmini(nmax,w)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1)
c 
c       this is the initialization subroutine for the subroutine
c       pnmeva (see) evaluating a collection of associated Legendre
c       functions at a point on the interval (-1,1).
c 
c                  input parameters:
c 
c  nmax - the highest degree of the Legendre polynomial that the
c          subroutine will be expected to evaluate
c 
c                  output parametrs:
c 
c  w - array containing various types of information to be used by
c          the subroutine pnmeva; the subroutine uses 6*nmax + 50
c          real *8 elements; these should not be changed between the
c          call to this subroutine and subsequent calls to the subroutine
c          pnmeva (see)
c 
c 
c        . . . allocate memory for the subroutine pnmini0
c 
        icmm=12
        lcmm=nmax+6
c 
        icmmp1=icmm+lcmm
        lcmmp1=nmax+6
c 
        iarrints=icmmp1+lcmmp1
        larrints=nmax+6
        larrints=larrints*2
c 
        iarrhalfs=iarrints+larrints
        larrhalfs=nmax+6
        larrhalfs=larrhalfs*2
c 
c        initialize the subroutine for the evaluation of
c        normalized assocated Legendre functions
c 
        call pnmini0(nmax+2,w(icmm),w(icmmp1),w(iarrints),
     1      w(iarrhalfs))
c 
        w(1)=icmm+0.1
        w(2)=icmmp1+0.1
        w(3)=iarrints+0.1
        w(4)=iarrhalfs+0.1
        w(5)=nmax+0.1
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine pnmeva(ier,x,n,m,rs,ders,w)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),rs(1),ders(1)
c 
c        this subroutine evaluates the functions R_k^m, with
c        k=m, m+1, . . . ,n; here R_k^m is the normalized version
c        of the associated Legendre function P_k^m. The functions
c        R_k^m are evaluated at the point x \in (-1,1) that is
c        supplied by the user. The subroutine also returns the
c        derivatives (with respect to x) of the values returned
c        in array rs
c 
c 
c                    Input parameters:
c 
c  x - the point where the functions R_k^m are to be evaluated
c  n - the largest value of k for which the functions
c        R_k^m are to be evaluated.
c  m - the value of the second parameter for which the functions
c        R_k^m are to be evaluated.
c  w - the array produced by the initializaion subroutine pnmini(see)
c 
c                    Output parameters:
c 
c  ier - error return code;
c              ier=0 means successful execution of the subroutine
c              ier=4 means that n > m; this is a fatal error
c              ier=16 means that n > nmax specified in the entry
c              pnmini; this is a fatal error
c  rs -  values at the point x of normalized associated Legendre
c        functions P_k^m, with k=m,m+1,. . . ,n. Note that the value
c        of the k-th Legendre function is found in the k-th element
c        of the array rs; thus, the first m-1 elements of rs are not
c        changed by the subroutine.
c  ders- the derivatives at the points x of the values rs
c 
c        . . . retrieve the locations of various arrays in w from
c              the beginning of w
c 
        icmm=w(1)
        icmmp1=w(2)
        iarrints=w(3)
        iarrhalfs=w(4)
        nmax=w(5)
c 
c        evaluate the requested normalized P_n^m with the
c        user-specified m
c 
        ier=0
        if(n .le. nmax) goto 1200
        ier=16
        return
 1200 continue
c 
        if(n .ge. m) goto 1400
        ier=4
        return
 1400 continue
c 
        call pnmeva0(x,m,n,w(icmm),w(icmmp1),w(iarrints),
     1      w(iarrhalfs),rs,ders)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine pnmini0(nmax,cmm,cmmp1,arrints,arrhalfs)
        implicit real *8 (a-h,o-z)
        save
        dimension cmm(0:1),cmmp1(0:1),arrints(0:1),
     1    arrhalfs(0:1)
c 
c       construct the arrays cmm, cmmp1
c 
        done=1
        half=done/2
c 
        cmm(0)=1
        cmmp1(0)=dsqrt(half)
c 
        do 1400 i=1,nmax
c 
        cmm(i)=-cmm(i-1)*dsqrt( (2*i-done)*2*i)/(i*2)
c 
        cmmp1(i)=
     1   -cmmp1(i-1)*dsqrt( (2*i+done)*(2*i+2) )/(i+1)/2
c 
 1400 continue
c 
        do 1600 i=0,nmax
c 
        cmm(i)=cmm(i)*dsqrt(i+half)
c 
        cmmp1(i)=cmmp1(i)*dsqrt( (2*i+2)*(i+3*half) )
 1600 continue
c 
c       construct the arrays of square roots of
c       integers and half-integers
c 
        do 1800 i=0,nmax*2
c 
        arrints(i)=dsqrt(i*done)
c 
        if(i .eq. 0) goto 1800
        arrhalfs(i)=dsqrt(i-half)
 1800 continue
c 
        return
        end
c 
c 
c 
c 
c 
  
        subroutine pnmeva0(x,m,n,cmm,cmmp1,arrints,
     1      arrhalfs,rs,ders)
        implicit real *8 (a-h,o-z)
        save
        dimension cmm(0:1),cmmp1(0:1),arrints(0:1),
     1    arrhalfs(0:1),rs(1),ders(1)
c 
c       calculate p_m^m and p_{m+1}^m
c 
        done=1
        half=done/2
c 
        rmm=(1-x**2)**(m*half) * cmm(m)
c 
        rmmp1=(1-x**2)**(m*half)*x*cmmp1(m)
c 
        dmm=-2*m*half*x * (1-x**2)**(m*half-1)
c 
        dmm=dmm*cmm(m)
c 
        dmmp1=(1-x**2)**(m*half) -2*m*half*x**2
     1     * (1-x**2)**(m*half-1)
c 
        dmmp1=dmmp1*cmmp1(m)
c 
c       conduct the whole recursion
c 
        rs(m)=rmm
        rs(m+1)=rmmp1
c 
        ders(m)=dmm
        ders(m+1)=dmmp1
c 
        do 2000 i=m+1,n
c 
        rat1=(2*i+1)* arrints(i-m+1)*arrhalfs(i+2) /
     1       (arrhalfs(i+1)*arrints(i+m+1)*(i-m+done) )
c 
        rat2=arrints(i+m) *arrints(i-m) *arrhalfs(i+2) /
     1      ( arrhalfs(i) * arrints(i-m+1)*arrints(i+m+1) )
c 
        rs(i+1)=rat1*x*rs(i)-rat2*rs(i-1)
c 
        ders(i+1)=rat1*(x*ders(i)+rs(i))-rat2*ders(i-1)
 2000 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine pnmeva_single(ier,x,n,m,w,rsn,rsnp1)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1)
c 
c        this subroutine evaluates the functions R_n^m, R_{n+1}^m, 
c        of the user-specified argument x. It is a castrated version
c        of the subroutine pnmeva (see)
c 
c                    Input parameters:
c 
c  x - the point where the functions R_k^m are to be evaluated
c  n - the largest value of k for which the functions
c        R_k^m are to be evaluated.
c  m - the value of the second parameter for which the functions
c        R_k^m are to be evaluated.
c  w - the array produced by the initializaion subroutine pnmini(see)
c 
c                    Output parameters:
c 
c  ier - error return code;
c              ier=0 means successful execution of the subroutine
c              ier=4 means that n > m; this is a fatal error
c              ier=16 means that n > nmax specified in the entry
c              pnmini; this is a fatal error
c  rsn -  value at the point x of normalized associated Legendre
c        functions R_n^m
c  rsnp1 -  value at the point x of normalized associated Legendre
c        functions R_{n+1}^m
c 
c        . . . retrieve the locations of various arrays in w from
c              the beginning of w
c 
        icmm=w(1)
        icmmp1=w(2)
        iarrints=w(3)
        iarrhalfs=w(4)
        nmax=w(5)
c 
c        evaluate the requested normalized P_n^m with the
c        user-specified m
c 
        ier=0
        if(n .le. nmax) goto 1200
        ier=16
        return
 1200 continue
c 
        if(n .ge. m) goto 1400
        ier=4
        return
 1400 continue
c 
        call pnmeva0_single(x,m,n,w(icmm),w(icmmp1),w(iarrints),
     1      w(iarrhalfs),rsn,rsnp1)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine pnmeva0_single(x,m,n,cmm,cmmp1,arrints,
     1      arrhalfs,rsn,rsnp1)
        implicit real *8 (a-h,o-z)
        save
        dimension cmm(0:1),cmmp1(0:1),arrints(0:1),
     1    arrhalfs(0:1)
c 
c       calculate p_m^m and p_{m+1}^m
c 
        done=1
        half=done/2
c 
        rmm=(1-x**2)**(m*half) * cmm(m)
c 
        rmmp1=(1-x**2)**(m*half)*x*cmmp1(m)
c 
        dmm=-2*m*half*x * (1-x**2)**(m*half-1)
c 
        dmm=dmm*cmm(m)
c 
        dmmp1=(1-x**2)**(m*half) -2*m*half*x**2
     1     * (1-x**2)**(m*half-1)
c 
        dmmp1=dmmp1*cmmp1(m)
c 
c       conduct the whole recursion
c 
        rsim1=rmm
        rsi=rmmp1
c 
        do 2000 i=m+1,n
c 
        rat1=(2*i+1)* arrints(i-m+1)*arrhalfs(i+2) /
     1       (arrhalfs(i+1)*arrints(i+m+1)*(i-m+done) )
c 
        rat2=arrints(i+m) *arrints(i-m) *arrhalfs(i+2) /
     1      ( arrhalfs(i) * arrints(i-m+1)*arrints(i+m+1) )
c 
        rsip1=rat1*x*rsi-rat2*rsim1
c
        rsim1=rsi
        rsi=rsip1
c
 2000 continue
c
        rsn=rsim1
        rsnp1=rsi
c 
        return
        end

