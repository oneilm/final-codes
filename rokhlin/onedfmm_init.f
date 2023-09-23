        implicit real *8 (a-h,o-z)
        dimension xs(10 000),ys(10 000),w(1000 000)
c
        complex *16 sigma(10 000),ima,
     1      vals(10000),vals2(10000),vals7(10000),
     2      diff(10000),sigma2(10000),diff2(10000),vals8(10000)
c
        data ima/(0.0d0,1.0d0)/

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
        itype=0
        call legeexps(itype,n,xs,u,v,whts)
c
        m=n/2
cccc        m=n-1
        call legeexps(itype,m,ys,u,v,whts)
c
        done=1
        do 1200 i=1,n
c
        sigma(i)=cos(done*i)+ima/i
        sigma2(i)=sin(sqrt(done*i))-ima/sqrt(done*i)
cccc        sigma2(i)=i
 1200 continue

cccc        call prin2('xs=*',xs,n)
cccc        call prin2('ys=*',ys,m)
cccc        call prin2('sigma=*',sigma,n*1)
cccc        call prin2('sigma2=*',sigma2,n*1)

        ndigits=13
        lenw=1000 000
        call onedfmm_init(ier,xs,n,ys,m,ndigits,w,lenw,keep)

        call prinf('before evaluation,n=*',n,1)

cc        PRINT *, 'ENTER dummy'
cc        READ *,dummy
c
C
        call onedfmm_eval2c(w,sigma,sigma2,vals,vals2,n,m)

        call prin2('after evaluation, vals=*',vals,m*2)
        call prin2('vals2=*',vals2,m*2)


        call opertest(xs,sigma,n,ys,m,vals7)
        call opertest(xs,sigma2,n,ys,m,vals8)


        call prin2('vals7=*',vals7,m)
        call prin2('vals=*',vals,m)

        dd=0
        ddd=0

        do 2200 i=1,m
c
        diff(i)=vals(i)-vals7(i)


        dd=dd+abs(diff(i))**2
        ddd=ddd+abs(vals(i))**2

        diff2(i)=vals2(i)-vals8(i)
 2200 continue
c
        call prin2('and diff=*',diff,m)
        call prin2('and diff2=*',diff2,m)
        call prinf('and keep=*',keep,1)
c
        errrel=sqrt(dd/ddd)
        call prin2('and errrel=*',errrel,1)
        stop
        end
c
c
c
c
c
        subroutine opertest(xs,sigma,n,ys,m,vals)
        implicit real *8 (a-h,o-z)
        save
        complex *16  sigma(1),vals(1),d
        dimension xs(1),ys(1)
c
        do 1400 i=1,m
        d=0
        do 1200 j=1,n
c
cccc        if(xs(j) .ge. ys(i)) goto 1200
c
        d=d-sigma(j)/(xs(j)-ys(i))
 1200 continue
c
        vals(i)=d
 1400 continue
c
        return
        end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        This is the end of the debugging code, and the beginning
c        of the one-dimensional FMM code proper.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c        This file contains five user-callable subroutines: onedfmm_init,
c        onedfmm_eval, onedfmm_eval2, onedfmm_evalc, onedfmm_eval2c. 
c        An adventurous user might also consider callable the entries 
c        onedfmm_eval_oneside, onedfmm_eval_oneside2, onedfmm_eval_onesidec,
c        onedfmm_eval_oneside2c of thw subroutine onedfmm_init.
c        Following is a brief description of the first five of these
c        subroutines.
c
c   onedfmm_init - prepares the various data for the rapid application 
c        of the one-sided FMM on the line, i.e. of the sums of the form
c
c           \sum_{j=1}^n_i  sigma(j)/(ys(i)-xs(j)),               (1)
c
c        for all j .leq. n_i. The actual evaluation of the sums
c        (1) is performed by the entries onedfmm_eval_oneside,
c        onedfmm_eval_oneside2, onedfmm_eval_onesidec,
c        onedfmm_eval_oneside2c of this subroutine (see below).
c
c   onedfmm_eval - applies the FMM operator on the line to the 
c        user-provided real vector. In other words, it evaluates 
c        the sum (1) for all pairs (xs(j),ys(i)) such that 
c        
c            xs(j) \neq ys(i).                                    (2)
c
c   onedfmm_eval2 - applies the FMM operator on the line to TWO 
c        user-provided real vectors. In other words,
c        it evaluates the sums
c
c           \sum_{j=1}^n sigma(j)/(ys(i)-xs(j)),                  (3)
c 
c
c           \sum_{j=1}^n sigma2(j)/(ys(i)-xs(j)),                 (4)
c 
c        subject to the condition (2)
c   
c   onedfmm_evalc - the same as onedfmm_eval, but with complex sigma
c   onedfmm_eval2c - the same as onedfmm_evalc, but with complex 
c        sigma, sigma2
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
        subroutine onedfmm_eval(w,sigma,vals,n,m)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),sigma(1),vals(1)
c
c        This subroutine applies the FMM operator on the line
c        to the  user-provided real vector. In other words,
c        it evaluates the sums
c
c
c           \sum_{j=1}^n sigma(j)/(ys(i)-xs(j)),                (1)
c 
c        omitting the pairs (i,j) such that 
c
c           xs(j)= ys(i)                                        (2)
c        
c        Please note that the subroutine oned_fmm_init (see)
c        should be called prior to the call to this subroutine.
c
c        IMPORTANT!!! THIS SUBROUTINE ASSUMES THAT THE ARRAYS 
C        XS, YS (SEE SUBROUTINE ONED_FMM_INIT FOR THE DEFINITION
C        OF THESE) SATISFY THE FOLLOWING CONDITIONS:
C
C         XS(I)=XS(N-I+1)                                       (3)
C
C        FOR ALL 1 \LEQ I \LEQ N,
C
C        AND
C
C         YS(I)=YS(M-I+1)                                       (4)
C
C        FOR ALL 1 \LEQ I \LEQ M.
C
c
c                     Input parameters:
c 
c  w - array containing various types of data, created by a preceding
c        call to the subroutine onedfmm_init (see)
c  sigma - the "charge" density in (1) above
c  n - the number of charges
c  m - the number of points where the sum (1) is to be
c      evaluated
c  
c                     Output parameters:
c
c  vals - the sum (1) devaluated at the m nodes ys(i)
C
c
c        . . . allocate memory
c
        isigmarev=w(9)
        ivalsrev=w(11)
c
c       . . . apply the FMM
c
        call onedfmm_eval10(w,sigma,vals,n,m,w(isigmarev),
     1      w(ivalsrev) )
c
        return
        end
c
c
c
c
c
        subroutine onedfmm_eval2(w,sigma,sigma2,vals,vals2,n,m)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),sigma(1),vals(1),sigma2(1),vals2(1)
c
c        This subroutine applies the FMM operator on the line
c        to TWO user-provided real vectors. In other words,
c        it evaluates the sums
c
c
c           \sum_{j=1}^n sigma(j)/(ys(i)-xs(j)),                (1)
c 
c
c           \sum_{j=1}^n sigma2(j)/(ys(i)-xs(j)),               (2)
c 
c        omitting the pairs (i,j) such that 
c
c           xs(j)= ys(i)                                        (3)
c        
c        Please note that the subroutine oned_fmm_init (see)
c        should be called prior to the call to this subroutine.
c
c        IMPORTANT!!! THIS SUBROUTINE ASSUMES THAT THE ARRAYS 
C        XS, YS (SEE SUBROUTINE ONED_FMM_INIT FOR THE DEFINITION
C        OF THESE) SATISFY THE FOLLOWING CONDITIONS:
C
C         XS(I)=XS(N-I+1)                                       (4)
C
C        FOR ALL 1 \LEQ I \LEQ N,
C
C        AND
C
C         YS(I)=YS(M-I+1)                                       (5)
C
C        FOR ALL 1 \LEQ I \LEQ M.
C
c
c                     Input parameters:
c 
c  w - array containing various types of data, created by a preceding
c        call to the subroutine onedfmm_init (see)
c  sigma, sigma2 - the "charge" densities in (1), (2) above
c  n - the number of charges
c  m - the number of points where the sums (1), (2) are to be
c      evaluated
c  
c                     Output parameters:
c
c  vals - the sum (1) devaluated at the m nodes ys(i)
c  vals2 - the sum (2) devaluated at the m nodes ys(i)
c
C       
c        . . . allocate memory
c
        isigmarev=w(9)
        isigma2rev=w(10)
        ivalsrev=w(11)
        ivals2rev=w(12)
c
c       . . . apply the FMM
c
        call onedfmm_eval20(w,sigma,sigma2,
     1      vals,vals2,n,m,w(isigmarev),w(isigma2rev),
     2      w(ivalsrev),w(ivals2rev) )
c
        return
        end
c
c
c
c
c
        subroutine onedfmm_evalc(w,sigma,vals,n,m)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1)
        complex *16 sigma(1),vals(1)
c
c        This subroutine applies the FMM operator on the line
c        to the  user-provided complex vector. In other words,
c        it evaluates the sums
c
c
c           \sum_{j=1}^n sigma(j)/(ys(i)-xs(j)),                (1)
c 
c        omitting the pairs (i,j) such that 
c
c           xs(j)= ys(i)                                        (2)
c        
c        Please note that the subroutine oned_fmm_init (see)
c        should be called prior to the call to this subroutine.
c
c        IMPORTANT!!! THIS SUBROUTINE ASSUMES THAT THE ARRAYS 
C        XS, YS (SEE SUBROUTINE ONED_FMM_INIT FOR THE DEFINITION
C        OF THESE) SATISFY THE FOLLOWING CONDITIONS:
C
C         XS(I)=XS(N-I+1)                                       (3)
C
C        FOR ALL 1 \LEQ I \LEQ N,
C
C        AND
C
C         YS(I)=YS(M-I+1)                                       (4)
C
C        FOR ALL 1 \LEQ I \LEQ M.
C
c
c                     Input parameters:
c 
c  w - array containing various types of data, created by a preceding
c        call to the subroutine onedfmm_init (see)
c  sigma - the "charge" density in (1) above
c  n - the number of charges
c  m - the number of points where the sum (1) is to be
c      evaluated
c  
c                     Output parameters:
c
c  vals - the sum (1) devaluated at the m nodes ys(i)
c
c        . . . allocate memory
c
        isigmarev=w(9)
        ivalsrev=w(11)
c
c       . . . apply the FMM
c
        call onedfmm_eval10c(w,sigma,vals,n,m,w(isigmarev),
     1      w(ivalsrev) )
c
        return
        end
c
c
c
c
c
        subroutine onedfmm_eval2c(w,sigma,sigma2,vals,vals2,n,m)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1)
        complex *16 sigma(1),vals(1),sigma2(1),vals2(1)
c
c        This subroutine applies the FMM operator on the line
c        to TWO user-provided complex vectors. In other words,
c        it evaluates the sums
c
c
c           \sum_{j=1}^n sigma(j)/(ys(i)-xs(j)),                (1)
c 
c
c           \sum_{j=1}^n sigma2(j)/(ys(i)-xs(j)),               (2)
c 
c        omitting the pairs (i,j) such that 
c
c           xs(j)= ys(i)                                        (3)
c        
c        Please note that the subroutine oned_fmm_init (see)
c        should be called prior to the call to this subroutine.
c
c        IMPORTANT!!! THIS SUBROUTINE ASSUMES THAT THE ARRAYS 
C        XS, YS (SEE SUBROUTINE ONED_FMM_INIT FOR THE DEFINITION
C        OF THESE) SATISFY THE FOLLOWING CONDITIONS:
C
C         XS(I)=XS(N-I+1)                                       (4)
C
C        FOR ALL 1 \LEQ I \LEQ N,
C
C        AND
C
C         YS(I)=YS(M-I+1)                                       (5)
C
C        FOR ALL 1 \LEQ I \LEQ M.
c
c                     Input parameters:
c 
c  w - array containing various types of data, created by a preceding
c        call to the subroutine onedfmm_init (see)
c  sigma, sigma2 - the "charge" densities in (1), (2) above
c  n - the number of charges
c  m - the number of points where the sums (1), (2) are to be
c      evaluated
c  
c                     Output parameters:
c
c  vals - the sum (1) devaluated at the m nodes ys(i)
c  vals2 - the sum (2) devaluated at the m nodes ys(i)
c
c
c        . . . allocate memory
c
        isigmarev=w(9)
        isigma2rev=w(10)
        ivalsrev=w(11)
        ivals2rev=w(12)
c
c       . . . apply the FMM
c
        call onedfmm_eval20c(w,sigma,sigma2,
     1      vals,vals2,n,m,w(isigmarev),w(isigma2rev),
     2      w(ivalsrev),w(ivals2rev) )
c
        return
        end

c
c
c
c
c
        subroutine onedfmm_init(ier,xs,n,ys,m,ndigits,
     1      w,lenw,keep)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),ys(1),sigma(1),rlams(100),ws(100),
     1      vals(1),w(1),sigma2(1),sigma3(1),sigma4(1),
     2      vals2(1),vals3(1),vals4(1)
c
c        This subroutine prepares the various data for the 
c        rapid application of the one-sided FMM on the line, 
c        i.e. of the sums of the form
c
c           \sum_{j=1}^n_i  sigma(j)/(ys(i)-xs(j)),               (1)
c 
c        where n_i is the larges integer such that 
c
c           xs(j) < ys(i)                                         (2)
c
c        for all j .leq. n_i. The actual evaluation of the sums
c        (1) is performed by the entries onedfmm_eval_oneside,
c        onedfmm_eval_oneside2, onedfmm_eval_onesidec,
c        onedfmm_eval_oneside2c of this subroutine (see below).
c
c        IMPORTANT: PLEASE NOTE THAT THE ELEMENTS IN BOTH ARRAYS 
C        XS, YS MUST BE IN INCREASING ORDER (NOT NECESSARILY 
C        STRICTLY INCREASING)
C
c                Input parameters:
c
c  xs - the points on the line where the charges are located
c  n - the numer of elements in xs
c  ys - the points on the line where the fiels is to be evaluated
c  m - the numer of elements in ys
c  ndigits - the L^2 accuracy (the number of decimal digits) with
c        which the result is to be evaluated. 
c  lenw - the length of the user-supplied array w, in real *8 elements
c
c                  Output parameters:
c
c  ier - error return code (not imploemented at this time) -10.27.03
c  w - array containingvarious types of data to be used by subroutine
c        onedfmm_eval2 and its relatives for the evaluation of expressions
c        of the form (1)
c  keep - the number of real *8 elements in array w that should
c        not be changed between the call to this subroutine and 
c        subsequent calls to onedfmm_eval2 and its relatives; please
c        note that this is also the number of elements in w actually
c        used by this subroutine
c
c        . . . retrieve the nodes and weights of the quadrature
c              from the subroutine onedfmm_retr
c
        call onedfmm_retr(ndigits,rlams,ws,nrlams,h)
c
c       . . . evaluate the operator
c
        n7=n
        m7=m
c
c        allocate memory
c
        iiends=21
        liends=m*2+4
c
        iiaddr=iiends+liends
        liaddr=n*2+4
c
        icoef=iiaddr+liaddr
        lcoef=nrlams*2+2
c
        iccs=icoef+lcoef
        lccs=lenw-iccs-2
c
        w(1)=nrlams +0.1
        w(2)=n+0.1
        w(3)=m+0.1
c
        w(4)=iiends+0.1
        w(5)=iiaddr+0.1
        w(6)=icoef+0.1
        w(7)=iccs+0.1
c
        call oned_fmm_init(xs,h,n,ys,m,w(iiends),
     1      rlams,ws,nrlams,w(iccs),w(iiaddr),keepccs)
c
        keep=keepccs+iccs
c
        isigmarev=keep+2
        lsigmarev=2*n+4
c
        isigma2rev=isigmarev+lsigmarev
        lsigma2rev=2*n+4
c
        ivalsrev=isigma2rev+lsigma2rev
        lvalsrev=m*2+4
c
        ivals2rev=ivalsrev+lvalsrev
        lvals2rev=m*2+4
c
        keep=ivals2rev+lvals2rev

        w(8)=keep+0.1
        w(9)=isigmarev+0.1
        w(10)=isigma2rev+0.1
        w(11)=ivalsrev+0.1
        w(12)=ivals2rev+0.1
c
        return
c
c
c
c
        entry  onedfmm_eval_oneside(w,sigma,vals)
c
        call onedfmm_wdecode(w,nrlams,n7,m7,iiends,
     1      iiaddr,icoef,iccs)
c
        call oned_fmm_eval0(sigma,n7,m7,w(iiends),
     1      nrlams,vals,w(iccs),w(iiaddr),w(icoef) )
c
        return
c
c
c
c
        entry  onedfmm_eval_oneside2(w,sigma,sigma2,vals,vals2)
c
        call onedfmm_wdecode(w,nrlams,n7,m7,iiends,
     1      iiaddr,icoef,iccs)
c
        call oned_fmm_eval02(sigma,sigma2,n7,m7,w(iiends),
     1      nrlams,vals,vals2,w(iccs),w(iiaddr),w(icoef) )
c
        return
c
c
c
c
        entry  onedfmm_eval_oneside4(w,sigma,sigma2,sigma3,
     1      sigma4,vals,vals2,vals3,vals4)
c
        call onedfmm_wdecode(w,nrlams,n7,m7,iiends,
     1      iiaddr,icoef,iccs)
c
        call oned_fmm_eval04(sigma,sigma2,sigma3,sigma4,
     1      n7,m7,w(iiends),nrlams,vals,vals2,vals3,vals4,
     2      w(iccs),w(iiaddr),w(icoef) )
c
        return
c
c
c
c
        entry  onedfmm_eval_onesidec(w,sigma,vals)
c
        call onedfmm_wdecode(w,nrlams,n7,m7,iiends,
     1      iiaddr,icoef,iccs)
c
        call oned_fmm_eval0c(sigma,
     1      n7,m7,w(iiends),nrlams,vals,
     2      w(iccs),w(iiaddr),w(icoef) )
c
        return
c
c
c
c
        entry  onedfmm_eval_oneside2c(w,sigma,sigma2,
     1      vals,vals2)
c
        call onedfmm_wdecode(w,nrlams,n7,m7,iiends,
     1      iiaddr,icoef,iccs)
c
        call oned_fmm_eval02c(sigma,sigma2,
     1      n7,m7,w(iiends),nrlams,vals,vals2,
     2      w(iccs),w(iiaddr),w(icoef) )
c
        return
c
c
c
c
        entry  onedfmm_eval_oneside4c(w,sigma,sigma2,sigma3,
     1      sigma4,vals,vals2,vals3,vals4)
c
        call onedfmm_wdecode(w,nrlams,n7,m7,iiends,
     1      iiaddr,icoef,iccs)
c
        call oned_fmm_eval04c(sigma,sigma2,sigma3,sigma4,
     1      n7,m7,w(iiends),nrlams,vals,vals2,vals3,vals4,
     2      w(iccs),w(iiaddr),w(icoef) )
c
        return
        end
c
c
c
c
c
        subroutine onedfmm_eval10(w,sigma,
     1      vals,n,m,sigmarev,valsrev)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),sigma(1),vals(1),
     1      sigmarev(1),valsrev(1)
c
c       flip the user-supplied arrays sigma,sigma2
c
        do 1200 i=1,n
c
        sigmarev(i)=sigma(n-i+1)
 1200 continue
c
c        apply the one-sided FMM to both original and 
c        flipped arrays
c
        call  onedfmm_eval_oneside2(w,sigma,
     1      sigmarev,vals,valsrev)
c
c       put them things back together
c
        do 1400 i=1,m
c
        vals(i)=vals(i)-valsrev(m-i+1)
 1400 continue
c
        return
        end
c
c
c
c
c
        subroutine onedfmm_eval10c(w,sigma,
     1      vals,n,m,sigmarev,valsrev)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1)
        complex *16 sigma(1),vals(1),
     1      sigmarev(1),valsrev(1)
c
c       flip the user-supplied array sigma
c
        do 1200 i=1,n
c
        sigmarev(i)=sigma(n-i+1)
 1200 continue
c
c        apply the one-sided FMM to both original and 
c        flipped arrays
c
        call  onedfmm_eval_oneside2c(w,sigma,
     1      sigmarev,vals,valsrev)
c
c       put them things back together
c
        do 1400 i=1,m
c
        vals(i)=vals(i)-valsrev(m-i+1)
 1400 continue
c
        return
        end
c
c
c
c
c
        subroutine onedfmm_eval20c(w,sigma,sigma2,
     1      vals,vals2,n,m,sigmarev,sigma2rev,
     2      valsrev,vals2rev)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1)
        complex *16 sigma(1),vals(1),sigma2(1),
     1      vals2(1),sigmarev(1),sigma2rev(1),
     2      valsrev(1),vals2rev(1)
c
c       flip the user-supplied arrays sigma,sigma2
c
        do 1200 i=1,n
c
        sigmarev(i)=sigma(n-i+1)
        sigma2rev(i)=sigma2(n-i+1)
 1200 continue
c
c        apply the one-sided FMM to both original and 
c        flipped arrays
c
        call  onedfmm_eval_oneside4c(w,sigma,sigma2,
     1      sigmarev,sigma2rev,vals,vals2,valsrev,vals2rev)
c
c       put them things back together
c
        do 1400 i=1,m
c
        vals(i)=vals(i)-valsrev(m-i+1)
        vals2(i)=vals2(i)-vals2rev(m-i+1)
 1400 continue
c
        return
        end
c
c
c
c
c
        subroutine oned_fmm_eval04c(sigma,sigma2,sigma3,sigma4,
     1      n,m,iends,nlams,vals,vals2,vals3,vals4,ccs,iaddr,
     2      coef)
        implicit real *8 (a-h,o-z)
        save
        complex *16 sigma(n),vals(1),d,d2,d3,d4,
     1      coef(1),sigma2(n),vals2(1),coef2(100),
     2      sigma3(1),sigma4(1),coef3(100),coef4(100),vals3(1),
     3      vals4(1)
c
        dimension iends(2,1),iaddr(2,1),ccs(1)
c
c        construct the right-going exponential expansions
c        for all points in array xs
c
        do 1100 i=1,m
c
        vals(i)=0
        vals2(i)=0
        vals3(i)=0
        vals4(i)=0
 1100 continue
c
        ii=0
        do 1200 i=1,nlams
c
        coef(i)=sigma(1)
        coef2(i)=sigma2(1)
        coef3(i)=sigma3(1)
        coef4(i)=sigma4(1)
 1200 continue
c
        do 2000 i=1,n
c
c       if this point in x interacts with the y array,
c       account for this interaction
c
        if(iaddr(1,i) .eq. 0) goto 1800
c
        do 1600 j=iaddr(1,i),iaddr(2,i)
c
        d=0
        d2=0
        d3=0
        d4=0
        do 1500 jj=1,nlams
c
        ii=ii+1
        d=d+ccs(ii)*coef(jj)
        d2=d2+ccs(ii)*coef2(jj)
        d3=d3+ccs(ii)*coef3(jj)
        d4=d4+ccs(ii)*coef4(jj)
 1500 continue
c
        vals(j)=d                
        vals2(j)=d2                
        vals3(j)=d3
        vals4(j)=d4
c
 1600 continue
c
 1800 continue
c
        if(i .eq. n) goto 2000
c
        do 1900 j=1,nlams
c
        ii=ii+1
        coef(j)=coef(j)*ccs(ii)+sigma(i+1)
        coef2(j)=coef2(j)*ccs(ii)+sigma2(i+1)
        coef3(j)=coef3(j)*ccs(ii)+sigma3(i+1)
        coef4(j)=coef4(j)*ccs(ii)+sigma4(i+1)
c
 1900 continue
c
 2000 continue
c
c        add the part of the sum that is to be 
c        evaluated directly
c
        do 3400 i=1,m
c
        i1=iends(1,i)+1
        i2=iends(2,i)
c
        if(i2 .lt. i1) goto 3400
c
        d=0
        d2=0
        d3=0
        d4=0
        do 3200 j=i1,i2
c
        ii=ii+1
        d=d+sigma(j)*ccs(ii)
        d2=d2+sigma2(j)*ccs(ii)
        d3=d3+sigma3(j)*ccs(ii)
        d4=d4+sigma4(j)*ccs(ii)
 3200 continue
c
        vals(i)=vals(i)+d
        vals2(i)=vals2(i)+d2
        vals3(i)=vals3(i)+d3
        vals4(i)=vals4(i)+d4
 3400 continue
c
        return
        end
c
c
c
c
c
        subroutine onedfmm_eval20(w,sigma,sigma2,
     1      vals,vals2,n,m,sigmarev,sigma2rev,
     2      valsrev,vals2rev)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),sigma(1),vals(1),sigma2(1),
     1      vals2(1),sigmarev(1),sigma2rev(1),
     2      valsrev(1),vals2rev(1)
c
c       flip the user-supplied arrays sigma,sigma2
c
        do 1200 i=1,n
c
        sigmarev(i)=sigma(n-i+1)
        sigma2rev(i)=sigma2(n-i+1)
 1200 continue
c
c        apply the one-sided FMM to both original and 
c        flipped arrays
c
        call  onedfmm_eval_oneside4(w,sigma,sigma2,
     1      sigmarev,sigma2rev,vals,vals2,valsrev,vals2rev)
c
c       put them things back together
c
        do 1400 i=1,m
c
        vals(i)=vals(i)-valsrev(m-i+1)
        vals2(i)=vals2(i)-vals2rev(m-i+1)
 1400 continue
c
        return
        end
c
c
c
c
c
        subroutine oned_fmm_eval0c(sigma,
     1      n,m,iends,nlams,vals,ccs,iaddr,
     2      coef)
        implicit real *8 (a-h,o-z)
        save
        complex *16 sigma(n),vals(1),d,d2,
     1      coef(1)
c
        dimension iends(2,1),iaddr(2,1),ccs(1)
c
c        construct the right-going exponential expansions
c        for all points in array xs
c
        do 1100 i=1,m
c
        vals(i)=0
 1100 continue
c
        ii=0
        do 1200 i=1,nlams
c
        coef(i)=sigma(1)
 1200 continue
c
        do 2000 i=1,n
c
c       if this point in x interacts with the y array,
c       account for this interaction
c
        if(iaddr(1,i) .eq. 0) goto 1800
c
        do 1600 j=iaddr(1,i),iaddr(2,i)
c
        d=0
        d2=0
        do 1500 jj=1,nlams
c
        ii=ii+1
        d=d+ccs(ii)*coef(jj)
 1500 continue
c
        vals(j)=d                
c
 1600 continue
c
 1800 continue
c
        if(i .eq. n) goto 2000
c
        do 1900 j=1,nlams
c
        ii=ii+1
        coef(j)=coef(j)*ccs(ii)+sigma(i+1)
c
 1900 continue
c
 2000 continue
c
c        add the part of the sum that is to be 
c        evaluated directly
c
        do 3400 i=1,m
c
        i1=iends(1,i)+1
        i2=iends(2,i)
c
        if(i2 .lt. i1) goto 3400
c
        d=0
        d2=0
        do 3200 j=i1,i2
c
        ii=ii+1
        d=d+sigma(j)*ccs(ii)
 3200 continue
c
        vals(i)=vals(i)+d
 3400 continue
c
        return
        end
c
c
c
c
c
        subroutine oned_fmm_eval02c(sigma,sigma2,
     1      n,m,iends,nlams,vals,vals2,ccs,iaddr,
     2      coef)
        implicit real *8 (a-h,o-z)
        save
        complex *16 sigma(n),vals(1),d,d2,
     1      coef(1),sigma2(n),vals2(1),coef2(100)
c
        dimension iends(2,1),iaddr(2,1),ccs(1)
c
c        construct the right-going exponential expansions
c        for all points in array xs
c
        do 1100 i=1,m
c
        vals(i)=0
        vals2(i)=0
 1100 continue
c
        ii=0
        do 1200 i=1,nlams
c
        coef(i)=sigma(1)
        coef2(i)=sigma2(1)
 1200 continue

        call prin2('coef=*',coef,nlams*2)
c
        do 2000 i=1,n
c
c       if this point in x interacts with the y array,
c       account for this interaction
c
        if(iaddr(1,i) .eq. 0) goto 1800
c
        do 1600 j=iaddr(1,i),iaddr(2,i)
c
        d=0
        d2=0
        do 1500 jj=1,nlams
c
        ii=ii+1
        d=d+ccs(ii)*coef(jj)
        d2=d2+ccs(ii)*coef2(jj)
 1500 continue
c
        vals(j)=d                
        vals2(j)=d2                
c
 1600 continue
c
 1800 continue
c
        if(i .eq. n) goto 2000
c
        do 1900 j=1,nlams
c
        ii=ii+1
        coef(j)=coef(j)*ccs(ii)+sigma(i+1)
        coef2(j)=coef2(j)*ccs(ii)+sigma2(i+1)
c
 1900 continue
c
 2000 continue
c
c        add the part of the sum that is to be 
c        evaluated directly
c
        do 3400 i=1,m
c
        i1=iends(1,i)+1
        i2=iends(2,i)
c
        if(i2 .lt. i1) goto 3400
c
        d=0
        d2=0
        do 3200 j=i1,i2
c
        ii=ii+1
        d=d+sigma(j)*ccs(ii)
        d2=d2+sigma2(j)*ccs(ii)
 3200 continue
c
        vals(i)=vals(i)+d
        vals2(i)=vals2(i)+d2
 3400 continue

        call prin2('vals=*',vals,m*2)

c
        return
        end
c
c
c
c
c
        subroutine oned_fmm_eval04(sigma,sigma2,sigma3,sigma4,
     1      n,m,iends,nlams,vals,vals2,vals3,vals4,ccs,iaddr,
     2      coef)
        implicit real *8 (a-h,o-z)
        save
        dimension sigma(n),iends(2,1),vals(1),ccs(1),
     1      iaddr(2,1),coef(1),sigma2(n),vals2(1),coef2(100),
     2      sigma3(1),sigma4(100),coef3(100),coef4(1),vals3(1),
     3      vals4(1)
c
c        construct the right-going exponential expansions
c        for all points in array xs
c
        do 1100 i=1,m
c
        vals(i)=0
        vals2(i)=0
        vals3(i)=0
        vals4(i)=0
 1100 continue
c
        ii=0
        do 1200 i=1,nlams
c
        coef(i)=sigma(1)
        coef2(i)=sigma2(1)
        coef3(i)=sigma3(1)
        coef4(i)=sigma4(1)
 1200 continue
c
        do 2000 i=1,n
c
c       if this point in x interacts with the y array,
c       account for this interaction
c

        if(iaddr(1,i) .eq. 0) goto 1800
c
        do 1600 j=iaddr(1,i),iaddr(2,i)
c
        d=0
        d2=0
        d3=0
        d4=0
        do 1500 jj=1,nlams
c
        ii=ii+1
        d=d+ccs(ii)*coef(jj)
        d2=d2+ccs(ii)*coef2(jj)
        d3=d3+ccs(ii)*coef3(jj)
        d4=d4+ccs(ii)*coef4(jj)
 1500 continue
c
        vals(j)=d                
        vals2(j)=d2                
        vals3(j)=d3
        vals4(j)=d4
c
 1600 continue
c
 1800 continue
c
        if(i .eq. n) goto 2000
c
        do 1900 j=1,nlams
c
        ii=ii+1
        coef(j)=coef(j)*ccs(ii)+sigma(i+1)
        coef2(j)=coef2(j)*ccs(ii)+sigma2(i+1)
        coef3(j)=coef3(j)*ccs(ii)+sigma3(i+1)
        coef4(j)=coef4(j)*ccs(ii)+sigma4(i+1)
c
 1900 continue
c
 2000 continue
c
c        add the part of the sum that is to be 
c        evaluated directly
c
        do 3400 i=1,m
c
        i1=iends(1,i)+1
        i2=iends(2,i)
c
        if(i2 .lt. i1) goto 3400
c
        d=0
        d2=0
        d3=0
        d4=0
        do 3200 j=i1,i2
c
        ii=ii+1
        d=d+sigma(j)*ccs(ii)
        d2=d2+sigma2(j)*ccs(ii)
        d3=d3+sigma3(j)*ccs(ii)
        d4=d4+sigma4(j)*ccs(ii)
 3200 continue
c
        vals(i)=vals(i)+d
        vals2(i)=vals2(i)+d2
        vals3(i)=vals3(i)+d3
        vals4(i)=vals4(i)+d4
 3400 continue
c
        return
        end
c
c
c
c
c
        subroutine oned_fmm_eval02(sigma,sigma2,n,m,iends,
     1      nlams,vals,vals2,ccs,iaddr,coef)
        implicit real *8 (a-h,o-z)
        save
        dimension sigma(n),iends(2,1),vals(1),ccs(1),
     1      iaddr(2,1),coef(1),sigma2(n),vals2(1),coef2(100)
c
c        construct the right-going exponential expansions
c        for all points in array xs
c
        do 1100 i=1,m
c
        vals(i)=0
        vals2(i)=0
 1100 continue
c
        ii=0
        do 1200 i=1,nlams
        coef(i)=sigma(1)
        coef2(i)=sigma2(1)
 1200 continue
c
        do 2000 i=1,n
c
c       if this point in x interacts with the y array,
c       account for this interaction
c
        if(iaddr(1,i) .eq. 0) goto 1800
c
        do 1600 j=iaddr(1,i),iaddr(2,i)
c
        d=0
        d2=0
        do 1500 jj=1,nlams
c
        ii=ii+1
        d=d+ccs(ii)*coef(jj)
        d2=d2+ccs(ii)*coef2(jj)
 1500 continue
c
        vals(j)=d                
        vals2(j)=d2                
c
 1600 continue
c
 1800 continue
c
        if(i .eq. n) goto 2000
c
        do 1900 j=1,nlams
c
        ii=ii+1
        coef(j)=coef(j)*ccs(ii)+sigma(i+1)
        coef2(j)=coef2(j)*ccs(ii)+sigma2(i+1)
c
 1900 continue
c
 2000 continue
c
c        add the part of the sum that is to be 
c        evaluated directly
c
        do 3400 i=1,m
c
        i1=iends(1,i)+1
        i2=iends(2,i)
c
        if(i2 .lt. i1) goto 3400
c
        d=0
        d2=0
        do 3200 j=i1,i2
c
        ii=ii+1
        d=d+sigma(j)*ccs(ii)
        d2=d2+sigma2(j)*ccs(ii)
 3200 continue
c
        vals(i)=vals(i)+d
        vals2(i)=vals2(i)+d2
 3400 continue
c
        return
        end
c
c
c
c
c
        subroutine oned_fmm_eval0(sigma,n,m,iends,
     1      nlams,vals,ccs,iaddr,coef)
        implicit real *8 (a-h,o-z)
        save
        dimension sigma(n),iends(2,1),vals(1),ccs(1),
     1      iaddr(2,1),coef(1)
c
c        construct the right-going exponential expansions
c        for all points in array xs
c
        do 1100 i=1,m
c
        vals(i)=0
 1100 continue
c
        ii=0
        do 1200 i=1,nlams
        coef(i)=sigma(1)
 1200 continue
c
        do 2000 i=1,n
c
c       if this point in x interacts with the y array,
c       account for this interaction
c
        if(iaddr(1,i) .eq. 0) goto 1800
c
        do 1600 j=iaddr(1,i),iaddr(2,i)
c
        d=0
        do 1500 jj=1,nlams
c
        ii=ii+1
        d=d+ccs(ii)*coef(jj)
 1500 continue
c
        vals(j)=d                
c
 1600 continue
c
 1800 continue
c
        if(i .eq. n) goto 2000
c
        do 1900 j=1,nlams
c
        ii=ii+1
        coef(j)=coef(j)*ccs(ii)+sigma(i+1)
c
 1900 continue
c
 2000 continue
c
c        add the part of the sum that is to be 
c        evaluated directly
c
        do 3400 i=1,m
c
        i1=iends(1,i)+1
        i2=iends(2,i)
c
        if(i2 .lt. i1) goto 3400
c
        d=0
        do 3200 j=i1,i2
c
        ii=ii+1
        d=d+sigma(j)*ccs(ii)
 3200 continue
c
        vals(i)=vals(i)+d
 3400 continue
c
        return
        end
c
c
c
c
c
        subroutine oned_fmm_init(xs,h,n,ys,m,iends,
     1      rlams,ws,nlams,ccs,iaddr,keepccs)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(n),ys(m),iends(2,1),
     1      rlams(1),ws(1),ccs(1),iaddr(2,1)
c
c       construct the control table
c
        call oned_fmm_bldgeo(xs,h,n,ys,m,iends)
        call oned_fmm_table0(iends,n,m,iaddr)
c
c        construct the right-going exponential expansions
c        for all points in array xs
c
        ii=0
        do 2000 i=1,n
c
c       if this point in x interacts with the y array,
c       account for this interaction
c
        if(iaddr(1,i) .eq. 0) goto 1800
c
        do 1600 j=iaddr(1,i),iaddr(2,i)
c
        i1=i
c        
        d=0
        do 1500 jj=1,nlams
c
        ii=ii+1
        ccs(ii)=exp(-(ys(j)-xs(i))*rlams(jj))*ws(jj)
 1500 continue
c
 1600 continue
c
 1800 continue
c
        if(i .eq. n) goto 2000
c
        do 1900 j=1,nlams
c
        ii=ii+1
        ccs(ii)=exp(-(xs(i+1)-xs(i)) *rlams(j) )
c
 1900 continue
c
 2000 continue

c        add the part of the sum that is to be 
c        evaluated directly
c
        do 3400 i=1,m
c
        i1=iends(1,i)+1
        i2=iends(2,i)
c
        if(i2 .lt. i1) goto 3400
c
        d=0
        do 3200 j=i1,i2
c
        ii=ii+1
        ccs(ii)=1/(ys(i)-xs(j))
 3200 continue
c
 3400 continue
c
        keepccs=ii
        return
        end
c
c
c
c
c
        subroutine oned_fmm_table0(iends,n,m,iaddr)
        implicit real *8 (a-h,o-z)
        save
        dimension iends(2,1),iaddr(2,1)
c
        do 1400 i=1,n
c
        iaddr(1,i)=0
        iaddr(2,i)=0
 1400 continue
c
c       for each point in x, construct the first and the last
c       points in y for which the expansion at the said point in 
c       x is to be used
c
        do 1600 i=1,m
c
        if(iends(1,i) .eq. 0) goto 1600
c
        j=iends(1,i)
c
        if(iaddr(1,j) .eq. 0) iaddr(1,j)=i
        iaddr(2,j)=i
 1600 continue
c
        return
        end
c
c
c
c
c
        subroutine oned_fmm_bldgeo(xs,h,n,ys,m,iends)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(n),ys(m),iends(2,m)
c
c        For the user-supplied monotonically increasing arrays 
c        xs, ys, and a positive real h, this subroutine returns 
c        the integer array iends(2,m). The element iends(1,i) 
c        contains the index j of the first element in xs such that
c
c                xs(j) \leq ys(i).
c
c        The element iends(2,i) contains the index j of the last 
c        element in x such that 
c
c                xs(j)+h < ys(i)
c        
c
c        . . . for each point in array ys, construct the 
c              control table
c
        do 2000 i=1,m
c
        iends(1,i)=0
        iends(2,i)=0
c
        do 1200 j=1,n
c
        if(xs(j)+h .le. ys(i)) iends(1,i)=j
        if(xs(j) .lt. ys(i)) iends(2,i)=j
c
 1200 continue
c
 2000 continue
c
        return
        end
c
c
c
c
c
        subroutine onedfmm_wdecode(w,nrlams,n7,m7,iiends,
     1      iiaddr,icoef,iccs)
        implicit real *8 (a-h,o-z)
        dimension w(1)
c
        nrlams=w(1)
        n7=w(2)
        m7=w(3)
c
        iiends=w(4)
        iiaddr=w(5)
        icoef=w(6)
        iccs=w(7)
c
        return
        end
c
c
c
c
c
        subroutine onedfmm_retr(ndigits,rl,ws,nrl,h)
        implicit real *8 (a-h,o-z)
        save
        dimension rl(13),ws(13),
     1      rl6(6),ws6(6),rl8(8),ws8(8),rl10(10),ws10(10),
     2      rl12(12),ws12(12),rl14(14),ws14(14),
     3      rl16(16),ws16(16),rl18(18),ws18(18),
     4      rl20(20),ws20(20),rl22(22),ws22(22),
     5      rl25(25),ws25(25),rl27(27),ws27(27),
     6      rl32(32),ws32(32),
     7      rl30(30),ws30(30)
c
        data rl6/
     1   .9910154337706231D+00, .9084071844138055D+01,
     2   .4747421717317850D+02, .1968789147108246D+03,
     3   .7122438107923539D+03, .2339728223490408D+04/
c        
        data ws6/
     1   .2907563180882740D+01, .1659805592378940D+02,
     2   .7198030139035295D+02, .2645679486948234D+03,
     3   .8778984157883630D+03, .2728700865513198D+04/
c
        data rl8/
     1   .6362917667647940D+00, .4531226927560700D+01,
     2   .1831440756876937D+02, .6076421269052143D+02,
     3   .1813865842334534D+03, .5044535357482889D+03,
     4   .1316550694135031D+04, .3332825423985843D+04/
c        
        data ws8/
     1   .1745629958618850D+01, .7031355038932648D+01,
     2   .2343196595700368D+02, .6889763047870008D+02,
     3   .1920250019153495D+03, .4982169872818122D+03,
     4   .1231249697523693D+04, .3143898056960176D+04/
c
        data rl10/
     1   .4752294642218484D+00, .2917512985686099D+01,
     2   .9655226623408764D+01, .2703534104519955D+02,
     3   .7003264926193296D+02, .1717850737661778D+03,
     4   .4029577529471144D+03, .9098896960860491D+03,
     5   .1995155959039442D+04, .4365610024629618D+04/
c        
        data ws10/
     1   .1261101956356473D+01, .3962406657099427D+01,
     2   .1050702541902272D+02, .2660342130786542D+02,
     3   .6459426609162857D+02, .1501011629872073D+03,
     4   .3354527820608259D+03, .7258987009994283D+03,
     5   .1548263458889049D+04, .3530738068442793D+04/
c
        data rl12/
     1   .3847353888490337D+00, .2230672677934408D+01,
     2   .6612754715381708D+01, .1625669949595771D+02,
     3   .3713775896583748D+02, .8136716415913352D+02,
     4   .1727457419682600D+03, .3570420754841424D+03,
     5   .7205048377333785D+03, .1423384054381622D+04,
     6   .2770079316040960D+04, .5451267042220382D+04/
c        
        data ws12/
     1   .1008371624995128D+01, .2840386996924615D+01,
     2   .6354124661571288D+01, .1388452163654476D+02,
     3   .2980797224763764D+02, .6245424842028146D+02,
     4   .1276365399042932D+03, .2548063999137857D+03,
     5   .4978244804903927D+03, .9559339532000802D+03,
     6   .1837897798437768D+04, .3861157114749360D+04/
c        
        data rl14/
     1   .3202872733708966D+00, .1797667042345338D+01,
     2   .5003626333533760D+01, .1128660450684209D+02,
     3   .2349743548212166D+02, .4709090509197091D+02,
     4   .9207035928182395D+02, .1763065915848103D+03,
     5   .3310449002798700D+03, .6105937927736252D+03,
     6   .1109800505866439D+04, .1994534311599102D+04,
     7   .3568390797019952D+04, .6518991935089340D+04/

        data ws14/
     1   .8336429572185068D+00, .2203188720993088D+01,
     2   .4422172642332030D+01, .8601795248031861D+01,
     3   .1669780400654837D+02, .3210263612364786D+02,
     4   .6074724676695964D+02, .1127808853412877D+03,
     5   .2054646573761684D+03, .3689637985633234D+03,
     6   .6561823721201911D+03, .1161075112309891D+04,
     7   .2086338299626188D+04, .4153478705227042D+04/
        data rl16/
     1   .2772538637975959D+00, .1529382828275024D+01,
     2   .4107906064038672D+01, .8765601880453531D+01,
     3   .1703447276037033D+02, .3173077009029402D+02,
     4   .5770704425858937D+02, .1031756851397945D+03,
     5   .1819822123094220D+03, .3171011711861000D+03,
     6   .5463449481991935D+03, .9316002695045712D+03,
     7   .1573990034726075D+04, .2641999855893926D+04,
     8   .4437364125250710D+04, .7643089690966770D+04/
c        
        data ws16/
     1   .7189324203990844D+00, .1835062458126015D+01,
     2   .3441737472221243D+01, .6120654948429199D+01,
     3   .1087474755836901D+02, .1930577798894369D+02,
     4   .3397757155777860D+02, .5921168118527332D+02,
     5   .1021386527216847D+03, .1742508098868775D+03,
     6   .2943144399190919D+03, .4926186856075870D+03,
     7   .8193322949627307D+03, .1364184608992976D+04,
     8   .2325176573295422D+04, .4426855212029718D+04/
c
        data rl18/
     1   .2433354395387177D+00, .1327301607790609D+01,
     2   .3485987505982816D+01, .7180435155425652D+01,
     3   .1331324915684459D+02, .2349061385261277D+02,
     4   .4039468889801780D+02, .6837203791767776D+02,
     5   .1143861425920434D+03, .1895066011075467D+03,
     6   .3112093402349025D+03, .5068946049675827D+03,
     7   .8192797608928415D+03, .1314821902706027D+04,
     8   .2097699808055681D+04, .3336558142685885D+04,
     9   .5330621250111263D+04, .8762701035750893D+04/
c        
        data ws18/
     1   .6294188274550986D+00, .1570784517951019D+01,
     2   .2820338185003186D+01, .4714694581603664D+01,
     3   .7812662388224537D+01, .1297970288039020D+02,
     4   .2153729629211060D+02, .3555391799671030D+02,
     5   .5828168445588339D+02, .9481481747356486D+02,
     6   .1530700520019215D+03, .2452845408305543D+03,
     7   .3903425206023612D+03, .6177040648292923D+03,
     8   .9752492525035492D+03, .1549483993863369D+04,
     9   .2536952798031864D+04, .4671277338439139D+04/
c
        data rl20/
     1  0.2170177283456096D+00,0.1174746269941844D+01,
     2  0.3038588014971929D+01,0.6106519063154826D+01,
     3  0.1094128682794228D+02,0.1852175882405450D+02,
     4  0.3044489361780420D+02,0.4921548645262547D+02,
     5  0.7869483353646393D+02,0.1247951233904915D+03,
     6  0.1965346659035602D+03,0.3076212311771787D+03,
     7  0.4787976066774670D+03,0.7412742780417782D+03,
     8  0.1141769161713824D+04,0.1750249651947201D+04,
     9  0.2673135693940089D+04,0.4079584646177636D+04,
     *  0.6267701180069594D+04,0.9920850864038967D+04/
c        
        data ws20/
     1  0.5603915249742358D+00,0.1377165299080497D+01,
     2  0.2398448193254439D+01,0.3828666617880930D+01,
     3  0.6000739403346859D+01,0.9421722661028075D+01,
     4  0.1483472202461460D+02,0.2333871269820330D+02,
     5  0.3658815589813275D+02,0.5708961663836639D+02,
     6  0.8863578237967738D+02,0.1369383568538534D+03,
     7  0.2105343746329764D+03,0.3220819438227263D+03,
     8  0.4903636271509675D+03,0.7439410711250566D+03,
     9  0.1129114610031477D+04,0.1730587791966919D+04,
     *  0.2743767448420053D+04,0.4910017945848916D+04/
c
        data rl22/
     1  0.1960560208802185D+00,0.1055782234290130D+01,
     2  0.2703249463235939D+01,0.5344984991832819D+01,
     3  0.9358368219827713D+01,0.1538369826608681D+02,
     4  0.2444563059786270D+02,0.3811899590983725D+02,
     5  0.5876985733720651D+02,0.8990509468408722D+02,
     6  0.1366783137535826D+03,0.2066397767568919D+03,
     7  0.3108490611196090D+03,0.4654600369683961D+03,
     8  0.6939326543791237D+03,0.1030240652061212D+04,
     9  0.1523701739546859D+04,0.2246421941912387D+04,
     *  0.3305936255598089D+04,0.4870854956721260D+04,
     1  0.7237416678335674D+04,0.1109050193236164D+05/
c        
        data ws22/
     1  0.5056719624275326D+00,0.1229779608855305D+01,
     2  0.2098532545879258D+01,0.3246289571928536D+01,
     3  0.4885205352673618D+01,0.7333549488551948D+01,
     4  0.1104884508210339D+02,0.1668649828522936D+02,
     5  0.2519148645428995D+02,0.3792573431209146D+02,
     6  0.5686077882216572D+02,0.8487927930498260D+02,
     7  0.1261955689076592D+03,0.1868845131180676D+03,
     8  0.2756443762742277D+03,0.4050721786581285D+03,
     9  0.5936720945791979D+03,0.8693102683326101D+03,
     *  0.1276912043149237D+04,0.1899361316169435D+04,
     1  0.2931852157985050D+04,0.5125655583569004D+04/
c        
        data rl25/
     1  0.1710023197238303D+00,0.9158520228462631D+00,
     2  0.2320353217523775D+01,0.4511763835308422D+01,
     3  0.7713017354924120D+01,0.1229048310619602D+02,
     4  0.1881746781234381D+02,0.2815310891204626D+02,
     5  0.4154497636805841D+02,0.6077392856106784D+02,
     6  0.8836294808886656D+02,0.1278725839711497D+03,
     7  0.1843151210609086D+03,0.2647290229979918D+03,
     8  0.3789734304010005D+03,0.5408254967515636D+03,
     9  0.7694973719208770D+03,0.1091746605188409D+04,
     *  0.1544861551557812D+04,0.2181041149073148D+04,
     1  0.3074292920551870D+04,0.4332640881446983D+04,
     2  0.6123710738601623D+04,0.8742550408140323D+04,
     *  0.1287901503616169D+05/
c        
        data ws25/
     1  0.4405070603487727D+00,0.1059601682287198D+01,
     2  0.1770085638476065D+01,0.2649093794953297D+01,
     3  0.3813410326224205D+01,0.5435777081518133D+01,
     4  0.7759718956988919D+01,0.1111713314392624D+02,
     5  0.1595997499628655D+02,0.2291358613923074D+02,
     6  0.3285106364649719D+02,0.4699416112293652D+02,
     7  0.6705093719656444D+02,0.9540204270704585D+02,
     8  0.1353580710502932D+03,0.1915135517051573D+03,
     9  0.2702411329363108D+03,0.3804014324069568D+03,
     *  0.5344188670611520D+03,0.7500868457862241D+03,
     1  0.1054043188833117D+04,0.1489666508150906D+04,
     2  0.2138637544396438D+04,0.3197566570469574D+04,
     *  0.5432840286100413D+04/
c
        data rl27/
     1  0.1575353168041662D+00,0.8415492372334136D+00,
     2  0.2121591236870803D+01,0.4093331013436899D+01,
     3  0.6920115664949059D+01,0.1086438046965609D+02,
     4  0.1633058975440849D+02,0.2391877457638511D+02,
     5  0.3449057783656910D+02,0.4925640830296323D+02,
     6  0.6989688699282495D+02,0.9873282917931064D+02,
     7  0.1389591925585043D+03,0.1949627791870424D+03,
     8  0.2727505571634212D+03,0.3805275873928339D+03,
     9  0.5294804963648602D+03,0.7348495726475511D+03,
     *  0.1017405228043853D+04,0.1405494047114722D+04,
     1  0.1937909191938474D+04,0.2668088570734987D+04,
     2  0.3670803017726510D+04,0.5054254438995545D+04,
     3  0.6985803095539778D+04,0.9760555602123270D+04,
     *  0.1407256148025137D+05/
c        
        data ws27/
     1  0.4055767002886141D+00,0.9705248998889890D+00,
     2  0.1605160062074552D+01,0.2365023429670761D+01,
     3  0.3331564803403754D+01,0.4623485456960240D+01,
     4  0.6407859739653024D+01,0.8910697290257594D+01,
     5  0.1243265431435451D+02,0.1737631118717147D+02,
     6  0.2428759742560504D+02,0.3391154455874941D+02,
     7  0.4726527131798085D+02,0.6573355591819771D+02,
     8  0.9119665320283983D+02,0.1262042411309599D+03,
     9  0.1742174063094966D+03,0.2399489289904696D+03,
     *  0.3298400309953051D+03,0.4527355808477365D+03,
     1  0.6208954790870368D+03,0.8517391726980110D+03,
     2  0.1171343995543024D+04,0.1622518648477326D+04,
     3  0.2286712251653505D+04,0.3361764751766289D+04,
     *  0.5623841237445988D+04/
c        
        data rl30/
     1  0.1411423972295074D+00,0.7518544578176630D+00,
     2  0.1885389916680853D+01,0.3607680322632006D+01,
     3  0.6027946667178653D+01,0.9316386503008391D+01,
     4  0.1372871621968145D+02,0.1963721622987381D+02,
     5  0.2756770412542085D+02,0.3824453401964708D+02,
     6  0.5264959920136206D+02,0.7210283556405059D+02,
     7  0.9837099848102204D+02,0.1338124264996038D+03,
     8  0.1815702961665052D+03,0.2458299530238740D+03,
     9  0.3321572979331380D+03,0.4479415496504093D+03,
     *  0.6029792405815838D+03,0.8102482967375834D+03,
     1  0.1086932459010325D+04,0.1455794716182759D+04,
     2  0.1947062186830640D+04,0.2601100350657114D+04,
     3  0.3472425979246390D+04,0.4636264262130348D+04,
     4  0.6200607260248036D+04,0.8332172605190495D+04,
     5  0.1132626743607845D+05,0.1588380873149953D+05/
c        
        data ws30/
     1  0.3631378106560880D+00,0.8640618468665494D+00,
     2  0.1413798428821516D+01,0.2048694644634513D+01,
     3  0.2819793879452562D+01,0.3799184973967351D+01,
     4  0.5086909177675927D+01,0.6816992798108090D+01,
     5  0.9163744213994908D+01,0.1235223841866342D+02,
     6  0.1667608382196640D+02,0.2252246674579401D+02,
     7  0.3040373307416541D+02,0.4099873475014893D+02,
     8  0.5520897436357335D+02,0.7423069770650116D+02,
     9  0.9964591770974953D+02,0.1335425660209947D+03,
     *  0.1786789575117243D+03,0.2387005954518417D+03,
     1  0.3184307472961597D+03,0.4242882359696212D+03,
     2  0.5649099220553905D+03,0.7521492355944543D+03,
     3  0.1002854514249030D+04,0.1342439088153352D+04,
     4  0.1813082986379242D+04,0.2496149216616223D+04,
     5  0.3591793302474308D+04,0.5891073903552565D+04/
c
        data rl32/
     1  0.1318864482433561D+00,0.7015471150019385D+00,
     2  0.1754534893633177D+01,0.3343522031312329D+01,
     3  0.5554343974051792D+01,0.8518449145098737D+01,
     4  0.1243025434952039D+02,0.1756938307899011D+02,
     5  0.2432720447341466D+02,0.3323811672286271D+02,
     6  0.4501840720688799D+02,0.6061732075956662D+02,
     7  0.8128541426822347D+02,0.1086654362365695D+03,
     8  0.1449119260349479D+03,0.1928474483019565D+03,
     9  0.2561654939436937D+03,0.3396939087701197D+03,
     *  0.4497366091659148D+03,0.5945167177892597D+03,
     1  0.7847532687947653D+03,0.1034413889203996D+04,
     2  0.1361705472517027D+04,0.1790396764460982D+04,
     3  0.2351625467672229D+04,0.3086471477506362D+04,
     4  0.4049855990904861D+04,0.5317007937050366D+04,
     5  0.6995542343428081D+04,0.9251719537000234D+04,
     6  0.1238128461286760D+05,0.1709020154236778D+05/
c        
        data ws32/
     1  0.3392122766739988D+00,0.8048245559986995D+00,
     2  0.1309785603152200D+01,0.1882290370980470D+01,
     3  0.2560971991528591D+01,0.3399268440854769D+01,
     4  0.4470503843042344D+01,0.5872467550273891D+01,
     5  0.7731605691559419D+01,0.1020876103536085D+02,
     6  0.1350890030615452D+02,0.1789585249869135D+02,
     7  0.2371215182678357D+02,0.3140440064496308D+02,
     8  0.4155565804990520D+02,0.5492662010173193D+02,
     9  0.7250846888948826D+02,0.9559170636032601D+02,
     *  0.1258546988070001D+03,0.1654796226397787D+03,
     1  0.2173051231998507D+03,0.2850288901950443D+03,
     2  0.3734862411053285D+03,0.4890442690134872D+03,
     3  0.6401973968250038D+03,0.8385465686303599D+03,
     4  0.1100575618227509D+04,0.1451278189828796D+04,
     5  0.1932536678828696D+04,0.2625950296252501D+04,
     6  0.3733546130586620D+04,0.6056153284886580D+04/
c
        h=0.001
c
        if(ndigits .eq. 2) 
     1      call onedfmm_arrmove(rl6,ws6,rl,ws,6,nrl)
c
        if(ndigits .eq. 3) 
     1      call onedfmm_arrmove(rl8,ws8,rl,ws,8,nrl)
c
        if(ndigits .eq. 4) 
     1      call onedfmm_arrmove(rl10,ws10,rl,ws,10,nrl)
c
        if(ndigits .eq. 5)
     1      call onedfmm_arrmove(rl12,ws12,rl,ws,12,nrl)
c
        if(ndigits .eq. 6)
     1      call onedfmm_arrmove(rl14,ws14,rl,ws,14,nrl)
c
        if(ndigits .eq. 7) 
     1      call onedfmm_arrmove(rl16,ws16,rl,ws,16,nrl)
c
        if(ndigits .eq. 8) 
     1      call onedfmm_arrmove(rl18,ws18,rl,ws,18,nrl)
c
        if(ndigits .eq. 9) 
     1      call onedfmm_arrmove(rl20,ws20,rl,ws,20,nrl)
c
        if(ndigits .eq. 10)
     1      call onedfmm_arrmove(rl22,ws22,rl,ws,22,nrl)
c
        if(ndigits .eq. 11)
     1      call onedfmm_arrmove(rl25,ws25,rl,ws,25,nrl)
c
        if(ndigits .eq. 12)
     1      call onedfmm_arrmove(rl27,ws27,rl,ws,27,nrl)
c
        if(ndigits .eq. 13)
     1      call onedfmm_arrmove(rl30,ws30,rl,ws,30,nrl)
c
        if(ndigits .eq. 14)
     1      call onedfmm_arrmove(rl32,ws32,rl,ws,32,nrl)
c
        return
        end
c 
c 
c 
c 
c 
        subroutine onedfmm_arrmove(xs1,ws1,xs2,ws2,nrl1,nrl2)
        implicit real *8 (a-h,o-z)
        save
        dimension xs1(1),ws1(1),xs2(1),ws2(1)
c 
        nrl2=nrl1
c        
        do 1200 i=1,nrl2
        xs2(i)=xs1(i)
        ws2(i)=ws1(i)
 1200 continue
        return
        end
  


