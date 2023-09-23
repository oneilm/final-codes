        implicit real *8 (a-h,o-z)
        real *8 a(3000 000),b(3000 000)
c 
        CALL PRINI(6,13)
C 
c       . . . construct the testing parameters
c 
         PRINT *, 'ENTER n'
         READ *,n
         CALL PRInf('n=*',n,1)
  
         na=n*10
cccc         na=n/2
         nb=n
         rlammax=1550
         rlammax=5000
         rlammax=500
  

        call mach_zero(zero_mach)

        call prin2('zero_mach=*',zero_mach,1)

cccc        stop         

         eps=1.0d-12
cccc         eps=1.0d-10
         call matrtest(n,a,na,b,nb,rlammax,eps)
  
  
         stop
         end
c 
c 
c 
c 
c 
        subroutine matrtest(n,a,na,b,nb,rlammax,eps)
        implicit real *8 (a-h,o-z)
        save
        dimension xs1(10000),xs2(10000),
     1      ws1(10000),rlams(10000),
     2      wsa(10000),ipivots(10000),errs(10000),
     3      freqs(10000),bads(10000),bsizes(10000),
     4      w(10 000 000),coefs(1000 000),
     5      rnorms(10000)

        complex *16 a(n,na),b(n,nb),ima
c
        data ima/(0.0d0,1.0d0)/ 
  
c 
c        construct the discretizations of two intervals
c 
        itype=1
        call legeexps(itype,nb,xs1,u,v,ws1)
c 
        itype=1
        call legeexps(itype,n,xs2,u,v,ws1)
c 
        do 1200 i=1,n
c 
        xs2(i)=xs2(i)+4
c 
 1200 continue
  
  
        itype=1
        call legeexps(itype,na,rlams,u,v,wsa)
  
c 
        do 1300 i=1,na
c 
        rlams(i)=(rlams(i)+1)*rlammax/2
 1300 continue
  
  
        call prin2('xs1=*',xs1,nb)
        call prin2('xs2=*',xs2,n)
        call prin2('rlams=*',rlams,na)
  
  
cccc        stop
  
c 
c        construct the test matrices a,b
c 
        done=1
        do 1600 i=1,nb
        do 1400 j=1,n
c 
           b(j,i)=done/(xs2(j)-xs1(i)) + xs2(i)*ima
           b(j,i)=done/(xs2(j)-xs1(i)) 
c 
 1400 continue
 1600 continue
c 
c        construct the test matrices a,b
c 
        done=1
        do 2600 i=1,na
        do 2400 j=1,n
c 
cccc        a(j,i)=exp(-(xs2(j)-3)*rlams(i)*ima)
        a(j,i)=exp(-(xs2(j)-3)*rlams(i)) + ima
cccc        a(j,i)=exp(-(xs2(j)-3)*rlams(i)) 

 2400 continue
 2600 continue
  
  
        iftest=1
        call cabskel(ier,n,a,na,b,nb,eps,iftest,ncols,ipivots,
     1      coefs,errs,rnorms,bads,bsizes,w)
  
        call prinf('after cabskel, ipivots=*',ipivots,ncols)
        call prin2('after cabskel, bads=*',bads,ncols)
        call prin2('bads=*',bads,ncols)
        call prin2('coefs=*',coefs,ncols*nb*2)
        call prinf('and ncols=*',ncols,1)
        call prinf('and ier=*',ier,1)
c 
c        extract the frequencies
c 
        do 3200 i=1,ncols
c 
        freqs(i)=rlams(ipivots(i))
 3200 continue
c 
        call prin2('and freqs=*',freqs,ncols)
        call prin2('and bsizes=*',bsizes,ncols+1)
        call prin2('and errs=*',errs,nb)
  
  
        return
        end
c 
c 
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code, and the beginning
c        of the skeletonization code proper. This file contains
c        two user-callable subroutines: cabskel and cabskel_basic.
c        Following is a brief description of these two subroutines.
c 
c   cabskel - skeletonizes the columns of the complex input matrix 
c        a IN TERMS OF THE COLUMNS OF THE INPUT MATRIX B; then, it 
c        uses the subroutine cabskel_basic to approximate all 
c        columns of b with the skeleton columns of a.
c 
c   cabskel_basic - skeletonizes the columns of the complex input 
c        matrix a IN TERMS OF THE COLUMNS OF THE complex INPUT 
c        MATRIX B.
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c 
        subroutine cabskel(ier,n,a,na,b,nb,eps,iftest,
     1      ncols,ipivots,coefs,errs,rnorms,bads,bsizes,w)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,na),b(n,nb),ipivots(1),errs(1),
     1      bads(1),bsizes(1),w(1),coefs(1),rnorms(1)
c 
c        This subroutine skeletonizes the columns of the
c        complex input matrix a IN TERMS OF THE COLUMNS OF 
C        THE COMPLEX INPUT MATRIX B; then, uses the subroutine 
c        ncleastsq to approximate all columns of b with the 
c        skeleton columns of a.
c 
c                     Input parameters:
c 
c  n - the number of rows in each of the matrices a, b
c  a - the matrix to be skeletonized; destroyed utterly
c  na - the number of columns in the matrix a
c  b - the matrix with respect to which the skeletonization
c        is to be performed; destroyed utterly
c  nb - the number of columns in the matrix b
c  eps - the precision to which the skeletonization is to be
c        performed
c  iftest - the integer parameter telling the subroutine
c        whether it should test the obtained approximation.
c        Setting iftest=1 will cause the subroutine to
c        perform the testing, and store the result in the
c        output array errs (see below). Setting iftest=0 will
c        cause the testing to be skipped
c 
c                     Output parameters:
c 
c  ier - error return code.
c     ier=0 means that the prescribed accuracy eps has been achieved
c     ier=4 means that the norm of the leading column in matrix A
c           became less than 10*(machine_eps) before the prescribed 
c           accuracy eps has been achieved. This is NOT a fatal 
c           error
c  ncols - the number of columns in the skeleton of a
c  ipivots - integer array containing the numbers of the
c        columns in the skeleton of a
c  coefs - (ncols \times nb) - array of coefficients, expressing
c        all columns in the array b as linear combinations of
c        the selected columns of a, via the formula
c 
c        b(i)=\sum_{j=1}^{ncols} coefs(j,i) * a(ipivots(j))
c 
c        (in the above, a(i) denotes the i-th column of a,
c        and b(i) denotes the i-th column of b)
c  errs - the i-th elemenst of the array errs is the L^2 error
c        in the approximation of the i-th columns of b
c  rnorms - the array of length ncols containing the norms of
c        the pivot columns of a at the time of piviting
c  bads - the coefficients describing the "badness" of the obtained
c        skeleton. Explanation: if any one of the ncols elements of
c        the array bads is (for example) 1000, the approximation
c        will lose about 3 digits. If all of the elements of the
c        array bads are reasonably small, it is very likely that the
c        approximation loses no digits. HOWEVER, IT IS POSSIBLE TO
C        CONSTRUCT UNPLEASANT (BUT HIGHLY UNLIKELY) COUNTEREXAMPLES
c  bsizes - the array describing the decay of the error of the
c        approximation to the matrix b as a function of the number
c        of columns of a used. Explanation: the k-the element of
c        bsizes is the magnitude (in L^2 sense) of the remainder
c        after the k-th step.
c 
c                     Work array:
c 
c  w - must be at least na*nb+2*na+nb+(na+nb)*n +2*n + 200
c        Complex *16 elements
c 
c        . . . create copies of the matrices a, b
c 
        ier=0
c
        ia2=1
        la2=na*n+10
        la2=la2*2
c 
        ib2=ia2+la2
        lb2=nb*n+10
        lb2=lb2*2
c 
        iarr=ib2+lb2
        larr=n+10
        larr=larr*2
c 
        iw=iarr+larr
c 
        call cleascop(b,w(ib2),n*nb)
        call cleascop(a,w(ia2),n*na)
c 
        call cabskel_basic(ier,n,w(ia2),na,w(ib2),nb,eps,
     1      ncols,ipivots,bads,bsizes,w(iw))
c 
c        use least squares to approximate all of the columns
c        of a by linear combinations of the selected columns
c        of a
c 
        call cabskel0(n,a,na,b,nb,ipivots,ncols,
     1      w(ia2),coefs,errs,rnorms,w(iw),w(iarr),iftest,eps)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cabskel0(n,a,na,b,nb,ipivots,ncols,
     1      c,coefs,errs,rnorms,w,arr,iftest,eps)
        implicit real *8 (a-h,o-z)
        save
        dimension ipivots(1),errs(1)
c 
        complex *16 a(n,na),b(n,nb),c(n,ncols),
     1      w(1),rnorms(1),coefs(ncols,nb),arr(1)
c 
c        extract from the array a the pivots stored
c        in the array ipivots
c 
        do 1400 i=1,ncols
c 
        ii=ipivots(i)
        do 1200 j=1,n
c 
        c(j,i)=a(j,ii)
 1200 continue
 1400 continue
c 
c        construct the decomposition of c
c 
        call prin2('eps=*',eps,1)
cccc        stop

        eps2=1.0d-30

cccc        call ncleastsq(c,n,ncols,eps*4,ncols2,rnorms,w)
        call ncleastsq(c,n,ncols,eps2,ncols2,rnorms,w)
c 
        call prin2('after nrleastsq, rnorms=*',rnorms,ncols)
c 
c        construct the approximations for all columns of b
c 
        do 1800 i=1,nb
c 
        call ncleasts2(w,b(1,i),coefs(1,i))
 1800 continue
c 
c       check the error of the approximation
c 
        if(iftest .eq. 0) return
c  
        do 2800 i=1,nb
c 
        do 2200 j=1,n
        arr(j)=0
 2200 continue
c 
        do 2600 j=1,ncols
        do 2400 k=1,n
c 
        arr(k)=arr(k)+coefs(j,i)*c(k,j)
 2400 continue
 2600 continue
c 
        err=0
        do 2700 j=1,n
c 
        err=err+(arr(j)-b(j,i))*conjg(arr(j)-b(j,i))
 2700 continue
c 
        err=sqrt(err)
        errs(i)=err
c 
 2800 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cabskel_basic(ier,n,a,na,b,nb,eps,
     1      ncols,ipivots,bads,bsizes,w)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,na),b(n,nb),ipivots(1),
     1      bads(1),bsizes(1),w(1)
c 
c        This subroutine skeletonizes the columns of the
c        input matrix a IN TERMS OF THE COLUMNS OF THE
C        INPUT MATRIX B. BOTH MATRICES ARE COMPLEX
c 
c                     Input parameters:
c 
c  n - the number of rows in each of the matrices a, b
c  a - the matrix to be skeletonized; destroyed utterly
c  na - the number of columns in the matrix a
c  b - the matrix with respect to which the skeletonization
c        is to be performed; destroyed utterly
c  nb - the number of columns in the matrix b
c  eps - the precision to which the skeletonization is to be
c        performed
c 
c                     Output parameters:
c 
c  ncols - the number of columns in the skeleton of a
c  ipivots - integer array containing the numbers of the
c        columns in the skeleton of a
c  bads - the coefficients describing the "badness" of the obtained
c        skeleton. Explanation: if any one of the ncols elements of
c        the array bads is (for example) 1000, the approximation
c        will lose about 3 digits. If all of the elements of the
c        array bads are reasonably small, it is very likely that the
c        approximation loses no digits. HOWEVER, IT IS POSSIBLE TO
C        CONSTRUCT UNPLEASANT (BUT HIGHLY UNLIKELY) COUNTEREXAMPLES
c  bsizes - the array describing the decay of the error of the
c        approximation to the matrix b as a function of the number
c        of columns of a used. Explanation: the k-the element of
c        bsizes is the magnitude (in L^2 sense) of the remainder
c        after the k-th step.
c 
c                     Work array:
c 
c  w - must be at least na*nb+2*na+nb+100 complex *16 elements
c 
c        . . . allocate memory
c 
        iprods=1
        lprods=na*nb+10
        lprods=lprods*2
c 
        isizesab=iprods+lprods
        lsizesab=na+10
        lsizesab=lsizesab*2
c 
        ialphas=isizesab+lsizesab
        lalphas=na+10
        lalphas=lalphas*2
c 
        ibetas=ialphas+lalphas
        lbetas=nb+10
        lbetas=lbetas*2
c 
        call mach_zero(d) 
        aeps=d*10
cccc        aeps=d*2
        call cabskel_basic0(ier,n,a,na,b,nb,w(iprods),
     1      eps,aeps,ipivots,ncols,bads,bsizes,
     2      w(isizesab),w(ialphas),w(ibetas))
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cabskel_basic0(ier,n,a,na,b,nb,prods,
     1      eps,aeps,ipivots,ncols,bads,bsizes,
     2      sizesab,alphas,betas)
        implicit real *8 (a-h,o-z)
        save
c
        complex *16 a(n,1),b(n,1),prods(nb,na),
     1      alphas(1),betas(1),
     2      cd,prod
c 
        dimension sizesab(1),
     1      ipivots(1),
     2      bads(1),bsizes(1)
c
        ier=0
c 
c        skeletonize a with respect to b
c 
c        . . . initialize the array of pivots
c 
        do 1200 i=1,na
        ipivots(i)=i
 1200 continue
c 
        call cleascap(a,a,n*na,cd)
c
        sizea0=cd
        sizea0=sqrt(sizea0)

        call cleascap(b,b,n*nb,cd)
        sizeb0=cd
        sizeb0=sqrt(sizeb0)
c
        bsizes(1)=sizeb0
        bsize_prev=sizeb0
c 
        epsb=eps*sizeb0
c 
c        construct the matrix of inner products of columns of a,b
c 
        call cabskel_fndprods(n,a,na,b,nb,prods,sizesab)
c 
c        construct the skeletonization of a with respect to b
c 
        do 3000 i=1,na
c
        call prinf('i=*',i,1)
c 
c        find the largest column of the matrix prods
c 
        imax=i
        dmax=sizesab(i)
        do 1400 j=i+1,na
c 
        if(sizesab(j) .le. dmax) goto 1400
            dmax=sizesab(j)
            imax=j
c 
 1400 continue
c 
c        put the found column of a in the i-th position
c 
        do 1500 j=1,n
c 
        cd=a(j,imax)
        a(j,imax)=a(j,i)
        a(j,i)=cd
c 
 1500 continue
c 
        do 1600 j=1,nb
c 
        cd=prods(j,imax)
        prods(j,imax)=prods(j,i)
        prods(j,i)=cd
c 
 1600 continue
c 
        iii=ipivots(i)
        ipivots(i)=ipivots(imax)
        ipivots(imax)=iii
c 
        d=sizesab(i)
        sizesab(i)=sizesab(imax)
        sizesab(imax)=d
c 
c       orthogonalize the i-th column of a to all preceding columns
c 
        if(i .eq. 1) goto 1700
c 
        do 1650 j=1,i-1
c 
        call cleascap(a(1,j),a(1,i),n,cd)
c 
        do 1630 k=1,n
        a(k,i)=a(k,i)-a(k,j)*cd
 1630 continue
 1650 continue
c 
 1700 continue
c 
c        normalize the i-th column of a
c 
        call cleascap(a(1,i),a(1,i),n,cd)
        d=cd

        call prin2('normalizing factor d=*',d,1)
        call prin2('and d/sizea0=*',d/sizea0,1)

        if(d/sizea0 .lt. aeps/n/na) then 
            ier=4
            goto 3200
        endif

cccc        if(d/sizea0 .lt. 1.0d-15) goto 3200
c 
        d=1/sqrt(d)
c 
        dnorm=d
        do 1800 j=1,n
c 
        a(j,i)=a(j,i)*d
 1800 continue
c 
c        orthogonalize all columns of b to the i-th
c        column of a
c 
        bsize=0
c 
        do 2400 j=1,nb
c 
        call cleascap(b(1,j),a(1,i),n,prod)
c
        betas(j)=prod
c 
        do 2200 k=1,n
        b(k,j)=b(k,j)-a(k,i)*prod
c 
        bsize=bsize+b(k,j)*conjg(b(k,j))
 2200 continue
 2400 continue
c 
        bsize=sqrt(bsize)
c 
        bads(i)=bsize*dnorm
        bsizes(i+1)=bsize
c 
        ncols=i
        if(bsize .lt. epsb) goto 3200
c 
c        orthogonalize all subsequent columns of a to the i-th
c        column of a
c 
        if(i .eq. na) goto 2900
c 
        do 2500 j=1,na
        alphas(j)=0
 2500 continue
c 
        do 2800 j=i+1,na
c 
        call cleascap(a(1,j),a(1,i),n,prod)
        alphas(j)=prod
c 
        do 2600 k=1,n
        a(k,j)=a(k,j)-a(k,i)*prod
 2600 continue
 2800 continue
c 
 2900 continue
c 
c       adjust the matrix prods
c 
        call cabskel_adjprods(n,a,na,b,nb,prods,sizesab,
     1      alphas,betas,i)
c
        if(bsize_prev*1.0d-4 .gt. bsizes(i+1) ) then
            call cabskel_fndprods(n,a,na,b,nb,prods,sizesab)
c
            bsize_prev=bsizes(i+1)
        endif
c
 3000 continue
c  
 3200 continue
c 
c       bubble them things
c 
cccc        call cabskel_bubble(ipivots,ncols)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cabskel_bubble(ipivots,ncols)
        save
        dimension ipivots(1)
c 
        do 1400 i=1,ncols
        do 1200 j=1,ncols-1
c 
        if(ipivots(j) .lt. ipivots(j+1)) goto 1200
c 
        k=ipivots(j)
        ipivots(j)=ipivots(j+1)
        ipivots(j+1)=k
 1200 continue
 1400 continue
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine cabskel_adjprods(n,a,na,b,nb,prods,sizesab,
     1      alphas,betas,iii)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,na),b(n,nb),prods(nb,na),
     1      alphas(1),betas(1)
c 
        dimension sizesab(1)
c 
c        adjust the matrix of inner products of columns of a,b
c 
        do 2400 i=iii+1,na
        do 2200 j=1,nb
c 
        prods(j,i)=prods(j,i)-conjg(betas(j)) * alphas(i)
c 
 2200 continue
c 
 2400 continue
c 
        do 2600 j=1,nb
c 
        prods(j,iii)=0
 2600 continue
c 
c        recalculate the vector sizesab
c 
        do 3400 i=1,na
c 
        sizesab(i)=0
        do 3200 j=1,nb
c 
        sizesab(i)=sizesab(i)+prods(j,i)*conjg(prods(j,i))
c 
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
        subroutine cabskel_fndprods(n,a,na,b,nb,prods,sizesab)
        implicit real *8 (a-h,o-z)
        save
        complex *16  a(n,na),b(n,nb),prods(nb,na)
        dimension sizesab(1)
c 
c        construct the matrix of inner products of columns of a,b
c 
  
        do 1400 i=1,na
c 
        sizesab(i)=0
        do 1200 j=1,nb
c 
        call cleascap(a(1,i),b(1,j),n,prods(j,i))
c 
        sizesab(i)=sizesab(i)+prods(j,i)*conjg(prods(j,i))
c 
 1200 continue
c 
 1400 continue
c 
        return
        end
