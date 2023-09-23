        implicit real *8 (a-h,o-z)
        dimension zs1(2,100 000),zs2(2,100 000),
     1      zs3(2,100 000),w(1000 000),
     2      ipivots(10000),fillin(1000 000),eval(1000 000),
     3      expand(1000 000),askel(1000 000),bskel(1000 000)
        real *8 ima,rk,a(100 000),b(100 000)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
         PRINT *, 'ENTER n'
         READ *,n
         CALL PRINF('n=*',n,1 )
c 
c        Construct and plot the test structure
c 
        rk=1
        eps=1.0d-15
c 
        n2=n
        n1=n/2
        n3=n*2
c 
        n2=n
        n1=n
        n3=n
c 
        call circles3_bld(n1,zs1,n2,zs2,n3,zs3)
  
        iw=21
        call zquaplot3(iw,zs1,n1,2,zs2,n2,2,zs3,n3,2,
     1      'test structure*')
c 
c        create the interaction matrices
c 
  
c 
c        . . . between the first and the second circles
c 
        diag=0
  
        call prin2('before matr_crea, rk=*',rk,2)
        call prinf('before matr_crea, n1=*',n1,1)
        call prinf('before matr_crea, n2=*',n2,1)
  
        call matr_crea(zs1,n1,zs2,n2,diag,a,rk)
c 
c        . . . between the second and the third circles
c 
        call matr_crea(zs1,n2,zs2,n3,diag,b,rk)
  
        eps=1.0d-13
  
c 
        lenw=1000 000
        call rskel_twomatr(ier,a,b,n1,n2,n3,eps,
     1      ncols,ipivots,fillin,eval,expand,w,lenw,lused)
c 
  
        call prinf('after matr_test, ncols=*',ncols,1)
        call prinf('after matr_test, lused=*',lused,1)
  
  
        call basic_test(a,b,n1,n2,n3,ncols,ipivots,
     1     askel,bskel,fillin,eval,expand)
  
  
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine basic_test(a,b,n,m,k,ncols,ipivots,
     1     askel,bskel,fillin,eval,expand)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension a(n,m),b(m,k),askel(n,ncols),bskel(ncols,k),
     1      ipivots(1),fillin(1),eval(1),expand(1)
c 
        dimension ab(1000 000),uu(1000 000),vv(1000 000),
     1      errs(1000000)
c 
c       extract the submatrices askel, bskel from the
c       matrices a, b, respectively
c 
        do 1600 i=1,ncols
c 
        ii=ipivots(i)
        do 1200 j=1,n
c 
        askel(j,i)=a(j,ii)
 1200 continue
c 
         do 1400 j=1,k
c 
         bskel(i,j)=b(ii,j)
 1400 continue
c 
 1600 continue
c 
c        multiply the matrices a,b
c 
        call rskel_matmult(a,b,ab,n,m,k)
c 
c       test the matrix fillin
c 
        call rskel_matmult(fillin,bskel,uu,ncols,ncols,k)
        call rskel_matmult(askel,uu,vv,n,ncols,k)
c 
c        calculate the difference between vv and ab
c 
        errmax=0
        dd=0
        ddd=0
        do 2200 i=1,n*k
c 
        d=abs(vv(i)-ab(i))
        if(d .gt. errmax) errmax=d
c 
        dd=dd+d**2
        ddd=ddd+ab(i)**2
c 
 2200 continue
c 
        errrel=sqrt(dd/ddd)
        call prin2('in basic_test, errmax=*',errmax,1)
        call prin2('in basic_test, errrel=*',errrel,1)
c 
c        test the matrix eval
c 
        call rskel_matmult(eval,bskel,uu,m,ncols,k)
        call rskel_matmult(a,uu,vv,n,m,k)
c 
c        calculate the difference between vv and ab
c 
        errmax=0
        dd=0
        ddd=0
        do 3200 i=1,n*k
c 
        d=abs(vv(i)-ab(i))
        errs(i)=d
  
        if(d .gt. errmax) errmax=d
c 
        dd=dd+d**2
        ddd=ddd+ab(i)**2
c 
 3200 continue
c 
        errrel=sqrt(dd/ddd)
        call prin2('in basic_test errmax for eval =*',errmax,1)
        call prin2('in basic_test errrel for eval =*',errrel,1)
cccc        call prin2('and errs for eval=*',errs,n*k)
c 
c        test the matrix expand
c 
        call rskel_matmult(expand,b,uu,ncols,m,k)
        call rskel_matmult(askel,uu,vv,n,ncols,k)
  
c 
c        calculate the difference between vv and ab
c 
        errmax=0
        dd=0
        ddd=0
        do 4200 i=1,n*k
c 
        d=abs(vv(i)-ab(i))
        if(d .gt. errmax) errmax=d
  
c 
        dd=dd+d**2
        ddd=ddd+ab(i)**2
c 
 4200 continue
c 
        errrel=sqrt(dd/ddd)
        call prin2('in basic_test errmax for expand =*',errmax,1)
        call prin2('in basic_test errrel for expand =*',errrel,1)
  
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine matr_crea(zs1,n1,zs2,n2,diag,a,rk)
        implicit real *8 (a-h,o-z)
c 
        save
        real *8 a(n1,n2),rk,z,diag
        dimension zs1(2,1),zs2(2,1)
c 
c        construct the matrix of interactions
c 
cccc        call prinf('in matr_crea, n1=*',n1,1)
cccc        call prinf('in matr_crea, n2=*',n2,1)
c 
        do 1400 i=1,n1
        do 1200 j=1,n2
c 
cccc        call prinf('i=*',i,1)
cccc        call prinf('j=*',j,1)
  
        d=(zs1(1,i)-zs2(1,j))**2+(zs1(2,i)-zs2(2,j))**2
c 
        if(d .lt. 1.0d-9) goto 1200
  
        d=sqrt(d)
c 
        ifexpon=1
        z=d*rk
c 
        call jy0rea(z,fj,fy)
        a(i,j)=fy
 1200 continue
c 
ccc        a(i,i)=diag
 1400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine circles3_bld(n1,zs1,n2,zs2,n3,zs3)
        implicit real *8 (a-h,o-z)
        save
        real *8 zs1(2,1),zs2(2,1),zs3(2,1)
c 
c       construct the three circles
c 
        coef=0.8
        r=coef
        done=1
        pi=atan(done)*4
        h=2*pi/n1
c 
        do 1600 i=1,n1
c 
        zs1(1,i)=r*cos((i-1)*h)-2
        zs1(2,i)=r*sin((i-1)*h)
c 
 1600 continue
c 
        h=2*pi/n2
c 
        do 2600 i=1,n2
c 
        zs2(1,i)=r*cos((i-1)*h)
        zs2(2,i)=r*sin((i-1)*h)
c 
 2600 continue
c 
        h=2*pi/n3
c 
        do 3600 i=1,n3
c 
        zs3(1,i)=r*cos((i-1)*h)+2
        zs3(2,i)=r*sin((i-1)*h)
c 
 3600 continue
c 
        return
        end
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code and the beginning
c        of the two-sided matrix compression code proper
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c        This file contains two user-callable subroutines:
c        rskel_twomatr and rskel_twomatr_basic. Following is
c        a brief description of these subroutines.
c 
c 
c   rskel_twomatr - starts with constructing a two-sided
c        skeletonization of the matrix
c 
c                c(n,k)=a(n,m)*b(m,k).                              (1)
c 
c        In other words, it finds the rank (to precision eps)
c        ncols of the matrix a*b and a subset of ncols coordinates
c        in R^m such that the matrix
c 
c        a      * b         = c                                     (2)
c         |skel    |skel       |skel
c 
c        has the smallest singular value that is reasonably close
c        to the smallest singular value of c.
c 
c        Then, this subroutine constructs three mappings
c        fillin(ncols,ncols), eval(m,ncols), expand(ncols,m)
c        such that
c 
c        a(n,ncols) * fillin (ncols,ncols) * b(ncols,k)  = a*b,    (3)
c         |skel                               |skel
c 
c        a(n,m) * eval(m,ncols) * b(ncols,k)  = a*b,               (4)
c                                  |skel
c 
c        a(n,m) * expand(ncols,m) * b(m,k)    = a*b,               (5)
c         |skel
c 
c        all of the above equalities to be understood to precision
c        eps.
c 
c   rskel_twomatr_basic - constructs a two-sided skeletonization
c        of the matrix
c 
c                c(n,k)=a(n,m)*b(m,k).                              (1)
c 
c        In other words, it finds the rank (to precision eps)
c        ncols of the matrix a*b and a subset of ncols coordinates
c        in R^m such that the matrix
c 
c        a      * b         = c                                     (2)
c         |skel    |skel       |skel
c 
c        has the smallest singular value that is reasonably close
c        to the smallest singular value of c.
c 
c        PLEASE NOTE THAT BOTH MATRICES A,B ARE UTTERLY DESTROYED
C        BY THIS SUBROUTINE!!!!
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine rskel_twomatr(ier,a,b,n,m,k,eps,
     1      ncols,ipivots,fillin,eval,expand,w,lenw,lused)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension a(n,m),b(m,k),fillin(1),eval(1),
     1      expand(1),ipivots(1),w(1)
c 
c        This subroutine starts with constructing a two-sided
c        skeletonization of the matrix
c 
c                c(n,k)=a(n,m)*b(m,k).                              (1)
c 
c        In other words, it finds the rank (to precision eps)
c        ncols of the matrix a*b and a subset of ncols coordinates
c        in R^m such that the matrix
c 
c        a      * b         = c                                     (2)
c         |skel    |skel       |skel
c 
c        has the smallest singular value that is reasonably close
c        to the smallest singular value of c.
c 
c        Then, this subroutine constructs three mappings
c        fillin(ncols,ncols), eval(m,ncols), expand(ncols,m)
c        such that
c 
c        a(n,ncols) * fillin (ncols,ncols) * b(ncols,k)  = a*b,    (3)
c         |skel                               |skel
c 
c        a(n,m) * eval(m,ncols) * b(ncols,k)  = a*b,               (4)
c                                  |skel
c 
c        a(n,m) * expand(ncols,m) * b(m,k)    = a*b,               (5)
c         |skel
c 
c        all of the above equalities to be understood to precision
c        eps.
c 
c                        Input parameters:
c 
c  a,b,n,m,k - parameters in (1)-(5)
c  eps - the precision to which the calculations are to be performed
c  lenw - the length of the work array w supplied by the user to
c        the subroutine
c 
c                        Output parameters:
c 
c  ier - error return code.
c    ier=0 means successful execution
c    ier > 0 means that the subroutine ran out of memory in array w
c  ipivots - integer array of length ncols containing the sequence
c        numbers of the columns of a and the rows of b constituting
c        the rows and columns of the matrices a     , b
c                                              |skel   |skel
c        respectively.
c  fillin, eval, expand - the matrices in (3)-(5)
c  lused - the amount of space in array w (in real *8 words)
c        actually used by the subroutine
c 
c 
c       copy the matrices a,b into the work array
c 
        ier=0
c 
        ia2=1
        la2=n*m+10
c 
        ib2=ia2+la2
        lb2=m*k+10
c 
        irnorms=ib2+lb2
        lrnorms=m+10
c 
        lused=irnorms+lrnorms
c 
        if(lused .gt. lenw) then
            ier=1024
            return
        endif
c 
        call rskel_arrmove(a,w(ia2),n*m)
        call rskel_arrmove(b,w(ib2),k*m)
c 
c        compress the matrix a*b
c 
        call rskel_twomatr_basic(w(ia2),w(ib2),n,m,k,eps,
     1      ipivots,ncols,w(irnorms) )
c 
c        construct the matrices fillin, eval, expand
c 
        iaskel=1
        laskel=n*ncols+10
c 
        ibskel=iaskel+laskel
        lbskel=k*ncols
c 
        iab=ibskel+lbskel
        lab=n*k
c 
        iw2=iab+lab
        lw2=lenw-iw2
c 
        if(lused .lt. iw2) lused=lw2
c 
        if(lw2 .lt. 1000) then
            ier=512
            return
        endif
c 
        eps2=eps*10
        call prinf('and lw2=*',lw2,1)
  
        call rskel_twomatr0(jer,a,b,n,m,k,ncols,ipivots,
     1     w(iaskel),w(ibskel),fillin,eval,expand,w(iab),
     2      w(iw2),lw2,eps2,ltot)
c 
        if(jer .ne. 0) ier=128
c 
        lused=lused+ltot
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine rskel_twomatr0(ier,a,b,n,m,k,ncols,ipivots,
     1     askel,bskel,fillin,eval,expand,ab,w,lw,eps,ltot)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension a(n,m),b(m,k),askel(n,ncols),bskel(ncols,k),
     1      ipivots(1),fillin(1),eval(1),expand(1)
c 
        dimension ab(1),w(1)
c 
c       extract the submatrices askel, bskel from the
c       matrices a, b, respectively
c 
        ier=0
c 
        do 1600 i=1,ncols
c 
        ii=ipivots(i)
        do 1200 j=1,n
c 
        askel(j,i)=a(j,ii)
 1200 continue
c 
         do 1400 j=1,k
c 
         bskel(i,j)=b(ii,j)
 1400 continue
c 
 1600 continue
c 
c        multiply the matrices a,b
c 
        call rskel_matmult(a,b,ab,n,m,k)
c 
c        construct the matrix fillin
c 
        delta=0
cccc        call leamatr(jer,askel,bskel,ab,n,ncols,ncols,k,eps,
        call rleamatr(jer,askel,bskel,ab,n,ncols,ncols,k,eps,
     1      delta,fillin,w,lw,ltot)
c 
        if(jer .ne. 0) then
            ier=64
            return
        endif
c 
c        construct the matrix eval
c 
cccc        call leamatr(jer,a,bskel,ab,n,m,ncols,k,eps,
        call rleamatr(jer,a,bskel,ab,n,m,ncols,k,eps,
     1      delta,eval,w,lw,ltot)
c 
        if(jer .ne. 0) then
            ier=32
            return
        endif
c 
c        construct the matrix expand
c 
cccc        call leamatr(jer,askel,b,ab,n,ncols,m,k,eps,
        call rleamatr(jer,askel,b,ab,n,ncols,m,k,eps,
     1      delta,expand,w,lw,ltot)
c 
        if(jer .ne. 0) ier=16
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine rskel_matmult(a,x,c,n,m,k)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),x(m,k),c(n,k)
c 
        do 1600 i=1,n
        do 1400 j=1,k
c 
        cd=0
        do 1200 ll=1,m
c 
        cd=cd+a(i,ll)*x(ll,j)
 1200 continue
c 
        c(i,j)=cd
 1400 continue
 1600 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine rskel_twomatr_basic(a,b,n,m,k,eps,
     1      ipivots,ncols,rnorms)
        implicit real *8 (a-h,o-z)
        save
        dimension ipivots(1),rnorms(1)
        real *8 a(n,m),b(m,k),cd
c 
c        This subroutine constructs a two-sided skeletonization
c        of the matrix
c 
c                c(n,k)=a(n,m)*b(m,k).                              (1)
c 
c        In other words, it finds the rank (to precision eps)
c        ncols of the matrix a*b and a subset of ncols coordinates
c        in R^m such that the matrix
c 
c        a      * b         = c                                     (2)
c         |skel    |skel       |skel
c 
c        has the smallest singular value that is reasonably close
c        to the smallest singular value of c.
c 
c        PLEASE NOTE THAT BOTH MATRICES A,B ARE UTTERLY DESTROYED
C        BY THIS SUBROUTINE!!!!
C 
c 
c                     Input parameters:
c 
c  a,b - the two matrices in (1)
c  n,m,k - the dimensionalities in (1)
c  eps - the precision with which the ranks are calculated
c 
c                     Output parameters:
c 
c  ipivots - integer array containing the sequence numbers of the
c        ncols coordinates in R^m spanning the whole matrix c
c        to precision eps
c  ncols - the rank of c to precision eps
c 
c                     Work arrays:
c 
c  rnorms - must be at least ncols real *8 locations long
c 
c        . . . initialize the array of pivots
c 
        do 1100 i=1,m
        ipivots(i)=i
 1100 continue
c 
c       . . . prepare the arrays of values of norms
c             of columns of a and rows of b
c 
        done=1
c 
        m0=1
        call rskel_rnorms_calc(a,b,n,m,k,rnorms,m0,tot)
c 
        thresh=tot*eps**2
cccc        thresh=thresh/1.0d6
  
  
  
        call prin2('thresh=*',thresh,1)
c 
c       . . . conduct screwed-up gram-schmidt iterations
c 
        do 4000 i=1,m
c 
c       find the pivot
c 
        ipivot=i
        rn=rnorms(i)
c 
        do 1500 j=i+1,m
        if(rnorms(j) .le. rn) goto 1500
        rn=rnorms(j)
        ipivot=j
 1500 continue
c 
c       put the column of a number ipivot in the i-th place,
c       and the row of b number ipivot likewise in the i-th
c       place
c 
        do 1600 j=1,n
        cd=a(j,i)
        a(j,i)=a(j,ipivot)
        a(j,ipivot)=cd
 1600 continue
c 
        do 1700 j=1,k
        cd=b(i,j)
        b(i,j)=b(ipivot,j)
        b(ipivot,j)=cd
 1700 continue
c 
        iijj=ipivots(i)
        ipivots(i)=ipivots(ipivot)
        ipivots(ipivot)=iijj
c 
        d=rnorms(ipivot)
        rnorms(ipivot)=rnorms(i)
        rnorms(i)=d
c 
c       orthogonalize the i-th column of a to all preceding ones
c 
        if(i .eq. 1) goto 2300
        do 2200 j=1,i-1
c 
        call rskel_rscap(a(1,i),a(1,j),n,cd)
c 
        do 2100 l=1,n
        a(l,i)=a(l,i)-a(l,j)*cd
 2100 continue
 2200 continue
 2300 continue
c 
c       orthogonalize the i-th row of b to all preceding ones
c 
        if(i .eq. 1) goto 2600
        do 2500 j=1,i-1
c 
        call rskel_rscap2(b,m,k,i,j,cd)
c 
        do 2400 l=1,k
        b(i,l)=b(i,l)-b(j,l)*cd
 2400 continue
 2500 continue
 2600 continue
c 
c       normalize the i-th column of a
c 
        call rskel_rscap(a(1,i),a(1,i),n,cd)
c 
        d=cd
        d=done/sqrt(d)
c 
        do 2700 j=1,n
        a(j,i)=a(j,i)*d
 2700 continue
c 
c       if the remaining part of the matrix is less than
c       thresh, get out
c 
        ncols=i
  
        call rskel_rscap2(b,m,k,i,i,cd)
        d=cd
c 
c       normalize the i-th row of b
c 
        d=done/sqrt(d)
  
        do 2800 j=1,k
        b(i,j)=b(i,j)*d
 2800 continue
c 
        if(i .eq. m) goto 3800
c 
c        orthogonalize all subsequent columns of a to the i-th
c        all subsequent columns of a
c 
        do 3200 j=i+1,m
c 
        call rskel_rscap(a(1,i),a(1,j),n,cd)
c 
        do 3000 l=1,n
        a(l,j)=a(l,j)-a(l,i)*cd
 3000 continue
 3200 continue
c 
c        orthogonalize all subsequent columns of a to
c        all subsequent columns of a
c 
        do 3600 j=i+1,m
c 
        if(rnorms(j) .lt. thresh/100) goto 3600
c 
        call rskel_rscap2(b,m,k,i,j,cd)
c 
        do 3400 l=1,k
        b(j,l)=b(j,l)-b(i,l)*cd
 3400 continue
 3600 continue
c 
 3800 continue
c 
c        recalculate the array rnorms
c 
        m0=i+1
        if(m0 .ne. m+1)
     1      call rskel_rnorms_calc(a,b,n,m,k,rnorms,m0,tot)
c 
        if(tot .lt. thresh) return
c 
 4000 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine rskel_rnorms_calc(a,b,n,m,k,rnorms,m0,tot)
        implicit real *8 (a-h,o-z)
        save
        dimension rnorms(1)
        real *8 a(n,m),b(m,k)
c 
c       . . . prepare the arrays of values of norms
c             of columns of a and rows of b
c 
        done=1
        tot=0
c 
        do 1400 i=m0,m
cccc        do 1400 i=1,m
c 
        d=0
        do 1200 j=1,n
c 
        d=d+a(j,i)**2
 1200 continue
c 
        da=d
c 
        d=0
        do 1300 j=1,k
c 
        d=d+b(i,j)**2
 1300 continue
c 
        db=d
c 
        rnorms(i)=sqrt(da*db)
        tot=rnorms(i)+tot
cccc        tot=rnorms(i)
 1400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine rskel_rscap2(b,m,k,ii,jj,prod)
        implicit real *8 (a-h,o-z)
        save
        dimension b(m,k)
c 
        prod=0
        do 1200 i=1,k
        prod=prod+b(ii,i)*b(jj,i)
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine rskel_rscap(x,y,n,prod)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1)
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
        subroutine rskel_arrmove(x,y,n)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1)
c 
        do 1200 i=1,n
        y(i)=x(i)
 1200 continue
        return
        end
  
  
  
