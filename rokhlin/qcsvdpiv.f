        implicit real *16 (a-h,o-z)
        complex *32 a(10 000),u(10 000),
     1      v(10 000),a2(100 000),w(100 000),
     2      a3(10 000)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
         PRINT *, 'ENTER n'
         READ *,n
         CALL PRINF('n=*',n,1 )
c 
         PRINT *, 'ENTER eps2'
         READ *,eps2
         CALL PRINQ('eps2=*',eps2,1 )
c 
c       test the gram-schmidt procedure
c 
        m=n/2+7
        m=n
        lw=300000
        call svdtes(a,n,m,u,v,w,lw,eps2,a2,a3)
  
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine svdtes(a,n,m,u,v,w,lw,eps,a2,a3)
        implicit real *16 (a-h,o-z)
        save
        complex *32 a(n,m),u(n,m),v(m,n),w(1),s(2000),
     1      a2(n,m),w1(10 000),w2(10 000),a3(n,m),ima
        real *16 tt(2000)
        data ima/(0.0d0,1.0d0)/
c 
c        construct the test matrix to be compressed
c 
        h=1
        h=h/n
        do 1100 i=1,2000
        tt(i)=i*h
 1100 continue
c 
        do 1400 i=1,m
        do 1200 j=1,n
ccc        a(j,i)=i**2-j**2
        a(j,i)=qexp(-tt(i)*tt(j))+qexp(-tt(i)/tt(j)) +ima
        a(j,i)=a(j,i)+a(j,i)**2*ima
 1200 continue
 1400 continue
c 
c      copy a into a2
c 
        do 1800 i=1,m
        do 1600 j=1,n
        a2(j,i)=a(j,i)
 1600 continue
 1800 continue
         nprint=100
         if(n**2*2 .lt. nprint) nprint=n**2*2
         call prinq('a as created*',a,nprint)
c 
c        perform SVD of a
c 
         lw=100 000
         call qcsvdpiv(ier,a,n,m,u,v,s,ncols,eps,
     1      w,lw,ltot)
          call prinf('after qcsvdpiv, ltot=*',ltot,1)
          call prinf('after qcsvdpiv, ncols=*',ncols,1)
c 
          call prinf('after qcsvdpiv, ncols=*',ncols,1)
c 
         nprint=100
         if(n*ncols*2 .lt. nprint) nprint=n*ncols*2
          call prinq('after qcsvdpiv, u=*',u,nprint)
          call prinq('after qcsvdpiv, s=*',s,ncols*2)
cccc          call prinq('after qcsvdpiv, v=*',v,ncols*m*2)
c 
c        finally, multiply them things together
c 
         call mulback(u,w1,v,n,ncols,m,a3,w2,s)
  
cccc         call prinq('after mulback, a3=*',a3,n*m*2)
cccc         call prinq('while a2=*',a2,n*m*2)
c 
c        now, check the error
c 
        dmax=0
        do 1860 i=1,n
        do 1840 j=1,m
c 
        dd=cqabs(a3(i,j)-a2(i,j))
        if(dd .gt. dmax) dmax=dd
 1840 continue
 1860 continue
         call prinq('and maximum error=*',dmax,1)
c 
c       test the orthogonality of columns of u
c 
cccc          call prinq('before orttes, u=*',u, n*n)
        call orttes(u,n,ncols,w1)
c 
       return
       end
c 
c 
c 
c 
c 
        subroutine mulback(u,t,w,n,ncols,m,c,work,s)
        implicit real *16 (a-h,o-z)
        save
        complex *32 u(n,ncols),t(ncols,ncols),w(m,ncols),
     1      work(n,ncols),c(n,m),s(1),d
c 
c       multiply u by t
c 
cccc         call prinq('inside mulback, u=*',u,n*ncols)
cccc         call prinq('inside mulback, w=*',w,m*ncols)
c 
        do 1050 i=1,ncols
        do 1040 j=1,ncols
        t(i,j)=0
 1040 continue
        t(i,i)=s(i)
 1050 continue
cccc        call prinq('inside mulback, t=*',t,ncols*ncols)
c 
        do 1600 i=1,n
        do 1400 j=1,ncols
        d=0
        do 1200 k=1,ncols
        d=d+u(i,k)*t(k,j)
ccc        d=d+u(i,k)*t(j,k)
 1200 continue
        work(i,j)=d
 1400 continue
 1600 continue
c 
c       multiply the result by w^*
c 
        do 2600 i=1,n
        do 2400 j=1,m
        d=0
        do 2200 k=1,ncols
        d=d+work(i,k)*qconjg(w(j,k))
 2200 continue
        c(i,j)=d
 2400 continue
 2600 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine orttes(u,n,ncols,w)
        implicit real *16 (a-h,o-z)
        save
        complex *32 u(n,ncols),w(ncols,ncols)
c 
c       test the orthogonality of columns of u
c 
        do 1880 i=1,ncols
        do 1870 j=1,ncols
        call qcsvdscapr(u(1,i),u(1,j),n,w(i,j))
 1870 continue
 1880 continue
cccc         call prinq('matrix of inner products*',w,ncols**2*2)
c 
c       now, check how non-orthogonal it is
c 
        done=1
        do 2400 i=1,ncols
        do 2200 j=1,ncols
        if( (i .ne. j) .and. (dd .lt. cqabs(w(i,j))) )
     1   dd=cqabs(w(i,j) )
c 
        if( (i .eq. j) .and. (dd .lt. cqabs(w(i,j)-done)) )
     1    dd=cqabs(w(i,j)-done)
 2200 continue
 2400 continue
        call prinq('and non-orthogonality of u is*',dd,1)
        return
        end
c 
c 
c 
c 
c 
        subroutine matmua(a,b,n,m,ncols,c)
        implicit real *16 (a-h,o-z)
        save
        complex *32 a(n,ncols),b(m,ncols),c(n,m),cd
c 
        do 1600 i=1,n
        do 1400 j=1,m
c 
        cd=0
        do 1200 k=1,ncols
        cd=cd+a(i,k)*qconjg(b(j,k))
 1200 continue
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c          this is the end of the debugging output and the beginning
c          of the actual SVD code
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c 
c 
        subroutine qcsvdpiv(ier,a,n,m,u,v,s,ncols,eps,
     1      w,lw,ltot)
        implicit real *16 (a-h,o-z)
        save
        complex *32 a(n,m),u(n,1),v(m,1),s(1),w(1)
c 
c        this subroutine constructs the singular value decomposition
c        of the usert-specified matrix a, under the assumption
c        that the rank of a is considerably smaller than its
c        dimensionalities n,m. the output of this subroutine
c        consists of two matrices u,v, and a vector s such that
c 
c                 a = u d v^*                                   (1)
c 
c        with t the diagonal matrix with the vector s on the
c        diagonal, and u, v two matrices whose columns are
c        orthonormal. the reason for the existence of this
c        subroutine is the fact that on exit, the matrices u,v
c        have dimensionalities (n,ncols), (m,ncols), respectively,
c        with ncols the numerical rank of a.
c 
c                   input parameters:
c 
c  a - the matrix to be svd-decomposed; note that it is NOT damaged
c       by this subroutine in any way
c  n,m - dimensionalities of a
c  eps - the accuracy with which the calculation will be performed;
c        should not be less than the machine precision; recommended
c        value: 1.0d-16
c  lw - the number of real *16 words provided by the user in the work
c        array w. must be at least m+n+10 (otherwise, the subroutine
c        can not even estimate the actual memory requirements).
c 
c                   output parameters:
c 
c  ier - error return code
c      ier=0 means successful execution
c      ier=4 means that the amount of space lw allocated supplied
c            in the work array w is insufficient
c      ier=1024 means that the amount of space lw allocated supplied
c            in the work array w is grossly insufficient
c  u - the matrix in (1). note that the matrix u has to be dimensioned
c        u(n,m), but on exit from this subroutine only the first
c        ncols columns of  b  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit b is dimensioned
c        b(n,ncols)
c  v - the second matrix in (1). it is dimensioned v(m,ncols) on
c        exit. on entry the parameter ncols is unknown,
c        so v should be dimensioned v(m,n) by the user. however,
c        on exit only the first ncols columns of v are meaningful
c  s - the array of singular values of a; note that this is a
c        COMPLEX *32 array; imaginary parts of its elements are
c        (hopefully) equal to zero.
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices b, v on exit
c  ltot - the total amount of space (in terms of complex *32 words)
c        needed by the subroutine in the array w
c 
c                   work arrays:
c 
c  w - must be sufficiently long. 3*n*m+ 2*(min(m,n))**2 +2*n+2*m+100
c        are always sufficient. however, most of the time it is much
c        less than that, and depends on the rank of a. the actual
c        amount of space available in w is communicated to the
c        subroutine via the parameter lw (see above), and if this
c        amount is insufficient, the output error code ier is set
c        accordingly, and the memory requirements of the subroutine
c        are returned in the parameter ltot (see above).
c 
c        . . . using gram-schmidt process with pivoting, decompose
c              the matrix a in the form
c 
c                    a=U  V^*,
c 
c              with  u an orthogonal matrix of minimum rank,
c              and v whatever it wants to be
c 
        if( (lw .gt. m) .and. (lw. gt. n)) goto 1200
        ier=1024
        return
 1200 continue
c 
        ifpivot=1
        call qcsvdpiv0(a,n,m,u,v,ncols,w,eps,ifpivot)
cccc         call prinq('after first qcsvdpiv0, rnorms=*',w,ncols)
c 
c       allocate memory for the subsequent operations
c 
        iw=1
        liw=m*ncols+10
c 
        iu2=iw+liw
        lu2=ncols**2+10
c 
        iv2=iu2+lu2
        lv2=ncols**2+10
c 
        it=iv2+lv2
        lt=lv2
c 
        irnorms=it+lt
        lrnorms=n
        if(m .gt. n) lrnorms=m
        lrnorms=lrnorms+10
c 
        iww=irnorms+lrnorms
        lww=n
        if(m .gt. n) lww=m
        lww=lww+10
c 
        ltot=iww+lww
        if(ltot .le. lw) goto 1400
        ier=4
        return
 1400 continue
c 
c       conclude the construction of the singular value decomposition
c 
        iffirst=1
        call qcsvdpiv1(a,u,v,w(iw),w(it),n,m,ncols,
     1    w(irnorms),eps,w(iv2),w(iu2),w(iww),iffirst)
cccc         call prinq('after qcsvdpiv1, w(it)=*',w(it),ncols)
c 
c       copy the right singular matrix into the array v
c 
        jj=0
        do 1800 i=1,ncols
        do 1600 j=1,m
        jj=jj+1
        v(j,i)=w(jj)
 1600 continue
 1800 continue
c 
c       now, copy the singular values into the array s
c 
        do 2000 i=1,ncols
        s(i)=w(it+i-1)
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine qcsvdpiv1(a,u,v,w,t,n,m,ncols,rnorms,
     1    eps2,v2,u2,ww,iffirst)
        implicit real *16 (a-h,o-z)
        save
        complex *32 a(n,m),u(n,1),v(m,1),w(m,1),t(1),
     1      v2(1),u2(1),ww(1)
        dimension rnorms(1)
c 
c        using gram-schmidt process with pivoting, decompose
c        the matrix a in the form
c 
c          a=U  V^*,
c 
c        with  u an orthogonal matrix of minimum rank,
c        and v whatever it wants to be
c 
        if(iffirst .eq. 0) goto 1200
        ifpivot=1
        call qcsvdpiv0(a,n,m,u,v,ncols,rnorms,eps2,ifpivot)
cccc         call prinq('after first qcsvdpiv0, rnorms=*',rnorms,ncols)
 1200 continue
c 
c        using gram-schmidt process without pivoting, decompose
c        the matrix v in the form
c 
c          v=w t^*,
c 
c        with  w an orthogonal matrix of minimum rank,
c        and t a triangular matrix of dimensionality ncols * ncols
  
c 
        ifpivot=0
        call qcsvdpiv0(v,m,ncols,w,t,ncols2,rnorms,eps2,ifpivot)
cccc         call prinq('after second qcsvdpiv0, rnorms=*',rnorms,ncols)
c 
c        use jacobi rotations to sv-decompose the matrix t
c 
        ifuv=1
c 
cccc        subroutine qd2mcudv(a,n,u,v,eps,ifuv,numrows)
  
cccc        call prinq('before qd2mcudv, t=*',t,ncols*ncols*2)
c 
        call qd2mcudv(t,ncols,u2,v2,eps2,ifuv,numrows)
cccc        call prinq('after d2mudv, t=*',t,ncols)
cccc        call prinf('after d2mudv, numrows=*',numrows,1)
c 
c       multiply the matrix u from the right by u2, and the
c       matrix w from the right by v2^*
c 
cccc         call prinq('after qd2mcudv, t=*',t,n*n*2)
        call qcsvdfinm(u,u2,w,v2,n,m,ncols,ww)
        return
        end
c 
c 
c 
c 
c 
        subroutine qcsvdfinm(u,u2,w,v2,n,m,ncols,ww)
        implicit real *16 (a-h,o-z)
        save
        complex *32 u(n,ncols),w(m,ncols),u2(ncols,ncols),
     1      v2(ncols,ncols),ww(1),d
c 
c       multiply u from the right by u2
c 
        do 1600 i=1,n
        do 1400 j=1,ncols
c 
        d=0
        do 1200 k=1,ncols
        d=d+u(i,k)*u2(k,j)
 1200 continue
        ww(j)=d
 1400 continue
        do 1500 j=1,ncols
        u(i,j)=ww(j)
 1500 continue
 1600 continue
c 
c         multiply w from the right by the adjoint of v2
c 
        zero=0
        do 2600 i=1,m
        do 2400 j=1,ncols
c 
        d=zero
        do 2200 k=1,ncols
        d=d+w(i,k)*qconjg(v2(j,k))
 2200 continue
        ww(j)=d
 2400 continue
        do 2500 j=1,ncols
        w(i,j)=ww(j)
 2500 continue
 2600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine qcsvdpiv0(a,n,m,b,v,ncols,rnorms,eps,ifpivot)
        implicit real *16 (a-h,o-z)
        save
        complex *32 a(n,m),b(n,m),v(m,1),prod
        dimension rnorms(1)
c 
c       this matrix applies the compressing gram-schmidt
c       process to the matrix a, obtaining its decomposition
c       in the form
c 
c             a=b v^T,                                        (1)
c 
c       with the matrices b, v having dimensionalities
c       b(n,ncols), v(m,ncols), respectively. the reason for the
c       existence of this subroutine is the hope that the
c       dimensionality ncols, determined by this subroutine,
c       is comsiderably lower than either n or m
c 
c                     input parameters:
c 
c  a - the matrix to be decomposed (note that a is NOT damaged
c       by this subroutine in any way)
c  n,m - the dimensionalities of a
c  eps - the accuracy to which the decomposition (1) will be computed.
c 
c                     output parameters:
c 
c  b - the matrix in (1). note that the matrix b has to be dimensioned
c        b(n,m), but on exit from this subroutine only the first
c        ncols columns of  b  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit b is dimensioned
c        b(n,ncols)
c  v - the second matrix in (1). it is dimensioned v(m,ncols) on
c        exit. on entry the parameter ncols is unknown,
c        so v should be dimensioned v(m,n) by the user. however,
c        on exit only the first ncols columns of v are meaningful
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices b, v on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *16 elements long.
c 
c        . . . copy the user-supplied matrix a into b
c 
        do 1400 i=1,m
        do 1200 j=1,n
        b(j,i)=a(j,i)
 1200 continue
 1400 continue
c 
c        apply the gram-schmidt proces (with pivoting) to b
c 
         call qcsvdgrm(b,n,m,rnorms,eps,ncols,ifpivot)
c 
c        project the original matrix on the obtained orthogonal
c        columns
c 
        do 2400 i=1,ncols
        do 2200 j=1,m
c 
        call qcsvdscapr(a(1,j),b(1,i),n,prod)
cccc        v(j,i)=prod
        v(j,i)=qconjg(prod)
 2200 continue
 2400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine qcsvdgrm(b,n,m,rnorms,eps,ncols,ifpivot)
        implicit real *16 (a-h,o-z)
        save
        complex *32 b(n,m),cd
        dimension rnorms(1)
c 
c       this subroutine applies a pivoted double gram-schmidt
c       procedure to the matrix b, resulting in a collection of
c       orthogonal vectors. the number of these vectors is
c       equal to the numerical rank of b, given the precision
c       of computations eps, which is supplied by the user.
c 
c                    input paramneters:
c 
c  b - the matrix to be gram-schmidt'ed. it is destroyed by this
c       subroutine
c  n,m - dimensionalities of the matrix b
c  eps - the machine accuracy
c 
c                     output parameters:
c 
c  b - the matrix of gram-schmidt vectors of the matrix a. note
c        that on exit from this subroutine only the first
c        ncols columns of  b  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit b is dimensioned
c        b(n,ncols)
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrix b on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *16 elements long.
c 
c       . . . prepare the array of values of norms
c             of columns
c 
        done=1
        dtot=0
        do 1400 i=1,m
        d=0
        do 1200 j=1,n
cccc        d=d+b(j,i)**2
        d=d+b(j,i)*qconjg(b(j,i))
 1200 continue
        rnorms(i)=d
        dtot=dtot+d
 1400 continue
c 
c       . . . conduct gram-schmidt iterations
c 
        thresh=dtot*eps**2
        do 4000 i=1,m
c 
c       find the pivot
c 
         if(ifpivot .eq. 0) goto 2700
        ipivot=i
        rn=rnorms(i)
c 
        do 2200 j=i+1,m
        if(rnorms(j) .le. rn) goto 2200
        rn=rnorms(j)
        ipivot=j
 2200 continue
 2400 continue
c 
c       put the column number ipivot in the i-th place
c 
        do 2600 j=1,n
        cd=b(j,i)
        b(j,i)=b(j,ipivot)
        b(j,ipivot)=cd
 2600 continue
c 
        d=rnorms(ipivot)
        rnorms(ipivot)=rnorms(i)
        rnorms(i)=d
 2700 continue
c 
c       orthogonalize the i-th column to all preceeding ones
c 
        if(i .eq. 1) goto 2790
        do 2780 j=1,i-1
c 
        call qcsvdscapr(b(1,i),b(1,j),n,cd)
c 
cccc        cd=qconjg(cd)
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*cd
 2770 continue
 2780 continue
 2790 continue
c 
c       normalize the i-th column
c 
        call qcsvdscapr(b(1,i),b(1,i),n,cd)
c 
c       if the norm of the remaining part of the matrix
c       is sufficiently small - terminate the process
c 
        d=cd
        if( (ifpivot .eq. 1) .and.
     1      (d .lt. thresh) ) return
        ncols=i
c 
        d=done/qsqrt(d)
        do 2800 j=1,n
        b(j,i)=b(j,i)*d
 2800 continue
c 
c        orthogonalize everything else to it
c 
        do 3200 j=i+1,m
c 
        if(rnorms(j) .lt. thresh) goto 3200
c 
        call qcsvdscapr(b(1,i),b(1,j),n,cd)
        cd=qconjg(cd)
c 
         rrn=0
        do 3000 l=1,n
        b(l,j)=b(l,j)-b(l,i)*cd
        rrn=rrn+b(l,j)*qconjg(b(l,j))
 3000 continue
        rnorms(j)=qsqrt(rrn)
 3200 continue
 3400 continue
 4000 continue
  
 4200 continue
  
  
        return
        end
c 
c 
c 
c 
c 
        subroutine qcsvdscapr(x,y,n,prod)
        implicit complex *32 (a-h,o-z)
        save
        dimension x(1),y(1)
c 
        prod=0
        do 1200 i=1,n
        prod=prod+x(i)*qconjg(y(i))
 1200 continue
        return
        end
c 
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c            this is the end of the debuging code and the beginning
c            of the actual SVD routines
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c 
        subroutine qd2mcudv(a,n,u,v,eps,ifuv,numrows)
        implicit complex *32 (a-h,o-z)
        save
        dimension a(n,n),u0(2,2),v0(2,2),u(n,n),v(n,n),
     1      b(2,2),rlam(2),u1(2,2),v1(2,2)
        real *16 rlam(2),eps,thresh
c 
c        this subroutine uses the classical jacobi iteration
c        to construct the singular value decomposition of a
c        general complex *32 matrix, of the form
c 
c          a = u s v,                                        (1)
c 
c        with u,v unitary matrices, and s diagonal, with
c        the diagonal elements monotonically decreasing.
c 
c                     input parameters:
c 
c  a - the n * n matrix to be sv-decomposed, is destroyed
c        by the subroutine
c  n - the dimensionality of a
c  eps - the accuracy to which the computations will be performed
c  ifuv - tells the subroutine whether the singular vectors are
c        to be computed, or only singular values
c 
c                     output parameters:
c 
c  a - the first column of a contains the (ordered) singular values ,
c         {\bf IN COMPLEX FORM}
c  u,v - orthogonal matrices in (1)
c  numrows - the number of jacobi iterations (i.e. single rotations)
c        that the algorithm has taken
c 
c       . . . one sweep after another, conduct jacobi
c             iterations
c 
        maxsweep=1000
        if(ifuv .eq. 0) goto 1150
        zero=0
        done=1
        do 1100 i=1,n
        do 1050 j=1,n
        u(j,i)=zero
        v(j,i)=zero
 1050 continue
        u(i,i)=done
        v(i,i)=done
 1100 continue
 1150 continue
        numrows=0
        thresh=qsqrt(eps)
        thresh=qsqrt(thresh)
        thresh=qsqrt(thresh)
        do 3400 iijjkk=1,4
c 
        do 3000 ijk=1,maxsweep
        ifany=0
c 
        do 2000 i=1,n
        do 1800 j=1,n
         if(j .eq. i) goto 1800
         if( (cqabs(a(i,j)) .lt. thresh) .and.
     1       (cqabs(a(j,i)) .lt. thresh) ) goto 1800
        ifany=1
        numrows=numrows+1
c 
c       construct the svd of the (i,j)-th 2 * 2 submatrix
c 
        b(1,1)=a(i,i)
        b(1,2)=a(i,j)
        b(2,1)=a(j,i)
        b(2,2)=a(j,j)
        call qconeuv(ier,b,rlam,u0,v0)
c 
c     transform the big matrix
c 
c       . . . the rows
c 
        do 1200 l=1,n
        di=u0(1,1)*a(i,l)+u0(1,2)*a(j,l)
        dj=u0(2,1)*a(i,l)+u0(2,2)*a(j,l)
        a(i,l)=di
        a(j,l)=dj
 1200 continue
c 
c       . . . and the columns
c 
        do 1400 l=1,n
        di=v0(1,1)*a(l,i)+v0(2,1)*a(l,j)
        dj=v0(1,2)*a(l,i)+v0(2,2)*a(l,j)
        a(l,i)=di
        a(l,j)=dj
 1400 continue
c 
c       if the user so requested, adjust the matrices u,v
c 
        if(ifuv .eq. 0) goto 1800
c 
        call qd2mcpon(u0,u1)
        call qd2mcpon(v0,v1)
        do 1600 l=1,n
c 
c       the matrix u
c 
        di=u1(1,1)*u(l,i)+u1(2,1)*u(l,j)
        dj=u1(1,2)*u(l,i)+u1(2,2)*u(l,j)
        u(l,i)=di
        u(l,j)=dj
c 
c       . . . and the columns
c 
        di=v1(1,1)*v(i,l)+v1(1,2)*v(j,l)
        dj=v1(2,1)*v(i,l)+v1(2,2)*v(j,l)
        v(i,l)=di
        v(j,l)=dj
 1600 continue
c 
 1800 continue
 2000 continue
c 
         call prinf('ifany=*',ifany,1)
      if(ifany .eq. 0) goto 3200
 3000 continue
 3200 continue
        thresh=thresh**2
 3400 continue
c 
c        now, reorder the singular values of a to put them
c        in increasing order
c 
        do 4200 i=1,n
        a(i,1)=a(i,i)
 4200 continue
c 
       call qd2mcsrt(a(1,1),n,u,v,a(1,2),ifuv)
        return
        end
c 
c 
c 
c 
c 
        subroutine qd2mcsrt(s,n,u,v,iw,ifuv)
        implicit real *16 (a-h,o-z)
        save
        complex *32 u(n,n),v(n,n),s(1),com,cd
        dimension rea(2),iw(1)
        equivalence (rea(1),com)
c 
c        if some of the alleged singular values
c        are negative - change their sign, and that
c        of the corresponding left eigenvectors
c 
        do 1100 i=1,n
        com=s(i)
        ds=rea(1)
        if(ds .ge. 0) goto 1100
        s(i)=-s(i)
        if(ifuv .eq. 0) goto 1100
        do 1050 j=1,n
        v(i,j)=-v(i,j)
 1050 continue
 1100 continue
  
c        using the bubble sort, reorder the elements
c        of array s
c 
        do 1200 i=1,n
        iw(i)=i
 1200 continue
c 
        do 2000 i=1,n-1
c 
        ifmoved=0
        do 1600 j=1,n-i
        com=s(j)
        dsj=rea(1)
        com=s(j+1)
        dsjp1=rea(1)
        if(dsj. ge. dsjp1) goto 1600
        jj=iw(j)
        iw(j)=iw(j+1)
        iw(j+1)=jj
c 
        cd=s(j)
        s(j)=s(j+1)
        s(j+1)=cd
 1600 continue
 2000 continue
        if(ifuv .eq. 0) return
c 
c       now, reorder the rows of u and columns of v
c 
        do 2400 i=1,n-1
        ii=iw(i)
c 
c       exchange the rows number i and ii in the matrix u
c 
        do 2200 j=1,n
        cd=v(i,j)
        v(i,j)=v(ii,j)
        v(ii,j)=cd
c 
c       exchange columns number i and ii in the matrix v
c 
        cd=u(j,i)
        u(j,i)=u(j,ii)
        u(j,ii)=cd
 2200 continue
c 
        do 2300 j=i+1,n
        if(iw(j) .ne. i) goto 2300
        iw(j)=ii
        goto 2350
 2300 continue
 2350 continue
c 
 2400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine qconeuv(ier,a0,rlam,u,v)
        implicit real *16 (a-h,o-z)
        save
        complex *32 z,b(2,2),u(2,2),a0(2,2),
     1      uu(2,2),v(2,2),w1(2,2),w2(2,2),
     2      w3(2,2),cd,clam1,clam2,clam,x(2),y(2)
        dimension rlam(2)
c 
c       multiply a by its adjoint
c 
        call qd2mcpon(a0,w1)
        call qd2mprodc(w1,a0,b)
c 
c       find the eigenvalues of a^*a
c 
        two=2
        four=4
        cd=(b(1,1)-b(2,2))**2+four*b(1,2)*b(2,1)
        cd=cqsqrt(cd)
        clam1= (b(1,1)+b(2,2)+cd )/two
        clam2= (b(1,1)+b(2,2)-cd )/two
c 
c        find the eigenvectors
c 
        clam=clam1
        dc2=clam2*qconjg(clam2)
        dc1=clam1*qconjg(clam1)
        if( dc2 .gt. dc1 ) clam=clam2
c 
        d11=(b(1,1)-clam)*qconjg((b(1,1)-clam))
        d22=(b(2,2)-clam)*qconjg((b(2,2)-clam))
        d12=b(1,2)*qconjg(b(2,2))
c 
        if( (d11 .lt. d22) .or. (d11. lt. d12) ) goto 1400
c 
        x(2)=1
        x(1)=-b(1,2)/(b(1,1)-clam)
        goto 1800
c 
 1400 continue
c 
        if(d22 .lt. d12) goto 1600
c 
        x(1)=1
        x(2)=-b(2,1)/(b(2,2)-clam)
        goto 1800
c 
 1600 continue
c 
        x(2)=1
        x(1)=(b(1,1)-clam)/b(1,2)
c 
 1800 continue
c 
        d=x(1)*qconjg(x(1))+x(2)*qconjg(x(2))
        d=qsqrt(d)
        x(1)=x(1)/d
        x(2)=x(2)/d
c 
        call qd2md2ort(x,y)
c 
        uu(1,1)=x(1)
        uu(2,1)=x(2)
        uu(1,2)=y(1)
        uu(2,2)=y(2)
c 
        call qd2mcpon(uu,w2)
c 
c        now, obtain v
c 
         call qd2mprodc(a0,uu,v)
c 
         d1=v(1,1)*qconjg(v(1,1))+v(2,1)*qconjg(v(2,1))
         d2=v(1,2)*qconjg(v(1,2))+v(2,2)*qconjg(v(2,2))
c 
         if(d2 .gt. d1) goto 3000
c 
         d=qsqrt(d1)
         call qd2md2ort(v(1,1),v(1,2))
         goto 3200
c 
 3000 continue
c 
         d=qsqrt(d2)
         call qd2md2ort(v(1,2),v(1,1))
 3200 continue
c 
c       obtain the diagonal matrix
c 
         call qd2mprodc(a0,uu,w1)
         call qd2mcpon(v,b)
         call qd2mprodc(b,w1,w3)
c 
        call qd2mcopy2(w2,u)
c 
c       finally, make sure that the diagonal elements are real
c 
        z=1
        dzero=0
        d=w3(1,1)*qconjg(w3(1,1))
        if(d .ne. dzero) z=w3(1,1)/cqabs(w3(1,1))
        rlam(1)=cqabs(w3(1,1))
        v(1,1)=v(1,1)*z
        v(2,1)=v(2,1)*z
c 
        z=1
        dzero=0
        d=w3(2,2)*qconjg(w3(2,2))
        if(d .ne. dzero) z=w3(2,2)/cqabs(w3(2,2))
        rlam(2)=cqabs(w3(2,2))
        v(1,2)=v(1,2)*z
        v(2,2)=v(2,2)*z
c 
c        invert u, v
c 
        call qd2mcopy2(v,w3)
        call qd2mcpon(u,v)
        call qd2mcpon(w3,u)
        return
        end
c 
c 
c 
c 
c 
        subroutine qd2mcopy2(a,b)
        implicit complex *32 (a-h,o-z)
        save
        dimension a(2,2),b(2,2)
c 
        b(1,1)=a(1,1)
        b(2,2)=a(2,2)
        b(1,2)=a(1,2)
        b(2,1)=a(2,1)
        return
c 
c 
c 
c 
        entry qd2mcpon(a,b)
c 
        b(1,1)=qconjg(a(1,1))
        b(2,2)=qconjg(a(2,2))
        b(2,1)=qconjg(a(1,2))
        b(1,2)=qconjg(a(2,1))
        return
        end
c 
c 
c 
c 
c 
        subroutine qd2md2ort(x,y)
        implicit real *16 (a-h,o-z)
        save
        complex *32 x(2),y(2)
c 
c       normalize x
c 
        d1=x(1)*qconjg(x(1))
        d2=x(2)*qconjg(x(2))
        d=d1+d2
        d=qsqrt(d)
        x(1)=x(1)/d
        x(2)=x(2)/d
c 
c       make y into a unit vector orthogonal to x
c 
c       . . . when |x(1)| > |x(2)|
c 
        if(d1 .lt. d2) goto 1200
        y(2)=x(1)
        y(1)=-qconjg(x(2))*x(1)/qconjg(x(1))
        return
 1200 continue
c 
c       . . . when |x(2)| > |x(1)|
c 
        y(1)=x(2)
        y(2)=-qconjg(x(1))*x(2)/qconjg(x(2))
        return
        end
c 
c 
c 
c 
c 
        subroutine qd2mprodc(a,b,c)
        implicit complex *32 (a-h,o-z)
        save
        dimension a(2,2),b(2,2),c(2,2)
c 
        c(1,1)=a(1,1)*b(1,1)+a(1,2)*b(2,1)
        c(1,2)=a(1,1)*b(1,2)+a(1,2)*b(2,2)
        c(2,1)=a(2,1)*b(1,1)+a(2,2)*b(2,1)
        c(2,2)=a(2,1)*b(1,2)+a(2,2)*b(2,2)
        return
        end
