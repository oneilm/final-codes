        implicit real *8 (a-h,o-z)
        dimension a(500 000),u(500 000),
     1      v(500 000),a2(500 000),w(1000 000),
     2      a3(500 000)
c 
cccc        dimension a(50 000),u(50 000),
cccc     1      v(50 000),a2(50 000),w(100 000),
cccc     2      a3(50 000)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
         PRINT *, 'ENTER n'
         READ *,n
         CALL PRINF('n=*',n,1 )
c 
         PRINT *, 'ENTER log10(eps)'
         READ *,d
         CALL PRIN2('log10(eps)=*',d,1 )
c 
         eps=1.0d0/10.0d0**d
c 
c       test the gram-schmidt procedure
c 
        m=n/2+7
        m=n
        lw=1000 000
        call svdtes(a,n,m,u,v,w,lw,eps,a2,a3)
  
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine svdtes(a,n,m,u,v,w,lw,eps,a2,a3)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension a(n,m),u(n,m),v(m,n),w(1),s(200 000),
     1      a2(n,m),w1(500 000),w2(500 000),a3(n,m),tt(200000)
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
        a(j,i)=i**2-j**2
        a(j,i)=dexp(-tt(i)*tt(j))+dexp(-tt(i)/tt(j))
 1200 continue
c 
        a(i,i)=10000
c 
 1400 continue
c 
c      copy a into a2
c 
        do 1800 i=1,m
        do 1600 j=1,n
        a2(j,i)=a(j,i)
 1600 continue
 1800 continue
c 
c        perform svd of the matrix a
c 
        call fvdpivot(ier,a,n,m,u,v,s,ncols,eps,
     1      w,lw,ltot)
        call prinf('after fvdpivot, ier=*',ier,1)
        call prinf('after fvdpivot, ncols=*',ncols,1)
        call prinf('after fvdpivot, ltot=*',ltot,1)
c 
c       test the accuracy
c 
        call mulback(u,w1,v,n,ncols,m,a3,w2,s)
ccc        call prin2('a3=*',a3,n*m)
ccc        call prin2('while a2=*',a2,n*m)
c 
        err=0
        do 2200 i=1,m
        do 2000 j=1,n
        d=d+(a2(j,i)-a3(j,i))**2
 2000 continue
 2200 continue
        call prin2('and the error is*',dsqrt(d),1)
c 
c        look at the inner products of the columns of
c        the matrix u
c 
        ii=0
        do 3400 i=1,ncols
        do 3200 j=1,ncols
c 
        ii=ii+1
        call fvdscapr(u(1,j),u(1,i),n,w1(ii))
 3200 continue
 3400 continue
  
        call prin2('and inner products of columns of u are*',
     1        w1,ii)
c 
        call prin2('and again, the error is*',dsqrt(d),1)
  
       return
       end
c 
c 
c 
c 
c 
        subroutine mulback(u,t,w,n,ncols,m,c,work,s)
        implicit real *8 (a-h,o-z)
        save
        dimension u(n,ncols),t(ncols,ncols),w(m,ncols),
     1      work(n,ncols),c(n,m),s(1)
c 
c       multiply u by t
c 
cccc         call prin2('inside mulback, u=*',u,n*ncols)
cccc         call prin2('inside mulback, w=*',w,m*ncols)
c 
        do 1050 i=1,ncols
        do 1040 j=1,ncols
        t(i,j)=0
 1040 continue
        t(i,i)=s(i)
 1050 continue
cccc        call prin2('inside mulback, t=*',t,ncols*ncols)
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
        d=d+work(i,k)*w(j,k)
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c          this is the end of the debugging output and the beginning
c          of the actual FVD code
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c 
c 
        subroutine fvdpivot(ier,a,n,m,u,v,s,ncols,eps,
     1      w,lw,ltot)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),u(n,1),v(m,1),s(1),w(1)
c 
c        this subroutine constructs a distorted singular value
c        decomposition (to be referred to as the FVD) of the
c        user-specified matrix A. This subroutine is particularly
c        efficient when the rank of a is considerably smaller than
c        its dimensionalities n,m. The output of this subroutine
c        consists of two matrices u,v, and a vector s such that
c 
c                 a = u d v^*                                   (1)
c 
c        with d the diagonal matrix with the vector s on the
c        diagonal, V a matrix all of whose columns are orthonormal,
c        and U a matrix all of whose singular values are reasonably
c        close to 1. The reason for the existence of this
c        subroutine is the fact that on exit, the matrices u,v
c        have dimensionalities (n,ncols), (m,ncols), respectively,
c        with ncols the numerical rank of a.
c 
c        Please note that the output of this subroutine is not
c        quite a singular value decomposition of A. On the other
c        hand, it is a perfectly acceptable substitute for the SVD
c        in the context of least squares problems, and is much
c        cheaper to evaluate than the "true" SVD.
c 
c                   input parameters:
c 
c  a - the matrix to be svd-decomposed; note that it is NOT damaged
c       by this subroutine in any way
c  n,m - dimensionalities of a
c  eps - the accuracy with which the calculation will be performed;
c        should not be less than the machine precision; recommended
c        value: 1.0d-16
c  lw - the number of real *8 words provided by the user in the work
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
c      ier=4096 means that the subroutine failed to separate two
c            singular values after 30 different randomization.
c            Something is seriously wrong.
c  u - the matrix in (1). note that the matrix u has to be dimensioned
c        u(n,m), but on exit from this subroutine only the first
c        ncols columns of  b  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit b is dimensioned
c        b(n,ncols)
c  v - the second matrix in (1). it is dimensioned v(m,ncols) on
c        exit. on entry the parameter ncols is unknown,
c        so v should be dimensioned v(m,n) by the user. however,
c        on exit only the first ncols columns of v are meaningful
c  s - the array of singular values of a
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices b, v on exit
c  ltot - the total amount of space (in terms of real *8 words)
c        needed by the subroutine in the array w
c 
c                   work arrays:
c 
c  w - must be sufficiently long. 3*n*m+ 2*(min(m,n))**2 +2*n+2*m+100
c        are always sufficient. however, must of the time it is much
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
        call fvdpiv0(a,n,m,u,v,ncols,w,eps,ifpivot)
cccc         call prin2('after first fvdpiv0, rnorms=*',w,ncols)
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
        lenww=lw-iww-6
c 
c       conclude the construction of the singular value decomposition
c 
        call fvdpiv1(ier,a,u,v,w(iw),w(it),n,m,ncols,
     1    w(irnorms),eps,w(iv2),w(iu2),w(iww),lenww,lused,s)
c 
        ltot=iww+lused
c 
        if(ier .ne. 0) return
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
        return
        end
c 
c 
c 
c 
c 
        subroutine fvdpiv1(ier,a,u,v,w,t,n,m,ncols,rnorms,
     1    eps2,v2,u2,ww,lenww,lused,rlams)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),u(n,1),v(m,1),w(m,1),t(1),rnorms(1),
     1      v2(1),u2(1),ww(1),rlams(1)
c 
c        using gram-schmidt process without pivoting, decompose
c        the matrix v in the form
c 
c          v=w t^*,
c 
c        with  w an orthogonal matrix of minimum rank,
c        and t a triangular matrix of dimensionality ncols * ncols
c 
        ier=0
        ifpivot=0
        call fvdpiv0(v,m,ncols,w,t,ncols2,rnorms,eps2,ifpivot)
c 
c        use inverse power to sv-decompose the matrix t
c 
        iaa=1
        laa=ncols+10
c 
        ibb=iaa+laa
        lbb=ncols+10
c 
        iaaold=ibb+lbb
        laaold=ncols+10
c 
        ibbold=iaaold+laaold
        lbbold=ncols+10
c 
        iyrand=ibbold+lbbold
        lyrand=ncols+10
c 
        iwww=iyrand+lyrand
        lenwww=lenww-iwww-5
c 
        call fvd_afterpiv(ier,t,ncols,u2,v2,ww(iaa),ww(ibb),
     1      eps2,v,ww(iwww),lenwww,lused,rlams,ww(iyrand),
     2      ww(iaaold),ww(ibbold) )
c 
        lused=iwww+lused
c 
        if(ier .ne. 0) return
c 
        call fvdfinm(u,u2,w,v2,n,m,ncols,ww)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine fvd_afterpiv(ier,t,ncols,u2,v2,aa,bb,eps,
     1      v3,www,lenwww,lused,rlams,yrand,aaold,bbold)
        implicit real *8 (a-h,o-z)
        save
        dimension v2(ncols,1),u2(ncols,1),t(ncols,ncols),
     1      aa(1),bb(1),www(1),v3(1),rlams(1),yrand(1),
     2      aaold(1),bbold(1)
c 
c       bidiagonalize the input triangular matrix
c 
        ier=0
c 
        call fvd_bidiag(t,ncols,u2,v2)
c 
        do 1200 i=1,ncols-1
c 
        aaold(i)=t(i,i)
        bbold(i)=t(i+1,i)
 1200 continue
c 
        aaold(ncols)=t(ncols,ncols)
c 
c        conduct outer iterations designed to reduce
c        the probability of degenerate singular values
c        of the bidiagonal matrix
c 
        do 3000 ijk=1,30
c 
c       multiply the rows of the bidiagonal matrix by random
c       numbers between 1 and 2, and scale the matrix u2
c       by the inverse of the said matrix
c 
cccc        call bidiag_rand(ncols,yrand)
        call bidiag_rand(ncols,yrand)
c 
        do 1400 i=1,ncols
c 
        yrand(i)=1+yrand(i)
cccc        yrand(i)=1
 1400 continue
c 
        do 1600 i=1,ncols-1
c 
        aa(i+1)=yrand(i+1)*aaold(i+1)
        bb(i)=yrand(i+1)*bbold(i)
 1600 continue
c 
        aa(1)=aaold(1)*yrand(1)
c 
        do 2400 i=1,ncols
        do 2200 j=1,ncols
c 
        u2(j,i)=u2(j,i)/yrand(i)
 2200 continue
 2400 continue
c 
        eps3=1.0d-13
c 
        call bidiag_svd(jer,aa,bb,ncols,eps3,t,v3,rlams,nlams,
     1      www,lenwww,lused)
c 
        if(jer .eq. 0) goto 3200
c 
        if(jer .gt. 4) then
            ier=32
            return
        endif
c 
 3000 continue
c 
        ier=4096
        return
 3200 continue
c 
c       multiply the obtained matrices of singular vectors by the
c       matrices u2, v2
c 
        call matmult(u2,t,www,ncols)
        call fvdcopy(www,u2,ncols**2)
        call fvd_matamult(v3,v2,www,ncols)
        call fvdcopy(www,v2,ncols**2)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine fvd_matamult(a,x,c,n)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),x(n,n),c(n,n)
c 
        do 1600 i=1,n
        do 1400 j=1,n
c 
        cd=0
        do 1200 ll=1,n
c 
        cd=cd+a(ll,i)*x(ll,j)
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
        subroutine fvdcopy(a,b,n)
        implicit real *8 (a-h,o-z)
        save
        dimension a(1),b(1)
c 
        do 1200 i=1,n
        b(i)=a(i)
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine fvd_bidiag(a,n,u,v)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),b(3),ub(2,2),u(n,n),v(n,n)
c 
c        This subroutine reduces the user-supplied lower
c        triangular matrix a to the lower bidiagonal form,
c        so that
c 
c                   A=U \circ B \circ V,                            (1)
c 
c        where A in the input matrix, B is the bidiagonalized
c        matrix, and U, V are orthogonal matrices
c 
c                      Input parameters:
c 
c  a - the triangular matrix to be rediced
c  n - the fdimensionality of a
c 
c                      Output parametders:
c 
c  a - the result of applying the bidiagonalization procedure to a
c  u, v - the orthogonal matrices in (1)
c 
c 
c        . . . initialize the matrices u,v
c 
        do 1400 i=1,n
        do 1200 j=1,n
c 
        u(j,i)=0
        v(j,i)=0
 1200 continue
c 
        u(i,i)=1
        v(i,i)=1
 1400 continue
c 
c        starting with the left lower corner, use Jacobi
c        rotations to bidiagonalize a
c 
        do 2000 i=n,3,-1
c 
        do 1800 j=1,i-2
c 
c       find the 2*2-jacobi matrix to be applied to the columns
c 
        b(1)=a(i,j+1)
        b(2)=a(i,j)
        call fvd_qrrotfnd(b,ub)
c 
c       transform the columns
c 
        call fvd_cols_rotate(ub,a,n,j+1,j)
c 
c        apply the inverse of the column transformation
c        to the rows of matrix v
c 
        call fvd_rows_rotate(ub,v,n,j+1,j)
c 
c       find the 2*2-jacobi matrix to be applied to the rows
c 
        b(1)=a(j+1,j+1)
        b(2)=a(j,j+1)
c 
        call fvd_qrrotfnd(b,ub)
        call fvd_rows_rotate(ub,a,n,j+1,j)
c 
c        apply the inverse of the row transformation
c        to the columns of matrix u
c 
        call fvd_cols_rotate(ub,u,n,j+1,j)
c 
 1800 continue
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
        subroutine fvd_cols_rotate(u,a,n,ii,jj)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(n,n),u(2,2)
c 
        do 1200 i=1,n
c 
        d1=u(1,1)*a(i,ii)+u(1,2)*a(i,jj)
        d2=u(2,1)*a(i,ii)+u(2,2)*a(i,jj)
c 
        a(i,ii)=d1
        a(i,jj)=d2
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine fvd_rows_rotate(u,a,n,ii,jj)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(n,n),u(2,2)
c 
        do 1200 i=1,n
c 
        d1=u(1,1)*a(ii,i)+u(1,2)*a(jj,i)
        d2=u(2,1)*a(ii,i)+u(2,2)*a(jj,i)
c 
        a(ii,i)=d1
        a(jj,i)=d2
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine fvd_qrrotfnd(a,u)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(2),u(2,2)
c 
        u22=-a(1)
        u12=a(2)
c 
        d=sqrt(u22**2+u12**2)
c 
        if(d .eq. 0) then
c 
            u(2,2)=1
            u(1,2)=0
            u(1,1)=1
            u(2,1)=0
            return
        endif
c 
        u(2,2)=u22/d
        u(1,2)=u12/d
        u(1,1)=-u(2,2)
        u(2,1)=u(1,2)
  
        do 1400 i=1,2
        do 1200 j=1,2
        u(j,i)=-u(j,i)
 1200 continue
 1400 continue
c 
        u(2,2)=u22/d
        u(2,1)=u12/d
c 
        u(1,2)=-u(2,1)
        u(1,1)=u(2,2)
c 
        do 2400 i=1,2
        do 2200 j=1,2
        u(j,i)=-u(j,i)
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
        subroutine fvdfinm(u,u2,w,v2,n,m,ncols,ww)
        implicit real *8 (a-h,o-z)
        save
        dimension u(n,ncols),w(m,ncols),u2(ncols,ncols),
     1      v2(ncols,ncols),ww(1)
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
        do 2600 i=1,m
        do 2400 j=1,ncols
c 
        d=0
        do 2200 k=1,ncols
        d=d+w(i,k)*v2(j,k)
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
        subroutine fvdpiv0(a,n,m,b,v,ncols,rnorms,eps,ifpivot)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),b(n,m),v(m,1),rnorms(1)
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
c        real *8 elements long.
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
         call fvdgrm(b,n,m,rnorms,eps,ncols,ifpivot)
c 
c        project the original matrix on the obtained orthogonal
c        columns
c 
        do 2400 i=1,ncols
        do 2200 j=1,m
c 
        call fvdscapr(a(1,j),b(1,i),n,prod)
        v(j,i)=prod
 2200 continue
 2400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine fvdgrm(b,n,m,rnorms,eps,ncols,ifpivot)
        implicit real *8 (a-h,o-z)
        save
        dimension b(n,m),rnorms(1)
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
c        real *8 elements long.
c 
c       . . . prepare the array of values of norms
c             of columns
c 
        done=1
        dtot=0
        do 1400 i=1,m
        d=0
        do 1200 j=1,n
        d=d+b(j,i)**2
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
        d=b(j,i)
        b(j,i)=b(j,ipivot)
        b(j,ipivot)=d
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
        call fvdscapr(b(1,i),b(1,j),n,d)
c 
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*d
 2770 continue
 2780 continue
 2790 continue
c 
c       normalize the i-th column
c 
        call fvdscapr(b(1,i),b(1,i),n,d)
c 
c       if the norm of the remaining part of the matrix
c       is sufficiently small - terminate the process
c 
        if( (ifpivot .eq. 1) .and.
     1      (d .lt. thresh) ) return
        ncols=i
c 
        d=done/dsqrt(d)
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
        call fvdscapr(b(1,i),b(1,j),n,d)
c 
         rrn=0
        do 3000 l=1,n
        b(l,j)=b(l,j)-b(l,i)*d
        rrn=rrn+b(l,j)**2
 3000 continue
        rnorms(j)=dsqrt(rrn)
 3200 continue
 3400 continue
 4000 continue
c 
 4200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine fvdscapr(x,y,n,prod)
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
