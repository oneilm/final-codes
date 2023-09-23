        implicit real *8 (a-h,o-z)
        complex *16 a(1000 000),u(1000 000),
     1      v(1000 000),a2(1000 000),w(2000 000),
     2      a3(1000 000)
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
        lw=2000 000
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
        save
        complex *16 a(n,m),u(n,m),v(m,n),w(1),s(20 000),
     1      a2(n,m),w1(1000 000),w2(1000 000),a3(n,m),ima
        real *8 tt(20 000),yrand(100000)
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
        call bidiag_rand(n,yrand)
  
  
  
        do 1400 i=1,m
        do 1200 j=1,n
        a(j,i)=i**2-j**2
        a(j,i)=dexp(-tt(i)*tt(j))+dexp(-tt(i)/tt(j))
        a(j,i)=a(j,i)+a(j,i)**2*ima
 1200 continue
c 
        a(i,i)=100000 *(1+yrand(i))
        a(i,i)=100000 *yrand(i)
        a(i,i)=100000
cccc        a(i,i)=i**3
  
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
         call prin2('a as created*',a,nprint)
c 
c        perform SVD of a
c 
         lw=2000 000
         call csvdpiv(ier,a,n,m,u,v,s,ncols,eps,
     1      w,lw,ltot)
          call prinf('after csvdpiv, ltot=*',ltot,1)
          call prinf('after csvdpiv, ncols=*',ncols,1)
c 
          call prinf('after csvdpiv, ncols=*',ncols,1)
c 
         nprint=100
         if(n*ncols*2 .lt. nprint) nprint=n*ncols*2
          call prin2('after csvdpiv, u=*',u,nprint)
          call prin2('after csvdpiv, s=*',s,ncols*2)
cccc          call prin2('after csvdpiv, v=*',v,ncols*m*2)
c 
c        finally, multiply them things together
c 
         call mulback(u,w1,v,n,ncols,m,a3,w2,s)
  
cccc         call prin2('after mulback, a3=*',a3,n*m*2)
cccc         call prin2('while a2=*',a2,n*m*2)
c 
c        now, check the error
c 
        dmax=0
        do 1860 i=1,n
        do 1840 j=1,m
c 
        dd=cdabs(a3(i,j)-a2(i,j))
        if(dd .gt. dmax) dmax=dd
 1840 continue
 1860 continue
         call prin2('and maximum error=*',dmax,1)
c 
c       test the orthogonality of columns of u
c 
cccc          call prin2('before orttes, u=*',u, n*n)
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
        implicit real *8 (a-h,o-z)
        save
        complex *16 u(n,ncols),t(ncols,ncols),w(m,ncols),
     1      work(n,ncols),c(n,m),s(1),d
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
        d=d+work(i,k)*dconjg(w(j,k))
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
        implicit real *8 (a-h,o-z)
        save
        complex *16 u(n,ncols),w(ncols,ncols)
c 
c       test the orthogonality of columns of u
c 
        do 1880 i=1,ncols
        do 1870 j=1,ncols
        call csvdscapr(u(1,i),u(1,j),n,w(i,j))
 1870 continue
 1880 continue
cccc         call prin2('matrix of inner products*',w,ncols**2*2)
c 
c       now, check how non-orthogonal it is
c 
        done=1
        do 2400 i=1,ncols
        do 2200 j=1,ncols
        if( (i .ne. j) .and. (dd .lt. cdabs(w(i,j))) )
     1   dd=cdabs(w(i,j) )
c 
        if( (i .eq. j) .and. (dd .lt. cdabs(w(i,j)-done)) )
     1    dd=cdabs(w(i,j)-done)
 2200 continue
 2400 continue
        call prin2('and non-orthogonality of u is*',dd,1)
        return
        end
c 
c 
c 
c 
c 
        subroutine matmua(a,b,n,m,ncols,c)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,ncols),b(m,ncols),c(n,m),cd
c 
        do 1600 i=1,n
        do 1400 j=1,m
c 
        cd=0
        do 1200 k=1,ncols
        cd=cd+a(i,k)*dconjg(b(j,k))
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
        subroutine bidiag_rand(n,y)
        implicit real *8 (a-h,o-z)
        save
        dimension y(1)
        data ifcalled/0/
c 
c       generate pseudo-random numbers
c 
        if(ifcalled .ne. 0) goto 1100
c 
        done=1
        pi=atan(done)*4
        add=sqrt(7*done)
        x=1.1010101010101
        ifcalled=1
 1100 continue
c 
        do 1200 i=1,n
c 
        phi=x*100000*pi+add
        j=phi
        phi=(phi-j)
        x=phi
c 
        y(i)=x
 1200 continue
c 
        return
c 
c 
c 
c 
        entry bidiag_rand_reset(x7)
c 
        done=1
        pi=atan(done)*4
        add=sqrt(7*done)
        x=x7
        ifcalled=1
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
        subroutine csvdpiv(ier,a,n,m,u,v,s,ncols,eps,
     1      w,lw,ltot)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,m),u(n,1),v(m,1),s(1),w(1)
c 
c        this subroutine constructs the singular value decomposition
c        of the usert-specified matrix a, under the assumption
c        that the rank of a is considerably smaller than its
c        dimensionalities n,m. the output of this subroutine
c        consists of two matrices u,v, and a vector s such that
c 
c                 a = u t v^*                                   (1)
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
c        COMPLEX *16 array; imaginary parts of its elements are
c        (hopefully) equal to zero.
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices b, v on exit
c  ltot - the total amount of space (in terms of complex *16 words)
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
        call csvdpiv0(a,n,m,u,v,ncols,w,eps,ifpivot)
c 
c       allocate memory for the subsequent operations
c 
        iw=1
        liw=m*ncols+10
c 
        iu2=iw+liw
        lu2=ncols**2+100
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
        call csvdpiv1(a,u,v,w(iw),w(it),n,m,ncols,
     1    w(irnorms),eps,w(iv2),w(iu2),w(iww),iffirst)
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
        subroutine csvdpiv1(a,u,v,w,t,n,m,ncols,rnorms,
     1    eps2,v2,u2,ww,iffirst)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,m),u(n,1),v(m,1),w(m,1),t(1),
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
        call csvdpiv0(a,n,m,u,v,ncols,rnorms,eps2,ifpivot)
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
        call csvdpiv0(v,m,ncols,w,t,ncols2,rnorms,eps2,ifpivot)
c 
c       reduce matrix t to a real bidiagonal form
c 
        call csvd_rbidiag(t,ncols,u2,v2)
c 
        call csvdfinm(u,u2,w,v2,n,m,ncols,ww)
c 
c        use jacobi rotations to sv-decompose the matrix t
c 
        iiaa=ncols**2/2+10
        call csvd_bidiag_real(t,ncols,u2(iiaa),u2,v2,eps2)
        call csvdrfinm(u,u2,w,v2,n,m,ncols,ww)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine csvd_bidiag_real(t,ncols,a,u3,v3,eps)
        implicit real *8 (a-h,o-z)
        save
        dimension u3(ncols,1),v3(ncols,1),a(ncols,ncols)
c 
        complex *16 t(ncols,ncols)
c 
c       construct the real bidiagonal matrix
c 
        do 1600 i=1,ncols
        do 1400 j=1,ncols
c 
        a(j,i)=0
 1400 continue
c 
        a(i,i)=t(i,i)
        if(i .lt. ncols) a(i+1,i)=t(i+1,i)
 1600 continue
c 
c        diagonalize the real matrix
c 
        ifuv=1
        call csvdrd2mudv(a,ncols,u3,v3,eps,ifuv,numrows)
c 
        call prinf('after csvdrd2mudv, numrows=*',numrows,1)
c 
        do 1800 i=1,ncols
c 
        t(i,1)=a(i,1)
 1800 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine csvdrfinm(u,u2,w,v2,n,m,ncols,ww)
        implicit real *8 (a-h,o-z)
        save
        complex *16 u(n,ncols),w(m,ncols),
     1      ww(1),d
c 
        real *8 u2(ncols,ncols),v2(ncols,ncols)
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
        subroutine csvdfinm(u,u2,w,v2,n,m,ncols,ww)
        implicit real *8 (a-h,o-z)
        save
        complex *16 u(n,ncols),w(m,ncols),u2(ncols,ncols),
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
        d=d+w(i,k)*dconjg(v2(j,k))
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
        subroutine csvdpiv0(a,n,m,b,v,ncols,rnorms,eps,ifpivot)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,m),b(n,m),v(m,1),prod
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
         call csvdgrm(b,n,m,rnorms,eps,ncols,ifpivot)
c 
c        project the original matrix on the obtained orthogonal
c        columns
c 
        do 2400 i=1,ncols
        do 2200 j=1,m
c 
        call csvdscapr(a(1,j),b(1,i),n,prod)
        v(j,i)=dconjg(prod)
 2200 continue
 2400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine csvdgrm(b,n,m,rnorms,eps,ncols,ifpivot)
        implicit real *8 (a-h,o-z)
        save
        complex *16 b(n,m),cd
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
        d=d+b(j,i)*dconjg(b(j,i))
 1200 continue
        rnorms(i)=sqrt(d)
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
        call csvdscapr(b(1,i),b(1,j),n,cd)
c 
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*cd
 2770 continue
 2780 continue
 2790 continue
c 
c       normalize the i-th column
c 
        call csvdscapr(b(1,i),b(1,i),n,cd)
c 
c       if the norm of the remaining part of the matrix
c       is sufficiently small - terminate the process
c 
        d=cd
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
        call csvdscapr(b(1,i),b(1,j),n,cd)
        cd=dconjg(cd)
c 
         rrn=0
        do 3000 l=1,n
        b(l,j)=b(l,j)-b(l,i)*cd
        rrn=rrn+b(l,j)*dconjg(b(l,j))
 3000 continue
        rnorms(j)=dsqrt(rrn)
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
        subroutine csvdscapr(x,y,n,prod)
        implicit complex *16 (a-h,o-z)
        save
        dimension x(1),y(1)
c 
        prod=0
        do 1200 i=1,n
        prod=prod+x(i)*dconjg(y(i))
 1200 continue
        return
        end
  
c 
c 
c 
c 
c 
        subroutine csvdrd2mudv(a,n,u,v,eps,ifuv,numrows)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),u0(2,2),v0(2,2),u(n,n),v(n,n),
     1      b(2,2),rlam(2)
c 
c        this subroutine uses the classical jacobi iteration
c        to construct the singular value decomposition of a
c        general real matrix, of the form
c 
c          a = u s v,                                        (1)
c 
c        with u,v orthogonal matrices, and s diagonal, with
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
c  a - the first column of a contains the (ordered) singular values
c  u,v - orthogonal matrices in (1)
c  numrows - the number of jacobi iterations (i.e. single rotations)
c        that the algorithm has taken
c 
c       . . . one sweep after another, conduct jacobi
c             iterations
c 
        maxsweep=10
        if(ifuv .eq. 0) goto 1150
        zero=0
        do 1100 i=1,n
        do 1050 j=1,n
        u(j,i)=zero
        v(j,i)=zero
 1050 continue
        u(i,i)=1
        v(i,i)=1
 1100 continue
 1150 continue
        numrows=0
        thresh=dsqrt(eps)
        thresh=dsqrt(thresh)
        thresh=dsqrt(thresh)
        do 3400 iijjkk=1,4
c 
        do 3000 ijk=1,maxsweep
        ifany=0
c 
        do 2000 i=1,n
        do 1800 j=1,n
         if(j .eq. i) goto 1800
         if( (dabs(a(i,j)) .lt. thresh) .and.
     1       (dabs(a(j,i)) .lt. thresh) ) goto 1800
        ifany=1
        numrows=numrows+1
c 
c       construct the svd of the (i,j)-th 2 * 2 submatrix
c 
        b(1,1)=a(i,i)
        b(1,2)=a(i,j)
        b(2,1)=a(j,i)
        b(2,2)=a(j,j)
        call csvdroneuv(b,rlam,u0,v0)
c 
c       transform the big matrix
c 
c       . . . the rows
c 
        do 1200 l=1,n
        di=u0(1,1)*a(i,l)-u0(1,2)*a(j,l)
        dj=-u0(2,1)*a(i,l)+u0(2,2)*a(j,l)
        a(i,l)=di
        a(j,l)=dj
 1200 continue
c 
c       . . . and the columns
c 
        do 1400 l=1,n
        di=v0(1,1)*a(l,i)+v0(1,2)*a(l,j)
        dj=v0(2,1)*a(l,i)+v0(2,2)*a(l,j)
        a(l,i)=di
        a(l,j)=dj
 1400 continue
c 
c 
c       if the user so requested, adjust the matrices u,v
c 
        if(ifuv .eq. 0) goto 1800
c 
        do 1600 l=1,n
c 
c       the matrix u
c 
        di=u0(1,1)*u(l,i)-u0(1,2)*u(l,j)
        dj=-u0(2,1)*u(l,i)+u0(2,2)*u(l,j)
        u(l,i)=di
        u(l,j)=dj
c 
c       . . . and the columns
c 
        di=v0(1,1)*v(i,l)+v0(1,2)*v(j,l)
        dj=v0(2,1)*v(i,l)+v0(2,2)*v(j,l)
        v(i,l)=di
        v(j,l)=dj
 1600 continue
c 
 1800 continue
 2000 continue
c 
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
        call csvdrd2mrt(a(1,1),n,u,v,a(1,2),ifuv)
        return
        end
c 
c 
c 
c 
c 
        subroutine csvdrd2mrt(s,n,u,v,iw,ifuv)
        implicit real *8 (a-h,o-z)
        save
        dimension u(n,n),v(n,n),s(1),iw(1)
c 
c        if some of the alleged singular values
c        are negative - change their sign, and that
c        of the corresponding left eigenvectors
c 
        do 1100 i=1,n
        if(s(i) .ge. 0) goto 1100
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
        if(s(j). ge. s(j+1)) goto 1600
        jj=iw(j)
        iw(j)=iw(j+1)
        iw(j+1)=jj
c 
        d=s(j)
        s(j)=s(j+1)
        s(j+1)=d
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
        d=v(i,j)
        v(i,j)=v(ii,j)
        v(ii,j)=d
c 
c       exchange columns number i and ii in the matrix v
c 
        d=u(j,i)
        u(j,i)=u(j,ii)
        u(j,ii)=d
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
        subroutine csvdroneuv(a,rlam,u,v)
        implicit real *8 (a-h,o-z)
        save
        dimension a(2,2),u(2,2),v(2,2),rlam(2),
     1      b(2,2),w(2,2),vstar(2,2),z(2,2),rlams(2,2)
        data eps/1.0d-10/,eps2/1.0d-1/,done/1.0d0/
c 
c       this subroutine produces the singular value decomposition
c       of a 2 * 2 real matrix a, so that on exit,
c 
c       a = u d v,
c 
c       with u and v orthogonal, and d a diagonal matrix with the
c       vector rlam=(rlam(1),rlam(2)) on the diagonal.
c 
c       . . . simmetrize a
c 
        den=a(2,2)+a(1,1)
        rnum=a(2,1)-a(1,2)
c 
c       if denominator is not too small, use exact formulae
c 
        if(dabs(den) .lt. eps*dabs(rnum)) goto 1200
        tgphi=rnum/den
c 
        if(dabs(tgphi) .lt. eps2) goto 1100
c 
        dd=done/(done+tgphi**2)
        alpha=dsqrt(dd)
        beta=dsqrt(done-dd)
        if(tgphi .lt. 0) beta=-beta
        goto 1400
 1100 continue
c 
        dd=tgphi**2/(1+tgphi**2)
c 
        beta=tgphi-tgphi**3/2
        alpha=dsqrt(done-beta**2)
c 
cccc        beta=dsqrt(dd)
cccc        alpha=dsqrt(done-dd)
cccc        if(tgphi .lt. 0) beta=-beta
c 
cccc        phi=datan(tgphi)
cccc        alpha=dcos(phi)
cccc        beta=dsin(phi)
        goto 1400
 1200 continue
c 
c       denominator is too small. use taylor series
c 
        alpha=den/rnum
        beta=1
 1400 continue
c 
c      calculate the simmetrizing matrix and the simmetric one
c 
        w(1,1)=alpha
        w(1,2)=beta
        w(2,1)=-beta
        w(2,2)=alpha
        call csvdrprod2(w,a,b)
c 
c       now, diagonalize the simmetrized matrix
c 
        den=b(2,2)-b(1,1)
        rnum=b(1,2)*2
c 
c       if denominator is not too small, use exact formulae
c 
        if(dabs(den) .lt. eps*dabs(rnum)) goto 1600
        tg2phi=-rnum/den
        if(dabs(tg2phi) .lt. 1.0d-3) goto 1500
c 
        dd=dsqrt(4/tg2phi**2+4)
        if(tg2phi .lt. 0) tgphi=(-2/tg2phi-dd)/2
        if(tg2phi .gt. 0) tgphi=(-2/tg2phi+dd)/2
c 
        phi=datan(tgphi)
        alpha=dcos(phi)
        beta=dsin(phi)
        goto 1800
 1500 continue
c 
        beta=tg2phi/2
        alpha=dsqrt(done-beta**2)
  
cccc        phi=datan(tg2phi)/2
cccc        alpha=dcos(phi)
cccc        beta=dsin(phi)
        goto 1800
 1600 continue
c 
c       denominator is too small. use taylor expansion
c 
        alpha=dsqrt((1-den/rnum)/2)
        beta=dsqrt((1+den/rnum)/2)
 1800 continue
c 
c       construct the diagonalizing matrix
c       and the resulting diagonal
c 
        v(1,1)=alpha
        v(1,2)=beta
        v(2,1)=-beta
        v(2,2)=alpha
c 
        vstar(1,1)=v(1,1)
        vstar(1,2)=v(2,1)
        vstar(2,1)=v(1,2)
        vstar(2,2)=v(2,2)
c 
        call csvdrprod2(v,w,u)
c 
c       finally, compute the resulting diagonal elements
c 
        call csvdrprod2(u,a,z)
        call csvdrprod2(z,vstar,rlams)
        rlam(1)=rlams(1,1)
        rlam(2)=rlams(2,2)
c 
        u(1,2)=-u(1,2)
        u(2,1)=-u(2,1)
        return
        end
c 
c 
c 
c 
c 
        subroutine csvdrprod2(a,b,c)
        implicit real *8 (a-h,o-z)
        save
        dimension a(2,2),b(2,2),c(2,2)
c 
        c(1,1)=a(1,1)*b(1,1)+a(1,2)*b(2,1)
        c(1,2)=a(1,1)*b(1,2)+a(1,2)*b(2,2)
        c(2,1)=a(2,1)*b(1,1)+a(2,2)*b(2,1)
        c(2,2)=a(2,1)*b(1,2)+a(2,2)*b(2,2)
        return
        end
  
c 
c 
c 
c 
c 
        subroutine csvd_rbidiag(a,n,u,v)
        implicit complex *16 (a-h,o-z)
        save
        complex *16 a(n,n),b(3),ub(2,2),u(n,n),v(n,n)
c 
c        This subroutine reduces the user-supplied COMPLEX
c        lower triangular matrix a to the REAL lower bidiagonal
c        form, so that
c 
c                   A=U \circ B \circ V,                            (1)
c 
c        where A in the input matrix, B is the bidiagonalized
c        matrix, and U, V are orthogonal matrices
c 
c                      Input parameters:
c 
c  a - the triangular matrix to be reduced
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
  
        call csvd_rotfnd(b,ub)
c 
c       transform the columns
c 
        call csvd_cols_rotate(ub,a,n,j+1,j)
c 
c        apply the inverse of the column transformation
c        to the rows of matrix v
c 
        call csvd_arrconj(ub,4)
        call csvd_rows_rotate(ub,v,n,j+1,j)
c 
c       find the 2*2-jacobi matrix to be applied to the rows
c 
        b(1)=a(j+1,j+1)
        b(2)=a(j,j+1)
c 
        call csvd_rotfnd(b,ub)
c 
        call csvd_rows_rotate(ub,a,n,j+1,j)
c 
c        apply the inverse of the row transformation
c        to the columns of matrix u
c 
        call csvd_arrconj(ub,4)
        call csvd_cols_rotate(ub,u,n,j+1,j)
c 
 1800 continue
c 
 2000 continue
c 
c       multiply the obtained bidiagonal matrix by diagonal
c       matrices from left and right, converting it into a
c       bidiagonal matrix with real elements
c 
        if(a(1,1) .eq. 0) goto 2200
c 
        cd=a(1,1)/abs(a(1,1))
        a(1,1)=a(1,1)*conjg(cd)
        call csvd_colmult(u,n,1,cd)
c 
 2200 continue
c 
        do 3000 i=2,n
c 
        if(a(i,i-1) .eq. 0) goto 2600
c 
        cd=a(i,i-1)/abs(a(i,i-1))
        a(i,i-1)=a(i,i-1)*conjg(cd)
        a(i,i)=a(i,i)*conjg(cd)
        call csvd_colmult(u,n,i,cd)
c 
 2600 continue
c 
        if(a(i,i) .eq. 0) goto 3000
c 
        cd=a(i,i)/abs(a(i,i))
        a(i,i)=a(i,i)*conjg(cd)
        call csvd_rowmult(v,n,i,cd)
c 
        if(i .eq. n) goto 3000
c 
        a(i+1,i)=a(i+1,i)*conjg(cd)
c 
 3000 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine csvd_rotfnd(a,u)
        implicit complex *16 (a-h,o-z)
        save
        dimension a(2),u(2,2)
c 
        u21=-a(2)
        u22=a(1)
c 
        d=sqrt(u22*conjg(u22)+u21*conjg(u21))
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
        u(2,1)=u21/d
c 
        u(1,1)=-conjg(u(2,2))
        u(1,2)=conjg(u(2,1))
        return
        end
c 
c 
c 
c 
c 
        subroutine csvd_colmult(a,n,ii,coef)
        implicit complex *16 (a-h,o-z)
        save
        dimension a(n,n)
c 
        do 1200 i=1,n
c 
        a(i,ii)=a(i,ii)*coef
 1200 continue
c 
        return
c 
c 
c 
c 
        entry csvd_rowmult(a,n,ii,coef)
c 
        do 2200 i=1,n
c 
        a(ii,i)=a(ii,i)*coef
 2200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine csvd_cols_rotate(u,a,n,ii,jj)
        implicit complex *16 (a-h,o-z)
        save
        dimension a(n,n),u(2,2)
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
        subroutine csvd_rows_rotate(u,a,n,ii,jj)
        implicit complex *16 (a-h,o-z)
        save
        dimension a(n,n),u(2,2)
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
        subroutine csvd_arrconj(a,n)
        implicit complex *16 (a-h,o-z)
        save
        dimension a(n)
c 
        do 1200 i=1,n
c 
        a(i)=conjg(a(i))
 1200 continue
c 
        return
        end
  
