        implicit real *8 (a-h,o-z)
        dimension a(50 000),u(50 000),
     1      v(50 000),a2(50 000),w(100 000),
     2      a3(50 000)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
         PRINT *, 'ENTER n'
         READ *,n
         CALL PRINF('n=*',n,1 )
C 
         PRINT *, 'ENTER m'
         READ *,m
         CALL PRINF('m=*',m,1 )
c 
c 
        PRINT *, 'ENTER -log(eps) base 10'
        READ *,db
        CALL PRIN2('db=*',db,1 )
c 
        eps=0.1**db
        call prin2('and eps=*',eps,1)
c 
c       test the gram-schmidt procedure
c 
        lw=300000
        call leastes(a,n,m,u,v,w,lw,eps,a2,a3)
  
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine leastes(a,n,m,u,v,w,lw,eps,a2,c)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,m),v(m,n),w(50 000),t(50 000),
     1      c(n,m),rnorms(1000),tt(1000),a2(n,m),cd,
     2      x(1000),y(1000),work(1000),x2(1000),ima
c 
        data ima/(0.0d0,1.0d0)/
c 
c        construct the test matrix to be compressed
c 
        h=1
        h=h/n
        do 1100 i=1,2000
        tt(i)=i*h+ima/10
cccc        tt(i)=i*h
 1100 continue
c 
        do 1400 i=1,m
        do 1200 j=1,n
        a(j,i)=cdexp(-tt(i)*tt(j))+cdexp(-tt(i)/tt(j))
cccc        a(n-j+1,i)=dexp(-tt(i)*tt(j))+dexp(-tt(i)/tt(j))
 1200 continue
 1400 continue
  
        call prin2('a as created*',a,n*m)
  
c 
c      copy a into a2
c 
        do 1800 i=1,m
        do 1600 j=1,n
        a2(j,i)=a(j,i)
 1600 continue
 1800 continue
c 
c        perform the QR of the matrix a
c 
  
        call cleastsq(a,u,w,t,n,m,ncols,rnorms,eps,v)
  
cccc        call cleast1(a,u,w,t,n,m,ncols,rnorms,
cccc     1    eps,v)
c 
        call prinf('after clieast1, ncols=*',ncols,1)
        call prin2('after clieast1, rnorms=*',rnorms,ncols)
        call prin2('after clieast1, a=*',a,n*m)
        call prin2('after clieast1, u=*',u,ncols*n)
        call prin2('after clieast1, w=*',w,ncols*m)
  
cccc           if(2 .ne. 3) stop
c 
c       test the accuracy
c 
        call cmulback(u,w,t,n,ncols,m,c,v)
c 
          call prin2('and after cmulback, c=*',c,n*m)
  
c 
cccc                if(2 .ne. 3) stop
        cd=0
        do 2200 i=1,m
        do 2000 j=1,n
        cd=cd+(a2(j,i)-c(j,i))**2
 2000 continue
 2200 continue
        call prin2('and the error is*',cdsqrt(cd),1)
c 
c        check the least squares inverter
c 
        do 3200 i=1,m
        x(i)=i
 3200 continue
c 
        call prin2('checking the inverter, x=*',x,m*2)
  
        call cmatve(a2,x,y,n,m)
  
  
        call prin2('checking the inverter, y=*',y,m*2)
  
c 
        call cleasts2(u,w,t,n,m,ncols,y,x2,work)
  
        call prin2('checking the inverter, x2=*',x2,m*2)
  
       return
       end
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c       this is the end of the debugging code and the beginning of the
c       least squares subroutines proper
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine cleastsq(a,u,w,t,n,m,ncols,rnorms,eps,v)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,m),u(n,1),v(m,1),w(m,1),t(1)
        dimension rnorms(1)
c 
c        this subroutine decomposes the input matrix a in the
c        form
c               A=U * T * W^*                                       (1)
c 
c        with the columns of each of the matrices U,W orthonormal,
c        and T a lower triangular matrix of dimensionality (ncols,ncols)
c 
c   NOTE: this subroutine uses the subroutine leastsq1 (see) to perform
c        almost all work. However, when m < n, it tends to be much
c        more efficient than leastsq1, since it performs the two
c        Gram-Schmidt processes in the optimal order (starting with
c        the longer gram-schmidt involving shorter vectors)
c 
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
c  u - the matrix in (1). note that the matrix u has to be dimensioned
c        u(n,m), but on exit from this subroutine only the first
c        ncols columns of  u  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit u is dimensioned
c        u(n,ncols)
c  w - the third matrix in (1). it is dimensioned w(m,ncols) on
c        exit. on entry the parameter ncols is unknown,
c        so w should be dimensioned w(m,n) by the user. however,
c        on exit only the first ncols columns of w are meaningful
c  t - the second matrix in (1); On exit, it is structured as an
c        ncols * ncols matrix; however, on entry ncols is not
c        known, so it should be at least m*n complex *16 locations long.
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices u, w on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c 
c 
c                      work arrays:
c 
c  v - must be at least n*m real *8 locations long
c 
c       . . . if m > n , construct the decomposition of the matrix as
c             specified by the user
c 
        if(m .lt. n-2) goto 2000
c 
        call cleast1(a,u,w,t,n,m,ncols,rnorms,eps,v)
cccc        call prin2('after first leastsq1, u=*',u,n*ncols)
c 
        return
c 
 2000 continue
c 
c       n is greater than m. transpose the matrix, and decompose
c       the transpose
c 
        call cleastra(a,n,m,v)
        call cleascop(v,a,n*m)
c 
        call cleast1(a,w,u,t,m,n,ncols,rnorms,eps,v)
c 
cccc        call prin2('after leastsq1, u=*',u,n*ncols)
c 
c        transpose back everything that needs to be transposed back
c 
        call cleastra(a,m,n,v)
        call cleascop(v,a,n*m)
c 
        call cleastra(t,ncols,ncols,v)
        call cleascop(v,t,ncols**2)
c 
        call cleasrer(t,ncols,ncols,v)
        call cleascop(v,t,ncols**2)
c 
        call cleasrec(u,n,ncols,v)
        call cleascop(v,u,n*ncols)
c 
c 
c 
        call cleasrec(w,m,ncols,v)
        call cleascop(v,w,m*ncols)
c 
        call cleasrec(t,ncols,ncols,v)
        call cleascop(v,t,ncols**2)
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine cleasts2(u,w,t,n,m,ncols,y,x,work,whts)
        implicit real *8 (a-h,o-z)
        save
        complex *16 u(n,ncols),w(m,ncols),t(ncols,ncols)
        dimension whts(1)
        complex *16 x(m),y(n),work(ncols),d
c 
c         this subroutine uses the QR-type (not exactly QR) decomposition
c         of a matrix to solve in the least squares sense the linear
c         system
c 
c                    A X = Y                                             (1)
c 
c         The expansion used must have been produced by a prior call
c         to the subroutine leastsq1 (see), and is of the form
c 
c               A=U * T * W^*                                            (2)
c 
c        with the columns of each of the matrices U,W orthonormal,
c        and T a lower triangular matrix of dimensionality (ncols,ncols)
c 
c                     input parameters:
c 
c  u - the matrix in (1). note that the matrix u has to be dimensioned
c        u(n,m), but on exit from this subroutine only the first
c        ncols columns of  u  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit u is dimensioned
c        u(n,ncols)
c  w - the third matrix in (1). it is dimensioned w(m,ncols) on
c        exit. on entry the parameter ncols is unknown,
c        so w should be dimensioned w(m,n) by the user. however,
c        on exit only the first ncols columns of w are meaningful
c  t - the second matrix in (1)
c  n,m - the dimensionalities of a
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices u, w on exit
c  y - the right-hand side in (1)
c 
c                     output parameters:
c 
c  x - the solution of the system (1) in the least squares sense
c 
c                      work arrays:
c 
c  work - must be at least ncols real *8 locations long
c 
c        . . . apply to the right-hand side the matrux  U^*
c 
        do 1400 i=1,ncols
        d=0
        do 1200 j=1,n
cccc        d=d+u(j,i)*y(j)*whts(j)
cccccc        d=d+u(j,i)*y(j)
        d=d+dconjg(u(j,i))*y(j)
 1200 continue
        work(i)=d
 1400 continue
c 
c       apply to the vector work the inverse of the matrix t
c 
        x(1)=work(1)/t(1,1)
        do 2000 i=2,ncols
c 
        d=0
        do 1600 j=1,i-1
        d=d+t(i,j)*x(j)
 1600 continue
c 
        x(i)=(work(i)-d)/t(i,i)
 2000 continue
c 
        do 2200 i=1,ncols
        work(i)=x(i)
 2200 continue
c 
c       multiply work by the matrix w
c 
        do 2600 i=1,m
        d=0
        do 2400 j=1,ncols
        d=d+w(i,j)*work(j)
 2400 continue
        x(i)=d
 2600 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine cleastra(a,n,m,b)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,m),b(m,n),x(1),y(1)
c 
c       transpose a
c 
        do 1400 i=1,n
        do 1200 j=1,m
cccc        b(j,i)=a(i,j)
        b(j,i)=dconjg(a(i,j))
 1200 continue
 1400 continue
        return
c 
c 
c 
c 
        entry cleascop(x,y,n)
c 
        do 2400 i=1,n
        y(i)=x(i)
 2400 continue
        end
c 
c 
c 
c 
c 
        subroutine cleast1(a,u,w,t,n,m,ncols,rnorms,
     1    eps,v)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,m),u(n,1),v(m,1),w(m,1),t(1),rnorms(1)
c 
c        this subroutine decomposes the input matrix a in the
c        form
c               A=U * T * W^*                                       (1)
c 
c        with the columns of each of the matrices U,W orthonormal,
c        and T a lower triangular matrix of dimensionality (ncols,ncols)
c 
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
c  u - the matrix in (1). note that the matrix u has to be dimensioned
c        u(n,m), but on exit from this subroutine only the first
c        ncols columns of  u  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit u is dimensioned
c        u(n,ncols)
c  w - the third matrix in (1). it is dimensioned w(m,ncols) on
c        exit. on entry the parameter ncols is unknown,
c        so w should be dimensioned w(m,n) by the user. however,
c        on exit only the first ncols columns of w are meaningful
c  t - the second matrix in (1)
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices u, w on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c 
c 
c                      work arrays:
c 
c  v - must be at least n*m real *8 locations long
c 
c        . . . using gram-schmidt process with pivoting, decompose
c              the matrix a in the form
c 
c          a=U  V^*,
c 
c        with  u an orthogonal matrix of minimum rank,
c        and v whatever it wants to be
c 
        ifpivot=1
        call cleaspiv(a,n,m,u,v,ncols,rnorms,eps,ifpivot)
         call prin2('after first leaspiv0, rnorms=*',rnorms,ncols)
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
        call cleaspiv(v,m,ncols,w,t,ncols2,rnorms,eps,ifpivot)
         call prin2('after second leaspiv0, rnorms=*',rnorms,ncols2)
  
  
        if(2 .ne. 3) return
c 
c       if the need be - restructure the matrix t
c 
        if(ncols2 .eq. ncols) return
c 
cccc        call leastres(t,v,t,ncols,ncols2)
        ncols=ncols2
  
        return
        end
c 
c 
c 
c 
c 
c 
        subroutine cmatve(a,x,y,n,m)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,m),x(n),y(m),d
c 
c       apply the matrix a to the vector x getting y
c 
        do 1400 i=1,n
        d=0
        do 1200 j=1,m
        d=d+a(i,j)*x(j)
 1200 continue
        y(i)=d
 1400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine cmulback(u,w,t,n,ncols,m,c,work)
        implicit real *8 (a-h,o-z)
        save
        complex *16 w(m,ncols),u(n,ncols),c(n,m),t(ncols,ncols),
     1      work(n,ncols),d
c 
c       multiply u by t
c 
        do 2600 i=1,n
        do 2400 j=1,ncols
        d=0
        do 2200 k=1,ncols
        d=d+u(i,k)*t(k,j)
 2200 continue
        work(i,j)=d
 2400 continue
 2600 continue
c 
c        multiply work by w ^*
c 
        do 3600 i=1,n
        do 3400 j=1,m
c 
        d=0
        do 3200 k=1,ncols
        d=d+work(i,k)*dconjg(w(j,k))
 3200 continue
        c(i,j)=d
 3400 continue
 3600 continue
c 
        return
        end
c 
c 
c 
c 
c 
c 
        subroutine cleaspiv(a,n,m,b,v,ncols,rnorms,eps,ifpivot)
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
         call cleasgrm(b,n,m,rnorms,eps,ncols,ifpivot)
c 
c        project the original matrix on the obtained orthogonal
c        columns
c 
        do 2400 i=1,ncols
        do 2200 j=1,m
c 
        call cleascap(a(1,j),b(1,i),n,prod)
cccc        v(j,i)=prod
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
        subroutine cleasgrm(b,n,m,rnorms,eps,ncols,ifpivot)
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
cccc        d=d+b(j,i)**2
        d=d+b(j,i)*dconjg(b(j,i))
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
        call cleascap(b(1,i),b(1,j),n,cd)
c 
cccc        cd=dconjg(cd)
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*cd
 2770 continue
 2780 continue
 2790 continue
c 
c       normalize the i-th column
c 
        call cleascap(b(1,i),b(1,i),n,cd)
c 
c       if the norm of the remaining part of the matrix
c       is sufficiently small - terminate the process
c 
        d=cd
        if(d .lt. thresh) return
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
        call cleascap(b(1,i),b(1,j),n,cd)
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
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cleascap(x,y,n,prod)
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
        subroutine cleasrec(a,n,m,b)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,m),b(n,m)
c 
c       reverse the columns of a
c 
        do 1400 i=1,m
        do 1200 j=1,n
        b(j,i)=a(j,m-i+1)
 1200 continue
 1400 continue
        return
c 
c 
c 
c 
        entry cleasrer(a,n,m,b)
c 
c       reverse the rows of a
c 
        do 2400 i=1,m
        do 2200 j=1,n
        b(j,i)=a(n-j+1,i)
 2200 continue
 2400 continue
        return
        end
