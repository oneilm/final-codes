  
c 
c 
c 
c 
c 
      SUBROUTINE clegeexe(X,VAL,PEXP,N)
      IMPLICIT REAL *8 (A-H,O-Z)
        save
      complex *16 PEXP(1),val
C 
C     This subroutine computes the value o a Legendre
c     expansion with coefficients PEXP at point X in interval [-1,1]
C 
c                input parameters:
c 
C     X = evaluation point
C     PEXP = expansion coefficients
C     N  = order of expansion
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c 
c                output parameters:
c 
C     VAL = computed value
C 
        done=1
        pjm2=1
        pjm1=x
c 
        val=pexp(1)*pjm2+pexp(2)*pjm1
c 
        DO 600 J = 2,N
c 
        pj= ( (2*j-1)*x*pjm1-(j-1)*pjm2 ) / j
        val=val+pexp(j+1)*pj
c 
        pjm2=pjm1
        pjm1=pj
 600   CONTINUE
c 
        RETURN
        END
c 
c 
c 
c 
c 
        subroutine scawrong(fldcoe,g,nphi,ntheta,fld)
        implicit real *8 (a-h,o-z)
        save
        complex *16 fldcoe(1),g(1),fld
c 
        fld=0
        do 1200 i=1,nphi*ntheta
        fld=fld+g(i)*fldcoe(i)
 1200 continue
        return
        end
c 
c 
c 
c 
c 
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
        subroutine leastsq(a,u,w,t,n,m,ncols,rnorms,eps,v)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),u(n,1),v(m,1),w(m,1),t(1),rnorms(1)
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
c       . . . if m > n , construct the decomposition of the matrix as
c             specified by the user
c 
        if(m .lt. n-2) goto 2000
cccc        if(n .lt. m) goto 2000
c 
        call leastsq1(a,u,w,t,n,m,ncols,rnorms,eps,v)
cccc        call prin2('after first leastsq1, u=*',u,n*ncols)
c 
        return
c 
 2000 continue
c 
c       n is greater than m. transpose the matrix, and decompose
c       the transpose
c 
        call leastran(a,n,m,v)
        call leascopy(v,a,n*m)
c 
        call leastsq1(a,w,u,t,m,n,ncols,rnorms,eps,v)
c 
cccc        call prin2('after leastsq1, u=*',u,n*ncols)
c 
c        transpose back everything that needs to be transposed back
c 
        call leastran(a,n,m,v)
        call leascopy(v,a,n*m)
c 
        call leastran(t,ncols,ncols,v)
        call leascopy(v,t,ncols**2)
c 
        call leasrevr(t,ncols,ncols,v)
        call leascopy(v,t,ncols**2)
c 
        call leasrevc(u,n,ncols,v)
        call leascopy(v,u,n*ncols)
c 
c 
c 
        call leasrevc(w,m,ncols,v)
        call leascopy(v,w,m*ncols)
c 
        call leasrevc(t,ncols,ncols,v)
        call leascopy(v,t,ncols**2)
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine leastran(a,n,m,b)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),b(m,n),x(1),y(1)
c 
c       transpose a
c 
        do 1400 i=1,n
        do 1200 j=1,m
        b(j,i)=a(i,j)
 1200 continue
 1400 continue
        return
c 
c 
c 
c 
        entry leascopy(x,y,n)
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
        subroutine leastsq2(u,w,t,n,m,ncols,y,x,work)
        implicit real *8 (a-h,o-z)
        save
        dimension u(n,ncols),w(m,ncols),t(ncols,ncols),
     1      x(m),y(n),work(ncols)
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
        d=d+u(j,i)*y(j)
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
        subroutine cleastsq2(u,w,t,n,m,ncols,y,x,work,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension u(n,ncols),w(m,ncols),t(ncols,ncols),whts(1)
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
        d=d+u(j,i)*y(j)*whts(j)
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
        subroutine leastsq1(a,u,w,t,n,m,ncols,rnorms,
     1    eps,v)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),u(n,1),v(m,1),w(m,1),t(1),rnorms(1)
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
        call leaspiv0(a,n,m,u,v,ncols,rnorms,eps,ifpivot)
cccc         call prin2('after first leaspiv0, rnorms=*',rnorms,ncols)
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
        call leaspiv0(v,m,ncols,w,t,ncols2,rnorms,eps,ifpivot)
cccc         call prin2('after second leaspiv0, rnorms=*',rnorms,ncols2)
c 
c       if the need be - restructure the matrix t
c 
        if(ncols2 .eq. ncols) return
c 
        call leastres(t,v,t,ncols,ncols2)
        ncols=ncols2
  
        return
        end
c 
c 
c 
c 
c 
        subroutine leastres(a,b,c,n1,n2)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n1,n1),b(n2,n2),c(n2,n2)
c 
        do 1400 i=1,n2
        do 1200 j=1,n2
        b(j,i)=a(j,i)
 1200 continue
 1400 continue
c 
        do 2400 i=1,n1
        do 2200 j=1,n1
        c(j,i)=b(j,i)
 2200 continue
 2400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine leaspiv0(a,n,m,b,v,ncols,rnorms,eps,ifpivot)
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
         call leasgrm(b,n,m,rnorms,eps,ncols,ifpivot)
c 
c        project the original matrix on the obtained orthogonal
c        columns
c 
        do 2400 i=1,ncols
        do 2200 j=1,m
c 
        call leascapr(a(1,j),b(1,i),n,prod)
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
        subroutine leasgrm(b,n,m,rnorms,eps,ncols,ifpivot)
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
        call leascapr(b(1,i),b(1,j),n,d)
c 
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*d
 2770 continue
 2780 continue
 2790 continue
c 
c       normalize the i-th column
c 
        call leascapr(b(1,i),b(1,i),n,d)
c 
c       if the norm of the remaining part of the matrix
c       is sufficiently small - terminate the process
c 
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
        call leascapr(b(1,i),b(1,j),n,d)
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
  
 4200 continue
  
  
        return
        end
c 
c 
c 
c 
c 
        subroutine leascapr(x,y,n,prod)
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
        subroutine leasrevc(a,n,m,b)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),b(n,m)
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
        entry leasrevr(a,n,m,b)
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
c 
c 
c 
c 
c 
        subroutine gauswhts(n,ts,whts,ifwhts)
        implicit real *8 (a-h,o-z)
        save
        dimension ts(1),whts(1)
c 
c        this subroutine constructs the nodes and the
c        weights of the n-point gaussian quadrature on
c        the interval [-1,1]
c 
c                input parameters:
c 
c  n - the number of nodes in the quadrature
c 
c                output parameters:
c 
c  ts - the nodes of the n-point gaussian quadrature
c  w - the weights of the n-point gaussian quadrature
c 
c       . . . construct the array of initial approximations
c             to the roots of the n-th legendre polynomial
c 
        eps=1.0d-14
        ZERO=0
        DONE=1
        pi=datan(done)*4
        h=pi/(2*n)
        do 1200 i=1,n
        t=(2*i-1)*h
        ts(n-i+1)=dcos(t)
1200  CONTINUE
c 
c         use newton to find all roots of the legendre polynomial
c 
        ts(n/2+1)=0
        do 2000 i=1,n/2
c 
        xk=ts(i)
        deltold=1
        do 1400 k=1,10
        call legpol(xk,n,pol,der)
        delta=-pol/der
cccc         call prin2('delta=*',delta,1)
        xk=xk+delta
        if(deltold .le. eps) goto 1600
        deltold=dabs(delta)
 1400 continue
 1600 continue
        ts(i)=xk
        ts(n-i+1)=-xk
 2000 continue
  
        call prinf('in gauswhts, ifwhts=*',ifwhts,1)
c 
        if(ifwhts .eq. 0) return
c 
c       now, use the explicit integral formulae
c       to obtain the weights
c 
        a=-1
        b=1
        do 2200 i=1,n/2+1
        call prodend(a,ts,n,i,fm)
        call prodend(b,ts,n,i,fp)
        whts(i)=fp-fm
        whts(n-i+1)=whts(i)
 2200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine legpol(x,n,pol,der)
        implicit real *8 (a-h,o-z)
c 
        save
        pkm1=1
        pk=x
c 
        pk=1
        pkp1=x
c 
c        if n=0 or n=1 - exit
c 
        if(n .ge. 2) goto 1200
        pol=1
        der=0
        if(n .eq. 0) return
c 
        pol=x
        der=1
        return
 1200 continue
c 
c       n is greater than 1. conduct recursion
c 
        do 2000 k=1,n-1
        pkm1=pk
        pk=pkp1
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
 2000 continue
c 
c        calculate the derivative
c 
        pol=pkp1
        der=n*(x*pkp1-pk)/(x**2-1)
        return
        end
c 
c 
c 
c 
c 
        subroutine prodend(x,xs,n,i,f)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1)
c 
c      evaluate the product
c 
        f=1
        dlarge=1.0d20
        dsmall=f/dlarge
c 
        large=0
        do 2000 j=1,n
        dd=dabs(f)
        if( dd .gt. dsmall) goto 1200
         f=f*10000
         large=large-1
 1200 continue
c 
        if( dd .lt. dlarge) goto 1400
        f=f/10000
        large=large+1
 1400 continue
        if(j .eq. i) goto 2000
        f=f*(x-xs(j))/(xs(i)-xs(j))
 2000 continue
        d10000=10000
        f=f*d10000**large
        f=f**2*(x-xs(i))
        return
        end
c 
c 
c 
c 
cccc        SUBROUTINE RSplot(Z,N,IW)
        SUBROUTINE RSPLOT(X,Y,N,IW,mes)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        DIMENSION X(1),Y(1),Z(2,1)
        character *1 mes(1)
C 
C        THIS SUBROUTINE PLOTS THE CURVE SPECIFIED BY ITS NODES
C        X,Y (CONCEPTUALLY, THE CURVE CONNECTS THE POINT (X(1),Y(1))
C        WITH THE POINT (X(N),Y(N)). IF THE OUTPUT FILE IS 'INCLUDE' D
C        INTO TD OR TD3 STREAM, A FILE NAMED "PLOT.IMP" IS CREATED.
C        THIS IS A TRUE, UNADULTERATED 'IMPRINT' FILE, WHICH CAN AND
C        SHOULD BE "IMPRINTED".
C     NOTE: NATURALLY, CERTAIN DECISIONS WERE MADE IN WRITING THIS
C        ROUTINE (FRAME, SIZE, AND OTHER PARAMETERS OF THE PLOT.
C        IF YOU DON'T LIKE THEM - YOU ARE FREE TO WRITE YOUR OWN
C        SUBROUTINE!!!
C 
C  INPUT PARAMETERS:
C   X,Y - THE COORDINATES OF THE CURVE TO PLOT
C   N - THE NUMBER OF ELEMENTS IN ARAYS X,Y.
C   IW - THE FORTRAN UNIT NUMBER ON WHICH THE OUTPUT DATA SET IS WRITTEN
C 
C  OUTPUT PARAMETERS: NONE
C 
c 
c       . . . print the label
c 
        IF(IW.EQ.0) RETURN
        call titlpr(MES,Iw)
C 
C 
C . . . FIND THE LIMITS TO BE SET
C 
        XMIN=1.0D20
        YMIN=1.0D20
        XMAX=-1.0D20
        YMAX=-1.0D20
        DO 1000 I=1,N
       IF(XMIN.GT.X(I))XMIN=X(I)
       IF(XMAX.LT.X(I))XMAX=X(I)
       IF(YMIN.GT.Y(I))YMIN=Y(I)
       IF(YMAX.LT.Y(I))YMAX=Y(I)
ccccc        IF(XMIN.GT.Z(1,I))XMIN=Z(1,I)
ccccc        IF(XMAX.LT.Z(1,I))XMAX=Z(1,I)
ccccc        IF(YMIN.GT.Z(2,I))YMIN=Z(2,I)
ccccc        IF(YMAX.LT.Z(2,I))YMAX=Z(2,I)
 1000 CONTINUE
        DX=(XMAX-XMIN)
        DY=YMAX-YMIN
        DD=DX
        IF(DY.GT.DX) DD=DY
        DMIN=XMIN
        IF(YMIN.LT.XMIN) DMIN=YMIN
        DMAX=XMAX
        IF(YMAX.GT.XMAX) DMAX=YMAX
        DMIN=DMIN-DD/10
        DMAX=DMAX+DD/10
 1200 FORMAT(1X,37H#SET DEVICE postscript FILE 'plot.ps',/,
c1200 FORMAT(1X,34H#SET DEVICE IMAGEN FILE 'PLOT.IMP',/,
C1200 FORMAT(1X,26HSET IMAGEN FILE 'PLOT.IMP',/,
ccc     A       1X,'SET DEVICE SUN',/,
C    1  1X,'SET DEV MC',/,
     2  1X,'SET AXIS TOP OFF',/,
     3  1X,'SET AXIS RIGHT OFF',/,
     4  1X,'SET LIMITS X ', 2(1X,E11.5),' Y ',2(1X,E11.5),/,
CCC  5  1X,'SET LIMITS Y ', 2(2X,E11.5),/,
     5  1X,'SET SCALE EQUAL')
          IF(IW.EQ.0) RETURN
          WRITE(IW,1200) DMIN,DMAX,DMIN,DMAX
 1400 FORMAT(2(2X,E15.7))
ccccc       WRITE(IW,1400) (Z(1,I),Z(2,I),I=1,N)
       WRITE(IW,1400) (X(I),Y(I),I=1,N)
cccc 1500 FORMAT(1X,'plot ''*'' ')
 1500 FORMAT(1X,'PLOT')
C1600 FORMAT(1X,'FRAME')
C1800 FORMAT(1X,'EXIT')
 1450 FORMAT(1X,'JOIN')
c1500 FORMAT(1X,'PLOT')
 1600 FORMAT(1X,'#FRAME')
 1800 FORMAT(1X,'#EXIT')
        WRITE(IW,1500)
        WRITE(IW,1600)
        WRITE(IW,1800)
        RETURN
C 
C 
C 
C 
        ENTRY RSGRAF(X,Y,N,IW,mes)
C 
C        THIS SUBROUTINE PLOTS THE Y AS A FUNCTION OF X.
C        IF THE OUTPUT FILE IS 'INCLUDE' D
C        INTO TD OR TD3 STREAM, A FILE NAMED "PLOT.IMP" IS CREATED.
C        THIS IS A TRUE, UNADULTERATED 'IMPRINT' FILE, WHICH CAN AND
C        SHOULD BE "IMPRINTED".
C     NOTE: NATURALLY, CERTAIN DECISIONS WERE MADE IN WRITING THIS
C        ROUTINE (FRAME, SIZE, AND OTHER PARAMETERS OF THE PLOT).
C        IF YOU DON'T LIKE THEM - YOU ARE FREE TO WRITE YOUR OWN
C        SUBROUTINE!!!
C 
C  INPUT PARAMETERS:
C   X,Y - THE COORDINATES OF THE CURVE TO PLOT
C   N - THE NUMBER OF ELEMENTS IN ARAYS X,Y.
C   IW - THE FORTRAN UNIT NUMBER ON WHICH THE OUTPUT DATA SET IS WRITTEN
C 
C  OUTPUT PARAMETERS: NONE
C 
        IF(IW.EQ.0) RETURN
        call titlpr(MES,Iw)
 2200 FORMAT(1X,37H#SET DEVICE postscript FILE 'plot.ps',/,
c2200 FORMAT(1X,26HSET IMAGEN FILE 'PLOT.IMP',/,
ccc     A       1X,'SET DEVICE SUN',/,
ccc  1  1X,'SET DEV MC',/,
     2  1X,'SET AXIS TOP OFF',/,
     3  1X,'SET AXIS RIGHT OFF')
CC   3  1X,'SET AXIS RIGHT OFF',/,
CC   4  1X,'SET SCALE EQUAL')
          WRITE(IW,2200)
          WRITE(IW,1400) (X(I),Y(I),I=1,N)
          WRITE(IW,1450)
          WRITE(IW,1600)
          WRITE(IW,1800)
  
        close(iw)
  
  
          RETURN
          END
c 
c 
c 
c 
c 
        SUBROUTINE titlpr(MES,IP)
        save
        CHARACTER *1 MES(1),AST,card(80),blank,title(5),quo
        DATA AST/'*'/,quo/''''/,blank/' '/
     1      title/'t','i','t','l','e'/
C 
C         DETERMINE THE LENGTH OF THE MESSAGE
C 
        I=0
        DO 1400 I=1,10000
        IF(MES(I).EQ.AST) GOTO 1600
        I1=I
 1400 CONTINUE
 1600 CONTINUE
 1800 FORMAT(1X, 80A1)
         do 2000 i=1,80
         card(i)=blank
 2000 continue
         do 2200 i=1,5
         card(i+1)=title(i)
 2200 continue
         card(8)=quo
c 
         do 2400 i=1,i1
         card(i+8)=mes(i)
 2400 continue
         card(8+i1+1)=quo
c 
         write(ip,1800) card
         return
         end
  
C 
C 
C 
C 
        SUBROUTINE PRINI(IP1,IQ1)
        save
        CHARACTER *1 MES(1), AA(1)
        REAL *4 A(1)
        REAL *8 A2(1)
        REAL *16 A4(1)
        INTEGER *4 IA(1)
        INTEGER *2 IA2(1)
        IP=IP1
        IQ=IQ1
        RETURN
  
C 
C 
C 
C 
C 
        ENTRY PRIN(MES,A,N)
        CALL  MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1200)(A(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1200)(A(J),J=1,N)
 1200 FORMAT(6(2X,E11.5))
         RETURN
C 
C 
C 
C 
        ENTRY PRIN2(MES,A2,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1400)(A2(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1400)(A2(J),J=1,N)
 1400 FORMAT(6(2X,E11.5))
        RETURN
C 
C 
C 
C 
        ENTRY PRINQ(MES,A4,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1500)(A4(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1500)(A4(J),J=1,N)
 1500 FORMAT(6(2X,e11.5))
        RETURN
C 
C 
C 
C 
        ENTRY PRINF(MES,IA,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1600)(IA(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1600)(IA(J),J=1,N)
 1600 FORMAT(10(1X,I7))
        RETURN
C 
C 
C 
C 
        ENTRY PRINF2(MES,IA2,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1600)(IA2(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1600)(IA2(J),J=1,N)
        RETURN
C 
C 
C 
C 
        ENTRY PRINA(MES,AA,N)
        CALL MESSPR(MES,IP,IQ)
 2000 FORMAT(1X,80A1)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,2000)(AA(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,2000)(AA(J),J=1,N)
        RETURN
        END
c 
c 
c 
c 
c 
        SUBROUTINE MESSPR(MES,IP,IQ)
        save
        CHARACTER *1 MES(1),AST
        DATA AST/'*'/
C 
C         DETERMINE THE LENGTH OF THE MESSAGE
C 
        I=0
        DO 1400 I=1,10000
        IF(MES(I).EQ.AST) GOTO 1600
        I1=I
 1400 CONTINUE
 1600 CONTINUE
         IF ( (I1.NE.0) .AND. (IP.NE.0) )
     1     WRITE(IP,1800) (MES(I),I=1,I1)
         IF ( (I1.NE.0) .AND. (IQ.NE.0) )
     1     WRITE(IQ,1800) (MES(I),I=1,I1)
 1800 FORMAT(1X,80A1)
         RETURN
         END
C 
C 
C 
C 
C 
        SUBROUTINE ZTIME(I)
        save
        J=1
        J=7-I+J
CCCC    I=MRUN(J)
        RETURN
        END
  
