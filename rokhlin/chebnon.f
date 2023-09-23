        implicit real *8 (a-h,o-z)
        dimension ts(2000),f(2000),coefs(40),
     1      errs(2000),w(40000),f2(2000),jjs(40)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
         PRINT *, 'ENTER n'
         READ *,n
         CALL PRINF('n=*',n,1 )
c 
         PRINT *, 'ENTER m'
         READ *,m
         CALL PRINf('m=*',m,1 )
c 
c        construct the nodes where the function to be interpolated
c        is tabulated, and the values of the function
c 
        do 1200 i=1,n
        ts(i)=i
        f(i)=1+ts(i)**2-7*(ts(i)*0.1d0)**3
 1200 continue
c 
c       construct an interpolating polynomial
c 
        call chebnon(ier,ts,f,n,m,
     1      coefs,u,v,cond,w,errs,errmax,jjs)
        call prinf('after chebnon, ier=*',ier,1)
        call prin2('after chebnon, errmax=*',errmax,1)
        call prin2('and cond=*',cond,1)
        call prin2('and coefs=*',coefs,m+1)
        call prinf('and jjs=*',jjs,m+1)
c 
c       calculate the expansion at the original points
c 
        do 2200 i=1,n
        call chebnev(ts(i),u,v,coefs,m,f2(i))
        w(i)=f2(i)-f(i)
 2200 continue
         call prin2('and recomputed f are*',f2,n)
         call prin2('while original f are*',f,n)
         call prin2('and the differences are*',w,n)
        stop
        end
c 
c 
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c          this is the end of the debugging code and the beginning
c          of the actual approximation routine
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c 
c 
        subroutine chebnon(ier,ts,f,n,m,
     1      coefs,u,v,cond,w,errs,errmax,jjs)
        implicit real *8 (a-h,o-z)
        save
        dimension ts(1),f(1),coefs(1),w(1),errmax(1),jjs(1)
c 
c        this subroutine constructs a polynomial approximation
c        to a function tabulated at a bunch of nodes. the nodes
c        do not have to be equispaced, and do not have to be in
c        an increasing order. subsequently to a call to this
c        subroutine, the approximation can be evaluated at an
c        arbitrary point by a call to the entry chebnev of this
c        subroutine (see below).
c 
c    explanation: the function f(t) is approximated by the
c    expression
c 
c        f(t)= \sum_{k=0}^{m-1} T_k (u \cdot t + v),              (1)
c 
c    with the coefficients u,v and coefs(k), k=1,2, ...,m-1
c    produced by this subroutine
c 
c                      input parameters:
c 
c  ts - the nodes where the user-specified function is tabulated
c  f - the values of the user-specified function at the points ts
c  n - the number of elements in each of the arrays ts, f
c  m - the order of the approximation
c 
c                      output parameters:
c 
c  ier - error return code.
c        ier=0 means successful conclusion
c        ier=16 means a truly terrible distribution of points ts.
c          in the opinion of the subroutine, an approximation like
c          this makes no sense.
c  coefs - the m+1 coefficients of the m-th order chebychev expansion
c         approximating the user-supplied function.
c  u,v - the coefficients in formula (1)
c  cond - the condition number of the (m+1) \times (m+1) linear
c          system that had to be solved to obtain the approximation (1)
c  errs - the array of errors made by the approximation (1) at the
c         nodes ts (n of them, obviously).
c  errmax - the maximal (in absolute value) element in errs
c  jjs - the indices in the arrays ts, f of the elements actually
c         used for the interpolation (m+1 of them, naturally)
c 
c                      work arrays:
c 
c  w - must be at least 2*m**2 + 10*m +20 real *8 elements long
c 
c 
c        . . . allocate memory for all work arrays
c 
        mp=m+1
c 
        ix=1
        lx=mp+1
c 
        iuvx=ix+lx
        luvx=mp+1
c 
        iamatr=iuvx+luvx
        lamatr=mp**2+mp+2
c 
        irhs=iamatr+lamatr
        lrhs=mp+1
c 
        ijjs=irhs+lrhs
        ljjs=mp+1
c 
        iwork=ijjs+ljjs
c 
c       find the coefficients of the chebychev expansion
c       approximating the user-supplied function
c 
        call chebnon0(ier,ts,n,f,mp,coefs,w(ix),w(iuvx),
cccc     1      w(ijjs),w(iamatr),w(irhs),w(iwork),errs,errmax,
     1      jjs,w(iamatr),w(irhs),w(iwork),errs,errmax,
     2      u1,v1,u2,v2,cond)
c 
        u=u2
        v=v2
        return
c 
c 
c 
c 
        entry chebnev(x,u,v,coefs,m,val)
c 
c        this entry evaluates at the point x the value of the  chebychev
c        expansion composed with a linear change of variables, given
c        by the formula
c 
c        val= \sum_{k=0}^{m-1} T_k (u \cdot x + v).              (1)
c 
c                  input parameters:
c 
c  x - the point on the line where the expansion is to be evaluated
c  u, v - the coefficients of a linear function, normally obtained
c        by calling the entry chebnon of this subroutine (see above).
c  coefs - the coefficients of the chebychev expansion, normally
c        obtained by calling the entry chebnon of this subroutine
c        (see above).
c  m - the order of the chebychev expansion to be evaluated
c 
c                  output parameters:
c 
c  val - the value of the expansion (1) at the point x
  
c        . . . evaluate the expansion at the user-specified point x
c 
        xx=u*x+v
        call CHEBNOCH(Xx,VAL,der,coefs,m)
        return
        end
c 
c 
c 
c 
c 
        subroutine chebnon0(ier,ts,n,f,m,coefs,x,uvx,jjs,
     1      amatr,rhs,work,errs,errmax,u1,v1,u2,v2,cond)
        implicit real *8 (a-h,o-z)
        save
        dimension x(20),jjs(20),ts(1),f(1),coefs(1),
     1      x(1),uvx(1),amatr(m,1),rhs(1),work(m,1),errs(1)
c 
c        find the chebychev nodes on the interval [-1,1]
c 
        ZERO=0
        DONE=1
        pi=datan(done)*4
        h=pi/(2*m)
        do 1200 i=1,m
        t=(2*i-1)*h
        x(m-i+1)=dcos(t)
1200  CONTINUE
cccc        call prin2('chebychev nodes as constructed*',x,m)
cccc        call prin2('and ts=*',ts,n)
c 
c       find the largest and the smallest among the nodes ts
c 
        dmin=1.0d20
        dmax=-1.0d20
        do 1300 i=1,n
        if(ts(i) .gt. dmin) goto 1240
        dmin=ts(i)
        imin=i
 1240 continue
c 
        if(ts(i) .lt. dmax) goto 1300
        dmax=ts(i)
        imax=i
 1300 continue
c 
c        find the linear mappings converting the interval [-1,1]
c        into the interval [ts(1),ts(n)], and back
c 
        u1=(dmax-dmin)/2
        v1=(dmax+dmin)/2
c 
c       . . . and back
c 
        u2=2/(dmax-dmin)
        v2=1-dmax*u2
c 
c       construct the chebychev nodes on the interval
c       [ts(1),ts(n)]
c 
        do 1400 i=1,m
        uvx(i)=u1*x(i)+v1
 1400 continue
c 
c       find the nodes in the array ts nearest to the images
c       of the nodes x(i)
c 
        do 1800 i=1,m
        ii=0
        rmin=1.0d20
        do 1600 j=1,n
        d=dabs(uvx(i)-ts(j))
        if(d .gt. rmin) goto 1600
        ii=j
        rmin=d
 1600 continue
c 
        jjs(i)=ii
 1800 continue
cccc        call prinf('array jjs as found*',jjs,m)
c 
c       make sure no two elements of the array jjs coincide
c 
        ier=0
        do 2000 i=1,m
        do 1900 j=1,m
        if(i .eq. j) goto 1900
        if(jjs(i) .eq. jjs(j)) ier=16
 1900 continue
 2000 continue
        if(ier .eq. 16) return
c 
c        find the matrix of values of chebychev polynomials
c        at the points ts(jjs(i))
c 
        do 2400 i=1,m
c 
        d=ts(jjs(i))
        xx=u2*d+v2
        amatr(i,1)=1
        amatr(i,2)=xx
c 
        do 2200 j=2,m-1
        amatr(i,j+1)=2*xx*amatr(i,j)-amatr(i,j-1)
 2200 continue
 2400 continue
cccc        call prin2('amatr as calculated*',amatr,m*m)
c 
c       construct the right-hand side of the linear system
c       to be solved
c 
        do 2600 i=1,m
        rhs(i)=f(jjs(i))
 2600 continue
cccc      call prin2('rhs as calculated*',rhs,m)
c 
c        solve the linear system
c 
        call chebnort(amatr,m,work,cond)
cccc         call prin2('after chebnort, cond=*',cond,1)
c 
        call chebnoap(amatr,m,rhs,coefs)
cccc        call prin2('coefs as calculated*',coefs,m)
c 
c       test the accuracy of the approximation
c 
        errmax=0
        do 3200 i=1,n
        xx=u2*ts(i)+v2
c 
        call CHEBNOCH(Xx,VAL,der,coefs,m-1)
        errs(i)=val-f(i)
        d=dabs(errs(i))
        if(d .gt. errmax) errmax=d
 3200 continue
cccc        call prin2('in chebnon0, errs=*',errs,n)
cccc        call prin2('in chebnon0, errmax=*',errmax,1)
        return
        end
c 
c 
c 
c 
c 
      SUBROUTINE CHEBNOCH(X,VAL,der,TEXP,N)
C 
C     This subroutine computes the value and the derivative
c     of a chebychev expansion with coefficients TEXP
C     at point X in interval [-1,1]
C 
c                input parameters:
c 
C     X = evaluation point
C     TEXP = expansion coefficients
C     N  = order of expansion (0 ,..., N-1)
c 
c                output parameters:
c 
C     VAL = computed value
C     der = computed value of the derivative
C 
      IMPLICIT REAL *8 (A-H,O-Z)
        save
      REAL *8 TEXP(1)
C 
        done=1
        tjm2=1
        tjm1=x
        derjm2=0
        derjm1=1
c 
        val=texp(1)*tjm2+texp(2)*tjm1
        der=texp(2)
c 
      DO 600 J = 2,N
c 
        tj=2*x*tjm1-tjm2
        val=val+texp(j+1)*tj
c 
        derj=2*tjm1+2*x*derjm1-derjm2
        der=der+texp(j+1)*derj
c 
        tjm2=tjm1
        tjm1=tj
        derjm2=derjm1
        derjm1=derj
 600   CONTINUE
c 
      RETURN
      END
c 
c 
c 
c 
c 
        subroutine chebnort(a,n,work,cond)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(1),work(1)
c 
c        this subroutine inverts a user-supplied real matrix
c        a by means of the gram-schmidt process. in fact, this
c        is a memory management routine for the routine chebnor1,
c        which performs the actual inversion
c 
c               input parameters:
c 
c  a - the n*n - matrix to be inverted (destroyed by the
c         subroutine)
c  n - the dimensionality of the matrix a
c 
c               output parameters:
c 
c  a - the inverse of a on input
c 
c               work array
c 
c  work - must be at least  n * (n+1)+2 real *8 locations long
c 
c 
c        . . . allocate memory for the matrix inverter
c 
        ib=1
        lb=n*n
c 
        iw=ib+lb+1
c 
c       invert the matrix
c 
        call chebnor1(a,work(ib),n,work(iw),cond)
        return
        end
c 
c 
c 
c 
c 
        subroutine chebnor1(a,b,n,work,cond)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(n,n),b(n,n),cd,zero,
     1      work(1)
c 
c       set the matrix b to unity
c 
        done=1
        zero=0
        do 1400 i=1,n
        do 1200 j=1,n
        b(j,i)=zero
 1200 continue
        b(i,i)=1
 1400 continue
c 
c        conduct the gram-schmidt process
c 
        dnormax=-1
        dnormin=1.0d20
        do 4000 i=1,n
c 
c        normalize the i-th row
c 
        d=0
        do 2200 j=1,n
        d=d+a(i,j)**2
 2200 continue
c 
        d=done/dsqrt(d)
        if(d .gt. dnormax) dnormax=d
        if(d .lt. dnormin) dnormin=d
        do 2400 j=1,n
        a(i,j)=a(i,j)*d
        b(i,j)=b(i,j)*d
 2400 continue
c 
c       orthogonalize all subsequent rows to the i-th one
c 
        if(i .eq. n) goto 4000
         do 3200 ijk=1,2
        do 3000 j=i+1,n
c 
        cd=0
        do 2600 k=1,n
        cd=cd+a(i,k)*a(j,k)
 2600 continue
c 
        do 2800 k=1,n
        a(j,k)=a(j,k)-cd*a(i,k)
        b(j,k)=b(j,k)-cd*b(i,k)
 2800 continue
 3000 continue
 3200 continue
 4000 continue
        cond=dnormax/dnormin
c 
c       now, multiply the abjoint of the resulting
c       orthogonal matrix by the triangular one,
c       obtaining the inverse of the original a
c 
        do 5000 i=1,n
        do 4800 j=1,n
        cd=0
        do 4200 k=1,n
        cd=cd+a(k,i)*b(k,j)
 4200 continue
        work(j)=cd
 4800 continue
        do 4900 j=1,n
        a(j,i)=work(j)
 4900 continue
 5000 continue
c 
c        now, transpose a
c 
        do 5400 i=2,n
        do 5200 j=1,i
        cd=a(i,j)
        a(i,j)=a(j,i)
        a(j,i)=cd
 5200 continue
 5400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine chebnoap(a,n,x,y)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(n,n),x(1),y(1),cd
c 
c        apply the matrix a to the vector x obtaining y
c 
        do 1400 i=1,n
        cd=0
        do 1200 j=1,n
        cd=cd+a(i,j)*x(j)
 1200 continue
        y(i)=cd
 1400 continue
        return
        end
  
  
  
  
  
