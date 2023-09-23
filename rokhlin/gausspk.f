c 
c enclosed is a test program, and then the adaptive gaussian quadrature
c   function itself.
c 
      program test
      real*8 t(2000),w(2000),stack(2000),dint,eps,defint,a,c
      integer n,status,idefint,junk
      external fun
c 
c 
c       create the input parameters: interval is (0,1)
c       # of gaussian weights n is user-specified
c       absolute error eps is user-specified
c 
      a = -1
      c = 1
c 
      print *, 'enter (integer) number of gaussian weights,',
     $     ' (real) absolute error'
      read *,n,eps
      write(6,19) n,eps
 19   format('# of weights= ',i4,', absolute error=',e8.2)
c 
c       construct gaussian quadratures: zero is always returned
c 
      junk=idefint(stack,t,w,n)
c 
c       construct the test function and integrate it
c 
      dint=defint(t,w,fun,a,c,eps,n,status,stack)
c 
      if (status .eq. 0) then
	write(6,99) dint
99	format ('Integral is ',e20.14)
      else
         print *,"defint returned error"
      endif
c 
      stop
      end
c
      real*8 function fun(t)
      real*8 pi,t
      integer init
      data init/0/
c 
      if (init.eq.0) then
         pi = 4 * datan(1.D0)
         init = 1
      endif
c 
      fun=1
c$$$          fun = 1.D0/dsqrt(dabs(t))
c$$$      fun = t**-.9
c$$$      fun = dsin(1/t)
c$$$      fun=dexp(-dabs(10.*t))
c$$$      fun = dlog(dabs(t))
c$$$      fun=15.5*(t**30)
c$$$      fun=dcos(201.*pi*t/2.)
      return
      end
c
c compute definite integral of a function using adaptive gaussian quadrature
c 
c  IMPORTANT NOTE:  Before defint is called, the entry point
c		idefint(stack,t,w,n)
c    must be called to initialize the vectors t and w, with n specifying
c    the number of gaussian weights desired, and with stack a temporary
c    vector of size n.  t and w only have to be initialized once, provided
c    that the value of n provided for idefint is the same value provided
c    for each call of defint.
c 
c    The entry point idefint is an integer function which always returns
c    the value 0.
c 
c    As is well-known, a gaussian quadrature rule using
c    n weights integrates all polynomials of order up to 2n-1 exactly.
c 
c  input parameters:
c 
c  t (real*8) = vector of quadrature roots (on interval -1,1): size n vector
c  w (real*8) = vector of weights: size n vector
c  n (integer) = # of pts in quadrature rule
c 
c  fun = real*8 function, called in the form: result=fun(x), with x
c    a real*8 point in the interval (a,c)
c 
c  a (real*8) = left end of interval of interest
c  c (real*8) = right end of interval
c 
c  eps (real*8) = desired accuracy (to absolute precision)
c 
c output parameters:
c  status (integer) =
c     returns -1 if an error occurred
c     otherwise, returns 0
c 
c Technical note follows, and can be ignored:
c  (-1 is returned either when some of the data got trashed, due
c   to not passing large enough arrays, or when the eps requested
c   is way, way smaller in magnitude than the actual integral,
c   by a factor of 10**-30)
c 
c temporaries required:
c  stack = real scratch array
c    when defint is called, stack should be size 700 at least
c    when idefint is called, stack should be size n at least,
c    where n is the number of pts in the quadrature rule
c 
      real*8 function defint(t,w,fun,a,c,eps,n,status,stack)
      external fun
      real*8 t(n),w(n),a,c,eps,sum,xc,xl,fun,qxc,qxl,qsum,
     $     stack(7,100),err,tsum
      integer n,i,status,top,rts(100),temp(2),idefint
c 
c technical note:
c    for each stack frame:
c      1 = sum1 (integral of left half)
c      2 = sum2 (integral of right half)
c      3 = sxl  (radius of each half)
c      4 = xc1  (center of left half)
c      5 = xc2  (center of right half)
c      6 = ret  (return value)
c      7 = temp
c 
c    also: rts array simulates return address
c 
      top=1
      qxc = (a+c)/2
      qxl = (c-a)/2
c 
      qsum=0
      status=0
  
      do 100 i=1,n
         qsum= qsum+ w(i)*fun(qxc+qxl*t(i))
 100  continue
      qsum=qxl*qsum
c 
c now simulate recursive function calls
c  rts = return type: 1 = called by left traverse, 2=called by
c    right traverse, 3 = initial call
c 
      rts(top)=3
      xc=qxc
      xl=qxl
      sum=qsum
c 
c start of "recursive" function
c 
 200  continue
c 
c$$$      print *,top,"----------"
c$$$      print *,top,"At top of rec function"
c$$$      print *,top,"xc=",xc,", xl=",xl,", sum=",sum,", top=",top
c 
      stack(3,top)= xl/2
      stack(4,top)= xc-stack(3,top)
      stack(5,top)= xc+stack(3,top)
c 
c$$$      print *,top,"sxl=",stack(3,top),", xc1=",stack(4,top),", xc2=",
c$$$     $     stack(5,top)
c 
      stack(1,top)=0
  
      do 210 i=1,n
         stack(1,top) = stack(1,top) + w(i)*fun(stack(4,top)+
     $        stack(3,top)*t(i))
 210  continue
      stack(1,top) = stack(3,top)*stack(1,top)
  
      stack(2,top)=0
  
      do 220 i=1,n
         stack(2,top) = stack(2,top) + w(i)*fun(stack(5,top)+
     $        stack(3,top)*t(i))
 220  continue
      stack(2,top) = stack(3,top)*stack(2,top)
  
      tsum=stack(1,top)+stack(2,top)
c 
c$$$      print *,top,"sum1=", stack(1,top),", sum2=",stack(2,top),
c$$$     $     ", tsum=",tsum
c 
      err=dabs(sum-stack(1,top)-stack(2,top))
c 
c$$$      print *,top,"err = ",err
c 
      if (err .lt. eps) then
         stack(6,top) = stack(1,top)+stack(2,top)
c 
c$$$         print *,top,"return val = ",stack(6,top),", rts=",rts(top)
c 
c simulate a return
c 
         if (rts(top).eq.1) goto 300
         if (rts(top).eq.2) goto 400
         if (rts(top).eq.3) goto 500
c 
c a technical error: bad stack return value
c 
         status=-1
         return
      else
c 
c simulate a call
c 
c$$$         print *,top,"recursive call: returning to 1"
c 
         rts(top+1)=1
         xc = stack(4,top)
         xl = stack(3,top)
         sum = stack(1,top)
         top=top+1
         if (top .gt. 100) then
            status=-1
            return
         endif
         goto 200
      endif
c 
c first return location
c 
 300  continue
      top=top-1
c 
c$$$      print *,top,"returning to 1"
c 
      stack(7,top)=stack(6,top+1)
c 
c$$$      print *,top,"temp = ",stack(7,top)
c 
c simulate a call
c 
c$$$      print *,top,"recursive call: returning to 2"
c 
      rts(top+1)=2
      xc=stack(5,top)
      xl=stack(3,top)
      sum = stack(2,top)
c 
      top=top+1
      if (top .gt. 100) then
c 
c stack overflow: set status to -1 and return
c 
         status=-1
         return
      endif
      goto 200
c 
c second return location
c 
 400  continue
      top = top-1
c 
c$$$      print *,top,"returning to 2"
c 
      stack(6,top) = stack(7,top) + stack(6,top+1)
c 
c$$$      print *,top,"return value = ",stack(6,top)
c 
c simulate a return
c 
      if (rts(top) .eq.1) goto 300
      if (rts(top) .eq.2) goto 400
      if (rts(top) .eq.3) goto 500
c 
c a technical error: bad stack return value
      status=-1
      return
c 
c final return location: we're done
c 
 500  continue
c 
c$$$      print *,top,"Final return"
c 
      top=top-1
      if (top .ne. 0) then
c 
c a technical error: rval = 3 (but not last elt on stack)
c 
         status=-1
         return
      endif
c 
      defint = stack(6,1)
      return
c 
c entry point to initialize gaussian quadratures
c 
      entry idefint(stack,t,w,n)
      i=1
      top=0
c 
c kind=1, alpha=beta=kpts=0
c 
      call gaussq(i,n,top,top,top,temp,stack,t,w)
      idefint=0
      return
      end
c
c          data set gaussq     at level 001 as of 07/23/82
c                                1/20/75
      subroutine gaussq(kind, n, alpha, beta, kpts, endpts, b, t, w)
c 
c           this set of routines computes the nodes t(j) and weights
c        w(j) for gaussian-type quadrature rules with pre-assigned
c        nodes.  these are used when one wishes to approximate
c 
c                 integral (from a to b)  f(x) w(x) dx
c 
c                              n
c        by                   sum w  f(t )
c                             j=1  j    j
c 
c        (note w(x) and w(j) have no connection with each other.)
c        here w(x) is one of six possible non-negative weight
c        functions (listed below), and f(x) is the
c        function to be integrated.  gaussian quadrature is particularly
c        useful on infinite intervals (with appropriate weight
c        functions), since then other techniques often fail.
c 
c           associated with each weight function w(x) is a set of
c        orthogonal polynomials.  the nodes t(j) are just the zeroes
c        of the proper n-th degree polynomial.
c 
c     input parameters (all real numbers are in double precision)
c 
c        kind     an integer between 1 and 6 giving the type of
c                 quadrature rule:
c 
c        kind = 1:  legendre quadrature, w(x) = 1 on (-1, 1)
c        kind = 2:  chebyshev quadrature of the first kind
c                   w(x) = 1/sqrt(1 - x*x) on (-1, +1)
c        kind = 3:  chebyshev quadrature of the second kind
c                   w(x) = sqrt(1 - x*x) on (-1, 1)
c        kind = 4:  hermite quadrature, w(x) = exp(-x*x) on
c                   (-infinity, +infinity)
c        kind = 5:  jacobi quadrature, w(x) = (1-x)**alpha * (1+x)**
c                   beta on (-1, 1), alpha, beta .gt. -1.
c                   note: kind=2 and 3 are a special case of this.
c        kind = 6:  generalized laguerre quadrature, w(x) = exp(-x)*
c                   x**alpha on (0, +infinity), alpha .gt. -1
c 
c        n        the number of points used for the quadrature rule
c        alpha    real parameter used only for gauss-jacobi and gauss-
c                 laguerre quadrature (otherwise use 0.d0).
c        beta     real parameter used only for gauss-jacobi quadrature--
c                 (otherwise use 0.d0)
c        kpts     (integer) normally 0, unless the left or right end-
c                 point (or both) of the interval is required to be a
c                 node (this is called gauss-radau or gauss-lobatto
c                 quadrature).  then kpts is the number of fixed
c                 endpoints (1 or 2).
c        endpts   real array of length 2.  contains the values of
c                 any fixed endpoints, if kpts = 1 or 2.
c        b        real scratch array of length n
c 
c     output parameters (both double precision arrays of length n)
c 
c        t        will contain the desired nodes.
c        w        will contain the desired weights w(j).
c 
c     subroutines required
c 
c        solve, class, and imtql2 are provided.  underflow may sometimes
c        occur, but it is harmless if the underflow interrupts are
c        turned off.  to do this, the first call of the main program
c        should be
c                  call traps (0, 0, 10000, 0, 0)    in watfiv
c        or
c                  call init                         in fortran g or h.
c 
c     accuracy
c 
c        the routine was tested up to n = 512 for legendre quadrature,
c        up to n = 136 for hermite, up to n = 68 for laguerre, and up
c        to n = 10 or 20 in other cases.  in all but two instances,
c        comparison with tables in ref. 3 showed 12 or more significant
c        digits of accuracy.  the two exceptions were the weights for
c        hermite and laguerre quadrature, where underflow caused some
c        very small weights to be set to zero.  this is, of course,
c        completely harmless.
c 
c     method
c 
c           the coefficients of the three-term recurrence relation
c        for the corresponding set of orthogonal polynomials are
c        used to form a symmetric tridiagonal matrix, whose
c        eigenvalues (determined by the implicit ql-method with
c        shifts) are just the desired nodes.  the first components of
c        the orthonormalized eigenvectors, when properly scaled,
c        yield the weights.  this technique is much faster than using a
c        root-finder to locate the zeroes of the orthogonal polynomial.
c        for further details, see ref. 1.  ref. 2 contains details of
c        gauss-radau and gauss-lobatto quadrature only.
c 
c     references
c 
c        1.  golub, g. h., and welsch, j. h., "calculation of gaussian
c            quadrature rules," mathematics of computation 23 (april,
c            1969), pp. 221-230.
c        2.  golub, g. h., "some modified matrix eigenvalue problems,"
c            siam review 15 (april, 1973), pp. 318-334 (section 7).
c        3.  stroud and secrest, gaussian quadrature formulas, prentice-
c            hall, englewood cliffs, n.j., 1966.
c 
c     ..................................................................
c 
        save
      real*8 b(n), t(n), w(n), endpts(2), muzero, t1,
     x gam, solve, dsqrt, alpha, beta
      integer kind,n,kpts,i,ierr
c 
      call class (kind, n, alpha, beta, b, t, muzero)
c 
c           the matrix of coefficients is assumed to be symmetric.
c           the array t contains the diagonal elements, the array
c           b the off-diagonal elements.
c           make appropriate changes in the lower right 2 by 2
c           submatrix.
c 
      if (kpts.eq.0)  go to 100
      if (kpts.eq.2)  go to  50
c 
c           if kpts=1, only t(n) must be changed
c 
      t(n) = solve(endpts(1), n, t, b)*b(n-1)**2 + endpts(1)
      go to 100
c 
c           if kpts=2, t(n) and b(n-1) must be recomputed
c 
   50 gam = solve(endpts(1), n, t, b)
      t1 = ((endpts(1) - endpts(2))/(solve(endpts(2), n, t, b) - gam))
      b(n-1) = dsqrt(t1)
      t(n) = endpts(1) + gam*t1
c 
c           note that the indices of the elements of b run from 1 to n-1
c           and thus the value of b(n) is arbitrary.
c           now compute the eigenvalues of the symmetric tridiagonal
c           matrix, which has been modified as necessary.
c           the method used is a ql-type method with origin shifting
c 
  100 w(1) = 1.0d0
      do 105 i = 2, n
  105    w(i) = 0.0d0
c 
      call imtql2 (n, t, b, w, ierr)
      do 110 i = 1, n
  110    w(i) = muzero * w(i) * w(i)
c 
      return
      end
c 
c 
c
      real*8 function solve(shift, n, a, b)
c 
c       this procedure performs elimination to solve for the
c       n-th component of the solution delta to the equation
c 
c             (jn - shift*identity) * delta  = en,
c 
c       where en is the vector of all zeroes except for 1 in
c       the n-th position.
c 
c       the matrix jn is symmetric tridiagonal, with diagonal
c       elements a(i), off-diagonal elements b(i).  this equation
c       must be solved to obtain the appropriate changes in the lower
c       2 by 2 submatrix of coefficients for orthogonal polynomials.
c 
c 
      real*8 shift, a(n), b(n), alpha
      integer n,i,nm1
c 
      alpha = a(1) - shift
      nm1 = n - 1
      do 10 i = 2, nm1
   10    alpha = a(i) - shift - b(i-1)**2/alpha
      solve = 1.0d0/alpha
      return
      end
c 
c 
c
      subroutine class(kind, n, alpha, beta, b, a, muzero)
c 
c           this procedure supplies the coefficients a(j), b(j) of the
c        recurrence relation
c 
c             b p (x) = (x - a ) p   (x) - b   p   (x)
c              j j            j   j-1       j-1 j-2
c 
c        for the various classical (normalized) orthogonal polynomials,
c        and the zero-th moment
c 
c             muzero = integral w(x) dx
c 
c        of the given polynomial's weight function w(x).  since the
c        polynomials are orthonormalized, the tridiagonal matrix is
c        guaranteed to be symmetric.
c 
c           the input parameter alpha is used only for laguerre and
c        jacobi polynomials, and the parameter beta is used only for
c        jacobi polynomials.  the laguerre and jacobi polynomials
c        require the gamma function.
c 
c     ..................................................................
c 
        save
      real*8 a(n), b(n), muzero, alpha, beta
      real*8 abi, a2b2, dgamma, pi, dsqrt, ab
      integer kind,n,nm1,i
      data pi / 3.141592653589793d0/
c 
      nm1 = n - 1
      go to (10, 20, 30, 40, 50, 60), kind
c 
c              kind = 1:  legendre polynomials p(x)
c              on (-1, +1), w(x) = 1.
c 
   10 muzero = 2.0d0
      do 11 i = 1, nm1
         a(i) = 0.0d0
         abi = i
   11    b(i) = abi/dsqrt(4*abi*abi - 1.0d0)
      a(n) = 0.0d0
      return
c 
c              kind = 2:  chebyshev polynomials of the first kind t(x)
c              on (-1, +1), w(x) = 1 / sqrt(1 - x*x)
c 
   20 muzero = pi
      do 21 i = 1, nm1
         a(i) = 0.0d0
   21    b(i) = 0.5d0
      b(1) = dsqrt(0.5d0)
      a(n) = 0.0d0
      return
c 
c              kind = 3:  chebyshev polynomials of the second kind u(x)
c              on (-1, +1), w(x) = sqrt(1 - x*x)
c 
   30 muzero = pi/2.0d0
      do 31 i = 1, nm1
         a(i) = 0.0d0
   31    b(i) = 0.5d0
      a(n) = 0.0d0
      return
c 
c              kind = 4:  hermite polynomials h(x) on (-infinity,
c              +infinity), w(x) = exp(-x**2)
c 
   40 muzero = dsqrt(pi)
      do 41 i = 1, nm1
         a(i) = 0.0d0
   41    b(i) = dsqrt(i/2.0d0)
      a(n) = 0.0d0
      return
c 
c              kind = 5:  jacobi polynomials p(alpha, beta)(x) on
c              (-1, +1), w(x) = (1-x)**alpha + (1+x)**beta, alpha and
c              beta greater than -1
c 
   50 ab = alpha + beta
      abi = 2.0d0 + ab
      muzero = 2.0d0 ** (ab + 1.0d0) * dgamma(alpha + 1.0d0) * dgamma(
     x beta + 1.0d0) / dgamma(abi)
      a(1) = (beta - alpha)/abi
      b(1) = dsqrt(4.0d0*(1.0d0 + alpha)*(1.0d0 + beta)/((abi + 1.0d0)*
     1  abi*abi))
      a2b2 = beta*beta - alpha*alpha
      do 51 i = 2, nm1
         abi = 2.0d0*i + ab
         a(i) = a2b2/((abi - 2.0d0)*abi)
   51    b(i) = dsqrt (4.0d0*i*(i + alpha)*(i + beta)*(i + ab)/
     1   ((abi*abi - 1)*abi*abi))
      abi = 2.0d0*n + ab
      a(n) = a2b2/((abi - 2.0d0)*abi)
      return
c 
c              kind = 6:  laguerre polynomials l(alpha)(x) on
c              (0, +infinity), w(x) = exp(-x) * x**alpha, alpha greater
c              than -1.
c 
   60 muzero = dgamma(alpha + 1.0d0)
      do 61 i = 1, nm1
         a(i) = 2.0d0*i - 1.0d0 + alpha
   61    b(i) = dsqrt(i*(i + alpha))
      a(n) = 2.0d0*n - 1 + alpha
      return
      end
c 
c 
c 
c
c 
        real*8 function dgamma(x)
        real*8 x
        dgamma=1
        return
        end
c
c 
      subroutine imtql2(n, d, e, z, ierr)
c 
c     this subroutine is a translation of the algol procedure imtql2,
c     num. math. 12, 377-383(1968) by martin and wilkinson,
c     as modified in num. math. 15, 450(1970) by dubrulle.
c     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
c     this is a modified version of the 'eispack' routine imtql2.
c 
c     this subroutine finds the eigenvalues and first components of the
c     eigenvectors of a symmetric tridiagonal matrix by the implicit ql
c     method.
c 
c     on input:
c 
c        n is the order of the matrix;
c 
c        d contains the diagonal elements of the input matrix;
c 
c        e contains the subdiagonal elements of the input matrix
c          in its first n-1 positions.  e(n) is arbitrary;
c 
c        z contains the first row of the identity matrix.
c 
c      on output:
c 
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1, 2, ..., ierr-1;
c 
c        e has been destroyed;
c 
c        z contains the first components of the orthonormal eigenvectors
c          of the symmetric tridiagonal matrix.  if an error exit is
c          made, z contains the eigenvectors associated with the stored
c          eigenvalues;
c 
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c 
c     ------------------------------------------------------------------
c 
        save
      integer i, j, k, l, m, n, ii, mml, ierr
      real*8 d(n), e(n), z(n), b, c, f, g, p, r, s, machep
      real*8 dsqrt, dabs, dsign
c 
c     :::::::::: machep is a machine dependent parameter specifying
c                the relative precision of floating point arithmetic.
c                machep = 16.0d0**(-13) for long form arithmetic
c                on s360 ::::::::::
cccc  data machep/z3410000000000000/
      data machep/1.0d-14/
c 
      ierr = 0
      if (n .eq. 1) go to 1001
c 
      e(n) = 0.0d0
      do 240 l = 1, n
         j = 0
c     :::::::::: look for small sub-diagonal element ::::::::::
  105    do 110 m = l, n
            if (m .eq. n) go to 120
            if (dabs(e(m)) .le. machep * (dabs(d(m)) + dabs(d(m+1))))
     x         go to 120
  110    continue
c 
  120    p = d(l)
         if (m .eq. l) go to 240
         if (j .eq. 30) go to 1000
         j = j + 1
c     :::::::::: form shift ::::::::::
         g = (d(l+1) - p) / (2.0d0 * e(l))
         r = dsqrt(g*g+1.0d0)
         g = d(m) - p + e(l) / (g + dsign(r, g))
         s = 1.0d0
         c = 1.0d0
         p = 0.0d0
         mml = m - l
c 
c     :::::::::: for i=m-1 step -1 until l do -- ::::::::::
         do 200 ii = 1, mml
            i = m - ii
            f = s * e(i)
            b = c * e(i)
            if (dabs(f) .lt. dabs(g)) go to 150
            c = g / f
            r = dsqrt(c*c+1.0d0)
            e(i+1) = f * r
            s = 1.0d0 / r
            c = c * s
            go to 160
  150       s = f / g
            r = dsqrt(s*s+1.0d0)
            e(i+1) = g * r
            c = 1.0d0 / r
            s = s * c
  160       g = d(i+1) - p
            r = (d(i) - g) * s + 2.0d0 * c * b
            p = s * r
            d(i+1) = g + p
            g = c * r - b
c     :::::::::: form first component of vector ::::::::::
            f = z(i+1)
            z(i+1) = s * z(i) + c * f
  200       z(i) = c * z(i) - s * f
c 
         d(l) = d(l) - p
         e(l) = g
         e(m) = 0.0d0
         go to 105
  240 continue
c 
c     :::::::::: order eigenvalues and eigenvectors ::::::::::
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
c 
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
c 
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
         p = z(i)
         z(i) = z(k)
         z(k) = p
  300 continue
c 
      go to 1001
c     :::::::::: set error -- no convergence to an
c                eigenvalue after 30 iterations ::::::::::
 1000 ierr = l
 1001 return
c     :::::::::: last card of imtql2 ::::::::::
      end
