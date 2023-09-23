        implicit real *8 (a-h,o-z)
        dimension f(1000),u(200 000),v(200 000),
     1      whts(1000),zeroes(10 000),tt(1000),
     2      w(100 000)
c 
        call prini(6,13)
c 
C 
C       SET ALL PARAMETERS
C 
C 
        PRINT *, 'ENTER n'
        READ *,n
        CALL PRINf('n=*',n,1 )
C 
c       construct the test function whose zeroes are to be found
c 
        a=0
        b=1
        large=100
c 
         itype=0
         call legeexps(itype,n,tt,u,v,whts)
c 
         do 1100 i=1,n
         tt(i)=(b-a)/2*tt(i)+(b+a)/2
 1100 continue
c 
         call prin2('tt as created*',tt,n)
c 
        hh=1
        hh=hh/10
        done=1
        pi=datan(done)*4
        do 1200 i=1,n
ccc        f(i)=(tt(i)-hh)*(tt(i)-hh*2)*(tt(i)-hh*3)
        f(i)=dcos(tt(i)*10*pi)
 1200 continue
          call prin2('f as created*',f,n)
c 
c 
c        now, test the properly packaged subroutine
c 
        ifconstr=1
        eps=1.0d-10
        lw=100000
        call legexrts(ier,a,b,f,n,eps,large,ifconstr,
     1      zeroes,nzeroes,w,lw)
c 
        call prinf('after first legexrts, ier=*',ier,1)
        call prin2('after legexrts, zeroes=*',zeroes,nzeroes)
c 
        iffirst=0
        call legexrts(ier,a,b,f,n,eps,large,iffirst,
     1      zeroes,nzeroes,w,lw)
c 
        call prinf('after second legexrts, ier=*',ier,1)
        call prin2('after second legexrts, zeroes=*',zeroes,nzeroes)
  
  
  
  
  
        stop
        end
c 
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c      this is the end of debugging code and the beginning
c      of the root-finding subroutine
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        subroutine legexrts(ier,a,b,f,n,eps,large,iffirst,
     1      zeroes,nzeroes,w,lw)
        implicit real *8 (a-h,o-z)
        save
        dimension f(1),zeroes(1),w(1)
c 
c        this subroutine attempts to find all roots of a
c        user-specified function located on the interval [a,b]
c        the function is specified by its values tabulated at the
c        gaussian nodes on that same interval [a,b]. evaluating this
c        function at these nodes is the user's responsibility; however,
c        the gaussian nodes on the interval [-1,1] can be obtained by
c        calling the subroutine legeexps (see).
c 
c                     input parameters:
c 
c  a,b - the ends of the interval on which the roots of the function
c        are to be looked for.
c  f - the table of values of the function (whose roots are to be found);
c        these values are supposed to be at the gaussian nodes on the
c        interval [a,b]
c  eps - the relative accuracy with which the roots will be found
c  large - the number of points at which the interval [a,b] will be sampled
c        in search of the zero-crossing. this number should be LARGE.
c  iffirst - the parameter telling the subroutine if this is the first
c        time it is being called with this n. if this is the first call
c        (iffirst=1), the subroutine will generate certain internally
c        used matrices, at a cost in CPU time. if this is not the first
c        call (iffirst=0), the CPU time will be saved, since these matrices
c        have already been generated. note that if the the array w has been
c        changed between the calls, the parameter iffirst haas to be reset
c        to 1.
c  lw - the length (in real *8 words) of the work array w (see below)
c 
c                     output parameters:
c 
c  zeroes - the array of roots of the function f on the interval [a,b]
c  nzeroes - the number of elements in array zeroes (the number of roots
c        found on the inetrval [a,b])
c 
c                     work arrays:
c 
c  w - must be at least max(2*n**2+5*n+20, n**2+large+5*n+20)
c        real *8 locations long. it is better to not change the contents
c        of w between calls to this subroutine with the same n, since
c        in this case the parameter iffirst can be set to zero for
c        the subsequent calls, saving some CPU time.
c 
c        . . . allocate memory for the search for the zeroes
c 
        ier=0
        iu=1
        lu=n**2+2
c 
        iv=iu+lu
        lv=n**2+2
        if(large .gt. n**2) lv=large+2
c 
        it=iv+lv
        lt=n+2
c 
        iwhts=it+lt
        lwhts=n+2
c 
        isuspect=iwhts+lwhts
        lsuspect=n+2
c 
        itexp=isuspect+lsuspect
        ltexp=n+2
c 
        ltot=itexp+ltexp
        if(ltot .le. lw) goto 1200
        ier=16
        return
 1200 continue
c 
c       find the roots of the user-supplied function f
c 
        call legexrt0(ier,f,n,iffirst,large,w(iu),w(iv),w(it),
     1     w(iwhts),w(isuspect),zeroes,nzeroes,a,b,w(itexp),eps)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine legexrt0(ier,f,n,ifconstr,large,u,v,t,whts,
     1      suspect,zeroes,ncross,a,b,texp,eps)
        implicit real *8 (a-h,o-z)
        save
        dimension f(1),u(1),v(1),t(1),whts(1),suspect(1),
     1      zeroes(1)
c 
c        this subroutine finds the roots of a polynomial
c        tabulated at gaussian nodes on the interval [a,b]
c 
        ier=0
c 
c        if needed, construct the matrix converting the
c        values of the legendre expansion into the coefficients
c 
         if(ifconstr .eq. 0) goto 1200
c 
         itype=2
         call legeexps(itype,n,t,u,v,whts)
 1200 continue
c 
c        . . . construct the legendre expansion of f
c 
         call legexmtv(u,f,texp,n)
cccc         call prin2('legendre expansion of f is*',texp,n)
c 
c        find all zero crossings
c 
         if(ifconstr .eq. 0) goto 1300
         itype=0
         call legeexps(itype,large,v,udummy,vdummy,whtdummy)
cccccc          call prin2('many nodes are*',v,large)
 1300 continue
c 
         call legeFDER(v(1),f2,der,texp,N)
c 
         ncross=0
         do 1400 i=2,large
c 
        f1=f2
c 
         call legeFDER(v(i),f2,der,texp,N)
c 
        d=f1*f2
        if(d .gt. 0) goto 1400
        ncross=ncross+1
        suspect(ncross)=(f2*v(i-1)-f1*v(i))/(f2-f1)
 1400 continue
cccc         call prin2('crossings are*',suspect,ncross)
c 
c       tighten the found zero crossings by newton
c 
        if(ncross .eq. 0) return
c 
        do 2000 i=1,ncross
c 
        tk=suspect(i)
        ifstop=0
c 
        do 1600 j=1,10
        call legeFDER(tk,fk,derk,texp,N)
        tkp1=tk-fk/derk
c 
cccc        call prin2('tkp1-tk=*',tkp1-tk,1)
        if(ifstop .eq. 1) goto 1800
        if(dabs(tkp1-tk) .lt. eps) ifstop=1
        tk=tkp1
 1600 continue
        ier=8
 1800 continue
        zeroes(i)=tkp1
 2000 continue
c 
c       shift the zeroes onto the proper interval
c 
        uu=(b-a)/2
        vv=(b+a)/2
        do 2200 i=1,ncross
        zeroes(i)=zeroes(i)*uu+vv
 2200 continue
  
        return
        end
c 
c 
c 
c 
c 
        subroutine legexmtv(a,x,y,n)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),x(1),y(1)
c 
c       apply the matrix a to vector x getting y
c 
        do 1400 i=1,n
        d=0
        do 1200 j=1,n
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
c 
c 
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        this is the end of the debugging code and the beginning
c        of the legendre expansion routine
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c 
c 
        subroutine legeexps(itype,n,x,u,v,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),whts(1),u(n,n),v(n,n)
  
c 
c         this subroutine constructs the gaussiaqn nodes
c         on the interval [-1,1], and the weights for the
c         corresponding order n quadrature. it also constructs
c         the matrix v converting the coefficients
c         of a legendre expansion into its values at the n
c         gaussian nodes, and its inverse u, converting the
c         values of a function at n gaussian nodes into the
c         coefficients of the corresponding legendre series.
c         no attempt has been made to make this code efficient,
c         but its speed is normally sufficient, and it is
c         mercifully short.
c 
c                 input parameters:
c 
c  itype - the type of the calculation to be performed
c          itype=0 means that only the gaussian nodes are
c                  to be constructed.
c          itype=1 means that only the nodes and the weights
c                  are to be constructed
c          itype=2 means that the nodes, the weights, and
c                  the matrices u, v are to be constructed
c  n - the number of gaussian nodes and weights to be generated
c 
c                 output parameters:
c 
c  x - the order n gaussian nodes - computed independently
c          of the value of itype.
c  u - the n*n matrix converting the  values at of a polynomial of order
c         n-1 at n legendre nodes into the coefficients of its
c         legendre expansion - computed only in itype=2
c  v - the n*n matrix converting the coefficients
c         of an n-term legendre expansion into its values at
c         n legendre nodes (note that v is the inverse of u)
c          - computed only in itype=2
c  whts - the corresponding quadrature weights - computed only
c         if itype .ge. 1
c 
c       . . . construct the nodes and the weights of the n-point gaussian
c             quadrature
c 
        ifwhts=0
        if(itype. gt. 0) ifwhts=1
        call legewhts(n,x,whts,ifwhts)
c 
c       construct the matrix of values of the legendre polynomials
c       at these nodes
c 
        if(itype .ne. 2) return
        do 1400 i=1,n
c 
        call legepols(x(i),n-1,u(1,i) )
 1400 continue
c 
        do 1800 i=1,n
        do 1600 j=1,n
        v(i,j)=u(j,i)
 1600 continue
 1800 continue
c 
c       now, v converts coefficients of a legendre expansion
c       into its values at the gaussian nodes. construct its
c       inverse u, converting the values of a function at
c       gaussian nodes into the coefficients of a legendre
c       expansion of that function
c 
        do 2800 i=1,n
        d=1
        d=d*(2*i-1)/2
        do 2600 j=1,n
        u(i,j)=v(j,i)*whts(j)*d
 2600 continue
 2800 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine legewhts(n,ts,whts,ifwhts)
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
        call legepol(xk,n,pol,der)
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
c 
c       now, use the explicit integral formulae
c       to obtain the weights
c 
        if(ifwhts .eq. 0) return
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
        subroutine legepol(x,n,pol,der)
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
c 
        subroutine legepols(x,n,pols)
        implicit real *8 (a-h,o-z)
        save
        dimension pols(1)
c 
        pkm1=1
        pk=x
c 
        pk=1
        pkp1=x
c 
c 
c        if n=0 or n=1 - exit
c 
        if(n .ge. 2) goto 1200
        pols(1)=1
        if(n .eq. 0) return
c 
        pols(2)=x
        return
 1200 continue
c 
        pols(1)=1
        pols(2)=x
c 
c       n is greater than 2. conduct recursion
c 
        do 2000 k=1,n-1
        pkm1=pk
        pk=pkp1
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
        pols(k+2)=pkp1
 2000 continue
c 
        return
        end
c 
c 
c 
c 
c 
      SUBROUTINE legeFDER(X,VAL,der,PEXP,N)
      IMPLICIT REAL *8 (A-H,O-Z)
        save
      REAL *8 PEXP(1)
C 
C     This subroutine computes the value and the derivative
c     of a gaussian expansion with coefficients PEXP
C     at point X in interval [-1,1]
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
C     der = computed value of the derivative
C 
C 
        done=1
        pjm2=1
        pjm1=x
        derjm2=0
        derjm1=1
c 
        val=pexp(1)*pjm2+pexp(2)*pjm1
        der=pexp(2)
c 
        dd=x**2-1
        dd=done/dd
        DO 600 J = 2,N
c 
        pj= ( (2*j-1)*x*pjm1-(j-1)*pjm2 ) / j
        val=val+pexp(j+1)*pj
c 
        derj=j*(x*pj-pjm1)*dd
        der=der+pexp(j+1)*derj
c 
        pjm2=pjm1
        pjm1=pj
        derjm2=derjm1
        derjm1=derj
 600   CONTINUE
c 
      RETURN
      END
  
