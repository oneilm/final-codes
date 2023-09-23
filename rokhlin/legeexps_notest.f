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
        subroutine legeinmt(n,ainte,adiff,x,whts,endinter,
     1      itype,w)
        implicit real *8 (a-h,o-z)
        save
        dimension ainte(1),w(1),x(1),whts(1),adiff(1),endinter(1)
c 
c 
c        for the user-specified n, this subroutine constructs
c        the matrices of spectral indefinite integration and/or
c        spectral differentiation on the n Gaussian nodes
c        on the interval [-1,1]. Actually, this is omnly a
c        memory management routine. All the actual work is done
c        by the subroutine legeinm0 (see)
c 
c                           input parameters:
c 
c  n - the number of Gaussian nodes on the interval [-1,1]
c  itype - the type of the calculation to be performed
c          EXPLANATION:
c       itype=1 means that only the matrix ainte will
c               be constructed
c       itype=2 means that only the matrix adiff will
c               be constructed
c       itype=3 means that both matrices ainte and adiff
c               will be constructed
c 
c                           output paramaters:
c 
c  ainte - the matrix of spectral indefinite integration on
c          the Gaussian nodes
c  adiff - the matrix of spectral differentiation on
c          the Gaussian nodes
c  x - the n Gaussian nodes on the intervl [-1,1]
c  whts - the n Gaussian weights on the interval [-1,1]
c  endinter - the interpolation coefficients converting the
c          values of a function at n Gaussian nodes into its
c          value at 1 (the right end of the interval)
c 
c                           work arrays:
c 
c  w - must be 3* n**2 + 2*n +50 *8 locations long
c 
c        . . . allocate memory for the construction of the integrating
c              matrix
c 
        ipolin=1
        lpolin=n+5
c 
        ipolout=ipolin+lpolin
        lpolout=n+5
c 
        iu=ipolout+lpolout
        lu=n**2+1
c 
        iv=iu+lu
        lv=n**2+1
c 
        iw=iv+lv
        lw=n**2+1
c 
        ltot=iw+lw
c 
c        construct the integrating matrix
c 
        call legeinm0(n,ainte,adiff,w(ipolin),w(ipolout),
     1      x,whts,w(iu),w(iv),w(iw),itype,endinter)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine legeinm0(n,ainte,adiff,polin,polout,
     1      x,whts,u,v,w,itype,endinter)
        implicit real *8 (a-h,o-z)
        save
        dimension ainte(n,n),u(n,n),v(n,n),w(n,n),
     1      endinter(1),x(n),whts(n),polin(n),polout(n),
     2      adiff(n,n)
c 
c        for the user-specified n, this subroutine constructs
c        the matrices of spectral indefinite integration and/or
c        spectral differentiation on the n Gaussian nodes
c        on the interval [-1,1]
c 
c                           input parameters:
c 
c  n - the number of Gaussian nodes on the interval [-1,1]
c  itype - the type of the calculation to be performed
c          EXPLANATION:
c       itype=1 means that only the matrix ainte will
c               be constructed
c       itype=2 means that only the matrix adiff will
c               be constructed
c       itype=3 means that both matrices ainte and adiff
c               will be constructed
c 
c                           output paramaters:
c 
c  ainte - the matrix of spectral indefinite integration on
c          the Gaussian nodes
c  adiff - the matrix of spectral differentiation on
c          the Gaussian nodes
c  x - the n Gaussian nodes on the intervl [-1,1]
c  whts - the n Gaussian weights on the interval [-1,1]
c 
c                           work arrays:
c 
c  polin, polout - must be n+3 real *8 locations each
c 
c  u, v, w - must be n**2+1 real *8 locations each
c 
c        . . . construct the matrices of the forward and inverse
c              Legendre transforms
c 
        itype2=2
        call legeexps(itype2,n,x,u,v,whts)
c 
cccc         call prin2('after legeexps, u=*',u,n*n)
c 
c        if the user so requested,
c        construct the matrix converting the coefficients of
c        the Legendre series of a function into the coefficients
c        of the indefinite integral of that function
c 
        if(itype. eq. 2) goto 2000
c 
        do 1600 i=1,n
c 
        do 1200 j=1,n+2
        polin(j)=0
 1200 continue
c 
        polin(i)=1
c 
        call legeinte(polin,n,polout)
c 
        do 1400 j=1,n
        ainte(j,i)=polout(j)
 1400 continue
c 
 1600 continue
c 
cccc         call prin2('ainte initially is*',ainte,n*n)
c 
c        multiply the three, obtaining the integrating matrix
c 
        call matmul(ainte,u,w,n)
        call matmul(v,w,ainte,n)
c 
 2000 continue
c 
c        if the user so requested,
c        construct the matrix converting the coefficients of
c        the Legendre series of a function into the coefficients
c        of the derivative of that function
c 
        if(itype. eq. 1) goto 3000
c 
        do 2600 i=1,n
c 
        do 2200 j=1,n
        polin(j)=0
 2200 continue
c 
        polin(i)=1
c 
        call legediff(polin,n,polout)
c 
        do 2400 j=1,n
        adiff(j,i)=polout(j)
cccc        ainte(i,j)=polout(j)
 2400 continue
c 
 2600 continue
c 
cccc         call prin2('adiff initially is*',adiff,n*n)
c 
c        multiply the three, obtaining the integrating matrix
c 
        call matmul(adiff,u,w,n)
        call matmul(v,w,adiff,n)
c 
 3000 continue
c 
c        construct the vector of interpolation coefficients
c        converting the values of a polynomial at the Gaussian
c        nodes into its value at the right end of the interval
c 
        do 3400 i=1,n
c 
        d=0
        do 3200 j=1,n
        d=d+u(j,i)
 3200 continue
        endinter(i)=d
 3400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine legeinte(polin,n,polout)
        implicit real *8 (a-h,o-z)
        save
        dimension polin(1),polout(1)
c 
c       this subroutine computes the indefinite integral of the
c       legendre expansion polin getting the expansion polout
c 
c 
c                       input parameters:
c 
c  polin - the legendre expansion to be integrated
c  n - the order of the expansion polin
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c         also nothe that the order of the integrated expansion is
c         n+1 (who could think!)
c 
c                       output parameters:
c 
c  polout - the legendre expansion of the integral of the function
c         represented by the expansion polin
c 
        do 1200 i=1,n+2
        polout(i)=0
 1200 continue
c 
        do 2000 k=2,n+1
        j=k-1
c 
        polout(k+1)=polin(k)/(2*j+1)+polout(k+1)
        polout(k-1)=-polin(k)/(2*j+1)+polout(k-1)
c 
 2000 continue
c 
        polout(2)=polin(1)+polout(2)
c 
        dd=0
        sss=-1
        do 2200 k=2,n+1
c 
        dd=dd+polout(k)*sss
        sss=-sss
 2200 continue
c 
        call prin2('dd=*',dd,1)
        polout(1)=-dd
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine legediff(polin,n,polout)
        implicit real *8 (a-h,o-z)
        save
        dimension polin(1),polout(1)
c 
c       this subroutine differentiates the legendre
c       expansion polin getting the expansion polout
c 
c 
c                       input parameters:
c 
c  polin - the legendre expansion to be differentiated
c  n - the order of the expansion polin
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c         also nothe that the order of the integrated expansion is
c         n+1 (who could think!)
c 
c                       output parameters:
c 
c  polout - the legendre expansion of the derivative of the function
c         represented by the expansion polin
c 
        do 1200 k=1,n+1
        polout(i)=0
 1200 continue
c 
        pk=polin(n+1)
        pkm1=polin(n)
        pkm2=0
        do 2000 k=n+1,2,-1
c 
        j=k-1
c 
        polout(k-1)=pk*(2*j-1)
        if(k .ge. 3) pkm2=polin(k-2)+pk
c 
        pk=pkm1
        pkm1=pkm2
c 
 2000 continue
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
C     at point X in interval [-1,1].
c 
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
        DO 600 J = 2,N
c 
        pj= ( (2*j-1)*x*pjm1-(j-1)*pjm2 ) / j
        val=val+pexp(j+1)*pj
c 
        derj=(2*j-1)*(pjm1+x*derjm1)-(j-1)*derjm2
  
ccc         call prin2('derj=*',derj,1)
  
  
        derj=derj/j
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
c 
c 
c 
c 
c 
      SUBROUTINE legeexev(X,VAL,PEXP,N)
      IMPLICIT REAL *8 (A-H,O-Z)
        save
      REAL *8 PEXP(1)
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
        der=pexp(2)
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
        subroutine matmul(a,b,c,n)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),b(n,n),c(n,n)
c 
        do 2000 i=1,n
        do 1800 j=1,n
        d=0
        do 1600 k=1,n
        d=d+a(i,k)*b(k,j)
 1600 continue
        c(i,j)=d
 1800 continue
 2000 continue
        return
c 
c 
c 
c 
        entry matmua(a,b,c,n)
ccc          call prin2('in matmua, a=*',a,n**2)
ccc          call prin2('in matmua, b=*',b,n**2)
        do 3000 i=1,n
        do 2800 j=1,n
        d=0
        do 2600 k=1,n
        d=d+a(i,k)*b(j,k)
 2600 continue
        c(i,j)=d
 2800 continue
 3000 continue
ccc          call prin2('exiting, c=*',c,n**2)
        return
        end
  
