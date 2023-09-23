        implicit real *8 (a-h,o-z)
        dimension t(1000),whts(1000),u(10000),t2(1000),
     1    c1(10000),c2(10000),f(100),g(100),tjs(100),
     2    v(10000),f2(100),valnew(100),polin(100),polout(100),
     3    polout2(100),ainte(10 000),adiff(10 000),w(40 000),
     4    aprod(40 000),endinter(100),amatrint(10 000),
     5      ts(100),xs(100),fts(100),fxs(100),fxs2(100)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER n '
        READ *,n
        CALL PRINF('n=*',n,1 )
c 
        PRINT *, 'ENTER m '
        READ *,m
        CALL PRINF('m=*',m,1 )
c 
        PRINT *, 'ENTER y'
        READ *,y
        CALL PRIN2('y=*',y,1 )
c 
c        test the weights of the chebychev quardature
c 
        itype=2
        call chebexps(itype,n,t,u,v,whts)
        call prin2('after chebexps, t=*',t,n)
        call prin2('after chebexps, whts=*',whts,n)
        call prin2('after chebexps, u=*',u,n**2)
c 
c 
        rint=0
        do 2200 i=1,n
        rint=rint+fun(t(i),m)*whts(i)
 2200 continue
        call prin2('and rint=*',rint,1)
c 
           done=1
c 
ccc        rint2=(done**(m+1)-(-done)**(m+1))/(m+1)
c 
ccc         rint2= ((-1)**(m+1)-done)/(m+1)-
ccc     1       ((-1)**(m-1)-done)/(m-1)
         rint2=-rint2/2
        call prin2('while rint2=*',rint2,1)
        call prin2('and the difference is*',rint-rint2,1)
c 
c        maltiply u by v its adjoint
c 
        call matmul_cheb(u,v,c1,n)
         call prin2('and u by v is*',c1,n**2)
c 
c        see what it does to a chebychev polynomial
c 
        do 3200 i=1,n
        f(i)=t(i)
        call tjcom(t(i),tjs(2),m+3)
        f(i)=tjs(m)
 3200 continue
        call matvec(u,f,g,n)
         call prin2('and after matvec, g=*',g,n)
c 
c       now, check the subroutine chfunder
c       . . . computation of the value of the expansion
c 
        do 3400 i=1,n
c 
        call CHFUNDER(t(i),f2(i),der,g,N)
 3400 continue
        call prin2('T_m at chebychev nodes is*',f,n)
        call prin2('and T_m computed by chfunder*',f2,n)
c 
c       . . . and of its derivative
c 
        h=0.001
ccc        h=0.01
        do 3600 i=1,n
        call CHFUNDER(t(i),fff,f(i),g,N)
        call CHFUNDER(t(i)+h,fp,fff,g,N)
        call CHFUNDER(t(i)-h,fm,fff,g,N)
        f2(i)=(fp-fm)/(2*h)
 3600 continue
         call prin2('derivative via expansion is*',f,n)
         call prin2('derivative via finite differences*',f2,n)
c 
c       check the evaluation of the expansion via the
c       new subroutine
c 
        do 3800 i=1,n
c 
        call chebexev(t(i),VALnew(i),g,N)
 3800 continue
        call prin2('value of the expansion via new subroutine*',
     1      valnew,n)
c 
c        test the subroutines for the integration and differentiation
c        of chebychev series
c 
        do 4000 i=1,100
        polin(i)=0
        polout(i)=-777
 4000 continue
        polin(m)=1
c 
        call chebinte(polin,n,polout)
c 
        call chebdiff(polout,n+1,polout2)
c 
        call prin2('testing chebdiff and chebinte, polin=*',polin,n+3)
        call prin2(' polout=*',polout,n+3)
        call prin2(' polout2=*',polout2,n+3)
c 
        do 4200 i=1,n
        polout2(i)=polout2(i)-polin(i)
 4200 continue
c 
        call prin2('and polout2-polin=*',polout2,n)
c 
c        evaluate the obtained integrated expansion at the point -1
c 
        dd=-1
c 
        call chebexev(dd,VAlm1,polout,N+1)
c 
        call prin2('and valm1=*',valm1,1)
c 
        dd=1
c 
        call chebexev(dd,VAl1,polout,N+1)
        call prin2('and val1=*',val1,1)
c 
c        test the routine for the evaluation of all
c        chebychev polynomials at a point
c 
        call chebpols(y,n,polin)
c 
        call prin2('chebychev polynomials via the subroutine*',
     1      polin,n+3)
c 
        do 4400 i=1,n+1
        polout(i)=dcos(dacos(y) * (i-1) )
 4400 continue
c 
        call prin2('and via direct formula*',polout,n+3)
c 
c        test the stand-alone subroutine for the evaluation of a
c        single Chebychev polynomial and its derivative
c 
        call chebpol(y,m,pol,der)
c 
        call prin2('testing chebpol, pol=*',pol,1)
        call prin2('testing chebpol, der=*',der,1)
c 
        hh=0.00001
        call chebpol(y-hh,m,polm,der)
        call chebpol(y+hh,m,polp,der)
c 
        call prin2('and derivative numerically,*',
     1      (polp-polm)/hh/2,1)
c 
c       test the subroutine generating the matrices
c       of spectral integration and differentiation on the
c       Chebychev nodes
c 
        itype=3
        call chebinmt(n,ainte,adiff,t,whts,endinter,
     1      itype,w)
c 
        call prin2('after chebinmt, whts=*',whts,n)
        call prin2('after chebinmt, t=*',t,n)
        call prin2('after chebinmt, endinter=*',endinter,n)
        call prin2('after chebinmt, ainte=*',ainte,n**2)
        call prin2('after chebinmt, adiff=*',adiff,n**2)
c 
c       apply the matrices ainte, adiff to a test vector and see what
c       happens
c 
        do 4600 i=1,n
        polin(i)=t(i)**(m-1)
cccc        polin(i)=1
 4600 continue
c 
c       . . . integrate
c 
       call prin2('polin before indefinite integration*',polin,n)
  
        call matvec(ainte,polin,polout,n)
  
       call prin2('polout after integration*',polout,n)
  
       do 4800 i=1,n-1
       f2(i)=(polout(i+1)-polout(i)) /(t(i+1)-t(i))
 4800 continue
        call prin2('numerical derivatives of polout*',f2,n-1)
c 
        call matvec(adiff,polout,polout2,n)
c 
       call prin2('polout2 after differentiation*',polout2,n)
c 
       do 5200 i=1,n
       polout2(i)=polout2(i)-polin(i)
 5200 continue
       call prin2('and polout2-polin=*',polout2,n)
c
c        test the subroutines generating interpolation 
c        vectors and matrices
c
        done=1
        hm=2*done/(m-1)
c
        do 5400 i=1,m
c
        xs(i)=-1+(i-1)*h
        fxs(i)=xs(i)**(m-1)
 5400 continue
c
        call prin2('xs=*',xs,m)        
c
        call chematrin(n,m,xs,amatrint,ts,w)

c
        do 5600 i=1,n
c
        fts(i)=ts(i)**(m-1)
 5600 continue
c
        call prin2('fts=*',fts,n)                

        call chematvec2(amatrint,m,n,fts,fxs2)

        call prin2('fxs=*',fxs,m)                
        call prin2('fxs2=*',fxs2,m)                



        stop
        end
c 
c 
c 
c 
c 
        subroutine fun2(t,f,y)
        implicit real *8 (a-h,o-z)
        save
        f=t/(1+y*t)
        return
        end
c 
c 
c 
c 
c 
        subroutine tjcom(x,tjs,n)
        implicit real *8 (a-h,o-z)
        save
        dimension tjs(1)
c 
        tjs(0)=1
        tjs(1)=x
        do 1200 i=2,n
        tjs(i)=2*x*tjs(i-1)-tjs(i-2)
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine matvec(a,x,y,n)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),x(1),y(1)
c 
        do 2000 i=1,n
        d=0
        do 1400 j=1,n
        d=d+a(i,j)*x(j)
 1400 continue
        y(i)=d
 2000 continue
        return
        end
  
c 
c 
c 
c 
c 
        function fun(x,m)
        implicit real *8 (a-h,o-z)
ccc        fun=x**m
        save
        fun=dcos(m*dacos(x))
        return
        end
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c               this is the end of the debugging code and the beginning
c               of the actual chebychev expansion routine
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        subroutine chebexps(itype,n,x,u,v,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),whts(1),u(n,n),v(n,n)
c 
c         this subroutine constructs the chebychev nodes
c         on the interval [-1,1], and the weights for the
c         corresponding order n quadrature. it also constructs
c         the matrix converting the coefficients
c         of a chebychev expansion into its values at the n
c         chebychev nodes. no attempt has been made to
c         make this code efficient, but its speed is normally
c         sufficient, and it is mercifully short.
c 
c                 input parameters:
c 
c  itype - the type of the calculation to be performed
c          itype=0 means that only the chebychev nodes are
c                  to be constructed.
c          itype=1 means that only the nodes and the weights
c                  are to be constructed
c          itype=2 means that the nodes, the weights, and
c                  the matrices u, v are to be constructed
c  n - the number of chebychev nodes and weights to be generated
c 
c                 output parameters:
c 
c  x - the order n chebychev nodes - computed independently
c          of the value of itype.
c  u - the n*n matrix converting the  values at of a polynomial of order
c         n-1 at n chebychev nodes into the coefficients of its
c         chebychev expansion - computed only in itype=2
c  v - the n*n matrix converting the coefficients
c         of an n-term chebychev expansion into its values at
c         n chebychev nodes (note that v is the inverse of u)
c          - computed only in itype=2
c  whts - the corresponding quadrature weights - computed only
c         if itype .ge. 1
c 
c       . . . construct the chebychev nodes on the interval [-1,1]
c 
        ZERO=0
        DONE=1
        pi=datan(done)*4
        h=pi/(2*n)
ccc        do 1200 i=n,1,-1
        do 1200 i=1,n
        t=(2*i-1)*h
        x(n-i+1)=dcos(t)
1200  CONTINUE
        call prin2('chebychev nodes as constructed*',x,n)
c 
        if(itype .eq. 0) return
c 
c        construct the weights of the quadrature
c        formula based on the chebychev nodes,
c        and also the matrix of the chebychev transform
c 
c        . . . construct the first two rows of the matrix
c 
         if(itype .le. 1) goto 1350
        do 1300 i=1,n
        u(1,i)=1
        u(2,i)=x(i)
 1300 continue
 1350 continue
c 
c       construct all quadrature weights and the rest of the rows
c 
        do 2000 i=1,n
c 
c       construct the weight for the i-th node
c 
        Tjm2=1
        Tjm1=x(i)
        whts(i)=2
c 
        ic=-1
        do 1400 j=2,n-1
c 
c       calculate the T_j(x(i))
c 
        Tj=2*x(i)*Tjm1-Tjm2
c 
        if(itype .eq. 2) u(j+1,i)=tj
c 
        tjm2=tjm1
        tjm1=tj
c 
c       calculate the contribution of this power to the
c       weight
c 
  
        ic=-ic
        if(ic .lt. 0) goto 1400
        rint=-2*(done/(j+1)-done/(j-1))
        whts(i)=whts(i)-rint*tj
ccc        whts(i)=whts(i)+rint*tj
 1400 continue
        whts(i)=whts(i)/n
 2000 continue
           if(itype .ne. 2) return
c 
c        now, normalize the matrix of the chebychev transform
c 
        do 3000 i=1,n
c 
        d=0
        do 2200 j=1,n
        d=d+u(i,j)**2
 2200 continue
        d=done/dsqrt(d)
        do 2400 j=1,n
        u(i,j)=u(i,j)*d
 2400 continue
 3000 continue
c 
c        now, rescale the matrix
c 
        ddd=2
        ddd=dsqrt(ddd)
        dd=n
        dd=done/dsqrt(dd/2)
        do 3400 i=1,n
        do 3200 j=1,n
        u(j,i)=u(j,i)*dd
 3200 continue
        u(1,i)=u(1,i)/ddd
 3400 continue
c 
c        finally, construct the matrix v, converting the values at the
c        chebychev nodes into the coefficients of the chebychev
c        expansion
c 
        dd=n
        dd=dd/2
        do 4000 i=1,n
        do 3800 j=1,n
        v(j,i)=u(i,j)*dd
 3800 continue
 4000 continue
c 
        do 4200 i=1,n
        v(i,1)=v(i,1)*2
 4200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine chebinmt(n,ainte,adiff,x,whts,endinter,
     1      itype,w)
        implicit real *8 (a-h,o-z)
        save
        dimension ainte(1),w(1),x(1),whts(1),adiff(1),endinter(1)
c 
c 
c        for the user-specified n, this subroutine constructs
c        the matrices of spectral indefinite integration and/or
c        spectral differentiation on the n Chebychev nodes
c        on the interval [-1,1]. Actually, this is only a
c        memory management routine. All the actual work is done
c        by the subroutine legeinm0 (see)
c 
c                           input parameters:
c 
c  n - the number of Chebychev nodes on the interval [-1,1]
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
c          the Chebychev nodes
c  adiff - the matrix of spectral differentiation on
c          the Chebychev nodes
c  x - the n Chebychev nodes on the intervl [-1,1]
c  whts - the n Chebychev weights on the interval [-1,1]
c  endinter - the interpolation coefficients converting the
c          values of a function at n Chebychev nodes into its
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
        call chebinm0(n,ainte,adiff,w(ipolin),w(ipolout),
     1      x,whts,w(iu),w(iv),w(iw),itype,endinter)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine chebinm0(n,ainte,adiff,polin,polout,
     1      x,whts,u,v,w,itype,endinter)
        implicit real *8 (a-h,o-z)
        save
        dimension ainte(n,n),u(n,n),v(n,n),w(n,n),
     1      endinter(1),x(n),whts(n),polin(n),polout(n),
     2      adiff(n,n)
c 
c        for the user-specified n, this subroutine constructs
c        the matrices of spectral indefinite integration and/or
c        spectral differentiation on the n Chebychev nodes
c        on the interval [-1,1]
c 
c                           input parameters:
c 
c  n - the number of Chebychev nodes on the interval [-1,1]
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
c          the Chebychev nodes
c  adiff - the matrix of spectral differentiation on
c          the Chebychev nodes
c  x - the n Chebychev nodes on the intervl [-1,1]
c  whts - the n Chebychev weights on the interval [-1,1]
c 
c                           work arrays:
c 
c  polin, polout - must be n+3 real *8 locations each
c 
c  u, v, w - must be n**2+1 real *8 locations each
c 
c        . . . construct the matrices of the forward and inverse
c              Chebychev transforms
c 
        itype2=2
        call chebexps(itype2,n,x,u,v,whts)
c 
cccc         call prin2('after chebexps, u=*',u,n*n)
c 
c        if the user so requested,
c        construct the matrix converting the coefficients of
c        the Chebychev series of a function into the coefficients
c        of the indefinite integral of that function
c 
        if(itype. eq. 2) goto 2000
c 
        do 1600 i=1,n
c 
        do 1200 j=1,n
        polin(j)=0
 1200 continue
c 
        polin(i)=1
        call chebinte(polin,n+1,polout)
c 
        do 1400 j=1,n
        ainte(j,i)=polout(j)
 1400 continue
c 
 1600 continue
c 
c        multiply the three, obtaining the integrating matrix
c 
        call matmul_cheb(ainte,u,w,n)
        call matmul_cheb(v,w,ainte,n)
c 
 2000 continue
c 
c        if the user so requested,
c        construct the matrix converting the coefficients of
c        the Chebychev series of a function into the coefficients
c        of the derivative of that function
c 
        if(itype. eq. 1) goto 3000
c 
        do 2600 i=1,n
c 
        do 2200 j=1,n+3
        polin(j)=0
 2200 continue
c 
        polin(i)=1
        call chebdiff(polin,n+1,polout)
c 
        do 2400 j=1,n
        adiff(j,i)=polout(j)
 2400 continue
c 
 2600 continue
c 
cccc         call prin2('adiff initially is*',adiff,n*n)
c 
c        multiply the three, obtaining the differentiating matrix
c 
        call matmul_cheb(adiff,u,w,n)
        call matmul_cheb(v,w,adiff,n)
c 
 3000 continue
c 
c        construct the vector of interpolation coefficients
c        converting the values of a polynomial at the Chebychev
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
        subroutine chebpol(x,n,pol,der)
        implicit real *8 (a-h,o-z)
c 
        save
        d=dacos(x)
        pol=dcos(n*d)
        der=dsin(n*d)*n/dsqrt(1-x**2)
        return
        end
c 
c 
c 
c 
c 
        subroutine chebpols(x,n,pols)
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
        pkp1=2*x*pk-pkm1
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
        subroutine chebinte(polin,n,polout)
        implicit real *8 (a-h,o-z)
        save
        dimension polin(1),polout(1)
c 
c       this subroutine computes the indefinite integral of the
c       Chebychev expansion polin getting the expansion polout
c 
c 
c                       input parameters:
c 
c  polin - the Chebychev expansion to be integrated
c  n - the order of the expansion polin
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c         also nothe that the order of the integrated expansion is
c         n+1 (who could think!)
c 
c                       output parameters:
c 
c  polout - the Chebychev expansion of the integral of the function
c         represented by the expansion polin
c 
        do 1200 i=1,n+2
        polout(i)=0
 1200 continue
c 
        polout(2)=polin(1)
        if(n .eq. 0) return
        polout(3)=polin(2)/4
        if(n .eq. 1) return
c 
        do 2000 k=3,n
c 
        polout(k+1)=polin(k)/(2*k)
        polout(k-1)=-polin(k)/(2*k-4)+polout(k-1)
c 
 2000 continue
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
        subroutine chebdiff(polin,n,polout)
        implicit real *8 (a-h,o-z)
        save
        dimension polin(1),polout(1)
c 
c       this subroutine differentiates the Chebychev
c       expansion polin getting the expansion polout
c 
c 
c                       input parameters:
c 
c  polin - the Chebychev expansion to be differentiated
c  n - the order of the expansion polin
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c         also nothe that the order of the integrated expansion is
c         n+1 (who could think!)
c 
c                       output parameters:
c 
c  polout - the Chebychev expansion of the derivative of the function
c         represented by the expansion polin
c 
        do 1200 i=1,n+1
        polout(i)=0
 1200 continue
c 
        do 2000 k=1,n-1
c 
        polout(n-k)=polout(n-k+2)+(n-k)*polin(n-k+1) *2
 2000 continue
        polout(1)=polout(1)/2
c 
        return
        end
c 
c 
c 
c 
c 
      SUBROUTINE chebexev(X,VAL,TEXP,N)
      IMPLICIT REAL *8 (A-H,O-Z)
        save
      REAL *8 TEXP(1)
C 
C     This subroutine computes the value o a Chebychev
c     expansion with coefficients TEXP at point X in interval [-1,1]
C 
c                input parameters:
c 
C     X = evaluation point
C     TEXP = expansion coefficients
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
        val=texp(1)*pjm2+texp(2)*pjm1
c 
        DO 600 J = 2,N
c 
        pj= 2*x*pjm1-pjm2
        val=val+texp(j+1)*pj
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
      SUBROUTINE CHFUNDER(X,VAL,der,TEXP,N)
C 
C     This subroutine computes the value and the derivative
c     of a chebychev expansion with coefficients TEXP
C     at point X in interval [-1,1]
C 
c                input parameters:
c 
C     X = evaluation point
C     TEXP = expansion coefficients
C     N  = order of expansion
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
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
      SUBROUTINE CHebcval(X,VAL,TEXP,N)
C 
C     This subroutine computes the value and the derivative
c     of a chebychev expansion with coefficients TEXP
C     at point X in interval [-1,1]
C 
c                input parameters:
c 
C     X = evaluation point
C     TEXP = expansion coefficients
C     N  = order of expansion
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c 
c                output parameters:
c 
C     VAL = computed value
C     der = computed value of the derivative
C 
      IMPLICIT REAL *8 (A-H,O-Z)
        save
      complex *16 TEXP(1),val
C 
        done=1
        tjm2=1
        tjm1=x
c 
        val=texp(1)*tjm2+texp(2)*tjm1
c 
      DO 600 J = 2,N
c 
        tj=2*x*tjm1-tjm2
        val=val+texp(j+1)*tj
c 
        tjm2=tjm1
        tjm1=tj
 600   CONTINUE
c 
      RETURN
      END
c 
c 
c 
c 
c 
      SUBROUTINE CHebval(X,VAL,TEXP,N)
C 
C     This subroutine computes the value and the derivative
c     of a chebychev expansion with coefficients TEXP
C     at point X in interval [-1,1]
C 
c                input parameters:
c 
C     X = evaluation point
C     TEXP = expansion coefficients
C     N  = order of expansion
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c 
c                output parameters:
c 
C     VAL = computed value
C     der = computed value of the derivative
C 
      IMPLICIT REAL *8 (A-H,O-Z)
        save
      real *8 TEXP(1),val
C 
        done=1
        tjm2=1
        tjm1=x
c 
        val=texp(1)*tjm2+texp(2)*tjm1
c 
      DO 600 J = 2,N
c 
        tj=2*x*tjm1-tjm2
        val=val+texp(j+1)*tj
c 
        tjm2=tjm1
        tjm1=tj
 600   CONTINUE
c 
      RETURN
      END
c 
c 
c 
c 
c 
      SUBROUTINE CHFUNcDER(X,VAL,der,TEXP,N)
C 
C     This subroutine computes the value and the derivative
c     of a chebychev expansion with coefficients TEXP
C     at point X in interval [-1,1]
C 
c                input parameters:
c 
C     X = evaluation point
C     TEXP = expansion coefficients
C     N  = order of expansion
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c 
c                output parameters:
c 
C     VAL = computed value
C     der = computed value of the derivative
C 
      IMPLICIT REAL *8 (A-H,O-Z)
        save
      complex *16 TEXP(1),val,der
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
        subroutine matmul_cheb(a,b,c,n)
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
        end
c 
c 
c 
c 
c 
        subroutine chematrin(n,m,xs,amatrint,ts,w)
        implicit real *8 (a-h,o-z)
        save
        dimension amatrint(m,n),xs(1),w(1),ts(1)
c 
c 
c        This subroutine constructs the matrix interpolating
c        functions from the n-point Chebychev grid on the interval [-1,1]
c        to an arbitrary m-point grid (the nodes of the latter are
c        user-provided)
c 
c                 Input parameters:
c 
c  n - the number of interpolation nodes
c  m - the number of nodes to which the functions will be interpolated
c  xs - the points at which the function is to be interpolated
c 
c                  Output parameters:
c 
c  amatrint - the m \times n matrix conerting the values of a function
c        at the n Chebychev nodes into its values at m user-specified
c        (arbitrary) nodes
c  ts - the n Chebychev nodes on the interval [-1,1]
c 
c                  Work arrays:
c 
c  w - must be at least 2*n**2+n + 100 real *8 locations long
c 
  
        icoefs=1
        lcoefs=n+2
c 
        iu=icoefs+lcoefs
        lu=n**2+10
c 
        iv=iu+lu
c 
        ifinit=1
        do 2000 i=1,m
c 
        call chevecin(n,xs(i),ts,w(iu),w(iv),w(icoefs),ifinit)
c 
        do 1400 j=1,n
        amatrint(i,j)=w(j)
 1400 continue
c 
        ifinit=0
 2000 continue
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine chevecin(n,x,ts,u,v,coefs,ifinit)
        implicit real *8 (a-h,o-z)
        save
        dimension u(n,n),v(n,n),ts(1),coefs(1)
c 
c        This subroutine constructs the coefficients of the
c        standard interpolation formula connecting the values of a
c        function at n Chebychev nodes on the interval [a,b] with
c        its value at the point x \in R^1
c 
c                 Input parameters:
c 
c  n - the number of interpolation nodes
c  x - the points at which the function is to be interpolated
c  ts - the n Chebychev nodes on the interval [-1,1]; please note that
c        it is an input parameter only if the parameter ifinit (see
c        below) has been set to 1; otherwise, it is an output parameter
c  u - the n*n matrix converting the  values at of a polynomial of order
c         n-1 at n Chebychev nodes into the coefficients of its
c        Chebychev expansion; please note that
c        it is an input parameter only if the parameter ifinit (see
c        below) has been set to 1; otherwise, it is an output parameter
c  ifinit - an integer parameter telling the subroutine whether it should
c        initialize the Chebychev expander;
c     ifinit=1 will cause the subroutine to perform the initialization
c     ifinit=0 will cause the subroutine to  skip the initialization
c 
c                  Output parameters:
c 
c  coefs - the interpolation coefficients
c 
c                 Work arrays:
c 
c  v - must be at least n*n real *8 locations long
c 
c       . . . construct the n Chebychev nodes on the interval [-1,1];
c             also the corresponding Chebychev expansion-evaluation
c             matrices
c 
        itype=2
        if(ifinit .ne.0) call chebexps(itype,n,ts,u,v,coefs)
c 
c       evaluate the n Chebychev polynomials at the point where the
c       functions will have to be interpolated
c 
        call chebpols(x,n+1,v)
c 
c       apply the interpolation matrix to the ector of values
c       of polynomials from the right
c 
        call chematvec(u,v,coefs,n)
        return
        end
c 
c 
c 
c 
c 
        subroutine chematvec(a,x,y,n)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),x(n),y(n)
c 
        do 1400 i=1,n
        d=0
        do 1200 j=1,n
        d=d+a(j,i)*x(j)
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
        subroutine chematvec2(a,n,m,x,y)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(n,m),x(m),y(n),cd
c 
c        apply the matrix a to the vector x obtaining y
c 
        do 1400 i=1,n
        cd=0
        do 1200 j=1,m
        cd=cd+a(i,j)*x(j)
 1200 continue
        y(i)=cd
 1400 continue
c
        return
        end
