       implicit real *16(a-h,o-z)
       dimension x(1000),u(10000),v(10000),whts(10000),
     1     f(1000),coefs(1000),f2(1000)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
         PRINT *, 'ENTER n'
         READ *,n
         CALL PRINF('n=*',n,1 )
C 
         PRINT *, 'ENTER a'
         READ *,a
         CALL PRINQ('a=*',a,1 )
c 
c       test the extended precision qatan
c 
        call foolcoi(ier)
c 
        call foolcosi(a,extcos,extsin)
        tga=extsin/extcos
        aa=qatan(tga)
        call prinq('aa=*',aa,1)
        call prinq('and aa-a=*',aa-a,1)
c 
c       initialize the chebychev decomposition routine
c 
        itype=2
        call qchebexp(itype,n,x,u,v,whts)
        call prinq('after qchebexps, whts=*',whts,n)
c 
c       construct a test function
c 
        do 1200 i=1,n
        f(i)=x(i)**5
        f(i)=qexp(x(i))
 1200 continue
         call prinq('f as construcrted*',f,n)
c 
c       construct the chebychev expansion of the test function
c 
        call matvec2(u,f,n,coefs)
c 
        call prinq('coefs are*',coefs,n)
c 
c       reconstruct the function from its chebychev expansion and
c       check the error
c 
        call matvec2(v,coefs,n,f2)
c 
        do 1400 i=1,n
        f2(i)=f2(i)-f(i)
 1400 continue
        call prinq('and f2-f=*',f2,n)
c 
c        check it at the point a
c 
        call QCHFUNDE(a,VAL,der,coefs,N-1)
         call prinq('val(a)=*',val,1)
         val2=qexp(a)
         call prinq('and analytically, val2=*',val2,1)
         call prinq('and the difference*',val2-val,1)
         call prinq('and the derivative*',der-val2,1)
  
  
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine matvec2(a,x,n,y)
        implicit real *16(a-h,o-z)
        save
        dimension a(n,n),x(1),y(1)
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c               this is the end of the debugging code and the beginning
c               of the actual chebychev expansion routine
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        subroutine qchebexp(itype,n,x,u,v,whts)
        implicit real *16(a-h,o-z)
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
        pi=qatan(done)*4
        h=pi/(2*n)
        call foolcoi(ier)
ccc        do 1200 i=n,1,-1
        do 1200 i=1,n
        t=(2*i-1)*h
cccc        x(n-i+1)=qcos(t)
        call foolcosi(t,x(n-i+1),extsin)
1200  CONTINUE
ccc        call prinq('chebychev nodes as constructed*',x,n)
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
        d=done/qsqrt(d)
        do 2400 j=1,n
        u(i,j)=u(i,j)*d
 2400 continue
 3000 continue
c 
c        now, rescale the matrix
c 
        ddd=2
        ddd=qsqrt(ddd)
        dd=n
        dd=done/qsqrt(dd/2)
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
      SUBROUTINE QCHFUNDE(X,VAL,der,TEXP,N)
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
      IMPLICIT REAL *16(A-H,O-Z)
        save
      REAL *16TEXP(1)
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
        subroutine foolcosi(x,extcos,extsin)
        implicit real *16 (a-h,o-z)
        save
        dimension factin(0:60)
c 
c        compute, slowly but correctly, the sine and cosine
c        of x with extended precision
c 
        done=1
        two=2
        half=done/2
        n=20
        small=0.01
c 
        nhalf=0
        xx=x
ccccc         if(2 .ne. 3) goto 1500
        do 1200 i=1,100
        if(qabs(xx) .lt. small) goto 1400
        xx=xx*half
        nhalf=nhalf+1
 1200 continue
 1400 continue
 1500 continue
cccc         call prinf('nhalf=*',nhalf,1)
cccc         call prinq('xx=*',xx,1)
c 
c       calculate sine and cosine of xx
c 
        extcos=done
        extsin=0
c 
        d1=xx
        d2=xx**2
        dd=d2
        do 1600 i=1,n
        extcos=extcos+d2*factin(2*i)
        extsin=extsin+d1*factin(2*i-1)
        d1=d1*dd
        d2=d2*dd
 1600 continue
        extsin=-extsin
cccc          if(2 .ne. 3) return
c 
c       now, square the thing nhalf times
c 
        if(nhalf .eq. 0) return
        do 1800 i=1,nhalf
cccc          call prinf('i=*',i,1)
cccc          call prinq('extcos=*',extcos,1)
cccc          call prinq('extsin=*',extsin,1)
        d2=extcos**2-extsin**2
        d1=two*extcos*extsin
        extsin=d1
        extcos=d2
 1800 continue
        return
c 
c 
c 
c 
        entry foolcoi(ier)
c 
c       this is initialization entry point
c 
        coe=1
        factin(0)=1
        do 2200 i=1,55
        coe=-coe
        d=i
        factin(i)=factin(i-1)/d*coe
 2200 continue
         call prinq('factin as created*',factin(0),55)
        return
        end
