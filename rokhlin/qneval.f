        implicit real *8 (a-h,o-z)
        dimension qfuns(1000),qfuns3(1000)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER n'
        READ *,n
        CALL PRInf('n=*',n,1 )
C 
        PRINT *, 'ENTER x'
        READ *,x
        CALL PRIn2('x=*',x,1 )
c 
        call prin2('and x-1=*',x-1,1)
c 
c        evaluate Q_n (x) via qlegfuns
c 
        call qlegfuns(x,n,qfuns(2))
c 
        call prin2('qfuns after qlegfuns is*',qfuns,n+1)
c 
c        now, evaluate the functions via the qneval
c 
        call qneval(x,n,qfuns3(2))
c 
        call prin2('and qfuns3 via qneval*',qfuns3,n+1)
c 
        stop
        end
c 
c 
c 
c 
c 
        subroutine qlegfuns(x,n,qfuns)
        implicit real *8 (a-h,o-z)
        save
        dimension qfuns(1)
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c 
c        construct Q_0(x),Q_1(x)
c 
        done=1
        i0=0
        delta=1.0d-20
c 
        d=log( (done+x+ima*delta)/(done-x+ima*delta) )
c 
        qfuns(i0)=d/2
        qfuns(1)=x/2*d-1
c 
c       recurse up
c 
        do 1200 i=1,n-1
c 
        qfuns(i+1)=( (2*i+1)*x*qfuns(i)-i*qfuns(i-1) ) /(i+1)
 1200 continue
c 
        return
        end
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c      This is the end of the testing code and the beginning of the
c      code for the evaluation of the Legendre functions Q_i
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine qneval(x,n,qfuns)
        implicit real *8 (a-h,o-z)
        save
        dimension qfuns(1)
c 
c        for user-specified real x and integer n, this subroutine
c        evaluates the Legendre functions Q_0, Q_1, . . . ,Q_n at
c        the point x anywhere on the real line (except |x|=1).
c 
c                        Input parameters:
c 
c  x - the point where the functions Q_i will be evaluated
c  n - maximum order of Q to be evaluated
c 
c                        Output parameters:
c 
c  qfuns - the values of the Legendre functions at the
c      point x (n+1 of them things)
c 
c           IMPORTANT NOTE:
c 
c      ARRAY QFUNS MUST HAVE ELEMENT NUMBER 0 !!!!!!
c 
c        . . . if x is inside the interval [-1,1], or sufficiently
c              close to 1 or -1, evaluate the Legenedre functions
c              via the simple recursion
c 
        if(dabs(x) .lt. 1) goto 1100
c 
        delta=dabs(x)-1
        b=1+delta+sqrt(2*delta+delta**2)
c 
        if( (n+1)*log(b) .gt. 2.3) goto 2000
c 
 1100 continue
c 
c        construct Q_0(x),Q_1(x)
c 
        done=1
        i0=0
        d=log( abs((done+x)/(done-x)) )
c 
        qfuns(i0)=d/2
        qfuns(1)=x/2*d-1
c 
c       recurse up
c 
        do 1200 i=1,n-1
c 
        qfuns(i+1)=( (2*i+1)*x*qfuns(i)-i*qfuns(i-1) ) /(i+1)
 1200 continue
c 
        return
c 
 2000 continue
c 
c       the point x is outside the interval [-1,1].
c       determine how far we will have to start to get it by
c       the combination of recursing up and scaling
c 
        eps=1.0d-20
c 
        delta=dabs(x)-1
        b=1+delta+sqrt(2*delta+delta**2)
        nn=-log(eps)/log(b)+1
c 
        call prinf('nn as obtained in qninte is*',nn,1)
c 
c       use recursion down coupled with scaling to get the
c       Legendre functions Q_i
c 
        done=1
        i0=0
        d=log( abs((done+x)/(done-x)) )
c 
        q0=d/2
c 
        do 2050 i=i0,n
        qfuns(i)=0
 2050 continue
c 
c        if nn is greater than n+1, do the preliminary recursion
c        to get the starting values for qfuns(n+1),qfuns(n)
c 
        fip1=0
        fi  =1
c 
        if(nn .le. n+1) goto 2150
c 
        do 2100 i=nn,n+1,-1
c 
         fim1=( (2*i+1)*x*fi-(i+1)*fip1 ) / i
c 
        fip1=fi
        fi=fim1
 2100 continue
c 
        nn=n
 2150 continue
c 
        qfuns(nn+1)=fip1
        qfuns(nn)=fi
c 
c       recurse down starting with i=n
c 
        do 2200 i=nn,1,-1
c 
         qfuns(i-1)=( (2*i+1)*x*qfuns(i)-(i+1)*qfuns(i+1) ) / i
 2200 continue
c 
c        scale the values of the recursed (unscaled) Q_i to
c        obtain the correct one
c 
        rat=q0/qfuns(i0)
c 
        call prin2('and rat=*',rat,1)
c 
        do 2400 i=i0,n+1
c 
        qfuns(i)=rat*qfuns(i)
 2400 continue
 3000 continue
c 
        return
        end
  
  
