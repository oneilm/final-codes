        implicit real *8 (a-h,o-z)
      dimension t(1000),whts(1000)
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
c        test the weights of the chebychev quardature
c 
        call chebwhts(n,t,whts)
        call prin2('after chewhts, t=*',t,n)
        call prin2('after chewhts, whts=*',whts,n)
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
         rint2= ((-1)**(m+1)-done)/(m+1)-
     1       ((-1)**(m-1)-done)/(m-1)
         rint2=-rint2/2
        call prin2('while rint2=*',rint2,1)
        call prin2('and the difference is*',rint-rint2,1)
c 
        stop
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
c 
c 
c 
        subroutine chebwhts(n,x,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),whts(1)
c 
c         this subroutine constructs the chebychev nodes
c         on the interval [-1,1], and the weights for the
c         corresponding order n quadrature. the code is
c         a slow one (requiering order n**2 operations),
c         but its speed is normally sufficient, and it is
c         mercifully short.
c 
c                 input parameters:
c 
c  n- the number of chebychev nodes and weights to be generated
c 
c                 output parameters:
c 
c  x - the order n chebychev nodes
c  whts - the corresponding quadrature weights
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
c        construct the weights of the quadrature
c        formula based on the chebychev nodes
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
        do 1400 j=2,n
c 
c       calculate the T_j(x(i))
c 
        Tj=2*x(i)*Tjm1-Tjm2
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
        return
        end
