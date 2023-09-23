        implicit real *8 (a-h,o-z)
        external fun
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
c 
        PRINT *, 'ENTER m '
        READ *,m
        CALL PRINF('m=*',m,1 )
c 
c        integrate
c 
        jtest=4
        a=0
        b=1
        eps=1.0d-15
c 
        call adapcheb(ier,a,b,fun,jtest,par2,m,eps,
     1      rint,maxrec,numint)
        call prin2('after adapcheb, rint=*',rint,1)
        call prinf('after adapcheb, ier=*',ier,1)
        call prinf('after adapcheb, maxrec=*',maxrec,1)
        call prinf('after adapcheb, numint=*',numint,1)
c 
c        evaluate the integral via trapezoidal rule
c 
        trap=0
        n=10000
        h=(b-a)/(n+1)
        do 2000 i=1,n
        t=a+i*h
        trap=trap+fun(t,jtest,dummy)
 2000 continue
        trap=trap*h
        call prin2('and trap=*',trap,1)
  
  
  
  
c 
        stop
        end
c 
c 
c 
c 
c 
c 
        function fun(x,jtest,dummy)
        implicit real *8 (a-h,o-z)
c 
        save
        fun=x**jtest
         half=1
         half=half/2
         done=1
         fun=dlog(x)+dlog(1-x)+dlog(dabs(x-half))
cccc         fun=dlog(dabs(x))+dlog(dabs(done-x))
cccc          fun=half/x
ccc         call prin2('inside fun, fun=*',fun,1)
        return
        end
c 
c 
c 
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        this is the end of the debugging code and the
c        start of the actual quadrature routines.
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c 
c 
        subroutine adapcheb(ier,a,b,fun,par1,par2,m,eps,
     1      rint,maxrec,numint)
        implicit real *8 (a-h,o-z)
        save
        dimension t(100),w(100),stack(400),vals(200),
     1      par1(1),par2(1)
        data m7/-2341034/
c 
c       this subroutine uses the adaptive chebychev quadrature
c       to evaluate the integral of the user-supplied function
c       fun on the user-specified interval [a,b]
c 
c                       input parameters:
c 
c  a,b - the ends of the interval on which the integral is to
c       be evaluated
c  fun - the user-supplied function to be integrated. the calling
c       sequence of fun must be
c 
c        fun(x,par1,par2).                            (1)
c 
c        in (1), x is a point on the interval [a,b] where
c        the function is to be evaluated, and par1, par2
c        are two parameters to be used by fun; they can be
c        variables or arrays, real or integer, as desired.
c        fun is assumed to be real *8.
c  par1, par2 - partameters to be used by the user-supplied
c       function fun (see above)
c  m - the order of the quadrature to me used on each subinterval
c  eps - the accuracy (absolute) to which the integral will be
c       evaluated
c 
c                       output parameters:
c 
c  ier - error return code.
c          ier=0 means normal conclusion
c          ier=8 means that at some point, one subinterval in the
c                subdivision was smaller than (b-a)/2**200. this
c                is a fatal error.
c          ier=16 means that the total number of subintervals in the
c                adaptive subdivision of [a,b] turned out to be greater
c                than 100000.  this is a fatal error.
c 
c  rint - the integral as evaluated
c  maxrec - the maximum depth to which the recursion went at its
c         deepest point. can not be greater than 200, since at that
c         point ier is set to 8 and the execution of the subroutine
c         terminated.
c  numint - the totla number of intervals in the subdivision. can not
c         be greater than 100000,  since at that
c         point ier is set to 16 and the execution of the subroutine
c         terminated.
c 
c 
c          . . . check if the chebychev quadarture has been
c                initialized at a preceeding call to this routine
c 
        if(m .eq. m7) goto 1200
        call chebwhts(m,t,w)
        m7=m
 1200 continue
c 
c        integrate the user-supplied function using the
c        cbebychev quadratures
c 
        nnmax=100000
        maxdepth=200
c 
        call adinrec(ier,stack,a,b,fun,par1,par2,t,w,m,
     1      vals,nnmax,eps,rint,maxdepth,maxrec,numint)
        return
        end
c 
c 
c 
c 
c 
        subroutine adinrec(ier,stack,a,b,fun,
     1      par1,par2,t,w,m,vals,nnmax,eps,
     2      rint,maxdepth,maxrec,numint)
        implicit real *8 (a-h,o-z)
        save
        dimension stack(2,1),t(1),w(1),vals(1),par1(1),par2(1)
c 
c       start the recursion
c 
        stack(1,1)=a
        stack(2,1)=b
        call oneint(a,b,fun,par1,par2,t,w,m,vals(1))
c 
c       recursively integrate the thing
c 
        j=1
        rint=0
        ier=0
        maxrec=0
        do 3000 i=1,nnmax
ccc        call prinf('i=*',i,1)
        numint=i
        if(j .gt. maxrec) maxrec=j
ccc        call prinf('j=*',j,1)
c 
c       subdivide the current subinterval
c 
         c=(stack(1,j)+stack(2,j))/2
        call oneint(stack(1,j),c,fun,
     1      par1,par2,t,w,m,value2)
c 
        call oneint(c,stack(2,j),fun,
     1      par1,par2,t,w,m,value3)
c 
        dd=dabs(value2+value3-vals(j))
cccc         call prin2('in adinrec, dd=*',dd,1)
        ifdone=0
        if(dd .le. eps) ifdone=1
c 
c       if the function on this subinterval has been
c       integrated with sufficient accuracy - add the
c       value to that of the global integral and move up
c       in the stack
c 
        if(ifdone  .eq. 0) goto 2000
c 
        rint=rint+value2+value3
        j=j-1
c 
c        if the whole thing has been integrated - return
c 
        if(j .eq. 0) return
        goto 3000
 2000 continue
c 
c       if the function on this subinterval has not been
c       integrated with sufficient accuracy - move
c       down the stack
c 
        stack(1,j+1)=stack(1,j)
        stack(2,j+1)=(stack(1,j)+stack(2,j))/2
        vals(j+1)=value2
c 
        stack(1,j)=(stack(1,j)+stack(2,j))/2
        vals(j)=value3
c 
        j=j+1
c 
c       if the depth of the recursion has become excessive - bomb
c 
        if(j .le. maxdepth) goto 3000
        ier=8
        return
 3000 continue
        ier=16
        return
        end
c 
c 
c 
c 
c 
        subroutine oneint(a,b,fun,par1,par2,t,w,m,rint)
        implicit real *8 (a-h,o-z)
        save
        dimension t(1),w(1),par1(1),par2(1)
c 
c       integrate the function fun on the interval [a,b]
c 
        rint=0
        u=(b-a)/2
        v=(b+a)/2
        do 1200 i=1,m
        tt=u*t(i)+v
        rint=rint+fun(tt,par1,par2)*w(i)
 1200 continue
        rint=rint*u
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
cccc        call prin2('chebychev nodes as constructed*',x,n)
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
 1400 continue
        whts(i)=whts(i)/n
 2000 continue
        return
        end
