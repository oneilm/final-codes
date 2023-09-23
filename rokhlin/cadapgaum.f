        implicit real *8 (a-h,o-z)
        external fun
        dimension work(200 000)
        complex *16 vals(100),traps(100),rints(100)
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
c 
        PRINT *, 'ENTER nn '
        READ *,nn
        CALL PRINF('nn=*',nn,1 )
c 
c        integrate
c 
  
        jtest=4
        a=0
        b=1
        eps=1.0d-15
c 
        call cadapgaum(ier,a,b,fun,nn,jtest,par2,m,eps,
     1      rints,maxrec,numint,work)
        call prin2('after cadapgaum, rints=*',rints,2*nn)
        call prinf('after cadapgaum, ier=*',ier,1)
        call prinf('after cadapgaum, maxrec=*',maxrec,1)
        call prinf('after cadapgaum, numint=*',numint,1)
  
  
ccc        if(2 .ne. 3) stop
c 
c        evaluate the integral via trapezoidal rule
c 
        do 2400 i=1,nn
        traps(i)=0
 2400 continue
c 
        trap=0
        n=100000
        h=(b-a)/(n+1)
        do 3000 i=1,n
        t=a+i*h
        call fun(t,nn,jtest,dummy,vals)
  
        do 2600 j=1,nn
        traps(j)=traps(j)+vals(j)
 2600 continue
c 
cccc        trap=trap+val
 3000 continue
c 
        do 3200 j=1,nn
        traps(j)=traps(j)*h
 3200 continue
  
cccc      trap=trap*h
  
        call prin2('and traps=*',traps,2*nn)
  
  
  
c 
        stop
        end
c 
c 
c 
c 
c 
c 
        subroutine fun(x,n,jtest,dummy,vals)
        implicit real *8 (a-h,o-z)
        save
        complex *16 vals(1),eye
c 
	data eye/(0.0d0,1.0d0)/
c 
         half=1
         half=half/2
         done=1
  
        do 1200 i=1,n
        vals(i)=x**(i-1)+eye*((1-x)**(i-1))
 1200 continue
  
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
        subroutine cadapgaum(ier,a,b,fun,n,par1,par2,m,eps,
     1      rints,maxrec,numint,work)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension t(100),w(100),stack(400),
     1      par1(1),par2(1),work(1)
	complex *16 rints(1)
c 
        data m7/-2341034/
c 
c       this subroutine uses the adaptive Gaussian quadrature
c       to evaluate the integral of the user-supplied vector-valued
c       function fun on the user-specified interval [a,b].
c       PLEASE NOTE THAT THIS IS A COMPLEX-VALUED VERSION OF THE
C       SUBROUTINE ADAPGAUM (SEE).
c 
c                       input parameters:
c 
c  a,b - the ends of the interval on which the integral is to
c       be evaluated
c  fun - the user-supplied subroutine evaluating the function
c        to be integrated. the calling sequence of fun must be
c 
c        fun(x,n,par1,par2,vals).                            (1)
c 
c        in (1), x is a point on the interval [a,b] where
c        the function is to be integrated, n is the dimensionality
c        of the vector-valued function to be integrated, and
c        par1, par2 are two parameters to be used by fun; they can
c        be  variables or arrays, real or integer, as desired.
c        fun is assumed to be real *8.
c  par1, par2 - parameters to be used by the user-supplied
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
c  rints - the vector of integrals as evaluated
c  maxrec - the maximum depth to which the recursion went at its
c         deepest point. can not be greater than 200, since at that
c         point ier is set to 8 and the execution of the subroutine
c         terminated.
c  numint - the total number of intervals in the subdivision. can not
c         be greater than 100000,  since at that
c         point ier is set to 16 and the execution of the subroutine
c         terminated.
c 
c                       work array:
c 
c  work - ideally, should be 202*n + 20 real *8 locations long.
c         However, the subroutine acutually uses only (numrec+2)*n + 20
c         elements. Thus, in most cases, the subroutine will use much
c         less space. However, caveat emptor!
c 
c 
c          . . . check if the gauss quadarture has been
c                initialized at a preceeding call to this routine
c 
        if(m .eq. m7) goto 1200
cccc        call gauswhmu(m,t,w)
c 
        ifwhts=1
        call legewhts(m,t,w,ifwhts)
c 
        m7=m
 1200 continue
c 
c        integrate the user-supplied function using the
c        adaptive gaussian quadratures
c 
        nnmax=100000
        maxdepth=200
c 
        ivalue2=1
        lvalue2=2*n+10
c 
        ivalue3=ivalue2+lvalue2
        lvalue3=2*n+10
c 
        ivaltmp=ivalue3+lvalue3
        lvaltmp=2*n+10
c 
        ivals=ivaltmp+lvaltmp
c 
        call adinrecm(ier,stack,a,b,fun,n,par1,par2,t,w,m,
     1      work(ivals),nnmax,eps,rints,maxdepth,maxrec,
     2      numint,work(ivalue2),work(ivalue3),work(ivaltmp) )
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine adinrecm(ier,stack,a,b,fun,n,
     1      par1,par2,t,w,m,vals,nnmax,eps,
     2      rints,maxdepth,maxrec,numint,value2,value3,valtmp)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension stack(2,1),t(1),w(1),par1(1),par2(1)
	complex *16 rints(1),vals(n,1),value2(1),value3(1),valtmp(1)
c 
c       start the recursion
c 
        stack(1,1)=a
        stack(2,1)=b
        call coneintmu(a,b,fun,n,par1,par2,t,w,m,vals(1,1),valtmp)
c 
c       recursively integrate the thing
c 
        j=1
        do 1200 i=1,n
        rints(i)=0
 1200 continue
c 
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
        call coneintmu(stack(1,j),c,fun,n,
     1      par1,par2,t,w,m,value2,valtmp)
c 
        call coneintmu(c,stack(2,j),fun,n,
     1      par1,par2,t,w,m,value3,valtmp)
c 
        dd=0
        do 1400 jj=1,n
        ddd=abs(value2(jj)+value3(jj)-vals(jj,j) )
        if(ddd .gt. dd) dd=ddd
 1400 continue
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
        do 1600 jj=1,n
        rints(jj)=rints(jj)+value2(jj)+value3(jj)
 1600 continue
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
        call carrmomom(value2,vals(1,j+1),n)
c 
        stack(1,j)=(stack(1,j)+stack(2,j))/2
        call carrmomom(value3,vals(1,j),n)
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
        subroutine carrmomom(x,y,n)
        implicit real *8 (a-h,o-z)
        save
        complex *16 x(1),y(1)
c 
        do 1200 i=1,n
        y(i)=x(i)
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine coneintmu(a,b,fun,n,par1,par2,t,w,m,rints,vals)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension t(1),w(1),par1(1),par2(1)
	complex *16 rints(1),vals(1)
c 
c       integrate the function fun on the interval [a,b]
c 
        do 1200 j=1,n
        rints(j)=0
 1200 continue
c 
        u=(b-a)/2
        v=(b+a)/2
c 
        do 1600 i=1,m
           tt=u*t(i)+v
           call fun(tt,n,par1,par2,vals)
c 
	   uw=u*w(i)
        do 1400 j=1,n
           rints(j)=rints(j)+vals(j)*uw
 1400 continue
c 
 1600 continue
c 
        return
        end
