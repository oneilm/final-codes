        implicit real *8 (a-h,o-z)
        dimension funjs(10000)
        external fun
c 
        CALL PRINI(6,13)
C 
c       . . . construct the testing parameters
c 
c 
         PRINT *, 'ENTER x'
         READ *,x
         CALL PRIN2('x=*',x,1)
c 
c       calculate the integral sine of x
c 
c 
c         . . . via the subroutine
c 
  
         done=1
         pi=atan(done)*4
         rlim=10
         a=pi
  
         a=5
  
  
         a=x
  
         call fourlog(rlim,a,rint)
  
         call prin2('after fourlog, rint=*',rint,1)
  
  
        aa=-rlim
        bb=rlim
        eps=1.0d-14
        m=20
  
        call adapgaus(ier,aa,bb,fun,a,par2,m,eps,
     1      rint2,maxrec,numint)
  
        call prin2('after adapgaus, rint2=*',rint2,1)
        call prin2('and rint2-rint=*',rint2-rint,1)
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine fourlog(rlim,a,rint)
        implicit real *8 (a-h,o-z)
        save
        dimension funjs(130)
c 
c        This subroutine evaluates the integral
c 
c        int_{-rlim}^{rlim} cos(a*x) * log(|x|) dx               (1)
c 
        d1=log(rlim)/a*sin(rlim*a)
        dd=a*rlim
        call sinint(dd,sinin,funjs)
        rint=d1-sinin/(rlim*a) *rlim
        rint=rint*2
c 
        return
        end
c 
c 
c 
c 
c 
        function fun(x,a,par2)
        implicit real *8(a-h,o-z)
        save
  
        fun=cos(a*x)*log(abs(x))
        return
        end
C 
C 
C 
C 
C 
        subroutine trap(x,rint)
        implicit real *8(a-h,o-z)
c 
c       this subroutine uses the end-point corrected
c       trapezoidal quadrature rule to calculate the
c       integral sine rint of the argument x. this is
c       a purely debugging/testing device, using
c       100000 (!!) nodes.
c 
        save
        n=100000
        h=x/(n-1)
        rint=1
        rint=rint/2
        do 1200 i=1,n-1
        d=i*h
        rint=rint+dsin(d)/d
 1200 continue
        rint=rint-dsin(d)/d/2
        rint=rint*h
cccc        call prin2('trapezoidal integral is*',rint,1)
c 
c        now, introduce fourth-order corrections
c 
        d1=-h
        d2=h
        d3=x-h
        d4=x+h
        corr=-dsin(d1)/d1+dsin(d2)/d2+dsin(d3)/d3-dsin(d4)/d4
        rint=rint+corr*h/24
        return
        end
c 
c 
c 
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c       this is the end of the debugging routines and the
c       beginning of the actual code for the calculation of the
c       integral sine function and of the Fourier transform of
c       the logarithm
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c 
        subroutine fourlog(rlim,a,rint)
        implicit real *8 (a-h,o-z)
        save
        dimension funjs(130)
c 
c        This subroutine evaluates the integral
c 
c        int_{-rlim}^{rlim} cos(a*x) * log(|x|) dx               (1)
c 
        d1=log(rlim)/a*sin(rlim*a)
        dd=a*rlim
        call sinint(dd,sinin,funjs)
        rint=d1-sinin/(rlim*a) *rlim
        rint=rint*2
c 
        return
        end
c 
c 
c 
c 
c 
        function fun(x,a,par2)
        implicit real *8(a-h,o-z)
        save
  
        fun=cos(a*x)*log(abs(x))
        return
        end
c 
c 
c 
c 
c 
        subroutine sinint(x,sinin,funjs)
        implicit real *8 (a-h,o-z)
        save
        dimension funjs(1)
c 
c        this subroutine evaluates the integral sine
c        of a real number, defined by the formula
c 
c        sinint(x)=\int_0^x sin(x)/x dx           (1)
c 
c            input parameters:
c 
c  x - the argument of the integral sine to be calculated
c         must be strictly greater than 0.
c 
c            output parameters:
c 
c  sinin - integral sine of x defined by formula (1) above
c 
c            work arrays:
c 
c  funjs - must be at least 120 real *8 locations long.
c 
c    remark:
c 
c       while this subroutine is not particularly fast,
c       it produces at least 13-digit accuracy in all
c       regimes.
c 
c        . . . if the argument is less than 120 - use
c              the expansion of the sine integral into bessel
c              functions
c 
        if(x .gt. 120) goto 2000
        call sinin1(x,funjs(2),sinin)
        return
 2000 continue
c 
c       the argument is greater than 120 - use
c       the asymptotic expansion
c 
        call sinin2(x,sinin)
        return
        end
c 
c 
c 
c 
c 
  
        subroutine sinin2(x,sinin)
        implicit real *8 (a-h,o-z)
        save
        dimension fact(10)
c 
c        this subroutine evaluates the integral sine
c        of a real number, defined by the formula
c 
c        sinint(x)=\int_0^x sin(x)/x dx           (1)
c 
c        this subroutine is to be used for values of
c        x greater than 120, since it uses the asymptotic expansion
c 
c            input parameters:
c 
c  x - the argument of the integral sine to be calculated
c 
c            output parameters:
c 
c  sinin - integral sine of x defined by formula (1) above
c 
c       . . . construct the table of factorials
c 
        fact(1)=1
        do 1200 i=2,10
        fact(i)=fact(i-1)*i
 1200 continue
c 
c       calculate the auxiliary functions f, g via
c       asymptotic expansions
c 
        z=x
        f=1-fact(2)/z**2+fact(4)/z**4-fact(6)/z**6
        g=1-fact(3)/z**2+fact(5)/z**4-fact(7)/z**6
c 
        f=f/z
        g=g/z**2
c 
c       now, calculate the integral sine
c 
        done=1
        pi=datan(done)*4
c 
        sinin=pi/2-f*dcos(z)-g*dsin(z)
        return
        end
c 
c 
c 
c 
c 
        subroutine sinin1(x,funjs,sinin)
        implicit real *8 (a-h,o-z)
        save
        dimension funjs(1)
c 
c        this subroutine evaluates the integral sine
c        of a real number, defined by the formula
c 
c        sinint(x)=\int_0^x sin(x)/x dx           (1)
c 
c        this subroutine should be used for values of
c        x that are .lt. 120 or so, since its cost is
c        proportional to the magnitude of x, and for
c        x .ge. 120, an asymptotic expansion is quite
c        satisfactory.
c 
c            input parameters:
c 
c  x - the argument of the integral sine to be calculated
c 
c            output parameters:
c 
c  sinin - integral sine of x defined by formula (1) above
c 
c            work arrays:
c 
c  funjs - must be at least 120 real *8 locations long, and
c          munst contain the element number 0.
c 
c        . . . find the index at which the bessel function
c              of x/2 becomes effectively zero
c 
        i0=0
        rlarge=1.0d20
        xhalf=x/2
        x1=1
        x1=x1/xhalf
        funjs(i0)=0
        funjs(1)=1
        do 1200 i=1,1000
        k=i
        funjs(i+1)=(2*i+1)*x1*funjs(i)-funjs(i-1)
        if(dabs(funjs(i+1)) .gt. rlarge) goto 1400
 1200 continue
 1400 continue
cccc          call prinf('in sinin1, k=*',k,1)
c 
c       construct the unscaled bessel functions of x/2
c 
        funjs(k)=0
        funjs(k-1)=1
        do 1600 i=k-1,1,-1
        funjs(i-1)=(2*i+1)*x1*funjs(i)-funjs(i+1)
 1600 continue
c 
c        determine the scaling coefficient
c 
        fj0=sin(xhalf)/xhalf
        fj1=sin(xhalf)/xhalf**2-cos(xhalf)/xhalf
        scale=(fj0**2+fj1**2)/(funjs(i0)**2+funjs(1)**2)
        scale=sqrt(scale)
c 
c        calculate the correct values of the bessel functions
c 
        do 1800 i=0,k
        funjs(i)=funjs(i)*scale
 1800 continue
cccc        call prin2('in sinin1, funjs=*',funjs(i0),k)
c 
c        now, calculate the integral sine of x
c 
        sinin=0
        do 2000 i=0,k-1
        sinin=sinin+funjs(i)**2
 2000 continue
         sinin=sinin*x
         return
        end
  
  
  
  
  
  
  
  
  
  
  
  
