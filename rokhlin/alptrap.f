        implicit real *8 (a-h,o-z)
        dimension t(10000),f(10000)
        external fun,funini
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
        PRINT *, 'ENTER k '
        READ *,k
        CALL PRINF('k=*',k,1 )
c 
        PRINT *, 'ENTER jtest '
        READ *,jtest
        CALL PRINF('jtest=*',jtest,1 )
c 
c        plot the smooth step function
c 
        done=1
        pi=datan(done)*4
        h=pi/2/(m-1)
c 
        h=done/(m-1)
        do 1200 i=1,m
        t(i)=(i-1)*h
 1200 continue
c 
c        verify alpert's end-point corrections
c 
c        . . . construct the values of the function
c              to be integrated
c 
        do 1400 i=1,m
        f(i)=t(i)**jtest
 1400 continue
c 
c        now, calculate the integral
c 
        call alptrap(ier,f,m,h,k,trap)
        call prinf('after alptrap, ier=*',ier,1)
        call prin2('after alptrap, trap=*',trap,1)
        call prin2('and trap-1/(jtest+1)=*',
     1      trap-done/(jtest+1),1)
c 
c       now, calculate the integral using the other subroutine
c 
        d=funini(jtest)
        a=0
        b=1
        call alptrap2(ier,fun,m,a,b,k,trap2)
        call prin2('after alptrap2, trap2=*',trap2,1)
        call prin2('and trap2-1/(jtest+1)=*',
     1      trap2-done/(jtest+1),1)
  
  
  
  
c 
        stop
        end
c 
c 
c 
c 
c 
        function fun(x)
        implicit real *8 (a-h,o-z)
c 
        save
        fun=x**jtest
        return
c 
c 
c 
c 
        entry funini(jtest7)
        jtest=jtest7
        return
        end
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c        this is the end of the debugging code and the
c        beginning of the actual quadrature coded
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        subroutine alptrap2(ier,fun,n,a,b,m,trap)
        implicit real *8 (a-h,o-z)
        save
        dimension weights(12),w(12),w1(12)
        data m7/-2341034/
c 
c       this subroutine uses an end-point corrected
c       trapezoidal rule to evaluate the integral of
c       user-supplied function fun on the interval [a,b].
c       the user supplies the function to be integrated
c       in the form of the function name fun with the calling
c       sequence fun(x)
c 
c                   input parameters:
c 
c  fun - the real *8 function that will be called by this
c       subroutine at n equispaced nodes on the interval
c       [a,b], on which the integral will be calculated.
c       note that the function will be evaluated at the ends
c       of this interval.
c  n - the number of nodes in the end-point corrected trapezoidal
c       rule that the subroutine will use.
c  a - the left end of the interval in which the function will be
c       integrated
c  b - the right end of the interval in which the function will be
c       integrated
c  m - the order of the quadrature rule. must be 2,4,6,8,10, or 12.
c       if any other order is specified, the error return code ier
c       is set to 8, and the execution of the subroutine is
c       terminated.
c 
c                   output parameters:
c 
c  ier - error return code.
c            ier=0 means normal conclusion.
c            ier=8 means that the order m of the quadrature
c                  specified is not a permitted one (i.e.
c                  m is not equal to 2,4,6,8,10, or 12.
c            ier=16 means that n < m.
c 
c        . . . if this is the first call to this subroutine
c              with this m - get the weights
c 
        h=(b-a)/(n-1)
        t0=a
        if(m .eq. m7) goto 1100
        call alptrap0(ier,m,weights)
        if(ier .ne. 0) return
        m7=m
 1100 continue
c 
c        check if the number of points in the user-supplied
c        discretization of the functiojn exceeds the order
c        requested by the user
c 
        if(n .ge. m) goto 1150
        ier=16
        return
 1150 continue
c 
c        calculate the trapezoidal approximation
c 
        trap=0
        do 1200 i=1,m-1
        t=t0+(i-1)*h
        d=fun(t)
        w(i)=d
        trap=trap+d
 1200 continue
c 
        do 1250 j=1,m-1
        i=n-j+1
        t=t0+(i-1)*h
        d=fun(t)
        w1(j)=d
        if(i .ge. m) trap=trap+d
 1250 continue
c 
        do 1300 i=m,n-m+1
        t=t0+(i-1)*h
        d=fun(t)
        trap=trap+d
 1300 continue
c 
        trap=trap-(w(1)+w1(1))/2
c 
c       now, apply the end-point corrections
c 
        do 1400 i=1,m-1
        trap=trap+(w(i)+w1(i))*weights(i)
 1400 continue
        trap=trap*h
        return
        end
c 
c 
c 
c 
c 
        subroutine alptrap(ier,f,n,h,m,trap)
        implicit real *8 (a-h,o-z)
        save
        dimension f(1),weights(12)
        data m7/-2341034/
c 
c       this subroutine uses an end-point corrected
c       trapezoidal rule to evaluate the integral of
c       user-supplied function. the user supplies the
c       function to be integrated in the form of a table
c       of equispaced values; it is assumed that the
c       end-point values are among those supplied.
c 
c                   input parameters:
c 
c  f - the table of values of the function to be integrated
c  n - the number of values in the array f
c  h - the sampling interval in the table f
c  m - the order of the quadrature rule. must be 2,4,6,8,10, or 12.
c       if any other order is specified, the error return code ier
c       is set to 8, and the execution of the subroutine is
c       terminated.
c 
c                   output parameters:
c 
c  ier - error return code.
c            ier=0 means normal conclusion.
c            ier=8 means that the order m of the quadrature
c                  specified is not a permitted one (i.e.
c                  m is not equal to 2,4,6,8,10, or 12.
c            ier=16 means that n < m.
c 
c        . . . if this is the first call to this subroutine
c              with this m - get the weights
c 
        if(m .eq. m7) goto 1100
        call alptrap0(ier,m,weights)
        if(ier .ne. 0) return
        m7=m
 1100 continue
c 
c        check if the number of points in the user-supplied
c        discretization of the functiojn exceeds the order
c        requested by the user
c 
        if(n .ge. m) goto 1150
        ier=16
        return
 1150 continue
c 
c        calculate the trapezoidal approximation
c 
        trap=0
        do 1200 i=1,n
        trap=trap+f(i)
 1200 continue
c 
        trap=trap-(f(1)+f(n))/2
c 
c       now, apply the end-point corrections
c 
        do 1400 i=0,m-2
        trap=trap+(f(1+i)+f(n-i))*weights(i+1)
 1400 continue
        trap=trap*h
        return
        end
c 
c 
c 
c 
c 
        subroutine alptrap0(ier,m,weights)
        implicit real *8 (a-h,o-z)
        save
        dimension weights(1)
        integer *4 c4(4),c6(6),c8(8),c10(10),c12(12),cc12(12)
        data c4/24,-3,4,-1/,c6/1440,-245,462,-336,146,-27/,
     1      c8/120960,-23681,55688,-66109,57024,
     2      -31523,9976,-1375/,
     3      c10/7257600,-1546047,4274870,-6996434,9005886,
     4      -8277760,5232322,-2161710,526154,-57281/
c 
            data c12/9580032,-2162543,6795432,-14129473,
     6       24158814,-31035790,29399424,-20232241,
     7       9845153,-3214558,632535,-56752/
c 
            data cc12/00,-35,84,-89,
     6       96,-86,00,-14,
     7       04,-11,16,-65/
c 
c       this subroutine returns to the user the end-point
c       corrections for the trapezoidal rule, resulting
c       in a rule of order m
c 
        ier=8
c 
c       . . . if m=2
c 
        if(m .ne. 2) goto 1400
c 
        ier=0
        weights(1)=0
        return
 1400 continue
c 
c       . . . if m=4
c 
        if(m .ne. 4) goto 2400
c 
        ier=0
        d=c4(1)
        do 2200 i=1,m-1
        weights(i)=c4(i+1)/d
 2200 continue
        return
 2400 continue
c 
c 
c       . . . if m=6
c 
        if(m .ne. 6) goto 3400
c 
        ier=0
        d=c6(1)
        do 3200 i=1,m-1
        weights(i)=c6(i+1)/d
 3200 continue
        return
 3400 continue
c 
c 
c       . . . if m=8
c 
        if(m .ne. 8) goto 4400
c 
        ier=0
        d=c8(1)
        do 4200 i=1,m-1
        weights(i)=c8(i+1)/d
 4200 continue
        return
 4400 continue
c 
c 
c       . . . if m=10
c 
        if(m .ne. 10) goto 5400
c 
        ier=0
        d=c10(1)
        do 5200 i=1,m-1
        weights(i)=c10(i+1)/d
 5200 continue
        return
 5400 continue
c 
c 
c       . . . if m=12
c 
        if(m .ne. 12) goto 6400
c 
        ier=0
        d=c12(1)*100+cc12(1)
        do 6200 i=1,m-1
        weights(i)=c12(i+1)/d*100+cc12(i+1)/d
 6200 continue
        return
 6400 continue
c 
        return
        end
