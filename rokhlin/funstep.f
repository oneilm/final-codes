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
c 
c        plot the smooth step function
c 
        done=1
        pi=datan(done)*4
        h=pi/2/(m-1)
c 
        do 1200 i=1,m
        t(i)=(i-1)*h
        call funstep(t(i),f(i))
        if(dabs(f(i)) .lt. 1.0d-20) f(i)=0
 1200 continue
c 
        iw=36
        call RSGRAF(t,f,m,IW,'step function*')
  
  
  
  
c 
        stop
        end
c 
c 
c 
c 
c 
        subroutine funstep(x,f)
        implicit real *8 (a-h,o-z)
        save
        data ifcalled/0/
c 
c       this subroutine evaluates an elementary function f
c       such that it is
c 
c   1. c^{\infty} on the real line.
c   2. f(x)=1 for all x .le. 0
c   3. f(x)=0 for all x .ge. \pi/2
c 
c       the said function is defined by the formulae
c 
c        d=dtan(x)
c        dd=1-dexp(-1/d)
c        f=dexp(-1/dd+1)
c 
c                  input parameters:
c 
c  x - the point on the interval [0, \pi/2] where the function
c       is to be evaluated
c 
c                  output parametedrs:
c 
c  f - the value of the function
c 
        if(ifcalled .eq. 1) goto 1100
        ifcalled=1
        done=1
        pi=datan(done)*4
 1100 continue
c 
c       if x .le. 0  set the step function to one
c 
        if(x .gt. 0) goto 1200
        f=1
        return
 1200 continue
c 
c       if x .gt pi/2 - set the function to zero
c 
        if(x .lt. pi/2) goto 1400
        f=0
        return
 1400 continue
c 
c       x is in the transition region - calculate the transition function
c 
        d=dtan(x)
        dd=done-dexp(-done/d)
        f=dexp(-done/dd+done)
        return
        end
