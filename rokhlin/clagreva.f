        implicit real *8 (a-h,o-z)
        complex *16 ima,zs(100),coefs(100),z,f,f2,fun,
     1      fs(100)
c 
        data ima/(0.0d0,1.0d0)/
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
         PRINT *, 'ENTER n'
         READ *,n
         CALL PRINF('n=*',n,1 )
C 
c        construct the nodes on which the interpolation is
c        to be based and the values to be interpolated
c 
         zs(1)=1
        zs(2)=2
        zs(3)=3
cccc        n=3
c 
        done=1
        pi=datan(done)*4
        h=2*pi/n
        do 1200 i=1,n
        t=(i-1)*h
        zs(i)=cdexp(ima*t)
 1200 continue
c 
        do 1400 i=1,n
        fs(i)=fun(zs(i))
 1400 continue
c 
c        initialize the interpolation routine
c 
        call clagrini(zs,n,coefs)
c 
c       interpolate the function fun at the pints z
c 
        z=0.7+ima*0.56
ccc        z=z*2
        call clagreva(zs,coefs,n,z,fs,f)
        call prin2('f from clagreva is*',f,2)
        f2=fun(z)
        call prin2('and f2 obtained directly*',f2,2)
        call prin2('and f2-f=*',f2-f,2)
  
  
c 
  
  
        stop
        end
c 
c 
c 
c 
c 
        function fun(x)
        save
        complex *16 fun,x
        fun=cdexp(x)
        return
        end
c 
c 
c 
c 
c 
        subroutine clagrini(t,n,coefs)
        implicit real *8 (a-h,o-z)
        save
        complex *16 t(1),coefs(1),f(1),cd,
     1      coesum,sum,snear,fun,x
c 
c       this subroutine initializes the lagrange interpolation
c       subroutine in the complex plane. the actual interpolation
c       is performed by the entry clagreva (see below).
c 
c                         input parameters:
c 
c  t - the nodes on which the inetrpolation is based
c  n - the number of elements in array t
c 
c                         output parameters:
c 
c  coefs -  an array of length n to be used by the entry clagreva
c 
c        . . . initialize the interpolation routine based on
c              user-specified nodes t(i)
c 
        done=1
        do 1400 i=1,n
        cd=1
        do 1200 j=1,n
        if(j.eq.i) goto 1200
        cd=cd*(t(i)-t(j))
 1200 continue
        coefs(i)=done/cd
 1400 continue
        call prin2('in clagrini, coefs=*',coefs,n*2)
        return
c 
c 
c 
c 
       entry clagreva(t,coefs,n,x,f,fun)
c 
c       this entry interpolates a function f tabulated at the
c       nodes t(i) in the complex plane to obtain  the value of f at
c       the point x.
c 
c                       input parameters:
c 
c  t - the nodes on which the interpolation is based.
c  coefs -  an array of length n produced by the entry clagrini
c  n - the number of elements in array t
c  x - the point in the complex plane where f is to be evaluated
c  f - the values of the function f at the nodes t(i)
c 
c                         output parameters:
c 
c  fun - the inetrpolated value of f at the point x
c 
c       . . . find the node nearest to this point
c 
        knear=1
        cd=x-t(1)
        dd=cd*dconjg(cd)
        dist=dd
        do 2400 i=2,n
        cd=x-t(i)
        dd=cd*dconjg(cd)
        if(dd .gt. dist) goto 2400
        dist=dd
        knear=i
 2400 continue
        call prinf('knear=*',knear,1)
c 
c       calculate the coefficient in front of the sum
c 
ccc          call prinf('knear as found is*',knear,1)
        coesum=1
        do 4200 i=1,n
        coesum=coesum*(x-t(i))
 4200 continue
        call prin2('coesum=*',coesum,2)
c 
c       calculate the sum excluding the term corresponding
c       to knear
c 
        sum=0
        do 4400 i=1,n
        if(i .eq. knear) goto 4400
c 
        sum=sum+coefs(i)*f(i)/(x-t(i))
 4400 continue
        sum=sum*coesum
c 
c       account for the term knear
c 
        snear=1
        do 4600 i=1,n
        if(i .eq. knear) goto 4600
        snear=snear*(x-t(i))
 4600 continue
        call prin2('snear=*',snear,2)
        snear=snear*coefs(knear)*f(knear)
        call prin2('snear=*',snear,2)
c 
        fun=sum+snear
        return
        end
