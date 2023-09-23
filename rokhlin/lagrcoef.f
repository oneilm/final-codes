c 
c 
c 
c 
c 
        subroutine lagrini(t,n,coefs)
        implicit real *8 (a-h,o-z)
        save
        real *8 t(1),coefs(1),coefsx(1),cd,
     1      coesum,snear,x
c 
c       this subroutine performs initialization for the entry
c       lagrcoef (see below). its output is the array coefs
c       to be used by lagrcoef to construct interpolation
c       coefficients for the interpolation on the real line.
c 
c                         input parameters:
c 
c  t - the nodes on which the interpolation is based.
c  n - the number of elements in array t
c 
c                         output parameters:
c 
c  coefs -  an array of length n to be used by the entry lagreva
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
cccc        call prin2('in lagrini, coefs=*',coefs,n)
        return
c 
c 
c 
c 
       entry lagrcoef(t,coefs,n,x,coefsx)
c 
c       this entry constructs interpolation coiefficients
c       connecting the values of a function f tabulated at the
c       nodes t(i) on the real line  with the value of f at
c       the real polint x
c 
c                       input parameters:
c 
c  t - the nodes on which the interpolation is based.
c  coefs -  an array of length n produced the entry lagrini (see)
c  n - the number of elements in arrays t, coefs
c  x - the point in the complex plane where f is to be evaluated
c 
c                         output parameters:
c 
c  coefsx - the array of interpolation coefficients
c 
c       . . . find the node nearest to this point
c 
        knear=1
        cd=x-t(1)
        dd=cd**2
        dist=dd
        do 2400 i=2,n
        cd=x-t(i)
        dd=cd**2
        if(dd .gt. dist) goto 2400
        dist=dd
        knear=i
 2400 continue
cccccc        call prinf('knear=*',knear,1)
c 
c       calculate the coefficient in front of the sum
c 
ccc          call prinf('knear as found is*',knear,1)
        coesum=1
        do 4200 i=1,n
        coesum=coesum*(x-t(i))
 4200 continue
cccccc        call prin2('coesum=*',coesum,1)
c 
c       calculate the sum excluding the term corresponding
c       to knear
c 
        do 4400 i=1,n
        if(i .eq. knear) goto 4400
c 
        coefsx(i)=coefs(i)*coesum/(x-t(i))
 4400 continue
c 
c       account for the term knear
c 
        snear=1
        do 4600 i=1,n
        if(i .eq. knear) goto 4600
        snear=snear*(x-t(i))
 4600 continue
cccccc        call prin2('snear=*',snear,1)
        coefsx(knear)=snear*coefs(knear)
cccccc        call prin2('snear=*',snear,1)
c 
        return
        end
