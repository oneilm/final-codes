        implicit real *8 (a-h,o-z)
        dimension ts(10 000),fs(10 000)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER n'
        READ *,n
        CALL PRInf('n=*',n,1 )
c 
        PRINT *, 'ENTER x'
        READ *,x
        CALL PRIn2('x=*',x,1 )
c 
c       construct the test nodes and function values at them
c 
        m=160
        b=20
  
c 
        done=1
        pi=datan(done)*4
        a=2*pi
c 
        h=a/m
        do 1200 i=1,m
        ts(i)=(i-1)*h
c 
        fs(i)=cos(b*ts(i))
 1200 continue
c 
c       interpolate the function at the point x
c 
        ifinit=1
  
        call perilagr(ts,fs,m,x,n,fout,ifinit)
  
        call prin2('after eqlagrev, fout=*',fout,1)
  
        d=cos(b*x)
  
        call prin2('and f(x) directly is=*',d,1)
        call prin2('and the difference is*',fout-d,1)
  
cccc        call prin2('and the differences*',errs,nnn)
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine perilagr(ts,fs,m,x,n,fout,ifinit)
        implicit real *8 (a-h,o-z)
        save
        real *8 ts(1),coefs(100),coefsx(100),fs(1),
     1      tsnew(100),fsnew(100)
c 
c       this subroutine interpolates (via standard Lagrange
c       interpolation) a user-specified perioudic function
c       from the user-provided equispaced grid to the user-supplied
c       point x on the line. The point is permitted to be outside the
c       interval where the discretization is located, but by no more
c       than the length of the interval on which the function is
c       originally defined (i.e. by no more than ts(m)-ts(1)). Please
c       note that the subpoutine assumes that ts(m) .neq. ts(1),
c       just like standard FFT routines do; conceptually,
c       ts(m+1)=ts(1).
c 
c 
c                    Input parameters:
c 
c  ts - the equispaced nodes at which the user-supplied function is
c       tabulated
c  fs - the values of the interpolant at the nodes ts
c  m - the number of nodes in the array ts
c  x - the point where the function is to be interpolated
c  ifinit - an integer parameter telling the subroutine whether the
c       Lagrange interpolation routine has to be initialized.
c          ifinit=1 will cause the initialization to be performed
c          ifinit=0 will cause the initialization to be skipped.
c       Please note that ifinit has to be set to one each time the
c       sampling distance in the array ts changes
c 
c                     Output parameters:
c 
c  fout - the interpolated value of the function at the point x
c 
c 
c 
c       if the user so requested - initialize the interpolation routine
c 
        if(ifinit .eq. 1) call lagrini(ts,n,coefs)
c 
c       use interpolation to find the value of the function
c       for the user-specified x
c 
c       . . . find on which subinterval x lives
c 
        ier=0
c 
        h=ts(2)-ts(1)
        iint=(x-ts(1))/h
        istart=iint-n/2+2
        iend=istart+n-1
c 
       do 2300 i=1,n
c 
       tsnew(i)=ts(1)+(istart+i-2)*h
       fsnew(i)=fs(iwrap(m,istart+i-1))
 2300 continue
c 
        call lagrcoef(tsnew,coefs,n,x,coefsx)
c 
        fout=0
        do 2400 i=1,n
        fout=fout+coefsx(i)*fsnew(i)
 2400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        function iwrap(m,i)
c 
        save
        if( (i .le. m) .and. (i .ge. 1)) then
           iwrap=i
           return
        endif
c 
        if(i .lt. 1) iwrap=m+i
        if(i .gt. m) iwrap=i-m
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine eqlagrev(ier,ts,fs,m,x,n,fout,ifinit)
        implicit real *8 (a-h,o-z)
        save
        real *8 ts(1),coefs(100),coefsx(100),fs(1)
c 
c       this subroutine interpolates (via standard Lagrange
c       interpolation) a user-specified function from the
c       user-provided equispaced grid to the user-supplied
c       point x on the line.
c 
c                    Input parameters:
c 
c  ts - the equispaced nodes at which the user-supplied function is
c       tabulated
c  fs - the values of the interpolant at the nodes ts
c  m - the number of nodes in the array ts
c  x - the point where the function is to be interpolated
c  ifinit - an integer parameter telling the subroutine whether the
c       Lagrange interpolation routine has to be initialized.
c          ifinit=1 will cause the initialization to be performed
c          ifinit=0 will cause the initialization to be skipped.
c       Please note that ifinit has to be set to one each time the
c       sampling distance in the array ts changes
c 
c                     Output parameters:
c 
c  ier - error return code.
c              ier=0 means no trouble
c              ier=4 means no trouble that x < ts(1)
c              ier=8 means no trouble that x > ts(1)
c       None of these conditions are fatal; however, when a non-zero
c       error code is returned, great care should be exercised (Runge
c       phenomenon and related trouble)
c  fout - the interpolated value of the function at the point x
c 
c 
c       if the user so requested - initialize the interpolation routine
c 
        if(ifinit .eq. 1) call lagrini(ts,n,coefs)
c 
c       use interpolation to find the value of the function
c       for the user-specified x
c 
c       . . . find on which subinterval x lives
c 
        ier=0
        if(x .lt. ts(1)) then
           ier=4
           istart=1
           goto 2200
        endif
c 
        if(x .gt. ts(m)) then
           ier=8
           istart=m-n+1
           goto 2200
        endif
c 
        h=ts(2)-ts(1)
        iint=(x-ts(1))/h
        istart=iint-n/2+2
        if(istart .lt. 1) istart=1
        iend=istart+n-1
        if(iend .gt. m) istart=m-n+1
c 
 2200 continue
c 
cccc        call prinf('and istart=*',istart,1)
c 
        call lagrcoef(ts(istart),coefs,n,x,coefsx)
c 
cccc        call prin2('and coefsx=*',coefsx,n)
cccc        call prin2('while ts(istart)=*',ts(istart),n)
  
        fout=0
        do 2400 i=1,n
        fout=fout+coefsx(i)*fs(istart+i-1)
 2400 continue
c 
        return
        end
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
c 
c       calculate the coefficient in front of the sum
c 
        coesum=1
        do 4200 i=1,n
        coesum=coesum*(x-t(i))
 4200 continue
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
        coefsx(knear)=snear*coefs(knear)
c 
        return
        end
