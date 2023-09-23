        implicit real *8 (a-h,o-z)
        complex *16 ima,zs(100),coefs(100),z,f,f2,fun,
     1      fs(100),coefsx(100),a(42000)
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
        call lagrcini(zs,n,coefs)
c 
c       interpolate the function fun at the points z
c 
        z=0.7+ima*0.56
cccc        z=z*2
        call lagrcoef(zs,coefs,n,z,coefsx)
         call prin2('after lagrcoef, coefsx=*',coefsx,n*2)
c 
        call expeva(coefsx,fs,n,f)
        call prin2('f numerically is*',f,2)
c 
        f2=fun(z)
        call prin2('and f2 obtained directly*',f2,2)
        call prin2('and f2-f=*',f2-f,2)
c 
c       construct the whole interpolation matrix
c 
        ir=103
        call mattest(a,n,ir)
  
  
c 
  
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine mattest(a,n,ir)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,n),x(1000),y(1000),w(10000),
     1      f1(1000),f2(1000),f3(1000),fun2
        dimension xwhts(1000),ywhts(1000)
c 
c       read from disk the arrays x, y
c 
        call pntsrd(ir,nroots,x,xwhts,
     1      xtheta,r)
c 
        call prin2('after pntsrd, x=*',x,nroots*2)
c 
        call pntsrd(ir,nroots,y,ywhts,
     1      ytheta,r)
        call prin2('after pntsrd, y=*',y,nroots*2)
c 
c       now, construct the matrix converting the values
c       of an analytical function on x into its values on y
c 
        call matinter(x,y,n,a,w)
        call prin2('after matinter, a=*',a,n*n*2)
c 
c       test the interpolation matrix
c 
c       . . . construct values of a function at the points x and y
c 
        do 2200 i=1,n
        f1(i)=fun2(x(i))
        f2(i)=fun2(y(i))
 2200 continue
c 
c        . . . interpolate
c 
        call cmatvec(a,f1,f3,n)
        call prin2('f1=*',f1,n*2)
        call prin2('f2=*',f2,n*2)
        call prin2('f3=*',f3,n*2)
  
        return
        end
c 
c 
c 
c 
c 
        function fun2(z)
        implicit real *8 (a-h,o-z)
        save
        complex *16 z,fun2
c 
        fun2=z**9
        return
        end
c 
c 
c 
c 
c 
        subroutine cmatvec(a,x,y,n)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,n),x(1),y(1),cd
c 
        do 1400 i=1,n
        cd=0
        do 1200 j=1,n
        cd=cd+a(i,j)*x(j)
 1200 continue
        y(i)=cd
 1400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine expeva(coefsx,fs,n,f)
        implicit real *8 (a-h,o-z)
        save
        complex *16 coefsx(1),fs(1),f
c 
        f=0
        do 1200 i=1,n
        f=f+coefsx(i)*fs(i)
 1200 continue
        return
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
        subroutine matinter(x,y,n,a,w)
        implicit real *8 (a-h,o-z)
        save
        complex *16 x(1),y(1),a(n,n),w(1)
c 
c        this subroutine constructs a matrix converting the values
c        of a complex analytic function at the points x into its
c        values at the points y. the choice of the points x that
c        would assure the stability of this procedure is the
c        user's responsibility.
c 
c                     input parameters:
c 
c  x - the points from which the interpolation is to be performed
c  y - the points onto which the interpolation is to be performed
c  n - the number of elements in each of the arrys x, y
c 
c                     output parameters:
c 
c  a - the n \times n matrix interpolating functions from the
c        points x onto points y
c 
c                     work arrays:
c 
c  w - must be at least 2*n+2 complex *16 locations long
c 
c       . . . initialize the lagrange interpolation subroutine
c 
        call lagrcini(x,n,w)
c 
c       one column after another, construct the interpolation
c       matrix
c 
        do 2000 i=1,n
c 
c       construct a row of the interpolation matrix
c 
        icoefsx=n+1
        call lagrcoef(x,w,n,y(i),w(icoefsx))
c 
c       store this row
c 
        do 1600 j=1,n
        a(i,j)=w(icoefsx+j-1)
 1600 continue
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine lagrcini(t,n,coefs)
        implicit real *8 (a-h,o-z)
        save
        complex *16 t(1),coefs(1),coefsx(1),cd,
     1      coesum,snear,x
c 
c       this subroutine performs initialization for the entry
c       lagrcoef (see below). its output is the array coefs
c       to be used by lagrcoef to construct interpolation
c       coefficients for the interpolation in the complex plane.
c 
c                         input parameters:
c 
c  t - the nodes on which the inetrpolation is based. they must
c      be in increasing order.
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
        call prin2('in lagrini, coefs=*',coefs,n*2)
        return
c 
c 
c 
c 
       entry lagrcoef(t,coefs,n,x,coefsx)
c 
c       this entry constructs interpolation coiefficients
c       connecting the values of a function f tabulated at the
c       nodes t(i) in the complex plane with the value of f at
c       the point x in the complex plane
c 
c                       input parameters:
c 
c  t - the nodes on which the interpolation is based.
c  coefs -  an array of length n produced the entry lagrcini (see)
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
        call prin2('snear=*',snear,2)
        coefsx(knear)=snear*coefs(knear)
        call prin2('snear=*',snear,2)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine pntswrt(iw,nroots,zroots,weights,
     1      theta,r)
        implicit real *8 (a-h,o-z)
        save
        dimension zroots(2,1),weights(1)
c 
c       write the header of the record
c 
 1200 format('   nroots= ',i6,'   theta= ',d22.16,
     1    '  r= ',d22.16)
        write(iw,1200) nroots, theta,r
c 
c        now, write the locations of the points in the plane
c 
 1400 format('         the nodes in the plane are')
        write(iw,1400)
 1600 format(2(2x,e22.16))
        write(iw,1600) (zroots(1,i),zroots(2,i),i=1,nroots)
c 
c       . . . and the weights corresponding to these roots
c 
 1800 format('         weights corresponding to above roots')
        write(iw,1800)
        write(iw,1600) (weights(i),i=1,nroots)
        return
c 
c 
c 
c 
        entry pntsrd(ir,nroots,zroots,weights,
     1      theta,r)
c 
c       read the header from disk
c 
ccc 2200 format(11x,i6,8x,d22.16,5x,d22.16)
 2200 format(11x,i6,10x,d22.16,5x,d22.16)
        read (ir,2200) nroots,theta,r
        call prinf('nroots as read*',nroots,1)
        call prin2('theta as read*',theta,1)
        call prin2('r as read*',r,1)
cccc          if(2 .ne. 3) stop
c 
c       read the locations of the points
c 
 2400 format(10x)
 2600 format(2(2x,e22.16))
        read(ir,2400)
        read(ir,2600) (zroots(1,i),zroots(2,i),i=1,nroots)
c 
c       read the weights
c 
        read(ir,2400)
        read(ir,2600) (weights(i),i=1,nroots)
c 
        return
        end
c 
c 
