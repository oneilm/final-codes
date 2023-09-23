         implicit real *8 (a-h,o-z)
         complex *16 ima,fjs(10000),z
         data ima/(0.0d0,1.0d0)/
c 
        CALL PRINI(6,13)
C 
c       . . . construct the testing parameters
c 
         PRINT *, 'ENTER x'
         READ *,x
         CALL PRIn2('x=*',x,1)
c 
c        test the spherical Bessel function computer
c 
        z=x*ima
        z=x
        eps=1.0d-15
        call jfuns3d(z,eps,fjs(2),n)
        call prin2('after jfuns3d, fjs=*',fjs,n*2)
  
  
         stop
         end
c 
c 
c 
c 
c 
        subroutine jfuns3d(z,eps,fjs,n)
        implicit real *8 (a-h,o-z)
        save
        complex *16 z,fjs(1),zinv,com,fj0,scale,fj1
        dimension rea(2)
        equivalence (rea(1),com)
c 
c        this subroutine evaluates all spherical Bessel functions
c        of the complex argument z that are greater (in absolute
c        value) than eps.
c 
c                        input parameters:
c 
c  z - the argument
c  eps - machine zero - all Besel functions that are smaller
c        than eps are declared to be zero.
c 
c                        output parameters:
c 
c  fjs - the array of Bessel functions.
c 
c      Attention: the subroutie assumes that this array has
c                 element number zero!!!!!!
c 
c  n - the number of elements in fjs that are greater (in
c      absolute value) than eps.
c 
c       . . . conduct the recursion up to determine how high
c             we should go
c 
        done=1
        i0=0
        fjs(i0)=0
        fjs(1)=1
        rlarge=done/eps*1.0d17
        zinv=done/z
c 
        do 1200 i=1,10000000
        nn=i
        fjs(i+1)=(2*i+done)*zinv*fjs(i)-fjs(i-1)
c 
        com=fjs(i+1)
        dd=dabs(rea(1))+dabs(rea(2))
        if(dd .gt. rlarge) goto 1400
 1200 continue
c 
 1400 continue
        call prinf('nn=*',nn,1)
c 
c       conduct recursion down
c 
        fjs(nn+1)=0
        fjs(nn)=1
        do 2200 i=nn,1,-1
c 
        fjs(i-1)=(2*i+done)*zinv*fjs(i)-fjs(i+1)
 2200 continue
         call prin2('fjs before scaling*',fjs(i0),nn*2+2)
c 
c       determine the scaling parameter
c 
        fj0=cdsin(z)*zinv
        fj1=fj0*zinv-cdcos(z)*zinv
c 
        com=fj0
        d0=dabs(rea(1))+dabs(rea(2))
        com=fj1
        d1=dabs(rea(1))+dabs(rea(2))
c 
        if(d1 .gt. d0) goto 2250
        scale=fj0/fjs(i0)
        goto 2300
 2250 continue
c 
        scale=fj1/fjs(1)
 2300 continue
        call prin2('scale=*',scale,2)
c 
c       . . . scale
c 
        do 2400 i=i0,nn
        fjs(i)=fjs(i)*scale
 2400 continue
c 
c       find the number of Bessel functions that are greater than eps
c 
        do 2600 i=nn,1,-1
        n=i
        com=fjs(i)
        dd=dabs(rea(1))+dabs(rea(2))
        if(dd .gt. eps) goto 2800
 2600 continue
 2800 continue
        return
        end
  
  
