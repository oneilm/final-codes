       implicit real *8 (a-h,o-z)
        complex *16 ima,fs(10 000),zs(10 000),diffs(10000),
     1      poles(100),fsbomb(100)
        dimension ts(10000),fsrea(10000),fsima(10000),
     1      deltas(100)
  
  
  
c 
        data ima/(0.0d0,1.0d0)/
c 
        call prini(6,13)
c 
        PRINT *, 'ENTER delta'
        READ *,delta
        CALL PRIN2('delta=*',delta,1 )
c 
c       test the accuracy of the approximation on the interval [-1,1]
c 
        c=1.5
        c=1
c 
        h=2
        n=100
        h=h/(n-1)
        do 1200 i=1,n
c 
        ts(i)=-1+(i-1)*h
c 
        zs(i)=ts(i)
c 
        call moderfun3(c,delta,zs(i),fs(i))
  
  
  
        diffs(i)=fs(i)-1
 1200 continue
c 
        call prin2('ts=*',ts,n)
        call prin2('fs=*',fs,n*2)
        call prin2('diffs=*',diffs,n*2)
  
  
  
  
  
  
  
  
  
        call moderfun_data(c,delta,poles)
  
  
        call prin2('after moderfun_data, poles=*',poles,16)
c 
        do 1400 i=1,8
c 
        call moderfun3(c,delta,poles(i),fsbomb(i))
  
 1400 continue
  
        call prin2('and fsbomb=*',fsbomb,16)
  
  
  
c 
c       plot the values of the moderating function on
c       the left wing
c 
        h=100
        h=20
        n=100
        h=h/(n-1)
        do 2200 i=1,n
c 
        ts(i)=(i-1)*h
c 
        zs(i)=-1+(i-1)*h * ima
c 
        call moderfun3(c,delta,zs(i),fs(i))
  
  
  
  
        fsrea(i)=fs(i)
        fsima(i)=-ima*fs(i)
  
c 
cccc        fsrea(i)=log10(abs(fsrea(i))+1.0d-15)
cccc        fsima(i)=log10(abs(fsima(i))+1.0d-15)
 2200 continue
c 
        call prin2('fs=*',fs,n*2)
        call prin2('zs=*',zs,n)
  
        iw=21
  
        call quagraph2(iw,ts,fsrea,n,3,ts,fsima,n,3,
     1      'moderator at -ima+y*')
  
  
  
  
  
c 
        stop
        end
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code, and the beginning
c        of the "moderator" function code proper.
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c        There are two user-callable subroutines in this file:
c        moderfun3 and moderfun_data. Moderfun3 evaluates the
c        moderating function at the user-supplied point;
c        moderfun_data returns to the user 8 poles of the
c        moderating function.
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        subroutine moderfun3(c,delta,z,f)
        implicit real *8 (a-h,o-z)
        save
        complex *16 z,f,f2,ima,cw
        data done/1.0d0/,ima/(0.0d0,1.0d0)/
c 
c        This subroutine returns to the user the value of a
c        "moderator" function at the user-supplied complex
c        point z. The "moderator" function is close to 1 on
c        the real axis, and decays rapidly in the upper half
c        -plane. The "moderator" function is determined by
c        the parameters c, delta in a fairly complex manner.
c        In the upper half-plane, it decays as exp(-c*y);
c        it has an infinite number of poles, all of which are
c        located in
c        the upper half-plane, in groups of 8. Within each group,
c        the poles sre on a vertical line; such vertical lines
c        are spaced with the interval 2*pi/c. The subroutine
c        moderfun_data returns one such group, with x=pi/c.
c 
c                Input parameters:
c 
c  c, delta - the parameters defining the functions. Recommended
c        values: c=1, delta=0.0002 give roughly 15-digit accuracy
c  z - the point in the complex plane where the moderator function
c        is to be evaluated
c 
c                Output parameters:
c 
c  f - the value of the moderator at the point z
c 
        cw=exp(-ima*c*z)
c 
        call moderfun2(c,delta,z,f,cw)
c 
        delta2=delta/(2.0d0**(done/3))
        call moderfun2(c,delta2,z,f2,cw)
        f=2*f2-f
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine moderfun_data(c,delta,poles)
        implicit real *8 (a-h,o-z)
        save
        complex *16 ima,poles(8)
        dimension deltas(8)
c 
        data ima/(0.0d0,1.0d0)/
c 
c        This subroutine returns to the user 8 poles of the
c        moderating function, corresponding to the parameters
c        c, delta. Please note that the moderating function
c        (evaluated by the subroutine moderfun3 (see)) has an
c        infinite number of poles. The poles are all located in
c        the upper half-plane, in groups of 8. Within each group,
c        the poles sre on a vertical line; such vertical lines
c        are spaced with the interval 2*pi/c. This subroutine
c        returns one such group; for this group, x=pi/c.
c 
c       construct the array of deltas
c 
        done=1
        deltas(1)=delta
        deltas(2)=delta/2
        deltas(3)=delta/sqrt(2.0d0)
        deltas(4)=delta/sqrt(8.0d0)
        deltas(5)=delta/ (2.0d0**(done/3))
        deltas(6)=delta/ (16.0d0**(done/3))
        deltas(7)=delta/ (2.0d0**(5*done/6) )
        deltas(8)=delta/ (2.0d0**(11*done/6))
c 
c        construct the array of corresponding poles
c 
        done=1
        pi=atan(done)*4
        u=pi/c
c 
        do 1200 i=1,8
c 
        v=log(1/deltas(i))/c
        poles(i)=u+ima*v
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine moderfun2(c,delta,z,f,cw)
        implicit real *8 (a-h,o-z)
        save
        complex *16 z,f,f2,cw
c 
        call moderfun1(c,delta,z,f,cw)
  
        delta2=delta/sqrt(2.0d0)
  
        call moderfun1(c,delta2,z,f2,cw)
  
        f=2.0d0*f2-f
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine moderfun1(c,delta,z,f,cw)
        implicit real *8 (a-h,o-z)
        save
        complex *16 z,f,cw
        data done/1.0d0/,half/0.5d0/,c15/1.5d0/
c 
        f=(done+(c15*delta)*cw)/
     1      ((done+(delta*half)*cw)*(done+delta*cw))
c 
        return
        end
