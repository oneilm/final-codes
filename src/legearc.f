c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       this is the end of the debugging code
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       the user-callable routines in this file which are contained
c       below are as follows:
c
c       legearcnodes - 
c
c       zlegearcnodes -
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
        subroutine zlegearcnodes(ier,n,czs,k,zsout,zwhtsout,
     1      tsout,rlout)
        implicit real *8 (a-h,o-z)
        real *8 tsout(1),cxs(0:10000),cys(0:10000),
     1      xsout(10000),ysout(10000),ws(10000)
            complex *16 czs(0:1),zsout(1),zwhtsout(1),ima,z1,z2,
     1      zsc
c
c       this routine is the same as legearcnodes, except meant for
c       curves parameterized in the complex plane (i.e. have one
c       complex legendre expansion instead of two real)
c
c       input:
c
c         n - the order of the legendre expansions describing the curve
c         czs - complex legendre expansion for z, indexed 0 to n,
c             assumed to be parameterized on the interval [-1,1] as usual
c         k - order of quadrature to construct on the curve
c
c       output:
c
c         ier - error code, 0 is good, 1 need more nodes for arclength,
c             16 newton did not converge, 16 is fatal
c         zsout - COMPLEX legendre nodes on the curve
c         zwhtsout - legendre weights (scaled by rlout), note that these
c             are COMPLEX, and depend on the endpoints of the curve
c         tsout - the values of t which yield legendre nodes in arclength
c         rlout - the total length of the curve
c
c
        done=1
        pi=4*atan(done)
        ima=(0,1)

c
c       convert to real and call legearcnodes
c
        do 1600 i=0,n
        cxs(i)=czs(i)
        cys(i)=-ima*czs(i)
 1600 continue
c
        call legearcnodes(ier,n,cxs,cys,k,xsout,ysout,ws,
     1      tsout,rlout)

c
c       calculate end points to scale the weights with
c
        t=-1
        call legeexev(t,x,cxs,n)
        call legeexev(t,y,cys,n)
        z1=x+ima*y
c
        t=1
        call legeexev(t,x,cxs,n)
        call legeexev(t,y,cys,n)
        z2=x+ima*y
c
        zsc=(z2-z1)/rlout
c
        do 1800 i=1,k
        zsout(i)=xsout(i)+ima*ysout(i)
        zwhtsout(i)=zsc*ws(i)
 1800 continue
c
        return
        end
c
c
c
c
c
        subroutine legearcnodes(ier,n,cxs,cys,k,xsout,ysout,whtsout,
     1      tsout,rlout)
        implicit real *8 (a-h,o-z)
        real *8 cxs(0:1),cys(0:1),xsout(1),ysout(1),tsout(1),
     1      whtsout(1)
c
c       this routine receives a parameterized 2d curve given by its
c       length n legendre expansions and returns k legendre nodes
c       along the curve in arclength
c
c       input:
c
c         n - the order of the legendre expansions describing the curve
c         cxs,cys - legendre expansion for x and y, indexed 0 to n,
c             assumed to be parameterized on the interval [-1,1] as usual
c         k - order of quadrature to construct on the curve
c
c       output:
c
c         ier - error code, 0 is good, 1 need more nodes for arclength,
c             16 newton did not converge, 16 is fatal
c         xsout,ysout - legendre nodes on the curve
c         whtsout - legendre weights (scaled by rlout)
c         tsout - the values of t which yield legendre nodes in arclength
c         rlout - the total length of the curve
c
c
        done=1
        pi=atan(done)
        ier=0

c
c       calculate the total arclength of the curve
c
        t=1
        call onearc(ier7,n,cxs,cys,t,rlout,dsdt)
c
        if (ier7 .ne. 0) then
          ier=1
          call prinf('in legearcnodes, ier=*',ier,1)
        endif
c
c       create the k nodes on [0,rlout]
c
        ifwhts=1
        call legewhts(k,tsout,whtsout,ifwhts)
c
        do 2000 i=1,k
        tsout(i)=(tsout(i)+1)/2*rlout
        whtsout(i)=whtsout(i)/2*rlout
 2000 continue

c
c       now solve for the t such that s(i)=tsout(i) using newton
c
        do 3000 j=1,k
c
        t1=tsout(j)/rlout*2-1
        thresh=1.0d-8
c
        do 2400 i=1,30
        t0=t1
        call onearc(ier7,n,cxs,cys,t0,s0,dsdt0)
        t1=t0-(s0-tsout(j))/dsdt0
        d=t1-t0
        if (abs(d) .le. thresh) goto 2600
 2400 continue
c
        ier=16
        call prinf('bomb! newton was not close! ier=*',ier,0)
        return
     

 2600 continue

c
c       do two more newton iterations for good measure
c
        t0=t1
        call onearc(ier,n,cxs,cys,t0,s0,dsdt0)
        t1=t0-(s0-tsout(j))/dsdt0
c
        t0=t1
        call onearc(ier,n,cxs,cys,t0,s0,dsdt0)
        t1=t0-(s0-tsout(j))/dsdt0
c
        d=abs(t1-t0)/abs(t1)
        thresh=1.0d-12
        if (d .gt. thresh) then
          ier=16
          call prinf('original iterations=*',i,1)
          call prinf('point=*',j,1)
          call prin2('error in newton=*',d,1)
          call prinf('bomb! newton did not converge! ier=*',ier,0)
          stop
          return
        endif
c
        tsout(j)=t1
c
 3000 continue

c
c       now evaluate the curve at the points
c
        do 3400 i=1,k
        call legeexev(tsout(i),xsout(i),cxs,n)
        call legeexev(tsout(i),ysout(i),cys,n)
 3400 continue
c
        return
        end
c
c
c
c
c
        subroutine onearc(ier,n,cxs,cys,t,s,dsdt)
        implicit real *8 (a-h,o-z)
        real *8 cxs(0:1),cys(0:1),ts(10 000),whts(10 000)
c
c       this routine calculates the arclength from -1 to t of
c       the curve described by legendre coefficients cxs,cys
c
c       input:
c
c         n - the order of the legendre expansions describing the curve
c         cxs,cys - legendre expansion for x and y, indexed 0 to n
c         t - upper integration limit
c
c       output:
c
c         s - the arclength from -1 to t
c         dsdt - the derivative of arclength with respect to the
c             parameterization t at the point t
c
c
        done=1
        pi=4*atan(done)
        ier=0
c
        kk=n+10
        u=t+1
c
        ifwhts=1
        call legewhts(kk,ts,whts,ifwhts)
c
        do 1200 i=1,kk
        ts(i)=(ts(i)+1)/2*u-1
        whts(i)=whts(i)/2*u
 1200 continue
c
        cd=0
        do 1600 i=1,kk
        call legefder(ts(i),valx,dvalx,cxs,n)
        call legefder(ts(i),valy,dvaly,cys,n)
        dd=sqrt(dvalx**2+dvaly**2)
        cd=cd+dd*whts(i)
 1600 continue
c
        s1=cd

c
c       calculate again with 15 more nodes to check
c
        kk=kk+15
c
        ifwhts=1
        call legewhts(kk,ts,whts,ifwhts)
c
        do 2200 i=1,kk
        ts(i)=(ts(i)+1)/2*u-1
        whts(i)=whts(i)/2*u
 2200 continue
c
        cd=0
        do 2600 i=1,kk
        call legefder(ts(i),valx,dvalx,cxs,n)
        call legefder(ts(i),valy,dvaly,cys,n)
        dd=sqrt(dvalx**2+dvaly**2)
        cd=cd+dd*whts(i)
 2600 continue
c
        s2=cd
c
        d=abs(s2-s1)/abs(s2)
        thresh=1.0d-14
        if (d .gt. thresh) then
          call prin2('arclength only calculated to rel prec*',d,1)
          ier=1
        endif
c
        s=s2

c
c       and calculate the derivative at t
c
        call legefder(t,valx,dvalx,cxs,n)
        call legefder(t,valy,dvaly,cys,n)
        dsdt=sqrt(dvalx**2+dvaly**2)
c
        return
        end
