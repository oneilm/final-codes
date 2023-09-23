        IMPLICIT REAL *8 (A-H,O-Z)
        DIMENSION X0(1000),Y0(1000),xs2(1000),ys2(1000),
     1      ww(10 000 000),xsder(1000),ysder(1000)
c 
        CALL PRINI(6,13)
c 
C      INPUT ALL PARAMETERS
C 
         PRINT *, 'ENTER N'
         READ *,N
         CALL PRINF('N=*',N,1 )
c 
c        retrieve the initial polygon from the subroutine
c 
        c=100
cccc        c=40
        sgain=1.1
  
        call polyget(X0,Y0,N0)
c 
        iw=21
        call quaplot2(iw,x0,y0,n0,2,x0,y0,n0,3,
     1      'initial polygon*')
c 
        time1=clotatim()
  
        lenww=10 000 000
        call curvefilt(ier,x0,y0,n0,c,sgain,n,
     1      rlen,xs2,ys2,xsder,ysder,ww,lenww,accur,sgainout)
  
  
        time2=clotatim()
  
  
        call prin2('and time for curvefilt is*',time2-time1,1)
  
  
  
        call prin2('and rlen=*',rlen,1)
        call prin2('and accur=*',accur,1)
        call prin2('and sgainout=*',sgainout,1)
  
        STOP
        END
C 
C 
C 
C 
C 
        SUBROUTINE polyget(X,Y,N)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        DIMENSION X(1),Y(1),izs(2,6)
c 
        data izs/1,1, 3,2, 5,3, 6,4, 7,6, 6,9/
c 
        n=6
        do 1200 i=1,n
c 
        x(i)=izs(1,i)
        y(i)=izs(2,i)
 1200 continue
c 
        return
        end
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c        This is the end of the edbugging code, and the beginning
c        of the resampling code proper
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine curvefilt(ier,x0,y0,n0,c,sgain,n,
     1      rlen,xs2,ys2,xsder,ysder,ww,lenww,accur,
     2      sgainout)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        DIMENSION X0(1),Y0(1),xsder(1),ysder(1),xs2(1),ys2(1),
     1      ww(1)
c 
c       This subroutine filters and resamples the user-supplied
c       curve. The curve is NOT assumed to be closed, and is
c       supplied by the user in the form of a collection of
c       n0 points in the plane, with the coordinates (x0(i),y0(i)).
c       The smoothing process consists of filtering the function
c       (\theta(l), where \theta is the angle between the tangent
c       and the x axis, and l is the arc-length. The filtering is
c       in the Slepian sense, with band-limit c and supergain sgain,
c       both user-supplied.
c 
c       The curve is returned as a collection of nodes (xs2,ys2)
c       tabulating the curve at Legendre nodes in arclength. The
c       subroutine also returns the derivatives of x,t with respect
c       to arclength, in the arrays xsder,ysder
c 
c                     Input parameters:
c 
c  x0,y0 - the coordinates of the points in R^2, defining the curve
c  n0 - the number of elements in each of the arrays x0,y0
c  c - the band-limit to be used in filtering
c  sgain - supergain to be used in filtering
c  n - the number of Gaussian nodes to be returned in the arrays
c       xs2,ys2
c  lenww - the amount of space (in terms of real *8 elements)
c       supplied to the subroutine in the array ww
c 
c                     Output parameters:
c 
c  ier - error return code.
c     ier = 0 means successful conclusion.
c     ier > 0 means that the amount of space provided in the work
c       array ww (see below) is insufficient
c  rlen - the length of the resumpled curve
c  xs2,ys2 - the coordinates of the n Gaussian nodes on the
c       resampled curve
c  xsder,ysder - the derivatives of x,y with respect to arc-length
c       at the points xs2,ys2
c  accur - the accuracy of the resampling.
c     EXPLANATION: accur is the sum of absolute values of the last
c                  5 elements in the Legendre series of x,y as a
c                  function of arclength
c  sgainout - the actual supergain of the angle between the tangent
c     to the curve and the real axis
c 
c 
c       . . . allocate memory
c 
        ixslege=1
        lxslege=n+2
c 
        iwslege=ixslege+lxslege
        lwslege=n+2
c 
        iuu=iwslege+lwslege
        luu=n*n+10
c 
        ivv=iuu+luu
        lvv=n*n+10
c 
        ifis=ivv+lvv
        lfis=n+2
c 
        ifis2=ifis+lfis
        lfis2=n+2
c 
        icoefs=ifis2+lfis2
        lcoefs=n+2
c 
        icoefs2=icoefs+lcoefs
        lcoefs2=n+2
c 
        iww=icoefs2+lcoefs2
c 
        lenww2=lenww-iww
        if(lenww2 .lt. 1000) then
            ier=1024
            return
        endif
c 
        call curvefilt0(jer,x0,y0,n0,c,sgain,n,
     1      ww(ixslege),ww(iwslege),ww(iuu),ww(ivv),
     2      ww(iww),lenww2,xs2,ys2,xsder,ysder,ww(ifis),
     3      ww(ifis2),ww(icoefs),ww(icoefs2),accur,
     4      rlen,sgainout)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine curvefilt0(ier,x0,y0,n0,c,sgain,n,
     1      xslege,wslege,uu,vv,ww,lenww2,xs2,ys2,
     2      xsder,ysder,fis,fis2,coefs,coefs2,accur,
     3      rlen,sgainout)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        DIMENSION X0(1),Y0(1),fis(1),xslege(1),wslege(1),
     1      xsder(1),ysder(1),xs2(1),ys2(1),fis2(1),ww(1),
     2      uu(1),vv(1),coefs(1),coefs2(1)
c 
c       construct the Legendre structure
c 
        ier=0
        itype=2
        call legeexps(itype,n,xslege,uu,vv,wslege)
c 
c        construct the original table of angles
c 
        call crudenodes(x0,y0,n0,n,fis,xslege,wslege,
     1      rlen,xs2,ys2,xsder)
c 
  
        ifinit=1
        call onefilt(x0,y0,rlen,c,xslege,wslege,n,fis,fis2,
     1      ifinit,sgain,ww,sgainout)
c 
c 
        do 1150 i=1,n
c 
        fis(i)=fis2(i)
 1150 continue
c 
c        conduct Newton iterations
c 
  
        ifreco=1
        delta=0.001
        eps=1.0d-10
  
        ifout=0
        do 2600 ijk=1,10
c 
        call oneiter(fis,n,wslege,rlen,x0,y0,n0,
     1        delta,err,ifreco,xslege,sgain,
     2        xs2,ys2,xsder,ysder,ww)
c 
        call prinf('after oneiter, ijk=*',ijk,1)
        call prin2('after oneiter, err=*',err,1)
  
        if(err .lt. eps) ifout=ifout+1
        if(ifout .ge. 2) goto 2800
        ifreco=0
c 
 2600 continue
c 
 2800 continue
c 
c       reconstruct the filtered curve
c 
        call onexys(x0,y0,rlen,uu,vv,n,
     1    fis,xs2,ys2,n0,xsder,ysder,coefs,
     2      coefs2,accur)
c 
        iw=24
        call quaplot2(iw,x0,y0,n0,3,xs2,ys2,n,3,
     1      'after Newtons, resampling vs. initial polygon*')
c 
        iw=25
        call quaplot(iw,xs2,ys2,n,3,
     1      'resampled curve alone*')
c 
c        calculate the second derivatives of x,y with respect
c        to arc length
c 
  
  
        Return
        END
c 
c 
c 
c 
c 
        subroutine onexys(x0,y0,rlen,uu,vv,n,
     1    fis2,xs2,ys2,n0,xsder,ysder,coefs,coefs2,
     2      accur)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        DIMENSION uu(n,n),vv(n,n),xsder(1),ysder(1),
     1      coefs(1),xs2(1),ys2(1),coefs2(1),
     2      fis2(1),x0(1),y0(1)
c 
        do 2600 i=1,n
c 
        xsder(i)=cos(fis2(i))
        ysder(i)=sin(fis2(i))
 2600 continue
c 
c        integrate the tangent vectors, getting the curve
c 
        call matvecr(uu,xsder,coefs,n)
c 
        call legeinte(coefs,n,coefs2)
c 
        accur=abs(coefs2(n))+abs(coefs2(n-1))+
     1      abs(coefs2(n-2))+abs(coefs2(n-3))+abs(coefs2(n-4))
c 
        call matvecr(vv,coefs2,xs2,n)
        call matvecr(uu,ysder,coefs,n)
c 
        call legeinte(coefs,n,coefs2)
c 
        accur=accur+abs(coefs2(n))+abs(coefs2(n-1))+
     1      abs(coefs2(n-2))+abs(coefs2(n-3))+abs(coefs2(n-4))
c 
        call matvecr(vv,coefs2,ys2,n)
c 
        do 2800 i=1,n
c 
        xs2(i)=xs2(i)*rlen/2 + x0(1)
        ys2(i)=ys2(i)*rlen/2 + y0(1)
c 
 2800 continue
  
  
        return
        end
c 
c 
c 
c 
c 
        subroutine onefilt(x0,y0,rlen,c,xslege,wslege,n,
     1    fis,fis2,ifinit,sgain,w,sgainout)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        DIMENSION xslege(1),wslege(1),w(1),
     1      fis2(1),x0(1),y0(1),fis(1)
c 
        complex *16 cder,cval
c 
c 
c        If the user so requested, initialize the prolate
c        filtering routine
c 
        if(ifinit .eq. 0) goto 1600
c 
        lenw=1000 000
        oversamp=1
c 
c       . . . allocate memory for the initialization
c 
        irlams=1
        lrlams=c*4+400
c 
        irlams2=irlams+lrlams
        lrlams2=c*2+200
c 
        ittexp=irlams2+lrlams2
        lttexp=c*10+400
c 
        iw=ittexp+lttexp
c 
        call slepapin(ier,c,xslege,wslege,n,w(iw),lenw,
     1      w(irlams),w(irlams2),nvects,nhigh,oversamp,
     2      w(ittexp),nexpons,keep)
c 
c 
        call prinf('after slepapin, ier=*',ier,1)
        call prinf('after slepapin, keep=*',keep,1)
        call prinf('after slepapin, nhigh=*',nhigh,1)
        call prinf('after slepapin, nvects=*',nvects,1)
        call prin2('after slepapin, w(irlams2)=*',w(irlams2),nvects)
c 
c       . .. collect garbage
c 
        iw2=21
        call prolarrm(w(iw),w(iw2),keep)
c 
        w(1)=iw2+0.1
        w(2)=keep+23
        w(3)=nvects+0.1
        w(4)=nhigh+0.1
        w(5)=nexpons+0.1
        w(6)=nhigh+0.1
c 
 1600 continue
c 
c        construct the smoothed version of the angle
c 
        nvects=w(3)
        nhigh=w(4)
        nexpons=w(5)
        nhigh=w(6)
c 
        iw=w(1)
        lw=w(2)
        ialphas=lw+iw
        lalphas=nvects*2+100
c 
        isrccoe=ialphas+lalphas
        lsrccoe=nhigh*2+100
c 
        ipattcoe=isrccoe+lsrccoe
        lpattcoe=nhigh*2+100
c 
        ipattexp=ipattcoe+lpattcoe
        lpattexp=nexpons*2+100
c 
        icfis=ipattexp+lpattexp+1000
        lcfis=n*2+10
c 
        j=0
        do 2250 i=1,n
c 
        w(icfis+j)=fis(i)
        j=j+1
        w(icfis+j)=0
        j=j+1
 2250 continue
c 
        call slepapev(ier,sgain,w(icfis),w(iw),
     1      w(ialphas),ncoefs,w(isrccoe),w(ipattcoe),nlege,
     2      w(ipattexp),errl2,sgainout)
c 
  
        call prin2('after slepapev, sgainout=*',sgainout,1)
  
        call prin2('after slepapev, errl2=*',errl2,1)
  
        call prinf('in onefilt after slepapev, nlege=*',nlege,1)
  
  
        do 2300 i=1,n
c 
        call legecFDE(xslege(i),cval,cder,w(iPattcoe),Nlege-1)
  
        fis2(i)=cval
 2300 continue
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine oneiter(thetas,n,wslege,rlen,x0,y0,n0,
     1        delta,err,ifreco,xslege,sgain,xs2,ys2,xs,ys,ww)
        implicit real *8 (a-h,o-z)
        save
        dimension thetas(1),x0(1),y0(1),ww(1),
     1      wslege(1),xs(1),ys(1),xs2(1),ys2(1)
c 
c 
c        if the user so requested, construct the corrections
c        to be added (after scaling) to the angle
c 
        if(ifreco .eq. 0) goto 1400
c 
        do 1200 i=1,n
        xs(i)=dcos(thetas(i))
        ys(i)=dsin(thetas(i))
 1200 continue
c 
        ifinit=0
c 
        call onefilt(x0,y0,rlen,c,xslege,wslege,n,
     1    xs,xs2,ifinit,sgain,ww,sgainout)
c 
        call onefilt(x0,y0,rlen,c,xslege,wslege,n,
     1    ys,ys2,ifinit,sgain,ww,sgainout)
c 
 1400 continue
c 
c        . . . calculate the discrepancy of the user-supplied
c              tangent angle
c 
        call comdev(thetas,n,wslege,rlen,x0,y0,cosint,sinint)
c 
        cosint=cosint-x0(n0)
        sinint=sinint-y0(n0)
c 
        err=cosint**2+sinint**2
        err=sqrt(err)
  
         call prin2('in oneiter, cosint=*',cosint,1)
         call prin2('in oneiter, sinint=*',sinint,1)
c 
c       obtain the matrix connecting the integrals with the
c       perturbations
c 
c 
c        . . . construct the perturbed tangent vectors
c              and the corresponding discrepancies
c 
        do 1600 i=1,n
        xs(i)=thetas(i)+delta*xs2(i)
 1600 continue
c 
        call comdev(xs,n,wslege,rlen,x0,y0,a11p,a12p)
c 
        do 1800 i=1,n
        xs(i)=thetas(i)-delta*xs2(i)
 1800 continue
c 
        call comdev(xs,n,wslege,rlen,x0,y0,a11m,a12m)
c 
        do 2000 i=1,n
        xs(i)=thetas(i)+delta*ys2(i)
 2000 continue
c 
        call comdev(xs,n,wslege,rlen,x0,y0,a21p,a22p)
c 
        do 2200 i=1,n
        xs(i)=thetas(i)-delta*ys2(i)
 2200 continue
c 
        call comdev(xs,n,wslege,rlen,x0,y0,a21m,a22m)
c 
c         calculate the matrix for the newton
c         via finite differences
c 
        a11=(a11p-a11m)/(2*delta)
        a21=(a12p-a12m)/(2*delta)
        a12=(a21p-a21m)/(2*delta)
        a22=(a22p-a22m)/(2*delta)
c 
c        solve the linear system
c 
        y1=cosint
        y2=sinint
c 
        det=a22*a11-a12*a21
        det1=y1*a22-y2*a12
        det2=y2*a11-y1*a21
c 
        delx=det1/det
        dely=det2/det
c 
c       adjust thetas
c 
        do 2400 i=1,n
        thetas(i)=thetas(i)-(delx*xs2(i)+dely*ys2(i))
 2400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine comdev(theta,n,wslege,rlen,x0,y0,dx,dy)
        implicit real *8 (a-h,o-z)
        save
        dimension theta(1),wslege(1),x0(1),y0(1)
c 
        dx=0
        dy=0
        do 1200 i=1,n
        ddx=dcos(theta(i))
        ddy=dsin(theta(i))
        dx=dx+ddx*wslege(i)
        dy=dy+ddy*wslege(i)
 1200 continue
c 
        dx=dx*rlen/2+x0(1)
        dy=dy*rlen/2+y0(1)
        return
        end
c 
c 
c 
c 
c 
        subroutine matvecr(a,x,y,n)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        DIMENSION a(n,n)
        real *8 x(1),y(1),cd
c 
        do 2000 i=1,n
c 
        cd=0
        do 1800 j=1,n
c 
        cd=cd+a(i,j)*x(j)
 1800 continue
c 
        y(i)=cd
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine matvec(a,x,y,n)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        DIMENSION a(n,n)
        complex *16 x(1),y(1),cd
c 
        do 2000 i=1,n
c 
        cd=0
        do 1800 j=1,n
c 
        cd=cd+a(i,j)*x(j)
 1800 continue
c 
        y(i)=cd
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine crudenodes(x0,y0,n0,n,fis,xslege,wslege,
     1      rlen,thetas,s,ss)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        DIMENSION X0(1),Y0(1),fis(1),thetas(1),s(1),
     1      xslege(1),wslege(1),ss(1)
C 
        nn=n
c 
c        determine the total length of the user-supplied polygon
c 
        rlen=0
        do 1200 i=1,n0-1
c 
        d=(x0(i+1)-x0(i))**2+(y0(i+1)-y0(i))**2
        rlen=rlen+sqrt(d)
 1200 continue
c 
c        construct the linear mapping, converting the interval
c              [-1,1] into the interval [0,rlen]
c 
        alpha=rlen/2
        beta=rlen/2
c 
        do 1300 i=1,nn
c 
        ss(i)=alpha*xslege(i)+beta
 1300 continue
c 
c        construct the resampling of the tangent angle
c        at Gaussian nodes
c 
        call REPOLY(IER,X0,Y0,N0,S,nn,H,fis,thetas,ss)
c 
        return
        end
C 
C 
C 
C 
C 
        SUBROUTINE REPOLY(IER,X,Y,N,S,M,H,fis,thetas,ss)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        DIMENSION X(1),Y(1),S(1),fis(1),thetas(1),ss(1)
C 
C       THIS SUBROUTINE FOR A USER-SPECIFIED POLYGONAL CURVE
C       (X,Y) CONSTRUCTS THE ARRAY S OF ITS CUMULATIVE ARC-LENGTHS,
c       snd the array fis  of its tangent angles
C       THE USER-PROVIDED POLYGON MUST BE CLOSED, I.E.
C       X(N)=X(1),
C       Y(N)=Y(1).
c       note also that m must be greater than n
C 
C                  INPUT PARAMETERS:
C 
C   X,Y - COORDINATS OF THE USER-SPECIFIED POLYGON.
C   N - THE NUMBER OF ELEMENTS IN ARRAYS X,Y
C   M - THE NUMBER OF NODES IN THE EQUISPACED RESAMPLING
C       OF THE CURVE
C 
C                  OUTPUT PARAMETERS:
C 
C   IER - ERROR RETURN CODE.
C      IER=0 MEANS SUCCESSFULL EXECUTION
C      IER=16 MEANS THAT THE USER-SUPPLIED POLYGON IS NOT
C         CLOSED. THIS IS A TERMINAL ERROR.
C   S - THE ARRAY OF CUMULATIVE POLYGONAL ARC-LENGTHS ON THE
C       USER-SUPPLIED POLYGON (NOTE THAT S(1)=0).
C   H - THE SAMPLING INTERVAL ON THE EQUISPACED-RESAMPLED
C       CURVE
c   fis - the tangent angles of the resampled curve
c 
c                  work arrays:
c 
c   thetas - must be at least m+2 real *8 locations long
c 
        IER=0
        EPS0=1.0D-13
cccc        CALL PRIN2('IN REPOLY, X=*',X,N)
cccc        CALL PRIN2('IN REPOLY, Y=*',Y,N)
CCC     CALL PRINF('IN REPOLY, N=*',N,1)
CCC     CALL PRINF('IN REPOLY, M=*',M,1)
C 
C       CONSTRUCT THE ARRAY OF USER-PROVIDED LENGTHS FOR ALL
C       NODES OF THE USER-PROVIDED POLYGON
C 
        S(1)=0
        DO 1200 I=2,N
        D= (X(I)-X(I-1))**2 + (Y(I)-Y(I-1))**2
        S(I)=S(I-1)+DSQRT(D)
 1200 CONTINUE
cccc         call prin2('inside repoly, s=*',s,n)
C 
        H=S(N)/(M-1)
ccc        H=S(N)/m
c 
c        construct the angles between the x axis and the
c        various sides of the polygon
c 
c       . . . construct the angles between the consecutive
c             sides of the polygon
c 
        call REANCR(X,Y,N,FIS)
ccc         call prin2('angles between sides of polygon are*',fis,n)
c 
c       find the angle between the first side of the polygon and
c       the x axis
c 
        wx0=x(1)-1
        wy0=y(1)
c 
        call REANGL(wX0,wY0,X(1),Y(1),X(2),Y(2),thetas(1))
c 
c        construct the angles
c 
         do 2200 i=2,n
         thetas(i)=thetas(i-1)+fis(i)
 2200 continue
ccc         call prin2('tangent angles of original polygon are*',
ccc     1    thetas,n)
c 
c       construct the angles between the x axis and the various
c       sides of the polygon, but tabulated at the nodes of
c       the finer discretization
c 
        j=1
ccc         call prin2('in repoly, h=*',h,1)
        do 2400 i=1,m
        if(ss(i) .gt.s(j+1)) j=j+1
        fis(i)=thetas(j)
 2400 continue
c 
        RETURN
        END
C 
C 
C 
C 
C 
        SUBROUTINE REANCR(X,Y,N,FIS)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        DIMENSION X(1),Y(1),FIS(1)
C 
C       FOR THE USER-SPECIFIED POLYGON, THIS SUBROUTINE
C       CONSTRUCTS THE ANGLES MADE BY THE CONSECUTIVE SIDES
C       OF THE POLYGON. THE POLYGON IS ASSUMED TO BE CLOSED.
C 
C                  INPUT PARAMETERS:
C 
C   X,Y - THE COORDINATES OF THE VERTICES OF THE USER-SUPPLIED
C         POLYGON
C   N - THE NUMBER OF ELEMENTS IN ARRAYS X,Y
C 
C                  OUTPUT PARAMETERS:
C 
C   FIS - THE ANGLES MADE BY CONSECUTIVE SIDES OF THE POLYGON
C 
C      REMARK:
C       IT IS ASSUMED THAT THE POLYGON IS CLOSED, I.E.
C       X(N)=X(1),
C       Y(N)=Y(1).
C 
C       . . . INITIALIZE THE ANGLE-COMPUTING ROUTINE AND
C             CALCULATE THE ANGLE BETWEEN THE LAST AND THE FIRST
C             SIDES. THIS ANGLE IS DECLARED TO BE AT THE
C             FIRST VERTEX.
C 
        CALL REANGI(IER)
        CALL REANGL(X(N-1),Y(N-1),X(N),Y(N),X(2),Y(2),FIS(1))
C 
C       CALCULATE THE REST OF THE ANGLES
C 
        DO 1200 I=2,N-1
        CALL REANGL(X(I-1),Y(I-1),X(I),Y(I),X(I+1),Y(I+1),FIS(I))
CCCC    CALL PRIN2('AFTER REANGL, FIS(I)  IS*',FIS(I),1)
 1200 CONTINUE
        FIS(N)=FIS(1)
        RETURN
        END
C 
C 
C 
C 
C 
        SUBROUTINE REANGL(X1,Y1,X2,Y2,X3,Y3,FI)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        DATA EPS/1.0D-30/
C 
C 
C       THIS SUBROUTINE FINDS THE ANGLE FI BETWEEN THE
C       LINE SEGMENTS  [(X1,Y1),(X2,Y2)] AND  [(X2,Y2),(X3,Y3)]
C       INITIALIZATION ENTRY POINT REANGI (SEE BELOW) HAS TO BE
C       CALLED PRIOR TO THE INVOCATION OF THIS MAIN ENTRY
C 
C        . . . FIND THE (UNSCALED) SIN AND COS OF THE ANGLE
C 
        UX=X2-X1
        UY=Y2-Y1
C 
        VX=X3-X2
        VY=Y3-Y2
C 
        SINFI=UX*VY-VX*UY
        COSFI=UX*VX+UY*VY
C 
C       IF THE ANGLE IS NOT RIGHT - FIND IT AND EXIT
C 
          IF(DABS(COSFI) .LT. EPS) GOTO 2000
        FI=DATAN(SINFI/COSFI)
CCCC    IF(COSFI .LE. 0) FI=FI+PI
        IF( (COSFI .LE. 0) .AND. (SINFI.GE.0) ) FI=FI+PI
        IF( (COSFI .LE. 0) .AND. (SINFI.LE.0) ) FI=FI-PI
        RETURN
 2000 CONTINUE
C 
C       THE ANGLE IS RIGHT. FIND IT
C 
CCCC    CALL PRIN2('IN REANGL, THE ANGLE IS RIGHT. COSFI IS*',COSFI,1)
        FI=PI2
        IF(SINFI .LT. 0.0D0) FI=-FI
CCC     IF(SINFI .GT. 0.0D0) FI=-FI
        RETURN
C 
C 
C 
C 
        ENTRY REANGI(IER)
C 
C       THIS IS INITIALIZATION ENTRY POINT. IT HAS TO BE CALLED
C       PRIOR TO ANY CALLS TO THE MAIN ENTRY REANGL (ABOVE).
c 
        DONE=1
        PI=DATAN(DONE)*4
        PI2=PI/2
        RETURN
        END
