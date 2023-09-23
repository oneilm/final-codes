        implicit real *8 (a-h,o-z)
        dimension xsides(1000),ysides(1000),
     1      wholemat(1000 000),xsnorms(1000),ysnorms(1000),
     2      whts3(1000),corners(2,500)
        dimension ww(4000 000)
c 
        external chnkpnt,interact2
c 
cccc        data corners/ 0,0, 2,1, 3,1, 3,3, 0,3/
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER k'
        READ *,k
        CALL PRINf('k=*',k,1 )
  
c 
c 
  
  
  
  
cccc        iw=22
cccc        n=3000
cccc        call o1rauwrt(iw,n,ww)
  
  
cccc        iw=23
cccc        call o3rwrall(iw,ww)
  
  
cccc        stop
c 
c        create the test curve
c 
        ncorners=5
  
        done=1
        pi=atan(done)*4
        phi=pi/40
        phi=pi/74
        phi=pi/10
  
        call corners4(phi,corners,ncorners)
  
        call boat36(corners,ncorners)
        call boat8(corners,ncorners)
  
        call gunboat2(corners,ncorners)
cccc        call gunboat(corners,ncorners)
  
        call prin2('corners as created*',corners,ncorners*2)
  
        call chnkinit(corners,ncorners)
  
c 
c       construct the matrix of interactions and solve the resulting
c       linear system
c 
        n=ncorners*20
c 
c 
        xsource=10
        ysource=30
c 
        xreceiv=1.5
        yreceiv=2
  
  
        xreceiv=2.99
        yreceiv=2
  
  
  
        xreceiv=0.7
        yreceiv=0.00
  
cccc        xreceiv=1.5
cccc        yreceiv=2
        xreceiv=0.7
        yreceiv=0.00
  
        xreceiv=10
        yreceiv=5
  
  
  
        ir1=23
  
  
        call wholmatr(ir1,xsides,ysides,xsnorms,ysnorms,
     1      ncorners,whts3,wholemat,
     2      chnkpnt,interact2)
  
  
        call prinf('after wholmatr, ncorners=*',ncorners,1)
        npts=ncorners*20
  
  
        iw=32
        call lotaplot(iw,xsides,ysides,npts,
     1      'polygon as constructed*')
  
  
        call prin2('after wholmatr, xsnorms=*',xsnorms,npts)
        call prin2('after wholmatr, ysnorms=*',ysnorms,npts)
  
  
        ir=22
  
        call dirsol2(ir,xsides,ysides,xsnorms,ysnorms,
     1      npts,whts3,wholemat,xsource,ysource,xreceiv,yreceiv,
     2      interact2)
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine dirsol2(ir,xsides,ysides,xsnorms,ysnorms,
     1      n,whts3,wholemat,xsource,ysource,xreceiv,yreceiv,
     2      interact)
        implicit real *8 (a-h,o-z)
        save
        dimension wholemat(n,n),xsides(1),ysides(1),
     1     xsnorms(1),ysnorms(1),whts3(1),work(100 000),
     2     rhs(1000),sol(1000)
  
        dimension uu(1000 000),vv(1000 000),ss(1000),
     1      ww(4000 000)
c 
        external interact,chnkpnt
c 
c        add the diagonal terms to the matrix of interactions
c 
        done=1
        pi=atan(done)*4
        do 2400 i=1,n
c 
        wholemat(i,i)=wholemat(i,i)+pi
cccc        wholemat(i,i)=wholemat(i,i)-pi
 2400 continue
  
        if(2 .ne. 3) goto 2500
  
        eps=1.0d-12
        lww=4000 000
        call svdpivot(ier,wholemat,n,n,uu,vv,ss,ncols,eps,
     1      ww,lww,ltot)
  
        call prin2('after svdpivot, ss=*',ss,n)
  
          stop
  
  
 2500 continue
c 
c       invert the matrix of interactions
c 
        call orthom(wholemat,n,work,cond)
  
        call prin2('after orthom, cond=*',cond,1)
        call prinf('after orthom, n=*',n,1)
c 
c        evaluate the right-hand side
c 
        do 1600 i=1,n
c 
        d=(xsides(i)-xsource)**2+(ysides(i)-ysource)**2
c 
        rhs(i)=log(d)
 1600 continue
  
        call prin2('rhs as created=*',rhs,n)
  
        call matvec8(wholemat,rhs,sol,n)
  
        call prin2('and sol=*',sol,n)
  
c 
c       evaluate the potential at the receiver via the
c       solution of the boundary value problem
c 
        call soleval(xsides,ysides,xsnorms,ysnorms,
     1    whts3,sol,n,interact,xreceiv,yreceiv,f1)
  
  
        call prin2('and f1=*',f1,1)
  
  
        dd=(xsource-xreceiv)**2+(ysource-yreceiv)**2
  
        f2=log(dd)
        call prin2('and f2=*',f2,1)
        call prin2('and f2-f1=*',f2-f1,1)
        call prin2('and f2/f1=*',f2/f1,1)
c 
c       evaluate the field via a grotesquely oversamppled
c       discretization of the boundary
c 
        call resaeva(ir,chnkpnt,interact,sol,n,
     1     xreceiv,yreceiv,f3)
  
  
        call prin2('and after resaeva, f3=*',f3,1)
        call prin2('and f3/f2=*',f3/f2,1)
  
        call prin2('and f3-f2=*',f3-f2,1)
  
  
        return
        end
  
c 
c 
c 
c 
c 
  
        subroutine chnkpnt(ichunk,t,xy,dxydt)
        implicit real *8 (a-h,o-z)
        save
        dimension xy(2),dxydt(2),pts(2,10 000),corners(2,1),
     1      xalphas(10 000),xbetas(10 000),yalphas(10 000),
     2      ybetas(10 000)
c 
c 
c 
c        This subroutine returns to the user a point in a "standard"
c        parametrization of a polygon. The polygon must have been
c        supplied to the subroutine by the user via a preceding call
c        to the entry chnkinit (see below) of this subroutine.
c 
        i=ichunk
        xy(1)=xalphas(i)*t+xbetas(i)
        xy(2)=yalphas(i)*t+ybetas(i)
c 
        dxydt(1)=xalphas(i)
        dxydt(2)=yalphas(i)
c 
        return
c 
c 
c 
c 
        entry chnkinit(corners,ncorners)
c 
        do 2200 i=1,ncorners
c 
        pts(1,i)=corners(1,i)
        pts(2,i)=corners(2,i)
 2200 continue
c 
        ncp1=ncorners+1
        pts(1,ncp1)=corners(1,1)
        pts(2,ncp1)=corners(2,1)
c 
c       construct the linear mappings parametrizing the sides
c       of the polygon by the interval [-1,1]
c 
        do 2400 i=1,ncorners
c 
        xalphas(i)=(pts(1,i+1)-pts(1,i))/2
        xbetas(i)=(pts(1,i+1)+pts(1,i))/2
c 
        yalphas(i)=(pts(2,i+1)-pts(2,i))/2
        ybetas(i)=(pts(2,i+1)+pts(2,i))/2
c 
 2400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine interact2(xy,xynorm,targ,targnorm,fout)
        implicit real *8 (a-h,o-z)
        save
        dimension xy(2),xynorm(2),targ(2),targnorm(2)
c 
        dx=targ(1)-xy(1)
        dy=targ(2)-xy(2)
c 
        dd=dx*xynorm(1)+dy*xynorm(2)
  
        fout=dd/(dx**2+dy**2)
  
  
cccc        return
  
  
        fout=fout+0.1
  
        return
  
        ddd=(dx**2+dy**2)
        d=( (targ(1)-xy(1))**2+(targ(2)-xy(2))**2)
  
        fout=fout+log(d)+xy(1)
c 
cccc        fout=fout+10*log(d)+xy(1)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine corners4(phi,corners,ncorners)
        implicit real *8 (a-h,o-z)
        save
        dimension corners(2,1)
c 
c        construct the parallelogram
c 
        ncorners=4
c 
        corners(1,1)=-1
        corners(2,1)=0
c 
        corners(1,2)=0
        corners(2,2)=-tan(phi)
c 
        corners(1,3)=1
        corners(2,3)=0
c 
        corners(1,4)=0
        corners(2,4)=tan(phi)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine boat18(corners,ncorners)
        implicit real *8 (a-h,o-z)
        save
        dimension corners0(2,13),corners(2,13)
c 
        data corners0/10,3, 15,3, 17,4.5, 19,6, 16,5, 12,5,
     1       12,7, 7,7, 7,5, 4,5, 1,6, 3,4.5, 5,3/
c 
        ncorners=13
        do 1200 i=1,13
c 
        corners(1,i)=corners0(1,i)
        corners(2,i)=corners0(2,i)
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine boat8(corners,ncorners)
        implicit real *8 (a-h,o-z)
        save
        dimension corners0(2,13),corners(2,13)
c 
        data corners0/10,4, 15,4, 17,6, 20,9, 16,6, 12,6,
     1       12,8, 8,8, 8,6, 4,6, 0,9, 3,6, 5,4/
c 
        ncorners=13
        do 1200 i=1,13
c 
        corners(1,i)=corners0(1,i)
        corners(2,i)=corners0(2,i)
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine boat36(corners,ncorners)
        implicit real *8 (a-h,o-z)
        save
        dimension corners0(2,9),corners(2,9)
c 
        data corners0/10,3, 16,3, 20,6, 13,6, 13,8, 7,8,
     1       7,6, 0,6, 4,3/
c 
        ncorners=9
        do 1200 i=1,9
c 
        corners(1,i)=corners0(1,i)
        corners(2,i)=corners0(2,i)
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine gunboat(corners,ncorners)
        implicit real *8 (a-h,o-z)
        save
        dimension corners0(2,17),corners(2,17)
c 
        data corners0/10,3, 16,3, 20,6, 13,6, 13,7, 17,8,
     1       21,9, 25,10, 25,11, 21,10, 17,9, 13,8, 13,9,
     2       7,9, 7,6, 0,6, 4,3/
c 
        ncorners=17
        do 1200 i=1,17
c 
        corners(1,i)=corners0(1,i)
        corners(2,i)=corners0(2,i)
 1200 continue
  
  
        corners(1,9)=corners(1,9)-0.4
        corners(2,9)=corners(2,9)-0.1
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine gunboat2(corners,ncorners)
        implicit real *8 (a-h,o-z)
        save
        dimension corners0(2,26),corners(2,17)
c 
        data corners0/12,3, 19,3, 22,6, 15,6, 15,7, 22,7.5,
     1       29,8, 29,9, 22,8.5, 15,8, 15,9, 9,9, 9,6, 7,6,
     2       7,7, 5,7, 5,6.75, 4,6.875, 3,7, 3,6.5, 4,6.375,
     3       5,6.25, 5,6, 4,6, 2,6, 5,3/
c 
        ncorners=26
        do 1200 i=1,26
c 
        corners(1,i)=corners0(1,i)
        corners(2,i)=corners0(2,i)
 1200 continue
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine matavec8(amatr,nn,nlarge,x,y)
        implicit real *8 (a-h,o-z)
        save
        real *8 x(1),y(1),amatr(nn,nlarge)
c 
        do 1400 i=1,nlarge
c 
        d=0
        do 1200 j=1,nn
c 
        d=d+amatr(j,i)*x(j)
 1200 continue
c 
        y(i)=d
 1400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine matvec8(a,x,y,n)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),x(n),y(n)
c 
        do 1400 i=1,n
        d=0
        do 1200 j=1,n
        d=d+a(i,j)*x(j)
 1200 continue
        y(i)=d
 1400 continue
        return
        end
c 
c 
c 
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c    This is the end of the debugging code and the beginning of the
c    matrix-construction code proper
c 
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c       This file contains machinery for the construction of
c       interactions of polygonal boundaries with themselves. The
c       only non-trivial limitation on the polygon is that no angle
c       can be less that 30 degrees. The file contains 5 user-callable
c       subroutines: o2rrd2, o2rwrt2, o2rmatb2, o2rappb2, o2rvecne,
c       o2rsubm, o1rmtot. Following is a description of these
c       subroutines:
c 
c 
c   wholmatr - constructs the matrix of interactions corresponding to
c       the user-specified boundary and interaction law.
c       Please note that the subroutine does NOT add in the diagonal
c       part of the matrix (the part that makes the equation a second
c       order one); the reason for this is that the diagonal part is
c       very much problem-specific: interior vs. exterior, etc.
c 
c   o2rsubm - constructs the matrix of interactions of a single 20-point
c       chunk with another 20-point chunk.
c 
c   o1rmtot - constructs the matrix of interactions on a single chunk
c       discretized by 20 points.
c 
c   o3rvecne - returns to the user the coefficients of the linear form
c       used for the evaluation of the potential created at the
c       user-supplied point targ by a "charge" (dipole, etc.)
c       distribution on the user-specified chunk in R^2  (with targ
c       OUTSIDE the chunk on which the charge is distributed)
c 
c   soleval - determines the field generated at a point by a charge
c       (dipole, etc) distriburion on a surface.
c 
c   resaeva - resamples a user-specified function (conceptually, a
c       solution of a boundary integral equation) to a "universal"
c       oversampled grid, and uses the resampled version to evaluate
c       the resulting field at a point.
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        subroutine wholmatr(ir,xsides,ysides,xsnorms,ysnorms,
     1      nsides,weights,wholemat,chnkpnt,interact)
        implicit real *8 (a-h,o-z)
        save
        dimension wholemat(nsides*20,nsides*20),xsides(1),ysides(1),
     1     xsnorms(1),ysnorms(1),weights(1),ts(21),whts(21)
c 
        dimension xs7(21),whts7(21),totmat(20,20),
     1      xys(21),whtsxys(21),xy(2),dxydt(2)
c 
        external interact
c 
        data ifcalled/0/
c 
c       This subroutine constructs the matrix of interactions
c       corresponding to the user-specified boundary and interaction
c       law. The boundary is supplied via the subroutine chnkpnt, and
c       the interaction law is supplied via the subroutine interact.
c       Please note that the subroutine does not add in the diagonal
c       part of the matrix (the part that makes trhe equation a second
c       order one); the reason for this is that the diagonal part is
c       very much problem-specific: interior vs. exterior, etc.
c 
c                      Input parameters:
c 
c  ir - the FORTRAN file number from which the various precomputed data
c       are to be read. Hopefully, these data have been stored there by
c       a prior call to the subroutine o3rwrall (see)
c 
c  nsides - the number of sides on the polygonal boundary
cccc---------------------------------------------------------------------
c 
c  chnkpnt - the subroutine providing the parametrizations of the
c        polygonal boundary on which the matrix is being created.
c        Conceptually,
c        chnkpnt implements nside mappings of the interval [-1,1] to R^2;
c        the output parameters xy, dxydt are the image of the point
c        t in R^2 and the derivative of the image with respect to
c        t. In other words, dxydt is a tangent vector to the curve
c        at the point xy, but the length of dxydt is not equal to 1
c        (generally speaking).
c 
c        The calling sequence of chnkpnt is
c 
c         chnkpnt(ichunk,t,xy,dxydt)
c 
c          Input parameters of chnkpnt:
c 
c  ichunk - the sequence number of the side on which the point lives
c  t - the point on the interval [-1,1] parametrizing the side
c        number ichunk
c 
c          Output parameters of chnkpnt:
c 
c  xy - the image in R^2 of the point t on the interval [-1,1]
c  dxydt - the derivative of xy with respect to t
c 
cccc-------------------------------------------------------------------
c 
cccc-------------------------------------------------------------------
c 
c  interact - the subroutine evaluating the interaction between two
c        points on the interval [-1,1]. The calling sequence of
c        interact is
c 
c         interact(xy,xynorm,targ,targnorm,fout)
c 
c                 Input parameters of interact:
c 
c  xy - the source point in R^2
c  xynorm - the (unit) normal at the source point
c  targ - the target point in R^2
c  targnorm - the (unit) notmal at the target
c 
c                 Output parameters of interact:
c 
c  fout - the potential (normal derivative, whatever) created at
c        the target by the source
c 
cccc-------------------------------------------------------------------
c 
c 
c                      Output parameters:
c 
c  xsides,ysides - the coordinates in R^2 of the nodes of the
c        obtained discretization of the boundary (nsides*20 of
c        them things)
c  xsnorms,ysnorms - the normals to the boundary at the nodes
c        (xsides,ysides)
c  weights - the weights of the quadrature formula based on
c        the nodse (xsides,ysides)
c  wholemat - the complete matrix of interactions (except for
c        diagonal terms); its dimensionality is (nsides*20) \times
c        (nsides*20)
c 
c       . . . construct the 20 support nodes on the interval [-1,1]
c 
        if(ifcalled .eq. 0) call o1rnodes(ts,whts,n20)
        ifcalled=1
  
        call prinf('in wholmatr, nsides=*',nsides,1)
  
c 
c       construct the discretization of the user-specified object
c 
        ii=0
        do 1600 i=1,nsides
c 
        rlen=0
        do 1400 j=1,20
c 
        ii=ii+1
c 
        call chnkpnt(i,ts(j),xy,dxydt)
c 
        xsides(ii)=xy(1)
        ysides(ii)=xy(2)
c 
c 
        d=sqrt(dxydt(1)**2+dxydt(2)**2)
c 
        rlen=rlen+d*whts(j)
        weights(ii)=whts(j)
c 
        xsnorms(ii)=-dxydt(2)/d
        ysnorms(ii)=dxydt(1)/d
 1400 continue
c 
        do 1500 j=1,20
c 
        weights(ii-j+1)=weights(ii-j+1)*rlen/2
 1500 continue
c 
 1600 continue
c 
        npts=ii
c 
c        construct the matrix of interactions, minus the diagonal
c        blocks
c 
        done=1
        pi=atan(done)*4
  
        ir11=ir
  
        call prinf('in wholmatr before wholmat0, nsides=*',nsides,1)
  
        call wholmat0(ir11,nsides,chnkpnt,interact,wholemat)
  
  
        call prinf('in wholmatr after wholmat0, nsides=*',nsides,1)
  
c 
c       construct the diagonal blocks
c 
        ir22=ir
  
        do 2200 i=1,nsides
c 
        call o1rmtot(ir22,i,xs7,whts7,totmat,
     1      chnkpnt,interact,xys,whtsxys)
c 
        ir22=0
c 
        ii0=(i-1)*20
c 
        do 2150 ii=1,20
        do 2100 jj=1,20
c 
        wholemat(ii0+ii,ii0+jj)=totmat(ii,jj)
 2100 continue
 2150 continue
 2200 continue
  
        call prinf('exiting wholmatr, nsides=*',nsides,1)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine wholmat0(ir,nchunks,chnkpnt,interact,wholemat)
        implicit real *8 (a-h,o-z)
        save
        dimension wholemat(nchunks*20,nchunks*20),amatr(20,20)
c 
        nc20=nchunks*20
        do 1400 i=1,nc20
        do 1200 j=1,nc20
        wholemat(j,i)=0
 1200 continue
 1400 continue
c 
c       construct the matrix discretizing the Dirichlet problem,
c       one 20 * 20 submatrix at a time
c 
        irr=ir
c 
        do 3000 i=1,nchunks
        do 2800 j=1,nchunks
c 
        if(i .eq. j) goto 2800
c 
        call o2rsubm(irr,i,j,chnkpnt,interact,
     1      amatr)
c 
        irr=0
c 
        ii0=(i-1)*20
        jj0=(j-1)*20
        do 2600 ii=1,20
        do 2400 jj=1,20
c 
        wholemat(ii0+ii,jj0+jj)=amatr(ii,jj)
 2400 continue
 2600 continue
c 
 2800 continue
 3000 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o2rsubm(ir,ichunk,jchunk,chnkpnt,interact,
     1      amatr)
        implicit real *8 (a-h,o-z)
        save
        dimension amatr(20,24),quadr(24),targ(2),targnorm(2),
     1      dtardt(2),xs7(21),whts7(21),xys7(50),whtsxys7(100)
c 
c 
c        This subroutine constructs the matrix of interactions
c        of a single 20-point chunk with another 20 \times 20 chunk.
c 
c                      Input parameters:
c 
c  ir - the FORTRAN unit number from which the data are to be read;
c        ir .leq. 0 will cause the subroutine to skip the reading. In
c        this case, the parameters xs, whts are input parameters.
c  ichunk - the sequence number of the target shunk, to be used by
c       the subroutine chnkpnt
c  jchunk - the sequence number of the source shunk, to be used by
c       the subroutine chnkpnt
c  chnkpnt - the subroutine providing the parametrization of the
c        chunk on which the submatrix is being created. Conceptually,
c        chnkpnt implements a mapping of the interval [-1,1] to R^2;
c        the output parameters xy, dxydt are the image of the point
c        t in R^2 and the derivative of the image with respect to
c        t. In other words, dxydt is a tangent vector to the curve
c        at the point xy, but the length of dxydt is not equal to 1
c        (generally speaking)
c 
c        The calling sequence of chnkpnt is
c 
c         chnkpnt(ichunk,t,xy,dxydt)
c 
c          Input parameters of chnkpnt:
c 
c  ichunk - the sequence number of the chunk being processed
c  t - the point on the interval [-1,1] parametrizing the chunk
c        number ichunk
c 
c          Output parameters of chnkpnt:
c 
c  xy - the image in R^2 of the point t on the interval [-1,1]
c  dxydt - the derivative of xy with respect to t
c 
cccc-------------------------------------------------------------------
c 
cccc-------------------------------------------------------------------
c 
c  interact - the subroutine evaluating the interaction between two
c        points on the interval [-1,1]. The calling sequence of
c        interact is
c 
c         interact(source,target,fout)
c 
cccc-------------------------------------------------------------------
c 
c                   Output parameters:
c 
c  amatr - the 20 \times 20 matrix of interactions between the chinks
c       with numbers ichunk, jchunk
c 
c 
c       . . . construct the matrix discretizing the Dirichlet
c       problem, one 20 * 20 submatrix at a time
c 
  
        if(ir .ne. 0) call o1rnodes(xs7,whts7,n20)
  
        irr=ir
        do 3000 i=1,20
c 
        call chnkpnt(ichunk,xs7(i),targ,dtardt)
c 
        call o3rvecne(ier,irr,jchunk,targ,targnorm,
     1      xs7,whts7,quadr,chnkpnt,interact,xys7,whtsxys7)
c 
        irr=0
c 
        do 1200 j=1,20
c 
         amatr(i,j)=quadr(j)
 1200 continue
 3000 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o3rvecne(ier,ir,ichunk,targ,targnorm,
     1      xs,whts,vectquad,chnkpnt,interact,xys,whtsxys)
        implicit real *8 (a-h,o-z)
        save
        real *8 xs(1),whts(1),vectquad(1)
        dimension xys(2,1),whtsxys(1)
c 
        dimension xy(2),targ(2),targnorm(2),dxydt(2),
     1      thresh(2,26),dxydts(2,100),xynorm(2)
c 
        dimension numspts(20),
     1      points1(24,20),weights1(24,20),coefslrg1(20,24,20),
     2      points2(27,21),weights2(27,21),coefslrg2(20,27,21),
     3      points3(31,21),weights3(31,21),coefslrg3(20,31,21),
     4      points4(38,21),weights4(38,21),coefslrg4(20,38,21),
     5      points5(47,42),weights5(47,42),coefslrg5(20,47,42)
c 
        dimension thresh2(2,50)
c 
        external interact,chnkpnt
c 
         data ifcalled/0/
c 
c        This subroutine returns to the user the coefficients vectquad
c        of the linear form used for the evaluation of the potential
c        created at the user-supplied point targ by a "charge"
c        (dipole, etc.) distribution on the user-specified chunk in
c        R^2. The chunk is not assumed to be a line (i.e. it can be
c        curved), and the function is assumed to be "appropriately-
c        behaved". More specifically, the subroutine returns to the
c        user the coefficients vectquad connecting the values of the
c        "charge" (dipole, etc.) distribution at the support nodes
c        on the user-supplied chunk with the value of the potential
c        (normal derivative of the potential, etc) at the user-specified
c        point targ. The point is assumed to be arbitrary, except that
c        it can not be "too close" to the chunk on which the charge is
c        distributed. If the chunk is straight, the quadratures are valid
c        outside the intersection of two half-cones with vertex angles
c        of 10^o (5^o on each side of the chunk), with vertices at
c        the ends of the chunk. The quadratures are also invalid within
c        1.0E-6 *len from each of the vertices, with len the length of
c        the chunk.
c 
c                        Input parameters:
c 
c  ir - the FORTRAN unit number from which the data are to be read;
c        ir .leq. 0 will cause the subroutine to skip the reading. In
c        this case, the parameters xs, whts are input parameters.
c  ichunk - the sequence number of the chunk being processed
c  targ - the locatrion in R^2 of the point where the potential is
c        to be evaluated
c  targnorm - the normal to the curve (if it is needed) at the point
c        targ
c 
c 
cccc---------------------------------------------------------------------
c 
c  chnkpnt - the subroutine providing the parametrization of the
c        chunk on which the submatrix is being created. Conceptually,
c        chnkpnt implements a mapping of the interval [-1,1] to R^2;
c        the output parameters xy, dxydt are the image of the point
c        t in R^2 and the derivative of the image with respect to
c        t. In other words, dxydt is a tangent vector to the curve
c        at the point xy, but the length of dxydt is not equal to 1
c        (generally speaking)
c 
c        The calling sequence of chnkpnt is
c 
c         chnkpnt(ichunk,t,xy,dxydt)
c 
c          Input parameters of chnkpnt:
c 
c  ichunk - the sequence number of the chunk being processed
c  t - the point on the interval [-1,1] parametrizing the chunk
c        number ichunk
c 
c          Output parameters of chnkpnt:
c 
c  xs - the image in R^2 of the point t on the interval [-1,1]
c  dxydt - the derivative of xy with respect to t
c 
cccc-------------------------------------------------------------------
c 
cccc-------------------------------------------------------------------
c 
c  interact - the subroutine evaluating the interaction between two
c        points on the interval [-1,1]. The calling sequence of
c        interact is
c 
c         interact(xy,xynorm,targ,targnorm,fout)
c 
c                 Input parameters of interact:
c 
c  xy - the source point in R^2
c  xynorm - the (unit) normal at the source point
c  targ - the target point in R^2
c  targnorm - the (unit) notmal at the target
c 
c                 Output parameters of interact:
c 
c  fout - the potential (normal derivative, whatever) created at
c        the target by the source
c 
cccc-------------------------------------------------------------------
c 
c 
c                 Output parameters:
c 
c  xs - the 20 support nodes on the interval [-1,1]; PLEASE NOTE THAT
C         THIS IS AN OUTPUT PARAMETER ONLY IF IR > 0; OTHERWISE, IT
C         IS AN INPUT PARAMETER!
c  whts - quadrature weights (pretty universal) corresponding to the
c        nodes xs;  PLEASE NOTE THAT THIS IS AN OUTPUT PARAMETER ONLY
c        IF IR > 0; OTHERWISE, IT IS AN INPUT PARAMETER!
c  vectquad - the coefficients of the linear form converting
c        the values of the charge" (dipole, etc.) distribution at the
c        support nodes on the user-supplied chunk into the value of
c        the potential (normal derivative of the potential, etc.) at
c        the user-specified point targ.
c  xys - the images of the 20 support nodes in R^2
c  whtsxys - the quadrature weights corresponding to the images xys of
c        the 20 support nodes in R^2
  
c 
c        . . . if this is the first call to this subroutine
c              - construct the array of thresholds
c 
        if(ifcalled .eq. 1) goto 1300
        ifcalled=1
c 
        thresh(1,1)=1.0d-6
        thresh(2,1)=thresh(1,1)*2
c 
        do 1200 i=1,20
        thresh(1,i+1)=thresh(2,i)
        thresh(2,i+1)=thresh(1,i+1)*2
 1200 continue
c 
        d=sqrt(2.0d0)
c 
        thresh2(1,1)=1.0d-6
        thresh2(2,1)=thresh2(1,1)*d
c 
        do 1250 i=1,41
        thresh2(1,i+1)=thresh2(2,i)
        thresh2(2,i+1)=thresh2(1,i+1)*d
 1250 continue
c 
        call prin2('thresh as created*',thresh,42)
        call prin2('thresh2 as created*',thresh2,84)
c 
 1300 continue
c 
c       . . . if necessary, retrieve from disk all of the data
c 
        ier=0
c 
        if(ir .gt. 0)
     a    call o3rrdall(ir,xs,whts,
     1      points1,weights1,numspts,coefslrg1,
     2      points2,weights2,coefslrg2,npts2,nregions2,
     3      points3,weights3,coefslrg3,npts3,nregions3,
     4      points4,weights4,coefslrg4,npts4,nregions4,
     5      points5,weights5,coefslrg5,npts5,nregions5)
c 
c       construct and return to the user the locations in R^2 of
c       the support nodes, and the corresponding weights; determine
c       the length of this chunk
c 
        rlen=0
        do 1400 i=1,20
c 
        call chnkpnt(ichunk,xs(i),xys(1,i),dxydt)
c 
        dxydts(1,i)=dxydt(1)
        dxydts(2,i)=dxydt(2)
c 
        whtsxys(i)=whts(i)*sqrt(dxydt(1)**2+dxydt(2)**2)
c 
        rlen=rlen+whtsxys(i)
 1400 continue
c 
c       determine which azimuthal set of data should be used
c 
        done=1
        call chnkpnt(ichunk,done,xy,dxydt)
c 
        xp=xy(1)
        yp=xy(2)
c 
        d2=(targ(1)-xy(1))**2+(targ(2)-xy(2))**2
c 
        call chnkpnt(ichunk,-done,xy,dxydt)
  
        xm=xy(1)
        ym=xy(2)
        d1=(targ(1)-xy(1))**2+(targ(2)-xy(2))**2
c 
c        if the target is sufficiently far from the chunk
c        - use the quadrature based on the support nodes
c 
         d3=d1
         if(d2 .lt. d1) d3=d2
c 
         rr=sqrt(d3)
         if(rr .lt. rlen/2) goto 1550
c 
        do 1500 i=1,20
c 
        xy(1)=xys(1,i)
        xy(2)=xys(2,i)
        dxydt(1)=dxydts(1,i)
        dxydt(2)=dxydts(2,i)
c 
        d=sqrt(dxydt(1)**2+dxydt(2)**2)
        xynorm(1)=-dxydt(2)/d
        xynorm(2)=dxydt(1)/d
c 
        call interact(xy,xynorm,targ,targnorm,fout)
c 
        vectquad(i)=fout*d*whts(i)
 1500 continue
c 
        return
 1550 continue
c 
        done=1
        pi=atan(done)*4
        theta=pi/6
c 
        izone=1
        call ifinside(theta,xm,ym,xp,yp,targ(1),targ(2),ifin)
c 
        if(ifin .eq. 0) goto 2400
c 
        theta=pi/12
c 
        izone=2
        call ifinside(theta,xm,ym,xp,yp,targ(1),targ(2),ifin)
c 
        if(ifin .eq. 0) goto 2400
        theta=pi/24
c 
        izone=3
        call ifinside(theta,xm,ym,xp,yp,targ(1),targ(2),ifin)
  
        if(ifin .eq. 0) goto 2400
c 
        theta=pi/36
        izone=4
        call ifinside(theta,xm,ym,xp,yp,targ(1),targ(2),ifin)
c 
        if(ifin .eq. 0) goto 2400
c 
        ier=4
 2400 continue
c 
c       find the end of the interval nearest the target
c 
        ifright=1
        rr=d2
        if(d2 .gt. d1) then
            ifright=0
            rr=d1
        endif
c 
        rr=sqrt(rr)
        rr=rr/rlen*2
c 
c        find the region in which the target is located,
c        given that the target is in the 7.5 degree zone or farther
c 
         if(izone .eq. 4) goto 3000
c 
        do 2600 i=1,21
c 
        ii=-7
        ii=i
        if( (rr .ge. thresh(1,i)) .and.
     1      (rr .le. thresh(2,i) ) ) goto 2800
 2600 continue
c 
        ier=4
        return
 2800 continue
  
        goto 4200
c 
 3000 continue
c 
c        find the region in which the target is located,
c        given that the target is in the 5 degree zone or farther
c 
         if(izone .ne. 4) goto 4000
c 
        do 3200 i=1,42
c 
        ii=-7
        ii=i
        if( (rr .ge. thresh2(1,i)) .and.
     1      (rr .le. thresh2(2,i) ) ) goto 3400
 3200 continue
c 
        ier=4
        return
 3400 continue
c 
 4000 continue
c 
 4200 continue
c 
        if(izone .eq. 1)
     1      call onevecqu(points2,weights2,coefslrg2,
     2          npts2,nregions2,ii,ifright,chnkpnt,interact,
     3          targ,targnorm,vectquad,ichunk)
c 
        if(izone .eq. 2)
     1      call onevecqu(points3,weights3,coefslrg3,
     2          npts3,nregions3,ii,ifright,chnkpnt,interact,
     3          targ,targnorm,vectquad,ichunk)
c 
        if(izone .eq. 3)
     1      call onevecqu(points4,weights4,coefslrg4,
     2          npts4,nregions4,ii,ifright,chnkpnt,interact,
     3          targ,targnorm,vectquad,ichunk)
c 
        if(izone .eq. 4)
     1      call onevecqu(points5,weights5,coefslrg5,
     2          npts5,nregions5,ii,ifright,chnkpnt,interact,
     3          targ,targnorm,vectquad,ichunk)
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine onevecqu(points,weights,coefslrg,
     1      npts,nregions,ii,ifright,chnkpnt,interact,
     2      targ,targnorm,vectquad,ichunk)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension xy(2),xynorm(2),targ(2),targnorm(2),dxydt(2),
     1      formint(100),vectquad(1)
c 
        dimension points(npts,nregions),weights(npts,nregions),
     1      coefslrg(20,npts,nregions)
c 
        external interact,chnkpnt
c 
c 
c        This subroutine returns to the user the coefficients vectquad
c        of the linear form used for the evaluation of the potential
c        created at the user-supplied point targ by a "charge"
c        (dipole, etc.) distribution on the user-specified chunk in
c        R^2. The chunk is not assumed to be a line (i.e. it can be
c        curved), and the function is assumed to be "appropriately-
c        behaved". More specifically, the subroutine returns to the
c        user the coefficients vectquad of the linear form converting
c        the values of the charge" (dipole, etc.) distribution at the
c        support nodes on the user-supplied chunk into the value of
c        the potential (normal derivative of the potential, etc) at
c        the user-specified point targ. Please note that this subroutine
c        is the computational module used by the subroutine o3rvecne
c        (see). It depends on o3rvecne for all of the bookkeeping,
c        logic, storage-retrieval, etc; this subroutine has no uses as
c        a stand-alone device.
c 
c                        Input parameters:
c 
c  points - the nodes of all quadrature formulae for all nregions
c        regions, corresponding to the zone the target point targ is
c  weights - the quadrature weights corresponding to the nodes points
c  coefslrg - the coefficients of interpolation formulae converting the
c        values of "appropriately behaved" functions at the support
c        nodes into their values at thg nodes points of auxiliary
c        quadrature formulae.
c  npts - the number of nodes each of the auxiliary quadrature formulae
c        has.
c  nregions - the number of regions corresponding to the zone in which
c        the target point targ is located. As a matter of fact, this
c        parameter is not used by this subroutine, and is provided
c        solely to amuze the user
c  ii - the region in which the target point targ is located. Should
c        be between 1 and nregions
c  ifright - the parameter telling the subroutine which of the ends
c        of the interval (right or left) is nearest to the point targ:
c        ifright=1 means that the right end is nearer
c        ifright=0 means that the left end is nearer
c 
cccc---------------------------------------------------------------------
c 
c  chnkpnt - the subroutine providing the parametrization of the
c        chunk on which the submatrix is being created. Conceptually,
c        chnkpnt implements a mapping of the interval [-1,1] to R^2;
c        the output parameters xy, dxydt are the image of the point
c        t in R^2 and the derivative of the image with respect to
c        t. In other words, dxydt is a tangent vector to the curve
c        at the point xy, but the length of dxydt is not equal to 1
c        (generally speaking)
c 
c        The calling sequence of chnkpnt is
c 
c         chnkpnt(ichunk,t,xy,dxydt)
c 
c          Input parameters of chnkpnt:
c 
c  ichunk - the sequence number of the chunk being processed
c  t - the point on the interval [-1,1] parametrizing the chunk
c        number ichunk
c 
c          Output parameters of chnkpnt:
c 
c  xs - the image in R^2 of the point t on the interval [-1,1]
c  dxydt - the derivative of xy with respect to t
c 
cccc-------------------------------------------------------------------
c 
cccc-------------------------------------------------------------------
c 
c  interact - the subroutine evaluating the interaction between two
c        points on the interval [-1,1]. The calling sequence of
c        interact is
c 
c         interact(xy,xynorm,targ,targnorm,fout)
c 
c                 Input parameters of interact:
c 
c  xy - the source point in R^2
c  xynorm - the (unit) normal at the source point
c  targ - the target point in R^2
c  targnorm - the (unit) notmal at the target
c 
c                 Output parameters of interact:
c 
c  fout - the potential (normal derivative, whatever) created at
c        the target by the source
c 
cccc-------------------------------------------------------------------
c 
c  ichunk - the sequence number of the chunk being processed, to be
c      used by the subroutine chnkpnt (see above)
C  targ - the locatrion in R^2 of the point where the potential is
c        to be evaluated
c  targnorm - the normal to the curve (if it is needed) at the point
c        targ
c  ichunk - the sequence number of the chunk being processed, to be
c      used by the subroutine chnkpnt (see above)
c 
c                 Output parameters:
c 
c  vectquad - the coefficients of the linear form converting
c        the values of the charge" (dipole, etc.) distribution at the
c        support nodes on the user-supplied chunk into the value of
c        the potential (normal derivative of the potential, etc.) at
c        the user-specified point targ.
c 
c 
c       . . . construct the array of values of potentials created
c             at the user-supplied target nodes by charges at the
c             auxiliary nodes
c 
        naux=npts
        do 2200 j=1,naux
c 
        dd=points(j,ii)
        if(ifright .eq. 0) dd=-dd
c 
        call chnkpnt(ichunk,dd,xy,dxydt)
c 
        d=sqrt(dxydt(1)**2+dxydt(2)**2)
        xynorm(1)=-dxydt(2)/d
        xynorm(2)=dxydt(1)/d
        call interact(xy,xynorm,targ,targnorm,fout)
c 
        formint(j)=fout *d
 2200 continue
c 
c       . . . constructing the column corresponding to the support
c             node number k
c 
        do 3400 i=1,20
c 
        cd=0
        do 3200 j=1,naux
c 
        cd=cd+weights(j,ii)*formint(j)*coefslrg(i,j,ii)
 3200 continue
c 
        vectquad(i)=cd
 3400 continue
c 
        if(ifright .eq. 1) return
        do 3600 i=1,10
c 
        d=vectquad(i)
        vectquad(i)=vectquad(20-i+1)
        vectquad(20-i+1)=d
 3600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o1rmtot(ir,ichunk,xs,whts,totmat,
     1      chnkpnt,interact,xys,whtsxys)
        implicit real *8 (a-h,o-z)
        save
        real *8 rintmat(24,20),totmat(20,20),xs(1),whts(1)
        dimension numspts(21),weights(24,20),points(24,20),
     1      coefslrg(20,24,20),xys(2,1),whtsxys(1)
c 
        dimension xy(2),xynorm(2),targ(2),targnorm(2),dxydt(2),
     1      dtardt(2)
c 
        external interact,chnkpnt
c 
c        This subroutine constructs the matrix of interactions
c        on a single chunk discretized by 20 points. The chunk
c        must be parametrized by the interval [-1,1], the
c        parametrization being provided by the user-supplied
c        subroutine chnkpnt (see description below), and is
c        discretized by 20 "support nodes", returned by this
c        subroutine.
c 
c                        Input parameters:
c 
c  ir - the FORTRAN unit number from which the data are to be read;
c        ir .leq. 0 will cause the subroutine to skip the reading. In
c        this case, the parameters xs, whts are input parameters.
c 
cccc---------------------------------------------------------------------
c 
c  chnkpnt - the subroutine providing the parametrization of the
c        chunk on which the submatrix is being created. Conceptually,
c        chnkpnt implements a mapping of the interval [-1,1] to R^2;
c        the output parameters xy, dxydt are the image of the point
c        t in R^2 and the derivative of the image with respect to
c        t. In other words, dxydt is a tangent vector to the curve
c        at the point xy, but the length of dxydt is not equal to 1
c        (generally speaking)
c 
c        The calling sequence of chnkpnt is
c 
c         chnkpnt(ichunk,t,xy,dxydt)
c 
c          Input parameters of chnkpnt:
c 
c  ichunk - the sequence number of the chunk being processed
c  t - the point on the interval [-1,1] parametrizing the chunk
c        number ichunk
c 
c          Output parameters of chnkpnt:
c 
c  xs - the image in R^2 of the point t on the interval [-1,1]
c  dxydt - the derivative of xy with respect to t
c 
cccc-------------------------------------------------------------------
c 
cccc-------------------------------------------------------------------
c 
c  interact - the subroutine evaluating the interaction between two
c        points on the interval [-1,1]. The calling sequence of
c        interact is
c 
c         interact(source,target,fout)
c 
cccc-------------------------------------------------------------------
c 
c                 Output parameters:
c 
c  xs - the 20 support nodes on the interval [-1,1]
c  whts - quadrature weights (pretty universal) corresponding to the
c        nodes xs
c  totmat - the 20 \times 20 matrix of interactions between the
c       support nodes
c  xys - the images of the 20 support nodes in R^2
c  xys - the quadrature weights corresponding to the images xys of
c        the 20 support nodes in R^2
c 
c       . . . if necessary, retrieve from disk all of the data
c 
        if(ir .gt. 0) then
c 
            rewind(ir)
            call o1rread(ir,xs,whts,points,weights,numspts,coefslrg)
        endif
c 
c       construct and return to the user the locations in R^2 of
c       the support nodes, and the corresponding weights
c 
        do 1400 i=1,20
c 
        call chnkpnt(ichunk,xs(i),xys(1,i),xynorm)
c 
        whtsxys(i)=whts(i)*sqrt(xynorm(1)**2+xynorm(2)**2)
 1400 continue
c 
c       construct the matrix of values of potentials created
c       at the support nodes by charges at the auxiliary nodes
c 
        n=20
        do 2400 i=1,n
c 
        call chnkpnt(ichunk,xs(i),targ,dtardt)
c 
        d=sqrt(dtardt(1)**2+dtardt(2)**2)
        targnorm(1)=-dtardt(2)/d
        targnorm(2)=dtardt(1)/d
c 
        naux=numspts(i)
        do 2200 j=1,naux
c 
        call chnkpnt(ichunk,points(j,i),xy,dxydt)
c 
        d=sqrt(dxydt(1)**2+dxydt(2)**2)
        xynorm(1)=-dxydt(2)/d
        xynorm(2)=dxydt(1)/d
  
        call interact(xy,xynorm,targ,targnorm,fout)
c 
        rintmat(j,i)=fout *d
 2200 continue
 2400 continue
c 
c        construct the matrix connecting the values of the charge
c        at the support nodes with the values of the potential at the
c        same nodes
c 
        do 3800 k=1,n
c 
c       . . . constructing the column corresponding to the support
c             node number k
c 
        nk=numspts(k)
        do 3600 j=1,n
c 
        cd=0
        do 3400 i=1,nk
c 
        cd=cd+weights(i,k)*rintmat(i,k)*coefslrg(j,i,k)
 3400 continue
c 
        totmat(k,j)=cd
 3600 continue
 3800 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine soleval(xsides,ysides,xsnorms,ysnorms,
     1    whts,sol,n,interact,xreceiv,yreceiv,fout)
        implicit real *8 (a-h,o-z)
        save
        dimension xsides(1),ysides(1),xsnorms(1),ysnorms(1),
     1      whts(1),sol(1000),xy(2),xynorm(2),targ(2)
c 
c       This subroutine determines the field generated at a
c       point by a charge (dipole, etc) distriburion on a surface.
c 
c                  Input parameters:
c 
c  xsides, ysides - the coordinates of the diecretization of the
c       surface at which the charge distribution is tabulated
c  xnorms,ynorms - the normals at the nodes (xsides,ysides)
c  whts - the quadrature weights corresponding to the discretization
c       (xsides,ysides)
c  sol - the charge distribution
c  n - the number of nodes in the discretization of the boundary of
c       the scatterer
c 
c 
cccc-------------------------------------------------------------------
c 
c  interact - the subroutine evaluating the interaction between two
c        points on the interval [-1,1]. The calling sequence of
c        interact is
c 
c         interact(xy,xynorm,targ,targnorm,fout)
c 
c                 Input parameters of interact:
c 
c  xy - the source point in R^2
c  xynorm - the (unit) normal at the source point
c  targ - the target point in R^2
c  targnorm - the (unit) notmal at the target
c 
c                 Output parameters of interact:
c 
c  fout - the potential (normal derivative, whatever) created at
c        the target by the source
c 
cccc-------------------------------------------------------------------
c 
c  (xreceiv,yreceiv) - the coordinates of the point at which the
c       potential is to be evaluated
c 
c                   Output parameters:
c 
c  fout - the potential at the point (xreceiv,yreceiv)
c 
c 
c       . . . evaluate the potential at the receiver via the
c             solution of the boundary value problem
c 
        fout=0
        do 1800 i=1,n
c 
        xy(1)=xsides(i)
        xy(2)=ysides(i)
c 
        xynorm(1)=xsnorms(i)
        xynorm(2)=ysnorms(i)
c 
        targ(1)=xreceiv
        targ(2)=yreceiv
  
        call interact(xy,xynorm,targ,targnorm,ff)
c 
        fout=fout+sol(i)*ff*whts(i)
 1800 continue
  
        return
        end
  
c 
c 
c 
c 
c 
        subroutine resaeva(ir,chnkpnt,interact,sigma,npts,
     1     xreceiv,yreceiv,f1)
        implicit real *8 (a-h,o-z)
        save
        real *8 whts(10 000),xs(10 000),points(21),weights(21),
     1      amatr(100 000),sigma(1)
        dimension xy(2),xynorm(2),targ(2),dxydt(2)
        dimension xys(2,100 000),xynorms(2,100 000),
     1      whtsxys(100 000),eta(100 000)
c 
        external chnkpnt
c 
c       This subroutine resamples a user-specified function
c       (conceptually, a solution of a boundary integral equation)
c       to a "universal" oversamppled grid, and uses the resampled
c       version to evaluate the resulting field at a point. The user
c       can skip the resampling by changing the parameter itype (see
c       below). PLEASE NOTE THAT THIS IS A VERY SLOW CODE; IT IS
C       MEANT MORE OR LESS EXCLUSIVELY AS A DEBUGGING AND VERIFICATION
C       TOOL, TO THE POINT THAT EVEN MEMORY MANAGEMENT HAS BEEN SKIPPED.
C 
c                     Input parameters:
c 
c 
c  ir - the FORTRAN unit number from which the "universal" discretization
c       of the interval [-1,1] is to be read, together with the related
c       interpolation matrix; hopefully, these data have been stored on
c       disk via a preceding call to the subroutine o1rauwrt (see)
c 
cccc-------------------------------------------------------------------
c 
c  interact - the subroutine evaluating the interaction between two
c        points on the interval [-1,1]. The calling sequence of
c        interact is
c 
c         interact(xy,xynorm,targ,targnorm,fout)
c 
c                 Input parameters of interact:
c 
c  xy - the source point in R^2
c  xynorm - the (unit) normal at the source point
c  targ - the target point in R^2
c  targnorm - the (unit) notmal at the target
c 
c                 Output parameters of interact:
c 
c  fout - the potential (normal derivative, whatever) created at
c        the target by the source
c 
cccc-------------------------------------------------------------------
c 
c  sigma - the charge distribution tabulated at the original nodes
c       (as produced by the subroutine wholmatr)
c  npts - the number of nodes in the discretization of the boundary of
c       the scatterer
c  (xreceiv,yreceiv) - the coordinates of the point at which the
c       potential is to be evaluated
c 
c                   Output parameters:
c 
c  f1 - the potential at the point (xreceiv,yreceiv)
c 
c        . . . read from disk the nodes of the unoiversal discretization
c              of the interval [-1,1], and related data
c 
        call o1raurd(ir,n,xs,whts,amatr,points,weights)
c 
        call prinf('in resaeva, n=*',n,1)
c 
c        construct the nodes of the fine resampling of the
c        user-specified boundary
c 
         nsides=npts/20
         ii=0
         do 1600 i=1,nsides
c 
         do 1400 j=1,n
c 
         call chnkpnt(i,xs(j),xy,dxydt)
c 
        d=sqrt(dxydt(1)**2+dxydt(2)**2)
        xynorm(1)=-dxydt(2)/d
        xynorm(2)=dxydt(1)/d
c 
        ii=ii+1
c 
        xys(1,ii)=xy(1)
        xys(2,ii)=xy(2)
c 
        xynorms(1,ii)=xynorm(1)
        xynorms(2,ii)=xynorm(2)
c 
        whtsxys(ii)=whts(j)*sqrt(dxydt(1)**2+dxydt(2)**2)
 1400 continue
 1600 continue
c 
c       interpolate the user-supplied charge density to the
c       grossly resampled nodes
c 
        do 2600 i=1,nsides
c 
        ii=(i-1)*20+1
        jj=(i-1)*n+1
c 
        call matavec8(amatr,20,n,sigma(ii),eta(jj) )
 2600 continue
c 
c       evaluate the potential of the oversampled charge distribution
c       at the user-supplied point
c 
        f1=0
        do 1800 i=1,n*nsides
  
        xy(1)=xys(1,i)
        xy(2)=xys(2,i)
c 
        xynorm(1)=xynorms(1,i)
        xynorm(2)=xynorms(2,i)
c 
        targ(1)=xreceiv
        targ(2)=yreceiv
  
        call interact(xy,xynorm,targ,targnorm,fout)
c 
        f1=f1+eta(i)*fout*whtsxys(i)
 1800 continue
c 
        call prin2('in resaeva, f1=*',f1,1)
c 
        return
        end
