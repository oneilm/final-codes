        implicit real *8 (a-h,o-z)
        complex *16 rk,sol(10000),
     1      rhs3(10 000),ima,fields(10000),rhs(10000)
        dimension zs(2,1000),ts(10000),
     1      ts3(10 000),vects3(2,10 000),fsrea(10 000),
     2      fsima(10 000),farmag(10 000),xstest(10000),
     3      ystest(10000),xslege(10000),whts(10000),
     4      wslege(10000),w(1000 000)
c 
c 
        data ima/(0.0d0,1.0d0)/
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
c 
        PRINT *, 'ENTER npts '
        READ *,npts
        CALL PRINF('npts=*',npts,1 )
  
c 
        PRINT *, 'ENTER nn '
        READ *,nn
        CALL PRINF('nn=*',nn,1 )
  
c 
c       set up the test
c 
        eps=1.0d-10
c 
        rk=3.14 + ima/1.0d10
c 
        x0=0
        y0=0
        delta=0.1
c 
c 
c       construct the collection of charges
c 
        ncharges=12
  
        done=1
        pi=atan(done)*4
        h=2*pi/ncharges
c 
        a=done/5.3
        a=done/2.5
        a=done
        do 1600 i=1,ncharges
c 
        zs(1,i)=(i-1)*a
        zs(2,i)=0
cccc        zs(2,i)=(-1)**i * 0.1
c 
 1600 continue
c 
        call prin2('zs as constructed*',zs,2*ncharges)
c 
c       construct the points on the horizon where the far-field
c       values of the potential being designed are to be enforced
c 
c 
c        . . . construct the Legendre nodes, etc.
c 
        n=npts
c 
        itype=1
        call legeexps(itype,n/2,xslege,u,v,wslege)
  
c 
c       construct the images of the Legendre nodes on the interval [0,2*pi]
c 
        done=1
        pi=atan(done)*4
  
  
  
c 
c       construct the right-hand side and the weights
c 
        do 2200 i=1,npts+1
c 
        rhs(i)=0
        whts(i)=1
c 
 2200 continue
c 
        rhs(npts+1)=1
c 
        do 2400 i=1,n/2
c 
        ts(i)=delta*xslege(i)
        ts(n/2+i)=ts(i)+pi
 2400 continue
c 
        ts(n+1)=pi/2
  
        call prin2('ts as created*',ts,npts+1)
c 
c 
c        construct the normal vectors at all n+1 points
c 
        ts(n+1)=pi/2
c 
c       construct the matrices to be used in the design the
c       weird Green's function
c 
  
        call blind_bld(rk,x0,y0,zs,ncharges,
     1      ts,rhs,whts,npts+1,sol,eps,errmax,w)
  
  
        call prin2('and after blind_bld, sol=*',sol,ncharges*2)
        call prin2('and after blind_bld, errmax=*',errmax,1)
  
  
c 
c       evaluate the far field at a large collection of points
c       and plot it
c 
        ntest=100
        h=2*pi/ntest
        do 3200 i=1,ntest
c 
        ts3(i)=(i-1)*h
        vects3(1,i)=cos(ts3(i))
        vects3(2,i)=sin(ts3(i))
 3200 continue
c 
        call prin2('vects3=*',vects3,ntest*2)
c 
        call farfld_eval(rk,x0,y0,ntest,ncharges,zs,
     1      vects3,sol,rhs3)
  
        call prin2('rhs3 as calculated before plotting*',
     1      rhs3,ntest*2)
  
  
        do 3400 i=1,ntest
c 
        fsrea(i)=rhs3(i)
        fsima(i)=-ima*rhs3(i)
  
        farmag(i)=log10(abs(rhs3(i))+1.0d-15)
  
 3400 continue
c 
        iw=23
        call quagraph2(iw,ts3,fsrea,ntest,3,ts3,fsima,ntest,3,
     1      'both real and imaginary parts of far field*')
  
        call prin2('and again, zs=*',zs,ncharges*2)
        call prin2('and sol=*',sol,ncharges*2)
  
  
        iw=24
        call quagraph(iw,ts3,farmag,ntest,3,
     1      'magnitude of far field*')
  
c 
c        evaluate and plot a bunch of values of the obtained
c        filed at points in R^3
c 
  
        rr=50
        done=1
        pi=atan(done)*4
c 
        hh=2*pi/ntest
c 
        do 4200 i=1,ntest
c 
        ts(i)=(i-1)*h
        xstest(i)=rr*cos(ts(i))
        ystest(i)=rr*sin(ts(i))
c 
        call charges_field(rk,zs,sol,ncharges,
     1      xstest(i),ystest(i),fields(i))
c 
        fsrea(i)=fields(i)
        fsima(i)=-ima*fields(i)
        farmag(i)=log10(abs(fields(i))+1.0d-15)
c 
 4200 continue
  
        call prin2('fields as calculated*',fields,ntest*2)
  
        iw=25
        call quagraph2(iw,ts3,fsrea,ntest,3,ts,fsima,ntest,3,
     1      'field in physical space*')
  
  
  
        iw=26
        call quagraph(iw,ts,farmag,ntest,3,
     1      'magnitude of field in physical space*')
  
  
        m=ncharges
  
cccc        call testmatr(w,nn,sol,m)
        call testconvol(nn,sol,m)
  
  
        stop
        end
  
c 
c 
c 
c 
c 
        subroutine testconvol(n,sol,m)
        implicit real *8 (a-h,o-z)
        save
        complex *16 wsave(100000),fs(10000),ima,sol(1)
        dimension rnums(10 000),rmags(10 000)
c 
        data ima/(0.0d0,1.0d0)/
c 
  
c 
c        initialize the FFT
c 
  
        call DCFFTI(N,WSAVE)
c 
c       construct the FFT of my matrix
c 
        do 1400 i=1,n
c 
        fs(i)=0
 1400 continue
c 
        do 1600 i=1,m
c 
        fs(i)=sol(i)
 1600 continue
  
  
        call prin2('before fft, fs=*',fs,n*2)
  
  
  
        call DCFFTb(N,fs,WSAVE)
  
        call prin2('after fft, fs=*',fs,n*2)
  
  
c 
c        plot the magnitude of singular values
c 
        do 2200 i=1,n
c 
        rnums(i)=i
        rmags(i)=abs(fs(i))
 2200 continue
  
  
        iw=27
        call quagraph(iw,rnums,rmags,n,3,
     1      'singular values*')
  
  
        return
        end
c 
c 
c 
c 
c 
        subroutine testmatr(a,n,sol,m)
        implicit real *8 (a-h,o-z)
        save
        dimension uu(1000 000),vv(1000 000),
     1      w(2000 000),rsol(20),rsol2(20),rnums(10000),
     2      rmags(10 000)
  
        complex *16 a(n,n),sol(20),ima,ss(2000)
  
c 
        data ima/(0.0d0,1.0d0)/
c 
  
c 
        do 1400 i=1,n
        do 1200 j=1,n
c 
        a(j,i)=0
 1200 continue
c 
        a(i,i)=1
 1400 continue
c 
        d=1
        do 1600 i=1,n-12
c 
        do 1500 j=1,m
c 
        a(i,i+j-1)=sol(j)
 1500 continue
c 
 1600 continue
  
cccc        call prin2('a as created*',a,n*n)
c 
c        SVD the obtained matrix
c 
  
        eps=1.0d-12
        lw=2000 000
  
        call csvdpiv(ier,a,n,n,uu,vv,ss,ncols,eps,
     1      w,lw,ltot)
  
  
        call prin2('after svdpivot, ss=*',ss,ncols*2)
        call prinf('after svdpivot, ncols=*',ncols,1)
        call prinf('after svdpivot, ier=*',ier,1)
  
        call prin2('after svdpivot, cond=*',ss(1)/ss(ncols-5),2)
  
c 
c        plot the magnitude of singular values
c 
        do 2200 i=1,n
c 
        rnums(i)=i
        rmags(i)=abs(ss(i))
 2200 continue
  
  
        iw=27
        call quagraph(iw,rnums,rmags,n,3,
     1      'singular values*')
  
  
  
        return
        end
  
  
c 
c 
c 
c 
c 
  
        subroutine blind_bld(rk,x0,y0,zs,nzs,
     1      ts,rhs,whts,npts,sol,eps,errmax,w)
        implicit real *8 (a-h,o-z)
        save
        dimension zs(2,1),whts(1),ts(1),w(1)
c 
        complex *16 rhs(1),sol(1),rk
c 
c        This subroutine constructs a Helmholtz potential
c        in the plane having a user-prescribed far-field
c        signiture in certain (user-prescribed) directions.
c        The potential to be constructed is generated by
c        a collection of charges located at user-prescribed
c        positions in the plane, and the subroutine uses a
c        least-squares scheme to fiind the appropriate
c        (complex) intensities. PLEASE NOTE THAT THIS
c        SUBROUTINE IS NOT DESIGNED TO BE VERY FAST OR TO
c        USE MEMORY VERY EFFICIENTLY
c 
c                        Input parameters:
c 
c  rk - Helmholtz coefficient (complex)
c  (x0,y0) - the center of the far-field expansion
c  zs - the locations in R^2 of the charges in R^2
c  nzs - the number of charges in the array zs
c  ts - the collection of points on the interval [0,2*pi] (on
c        the horizon) at which the user has specified the values
c        of the far-field
c  rhs - the value of the far-field (complex) specified at the
c        points in array ts
c  whts - the weights (real) to be used in the least-squares
c        process. Please note that at this time, the only value
c        of the weights to have been tested is 1 (5.24.02)
c  npts - the number of points on the interval [0,2*pi] at which
c        the function is to be approximated (in other words, the
c        number of elements in each of the arrays zs, rhs, npts
c  eps - the accuracy to be used in the least-squares process
c        (the underlying code used by this subroutine is CLEASTSQ)
c 
c                        Output parameters:
c 
c  sol - the (complex) charges located at the points zs and generating
c        the required Helmholtz potential
c  errmax - the maximum deviation of the obtained far-field from the
c        user's specifications
c 
c                        Work arrays:
c 
c  w - must be at least 10*npts*nzs + 12*nzs+10*npts+1000
c        real *8 elements long.
c 
c 
c       . . . allocate memory for bld0
c 
        iamatr=1
        lamatr=(npts+1)*nzs*2 + 10
c 
        iuu=iamatr+lamatr
        luu=(npts+1)*nzs*2+10
c 
        ivv=iuu+luu
        lvv=(npts+1)*nzs*2+10
c 
        itt=ivv+lvv
        ltt=(npts+1)*nzs*2+10
c 
        iww=itt+ltt
        lww=(npts+1)*nzs*2+10
c 
        ivects=iww+lww
        lvects=npts*2+10
c 
        irnorms=ivects+lvects
        lrnorms=(npts+1+nzs)*2
c 
        irhs3=irnorms+lrnorms
        lrhs3=npts*2+10
c 
        iff=irhs3+lrhs3
        lff=npts*2+10
c 
        call blind_bld0(rk,npts-1,zs,nzs,w(iamatr),
     1      w(ivects),x0,y0,sol,errmax,eps,ts,rhs,whts,
     2      w(iuu),w(ivv),w(iww),w(itt), w(irnorms),w(irhs3),
     3      w(iff))
c 
        return
        end
c 
c 
c 
c 
c 
  
        subroutine blind_bld0(rk,npts,zs,ncharges,amatr,
     1      vects,x0,y0,sol,errmax,eps,ts,rhs,whts,
     2      uu,vv,ww,tt, rnorms,rhs3,f)
        implicit real *8 (a-h,o-z)
        save
        complex *16 rk,amatr(npts+1,1),rhs(1)
        dimension vects(2,1),zs(2,1),whts(1),ts(1)
c 
        complex *16 uu(1),ww(1),tt(1),
     1      rnorms(1000),vv(1),sol(1),
     2      rhs3(10 000),f(1000)
c 
c       construct the matrices to be used in the design of the
c       weird Green's function
c 
        call matrices_bld(rk,npts,zs,ncharges,amatr,vects,
     1      x0,y0,ts,f)
c 
c       use the least square procedure to construct the
c       beam-shaped green's function
c 
        call cleastsq(amatr,uu,ww,tt,npts+1,ncharges,ncols,
     1      rnorms,eps,vv)
  
        call prinf('in blindbld after cleastsq, ncols=*',ncols,1)
        call prin2('in blinbld fter cleastsq, rnorms=*',
     1      rnorms,ncols)
c 
c       solve the least squares problem
c 
        call cleasts2(uu,ww,tt,npts+1,ncharges,
     1      ncols,rhs,sol,vv,whts)
  
        call prin2('sol=*',sol,ncharges*2)
c 
c       test the result
c 
c 
        call farfld_eval(rk,x0,y0,npts,ncharges,zs,
     1      vects,sol,rhs3)
c 
        call prin2('rhs3=*',rhs3,npts*2+2)
  
        errmax=0
        do 2250 i=1,npts+1
c 
        d=abs(rhs3(i)-rhs(i))
        if(errmax .lt. d) errmax=d
 2250 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine matrices_bld(rk,npts,zs,ncharges,amatr,
     1      vects,x0,y0,ts,f)
        implicit real *8 (a-h,o-z)
        save
        complex *16 rk,f(1),amatr(npts+1,1)
        dimension ts(1),vects(2,1),zs(2,1)
c 
c 
        n=npts
c 
c        construct the normal vectors at all n+1 points
c 
        do 2600 i=1,n+1
c 
        vects(1,i)=cos(ts(i))
        vects(2,i)=sin(ts(i))
 2600 continue
c 
c       construct the matrix of the far-fields generated by
c       the user-supplied charges
c 
        do 3600 i=1,ncharges
c 
        call hexp2d(x0,y0,zs(1,i),zs(2,i),vects,
     1    rk,npts+1,f)
c 
        do 3200 j=1,npts+1
c 
        amatr(j,i)=f(j)
 3200 continue
c 
 3600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine hexp2d(x0,y0,x,y,vectors,
     1    rk,nphi,f)
        implicit real *8 (a-h,o-z)
        save
        complex *16 rk,f(1),ima
        dimension vectors(2,nphi)
         data ima/(0.0d0,1.0d0)/
c 
c       this subroutine constructs the far-field
c       representation of the potential created at a
c       unity charge located at the point (x,y,z).
c 
c               input parameters:
c 
c  x0,y0 - coordinates of the center of the expansion
c  x,y - coordinates of the charge whose potential
c           is being expanded
c  vectors - the normals to the unity circle
c  rk - the helmholtz coefficient
c  nphi - the number of elements in the discretization
c        of the circle
c 
c               output parameters:
c 
c  f - the far-field representation of the potential of
c       the charge being expanded
c 
        dx=x-x0
        dy=y-y0
c 
        do 1800 j=1,nphi
c 
c       calculate the projection of the location of the
c       charge on the direction in question
c 
        proj=vectors(1,j)*dx+vectors(2,j)*dy
c 
c       calculate the far-field expansion
c 
        f(j)=cdexp(-ima*rk*proj) /ima
 1800 continue
        return
        end
c 
c 
c 
c 
c 
  
        subroutine farfld_eval(rk,x0,y0,npts,ncharges,zs,
     1      vects,sol,rhs2)
        save
  
        implicit real *8 (a-h,o-z)
        complex *16 rk,sol(1),rhs2(1),f(10000)
        dimension zs(2,1),vects(2,1)
  
c 
        do 3200 i=1,npts+1
c 
        rhs2(i)=0
 3200 continue
c 
        do 3600 i=1,ncharges
c 
        call hexp2d(x0,y0,zs(1,i),zs(2,i),vects,
     1    rk,npts+1,f)
  
        do 3400 j=1,npts+1
c 
        rhs2(j)=rhs2(j)+f(j)*sol(i)
 3400 continue
 3600 continue
c 
        return
        end
c 
c 
c 
c 
c 
  
        subroutine charges_field(alpha,zs,coefs,nchgs,x,y,field)
        implicit real *8 (a-h,o-z)
        save
        complex *16 field,coefs(1),cd,h0,h1,alpha
        dimension zs(2,1)
c 
c       This subroutine evaluates the potential at the point x on
c       the real line of a collection of nchgs charges; the potential
c       is (obviously) given by the formula
c 
c        \sum_{j=1}^{nchgs} coefs(j) *
c          H_0( alpha*sqrt( (x-x_j)^2 + y_j^2) * exp(ima*beta*x)       (1)
c 
c              Input parameters:
c 
c  alpha - the Helmholtz coefficient (complex)
c  zs - the locations of the charges whose potential is to be evaluated
c  coefs - the (complex) intensitiesof the charges whose potential is
c       to be evaluated
c  nchgs - the number of said charges (the number of nodes in each of
c       the arrays zs, coefs)
c  (x,y) - the point on the real line where the potential is to be
c       evaluated
c 
c 
c       . . . construct the combined field of the nchgs charges
c 
c 
        field=0
        do 1200 i=1,nchgs
c 
        x0=zs(1,i)
        y0=zs(2,i)
c 
        d=sqrt((x-x0)**2+(y-y0)**2)
        cd=d*alpha
c 
        ifexpon=1
c 
        call hank103(cd,h0,h1,ifexpon)
c 
        field=field+h0*coefs(i)
 1200 continue
c 
        return
        end
