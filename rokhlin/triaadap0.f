      implicit real *8 (a-h,o-z)
       dimension xy1(2),xy2(2),xy3(2),xs(10),ys(10),
     1     w(40 000)
       external fun1
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
C 
        PRINT *, 'ENTER m'
        READ *,m
        CALL PRINF('m=*',m,1 )
c 
c 
        xy1(1)=1
        xy1(2)=2
        xy2(1)=-3
        xy2(2)=4
        xy3(1)=5.333
        xy3(2)=11
c 
        call prin2('xy1 as created*',xy1,2)
        call prin2('xy2 as created*',xy2,2)
        call prin2('xy3 as created*',xy3,2)
c 
        xs(1)=xy1(1)
        xs(2)=xy2(1)
        xs(3)=xy3(1)
        xs(4)=xy1(1)
c 
        ys(1)=xy1(2)
        ys(2)=xy2(2)
        ys(3)=xy3(2)
        ys(4)=xy1(2)
c 
        iw=21
        ifjoin=1
        call RSPLOT(Xs,Ys,4,IW,'original triangle*',ifjoin)
c 
c       recursively evaluate the integral of the test function
c 
        nnmax=1000 000
        maxdepth=200
        eps=1.0d-8
        iw=21
c 
        call triaadap(ier,xy1,xy2,xy3,fun1,
     1      par1,par2,m,eps,rint,maxrec,numfunev,w)
c 
        call prin2('after triaadap, rint=*',rint,1)
        call prinf('after triaadap, numfunev=*',numfunev,1)
  
  
  
  
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine fun1(x,y,par1,par2,f)
        implicit real *8 (a-h,o-z)
c 
        save
        f=x**8+y**8
c 
        f=dsin(dsqrt(x**2+y**2))
  
        f=f*dlog( (x-1)**2+(y-2)**2)
cccc        f=dlog( (x-1)**2+(y-2)**2)
        return
        end
c 
c 
c 
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        this is the end of the debugging code and the
c        start of the actual quadrature routines.
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c 
c 
        subroutine triaadap(ier,vert1,vert2,vert3,fun,
     1      par1,par2,m,eps,rint,maxrec,numfunev,w)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),par1(1),par2(1),vert1(2),vert2(2),vert3(2)
c 
c       this subroutine uses the adaptive tensor product gaussian
c       quadrature to integrate a user-supplied function on a
c       triangle in the plane (also user-supplied)
c 
c                       input parameters:
c 
c  vert1,vert2,vert3 - the vertices in the plane of the triangle
c       over which the function is to be integrated
c  fun - the user-supplied function to be integrated. the calling
c       sequence of fun must be
c 
c        call fun(x,y,par1,par2,rint).                            (1)
c 
c        in (1), (x,y) is a point in the plane where
c        the function is to be evaluated, and par1, par2
c        are two parameters to be used by fun; they can be
c        variables or arrays, real or integer, as desired.
c        rint is assumed to be real *8.
c  par1, par2 - parameters to be used by the user-supplied
c       subroutine fun (see above)
c  m - the order of the quadrature to me used on each subinterval
c  eps - the accuracy (absolute) to which the integral will be
c       evaluated
c 
c                       output parameters:
c 
c  ier - error return code.
c          ier=0 means normal conclusion
c          ier=8 means that at some point, the depth of recursion
c                reached 200. this is a fatal error.
c          ier=16 means that the total number of subtriangles in the
c                adaptive subdivision of [a,b] turned out to be greater
c                than nnmax*4.  this is a fatal error.
c 
c  rint - the integral as evaluated
c  maxrec - the maximum depth to which the recursion went at its
c         deepest point. can not be greater than 200, since at that
c         point ier is set to 8 and the execution of the subroutine
c         terminated.
c  numfunev - the total number of function evaluations (calls to the
c         user-supplied subroutine fun) that has occured in the
c         calculation
c 
c                         work arrays:
c 
c  w - must be at least 3*m**2+4*m+1500 real *8 elements long
c 
c        . . . integrate the user-supplied function using the
c              adaptive gaussian quadratures
c 
        nnmax=100000
        maxdepth=200
        iiiw=21
c 
c        allocate memory for the subroutine trianrec
c 
        istack=1
        lstack=1207
c 
        iw=istack+lstack
        lw=3*m**2+50
c 
        ivals=iw+lw
        lvals=207
c 
        iww=ivals+lvals
        lww=4*m+10
c 
c       . . . integrate
c 
        call trianrec(ier,w(istack),vert1,vert2,vert3,fun,
     1      par1,par2,w(iw),m,w(ivals),nnmax,eps,
     2      rint,maxdepth,maxrec,numint,iiiw,numfunev,w(iww))
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine triaarrm(x,y,n)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1)
c 
        do 1200 i=1,n
        y(i)=x(i)
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine trianrec(ier,stack,z1,z2,z3,fun,
     1      par1,par2,w,m,vals,nnmax,eps,
     2      rint,maxdepth,maxrec,numint,iw,numfunev,ww)
        implicit real *8 (a-h,o-z)
        save
        dimension stack(2,3,1),w(1),vals(1),par1(1),par2(1)
        dimension z1(2),z2(2),z3(2),zout(2,3,10),ww(1)
c 
c       this subroutine uses the adaptive tensor product gaussian
c       quadrature to integrate a user-supplied function on a
c       triangle in the plane (also user-supplied)
c 
c                       input parameters:
c 
c  z1,z2,z3 - the vertices in the plane of the triangle
c       over which the function is to be integrated
c  fun - the user-supplied function to be integrated. the calling
c       sequence of fun must be
c 
c        call fun(x,y,par1,par2,rint).                            (1)
c 
c        in (1), (x,y) is a point in the plane where
c        the function is to be evaluated, and par1, par2
c        are two parameters to be used by fun; they can be
c        variables or arrays, real or integer, as desired.
c        rint is assumed to be real *8.
c  par1, par2 - parameters to be used by the user-supplied
c       subroutine fun (see above)
c  m - the order of the quadrature to be used on each subinterval
c  nnmax - maximum permitted number of subdivisions
c  eps - the accuracy (absolute) to which the integral will be
c       evaluated
c  maxrec - maximum permitted depth of recursion
c  iw - the FORTRAN unit number on which the subroutine will
c       plot the subdivision of the triangle used for the integration.
c       setting ls .leq. 0 suppresses the plotting
c 
c                       output parameters:
c 
c  ier - error return code.
c          ier=0 means normal conclusion
c          ier=8 means that at some point, the depth of recursion
c                reached 200. this is a fatal error.
c          ier=16 means that the total number of subtriangles in the
c                adaptive subdivision of the user-supplied one turned
c                out to be greater than nnmax*4.  this is a fatal error.
c  rint - the integral as evaluated
c  maxrec - the maximum depth to which the recursion went at its
c         deepest point. can not be greater than 200, since at that
c         point ier is set to 8 and the execution of the subroutine
c         terminated.
c  numint - the total number of triangles in the subdivision divided by
c         four. can not be greater than nnmax,  since at that
c         point ier is set to 16 and the execution of the subroutine
c         terminated.
c  numfunev - the total number of calls to the subroutine fun that
c         have been performed ("the number of function evaluations")
c 
c                         work arrays:
c 
c  stack - must be at least 1206 real *8 elements long
c  w - must be at least 3*m**2+50 real *8 elements long
c  vals - must be at least 200 real *8 elements long
c  ww - must be at leats m*4+10 real *8 elements long
c 
c       . . . start the recursion
c 
        numfunev=0
c 
        call triaarrm(z1,stack(1,1,1),2)
        call triaarrm(z2,stack(1,2,1),2)
        call triaarrm(z3,stack(1,3,1),2)
c 
        numfunev=numfunev+m**2
c 
        call triaint1(m,z1,z2,z3,fun,
     1      par1,par2,vals(1),w,ww)
c 
c       recursively integrate the thing
c 
        j=1
        rint=0
        ier=0
        maxrec=0
        do 3000 i=1,nnmax
ccc        call prinf('i=*',i,1)
        numint=i
        if(j .gt. maxrec) maxrec=j
cccc        call prinf('j=*',j,1)
c 
c       subdivide the current triangle
c 
        call triadiv(stack(1,1,j),stack(1,2,j),stack(1,3,j),zout)
c 
        numfunev=numfunev+m**2
        call triaint1(m,zout(1,1,1),zout(1,2,1),zout(1,3,1),fun,
     1      par1,par2,value1,w,ww)
c 
        numfunev=numfunev+m**2
        call triaint1(m,zout(1,1,2),zout(1,2,2),zout(1,3,2),fun,
     1      par1,par2,value2,w,ww)
c 
        numfunev=numfunev+m**2
        call triaint1(m,zout(1,1,3),zout(1,2,3),zout(1,3,3),fun,
     1      par1,par2,value3,w,ww)
c 
        numfunev=numfunev+m**2
        call triaint1(m,zout(1,1,4),zout(1,2,4),zout(1,3,4),fun,
     1      par1,par2,value4,w,ww)
c 
        dd=dabs(value1+value2+value3+value4-vals(j))
c 
cccc         call prin2('in trianrec, dd=*',dd,1)
c 
        ifdone=0
        if(dd .le. eps) ifdone=1
c 
c       if the function on this subinterval has been
c       integrated with sufficient accuracy - add the
c       value to that of the global integral and move up
c       in the stack
c 
        if(ifdone  .eq. 0) goto 2000
c 
        rint=rint+value1+value2+value3+value4
        j=j-1
c 
c        if the whole thing has been integrated - return
c 
        if(j .eq. 0) return
        goto 3000
 2000 continue
c 
c       if the function on this subinterval has not been
c       integrated with sufficient accuracy - move
c       down the stack
c 
        call triaarrm(zout(1,1,1),stack(1,1,j),6)
        call triaarrm(zout(1,1,2),stack(1,1,j+1),6)
        call triaarrm(zout(1,1,3),stack(1,1,j+2),6)
        call triaarrm(zout(1,1,4),stack(1,1,j+3),6)
c 
        vals(j)=value1
        vals(j+1)=value2
        vals(j+2)=value3
        vals(j+3)=value4
c 
        j=j+3
c 
c        plot the newly created triangles
c 
        if(iw .le. 0) goto 2200
  
        call triaplot(iw,zout(1,1,1) )
c 
        call triaplot(iw,zout(1,1,2) )
c 
        call triaplot(iw,zout(1,1,3) )
c 
        call triaplot(iw,zout(1,1,4) )
c 
 2200 continue
c 
c       if the depth of the recursion has become excessive - bomb
c 
        if(j .le. maxdepth) goto 3000
        ier=8
        return
 3000 continue
        ier=16
        return
        end
c 
c 
c 
c 
  
        subroutine triaplot(iw,z)
        implicit real *8 (a-h,o-z)
        save
        dimension z(2,3)
c 
 1200 format(2x,e11.5,2x,e11.5)
c 
        write(iw,1200) z(1,1),z(2,1)
c 
        write(iw,1200) z(1,2),z(2,2)
c 
        write(iw,1200) z(1,3),z(2,3)
c 
        write(iw,1200) z(1,1),z(2,1)
c 
 1400 format(2x, 'join')
        write(iw,1400)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine triadiv(z1,z2,z3,zout)
        implicit real *8 (a-h,o-z)
        save
        dimension z1(2),z2(2),z3(2),zout(2,3,4),w12(2),
     1      w13(2),w23(2)
c 
c        this subroutine subdivides the input triangle
c        with the nodes z1, z2, z3 into four triangles
c 
c        . . . construct the middles of the sides
c 
        w12(1)=(z1(1)+z2(1))/2
        w12(2)=(z1(2)+z2(2))/2
c 
        w13(1)=(z1(1)+z3(1))/2
        w13(2)=(z1(2)+z3(2))/2
c 
        w23(1)=(z2(1)+z3(1))/2
        w23(2)=(z2(2)+z3(2))/2
c 
c       . . . construct the new triangles
c 
        zout(1,1,1)=w12(1)
        zout(2,1,1)=w12(2)
c 
        zout(1,2,1)=z1(1)
        zout(2,2,1)=z1(2)
c 
        zout(1,3,1)=w13(1)
        zout(2,3,1)=w13(2)
c 
c 
        zout(1,1,2)=w13(1)
        zout(2,1,2)=w13(2)
c 
        zout(1,2,2)=z3(1)
        zout(2,2,2)=z3(2)
c 
        zout(1,3,2)=w23(1)
        zout(2,3,2)=w23(2)
c 
c 
        zout(1,1,3)=w23(1)
        zout(2,1,3)=w23(2)
c 
        zout(1,2,3)=z2(1)
        zout(2,2,3)=z2(2)
c 
        zout(1,3,3)=w12(1)
        zout(2,3,3)=w12(2)
c 
c 
        zout(1,1,4)=w23(1)
        zout(2,1,4)=w23(2)
c 
        zout(1,2,4)=w12(1)
        zout(2,2,4)=w12(2)
c 
        zout(1,3,4)=w13(1)
        zout(2,3,4)=w13(2)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine triaint1(n,vert1,vert2,vert3,fun,
     1      par1,par2,rint,w,work)
        implicit real *8 (a-h,o-z)
        save
        dimension vert1(1),vert2(1),vert3(1),w(1),work(450)
        data nold/-10/
c 
c        construct the quadrature formula on this triangle
c 
        ifinit=1
        if(n .eq. nold) ifinit=0
c 
        irnodes=1
        lrnodes=n**2 *2 +10
c 
        iweights=irnodes+lrnodes
        lweights=n**2+5
c 
        call triagauc(n,vert1,vert2,vert3,w(irnodes),
     1      w(iweights),ifinit,work)
c 
c       integrate the user-specified function
c 
        call triaint0(n,w(irnodes),w(iweights),fun,par1,par2,rint)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine triaint0(n,rnodes,weights,fun,par1,
     1    par2,rint)
        implicit real *8 (a-h,o-z)
        save
        dimension rnodes(2,n,n),weights(n,n),par1(1),par2(1)
c 
c       this subroutine integrates the user-specified function
c       fun over a triangle in the plane; it uses the quadrature
c       formula with n**2 nodes rnodes and weights weights;
c       presumably, these have been created by a prior call to
c       the subroutine triagauc (see)
c 
        rint=0
        do 1400 i=1,n
        do 1200 j=1,n
c 
        call fun(rnodes(1,j,i),rnodes(2,j,i),par1,par2,d)
  
cccc        call prin2('in triinte, d=*',d,1)
c 
        rint=rint+d*weights(j,i)
 1200 continue
 1400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine triagauc(n,vert1,vert2,vert3,rnodes,
     1      weights,ifinit,w)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),vert1(1),
     1      vert3(1),vert2(1),rnodes(1),weights(1)
c 
c       this subroutine constructs a tensor product Gaussian
c       quadrature formula on the triangle in the plane
c       specified by its vertices. Note that the one-dimensional
c       Gaussian quadratures are lined up parallel to the side of
c       the triangle opposite the vertex vert1; the total number
c       of nodes and weights to be created is n**2. Actually, this
c       is simply a memory manager for the subroutine triagau0,
c       which does all the work.
c 
c                              input parameters:
c 
c  n - the number of nodes in each direction. note that the
c       total number of nodes to be created is n**2, and the
c       algebraic order of the quadrature is (n-1)*2
c  vert1,vert2,vert3 - the vertices of the triangle on which the
c       quadrature rule is to be constructed
c 
c                              output parameters:
c 
c  rnodes - the array of n**2 nodes in the plane (all inside the
c       user-supplied traiangle)
c  weights - the quadrature weights corresponding to the nodes rnodes
c 
c                              work arrays:
c 
c  w - must be at least 4*n+5 real *8 locations
c 
c 
c        . . . allocate memory in the work array
c 
        it=1
        lt=n+1
c 
        iwhts=it+lt
        lwhts=n+2
c 
        itvert=iwhts+lwhts
        ltvert=n+2
c 
        iwhtsver=itvert+ltvert
        lwhtsver=n+1
c 
c       construct the quadrature formula on the triangle
c 
        call triagau0(n,vert1,vert2,vert3,rnodes,weights,
     1      ifinit,w(it),w(iwhts),w(itvert),w(iwhtsver) )
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine triagau0(n,vert1,vert2,vert3,rnodes,weights,
     1      ifinit,t,whts,tvert,whtsvert)
        implicit real *8 (a-h,o-z)
        save
        dimension w(12),z1(2),z2(2),z3(2),vert1(1),
     1      vert3(1),vert2(1),rnodes(2,n,n),weights(n,n)
c 
        dimension t(1),whts(1),tvert(1),whtsvert(1)
c 
c       this subroutine constructs a tensor product Gaussian
c       quadrature formula on the triangle in the plane
c       specified by its vertices. Note that the one-dimensional
c       Gaussian quadratures are lined up parallel to the side of
c       the triangle opposite the vertex vert1; the total number
c       of nodes and weights to be created is n**2.
c 
c                              input parameters:
c 
c  n - the number of nodes in each direction. note that the
c       total number of nodes to be created is n**2, and the
c       algebraic order of the quadrature is (n-1)*2
c  vert1,vert2,vert3 - the vertices of the triangle on which the
c       quadrature rule is to be constructed
c 
c                              output parameters:
c 
c  rnodes - the array of n**2 nodes in the plane (all inside the
c       user-supplied traiangle)
c  weights - the quadrature weights corresponding to the nodes rnodes
c 
c                              work arrays:
c 
c  t,whts,tvert,whtsvert - must be n+1 real *8 locations each
c 
c 
c       . . . if required, construct the Gaussian nodes and
c             weights on the interval [-1,1]
c 
        if(ifinit .eq. 1) call gauswhts(n,t,whts)
c 
c        construct affine transformation putting one side
c       of the user-specified triangle on the real axis. The side
c       to be put on the real axis is the one opposite the vertex vert1.
c 
        call trianini(vert1,vert2,vert3,w)
c 
c       construct the discretization of the version of the triangle
c       moved so that its side is on the real axis
c 
        call trianfor(w,vert1,z1)
        call trianfor(w,vert2,z2)
        call trianfor(w,vert3,z3)
c 
c       . . . construct the gaussian discretization of the
c             hight of the triangle
c 
        b=z1(2)
        u=b/2
        v=b/2
        d=dabs(b/2)
        do 1200 i=1,n
        tvert(i)=u*t(i)+v
        whtsvert(i)=whts(i)*d
 1200 continue
c 
c        find the equations of the two non-horizontal sides of
c        the shifted triangle
c 
        call trialine(z1(1),z2(1),z1(2),z2(2),a12,b12,c12)
c 
        call trialine(z1(1),z3(1),z1(2),z3(2),a13,b13,c13)
c 
c        one after another, construct the horizontal discretizations
c 
        do 2000 i=1,n
c 
c       find the current this horizontal section
c 
        y=tvert(i)
        x12=-(b12*y+c12)/a12
        x13=-(b13*y+c13)/a13
  
cccc        call prin2('x12=*',x12,1)
cccc        call prin2('x13=*',x13,1)
c 
        d=dabs(x12-x13)/2
c 
        u=(x13-x12)/2
        v=(x13+x12)/2
c 
cccc        call prin2('second u=*',u,1)
cccc        call prin2('second v=*',v,1)
c 
        do 1400 j=1,n
c 
        thor=u*t(j)+v
c 
        rnodes(1,j,i)=thor
        rnodes(2,j,i)=y
c 
        weights(j,i)=d*whts(j)*whtsvert(i)
c 
 1400 continue
 2000 continue
c 
c       now, move the nodes in the plane back to the original triangle
c 
        do 2400 i=1,n
        do 2200 j=1,n
c 
        call trianbak(w,rnodes(1,j,i),z3)
        rnodes(1,j,i)=z3(1)
        rnodes(2,j,i)=z3(2)
 2200 continue
 2400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine trialine(x1,x2,y1,y2,a,b,c)
        implicit real *8 (a-h,o-z)
c 
c        construct the equation of the line passing through the
c        two points (x1,y1), (x2,y2)
c 
        save
        dx=x2-x1
        dy=y2-y1
c 
        if(dabs(dx) .lt. dabs(dy)) goto 1200
c 
        b=1
        a=-dy/dx
        c=-(a*x1+b*y1)
c 
        goto 1400
c 
 1200 continue
c 
        a=1
        b=-dx/dy
        c=-(a*x1+b*y1)
c 
 1400 continue
c 
c       . . . normalize the coefficients
c 
        d=a**2+b**2
        d=dsqrt(d)
c 
        a=a/d
        b=b/d
        c=c/d
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine trianini(vert1,vert2,vert3,w)
        implicit real *8 (a-h,o-z)
        save
        dimension vert1(2),vert2(2),vert3(2),w(1),zin(1),zout(1)
        data ixshift/1/,iyshift/2/,iumat/3/,ivmat/7/
c 
c       this subroutine constructs an affine transformation putting one side
c       of the user-specified triangle on the real axis. The side
c       to be put on the real axis is the one opposite the vertex vert1.
c       this subroutine also constructs the transformation inverse to the
c       above one. The actual transformations are applied to particular
c       user-specified points by the entries trianfor and trianbak of this
c       subroutine
c 
c                              input parameters:
c 
c  vert1,vert2,vert3 - the vertices of the triangle one of whose
c       sides is to be put on the x axis
c 
c                              output parameters:
c 
c  w - the real *8 array of length 10 containing the transformations;
c       it is to be used by the entries trianfor, trianbak (see below)
c 
        call triambld(vert1(1),vert1(2),vert2(1),vert2(2),vert3(1),
     1      vert3(2),w(ixshift),w(iyshift),w(iumat),w(ivmat) )
c 
        return
c 
c 
c 
c 
        entry trianfor(w,zin,zout)
c 
c       this entry applies to the user-specified point zin \in R^2
c       the first of the transformations constructed by the entry
c       trianini (see above), obtaining the point zout \in R^2.
c 
c 
        call triarotf(zin(1),zin(2),w(ixshift),w(iyshift),
     1      w(iumat),zout(1),zout(2) )
c 
cccc        call prin2('in trianfor, zout=*',zout,2)
  
        return
c 
c 
c 
c 
        entry trianbak(w,zin,zout)
c 
c       this entry applies to the user-specified point zin \in R^2
c       the second of the transformations constructed by the entry
c       trianini (see above), obtaining the point zout \in R^2.
c 
        call triarotb(zin(1),zin(2),w(ixshift),w(iyshift),
     1      w(ivmat),zout(1),zout(2) )
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine triarotf(x,y,xshift,yshift,umat,xout,yout)
        implicit real *8 (a-h,o-z)
        save
        dimension umat(2,2)
c 
c        apply the shift
c 
        xx=x+xshift
        yy=y+yshift
  
cccc         call prin2('in triarotf, xx=*',xx,1)
cccc         call prin2('in triarotf, yy=*',yy,1)
  
c 
c        apply the rotation
c 
        xout=umat(1,1)*xx+umat(1,2)*yy
        yout=umat(2,1)*xx+umat(2,2)*yy
        return
        end
c 
c 
c 
c 
c 
        subroutine triarotb(x,y,xshift,yshift,vmat,xout,yout)
        implicit real *8 (a-h,o-z)
        save
        dimension vmat(2,2)
  
c 
c        apply the rotation
c 
        xx=vmat(1,1)*x+vmat(1,2)*y
        yy=vmat(2,1)*x+vmat(2,2)*y
c 
cccc         call prin2('in triarotb, xx=*',xx,1)
cccc         call prin2('in triarotb, yy=*',yy,1)
c 
c        apply the shift
c 
        xout=xx-xshift
        yout=yy-yshift
        return
        end
c 
c 
c 
c 
c 
        subroutine triambld(x1,y1,x2,y2,x3,y3,
     1      xshift,yshift,umat,vmat)
        implicit real *8 (a-h,o-z)
        save
        dimension umat(2,2),vmat(2,2)
c 
c        construct a motion in the plane that will put the side of
c        the triangle opposite the vertex (x1,y1) on th2 x axis
c 
c        . . . the translation
c 
        xshift=-x2
        yshift=-y2
c 
        xx1=x1+xshift
        xx2=x2+xshift
        xx3=x3+xshift
c 
        yy1=y1+yshift
        yy2=y2+yshift
        yy3=y3+yshift
c 
c        . . . rotation
c 
        u=(x3-x2)
        v=y3-y2
        d=dsqrt(u**2+v**2)
        u=u/d
        v=v/d
c 
        umat(1,1)=u
        umat(1,2)=v
        umat(2,1)=-v
        umat(2,2)=u
c 
        vmat(1,1)=umat(1,1)
        vmat(2,2)=umat(2,2)
        vmat(1,2)=umat(2,1)
        vmat(2,1)=umat(1,2)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine gauswhts(n,ts,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension ts(1),whts(1)
c 
c        this subroutine constructs the nodes and the
c        weights of the n-point gaussian quadrature on
c        the interval [-1,1]
c 
c                input parameters:
c 
c  n - the number of nodes in the quadrature
c 
c                output parameters:
c 
c  ts - the nodes of the n-point gaussian quadrature
c  w - the weights of the n-point gaussian quadrature
c 
c       . . . construct the array of initial approximations
c             to the roots of the n-th legendre polynomial
c 
        eps=1.0d-14
        ZERO=0
        DONE=1
        pi=datan(done)*4
        h=pi/(2*n)
        do 1200 i=1,n
        t=(2*i-1)*h
        ts(n-i+1)=dcos(t)
1200  CONTINUE
c 
c         use newton to find all roots of the legendre polynomial
c 
        ts(n/2+1)=0
        do 2000 i=1,n/2
c 
        xk=ts(i)
        deltold=1
        do 1400 k=1,10
        call legpol(xk,n,pol,der)
        delta=-pol/der
cccc         call prin2('delta=*',delta,1)
        xk=xk+delta
        if(deltold .le. eps) goto 1600
        deltold=dabs(delta)
 1400 continue
 1600 continue
        ts(i)=xk
        ts(n-i+1)=-xk
 2000 continue
c 
c       now, use the explicit integral formulae
c       to obtain the weights
c 
        a=-1
        b=1
        do 2200 i=1,n/2+1
        call prodend(a,ts,n,i,fm)
        call prodend(b,ts,n,i,fp)
        whts(i)=fp-fm
        whts(n-i+1)=whts(i)
 2200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine legpol(x,n,pol,der)
        implicit real *8 (a-h,o-z)
c 
        save
        pkm1=1
        pk=x
c 
        pk=1
        pkp1=x
c 
c        if n=0 or n=1 - exit
c 
        if(n .ge. 2) goto 1200
        pol=1
        der=0
        if(n .eq. 0) return
c 
        pol=x
        der=1
        return
 1200 continue
c 
c       n is greater than 1. conduct recursion
c 
        do 2000 k=1,n-1
        pkm1=pk
        pk=pkp1
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
 2000 continue
c 
c        calculate the derivative
c 
        pol=pkp1
        der=n*(x*pkp1-pk)/(x**2-1)
        return
        end
c 
c 
c 
c 
c 
        subroutine prodend(x,xs,n,i,f)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1)
c 
c      evaluate the product
c 
        f=1
        dlarge=1.0d20
        dsmall=f/dlarge
c 
        large=0
        do 2000 j=1,n
        dd=dabs(f)
        if( dd .gt. dsmall) goto 1200
         f=f*10000
         large=large-1
 1200 continue
c 
        if( dd .lt. dlarge) goto 1400
        f=f/10000
        large=large+1
 1400 continue
        if(j .eq. i) goto 2000
        f=f*(x-xs(j))/(xs(i)-xs(j))
 2000 continue
        d10000=10000
        f=f*d10000**large
        f=f**2*(x-xs(i))
        return
        end
  
  
  
