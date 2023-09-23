       implicit real *8 (a-h,o-z)
       dimension pin(2,3),w(1000 000)
c
       complex *16 cint
c
       external fun1
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
         PRINT *, 'ENTER n'
         READ *,n
         CALL PRINF('n=*',n,1 )
c
c        test the basic tensor product quadrature
c
         x1=1
         x2=2
c
         y1=3
         y2=5
c
         x1=0
         x2=1
c
         y1=0
         y2=1
c
        pin(1,1)=x1
        pin(2,1)=x2
c
        pin(1,2)=y1
        pin(2,2)=y2
c
c        integrate the function recursively
c
        eps=1.0d-10
        lenw=1000 000
        call adapgaus2d(ier,pin,n,fun1,par1,par2,eps,
     1      cint,numfunev,w,lenw,lused)


        call prin2('after adapgaus2d_rec, cint=*',cint,2)
        call prinf('after adapgaus2d_rec, numfunev=*',numfunev,1)



        stop
        end
c 
c 
c 
c 
c 
        subroutine fun1(x,y,par1,par2,f)
        implicit real *8 (a-h,o-z)
        complex *16 f,cd,ima
c
        save
        data ima/(0.0d0,1.0d0)/
c
        f=1
        cd=-1
        f=1-sqrt(cd+0.0000001)


        f=cos(x)+sin(y)+ima*(cos(2*x)+sin(3*y) )


cccc        f=cos(x)
cc        f=f*log(x**2+y**2) *sqrt(x**2+y**2)

        f=f/sqrt(x**2+y**2) *x**2

cccc        f=sin(x)+1

        
cc        call prin2('x =*',x,1)
cc        call prin2('y =*',y,1)
cc        call prin2('f as created*',f,2)

        f=log(x**2+y**2)


        return
        end
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c 
c        This is the end of debugging code, and the beginning of the
c        quadrature code proper
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c 
        subroutine adapgaus2d(ier,pin,n,fun,par1,par2,eps,
     1      cint,numfunev,w,lenw,lused)
        implicit real *8 (a-h,o-z)
        save
        dimension pin(2,3),par1(1),par2(1),w(1)
        complex *16 cint
c 
c        this subroutine uses the adaptive tensor product gaussian
c        quadrature to integrate a user-supplied function on a
c        rectangle in R^2 (also user-supplied). The edges of the 
c        said rectangle are parallel to the corresponding axes.        
c 
c                       input parameters:
c 
c  pin - array containing three pairs of numbers (x1,x2), (y1,y2),
c        (z1,z2); the domain of integration is bounded by the six
c        planes thus defined
c  m - the number of nodes in the quadrature to me used in each
c        direction on each cube. Recommended value: 16
c  fun - the user-supplied function to be integrated. the calling
c       sequence of fun must be
c 
c        call fun(x,y,par1,par2,f).                              (1)
c 
c        in (1), (x,y,z) is a point in R^3 where the function is 
c        to be evaluated, and par1, par2 are two parameters to be 
c        used by fun; they can be variables or arrays, real or 
c        integer, as desired. f is assumed to be complex *16
c  par1, par2 - parameters to be used by the user-supplied
c       subroutine fun (see above)
c  eps - the accuracy (absolute) to which the integral will be
c       evaluated
c  lenw - the length of the user-provided work array w (in real *8
c       locations
c 
c                       output parameters:
c 
c  ier - error return code.
c          ier=0 means normal conclusion
c          ier=2048 means that the user-supplied work array w is 
c                   too small
c          ier = 16 means that the maximum permitted recursion depth
c                   has been exceeded
c          ier = 32 means that the maximum permitted number of steps
c                   has been exceeded
c 
c  cint - the integral as evaluated
c  numfunev - the total number of function evaluations (calls to the
c         user-supplied subroutine fun) that has occured in the
c         calculation
c  lused - the actual number of elements in array w used by the 
c         subroutine (in real *8 locations)
c
c                         work arrays:
c 
c  w - must be sufficiently long
c 
c        . . . integrate the user-supplied function using the
c              adaptive gaussian quadratures
c
c        allocate memory
c
        ier=0
        istack=1
        lstack=2000
c
        irnodes=istack+lstack
        lrnodes=3*n**2+20 
c
        icoefs=irnodes+lrnodes
        lcoefs=n**2+20
c
        ivals=icoefs+lcoefs
        lvals=2000
c
        lused=ivals+lvals
        call prinf('and lused=*',lused,1)
c
        if(lused .gt. lenw) then
            ier=2048
            return
        endif
c
c        perform adaptive integration
c
        call adapgaus2d_rec(ier,pin,n,fun,par1,par2,eps,
     1      cint,numfunev,w(istack),w(irnodes),w(icoefs),w(ivals))
c
        return
        end
c 
c 
c 
c 
c 
        subroutine adapgaus2d_rec(ier,pin,n,fun,par1,par2,eps,
     1      cint,numfunev,stack,rnodes,coefs,vals)
        implicit real *8 (a-h,o-z)
        save
        dimension stack(2,2,1),pin(2,3),rnodes(1),
     1      coefs(1),zout(2,2,10),par1(1),par2(1)
        complex *16 vals(1),cint,cd,vvs(10)
c
c       . . . start the recursion
c 
        call adapgaus2d_arrm(pin,stack(1,1,1),4)
c 
        call adapgaus2d_onequadr(pin,n,fun,par1,par2,
     1      rnodes,coefs,vals(1))
c 
c       recursively integrate the thing
c 
        maxdepth=300
c        
        j=1
        numfunev=n**2
        nnmax=1 000 000
        cint=0
        ier=0
        maxrec=0
        do 3000 i=1,nnmax
        call prinf('i=*',i,1)
c
        if(j .gt. maxrec) maxrec=j
c 
c       subdivide the current cube
c 
        call adapgaus2d_squarediv(stack(1,1,j),zout)
c
c        calculate all four integrals
c
        cd=0
        do 1200 ll=1,4
c
        call adapgaus2d_onequadr(zout(1,1,ll),n,fun,par1,par2,
     1      rnodes,coefs,vvs(ll) )
c
        numfunev=numfunev+n**2
c
        cd=cd+vvs(ll)
 1200 continue
c
        dd=abs(cd-vals(j))
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
        cint=cint+cd
        j=j-1
c 
c        if the whole thing has been integrated - return
c 
        if(j .eq. 0) then
            call prin2('victory; integral=*',cint,2)
            call prinf('numfunev=*',numfunev,1)
            return
        endif
c
        goto 3000
 2000 continue
c 
c       if the function on this subinterval has not been
c       integrated with sufficient accuracy - move
c       down the stack
c 
        call adapgaus2d_arrm(zout,stack(1,1,j),48)
        call adapgaus2d_arrm(vvs,vals(j),8)
c 
        j=j+3
        if(j .gt. maxdepth) then 
            ier=32
            return
        endif
c
 3000 continue
        ier=64
        return
        end
c 
c 
c 
c 
c 
        subroutine adapgaus2d_arrm(x,y,n)
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
        subroutine adapgaus2d_squarediv(pin,pout)
        implicit real *8 (a-h,o-z)
        save
        dimension pin(2,3),pout(2,2,8)
c 
c        this subroutine subdivides the input cube
c        into 8 new ones
c 
        x1=pin(1,1)
        x3=pin(2,1)
        x2=(x1+x3)/2
c 
        y1=pin(1,2)
        y3=pin(2,2)
        y2=(y1+y3)/2
c 
c
        pout(1,1,1)=x1
        pout(2,1,1)=x2
c
        pout(1,2,1)=y1
        pout(2,2,1)=y2
c
c
        pout(1,1,2)=x2
        pout(2,1,2)=x3
c
        pout(1,2,2)=y1
        pout(2,2,2)=y2
c
c
        pout(1,1,3)=x1
        pout(2,1,3)=x2
c
        pout(1,2,3)=y2
        pout(2,2,3)=y3
c
        pout(1,1,4)=x2
        pout(2,1,4)=x3
c
        pout(1,2,4)=y2
        pout(2,2,4)=y3
c
        return
        end
c 
c 
c 
c 
c 
        subroutine adapgaus2d_onequadr(p,n,fun,par1,par2,
     1      rnodes,coefs,cint)
        implicit real *8 (a-h,o-z)
        save
        dimension ts(100),ws(100),rnodes(2,1),coefs(1),p(2,2),
     1      xs(100),ys(100),wsx(100),wsy(100)
c
        complex *16 cint,f,cd
c
        data nold/-1/
c
c       if the user so requested, initialize the quadratures
c
        if(n .ne. nold) then
            itype=1
            call legeexps(itype,n,ts,u,v,ws)
            nold=n
        endif 
c
c
        x1=p(1,1)
        x2=p(2,1)
c
        y1=p(1,2)
        y2=p(2,2)
c
        vx=(x2+x1)/2
        ux=(x2-x1)/2       
c
        vy=(y2+y1)/2
        uy=(y2-y1)/2       
c
        do 1100 i=1,n
c
        xs(i)=ux*ts(i)+vx
        wsx(i)=ws(i)*ux
c
        ys(i)=uy*ts(i)+vy
        wsy(i)=ws(i)*uy
c
 1100 continue
c
        ii=0
        do 1600 i=1,n
        do 1400 j=1,n
c
        ii=ii+1
        rnodes(1,ii)=xs(i)
        rnodes(2,ii)=ys(j)
        coefs(ii)=wsx(i)*wsy(j)
 1400 continue
 1600 continue
c
c       evaluate the integral via the tensor product quadrature
c
        cd=0
        do 2200 i=1,n**2
c
        call fun(rnodes(1,i),rnodes(2,i),par1,par2,f)
c
        cd=cd+coefs(i)*f
 2200 continue
c 
        cint=cd
        return
        end
