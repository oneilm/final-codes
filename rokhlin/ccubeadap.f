       implicit real *8 (a-h,o-z)
       dimension cent(3)
       dimension w(1000 000)
c 
       complex *16 cint
c 
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
c        integrate the test function adaptively
c 
        eps=1.0d-12
cccc        eps=1.0d-5
        a=2
        cent(1)=0
        cent(2)=0
        cent(3)=0
c 
        call ccubeadap(ier,cent,a,fun1,
     1      par1,par2,m,eps,
     2      cint,numint,numfunev,w)
  
  
  
        call prinf('after ccubeadap, ier=*',ier,1)
        call prinf('after ccubeadap, numint=*',numint,1)
        call prinf('after ccubeadap, numfunev=*',numfunev,1)
        call prin2('after ccubeadap, cint=*',cint,2)
  
  
c 
c 
        stop
        end
c 
c 
c 
c 
c 
        subroutine fun1(x,y,z,par1,par2,f)
        implicit real *8 (a-h,o-z)
c 
        save
        complex *16 f,ima
c 
        data ima/(0.0d0,1.0d0)/
  
c 
        f=x**3+y**4+z**5+cos(x)-sin(y)+1.7d0**z
  
  
cccc        f=f/sqrt(x**2+y**2+z**2)
  
        f=f+ima*f*(x**2+y**2+z**2+1)
  
cccc        f=f*sqrt(2.0d0)
  
  
cccc        call prin2('f=*',f,1)
  
cccc        f=1
c 
        return
        end
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code, and the beginning of the
c        quadrature code proper
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine ccubeadap(ier,cent,side,fun,par1,par2,m,eps,
     1      cint,numint,numfunev,w)
        implicit real *8 (a-h,o-z)
c 
        save
        external fun
        dimension par1(1),par2(1),cent(1),w(1)
        complex *16 cint
c 
c       this subroutine uses the adaptive tensor product gaussian
c       quadrature to integrate a user-supplied COMPLEX function on a
c       cube in R^3 (also user-supplied)
c 
c                       input parameters:
c 
c  cent - the center of the cube on which the function is to be
c       integrated
c  side - the length of the side of the said cube
c  fun - the user-supplied function to be integrated. The calling
c       sequence of fun must be
c 
c        call fun(x,y,z,par1,par2,f).                            (1)
c 
c        in (1), (x,y,z) is a point in the plane where the
c        function is to be evaluated, and par1, par2 are two
c        parameters to be used by fun; they can be variables
c        or arrays, real or integer, as desired.
c  par1, par2 - parameters to be used by the user-supplied
c       subroutine fun (see above)
c  m - the order of the quadrature to me used on each subcube
c  eps - the accuracy (absolute) to which the integral will be
c       evaluated
c 
c                       output parameters:
c 
c  ier - error return code.
c          ier=0 means normal conclusion
c          ier=8 means that at some point, the depth of recursion
c                reached 200. this is a fatal error.
c          ier=16 means that the total number of subcubes in the
c                adaptive subdivision of the user-suppplied cube
c                turned out to be greater than nnmax*8. This is a
c                fatal error.
c 
c  cint - the integral as evaluated
c  numint - the total number of subcubes examined
c  numfunev - the total number of function evaluations (calls to the
c         user-supplied subroutine fun) that has occured in the
c         calculation
c 
c                         work arrays:
c 
c  w - must be at least 4*m**3+5100 real *8 elements long
c 
c 
c        . . . allocate memory
c 
        izs=1
        lzs=m**3*3+10
c 
        iwhts=izs+lzs
        lwhts=m**3+10
c 
        istack=iwhts+lwhts
        lstack=4000
c 
        ivals=istack+lstack
        lvals=1000*2
c 
        nnmax=100000
        maxdepth=1000
c 
        call ccuberec(ier,cent,side,fun,
     1      par1,par2,m,nnmax,eps,
     2      cint,maxdepth,maxrec,numint,numfunev,w(izs),w(iwhts),
     3      w(istack),w(ivals) )
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine ccuberec(ier,cent,side,fun,
     1      par1,par2,m,nnmax,eps,
     2      cint,maxdepth,maxrec,numint,numfunev,zs,whts,
     3      stack,vals)
        implicit real *8 (a-h,o-z)
c 
        save
        external fun
c 
        dimension stack(4,1),par1(1),par2(1),cent(1),
     1      centsout(3,10),zs(1),whts(1)
c 
        complex *16 cint,vals(1),values(10),cddd
c 
c        initialize the single cube integrator
c 
        call cubegau0(m,zs,whts)
c 
c       . . . start the recursion
c 
        numfunev=0
c 
        stack(1,1)=cent(1)
        stack(2,1)=cent(1)
        stack(3,1)=cent(3)
        stack(4,1)=side
        a=side
c 
        numfunev=numfunev+m**3
c 
        call ccubegau1(m,zs,whts,cent,a,fun,par1,par2,vals(1) )
c 
c       recursively integrate the thing
c 
        j=1
        cint=0
        ier=0
        maxrec=0
        do 3000 i=1,nnmax
c 
        numint=i
        if(j .gt. maxrec) maxrec=j
c 
cccc        call prinf('j=*',j,1)
c 
c       subdivide the current cube
c 
        call cubediv(stack(1,j),stack(4,j),centsout,aout)
c 
        cddd=0
        do 1600 ll=1,8
c 
        call ccubegau1(m,zs,whts,centsout(1,ll),aout,
     1      fun,par1,par2,values(ll) )
        cddd=cddd+values(ll)
 1600 continue
c 
        numfunev=numfunev+8*m**3
c 
        dd=abs(cddd-vals(j))
c 
cccc         call prin2('in ccuberec, dd=*',dd,1)
c 
        ifdone=0
        if(dd .le. eps) ifdone=1
c 
c       if the function on this subinterval has been
c       integrated with sufficient accuracy - add the
c       value to that of the global integral and move up
c       in the stack
c 
        if(ifdone .eq. 0) goto 2000
c 
        cint=cint+cddd
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
        do 2100 ll=1,8
c 
        stack(1,j+ll-1)=centsout(1,ll)
        stack(2,j+ll-1)=centsout(2,ll)
        stack(3,j+ll-1)=centsout(3,ll)
c 
        stack(4,j+ll-1)=aout
c 
        vals(j+ll-1)=values(ll)
c 
 2100 continue
c 
        j=j+7
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
c 
        subroutine ccubegau1(n,zs,whts,cent,a,fun,par1,par2,cint)
        implicit real *8 (a-h,o-z)
        save
        dimension cent(3),zs(3,n,n,n),whts(n,n,n)
c 
        complex *16 f,cint
c 
c       convert the discretization of the standard cube
c       into the discretization of cube with side a and
c       center cent
c 
        cc=a/2
        ccc=cc**3
c 
        cint=0
        do 2000 i=1,n
        do 1800 j=1,n
        do 1600 k=1,n
c 
        x=zs(1,i,j,k)*cc+cent(1)
        y=zs(2,i,j,k)*cc+cent(2)
        z=zs(3,i,j,k)*cc+cent(3)
c 
        w=whts(i,j,k)*ccc
c 
        call fun(x,y,z,par1,par2,f)
c 
        cint=cint+w*f
c 
 1600 continue
 1800 continue
 2000 continue
c 
        return
        end
