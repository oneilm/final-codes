       implicit real *8 (a-h,o-z)
       dimension cent(3)
       dimension w(1000 000)
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
        call cubeadap(ier,cent,a,fun1,
     1      par1,par2,m,eps,
     2      rint,numint,numfunev,w)
  
  
  
        call prinf('after cubeadap, ier=*',ier,1)
        call prinf('after cubeadap, numint=*',numint,1)
        call prinf('after cubeadap, numfunev=*',numfunev,1)
        call prin2('after cubeadap, rint=*',rint,1)
  
  
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
        f=x**3+y**4+z**5+cos(x)-sin(y)+1.7d0**z
  
  
cccc        f=f/sqrt(x**2+y**2+z**2)
  
cccc        f=f*(x**2+y**2+z**2+1)
  
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
        subroutine cubeadap(ier,cent,side,fun,par1,par2,m,eps,
     2      rint,numint,numfunev,w)
        implicit real *8 (a-h,o-z)
c 
        save
        external fun
        dimension par1(1),par2(1),cent(1),w(1)
c 
c       this subroutine uses the adaptive tensor product gaussian
c       quadrature to integrate a user-supplied function on a
c       cube in R^2 (also user-supplied)
c 
c                       input parameters:
c 
c  cent - the center of the cube on which the function is to be
c       integrated
c  side - the length of the side of the said cube
c  fun - the user-supplied function to be integrated. the calling
c       sequence of fun must be
c 
c        call fun(x,y,z,par1,par2,rint).                            (1)
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
c  rint - the integral as evaluated
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
        lvals=1000
c 
        nnmax=100000
        maxdepth=1000
c 
        call cuberec(ier,cent,side,fun,
     1      par1,par2,m,nnmax,eps,
     2      rint,maxdepth,maxrec,numint,numfunev,w(izs),w(iwhts),
     3      w(istack),w(ivals) )
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cuberec(ier,cent,side,fun,
     1      par1,par2,m,nnmax,eps,
     2      rint,maxdepth,maxrec,numint,numfunev,zs,whts,
     3      stack,vals)
        implicit real *8 (a-h,o-z)
c 
        save
        external fun
c 
        dimension stack(4,1),vals(1),par1(1),par2(1),cent(1),
     1      centsout(3,10),zs(1),whts(1),values(10)
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
        call cubegau1(m,zs,whts,cent,a,fun,par1,par2,vals(1) )
c 
c       recursively integrate the thing
c 
        j=1
        rint=0
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
        ddd=0
        do 1600 ll=1,8
c 
        call cubegau1(m,zs,whts,centsout(1,ll),aout,
     1      fun,par1,par2,values(ll) )
        ddd=ddd+values(ll)
 1600 continue
c 
        numfunev=numfunev+8*m**3
c 
        dd=dabs(ddd-vals(j))
c 
cccc         call prin2('in cuberec, dd=*',dd,1)
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
        rint=rint+ddd
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
        subroutine cubediv(cent,a,centout,aout)
        implicit real *8 (a-h,o-z)
        save
        dimension cent(3),centout(3,8)
c 
c        this subroutine subdivides the input cube
c        with the center cent and side a eight cubes
c 
         aout=a/2
         b=aout/2
c 
         centout(1,1)=cent(1)+b
         centout(2,1)=cent(2)+b
         centout(3,1)=cent(3)+b
c 
         centout(1,2)=cent(1)-b
         centout(2,2)=cent(2)+b
         centout(3,2)=cent(3)+b
c 
         centout(1,3)=cent(1)+b
         centout(2,3)=cent(2)-b
         centout(3,3)=cent(3)+b
c 
         centout(1,4)=cent(1)+b
         centout(2,4)=cent(2)+b
         centout(3,4)=cent(3)-b
c 
c 
         centout(1,5)=cent(1)-b
         centout(2,5)=cent(2)-b
         centout(3,5)=cent(3)-b
c 
         centout(1,6)=cent(1)+b
         centout(2,6)=cent(2)-b
         centout(3,6)=cent(3)-b
c 
         centout(1,7)=cent(1)-b
         centout(2,7)=cent(2)+b
         centout(3,7)=cent(3)-b
c 
         centout(1,8)=cent(1)-b
         centout(2,8)=cent(2)-b
         centout(3,8)=cent(3)+b
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cubegau1(n,zs,whts,cent,a,fun,par1,par2,rint)
        implicit real *8 (a-h,o-z)
        save
        dimension cent(3),zs(3,n,n,n),whts(n,n,n)
c 
c       convert the discretization of the standard cube
c       into the discretization of cube with side a and
c       center cent
c 
        cc=a/2
        ccc=cc**3
c 
        rint=0
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
        rint=rint+w*f
c 
 1600 continue
 1800 continue
 2000 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cubegau0(n,rnodes,weights)
        implicit real *8 (a-h,o-z)
        save
        dimension ts(100),ws(100),
     1      rnodes(3,n,n,n),weights(n,n,n)
c 
c       construct the Gaussian nodes and weights on the
c       interval [-1,1]
c 
        call cubegaus0(n,ts,ws)
c 
c       construct the tensor product discretization of the cube
c 
        do 2000 i=1,n
        do 1800 j=1,n
        do 1600 k=1,n
c 
        rnodes(1,i,j,k)=ts(i)
        rnodes(2,i,j,k)=ts(j)
        rnodes(3,i,j,k)=ts(k)
c 
        weights(i,j,k)=ws(i)*ws(j)*ws(k)
 1600 continue
 1800 continue
 2000 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cubegaus0(n,ts,whts)
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
        ifout=0
c 
        do 1400 k=1,10
        call cubeadap_legepol(xk,n,pol,der,sum)
        delta=-pol/der
c 
        xk=xk+delta
c 
        if(delta .lt. eps) ifout=ifout+1
        if(ifout .eq. 3) goto 1600
c 
 1400 continue
 1600 continue
        ts(i)=xk
        ts(n-i+1)=-xk
 2000 continue
c 
c        construct the weights via the orthogonality relation
c 
        do 2400 i=1,(n+1)/2
        call cubeadap_legepol(ts(i),n,pol,der,sum)
        whts(i)=1/sum
        whts(n-i+1)=whts(i)
 2400 continue
  
        return
        end
c 
c 
c 
c 
c 
        subroutine cubeadap_legepol(x,n,pol,der,sum)
        implicit real *8 (a-h,o-z)
c 
        save
        done=1
        sum=0
c 
        pkm1=1
        pk=x
        sum=sum+pkm1**2 /2
        sum=sum+pk**2 *(1+done/2)
c 
        pk=1
        pkp1=x
c 
c        if n=0 or n=1 - exit
c 
        if(n .ge. 2) goto 1200
  
        sum=0
c 
        pol=1
        der=0
        sum=sum+pol**2 /2
        if(n .eq. 0) return
c 
        pol=x
        der=1
        sum=sum+pol**2*(1+done/2)
        return
 1200 continue
c 
c       n is greater than 1. conduct recursion
c 
        do 2000 k=1,n-1
        pkm1=pk
        pk=pkp1
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
        sum=sum+pkp1**2*(k+1+done/2)
 2000 continue
c 
c        calculate the derivative
c 
        pol=pkp1
        der=n*(x*pkp1-pk)/(x**2-1)
        return
        end
