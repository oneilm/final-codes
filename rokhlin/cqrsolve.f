        implicit real *8 (a-h,o-z)
        dimension a(1000 000),b(1000 000),a2(1000 000)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
         PRINT *, 'ENTER n'
         READ *,n
         CALL PRINF('n=*',n,1 )
c 
        call testit(a,n,b,a2)
  
c 
        stop
        end
c 
c 
c 
c 
c 
        subroutine testit(a,n,b,a2)
        implicit complex *16 (a-h,o-z)
        save
        complex *16 a(n,n),b(n,n),ima,
     1      x(1000),y(1000),x2(1000),a2(1),diffs(1000),
     2      y2(1000),xs7(10000),y7(10000),w(1000 000),
     3      uu(1000 000),vv(1000 000),s(1000),
     4      ainv(1000000)
c 
        real *8 clotatim,t1,t2,t3,t4
c 
        data ima/(0.0d0,1.0d0)/
c 
c        create the test matrix
c 
        do 1400 i=1,n
        do 1200 j=1,n
c 
        a(j,i)=i**2-j*ima -j**3
cccc        a(j,i)=i**2-j**3
  
cccc        a(j,i)=0
 1200 continue
  
        a(i,i)=a(i,i)+i**3 + ima*10000
  
cccc        a(i,i)=1
  
 1400 continue
  
        call arrmove(a,a2,n*n*2)
c 
        call prin2('test matrix as created *',a,n*n*2)
c 
c       construct its QR decomposition
c 
  
  
        t1=clotatim()
  
        call cqrdecom(a,n,b,rcond)
  
  
        t2=clotatim()
  
        call prin2('cpu time for qrdecom=*',t2-t1,1)
  
        call prin2('after cqrdecom, rcond=*',rcond,1)
c 
c        construct the inverse of a
c 
  
        call cqrinv(n,b,ainv)
  
  
        call prin2('after cqrinv, ainv=*',ainv,n*n*2)
        call prinf('and n=*',n,1)
  
        call matmult(a,ainv,vv,n)
  
        call prin2('and ainv a=*',vv,n*n*2)
  
        stop
        end
  
c 
c 
c 
c 
c 
        subroutine matmult(a,x,c,n)
        implicit complex *16 (a-h,o-z)
        save
        dimension a(n,n),x(n,n),c(n,n)
c 
        do 1600 i=1,n
        do 1400 j=1,n
c 
        cd=0
        do 1200 ll=1,n
c 
        cd=cd+a(i,ll)*x(ll,j)
cccc        cd=cd+a(i,ll)*x(j,ll)
 1200 continue
c 
        c(i,j)=cd
 1400 continue
 1600 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine arrmove(x,y,n)
        implicit complex *16 (a-h,o-z)
c 
        save
        dimension x(1),y(1)
c 
        do 1200 i=1,n
c 
        y(i)=x(i)
 1200 continue
c 
        return
        end
c 
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code and the beginning of the
c        linear solver proper
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c        This file contains three user-callable subroutines: cqrsolv,
c        cqrdecom, cqrsolve. Following is a brief description of these
c        subroutines.
c 
c   cqrsolv - uses a version of QR-decomposition to solve the user-supplied
c       system of complex linear algebraic equations, with the matrix a and
c       right-hand side rhs. Both the matrix a and the right-hand side rhs
c       are destroyed in the process. This is a primitive, short, and
c       reasonably efficient subroutine. IMPORTANT NOTE: THIS IS A
c       PRIMITIVE ROUTINE IN THAT IT SOLVES THE SYSTEM WITH A SINGLE
c       RIGHT-HAND SIDE, AND WILL PERFORM THE QR DECOMPOSITION AGAIN AND
c       AGAIN FOR EACH NEW RIGHT-HAND SIDE. IT SHOULD NOT BE USED IN THIS
c       REGIME! Instead, the subroutines qrdecom, qrsolve (see) should be
c       used! On the other hand, for a single right-hand side, this routine
c       is about twice more efficient than the subroutine qrdecom.
c 
c   cqrdecom - subroutine constructs a QR-decomposition of the real
c       user-supplied matrix a. It is expected that this decomposition
c       will be used by the subroutine qrsolve (see) for the solution
c       of linear systems with the matrix a; this subroutine has no
c       known uses as a stand-alone device.
c 
c   cqrsolve  - uses a version of QR-decomposition to solve the
c       user-supplied system of linear algebraic equations. The
c       QR-decomposition is prepared by a prior call to the subroutine
c       qrdecom (see); this subroutine has no known uses as a stand-alone
c       device.
c 
c   cqrinv - uses the QR decomposition obtained via a preceding call
c       to the subroutine cqrdecom (see) to construct the inverse of
c       the user-supplied matrix A.
c 
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        subroutine cqrinv(n,b,ainv)
        implicit complex *16 (a-h,o-z)
        save
        dimension b(1),ainv(1)
c 
c        This subroutine uses the QR decomposition obtained via
c        a preceding call to the subroutine cqrdecom (see)
c        to construct the inverse of the user-supplied matrix A.
c        Please note that the only inputs used by this subroutine
c        are the dimensionality n of the matrix A, and the array
c        b thas has been produced during a preceding call to
c        cqrdecom.
c 
c 
c 
c                      Input parameters:
c 
c  n - the dimensionality of the problem
c  b - the factorization of the matrix a (hopefully) produced via
c       a preceding call to the subroutine qrdecom (see)
c 
c                      Output parameters:
c 
c  ainv - the inverse of the matrix A
c 
c 
c        . . . construct the memory map
c 
        ia=1
        la=n*n+2
c 
        ib=ia+la
        lb=n*n+2
c 
        iz=ib+lb
        lz=n+2
c 
c       . . . apply the inverse of the factored matrix to the
c             right-hand side
c 
        call cqrinv0(n,b(ia),b(ib),ainv,b(iz) )
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine cqrinv0(n,u,b,ainv,z)
        implicit complex *16 (a-h,o-z)
        save
        dimension u(n,n),b(n,n),ainv(n,n),z(1)
c 
c        construct the inverse of a
c 
        do 2000 i=1,n
c 
        do 1200 j=1,n
        z(j)=b(i,j)
 1200 continue
c 
        call cqrtrin(u,z,n)
c 
        do 1400 j=1,n
        ainv(i,j)=z(j)
 1400 continue
c 
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine cqrcond(n,a,b,rcond,w)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,n),b(n,n),w(1)
c 
c       This subroutine uses a version of QR-decomposition
c       (obtained via a preceding call to the subroutine
c       cqrdecom) to estimate the condition number of the
c       complex matrix a
c 
c                      Input parameters:
c 
c  n - the dimensionality of the problem
c  a - the matrix whose condition number; must be the same as that
c       from which the subroutine cqrdecom has produced the paramater
c       b in the calling sequence of this subroutine
c  b - the factorization of the matrix a (hopefully) produced via
c       a preceding call to the subroutine qrdecom (see)
c 
c                      Output parameters:
c 
c  rcond - a (fairly good) estimate of the condition number of the
c       matrix a
c 
c                      Work arrays:
c 
c  w - must be at least 4*n+4 real *8 locations in length
c 
c 
c        . . . initialize the array to be used for forward iteration
c 
        call cqrsolve_rand(n*2,w)
        call cqrsolve_rand(n*2,w)
c 
c       estimate the maximum singular value of the matrix a
c 
        numit=10
        do 1800 i=1,numit
c 
        iy=n+1
        call cqrmatve(a,n,w,w(iy))
        call cqrmatve_adj(a,n,w(iy),w)
c 
        d=0
        do 1200 j=1,n
c 
        d=d+w(j)*conjg(w(j))
 1200 continue
c 
        d=sqrt(d)
c 
        do 1400 j=1,n
c 
        w(j)=w(j)/d
 1400 continue
c 
 1800 continue
c 
 2000 continue
c 
        rlammax=sqrt(d)
c 
c       estimate the minimum singular value of the matrix a
c 
        call cqrsolve_rand(n*2,w)
c 
        do 2800 i=1,numit
c 
        iy=n+1
c 
        call cqrsolve(n,b,w,w(iy))
        call cqrsolve_adj(n,b,w(iy),w)
c 
        d=0
        do 2200 j=1,n
c 
        d=d+w(j)*conjg(w(j))
 2200 continue
c 
        d=sqrt(d)
        do 2400 j=1,n
c 
        w(j)=w(j)/d
 2400 continue
c 
 2800 continue
c 
        rlammin=1/sqrt(d)
        rcond=rlammax/rlammin
c 
        return
        end
c 
c 
c 
c 
c 
c 
        subroutine cqrsolv(a,n,rhs,rcond)
        implicit complex *16 (a-h,o-z)
        save
        complex *16 a(n,n),u(2,2),aa(2),rhs(1)
c 
        real *8 dmax,dmin,rcond,size22
c 
c       This subroutine uses a version of QR-decomposition to solve
c       the user-supplied system of linear algebraic equations, with
c       the complex matrix a and right-hand side rhs. Both the matrix
c       a and the right-hand side rhs are destroyed in the process.
c 
c       IMPORTANT NOTE:
C        THIS IS A PRIMITIVE ROUTINE IN THAT IT SOLVES THE SYSTEM
C        WITH A SINGLE RIGHT-HAND SIDE, AND WILL PERFORM THE
C        QR DECOMPOSITION AGAIN AND AGAIN FOR EACH NEW RIGHT-HAND
C        SIDE. IT SHOULD NOT BE USED IN THIS REGIME! Instead, the
c        subroutines cqrdecom, cqrsolve (see) should be used
c 
c                    Input parameters:
c 
c  n - the dimensionality of the system being solved
c  b - the array containing the QR-decomposition of the matrix of the
c       system (hopefully, it has been produced by a prior call to the
c       subroutine qrdecom)
c  y - the right-hand side of the system to be solved; not damaged by this
c       subroutine in any way
c 
c                    Output parameters:
c 
c  x - the solution of the system
c  rcond - a fairly crude estimate of the condition number of a
c 
c        . . . transpose the input matrix a
c 
        size22=0
        do 1400 i=1,n
        do 1200 j=1,i
c 
        d=a(j,i)
        a(j,i)=a(i,j)
        a(i,j)=d
c 
        size22=size22+a(j,i)*conjg(a(j,i))
        size22=size22+a(i,j)*conjg(a(i,j))
 1200 continue
 1400 continue
c 
c       eliminate the upper right triangle
c 
        do 2000 i=1,n-1
c 
        do 1600 j=n,i+1,-1
c 
        aa(1)=a(i,j-1)
        aa(2)=a(i,j)
        call cqrrotfn2(aa,u,size22)
c 
        call cqrrotat(u,a(1,j-1),a(1,j),n,i)
c 
        d1=u(1,1)*rhs(j-1)+u(1,2)*rhs(j)
        d2=u(2,1)*rhs(j-1)+u(2,2)*rhs(j)
        rhs(j-1)=d1
        rhs(j)=d2
c 
 1600 continue
 2000 continue
c 
c       estimate the condition number
c 
        dmax=-1
        dmin=abs(a(1,1))
        do 2200 i=1,n
c 
        if(dmax .lt. abs(a(i,i)) ) dmax=abs(a(i,i))
        if(dmin .gt. abs(a(i,i)) ) dmin=abs(a(i,i))
 2200 continue
c 
        if(dmin .eq. 0) then
            rcond=-1
            return
        endif
c 
        rcond=dmax/dmin
c 
        call cqrtrin(a,rhs,n)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cqrsolve(n,b,y,x)
        implicit complex *16 (a-h,o-z)
        save
        dimension b(1),x(1),y(1)
c 
c       This subroutine uses a version of QR-decomposition to solve
c       the user-supplied system of complex linear algebraic equations.
c       The QR-decomposition is prepared by a prior call to the subroutine
c       cqrdecom (see); this subroutine has no known uses as a stand-alone
c       device.
c 
c                    Input parameters:
c 
c  n - the dimensionality of the system being solved
c  b - the array containing the QR-decomposition of the matrix of the
c       system (hopefully, it has been produced by a prior call to the
c       subroutine cqrdecom)
c  y - the right-hand side of the system to be solved; not damaged by this
c       subroutine in any way
c 
c                    Output parameters:
c 
c  x - the solution of the system
c 
c 
c        . . . construct the memory map
c 
        ia=1
        la=n*n+2
c 
        ib=ia+la
        lb=n*n+2
c 
        iz=ib+lb
        lz=n+2
c 
c       . . . apply the inverse of the factored matrix to the
c             right-hand side
c 
        call cqrtrinv(b(ia),y,b(iz),n)
        call cqrmatve(b(ib),n,b(iz),x)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cqrsolve_adj(n,b,y,x)
        implicit real *8 (a-h,o-z)
        save
        complex *16 b(1),x(1),y(1)
c 
c       This subroutine uses a version of QR-decomposition to
c       solve the user-supplied system of complex linear
c       algebraic equations
c 
c                   A^* (x)=y.                                         (1)
c 
c       Please note that A^* in (1) above denotes the adjoint of
c       the  matrix A. Also, please note that the QR-decomposition
c       is prepared by a prior call to the subroutine cqrdecom (see);
c       this subroutine has no known uses as a stand-alone device.
c       Also, please note that this subroutine is the companion to
c       the subroutine cqrsolve (see), solving the system of linear
c       equations
c 
c                   A(x)=y.                                            (2)
c 
c                    Input parameters:
c 
c  n - the dimensionality of the system being solved
c  b - the array containing the QR-decomposition of the matrix of the
c       system (hopefully, it has been produced by a prior call to the
c       subroutine cqrdecom)
c  y - the right-hand side of the system to be solved; not damaged by this
c       subroutine in any way
c 
c                    Output parameters:
c 
c  x - the solution of the system
c 
c 
c        . . . construct the memory map
c 
        ia=1
        la=n*n+2
c 
        ib=ia+la
        lb=n*n+2
c 
c        estimate the condition number of a
c 
        call cqrmatve_adj(b(ib),n,y,x)
c 
        call cqrtrin_conjg(b(ia),x,n)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cqrdecom(a,n,b,rcond)
        implicit complex *16 (a-h,o-z)
        save
        dimension a(1),b(1)
c 
c       This subroutine constructs a QR-decomposition of the complex
c       user-supplied matrix a. It is expected that this decomposition
c       will be used by the subroutine cqrsolve (see) for the solution
c       of linear systems with the matrix a; this subroutine has no
c       known uses as a stand-alone device.
c 
c                    Input parameters:
c 
c  a - the matrix to be QR-decomposed; not damaged by the subroutine
c       in any way
c  n - the dimensionality of the matrix a
c 
c                    Output parameters:
c 
c  b - the array containing the QR-decomposition of a; must be
c       2*n**2+n+10 real *8 locations long
c  rcond - a fairly crude estimate of the condition number of a
c 
c        . . . allocate memory
c 
        ia=1
        la=n*n+2
c 
        ib=ia+la
        lb=n*n+2
c 
        do 1200 i=1,n*n
c 
        b(i)=a(i)
 1200 continue
c 
c       construct the QR-decomposition
c 
        call cqrelim(b,n,b(ib),rcond)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cqrtrin_conjg(a,y,n)
        implicit complex *16 (a-h,o-z)
        save
        complex *16 a(n,n),y(1)
c 
c       apply the inverse of the trizngular matrix a to y
c 
        y(n)=y(n)/conjg(a(n,n))
        do 1400 i=n-1,1,-1
c 
        d=0
        do 1200 j=n,i+1,-1
c 
        d=d+conjg(a(j,i))*y(j)
 1200 continue
c 
        y(i)=(y(i)-d)/conjg(a(i,i))
 1400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine cqrtrin(a,y,n)
        implicit complex *16 (a-h,o-z)
        save
        complex *16 a(n,n),y(1)
c 
c       apply the inverse of the trizngular matrix a to y
c 
        y(n)=y(n)/a(n,n)
        do 1400 i=n-1,1,-1
c 
        d=0
        do 1200 j=n,i+1,-1
c 
        d=d+a(j,i)*y(j)
 1200 continue
c 
        y(i)=(y(i)-d)/a(i,i)
 1400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine cqrtrinv(a,y,x,n)
        implicit complex *16 (a-h,o-z)
        save
        dimension a(n,n),x(1),y(1)
c 
c       start the process
c 
        x(1)=y(1)/a(1,1)
        do 1400 i=2,n
c 
        d=0
        do 1200 j=1,i-1
c 
        d=d+a(i,j)*x(j)
 1200 continue
c 
        x(i)=(y(i)-d)/a(i,i)
 1400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine cqrelim(a,n,b,rcond)
        implicit complex *16 (a-h,o-z)
        save
        complex *16 a(n,n),u(2,2),aa(2),b(n,n)
c 
        real *8 dmax,dmin,rcond,size22
c 
c       construct the unity matrix
c 
        size22=0
        do 1400 i=1,n
        do 1200 j=1,n
c 
        b(j,i)=0
c 
        size22=size22+a(j,i)*conjg(a(j,i))
 1200 continue
c 
        b(i,i)=1
 1400 continue
c 
c       eliminate the upper right triangle
c 
        do 2000 i=1,n-1
c 
        do 1600 j=n,i+1,-1
c 
        aa(1)=a(i,j-1)
        aa(2)=a(i,j)
cccc        call cqrrotfn2(aa,u,size22)
        call cqrrotfn3(aa,u,size22)
  
  
cccc        call prin2('u=*',u,8)
c 
        call cqrrotat2(u,a(1,j-1),a(1,j),n,i)
c 
        ii=j-i
        call cqrrotat2(u,b(1,j-1),b(1,j),n,ii)
 1600 continue
 2000 continue
c 
c       estimate the condition number
c 
        dmax=-1
        dmin=abs(a(1,1))
        do 2200 i=1,n
c 
        if(dmax .lt. abs(a(i,i)) ) dmax=abs(a(i,i))
        if(dmin .gt. abs(a(i,i)) ) dmin=abs(a(i,i))
 2200 continue
c 
        if(dmin .eq. 0) then
            rcond=-1
            return
        endif
c 
        rcond=dmax/dmin
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cqrrotat2(a,x,y,n,n0)
        implicit complex *16 (a-h,o-z)
        save
        complex *16 a(2,2),x(1),y(1),d1,d2
        real *8 a11,a22
c 
        a11=a(1,1)
        a22=a(2,2)
c 
        do 1200 i=n0,n
c 
        d1=a11*x(i)+a(1,2)*y(i)
        d2=a(2,1)*x(i)+a22*y(i)
c 
cccc        d1=a(1,1)*x(i)+a(1,2)*y(i)
cccc        d2=a(2,1)*x(i)+a(2,2)*y(i)
c 
        x(i)=d1
        y(i)=d2
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cqrrotat(a,x,y,n,n0)
        implicit complex *16 (a-h,o-z)
        save
        complex *16 a(2,2),x(1),y(1),d1,d2
c 
        do 1200 i=n0,n
c 
        d1=a(1,1)*x(i)+a(1,2)*y(i)
        d2=a(2,1)*x(i)+a(2,2)*y(i)
c 
        x(i)=d1
        y(i)=d2
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cqrrotfn3(a,u,size22)
        implicit complex *16 (a-h,o-z)
        save
        dimension a(2),u(2,2)
        real *8 size22,d
c 
        u21=-a(2)
        u22=a(1)
c 
        d=u22*conjg(u22)+u21*conjg(u21)
c 
        if(d .lt. size22*1.0d-66) then
c 
            u(2,2)=1
            u(1,2)=0
            u(1,1)=1
            u(2,1)=0
            return
        endif
c 
        d=sqrt(d)
        u(2,2)=u22/d
        u(2,1)=u21/d
c 
        u(1,1)=-conjg(u(2,2))
        u(1,2)=conjg(u(2,1))
  
c 
        if(u(1,1) .ne. 0) then
            cd=abs(u(1,1))/u(1,1)
            u(1,1)=u(1,1)*cd
            u(1,2)=u(1,2)*cd
        endif
c 
c 
        if(u(2,2) .ne. 0) then
            cd=abs(u(2,2))/u(2,2)
            u(2,1)=u(2,1)*cd
            u(2,2)=u(2,2)*cd
        endif
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cqrrotfn2(a,u,size22)
        implicit complex *16 (a-h,o-z)
        save
        dimension a(2),u(2,2)
        real *8 size22,d
c 
        u21=-a(2)
        u22=a(1)
c 
        d=u22*conjg(u22)+u21*conjg(u21)
c 
        if(d .lt. size22*1.0d-66) then
c 
            u(2,2)=1
            u(1,2)=0
            u(1,1)=1
            u(2,1)=0
            return
        endif
c 
        d=sqrt(d)
        u(2,2)=u22/d
        u(2,1)=u21/d
c 
        u(1,1)=-conjg(u(2,2))
        u(1,2)=conjg(u(2,1))
        return
        end
c 
c 
c 
c 
c 
        subroutine cqrmatve(a,n,x,y)
        implicit complex *16 (a-h,o-z)
        save
        dimension a(n,n),x(1),y(1)
c 
c        apply the matrix a to the vector x obtaining y
c 
        do 1400 i=1,n
        cd=0
        do 1200 j=1,n
        cd=cd+a(i,j)*x(j)
 1200 continue
        y(i)=cd
 1400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine cqrmatve_adj(a,n,x,y)
        implicit complex *16 (a-h,o-z)
        save
        dimension a(n,n),x(1),y(1)
c 
c        apply the matrix a to the vector x obtaining y
c 
        do 1400 i=1,n
        cd=0
        do 1200 j=1,n
        cd=cd+conjg(a(j,i))*x(j)
 1200 continue
        y(i)=cd
 1400 continue
        return
        end
  
  
c 
c 
c 
c 
c 
        subroutine cqrsolve_rand(n,y)
        implicit real *8 (a-h,o-z)
        save
        dimension y(1)
        data ifcalled/0/
c 
c       generate pseudo-random numbers
c 
        if(ifcalled .ne. 0) goto 1100
c 
        done=1
        pi=atan(done)*4
        add=sqrt(7*done)
        x=1.1010101010101
        ifcalled=1
 1100 continue
c 
        do 1200 i=1,n
c 
        phi=x*100000*pi+add
        j=phi
        phi=(phi-j)
        x=phi
c 
        y(i)=x
 1200 continue
c 
        return
c 
c 
c 
c 
        entry qrsolve_rand_reset(x7)
c 
        done=1
        pi=atan(done)*4
        add=sqrt(7*done)
        x=x7
        ifcalled=1
        return
        end
