        IMPLICIT REAL *16 (A-H,O-Z)
         complex *32 a(2,2),ima,u(2,2),v(2,2),
     1       rlams(2,2),w1(2,2),w2(2,2),cb(20000),
     2       uu(20000),vv(20000),cb2(20000),cb3(20000),
     3       ww(20000)
         dimension rlam(2),b(20000),b2(20000)
         data ima/(0.0d0,1.0d0)/
c 
         CALL PRINI(6,13)
C 
C 
C       SET ALL PARAMETERS
C 
         PRINT *, 'ENTER n'
         READ *,n
         CALL PRINF('n=*',n,1 )
cccc
c 
cccc         ima=0
         a(1,1)=-(1+ima)
        a(2,2)=2-ima/3
        a(1,2)=3+ima*2
        a(2,1)=-(4-ima*3)
c 
ccc        a(2,2)=1.0d-17
        a(1,2)=1.0d-17
         call prinq('a as created*',a,8)
c 
         eps=1.0d-30
         call qoneuv(ier,a,rlam,u,v)
cccc         call prinq('after qoneuv, u=*',u,8)
cccc         call prinq('after qoneuv, v=*',v,8)
cccc         call prinq('after qoneuv, rlam=*',rlam,2)
c 
c       check the diagonalizing matrices
c 
        rlams(1,1)=rlam(1)
        rlams(2,2)=rlam(2)
        rlams(1,2)=0
        rlams(2,1)=0
        call d2mproqc(u,a,w1)
        call d2mproqc(w1,v,w2)
        call prinq('u a v =*',w2,8)
c 
c 
c       construct the big test matrix
c 
        call matrcre(b,n)
        call RIRAND(N**2*3,b)
        do 1200 i=1,n**2
        cb(i)=b(i)+ima*b(i*2+1)
cccc        cb(i)=b(i)
cccc        if(i .lt. n**2/2) cb(i)=0
        cb2(i)=cb(i)
 1200 continue
        call prinq('matrix as created*',cb,n*n*2)
c 
c       . . . and construct its svd
c 
        eps=1.0d-30
cccc        t1= second(tt)
        ifuv=1
        call qd2mcudv(cb,n,uu,vv,eps,ifuv,numrows)
ccc        call d2mudv(b,n,uuu,vvv,eps,ifuv,numrows)
cccc        t2= second(tt)
         call prinq('and CPU time is*',t2-t1,1)
         call prinf('and numrows=*',numrows,1)
        call prinq('singular values after qd2mcudv*',cb,n*2)
cccc        call prinq('matrix after d2mudv*',cb,n*n*2)
cccc        call prinq('uu after d2mudv*',uu,n*n*2)
cccc        call prinq('vv after d2mudv*',vv,n*n*2)
c 
c       now, test the svd
c 
         call svdtes(cb,uu,vv,ww,cb2,cb3,n)
  
         stop
         end
c 
c 
c 
c 
c 
        subroutine svdtes(s,u,v,w,b,w2,n)
        implicit real *16 (a-h,o-z)
        save
        complex *32 s(1),u(n,n),v(n,n),w(n,n),b(n,n),
     1      w2(n,n)
c 
c       construct the matrix of singular values
c 
        do 1400 i=1,n
        do 1200 j=1,n
        w(j,i)=0
 1200 continue
        w(i,i)=s(i)
 1400 continue
cccc         call prinq('in svdtes, w=*',w,n*n*2)
c 
c       construct all proper products
c 
        call matmul(u,w,w2,n)
         call prinq('after first matmul, w2=*',w2,8)
        call matmul(w2,v,w,n)
cccc          call prinq('and b recomputed from svd*',w,n*n*2)
cccc          call prinq('while b *',b,n*n*2)
c 
c        now, calculate the error of the decomposition
c 
        dd=0
        do 1800 i=1,n
        do 1600 j=1,n
        dd=dd+(b(j,i)-w(j,i))*qconjg(b(j,i)-w(j,i))
 1600 continue
 1800 continue
        dd=qsqrt(dd)
         call prinq('and l2-error is*',dd,1)
        return
        end
c 
c 
c 
c 
c 
        SUBROUTINE RIRAND(N,Y)
        IMPLICIT REAL *16 (A-H,O-Z)
        save
        DIMENSION Y(1)
        DATA IFCALL /0/
C 
C       CONDUCT PRELIMINARY RANDOMIZATION
C 
        IF (IFCALL.EQ.1 ) GOTO  1300
        IFCALL =1
CCCC    LAMBDA=13
        LAMBDA=2**7 +9
        MU=510511
        MU=1
CCC     IP=2**10
        IP=2**20
        D=1
        D=D/IP
        M=17
        DO 1200 I=1,100
        M1=M*LAMBDA  +MU
        J=M1/IP
        M=M1-J*IP
CCCC    CALL PRINF('M=*',M,1)
 1200 CONTINUE
CCC     CALL PRINF('PRELIMINARY RANDOMIZATION DONE, M=*',M,1)
 1300 CONTINUE
C 
C       GENERATE PSEUDO-RANDOM NUMBERS
C 
        DO 1400 I=1,N
        M1=M*LAMBDA +MU
        J=M1/IP
        M=M1-J*IP
        Y(I)=M*D
 1400 CONTINUE
        RETURN
c 
c 
c 
c 
        entry rirandin(ifcall7)
        ifcall=0
        return
        END
c 
c 
c 
c 
c 
        subroutine matmul(a,b,c,n)
        implicit real *16 (a-h,o-z)
        save
        complex *32 a(n,n),b(n,n),c(n,n),d
        do 2000 i=1,n
        do 1800 j=1,n
        d=0
        do 1400 k=1,n
        d=a(i,k)*b(k,j)+d
 1400 continue
        c(i,j)=d
 1800 continue
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine matrcre(a,n)
        implicit real *16 (a-h,o-z)
        save
        dimension a(n,n)
c 
        done=1
        do 1400 i=1,n
        do 1200 j=1,n
        a(j,i)=2*i**2+j**2-7*sin(1.2*j)
        a(i,j)=done/(i+j)
 1200 continue
 1400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine matrcre2(a,n)
        implicit real *16 (a-h,o-z)
        save
        dimension a(n,n)
c 
        do 1400 i=1,n
        do 1200 j=1,n
        a(j,i)=0
 1200 continue
        a(i,i)=1
 1400 continue
c 
        do 1600 i=1,n-1
        a(i,i+1)=1
 1600 continue
        return
        end
c 
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c            this is the end of the debuging code and the beginning
c            of the actual SVD routines
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c 
        subroutine qd2mcudv(a,n,u,v,eps,ifuv,numrows)
        implicit complex *32 (a-h,o-z)
        save
        dimension a(n,n),u0(2,2),v0(2,2),u(n,n),v(n,n),
     1      b(2,2),rlam(2),u1(2,2),v1(2,2)
        real *16 rlam(2),eps,thresh
c 
c        this subroutine uses the classical jacobi iteration
c        to construct the singular value decomposition of a
c        general complex *32 matrix, of the form
c 
c          a = u s v,                                        (1)
c 
c        with u,v unitary matrices, and s diagonal, with
c        the diagonal elements monotonically decreasing.
c 
c                     input parameters:
c 
c  a - the n * n matrix to be sv-decomposed, is destroyed
c        by the subroutine
c  n - the dimensionality of a
c  eps - the accuracy to which the computations will be performed
c  ifuv - tells the subroutine whether the singular vectors are
c        to be computed, or only singular values
c 
c                     output parameters:
c 
c  a - the first column of a contains the (ordered) singular values ,
c         {\bf IN COMPLEX FORM}
c  u,v - unitary matrices in (1)
c  numrows - the number of jacobi iterations (i.e. single rotations)
c        that the algorithm has taken
c 
c       . . . one sweep after another, conduct jacobi
c             iterations
c 
        maxsweep=1000
        if(ifuv .eq. 0) goto 1150
        zero=0
        done=1
        do 1100 i=1,n
        do 1050 j=1,n
        u(j,i)=zero
        v(j,i)=zero
 1050 continue
        u(i,i)=done
        v(i,i)=done
 1100 continue
 1150 continue
        numrows=0
        thresh=qsqrt(eps)
        thresh=qsqrt(thresh)
        thresh=qsqrt(thresh)
        do 3400 iijjkk=1,4
c 
        do 3000 ijk=1,maxsweep
        ifany=0
c 
        do 2000 i=1,n
        do 1800 j=1,n
         if(j .eq. i) goto 1800
         if( (cqabs(a(i,j)) .lt. thresh) .and.
     1       (cqabs(a(j,i)) .lt. thresh) ) goto 1800
        ifany=1
        numrows=numrows+1
c 
c       construct the svd of the (i,j)-th 2 * 2 submatrix
c 
        b(1,1)=a(i,i)
        b(1,2)=a(i,j)
        b(2,1)=a(j,i)
        b(2,2)=a(j,j)
        call qoneuv(ier,b,rlam,u0,v0)
c 
c     transform the big matrix
c 
c       . . . the rows
c 
        do 1200 l=1,n
        di=u0(1,1)*a(i,l)+u0(1,2)*a(j,l)
        dj=u0(2,1)*a(i,l)+u0(2,2)*a(j,l)
        a(i,l)=di
        a(j,l)=dj
 1200 continue
c 
c       . . . and the columns
c 
        do 1400 l=1,n
        di=v0(1,1)*a(l,i)+v0(2,1)*a(l,j)
        dj=v0(1,2)*a(l,i)+v0(2,2)*a(l,j)
        a(l,i)=di
        a(l,j)=dj
 1400 continue
c 
c       if the user so requested, adjust the matrices u,v
c 
        if(ifuv .eq. 0) goto 1800
c 
        call d2mcpqon(u0,u1)
        call d2mcpqon(v0,v1)
        do 1600 l=1,n
c 
c       the matrix u
c 
        di=u1(1,1)*u(l,i)+u1(2,1)*u(l,j)
        dj=u1(1,2)*u(l,i)+u1(2,2)*u(l,j)
        u(l,i)=di
        u(l,j)=dj
c 
c       . . . and the columns
c 
        di=v1(1,1)*v(i,l)+v1(1,2)*v(j,l)
        dj=v1(2,1)*v(i,l)+v1(2,2)*v(j,l)
        v(i,l)=di
        v(j,l)=dj
 1600 continue
c 
 1800 continue
 2000 continue
c 
         call prinf('ifany=*',ifany,1)
      if(ifany .eq. 0) goto 3200
 3000 continue
 3200 continue
        thresh=thresh**2
 3400 continue
c 
c        now, reorder the singular values of a to put them
c        in increasing order
c 
        do 4200 i=1,n
        a(i,1)=a(i,i)
 4200 continue
c 
       call qd2mcsrt(a(1,1),n,u,v,a(1,2),ifuv)
        return
        end
c 
c 
c 
c 
c 
        subroutine qd2mcsrt(s,n,u,v,iw,ifuv)
        implicit real *16 (a-h,o-z)
        save
        complex *32 u(n,n),v(n,n),s(1),com,cd
        dimension rea(2),iw(1)
        equivalence (rea(1),com)
c 
c        if some of the alleged singular values
c        are negative - change their sign, and that
c        of the corresponding left eigenvectors
c 
        do 1100 i=1,n
        com=s(i)
        ds=rea(1)
        if(ds .ge. 0) goto 1100
        s(i)=-s(i)
        if(ifuv .eq. 0) goto 1100
        do 1050 j=1,n
        v(i,j)=-v(i,j)
 1050 continue
 1100 continue
  
c        using the bubble sort, reorder the elements
c        of array s
c 
        do 1200 i=1,n
        iw(i)=i
 1200 continue
c 
        do 2000 i=1,n-1
c 
        ifmoved=0
        do 1600 j=1,n-i
        com=s(j)
        dsj=rea(1)
        com=s(j+1)
        dsjp1=rea(1)
        if(dsj. ge. dsjp1) goto 1600
        jj=iw(j)
        iw(j)=iw(j+1)
        iw(j+1)=jj
c 
        cd=s(j)
        s(j)=s(j+1)
        s(j+1)=cd
 1600 continue
 2000 continue
        if(ifuv .eq. 0) return
c 
c       now, reorder the rows of u and columns of v
c 
        do 2400 i=1,n-1
        ii=iw(i)
c 
c       exchange the rows number i and ii in the matrix u
c 
        do 2200 j=1,n
        cd=v(i,j)
        v(i,j)=v(ii,j)
        v(ii,j)=cd
c 
c       exchange columns number i and ii in the matrix v
c 
        cd=u(j,i)
        u(j,i)=u(j,ii)
        u(j,ii)=cd
 2200 continue
c 
        do 2300 j=i+1,n
        if(iw(j) .ne. i) goto 2300
        iw(j)=ii
        goto 2350
 2300 continue
 2350 continue
c 
 2400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine qoneuv(ier,a0,rlam,u,v)
        implicit real *16 (a-h,o-z)
        save
        complex *32 z,b(2,2),u(2,2),a0(2,2),
     1      uu(2,2),v(2,2),w1(2,2),w2(2,2),
     2      w3(2,2),cd,clam1,clam2,clam,x(2),y(2)
        dimension rlam(2)
c 
c       multiply a by its adjoint
c 
        call d2mcpqon(a0,w1)
        call d2mproqc(w1,a0,b)
c 
c       find the eigenvalues of a^*a
c 
        two=2
        four=4
        cd=(b(1,1)-b(2,2))**2+four*b(1,2)*b(2,1)
        cd=cqsqrt(cd)
        clam1= (b(1,1)+b(2,2)+cd )/two
        clam2= (b(1,1)+b(2,2)-cd )/two
c 
c        find the eigenvectors
c 
        clam=clam1
        dc2=clam2*qconjg(clam2)
        dc1=clam1*qconjg(clam1)
        if( dc2 .gt. dc1 ) clam=clam2
c 
        d11=(b(1,1)-clam)*qconjg((b(1,1)-clam))
        d22=(b(2,2)-clam)*qconjg((b(2,2)-clam))
        d12=b(1,2)*qconjg(b(2,2))
c 
        if( (d11 .lt. d22) .or. (d11. lt. d12) ) goto 1400
c 
        x(2)=1
        x(1)=-b(1,2)/(b(1,1)-clam)
        goto 1800
c 
 1400 continue
c 
        if(d22 .lt. d12) goto 1600
c 
        x(1)=1
        x(2)=-b(2,1)/(b(2,2)-clam)
        goto 1800
c 
 1600 continue
c 
        x(2)=1
        x(1)=(b(1,1)-clam)/b(1,2)
c 
 1800 continue
c 
        d=x(1)*qconjg(x(1))+x(2)*qconjg(x(2))
        d=qsqrt(d)
        x(1)=x(1)/d
        x(2)=x(2)/d
c 
        call d2mqort(x,y)
c 
        uu(1,1)=x(1)
        uu(2,1)=x(2)
        uu(1,2)=y(1)
        uu(2,2)=y(2)
c 
        call d2mcpqon(uu,w2)
c 
c        now, obtain v
c 
         call d2mproqc(a0,uu,v)
c 
         d1=v(1,1)*qconjg(v(1,1))+v(2,1)*qconjg(v(2,1))
         d2=v(1,2)*qconjg(v(1,2))+v(2,2)*qconjg(v(2,2))
c 
         if(d2 .gt. d1) goto 3000
c 
         d=qsqrt(d1)
         call d2mqort(v(1,1),v(1,2))
         goto 3200
c 
 3000 continue
c 
         d=qsqrt(d2)
         call d2mqort(v(1,2),v(1,1))
 3200 continue
c 
c       obtain the diagonal matrix
c 
         call d2mproqc(a0,uu,w1)
         call d2mcpqon(v,b)
         call d2mproqc(b,w1,w3)
c 
        call d2mcopq2(w2,u)
c 
c       finally, make sure that the diagonal elements are real
c 
        z=1
        dzero=0
        d=w3(1,1)*qconjg(w3(1,1))
        if(d .ne. dzero) z=w3(1,1)/cqabs(w3(1,1))
        rlam(1)=cqabs(w3(1,1))
        v(1,1)=v(1,1)*z
        v(2,1)=v(2,1)*z
c 
        z=1
        dzero=0
        d=w3(2,2)*qconjg(w3(2,2))
        if(d .ne. dzero) z=w3(2,2)/cqabs(w3(2,2))
        rlam(2)=cqabs(w3(2,2))
        v(1,2)=v(1,2)*z
        v(2,2)=v(2,2)*z
c 
c        invert u, v
c 
        call d2mcopq2(v,w3)
        call d2mcpqon(u,v)
        call d2mcpqon(w3,u)
        return
        end
c 
c 
c 
c 
c 
        subroutine d2mcopq2(a,b)
        implicit complex *32 (a-h,o-z)
        save
        dimension a(2,2),b(2,2)
c 
        b(1,1)=a(1,1)
        b(2,2)=a(2,2)
        b(1,2)=a(1,2)
        b(2,1)=a(2,1)
        return
c 
c 
c 
c 
        entry d2mcpqon(a,b)
c 
        b(1,1)=qconjg(a(1,1))
        b(2,2)=qconjg(a(2,2))
        b(2,1)=qconjg(a(1,2))
        b(1,2)=qconjg(a(2,1))
        return
        end
c 
c 
c 
c 
c 
        subroutine d2mqort(x,y)
        implicit real *16 (a-h,o-z)
        save
        complex *32 x(2),y(2)
c 
c       normalize x
c 
        d1=x(1)*qconjg(x(1))
        d2=x(2)*qconjg(x(2))
        d=d1+d2
        d=qsqrt(d)
        x(1)=x(1)/d
        x(2)=x(2)/d
c 
c       make y into a unit vector orthogonal to x
c 
c       . . . when |x(1)| > |x(2)|
c 
        if(d1 .lt. d2) goto 1200
        y(2)=x(1)
        y(1)=-qconjg(x(2))*x(1)/qconjg(x(1))
        return
 1200 continue
c 
c       . . . when |x(2)| > |x(1)|
c 
        y(1)=x(2)
        y(2)=-qconjg(x(1))*x(2)/qconjg(x(2))
        return
        end
c 
c 
c 
c 
c 
        subroutine d2mproqc(a,b,c)
        implicit complex *32 (a-h,o-z)
        save
        dimension a(2,2),b(2,2),c(2,2)
c 
        c(1,1)=a(1,1)*b(1,1)+a(1,2)*b(2,1)
        c(1,2)=a(1,1)*b(1,2)+a(1,2)*b(2,2)
        c(2,1)=a(2,1)*b(1,1)+a(2,2)*b(2,1)
        c(2,2)=a(2,1)*b(1,2)+a(2,2)*b(2,2)
        return
        end
