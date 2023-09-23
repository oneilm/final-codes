        implicit real *8 (a-h,o-z)
        dimension a(2,2),rlam(2),u(2,2),v(2,2),
     1      rlams(2,2),w1(2,2),w2(2,2),b(20000),
     2      uu(20000),vv(20000),ww(20000),c(20000),
     2      b2(20000)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
         PRINT *, 'ENTER n'
         READ *,n
         CALL PRINF('n=*',n,1 )
c 
        a(1,1)=-3
        a(1,2)=7
        a(2,1)=-11
        a(2,2)=1-1.0d-11
ccc        call d2muv0(a,rlam,u,v)
        call d2moneuv(a,rlam,u,v)
         call prin2('u as computed*',u,4)
         call prin2('v as computed*',v,4)
c 
c       now, test the calculation
c 
        rlams(1,1)=rlam(1)
        rlams(2,1)=0
        rlams(1,2)=0
        rlams(2,2)=rlam(2)
c 
        call d2mprod2(u,rlams,w1)
        call d2mprod2(w1,v,w2)
         call prin2('w2=*',w2,4)
         call prin2('while a=*',a,4)
c 
        w1(1,1)=w2(1,1)-a(1,1)
        w1(2,1)=w2(2,1)-a(2,1)
        w1(1,2)=w2(1,2)-a(1,2)
        w1(2,2)=w2(2,2)-a(2,2)
         call prin2('and errors are*',w1,4)
c 
c       construct the big test matrix
c 
        call matrcre(b,n)
        call RIRAND(N**2,b)
        do 1200 i=1,n**2/2
        b(i)=0
 1200 continue
        call prin2('matrix as created*',b,n*n)
c 
c       . . . and construct its svd
c 
        eps=1.0d-16
        t1= second(tt)
        ifuv=1
        call d2mudv(b,n,uu,vv,eps,ifuv,numrows)
ccc        call d2mudv(b,n,uuu,vvv,eps,ifuv,numrows)
        t2= second(tt)
         call prin2('and CPU time is*',t2-t1,1)
         call prinf('and numrows=*',numrows,1)
        call prin2('singular values after d2mudv*',b,n)
cccc        call prin2('matrix after d2mudv*',b,n*n)
cccc        call prin2('uu after d2mudv*',uu,n*n)
cccc        call prin2('vv after d2mudv*',vv,n*n)
c 
c       now, test the svd
c 
        call matrcre(b2,n)
        call RIRANDin(i)
        call RIRAND(N**2,b2)
        do 2200 i=1,n**2/2
        b2(i)=0
 2200 continue
        call svdtes(b,uu,vv,ww,b2,c,n)
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine svdtes(s,u,v,w,b,w2,n)
        implicit real *8 (a-h,o-z)
        save
        dimension s(1),u(n,n),v(n,n),w(n,n),b(n,n),
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
cccc         call prin2('in svdtes, w=*',w,n*n)
c 
c       construct all proper products
c 
        call matmul(u,w,w2,n)
        call matmul(w2,v,w,n)
cccc          call prin2('and b recomputed from svd*',w,n*n)
cccc          call prin2('while b *',b,n*n)
c 
c        now, calculate the error of the decomposition
c 
        dd=0
        do 1800 i=1,n
        do 1600 j=1,n
        dd=dd+(b(j,i)-w(j,i))**2
 1600 continue
 1800 continue
        dd=dsqrt(dd)
         call prin2('and l2-error is*',dd,1)
        return
        end
c 
c 
c 
c 
c 
        SUBROUTINE RIRAND(N,Y)
        IMPLICIT REAL *8 (A-H,O-Z)
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
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),b(n,n),c(n,n)
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
        implicit real *8 (a-h,o-z)
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
        implicit real *8 (a-h,o-z)
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
        subroutine d2mudv(a,n,u,v,eps,ifuv,numrows)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),u0(2,2),v0(2,2),u(n,n),v(n,n),
     1      b(2,2),rlam(2)
c 
c        this subroutine uses the classical jacobi iteration
c        to construct the singular value decomposition of a
c        general real matrix, of the form
c 
c          a = u s v,                                        (1)
c 
c        with u,v orthogonal matrices, and s diagonal, with
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
c  a - the first column of a contains the (ordered) singular values
c  u,v - orthogonal matrices in (1)
c  numrows - the number of jacobi iterations (i.e. single rotations)
c        that the algorithm has taken
c 
c       . . . one sweep after another, conduct jacobi
c             iterations
c 
        maxsweep=10
        if(ifuv .eq. 0) goto 1150
        zero=0
        do 1100 i=1,n
        do 1050 j=1,n
        u(j,i)=zero
        v(j,i)=zero
 1050 continue
        u(i,i)=1
        v(i,i)=1
 1100 continue
 1150 continue
        numrows=0
        thresh=dsqrt(eps)
        thresh=dsqrt(thresh)
        thresh=dsqrt(thresh)
        do 3400 iijjkk=1,4
c 
        do 3000 ijk=1,maxsweep
        ifany=0
c 
        do 2000 i=1,n
        do 1800 j=1,n
         if(j .eq. i) goto 1800
         if( (dabs(a(i,j)) .lt. thresh) .and.
     1       (dabs(a(j,i)) .lt. thresh) ) goto 1800
        ifany=1
        numrows=numrows+1
c 
c       construct the svd of the (i,j)-th 2 * 2 submatrix
c 
        b(1,1)=a(i,i)
        b(1,2)=a(i,j)
        b(2,1)=a(j,i)
        b(2,2)=a(j,j)
ccc        call d2muv0(b,rlam,u0,v0)
        call d2moneuv(b,rlam,u0,v0)
c 
c       transform the big matrix
c 
c       . . . the rows
c 
        do 1200 l=1,n
        di=u0(1,1)*a(i,l)-u0(1,2)*a(j,l)
        dj=-u0(2,1)*a(i,l)+u0(2,2)*a(j,l)
        a(i,l)=di
        a(j,l)=dj
 1200 continue
c 
c       . . . and the columns
c 
        do 1400 l=1,n
        di=v0(1,1)*a(l,i)+v0(1,2)*a(l,j)
        dj=v0(2,1)*a(l,i)+v0(2,2)*a(l,j)
        a(l,i)=di
        a(l,j)=dj
 1400 continue
c 
c 
c       if the user so requested, adjust the matrices u,v
c 
        if(ifuv .eq. 0) goto 1800
c 
        do 1600 l=1,n
c 
c       the matrix u
c 
        di=u0(1,1)*u(l,i)-u0(1,2)*u(l,j)
        dj=-u0(2,1)*u(l,i)+u0(2,2)*u(l,j)
        u(l,i)=di
        u(l,j)=dj
c 
c       . . . and the columns
c 
        di=v0(1,1)*v(i,l)+v0(1,2)*v(j,l)
        dj=v0(2,1)*v(i,l)+v0(2,2)*v(j,l)
        v(i,l)=di
        v(j,l)=dj
 1600 continue
c 
 1800 continue
 2000 continue
c 
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
       call d2msvsrt(a(1,1),n,u,v,a(1,2),ifuv)
        return
        end
c 
c 
c 
c 
c 
        subroutine d2msvsrt(s,n,u,v,iw,ifuv)
        implicit real *8 (a-h,o-z)
        save
        dimension u(n,n),v(n,n),s(1),iw(1)
c 
c        if some of the alleged singular values
c        are negative - change their sign, and that
c        of the corresponding left eigenvectors
c 
        do 1100 i=1,n
        if(s(i) .ge. 0) goto 1100
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
        if(s(j). ge. s(j+1)) goto 1600
        jj=iw(j)
        iw(j)=iw(j+1)
        iw(j+1)=jj
c 
        d=s(j)
        s(j)=s(j+1)
        s(j+1)=d
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
        d=v(i,j)
        v(i,j)=v(ii,j)
        v(ii,j)=d
c 
c       exchange columns number i and ii in the matrix v
c 
        d=u(j,i)
        u(j,i)=u(j,ii)
        u(j,ii)=d
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
c 
c 
c 
c 
c 
        subroutine d2moneuv(a7,rlam,u,v)
        implicit real *8 (a-h,o-z)
        save
        dimension u(2,2),v(2,2),u0(2,2),v0(2,2),wl(2,2),
     1      uk(2,2),vk(2,2),u2(2,2),v2(2,2),a(2,2),a7(2,2),
     2      rlam(2),a2(2,2),wr(2,2)
        data eps/1.0d-15/
c 
c        initialize the matrices u,v
c 
        uk(1,1)=1
        uk(2,2)=1
        uk(1,2)=0
        uk(2,1)=0
c 
        call d2mcopy2(uk,vk)
        call d2mcopy2(a7,a)
c 
        do 2000 i=1,10
c 
c       attempt to construct the SVD of a
c 
        call d2muv0(ier,a,rlam,u0,v0,eps)
        if(ier .eq. 0) goto 2200
c 
c        the attempt failed. rotate the matrix a a little from both sides
c 
c        . . . construct the perturbing matrices
c 
        if(i .ne. 1) goto 1600
        d=1
        d=d/10
        alpha=dcos(d)
        beta=dsin(d)
c 
        wl(1,1)=alpha
        wl(2,2)=alpha
        wl(1,2)=beta
        wl(2,1)=-beta
c 
        d=1
        d=d/7
        alpha=dcos(d)
        beta=dsin(d)
c 
        wr(1,1)=alpha
        wr(2,2)=alpha
        wr(1,2)=beta
        wr(2,1)=-beta
 1600 continue
c 
c        . . . perturb the matrix from the left
c 
        call d2mprod2(wl,a,a2)
        wl(1,2)=-wl(1,2)
        wl(2,1)=-wl(2,1)
        call d2mprod2(uk,wl,u2)
        wl(1,2)=-wl(1,2)
        wl(2,1)=-wl(2,1)
c 
c        . . . perturb the matrix from the right
c 
        call d2mprod2(a2,wr,a)
        wr(1,2)=-wr(1,2)
        wr(2,1)=-wr(2,1)
        call d2mprod2(wr,vk,v2)
        wr(1,2)=-wr(1,2)
        wr(2,1)=-wr(2,1)
c 
        call d2mcopy2(u2,uk)
        call d2mcopy2(v2,vk)
 2000 continue
 2200 continue
        call d2mprod2(u0,uk,u)
        call d2mprod2(vk,v0,v)
  
        return
        end
c 
c 
c 
c 
        subroutine d2mcopy2(a,b)
        implicit real *8 (a-h,o-z)
        save
        dimension a(4),b(4)
c 
        b(1)=a(1)
        b(2)=a(2)
        b(3)=a(3)
        b(4)=a(4)
        return
        end
  
  
c 
c 
c 
c 
c 
        subroutine d2muv0(ier,a,rlam,u,v,eps2)
        implicit real *8 (a-h,o-z)
        save
        dimension a(2,2),u(2,2),v(2,2),rlam(2),
     1      b(2,2),w(2,2),vstar(2,2),z(2,2),rlams(2,2)
        data eps/1.0d-10/
c 
c       this subroutine produces the singular value decomposition
c       of a 2 * 2 real matrix a, so that on exit,
c 
c       a = u d v,
c 
c       with u and v orthogonal, and d a diagonal matrix with the
c       vector rlam=(rlam(1),rlam(2)) on the diagonal.
c 
        ier=0
c 
c       . . . simmetrize a
c 
        den=a(2,2)+a(1,1)
        rnum=a(2,1)-a(1,2)
        dd=dabs(a(1,1))+dabs(a(2,1))+dabs(a(1,2))+dabs(a(2,2))
        if(dabs(rnum) .gt. eps2*dd) goto 1100
        alpha=1
        beta=0
        goto 1400
 1100 continue
c 
c       if denominator is not too small, use exact formulae
c 
        if(dabs(den) .le. eps*dabs(rnum)) goto 1200
        tgphi=rnum/den
        phi=datan(tgphi)
        alpha=dcos(phi)
        beta=dsin(phi)
        goto 1400
 1200 continue
c 
c       denominator is too small. bomb.
c 
        ier=4
        return
 1400 continue
c 
c      calculate the simmetrizing matrix and the simmetric one
c 
        w(1,1)=alpha
        w(1,2)=beta
        w(2,1)=-beta
        w(2,2)=alpha
        call d2mprod2(w,a,b)
c 
c       now, diagonalize the simmetrized matrix
c 
        den=b(2,2)-b(1,1)
        rnum=b(1,2)*2
        if(dabs(rnum) .gt. eps2*dd) goto 1500
        alpha=1
        beta=0
        goto 1800
 1500 continue
c 
c       if denominator is not too small, use exact formulae
c 
        if(dabs(den) .le. eps*dabs(rnum)) goto 1600
        tg2phi=-rnum/den
        phi=datan(tg2phi)/2
        alpha=dcos(phi)
        beta=dsin(phi)
        goto 1800
 1600 continue
c 
c       denominator is too small. bomb.
c 
        ier=8
        return
 1800 continue
c 
c       construct the diagonalizing matrix
c       and the resulting diagonal
c 
        v(1,1)=alpha
        v(1,2)=beta
        v(2,1)=-beta
        v(2,2)=alpha
c 
        vstar(1,1)=v(1,1)
        vstar(1,2)=v(2,1)
        vstar(2,1)=v(1,2)
        vstar(2,2)=v(2,2)
c 
        call d2mprod2(v,w,u)
c 
c       finally, compute the resulting diagonal elements
c 
        call d2mprod2(u,a,z)
        call d2mprod2(z,vstar,rlams)
        rlam(1)=rlams(1,1)
        rlam(2)=rlams(2,2)
c 
        u(1,2)=-u(1,2)
        u(2,1)=-u(2,1)
        return
        end
c 
c 
c 
c 
c 
        subroutine d2mprod2(a,b,c)
        implicit real *8 (a-h,o-z)
        save
        dimension a(2,2),b(2,2),c(2,2)
c 
        c(1,1)=a(1,1)*b(1,1)+a(1,2)*b(2,1)
        c(1,2)=a(1,1)*b(1,2)+a(1,2)*b(2,2)
        c(2,1)=a(2,1)*b(1,1)+a(2,2)*b(2,1)
        c(2,2)=a(2,1)*b(1,2)+a(2,2)*b(2,2)
        return
        end
