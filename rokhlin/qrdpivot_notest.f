  
c 
c 
c 
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
c 
        subroutine qrdpivot(a,n,m,ncols,eps,w,s)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),w(n,1),s(1)
c 
c       this subroutine constructs the left singular values and the
c       singular vectors of the user-specified rectangular matrix
c       a(n,m). In the process, the matrix A is destroyed, The
c       singular vectors are returned in the first ncols columns
c       of the matrix a.
c 
c                            Input parameters:
c 
c  a - the matrix whose singular vectors and singular values are
c       to be determined
c  n,m - dimensionalities of the matrix a
c  eps - the accuracy to which the calculations are to be performed.
c       specifically, all singular values that are less than eps
c       will be replaced with zero
c 
c                            Output parameters:
c 
c  a - the first ncols columns of a will contain the first ncols
c       singular vectors of a
c  ncols - the number of singular values of a that are greater
c       than eps
c  s - the singular values of a
c 
c                            Work array:
c 
c  w - must be at least (n*m)*9/8+n+10 real *8 locations long
  
c 
c       conduct preliminary QR-decomposition on the columns of a
c 
        call qrpiv0(a,w,n,m,ncols,eps,w(1,m+1) )
c 
        call prinf('after qrpiv0, ncols=*',ncols,1)
  
c 
c       transpose the matrix a "in place"
c 
        call transmat(n,m,a,w(1,m+1))
c 
c        convert a to a triangular form by Jacobi rotations
c        from the left; the rotations are discarded
c 
        call qrdmat(a,m,ncols)
c 
c        reformat a to the form of a matrix ncols*ncols
c 
        call qrdref(a,m,ncols,a)
c 
c       finally, run a through the Jacobi SVD routine
c 
        call qrmudv(a,ncols,w,eps,numrows,n)
c 
        call prinf('after qrmudv, numrows=*',numrows,1)
c 
c       store the left sungular vectors of the matrix A
c       in the first columns of a, while returning the
c       sigular values in the array ww
c 
        do 2200 i=1,ncols
        s(i)=a(i,1)
 2200 continue
c 
        do 2600 i=1,ncols
        do 2400 j=1,n
        a(j,i)=w(j,i)
 2400 continue
 2600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine qrdref(a,m,ncols,b)
        implicit real *8 (a-h,o-z)
        save
        dimension a(m,ncols),b(ncols,ncols)
c 
        do 1400 i=1,ncols
        do 1200 j=1,ncols
c 
        b(j,i)=a(j,i)
 1200 continue
 1400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine qrdmat(a,n,m)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m)
c 
c        eliminate the lower "triangle" in the matrix a
c        by Jacobi rotations
c 
        do 2000 i=1,m
c 
        do 1800 j=i+1,n
c 
c       find the rotation eliminating the element a(j,i) by
c       combining rows i and j
c 
c       . . . if |a(j,i)| > |a(i,i)|, transpose the rows
c 
        if(abs(a(j,i)) .lt. eps) goto 1800
        if(abs(a(j,i)) .le. abs(a(i,i)) ) goto 1600
c 
        do 1400 k=i,m
c 
        d=a(i,k)
        a(i,k)=a(j,k)
        a(j,k)=d
 1400 continue
c 
 1600 continue
c 
        phi=atan(a(j,i)/a(i,i))
        u11=cos(phi)
        u12=sin(phi)
        u21=-u12
        u22=u11
c 
        do 1700 k=i,m
c 
        d1=u11*a(i,k)+u12*a(j,k)
        d2=u21*a(i,k)+u22*a(j,k)
        a(i,k)=d1
        a(j,k)=d2
 1700 continue
c 
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
        subroutine mulmata(a,n,k,b,m,c)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,k),b(m,k),c(n,m)
c 
        do 2000 i=1,n
        do 1800 j=1,m
c 
        d=0
        do 1600 l=1,k
        d=d+a(i,l)*b(j,l)
 1600 continue
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
        subroutine qrpiv0(a,b,n,m,ncols,eps,ww)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),b(n,m),ww(1)
c 
c        . . . copy the user-supplied matrix a into b
c 
        do 1400 i=1,m
        do 1200 j=1,n
        b(j,i)=a(j,i)
 1200 continue
 1400 continue
c 
c        apply the gram-schmidt proces (with pivoting) to b
c 
         ifpivot=1
         call qrdgrm(b,n,m,ww,eps,ncols,ifpivot)
c 
c        project the original matrix on the obtained orthogonal
c        columns
c 
        do 2400 j=1,m
        do 2200 i=1,ncols
c 
        call qrdscapr(a(1,j),b(1,i),n,prod)
        ww(i)=prod
 2200 continue
c 
        do 2300 i=1,ncols
        a(i,j)=ww(i)
 2300 continue
c 
 2400 continue
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
        subroutine qrmudv(a,n,u,eps,numrows,nn)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),u0(2,2),v0(2,2),u(nn,n),
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
        call qrdoneuv(b,rlam,u0,v0)
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
        do 1600 l=1,nn
c 
c       the matrix u
c 
        di=v0(1,1)*u(l,i)+v0(1,2)*u(l,j)
        dj=v0(2,1)*u(l,i)+v0(2,2)*u(l,j)
c 
        u(l,i)=di
        u(l,j)=dj
c 
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
       call qrmsvsrt(a(1,1),n,u,v,a(1,2),ifuv)
        return
        end
c 
c 
c 
c 
c 
        subroutine qrmsvsrt(s,n,u,v,iw,ifuv)
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
        subroutine transmat(n,m,a,index)
        implicit real *8 (a-h,o-z)
        save
        logical *1 index(1)
        dimension a(1)
c 
c        this subroutine transposes a real matrix a(n,m)
c        "in place" returning the real matrix a(m,n). it
c        uses a work array index that has to be at least
c        m*n logical *1 positions long (or, equivalently,
c        at least n*m/32 real *8 positions long
c 
c        . . . set the index array to .false.
c 
        ntran=0
        nn=n*m
        do 1200 i=1,nn
c 
        index(i)=.false.
 1200 continue
c 
c       one cycle after another, conduct the transposition
c 
        iijj=0
c 
        do 3000 iijjkk=1,nn
c 
        iijj=iijj+1
c 
c       if this element has already been transposed - skip it
c 
        if(index(iijj) .eq. .true.) goto 2800
c 
c       this element has not been transposed; start a cycle
c 
        ij=iijj
        dold=a(ij)
c 
        do 1600 k=1,nn
c 
        ntran=ntran+1
c 
        if(index(ij) .eq. .true.) goto 2800
c 
        call lintomat(n,ij,i,j)
c 
        call mattolin(m,j,i,ji)
c 
c       if the element b(j,i) has been transposed already,
c       terminate this cycle
c 
        if(ji .eq. ij) goto 2800
c 
c        we are putting a(i,j) into b(j,i) and marking a(ij) as
c        having been transposed
c 
        dnew=a(ji)
        a(ji)=dold
        dold=dnew
        index(ij)=.true.
c 
        ij=ji
 1600 continue
c 
 2800 continue
 3000 continue
  
        call prinf('and number of elementary transpositions is*',
     1      ntran,1)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine mattolin(n,i,j,ij)
c 
        save
        ij=(i-1)*n + j
        return
        end
c 
        subroutine lintomat(n,ij,i,j)
c 
        save
        i=(ij-1)/n+1
        j=ij-(i-1)*n
        return
        end
c 
c 
c 
c 
c 
        subroutine qrdgrm(b,n,m,rnorms,eps,ncols,ifpivot)
        implicit real *8 (a-h,o-z)
        save
        dimension b(n,m),rnorms(1)
c 
c       this subroutine applies a pivoted double gram-schmidt
c       procedure to the matrix b, resulting in a collection of
c       orthogonal vectors. the number of these vectors is
c       equal to the numerical rank of b, given the precision
c       of computations eps, which is supplied by the user.
c 
c                    input paramneters:
c 
c  b - the matrix to be gram-schmidt'ed. it is destroyed by this
c       subroutine
c  n,m - dimensionalities of the matrix b
c  eps - the machine accuracy
c 
c                     output parameters:
c 
c  b - the matrix of gram-schmidt vectors of the matrix a. note
c        that on exit from this subroutine only the first
c        ncols columns of  b  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit b is dimensioned
c        b(n,ncols)
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrix b on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c 
c       . . . prepare the array of values of norms
c             of columns
c 
        done=1
        dtot=0
        do 1400 i=1,m
        d=0
        do 1200 j=1,n
        d=d+b(j,i)**2
 1200 continue
        rnorms(i)=d
        dtot=dtot+d
 1400 continue
c 
c       . . . conduct gram-schmidt iterations
c 
        thresh=dtot*eps**2
        do 4000 i=1,m
c 
c       find the pivot
c 
         if(ifpivot .eq. 0) goto 2700
        ipivot=i
        rn=rnorms(i)
c 
        do 2200 j=i+1,m
        if(rnorms(j) .le. rn) goto 2200
        rn=rnorms(j)
        ipivot=j
 2200 continue
 2400 continue
c 
c       put the column number ipivot in the i-th place
c 
        do 2600 j=1,n
        d=b(j,i)
        b(j,i)=b(j,ipivot)
        b(j,ipivot)=d
 2600 continue
c 
        d=rnorms(ipivot)
        rnorms(ipivot)=rnorms(i)
        rnorms(i)=d
 2700 continue
c 
c       orthogonalize the i-th column to all preceeding ones
c 
        if(i .eq. 1) goto 2790
        do 2780 j=1,i-1
c 
        call qrdscapr(b(1,i),b(1,j),n,d)
c 
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*d
 2770 continue
 2780 continue
 2790 continue
c 
c       normalize the i-th column
c 
        call qrdscapr(b(1,i),b(1,i),n,d)
c 
c       if the norm of the remaining part of the matrix
c       is sufficiently small - terminate the process
c 
        if(d .lt. thresh) return
        ncols=i
c 
        d=done/dsqrt(d)
        do 2800 j=1,n
        b(j,i)=b(j,i)*d
 2800 continue
c 
c        orthogonalize everything else to it
c 
        do 3200 j=i+1,m
c 
        if(rnorms(j) .lt. thresh) goto 3200
c 
        call qrdscapr(b(1,i),b(1,j),n,d)
c 
         rrn=0
        do 3000 l=1,n
        b(l,j)=b(l,j)-b(l,i)*d
        rrn=rrn+b(l,j)**2
 3000 continue
        rnorms(j)=dsqrt(rrn)
 3200 continue
 3400 continue
 4000 continue
  
 4200 continue
  
  
        return
        end
c 
c 
c 
c 
c 
        subroutine qrdscapr(x,y,n,prod)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1)
c 
        prod=0
        do 1200 i=1,n
        prod=prod+x(i)*y(i)
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine qrdoneuv(a,rlam,u,v)
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
c       . . . simmetrize a
c 
        den=a(2,2)+a(1,1)
        rnum=a(2,1)-a(1,2)
c 
c       if denominator is not too small, use exact formulae
c 
        if(dabs(den) .lt. eps*dabs(rnum)) goto 1200
        tgphi=rnum/den
        phi=datan(tgphi)
        alpha=dcos(phi)
        beta=dsin(phi)
        goto 1400
 1200 continue
c 
c       denominator is too small. use taylor series
c 
        alpha=den/rnum
        beta=1
 1400 continue
c 
c      calculate the simmetrizing matrix and the simmetric one
c 
        w(1,1)=alpha
        w(1,2)=beta
        w(2,1)=-beta
        w(2,2)=alpha
        call qrdprod2(w,a,b)
c 
c       now, diagonalize the simmetrized matrix
c 
        den=b(2,2)-b(1,1)
        rnum=b(1,2)*2
c 
c       if denominator is not too small, use exact formulae
c 
        if(dabs(den) .lt. eps*dabs(rnum)) goto 1600
        tg2phi=-rnum/den
        phi=datan(tg2phi)/2
        alpha=dcos(phi)
        beta=dsin(phi)
        goto 1800
 1600 continue
c 
c       denominator is too small. use taylor expansion
c 
        alpha=dsqrt((1-den/rnum)/2)
        beta=dsqrt((1+den/rnum)/2)
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
        call qrdprod2(v,w,u)
c 
c       finally, compute the resulting diagonal elements
c 
        call qrdprod2(u,a,z)
        call qrdprod2(z,vstar,rlams)
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
        subroutine qrdprod2(a,b,c)
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
