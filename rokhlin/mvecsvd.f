        implicit real *8 (a-h,o-z)
        dimension w(1000 000),rlams(10 000),
     1      a2(1000 000)
c 
        external matfun2
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER n'
        READ *,n
        CALL PRINf('n=*',n,1 )
C 
        PRINT *, 'ENTER nvects'
        READ *,nvects
        CALL PRINf('nvects=*',nvects,1 )
c 
c       construct the compressed matrix
c 
        eps=1.0d-14
c 
c       construct the SVD
c 
        lenw=1 000 000
        lenrlams=10 000
        ifsvd=0
        call mvecsvd(ier,matfun2,par1,par2,par3,par4,par5,n,
     1    nvects,ifsvd,eps,ncols,rlams,lenrlams,
     2    w,lenw,lused)
c 
        call prin2('after mvecsvd, rlams=*',rlams,ncols)
        call prinf('after mvecsvd, lused=*',lused,1)
  
  
c 
c       now, compute the SVD explicitly
c 
        call explic(a2,n,nvects,eps)
  
c 
        call diffcomp(w,a2,n,ncols,rlams)
  
  
        stop
        end
c 
c 
c 
c 
c 
  
        subroutine matfun2(i,j,f,par1,par2,par3,par4,par5)
        implicit real *8 (a-h,o-z)
c 
        save
        done=1
c 
        di=i
        dj=j
c 
        f=done/(di**2+dj)**2 *(-1)**i +(-1)**j
  
        par5=1
  
        return
        end
  
c 
c 
c 
c 
c 
        subroutine diffcomp(a,b,n,ncols,rlams)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,ncols),b(n,ncols),diffs(10 000),
     1      rlams(1)
c 
c       one column after another, compute the difference
c       between the matrices a, b. If for some column it
c       appears large, try replacing it by the sum
c 
       do 1400 i=1,ncols
       d1=0
       d2=0
       do 1200 j=1,n
       d1=d1+(b(j,i)-a(j,i))**2
       d2=d2+(b(j,i)+a(j,i))**2
 1200 continue
c 
       diffs(i)=sqrt(d1)
       if(d2 .lt. d1) diffs(i)=sqrt(d2)
 1400 continue
c 
       call prin2('and in diffcomp, diffs=*',diffs,ncols)
       call prin2('and,remember, rlams=*',rlams,ncols)
       return
       end
c 
c 
c 
c 
c 
        subroutine explic(a,n,nvects,eps)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,nvects),w(300 000),s(10000)
c 
c       construct the matrix explicitly
c 
        do 1400 i=1,nvects
        do 1200 j=1,n
c 
        call matfun2(j,i,a(j,i),par1,par2,par3,par4,par5)
 1200 continue
 1400 continue
c 
        call qrdpivot(a,n,nvects,ncols,eps,w,s)
c 
        call prin2('and singular values computed explicitly*',
     1      s,ncols)
c 
  
cccc        call prin2('and a explicitly*',a,n*ncols)
        call prinf('in the explicit computation, ncols=*',ncols,1)
  
  
        return
        end
c 
c 
c 
c 
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c       This is the end of the testing code and the beginning of the
c       SVD code proper.
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c 
c 
        subroutine mvecsvd(ier,matfun,par1,par2,par3,par4,par5,n,
     1    nvects,ifsvd,eps,ncols,rlams,lenrlams,w,lenw,lused)
        implicit real *8 (a-h,o-z)
        save
        dimension rlams(1),
     1      par1(1),par2(1),par3(1),w(1),par4(1),par5(1)
c 
        external matfun
c 
c       this subroutine constructs the left singular vectors and the
c       singular values of the user-specified rectangular matrix.
c       The matrix is provided via the user-supplied subroutine
c       matfun (see below). The principal output of the subroutine
c       are the compressed version of the matrix (returned in the
c       array a, and the vector of its singular values. The principal
c       anticipated uses for this subroutine are in the Least squares
c       and quadrature environments.
c 
c                    Input parameters:
c 
c  matfun - the subroutine evaluating the entries of the matrix whose
c       singular values and left singular vectors are to be found.
c       The calling sequence of matfun is
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c                   matfun(j,i,aji,par1,par2,par3,par4,par5)
c 
c       Description of the calling sequence of matfun:
c 
c                        Input:
c 
c  j,i - the index of the element of the matrix to be created
c  par1, par2, par3, par4, par5 - parameters to be used by
c        matfun (real, integer, variable, array - whatever)
c 
c                        Output:
c 
c  aji - the value of the element of the matrix with indices (j,i)
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c  n - the number of rows of the matrix to be SVD'ed
c  nvects - the number of columns of the matrix to be SVD'ed
c  ifsvd - an integer parameter telling the subroutine whether it
c       should return to the user the true singular vectors, or
c       the QR vectors. Setting ifsvd = 1 will cause the subroutine
c       to return the SVD. Setting ifsvd=0 will cause the return of
c       QR vectors. The latter option is sometimes slightly less
c       expensive, and sometimes drastically less expensiove.
c  eps - the accuracy to which the calculations are to be performed.
c       specifically, all singular values that are less than eps
c       will be replaced with zero
c  lenrlams - the amount of space (in real *8 words) provided in the
c        output array rlams. It is used by the subroutine to
c        bomb if there is not enough space allocated. See the description
c        of the output parameter ier below for details. Also, please note
c        that these are small arrays compared to the arrays w, a; one
c        should not skimp on their size.
c  lenw - the amount of space provided to the subroutine in the work
c        array w; should be large.
c 
c                            Output parameters:
c 
c  w - on exit, the first ncols*n elements of w are formatted as
c        an array (n,ncols). Its ncols columns will contain the first
c        ncols singular vectors of the input matrix
c  ncols - the number of singular values of a that are greater
c       than eps
c  rlams - the singular values of a (ncols of them)
c 
c                            Work array:
c 
c  w - must be big
c 
c        subdivide the nvects input vectors into blocks
c 
        ier=0
c 
        nopt=100
        if(nopt .ge. nvects/2) nopt=nvects/2-1
c 
        nblocks=nvects/nopt
        ileft=nvects-nblocks*nopt
c 
        if(ileft .lt. nopt/2) goto 1200
c 
c       the leftover chunk is large - make it into a separate matrix
c 
        nblocks=nblocks+1
        nlast=ileft
        goto 1400
c 
 1200 continue
c 
c       the leftover chunk is small - add it to the last
c       matrix already to be created
c 
        nlast=nopt+ileft
c 
 1400 continue
c 
        maxblock=nopt
        if(nlast .gt. maxblock) maxblock=nlast
c 
         call mveccini(nopt,nblocks,nlast)
c 
c        compress the user-supplied set of vectors
c 
        irhs=1
        lrhs=lenrlams+100
c 
        iw=irhs+lrhs
        lenw2=lenw-lrhs-2
c 
        call mveccomp(jer,matfun,par1,par2,par3,par4,par5,n,
     1    nblocks,maxblock,eps,rlams,w(irhs),lenrlams,
     2      ncols,w(iw),lenw2,lused)
c 
        call prinf('in mvecsvd after mveccomp, ncols=*',ncols,1)
c 
        lused=lused+iw
c 
c       if the user so requested, return to him the QR-vectors,
c       as opposed to singular vectors
c 
        if(jer .ne. 0) then
            ier=16
            call prinf('bombing after mveccomp, with jer=*',jer,1)
            return
        endif
c 
        call mvecmove(w(iw),w,ncols*n+2)
c 
        do 1500 i=1,ncols
c 
        rlams(i)=sqrt(rlams(i))
 1500 continue
c 
        if(ifsvd .eq. 0) return
c 
c 
c        project all of the user-specified vectors on the obtained
c        basis, and use the Jacobi rotations to construct the QR
c        decomposition of the projected matrix
c 
        ia=1
        la=ncols*n+2
c 
        ib=1+la
        lb=ncols**2+10
c 
        icol=ib+lb
        lcol=n+2
c 
        icol2=icol+lcol
        lcol2=ncols+2
c 
        lused2=icol2+lcol2
        if(lused2 .gt. lused) lused=lused2
c 
        call smallmat(w(ia),n,ncols,nvects,w(ib),par1,par2,par3,
     1      par4,par5,w(icol),w(icol2),matfun,eps)
  
  
cccc        call prin2('after smallmat, b=*',w(ib),ncols**2)
c 
c       construct the SVD of the matrix b
c 
        call mvecudv(w(ib),ncols,w(ia),eps,numrows,n)
  
        call prinf('after mvecudv, numrows=*',numrows,1)
c 
c       copy the singular values into the output array
c 
        do 2200 i=1,ncols
c 
        rlams(i)=w(ib+i-1)
 2200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine smallmat(phis,n,ncols,m,b,par1,par2,par3,
     1      par4,par5,col,col2,matfun,eps)
        implicit real *8 (a-h,o-z)
        save
        dimension phis(n,1),b(ncols,ncols),par1(1),par2(1),
     1      par3(1),par4(1),par5(1)
        dimension col(1),col2(1),aaa(3),uuu(2,2)
c 
        external matfun
c 
c        one column after another, construct the columns
c        of the original matrix, expand them into the
c        user-supplied orthogonal basis; use Jacobi rotations
c        to compress the whole mess of them
c 
c       . . . process the first ncols columns
c 
        do 2000 i=1,ncols
c 
        do 1400 j=1,n
c 
        call matfun(j,i,col(j),par1,par2,par3,par4,par5)
 1400 continue
c 
        do 1800 j=1,ncols
c 
        call mvecscap(col,phis(1,j),n,prod)
c 
        b(j,i)=prod
 1800 continue
c 
 2000 continue
c 
c        eliminate the lower "triangle" in the matrix b
c        by Jacobi rotations
c 
        do 3000 i=1,ncols-1
c 
        do 2800 j=i+1,ncols
c 
c       find the rotation eliminating the element b(j,i) by
c       combining rows i and j
c 
        aaa(1)=b(i,i)
        aaa(2)=b(i,j)
	call mverotfn(aaa,uuu)
c 
        do 2700 k=i,ncols
c 
        d1=uuu(1,1)*b(k,i)+uuu(1,2)*b(k,j)
        d2=uuu(2,1)*b(k,i)+uuu(2,2)*b(k,j)
        b(k,i)=d1
        b(k,j)=d2
 2700 continue
c 
 2800 continue
 3000 continue
c 
c        process columns with numbers ncols+1,...m
c 
        do 5000 k=ncols+1,m
c 
c       create column number k
c 
        do 3400 j=1,n
c 
        call matfun(j,k,col(j),par1,par2,par3,par4,par5)
 3400 continue
c 
c       expand column number k into the ncols basis functions
c 
        do 3800 j=1,ncols
c 
        call mvecscap(col,phis(1,j),n,prod)
c 
        col2(j)=prod
 3800 continue
c 
c        use Jacobi rotations to eliminate the newly
c        obtained column
c 
        do 4900 i=1,ncols
c 
c       . . . find the rotation
c 
        aaa(1)=b(i,i)
        aaa(2)=col2(i)
	call mverotfn(aaa,uuu)
c 
        do 4700 j=i,ncols
c 
        d1=uuu(1,1)*b(j,i)+uuu(1,2)*col2(j)
        d2=uuu(2,1)*b(j,i)+uuu(2,2)*col2(j)
        b(j,i)=d1
        col2(j)=d2
 4700 continue
c 
 4900 continue
c 
 5000 continue
c 
c       transpose b in place
c 
        do 5400 i=1,ncols
        do 5200 j=1,i
c 
        d=b(i,j)
        b(i,j)=b(j,i)
        b(j,i)=d
 5200 continue
 5400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine mveccrea(nummat,a,rhs,n,m,matfun,
     1      par1,par2,par3,par4,par5)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,1),rhs(1),par2(1),par3(1),par4(1),par5(1)
c 
c       construct the submatrix number nummat of the big matrix
c 
        imatr=nummat
        istart=(imatr-1)*nopt
c 
        m=nopt
        if(imatr .eq. nblocks) m=nlast
c 
        do 2000 i=1,m
        do 1800 j=1,n
c 
        call matfun(j,i+istart,d,par1,par2,par3,par4,par5)
c 
        a(j,i)=d
 1800 continue
c 
 2000 continue
c 
cccc        call prin2('exiting mveccrea, a=*',a,n*m)
  
        do 2200 i=1,m
        rhs(i)=0
 2200 continue
c 
        return
c 
c 
c 
c 
        entry mveccini(nopt7,nblocks7,nlast7)
c 
        nopt=nopt7
        nblocks=nblocks7
        nlast=nlast7
        return
        end
c 
c 
c 
c 
c 
        subroutine mveccomp(ier,matfun,par1,par2,par3,par4,
     1      par5,n,nblocks,maxblock,eps,rnorms,rhs,lenrhs,
     2      ncols,w,lenw,lwused)
        implicit real *8 (a-h,o-z)
        save
        dimension rnorms(1),rhs(1),w(1),par1(1),par2(1),
     1      par3(1),par4(1),par5(1)
        external matfun
c 
c       this subroutine compresses the user-provided matrix via a
c       nested pivoted gram-schmidt process, applying the same
c       linear mapping to the user provided vector. Both the matrix
c       and the vector are produced by the user-provided subroutine
c       matfun. The principal output of the subroutine are the
c       compressed version of the matrix (returned in the array
c       a), and the vector rhs after the gram-schmidt operator has
c       been applied to it. The principal anticipated uses for this
c       subroutine are in the Least squares and quadrature environments.
c       Please note that this is almost entirely a memory management
c       subroutine. Virtually all work is performed by the subroutine
c       mvecrecu (see).
c 
c                    Input parameters:
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c  matfun - the subroutine creating blocks of the matrix to be
c       compressed, and the corresponding parts of the right-hand side.
c       The calling sequence of matfun is
c 
c        matfun(nummat,a,rhs,n,m,par1,par2,par3,par4,par5)
c 
c        The input parameters for matfun:
c 
c    nummat - the sequence number of the submatrix to be generated
c    n - the number of rows in each submatrix (and also in the whole
c        matrix)
c    par1, par2, par3, par4, par5  - five parameters to be used by
c        matfun (real, integer, whatever)
c 
c        The output parameters for matfun:
c 
c    a - the n \times m submatrix of the big mtrix
c    rhs - the part of the right-hand side corresponding to the part
c        of the matrix returned in a
c    m - the number of columns in a
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c  par1,par2,par3,par4,par5 - five parameters to be used by matfun. They
c         can be real, integer, variable, array, whatever.
c  n - the number of points at which each of the functions is tabulated
c  nblocks - the number of blocks to be created by calls to the
c        subroutine matfun.
c  maxblock - the maximum possible number of columns in any of the
c        submatrices to be created by the subroutine matfun
c  eps - the relative precision to which the compression will be
c        performed
c  boxes - the box structure corresponding to the total number of
c        blocks to be created. It has been (hopefully) constructed
c        by a prior call to the subroutine mvecboxs (see)
c  lenrhs - the amount of space (in real *8 words) provided in the
c        output arrays rnorms, rhs. it is used by the subroutine to
c        bomb if there is not enough space allocated. See the description
c        of the output parameter ier below for details. Also, please note
c        that these are small arrays compared to the arrays w, a; one
c        should not skimp on their size.
c  lenw - the amount of space provided to the subroutine in the work
c        array w; should be large.
c 
c                         Output parameters:
c 
c  ier - error return code:
c        ier=0 means successful conclusion
c        ier=8 means that the amount of space provided in array w is
c                insufficient
c        ier=16 means that the amount of space provided in array a is
c                insufficient
c        ier=32 means that the amount of space provided in arrays rhs,
c               rnorms is insufficient
c 
c  w - the first n*ncols elements of w are formatted as an n \times ncols
c        matrix containing the compressed version of the matrix generated
c        in pieces by the subroutine matfun (see description above). The
c        proper way to view it as a collection of ncols vectors.
c  rnorms - the array of norms of pivot columns procduced in the
c        gram-schmidt process; ncols of them
c  rhs - the compressed version of the right-hand side constructed
c        in pieces by matfun together with the matrix compressed in a
c  ncols - the number of columns in a and of elements in rhs, rnorms
c  lwused - the maximum number of elements in the work array w used
c        by the subroutine at any time.
c 
c       . . . allocate space for boxes and build the box structure
c 
        iboxes=1
        ninire=2
        lboxes=nblocks*20/ninire + 20
c 
        call mvecboxs(nblocks,w(iboxes),nboxes)
c 
c       compress the user-supplied matrix
c 
        iw=iboxes+lboxes+1
        lenw2=lenw-iw-1
c 
        call mvecrecu(ier,matfun,par1,par2,par3,par4,par5,
     1      maxblock,eps,w(iboxes),n,rnorms,rhs,lenrhs,n
     2      cols,w(iw),lenw2,lwused)
c 
        call mvecmove(w(iw),w,ncols*n+2)
c 
        lwused=lwused+lboxes+1
        return
        end
c 
c 
c 
c 
c 
        subroutine mvecrecu(ier,matfun,par1,par2,par3,
     1    par4,par5,mmax,eps,boxes,n,rnorms,rhs,lenrhs,
     2      ncols,w,lenw,lwused)
        implicit real *8 (a-h,o-z)
        save
        dimension rnorms(1),rhs(1),w(1)
        integer *4 boxes(10,1),stack(100),par1(1),par2(1),
     1      par3(1),par4(1),par5(1)
c 
        external matfun
c 
c       this subroutine compresses the user-provided matrix via a
c       nested pivoted gram-schmidt process, applying the same
c       linear mapping to the user provided vector. Both the matrix
c       and the vector are produced by the user-provided subroutine
c       matfun. The principal output of the subroutine are the
c       compressed version of the matrix (returned in the array
c       a), and the vector rhs after the gram-schmidt operator has
c       been applied to it. The principal anticipated uses for this
c       subroutine are in the Least squares and quadrature environments.
c 
c                    Input parameters:
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c  matfun - the subroutine creating blocks of the matrix to be
c       compressed, and the corresponding parts of the right-hand side.
c       The calling sequence of matfun is
c 
c        matfun(nummat,a,rhs,n,m,par1,par2,par3,par4,par5)
c 
c        The input parameters for matfun:
c 
c    nummat - the sequence number of the submatrix to be generated
c    n - the number of rows in each submatrix (and also in the whole
c        matrix)
c    par1, par2, par3, par4, par5  - five parameters to be used by
c        matfun (real integer, whatever)
c 
c        The output parameters for matfun:
c 
c    a - the n \times m submatrix of the big mtrix
c    rhs - the part of the right-hand side corresponding to the part
c        of the matrix returned in a
c    m - the number of columns in a
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c  par1,par2,par3,par4,par5 - five parameters to be used by matfun. They
c         can be real, integer, variable, array, whatever.
c  mmax - the maximum possible number of columns in any of the
c        submatrix
c  eps - the relative precision to which the compression will be
c        performed
c  boxes - the box structure corresponding to the total number of
c        blocks to be created. It has been (hopefully) constructed
c        by a prior call to the subroutine mvecboxs (see)
c  n - the number of points at which each of the functions is tabulated
c  lenrhs - the amount of space (in real *8 words) provided in the
c        output arrays rnorms, rhs. it is used by the subroutine to
c        bomb if there is not enough space allocated. See the description
c        of the output parameter ier below for details. Also, please note
c        that these are small arrays compared to the arrays w, a; one
c        should not skimp on their size.
c  lenw - the amount of space provided to the subroutine in the work
c        array w; should be large.
c 
c                         Output parameters:
c 
c  ier - error return code:
c        ier=0 means successful conclusion
c        ier=8 means that the amount of space provided in array w is
c                insufficient
c        ier=32 means that the amount of space provided in arrays rhs,
c               rnorms is insufficient
c 
c  w - the first n*ncols elements of w are formatted as an n \times ncols
c        matrix containing the compressed version of the matrix generated
c        in pieces by the subroutine matfun (see description above).
c        The proper way to view it as a collection of ncols vectors.
c  rnorms - the array of norms of pivot columns procduced in the
c        gram-schmidt process; ncols of them
c  rhs - the compressed version of the right-hand side constructed
c        in pieces by matfun together with the matrix compressed in a
c  ncols - the number of columns in a and of elements in rhs, rnorms
c  lwused - the maximum number of elements in the work array w used
c        by the subroutine at any time.
c 
c       . . . initialize the recursion
c 
        lena=1000 0000
        ier=0
        lwused=0
        istack=1
        stack(1)=1
        ibox=1
        call mvecinit(ier)
c 
c       conduct the recursion
c 
        do 4000 ijk=1,1 000 000
c 
        ibox=stack(istack)
c 
c        if this is a leaf box - build and compress the corresponding
c        matrix; store the result in array w
c 
        if(boxes(6,ibox) .gt. 0) goto 2200
c 
        nummat=boxes(3,ibox)
c 
        call mvecinfo(inuse)
c 
        lbuff=mmax*2+mmax*n+100
c 
        irhs=inuse +2  +lbuff
        lrhs=mmax+2
c 
        if(lenrhs .ge. lrhs) goto 1200
        ier=32
        return
 1200 continue
c 
        irnorms=irhs+lrhs
        lrnorms=mmax+2
c 
        ia=irnorms+lrnorms
        la=mmax*n +10
c 
        ltot=ia+la
        if(ltot .gt. lwused) lwused=ltot
        if(ltot .lt. lenw) goto 1400
        ier=8
        return
 1400 continue
c 
        call mvecreac(matfun,par1,par2,par3,par4,par5,nummat,
     1    eps,w(ia),w(irhs),n,m,w(irnorms),ncols)
c 
        index=boxes(1,ibox)
        call mvecstor(ier,index,w(ia),w(irhs),n,ncols,w(irnorms),
     1      w,lenw)
c 
        if(ier .ne. 0) return
c 
        boxes(9,ibox)=1
        istack=istack-1
        goto 4000
c 
 2200 continue
c 
c       this is not a leaf box. if the compressed matrix has not been
c       constructed for its first son - construct and store it
c 
        ison1=boxes(6,ibox)
        ifdone=boxes(9,ison1)
        if(ifdone .gt. 0) goto 2400
c 
        istack=istack+1
        stack(istack)=ison1
        goto 4000
c 
 2400 continue
c 
c       if the compressed matrix has not been constructed
c       for its second son - construct and store it
c 
        ison2=boxes(7,ibox)
        ifdone=boxes(9,ison2)
        if(ifdone .gt. 0) goto 2600
c 
        istack=istack+1
        stack(istack)=ison2
        goto 4000
c 
 2600 continue
c 
c       the compressed matrices have been constructed for both sonnies
c       of this box. retrieve and merge them. store the result on disk
c 
        ison1=boxes(6,ibox)
        ifdone1=boxes(9,ison1)
c 
        ison2=boxes(7,ibox)
        ifdone2=boxes(9,ison2)
c 
        if( (ifdone1. gt. 0) .and. (ifdone2. gt. 0) ) goto 2800
c 
        call prinf('disaster happened at istack=*',istack,1)
        stop
 2800 continue
c 
        index1=boxes(1,ison1)
        call mvecpnt(ier,index1,ia1,irhs1,n1,ncols1,irnorms1,w,inuse)
c 
        index2=boxes(1,ison2)
        call mvecpnt(ier,index2,ia2,irhs2,n,ncols2,irnorms2,w,inuse)
c 
        call prinf('before mvecmrg, inuse=*',inuse,1)
c 
        lbuff=mmax*2+mmax*n+100
c 
        irhs=inuse +2 +lbuff
        lrhs=ncols1+ncols2+2
c 
        irnorms=irhs+lrhs
        lrnorms=ncols1+ncols2+2
c 
        ia=irnorms+lrnorms
        la=(ncols1+ncols2)*n+10
c 
        ltot=ia+la
        if(ltot .gt. lwused) lwused=ltot
        if(ltot .le. lenw) goto 3000
        ier=8
        return
 3000 continue
c 
        call mvecmrg(w(ia1),ncols1,w(irnorms1),w(irhs1),
     1      w(ia2),ncols2,w(irnorms2),w(irhs2),n,eps,w(ia),
     2      ncols,w(irnorms),w(irhs) )
c 
        if(ncols+2 .le. lenrhs) goto 3400
        ier=32
        return
 3400 continue
c 
        if(istack .ne. 1) goto 3600
c 
        call mvecmove(w(irnorms),rnorms,ncols+1)
        call mvecmove(w(irhs),rhs,ncols+1)
        call mvecmove(w(ia),w,ncols*n+2)
c 
        return
c 
 3600 continue
c 
        call mvecrem(ier)
        call mvecrem(ier)
c 
        index=boxes(1,ibox)
        call mvecstor(ier,index,w(ia),w(irhs),n,ncols,w(irnorms),
     1      w,lenw)
c 
        if(ier .ne. 0) return
c 
        boxes(9,ibox)=1
        istack=istack-1
        goto 4000
c 
 4000 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine mvecstor(ier,index,a,rhs,n,ncols,rnorms,w,lenw)
        implicit real *8 (a-h,o-z)
        save
        dimension a(1),rhs(1),rnorms(1),w(1)
        integer *4 map(5,100)
c 
c        this entry stores in tha array w the three arrays a,rhs,
c        rnorms of lengths ncols*n, ncols, ncols, respectively.
c        it also records the parameters n, ncols in the internal
c        array map, indexed by the key index.
c 
c        description of entries in array map
c 
c       map(1,i) - the key indexing the data stored in the i-th
c                  location
c       map(2,i) - the location in array w where the information
c                  corresponding to the i-th chunk begins
c       map(3,i) - the length in the array w of the area allocated
c                  to the i-th chunk
c       map(4,i) - equal to the parameter n for this chunk of data
c       map(5,i) - equal to the parameter ncols for this chunk of data
c 
c        find where in the array w the data are to be stored
c 
        ier=0
        ii=map(2,nstored)+map(3,nstored)+1
c 
        irnorms=ii+1
        lrnorms=ncols+1
c 
        irhs=irnorms+lrnorms
        lrhs=ncols+1
c 
        ia=irhs+lrhs
        la=n*ncols+2
c 
c       if the storage array w is too short - bomb
c 
        ltot=ia+la
        if(ltot .le. lenw) goto 1400
        ier=8
        return
 1400 continue
c 
c        store the junk in array w
c 
        call mvecmove(rnorms,w(irnorms),lrnorms)
        call mvecmove(rhs,w(irhs),lrhs)
        call mvecmove(a,w(ia),la)
c 
c       update the array map
c 
        nstored=nstored+1
        map(1,nstored)=index
        map(2,nstored)=ii
        map(3,nstored)=ia+la-ii+1
        map(4,nstored)=n
        map(5,nstored)=ncols
c 
        return
c 
c 
c 
c 
        entry mvecretr(ier,index,a,rhs,n,ncols,rnorms,w)
c 
c        this entry retrieves from the array w the three arrays
c        a,rhs, rnorms of lengths ncols*n, ncols, ncols, respectively;
c        it finds them by the key index. it also retrieves the
c        parameters n, ncols from the internal array map
c 
c        . . . find in the array map the entry with the key index
c 
        if(nstored .gt. 1) goto 2100
        ier=8
        return
 2100 continue
c 
        do 2200 i=nstored,1,-1
c 
        ier=0
        if(index .ne. map(1,i) ) goto 2200
        k=i
        goto 2400
 2200 continue
        ier=4
        return
 2400 continue
c 
        ii=map(2,k)
c 
c       unpack the beginning of the part of array w where the data
c       for this chunk are stored
c 
        n=map(4,k)
        ncols=map(5,k)
        ii=map(2,k)
c 
        irnorms=ii+1
        lrnorms=ncols+1
c 
        irhs=irnorms+lrnorms
        lrhs=ncols+1
c 
        ia=irhs+lrhs
        la=n*ncols+2
c 
c        retrieve the junk from array w
c 
        call mvecmove(w(irnorms),rnorms,lrnorms)
        call mvecmove(w(irhs),rhs,lrhs)
        call mvecmove(w(ia),a,la)
c 
        return
c 
c 
c 
c 
        entry mvecinfo(inuse)
c 
        if(nstored .le. 1) inuse=2
        if(nstored .gt. 1) inuse=map(2,nstored)+map(3,nstored)
  
cccc        call prinf('in mvecinfo,map=*',map,nstored*10)
  
        return
c 
c 
c 
c 
        entry mvecpnt(ier,index,i7a,i7rhs,n,ncols,i7rnorms,w,inuse)
c 
c        this entry returns to the user the addresses in the array w
c        of the three arrays a,rhs, rnorms, (hopefully) stored there
c        previously by a call to the entry mvecstor (see) of this
c        subroutine. it finds them by the key index. it also retrieves
c        the parameters n, ncols from the internal array map
c 
c        . . . find in the array map the entry with the key index
c 
        if(nstored .gt. 1)
     1      inuse=map(2,nstored)+map(3,nstored)
c 
        if(nstored .le. 1) inuse=2
c 
        if(nstored .gt. 1) goto 3100
        ier=8
        return
 3100 continue
c 
        do 3200 i=nstored,1,-1
c 
        ier=0
        if(index .ne. map(1,i) ) goto 3200
        k=i
        goto 3400
 3200 continue
        ier=4
        return
 3400 continue
c 
        ii=map(2,k)
c 
c       unpack the beginning of the part of array w where the data
c       for this chunk are stored
c 
        n=map(4,k)
        ncols=map(5,k)
        ii=map(2,k)
c 
        irnorms=ii+1
        lrnorms=ncols+1
c 
        irhs=irnorms+lrnorms
        lrhs=ncols+1
c 
        ia=irhs+lrhs
        la=n*ncols+2
c 
        i7a=ia
        i7rnorms=irnorms
        i7rhs=irhs
c 
        return
c 
c 
c 
c 
        entry mvecinit(ier)
        nstored=1
        map(1,1)=-1
        map(2,1)=1
        map(3,1)=0
        return
c 
c 
c 
c 
        entry mvecrem(ier)
        nstored=nstored-1
        return
        end
c 
c 
c 
c 
c 
        subroutine mvecmove(x,y,n)
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
        subroutine mvecreac(matfun,par1,par2,par3,par4,par5,
     1    nummat,eps,a,rhs,n,m,rnorms,ncols)
        implicit real *8 (a-h,o-z)
        save
        dimension a(1),rhs(1),rnorms(1),par1(1),par2(1),par3(1),
     1    par4(1),par5(1)
c 
        external matfun
c 
c        construct the matrix number nummat
c 
        call mveccrea(nummat,a,rhs,n,m,matfun,par1,par2,par3,
     1      par4,par5)
c 
c       compress the matrix via gram-schmidt procedure,
c       while applying the same transformation to the
c       right-hand side
c 
        call mvecpiv(a,n,m,rhs,rnorms,eps,ncols)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine mvecboxs(n,boxes,nboxes)
        implicit real *8 (a-h,o-z)
        save
        integer *4 boxes(10,1),ns(100),laddr(100)
c 
c        this subroutine constructs the box structure
c        corresponding to the user-supplied number n.
c        For each box, the meaning of the entries is as
c        follows:
c 
c     box(1) - the sequence number of the box
c     box(2) - the level of subdivision of this box
c     box(3) - the first point living in this box
c     box(4) - the last point living in this box
c     box(5) - the address in the array boxes of this box's dad
c     box(6) - the address in the array boxes of this box's first son
c     box(7) - the address in the array boxes of this box's second son
c     box(8) - the number of points inside this box
c     box(9) - the length in the master array of the data corresponding
c        to this box
c     box(10) - reserved for future expansion
c 
c       . . . create the box on the level of subdivision 0
c 
        boxes(1,1)=1
        boxes(2,1)=0
        boxes(3,1)=1
        boxes(4,1)=n
        boxes(5,1)=-1
        boxes(6,1)=0
        boxes(7,1)=0
        boxes(8,1)=n
c 
        ns(1)=1
        icurr=2
        laddr(1)=1
c 
c       one level after another, subdivide all boxes
c 
        do 2000 level=0,100
c 
cccc        call prinf('level=*',level,1)
cccc        call prinf('number of boxes on this level=*',
cccc     1      ns(level+1),1)
  
cccc        call prinf('laddr(level+1)=*',laddr(level+1),1)
  
c 
        nblev=0
        do 1800 i=1,ns(level+1)
  
cccc        call prinf('i=*',i,1)
  
  
c 
c       subdivide box number i on level level
c 
        ibb=laddr(level+1)+i-1
        call mvecbdiv(boxes(1,ibb),boxes(1,icurr),boxes(1,icurr+1),
     1      ifdiv,icurr)
c 
c       if this box has been subdivided - act accordingly
c 
        if(ifdiv .eq. 0) goto 1800
c 
        nblev=nblev+2
 1800 continue
c 
        if(nblev .eq. 0) goto 2200
c 
        ns(level+2)=nblev
        laddr(level+2)=laddr(level+1)+ns(level+1)
 2000 continue
c 
 2200 continue
c 
        nboxes=icurr-1
        return
        end
c 
c 
c 
c 
c 
        subroutine mvecbdiv(dad,son1,son2,ifdiv,icurr)
        implicit real *8 (a-h,o-z)
        save
        integer *4 dad(10),son1(10),son2(10)
c 
c       check if this box needs to be subdivided
c 
cccc        call prinf('in mvecbdiv, dad=*',dad,10)
  
        n=dad(4)-dad(3)+1
        if(n .ge. 2) goto 1200
        ifdiv=0
        return
 1200 continue
c 
c       this box needs to be subdivided. do so
c 
        do 1400 i=1,10
        son1(i)=0
        son2(i)=0
 1400 continue
c 
        son1(1)=icurr
        son2(1)=icurr+1
c 
        son1(2)=dad(2)+1
        son2(2)=dad(2)+1
c 
        son1(3)=dad(3)
c 
        d=dad(4)+dad(3)
        d=d/2+1.0d-8
        son1(4)=d
c 
        son2(3)=son1(4)+1
        son2(4)=dad(4)
c 
        son1(5)=dad(1)
        son2(5)=dad(1)
c 
        dad(6)=son1(1)
        dad(7)=son2(1)
c 
        son1(8)=son1(4)-son1(3)+1
        son2(8)=son2(4)-son2(3)+1
c 
        icurr=icurr+2
        ifdiv=1
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine mvecmrg(a,ncolsa,rnormsa,rhsa,
     1      b,ncolsb,rnormsb,rhsb,n,eps,c,ncolsc,
     2    rnormsc,rhsc)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,ncolsa),b(n,ncolsb),rnormsa(ncolsa),
     1      rnormsb(ncolsb),rnormsc(1),rhsa(ncolsa),rhsb(ncolsb),
     2      rhsc(1),c(n,1)
  
c 
c       pack the matrices a,b into the larger matrix c
c 
        do 1400 i=1,ncolsa
c 
        d=sqrt(rnormsa(i))
cccc        d=rnormsa(i)
c 
        do 1200 j=1,n
           c(j,i)=a(j,i)*d
 1200 continue
c 
        rhsc(i)=rhsa(i)*d
c 
 1400 continue
c 
        do 1800 i=1,ncolsb
c 
        d=sqrt(rnormsb(i))
cccc        d=rnormsb(i)
c 
        do 1600 j=1,n
        c(j,i+ncolsa)=b(j,i)*d
 1600 continue
c 
        rhsc(i+ncolsa)=rhsb(i)*d
c 
 1800 continue
c 
c         feed the resulting matrix c into the Gram-schmidt routine
c 
        ncc=ncolsa+ncolsb
  
cccc        call prinf('in mvecmrg before mvecpiv, ncolsa=*',ncolsa,1)
cccc        call prinf('in mvecmrg before mvecpiv, ncolsb=*',ncolsb,1)
cccc        call prinf('in mvecmrg before mvecpiv, ncc=*',ncc,1)
  
cccc        call prin2('c=*',c,n*ncc)
  
cccc        call prin2('rnormsa=*',rnormsa,ncolsa)
cccc        call prin2('rnormsb=*',rnormsb,ncolsb)
  
  
  
        call mvecpiv(c,n,ncc,rhsc,rnormsc,eps,ncolsc)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine amatv(a,n,m,x,y)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),x(n),y(m)
c 
        do 1400 i=1,n
        d=0
        do 1200 j=1,m
        d=d+a(i,j)*x(j)
 1200 continue
c 
        y(i)=d
 1400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine mvecpiv(b,n,m,rhs,rnorms,eps,ncols)
        implicit real *8 (a-h,o-z)
        save
        dimension b(n,m),rnorms(1),rhs(m)
c 
c       this subroutine applies a pivoted double gram-schmidt
c       procedure to the matrix b, resulting in a collection of
c       orthogonal vectors. the number of these vectors is
c       equal to the numerical rank of b, given the precision
c       of computations eps, which is supplied by the user.
c       In addition, this subroutine applies to the user-supplied
c       vector RHS the same set of transformations that it applies
c       to the matrix b in order to make the latter orthogonal.
c       The applications of this subroutine in the least squares
c       environment are obvious.
c 
c 
c                    input parameters:
c 
c  b - the matrix to be gram-schmidt'ed. it is destroyed by this
c       subroutine
c  n,m - dimensionalities of the matrix b
c  rhs - the vector (user-supplied) to which the subroutine will
c       apply the same transformation as to the matrix b
c  eps - the machine accuracy
c 
c                     output parameters:
c 
c  b - the matrix of gram-schmidt vectors of the matrix a. note
c        that on exit from this subroutine only the first
c        ncols columns of  b  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit b is dimensioned
c        b(n,ncols)
c  rhs - the vector (user-supplied) to which the subroutine has
c        applied the same transformation as to the matrix b
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
        ifpivot=1
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
        d=rhs(ipivot)
        rhs(ipivot)=rhs(i)
        rhs(i)=d
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
        call mvecscap(b(1,i),b(1,j),n,d)
c 
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*d
 2770 continue
c 
        rhs(i)=rhs(i)-d*rhs(j)
c 
 2780 continue
 2790 continue
c 
c       normalize the i-th column
c 
        call mvecscap(b(1,i),b(1,i),n,d)
c 
c       if the norm of the remaining part of the matrix
c       is sufficiently small - terminate the process
c 
        if( (ifpivot .eq. 1) .and.
     1      (d .lt. thresh) ) return
c 
        ncols=i
c 
        d=done/dsqrt(d)
        do 2800 j=1,n
        b(j,i)=b(j,i)*d
 2800 continue
c 
        rhs(i)=rhs(i)*d
c 
c        orthogonalize everything else to it
c 
        do 3200 j=i+1,m
c 
        if(rnorms(j) .lt. thresh) goto 3200
c 
        call mvecscap(b(1,i),b(1,j),n,d)
c 
         rrn=0
        do 3000 l=1,n
        b(l,j)=b(l,j)-b(l,i)*d
        rrn=rrn+b(l,j)**2
 3000 continue
c 
        rhs(j)=rhs(j)-d*rhs(i)
c 
        rnorms(j)=rrn
 3200 continue
 3400 continue
 4000 continue
c 
 4200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine mvecscap(x,y,n,prod)
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
        subroutine mvecudv(a,n,u,eps,numrows,nn)
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
        call mveoneuv(b,rlam,u0,v0)
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
       call mvesvsrt(a(1,1),n,u,v,a(1,2),ifuv,nn)
        return
        end
c 
c 
c 
c 
c 
        subroutine mvesvsrt(s,n,u,v,iw,ifuv,nn)
        implicit real *8 (a-h,o-z)
        save
        dimension u(nn,n),v(n,n),s(1),iw(1)
c 
c        if some of the alleged singular values
c        are negative - change their sign, and that
c        of the corresponding left eigenvectors
c 
        do 1100 i=1,n
        if(s(i) .ge. 0) goto 1100
        s(i)=-s(i)
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
c 
c       now, reorder the rows of u and columns of v
c 
        do 2400 i=1,n-1
        ii=iw(i)
c 
c 
c       exchange the rows number i and ii in the matrix u
c 
        do 2200 j=1,nn
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
        subroutine mveoneuv(a,rlam,u,v)
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
        call mveprod2(w,a,b)
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
        call mveprod2(v,w,u)
c 
c       finally, compute the resulting diagonal elements
c 
        call mveprod2(u,a,z)
        call mveprod2(z,vstar,rlams)
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
        subroutine mveprod2(a,b,c)
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
c 
c 
c 
c 
c 
	subroutine mverotfn(a,u)
	implicit real *8 (a-h,o-z)
	real *8 a(2),u(2,2)
c 
	u22=-a(1)
	u12=a(2)
c 
	d=sqrt(u22**2+u12**2)
c 
	if(d .eq. 0) then
c 
	    u(2,2)=1
	    u(1,2)=0
	    u(1,1)=1
	    u(2,1)=0
	    return
	endif
c 
	u(2,2)=u22/d
	u(1,2)=u12/d
	u(1,1)=-u(2,2)
	u(2,1)=u(1,2)
	return
	end
  
