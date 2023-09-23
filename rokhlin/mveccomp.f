        implicit real *8 (a-h,o-z)
        dimension a(100 000),rhs(30 000),rnorms(20 000),
     1      sol(30 000),w(1 000 000),t(10 000),rints(10 000),
     2      rints2(10 000),diff(10 000),sss(10 000)
c 
        external creamat,creamat2
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER n'
        READ *,n
        CALL PRINf('n=*',n,1 )
C 
        PRINT *, 'ENTER nblocks'
        READ *,nblocks
        CALL PRINf('nblocks=*',nblocks,1 )
c 
c       construct the compressed matrix
c 
        eps=1.0d-10
        maxblock=20
        lena=1000 000
        lenw=1000 000
        lenrhs=10 000
  
        ncolblck=10
c 
c       construct the array of sigularities
c 
        width=1.5
        thrleft=-0.75
  
        h=width/nblocks
        do 1400 i=1,nblocks
        sss(i)=(i-1)*h+thrleft
 1400 continue
  
        call prin2('sss as created*',sss,nblocks)
  
  
c 
ccc        call mveccomp(ier,creamat,ncolblck,par2,par3,n,
        call mveccomp(ier,creamat2,ncolblck,sss,par3,n,
     1    nblocks,maxblock,eps,a,lena,rnorms,rhs,lenrhs,
     2      ncols,w,lenw,lwused)
  
  
cccc        call mvecrecu(ier,creamat,mmax,eps,boxes,n,
cccc     1    a,la,rnorms,rhs,lrhs,ncols,w,lw,lwused)
  
c 
        call prinf('after mvecrecu, ier=*',ier,1)
        if(ier .ne. 0) stop
        call prinf('after mvecrecu, ncols=*',ncols,1)
        call prin2('after mvecrecu, rhs=*',rhs,ncols)
cccc        call prin2('after mvecrecu, a=*',a,ncols*n)
c 
c       obtain the quadrature weights
c 
        call amatv(a,n,ncols,rhs,sol)
  
        call prin2('after adjmatv, sol=*',sol,n)
c 
        rr=0
        rrr=0
        do 2200 i=1,n
        rr=rr+sol(i)
        rrr=rrr+abs(sol(i))
 2200 continue
  
        call prin2('and sum of weights*',rr,1)
        call prin2('and sum absolute values of of weights*',
     1      rrr,1)
c 
c       test the obtained quadratures
c 
c 
        itype=0
        call legeexps(itype,n,t,u,v,whts)
  
c 
        do 2300 i=1,n
        t(i)=(t(i)+1)/2
 2300 continue
  
        call prin2('in the main program, t=*',t,n)
c 
        ii=0
        done=1
        do 2800 i=1,ncolblck
        do 2600 j=1,nblocks
c 
        d=0
        dd=sss(j)
        do 2400 k=1,n
        d=d+sol(k)*t(k)**((i-1)+dd)
 2400 continue
c 
        ii=ii+1
        rints(ii)=d
        rints2(ii)=1/(i+dd )
c 
        diff(ii)=rints(ii)-rints2(ii)
  
 2600 continue
 2800 continue
  
        call prin2('rints=*',rints,ii)
        call prin2('and rints2=*',rints2,ii)
        call prin2('and diff=*',diff,ii)
  
        call prinf('and lwused=*',lwused,1)
  
  
         call prin2('and again, sum of absolute values of weights*',
     1      rrr,1)
  
  
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine creamat(nummat,a,rhs,n,m,ncolsblk,par2,par3)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,1),rhs(1),par2(1),par3(1),
     1      t(30 000),pols(1000)
c 
        m=ncolsblk
  
        call prinf('in creamat, m=*',m,1)
c 
c        construct a discretization of the interval [0,1]
c 
        done=1
        h=done/n
        do 1200 i=1,n
        t(i)=i*h
 1200 continue
c 
        itype=0
        call legeexps(itype,n,t,u,v,whts)
c 
        do 1400 i=1,n
        t(i)=(t(i)+1)/2
 1400 continue
c 
ccccc        call prin2('in creamat, t as constructed*',t,n)
c 
c       construct the matrices of values of functions at the nodes
c 
        dd=done/nummat
  
        do 2000 i=1,n
c 
cccc        call legepols(t(i),m,pols)
c 
  
c 
        do 1600 j=1,m
c 
cccc        a(i,j)=pols(i)* t(i)**d
        a(i,j)=t(i)**((j-1)+dd)
 1600 continue
c 
 2000 continue
  
c 
        do 2200 i=1,m
        rhs(i)=1/(i+dd )
 2200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine creamat2(nummat,a,rhs,n,m,ncolsblk,sings,par3)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,1),rhs(1),par2(1),par3(1),sings(1),
     1      t(30 000),pols(1000),work(30 000)
c 
        external funpol
c 
c        construct a discretization of the interval [0,1]
c 
        m=ncolsblk
c 
        done=1
        h=done/n
        do 1200 i=1,n
        t(i)=i*h
 1200 continue
c 
        itype=0
        call legeexps(itype,n,t,u,v,whts)
c 
        do 1400 i=1,n
        t(i)=(t(i)+1)/2
 1400 continue
c 
cccc        call prin2('in creamat2, t as constructed*',t,n)
c 
c       construct the matrices of values of functions at the nodes
c 
        do 2000 i=1,n
c 
cccc        call prinf('in creamat2, i=*',i,1)
  
        d=2*t(i)-1
  
cccc        call legepols(2*t(i)-1,ncolsblk,pols)
        call legepols(d,ncolsblk,pols)
c 
        do 1600 j=1,ncolsblk
c 
        a(i,j)=pols(j)* t(i)**sings(nummat)
 1600 continue
c 
 2000 continue
  
  
cccc        call prin2('in creamat2, a=*',a,n*ncolsblk)
  
  
c 
c       construct the vector of integrals of the functions
c 
        aa=0
        bb=1
        eps=1.0d-12
        mm=18
  
  
        call prin2('before adapgaum, sings(nummat)=*',
     1      sings(nummat),1)
  
        d=sings(nummat)
  
cccc        call adapgaum(ier,aa,bb,funpol,n,sings(nummat),par2,
        call adapgaum(ier,aa,bb,funpol,ncolsblk,d,par2,
     1       mm,eps,rhs,maxrec,numint,work)
c 
  
cccc        call prinf('after adapgaum, ier=*',ier,1)
cccc        call prin2('after adapgaum, rhs=*',rhs,ncolsblk)
  
  
  
        return
        end
c 
c 
c 
c        fun(x,n,par1,par2,vals).                            (1)
c 
c 
        subroutine funpol(x,n,s,par2,vals)
        implicit real *8 (a-h,o-z)
        save
        dimension vals(1)
c 
cccc        call prinf('in funpol, n=*',n,1)
cccc        call prin2('in funpol, s=*',s,1)
cccccc        call prin2('in funpol, x=*',x,1)
  
        call legepols(2*x-1,n,vals)
c 
  
cccc        call prin2('after legepols, vals=*',vals,n)
  
        d=x**s
        do 1200 i=1,n
        vals(i)=vals(i)*d
 1200 continue
  
  
cccc        call prin2('exiting funpol, vals=*',vals,n)
  
        return
        end
  
c 
c 
c 
c 
c 
        subroutine mveccomp(ier,creamat,par1,par2,par3,n,
     1    nblocks,maxblock,eps,a,lena,rnorms,rhs,lenrhs,
     2      ncols,w,lenw,lwused)
        implicit real *8 (a-h,o-z)
        save
        dimension a(1),rnorms(1),rhs(1),w(1),par1(1),par2(1),
     1      par3(1)
        external creamat
c 
c       this subroutine compresses the user-provided matrix via a
c       nested pivoted gram-schmidt process, applying the same
c       linear mapping to the user provided vector. Both the matrix
c       and the vector are produced by the user-provided subroutine
c       creamat. The principal output of the subroutine are the
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
c  creamat - the subroutine creating blocks of the matrix to be
c       compressed, and the corresponding parts of the right-hand side.
c       The calling sequence of creamat is
c 
c        creamat(nummat,a,rhs,n,m,par1,par2,par3)
c 
c        The input parameters for creamat:
c 
c    nummat - the sequence number of the submatrix to be generated
c    n - the number of rows in each submatrix (and also in the whole
c        matrix)
c    par1, par2, par3 - three parameters to be used by creamat (real,
c        integer, whatever)
c 
c        The input parameters for creamat:
c 
c    a - the n \times m submatrix of the big mtrix
c    rhs - the part of the right-hand side corresponding to the part
c        of the matrix returned in a
c    m - the number of columns in a
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c  par1,par2,par3 - three parameters to be used by creamat. They
c         can be real, integer, variable, array, whatever.
c  n - the number of points at which each of the functions is tabulated
c  nblocks - the number of blocks to be created by calls to the
c        subroutine creamat.
c  maxblock - the maximum possible number of columns in any of the
c        submatrices to be created by the subroutine creamat
c  eps - the relative precision to which the compression will be
c        performed
c  boxes - the box structure corresponding to the total number of
c        blocks to be created. It has been (hopefully) constructed
c        by a prior call to the subroutine mvecboxs (see)
c  lena - the amount of space (in real *8 words) provided in the
c        output array a; it is used by the subroutine to bomb if there
c        is not enough space allocated. See the description of the
c        output parameter ier below for details.
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
c  a - the n \times ncols matrix containing the cempressed version of
c        the matrix generated in pieces by the subroutine creamat
c        (see description above). The proper way to view it as a
c        collection of ncols vectors.
c  rnorms - the array of norms of pivot columns procduced in the
c        gram-schmidt process; ncols of them
c  rhs - the compressed version of the right-hand side constructed
c        in pieces by creamat together with the matrix compressed in a
c  ncols - the number of columns in a and of elements in rhs, rnorms
c  lwused - the maximum number of elements in the work array w used
c        by the subroutine at any time.
c 
  
  
  
c 
c       allocate space for boxes and build the box structure
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
        call mvecrecu(ier,creamat,par1,par2,par3,maxblock,eps,
     1    w(iboxes),n,a,lena,rnorms,rhs,lenrhs,ncols,w(iw),
     2      lenw2,lwused)
c 
        lwused=lwused+lboxes+1
        return
        end
c 
c 
c 
c 
c 
        subroutine mvecrecu(ier,creamat,par1,par2,par3,mmax,eps,
     1    boxes,n,a,lena,rnorms,rhs,lenrhs,ncols,w,lenw,lwused)
        implicit real *8 (a-h,o-z)
        save
        dimension a(1),rnorms(1),rhs(1),w(1)
        integer *4 boxes(10,1),stack(100)
        external creamat
c 
c       this subroutine compresses the user-provided matrix via a
c       nested pivoted gram-schmidt process, applying the same
c       linear mapping to the user provided vector. Both the matrix
c       and the vector are produced by the user-provided subroutine
c       creamat. The principal output of the subroutine are the
c       compressed version of the matrix (returned in the array
c       a), and the vector rhs after the gram-schmidt operator has
c       been applied to it. The principal anticipated uses for this
c       subroutine are in the Least squares and quadrature environments.
c 
c                    Input parameters:
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c  creamat - the subroutine creating blocks of the matrix to be
c       compressed, and the corresponding parts of the right-hand side.
c       The calling sequence of creamat is
c 
c        creamat(nummat,a,rhs,n,m,par1,par2,par3)
c 
c        The input parameters for creamat:
c 
c    nummat - the sequence number of the submatrix to be generated
c    n - the number of rows in each submatrix (and also in the whole
c        matrix)
c    par1, par2, par3 - three parameters to be used by creamat (real,
c        integer, whatever)
c 
c        The input parameters for creamat:
c 
c    a - the n \times m submatrix of the big mtrix
c    rhs - the part of the right-hand side corresponding to the part
c        of the matrix returned in a
c    m - the number of columns in a
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c  par1,par2,par3 - three parameters to be used by creamat. They
c         can be real, integer, variable, array, whatever.
c  mmax - the maximum possible number of columns in any of the
c        submatrix
c  eps - the relative precision to which the compression will be
c        performed
c  boxes - the box structure corresponding to the total number of
c        blocks to be created. It has been (hopefully) constructed
c        by a prior call to the subroutine mvecboxs (see)
c  n - the number of points at which each of the functions is tabulated
c  lena - the amount of space (in real *8 words) provided in the
c        output array a; it is used by the subroutine to bomb if there
c        is not enough space allocated. See the description of the
c        output parameter ier below for details.
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
c  a - the n \times ncols matrix containing the cempressed version of
c        the matrix generated in pieces by the subroutine creamat
c        (see description above). The proper way to view it as a
c        collection of ncols vectors.
c  rnorms - the array of norms of pivot columns procduced in the
c        gram-schmidt process; ncols of them
c  rhs - the compressed version of the right-hand side constructed
c        in pieces by creamat together with the matrix compressed in a
c  ncols - the number of columns in a and of elements in rhs, rnorms
c  lwused - the maximum number of elements in the work array w used
c        by the subroutine at any time.
c 
c       . . . initialize the recursion
c 
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
  
ccc        call prinf('after mvecinfo, inuse=*',inuse,1)
c 
        irhs=inuse +2
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
        call mvecreac(creamat,par1,par2,par3,nummat,eps,
     1    w(ia),w(irhs),n,m,w(irnorms),ncols)
  
        call prinf('after mvecreac, ncols=*',ncols,1)
  
c 
        if(ncols*n+3 .le. lena) goto 1600
        ier=16
        return
 1600 continue
c 
        call mvecmove(w(irnorms),rnorms,ncols+1)
        call mvecmove(w(irhs),rhs,ncols+1)
        call mvecmove(w(ia),a,ncols*n+2)
c 
        index=boxes(1,ibox)
        call mvecstor(ier,index,a,rhs,n,ncols,rnorms,w,lenw)
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
c       constructed for its first son - onstruct and store it
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
        irhs=inuse +2
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
     1      w(ia2),ncols2,w(irnorms2),w(irhs2),n,eps,w(ia),ncols,
     2    w(irnorms),w(irhs) )
c 
        if(ncols*n+2 .le. lena) goto 3200
        ier=16
        return
 3200 continue
c 
c 
        if(ncols+2 .le. lenrhs) goto 3400
        ier=32
        return
 3400 continue
c 
        call mvecmove(w(irnorms),rnorms,ncols+1)
        call mvecmove(w(irhs),rhs,ncols+1)
        call mvecmove(w(ia),a,ncols*n+2)
c 
        if(istack .eq. 1) return
c 
        call mvecrem(ier)
        call mvecrem(ier)
c 
        index=boxes(1,ibox)
        call mvecstor(ier,index,a,rhs,n,ncols,rnorms,w,lenw)
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
cccc        call prinf('stored the junk; map=*',map,nstored*10)
  
  
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
  
cccc         call prinf('in mvecretr, nstored=*',nstored,1)
  
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
  
cccc         call prinf('in mvecpnt, nstored=*',nstored,1)
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
        subroutine mvecreac(creamat,par1,par2,par3,nummat,
     1    eps,a,rhs,n,m,rnorms,ncols)
        implicit real *8 (a-h,o-z)
        save
        dimension a(1),rhs(1),rnorms(1),par1(1),par2(1),par3(1)
        external creamat
c 
c        construct the matrix number nummat
c 
        call creamat(nummat,a,rhs,n,m,par1,par2,par3)
c 
cccc        call prin2('in mvecreac after creamat, rhs=*',rhs,m)
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
