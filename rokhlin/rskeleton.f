        implicit real *8 (a-h,o-z)
        dimension a(500 000)
c 
         call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
         PRINT *, 'ENTER n'
         READ *,n
         CALL PRINF('n=*',n,1 )
C 
         PRINT *, 'ENTER m'
         READ *,m
         CALL PRINF('m=*',m,1 )
c 
c       test the gram-schmidt procedure
c 
        eps=1.0d-14
  
  
        lw=300000
        call leastes(a,n,m,eps)
  
  
  
  
        stop
        end
c 
c 
c 
c 
c 
  
        subroutine leastes(a,n,m,eps)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(n,m),
     1      tt(10000),cd,rowsout(200 000),colsout(200 000),
     2      aout(200 000),ima,w(1000 000),
     3      expand(500 000),eval(500 000)
        dimension irows(100 000),icols(100 000),iwhich(10)
c 
        data ima/(0.0d0,1.0d0)/
c 
c        construct the test matrices a, c
c 
        done=1
        h=1
        h=h/n
        do 1100 i=1,5000
        tt(i)=i*h+ima/10
        tt(i)=i*h
 1100 continue
c 
        amax=0
        do 1400 i=1,m
        do 1200 j=1,n
  
        a(j,i)=exp(-tt(i)*tt(j))+exp(-tt(i)/tt(j))
cccc        a(j,i)=exp(-tt(i)*tt(j))
c 
 1200 continue
 1400 continue
c 
        eps1=eps
  
  
cccc        call prin2('a as created*',a,n*m)
  
        iwhich(1)=1
        iwhich(2)=1
  
        iwhich(3)=1
  
        iwhich(4)=1
        iwhich(5)=1
        iwhich(6)=1
  
  
  
        call rskeleton3(a,n,m,eps,iwhich,irows,icols,ncols,
     1      aout,expand,eval,rowsout,colsout,errout,w)
  
cccc        call rskeleton(a,n,m,eps,irows,icols,ncols,
cccc     1      aout,expand,eval,rowsout,colsout,errout,w)
  
        call prinf('after rskeleton, irows=*',irows,ncols)
        call prinf('after rskeleton, icols=*',icols,ncols)
        call prin2('errout=*',errout,1)
  
c 
c        multiply colsout (n,ncols) by expand(ncols,m), and compare
c        the result with a.
c 
        call nrleamult(colsout,expand,w,n,ncols,m)
        call rleasubt(a,w,w,n*m)
        call rleascap(w,w,n*m,cd)
        d=cd
        d=sqrt(d)
c 
        call rleascap(a,a,n*m,cd)
        d2=cd
        d2=sqrt(d2)
c 
        errleft=d/d2
c 
        call prin2('after rskeleton, errleft=*',errleft,1)
c 
c        if the user so requested, multiply the obtained
c        eval(n,ncols) by rowsout(ncols,m), and compare the
c        result with a.
c 
        call nrleamult(eval,rowsout,w,n,ncols,m)
        call rleasubt(a,w,w,n*m)
        call rleascap(w,w,n*m,cd)
c 
        d=cd
        d=sqrt(d)
c 
        errright=d/d2
  
        call prin2('after rskeleton, errright=*',errright,1)
  
        return
        end
  
  
  
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c        This is the end of the debugging (testing) code, and the
c        beginning of the skeletonization code proper.
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c        This file contains four user-callable subroutines:
c        rskeleton, rskeleton3, rskel_basic, rskel_utv. Following
c        is a brief description of these four subroutines
c 
c   rskeleton - for the user-specified real n \times m - matrix a,
c       this subroutine finds a subset of ncols rows and a subset
c       of ncols columns, such that each of the three submatrices
c       dimensined rowsout(ncols,m), colsout(n,ncols),
c       aout(ncols,ncols) provides a good approximation to the
c       original matrix a. Then, it constructs the matrices
c       expand, eval, (approximately) solving the equations
c 
c        aout(ncols,ncols) * expand(ncols,m) = rowsout(ncols,m),     (1)
c 
c        colsout(n,ncols) * expand(ncols,m) = a(n,m),                (2)
c 
c        and
c 
c        eval(n,ncols) * aout(ncols,ncols)  = colsout(n,ncols),      (3)
c 
c        eval(n,k) * colsout(k,m) = a,                               (4)
c 
c        respectively
c 
c   rskeleton3 - has the purpose in life identical to that of
c        rskeleton; has an additional input parameter iwhich,
c        permitting the user to suppress some of the calculations
c        that are always performed by rskeleton. Such suppression
c        can be used to save both CPU time and storage.
c 
c   rskel_basic - this is a simplified and accelerated version
c        of the subroutines rskeleton3, rskeleton. the integer
c        vectors irows, icols are its only output: it does not
c        create any of the matrices rowsout, colsout, eval, expand,
c        whatever.
c 
c   rskel_utv - decomposes the user-supplied matrix a(n,m) in
c        the form
c 
c        a(n,m)=u(n,ncols) * t(ncols,ncols) * ( v(m,ncols) )^*,       (1)
c 
c       with the columns of the matrices u, v orthonormal, and
c       t a lower triangular matrix
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        subroutine rskel_basic(b,n,m,eps,irows,icols,ncols,w)
        implicit real *8 (a-h,o-z)
        save
        dimension rnorms(1),irows(1),icols(1)
        real *8 b(n,m),w(1)
c 
c       For the user-specified n \times m - matrix a, this
c       subroutine finds a subset of ncols rows and a subset
c       of ncols columns, such that each of the three submatrices
c       dimensined rowsout(ncols,m), colsout(n,ncols),
c       aout(ncols,ncols) provides a good approximation to the
c       original matrix a.
c 
c    PLEASE NOTE THAT THIS IS A BASIC VERSION OF THE SUBROUTINE
c    RSKELETON3 (SEE). IT IS FASTER AND SIMPLER THAN EITHER
C    RSKELETON3 OR RSKELETON, BUT THE INTEGER VECTORS IROWS,
C    ICOLS ARE ITS ONLY OUTPUT: IT DOES NOT CREATE ANY OF THE
C    MATRICES ROWSOUT, COLSOUT, EVAL, EXPAND, WHATEVER.
c 
c                        Input parameters:
c 
c  a - the original matrix being approximated
c  n,m - dimensions of a
c  eps - the precision to which the decomposition will be
c       constructed
c 
c                        Output parameters:
c 
c  irows - an integer array of length ncols, containing the sequence
c       numbers of "significant" rows of a
c  icols - an integer array of length ncols, containing the sequence
c       numbers of "significant" columns of a
c  ncols - the dimensions of the matrices aout, rowsout, colsout in
c 
c                        Work arrays:
c 
c  w - must be at least 2*n*m real *8 elements long
c 
c        perform the processing in the case when n .leq. m
c 
        if(n .gt. m) goto 2000
c 
c        . . . gram-schmidt the columns
c 
        call rskel_piv2(b,n,m,eps,w,icols,ncols)
  
cccc        call prin2('in rskel_basic, first rnorms=*',w,ncols)
c 
c        transpose the compressed matrix
c 
        call rleastra(b,n,ncols,w)
c 
c        . . . gram-schmidt the columns of the transposed matrix
c 
        call rskel_piv(w,ncols,n,b,irows)
c 
cccc        call prin2('in rskel_basic, second rnorms=*',b,ncols)
c 
        return
c 
 2000 continue
c 
c        perform the processing in the case when n .gt. m
c 
c        . . . transpose the matrix
c 
        call rleastra(b,n,m,w)
c 
c        . . . gram-schmidt the columns of the transposed matrix
c 
        call rskel_piv2(w,m,n,eps,b,irows,ncols)
  
cccc        call prin2('in rskel_basic, first rnorms=*',b,ncols)
c 
c        transpose the compressed matrix
c 
        call rleastra(w,m,ncols,b)
c 
        call rskel_piv(b,ncols,m,w,icols)
  
cccc        call prin2('in rskel_basic, second rnorms=*',w,ncols)
c 
        return
        end
  
  
c 
c 
c 
c 
c 
        subroutine rskel_piv2(b,n,m,eps,rnorms,ipivots,ncols)
        implicit real *8 (a-h,o-z)
        save
        dimension rnorms(1),ipivots(1)
        real *8 b(n,m),cd
c 
c        . . . initialize the array of pivots
c 
        do 1100 i=1,m
        ipivots(i)=i
 1100 continue
c 
c       . . . prepare the array of values of norms
c             of columns
c 
        done=1
        dtot=0
        do 1400 i=1,m
c 
        d=0
        do 1200 j=1,n
        d=d+b(j,i)*(b(j,i))
 1200 continue
        rnorms(i)=sqrt(d)
        dtot=dtot+d
 1400 continue
c 
        thresh=dtot*eps**2
  
        thresh=sqrt(thresh)
c 
c       . . . conduct gram-schmidt iterations
c 
        do 4000 i=1,n
c 
c       find the pivot
c 
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
        cd=b(j,i)
        b(j,i)=b(j,ipivot)
        b(j,ipivot)=cd
 2600 continue
c 
        iijj=ipivots(i)
        ipivots(i)=ipivots(ipivot)
        ipivots(ipivot)=iijj
c 
        d=rnorms(ipivot)
        rnorms(ipivot)=rnorms(i)
        rnorms(i)=d
c 
c       orthogonalize the i-th column to all preceding ones
c 
        if(i .eq. 1) goto 2790
        do 2780 j=1,i-1
c 
        call rleascap(b(1,i),b(1,j),n,cd)
c 
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*cd
 2770 continue
 2780 continue
 2790 continue
c 
c       normalize the i-th column
c 
        call rleascap(b(1,i),b(1,i),n,cd)
c 
        d=cd
        if(d .lt. thresh**2 ) return
c 
        ncols=i
c 
        d=done/dsqrt(d)
        do 2800 j=1,n
        b(j,i)=b(j,i)*d
 2800 continue
c 
        if(i .eq. m) goto 3400
c 
c        orthogonalize everything else to it
c 
        do 3200 j=i+1,m
c 
        if(rnorms(j) .lt. thresh) goto 3200
c 
        call rleascap(b(1,i),b(1,j),n,cd)
c 
        cd=(cd)
c 
        rrn=0
        do 3000 l=1,n
        b(l,j)=b(l,j)-b(l,i)*cd
        rrn=rrn+b(l,j)*(b(l,j))
 3000 continue
        rnorms(j)=dsqrt(rrn)
 3200 continue
 3400 continue
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
        subroutine rskeleton3(a,n,m,eps,iwhich,irows,
     1      icols,ncols,aout,expand,eval,rowsout,colsout,
     2      errout,w)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(n,m),rowsout(1),colsout(1),aout(1),
     1      expand(1),eval(1),w(1)
        dimension irows(1),icols(1),ifmatr(4),iwhich(10)
c 
c       For the user-specified n \times m - matrix a, this
c       subroutine finds a subset of ncols rows and a subset
c       of ncols columns, such that each of the three submatrices
c       dimensined rowsout(ncols,m), colsout(n,ncols),
c       aout(ncols,ncols) provides a good approximation to the
c       original matrix a. Then, it constructs the matrices
c       expand, eval, (approximately) solving the equations
c 
c        aout(ncols,ncols) * expand(ncols,m) = rowsout(ncols,m),     (1)
c 
c        colsout(n,ncols) * expand(ncols,m) = a(n,m),                (2)
c 
c        and
c 
c        eval(n,ncols) * aout(ncols,ncols)  = colsout(n,ncols),      (3)
c 
c        eval(n,ncols) * rowsout(ncols,m) = a,                       (4)
c 
c        respectively.
c 
c    PLEASE NOTE THAT THIS IS A MORE EFFICIENT VERSION OF THE
C    SUBROUTINE RSKELETON (SEE). IT IS DIFFERENT FROM RSKELETON
C    IN THAT IT HAS THE INPUT PARAMETER IWHICH. IF THE EFFICIENCY
c    (EITHER IN THE USE OF CPU TIME OR OF MEMORY) IS NOT A
c    CONSIDERATION, THE USER MIGHT CONSIDER USING RSKELETON.
c 
c                        Input parameters:
c 
c  a - the original matrix being approximated
c  n,m - dimensions of a
c  eps - the precision to which the decomposition will be
c       constructed
c  iwhich - an integer array of length 6 controlling several aspects
c       of the operation of the subroutine. Each entry can be either
c       0 or 1, to be interpreted as follows:
c 
c    iwhich(1)=1 means that expand will be returned
c    iwhich(2)=1 means that eval will be returned
c    iwhich(3)=1 means that errout will be returned
c    iwhich(4)=1 means that rowsout will be returned
c    iwhich(5)=1 means that colsout will be returned
c    iwhich(6)=1 means that aout will be returned
c 
c 
c                        Output parameters:
c 
c  irows - an integer array of length ncols, containing the sequence
c       numbers of "significant" rows of a
c  icols - an integer array of length ncols, containing the sequence
c       numbers of "significant" columns of a
c  ncols - the dimensions of the matrices aout, rowsout, colsout in
c       the formulae (1)-(4) above
c  aout - the (ncols,ncols) "significant" submatrix of the matrix a,
c       consisting of intersections of rows in the array rowsout and
c       columns in the array colsout. The most compressed form of a;
c       returned only if the inwhich(5) had been set to 1
c  rowsout - the (ncols,m) "significant" submatrix of the matrix a,
c       consisting of ncols rows of the latter; produced via a prior
c       call to rskel_extract; returned only if the inwhich(4) had
c       been set to 1;
c  colsout - the (ncols,n) "significant" submatrix of the matrix a;
c       returned only if the inwhich(5) had been set to 1
c  expand - the (ncols,m) matrix solving equations (1), (2) above;
c       returned only if the inwhich(1) had been set to 1
c  eval - the (n,ncols) matrix solving equations (3), (4) above;
c       returned only if the inwhich(2) had been set to 1
c  errout - the accuracy of the obtained approximation
c       returned only if the inwhich(3) had been set to 1;
c 
c                        Work arrays:
c 
c  w - must be at least 14*n*m+4*(m+n)+1000 real *8 elements long.
c       PLEASE NOTE THAT THIS IS A VERY CONSERVATIVE ESTIMATE
c 
c 
c        . . . construct the UTV decomposition of the input matrix a
c 
        iu=1
        lu=n*m+20
c 
        iv=iu+lu
        lv=n*m+20
c 
        it=iv+lv
        lt=n*m+20
c 
        irnorms=it+lt
        lrnorms=n+m
c 
        iw=irnorms+lrnorms
c 
        call rskel_utv(a,n,m,w(iu),w(it),w(iv),ncols,
     1      w(irnorms),eps,w(iw) )
c 
c        extract the "significant" coordinates from both
c        source and target spaces
c 
        ifmatr(1)=1
        ifmatr(2)=1
        ifmatr(3)=1
c 
c       . . . compress u,v
c 
        lu=n*ncols+20
c 
        iv2=iu+lu
        lv=m*ncols*2+20
c 
        call rleascop(w(iv),w(iv2),lv)
        iv=iv2
c 
        irowsnorms=iv+lv
        lrowsnorms=n+m
c 
        icolsnorms=irowsnorms+lrowsnorms
        lcolsnorms=n+m
c 
        irowsout=icolsnorms+lcolsnorms
        lrowsout=ncols*m+100
c 
        icolsout=irowsout+lrowsout
        lcolsout=ncols*n+100
c 
        iaout=icolsout+lcolsout
        laout=ncols**2+100
c 
        iw=iaout+laout
  
  
        call rskel_extract(w(iu),w(iv),n,m,ncols,irows,
     1      w(irowsnorms),icols,w(icolsnorms),ifmatr,a,
     2      w(irowsout),w(icolsout),w(iaout),w(iw) )
  
c 
c       collect garbage
c 
        irowsout2=1
        call rleascop(w(irowsout),w(irowsout2),lrowsout)
        irowsout=irowsout2
c 
        icolsout2=irowsout+lrowsout
        call rleascop(w(icolsout),w(icolsout2),lcolsout)
        icolsout=icolsout2
c 
        iaout2=icolsout+lcolsout
        call rleascop(w(iaout),w(iaout2),laout)
        iaout=iaout2
c 
c        construct the matrices expand, eval
c 
        irnorms=iaout+laout
        lrnorms=n+m
c 
        iw=irnorms+lrnorms
c 
        ifcheck =iwhich(3)
        errleft=0
        errright=0
c 
        call rskel_evalexp(a,w(iaout),ncols,n,m,w(irowsout),
     1      w(icolsout),eps,expand,eval,w(iw),w(irnorms),ifcheck,
     2      errleft,errright,iwhich)
c 
        call prin2('errleft=*',errleft,1)
        call prin2('errright=*',errright,1)
c 
        errout=errright
        if(errout .lt. errleft) errout=errleft
c 
c        if the user so requested, return the matrices
c        rowsout, colsout, aout
c 
        if(iwhich(4) .eq. 1)
     1      call rleascop(w(icolsout),colsout,ncols*n)
c 
        if(iwhich(5) .eq. 1)
     1      call rleascop(w(irowsout),rowsout,ncols*m)
c 
        if(iwhich(6) .eq. 1)
     1      call rleascop(w(iaout),aout,ncols**2)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine rskeleton(a,n,m,eps,irows,icols,ncols,
     1      aout,expand,eval,rowsout,colsout,errout,w)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(n,m),rowsout(1),colsout(1),aout(1),
     1      expand(1),eval(1),w(1)
        dimension irows(1),icols(1),ifmatr(4),iwhich(10)
c 
c       For the user-specified n \times m - matrix a, this
c       subroutine finds a subset of ncols rows and a subset
c       of ncols columns, such that each of the three submatrices
c       dimensined rowsout(ncols,m), colsout(n,ncols),
c       aout(ncols,ncols) provides a good approximation to the
c       original matrix a. Then, it constructs the matrices
c       expand, eval, (approximately) solving the equations
c 
c        aout(ncols,ncols) * expand(ncols,m) = rowsout(ncols,m),     (1)
c 
c        colsout(n,ncols) * expand(ncols,m) = a(n,m),                (2)
c 
c        and
c 
c        eval(n,ncols) * aout(ncols,ncols)  = colsout(n,ncols),      (3)
c 
c        eval(n,k) * colsout(k,m) = a,                               (4)
c 
c        respectively
c 
c 
c                        Input parameters:
c 
c  a - the original matrix being approximated
c  n,m - dimensions of a
c  eps - the precision to which the decomposition will be
c       constructed
c 
c                        Output parameters:
c 
c  irows - an integer array of length ncols, containing the sequence
c       numbers of "significant" rows of a
c  icols - an integer array of length ncols, containing the sequence
c       numbers of "significant" columns of a
c  ncols - the dimensions of the matrices aout, rowsout, colsout in
c       the formulae (1)-(4) above
c  aout - the (ncols,ncols) "significant" submatrix of the matrix a,
c       consisting of intersections of rows in the array rowsout and
c       columns in the array colsout. The most compressed form of a;
c       produced via a prior call to rskel_extract
c  rowsout - the (ncols,m) "significant" submatrix of the matrix a,
c       consisting of ncols rows of the latter; produced via a prior
c       call to rskel_extract
c  colsout - the (ncols,n) "significant" submatrix of the matrix a,
c  expand - the (ncols,m) matrix solving equations (1), (2) above
c  eval - the (n,ncols) matrix solving equations (3), (4) above
c  errout - the accuracy of the obtained approximation
c 
c                        Work arrays:
c 
c  w - must be at least max( 8*ncols**2+n+m+700, 2*n*m+20+n+m+700)
c        real *8 elements long
c 
c 
c        . . . construct the UTV decomposition of the input matrix a
c 
        iu=1
        lu=n*m+20
c 
        iv=iu+lu
        lv=n*m+20
c 
        it=iv+lv
        lt=n*m+20
c 
        irnorms=it+lt
        lrnorms=n+m
c 
        iw=irnorms+lrnorms
c 
        call rskel_utv(a,n,m,w(iu),w(it),w(iv),ncols,
     1      w(irnorms),eps,w(iw) )
c 
c        extract the "significant" coordinates from both
c        source and target spaces
c 
        ifmatr(1)=1
        ifmatr(2)=1
        ifmatr(3)=1
c 
c       . . . compress u,v
c 
        lu=n*ncols+20
c 
        iv2=iu+lu
        lv=m*ncols*2+20
c 
        call rleascop(w(iv),w(iv2),lv)
        iv=iv2
c 
        irowsnorms=iv+lv
        lrowsnorms=n+m
c 
        icolsnorms=irowsnorms+lrowsnorms
        lcolsnorms=n+m
c 
        call rskel_extract(w(iu),w(iv),n,m,ncols,irows,
     1         w(irowsnorms),icols,w(icolsnorms),ifmatr,a,
     2      rowsout,colsout,aout,w(iw) )
c 
c        construct the matrices expand, eval
c 
        irnorms=1
        lrnorms=n+m
c 
        iw=irnorms+lrnorms
  
  
        ifcheck =1
  
  
        iwhich(1)=1
        iwhich(2)=1
  
  
        call rskel_evalexp(a,aout,ncols,n,m,rowsout,
     1      colsout,eps,expand,eval,w(iw),w(irnorms),ifcheck,
     2      errleft,errright,iwhich)
c 
        call prin2('errleft=*',errleft,1)
        call prin2('errright=*',errright,1)
  
c 
        errout=errright
        if(errout .lt. errleft) errout=errleft
  
        return
        end
c 
c 
c 
c 
c 
        subroutine rskel_evalexp(a,aout,ncols,n,m,rowsout,
     1      colsout,eps,expand,eval,w,rnorms,ifcheck,
     2      errleft,errright,iwhich)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(n,m),
     1      cd,rowsout(ncols,m),colsout(n,ncols),
     2      aout(ncols,ncols),w(1),expand(1),eval(1)
        dimension rnorms(1),iwhich(2)
c 
c        This subroutine uses least squares to construct the
c        matrices expand, eval solving the equations
c 
c        aout(ncols,ncols) * expand(ncols,m) = rowsout(ncols,m),     (1)
c 
c        eval(n,ncols) * aout(ncols,ncols)  = colsout(n,ncols).      (2)
c 
c        According to the so-called thery, the solution expand
c        of the equation (1) should also solve the equation
c 
c        colsout(n,ncols) * expand(ncols,m) = a(n,m),                (3)
c 
c        and the solution eval of the equation (2) should also
c        solve the equation
c 
c        eval(n,k) * colsout(k,m) = a;                               (4)
c 
c        depending on the value of the parameter ifcheck, this
c        subroutine will also test the obtained solution by
c        substituting the obtained expand, eval into (3), (4),
c        respectively.
c 
c        This subroutine uses the results of prior calls to the
c        subroutines rskel_utv, rskel_extract; it has no known
c        uses as a stand-alone device.
c 
c                        Input parameters:
c 
c  a - the original matrix being approximated
c  aout - the (ncols,ncols) "significant" submatrix of the matrix a,
c       consisting of intersections of rows in the array rowsout and
c       columns in the array colsout. The most compressed form of a;
c       produced via a prior call to rskel_extract
c  ncols - the dimensions of the matrices u,v in (1); produced
c       via a prior call to rskel_utv
c  n,m - dimensions of a and other matrices
c  rowsout - the (ncols,m) "significant" submatrix of the matrix a,
c       consisting of ncols rows of the latter; produced via a prior
c       call to rskel_extract
c  colsout - the (ncols,n) "significant" submatrix of the matrix a,
c       consisting of ncols columns of the latter; produced via a prior
c       call to rskel_extract
c  eps - the precision to which the decomposition will be
c       constructed
c  ifcheck - integer parameter telling the subroutine whether to
c       test the accuracy of the obtained approximation;
c    ifcheck=0 will cause the subroutine to skip the checking
c    ifcheck=1 will cause the subroutine to perform the checking
c 
c                        Output parameters:
c 
c  expand - the (ncols,m) matrix solving the equation (1) above
c  eval - the (n,ncols) matrix solving the equation (2) above
c  errleft - the accuracy to which the matrix expand solves the
c            equation (3). Only calculated if the parameter ifcheck
c            had been set to 1
c  errright - the accuracy to which the matrix expand solves the
c            equation (4). Only calculated if the parameter ifcheck
c            had been set to 1
c 
c                        Work arrays:
c 
c  w - must be at least max( 8*ncols**2 +500, 2*n*m+20 ) real *8
c            elements long
c  rnorms - must be at least ncols+1 real *8 elements long
c 
c 
c       . . . if the user so requested, solve in the least
c             squares sense the equation
c 
c        aout(ncols,ncols) * expand(ncols,m) = rowsout(ncols,m)
c 
  
  
cccc        call prinf('iwhich=*',iwhich,2)
  
        if(iwhich(1) .eq. 0) goto 3000
c 
        ifcheck2=0
        call nrleamatll(aout,ncols,ncols,m,rowsout,expand,
     1      eps,ncols2,rnorms,w,ifcheck2,errl2l,errmaxl)
c 
c        if the user so requested, multiply colsout (n,ncols)
c        by the the obtained, expand(ncols,m), and compare
c        the result with a.
c 
        if(ifcheck .eq. 0) goto 3000
c 
        call nrleamult(colsout,expand,w,n,ncols,m)
        call rleasubt(a,w,w,n*m)
        call rleascap(w,w,n*m,cd)
        d=cd
        d=sqrt(d)
c 
        call rleascap(a,a,n*m,cd)
        d2=cd
        d2=sqrt(d2)
c 
        errleft=d/d2
c 
 3000 continue
c 
c       if the user so requested, solve in the least squares
c       sense the equation
c 
c        eval(n,ncols) * aout(ncols,ncols)  = colsout(n,ncols)
c 
        if(iwhich(2) .eq. 0) return
c 
        ifcheck2=0
        call nrleamatrr(aout,n,ncols,ncols,colsout,eval,
     1    eps,ncols2,rnorms,w,ifcheck2,errl2r,errmaxr)
c 
c        if the user so requested, multiply the obtained
c        eval(n,ncols) by rowsout(ncols,m), and compare the
c        result with a.
c 
        if(ifcheck .eq. 0) return
c 
        call nrleamult(eval,rowsout,w,n,ncols,m)
        call rleasubt(a,w,w,n*m)
        call rleascap(w,w,n*m,cd)
c 
        d=cd
        d=sqrt(d)
c 
        call rleascap(a,a,n*m,cd)
        d2=cd
        d2=sqrt(d2)
  
        errright=d/d2
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine rskel_extract(u,v,n,m,ncols,irows,rowsnorms,
     1         icols,colsnorms,ifmatr,a,rowsout,colsout,aout,w)
        implicit real *8 (a-h,o-z)
        save
        real *8 u(n,ncols),v(m,ncols),rowsout(ncols,m),a(n,m),
     1      colsout(n,ncols),aout(ncols,ncols),w(1)
        dimension irows(1),icols(1),rowsnorms(1),
     1      colsnorms(1),ifmatr(3)
c 
c       For the user-specified n \times m - matrix a, this
c       subroutine finds a subset of ncols rows and a subset
c       of ncols columns, such that each of the three submatrices
c       dimensined rowsout(ncols,m), colsout(n,ncols),
c       aout(ncols,ncols) provides a good approximation to the
c       original matrix a. As its input, this subroutine uses the
c       output of the subroutine rskel_utv (see). Depending on the
c       value of the user-supplied parameter ifmatr (see), this
c       subroutine will either return the sequence numbers of
c       such "significant" rows and columns (returned in the arrays
c       irows, icols, respectively), or will also return some or all
c       of the submatrices rowsout, colsout, aout.
c 
c                       Input parameters:
c 
c 
c  u,v - the matrices in the expansion of the form
c 
c        a(n,m)=u(n,ncols) * t(ncols,ncols) * ( v(m,ncols) )^*,       (1)
c 
c       produced via a prior call to the subroutine rskel_utv.
c  n,m - dimensions of a; must be the same as during the preceding
c       call to rskel_utv
C  ncols - the dimensions of the matrices u,v in (1); produced
c       via a prior call to rskel_utv
c  ifmatr - an integer 3-element array telling the subroutine which
c       of the submatrices rowsout, colsout, aout are to be returned
c       (in addition to the integer arrays irows, icols containing
c       the sequence numbers of the rows and columns of which the
c       said submatrices consist). Specifically,
c    ifmatr(1)=0 means that the array rowsout will NOT be constructed
c    ifmatr(2)=0 means that the array colsout will NOT be constructed
c    ifmatr(3)=0 means that the array aout will NOT be constructed
c       Any non-zero value of ifmatr(1),ifmatr(2),ifmatr(3) will
c       cause the corresponding matrix to be constructed.
c  a - the original matrix (see (1) above) about which the whole fuss
c       is happenning. PLEASE NOTE THAT IF ALL THREE ELEMENTS OF THE
C       ARRAY IFMATR HAVE BEEN SET TO ZERO, THIS PARAMETER IS NOT USED.
C       IN OTHER WORDS, THIS PARAMETER IS ONLY USED FOR THE CONSTRUCTION
C       OF THE SUBMATRICES ROWSOUT, COLSOUT, AOUT, BUT NOT FOR THE
C       CONSTRUCTION OF THE PARAMETERS IROWS, ICOLS (SEE BELOW).
c 
c                       Output parameters:
c 
c  irows - an integer array of length ncols, containing the sequence
c       numbers of "significant" rows of a
c  icols - an integer array of length ncols, containing the sequence
c       numbers of "significant" columns of a
c  rowsnorms - the array "rnorms" produced by the pivoted Gram-Schmidt
c       procedure while it was choosing the elements of the array
c       irows
c  colsnorms - the array "rnorms" produced by the pivoted Gram-Schmidt
c       procedure while it was choosing the elements of the array
c       icols
c  rowsout - the (ncols,m) "significant" submatrix of the matrix a,
c       consisting of ncols rows of the latter; the sequence numbers
c       of these rows are returned in the array irows. Only returned
c       if ifmatr(2) had been set to 1
c  colsout - the (ncols,n) "significant" submatrix of the matrix a,
c       consisting of ncols columns of the latter; the sequence numbers
c       of these columns are returned in the array icols. Only returned
c       if ifmatr(2) had been set to 1
  
c  aout - the (ncols,ncols) "significant" submatrix of the matrix a,
c       consisting of intersections of rows in the array rowsout and
c       columns in the array colsout. Only returned if ifmatr(3) had
c       been set to 1
c 
c                       Work arrays:
c 
c  w - must be at least n*m+2*(n+m) real *8 locations long
c 
c        extract the "significant" coordinates from both
c        source and target spaces
c 
        call rleastra(u,n,ncols,w)
c 
        iw=1
        lw=ncols*n+100
        if(m .gt. n) lw=ncols*m+100
c 
        irowsnorms=iw+lw
        lrowsnorms=n+m+2
c 
        iirows=irowsnorms+lrowsnorms
c 
        call rskel_piv(w,ncols,n,w(irowsnorms),w(iirows) )
c 
        call rleascop(w(irowsnorms),rowsnorms,ncols)
        call rskel_intcopy(w(iirows),irows,ncols)
  
        icolsnorms=iw+lw
        lcolsnorms=n+m+2
c 
        iicols=icolsnorms+lcolsnorms
c 
        call rleastra(v,m,ncols,w)
        call rskel_piv(w,ncols,m,w(icolsnorms),w(iicols) )
  
        call rleascop(w(icolsnorms),colsnorms,ncols)
        call rskel_intcopy(w(iicols),icols,ncols)
c 
c        if the user so requested, return to him the parts
c        of the matrix defined by the vectors irows,
c        icols
c 
        if(ifmatr(1) .eq. 0) goto 2600
c 
        do 2400 i=1,m
        do 2200 j=1,ncols
c 
        jj=irows(j)
        rowsout(j,i)=a(jj,i)
 2200 continue
 2400 continue
c 
 2600 continue
c 
        if(ifmatr(2) .eq. 0) goto 3600
c 
        do 3400 i=1,n
        do 3200 j=1,ncols
c 
        jj=icols(j)
        colsout(i,j)=a(i,jj)
 3200 continue
 3400 continue
c 
 3600 continue
c 
        if(ifmatr(3) .eq. 0) return
c 
        do 4400 i=1,ncols
        do 4200 j=1,ncols
c 
        jj=irows(j)
        ii=icols(i)
        aout(j,i)=a(jj,ii)
 4200 continue
 4400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine rskel_intcopy(ia,ib,n)
        implicit real *8 (a-h,o-z)
        save
        dimension ia(1),ib(1)
c 
        do 1200 i=1,n
c 
        ib(i)=ia(i)
 1200 continue
  
        return
        end
c 
c 
c 
c 
c 
        subroutine rskel_utv(a,n,m,u,t,v,ncols,rnorms,eps,w)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(n,m),u(n,m),v(m,n),w(1),t(1)
        dimension rnorms(1)
c 
c        This subroutine decomposes the user-supplied matrix
c        a(n,m) in the form
c 
c        a(n,m)=u(n,ncols) * t(ncols,ncols) * ( v(m,ncols) )^*,       (1)
c 
c       with the columns of the matrices u, v orthonormal, and
c       t a lower triangular matrix
c 
c 
c                 Input parameters:
c 
c  a - the matrix to be decomposed
c  n,m - dimensions of a
c  eps - the precision to which the decomposition will be
c       constructed
c 
c                 Output parameters:
c 
c  u,t,v - the matrices in (1) above. PLEASE NOTE THAT SINCE THE
C          PARAMETER NCOLS IS NOT KNOWN APRIORI, THE ARRAYS U,T,V
C          SHOULD BE OF LENGTH AT LEAST N*M*2 REAL *8 ELEMENTS.
C 
C  ncols - the dimensions of the matrices u,t,v in (1)
c 
c                 Work arrays:
c 
c  w - must be at least n*m*2 real *8 elements long.
c 
c 
c        . . . construct the decomposition
c 
c        A(n,m)= u(n,ncols) * ( w(m,ncols) )^*
c 
        ifpivot=1
        call rleaspiv(a,n,m,u,w,ncols1,v,eps,ifpivot)
  
cccc        call prin2('and rnorms=*',v,ncols1)
c 
c        construct the decomposition
c 
c        ( w(m,ncols1) ) = v(m,ncols) * ( t(ncols,ncols) )
c 
        ifpivot=0
        call rleaspiv(w,m,ncols1,v,t,ncols,rnorms,eps,ifpivot)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine rskel_piv(b,n,m,rnorms,ipivots)
        implicit real *8 (a-h,o-z)
        save
        dimension rnorms(1),ipivots(1)
        real *8 b(n,m),cd
c 
c       This subroutine extracts from the rectangular matrix b
c       (user-supplied) a square full-rank submatrix. The number
c       of columns m of b is expected to be greater than the number
c       n of its columns. The subroutine also returns the integer
c       array ipivots, containing the sequence numbers of columns
c       selected, and the real array rnorms, containing the
c       denominators by which the selected columns had to be
c       divided during the Gram-Schmidt process (a decent measure
c       of stability of the process).
c 
c                    input paramneters:
c 
c  b - the matrix to be gram-schmidt'ed. it is destroyed by this
c       subroutine
c  n,m - dimensionalities of the matrix b
c 
c                     output parameters:
c 
c  b - the matrix of gram-schmidt vectors of the matrix a. note
c        that on exit from this subroutine only the first
c        n columns of  b  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit b is dimensioned
c        b(n,n)
c  rnorms - the normalizing factors in the gram-schmidt process.
c        Only the first n of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c 
c 
c        . . . initialize the array of pivots
c 
        do 1100 i=1,m
        ipivots(i)=i
 1100 continue
c 
c       . . . prepare the array of values of norms
c             of columns
c 
        done=1
        do 1400 i=1,m
c 
        d=0
        do 1200 j=1,n
        d=d+b(j,i)*(b(j,i))
 1200 continue
        rnorms(i)=sqrt(d)
 1400 continue
c 
c       . . . conduct gram-schmidt iterations
c 
        do 4000 i=1,n
c 
c       find the pivot
c 
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
        cd=b(j,i)
        b(j,i)=b(j,ipivot)
        b(j,ipivot)=cd
 2600 continue
c 
        iijj=ipivots(i)
        ipivots(i)=ipivots(ipivot)
        ipivots(ipivot)=iijj
c 
        d=rnorms(ipivot)
        rnorms(ipivot)=rnorms(i)
        rnorms(i)=d
c 
c       orthogonalize the i-th column to all preceding ones
c 
        if(i .eq. 1) goto 2790
        do 2780 j=1,i-1
c 
        call rleascap(b(1,i),b(1,j),n,cd)
c 
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*cd
 2770 continue
 2780 continue
 2790 continue
c 
c       normalize the i-th column
c 
        call rleascap(b(1,i),b(1,i),n,cd)
c 
        d=cd
cccc        ncols=i
c 
        d=done/dsqrt(d)
        do 2800 j=1,n
        b(j,i)=b(j,i)*d
 2800 continue
c 
        if(i .eq. m) goto 3400
c 
c        orthogonalize everything else to it
c 
        do 3200 j=i+1,m
c 
        call rleascap(b(1,i),b(1,j),n,cd)
c 
        cd=(cd)
c 
        rrn=0
        do 3000 l=1,n
        b(l,j)=b(l,j)-b(l,i)*cd
        rrn=rrn+b(l,j)*(b(l,j))
 3000 continue
        rnorms(j)=dsqrt(rrn)
 3200 continue
 3400 continue
c 
 4000 continue
c 
        return
        end
