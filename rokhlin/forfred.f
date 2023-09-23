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
        complex *16 a(n,m),
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
cccc        tt(i)=i*h
 1100 continue
c 
        amax=0
        do 1400 i=1,m
        do 1200 j=1,n
  
  
        a(j,i)=cdexp(-tt(i)*tt(j))+cdexp(-tt(i)/tt(j))
cccc        a(j,i)=cdexp(-tt(i)*tt(j))
c 
 1200 continue
 1400 continue
c 
        eps1=eps
  
  
cccc        call prin2('a as created*',a,n*m*2)
  
        iwhich(1)=1
        iwhich(2)=1
  
        iwhich(3)=1
  
        iwhich(4)=1
        iwhich(5)=1
        iwhich(6)=1
  
  
  
        call cskeleton3(a,n,m,eps,iwhich,irows,icols,ncols,
     1      aout,expand,eval,rowsout,colsout,errout,w)
  
cccc        call cskeleton(a,n,m,eps,irows,icols,ncols,
cccc     1      aout,expand,eval,rowsout,colsout,errout,w)
  
        call prinf('after cskeleton, irows=*',irows,ncols)
        call prinf('after cskeleton, icols=*',icols,ncols)
        call prin2('errout=*',errout,1)
  
c 
c        multiply colsout (n,ncols) by expand(ncols,m), and compare
c        the result with a.
c 
        call ncleamult(colsout,expand,w,n,ncols,m)
        call cleasubt(a,w,w,n*m)
        call cleascap(w,w,n*m,cd)
        d=cd
        d=sqrt(d)
c 
        call cleascap(a,a,n*m,cd)
        d2=cd
        d2=sqrt(d2)
c 
        errleft=d/d2
c 
        call prin2('after cskeleton, errleft=*',errleft,1)
c 
c        if the user so requested, multiply the obtained
c        eval(n,ncols) by rowsout(ncols,m), and compare the
c        result with a.
c 
        call ncleamult(eval,rowsout,w,n,ncols,m)
        call cleasubt(a,w,w,n*m)
        call cleascap(w,w,n*m,cd)
c 
        d=cd
        d=sqrt(d)
c 
        errright=d/d2
  
        call prin2('after cskeleton, errright=*',errright,1)
  
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
c        This file contains three user-callable subroutines:
c        cskeleton, cskeleton3, and cskel_utv. Following is a
c        brief description of these two subroutines
c 
c   cskeleton - for the user-specified n \times m - matrix a, this
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
c   cskeleton3 - has the purpose in life identical to that of
c        cskeleton; has an additional input parameter iwhich,
c        permitting the user to suppress some of the calculations
c        that are always performed by cskeleton. Such suppression
c        can be used to save both CPU time and storage.
c 
c   cskel_utv - decomposes the user-supplied matrix a(n,m) in
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
c 
c 
        subroutine cskeleton3(a,n,m,eps,iwhich,irows,
     1      icols,ncols,aout,expand,eval,rowsout,colsout,
     2      errout,w)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,m),rowsout(1),colsout(1),aout(1),
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
C    SUBROUTINE CSKELETON (SEE). IT IS DIFFERENT FROM CSKELETON
C    IN THAT IT HAS THE INPUT PARAMETER IWHICH. IF THE EFFICIENCY
c    (EITHER IN THE USE OF CPU TIME OR OF MEMORY) IS NOT A
c    CONSIDERATION, THE USER MIGHT CONSIDER USING CSKELETON.
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
c       call to cskel_extract; returned only if the inwhich(4) had
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
        call cskel_utv(a,n,m,w(iu),w(it),w(iv),ncols,
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
        call cleascop(w(iv),w(iv2),lv)
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
  
  
        call cskel_extract(w(iu),w(iv),n,m,ncols,irows,
     1      w(irowsnorms),icols,w(icolsnorms),ifmatr,a,
     2      w(irowsout),w(icolsout),w(iaout),w(iw) )
  
c 
c       collect garbage
c 
        irowsout2=1
        call cleascop(w(irowsout),w(irowsout2),lrowsout)
        irowsout=irowsout2
c 
        icolsout2=irowsout+lrowsout
        call cleascop(w(icolsout),w(icolsout2),lcolsout)
        icolsout=icolsout2
c 
        iaout2=icolsout+lcolsout
        call cleascop(w(iaout),w(iaout2),laout)
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
        call cskel_evalexp(a,w(iaout),ncols,n,m,w(irowsout),
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
     1      call cleascop(w(icolsout),colsout,ncols*n)
c 
        if(iwhich(5) .eq. 1)
     1      call cleascop(w(irowsout),rowsout,ncols*m)
c 
        if(iwhich(6) .eq. 1)
     1      call cleascop(w(iaout),aout,ncols**2)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cskeleton(a,n,m,eps,irows,icols,ncols,
     1      aout,expand,eval,rowsout,colsout,errout,w)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,m),rowsout(1),colsout(1),aout(1),
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
c       produced via a prior call to cskel_extract
c  rowsout - the (ncols,m) "significant" submatrix of the matrix a,
c       consisting of ncols rows of the latter; produced via a prior
c       call to cskel_extract
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
        call cskel_utv(a,n,m,w(iu),w(it),w(iv),ncols,
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
        call cleascop(w(iv),w(iv2),lv)
        iv=iv2
c 
        irowsnorms=iv+lv
        lrowsnorms=n+m
c 
        icolsnorms=irowsnorms+lrowsnorms
        lcolsnorms=n+m
c 
        call cskel_extract(w(iu),w(iv),n,m,ncols,irows,
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
  
  
        call cskel_evalexp(a,aout,ncols,n,m,rowsout,
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
        subroutine cskel_evalexp(a,aout,ncols,n,m,rowsout,
     1      colsout,eps,expand,eval,w,rnorms,ifcheck,
     2      errleft,errright,iwhich)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,m),
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
c        subroutines cskel_utv, cskel_extract; it has no known
c        uses as a stand-alone device.
c 
c                        Input parameters:
c 
c  a - the original matrix being approximated
c  aout - the (ncols,ncols) "significant" submatrix of the matrix a,
c       consisting of intersections of rows in the array rowsout and
c       columns in the array colsout. The most compressed form of a;
c       produced via a prior call to cskel_extract
c  ncols - the dimensions of the matrices u,v in (1); produced
c       via a prior call to cskel_utv
c  n,m - dimensions of a and other matrices
c  rowsout - the (ncols,m) "significant" submatrix of the matrix a,
c       consisting of ncols rows of the latter; produced via a prior
c       call to cskel_extract
c  colsout - the (ncols,n) "significant" submatrix of the matrix a,
c       consisting of ncols columns of the latter; produced via a prior
c       call to cskel_extract
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
        call ncleamatll(aout,ncols,ncols,m,rowsout,expand,
     1      eps,ncols2,rnorms,w,ifcheck2,errl2l,errmaxl)
c 
c        if the user so requested, multiply colsout (n,ncols)
c        by the the obtained, expand(ncols,m), and compare
c        the result with a.
c 
        if(ifcheck .eq. 0) goto 3000
c 
        call ncleamult(colsout,expand,w,n,ncols,m)
        call cleasubt(a,w,w,n*m)
        call cleascap(w,w,n*m,cd)
        d=cd
        d=sqrt(d)
c 
        call cleascap(a,a,n*m,cd)
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
        call ncleamatrr(aout,n,ncols,ncols,colsout,eval,
     1    eps,ncols2,rnorms,w,ifcheck2,errl2r,errmaxr)
c 
c        if the user so requested, multiply the obtained
c        eval(n,ncols) by rowsout(ncols,m), and compare the
c        result with a.
c 
        if(ifcheck .eq. 0) return
c 
        call ncleamult(eval,rowsout,w,n,ncols,m)
        call cleasubt(a,w,w,n*m)
        call cleascap(w,w,n*m,cd)
c 
        d=cd
        d=sqrt(d)
c 
        call cleascap(a,a,n*m,cd)
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
        subroutine cskel_extract(u,v,n,m,ncols,irows,rowsnorms,
     1         icols,colsnorms,ifmatr,a,rowsout,colsout,aout,w)
        implicit real *8 (a-h,o-z)
        save
        complex *16 u(n,ncols),v(m,ncols),rowsout(ncols,m),a(n,m),
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
c       output of the subroutine cskel_utv (see). Depending on the
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
c       produced via a prior call to the subroutine cskel_utv.
c  n,m - dimensions of a; must be the same as during the preceding
c       call to cskel_utv
C  ncols - the dimensions of the matrices u,v in (1); produced
c       via a prior call to cskel_utv
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
c  w - must be at least n*m+2*(n+m) complex *16 locations long
c 
c        extract the "significant" coordinates from both
c        source and target spaces
c 
        call cleastra(u,n,ncols,w)
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
        call cskel_piv(w,ncols,n,w(irowsnorms),w(iirows) )
c 
        call cleascop(w(irowsnorms),rowsnorms,ncols)
        call cskel_intcopy(w(iirows),irows,ncols)
  
        icolsnorms=iw+lw
        lcolsnorms=n+m+2
c 
        iicols=icolsnorms+lcolsnorms
c 
        call cleastra(v,m,ncols,w)
        call cskel_piv(w,ncols,m,w(icolsnorms),w(iicols) )
  
        call cleascop(w(icolsnorms),colsnorms,ncols)
        call cskel_intcopy(w(iicols),icols,ncols)
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
        subroutine cskel_intcopy(ia,ib,n)
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
        subroutine cskel_utv(a,n,m,u,t,v,ncols,rnorms,eps,w)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,m),u(n,m),v(m,n),w(1),t(1)
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
        call cleaspiv(a,n,m,u,w,ncols1,v,eps,ifpivot)
  
cccc        call prin2('and rnorms=*',v,ncols1)
c 
c        construct the decomposition
c 
c        ( w(m,ncols1) ) = v(m,ncols) * ( t(ncols,ncols) )
c 
        ifpivot=0
        call cleaspiv(w,m,ncols1,v,t,ncols,rnorms,eps,ifpivot)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cskel_piv(b,n,m,rnorms,ipivots)
        implicit real *8 (a-h,o-z)
        save
        dimension rnorms(1),ipivots(1)
        complex *16 b(n,m),cd
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
        d=d+b(j,i)*dconjg(b(j,i))
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
        call cleascap(b(1,i),b(1,j),n,cd)
c 
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*cd
 2770 continue
 2780 continue
 2790 continue
c 
c       normalize the i-th column
c 
        call cleascap(b(1,i),b(1,i),n,cd)
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
        call cleascap(b(1,i),b(1,j),n,cd)
  
cc        subroutine cleascap(x,y,n,prod)
  
c 
        cd=conjg(cd)
c 
        rrn=0
        do 3000 l=1,n
        b(l,j)=b(l,j)-b(l,i)*cd
        rrn=rrn+b(l,j)*conjg(b(l,j))
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c       this is the end of the debugging code and the beginning of the
c       least squares subroutines proper
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c         This file contains 6 user-callable subroutines: ncleastsq,
c         ncleastsq2, cleamatr, ncleamatrr, ncleamatll, ncleamatrl.
c         Following is a brief description of these subroutines.
c 
c    NOTE: PLEASE NOTE THAT THE SUBROUTINE ncleamatrl IS ALMOST NOT
c         A USER-CALLABLE ONE. WHILE IT IS CONCIEVABLE THAT SOME
c         USERS WILL FIND IT USEFUL, MOST ARE EXPECTED TO BE HAPPY
c         WITH THE SUBROUTINES ncleamatrr, ncleamatll, TO THE EXTENT
c         THAT THEY ARE HAPPY WITH ANYTHING IN THIS FILE.
c 
c   ncleastsq - constructs a QR-type (not exactly QR) decomposition
c         of the input matrix A, to be used by the subroutine
c         ncleastsq2 for the solution of least squares problems of
c         the form
c 
c               A X \sim Y,                                             (1)
c 
c        and by the subroutine ncleamatrl for the least squares
c        solution of the matrix equation
c 
c         A(n,m) * X(m,k) = C(n,k).                                     (2)
c 
c   ncleastsq2 - uses the QR-type (not exactly QR) decomposition
c         of the matrix A produced by a prior call to nsleastsq (see)
c         to solve in the least squares sense the linear system
c 
c                    A X = RHS.                                         (3)
c 
c   cleamatr - solves in the least squares sense the matrix
c        equation
c 
c                     A(k,l) * X(l,m) * B(m,n) = C(k,n),                (3a)
c 
c       where A, B, C are user-specified complex matrices, and X
c       is the matrix to be found. Note that the dimensionalities
c       of the matrices A,B,C, X are as general as could be
c 
c   nclearmatll - solves in the least squares sense the matrix
c         equation
c 
c                    A(n,m) * X(m,k) = C(n,k),                          (4)
c 
c         without using any additional data (i.e. it performs all
c         factorizations itself).
c 
c   nclearmatrr - solves in the least squares sense the matrix
c         equation
c 
c                    X(n,m) * A(m,k) = C(n,k),                          (5)
c 
c         without using any additional data (i.e. it performs all
c         factorizations itself).
c 
c   ncleamatrl - solves in the least squares sense the
c        matrix equation
c 
c                    A(n,m) * X(m,k) = C(n,k),                          (6)
c 
c        using as input the matrix c and the array w obtained
c        via a preceding call to the subroutine ncleastsq
c        (see).
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        subroutine cleamatr(ier,a,b,c,k,l,m,n,eps,delta,x,
     1      w,lw,ltot)
        implicit real *8 (a-h,o-z)
c 
        save
        complex *16  a(k,l),x(l,m),b(m,n),c(k,n),w(1)
c 
c        this subroutine solves in the least squares sense the
c        equation
c 
c                     A(k,l) * X(l,m) * B(m,n) = C(k,n),               (1)
c 
c       where A, B, C are user-specified complex matrices, and X
c       is the matrix to be found. Note that the dimensionalities
c       of the matrices A,B,C, X are as general as could be
c 
c 
c                          input parameters:
c 
c  a,b,c - matrices in (1)
c  k,l,m,n - dimensionalities in (1)
c  eps - the parameter telling the subroutine at which point a
c       a singular value of the matrix is to be declared to be zero
c       and ignored (see subroutine cleasas for details)
c  delta - the "Tikhonov constant" (see subroutine cleasas for details)
c  lw - the length of the user-supplied work array
c 
c                          output parameters:
c 
c  ier - error return code
c      ier=0 means successful execution
c      ier .ne. 0 means that the amount of space lw allocated supplied
c            in the work array w is insufficient
c  x - the unknown matrix in (1)
c  ltot - the length of the work array w actually used by the subroutine
c 
c                          work array:
c 
c  w - must be sufficiently long
c 
c 
c        . . . construct the SVD of A
c 
        ier=0
        isa=1
        lsa=min(k,l)+2
c 
        iua=isa+lsa
        lua=k*l+10
c 
        iva=iua+lua
        lva=k*l+10
c 
        iwork=iva+lva
c 
        eps2=1.0d-13
        lw1=lw-iwork-2
c 
        call prinf('in cleamatr, lw1=*',lw1,1)
c 
  
cccc        subroutine csvdpiv(ier,a,n,m,u,v,s,ncols,eps,
cccc     1      w,lw,ltot)
  
  
        call csvdpiv(ier,a,k,l,w(iua),w(iva),w(isa),ncolsa,eps2,
     1      w(iwork),lw1,ltot1)
c 
        ltot=ltot1+iwork
c 
        if(ier .ne. 0) return
c 
        call prinf('after first svdpivot, ncolsa=*',ncolsa,1)
        call prin2('after first svdpivot, sa=*',w(isa),ncolsa)
c 
c        perform garbage collection
c 
        lua2=k*ncolsa+2
        lva2=l*ncolsa+2
        iva2=iua+lua2
        lua=lua2
c 
        call cleamem(w(iva),w(iva2),lva2)
c 
        iva=iva2
        lva=lva2
c 
c        . . . construct the SVD of B
c 
        isb=iva+lva
        lsb=min(m,n)+2
c 
        iub=isb+lsb
        lub=m*n+10
c 
        ivb=iub+lub
        lvb=m*n+10
c 
        iwork=ivb+lvb
        lw2=lw-iwork
cccc        call prinf('in cleamatr, lw2=*',lw2,1)
c 
        eps2=1.0d-13
c 
        call csvdpiv(ier,b,m,n,w(iub),w(ivb),w(isb),ncolsb,eps2,
     1      w(iwork),lw2,ltot2)
c 
        ii=iwork+ltot2
        if(ii .gt. ltot) ltot=ii
c 
cccc        call prinf('after second svdpivot, ltot=*',ltot,1)
c 
        if(ier .ne. 0) return
c 
        call prinf('after second svdpivot, ncolsb=*',ncolsb,1)
        call prin2('after second svdpivot, sb=*',w(isb),ncolsb)
c 
c        perform garbage collection
c 
        lub2=m*ncolsb+2
        lvb2=n*ncolsb+2
        ivb2=iub+lub2
c 
        call cleamem(w(ivb),w(ivb2),lvb2)
c 
        ivb=ivb2
        lvb=lvb2
c 
        iwork=ivb+lvb
        lwork=ncolsa*n+10
c 
        ii=iwork+lwork
        call prinf('in cleamatr, final ii=*',ii,1)
        if(ii .gt. ltot) ltot=ii
        if(ltot .lt. lw) goto 2200
        ier=2
        return
 2200 continue
c 
c         perform the remainder of the solution
c 
        call cleamat0(c,k,l,m,n,eps,delta,x,
     1      w(isa),w(iua),w(iva),w(isb),w(iub),w(ivb),
     2      w(iwork),lw,ltot,ncolsa,ncolsb)
  
  
        return
        end
  
c 
c 
c 
c 
c 
        subroutine ncleamatrr(a,n,m,k,c,x,eps,ncols,
     1    rnorms,w,ifcheck,errl2,errmax)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1)
        complex *16 x(n,m),c(n,k),a(m,k)
c 
c        This subroutine solves in the least squares sense the
c        matrix equation
c 
c        X(n,m) * A(m,k) = C(n,k).                              (1)
c 
c                     input parameters:
c 
c  a - the matrix to be decomposed (note that it is NOT damaged
c       by this subroutine in any way)
c  n, m, k - the dimensions in (1) above
c  c - the right-hand side above NOT damaged by this subroutine
c  eps - the accuracy to which the decomposition (1) will be computed.
c  ifcheck - integer parameter telling the subroutine whether the
c      error of the approximations should be calculated.
c   ifcheck=1 will cause the subroutine to evaluate (and return)
c      the errors - both maximum and relative l^2
c   ifcheck=0 will cause the subroutine skip the error evaluation.
c      This is a CPU time saving feature.
c 
c                     output parameters:
c 
c  x - the solution of (1) above in the least squares sense
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices u, w on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c  w - the first n*k complex *16 elements of w contain the
c        discrepancies in the approximations of the n*k elements
c        of the matrix c
c  errl2 - the l^2 error of the obtained approximation; in other
c        words, errl2=sqrt(sum_{i=1}^{n*k} |w_i|^2)
c 
c                    work arrays:
c 
c  w - must be at least n*m*8 + 500 real *8 elements long
c 
c 
c       . . . transpose both the matrices a, c
c 
        call cleastra(a,m,k,w)
        call cleascop(w,a,m*k)
c 
        call cleastra(c,n,k,w)
        call cleascop(w,c,n*k)
c 
c       solve least-squares problem
c 
c       A^* (k,m) * x^* (m,n) = C^* (k,n)
c 
        call ncleamatll(a,k,m,n,c,x,eps,ncols,
     1    rnorms,w,ifcheck,errl2,errmax)
c 
c       transpose everything back
c 
        iaia=n*m
        if(k*m .gt. iaia) iaia=k*m
        if(k*n .gt. iaia) iaia=k*n
  
        iaia=iaia*2+20
  
        call cleascop(w,w(iaia),n*k)
c 
        call cleastra(a,k,m,w)
        call cleascop(w,a,m*k)
c 
        call cleastra(c,k,n,w)
        call cleascop(w,c,n*k)
c 
        call cleastra(x,m,n,w)
        call cleascop(w,x,n*m)
c 
        call cleastra(w(iaia),k,n,w)
  
        return
        end
c 
c 
c 
c 
        subroutine ncleamatll(a,n,m,k,c,x,eps,ncols,
     1    rnorms,w,ifcheck,errl2,errmax)
        implicit real *8 (a-h,o-z)
        save
        complex *16 x(m,k),c(n,k),a(n,m),cd,w(1),cd2
c 
c        This subroutine solves in the least squares sense the
c        matrix equation
c 
c         A(n,m) * X(m,k) = C(n,k).                            (1)
c 
c                     input parameters:
c 
c  a - the matrix to be decomposed (note that it is NOT damaged
c       by this subroutine in any way)
c  n, m, k - the dimensions in (1) above
c  c - the right-hand side above, NOT damaged by this subroutine
c  eps - the accuracy to which the decomposition (1) will be computed.
c  ifcheck - integer parameter telling the subroutine whether the
c       error of the approximations should be calculated.
c    ifcheck=1 will cause the subroutine to evaluate (and return)
c       the errors - both maximum and relative l^2
c    ifcheck=0 will cause the subroutine skip the error evaluation.
c       This is a CPU time saving feature.
c 
c                     output parameters:
c 
c  x - the solution of (1) above in the least squares sense
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices u, w on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c  w - the first n*k complex *16 elements of w contain the
c        discrepancies in the approximations of the n*k elements
c        of the matrix c
c  errl2 - the l^2 error of the obtained approximation; in other
c        words, errl2=sqrt(sum_{i=1}^{n*k} |w_i|^2)
c 
c                    work arrays:
c 
c  w - must be at least n*m*8 + 500 real *8 elements long
c 
c 
c        . . . decompose the matrix a
c 
        call ncleastsq(a,n,m,eps,ncols,rnorms,w)
c 
c       obtain the matrix x
c 
        call ncleamatrl(w,c,x,n,m,k)
c 
c        multiply A by X and check how close the result is to C
c 
        if(ifcheck .eq. 0) return
  
        call ncleamult(a,x,w,n,m,k)
        call cleasubt(w,c,w,n*k)
c 
        call cleascap(w,w,n*k,cd)
        call cleascap(c,c,n*k,cd2)
c 
        errl2=sqrt(cd/cd2)
  
        errmax=0
        do 2200 i=1,n*k
c 
        d=w(i)*conjg(w(i))
        if(errmax .lt. d) errmax=d
 2200 continue
c 
        errmax=sqrt(errmax)
c 
        return
        end
  
  
  
c 
c 
c 
c 
c 
        subroutine ncleamatrl(w,c,x,n,m,k)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1)
        complex *16 x(m,k),c(n,k)
c 
c        This subroutine solves in the least squares sense the
c        matrix equation
c 
c         A(n,m) * X(m,k) = C(n,k),                            (1)
c 
c        using as input the matrix c and the array w obtained
c        via a preceding call to the subroutine ncleastsq
c        (see).
c 
c                     input parameters:
c 
c  w - the array containing the decomposition of the matrix a,
c        obtained via a preceding call to ncleastsq  (see)
c  c - the right-hand side above
c  n, m, k - the dimensions in (1) above
c 
c                     output parameters:
c 
c  x - the solution of (1) above in the least squares sense
c 
c 
c        . . . solve the least squares problems for columns of x
c              one after another
c 
        do 1400 i=1,k
c 
        call ncleasts2(w,c(1,i),x(1,i) )
 1400 continue
  
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine ncleastsq(a,n,m,eps,ncols,rnorms,w)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,m)
        dimension rnorms(1),w(1)
c 
c        This subroutine constructs the decomposition of the input
c        matrix a, to be used by the subroutine ncleastsq2 for the
c        solution of least squares problems of the form
c 
c               A X \sim Y,                                          (1)
c 
c        and by the subroutine ncleamatrl for the least squares
c        solution of the matrix equation
c 
c 
c         A(n,m) * X(m,k) = C(n,k).                                  (2)
c 
c        The decomposition is stored in the array w; the first
c        n*m*8 + 500 real *8 elements of the latter should not be
c        changed between the call to this subroutine and the
c        subsequent calls to ncleastsq2.
c 
c   NOTE: this subroutine uses the subroutine cleastsq to perform
c        all of the work; this is simply a memory manager for
c        cleastsq.
c 
c                     input parameters:
c 
c  a - the matrix to be decomposed (note that a is NOT damaged
c       by this subroutine in any way)
c  n, m - the dimensions of a
c  eps - the accuracy to which the decomposition (1) will be computed.
c 
c                     output parameters:
c 
c  w - the array containing the decomposition to be used by ncleastsq2.
c         Must be at least n*m*8 + 500 real *8 elements long
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c 
c 
c        . . . allocate memory for the decomposition of the matrix
c 
        iu=20
        lu=(n*m+10)*2
c 
        iw=iu+lu
        lw=(n*m+10)*2
c 
        it=iw+lw
        lt=(n*m+10)*2
c 
        iv=it+lt
        lv=(n*m+10)*2
c 
c       decompose the matrix
c 
        call cleastsq(a,w(iu),w(iw),w(it),n,m,ncols,
     1      rnorms,eps,w(iv))
c 
c       store various data at the beginning of the array w
c 
        w(1)=n+0.1
        w(2)=m+0.1
        w(3)=ncols+0.1
        w(4)=eps
        w(5)=iu+0.1
        w(6)=iw+0.1
        w(7)=it+0.1
        w(8)=iv+0.1
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine ncleasts2(w,rhs,x)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1)
        complex *16 x(1),rhs(1)
c 
c         This subroutine uses the QR-type (not exactly QR) decomposition
c         of a matrix to solve in the least squares sense the linear
c         system
c 
c                    A X = RHS                                             (1)
c 
c         The expansion used must have been produced by a prior call
c         to the subroutine ncleastsq (see), and is supplied to
c         this subroutine via the input array w.
c 
c                     input parameters:
c 
c  w - the array produced via a preceding call to the subroutine
c         ncleastsq (see); please note that the first  n*m*8 + 500
c         real *8 elements of this array should not be changed
c         between the call to this subroutine and the preceding call
c         to ncleastsq
c  rhs - the right-hand side in (1) (of length n)
c 
c                     output parameters:
c 
c  x - the solution (of length m) of the system (1) in the
c         least squares sense
c 
c 
c        . . . retrieve from the beginning of the array w all of the
c              relevant integer data
c 
        n=w(1)
        m=w(2)
        ncols=w(3)
        iu=w(5)
        iw=w(6)
        it=w(7)
        iv=w(8)
c 
c       solve the least squares problem
c 
        call cleasts2(w(iu),w(iw),w(it),n,m,ncols,rhs,
     1       x,w(iv) )
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cleastsq(a,u,w,t,n,m,ncols,rnorms,eps,v)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,m),u(n,1),v(m,1),w(m,1),t(1)
        dimension rnorms(1)
c 
c        this subroutine decomposes the input matrix a in the
c        form
c               A=U * T * W^*                                       (1)
c 
c        with the columns of each of the matrices U,W orthonormal,
c        and T a lower triangular matrix of dimensionality (ncols,ncols)
c 
c   NOTE: this subroutine uses the subroutine leastsq1 (see) to perform
c        almost all work. However, when m < n, it tends to be much
c        more efficient than leastsq1, since it performs the two
c        Gram-Schmidt processes in the optimal order (starting with
c        the longer gram-schmidt involving shorter vectors)
c 
c 
c                     input parameters:
c 
c  a - the matrix to be decomposed (note that a is NOT damaged
c       by this subroutine in any way)
c  n,m - the dimensionalities of a
c  eps - the accuracy to which the decomposition (1) will be computed.
c 
c                     output parameters:
c 
c  u - the matrix in (1). note that the matrix u has to be dimensioned
c        u(n,m), but on exit from this subroutine only the first
c        ncols columns of  u  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit u is dimensioned
c        u(n,ncols)
c  w - the third matrix in (1). it is dimensioned w(m,ncols) on
c        exit. on entry the parameter ncols is unknown,
c        so w should be dimensioned w(m,n) by the user. however,
c        on exit only the first ncols columns of w are meaningful
c  t - the second matrix in (1); On exit, it is structured as an
c        ncols * ncols matrix; however, on entry ncols is not
c        known, so it should be at least m*n complex *16 locations long.
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices u, w on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c 
c 
c                      work arrays:
c 
c  v - must be at least n*m real *8 locations long
c 
c       . . . if m > n , construct the decomposition of the matrix as
c             specified by the user
c 
        if(m .lt. n-2) goto 2000
c 
        call cleast1(a,u,w,t,n,m,ncols,rnorms,eps,v)
cccc        call prin2('after first leastsq1, u=*',u,n*ncols)
c 
        return
c 
 2000 continue
c 
c       n is greater than m. transpose the matrix, and decompose
c       the transpose
c 
        call cleastra(a,n,m,v)
        call cleascop(v,a,n*m)
c 
        call cleast1(a,w,u,t,m,n,ncols,rnorms,eps,v)
c 
cccc        call prin2('after leastsq1, u=*',u,n*ncols)
c 
c        transpose back everything that needs to be transposed back
c 
        call cleastra(a,m,n,v)
        call cleascop(v,a,n*m)
c 
        call cleastra(t,ncols,ncols,v)
        call cleascop(v,t,ncols**2)
c 
        call cleasrer(t,ncols,ncols,v)
        call cleascop(v,t,ncols**2)
c 
        call cleasrec(u,n,ncols,v)
        call cleascop(v,u,n*ncols)
c 
c 
c 
        call cleasrec(w,m,ncols,v)
        call cleascop(v,w,m*ncols)
c 
        call cleasrec(t,ncols,ncols,v)
        call cleascop(v,t,ncols**2)
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine cleasts2(u,w,t,n,m,ncols,y,x,work)
        implicit real *8 (a-h,o-z)
        save
        complex *16 u(n,ncols),w(m,ncols),t(ncols,ncols)
        complex *16 x(m),y(n),work(ncols),d
c 
c         this subroutine uses the QR-type (not exactly QR) decomposition
c         of a matrix to solve in the least squares sense the linear
c         system
c 
c                    A X = Y                                             (1)
c 
c         The expansion used must have been produced by a prior call
c         to the subroutine leastsq1 (see), and is of the form
c 
c               A=U * T * W^*                                            (2)
c 
c        with the columns of each of the matrices U,W orthonormal,
c        and T a lower triangular matrix of dimensionality (ncols,ncols)
c 
c                     input parameters:
c 
c  u - the matrix in (1). note that the matrix u has to be dimensioned
c        u(n,m), but on exit from this subroutine only the first
c        ncols columns of  u  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit u is dimensioned
c        u(n,ncols)
c  w - the third matrix in (1). it is dimensioned w(m,ncols) on
c        exit. on entry the parameter ncols is unknown,
c        so w should be dimensioned w(m,n) by the user. however,
c        on exit only the first ncols columns of w are meaningful
c  t - the second matrix in (1)
c  n,m - the dimensionalities of a
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices u, w on exit
c  y - the right-hand side in (1)
c 
c                     output parameters:
c 
c  x - the solution of the system (1) in the least squares sense
c 
c                      work arrays:
c 
c  work - must be at least ncols real *8 locations long
c 
c        . . . apply to the right-hand side the matrux  U^*
c 
        do 1400 i=1,ncols
        d=0
        do 1200 j=1,n
cccccc        d=d+u(j,i)*y(j)
        d=d+dconjg(u(j,i))*y(j)
 1200 continue
        work(i)=d
 1400 continue
c 
c       apply to the vector work the inverse of the matrix t
c 
        x(1)=work(1)/t(1,1)
        do 2000 i=2,ncols
c 
        d=0
        do 1600 j=1,i-1
        d=d+t(i,j)*x(j)
 1600 continue
c 
        x(i)=(work(i)-d)/t(i,i)
 2000 continue
c 
        do 2200 i=1,ncols
        work(i)=x(i)
 2200 continue
c 
c       multiply work by the matrix w
c 
        do 2600 i=1,m
        d=0
        do 2400 j=1,ncols
        d=d+w(i,j)*work(j)
 2400 continue
        x(i)=d
 2600 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine cleastra(a,n,m,b)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,m),b(m,n),x(1),y(1),z(1)
c 
c       transpose a
c 
        do 1400 i=1,n
        do 1200 j=1,m
cccc        b(j,i)=a(i,j)
        b(j,i)=dconjg(a(i,j))
 1200 continue
 1400 continue
        return
c 
c 
c 
c 
        entry cleascop(x,y,n)
c 
        do 2400 i=1,n
        y(i)=x(i)
 2400 continue
        return
c 
c 
c 
c 
        entry cleasubt(x,y,z,n)
c 
        do 3400 i=1,n
        z(i)=x(i)-y(i)
 3400 continue
        return
c 
c 
c 
c 
        entry cleaadd(x,y,z,n)
c 
        do 4400 i=1,n
        z(i)=x(i)+y(i)
 4400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine cleast1(a,u,w,t,n,m,ncols,rnorms,
     1    eps,v)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,m),u(n,1),v(m,1),w(m,1),t(1),rnorms(1)
c 
c        this subroutine decomposes the input matrix a in the
c        form
c               A=U * T * W^*                                       (1)
c 
c        with the columns of each of the matrices U,W orthonormal,
c        and T a lower triangular matrix of dimensionality (ncols,ncols)
c 
c 
c                     input parameters:
c 
c  a - the matrix to be decomposed (note that a is NOT damaged
c       by this subroutine in any way)
c  n,m - the dimensionalities of a
c  eps - the accuracy to which the decomposition (1) will be computed.
c 
c                     output parameters:
c 
c  u - the matrix in (1). note that the matrix u has to be dimensioned
c        u(n,m), but on exit from this subroutine only the first
c        ncols columns of  u  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit u is dimensioned
c        u(n,ncols)
c  w - the third matrix in (1). it is dimensioned w(m,ncols) on
c        exit. on entry the parameter ncols is unknown,
c        so w should be dimensioned w(m,n) by the user. however,
c        on exit only the first ncols columns of w are meaningful
c  t - the second matrix in (1)
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices u, w on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c 
c 
c                      work arrays:
c 
c  v - must be at least n*m real *8 locations long
c 
c        . . . using gram-schmidt process with pivoting, decompose
c              the matrix a in the form
c 
c          a=U  V^*,
c 
c        with  u an orthogonal matrix of minimum rank,
c        and v whatever it wants to be
c 
        ifpivot=1
        call cleaspiv(a,n,m,u,v,ncols,rnorms,eps,ifpivot)
c 
c        using gram-schmidt process without pivoting, decompose
c        the matrix v in the form
c 
c          v=w t^*,
c 
c        with  w an orthogonal matrix of minimum rank,
c        and t a triangular matrix of dimensionality ncols * ncols
c 
        ifpivot=0
        call cleaspiv(v,m,ncols,w,t,ncols2,rnorms,eps,ifpivot)
c 
        return
        end
c 
c 
c 
c 
c 
c 
        subroutine cleaspiv(a,n,m,b,v,ncols,rnorms,eps,ifpivot)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,m),b(n,m),v(m,1),prod
        dimension rnorms(1)
c 
c       this matrix applies the compressing gram-schmidt
c       process to the matrix a, obtaining its decomposition
c       in the form
c 
c             a=b v^T,                                        (1)
c 
c       with the matrices b, v having dimensionalities
c       b(n,ncols), v(m,ncols), respectively. the reason for the
c       existence of this subroutine is the hope that the
c       dimensionality ncols, determined by this subroutine,
c       is comsiderably lower than either n or m
c 
c                     input parameters:
c 
c  a - the matrix to be decomposed (note that a is NOT damaged
c       by this subroutine in any way)
c  n,m - the dimensionalities of a
c  eps - the accuracy to which the decomposition (1) will be computed.
c 
c                     output parameters:
c 
c  b - the matrix in (1). note that the matrix b has to be dimensioned
c        b(n,m), but on exit from this subroutine only the first
c        ncols columns of  b  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit b is dimensioned
c        b(n,ncols)
c  v - the second matrix in (1). it is dimensioned v(m,ncols) on
c        exit. on entry the parameter ncols is unknown,
c        so v should be dimensioned v(m,n) by the user. however,
c        on exit only the first ncols columns of v are meaningful
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices b, v on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
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
         call cleasgrm(b,n,m,rnorms,eps,ncols,ifpivot)
c 
c        project the original matrix on the obtained orthogonal
c        columns
c 
        do 2400 i=1,ncols
        do 2200 j=1,m
c 
        call cleascap(a(1,j),b(1,i),n,prod)
cccc        v(j,i)=prod
        v(j,i)=dconjg(prod)
 2200 continue
 2400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine cleasgrm(b,n,m,rnorms,eps,ncols,ifpivot)
        implicit real *8 (a-h,o-z)
        save
        complex *16 b(n,m),cd
        dimension rnorms(1)
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
        d=d+b(j,i)*dconjg(b(j,i))
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
        cd=b(j,i)
        b(j,i)=b(j,ipivot)
        b(j,ipivot)=cd
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
        call cleascap(b(1,i),b(1,j),n,cd)
c 
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*cd
 2770 continue
 2780 continue
 2790 continue
c 
c       normalize the i-th column
c 
        call cleascap(b(1,i),b(1,i),n,cd)
c 
c       if the norm of the remaining part of the matrix
c       is sufficiently small - terminate the process
c 
        d=cd
        if( (ifpivot .ne. 0) .and. (d .lt. thresh) ) return
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
        if(rnorms(j) .lt. thresh/100) goto 3200
c 
        call cleascap(b(1,i),b(1,j),n,cd)
        cd=dconjg(cd)
c 
         rrn=0
        do 3000 l=1,n
        b(l,j)=b(l,j)-b(l,i)*cd
        rrn=rrn+b(l,j)*dconjg(b(l,j))
 3000 continue
        rnorms(j)=dsqrt(rrn)
 3200 continue
 3400 continue
 4000 continue
  
 4200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cleascap(x,y,n,prod)
        implicit complex *16 (a-h,o-z)
        save
        dimension x(1),y(1)
c 
        prod=0
        do 1200 i=1,n
        prod=prod+x(i)*dconjg(y(i))
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine cleasrec(a,n,m,b)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,m),b(n,m)
c 
c       reverse the columns of a
c 
        do 1400 i=1,m
        do 1200 j=1,n
        b(j,i)=a(j,m-i+1)
 1200 continue
 1400 continue
        return
c 
c 
c 
c 
        entry cleasrer(a,n,m,b)
c 
c       reverse the rows of a
c 
        do 2400 i=1,m
        do 2200 j=1,n
        b(j,i)=a(n-j+1,i)
 2200 continue
 2400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine ncleamult(a,x,c,n,m,k)
        implicit complex *16 (a-h,o-z)
        save
        dimension a(n,m),x(m,k),c(n,k)
c 
        do 1600 i=1,n
        do 1400 j=1,k
c 
        cd=0
        do 1200 ll=1,m
c 
        cd=cd+a(i,ll)*x(ll,j)
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
        subroutine cleamat0(c,k,l,m,n,eps,delta,x,
     1      sa,ua,va,sb,ub,vb,work,lw,ltot,ncolsa,ncolsb)
        implicit real *8 (a-h,o-z)
c 
        save
        complex *16 x(l,m),c(k,n),work(1),
     1      ua(k,1),va(l,1),ub(m,1),vb(n,1),sa(1),sb(1)
c 
c        this subroutine solves in the least squares sense the
c        equation
c 
c                     A(k,l) * X(l,m) * B(m,n) = C(k,n),               (1)
c 
c       where A, B, C are user-specified matrices, and X is the
c       matrix to be found. Note that the dimensionalities of
c       the matrices A,B,C, X are as general as could be.
c       Actually, this subroutine uses as its input the SVDs of
c       the matrices A,B, provided by the calling subroutine cleamatr
c       (see), so that the equation (1) has the form
c 
c                    UA * SA * VA^* * X * UB * SB * VB^* = C       (2)
c 
c       . . . multiply the matrix c from the left by UA^* ;
c             note that the resulting matrix UAC is dimensioned
c             uac(ncolsa,n)
c 
        call cleamatl(ier,ua,k,ncolsa,c,k,n,work)
c 
c       . . . multiply the matrix UAC by VB from the right;
c             note that the resulting matrix uacvb is
c             dimensioned uacvb(ncolsa,ncolsb)
c 
        call cleamat(ier,work,ncolsa,n,vb,n,ncolsb,x)
c 
c        now, the equation (2) has assumed the form
c 
c                   SA * VA^* * X * UB * SB  = UACVB              (3)
c 
c       . . . multiply (3) by SA^{-1} from the left and
c             by SB^{-1} from the right in the appropriate least
c             squares sense
c 
         call cleasas(sa,ncolsa,sb,ncolsb,x,eps,delta)
c 
c        now, the equation (3) has assumed the form
c 
c                  VA^* * X * UB  = UACVB                     (4)
c 
c       . . . multiply (4) by VA from the left and
c             by UB^* from the right, obtaining the
c             solution X
c 
        call cleamat(ier,va,l,ncolsa,x,ncolsa,ncolsb,work)
c 
        call cleamar(ier,work,l,ncolsb,ub,m,ncolsb,x)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cleamat(ier,a,k,l,b,m,n,c)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(k,l),b(m,n),c(k,n),d
c 
c        this subroutine multiplies the matrix a by the matrix b,
c        getting the matrix c
c 
c        . . . check that the dimensionalities of the input
c              matrices agree
c 
        ier=4
        if(l .eq. m) ier=0
        if(ier .ne. 0) return
c 
c       multiply the matrix a by the matrix b getting c
c 
        do 2000 i=1,k
        do 1800 j=1,n
        d=0
        do 1600 kk=1,l
        d=d+a(i,kk)*b(kk,j)
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
        subroutine cleamatl(ier,a,k,l,b,m,n,c)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(k,l),b(m,n),c(l,n),d
c 
c        this subroutine multiplies the matrix a^* by the matrix b,
c        getting the matrix c
c 
c        check that the dimensionalities of the input matrices agree
c 
        ier=4
        if(k .eq. m) ier=0
        if(ier .ne. 0) return
c 
c       multiply the adjoint of the matrix a by the matrix b getting c
c 
        do 2000 i=1,l
        do 1800 j=1,n
        d=0
        do 1600 kk=1,k
cccc        d=d+a(kk,i)*b(kk,j)
        d=d+dconjg(a(kk,i))*b(kk,j)
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
        subroutine cleamar(ier,a,k,l,b,m,n,c)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(k,l),b(m,n),c(k,m),d
c 
c        this subroutine multiplies the matrix a by the matrix b^*,
c        getting the matrix c
c 
c        check that the dimensionalities of the input matrices agree
c 
        ier=4
        if(l .eq. n) ier=0
        if(ier .ne. 0) return
c 
c       multiply the matrix a by the the adjoint of the matrix b getting c
c 
        do 2000 i=1,k
        do 1800 j=1,m
        d=0
        do 1600 kk=1,l
cccc        d=d+a(i,kk)*b(j,kk)
        d=d+a(i,kk)*dconjg(b(j,kk))
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
        subroutine cleasas(sa,ncolsa,sb,ncolsb,uacvb,eps,delta)
        implicit real *8 (a-h,o-z)
        save
        complex *16 sa(1),sb(1),uacvb(ncolsa,ncolsb)
c 
c       this subroutine multiplies the user-supplied matrix UACVB
c       by SA^{-1} from the left and by SB^{-1} from the right in
c       the appropriate least squares sense
c 
        do 1800 i=1,ncolsb
        do 1600 j=1,ncolsa
c 
        d=abs(sa(j)*sb(i))
        if(d .gt. eps) goto 1400
        uacvb(j,i)=0
        goto 1600
c 
 1400 continue
c 
        uacvb(j,i)=uacvb(j,i)/(delta+d)
 1600 continue
 1800 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine cleamem(x,y,n)
        implicit real *8 (a-h,o-z)
        save
        complex *16 x(1),y(1)
c 
        do 1200 i=1,n
        y(i)=x(i)
 1200 continue
        return
        end
C 
C 
C 
C 
        SUBROUTINE PRINI(IP1,IQ1)
        save
        CHARACTER *1 MES(1), AA(1)
        REAL *4 A(1)
        REAL *8 A2(1)
        REAL *8 A4(1)
        INTEGER *4 IA(1)
        INTEGER *2 IA2(1)
        IP=IP1
        IQ=IQ1
        RETURN
  
C 
C 
C 
C 
C 
        ENTRY PRIN(MES,A,N)
        CALL  MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1200)(A(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1200)(A(J),J=1,N)
 1200 FORMAT(6(2X,E11.5))
         RETURN
C 
C 
C 
C 
        ENTRY PRIN2(MES,A2,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1400)(A2(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1400)(A2(J),J=1,N)
 1400 FORMAT(6(2X,E11.5))
        RETURN
C 
C 
C 
C 
        ENTRY PRINQ(MES,A4,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1500)(A4(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1500)(A4(J),J=1,N)
 1500 FORMAT(6(2X,e11.5))
        RETURN
C 
C 
C 
C 
        ENTRY PRINF(MES,IA,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1600)(IA(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1600)(IA(J),J=1,N)
 1600 FORMAT(10(1X,I7))
        RETURN
C 
C 
C 
C 
        ENTRY PRINF2(MES,IA2,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1600)(IA2(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1600)(IA2(J),J=1,N)
        RETURN
C 
C 
C 
C 
        ENTRY PRINA(MES,AA,N)
        CALL MESSPR(MES,IP,IQ)
 2000 FORMAT(1X,80A1)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,2000)(AA(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,2000)(AA(J),J=1,N)
        RETURN
        END
c 
c 
c 
c 
c 
        SUBROUTINE MESSPR(MES,IP,IQ)
        save
        CHARACTER *1 MES(1),AST
        DATA AST/'*'/
C 
C         DETERMINE THE LENGTH OF THE MESSAGE
C 
        I=0
        DO 1400 I=1,10000
        IF(MES(I).EQ.AST) GOTO 1600
        I1=I
 1400 CONTINUE
 1600 CONTINUE
         IF ( (I1.NE.0) .AND. (IP.NE.0) )
     1     WRITE(IP,1800) (MES(I),I=1,I1)
         IF ( (I1.NE.0) .AND. (IQ.NE.0) )
     1     WRITE(IQ,1800) (MES(I),I=1,I1)
 1800 FORMAT(1X,80A1)
         RETURN
         END
C 
C 
C 
C 
C 
        SUBROUTINE ZTIME(I)
        save
        J=1
        J=7-I+J
CCCC    I=MRUN(J)
        RETURN
        END
c 
c 
c 
c 
c 
        subroutine msgmerge(a,b,c)
        save
        character *1 a(1),b(1),c(1),ast
        data ast/'*'/
c 
        do 1200 i=1,1000
c 
        if(a(i) .eq. ast) goto 1400
        c(i)=a(i)
        iadd=i
 1200 continue
c 
 1400 continue
c 
        do 1800 i=1,1000
c 
        c(iadd+i)=b(i)
        if(b(i) .eq. ast) return
 1800 continue
        return
        end
c 
c 
c 
c 
c 
c 
        subroutine fileflush(iw)
        implicit real *8 (a-h,o-z)
c 
        save
        close(iw)
        open(iw,status='old')
        do 1400 i=1,1000000
c 
        read(iw,1200,end=1600)
 1200 format(1a1)
 1400 continue
 1600 continue
c 
        return
        end
  
  
c 
c 
c 
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c          this is the end of the debugging output and the beginning
c          of the actual SVD code
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c 
c 
        subroutine csvdpiv(ier,a,n,m,u,v,s,ncols,eps,
     1      w,lw,ltot)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,m),u(n,1),v(m,1),s(1),w(1)
c 
c        this subroutine constructs the singular value decomposition
c        of the usert-specified matrix a, under the assumption
c        that the rank of a is considerably smaller than its
c        dimensionalities n,m. the output of this subroutine
c        consists of two matrices u,v, and a vector s such that
c 
c                 a = u d v^*                                   (1)
c 
c        with t the diagonal matrix with the vector s on the
c        diagonal, and u, v two matrices whose columns are
c        orthonormal. the reason for the existence of this
c        subroutine is the fact that on exit, the matrices u,v
c        have dimensionalities (n,ncols), (m,ncols), respectively,
c        with ncols the numerical rank of a.
c 
c                   input parameters:
c 
c  a - the matrix to be svd-decomposed; note that it is NOT damaged
c       by this subroutine in any way
c  n,m - dimensionalities of a
c  eps - the accuracy with which the calculation will be performed;
c        should not be less than the machine precision; recommended
c        value: 1.0d-16
c  lw - the number of real *8 words provided by the user in the work
c        array w. must be at least m+n+10 (otherwise, the subroutine
c        can not even estimate the actual memory requirements).
c 
c                   output parameters:
c 
c  ier - error return code
c      ier=0 means successful execution
c      ier=4 means that the amount of space lw allocated supplied
c            in the work array w is insufficient
c      ier=1024 means that the amount of space lw allocated supplied
c            in the work array w is grossly insufficient
c  u - the matrix in (1). note that the matrix u has to be dimensioned
c        u(n,m), but on exit from this subroutine only the first
c        ncols columns of  b  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit b is dimensioned
c        b(n,ncols)
c  v - the second matrix in (1). it is dimensioned v(m,ncols) on
c        exit. on entry the parameter ncols is unknown,
c        so v should be dimensioned v(m,n) by the user. however,
c        on exit only the first ncols columns of v are meaningful
c  s - the array of singular values of a; note that this is a
c        COMPLEX *16 array; imaginary parts of its elements are
c        (hopefully) equal to zero.
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices b, v on exit
c  ltot - the total amount of space (in terms of complex *16 words)
c        needed by the subroutine in the array w
c 
c                   work arrays:
c 
c  w - must be sufficiently long. 3*n*m+ 2*(min(m,n))**2 +2*n+2*m+100
c        are always sufficient. however, most of the time it is much
c        less than that, and depends on the rank of a. the actual
c        amount of space available in w is communicated to the
c        subroutine via the parameter lw (see above), and if this
c        amount is insufficient, the output error code ier is set
c        accordingly, and the memory requirements of the subroutine
c        are returned in the parameter ltot (see above).
c 
c        . . . using gram-schmidt process with pivoting, decompose
c              the matrix a in the form
c 
c                    a=U  V^*,
c 
c              with  u an orthogonal matrix of minimum rank,
c              and v whatever it wants to be
c 
        if( (lw .gt. m) .and. (lw. gt. n)) goto 1200
        ier=1024
        return
 1200 continue
c 
        ifpivot=1
        call csvdpiv0(a,n,m,u,v,ncols,w,eps,ifpivot)
cccc         call prin2('after first csvdpiv0, rnorms=*',w,ncols)
c 
c       allocate memory for the subsequent operations
c 
        iw=1
        liw=m*ncols+10
c 
        iu2=iw+liw
        lu2=ncols**2+10
c 
        iv2=iu2+lu2
        lv2=ncols**2+10
c 
        it=iv2+lv2
        lt=lv2
c 
        irnorms=it+lt
        lrnorms=n
        if(m .gt. n) lrnorms=m
        lrnorms=lrnorms+10
c 
        iww=irnorms+lrnorms
        lww=n
        if(m .gt. n) lww=m
        lww=lww+10
c 
        ltot=iww+lww
        if(ltot .le. lw) goto 1400
        ier=4
        return
 1400 continue
c 
c       conclude the construction of the singular value decomposition
c 
        iffirst=1
        call csvdpiv1(a,u,v,w(iw),w(it),n,m,ncols,
     1    w(irnorms),eps,w(iv2),w(iu2),w(iww),iffirst)
cccc         call prin2('after csvdpiv1, w(it)=*',w(it),ncols)
c 
c       copy the right singular matrix into the array v
c 
        jj=0
        do 1800 i=1,ncols
        do 1600 j=1,m
        jj=jj+1
        v(j,i)=w(jj)
 1600 continue
 1800 continue
c 
c       now, copy the singular values into the array s
c 
        do 2000 i=1,ncols
        s(i)=w(it+i-1)
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine csvdpiv1(a,u,v,w,t,n,m,ncols,rnorms,
     1    eps2,v2,u2,ww,iffirst)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,m),u(n,1),v(m,1),w(m,1),t(1),
     1      v2(1),u2(1),ww(1)
        dimension rnorms(1)
c 
c        using gram-schmidt process with pivoting, decompose
c        the matrix a in the form
c 
c          a=U  V^*,
c 
c        with  u an orthogonal matrix of minimum rank,
c        and v whatever it wants to be
c 
        if(iffirst .eq. 0) goto 1200
        ifpivot=1
        call csvdpiv0(a,n,m,u,v,ncols,rnorms,eps2,ifpivot)
cccc         call prin2('after first csvdpiv0, rnorms=*',rnorms,ncols)
 1200 continue
c 
c        using gram-schmidt process without pivoting, decompose
c        the matrix v in the form
c 
c          v=w t^*,
c 
c        with  w an orthogonal matrix of minimum rank,
c        and t a triangular matrix of dimensionality ncols * ncols
  
c 
        ifpivot=0
        call csvdpiv0(v,m,ncols,w,t,ncols2,rnorms,eps2,ifpivot)
cccc         call prin2('after second csvdpiv0, rnorms=*',rnorms,ncols)
c 
c        use jacobi rotations to sv-decompose the matrix t
c 
        ifuv=1
c 
cccc        subroutine d2mcudv(a,n,u,v,eps,ifuv,numrows)
  
cccc        call prin2('before d2mcudv, t=*',t,ncols*ncols*2)
c 
        call d2mcudv(t,ncols,u2,v2,eps2,ifuv,numrows)
cccc        call prin2('after d2mudv, t=*',t,ncols)
cccc        call prinf('after d2mudv, numrows=*',numrows,1)
c 
c       multiply the matrix u from the right by u2, and the
c       matrix w from the right by v2^*
c 
cccc         call prin2('after d2mcudv, t=*',t,n*n*2)
        call csvdfinm(u,u2,w,v2,n,m,ncols,ww)
        return
        end
c 
c 
c 
c 
c 
        subroutine csvdfinm(u,u2,w,v2,n,m,ncols,ww)
        implicit real *8 (a-h,o-z)
        save
        complex *16 u(n,ncols),w(m,ncols),u2(ncols,ncols),
     1      v2(ncols,ncols),ww(1),d
c 
c       multiply u from the right by u2
c 
        do 1600 i=1,n
        do 1400 j=1,ncols
c 
        d=0
        do 1200 k=1,ncols
        d=d+u(i,k)*u2(k,j)
 1200 continue
        ww(j)=d
 1400 continue
        do 1500 j=1,ncols
        u(i,j)=ww(j)
 1500 continue
 1600 continue
c 
c         multiply w from the right by the adjoint of v2
c 
        zero=0
        do 2600 i=1,m
        do 2400 j=1,ncols
c 
        d=zero
        do 2200 k=1,ncols
        d=d+w(i,k)*dconjg(v2(j,k))
 2200 continue
        ww(j)=d
 2400 continue
        do 2500 j=1,ncols
        w(i,j)=ww(j)
 2500 continue
 2600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine csvdpiv0(a,n,m,b,v,ncols,rnorms,eps,ifpivot)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,m),b(n,m),v(m,1),prod
        dimension rnorms(1)
c 
c       this matrix applies the compressing gram-schmidt
c       process to the matrix a, obtaining its decomposition
c       in the form
c 
c             a=b v^T,                                        (1)
c 
c       with the matrices b, v having dimensionalities
c       b(n,ncols), v(m,ncols), respectively. the reason for the
c       existence of this subroutine is the hope that the
c       dimensionality ncols, determined by this subroutine,
c       is comsiderably lower than either n or m
c 
c                     input parameters:
c 
c  a - the matrix to be decomposed (note that a is NOT damaged
c       by this subroutine in any way)
c  n,m - the dimensionalities of a
c  eps - the accuracy to which the decomposition (1) will be computed.
c 
c                     output parameters:
c 
c  b - the matrix in (1). note that the matrix b has to be dimensioned
c        b(n,m), but on exit from this subroutine only the first
c        ncols columns of  b  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit b is dimensioned
c        b(n,ncols)
c  v - the second matrix in (1). it is dimensioned v(m,ncols) on
c        exit. on entry the parameter ncols is unknown,
c        so v should be dimensioned v(m,n) by the user. however,
c        on exit only the first ncols columns of v are meaningful
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices b, v on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
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
         call csvdgrm(b,n,m,rnorms,eps,ncols,ifpivot)
c 
c        project the original matrix on the obtained orthogonal
c        columns
c 
        do 2400 i=1,ncols
        do 2200 j=1,m
c 
        call svdscapr(a(1,j),b(1,i),n,prod)
cccc        v(j,i)=prod
        v(j,i)=dconjg(prod)
 2200 continue
 2400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine csvdgrm(b,n,m,rnorms,eps,ncols,ifpivot)
        implicit real *8 (a-h,o-z)
        save
        complex *16 b(n,m),cd
        dimension rnorms(1)
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
cccc        d=d+b(j,i)**2
        d=d+b(j,i)*dconjg(b(j,i))
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
        cd=b(j,i)
        b(j,i)=b(j,ipivot)
        b(j,ipivot)=cd
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
        call svdscapr(b(1,i),b(1,j),n,cd)
c 
cccc        cd=dconjg(cd)
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*cd
 2770 continue
 2780 continue
 2790 continue
c 
c       normalize the i-th column
c 
        call svdscapr(b(1,i),b(1,i),n,cd)
c 
c       if the norm of the remaining part of the matrix
c       is sufficiently small - terminate the process
c 
        d=cd
        if( (ifpivot .eq. 1) .and.
     1      (d .lt. thresh) ) return
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
        call svdscapr(b(1,i),b(1,j),n,cd)
        cd=dconjg(cd)
c 
         rrn=0
        do 3000 l=1,n
        b(l,j)=b(l,j)-b(l,i)*cd
        rrn=rrn+b(l,j)*dconjg(b(l,j))
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
        subroutine svdscapr(x,y,n,prod)
        implicit complex *16 (a-h,o-z)
        save
        dimension x(1),y(1)
c 
        prod=0
        do 1200 i=1,n
        prod=prod+x(i)*dconjg(y(i))
 1200 continue
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
        subroutine d2mcudv(a,n,u,v,eps,ifuv,numrows)
        implicit complex *16 (a-h,o-z)
        save
        dimension a(n,n),u0(2,2),v0(2,2),u(n,n),v(n,n),
     1      b(2,2),rlam(2),u1(2,2),v1(2,2)
        real *8 eps,thresh
c 
c        this subroutine uses the classical jacobi iteration
c        to construct the singular value decomposition of a
c        general complex *16 matrix, of the form
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
c  u,v - orthogonal matrices in (1)
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
         if( (cdabs(a(i,j)) .lt. thresh) .and.
     1       (cdabs(a(j,i)) .lt. thresh) ) goto 1800
        ifany=1
        numrows=numrows+1
c 
c       construct the svd of the (i,j)-th 2 * 2 submatrix
c 
        b(1,1)=a(i,i)
        b(1,2)=a(i,j)
        b(2,1)=a(j,i)
        b(2,2)=a(j,j)
        call oneuv(ier,b,rlam,u0,v0)
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
        call d2mcpcon(u0,u1)
        call d2mcpcon(v0,v1)
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
       call d2mcsrt(a(1,1),n,u,v,a(1,2),ifuv)
        return
        end
c 
c 
c 
c 
c 
        subroutine d2mcsrt(s,n,u,v,iw,ifuv)
        implicit real *8 (a-h,o-z)
        save
        complex *16 u(n,n),v(n,n),s(1),com,cd
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
        subroutine oneuv(ier,a0,rlam,u,v)
        implicit real *8 (a-h,o-z)
        save
        complex *16 z,b(2,2),u(2,2),a0(2,2),
     1      uu(2,2),v(2,2),w1(2,2),w2(2,2),
     2      w3(2,2),cd,clam1,clam2,clam,x(2),y(2)
        dimension rlam(2)
c 
c       multiply a by its adjoint
c 
        call d2mcpcon(a0,w1)
        call d2mprodc(w1,a0,b)
c 
c       find the eigenvalues of a^*a
c 
        two=2
        four=4
        cd=(b(1,1)-b(2,2))**2+four*b(1,2)*b(2,1)
        cd=cdsqrt(cd)
        clam1= (b(1,1)+b(2,2)+cd )/two
        clam2= (b(1,1)+b(2,2)-cd )/two
c 
c        find the eigenvectors
c 
        clam=clam1
        dc2=clam2*dconjg(clam2)
        dc1=clam1*dconjg(clam1)
        if( dc2 .gt. dc1 ) clam=clam2
c 
        d11=(b(1,1)-clam)*dconjg((b(1,1)-clam))
        d22=(b(2,2)-clam)*dconjg((b(2,2)-clam))
        d12=b(1,2)*dconjg(b(2,2))
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
        d=x(1)*dconjg(x(1))+x(2)*dconjg(x(2))
        d=dsqrt(d)
        x(1)=x(1)/d
        x(2)=x(2)/d
c 
        call d2md2ort(x,y)
c 
        uu(1,1)=x(1)
        uu(2,1)=x(2)
        uu(1,2)=y(1)
        uu(2,2)=y(2)
c 
        call d2mcpcon(uu,w2)
c 
c        now, obtain v
c 
         call d2mprodc(a0,uu,v)
c 
         d1=v(1,1)*dconjg(v(1,1))+v(2,1)*dconjg(v(2,1))
         d2=v(1,2)*dconjg(v(1,2))+v(2,2)*dconjg(v(2,2))
c 
         if(d2 .gt. d1) goto 3000
c 
         d=dsqrt(d1)
         call d2md2ort(v(1,1),v(1,2))
         goto 3200
c 
 3000 continue
c 
         d=dsqrt(d2)
         call d2md2ort(v(1,2),v(1,1))
 3200 continue
c 
c       obtain the diagonal matrix
c 
         call d2mprodc(a0,uu,w1)
         call d2mcpcon(v,b)
         call d2mprodc(b,w1,w3)
c 
        call d2mcopy2(w2,u)
c 
c       finally, make sure that the diagonal elements are real
c 
        z=1
        dzero=0
        d=w3(1,1)*dconjg(w3(1,1))
        if(d .ne. dzero) z=w3(1,1)/cdabs(w3(1,1))
        rlam(1)=cdabs(w3(1,1))
        v(1,1)=v(1,1)*z
        v(2,1)=v(2,1)*z
c 
        z=1
        dzero=0
        d=w3(2,2)*dconjg(w3(2,2))
        if(d .ne. dzero) z=w3(2,2)/cdabs(w3(2,2))
        rlam(2)=cdabs(w3(2,2))
        v(1,2)=v(1,2)*z
        v(2,2)=v(2,2)*z
c 
c        invert u, v
c 
        call d2mcopy2(v,w3)
        call d2mcpcon(u,v)
        call d2mcpcon(w3,u)
        return
        end
c 
c 
c 
c 
c 
        subroutine d2mcopy2(a,b)
        implicit complex *16 (a-h,o-z)
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
        entry d2mcpcon(a,b)
c 
        b(1,1)=dconjg(a(1,1))
        b(2,2)=dconjg(a(2,2))
        b(2,1)=dconjg(a(1,2))
        b(1,2)=dconjg(a(2,1))
        return
        end
c 
c 
c 
c 
c 
        subroutine d2md2ort(x,y)
        implicit real *8 (a-h,o-z)
        save
        complex *16 x(2),y(2)
c 
c       normalize x
c 
        d1=x(1)*dconjg(x(1))
        d2=x(2)*dconjg(x(2))
        d=d1+d2
        d=dsqrt(d)
        x(1)=x(1)/d
        x(2)=x(2)/d
c 
c       make y into a unit vector orthogonal to x
c 
c       . . . when |x(1)| > |x(2)|
c 
        if(d1 .lt. d2) goto 1200
        y(2)=x(1)
        y(1)=-dconjg(x(2))*x(1)/dconjg(x(1))
        return
 1200 continue
c 
c       . . . when |x(2)| > |x(1)|
c 
        y(1)=x(2)
        y(2)=-dconjg(x(1))*x(2)/dconjg(x(2))
        return
        end
c 
c 
c 
c 
c 
        subroutine d2mprodc(a,b,c)
        implicit complex *16 (a-h,o-z)
        save
        dimension a(2,2),b(2,2),c(2,2)
c 
        c(1,1)=a(1,1)*b(1,1)+a(1,2)*b(2,1)
        c(1,2)=a(1,1)*b(1,2)+a(1,2)*b(2,2)
        c(2,1)=a(2,1)*b(1,1)+a(2,2)*b(2,1)
        c(2,2)=a(2,1)*b(1,2)+a(2,2)*b(2,2)
        return
        end
