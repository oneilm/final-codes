        implicit real *8 (a-h,o-z)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER k '
        READ *,k
        CALL PRINF('k=*',k,1 )
C 
        PRINT *, 'ENTER nsings '
        READ *,nsings
        CALL PRINF('nsings=*',nsings,1 )
  
  
cccc        call testex(k,nsings)
  
        call testexs(k,nsings)
  
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine testexs(k,nsings)
        implicit real *8 (a-h,o-z)
        save
        dimension par1(100),ab(20 000),x(10000),whts(10000),
     1      rlams(10 000),rints(10 000),
     2      xsout(1000),uu(100 000),vv(100 000)
c 
        complex *16 w(4 000 000),fs2(1000 000),prods2(1000 000),
     1      ftest(10000),der,coefs(10000),errs(10000),ftest2(10000),
     2      whtsout(1000),forminte(1000),work(1000)
  
        external fun1
c 
c      construct all parameters
c 
        call fun1ini(nsings)
  
        a=1.0d-6
        a=1.0d-2
        a=0
ccc        a=1.0d-10
        b=0.5
cccc        b=1-a
  
        eps=1.0d-13
        lab=10 000
c 
  
  
        nfuns=nsings
  
        lw=4 000 000
        ifab=1
        lrlams=10000
        igraph=0
  
        call calsvbld(ier,a,b,k,fun1,par1,par2,
     1      nfuns,eps,ab,lab,ifab,nn,x,whts,n,rlams,lrlams,
     2      w,lw,ncols,rints,lused,igraph)
  
  
        call prinf('after calsvbld, nn=*',nn,1)
        call prin2('after calsvbld, rlams=*',rlams,ncols*2)
        call prin2('after calsvbld, rints=*',rints,ncols*2)
c 
c       check the obtained singular functions
c 
        call allsvchk(ier,w,ncols,eps,x,whts,n,rlams,
     1      fs2,prods2,ab,nn,k)
c 
c 
c       find a collection of points on the interval [a,b]
c       such that the matrix of values of singular vectors at these
c       points is well-conditioned
c 
        ifuv=1
        ncols2=ncols-5
        ncols2=ncols
  
        call calspnts(ier,ncols2,x,whts,ab,k,nn,w,
     1      xsout,whtsout,rints,uu,vv,ifuv,rcond,w(4*n*ncols+10) )
c 
        call prin2('after calspnts, rcond=*',rcond,1)
        call prin2('after calspnts, xsout=*',xsout,ncols2)
        call prin2('after calspnts, whtsout=*',whtsout,ncols2)
cccc        call prin2('after calspnts, uu=*',uu,ncols2**2*2)
c 
c       now, test the obtained matrices uu, vv
c 
  
        i=10
        do 2600 j=1,ncols2
c 
        call cnesteva(ier,w,ab,nn,k,xsout(j),i,ftest(j),der)
 2600 continue
c 
        call calqrma4(vv,ftest,coefs,ncols2)
  
        call prin2('coefs=*',coefs,ncols2*2)
  
        call calqrma4(uu,coefs,ftest2,ncols2)
  
        do 2200 i=1,ncols2
c 
        errs(i)=ftest(i)-ftest2(i)
 2200 continue
  
        call prin2('and ftest=*',ftest,ncols2*2)
        call prin2('and ftest2=*',ftest2,ncols2*2)
        call prin2('and errs=*',errs,ncols2*2)
c 
c       test the interpolation routine
c 
        xx=xsout(16)
        xx=0.3
        call cnesinte(ab,nn,k,w,vv,ncols,xx,
     1      forminte,work)
  
        call prin2('after cnesinte, forminte=*',forminte,ncols*2)
        call prinf('while ncols=*',ncols,1)
  
  
        itest=100
        ntest=20
        call apprtest(xsout,ncols,vv,itest,ntest,a,b,
     1      ab,nn,k,w)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine apprtest(xsout,ncols,v,itest,ntest,a,b,
     1      ab,nn,k,coefs)
        implicit real *8 (a-h,o-z)
        save
        dimension xsout(1),xstest(10000),v(1),ab(2,1)
c 
        complex *16 fsout(1000),fstest(1000),forminte(1000),cd,
     1      errs(1000),work(1000)
c 
c       evaluate the test function at the support nodes
c 
        do 1100 i=1,ncols
c 
        call fun1(xsout(i),itest,par1,par2,fsout(i))
 1100 continue
c 
        call prin2('in apprtest, xsout=*',xsout,ncols)
        call prin2('in apprtest, fsout=*',fsout,ncols*2)
c 
c       create the test nodes and the test function at them
c 
        h=(b-a)/ntest
        do 1200 i=1,ntest
c 
        xstest(i)=a+(i-1)*h+h/2
        call fun1(xstest(i),itest,par1,par2,fstest(i))
 1200 continue
c 
        call prin2('in apprtest, xstest=*',xstest,ntest)
        call prin2('in apprtest, fstest=*',fstest,ntest*2)
c 
c       for each of the test points, interpolate the function
c       from the support nodes to the test point
c 
        do 1600 i=1,ntest
c 
        call cnesinte(ab,nn,k,coefs,v,ncols,xstest(i),
     1      forminte,work)
c 
cccc        call prinf('i=*',i,1)
cccc        call prin2('and forminte=*',forminte,ncols*2)
c 
        cd=0
        do 1400 j=1,ncols
c 
        cd=cd+forminte(j)*fsout(j)
 1400 continue
c 
cccc        call prin2('and cd=*',cd,2)
cccc        call prin2('while fstest(i)=*',fstest(i),2)
cccc        call prin2('and fstest(i)-cd=*',fstest(i)-cd,2)
c 
        errs(i)=fstest(i)-cd
 1600 continue
  
        call prin2('and errs=*',errs,ntest*2)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine fun1(x,i,par1,par2,f)
        implicit real *8 (a-h,o-z)
        save
        dimension par1(1),par2(1),sings(10000)
c 
        complex *16 f,ima
c 
        data ima/(0.0d0,1.0d0)/
c 
        f=x**sings(i)+ima*(1-x)**sings(i)
ccccc        f=x**sings(i)
cccc        f=x**(sings(i)+1)+ima*(1-x)**sings(i)
  
        return
c 
c 
c 
        entry fun1ini(nsings7)
c 
        nsings=nsings7
        done=1
        h=10.3*done/nsings
        do 2200 j=1,nsings
        sings(j)=j*h
 2200 continue
c 
        call prin2('in fun1ini, sings=*',sings,nsings)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine allsvchk(ier,coefs,ncols,eps,x,whts,n,rlams,
     1      fs,prods,ab,nn,k)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),whts(1),rlams(1),ab(2)
        complex *16 fs(n,ncols),prods(ncols,ncols),cd,f,cds(1000)
c 
c       evaluate all singular functions at all nodes in array x
c 
        do 1400 i=1,ncols
c 
        do 1200 j=1,n
c 
        call cnestev2(ier,coefs,ab,nn,k,x(j),i,fs(j,i))
 1200 continue
 1400 continue
  
cccc         call prin2('in allsvchk, fs=*',fs,n*ncols*2)
cccc         call prin2('in allsvchk, whts=*',whts,n)
  
cccc           stop
  
  
c 
c       evaluate all inner products of all functions
c 
        do 2400 i=1,ncols
        do 2200 j=1,ncols
c 
        cd=0
        do 1800 ll=1,n
c 
        cd=cd+fs(ll,i)*conjg(fs(ll,j))*whts(ll)
 1800 continue
  
cccc         call prinf('in allsvchk, i=*',i,1)
cccc         call prinf('in allsvchk, j=*',j,1)
cccc         call prin2('in allsvchk, cd=*',cd,2)
c 
         prods(j,i)=cd
 2200 continue
 2400 continue
c 
cccc        call prin2('prods in allsvchk are*',prods,ncols**2*2)
  
  
        return
        end
c 
c 
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c        this is the end of the debugging code and the beginning
c        of the actual SVD routine
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c        This file contains a collection of subroutines for the
c        construction of singular value decompositions (or rather,
c        of left singular vectors) of families of complex functions,
c        of a real variable, and for the performance of certain
c        associated operations. The assumption is that the number of
c        functions to be compressed is too great for them to be
c        discretized and evaluated simultaneously; the user specifies
c        them things via subroutines that evaluate them. The file
c        contains 3 user-accesible subroutines: calsvbld, calspnts,
c        cnesinte. Following is a brief description of these subroutines.
c 
c   calsvbld - constructs the left singular functions of the set of
c       user-specified functions on the interval [a,b], and returns
c       to the user nested Legendre expansions for the first ncols
c       singular functions; all elements with higher number than
c       ncols have weights associated with them that are less than
c       eps (user-specified). The subroutine also returns to the
c       user the singular values associated with the said singular
c       functions. Once this subroutine has been called, any one of
c       the obtained singular functions (together with its derivative)
c       can be evaluated by a call to the subroutines cnesteva,
c       cnestev2, cnestem (see).
c 
c   calspnts - uses the output of the subroutine calsvbld to construct
c       the points xsout on the interval [a,b] such that the matrix of
c       values of left singular vectors (obtained by calsvbld) at the
c       points xsout is well-conditioned. Depending on the value of the
c       parameter ifuv, the subroutine will also produce the said matrix,
c       its inverse, and the quadrature weights whtsout corresponding
c       to the nodes xsout
c 
c   cnesinte - uses the output of the subroutines calsvbld,
c       calspnts (or of calsvcmb instead of calsvbld) to construct
c       in interpolation formula connecting the values of a function
c       at the nodes xsout with the value at the user-specified
c       point x on the interval [a,b].
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine calspnts(ier,ncols,xs,whts,ab,k,nn,coefs,
     1      xsout,whtsout,rints,uu,vv,ifuv,rcond,w)
        implicit real *8 (a-h,o-z)
        save
        dimension whts(1),xs(1),w(1),xsout(1),
     1      ab(1),coefs(1)
c 
        complex *16 uu(ncols,ncols),vv(ncols,ncols),whtsout(1),
     1      rints(1)
c 
c       This subroutine is a complex version of the subroutine allspnts
c       (see). It uses the output of the subroutine calsvbld
c       (or of calsvcmb) to construct the points xsout on the interval
c       [a,b] such that the matrix of values of left singular vectors
c       (obtained by allsvbld or allsvcmb) at the points xsout is
c       well-conditioned. Depending on the value of the parameter ifuv,
c       the subroutine will also produce the said matrix and its inverse.
c       It also returns the array whtsout of weights corresponding to
c       the nodes xsout, such that all of the functions supplied by the
c       user are integrated with accuracy eps.
c     PLEASE NOTE THAT THIS SUBROUTINE CAN ONLY BE USED AFTER ONE OF
c     SUBROUTINES CALSVBLD, CALSVCMB HAS BEEN CALLED. IT HAS NO USES
c     WHATSOEVER AS A STAND-ALONE DEVICE!!!!!!!!!!!!!!!!
c 
c                       Input parameters:
c 
c  ncols - the number of left singular vectors to be used (also the
c       number of nodes to construct). Must be less than or equal to
c       the number ncols returned by the preceding call to allsvbld or
c       allsvcmb.
c  xs - the discretization of the interval [a,b] consisting of nn
c        subintervals with k Legendre nodes on each. Must have been
c       produced by a preceding call to allsvbld or allsvcmb.
c  whts - the weights associated with the discretization x
c  ab - the array of ends of subintervals created. For each i=1,2,...,nn,
c       ab(1,i) is the left end of the i-th subinterval; ab(2,i) is
c       the right end of the i-th subinterval. Must have been produced
c       by a preceding call to allsvbld or allsvcmb.
c  k - the number of Gaussian nodes on each subinterval; must be the same
c       as the one given earlier to allsvbld or allsvcmb
c  nn - the number of subintervals in the nested discretization of
c       the interval [a,b]. Must have been produced by a preceding call
c       to allsvbld or allsvcmb.
c  coefs - contains the coefficients (nn*k*ncols of them) of
c        the nested Legendre expansions of the singular functions. Must
c        have been produced by a preceding call to calsvbld or calsvcmb.
c  ifuv - tells the subroutine whether to construct the matrix uu of values
c        of the singular function at the selected points xsout, and the
c        inverse vv of the matrix uu.
c 
c                        Output parameters:
c 
c  xsout - a collection of ncols points on the interval [a,b] such that
c        the matrix of values of the first ncols singular vectors at these
c        points is well-conditioned
c  whtsout - the quadrature weights corresponding to the nodes xsout
c  uu - the  matrix of values of the first ncols singular vectors at the
c        points xsout
c  vv - the  inverse of the matrix uu
c  rcond - an estimate of the condition number of the matrix uu; due to
c        poor scaling, often understates the accuracy.
c 
c                  Work arrays:
c 
c  w - should contain at least nn*k*ncols+nn*k*2 +100 real *8 elements
c 
c 
c        . . . allocate memory
c 
        n=nn*k
c 
        iuuadj=1
        luuadj=n*ncols+10
c 
        luuadj=luuadj*2
c 
        iipivots=iuuadj+luuadj
        lipivots=n+2
c 
        irnorms=iipivots+lipivots
        lrnorms=n+2
        lrnorms=lrnorms*2
c 
c        find a collection of points such that the matrix of values of
c        singular vectors at these points is well-conditioned
c 
        call calspnt0(ier,uu,n,ncols,w(iuuadj),xs,whts,
     1      xsout,whtsout,rints,w(iipivots),w(irnorms),coefs,rcond,k,
     2      ab,nn,vv,ifuv)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine calsvbld(ier,a,b,k,fun,par1,par2,
     1      nfuns,eps,ab,lab,ifab,nn,x,whts,n,rlams,lrlams,
     2      w,lw,ncols,rints,lused,igraph)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),rlams(1),w(1),ab(2,1),par1(1),par2(1),
     1      whts(1),rints(1)
c 
c       This subroutine is a complex version of the subroutine
c       allsvbld. It constructs the left singular functions
c       of the set of of user-specified functions on the interval
c       [a,b], and returns to the user nested Legendre expansions
c       for the first ncols singular functions; all elements with
c       higher number than ncols have weights associated with them
c       that are less than eps (user-specified). The subroutine also
c       returns to the user the singular values associated with
c       the said singular functions. Once this subroutine has
c       been called, any one of the obtained singular functions
c       (together with its derivative) can be evaluated by a call
c       to one of the subroutines cnesteva, cnestev2, cnestem (see).
c 
c       To this effect, the subroutine finds a subdivision of the
c       interval [a,b], such that on every subinterval, every one
c       of the user-provided collection of functions (given by the
c       subroutine fun (see description below) is interpolated to
c       precision eps by a k-point Gaussian quadrature.
c 
c       Then, the subroutine constructs such a nested Legendre
c       discretization of the interval [a,b] with the corresponding
c       set of weights, and applies the Gram-Schmidt process to the
c       resulting finite-dimensional object.
c 
c            Input parameters:
c 
c  a,b - the ends of the interval to be subdivided
c  k - the number of Gaussian nodes on each subinterval
c  fun - the user-supplied function to be integrated
c      The calling sequence of fun must be:
c 
c      fun(x,i,par1,par2,f),
c 
c       where x is the input argument where the function is to be
c       evaluated, i the sequence number of function, f (output)
c       is the value of the i-th function at the point x, and par
c       is an input parameter, carrying whatever information fun
c       might need.
c 
c  par1,par2 - whatever the user-supplied subroutine fun needs
c  nfuns - the number of functions to be compressed
c  eps - the accuracy to which the legendre expansion must represent
c      the user-supplied function for the subroutine to decide that
c      the representation is accurate
c  ab - the array of ends of subintervals created by a preceding
c      call to this subroutine.  For each i=1,2,...,nn,
c      ab(1,i) is the left end of the i-th subinterval; ab(2,i) is
c      the right end of the i-th subinterval. Note that this is an
c      input parameter only if the input parameter ifab (see below)
c      had been set to 0 by the user.
c  lab - the length of the user-provided array ab (see below)
c  ifab - the parameter telling the subroutine whether it should
c      create the subdivision ab of the interval [a,b] and
c      attendant parameters (x,whts,nn), or to view them as input
c      parameters. ifab=1 means that the subroutine will create
c      the subdivision; ifab=0 means that it will view the parameters
c      ab,nn,x,whts as input parameters.
c  nn - the number of subintervals in the nested discretization of
c      the interval [a,b], created by a preceding call to this
c      subroutine.  Note that this is an input parameter only
c      if the input parameter ifab (see above) had been set to 0 by
c      the user.
c  lrlams - the amount of space (in real *8 words) in the array rlams
c      in which the subroutine will return the singular values of the
c      matrix. This parameter is used by the subroutine to bomb if the
c      number of singular values that are greater than eps turns out
c      to be greater than the length of array rlams
c 
c  x - the discretization of the interval [a,b] consisting of nn
c      subintervals with k Legendre nodes on each,  created by a
c      preceding call to this subroutine.  Note that this is an
c      input parameter only if the input parameter ifab (see above)
c      had been set to 0 by the user.
c  whts - the weights associated with the discretization x, created
c      by a preceding call to this subroutine.   Note that this is
c      an in parameter only if the input parameter ifab (see
c      above) had been set to 0 by the user.
c  lw - the length of the user-provided array w. It has to be
c      large. If it is insufficient, ier os set to 64, and the
c      execution of the subroutine terminated
c  igraph - the integer parameter telling the subroutine whether
c      it should plot the obtained singular functions, and which
c      files to plot them on. For example, if igraph is set to
c      100, the first singular function will be plotted on the
c      file gn101, the  second singular function will be plotted
c      on the file gn102, the  third singular function will be
c      plotted on the file gn103, etc. Setting igraph .le. 0
c      will cause plotting to be suppressed. Incidentally, the
c      files gn101, gn102, ... are GNUPLOT-readable.
c 
c                       Output parameters:
c 
c  ier -error return code;
c          ier=0 means normal return
c          ier=4 means that the amount of space in the array ab
c                 (given by the parameter lab) is insufficient. very
c                 often, the real cause of the problem is eps that
c                 is too small.
c          ier=64 means that the amount of space in the user-provided
c                 array w (given by the parameter lw) is insufficient.
c  ab - the array of ends of subintervals created. For each i=1,2,...,nn,
c     ab(1,i) is the left end of the i-th subinterval; ab(2,i) is
c     the right end of the i-th subinterval. Note that this is an
c     output parameter only if the input parameter ifab (see above)
c     had been set to 1 by the user.
c  nn - the number of subintervals in the nested discretization of
c      the interval [a,b].  Note that this is an output parameter only
c      if the input parameter ifab (see above) had been set to 1 by
c      the user.
c  w - the first n*ncols elements of array w contain the coefficients of
c        the nested Legendre expansions of the singular functions;
c 
c        IMPORTANT NOTE: the first n*ncols elements of array w should
c                        not be changed between the call to this
c                        subroutine and subsequent calls to the
c                        subroutine nesteva.
c 
c  x - the discretization of the interval [a,b] consisting of nn
c        subintervals with k Legendre nodes on each.  Note that this is an
c        output parameter only if the input parameter ifab (see above)
c        had been set to 1 by the user.
c  whts - the weights associated with the discretization x.  Note that
c     this is an output parameter only if the input parameter ifab
c     (see above) had been set to 1 by the user.
c  n - the number of nodes in the discretization x of the interval [a,b];
c        equal to nn*k
c  rlams - the singular values of the matrix of functions given by the
c        function fun. Note that though only ncols elements of array rlams
c        are meaningful, nfuns elements have to be allocated by the
c        user.
c  ncols - the number of singular functions in the singular value
c      decomposition of the
c      nfuns functions given by the subroutine fun (see above) whose
c      corresponding weights are greater than eps (see above). It is
c      also the number of meaningful columns of arrays fs, coefs on
c      exit from this subroutine.
c  rints - the array of length ncols containing the integrals of the
c      first ncols singular functions of the operator
c  lused - the largest piece of the array w used by the subroutine at
c      any time during the execution
c 
c                  Work arrays:
c 
c  w - treated by the subroutine as both output parameter and a work
c      array. The subroutine utilizes the first lused elements of
c      w as a work space; after that, it returns in the first n*ncols
c      elements of w the coefficients of the nested legendre expansions
c      first ncols singular functions of the operator
c 
c        . . . construct the nested subdivision of the interval [a,b],
c              such that on each subinterval, the pair-wise products
c              of all user-provided functions are integrated accurately
c              with k nodes
c 
         ier=0
c 
         if(ifab .eq. 0) goto 1200
c 
         call cmaninco(jer,a,b,k,fun,nfuns,par1,par2,
     1    eps,ab,lab,nn,w,lw)
c 
         if(jer .ne. 0) ier=16
         if(ier .ne. 0) return
c 
 1200 continue
c 
c       . . . construct the all the nested Legendre nodes
c 
        ier=0
        call cnodewht(ab,nn,k,x,whts)
c 
c       construct the QR decompositions of the matrix of values of
c       the functions to be compressed at the nodes of the adaptive
c       discretization; also the nested Legendre expansions of the
c       first ncols such functions
c 
        n=nn*k
c 
c       . . . allocate memory
c 
        iu=1
        lu=k**2+2
c 
        iv=iu+lu
        lv=k**2+2
c 
        iww=iv+lv
        lww=n+2
        lww=lww*2
c 
        iww2=iww+lww
        lww2=n+2
c 
        icoefs=iww2+lww2
        lencoefs=lw-icoefs-10
c 
        lused=iww+lww
c 
        call calffuns(jer,x,whts,n,nfuns,fun,par1,par2,
     1      eps,w(icoefs),lencoefs,ncols,rlams,nn,k,w(iu),w(iv),
     2      w(iww),rints,lused2,lrlams,igraph,w(iww),w(iww2) )
c 
        lused=lused+lused2
c 
        if(jer .ne. 0) ier=8
        if(jer .ne. 0) return
c 
        call cmvecmov(w(icoefs),w,ncols*n*2)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cnesinte(ab,nn,k,coefs,vv,ncols,x,
     1      forminte,vals)
        implicit real *8 (a-h,o-z)
        save
        dimension ab(2,1),coefs(1)
c 
        complex *16 forminte(1),vals(1),vv(ncols,ncols)
c 
c       This subroutine uses the output of the subroutines calsvbld,
c       calspnts (or of calsvcmb instead of calsvbld) to construct
c       in interpolation formula connecting the values of a function
c       at the nodes xsout with the value at the user-specified
c       point x on the interval [a,b].
c 
c                          Input parameters:
c 
c  ab - the array of ends of subintervals created. For each i=1,2,...,nn,
c       ab(1,i) is the left end of the i-th subinterval; ab(2,i) is
c       the right end of the i-th subinterval. Must have been produced
c       by a preceding call to allsvbld or allsvcmb.
c  nn - the number of subintervals in the nested discretization of
c       the interval [a,b]. Must have been produced by a preceding call
c       to allsvbld or allsvcmb.
c  k - the number of Gaussian nodes on each subinterval; must be the same
c       as the one given earlier to allsvbld or allsvcmb
c  coefs - contains the coefficients (nn*k*ncols of them) of the nested
c       Legendre expansions of the singular functions. Must have been
c       produced by a preceding call to calsvbld or calsvcmb.
c  vv - the  inverse of the matrix uu of values of the first ncols
c        singular vectors at the points xsout (produced by subroutine
c        calsvbld or calsvcmb)
c  x - the point on the interval [a,b] to which the function is to be
c        interpolated
c 
c                          Output parameters:
c 
c  forminte - the coefficients of the linear form connection the value
c        of a function at the point x with its values at the points xsout
c 
c 
c                          Work arrays:
c 
c  vals - must be at least ncols complex *16 locations long
c 
c 
c       . . . evaluate all singular functions at the point x
c 
        call cnestem(ier,coefs,ab,nn,k,ncols,x,vals)
c 
c       apply the matrix vv to the vector vals, obtaining the
c       interpolation coefficients
c 
        call calqrma5(vv,vals,forminte,ncols)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine calspnt0(ier,uu,n,ncols,uuadj,xs,whts,
     1      xsout,whtsout,rints,ipivots,rnorms,w,rcond,
     2      k,ab,nn,vv,ifuv)
        implicit real *8 (a-h,o-z)
        save
        dimension ipivots(1),whts(1),xs(1),w(1),xsout(1),
     1      ab(1)
c 
        complex *16 uuadj(ncols,n),rnorms(1),uu(ncols,ncols),
     1      vv(ncols,ncols),der,rints(1),whtsout(1)
c 
c       evaluate all of the singular vectors at the nodes of
c       the discretization xs of the interval of definition
c 
        do 1450 j=1,n
c 
        call cnestem(ier,w,ab,nn,k,ncols,xs(j),uuadj(1,j) )
c 
        d=sqrt(whts(j))
        do 1430 i=1,ncols
c 
        uuadj(i,j)=uuadj(i,j)*d
 1430 continue
 1450 continue
  
cccc        call prin2('and whts=*',whts,n)
  
        eps2=1.0d-30
c 
        call calsppiv(uuadj,ncols,n,rnorms,eps2,ncols3,ipivots)
c 
c       sort the array ipivots
c 
        call calsbubb(ipivots,ncols3)
c 
c       pick out of the array xs the points at which the matrix
c       of values of the singular vectors is well-conditioned
c 
        do 1500 i=1,ncols
c 
        ii=ipivots(i)
        xsout(i)=xs(ii)
 1500 continue
c 
        if(ifuv .eq. 0) return
c 
c       construct the matrix of values of singular vectors at the
c       selected nodes
c 
        do 1800 i=1,ncols
c 
        do 1600 j=1,ncols
c 
        call cnesteva(ier,w,ab,nn,k,xsout(j),i,uu(j,i),der)
c 
        vv(j,i)=uu(j,i)
 1600 continue
 1800 continue
c 
c       invert the matrix uu
c 
        call corthom(vv,ncols,uuadj,rcond)
c 
c       evaluate the weights corresponding to the nodes xsout
c 
        call calqrma5(vv,rints,whtsout,ncols)
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine calsbubb(ipivots,ncols)
        save
        dimension ipivots(1)
c 
        do 1400 i=1,ncols
        do 1200 j=1,ncols-1
c 
        if(ipivots(j) .lt. ipivots(j+1)) goto 1200
c 
        k=ipivots(j)
        ipivots(j)=ipivots(j+1)
        ipivots(j+1)=k
 1200 continue
 1400 continue
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine calffuns(ier,x,whts,n,nfuns,fun,par1,par2,
     1      eps,fs,lfs,ncols,sss,nn,k,u,v,ww,rints,lused,lrlams,
     2      igraph,cww,ww2)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),whts(1),ww(1),par1(1),par2(1),
     1      u(1),v(1),xlege(100),whtslege(100),
     2      sss(1),ww2(1)
c 
        complex *16 fs(n,1),cww(1),ima,rints(1)
c 
        character *100 adummy,mes(100)
c 
        external calsvmat
c 
        data ima/(0.0d0,1.0d0)/
c 
c        . . . constuct the square roots of weights
c 
        ier=0
c 
        do 1100 i=1,n
        ww(i)=sqrt(whts(i))
 1100 continue
c 
c       conduct the SVD on the constructed functions
c 
        call cmvecsvd(jer,calsvmat,fun,par1,par2,x,ww,n,
     1    nfuns,eps,ncols,sss,lrlams,fs,lfs,lused)
c 
        if(jer .ne. 0) ier=8
        if(ier .ne. 0) return
c 
        call prin2('in alffuns after mvecsvd, sss=*',sss,ncols*2)
c 
c       compensate the vectors in fs for the multiplication by the
c       square roots of weights; at the same time, evaluate the
c       integrals of the QR-functions
c 
        do 1800 i=1,ncols
c 
        rints(i)=0
        do 1600 j=1,n
c 
        fs(j,i)=fs(j,i)/ww(j)
c 
        rints(i)=rints(i)+fs(j,i)*whts(j)
 1600 continue
 1800 continue
c 
        call prin2('integrals of fs evaluated in calffuns*',
     1      rints,ncols)
c 
c        if the user so requested, plot the obtained singular vectors
c 
        if(igraph .le. 0) goto 2100
c 
        do 2000 i=1,ncols
c 
        iw=igraph+i
        do 1850 j=1,n
c 
        ww(j)=fs(j,i)
        ww2(j)=-fs(j,i)*ima
 1850 continue
c 
 1860 format(i3,'*')
 1870 format(100a1)
c 
        write(adummy,1860) i
        call msgmerge(
     1      'real part of i-th singular function with i=*',
     2      adummy,mes)
c 
        call lotagraph2(iw,x,ww,n,x,ww2,n,mes)
c 
 2000 continue
  
 2100 continue
c 
c        construct Legendre expansions of all QR-functions on all
c        elementary intervals
c 
        itype=2
        call legeexps(itype,k,xlege,u,v,whtslege)
c 
        do 2400 i=1,ncols
        do 2200 j=1,nn
c 
        jj=(j-1)*k+1
        call calqrmat(u,fs(jj,i),cww(jj),k)
c 
 2200 continue
c 
        do 2300 j=1,n
        fs(j,i)=cww(j)
 2300 continue
c 
 2400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine calsppiv(b,n,m,rnorms,eps,ncols,ipivots)
        implicit real *8 (a-h,o-z)
        save
        dimension rnorms(1),ipivots(1)
        complex *16 b(n,m),cd
c 
c       This subroutine applies a pivoted double gram-schmidt
c       procedure to the matrix b, resulting in a collection of
c       orthogonal vectors. The number of these vectors is
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
        d=0
        do 1200 j=1,n
        d=d+b(j,i)*conjg(b(j,i))
 1200 continue
        rnorms(i)=d
        dtot=dtot+d
 1400 continue
c 
c       . . . conduct gram-schmidt iterations
c 
        thresh=dtot*eps**2
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
        d=b(j,i)
        b(j,i)=b(j,ipivot)
        b(j,ipivot)=d
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
        call cmvesca(b(1,i),b(1,j),n,cd)
c 
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*cd
 2770 continue
 2780 continue
 2790 continue
c 
c       normalize the i-th column
c 
        call cmvesca(b(1,i),b(1,i),n,cd)
c 
c       if the norm of the remaining part of the matrix
c       is sufficiently small - terminate the process
c 
        d=cd
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
        call cmvesca(b(1,i),b(1,j),n,cd)
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
        subroutine calqrmat(a,x,y,n)
        implicit complex *16 (a-h,o-z)
        save
        dimension x(1),y(1)
        real *8 a(n,n)
c 
        do 2000 i=1,n
        d=0
        do 1400 j=1,n
        d=d+a(i,j)*x(j)
 1400 continue
        y(i)=d
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine calqrma5(a,x,y,n)
        implicit complex *16 (a-h,o-z)
        save
        complex *16 a(n,n),x(1),y(1)
c 
        do 2000 i=1,n
        d=0
        do 1400 j=1,n
        d=d+a(j,i)*x(j)
 1400 continue
        y(i)=d
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine calqrma4(a,x,y,n)
        implicit complex *16 (a-h,o-z)
        save
        complex *16 a(n,n),x(1),y(1)
c 
        do 2000 i=1,n
        d=0
        do 1400 j=1,n
        d=d+a(i,j)*x(j)
 1400 continue
        y(i)=d
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine calsvmat(i,j,f,fun,par1,par2,xs,ww)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),ww(1)
        complex *16 f
c 
        call fun(xs(i),j,par1,par2,f)
c 
        f=f*ww(i)
        return
        end
