        implicit real *8 (a-h,o-z)
        dimension par1(100),ab(20 000),
     1      x(10 000),whts(10 000),rlams(10 000),
     2      rints(10 000),xtest(10 000),ftest(10 000),
     3      w(4000 000),w7(1000000)
c 
        dimension uu(1000 000),wsout(1000),
     1      xsout(10 000),vv(1000 000),
     2      coefs(10 000),ftest2(10 000),errs(10 000)
  
        external fun1,fun2
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
C 
        PRINT *, 'ENTER npols '
        READ *,npols
        CALL PRINF('npols=*',npols,1 )
c 
        a=0
        b=0.5
        eps=1.0d-13
c 
c       construct the SVD of the user-specified functions
c 
        lab=20 000
        lw=1000 000
        lrlams=10 000
        igraph=0
  
        call fun1ini(nsings)
  
        call allsvcmb(ier,a,b,k,fun1,par1,par2,fun2,
     1      nsings,npols,eps,
     2      rlams,lrlams,w,lw,ncols,nn,rints,lused,keep,igraph)
c 
        n=nn*k
        call prinf('after allsvcmb, ier=*',ier,1)
        call prinf('after allsvcmb, lused=*',lused,1)
        call prinf('after allsvcmb, ncols=*',ncols,1)
        call prin2('after allsvcmb, rlams=*',rlams,ncols)
c 
c       calculate the singular functions at test points and plot them
c 
        ntest=1000
        h=(b-a)/ntest
        do 1400 i=1,ntest
        xtest(i)=a+h/2+(i-1)*h
 1400 continue
c 
        do 2000 i=1,ncols
c 
        do 1600 j=1,ntest
c 
        call nesteva(ier,w,xtest(j),i,ftest(j),der)
 1600 continue
c 
        iw=200+i
        call lotagraph(iw,xtest,ftest,ntest,
     1      'QR-function from coefs*')
c 
 2000 continue
c 
c       find a collection of points on the interval [a,b]
c       such that the matrix of values of singular vectors at these
c       points is well-conditioned
c 
        ifuv=1
        ncols2=ncols-2
c 
        call allspnts(ier,ncols2,w,
     1      xsout,wsout,uu,vv,ifuv,rcond,w(keep+1000 000) )
c 
        call prin2('after allspnts, rcond=*',rcond,1)
        call prin2('after allspnts, xsout=*',xsout,ncols2)
        call prin2('after allspnts, wsout=*',wsout,ncols2)
c 
c       now, test the obtained matrices uu, vv
c 
  
        i=7
        do 2600 j=1,ncols2
c 
        call nesteva(ier,w,xsout(j),i,ftest(j),der)
 2600 continue
c 
        call allqrmat(vv,ftest,coefs,ncols2)
  
        call prin2('coefs=*',coefs,ncols2)
  
        call allqrmat(uu,coefs,ftest2,ncols2)
  
        do 2200 i=1,ncols2
c 
        errs(i)=ftest(i)-ftest2(i)
 2200 continue
  
        call prin2('and ftest=*',ftest,ncols2)
        call prin2('and ftest2=*',ftest2,ncols2)
        call prin2('and errs=*',errs,ncols2)
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine fun2(x,i,f)
        implicit real *8 (a-h,o-z)
c 
        save
        call legepol(x,i-1,f,der)
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
cccc        complex *16 f
  
c 
c 
        f=x**sings(i)
  
  
  
cccc        f=pol
c 
        return
c 
c 
c 
c 
        entry fun1ini(nsings7)
c 
        nsings=nsings7
        done=1
        h=done/nsings
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
c        of left singular vectors) of families of functions, and
c        for the performance of certain associated operations. The
c        assumption is that the number of functions to be compressed
c        is too great for them to be discretized and evaluated simul-
c        taneously; the user specifies them things via subroutines
c        that evaluate them. The file contains 5 user-accesible
c        subroutines: allsvbld, allsvcmb, nesteva, nestev2, allspnts.
c        Following is a brief description of these subroutines.
c 
c   allsvbld - constructs the left singular functions
c       of the set of of user-specified functions on the interval
c       [a,b], and returns to the user nested Legendre expansions
c       for the first ncols singular functions; all elements with
c       higher number than ncols have weights associated with them
c       that are less than eps (user-specified). The subroutine also
c       returns to the user the singular values associated with
c       the said singular functions. Once this subroutine has
c       been called, any one of the obtained singular functions
c       (together with its derivative) can be evaluated by a call
c       to the subroutines nesteva, nesteva2
c 
c   allsvcmb - this is a highly specialized version of the
c       subroutine allsvbld (see above). In most situations where
c       one would think of using this subroutine, one should
c       actually use allsvbld!!! This subroutine should be used
c       when the input functions \phi_i to be SVD'ed are naturally
c       represented as products
c 
c            \phi_i (x) = \psi_i (x) + \khi_i (x),                    (1)
c 
c       and the number of functions \psi_i is much greater than
c       the (numerical) rank of their set. The subroutine starts
c       with constructing the SVD of the collection of functions
c       \psi_i; then it multiplies the obtained singular vectors
c       by the functions \khi_i, and constructs the SVDs of the
c       product (scaled by the corresponding singular values, of
c       course). When the number of functions \psi_i is much
c       greater than the (numerical) rank of their set, this
c       subroutine tends to be considerably faster than allsvbld.
c 
c   alldecod - this is the decoder subroutine: for a curious user,
c       it will retrieve from tha array w (produced either by
c       allsvbld or allsvcmb) various types of information that the
c       curious user might be curious about. The curious user is
c       urged to remember the fate of the curious cat.
c 
c   nesteva - evaluates at the point x \in [a,b]
c        one of singular functions specified by their nested Legendre
c        expansions, whose coefficients are stored in the input array
c        coefs by a prior call to one of the subroutines allsvbld,
c        allsvcmb. It also evaluates the derivative of the singular
c        function.
c 
c   nestev2 - identical to nesteva, with the exception that it does
c        not evaluate the derivative of the singular function
c 
c   allspnts - uses the output of the subroutine allsvbld
c       (or of allsvcmb) to construct the points xsout on the interval
c       [a,b] such that the matrix of values of left singular vectors
c       (obtained by allsvbld or allsvcmb) at the points xsout is
c       well-conditioned. Depending on the value of the parameter ifuv,
c       the subroutine will also produce the said matrix and its inverse.
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c 
        subroutine allsvcmb(ier,a,b,k,fun1,par1,par2,fun2,
     1      nsings,npols,eps,rlams,lrlams,w,lw,ncols,nn,
     2      rints,lused,keep,igraph)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),rlams(1),w(1),par1(1),par2(1),
     1      whts(1),rints(1)
c 
        external fun1,allfncmb,fun2
c 
C       ATTENTION!!!
C 
c       This subroutine is a highly specialized version of the
c       subroutine allsvbld (see). In most situations where
c       one would think of using this subroutine, one should
c       actually use allsvbld!!! This subroutine should be used
c       when the input functions \phi_i to be SVD'ed are naturally
c       represented as products
c 
c            \phi_i (x) = \psi_i (x) + \khi_i (x),                    (1)
c 
c       and the number of functions \psi_i is much greater than
c       the (numerical) rank of their set. The subroutine starts
c       with constructing the SVD of the collection of functions
c       \psi_i; then it multiplies the obtained singular vectors
c       by the functions \khi_i, and constructs the SVDs of the
c       product (scaled by the corresponding singular values, of
c       course). When the number of functions \psi_i is much
c       greater than the (numerical) rank of their set, this
c       subroutine tends to be considerably faster than allsvbld.
c 
c       This subroutine constructs the left singular functions
c       of the set of of user-specified functions on the interval
c       [a,b], and returns to the user nested Legendre expansions
c       for the first ncols singular functions; all elements with
c       higher number than ncols have weights associated with them
c       that are less than eps (user-specified). The subroutine also
c       returns to the user the singular values associated with
c       the said singular functions. Once this subroutine has
c       been called, any one of the obtained singular functions
c       (together with its derivative) can be evaluated by a call
c       to the subroutine nesteva (see).
c 
c       To this effect, the subroutine finds a subdivision of the
c       interval [a,b], such that on every subinterval, every one
c       of the user-provided collection of functions (given by the
c       subroutine fun (see description below) is integrated to
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
c  fun1 - the user-supplied subroutine evaluationg the functions
c      \psi_1, \psi_2, \cdot \psi_{nsings} in (1) above. The calling
c      sequence of fun1 must be:
c 
c      fun1(x,i,par1,par2,f),
c 
c       where x is the input argument where the function is to be
c       evaluated, i the sequence number of function, f (output)
c       is the value of the i-th function at the point x, and par
c       is an input parameter, carrying whatever information fun
c       might need.
c 
c  par1,par2 - whatever the user-supplied subroutine fun1 needs
c  fun2 - the user-supplied subroutine evaluationg the functions
c      \khi_1, \khi_2, \cdot \khi_{nsings} in (1) above. The calling
c      sequence of fun2 must be:
c 
c      fun2(x,i,f),
c 
c       where x is the input argument where the function is to be
c       evaluated, i the sequence number of function, and f (output)
c       is the value of the i-th function at the point x. Please note
c       that fun2 is different from fun1 in that it does not have any
c       additional parameters that can be used to pass to it auxiliary
c       information. Sorry, Jack, you will have to pass them through
c       secondary ENTRY statements (or, God forbid, a COMMON statement);
c       in the latter case, you will get what you deserve.
c  nsings - the number of functions \psi_i in (1) above
c  npols - the number of functions \khi_i in (1) above
c  eps - the accuracy to which the legendre expansion must represent
c      the user-supplied function for the subroutine to decide that
c      the representation is accurate
c  lrlams - the amount of space (in real *8 words) in the array rlams
c      in which the subroutine will return the singular values of the
c      matrix. This parameter is used by the subroutine to bomb if the
c      number of singular values that are greater than eps turns out
c      to be greater than the length of array rlams. Unless the user
c      knows what he is doing, lrlams should be set to nsings*npols
c 
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
c  nn - the number of subintervals in the nested discretization of
c      the interval [a,b].
c  w - the first n*ncols elements of array w contain the coefficients of
c        the nested Legendre expansions of the singular functions;
c 
c        IMPORTANT NOTE: the first n*ncols elements of array w should
c                        not be changed between the call to this
c                        subroutine and subsequent calls to the
c                        subroutine nesteva.
c  rlams - the singular values of the matrix of functions given by the
c        function fun. Note that though only ncols elements of array rlams
c        are meaningful, nfuns elements have to be allocated by the
c        user.
c  ncols - the number of QR-functions in the QR-decomposition of the
c      nfuns functions gven by the subroutine fun (see above) whose
c      corresponding weights are greater than eps (see above). It is
c      also the number of meaningful columns of arrays fs, coefs on
c      exit from this subroutine.
c  rints - the array of length ncols containing the integrals of the
c      first ncols singular functions of the operator
c  lused - the largest piece of the array w used by the subroutine at
c      any time during the execution
c  keep - the number of elements in array w that should not be changed
c      between the call to this subroutine and subsequent calls to the
c      subroutines nesteval, nesteva2, etc.
c 
c                  Work arrays:
c 
c  w - treated by the subroutine as both output parameter and a work
c      array. The subroutine utilizes the first lused elements of
c      w as a work space; after that, it returns in the first n*ncols
c      elements of w the coefficients of the nested legendre expansions
c      first ncols singular functions of the operator
c 
c       . . . construct the SVD of the user-specified functions,
c             without the multiplyers
c 
        ifsvd=1
        call allsvbld(ier,a,b,k,fun1,par1,par2,ifsvd,nsings,eps,
     1      rlams,lrlams,w,lw,ncols,nn,
     2      rints,lused,keep,igraph)
c 
        call prinf('after first allsvbld, ier=*',ier,1)
c 
c       initialize the compound svd routine
c 
        call allfnini(npols,ncols,rlams,nn,k)
c 
c       calculate the svd via the compound routine
c 
        nfuns=npols*ncols
        iw2=keep+1
        lw2=lw-iw2-10
c 
        call allsvbld(ier,a,b,k,allfncmb,w,fun2,ifsvd,nfuns,eps,
     1      rlams,lrlams,w(iw2),lw2,ncols,nn,
     2      rints,lused2,keep2,igraph)
c 
c        copy the coefficients of Legendre expansions of the
c        new singular functions to the beginning of array w
c 
        do 2200 i=1,keep2
        w(i)=w(iw2+i-1)
 2200 continue
c 
        keep=keep2
        call prinf('after second allsvbld, ier=*',ier,1)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine allsvbld(ier,a,b,k,fun,par1,par2,ifsvd,
     1      nfuns,eps,rlams,lrlams,w,lw,ncols,nn,
     2      rints,lused,keep,igraph)
        implicit real *8 (a-h,o-z)
        save
        dimension rlams(1),w(1),par1(1),par2(1),rints(1)
c 
        external fun
c 
c       This subroutine constructs the left singular functions
c       of the set of of user-specified functions on the interval
c       [a,b], and returns to the user nested Legendre expansions
c       for the first ncols singular functions; all elements with
c       higher number than ncols have weights associated with them
c       that are less than eps (user-specified). The subroutine also
c       returns to the user the singular values associated with
c       the said singular functions. Once this subroutine has
c       been called, any one of the obtained singular functions
c       (together with its derivative) can be evaluated by a call
c       to the subroutine nesteva (see).
c 
c       To this effect, the subroutine finds a subdivision of the
c       interval [a,b], such that on every subinterval, every one
c       of the user-provided collection of functions (given by the
c       subroutine fun (see description below) is integrated to
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
c  ifsvd - an integer parameter telling the subroutine whether it
c       should return to the user the true singular vectors, or
c       the QR vectors. Setting ifsvd = 1 will cause the subroutine
c       to return the SVD. Setting ifsvd=0 will cause the return of
c       QR vectors. The latter option is sometimes slightly less
c       expensive, and sometimes drastically less expensiove.
c  nfuns - the number of functions to be compressed
c  eps - the accuracy to which the legendre expansion must represent
c      the user-supplied function for the subroutine to decide that
c      the representation is accurate
c  lrlams - the amount of space (in real *8 words) in the array rlams
c      in which the subroutine will return the singular values of the
c      matrix. This parameter is used by the subroutine to bomb if the
c      number of singular values that are greater than eps turns out
c      to be greater than the length of array rlams
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
c  w - the first n*ncols elements of array w contain the coefficients of
c        the nested Legendre expansions of the singular functions;
c 
c        IMPORTANT NOTE: the first n*ncols elements of array w should
c                        not be changed between the call to this
c                        subroutine and subsequent calls to the
c                        subroutine nesteva.
c 
c  rlams - the singular values of the matrix of functions given by the
c        function fun. Note that though only ncols elements of array rlams
c        are meaningful, nfuns elements have to be allocated by the
c        user.
c  ncols - the number of QR-functions in the QR-decomposition of the
c      nfuns functions gven by the subroutine fun (see above) whose
c      corresponding weights are greater than eps (see above). It is
c      also the number of meaningful columns of arrays fs, coefs on
c      exit from this subroutine.
c  nn - the number of subintervals in the nested discretization of
c      the interval [a,b].
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
c       allocate memory for the subroutine manincon
c 
        ier=0
c 
        iiww=1
        llww=10*k+8*k**2+30
        llww=llww*2
c 
        iab=iiww+llww
        llab=lw-iab-2
c 
         call manincon(jer,a,b,k,fun,nfuns,par1,par2,
     1    eps,w(iab),llab,nn,w(iiww),llww)
c 
         if(jer .ne. 0) ier=16
         if(ier .ne. 0) return
c 
c       perform garbage collection
c 
        iab2=21
        call mvecmove(w(iab),w(iab2),nn*2)
        iab=iab2
        lab=nn*2+2
c 
        ix=iab+lab
        lx=k*nn+2
c 
        iwhts=ix+lx
        lwhts=k*nn+2
c 
c       . . . construct all the nested Legendre nodes
c 
        ier=0
        call nodewhts(w(iab),nn,k,w(ix),w(iwhts) )
c 
        n=nn*k
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
        iu=iwhts+lwhts
        lu=k**2+2
c 
        iv=iu+lu
        lv=k**2+2
c 
        iww=iv+lv
        lww=n+2
c 
        icoefs=iww+lww
        lencoefs=lw-icoefs-10
c 
        lused=iww+lww
c 
        call allffuns(jer,w(ix),w(iwhts),n,nfuns,fun,par1,par2,
     1      eps,w(icoefs),lencoefs,ncols,rlams,nn,k,w(iu),w(iv),
     2      w(iww),rints,lused2,lrlams,igraph,ifsvd)
c 
        lused=lused+lused2
c 
        if(jer .ne. 0) ier=8
        if(jer .ne. 0) return
c 
c       compress the memory
c 
        icoefs2=iwhts+lwhts
        call mvecmove(w(icoefs),w(icoefs2),ncols*n)
        icoefs=icoefs2
c 
c       store at the beginning of array w the map of the said array w
c 
        w(1)=iab+0.1
        w(2)=ix+0.1
        w(3)=iwhts+0.1
        w(4)=icoefs+0.1
        w(5)=nn+0.1
        w(6)=k+0.1
        w(7)=n+0.1
        w(8)=icoefs+0.1
c 
        keep=icoefs+ncols*n+200
c 
        w(9)=keep+0.1
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine alldecod(w,k,nn,n,ab,x,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),x(1),whts(1),ab(2,1)
c 
c       This is the decoder subroutine: for a curious user, it will
c       retrieve from tha array w (produced either by allsvbld or
c       allsvcmb) various types of information that the curious
c       user might be curious about. The user is urged to remember
c       the fate of the curious cat.
c 
c                            Input parameters:
c 
c  w - the large array produced either by allsvbld or allsvcmb)
c 
c                            Output parameters:
c 
c 
c  k - the number of Gaussian nodes on each subinterval
c  nn - the number of subintervals in the nested discretization of
c      the interval [a,b].  Note that this is an output parameter only
c      if the input parameter ifab (see above) had been set to 1 by
c      the user.
c  n - the number of nodes in the discretization x of the interval [a,b];
c        equal to nn*k
c  ab - the array of ends of subintervals created. For each i=1,2,...,nn,
c     ab(1,i) is the left end of the i-th subinterval; ab(2,i) is
c     the right end of the i-th subinterval.
c  x - the discretization of the interval [a,b] consisting of nn
c      subintervals with k Legendre nodes on each,  created by a
c      preceding call to this subroutine.
c  whts - the weights associated with the discretization x, created
c      by a preceding call to this subroutine.
c 
        k=w(6)
        nn=w(5)
        iab=w(1)
        ix=w(2)
        iwhts=w(3)
        n=w(7)
c 
        call mvecmove(w(ix),x,n)
        call mvecmove(w(iwhts),whts,n)
        call mvecmove(w(iab),ab,nn*2)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine allffuns(ier,x,whts,n,nfuns,fun,par1,par2,
     1      eps,fs,lfs,ncols,sss,nn,k,u,v,ww,rints,lused,lrlams,
     2      igraph,ifsvd)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),fs(n,1),whts(1),ww(1),par1(1),par2(1),
     1      u(1),v(1),xlege(2000),whtslege(2000),
     2      rints(1),sss(1)
c 
        external allsvmat
c 
c        this subroutine constructs the left singular vectors of
c        a set of user-specified functions on the interval [a,b],
c        and nested Legendre expansions for the first ncols of the
c         said singular functions all singular functions with higher
c        number than ncols have the singular values associated with
c        them that are less than eps (user-specified).
c 
c            Input parameters:
c 
c  x - the discretization of the interval [a,b] (normally produced
c        by the subroutine manincon (see)
c  whts - the weights associated with the discretization x; (normally
c        produced by the subroutine manincon (see)
c  n - the number of nodes in the discretization x of the interval [a,b]
c  nfuns - the number of functions to be SVD'ed
c  fun - the subroutine evaluating the functions to be SVD'ed at
c        arbitrary points on the interval [a,b]. The calling sequence
c        of fun must be:
c 
c      fun(x,i,par1,par2,f),
c 
c      where x is the input argument where the function is to be
c      evaluated, i the sequence number of function, f (output)
c      is the value of the i-th function at the point x, and par
c      is an input parameter, carrying whatever information fun
c      might need.
c 
c  par1,par2 - whatever the user-supplied subroutine fun needs
c  eps - the accuracy to which the calculations will be performed
c  lfs - the length of the array fs (in real *8 words). If the length
c      is insufficient, the error return code ier is set to 8, and
c      the execution of the subroutine is terminated
c  nn - the number of subintervals in the nested discretization of
c      the interval [a,b] (normally produced by the subroutine
c      manincon (see)
c  k - the number of Gaussian nodes on each of nn subintervals in the
c      nested subdivision of the interval [a,b]; note that n=nn*k
c      - it better be!!!
c 
c           Output parameters:
c 
c  ier - error return code. ier=0 means successful execution
c                           ier=8 means that the length of the user-
c      provided array fs (given by the parameter lfs) is insufficient.
c  fs - nested Legendre expansions of the the first ncols of the
c      singular functions. Note that the array fs is used by this
c      subroutine both as the master work array and as the principal
c      output array. Its available length is communicated to the
c      subroutine via the parameter lfs (see above). Also, note
c      that the coefficients of the nested Legendre expansions
c      occupy the first n*ncols elements of the array fs.
c  ncols - the number of QR-functions in the QR-decomposition of the
c      nfuns functions gven by the subroutine fun (see above) whose
c      corresponding weights are greater than eps (see above). It is
c      also the number of meaningful columns of arrays fs, coefs on
c      exit from this subroutine.
c  sss - the singular values corresponding to the singular vectors
c      in fs.
c 
c              Work arrays:
c 
c  u,v - must be at least k**2+2 real *8 elements long each.
c  ww - must be at least n real *8 elements each
c 
c        . . . constuct the square roots of weights
c 
        ier=0
c 
        do 1100 i=1,n
        ww(i)=sqrt(whts(i))
 1100 continue
c 
c       conduct the Gram-Schmidt procedure on the constructed functions
c 
        call mvecsvd(jer,allsvmat,fun,par1,par2,x,ww,n,
     1    nfuns,ifsvd,eps,ncols,sss,lrlams,fs,lfs,lused)
c 
        if(jer .ne. 0) ier=8
        if(ier .ne. 0) return
c 
        call prin2('in alffuns after mvecsvd, sss=*',sss,ncols)
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
        call prin2('integrals of fs evaluated in allffuns*',
     1      rints,ncols)
c 
c        if the user so requested, plot the obtained QR-vectors
c 
        if(igraph .le. 0) goto 2100
  
        do 2000 i=1,ncols
c 
        iw=igraph+i
c 
        call lotagraph(iw,x,fs(1,i),n,'QR-function*')
 2000 continue
c 
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
        call allqrmat(u,fs(jj,i),ww(jj),k)
c 
 2200 continue
c 
        do 2300 j=1,n
        fs(j,i)=ww(j)
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
        subroutine nestem(ier,coefs,nfuns,x,vals)
        implicit real *8 (a-h,o-z)
        save
        dimension pols(300),coefs(1)
        real *8 vals(1)
c 
c        this subroutine evaluates at the point x \in [a,b] the
c        nfuns functions specified by their nested Legendre
c        expansions, whose coefficients are supplied in the input
c        array coefs. Normally, the latter has been produced by a
c        prior call to the subroutine allsvbld, or one of its
c        relatives.
c 
c                   Input parameters:
c 
c  coefs - the nested Legendre expansions of the functions stored
c      in the array fs.
c  nfuns - the number of functions whose expansions are stored in
c      the array coefs; also the number of functions whose values
c      will be returned in array vals
c  x - the point where the functions are to be evaluated; must be
c      on the interval [a,b]
c 
c                    Output parameters:
c 
c  ier - error return code. ier = 0 means successful execution
c                           ier > 0 means that the point x is
c                           outside the interval [a,b]
c  vals - values of the nfuns functions evaluated at the point x
c 
c       . . . find the subinterval in which the point x lives
c 
c 
c       decode the beginning of the array coefs
c 
        ier=0
c 
        k=coefs(6)
        nn=coefs(5)
        iab=coefs(1)
        n=coefs(7)
        icoefs=coefs(8)
c 
c       . . . find the subinterval in which the point x lives
c 
        call findinte(ier,x,coefs(iab),nn,intnum)
        if(ier .ne. 0) return
c 
c       the point x lives in the interval number intnum.
c       evaluate the expansion at this point
c 
        iii=iab+intnum*2-1
        u=2/(coefs(iii)-coefs(iii-1))
        v=1-coefs(iii)*u
c 
        t=u*x+v
c 
        jj=(intnum-1)*k+1
c 
c        evaluate all Legendre polynomials of t
c 
        call legepls2(t,k,pols)
c 
c       evaluate all expansions
c 
        do 1200 i=1,nfuns
c 
        vals(i)=0
 1200 continue
c 
        jj0=(intnum-1)*k
c 
        call nestmatv(coefs(icoefs),k,nn,nfuns,pols,vals,jj0)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine nestmatv(coefs,k,nn,nfuns,pols,vals,jj0)
        implicit real *8 (a-h,o-z)
        save
        dimension pols(1)
        real *8 coefs(nn*k,nfuns),vals(1)
c 
        do 1600 j=1,k
c 
        jj=jj0+j
c 
        do 1400 i=1,nfuns
c 
        vals(i)=vals(i)+coefs(jj,i)*pols(j)
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
        subroutine nesteva(ier,coefs,x,i,val,der)
        implicit real *8 (a-h,o-z)
        save
        dimension coefs(1),pjcoefs1(300),pjcoefs2(300)
        data ifcalled/0/
c 
c       decode the beginning of the array coefs
c 
        k=coefs(6)
        nn=coefs(5)
        iab=coefs(1)
        n=coefs(7)
        icoefs=coefs(8)
c 
c       . . . find the subinterval in which the point x lives
c 
        ier=0
        call findinte(ier,x,coefs(iab),nn,intnum)
        if(ier .ne. 0) return
c 
c       the point x lives in the interval number intnum.
c       evaluate the expansion at this point
c 
        iii=iab+intnum*2-1
        u=2/(coefs(iii)-coefs(iii-1))
        v=1-coefs(iii)*u
c 
        t=u*x+v
c 
        jj=(intnum-1)*k+1
        ninit=290
        if(ifcalled .eq. 1) ninit=0
        ifcalled=1
c 
        iijj=(i-1)*nn*k+jj +icoefs -1
c 
        call legeFDE2(t,VAL,der,coefs(iijj),k-1,
     1    pjcoefs1,pjcoefs2,ninit)
c 
        der=der*u
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine nestev2(ier,coefs,x,i,val)
        implicit real *8 (a-h,o-z)
        save
        dimension coefs(1),pjcoefs1(300),pjcoefs2(300)
        data ifcalled/0/
c 
c       decode the beginning of the array coefs
c 
        k=coefs(6)
        nn=coefs(5)
        iab=coefs(1)
        n=coefs(7)
        icoefs=coefs(8)
c 
c       . . . find the subinterval in which the point x lives
c 
        ier=0
        call findinte(ier,x,coefs(iab),nn,intnum)
        if(ier .ne. 0) return
c 
c       the point x lives in the interval number intnum.
c       evaluate the expansion at this point
c 
        iii=iab+intnum*2-1
        u=2/(coefs(iii)-coefs(iii-1))
        v=1-coefs(iii)*u
c 
        t=u*x+v
c 
        jj=(intnum-1)*k+1
c 
        ninit=290
        if(ifcalled .eq. 1) ninit=0
        ifcalled=1
c 
        iijj=(i-1)*nn*k+jj +icoefs -1
        call legeexe2(t,VAL,coefs(iijj),k-1,
     1      pjcoefs1,pjcoefs2,ninit)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine findinte(ier,x,ab,nn,intnum)
        implicit real *8 (a-h,o-z)
        save
        dimension ab(2,nn)
c 
        data intold/-10/,ithresh/10/
c 
c       check if the point is on the subinterval as the preceding one
c 
        ier=0
        if(intold .lt. 0) goto 2000
c 
        intnum=intold
        if( (x .ge. ab(1,intnum) ) .and. (x .le. ab(2,intnum) ) ) return
c 
 2000 continue
c 
c      the point is not on the same subinterval as the preceding one.
c      if nn is less than ithresh, use direct scan to find the proper
c      interval
c 
       if(nn .gt. ithresh) goto 3000
c 
       if(x .lt. ab(1,1)) then
           ier=4
           return
       endif
c 
        do 2200 j=1,nn
c 
           intnum=j
c 
        if(ab(2,j) .ge. x) goto 2400
 2200 continue
c 
        ier=4
        return
 2400 continue
c 
        intold=intnum
        return
c 
 3000 continue
c 
c      The point is not on the same subinterval as the preceding one,
c      and nn is greater than ithresh; use bisection to find the proper
c      interval
c 
       i1=1
       i2=nn
       i3=(i1+i2)/2
c 
cccc       nsteps=0
       do 3400 i=1,100
c 
       if(x .ge. ab(1,i3)) i1=i3
       if(x .le. ab(2,i3)) i2=i3
c 
       if(i2 .eq. i1) goto 3600
c 
       i3=(i1+i2)/2
 3400 continue
c 
 3600 continue
  
       if(x .lt. ab(1,i3)) i3=i3-1
       if(x .gt. ab(2,i3)) i3=i3+1
  
       intnum=i3
       intold=intnum
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine allsvmat(i,j,f,fun,par1,par2,xs,ww)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),ww(1)
c 
        call fun(xs(i),j,par1,par2,f)
c 
        f=f*ww(i)
        return
        end
c 
c 
c 
c 
c 
        subroutine manincon(ier,a,b,k,fun,nfuns,par1,par2,
     1    eps,ab,lab,nn,w,lenw)
        implicit real *8 (a-h,o-z)
        save
        dimension par1(1),ab(2,1),w(1),par2(1)
c 
c       this subroutine finds a subdivision of the interval [a,b],
c       such that on every subinterval, the every one of the
c       user-provided collection of functions given by the function
c       fun (see description below) is integrated to precision eps
c       by a k-point Gaussian quadrature.
c 
c                          input parameters:
c 
c  a,b - the ends of the interval to be subdivided
c  k - the number of Gaussian nodes on each subinterval
c  fun - the user-supplied function to be integrated
c      The calling sequence of fun must be:
c 
c      fun(x,i,par1,par2,f),
c 
c      where x is the input argument where the function is to be
c      evaluated, i the sequence number of function, f (output)
c      is the value of the i-th function at the point x, and par
c      is an input parameter, carrying whatever information fun
c      might need.
c 
c  par1,par2 - whatever the user-supplied subroutine fun needs
c  eps - the accuracy to which the legendre expansion must represent
c     the user-supplied function for the subroutine to decide that
c     the representation is accurate
c  lab - the length of the user-provided array ab (see below)
c  w - the array containing the data to be used by the subroutine,
c     created by a preceding call to it. It must be at
c     least 10*k+8*k**2+30 + 20 000 real *8 locations long.
c 
c                     output parameters:
c 
c  ier - error return code.
c         ier=0 means successful execution
c         ier=4 means that the space allocated in the array ab is
c         insufficient
c  ab - the array of ends of subintervals created. For each i=1,2,...,nn,
c     ab(1,i) is the left end of the i-th subinterval; ab(2,i) is
c     the right end of the i-th subinterval
c  nn - the number of subintervals created on the interval [a,b];
c     also the number of entries in the array ab (see above)
c  labused - the total amount of space used by the subroutine in the
c     array ab; note that normally, it  is somewhat greater than
c     the expected nn*2; The difference is not large, but should
c     be taken into account when allocating memory.
c  w - the array containing the data used by the subroutine,
c     It is only an output parameter if the user has set the
c     parameter ifinit (see above) to 1. Otherwise, it is an input
c     parameter, and must have been created  by a prior call to
c     this subroutine. Note that it should not be altered between
c     consecutive calls to this subroutine; the piece
c     used by this subroutine is 10*k+8*k**2+30 real *8 locations long.
c 
c 
c        one after another, construct the subdivisions
c        of the interval [a,b] corresponding to the
c        user-specified functions to be integrated, and merge
c        all these intervals
c 
c        . . . initialize the merging process
c 
        ier=0
        ipar1=1
        ifinit=1
c 
        iw=1
        lw=10*k+8*k**2+30
c 
        iab1=iw+lw
        lab1=lenw-iab1
c 
         call allincon(ier,a,b,k,fun,ipar1,par1,par2,eps,
     1    w(iab1),lab1,nn1,labused,w(iw),ifinit)
c 
        if(ier .ne. 0) return
c 
c        . . . merge
c 
        do 2000 i=2,nfuns
c 
        ipar1=i
        ifinit=0
c 
        iab2=iab1+nn1*2+2
        lab2=lenw-iab2
c 
         call allincon(ier,a,b,k,fun,ipar1,par1,par2,eps,
     1    w(iab2),lab2,nn2,labused,w(iw),ifinit)
c 
        if(ier .ne. 0) return
c 
        iab3=iab2+nn2*2+2
        lab3=(nn1+nn2)*2+4
c 
        iww=iab3+lab3
        lleft=lenw-iww
        if(lleft .gt. (nn1+nn2)*2+4) goto 1600
        ier=4
        return
 1600 continue
c 
        call  mrgstruc(w(iab1),nn1,w(iab2),nn2,
     1      w(iab3),nn3,eps,w(iww))
c 
        call arrmove(w(iab3),w(iab1),nn3*2)
        nn1=nn3
c 
 2000 continue
c 
        call arrmove(w(iab3),ab,nn3*2)
        nn=nn3
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine arrmove(x,y,n)
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
        subroutine nodewhts(ab,nn,k,x,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension ab(2,1),x(1),whts(1),t(2000),w(2000)
c 
c 
        data kold/-10/
c 
c       construct the nodes and weights of the k-point Gaussian
c       quadrature on the interval [-1,1]
c 
        itype=1
c 
        if(k .ne. kold) call legeexps(itype,k,t,u,v,w)
c 
        kold=k
c 
c       one subinterval after another, construct nodes and weights
c       of the quadrature on the big interval
c 
        do 2000 i=1,nn
        u=(ab(2,i)-ab(1,i))/2
        v=(ab(2,i)+ab(1,i))/2
        ik=(i-1)*k
c 
        do 1200 j=1,k
        x(ik+j)=t(j)*u+v
        whts(ik+j)=w(j)*u
 1200 continue
 2000 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine mrgstruc(ab1,nn1,ab2,nn2,ab3,nn3,eps,w)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),ab1(1),ab2(1),ab3(1)
c 
c        this subroutine merges the user-provided subdivisions to
c        get one single subdivision that is as dense as either of
c        the original ones at any point on the (common) interval
c        of their definition
c 
c                        input parameters:
c 
c  ab1 - the left and right end of subintervals in the
c        first subdivision
c  nn1 - the number of subintervals in the array ab1
c  ab2 - the left and right end of subintervals in the
c        second subdivision
c  nn2 - the number of subintervals in the array ab2
c  eps - tells the subroutine how short a subinterval must
c        be to be considered of length zero
c 
c                        input parameters:
c 
c  ab3 - the left and right end of subintervals in the
c        merged subdivision
c  nn3 - the number of subintervals in the array ab3
c 
c                        work arrays:
c 
c  w - must be at least (nn1+nn2)*2+2 real *8 elements long
c 
c 
c        . . . merge the arrays ab1,ab2 (both left and right ends)
c              into one long array with non-decreasing elements
c 
        call twomrge(ab1,nn1*2,ab2,nn2*2,w)
c 
c        now, scan the obtained merged array and discard repeating
c        elements
c 
        nnn=nn1*2+nn2*2
        ii=1
        ab3(ii)=w(1)
        do 1400 i=2,nnn
        d=w(i)-w(i-1)
        if(d .lt. eps) goto 1400
        ii=ii+1
        ab3(ii)=w(i)
 1400 continue
c 
c       now, duplicate all interior elements in array ab3
c 
        w(1)=ab3(1)
        iii=1
        do 1600 i=2,ii-1
        iii=iii+1
        w(iii)=ab3(i)
        iii=iii+1
        w(iii)=ab3(i)
 1600 continue
        w(iii+1)=ab3(ii)
        iii=iii+1
c 
c       copy the final array to where it belongs
c 
        nn3=iii/2
        do 1800 i=1,iii
        ab3(i)=w(i)
 1800 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine twomrge(a,na,b,nb,c)
        implicit real *8 (a-h,o-z)
        save
        real *8  a(1),b(1),c(1)
c 
c       merge the sorted arrays ia, ib obtaining the sorted array ic
c 
        ia=1
        ib=1
        ic=1
c 
        do 1600 i=1,na+nb
        ic=i
c 
        if(a(ia) .gt. b(ib) ) goto 1200
        c(i)=a(ia)
        ia=ia+1
        goto 1400
c 
 1200 continue
        c(i)=b(ib)
        ib=ib+1
 1400 continue
c 
        if(ib .gt. nb) goto 1800
        if(ia .gt. na) goto 2200
 1600 continue
c 
 1800 continue
c 
c       we have run out of elemnets in b. finish off the array a
c 
        do 2000 i=1,na-ia+1
        c(ic+i)=a(ia+i-1)
 2000 continue
        return
c 
 2200 continue
c 
c       we have run out of elements in a. finish off the array b
c 
        do 2400 i=1,nb-ib+1
        c(ic+i)=b(ib+i-1)
 2400 continue
        return
        end
c 
c 
c 
c 
c 
         subroutine allincon(ier,a,b,k,fun,ipar1,par1,
     1    par2,eps,ab,lab,nn,labused,w,ifinit)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),par1(1),par2(1),ab(2,1)
c 
c       this subroutine finds a subdivision of the interval [a,b],
c       such that on every subinterval, the user-provided function
c       fun is integrated to precision eps by a k-point Gaussian
c       quadrature. Actually, this is a memory management routine;
c       all the work is done by the routine allinter (see).
c 
c                          input parameters:
c 
c  a,b - the ends of the interval to be subdivided
c  k - the number of Gaussian nodes on each subinterval
c  fun - the user-supplied function to be integrated
c      The calling sequence of funget must be:
c 
c      fun(x,i,par1,par2,f),
c 
c      where x is the input argument where the function is to be
c      evaluated, i is the sequence  number of the function to be
c      evaluated, f  (output) is the function, and par1,par2,apr3
c      are input parameters, carrying whatever information fun
c      might need.
c 
c  par1,par2 - whatever the user-supplied subroutine fun needs
c  eps - the accuracy to which the legendre expansion must represent
c     the user-supplied function for the subroutine to decide that
c     the representation is accurate
c  lab - the length of the user-provided array ab (see below)
c  w - the array containing the data to be used by the subroutine,
c     created by a preceding call to it. It is only an input parameter
c     if the user has set the parameter ifinit (see below) to 0.
c     Otherwise, it is an output parameter, and must be at
c     least 10*k+8*k**2+30 real *8 locations long.
c  ifinit - the index telling the subroutine whether it should construct
c     the parameter w (see above), or view it as an input parameter.
c     In the latter case, this  must have been created by a prior call
c     to this subroutine with the same parameter k.
c     ifinit=1 means that the parameter w will be created.
c     ifinit=0 means that the parameter w will be viewed
c     as an input parameter.
c 
c                     output parameters:
c 
c  ier - error return code.
c         ier=0 means successful execution
c         ier=4 means that the space allocated in the array ab is
c         insufficient
c  ab - the array of ends of subintervals created. For each i=1,2,...,nn,
c     ab(1,i) is the left end of the i-th subinterval; ab(2,i) is
c     the right end of the i-th subinterval
c  nn - the number of subintervals created on the interval [a,b];
c     also the number of entries in the array ab (see above)
c  labused - the total amount of space used by the subroutine in the
c     array ab; note that normally, it  is somewhat greater than
c     the expected nn*2; The difference is not large, but should
c     be taken into account when allocating memory.
c  w - the array containing the data used by the subroutine,
c     It is only an output parameter if the user has set the
c     parameter ifinit (see above) to 1. Otherwise, it is an input
c     parameter, and must have been created  by a prior call to
c     this subroutine. Note that it should not be altered between
c     consecutive calls to this subroutine; the piece
c     used by this subroutine is 10*k+8*k**2+30 real *8 locations long.
c 
c       . . . allocate memory for the subroutine allinter
c 
        ier=0
c 
        it=1
        lt=k*2+2
c 
        ix=it+lt
        lx=k*2+2
c 
        if=ix+lx
        lf=k*2+2
c 
        icoefs=if+lf
        lcoefs=k*2+2
c 
        iwhts=icoefs+lcoefs
        lwhts=k*2+1
c 
        iu=iwhts+lwhts
        lu=k**2*4 +10
c 
        iv=iu+lu
        lv=k**2*4 +10
c 
c       recursively subdivide the interval [a,b] into subintervals
c       such that on each of them, the k-point Gaussian rule
c       will obtain the accuracy epc
c 
         call allinter(ier,a,b,fun,ipar1,par1,par2,k,ab,lab,
     1    eps,w(it),w(iu),w(iv),w(iwhts),w(ix),w(if),
     2       w(icoefs),nn,labused,ifinit)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine bubble(ab,nn)
        implicit real *8 (a-h,o-z)
        save
        dimension ab(2,1)
c 
c       sort the array ab with respect to the first coordinate
c 
        do 2000 i=1,nn
        do 1200 j=1,i-1
        if(ab(1,i) .ge. ab(1,j) ) goto 1200
        d=ab(1,i)
        ab(1,i)=ab(1,j)
        ab(1,j)=d
c 
        d=ab(2,i)
        ab(2,i)=ab(2,j)
        ab(2,j)=d
 1200 continue
 2000 continue
        return
        end
c 
c 
c 
c 
c 
         subroutine allinter(ier,a,b,fun,ipar1,par1,par2,
     1    k,ab,lab,eps,t,u,v,whts,x,f,coefs,nn,lused,ifinit7)
        implicit real *8 (a-h,o-z)
        save
        dimension t(1),u(2*k,2*k),v(2*k,2*k),whts(1),x(1),
     1      f(1),coefs(1),par1(1),par2(1),ab(2,1)
        equivalence (rea,intint)
c 
c       this subroutine finds a subdivision of the interval [a,b],
c       such that on every subinterval, the user-provided function
c       fun is integrated to precision eps by a k-point Gaussian
c       quadrature.
c 
c                          input parameters:
c 
c  a,b - the ends of the interval to be subdivided
c  fun - the user-supplied function to be integrated
c      The calling sequence of funget must be:
c 
c      funget(x,par1,par2,f),
c 
c      where x is the input argument where the function is to be
c      evaluated, f  (output) is the function, and par1,par2 are
c      both input parameters, carrying whatever information funget
c      might need.
c 
c  par1,par2 - whatever the user-supplied subroutine funget needs
c  k - the number of Gaussian nodes on each subinterval
c  lab - the length of the user-provided array ab (see below)
c  eps - the accuracy to which the legendre expansion must represent
c     the user-supplied function for the subroutine to decide that
c     the representation is accurate
c  t - the Gaussian nodes on the interval [-1,1]. This is only an
c     input parameter is ifinit7 (see below) has been set to 0;
c     otherwise, it is an output parameter.
c  u - the 2*k * 2*k  matrix converting the  values at of a polynomial
c     of order 2*k-1 at 2*k legendre nodes into the coefficients of its
c     legendre expansion. This is only an input parameter is ifinit7
c     (see below) has been set to 0; otherwise, it is an output parameter.
c  v - the 2*k * 2*k matrix converting the coefficients
c     of an 2*k-term legendre expansion into its values at 2*k legendre
c     nodes (note that v is the inverse of u). This is only an input
c     parameter is ifinit7 (see below) has been set to 0; otherwise, it
c     is an output parameter.
c  whts - the Gaussian quadrature weights (2*k of them).  This is only
c     an input parameter is ifinit7 (see below) has been set to 0;
c     otherwise, its an output parameter.
c  ifinit7 - the index telling the subroutine whether it should construct
c     the parameters t,u,v,whts (by calling the subroutine legeexps),
c     or view them as input parameters. In the latter case, these
c     must have been created by a prior call to this subroutine (or to
c     subroutine legeexps, for an adventurous and creative user).
c     ifinit7=1 means that the parameters t,u,v,whts will be created.
c     ifinit7=0 means that the parameters t,u,v,whts will be viewed
c     as input parameters.
c 
c                     output parameters:
c 
c  ier - error return code.
c         ier=0 means successful execution
c         ier=4 means that the space allocated in the array ab is
c         insufficient
c  ab - the array of ends of subintervals created. For each i=1,2,...,nn,
c     ab(1,i) is the left end of the i-th subinterval; ab(2,i) is
c     the right end of the i-th subinterval
c  t,u,v,whts - are only output is the parameter ifinit (see above)
c     has been set by the user to 1.
c  x - the Gaussian nodes on the interval [a,b]
c  f - the values of the user-specified function fun at
c     the nodes x
c  coefs - the coefficients of the Legendre expansion of the
c     function fun on the interval [a,b]
c  nn - the number of subintervals created on the interval [a,b];
c     also the number of entries in the array ab (see above)
c  lused - the total amount of space used by the subroutine in the
c     array ab; note that normally, it  is somewhat greater than
c     the expected nn*2; The difference is not large, but should
c     be taken into account when allocating memory.
c 
c       . . . start the recursive subdivision process
c 
        ab(1,1)=a
        ab(2,1)=b
        nn=0
        ifinit=ifinit7
        istore=1
        ier=0
c 
c       recursively subdivide the intervals till none need
c       to be subdivided
c 
        lused=0
        do 2000 i=1,1 000 000
c 
        nn=nn+1
c 
c       check if any more intervals need to be subdivided
c 
        if(nn .gt. istore) goto 2200
c 
c       check if this interval needs to be subdivided
c 
        d=ab(2,nn)-ab(1,nn)
c 
        call ifsplit(ier,ab(1,nn),ab(2,nn),k,fun,ipar1,par1,
     1    par2,eps/d,t,u,v,whts,x,ifinit,f,coefs)
c 
        ifinit=0
c 
        if(ier .eq. 0) goto 2000
c 
c       this interval needs to be subdivided. act accordingly
c 
        if(lused .lt. istore*2+4) lused=istore*2+4
        if(istore*2 +4 .lt. lab) goto 1800
        ier=4
        return
 1800 continue
c 
        istore=istore+1
        ab(1,istore)=ab(1,nn)
        ab(2,istore)=(ab(1,nn)+ab(2,nn))/2
c 
        istore=istore+1
        ab(1,istore)=ab(2,istore-1)
        ab(2,istore)=ab(2,nn)
c 
        intint=9654236
        ab(1,nn)=rea
        ab(2,nn)=1.0d35
 2000 continue
 2200 continue
c 
        nn=nn-1
c 
c       scan the intervals we have created, discarding those that have kids
c 
        ii=0
        do 2400 i=1,nn
        rea=ab(1,i)
c 
        if( (ab(2,i) .gt. 1.0d34 ) .and.
     1      (intint .eq. 9654236) ) goto 2400
        ii=ii+1
        ab(1,ii)=ab(1,i)
        ab(2,ii)=ab(2,i)
 2400 continue
        nn=ii
c 
c        bubble-sort the array ab
c 
        call bubble(ab,nn)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine ifsplit(ier,a,b,k,fun,ipar1,par1,
     1    par2,eps,t,u,v,whts,x,ifinit,f,coefs)
        implicit real *8 (a-h,o-z)
        save
        dimension t(1),u(2*k,2*k),v(2*k,2*k),whts(1),x(1),
     1      f(1),coefs(1),par1(1),par2(1)
c 
c       this subroutine determines whether a k-point legendre
c       expansion of the user-supplied function f on the interval
c       [a,b] represents it with accuracy eps.
c 
c                          input parameters:
c 
c  a,b - the ends of the interval on which the approximation is checked
c  k - the order of the approximation for which the accuracy is checked
c  funget - the user-supplied function for which the accuracy is checked
c      The calling sequence of funget must be:
c 
c      fun(x,i,par1,par2,f),
c 
c      where x is the input argument where the function is to be
c      evaluated, f  (output) is the function, and par1,par2 are
c      both input parameters, carrying whatever information funget
c      might need.
c 
c  par1,par2 - whatever the user-supplied subroutine funget needs
c  eps - the accuracy to which the legendre expansion must represent
c     the user-supplied function for the subroutine to decide that
c     the representation is accurate
c  t - the Gaussian nodes on the interval [-1,1]. This is only an
c     input parameter is ifinit (see below) has been set to 0;
c     otherwise, it is an output parameter.
c  u - the 2*k * 2*k  matrix converting the  values at of a polynomial
c     of order 2*k-1 at 2*k legendre nodes into the coefficients of its
c     legendre expansion. This is only an input parameter is ifinit
c     (see below) has been set to 0; otherwise, it is an output parameter.
c  v - the 2*k * 2*k matrix converting the coefficients
c     of an 2*k-term legendre expansion into its values at 2*k legendre
c     nodes (note that v is the inverse of u). This is only an input
c     parameter is ifinit (see below) has been set to 0; otherwise, it
c     is an output parameter.
c  whts - the Gaussian quadrature weights (2*k of them).  This is only
c     an input parameter is ifinit (see below) has been set to 0;
c     otherwise, its an output parameter.
c  ifinit - the index telling the subroutine whether it sholud construct
c     the parameters t,u,v,whts (by calling the subroutine legeexps),
c     or view them as input parameters. In the latter case, these
c     must have been created by a prior call to this subroutine (or to
c     subroutine legeexps, for an adventurous and creative user).
c     ifinit=1 means that the parameters t,u,v,whts will be created.
c     ifinit=0 means that the parameters t,u,v,whts will be viewed
c     as input parameters.
c 
c                     output parameters:
c 
c  t,u,v,whts - are only output is the parameter ifinit (see above)
c     has been set by the user to 1.
c  x - the Gaussian nodes on the interval [a,b]
c  f - the values of the user-specified function funget at
c     the nodes x
c  coefs - the coefficients of the Legendre expansion of the
c     function f on the interval [a,b]
c 
c       . . . if the need be, construct the gaussian nodes on the
c             interval [-1,1], and the matrix conevrting the values
c             of the function at these nodes into coefficients of
c             its legendre expansion
c 
        if(ifinit .eq. 0) goto 1200
c 
        itype=2
        call legeexps(itype,k*2,t,u,v,whts)
c 
 1200 continue
c 
c        discretize the interval and evaluate the
c        user-supplied function at its nodes
c 
        alpha=(b-a)/2
        beta=(b+a)/2
c 
        do 1400 i=1,k*2
        x(i)=alpha*t(i)+beta
cccc        call fun(x(i),ipar1,par1,par2,par3,f(i) )
        call fun(x(i),ipar1,par1,par2,f(i) )
 1400 continue
cccc        call prin2('in ifsplit, x=*',x,k*2)
cccc        call prin2('in ifsplit, f=*',f,k*2)
c 
c       construct the legendre expansion of the function on
c       the interval
c 
cccc        call allqrmat(u,f,coefs,k*2)
  
        call allqrma2(u,f,coefs,k*2,k+1)
c 
c       scan the coefficients and see if the last k
c       of them are small
c 
        ier=0
        do 1800 i=1,k
        d=dabs(coefs(k+i))
        if(d .gt. eps) ier=4
 1800 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine allqrma2(a,x,y,n,n0)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),x(1),y(1)
c 
        do 2000 i=n0,n
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
        subroutine allqrmat(a,x,y,n)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),x(1),y(1)
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
        subroutine allqrmata(a,x,y,n)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),x(1),y(1)
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
        subroutine allspnts(ier,ncols,coefs,
     1      xsout,wsout,uu,vv,ifuv,rcond,w)
        implicit real *8 (a-h,o-z)
        save
        dimension uu(ncols,ncols),w(1),xsout(1),wsout(1),
     1      vv(ncols,ncols),coefs(1)
c 
c       This subroutine uses the output of the subroutine allsvbld
c       (or of allsvcmb) to construct the points xsout on the interval
c       [a,b] such that the matrix of values of left singular vectors
c       (obtained by allsvbld or allsvcmb) at the points xsout is
c       well-conditioned. Depending on the value of the parameter ifuv,
c       the subroutine will also produce the said matrix and its inverse.
c 
c     PLEASE NOTE THAT THIS SUBROUTINE CAN ONLY BE USED AFTER ONE OF
c     SUBROUTINES ALLSVBLD, ALLSVCMB HAS BEEN CALLED. IT HAS NO USES
c     WHATSOEVER AS A STAND-ALONE DEVICE!!!!!!!!!!!!!!!!
c 
c                       Input parameters:
c 
c  ncols - the number of left singular vectors to be used (also the
c       number of nodes to construct). Must be less than or equal to
c       the number ncols returned by the preceding call to allsvbld or
c       allsvcmb.
c  coefs - the output of a preceding call to allsvbld or allsvcmb
c  ifuv - tells the subroutine whether to construct the matrix uu of values
c        of the singular function at the selected points xsout, and the
c        inverse vv of the matrix uu.
c 
c                        Output parameters:
c 
c  xsout - a collection of ncols points on the interval [a,b] such that
c        the matrix of values of the first ncols singular vectors at these
c        points is well-conditioned
c  wsout - the weights corresponding to the nodes xsout; only produced
c        if the input parameter ifuv has been set to 1
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
        keep=coefs(9)
        n=coefs(7)
        ix=coefs(2)
        iwhts=coefs(3)
c 
        iuuadj=1
        luuadj=n*ncols+10
c 
        iipivots=iuuadj+luuadj
        lipivots=n+2
c 
        irnorms=iipivots+lipivots
        lrnorms=n+2
c 
c        find a collection of points such that the matrix of values of
c        singular vectors at these points is well-conditioned
c 
        call prinf('in allspnts before allspnts0, n=*',n,1)
  
        call allspnt0(ier,uu,n,ncols,w(iuuadj),coefs(ix),coefs(iwhts),
     1      xsout,wsout,w(iipivots),w(irnorms),coefs,rcond,vv,ifuv)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine allspnt0(ier,uu,n,ncols,uuadj,xs,whts,
     1      xsout,wsout,ipivots,rnorms,w,rcond,vv,ifuv)
        implicit real *8 (a-h,o-z)
        save
        dimension uu(ncols,ncols),uuadj(ncols,n),rnorms(1),
     1      ipivots(1),whts(1),xs(1),w(1),xsout(1),
     2      vv(ncols,ncols),rints(1000),wsout(1)
c 
c       evaluate all of the singular vectors at the nodes of
c       the discretization xs of the interval of definition
c 
        do 1100 i=1,ncols
c 
        rints(i)=0
 1100 continue
  
cccc         call prin2('in allspnt03, rints zeroed=*',rints,ncols)
  
c 
        do 1400 i=1,ncols
        do 1200 j=1,n
c 
        call nesteva(ier,w,xs(j),i,uuadj(i,j),der)
 1200 continue
 1400 continue
c 
        do 1450 j=1,n
        do 1430 i=1,ncols
c 
        rints(i)=rints(i)+uuadj(i,j)*whts(j)
        uuadj(i,j)=uuadj(i,j)*sqrt(whts(j))
c 
 1430 continue
 1450 continue
  
         call prin2('in allspnt03, rints=*',rints,ncols)
c 
        eps2=1.0d-30
        call allsppiv(uuadj,ncols,n,rnorms,eps2,ncols3,ipivots)
c 
c       sort the array ipivots
c 
        call allsbubb(ipivots,ncols3)
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
        call nesteva(ier,w,xsout(j),i,uu(j,i),der)
c 
        vv(j,i)=uu(j,i)
 1600 continue
 1800 continue
c 
c       invert the matrix uu
c 
        call orthom(vv,ncols,uuadj,rcond)
c 
c        evaluate the weights corresponding to the nodes xsout
c 
        call allqrma5(vv,rints,wsout,ncols)
  
cccc        call prin2('and wsout=*',wsout,ncols)
  
  
        return
        end
c 
c 
c 
c 
c 
c 
        subroutine allqrma5(a,x,y,n)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),x(1),y(1)
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
        subroutine allsppiv(b,n,m,rnorms,eps,ncols,ipivots)
        implicit real *8 (a-h,o-z)
        save
        dimension b(n,m),rnorms(1),ipivots(1)
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
c  ipivots - the sequence numbers of the pivots used, in the order in
c        which they were used
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
        call mvecscap(b(1,i),b(1,j),n,d)
c 
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*d
 2770 continue
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
        call mvecscap(b(1,i),b(1,j),n,d)
  
cccc        subroutine mvecscap(x,y,n,prod)
  
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
        subroutine allsbubb(ipivots,ncols)
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
        subroutine allfncmb(x,i,coefs,fun2,f)
        implicit real *8 (a-h,o-z)
        save
        dimension singvals(10000),singval7(1),coefs(1)
c 
        external fun2
c 
c        This subroutine evaluates a composite function
c        number i at the point x. The function is composite
c        because it is a product of two functions: the first
c        one is evaluated by the user-provided subroutine fun2,
c        and the second one is evaluates by the subroutine
c        nestev2, using a nested Legendre expansion of some
c        other functions constructed earlier.
c 
        iord=(i-1)/nsings+1
        ising=i-(iord-1)*nsings
c 
        call fun2(x,iord,pol)
c 
        call nestev2(ier,coefs,x,ising,f)
c 
        f=f*singvals(ising)*pol
c 
        return
c 
c 
c 
c 
        entry allfnini(npols7,nsings7,singval7,nn7,k7)
c 
        nn=nn7
        k=k7
        npols=npols7
        nsings=nsings7
c 
        do 2200 j=1,nsings
        singvals(j)=singval7(j)
 2200 continue
c 
        return
        end
