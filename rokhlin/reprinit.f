        implicit real *8 (a-h,o-z)
        dimension w(12 000 000)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER nfuns'
        READ *,nfuns
        CALL PRINf('nfuns=*',nfuns,1 )
c 
        eps=1.0d-14
  
c 
        y0=0
        y1=6
c 
cccc        y0=100
        y1=y0+12
        rk=1
c 
        rlam0=0
        rlam1=rk-0.1d-14
        rlam1=8
  
  
  
  
  
        lenw=12 000 000
  
        call testit(nfuns,eps,y0,y1,rk,rlam0,rlam1,w,lenw)
c 
  
  
  
  
        stop
        end
  
c 
c 
c 
c 
c 
        subroutine testit(nfuns,eps,y0,y1,rk,rlam0,rlam1,w,lenw)
        implicit real *8 (a-h,o-z)
        save
        dimension ys7(100 000),ipivots(100 000),w(1),
     1      rlams(100 000),whts(100 000),xsout(1000)
c 
        external fun1
c 
c       initialize the function evaluator
c 
        call fun1in(nfuns,y0,y1,rk,ys7)
c 
c       construct the function approximator
c 
  
        call reprinit(ier,fun1,par1,par2,eps,ipivots,ncols,
     1      nfuns,rlam0,rlam1,rlams,whts,n,w,lenw,keep,lused,
     2      xsout)
c 
        call prinf('after reprinit, ipivots=*',ipivots,ncols)
        call prinf('after reprinit, keep=*',keep,1)
        call prinf('after reprinit, lused=*',lused,1)
  
cccc        call repretst2(fun1,par1,par2,ipivots,ncols,
cccc     1      nfuns,rlams,n,w,xsout)
  
        call repretst(fun1,par1,par2,ipivots,ncols,
     1      nfuns,rlams,n,w,xsout)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine repretst2(fun,par1,par2,ipivots,ncols,
     1      nfuns,xs,n,w,xsout)
        implicit real *8 (a-h,o-z)
        save
        dimension par1(1),par2(1),ipivots(1),xs(1),w(1),
     1      xsout(1),ftest(100 000),ftest2(100 000),ftest3(100 000),
     2      diffs(100 000),errs(100 000),rnorms(100 000),coefs(1000)
c 
        external fun
c 
c        test the accuracy of the obtained approximation
c        for all of the user-specified functions and all
c        of the nodes constructed by the subroutine reperegt0
c 
        ii=78
  
  
        do 3000 ii=1,nfuns
  
  
        call prinf('ii=*',ii,1)
  
        do 2200 i=1,n
c 
        call fun(xs(i),ii,par1,par2,ftest(i))
 2200 continue
c 
cccc        call prin2('ftest in repretst =*',ftest,ncols)
cccc        call prin2('xsout in repretst =*',xsout,ncols)
c 
        call reprcon1(ftest,w,coefs)
  
        call prin2('coefs2 in repretst =*',coefs,ncols)
c 
        d=0
        dd=0
  
        do 2400 i=1,n
c 
        call fun(xs(i),ii,par1,par2,ftest2(i))
c 
c 
        call repreval(fun,par1,par2,ipivots,coefs,ncols,
cccc     1      xs(i),ftest3(i))
     1      xs(i),ftest3(i))
c 
        diffs(i)=ftest3(i)-ftest2(i)
  
        d=d+diffs(i)**2
        dd=dd+ftest2(i)**2
  
 2400 continue
  
  
        err=sqrt(d/dd)
  
cccc        call prin2('and diffs=*',diffs,30)
  
  
        call prin2('and ftest2=*',ftest2,30)
        call prin2('and ftest3=*',ftest3,30)
        call prin2('and xs=*',xs,30)
        call prin2('and err=*',err,1)
  
        errs(ii)=err
        rnorms(ii)=sqrt(dd/n)
 3000 continue
  
  
        call prin2('errs=*',errs,nfuns)
        call prin2('and rnorms=*',rnorms,nfuns)
  
        return
        end
  
  
c 
c 
c 
c 
c 
        subroutine repretst(fun,par1,par2,ipivots,ncols,
     1      nfuns,xs,n,w,xsout)
        implicit real *8 (a-h,o-z)
        save
        dimension par1(1),par2(1),ipivots(1),xs(1),w(1),
     1      xsout(1),ftest(100 000),ftest2(100 000),
     2      diffs(100 000),errs(100 000),rnorms(100 000),
     3      ftest3(100 000),coefs(100 000)
c 
        external fun
c 
c        test the accuracy of the obtained approximation
c        for all of the user-specified functions and all
c        of the nodes constructed by the subroutine reperegt0
c 
        ii=78
  
  
        do 3000 ii=1,nfuns
  
  
        call prinf('ii=*',ii,1)
  
        do 2200 i=1,ncols
c 
        call fun(xsout(i),ii,par1,par2,ftest(i))
 2200 continue
c 
cccc        call prin2('ftest in repretst =*',ftest,ncols)
cccc        call prin2('xsout in repretst =*',xsout,ncols)
c 
        call reprcon2(ftest,w,coefs)
  
cccc        call prin2('coefs in repretst =*',coefs,ncols)
c 
        d=0
        dd=0
  
        do 2400 i=1,n
c 
        call fun(xs(i),ii,par1,par2,ftest2(i))
c 
c 
        call repreval(fun,par1,par2,ipivots,coefs,ncols,
     1      xs(i),ftest3(i))
c 
        diffs(i)=ftest3(i)-ftest2(i)
  
        d=d+diffs(i)**2
        dd=dd+ftest2(i)**2
  
 2400 continue
  
  
        err=sqrt(d/dd)
  
cccc        call prin2('and diffs=*',diffs,30)
  
  
cccc        call prin2('and ftest2=*',ftest2,30)
cccc        call prin2('and ftest3=*',ftest3,30)
cccc        call prin2('and xs=*',xs,30)
        call prin2('and err=*',err,1)
  
        errs(ii)=err
        rnorms(ii)=sqrt(dd/n)
 3000 continue
  
  
        call prin2('errs=*',errs,nfuns)
        call prin2('and rnorms=*',rnorms,nfuns)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine fun1(rlam,ii,par1,par2,f)
        implicit real *8 (a-h,o-z)
        save
        dimension ys(100 000),ys7(1),par1(1),par2(1)
c 
        complex *16 ima,cd,rk
c 
        data ima/(0.0d0,1.0d0)/
c 
c       evaluate the test function
c 
cccc        call prinf('in fun1, ii=*',ii,1)
cccc        call prin2('in fun1, rlam=*',rlam,1)
  
        if(ii .gt. nys) goto 1400
c 
        cd=cdexp(ima*ys(ii)*sqrt(rk**2-rlam**2))
        f=cd
  
cccc        call prin2('in fun1, f of first type =*',f,1)
  
        return
c 
 1400 continue
c 
        cd=cdexp(ima*ys(ii-nys)*sqrt(rk**2-rlam**2))
        f=-cd*ima
  
cccc        call prin2('in fun1, f of second type =*',f,1)
  
        return
c 
c 
c 
c 
        entry fun1in(nfuns7,y0,y1,rk7,ys7)
c 
c       initialize the array of values of y
c 
        delta=1.0d-13
        rk=rk7+delta*ima
        nfuns=nfuns7
        nys=nfuns/2
        hys=(y1-y0)/(nys-1)
  
        call prinf('in fun1in, nfuns7=*',nfuns7,1)
        call prinf('in fun1in, nys=*',nys,1)
        call prin2('in fun1in, y0=*',y0,1)
        call prin2('in fun1in, y1=*',y1,1)
        call prin2('in fun1in, rk=*',rk,2)
  
        call prin2('in fun1in, hys=*',hys,1)
  
        do 2200 i=1,nys
c 
ccccc        call prinf('in fun1in, i=*',i,1)
  
        ys(i)=(i-1)*hys + y0
        ys7(i)=ys(i)
 2200 continue
c 
        call prin2('in fun1in, ys=*',ys,nys)
c 
        return
        end
c 
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c       This is the end of the debugging code and the ebginning of the
c       expansion code proper
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c        This file contains a set of subroutines that, for a
c        user-specified collection of functions
c 
c        $f_1, f_2, ... f_{nfuns} : [a,b] \to R^1,                      (1)
c 
c        choose a subset of functions
c 
c        $\phi_1, \phi_2, ... \phi_{ncols} : [a,b] \to R^1,             (2)
c 
c        such that (hopefully) ncols is not too large and that the
c        functions (2) span (to the precision eps) all of the
c        functions (1). Then, a set of nodes
c 
c        $xout_1, xout_2, . . . , xout_{ncols}$                         (2a)
c 
c        is chosen on the interval [a,b] such that the
c        $ncols \times \ncols$-matrix $a$ with the coefficients
c 
c        a_{i,j} = \phi_i{x_j}                                          (3)
c 
c        has a generalized inverse that is stable in the appropriate
c        sense. Subsequently, any linear combination of the form
c 
c        F(x) = \sum_{i=1}^{nfuns} \alpha_i * f_i(x)                    (4)
c 
c        can be efficiently expressed in the form
c 
c        F(x) \sim \sum_{i=1}^{ncols} \coefs_i * \phi_i(x);             (5)
c 
c        the only data required by the procedure are the values of the
c        function F at the nodes $xout_1, xout_2, . . . ,xout_{ncols}$.
c 
c        The user-callable subroutines in this file are: reprinit,
c        reprcon1, reprcon2, repreval. Following is a brief description
c        of all four of them things.
c 
c  reprinit - the principal subroutine in this package. For the
c        user-specified collection of functions (1), it chooses a
c        of functions (2) such that (hopefully) ncols is not too
c        large and that the functions (2) span (to the precision eps)
c        all of the functions (1). Then, it chooses a set of nodes
c        (2a) on the interval [a,b] such that the
c        $ncols \times \ncols$-matrix $a$ with the coefficients (3)
c        has a generalized inverse that is stable in the appropriate
c        sense. It also constructs all of the machinery necessary for
c        expanding linear combinations of the form (4) in the form (5),
c        and returns the obtained data in the array w, to be subsequently
c        used by the subroutines reprcon1, reprcon2.
c 
c  reprcon2 - uses user-supplied values of a function of form (4) at the
c        nodes (2a) to return to the user the expansion of F in the form
c        (5). The latter can be subsequently evaluated at arbitrary nodes
c        on the interval [a,b] via the subroutine repreval.
c 
c  reprcon1 - uses user-supplied values of a function of form (4) at the
c        nodes xs to return to the user the expansion of F in the form
c        (5). The latter can be subsequently evaluated at arbitrary nodes
c        on the interval [a,b] via the subroutine repreval. PLEASE NOTE THAT
C        THIS SUBROUTINE PRODUCES THE SAME RESULT AS THE SUBROUTINE REPRCON2,
C        BUT REQUIRES THE ARRAY OF VALUES OF THE FUNCTION (4) AT N POINTS
C        XS; NORMALLY, N IS MUCH GREATER THAN NCOLS, AND THIS SUBROUTINE IS
C        NOT EXPECTED TO SEE MUCH USE. ON THE OTHER HAND, IT PRODUCES
C        SLIGHTLY MORE ACCURATE RESULTS THAN THOSE PRODUCED BY REPRCON2.
C        SOMETIMES,THE IMPROVEMENT MIGHT BE WORTH THE TROUBLE (BUT DO NOT
C        HOLD YOUR BREATH).
c 
c  repreval - evaluates an expansion produced by one of the subroutines
c        reprcon1, reprcon2 at an arbitrary points on the interval [a,b]
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        subroutine repreval(fun,par1,par2,ipivots,coefs,ncols,
     1      x,fout)
        implicit real *8 (a-h,o-z)
        save
        dimension par1(1),par2(1),ipivots(1),coefs(1)
c 
        external fun
c 
c       This subroutine evaluates a user-supplied "general" expansion
c       at the user-supplied point. The expansion has the form
c 
c        F(x) \sim \sum_{i=1}^{ncols} \coefs(i) * \phi_{ipivots(i)}(x);   (1)
c 
c       It uses the array ipivots, produced by the subroutine reprinit
c       (see) and the array coefs, normally produced by either the
c       subroutine reprcon2 (see), or the subroutine reprcon1 (see).
c       This subroutine has no known uses as a stand-alone device.
c 
c                  Input parameters:
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c  fun - the subroutine evaluationg the functions to be compressed.
c       The calling sequence for fun is
c 
c            call fun(x,i,par1,par2,fout),
c 
c                 Input parameters of the subroutine fun:
c 
c  x - the argument (must be on the interval [a,b])
c  i - the sequence number of the function to be evaluated (between
c      1 and nfuns)
c  par1, par2 - whatever parameters the subroutine fun needs
c 
c                 Output parameters of the subroutine fun:
c 
c  fout - the value of the i-th function at the point x
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c  par1, par2 - whatever the user-supplied subroutine fun needs
c  eps - the accuracy to which the functions are to be approximated
c  ipivots - the list of sequence numbers of input functions (defined
c        by the subroutine fun) that are used to approximate all other
c        functions in the family. There are ncols eleents returned in
c        this array; however, before this subroutine is called, the
c        value of ncols is not known, so ipivots should contain at least
c        nfuns elements in it
c  coefs - the coefficients in expansion (1) above
c  ncols - the number of elements in the sum (1) above; also the
c        number of elements in arrays ipivots, coefs.
c  x - the point at which the expansion (1) above is to be evaluated
c 
c                     Output parameters:
c 
c  fout - the value of the expansion (1) evaluated at the input point x
c 
c       . . . evaluate the "general function" expansion
c 
        fout=0
        do 1200 i=1,ncols
c 
        j=ipivots(i)
        call fun(x,j,par1,par2,d)
c 
        fout=fout+coefs(i)*d
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine reprcon2(f,w,coefs)
        implicit real *8 (a-h,o-z)
        save
        dimension coefs(1),f(1),w(1)
c 
c       This subroutine produces the coefficients of the "general"
c       expansion of the user-specified function f; the function is
c       specified by its values at the nodes xsout, returned by the
c       subroutine reprinit (see). It also uses the array w, also
c       produced by the subroutine reprinit. Please note that this
c       is supposed to be the principal means for the construction
c       of the array coefs, the secondary means being the subroutine
c       reprcon1 (see). This subroutine requires the values of the
c       function to be expanded tabulated at ncols nodes, returned
c       in the array xsout by the subroutine reprinit; the subroutine
c       reprcon1 requires the values of the function to be expanded
c       tabulated at the n nodes xs, also returned by the subroutine
c       reprinit, and n tends to be muuuch greater than ncols. On the
c       other hand, reprcon1 produces marginally more accurate results.
c 
c                      Input parameters:
c 
c  f - the function to be decomposed, tabulated at the nodes xsout
c       (ncols of them things, returned by the subroutine reprinit)
c  w - the array of data returned by the same subroutine reprinit
c 
c                      Output parameters:
c 
c  coefs - the ncols coefficients of the "general" expansion of the
c       function f
c 
c       . . . construct the memory map for the evaluation of the
c             coefficients of the decomposition of the user-supplied
c             function into the basis elements
c 
        iu=w(1)
        iv=w(2)
        is=w(3)
        iw1=w(4)
        iuuadj=w(7)
c 
        ncols=w(5)
        n=w(6)
c 
c       evaluate the coefficients of the decomposition of the
c       user-supplied function into the basis elements
c 
        call reprcon4(w(iv),w(is),ncols,f,
     1      coefs,w(iw1),w(iuuadj))
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine reprcon4(vv,s,ncols,f,coefs,w1,uuadj)
        implicit real *8 (a-h,o-z)
        save
        dimension vv(1),s(1),coefs(1),w1(1),f(1),
     1      uuadj(ncols,ncols)
c 
c       evaluate the coefficients of the decomposition of the
c       user-supplied function into the basis elements
c 
        call repmatve(uuadj,ncols,ncols,f,w1)
c 
        do 1600 j=1,ncols
c 
        w1(j)=w1(j)/s(j)
 1600 continue
c 
        call repmatve(vv,ncols,ncols,w1,coefs)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine reprcon1(f,w,coefs)
        implicit real *8 (a-h,o-z)
        save
        dimension coefs(1),f(1),w(1)
c 
c       This subroutine uses the data supplied in the array w and
c       (hopefully) created via a prior call to the subroutine
c       reprinit (see) to construct the decomposition of the function
c       whose values (n of them) at the nodes xs (returned by the
c       subroutine reprinit) are supplied to this subroutine in the
c       array f.
c 
c                           Input parameters:
c 
c  f - the values of the function to be approximated (n of them things)
c       tabulated at the nodes xs (returned by the subroutine reprinit)
c  w - the array returned by the subroutine reprinit (see)
c 
c                           Output parameters:
c 
c  coefs - coefficients of the expansion (ncols of them things)
c 
c 
c       . . . construct the memory map for the evaluation of the
c             coefficients of the decomposition of the user-supplied
c             function into the basis elements
c 
        iu=w(1)
        iv=w(2)
        is=w(3)
        iw1=w(4)
c 
        ncols=w(5)
        n=w(6)
c 
c       evaluate the coefficients of the decomposition of the
c       user-supplied function into the basis elements
c 
        call reprcon3(w(iu),w(iv),w(is),n,ncols,f,coefs,w(iw1))
  
  
        return
        end
c 
c 
c 
c 
c 
        subroutine reprcon3(uu,vv,s,n,ncols,f,coefs,w1)
        implicit real *8 (a-h,o-z)
        save
        dimension uu(1),vv(1),s(1),coefs(1),w1(1),f(1)
c 
c       evaluate the coefficients of the decomposition of the
c       user-supplied function into the basis elements
c 
        call repmatav(uu,ncols,n,f,w1)
  
        do 1600 j=1,ncols
c 
        w1(j)=w1(j)/s(j)
 1600 continue
c 
        call repmatve(vv,ncols,ncols,w1,coefs)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine reprinit(ier,fun,par1,par2,eps,ipivots,ncols,
     1      nfuns,a,b,xs,whts,n,w,lenw,keep,lused,xsout)
        implicit real *8 (a-h,o-z)
        save
        dimension par1(1),par2(1),ipivots(1),xs(1),whts(1),w(1),
     1      xsout(1)
c 
        external fun
c 
c       For nfuns user-specified functions on the
c       interval [a,b], this subroutine selects a subset of these
c       functions such that all nfuns of them things are linear
c       combinations of the selected subset, to the precision eps.
c       The required number of such functions is returned to the
c       user in the parameter ncols, and the sequence numbers of
c       the selected functions is returned in the array ipivots;
c       needless to say, ncols might turn out to be equal to nfuns,
c       if all nfuns input functions are strongly  linearly independent,
c       or if eps is too small. The subroutine also returns in the
c       array w various types of data to be used by the subroutine
c       coefsget (see) to evaluate the coefficients of expansions
c       of various user-supplied functions into functions defined by
c       the subroutine fun and the indices ipivots(1), ipivots(2), ...
c       ipivots(ncols). PLEASE NOTE THAT THIS SUBROUTINE HAS NO
c       USES AS A STAND-ALONE DEVICE!!
c 
c 
c                       Input parameters:
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c  fun - the subroutine evaluationg the functions to be compressed.
c       The calling sequence for fun is
c 
c            call fun(x,i,par1,par2,fout),
c 
c                 Input parameters of the subroutine fun:
c 
c  x - the argument (must be on the interval [a,b])
c  i - the sequence number of the function to be evaluated (between
c      1 and nfuns)
c  par1, par2 - whatever parameters the subroutine fun needs
c 
c                 Output parameters of the subroutine fun:
c 
c  fout - the value of the i-th function at the point x
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c  par1, par2 - whatever the user-supplied subroutine fun needs
c  eps - the accuracy to which the functions are to be approximated
c         REMARK: Please note that the subroutine tends to interpret
c                 this parameter conservatively, so that the obtained
c                 precision is usually better than eps (unless the
c                 machine zero gets in the way); specifying eps that is
c                 too small tends to lead to various unpleasant disasters.
c 
c 
c 
c  nfuns - the number of functions to be approximated
c  a,b - the interval on which the functions are to be approximated
c  lenw - the size (in real *8 elements) of the work array w
c 
c                       Output parameters:
c 
c  ier - error return code
c  ipivots - the list of sequence numbers of input functions (defined
c        by the subroutine fun) that are used to approximate all other
c        functions in the family. There are ncols eleents returned in
c        this array; however, before this subroutine is called, the
c        value of ncols is not known, so ipivots should contain at least
c        nfuns elements in it
c  ncols - the number of elements returned in ipivots
c  xs - the nodes on the interval [a,b] used by this subroutine to
c        tabulate the functions being analyzed; please note that these
c        nodes will be needed by the user to construct the data to be used
c        by the subroutine coefsget (see); it is hard to predict the
c        required size of this array, so make it biggish
c  whts - the weights corresponding to the nodes in the array xs
c  n - the number of elements in each of the arrays xs, whts
c  w - the first keep elements of this array contain the information to
c        be used by the subroutine coefsget (see) to construct expansions
c        of user-supplied functions.
c  keep - the number of elements in array w that has to be unchanged
c        between the call to this subroutine and subsequent call(s) to
c        subroutine coefsget (see)
c  lused - the maximum number of elements in array w used by this
c        subroutine at any time
c 
c 
c        . . . construct the subdivision of the interval that will
c              representall of the user-specified functions with
c              precision eps
c 
        k=20
        ier=0
c 
        iw=1
        lw=10*k+8*k**2+30
        lw=lw*2+20 000
c 
        iab=iw+lw
        lab=lenw-iw-10
c 
        call manincon(jer,a,b,k,fun,nfuns,par1,par2,
     1    eps,w(iab),lab,nn,w(iw),lw)
c 
        if(jer .ne. 0) ier=16
        if(jer .ne. 0) return
c 
        call prinf('after manincon, jer=*',jer,1)
        call prinf('after manincon, nn=*',nn,1)
c 
        call mvecmove(w(iab),w(iw),nn*2)
        iab=iw
        iw=iab+nn*2+4
        lwleft=lenw-iw-4
c 
c       allocate memory for the construction of the low-rank
c       representation of all of the user-supplied functions
c 
        n=nn*k
        ifs=iw
c 
        lfs=n*nfuns+10
c 
        irnorms=ifs+lfs
        lrnorms=nfuns+2
c 
        iuu=irnorms+lrnorms
        luu=n*nfuns+10
c 
        ivv=iuu+luu
        lvv=n*nfuns+10
c 
        is=ivv+lvv
        ls=nfuns+2
c 
        iww=is+ls
        lwleft=lenw-iww
c 
        call prinf('iww=*',iww,1)
c 
c       construct the low-rank representation of all of the
c       user-supplied functions
c 
        call reprini0(ier,xs,whts,fun,par1,par2,w(ifs),n,nfuns,
     1      nn,k,eps,ipivots,ncols,w(iab),w(iuu),w(ivv),
     2      w(is),w(irnorms),w(iww),lwleft,ltot)
c 
        lused=iww+ltot+5
c 
c       compress to obtained data, and store them in the array w for
c       use by the subroutine reprcon3 (see)
c 
        iiuu=10
        lluu=ncols*n+10
c 
        iivv=iiuu+lluu
        llvv=ncols**2+10
c 
        iiss=iivv+llvv
        llss=ncols+2
c 
        call mvecmove(w(iuu),w(iiuu),n*ncols)
        call mvecmove(w(ivv),w(iivv),ncols**2)
        call mvecmove(w(is),w(iiss),ncols)
c 
        w(1)=iiuu+0.1
        w(2)=iivv+0.1
        w(3)=iiss+0.1
c 
        iw1=iiss+ncols+2
        lw1=n+2
        w(4)=iw1+0.1
        keep=iw1+lw1
c 
        w(5)=ncols+0.1
        w(6)=n+0.1
c 
c       allocate memory for the construction of the representation
c       of the user's functions that are com,pressed from both sides
c 
        iipivots=keep+1
        lipivots=n+2
c 
        irnorms=iipivots+lipivots
        lrnorms=n+2
c 
        iww=irnorms+lrnorms
        lww=ncols**2+3*ncols+100
c 
        iuuadj=iww+lww
        luuadj=ncols*ncols+2
c 
        call reprpnts(ier,w(iiuu),n,ncols,w(iuuadj),xs,whts,
     1      xsout,w(iipivots),w(irnorms),w(iww),rcond)
c 
c        perform the final garbage collection
c 
        iuuadj2=keep+1
        luuadj2=ncols**2
c 
        call mvecmove(w(iuuadj),w(iuuadj2),ncols**2)
c 
        iuuadj=iuuadj2
        w(7)=iuuadj+0.1
        keep=iuuadj2+luuadj2
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine reprpnts(ier,uu,n,ncols,uuadj,xs,whts,
     1      xsout,ipivots,rnorms,w,rcond)
        implicit real *8 (a-h,o-z)
        save
        dimension uu(n,ncols),uuadj(ncols,n),rnorms(1),
     1      ipivots(1),whts(1),xs(1),w(1),xsout(1)
c 
c       apply the pivoted Gram-Schmidt to the adjoint of uu
c 
        do 1400 i=1,ncols
        do 1200 j=1,n
c 
        uuadj(i,j)=uu(j,i)/sqrt(whts(j))
 1200 continue
 1400 continue
c 
        eps2=1.0d-30
        call reprpiv(uuadj,ncols,n,rnorms,eps2,ncols3,ipivots)
c 
c       sort the array ipivots
c 
        call reprbubb(ipivots,ncols3)
c 
        call prinf('after reprbubb in reprpnts, ipivots=*',
     1      ipivots,ncols3)
c 
c       put the ncols columns of the original matrix (b) whose
c       numbers have been marked by the subroutine reprpiv in the
c       first ncols columns of the matrix b
c 
        do 1800 i=1,ncols
c 
        ii=ipivots(i)
        do 1600 j=1,ncols
c 
        uuadj(i,j)=uu(ii,j)/whts(ii)
 1600 continue
c 
        xsout(i)=xs(ii)
 1800 continue
c 
c       invert the matrix uuadj
c 
        call orthom(uuadj,ncols,w,cond)
c 
        call prin2('after orthom, cond=*',cond,1)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine reprini0(ier,xs,whts,fun,par1,par2,fs,n,nfuns,
     1      nn,k,eps,ipivots,ncols,ab,uu,vv,s,rnorms,w,lw,ltot)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension w(1),s(1),uu(1),vv(1),xs(n),whts(n),rnorms(1),
     1      ipivots(1),fs(n,1),ab(2,1)
c 
        external fun
c 
c        construct the nodes on the interval [a,b] corresponding to the
c        subdivision of the interval produced by manincon
c 
        ier=0
c 
        call nodewhts(ab,nn,k,xs,whts)
c 
c       construct the values of the functions
c 
        n=nn*k
        do 1400 i=1,nfuns
        do 1200 j=1,n
c 
        call fun(xs(j),i,par1,par2,fs(j,i))
c 
        fs(j,i)=fs(j,i)*sqrt(whts(j))
 1200 continue
 1400 continue
c 
c       construct the low-rank representation of all of the
c       user-supplied functions
c 
        call reprfuns(ier,fs,n,nfuns,eps,ncols,uu,vv,s,w,lw,
     1      rnorms,ipivots,fun,par1,par2,xs,whts,ltot)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine reprfuns(ier,b,n,m,eps,ncols,uu,vv,s,w,lw,
     1      rnorms,ipivots,fun,par1,par2,xs,whts,ltot)
        implicit real *8 (a-h,o-z)
        save
        dimension b(n,m),rnorms(1),ipivots(1),whts(1),
     1      w(1),s(1),uu(n,1),vv(m,1),xs(1)
c 
        external fun
c 
c       apply the pivoted Gram-Schmidt to b
c 
        ier=0
c 
        call reprpiv(b,n,m,rnorms,eps,ncols,ipivots)
c 
c       sort the array ipivots
c 
        call reprbubb(ipivots,ncols)
  
        call prinf('after reprbubb, ipivots=*',ipivots,ncols)
c 
c       put the ncols columns of the original matrix (b) whose
c       numbers have been marked by the subroutine reprpiv in the
c       first ncols columns of the matrix b
c 
        do 1800 i=1,ncols
c 
        ii=ipivots(i)
        do 1600 j=1,n
c 
        call fun(xs(j),ii,par1,par2,b(j,i))
c 
        b(j,i)=b(j,i)*sqrt(whts(j))
 1600 continue
c 
 1800 continue
c 
c       now, SVD the first ncols of the reordered input matrix
c 
        eps2=eps/1000
        call svdpivot(jer,b,n,ncols,uu,vv,s,ncols2,eps2,
     1      w,lw,ltot)
c 
        if(jer .ne. 0) then
            ier=8
            return
        endif
c 
        call prinf('after svdpivot, ier=*',ier,1)
        call prinf('after svdpivot, ncols2=*',ncols2,1)
  
        call prin2('after svdpivot, s=*',s,ncols2)
  
        do 2800 i=1,ncols
        do 2600 j=1,n
c 
        uu(j,i)=uu(j,i)*sqrt(whts(j))
 2600 continue
c 
 2800 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine reprbubb(ipivots,ncols)
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
        subroutine repmatav(a,n,m,x,y)
        implicit real *8 (a-h,o-z)
        save
        dimension a(m,n),x(m),y(n)
c 
        do 1400 i=1,n
        d=0
        do 1200 j=1,m
        d=d+a(j,i)*x(j)
 1200 continue
        y(i)=d
 1400 continue
	return
	end
c 
c 
c 
c 
c 
        subroutine repmatve(a,n,m,x,y)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),x(m),y(n)
c 
        do 1400 i=1,n
        d=0
        do 1200 j=1,m
        d=d+a(i,j)*x(j)
 1200 continue
        y(i)=d
 1400 continue
	return
	end
c 
c 
c 
c 
c 
        subroutine reprpiv(b,n,m,rnorms,eps,ncols,ipivots)
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
        call svdscapr(b(1,i),b(1,j),n,d)
c 
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*d
 2770 continue
 2780 continue
 2790 continue
c 
c       normalize the i-th column
c 
        call svdscapr(b(1,i),b(1,i),n,d)
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
        call svdscapr(b(1,i),b(1,j),n,d)
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
