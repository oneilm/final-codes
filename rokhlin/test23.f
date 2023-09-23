        implicit real *8 (a-h,o-z)
        dimension xsall(32,30),wsall(32,30),
     1      xsout(100),wsout(100),bigmatr(30,32,30),
     2      xsall2(32,30),wsall2(32,30),xs2(100),ws2(100),
     3      fs2(100),fslrg(32,30),fs3(100),diffs(100)
c 
        external funtes
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER npols '
        READ *,npols
        CALL PRINF('npols=*',npols,1 )
c 
c        retrieve from the subroutine the nodes and weights
c        for all 30 quadratures
c 
        call univquadrs(xsall,wsall)
  
        call prin2('after univquadrs, xsall=*',xsall,32)
        call prin2('after univquadrs, wsall=*',wsall,32)
c 
c       retrieve the support nodes and weights
c 
        call xswsget(xsout,wsout,npts)
  
        call prin2('after xswsget, xsout=*',xsout,30)
        call prin2('after xswsget, wsout=*',wsout,30)
c 
c       integrate the test function via the retrieved quadrature
c 
        inode=15
        inode=3
cccc        inode=11
  
        x0=xsout(inode)
c 
        rint=0
        do 1400 i=1,32
c 
        call funtest(x0,xsall(i,inode),f)
c 
        rint=rint+wsall(i,inode)*f
 1400 continue
  
        call prin2('and rint=*',rint,1)
c 
c       evaluate the same integral via adapgaus
c 
        a=-1
        b=1
        m=24
        eps=1.0d-15
  
  
        call adapgaus(ier,a,b,funtes,x0,par2,m,eps,
     1      rint2,maxrec,numint)
  
        call prin2('and rint2=*',rint2,1)
  
        call prin2('and rint2-rint=*',rint2-rint,1)
c 
c       construct and test the matrix of all interpolation
c       coefficients, interpolation from the 30 support nodes
c       to all auxiliary nodes (all 32*30 of them things!)
c 
        call univintall(bigmatr,xs2,ws2,xsall2,wsall2)
  
cccc        call prin2('after univintall, bigmatr=*',bigmatr,30*30*32)
  
        call prin2('after univintall, xs2=*',xs2,30)
        call prin2('after univintall, ws2=*',ws2,30)
c 
c       test the large interpolation matrix bigmatr by applying it to
c       the values of the test function at the support nodes, and
c       comparing the result with the values of the test function at
c       the auxiliary nodes
c 
        do 2200 i=1,30
  
        call funtest2(xs2(i),fs2(i))
 2200 continue
c 
        call prin2('fs2 as created*',fs2,30)
  
  
  
  
        call univallapp(bigmatr,fs2,fslrg)
  
  
        call prin2('first column of fslrg is*',fslrg,32)
        call prin2('and first column of xsall is*',xsall,32)
  
  
  
c 
c        conduct massive testing
c 
       dd=0
       ddd=0
       do 2800 i=1,30
c 
c       . . . evaluate directly the test function at the auxiliary
c             nodes corresponding to the i-th support node
c 
        do 2400 j=1,32
c 
        call funtest2(xsall(j,i),fs3(j))
c 
        diffs(j)=fslrg(j,i)-fs3(j)
c 
        dd=dd+fs3(j)**2
        ddd=ddd+diffs(j)**2
 2400 continue
c 
        call prinf('i=*',i,1)
  
        call prin2('and i-th column directly*',fs3,32)
  
        call prin2('and i-th column via bigmatr*',fslrg(1,i),32)
        call prin2('and the differences are *',diffs,32)
  
  
 2800 continue
  
        dddd=sqrt(ddd/dd)
  
        call prin2('and relative l2-error is*',dddd,1)
  
        stop
        end
  
c 
c 
c 
c 
c 
        subroutine funtest2(x,f)
        implicit real *8 (a-h,o-z)
c 
        save
  
  
        n=7
        delta=1.0d-16
        call legepol(x,n,pol,der)
  
        f=pol*(1+x)*log(1+x)
        f=pol*(1-x)*log(1-x)
cccc        f=pol
        return
        end
c 
c 
c 
c 
c 
        function funtes(x,x0,par2)
        implicit real *8 (a-h,o-z)
        save
  
        call funtest(x0,x,f)
        funtes=f
  
        return
        end
c 
c 
c 
c 
c 
        subroutine funtest(x0,x,f)
        implicit real *8 (a-h,o-z)
c 
        save
  
  
        n=11
        delta=1.0d-16
        call legepol(x,n,pol,der)
  
        f=log( (x-x0)**2+delta**2) * pol
        return
        end
  
c 
c 
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the testing code and the beginning of the
c        interpolation code proper
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c        The user-callable subroutines in this file are: univintall,
c        univinterp, ufuneval, matrret3, xswsget, uoutget, voutget,
c        vmatrget. Actually, the last four are entries in the subroutine
c        matrret3; each of them returns one of the parameters that are
c        all returned by matrret3 at once.
c 
c   univintall - returns to the user a set of 30 matrices.
c       Each of these matrices converts the values of a function at
c       the support nodes (30 of them things) into the values of that
c       function at the quadrature nodes corresponding to one support
c       node. Since there are 20 support nodes, the number of matrices
c       returned by this subroutine is also 20. On the other hand, there
c       is an anomaly associated with the matrices returned by this
c       subroutine:
c 
c   univinterp -   returns to the user the coefficients coefs
c        of the interpolation formula expressing the
c        values of a function at the user-specified point x
c        on the interval [-1,1] via that function's values at
c        the 30 support nodes.
c 
c        The function to be interpolated is assumed to be a linear
c        combination of functions
c 
c        P_n(x),
c 
c        P_n (x) \cdot (1-x) \cdot log(1-x),
c                                                                (1)
c        P_n (x) \cdot  \cdot(1+x) log(1+x),
c 
c        with n=0,1,2,...,9.
c 
c        The interpolation is exact to double precision.
c 
c   ufuneval - evaluates at the user-supplied point x on the
c       interval [-1,1] the 30 elements of the "standard"
c       orthonormal basis in the space spanned by the functions (1)
c 
c   matrret3 - returns to the user several pieces of data
c        to be used for the interpolation and integration on the
c        interval [-1,1] of functions of the form. Specifically,
c        returned are the 30 nodes xs and corresponding quadrature
c        weights ws, the matrix uout of values of elements of the
c        "standard" basis at the nodes xs, the inverse vout of uout,
c        and the matrix vmatrout that is used later (in the subroutine
c        ufuneval) for the evaluation of the 30 elements of the basis
c        at arbitrary points on the interval [-1,1]
c 
c   xswsget, uoutget, voutget, vmatrget - each of these entries returns
c        a part of the data returned by the principal entry matrret3;
c        which entry returns which part of data is obvious from the names
c        of parameters.
c 
c   univquadrs - returns to the user 32-point special-purpose
c        quadratures (30 of them things). The quadrature number k
c        (k=1,2,...,30) integrates functions of the form
c 
c        f(x)=phi(x)+psi(x)*log(|x-x_k|),                           (1)
c 
c        with phi and psi polynomials of order 19, and x_k one of the
c        support nodes returned by the subroutine matrret3 (or its entry
c        xswsget).
c 
c 
c 
c 
        subroutine univintall(bigmatr,xs,ws,xsall,wsall)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(31),ws(31),vmatr(30,30),v(30,30),
     1      xsall(32,30),wsall(32,30),bigmatr(30,32,30)
c 
c       This subroutine returns to the user a set of 30 matrices.
c       Each of these matrices converts the values of a function at
c       the support nodes (30 of them things) into the values of that
c       function at the quadrature nodes corresponding to one support
c       node. Since there are 20 support nodes, the number of matrices
c       returned by this subroutine is also 20. On the other hand, there
c       is an anomaly associated with the matrices returned by this
c       subroutine:
c 
c   1. The matrices returned by this subroutine are actually adjoints
c       of the matrices converting the values of a function at the
c       support nodes into its values at the quadrature nodes
c 
c       More specifically, the value f_{k,i} of the function f at the
c       i-th quadrature node corresponding to the k-th support node is
c       given by the formula
c 
c            f_{k,i}=\sum_{j=1}^{30} bigmatr(j,i,k)*f(j),             (1)
c 
c       with k=1,2,...,30, and i=1,2,...,32.
c 
c       In addition, the subroutine returns a number of other objects
c       that might be useful in the environment for which it has been
c       designed: support nodes xs and corresponding weights ws, the
c       quadrature nodes and weights corresponding to each of the
c       support nodes (all of these are returned in the arrays xsall,
c       wsall, respectively).
c 
c 
c                 Input parameters: None
c 
c                 Output parameters:
c 
c  bigmatr - a 30 \times 32 \times 30 - array, containing the set of 30
c       matrices described above.
c  xs - the 30 support nodes on the interval [-1,1]
c  ws - quadrature weights (pretty universal) corresponding to the
c       nodes xs.
c  xsall - the array containing the quadrature nodes corresponding to
c       all support nodes; the quadrature nodes corresponding to the
c       i-th support nodes are returned in the i-th column of the array
c       points. The number of quadrature nodes for each support node
c       is 32.
c  wsall - the array containing the quadrature weights corresponding to
c       the nodes return in the array points. One weight correponds to
c       each of the nodes in the array points, so that the structures of
c       these two arrays are identical.
c 
c        . . . retrieve the support nodes and corresponding weights,
c              the matrix vmatr used in the evaluation of the basis
c              functions, and the matrix v, converting values of a
c              function at the support nodes into the coefficints
c              of its expansion
c 
        call xswsget(xs,ws,npts)
c 
        call vmatrget(vmatr)
c 
        call voutget(v)
c 
c       retrieve all auxiliary quadrature nodes for all support nodes
c 
        call univquadrs(xsall,wsall)
c 
c       construct the matrices (all 32 of them things) converting
c       the values of a function at the 30 support nodes into its
c       values at all quadrature nodes corresponding to all support
c       nodes
c 
        do 2000 i=1,30
c 
c       for the i-th support node xs(i), construct the interpolation
c       formula for each of the corresponding auxiliary quadrature
c       nodes
c 
        do 1600 j=1,32
c 
        xx=xsall(j,i)
        call univinterp(xx,bigmatr(1,j,i))
c 
 1400 continue
 1600 continue
 2000 continue
        return
        end
  
c 
c 
c 
c 
c 
        subroutine univallapp(coefslrg,fs,fslrg)
        implicit real *8 (a-h,o-z)
        save
        dimension coefslrg(30,32,30),fs(30),
     1      fslrg(32,30)
c 
c       This subroutine applies to the user-specified vector fs of
c       length 30 (representing the values at the 30 support nodes
c       of an appropriately-behaved function f) the 30 interpolation
c       matrices produced by the subroutine univintall (see). The result
c       is the collection fslrg of (approximate) values  of the function
c       f at the nodes xsall (also produced by the subroutine univintall).
c       This subroutine is principally a tool for the debugging of the
c       subroutine univintall, since I am not aware of any other obvious
c       uses for this subroutine.
c 
c                          Input parameters:
c 
c  coefslrg - a 30 \times 32 \times 32 - array, containing the set of 30
c       matrices described above.
c  fs - the values of the function f to be interpolated at the 30 support
c       nodes on the interval [-1,1]
c 
c                          Output parameters:
c 
c  fslrg - the interpolated values of the function f at the quadrature
c       nodes xsall (the latter are returned by univintall).
c 
c        . . . given a user-supplied table of function values at the
c              support nodes, evaluate that function at the nodes of all
c              quadratures for all support nodes
c 
        do 2000 k=1,30
c 
        do 1800 i=1,32
        d=0
        do 1600 j=1,30
c 
        d=d+coefslrg(j,i,k)*fs(j)
 1600 continue
c 
        fslrg(i,k)=d
 1800 continue
 2000 continue
c 
        return
        end
c 
c 
        subroutine univinterp(x,coefs)
        implicit real *8 (a-h,o-z)
        save
        dimension usout(32),vmatr(901),coefs(1),v(30,30)
c 
        data ifcalled/0/
c 
c       This subroutine returns to the user the coefficients
c       forminte of the interpolation formula expressing the
c       values of a function at the user-specified point x
c       on the interval [-1,1] via that function's values at
c       the 30 support nodes.
c 
c        The function to be interpolated is assumed to be a linear
c        combination of functions
c 
c        P_n(x),
c 
c        P_n (x) \cdot (1-x) \cdot log(1-x),
c                                                                (1)
c        P_n (x) \cdot  \cdot(1+x) log(1+x),
c 
c        with n=0,1,2,...,9.
c 
c        The interpolation is exact to double precision.
c 
c                   Input parameters:
c 
c  x - the point on the interval [-1,1] at which the function is
c        to be interpolated
c 
c                   Output parameters:
c 
c  coefs - the coefficients of the interpolation formula
c        (30 of them things)
c 
c       . . . Evaluate all basis functions at the point x
c 
         call ufuneval(x,usout)
c 
c        if this is the first call to this subroutine, retrieve
c        the matrix v from the subroutine voutget
c 
         if(ifcalled .eq. 0) then
             call voutget(v)
             ifcalled=0
         endif
c 
c        apply the adjoint of the matrix v to the vector
c        of values of basis functions at the user-specified
c        point x, obtaining the interpolation formula
c 
        n=30
        do 1400 i=1,n
c 
        d=0
        do 1200 j=1,n
c 
        d=d+v(j,i)*usout(j)
 1200 continue
c 
        coefs(i)=d
 1400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine ufuneval(x,usout)
        implicit real *8 (a-h,o-z)
        save
        dimension vmatr(30,30),usout(1),fs(32)
c 
        data ifcalled/0/
c 
c       This subroutine evaluates at the user-supplied point x
c       on the interval [-1,1] the 30 elements of the "standard"
c       orthonormal basis in the space spanned by the functions
c 
c        P_n(x),
c 
c        P_n (x) \cdot (1-x) \cdot log(1-x),
c                                                                    (1)
c        P_n (x) \cdot  \cdot(1+x) log(1+x),
c 
c        with n=0,1,2,...,9.
c 
c                   Input parameters:
c 
c  x - the poiint on the interval [-1,1] at which the functions
c        are to be evaluated
c 
c                   Output parameters:
c 
c  usout - the 30 basis functions evaluated at the point x
c 
c       . . . if this is the first call to this subroutine - retrieve
c             the matrix vmatr from the subroutine vmatrget
c 
        if(ifcalled .eq. 0) then
            call vmatrget(vmatr)
            ifcalled=1
        endif
c 
c       evaluate all of the original functions at the user-specified
c       point x
c 
        nfuns=30
        ncols=30
c 
        call funs1(x,fs)
c 
c       apply the matrix vmatr to the obtained values
c 
        do 1600 i=1,ncols
c 
        d=0
        do 1400 j=1,nfuns
c 
        d=d+vmatr(j,i)*fs(j)
 1400 continue
c 
        usout(i)=d
 1600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine funs1(x,pols)
        implicit real *8 (a-h,o-z)
        save
        dimension pols(1)
c 
        n=10
        call legepols(x,n,pols)
c 
        dp=log(1+x)*(1+x)
        dm=log(1-x)*(1-x)
        n2=n*2
c 
        do 1200 i=1,10
c 
cccc        pols(n+i)=pols(i)*dp
cccc        pols(n2+i)=pols(i)*dm
        pols(10+i)=pols(i)*dp
        pols(20+i)=pols(i)*dm
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine matrret3(xsout,wsout,vmatrout,uout,vout,npts)
        implicit real *8 (a-h,o-z)
        save
        dimension uout(900),vout(900),vmatrout(900),
     1      u(900),v(900),vmatr(900),xs0(15),ws0(15),
     2      xsout(30),wsout(30)
c 
        dimension vmatr1(180),vmatr2(180),vmatr3(180),
     1      vmatr4(180),vmatr5(180)
c 
        dimension v1(180),v2(180),v3(180),v4(180),v5(180),
     1      u1(180),u2(180),u3(180),u4(180),u5(180)
c 
c        This subroutine returns to the user several pieces of data
c        to be used for the interpolation and integration on the
c        interval [-1,1] of functions of the form
c 
c        P_n(x),
c 
c        P_n (x) \cdot (1-x) \cdot log(1-x),
c                                                                    (1)
c        P_n (x) \cdot  \cdot(1+x) log(1+x),
c 
c        with n=0,1,2,...,9.  The interpolation is exact in exact
c        arithmetic; in practice, 14-15 digits are obtained with
c        double precision computations.
c 
c        Specifically, this subroutine returns to the user
c        the nodes xsout(i) and weights wsout(i) of the 30-point
c        Generalized Gaussian quadrature on the interval [-1,1] for
c        functions of the form
c 
c        P_n(x), P_n (x) \cdot log(1-x)\cdot (1-x),
c        P_n (x) \cdot log(1+x) \cdot(1+x), P_n(x),                 (2)
c        P_n (x) \cdot log(1-x)^2, P_n (x) \cdot log(1+x)^2,
c 
c        with n running from 0 to 18. It is easily seen that an
c        interpolation formula for functions of the form (1) based
c        on such nodes xsout(i) is exact in exact arithmetic, and is
c        numerically stable.
c 
c        This subroutine also returns to the user the matrix uout of
c        values at the nodes xsout(i) of a set of orthonormal basis
c        functions in the space spanned by the functions (1),
c        and the matrix vout, which is the inverse of uout.
c 
c        Finally, the subroutine returns to the user the matrix vmatrout,
c        used by the subroutine ufuneval (see) to evaluate the basis
c        functions at (arbitrary) user-specified points on [-1,1]
c 
c        PLEASE NOTE THAT EACH OF THE ABOVE PIECES OF INFORMATION CAN
c        BE OBTAINED SEPARATELY BY CALLING ONE OF THE ENTRIES xswsget,
c        uoutget, voutget, vmatrget (see below)
c 
c                          Input parameters: none
c 
c                          Output parameters:
c 
c  xsout - roots of the 30-node Generalized Gaussian quadrature
c       for functions (2) above
c  wsout - the weights corresponding to the roots xsout
c  vmatrout - a 30 \times 30 - matrix used by the subroutine ufuneval
c        (see) to evaluate the basis functions at (arbitrary)
c        user-specified points on [-1,1]
c  uout - the matrix of values at the nodes xsout(i) of a set of
c        orthonormal basis functions in the space spanned by the
c        functions (1)
c  vout - the inverse of uout.
c  npts - the number of nodes xsout returned by the subroutine. It
c        is equal to 30, and is included for historical reasons
c 
        equivalence (u1(1),u(1)),(u2(1),u(181)),(u3(1),u(361)),
     1      (u4(1),u(541)),(u5(1),u(721))
c 
        equivalence (v1(1),v(1)),(v2(1),v(181)),(v3(1),v(361)),
     1      (v4(1),v(541)),(v5(1),v(721))
c 
        equivalence (vmatr1(1),vmatr(1)),(vmatr2(1),vmatr(181)),
     1      (vmatr3(1),vmatr(361)),(vmatr4(1),vmatr(541)),
     2      (vmatr5(1),vmatr(721))
c 
c       following are the nodes of the spectral discretization,
c       and the corresponding weights (only half of them things,
c       where the nodes are positive)
c 
        data xs0/
     1  0.8073141183800038476D-01,0.2393822643556075629D+00,
     2  0.3898227981812752426D+00,0.5273106402545079198D+00,
     3  0.6481288501291113495D+00,0.7498587278001463136D+00,
     4  0.8315271563892494172D+00,0.8936088197087079388D+00,
     5  0.9378752307252273966D+00,0.9670976881178758543D+00,
     6  0.9846319772006163863D+00,0.9939396001308001918D+00,
     7  0.9981286463685694730D+00,0.9996124417392792462D+00,
     8  0.9999659071661138298D+00/
c 
        data ws0/
     1  0.1609903517461948937D+00,0.1554043590434651670D+00,
     2  0.1446800693902384301D+00,0.1296715084226410200D+00,
     3  0.1115592912437890501D+00,0.9173822045171811349D-01,
     4  0.7168028160282363210D-01,0.5278331126950577973D-01,
     5  0.3621844460824368490D-01,0.2279394975441856477D-01,
     6  0.1285944529528482946D-01,0.6279713696535565776D-02,
     7  0.2503053252639005211D-02,0.7256268032426555665D-03,
     8  0.1123734192596079203D-03/
c 
c        following is the matrix vmatr, with nfuns= 30, and ncols= 30
c 
c        following is block number   1 of length 180 elements
c 
c 
        data vmatr1/
     1  0.5190860690666699D+00,0.9418766963124588D-31,
     2  0.6287070065871464D-01,0.5579139876028903D-32,
     3  0.3117986978260758D-02,-.7025376203556077D-32,
     4  0.1932151108622572D-03,-.3945034172036193D-32,
     5  0.3700684663268723D-04,-.1765365050586553D-32,
     6  0.1266153232679384D+00,0.1814180236364119D+00,
     7  0.6451304415569473D-01,0.1417499885873664D-01,
     8  0.1200833719649023D-01,-.1887735987964728D-02,
     9  0.1823102939977707D-02,-.6871335118626719D-03,
     *  0.4640709214798913D-03,-.2433268512282166D-03,
     1  0.1266153232679384D+00,-.1814180236364119D+00,
     2  0.6451304415569473D-01,-.1417499885873664D-01,
     3  0.1200833719649023D-01,0.1887735987964728D-02,
     4  0.1823102939977707D-02,0.6871335118626719D-03,
     5  0.4640709214798913D-03,0.2433268512282166D-03,
     6  -.2252881812368924D-30,0.4072258429813258D+00,
     7  -.9460223601912072D-33,0.2964229861617185D-01,
     8  0.3845175501091084D-31,0.1647551556846964D-03,
     9  0.1670960804028539D-31,-.2816188107075611D-03,
     *  -.2516986805236672D-31,-.1095031086860236D-03,
     1  0.3466793552377989D+00,0.1545965156443364D+00,
     2  0.1419677952554743D+00,0.5693179120507344D-01,
     3  -.4053376351995485D-04,0.7905828875493719D-02,
     4  -.2775304162386797D-02,0.1208440222083929D-02,
     5  -.8692748046910362D-03,0.3477712098290506D-03,
     6  -.3466793552377989D+00,0.1545965156443364D+00,
     7  -.1419677952554743D+00,0.5693179120507344D-01,
     8  0.4053376351995485D-04,0.7905828875493719D-02,
     9  0.2775304162386797D-02,0.1208440222083929D-02,
     *  0.8692748046910362D-03,0.3477712098290506D-03,
     1  -.4076033461275789D+00,-.2167616198041076D-31,
     2  0.6573542920420847D+00,-.6160090465600240D-30,
     3  0.1739252745252064D+00,-.4390594520324175D-30,
     4  0.2284345452759945D-01,-.1803859753681701D-30,
     5  0.2360579354412650D-02,-.1102245342696692D-30,
     6  0.2042369551119511D+00,0.2129962011621904D+00,
     7  0.2217057380643719D+00,0.3087728888462040D+00,
     8  0.1416579790271640D+00,0.5469954725733388D-01,
     9  0.3886306909188649D-01,0.6366308725552560D-03,
     *  0.8020293820943368D-02,-.1842074849559431D-02,
     1  0.2042369551119511D+00,-.2129962011621904D+00,
     2  0.2217057380643719D+00,-.3087728888462040D+00,
     3  0.1416579790271640D+00,-.5469954725733388D-01,
     4  0.3886306909188649D-01,-.6366308725552560D-03,
     5  0.8020293820943368D-02,0.1842074849559431D-02,
     6  -.8653123231157054D-30,-.1299061181859290D+00,
     7  0.1572340906518505D-29,0.7635677339690197D+00,
     8  -.5294846676255314D-30,0.2638122982480983D+00,
     9  -.2008908264210305D-31,0.4674163765372556D-01,
     *  0.3640554554227736D-30,0.6522422403413034D-02,
     1  -.2081849910444527D+00,0.1709055525384973D+00,
     2  0.3052696539075312D+00,0.2845515894759603D+00,
     3  0.3896346674713111D+00,0.1940873290375675D+00,
     4  0.9400393091593424D-01,0.6266716930974585D-01,
     5  0.7218951603619598D-02,0.1486873762291096D-01,
     6  0.2081849910444527D+00,0.1709055525384973D+00,
     7  -.3052696539075312D+00,0.2845515894759603D+00,
     8  -.3896346674713111D+00,0.1940873290375675D+00,
     9  -.9400393091593424D-01,0.6266716930974585D-01,
     *  -.7218951603619598D-02,0.1486873762291096D-01,
     1  0.1976973965144306D+00,-.3734516656619291D-30,
     2  -.4025181000872247D+00,0.7385708225461805D-30,
     3  0.8100739456194036D+00,0.1090994784481838D-29,
     4  0.4418491660986742D+00,-.6483743041708585D-30,
     5  0.1064606297842446D+00,-.4854474006123024D-30,
     6  -.8143707401330786D-01,-.2267885498074705D+00,
     7  0.8693910891766059D-01,0.2020997202028233D+00,
     8  0.2952780954541649D+00,0.4947435840065840D+00,
     9  0.2674750433832420D+00,0.1861353925378834D+00,
     *  0.1088130329507768D+00,0.2500266821359668D-01,
     1  -.8143707401330786D-01,0.2267885498074705D+00,
     2  0.8693910891766059D-01,-.2020997202028233D+00,
     3  0.2952780954541649D+00,-.4947435840065840D+00,
     4  0.2674750433832420D+00,-.1861353925378834D+00,
     5  0.1088130329507768D+00,-.2500266821359668D-01,
     6  -.3838631572961684D-30,0.8820480999196980D-01,
     7  0.9197313521454902D-30,-.5131991383365538D+00,
     8  0.2613359502300411D-30,0.8120943001485559D+00,
     9  -.1432996375758046D-31,0.5732993300753925D+00,
     *  0.5434791453684320D-31,0.1505731788257444D+00,
     1  0.1080134465829285D+00,-.7219758879529641D-01,
     2  -.2801504917427863D+00,0.3739565762050579D-01,
     3  0.1434816651428419D+00,0.2998591406238608D+00,
     4  0.5576080021409004D+00,0.3170734794358965D+00,
     5  0.2559894983260003D+00,0.1423829939065980D+00,
     6  -.1080134465829285D+00,-.7219758879529641D-01,
     7  0.2801504917427863D+00,0.3739565762050579D-01,
     8  -.1434816651428419D+00,0.2998591406238608D+00,
     9  -.5576080021409004D+00,0.3170734794358965D+00,
     *  -.2559894983260003D+00,0.1423829939065980D+00/
c 
c        following is block number   2 of length 180 elements
c 
        data vmatr2/
     1  -.1198040645353958D+00,0.4220652513023511D-30,
     2  0.2560822464580096D+00,-.1249156591947122D-29,
     3  -.7416510591564669D+00,-.7475118592294685D-30,
     4  0.8167028668431188D+00,0.5777665005746207D-30,
     5  0.6821607853010644D+00,-.1022853586171253D-29,
     6  0.6405732994386110D-01,0.1369491242497761D+00,
     7  -.6800314265845261D-01,-.3079078748740727D+00,
     8  -.1651003785766686D-01,0.4293541817932609D-01,
     9  0.2844840008594223D+00,0.6178949867879765D+00,
     *  0.3402682309297377D+00,0.2781443007904754D+00,
     1  0.6405732994386110D-01,-.1369491242497761D+00,
     2  -.6800314265845261D-01,0.3079078748740727D+00,
     3  -.1651003785766686D-01,-.4293541817932609D-01,
     4  0.2844840008594223D+00,-.6178949867879765D+00,
     5  0.3402682309297377D+00,-.2781443007904754D+00,
     6  0.5395170843386642D-30,-.5876913828893095D-01,
     7  -.1170348414832198D-29,0.3803856089448061D+00,
     8  0.3368407547068277D-29,-.9644786984432735D+00,
     9  -.1412382927883208D-29,0.9414960567140432D+00,
     *  0.3224914969999215D-29,0.6630706457856801D+00,
     1  -.8224130524013937D-01,0.6311776819616217D-01,
     2  0.2070633854067871D+00,-.6262439286234566D-01,
     3  -.3537694455591381D+00,-.4815120878303709D-01,
     4  0.3996824581853064D-02,0.2819479484315209D+00,
     5  0.6690021371670097D+00,0.3380117870151612D+00,
     6  0.8224130524013937D-01,0.6311776819616217D-01,
     7  -.2070633854067871D+00,-.6262439286234566D-01,
     8  0.3537694455591381D+00,-.4815120878303709D-01,
     9  -.3996824581853064D-02,0.2819479484315209D+00,
     *  -.6690021371670097D+00,0.3380117870151612D+00,
     1  0.7121189922668154D-01,0.3923421559953240D-30,
     2  -.1594776207358413D+00,-.7533965907819838D-30,
     3  0.5425198455128847D+00,0.8909283873330296D-30,
     4  -.1185560705483878D+01,0.2151554320989026D-29,
     5  0.1351865849203226D+01,-.1095602783547685D-28,
     6  -.3359793190534303D-01,-.9474949851481089D-01,
     7  0.5999989178103426D-01,0.2417783220602193D+00,
     8  -.5675114864565171D-01,-.3929147037708558D+00,
     9  -.2508645673668621D-01,0.1006533356820295D+00,
     *  0.3148342564336787D+00,0.7057423278126757D+00,
     1  -.3359793190534303D-01,0.9474949851481089D-01,
     2  0.5999989178103426D-01,-.2417783220602193D+00,
     3  -.5675114864565171D-01,0.3929147037708558D+00,
     4  -.2508645673668621D-01,-.1006533356820295D+00,
     5  0.3148342564336787D+00,-.7057423278126757D+00,
     6  -.2829056576745307D-31,0.2810635416942685D-01,
     7  -.8551473112821602D-31,-.1769394247456111D+00,
     8  0.3042440329352344D-29,0.5644107361496500D+00,
     9  -.8375993365705359D-29,-.1327073755461514D+01,
     *  0.1468052393286970D-28,0.1996519818156844D+01,
     1  0.3335824121534016D-01,-.2285378185452265D-01,
     2  -.1046645875272570D+00,0.5097782636491152D-01,
     3  0.2354000901884921D+00,-.6374873906808757D-01,
     4  -.4633004704772744D+00,0.4827302474897111D-01,
     5  0.3257180114579109D+00,0.4409716793128581D+00,
     6  -.3335824121534016D-01,-.2285378185452265D-01,
     7  0.1046645875272570D+00,0.5097782636491152D-01,
     8  -.2354000901884921D+00,-.6374873906808757D-01,
     9  0.4633004704772744D+00,0.4827302474897111D-01,
     *  -.3257180114579109D+00,0.4409716793128581D+00,
     1  0.2230107950197852D-01,0.1041447314065674D-28,
     2  -.4933688660903728D-01,-.1304878227249390D-28,
     3  0.2187081621117584D+00,-.1381363052730761D-28,
     4  -.7584768037400983D+00,0.1797988630676581D-28,
     5  0.2553714555992971D+01,-.3103882463180010D-28,
     6  -.3224115411848906D-01,-.6630719534024553D-02,
     7  -.1660974474383666D-02,0.1436741236846466D+00,
     8  -.1028419477952570D+00,-.2747247152704879D+00,
     9  0.8302203323217589D-01,0.1275035531432154D+01,
     *  -.1087631718479918D+01,-.3080109715390725D+01,
     1  -.3224115411848906D-01,0.6630719534024553D-02,
     2  -.1660974474383666D-02,-.1436741236846466D+00,
     3  -.1028419477952570D+00,0.2747247152704879D+00,
     4  0.8302203323217589D-01,-.1275035531432154D+01,
     5  -.1087631718479918D+01,0.3080109715390725D+01,
     6  0.3547931773424427D-29,0.4302437322549794D-01,
     7  -.7170386589090024D-29,-.1937196893899313D-01,
     8  0.1729916489715483D-28,-.2111570105203402D+00,
     9  0.1939101902103566D-28,0.9286878233600557D-01,
     *  -.8357448633326348D-28,-.2512540262892322D+01,
     1  -.8198278035428782D-01,0.1371287078281486D+00,
     2  -.1328375416635606D+00,0.1342700101826260D+00,
     3  -.3336409696972651D+00,0.2581452154170117D+00,
     4  -.4809447534852646D+00,0.5633082767167574D+00,
     5  -.3417423775241333D+01,0.8096391339643235D+01,
     6  0.8198278035428782D-01,0.1371287078281486D+00,
     7  0.1328375416635606D+00,0.1342700101826260D+00,
     8  0.3336409696972651D+00,0.2581452154170117D+00,
     9  0.4809447534852646D+00,0.5633082767167574D+00,
     *  0.3417423775241333D+01,0.8096391339643235D+01/
c 
c        following is block number   3 of length 180 elements
c 
        data vmatr3/
     1  -.4519445222427587D-01,0.1187090348844600D-29,
     2  -.2754725935377573D+00,-.4804651814267573D-29,
     3  -.4393208611108050D+00,0.1088715811175366D-28,
     4  -.1353218323303171D+01,0.7264005091941596D-29,
     5  0.7732404133904157D+00,0.1375820306529475D-30,
     6  0.2242930945767874D+01,-.2633177203243481D+01,
     7  0.2579775498988509D+01,-.3518066061308811D+01,
     8  0.3681499856344962D+01,-.6308744632737693D+01,
     9  0.7813599267489793D+01,-.1398585072542893D+02,
     *  0.3126158051967833D+02,-.1588778521945739D+02,
     1  0.2242930945767874D+01,0.2633177203243481D+01,
     2  0.2579775498988509D+01,0.3518066061308811D+01,
     3  0.3681499856344962D+01,0.6308744632737693D+01,
     4  0.7813599267489793D+01,0.1398585072542893D+02,
     5  0.3126158051967833D+02,0.1588778521945739D+02,
     6  0.2677851030426297D-27,0.1505128537819797D+02,
     7  -.6301460405949597D-27,-.4876487683850678D+01,
     8  0.4889154821719518D-27,-.3600590997230667D+01,
     9  -.2069838079386056D-26,0.2496802246294824D+02,
     *  0.6748499011520918D-26,0.3774947449370860D+02,
     1  -.2111967785437550D+02,0.3995069142139690D+02,
     2  -.3605332767171667D+02,0.4267969943861913D+02,
     3  -.5429492840730098D+02,0.6230005236237095D+02,
     4  -.6414501912983671D+02,0.1039059987488881D+03,
     5  -.7392667656214123D+02,-.4303574740827587D+02,
     6  0.2111967785437550D+02,0.3995069142139690D+02,
     7  0.3605332767171667D+02,0.4267969943861913D+02,
     8  0.5429492840730098D+02,0.6230005236237095D+02,
     9  0.6414501912983671D+02,0.1039059987488881D+03,
     *  0.7392667656214123D+02,-.4303574740827587D+02,
     1  -.2423656558843612D+01,-.3745406451321677D-26,
     2  -.6315652883778623D+01,-.4358973208794650D-25,
     3  -.2298346172796427D+01,0.1018413922037642D-24,
     4  0.3987297578708171D+02,0.1176766553475021D-24,
     5  0.4425689810141042D+02,-.3770619039754974D-27,
     6  0.1629259914240558D+03,-.1731857995079773D+03,
     7  0.1731167127402070D+03,-.1892644463281696D+03,
     8  0.2012150407664095D+03,-.1876086695386810D+03,
     9  0.2171417468807008D+03,-.9912724077499929D+02,
     *  -.1975157281919880D+03,0.4772227854282930D+02,
     1  0.1629259914240558D+03,0.1731857995079773D+03,
     2  0.1731167127402070D+03,0.1892644463281696D+03,
     3  0.2012150407664095D+03,0.1876086695386810D+03,
     4  0.2171417468807008D+03,0.9912724077499929D+02,
     5  -.1975157281919880D+03,-.4772227854282930D+02,
     6  -.2043317493974032D-25,0.3275676463580989D+03,
     7  0.5145746467745036D-25,-.8949145552089381D+02,
     8  -.1112487625575881D-24,0.1625836402678505D+03,
     9  0.3072725516349378D-25,-.3926223857616492D+03,
     *  0.2262482902486713D-23,-.6676189137635828D+03,
     1  -.4009305803004565D+03,0.7800229523895357D+03,
     2  -.6112615222580785D+03,0.6625558672441790D+03,
     3  -.4961337014545619D+03,0.3595625257195584D+03,
     4  -.1969625368478577D+03,-.6443661157949528D+03,
     5  0.8226559880453132D+03,0.2099754518910125D+03,
     6  0.4009305803004565D+03,0.7800229523895357D+03,
     7  0.6112615222580785D+03,0.6625558672441790D+03,
     8  0.4961337014545619D+03,0.3595625257195584D+03,
     9  0.1969625368478577D+03,-.6443661157949528D+03,
     *  -.8226559880453132D+03,0.2099754518910125D+03,
     1  0.2133258184770809D+02,-.4866900007591294D-25,
     2  -.2685230268850312D+02,-.5075651613286407D-24,
     3  0.2207144494388403D+03,0.3059820769561111D-24,
     4  -.1594654806917585D+04,0.2788435024611937D-23,
     5  -.2250310168050241D+04,-.1615558619327445D-23,
     6  0.3116230546179615D+04,-.2923223672204792D+04,
     7  0.2475370940752951D+04,-.1528043844700531D+04,
     8  0.2658732937912536D+03,0.6193849747192275D+03,
     9  -.4207041394113960D+04,0.3343737706065144D+04,
     *  0.1680964082085384D+04,-.2142585512900544D+03,
     1  0.3116230546179615D+04,0.2923223672204792D+04,
     2  0.2475370940752951D+04,0.1528043844700531D+04,
     3  0.2658732937912536D+03,-.6193849747192275D+03,
     4  -.4207041394113960D+04,-.3343737706065144D+04,
     5  0.1680964082085384D+04,0.2142585512900544D+03,
     6  0.9954232562505281D-24,0.6120408280052088D+04,
     7  -.2587615500925754D-23,-.5766495357041915D+03,
     8  0.5585067637987079D-23,-.5787521128479725D+04,
     9  -.3076633976117882D-24,0.4828683552455679D+04,
     *  -.9159549792132871D-22,0.1552423796313185D+05,
     1  -.6911525715357361D+04,0.1274143624150054D+05,
     2  -.5829939884436477D+04,0.2513254319636186D+04,
     3  -.1205281909311690D+03,-.1084256392954625D+05,
     4  0.1456834217936439D+05,-.7069893176747019D+04,
     5  -.1202894570941072D+05,-.1513731745525924D+04,
     6  0.6911525715357361D+04,0.1274143624150054D+05,
     7  0.5829939884436477D+04,0.2513254319636186D+04,
     8  0.1205281909311690D+03,-.1084256392954625D+05,
     9  -.1456834217936439D+05,-.7069893176747019D+04,
     *  0.1202894570941072D+05,-.1513731745525924D+04/
c 
c        following is block number   4 of length 180 elements
c 
        data vmatr4/
     1  0.1424824901985506D+03,-.8299772428024136D-23,
     2  0.3474965926765900D+04,0.2288617633814403D-23,
     3  -.4593252962398307D+05,0.1189892008140665D-23,
     4  0.3199012764743287D+05,-.1907819301923842D-23,
     5  0.1053844041255780D+06,0.1325901537781713D-22,
     6  0.7923320141153010D+05,-.5347220252730389D+05,
     7  0.9840999138457713D+04,0.1163639916369922D+05,
     8  -.8916973725773980D+05,0.9263530410745086D+05,
     9  -.1025831721553501D+05,-.9126689662174181D+05,
     *  -.1902609437798623D+05,0.1289490454886533D+04,
     1  0.7923320141153010D+05,0.5347220252730389D+05,
     2  0.9840999138457713D+04,-.1163639916369922D+05,
     3  -.8916973725773980D+05,-.9263530410745086D+05,
     4  -.1025831721553501D+05,0.9126689662174181D+05,
     5  -.1902609437798623D+05,-.1289490454886533D+04,
     6  0.3484410709558824D-22,-.1814247055444161D+06,
     7  -.4383029845660981D-22,0.1233999940665918D+06,
     8  -.1401208748936001D-22,0.1726528020512690D+05,
     9  0.3088025800607170D-21,0.1353612251026699D+06,
     *  -.2517882944319958D-21,0.2861586139984804D+06,
     1  0.1611605665393441D+06,-.2630116860439235D+06,
     2  0.5345580934168737D+05,0.1798859294757339D+06,
     3  -.2618866036112362D+06,0.3356813151607411D+06,
     4  0.8608902646739449D+05,-.3840542798691744D+06,
     5  -.1693243867582155D+06,-.1267036770893636D+05,
     6  -.1611605665393441D+06,-.2630116860439235D+06,
     7  -.5345580934168737D+05,0.1798859294757339D+06,
     8  0.2618866036112362D+06,0.3356813151607411D+06,
     9  -.8608902646739449D+05,-.3840542798691744D+06,
     *  0.1693243867582155D+06,-.1267036770893636D+05,
     1  0.2585720874744819D+06,0.5459284802228740D-22,
     2  -.7610314074829316D+06,-.5750625729915210D-20,
     3  -.2102263081059073D+06,0.4464590282926914D-21,
     4  -.1062344102998029D+07,-.2700461137453328D-21,
     5  -.3263030405419245D+07,0.3047510912850490D-21,
     6  0.2490315848921509D+07,-.1508842275460512D+07,
     7  -.1476948143035753D+07,0.2696009764290626D+07,
     8  -.2916641604345258D+07,-.1477116553568333D+07,
     9  0.3486706516424937D+07,0.2106520874272544D+07,
     *  0.2437250076554367D+06,-.9521704832266498D+04,
     1  0.2490315848921509D+07,0.1508842275460512D+07,
     2  -.1476948143035753D+07,-.2696009764290626D+07,
     3  -.2916641604345258D+07,0.1477116553568333D+07,
     4  0.3486706516424937D+07,-.2106520874272544D+07,
     5  0.2437250076554367D+06,0.9521704832266498D+04,
     6  -.6183626316817279D-19,-.8069528694184526D+07,
     7  -.1557286415829250D-19,0.8559191584582117D+07,
     8  0.2404443149661216D-19,-.3448017069999335D+07,
     9  -.5314214552224225D-19,-.1421829369381912D+08,
     *  -.3914495604440262D-19,-.4759162531158232D+07,
     1  0.7092465804739067D+07,-.7008016888987293D+07,
     2  -.4378932388619389D+07,0.1513256687838507D+08,
     3  -.7323674304092003D+07,-.1054578689599380D+08,
     4  0.1019176664427810D+08,0.1020283628404903D+08,
     5  0.2343439308921740D+07,0.1166612639826824D+06,
     6  -.7092465804739067D+07,-.7008016888987293D+07,
     7  0.4378932388619389D+07,0.1513256687838507D+08,
     8  0.7323674304092003D+07,-.1054578689599380D+08,
     9  -.1019176664427810D+08,0.1020283628404903D+08,
     *  -.2343439308921740D+07,0.1166612639826824D+06,
     1  0.3692909155712698D+07,-.2213507733391912D-15,
     2  -.7438395494254816D+08,0.1096766147576451D-15,
     3  0.1081747193901678D+08,0.2565590001518817D-16,
     4  0.2155499329586463D+09,0.4558724938589720D-15,
     5  0.8198551718812234D+08,-.2031757983289935D-18,
     6  0.1352962978136282D+09,-.3440593931468877D+08,
     7  -.2025040312147361D+09,0.1308729342843078D+09,
     8  0.1373039028540160D+09,-.1386618541383676D+09,
     9  -.1525281051975389D+09,-.4362251224079190D+08,
     *  -.3267817985326198D+07,0.8033815240908052D+05,
     1  0.1352962978136282D+09,0.3440593931468877D+08,
     2  -.2025040312147361D+09,-.1308729342843078D+09,
     3  0.1373039028540160D+09,0.1386618541383676D+09,
     4  -.1525281051975389D+09,0.4362251224079190D+08,
     5  -.3267817985326198D+07,-.8033815240908052D+05,
     6  0.5135016724576667D-16,-.3427997228073664D+09,
     7  0.9378755662570491D-14,0.1702259670169568D+09,
     8  -.1394303089471572D-14,0.6803058691808170D+08,
     9  -.2757330567588799D-13,0.1021857476660534D+10,
     *  -.1064115072085929D-13,0.8595927611757958D+08,
     1  0.4552323371997043D+09,-.1118603702440569D+09,
     2  -.7018938279008217D+09,0.5057204979445499D+09,
     3  0.6201354823569313D+09,-.4908116552628774D+09,
     4  -.6988917476593455D+09,-.2636704530243885D+09,
     5  -.3641288940735648D+08,-.1256305405740580D+07,
     6  -.4552323371997043D+09,-.1118603702440569D+09,
     7  0.7018938279008217D+09,0.5057204979445499D+09,
     8  -.6201354823569313D+09,-.4908116552628774D+09,
     9  0.6988917476593455D+09,-.2636704530243885D+09,
     *  0.3641288940735648D+08,-.1256305405740580D+07/
c 
c        following is block number   5 of length 180 elements
c 
        data vmatr5/
     1  -.5313832747103602D+10,-.6256936957126857D-13,
     2  -.9258410749526798D+09,0.2233105013322759D-14,
     3  -.4076296496630839D+10,-.1072941044392323D-13,
     4  -.2225085869348745D+11,0.2638302354053794D-13,
     5  -.2098724371212300D+10,-.3940378662951591D-13,
     6  0.1136544619851658D+11,0.5963742684169094D+10,
     7  -.1724600005618587D+11,-.8863805245501890D+10,
     8  0.1262262847628579D+11,0.1444808781634454D+11,
     9  0.5711804248305848D+10,0.9556921377089287D+09,
     *  0.4910847967998103D+08,-.7937405679736574D+06,
     1  0.1136544619851658D+11,-.5963742684169094D+10,
     2  -.1724600005618587D+11,0.8863805245501890D+10,
     3  0.1262262847628579D+11,-.1444808781634454D+11,
     4  0.5711804248305848D+10,-.9556921377089287D+09,
     5  0.4910847967998103D+08,0.7937405679736574D+06,
     6  0.2700429045953624D-12,0.3944904513140195D+11,
     7  0.1274227155315915D-11,-.3725068677436177D+11,
     8  0.1684122847624703D-12,0.9429765074072301D+11,
     9  0.8999486448561353D-12,0.4892629178611181D+11,
     *  0.4894682709718860D-12,0.1563391590273382D+10,
     1  -.4307728889720394D+11,-.2478735303123212D+11,
     2  0.6607381644423791D+11,0.4083140365223793D+11,
     3  -.4621282674046010D+11,-.6270404906245213D+11,
     4  -.2920897981392729D+11,-.6339047357300685D+10,
     5  -.5886431448702278D+09,-.1479683400721104D+08,
     6  0.4307728889720394D+11,-.2478735303123212D+11,
     7  -.6607381644423791D+11,0.4083140365223793D+11,
     8  0.4621282674046010D+11,-.6270404906245213D+11,
     9  0.2920897981392729D+11,-.6339047357300685D+10,
     *  0.5886431448702278D+09,-.1479683400721104D+08,
     1  0.1722894448092670D+13,-.4183414949761398D-08,
     2  0.9119114921863091D+12,0.6471589749444142D-08,
     3  -.4562633269620064D+13,0.6501889321680611D-08,
     4  -.1466671765383033D+13,-.1726826817850891D-08,
     5  -.5224650084054552D+11,-.7595772850315042D-09,
     6  -.1478804482897373D+13,-.2049625396904567D+13,
     7  0.5123643269597118D+12,0.2430788727158754D+13,
     8  0.2019261109880588D+13,0.8409434332905432D+12,
     9  0.1895453512098165D+12,0.2105272979977230D+11,
     *  0.7843378025334641D+09,-.8747359845433506D+07,
     1  -.1478804482897373D+13,0.2049625396904567D+13,
     2  0.5123643269597118D+12,-.2430788727158754D+13,
     3  0.2019261109880588D+13,-.8409434332905432D+12,
     4  0.1895453512098165D+12,-.2105272979977230D+11,
     5  0.7843378025334641D+09,0.8747359845433505D+07,
     6  0.3004102779213300D-06,0.1283108049832896D+14,
     7  0.1245212240954499D-06,-.1227665093972232D+14,
     8  -.6446429148848069D-06,-.1507546222146058D+14,
     9  -.2075109816038374D-06,-.2035019809159856D+13,
     *  -.2389126006650022D-07,-.2960498549830592D+11,
     1  -.5792313795876348D+13,-.8232157928529096D+13,
     2  0.1790481679273997D+13,0.9971549627144425D+13,
     3  0.8905344261775276D+13,0.4089907204765451D+13,
     4  0.1068388906214385D+13,0.1525212723264332D+12,
     5  0.1011032659372751D+11,0.1913282544905190D+09,
     6  0.5792313795876348D+13,-.8232157928529096D+13,
     7  -.1790481679273997D+13,0.9971549627144425D+13,
     8  -.8905344261775276D+13,0.4089907204765451D+13,
     9  -.1068388906214385D+13,0.1525212723264332D+12,
     *  -.1011032659372751D+11,0.1913282544905190D+09,
     1  -.6805784028763431D+15,0.2519269814099966D-02,
     2  -.1883189060463907D+16,-.1012905274318669D-02,
     3  -.7588934272958729D+15,0.1008335106090694D-02,
     4  -.7653883508973891D+14,0.2989389641500426D-02,
     5  -.1278091283020069D+13,0.9649054069337576D-03,
     6  0.3407054710018120D+15,0.7694760944023110D+15,
     7  0.7213288699002137D+15,0.4174654402674320D+15,
     8  0.1585485594128051D+15,0.3905517614348636D+14,
     9  0.5867079384199861D+13,0.4665152687147462D+12,
     *  0.1314119363894102D+11,-.1053487662587143D+09,
     1  0.3407054710018120D+15,-.7694760944023110D+15,
     2  0.7213288699002137D+15,-.4174654402674320D+15,
     3  0.1585485594128051D+15,-.3905517614348636D+14,
     4  0.5867079384199865D+13,-.4665152687147504D+12,
     5  0.1314119363894200D+11,0.1053487662584902D+09,
     6  0.9983725143706265D-01,0.7076240504077969D+16,
     7  -.2249922322789175D+01,0.6152321328004941D+16,
     8  0.5136189508232793D+00,0.1375332909048575D+16,
     9  0.6747939648536814D+01,0.7569126407920718D+14,
     *  0.2554319442904379D+01,0.5724730101898438D+12,
     1  -.1406115098687702D+16,-.3216105822921871D+16,
     2  -.3095612801469361D+16,-.1869920381156787D+16,
     3  -.7572005922381434D+15,-.2050865609742830D+15,
     4  -.3563922284701381D+14,-.3632547581748715D+13,
     5  -.1802438520381346D+12,-.2646267699017172D+10,
     6  0.1406115098687710D+16,-.3216105822921869D+16,
     7  0.3095612801469348D+16,-.1869920381156795D+16,
     8  0.7572005922381518D+15,-.2050865609742744D+15,
     9  0.3563922284700432D+14,-.3632547581746001D+13,
     *  0.1802438520379303D+12,-.2646267699023298D+10/
c 
c        following is the matrix u, with npts= 30
c 
c        following is block number   1 of length 180 elements
c 
        data u1/
     1  0.1137263252338031D+01,0.1136629219814803D+01,
     2  0.1134008240454126D+01,0.1126738424888053D+01,
     3  0.1111000059960315D+01,0.1082527123314642D+01,
     4  0.1037872154425957D+01,0.9758101829134152D+00,
     5  0.8983181239827154D+00,0.8107013760733012D+00,
     6  0.7207139603406841D+00,0.6369536792997493D+00,
     7  0.5672586567990836D+00,0.5176711418311745D+00,
     8  0.4920142698239177D+00,0.4920142698239177D+00,
     9  0.5176711418311745D+00,0.5672586567990836D+00,
     *  0.6369536792997493D+00,0.7207139603406841D+00,
     1  0.8107013760733012D+00,0.8983181239827154D+00,
     2  0.9758101829134152D+00,0.1037872154425957D+01,
     3  0.1082527123314642D+01,0.1111000059960315D+01,
     4  0.1126738424888053D+01,0.1134008240454126D+01,
     5  0.1136629219814803D+01,0.1137263252338031D+01,
     6  -.1415284549701814D+01,-.1414880769487940D+01,
     7  -.1412513733688901D+01,-.1404515742147085D+01,
     8  -.1384782748785232D+01,-.1345545308780060D+01,
     9  -.1279314010052062D+01,-.1181382768635140D+01,
     *  -.1051948972538377D+01,-.8969740721905299D+00,
     1  -.7269553577364316D+00,-.5534514652718024D+00,
     2  -.3851230325659789D+00,-.2257224381256085D+00,
     3  -.7419936292064160D-01,0.7419936292064160D-01,
     4  0.2257224381256085D+00,0.3851230325659789D+00,
     5  0.5534514652718024D+00,0.7269553577364316D+00,
     6  0.8969740721905299D+00,0.1051948972538377D+01,
     7  0.1181382768635140D+01,0.1279314010052062D+01,
     8  0.1345545308780060D+01,0.1384782748785232D+01,
     9  0.1404515742147085D+01,0.1412513733688901D+01,
     *  0.1414880769487940D+01,0.1415284549701814D+01,
     1  0.2097708028616502D+01,0.2092441333952270D+01,
     2  0.2070544989788420D+01,0.2009878118238820D+01,
     3  0.1880214459819651D+01,0.1653148295430204D+01,
     4  0.1318651370814090D+01,0.8996507565036436D+00,
     5  0.4517851491658301D+00,0.4316484025786219D-01,
     6  -.2747109938129373D+00,-.4850468824248741D+00,
     7  -.6025728291171553D+00,-.6567622653688939D+00,
     8  -.6756090077537490D+00,-.6756090077537490D+00,
     9  -.6567622653688939D+00,-.6025728291171553D+00,
     *  -.4850468824248741D+00,-.2747109938129373D+00,
     1  0.4316484025786219D-01,0.4517851491658301D+00,
     2  0.8996507565036436D+00,0.1318651370814090D+01,
     3  0.1653148295430204D+01,0.1880214459819651D+01,
     4  0.2009878118238820D+01,0.2070544989788420D+01,
     5  0.2092441333952270D+01,0.2097708028616502D+01,
     6  -.2772606584822725D+01,-.2760938489555236D+01,
     7  -.2712672316657737D+01,-.2580181351490750D+01,
     8  -.2302441036627445D+01,-.1834151228896078D+01,
     9  -.1189633625390545D+01,-.4693787745584557D+00,
     *  0.1694551638791972D+00,0.5929830627516095D+00,
     1  0.7577617542809217D+00,0.7117682837807237D+00,
     2  0.5436326041684152D+00,0.3296785910052084D+00,
     3  0.1093501624243294D+00,-.1093501624243294D+00,
     4  -.3296785910052084D+00,-.5436326041684152D+00,
     5  -.7117682837807237D+00,-.7577617542809217D+00,
     6  -.5929830627516095D+00,-.1694551638791972D+00,
     7  0.4693787745584557D+00,0.1189633625390545D+01,
     8  0.1834151228896078D+01,0.2302441036627445D+01,
     9  0.2580181351490750D+01,0.2712672316657737D+01,
     *  0.2760938489555236D+01,0.2772606584822725D+01,
     1  0.3034454705281747D+01,0.3013159796395392D+01,
     2  0.2924657069329268D+01,0.2682552880631545D+01,
     3  0.2184196749450690D+01,0.1381483633999642D+01,
     4  0.3793623204652597D+00,-.5335589896017575D+00,
     5  -.1026310021327776D+01,-.9744586992039805D+00,
     6  -.5586247939732860D+00,-.8085342903334790D-01,
     7  0.2759411097016343D+00,0.4867352314135936D+00,
     8  0.5822576071855785D+00,0.5822576071855785D+00,
     9  0.4867352314135936D+00,0.2759411097016343D+00,
     *  -.8085342903334790D-01,-.5586247939732860D+00,
     1  -.9744586992039805D+00,-.1026310021327776D+01,
     2  -.5335589896017575D+00,0.3793623204652597D+00,
     3  0.1381483633999642D+01,0.2184196749450690D+01,
     4  0.2682552880631545D+01,0.2924657069329268D+01,
     5  0.3013159796395392D+01,0.3034454705281747D+01,
     6  -.3200406508843042D+01,-.3168136311997931D+01,
     7  -.3034359569460589D+01,-.2672142501183400D+01,
     8  -.1945957073225294D+01,-.8421469170236391D+00,
     9  0.3751993250788516D+00,0.1187712308769530D+01,
     *  0.1184065877448746D+01,0.4928882818465321D+00,
     1  -.2836093918294196D+00,-.6611426976948118D+00,
     2  -.6364732442145117D+00,-.4283215334535179D+00,
     3  -.1516739642253599D+00,0.1516739642253599D+00,
     4  0.4283215334535179D+00,0.6364732442145117D+00,
     5  0.6611426976948118D+00,0.2836093918294196D+00,
     6  -.4928882818465321D+00,-.1184065877448746D+01,
     7  -.1187712308769530D+01,-.3751993250788516D+00,
     8  0.8421469170236391D+00,0.1945957073225294D+01,
     9  0.2672142501183400D+01,0.3034359569460589D+01,
     *  0.3168136311997931D+01,0.3200406508843042D+01/
c 
c        following is block number   2 of length 180 elements
c 
        data u2/
     1  0.2792500730066501D+01,0.2757012794719952D+01,
     2  0.2609891854143673D+01,0.2213426616398066D+01,
     3  0.1431205580241231D+01,0.2921104800970631D+00,
     4  -.8250126884650850D+00,-.1275211134571350D+01,
     5  -.7058747365100495D+00,0.3839425839994380D+00,
     6  0.9933009326662305D+00,0.7062059386585097D+00,
     7  0.3096096453776940D-01,-.4282481508880927D+00,
     8  -.5813122188021777D+00,-.5813122188021777D+00,
     9  -.4282481508880927D+00,0.3096096453776940D-01,
     *  0.7062059386585097D+00,0.9933009326662305D+00,
     1  0.3839425839994380D+00,-.7058747365100495D+00,
     2  -.1275211134571350D+01,-.8250126884650850D+00,
     3  0.2921104800970631D+00,0.1431205580241231D+01,
     4  0.2213426616398066D+01,0.2609891854143673D+01,
     5  0.2757012794719952D+01,0.2792500730066501D+01,
     6  -.2367021160974388D+01,-.2329469359382649D+01,
     7  -.2175188675087089D+01,-.1765927361885092D+01,
     8  -.9840793007823477D+00,0.7832360240356199D-01,
     9  0.9484578748038819D+00,0.9765369818627541D+00,
     *  0.5140116597892515D-01,-.9343123827936188D+00,
     1  -.8919908465983105D+00,0.8120480666548857D-01,
     2  0.8336840231788015D+00,0.7503863054549613D+00,
     3  0.2559517274037907D+00,-.2559517274037907D+00,
     4  -.7503863054549613D+00,-.8336840231788015D+00,
     5  -.8120480666548857D-01,0.8919908465983105D+00,
     6  0.9343123827936188D+00,-.5140116597892515D-01,
     7  -.9765369818627541D+00,-.9484578748038819D+00,
     8  -.7832360240356199D-01,0.9840793007823477D+00,
     9  0.1765927361885092D+01,0.2175188675087089D+01,
     *  0.2329469359382649D+01,0.2367021160974388D+01,
     1  0.1754361614222790D+01,0.1724832315131352D+01,
     2  0.1602480279020129D+01,0.1276364290731360D+01,
     3  0.6550000291620660D+00,-.1695686364067352D+00,
     4  -.7746382388874168D+00,-.6196800601458804D+00,
     5  0.2528652891784685D+00,0.8607799182690834D+00,
     6  0.3233556975673544D+00,-.7535040129675802D+00,
     7  -.9472584605248116D+00,-.1060187731537234D-01,
     8  0.9489734704130381D+00,0.9489734704130381D+00,
     9  -.1060187731537234D-01,-.9472584605248116D+00,
     *  -.7535040129675802D+00,0.3233556975673544D+00,
     1  0.8607799182690834D+00,0.2528652891784685D+00,
     2  -.6196800601458804D+00,-.7746382388874168D+00,
     3  -.1695686364067352D+00,0.6550000291620660D+00,
     4  0.1276364290731360D+01,0.1602480279020129D+01,
     5  0.1724832315131352D+01,0.1754361614222790D+01,
     6  -.1747240825062306D+01,-.1712717108323309D+01,
     7  -.1572391334085862D+01,-.1207443525985710D+01,
     8  -.5392275158860430D+00,0.2836341797411758D+00,
     9  0.7714781638820850D+00,0.4374223438736767D+00,
     *  -.4330899612700505D+00,-.6971985184206539D+00,
     1  0.1802719730709580D+00,0.9359817375289874D+00,
     2  0.1884225599144946D+00,-.1061440800056436D+01,
     3  -.7454790223424396D+00,0.7454790223424396D+00,
     4  0.1061440800056436D+01,-.1884225599144946D+00,
     5  -.9359817375289874D+00,-.1802719730709580D+00,
     6  0.6971985184206539D+00,0.4330899612700505D+00,
     7  -.4374223438736767D+00,-.7714781638820850D+00,
     8  -.2836341797411758D+00,0.5392275158860430D+00,
     9  0.1207443525985710D+01,0.1572391334085862D+01,
     *  0.1712717108323309D+01,0.1747240825062306D+01,
     1  -.2284503357658974D+01,-.2241726677294864D+01,
     2  -.2063819530790370D+01,-.1591136370781003D+01,
     3  -.7082787786521192D+00,0.3902000566149867D+00,
     4  0.1012825413247230D+01,0.4867180527250674D+00,
     5  -.6418492104323924D+00,-.6945759948179559D+00,
     6  0.5211257376155955D+00,0.7475958274558100D+00,
     7  -.6092319242815857D+00,-.7725764714612253D+00,
     8  0.7334880241127122D+00,0.7334880241127122D+00,
     9  -.7725764714612253D+00,-.6092319242815857D+00,
     *  0.7475958274558100D+00,0.5211257376155955D+00,
     1  -.6945759948179559D+00,-.6418492104323924D+00,
     2  0.4867180527250674D+00,0.1012825413247230D+01,
     3  0.3902000566149867D+00,-.7082787786521192D+00,
     4  -.1591136370781003D+01,-.2063819530790370D+01,
     5  -.2241726677294864D+01,-.2284503357658974D+01,
     6  -.3954311978914984D+01,-.3834481314989599D+01,
     7  -.3383951574783495D+01,-.2311685805166204D+01,
     8  -.5893845057410609D+00,0.1062783359843805D+01,
     9  0.1355031917333477D+01,0.1130459171328051D-01,
     *  -.1053404128929248D+01,-.1743412638907552D+00,
     1  0.8745823120736635D+00,-.1170686563035750D-01,
     2  -.7763101755067945D+00,0.3161437194117501D+00,
     3  0.6014107786640433D+00,-.6014107786640433D+00,
     4  -.3161437194117501D+00,0.7763101755067945D+00,
     5  0.1170686563035750D-01,-.8745823120736635D+00,
     6  0.1743412638907552D+00,0.1053404128929248D+01,
     7  -.1130459171328051D-01,-.1355031917333477D+01,
     8  -.1062783359843805D+01,0.5893845057410609D+00,
     9  0.2311685805166204D+01,0.3383951574783495D+01,
     *  0.3834481314989599D+01,0.3954311978914984D+01/
c 
c        following is block number   3 of length 180 elements
c 
        data u3/
     1  0.5895440818793325D+01,0.5601665095955896D+01,
     2  0.4642372223311012D+01,0.2676703474076314D+01,
     3  0.6314805851188785D-01,-.1718490874049116D+01,
     4  -.1207685845495723D+01,0.6410424046225884D+00,
     5  0.9365189470490469D+00,-.4807159607540348D+00,
     6  -.5626636833545676D+00,0.5827988171680165D+00,
     7  0.1317337474243357D+00,-.5747221174160635D+00,
     8  0.3076146741044873D+00,0.3076146741044873D+00,
     9  -.5747221174160635D+00,0.1317337474243357D+00,
     *  0.5827988171680165D+00,-.5626636833545676D+00,
     1  -.4807159607540348D+00,0.9365189470490469D+00,
     2  0.6410424046225884D+00,-.1207685845495723D+01,
     3  -.1718490874049116D+01,0.6314805851188785D-01,
     4  0.2676703474076314D+01,0.4642372223311012D+01,
     5  0.5601665095955896D+01,0.5895440818793325D+01,
     6  -.8536226134455404D+01,-.7720698248335413D+01,
     7  -.5530152321222323D+01,-.2016804918325086D+01,
     8  0.1246418721234055D+01,0.1990448037025180D+01,
     9  0.1340568899731869D+00,-.1159417725844822D+01,
     *  -.2977412220798893D-01,0.7802670398487491D+00,
     1  -.3074608547618163D+00,-.4231605241406142D+00,
     2  0.5374866323029744D+00,-.4997463944104848D-01,
     3  -.4507512277522176D+00,0.4507512277522176D+00,
     4  0.4997463944104848D-01,-.5374866323029744D+00,
     5  0.4231605241406142D+00,0.3074608547618163D+00,
     6  -.7802670398487491D+00,0.2977412220798893D-01,
     7  0.1159417725844822D+01,-.1340568899731869D+00,
     8  -.1990448037025180D+01,-.1246418721234055D+01,
     9  0.2016804918325086D+01,0.5530152321222323D+01,
     *  0.7720698248335413D+01,0.8536226134455404D+01,
     1  0.1103141798908718D+02,0.9429328895428397D+01,
     2  0.5663963957677660D+01,0.7434021864272668D+00,
     3  -.2295171659869548D+01,-.1433512818174175D+01,
     4  0.9138033825965245D+00,0.7873488605138831D+00,
     5  -.7417309981958154D+00,-.2119327385193129D+00,
     6  0.6443084123780363D+00,-.2918478162354762D+00,
     7  -.2790713808008797D+00,0.5018341421081805D+00,
     8  -.2387120395172781D+00,-.2387120395172781D+00,
     9  0.5018341421081805D+00,-.2790713808008797D+00,
     *  -.2918478162354762D+00,0.6443084123780363D+00,
     1  -.2119327385193129D+00,-.7417309981958154D+00,
     2  0.7873488605138831D+00,0.9138033825965245D+00,
     3  -.1433512818174175D+01,-.2295171659869548D+01,
     4  0.7434021864272668D+00,0.5663963957677660D+01,
     5  0.9429328895428397D+01,0.1103141798908718D+02,
     6  -.1248433956741136D+02,-.9905243398881756D+01,
     7  -.4551734282146133D+01,0.1015569430488300D+01,
     8  0.2612538504966678D+01,0.1010098234840112D+00,
     9  -.1398381858102992D+01,0.2534646616832638D+00,
     *  0.7855118960036487D+00,-.6196229125747938D+00,
     1  -.1314044877691464D+00,0.5874713006713431D+00,
     2  -.4396313675285939D+00,-.6395397323694002D-01,
     3  0.4729463912118726D+00,-.4729463912118726D+00,
     4  0.6395397323694002D-01,0.4396313675285939D+00,
     5  -.5874713006713431D+00,0.1314044877691464D+00,
     6  0.6196229125747938D+00,-.7855118960036487D+00,
     7  -.2534646616832638D+00,0.1398381858102992D+01,
     8  -.1010098234840112D+00,-.2612538504966678D+01,
     9  -.1015569430488300D+01,0.4551734282146133D+01,
     *  0.9905243398881756D+01,0.1248433956741136D+02,
     1  0.1405503803813321D+02,0.1015588891905852D+02,
     2  0.3058957753098879D+01,-.2513306070705410D+01,
     3  -.2068971809302235D+01,0.1176225107087802D+01,
     4  0.8902444808437068D+00,-.1041972602627939D+01,
     5  0.1404656527246642D-01,0.6972804829906544D+00,
     6  -.5949859102354189D+00,0.3179307695229410D-01,
     7  0.4509726853929861D+00,-.5378151058859724D+00,
     8  0.2319354724081329D+00,0.2319354724081329D+00,
     9  -.5378151058859724D+00,0.4509726853929861D+00,
     *  0.3179307695229410D-01,-.5949859102354189D+00,
     1  0.6972804829906544D+00,0.1404656527246642D-01,
     2  -.1041972602627939D+01,0.8902444808437068D+00,
     3  0.1176225107087802D+01,-.2068971809302235D+01,
     4  -.2513306070705410D+01,0.3058957753098879D+01,
     5  0.1015588891905852D+02,0.1405503803813321D+02,
     6  -.1503407479953177D+02,-.9835236301676876D+01,
     7  -.1448845152396745D+01,0.3313033816630533D+01,
     8  0.9463206827365319D+00,-.1775205469040251D+01,
     9  0.1497264823878356D+00,0.9769720564241459D+00,
     *  -.8016423172901601D+00,0.2040394971479546D-01,
     1  0.5670293876333268D+00,-.6346072864621087D+00,
     2  0.2821724981337066D+00,0.1985321589510753D+00,
     3  -.5218979380394439D+00,0.5218979380394439D+00,
     4  -.1985321589510753D+00,-.2821724981337066D+00,
     5  0.6346072864621087D+00,-.5670293876333268D+00,
     6  -.2040394971479546D-01,0.8016423172901601D+00,
     7  -.9769720564241459D+00,-.1497264823878356D+00,
     8  0.1775205469040251D+01,-.9463206827365319D+00,
     9  -.3313033816630533D+01,0.1448845152396745D+01,
     *  0.9835236301676876D+01,0.1503407479953177D+02/
c 
c        following is block number   4 of length 180 elements
c 
        data u4/
     1  0.1655130527542360D+02,0.9450767760517566D+01,
     2  -.3749987235770745D+00,-.3590384795527415D+01,
     3  0.5001765333574129D+00,0.1573937740717607D+01,
     4  -.1131666042007823D+01,-.1239389574083093D+00,
     5  0.8319248661696513D+00,-.7428216449057439D+00,
     6  0.2155291348139220D+00,0.3154896462551432D+00,
     7  -.5831878447088045D+00,0.5198139571250358D+00,
     8  -.2041057083483440D+00,-.2041057083483440D+00,
     9  0.5198139571250358D+00,-.5831878447088045D+00,
     *  0.3154896462551432D+00,0.2155291348139220D+00,
     1  -.7428216449057439D+00,0.8319248661696513D+00,
     2  -.1239389574083093D+00,-.1131666042007823D+01,
     3  0.1573937740717607D+01,0.5001765333574129D+00,
     4  -.3590384795527415D+01,-.3749987235770745D+00,
     5  0.9450767760517566D+01,0.1655130527542360D+02,
     6  0.1753292360601324D+02,0.8750739333342171D+01,
     7  -.1838490132337566D+01,-.3214786597448908D+01,
     8  0.1579282614376747D+01,0.7575556574213355D+00,
     9  -.1385579432765773D+01,0.7528040255096071D+00,
     *  0.1260594948390180D+00,-.6586399381256727D+00,
     1  0.7229457346489030D+00,-.4494654653830360D+00,
     2  0.3823027345047653D-01,0.3349702365225073D+00,
     3  -.5508083225254948D+00,0.5508083225254948D+00,
     4  -.3349702365225073D+00,-.3823027345047653D-01,
     5  0.4494654653830360D+00,-.7229457346489030D+00,
     6  0.6586399381256727D+00,-.1260594948390180D+00,
     7  -.7528040255096071D+00,0.1385579432765773D+01,
     8  -.7575556574213355D+00,-.1579282614376747D+01,
     9  0.3214786597448908D+01,0.1838490132337566D+01,
     *  -.8750739333342171D+01,-.1753292360601324D+02,
     1  0.1878687182952331D+02,0.7666002689136455D+01,
     2  -.3404477426512804D+01,-.2209856730259832D+01,
     3  0.2323191365675765D+01,-.4837381241374515D+00,
     4  -.8373054020058460D+00,0.1125423782265202D+01,
     5  -.7425714367097070D+00,0.1553013097638267D+00,
     6  0.3329763413075318D+00,-.5946548273262859D+00,
     7  0.6179486540399921D+00,-.4503501119583696D+00,
     8  0.1635133695097584D+00,0.1635133695097584D+00,
     9  -.4503501119583696D+00,0.6179486540399921D+00,
     *  -.5946548273262859D+00,0.3329763413075318D+00,
     1  0.1553013097638267D+00,-.7425714367097070D+00,
     2  0.1125423782265202D+01,-.8373054020058460D+00,
     3  -.4837381241374515D+00,0.2323191365675765D+01,
     4  -.2209856730259832D+01,-.3404477426512804D+01,
     5  0.7666002689136455D+01,0.1878687182952331D+02,
     6  0.1971963569121103D+02,0.6563981433876758D+01,
     7  -.4405018736560776D+01,-.1050082130209461D+01,
     8  0.2342386818353591D+01,-.1362546494666642D+01,
     9  0.8744901545309482D-01,0.6939666994407514D+00,
     *  -.9197098985968038D+00,0.7724717214119335D+00,
     1  -.4480865237543868D+00,0.8662871779375656D-01,
     2  0.2287204902155826D+00,-.4546432850016119D+00,
     3  0.5712724672858390D+00,-.5712724672858390D+00,
     4  0.4546432850016119D+00,-.2287204902155826D+00,
     5  -.8662871779375656D-01,0.4480865237543868D+00,
     6  -.7724717214119335D+00,0.9197098985968038D+00,
     7  -.6939666994407514D+00,-.8744901545309482D-01,
     8  0.1362546494666642D+01,-.2342386818353591D+01,
     9  0.1050082130209461D+01,0.4405018736560776D+01,
     *  -.6563981433876758D+01,-.1971963569121103D+02,
     1  0.2067942779716816D+02,0.4949859019806576D+01,
     2  -.5181702549878119D+01,0.5043419076213845D+00,
     3  0.1669614825880018D+01,-.1773969848907003D+01,
     4  0.1065635192240751D+01,-.2880327125407381D+00,
     5  -.2804882023063036D+00,0.5972662349363789D+00,
     6  -.7058237790519812D+00,0.6663862791023850D+00,
     7  -.5309772322654065D+00,0.3382810368495153D+00,
     8  -.1158130899732296D+00,-.1158130899732296D+00,
     9  0.3382810368495153D+00,-.5309772322654065D+00,
     *  0.6663862791023850D+00,-.7058237790519812D+00,
     1  0.5972662349363789D+00,-.2804882023063036D+00,
     2  -.2880327125407381D+00,0.1065635192240751D+01,
     3  -.1773969848907003D+01,0.1669614825880018D+01,
     4  0.5043419076213845D+00,-.5181702549878119D+01,
     5  0.4949859019806576D+01,0.2067942779716816D+02,
     6  0.2129335618976210D+02,0.3494549488204937D+01,
     7  -.5397735247996340D+01,0.1690117919138163D+01,
     8  0.6921184395782759D+00,-.1470011448227040D+01,
     9  0.1384216112608755D+01,-.9856026933682936D+00,
     *  0.5464343239668447D+00,-.1702811919700702D+00,
     1  -.1183573202158310D+00,0.3253759090346138D+00,
     2  -.4648120382631614D+00,0.5499156420394142D+00,
     3  -.5900889803954503D+00,0.5900889803954503D+00,
     4  -.5499156420394142D+00,0.4648120382631614D+00,
     5  -.3253759090346138D+00,0.1183573202158310D+00,
     6  0.1702811919700702D+00,-.5464343239668447D+00,
     7  0.9856026933682936D+00,-.1384216112608755D+01,
     8  0.1470011448227040D+01,-.6921184395782759D+00,
     9  -.1690117919138163D+01,0.5397735247996340D+01,
     *  -.3494549488204937D+01,-.2129335618976210D+02/
c 
c        following is block number   5 of length 180 elements
c 
        data u5/
     1  0.2196092609771872D+02,0.1595478314491361D+01,
     2  -.5168300141959713D+01,0.2788107634919375D+01,
     3  -.6391477270421260D+00,-.5121573059338836D+00,
     4  0.9623228592445058D+00,-.1041801343785127D+01,
     5  0.9542823305898188D+00,-.8044861619839222D+00,
     6  0.6405308867560212D+00,-.4817366692518809D+00,
     7  0.3336817868268830D+00,-.1957038885801565D+00,
     8  0.6445673205078993D-01,0.6445673205078993D-01,
     9  -.1957038885801565D+00,0.3336817868268830D+00,
     *  -.4817366692518809D+00,0.6405308867560212D+00,
     1  -.8044861619839222D+00,0.9542823305898188D+00,
     2  -.1041801343785127D+01,0.9623228592445058D+00,
     3  -.5121573059338836D+00,-.6391477270421260D+00,
     4  0.2788107634919375D+01,-.5168300141959713D+01,
     5  0.1595478314491361D+01,0.2196092609771872D+02,
     6  -.2240653442859912D+02,-.1624413348135317D-01,
     7  0.4596754463601257D+01,-.3277230758630067D+01,
     8  0.1597186459567463D+01,-.4908179104346513D+00,
     9  -.1167643893425987D+00,0.4195200792755927D+00,
     *  -.5570893952640744D+00,0.6108369163701363D+00,
     1  -.6246559240100440D+00,0.6214386525847915D+00,
     2  -.6129066959421493D+00,0.6049511836631550D+00,
     3  -.6004046030991383D+00,0.6004046030991383D+00,
     4  -.6049511836631550D+00,0.6129066959421493D+00,
     5  -.6214386525847915D+00,0.6246559240100440D+00,
     6  -.6108369163701363D+00,0.5570893952640744D+00,
     7  -.4195200792755927D+00,0.1167643893425987D+00,
     8  0.4908179104346513D+00,-.1597186459567463D+01,
     9  0.3277230758630067D+01,-.4596754463601257D+01,
     *  0.1624413348135317D-01,0.2240653442859912D+02,
     1  -.2272621618186586D+02,0.1924096986934617D+01,
     2  0.3475177950876919D+01,-.3279096462207786D+01,
     3  0.2265412581566559D+01,-.1439748119036076D+01,
     4  0.8968139169445869D+00,-.5605190058159833D+00,
     5  0.3544316779078104D+00,-.2268336723139394D+00,
     6  0.1461032623641629D+00,-.9339917493413963D-01,
     7  0.5740824843508432D-01,-.3111014009906238D-01,
     8  0.9849853253947217D-02,0.9849853253947217D-02,
     9  -.3111014009906238D-01,0.5740824843508432D-01,
     *  -.9339917493413963D-01,0.1461032623641629D+00,
     1  -.2268336723139394D+00,0.3544316779078104D+00,
     2  -.5605190058159833D+00,0.8968139169445869D+00,
     3  -.1439748119036076D+01,0.2265412581566559D+01,
     4  -.3279096462207786D+01,0.3475177950876919D+01,
     5  0.1924096986934617D+01,-.2272621618186586D+02,
     6  -.2288123748340620D+02,0.3446943408052747D+01,
     7  0.2303133848581201D+01,-.2805768385552570D+01,
     8  0.2290447157095413D+01,-.1750653652775492D+01,
     9  0.1355593781124386D+01,-.1090153212630973D+01,
     *  0.9149789343290634D+00,-.7988725612694772D+00,
     1  0.7212770454411182D+00,-.6694181186395794D+00,
     2  0.6355425237637126D+00,-.6150518067727117D+00,
     3  0.6053803974276134D+00,-.6053803974276134D+00,
     4  0.6150518067727117D+00,-.6355425237637126D+00,
     5  0.6694181186395794D+00,-.7212770454411182D+00,
     6  0.7988725612694772D+00,-.9149789343290634D+00,
     7  0.1090153212630973D+01,-.1355593781124386D+01,
     8  0.1750653652775492D+01,-.2290447157095413D+01,
     9  0.2805768385552570D+01,-.2303133848581201D+01,
     *  -.3446943408052747D+01,0.2288123748340620D+02,
     1  0.2286335935781852D+02,-.5173667794137898D+01,
     2  -.6846367071036010D+00,0.1759952982968387D+01,
     3  -.1686305048351041D+01,0.1405156260732441D+01,
     4  -.1136469634391856D+01,0.9151130299926191D+00,
     5  -.7350900453818543D+00,0.5852144412725665D+00,
     6  -.4564126994689750D+00,0.3422182812650961D+00,
     7  -.2380211350894960D+00,0.1403510271972194D+00,
     8  -.4638320310944373D-01,-.4638320310944373D-01,
     9  0.1403510271972194D+00,-.2380211350894960D+00,
     *  0.3422182812650961D+00,-.4564126994689750D+00,
     1  0.5852144412725665D+00,-.7350900453818543D+00,
     2  0.9151130299926191D+00,-.1136469634391856D+01,
     3  0.1405156260732441D+01,-.1686305048351041D+01,
     4  0.1759952982968387D+01,-.6846367071036010D+00,
     5  -.5173667794137898D+01,0.2286335935781852D+02,
     6  -.2272004223821226D+02,0.6434878818404049D+01,
     7  -.6948280910718776D+00,-.6364897713112289D+00,
     8  0.7789786876735946D+00,-.6187466356764911D+00,
     9  0.3979075051947299D+00,-.1811348750715291D+00,
     *  -.1424403139987590D-01,0.1832315202911722D+00,
     1  -.3241642497852856D+00,0.4365234011572063D+00,
     2  -.5203420605323301D+00,0.5759200402163964D+00,
     3  -.6036004309760162D+00,0.6036004309760162D+00,
     4  -.5759200402163964D+00,0.5203420605323301D+00,
     5  -.4365234011572063D+00,0.3241642497852856D+00,
     6  -.1832315202911722D+00,0.1424403139987591D-01,
     7  0.1811348750715291D+00,-.3979075051947299D+00,
     8  0.6187466356764911D+00,-.7789786876735946D+00,
     9  0.6364897713112289D+00,0.6948280910718776D+00,
     *  -.6434878818404049D+01,0.2272004223821226D+02/
c 
c        following is the matrix v, with npts= 30
c 
c        following is block number   1 of length 180 elements
c 
        data v1/
     1  0.1277981643769441D-03,-.1590403705157395D-03,
     2  0.2357266525954249D-03,-.3115675543867951D-03,
     3  0.3409916607860436D-03,-.3596432632471868D-03,
     4  0.3137985029897685D-03,-.2659965904367070D-03,
     5  0.1971325422042324D-03,-.1963517237356980D-03,
     6  -.2566691641333348D-03,-.4445084307068533D-03,
     7  0.6627392214703927D-03,-.9584421605430728D-03,
     8  0.1238891490679367D-02,-.1406858727247356D-02,
     9  0.1582819038619263D-02,-.1661299816598187D-02,
     *  0.1836973764490119D-02,0.2183846046660424D-02,
     1  0.2344467444098586D-02,0.9098968094865845D-03,
     2  0.2756993572908802D-03,0.7157419900124440D-02,
     3  0.1256458291287097D-01,0.4781905323137212D-02,
     4  0.2140491517215169D-01,-.1489218492292852D-02,
     5  0.2567726503657701D-01,-.1358895281995831D-01,
     6  0.8247686173215908D-03,-.1026675394024209D-02,
     7  0.1518331446532427D-02,-.2003410305895294D-02,
     8  0.2186430452896011D-02,-.2298878180297393D-02,
     9  0.2000572893856299D-02,-.1690309961243810D-02,
     *  0.1251611300294515D-02,-.1242773193864306D-02,
     1  -.1626773602864815D-02,-.2782039146485499D-02,
     2  0.4064118392923898D-02,-.5604304300275790D-02,
     3  0.6843977231771675D-02,-.7177870542817420D-02,
     4  0.7361159627911900D-02,-.7205358447418939D-02,
     5  0.6913101513073936D-02,0.5828004810796818D-02,
     6  0.5000525772643877D-02,0.7965613920361793D-02,
     7  0.8543374950145372D-02,-.9257537319064509D-02,
     8  -.2346465430364462D-01,-.1842008895971291D-01,
     9  -.5776583792510076D-01,0.1665830680902771D-03,
     *  -.6157043004789963D-01,0.3293878841630409D-01,
     1  0.2838483026435497D-02,-.3535597114347611D-02,
     2  0.5182684453278354D-02,-.6789964062947048D-02,
     3  0.7320571282077490D-02,-.7595171324802184D-02,
     4  0.6532685930425193D-02,-.5444631623196677D-02,
     5  0.4011062027006479D-02,-.3935803540506831D-02,
     6  -.5165713021008028D-02,-.8470647035007139D-02,
     7  0.1162081042510035D-01,-.1383991501886005D-01,
     8  0.1417508274083656D-01,-.1140480141852545D-01,
     9  0.7666402346705793D-02,-.3544128816552965D-02,
     *  -.1003599395885032D-02,-.3974368222682438D-02,
     1  -.7863811993459364D-02,-.1490606853872479D-01,
     2  -.1880226987350207D-01,0.1037431204283005D-02,
     3  0.1655141755150430D-01,0.3506591962474145D-01,
     4  0.8133518465850812D-01,0.7635141351487709D-02,
     5  0.7138483378166335D-01,-.3902894412084513D-01,
     6  0.7075594708288464D-02,-.8819956724756269D-02,
     7  0.1262145907068896D-01,-.1620279940107169D-01,
     8  0.1684566510577066D-01,-.1678028238230447D-01,
     9  0.1389969703691252D-01,-.1108950031595130D-01,
     *  0.8015231814014897D-02,-.7582376152319207D-02,
     1  -.9992009515098929D-02,-.1451630325097083D-01,
     2  0.1680826972531991D-01,-.1266723079436264D-01,
     3  0.4670341707466538D-02,0.6388671387709764D-02,
     4  -.1579190422765766D-01,0.2072518297778072D-01,
     5  -.2248586722536560D-01,-.2079643297753054D-01,
     6  -.1448989769691356D-01,-.2793000538933416D-02,
     7  0.8646669561865471D-02,-.4004233971800675D-02,
     8  -.1084507422347944D-01,-.4555585565686930D-01,
     9  -.9302564848783165D-01,-.1755236535753770D-01,
     *  -.6562340549213800D-01,0.3706875632561275D-01,
     1  0.1428684450328887D-01,-.1780753801989289D-01,
     2  0.2417851505432073D-01,-.2960811523658210D-01,
     3  0.2808755773536034D-01,-.2502393512764448D-01,
     4  0.1840450006810487D-01,-.1265472973967272D-01,
     5  0.8422912128886412D-02,-.6934187455311209D-02,
     6  -.9107963535620145D-02,-.7579529745586635D-02,
     7  0.8126079714785931D-03,0.1603025762556734D-01,
     8  -.2951631410400395D-01,0.3358593503001794D-01,
     9  -.2659818213728298D-01,0.1223947215004901D-01,
     *  0.6380915774705383D-02,0.2084699461239980D-01,
     1  0.3038781630887894D-01,0.2671662892523829D-01,
     2  0.1683533975418074D-01,0.2241907894137389D-01,
     3  0.1646767374470205D-01,0.4531204286739569D-01,
     4  0.9530738248909595D-01,0.2619235276785705D-01,
     5  0.5402352758445220D-01,-.3193627337679837D-01,
     6  0.2467506884925503D-01,-.3067029214694104D-01,
     7  0.3768177913028111D-01,-.4180755037187636D-01,
     8  0.3148946924899179D-01,-.1919584887511735D-01,
     9  0.6658359513788769D-02,0.1785317760023479D-02,
     *  -.3865118871165516D-02,0.6465160933906952D-02,
     1  0.8894112781636309D-02,0.2422534831862316D-01,
     2  -.3917164575297387D-01,0.4536846042468085D-01,
     3  -.3267406374548426D-01,0.2310835164972578D-02,
     4  0.2680464986659163D-01,-.4052403322055983D-01,
     5  0.3591728485332131D-01,0.1680611492501035D-01,
     6  -.1143588889369994D-01,-.2809646612235629D-01,
     7  -.3669285370458688D-01,-.4570319258316545D-01,
     8  -.3227612341126312D-01,-.3546381570608098D-01,
     9  -.9132596860336032D-01,-.3228104322811189D-01,
     *  -.4229676367361372D-01,0.2655113519066795D-01/
  
c 
c        following is block number   2 of length 180 elements
c 
        data v2/
     1  0.3759011514134179D-01,-.4633476362128677D-01,
     2  0.4775950167294628D-01,-.4308668006177719D-01,
     3  0.1373991262770181D-01,0.1358913116157937D-01,
     4  -.2988068263075061D-01,0.3435165748162366D-01,
     5  -.2805620809376530D-01,0.2794172405243956D-01,
     6  0.3668303068642019D-01,0.4907687721190121D-01,
     7  -.4374014510282609D-01,0.4856792731414742D-02,
     8  0.3309546233747531D-01,-.5065440008570817D-01,
     9  0.3224815613659058D-01,0.5474145313208826D-02,
     *  -.4101956151477282D-01,-.4978820619103953D-01,
     1  -.3000474821182917D-01,0.5931266584800616D-03,
     2  0.3562852732654711D-01,0.6115844197805015D-01,
     3  0.5176244865963128D-01,0.1977277975181757D-01,
     4  0.8365918924846927D-01,0.3598954553278815D-01,
     5  0.3244123761080602D-01,-.2198553632185234D-01,
     6  0.5150649262009227D-01,-.6235729439521827D-01,
     7  0.4748654588161240D-01,-.2477536553127891D-01,
     8  -.2816300978492590D-01,0.6269139265345844D-01,
     9  -.6730986130187222D-01,0.5154486544191580D-01,
     *  -.3270875292243344D-01,0.2308861278206864D-01,
     1  0.2569053560650924D-01,0.5969283888939678D-03,
     2  0.3383605852782787D-01,-.6119916975609832D-01,
     3  0.4155972777866541D-01,0.1338491509007619D-01,
     4  -.5500261467672810D-01,0.5152345876385648D-01,
     5  -.6516493887235375D-02,0.3939235435133772D-01,
     6  0.5915322822169327D-01,0.3889702924580903D-01,
     7  -.1286666264614286D-01,-.6212142800628707D-01,
     8  -.6877128748581333D-01,-.1914109243114016D-02,
     9  -.7411332124591922D-01,-.3796093073083818D-01,
     *  -.2471106776122992D-01,0.1842888443833259D-01,
     1  0.6439169609959143D-01,-.7540399859223750D-01,
     2  0.3238408674195226D-01,0.1214659348690137D-01,
     3  -.7356619168975347D-01,0.8487417186130280D-01,
     4  -.5059730378164611D-01,0.3684441266252483D-02,
     5  0.1812544524234618D-01,-.3104402188988836D-01,
     6  -.4600788899595975D-01,-.7550851156359015D-01,
     7  0.6713016384080006D-01,-.2133103616767166D-02,
     8  -.5316815352857018D-01,0.5630023764516042D-01,
     9  0.1009887440521289D-02,-.5742281403496506D-01,
     *  0.5961271405712012D-01,0.9339648491307860D-02,
     1  -.5303308299029138D-01,-.6795914292008272D-01,
     2  -.2193821210977724D-01,0.4857202512175916D-01,
     3  0.7956694793499375D-01,-.1554267987577986D-01,
     4  0.6381047180666493D-01,0.3882441543830057D-01,
     5  0.1877055548441293D-01,-.1575968066949234D-01,
     6  0.7437230155592821D-01,-.8228680516608744D-01,
     7  0.3959865611212229D-02,0.5439921127648779D-01,
     8  -.8939510669593481D-01,0.4521669716586792D-01,
     9  0.3522221245148584D-01,-.8571214742749667D-01,
     *  0.7896642562835262D-01,-.6395974101535351D-01,
     1  -.6371919943211609D-01,-.1599357078280348D-01,
     2  -.4410020013880195D-01,0.7157930618879137D-01,
     3  -.1944181169902415D-01,-.5683817039548276D-01,
     4  0.6396490968119110D-01,0.1836556630047985D-02,
     5  -.6812966556135027D-01,-.6069689851423918D-01,
     6  0.1409668993193389D-01,0.7272619233307375D-01,
     7  0.5621845250063939D-01,-.2451882894236281D-01,
     8  -.8274408587698975D-01,0.3114072945810003D-01,
     9  -.5337448039697261D-01,-.3905018042768453D-01,
     *  -.1418327906099007D-01,0.1379268003155014D-01,
     1  0.8040233860725186D-01,-.8109862448231184D-01,
     2  -.3064656375119694D-01,0.8453536392479195D-01,
     3  -.6231978629627313D-01,-.3163926578620179D-01,
     4  0.1108119457063973D+00,-.9950987392509190D-01,
     5  0.3607332650959111D-01,0.2011100401121371D-01,
     6  0.5813644380640645D-01,0.9756761123509447D-01,
     7  -.6277022864202332D-01,-.3429919184455853D-01,
     8  0.7187819021827268D-01,-.1466392887220320D-01,
     9  -.6637439515794320D-01,0.6328984509054264D-01,
     *  0.2403244473242316D-01,0.8090445335313889D-01,
     1  0.3726095008243812D-01,-.5172288588194491D-01,
     2  -.7983208610956540D-01,-.4666989086039772D-02,
     3  0.7846609328558105D-01,-.4421026255450744D-01,
     4  0.4310803450869109D-01,0.3895043056453008D-01,
     5  0.1057001357035566D-01,-.1236273569537638D-01,
     6  0.8259474438856641D-01,-.7176688633361256D-01,
     7  -.6289676091118193D-01,0.9229606730028760D-01,
     8  -.1048438594871452D-01,-.8573136803732216D-01,
     9  0.9157479105208666D-01,0.1052995660270772D-01,
     *  -.9770799756973687D-01,0.1213701727052890D+00,
     1  0.9694185946902903D-01,-.1517885898007717D-02,
     2  0.7557230311984693D-01,-.5487272941993402D-01,
     3  -.3784405048298083D-01,0.7618254509889544D-01,
     4  0.4121314838338816D-02,-.8232092440753820D-01,
     5  0.4091875992352944D-01,-.5852080895148512D-01,
     6  -.7719384362895191D-01,0.1287777297029540D-01,
     7  0.8721675216697246D-01,0.3390596824704145D-01,
     8  -.6774595791691759D-01,0.5454194700986241D-01,
     9  -.3312279508604611D-01,-.3871971118643502D-01,
     *  -.7637467957655532D-02,0.1134376519389064D-01/
  
c 
c        following is block number   3 of length 180 elements
c 
        data v3/
     1  0.8207102182899553D-01,-.5571962708201524D-01,
     2  -.8718027872144364D-01,0.7865280261209546D-01,
     3  0.3992317879218299D-01,-.9208499586463285D-01,
     4  0.4479433304722521D-02,0.1206174558005112D+00,
     5  -.1370494228307044D+00,0.2726098050396384D-01,
     6  -.8814370385656740D-01,-.1123167638349943D+00,
     7  0.1905931574788104D-01,0.7776443006050010D-01,
     8  -.4037627099058021D-01,-.6360995932259201D-01,
     9  0.6524768447943799D-01,0.4085380338793891D-01,
     *  -.8438167495175443D-01,0.5758705481122377D-02,
     1  0.8946248785097568D-01,0.3150850180384178D-01,
     2  -.7737580704436400D-01,-.5913041638512950D-01,
     3  0.5196207560766796D-01,-.6215267822648939D-01,
     4  0.2342294692883558D-01,0.3847441446879099D-01,
     5  0.5166238974098463D-02,-.1064647345079029D-01,
     6  0.8044835199093310D-01,-.3507825081224286D-01,
     7  -.1020637188982216D+00,0.5123349039903803D-01,
     8  0.7564077672452713D-01,-.6656303072446115D-01,
     9  -.6655162870033364D-01,0.1166133091643111D+00,
     *  -.1647576170208659D-02,-.1649525189084934D+00,
     1  -.1200617591179100D+00,0.4913026133397624D-01,
     2  -.8931436220046571D-01,-.7767079209456311D-02,
     3  0.7798733303836547D-01,-.9934783866859496D-02,
     4  -.8357935444241633D-01,0.3082455252923807D-01,
     5  0.8078487993159309D-01,0.5183477212407461D-01,
     6  -.7002004542307683D-01,-.6910922353180870D-01,
     7  0.5289467899573298D-01,0.7744285255245627D-01,
     8  -.3259209677784270D-01,0.6714103899976156D-01,
     9  -.1395475290702212D-01,-.3828158163510529D-01,
     *  -.2990499853115919D-02,0.1021134091877940D-01,
     1  0.7920955036331044D-01,-.1194538154223972D-01,
     2  -.1087665317996457D+00,0.1760432084264359D-01,
     3  0.9373785696698937D-01,-.2441804745925145D-01,
     4  -.9358565881011606D-01,0.4120575238708103D-01,
     5  0.1527755722131697D+00,-.1200149381920604D+00,
     6  0.1180844975615403D+00,0.9682118571830679D-01,
     7  0.4952300774689208D-01,-.7256580802277535D-01,
     8  -.3843037471922214D-01,0.7613592123735488D-01,
     9  0.3733955217307740D-01,-.8399273341321167D-01,
     *  -.3286021116225620D-01,-.8845687632513103D-01,
     1  0.2633515448066852D-01,0.9044349799882367D-01,
     2  -.1875160312001496D-01,-.8703056214407745D-01,
     3  0.1109804029790578D-01,-.6960614701027802D-01,
     4  0.4634679839850164D-02,0.3817731874426382D-01,
     5  0.9792418042771173D-03,-.1000226957519711D-01,
     6  0.7920955036331044D-01,0.1194538154223972D-01,
     7  -.1087665317996457D+00,-.1760432084264359D-01,
     8  0.9373785696698937D-01,0.2441804745925145D-01,
     9  -.9358565881011606D-01,-.4120575238708103D-01,
     *  0.1527755722131697D+00,0.1200149381920604D+00,
     1  0.1180844975615403D+00,-.9682118571830679D-01,
     2  0.4952300774689208D-01,0.7256580802277535D-01,
     3  -.3843037471922214D-01,-.7613592123735488D-01,
     4  0.3733955217307740D-01,0.8399273341321167D-01,
     5  -.3286021116225620D-01,0.8845687632513103D-01,
     6  0.2633515448066852D-01,-.9044349799882367D-01,
     7  -.1875160312001496D-01,0.8703056214407745D-01,
     8  0.1109804029790578D-01,0.6960614701027802D-01,
     9  0.4634679839850164D-02,-.3817731874426382D-01,
     *  0.9792418042771174D-03,0.1000226957519711D-01,
     1  0.8044835199093310D-01,0.3507825081224286D-01,
     2  -.1020637188982216D+00,-.5123349039903803D-01,
     3  0.7564077672452713D-01,0.6656303072446115D-01,
     4  -.6655162870033364D-01,-.1166133091643111D+00,
     5  -.1647576170208659D-02,0.1649525189084934D+00,
     6  -.1200617591179100D+00,-.4913026133397624D-01,
     7  -.8931436220046571D-01,0.7767079209456311D-02,
     8  0.7798733303836547D-01,0.9934783866859496D-02,
     9  -.8357935444241633D-01,-.3082455252923807D-01,
     *  0.8078487993159309D-01,-.5183477212407461D-01,
     1  -.7002004542307683D-01,0.6910922353180870D-01,
     2  0.5289467899573298D-01,-.7744285255245627D-01,
     3  -.3259209677784270D-01,-.6714103899976156D-01,
     4  -.1395475290702212D-01,0.3828158163510529D-01,
     5  -.2990499853115920D-02,-.1021134091877940D-01,
     6  0.8207102182899553D-01,0.5571962708201524D-01,
     7  -.8718027872144364D-01,-.7865280261209546D-01,
     8  0.3992317879218299D-01,0.9208499586463285D-01,
     9  0.4479433304722521D-02,-.1206174558005112D+00,
     *  -.1370494228307044D+00,-.2726098050396384D-01,
     1  -.8814370385656740D-01,0.1123167638349943D+00,
     2  0.1905931574788104D-01,-.7776443006050010D-01,
     3  -.4037627099058021D-01,0.6360995932259201D-01,
     4  0.6524768447943799D-01,-.4085380338793891D-01,
     5  -.8438167495175443D-01,-.5758705481122377D-02,
     6  0.8946248785097568D-01,-.3150850180384178D-01,
     7  -.7737580704436400D-01,0.5913041638512950D-01,
     8  0.5196207560766796D-01,0.6215267822648939D-01,
     9  0.2342294692883558D-01,-.3847441446879099D-01,
     *  0.5166238974098463D-02,0.1064647345079029D-01/
c 
c        following is block number   4 of length 180 elements
c 
        data v4/
     1  0.8259474438856641D-01,0.7176688633361256D-01,
     2  -.6289676091118193D-01,-.9229606730028760D-01,
     3  -.1048438594871452D-01,0.8573136803732216D-01,
     4  0.9157479105208666D-01,-.1052995660270772D-01,
     5  -.9770799756973687D-01,-.1213701727052890D+00,
     6  0.9694185946902903D-01,0.1517885898007717D-02,
     7  0.7557230311984693D-01,0.5487272941993402D-01,
     8  -.3784405048298083D-01,-.7618254509889544D-01,
     9  0.4121314838338816D-02,0.8232092440753820D-01,
     *  0.4091875992352944D-01,0.5852080895148512D-01,
     1  -.7719384362895191D-01,-.1287777297029540D-01,
     2  0.8721675216697246D-01,-.3390596824704145D-01,
     3  -.6774595791691759D-01,-.5454194700986241D-01,
     4  -.3312279508604611D-01,0.3871971118643502D-01,
     5  -.7637467957655532D-02,-.1134376519389064D-01,
     6  0.8040233860725186D-01,0.8109862448231184D-01,
     7  -.3064656375119694D-01,-.8453536392479195D-01,
     8  -.6231978629627313D-01,0.3163926578620179D-01,
     9  0.1108119457063973D+00,0.9950987392509190D-01,
     *  0.3607332650959111D-01,-.2011100401121371D-01,
     1  0.5813644380640645D-01,-.9756761123509447D-01,
     2  -.6277022864202332D-01,0.3429919184455853D-01,
     3  0.7187819021827268D-01,0.1466392887220320D-01,
     4  -.6637439515794320D-01,-.6328984509054264D-01,
     5  0.2403244473242316D-01,-.8090445335313889D-01,
     6  0.3726095008243812D-01,0.5172288588194491D-01,
     7  -.7983208610956540D-01,0.4666989086039772D-02,
     8  0.7846609328558105D-01,0.4421026255450744D-01,
     9  0.4310803450869109D-01,-.3895043056453008D-01,
     *  0.1057001357035566D-01,0.1236273569537638D-01,
     1  0.7437230155592821D-01,0.8228680516608744D-01,
     2  0.3959865611212229D-02,-.5439921127648779D-01,
     3  -.8939510669593481D-01,-.4521669716586792D-01,
     4  0.3522221245148584D-01,0.8571214742749667D-01,
     5  0.7896642562835262D-01,0.6395974101535351D-01,
     6  -.6371919943211609D-01,0.1599357078280348D-01,
     7  -.4410020013880195D-01,-.7157930618879137D-01,
     8  -.1944181169902415D-01,0.5683817039548276D-01,
     9  0.6396490968119110D-01,-.1836556630047985D-02,
     *  -.6812966556135027D-01,0.6069689851423918D-01,
     1  0.1409668993193389D-01,-.7272619233307375D-01,
     2  0.5621845250063939D-01,0.2451882894236281D-01,
     3  -.8274408587698975D-01,-.3114072945810003D-01,
     4  -.5337448039697261D-01,0.3905018042768453D-01,
     5  -.1418327906099007D-01,-.1379268003155014D-01,
     6  0.6439169609959143D-01,0.7540399859223750D-01,
     7  0.3238408674195226D-01,-.1214659348690137D-01,
     8  -.7356619168975347D-01,-.8487417186130280D-01,
     9  -.5059730378164611D-01,-.3684441266252483D-02,
     *  0.1812544524234618D-01,0.3104402188988836D-01,
     1  -.4600788899595975D-01,0.7550851156359015D-01,
     2  0.6713016384080006D-01,0.2133103616767166D-02,
     3  -.5316815352857018D-01,-.5630023764516042D-01,
     4  0.1009887440521289D-02,0.5742281403496506D-01,
     5  0.5961271405712012D-01,-.9339648491307860D-02,
     6  -.5303308299029138D-01,0.6795914292008272D-01,
     7  -.2193821210977724D-01,-.4857202512175916D-01,
     8  0.7956694793499375D-01,0.1554267987577986D-01,
     9  0.6381047180666493D-01,-.3882441543830057D-01,
     *  0.1877055548441293D-01,0.1575968066949234D-01,
     1  0.5150649262009227D-01,0.6235729439521827D-01,
     2  0.4748654588161240D-01,0.2477536553127891D-01,
     3  -.2816300978492590D-01,-.6269139265345844D-01,
     4  -.6730986130187222D-01,-.5154486544191580D-01,
     5  -.3270875292243344D-01,-.2308861278206864D-01,
     6  0.2569053560650924D-01,-.5969283888939678D-03,
     7  0.3383605852782787D-01,0.6119916975609832D-01,
     8  0.4155972777866541D-01,-.1338491509007619D-01,
     9  -.5500261467672810D-01,-.5152345876385648D-01,
     *  -.6516493887235375D-02,-.3939235435133772D-01,
     1  0.5915322822169327D-01,-.3889702924580903D-01,
     2  -.1286666264614286D-01,0.6212142800628707D-01,
     3  -.6877128748581333D-01,0.1914109243114016D-02,
     4  -.7411332124591922D-01,0.3796093073083818D-01,
     5  -.2471106776122992D-01,-.1842888443833259D-01,
     6  0.3759011514134179D-01,0.4633476362128677D-01,
     7  0.4775950167294628D-01,0.4308668006177719D-01,
     8  0.1373991262770181D-01,-.1358913116157937D-01,
     9  -.2988068263075061D-01,-.3435165748162366D-01,
     *  -.2805620809376530D-01,-.2794172405243956D-01,
     1  0.3668303068642019D-01,-.4907687721190121D-01,
     2  -.4374014510282609D-01,-.4856792731414742D-02,
     3  0.3309546233747531D-01,0.5065440008570817D-01,
     4  0.3224815613659058D-01,-.5474145313208826D-02,
     5  -.4101956151477282D-01,0.4978820619103953D-01,
     6  -.3000474821182917D-01,-.5931266584800616D-03,
     7  0.3562852732654711D-01,-.6115844197805015D-01,
     8  0.5176244865963128D-01,-.1977277975181757D-01,
     9  0.8365918924846927D-01,-.3598954553278815D-01,
     *  0.3244123761080602D-01,0.2198553632185234D-01/
c 
c        following is block number   5 of length 180 elements
c 
        data v5/
     1  0.2467506884925503D-01,0.3067029214694104D-01,
     2  0.3768177913028111D-01,0.4180755037187636D-01,
     3  0.3148946924899179D-01,0.1919584887511735D-01,
     4  0.6658359513788769D-02,-.1785317760023479D-02,
     5  -.3865118871165516D-02,-.6465160933906952D-02,
     6  0.8894112781636309D-02,-.2422534831862316D-01,
     7  -.3917164575297387D-01,-.4536846042468085D-01,
     8  -.3267406374548426D-01,-.2310835164972578D-02,
     9  0.2680464986659163D-01,0.4052403322055983D-01,
     *  0.3591728485332131D-01,-.1680611492501035D-01,
     1  -.1143588889369994D-01,0.2809646612235629D-01,
     2  -.3669285370458688D-01,0.4570319258316545D-01,
     3  -.3227612341126312D-01,0.3546381570608098D-01,
     4  -.9132596860336032D-01,0.3228104322811189D-01,
     5  -.4229676367361372D-01,-.2655113519066795D-01,
     6  0.1428684450328887D-01,0.1780753801989289D-01,
     7  0.2417851505432073D-01,0.2960811523658210D-01,
     8  0.2808755773536034D-01,0.2502393512764448D-01,
     9  0.1840450006810487D-01,0.1265472973967272D-01,
     *  0.8422912128886412D-02,0.6934187455311209D-02,
     1  -.9107963535620145D-02,0.7579529745586635D-02,
     2  0.8126079714785931D-03,-.1603025762556734D-01,
     3  -.2951631410400395D-01,-.3358593503001794D-01,
     4  -.2659818213728298D-01,-.1223947215004901D-01,
     5  0.6380915774705383D-02,-.2084699461239980D-01,
     6  0.3038781630887894D-01,-.2671662892523829D-01,
     7  0.1683533975418074D-01,-.2241907894137389D-01,
     8  0.1646767374470205D-01,-.4531204286739569D-01,
     9  0.9530738248909595D-01,-.2619235276785705D-01,
     *  0.5402352758445220D-01,0.3193627337679837D-01,
     1  0.7075594708288464D-02,0.8819956724756269D-02,
     2  0.1262145907068896D-01,0.1620279940107169D-01,
     3  0.1684566510577066D-01,0.1678028238230447D-01,
     4  0.1389969703691252D-01,0.1108950031595130D-01,
     5  0.8015231814014897D-02,0.7582376152319207D-02,
     6  -.9992009515098929D-02,0.1451630325097083D-01,
     7  0.1680826972531991D-01,0.1266723079436264D-01,
     8  0.4670341707466538D-02,-.6388671387709764D-02,
     9  -.1579190422765766D-01,-.2072518297778072D-01,
     *  -.2248586722536560D-01,0.2079643297753054D-01,
     1  -.1448989769691356D-01,0.2793000538933416D-02,
     2  0.8646669561865471D-02,0.4004233971800675D-02,
     3  -.1084507422347944D-01,0.4555585565686930D-01,
     4  -.9302564848783165D-01,0.1755236535753770D-01,
     5  -.6562340549213800D-01,-.3706875632561275D-01,
     6  0.2838483026435497D-02,0.3535597114347611D-02,
     7  0.5182684453278354D-02,0.6789964062947048D-02,
     8  0.7320571282077490D-02,0.7595171324802184D-02,
     9  0.6532685930425193D-02,0.5444631623196677D-02,
     *  0.4011062027006479D-02,0.3935803540506831D-02,
     1  -.5165713021008028D-02,0.8470647035007139D-02,
     2  0.1162081042510035D-01,0.1383991501886005D-01,
     3  0.1417508274083656D-01,0.1140480141852545D-01,
     4  0.7666402346705793D-02,0.3544128816552965D-02,
     5  -.1003599395885032D-02,0.3974368222682438D-02,
     6  -.7863811993459364D-02,0.1490606853872479D-01,
     7  -.1880226987350207D-01,-.1037431204283005D-02,
     8  0.1655141755150430D-01,-.3506591962474145D-01,
     9  0.8133518465850812D-01,-.7635141351487709D-02,
     *  0.7138483378166335D-01,0.3902894412084513D-01,
     1  0.8247686173215908D-03,0.1026675394024209D-02,
     2  0.1518331446532427D-02,0.2003410305895294D-02,
     3  0.2186430452896011D-02,0.2298878180297393D-02,
     4  0.2000572893856299D-02,0.1690309961243810D-02,
     5  0.1251611300294515D-02,0.1242773193864306D-02,
     6  -.1626773602864815D-02,0.2782039146485499D-02,
     7  0.4064118392923898D-02,0.5604304300275790D-02,
     8  0.6843977231771675D-02,0.7177870542817420D-02,
     9  0.7361159627911900D-02,0.7205358447418939D-02,
     *  0.6913101513073936D-02,-.5828004810796818D-02,
     1  0.5000525772643877D-02,-.7965613920361793D-02,
     2  0.8543374950145372D-02,0.9257537319064509D-02,
     3  -.2346465430364462D-01,0.1842008895971291D-01,
     4  -.5776583792510076D-01,-.1665830680902771D-03,
     5  -.6157043004789963D-01,-.3293878841630409D-01,
     6  0.1277981643769441D-03,0.1590403705157395D-03,
     7  0.2357266525954249D-03,0.3115675543867951D-03,
     8  0.3409916607860436D-03,0.3596432632471868D-03,
     9  0.3137985029897685D-03,0.2659965904367070D-03,
     *  0.1971325422042324D-03,0.1963517237356980D-03,
     1  -.2566691641333348D-03,0.4445084307068533D-03,
     2  0.6627392214703927D-03,0.9584421605430728D-03,
     3  0.1238891490679367D-02,0.1406858727247356D-02,
     4  0.1582819038619263D-02,0.1661299816598187D-02,
     5  0.1836973764490119D-02,-.2183846046660424D-02,
     6  0.2344467444098586D-02,-.9098968094865845D-03,
     7  0.2756993572908802D-03,-.7157419900124440D-02,
     8  0.1256458291287097D-01,-.4781905323137212D-02,
     9  0.2140491517215169D-01,0.1489218492292852D-02,
     *  0.2567726503657701D-01,0.1358895281995831D-01/
c 
c        return to the user the matrices u,v,vmatr
c 
        do 1200 i=1,900
c 
        uout(i)=u(i)
        vout(i)=v(i)
        vmatrout(i)=vmatr(i)
 1200 continue
c 
c       return to the user the nodes and weights
c 
        npts=15
        do 1400 i=1,npts
c 
        xsout(i)=-xs0(npts-i+1)
        wsout(i)=ws0(npts-i+1)
c 
        xsout(npts+i)=xs0(i)
        wsout(npts+i)=ws0(i)
 1400 continue
c 
        npts=30
c 
        return
c 
c 
c 
c 
        entry xswsget(xsout,wsout,npts)
c 
        npts=15
        do 2200 i=1,npts
c 
        xsout(i)=-xs0(npts-i+1)
        wsout(i)=ws0(npts-i+1)
c 
        xsout(npts+i)=xs0(i)
        wsout(npts+i)=ws0(i)
 2200 continue
c 
        npts=30
c 
        return
c 
c 
c 
c 
        entry uoutget(uout)
c 
        do 2400 i=1,900
c 
        uout(i)=u(i)
 2400 continue
        return
c 
c 
c 
c 
        entry voutget(vout)
c 
        do 2500 i=1,900
c 
        vout(i)=v(i)
 2500 continue
        return
c 
c 
c 
c 
        entry vmatrget(vmatrout)
  
        do 2600 i=1,900
c 
        vmatrout(i)=vmatr(i)
 2600 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine univquadrs(xs,ws)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension xs1(32),ws1(32),xs2(32),ws2(32),xs3(32),ws3(32),
     1      xs4(32),ws4(32),xs5(32),ws5(32),xs6(32),ws6(32),
     2      xs7(32),ws7(32),xs8(32),ws8(32),xs9(32),ws9(32),
     3      xs10(32),ws10(32),xs11(32),ws11(32),xs12(32),ws12(32),
     4      xs13(32),ws13(32),xs14(32),ws14(32),xs15(32),ws15(32)
c 
        dimension xs(32,30),ws(32,30)
c 
c        This subroutine returns to the user 32-point special-purpose
c        quadratures (30 of them things). The quadrature number k
c        (k=1,2,...,30) integrates functions of the form
c 
c        f(x)=phi(x)+psi(x)*log(|x-x_k|),                           (1)
c 
c        with phi and psi polynomials of order 19, and x_k one of the
c        support nodes returned by the subroutine matrret3 (or its entry
c        xswsget).
c 
c                           Input parameters: none
c 
c                           Output parameters:
c 
c  xs - a 32 \times 30 array containing nodes of all 30 quadratures
c  ws - a 32 \times 30 array containing weights of all 30 quadratures
c 
c     . . . In this run, istart=1, and x0=-.9999659071661138D+00
c 
        data xs1/
     1  -.9999926842982632125D+00,-.9999766107184137465D+00,
     2  -.9999321725494536845D+00,-.9994134334906660373D+00,
     3  -.9971251084617104933D+00,-.9905239941429529562D+00,
     4  -.9759265980073636367D+00,-.9491146357724889995D+00,
     5  -.9062376613961078048D+00,-.8448343418335435955D+00,
     6  -.7650219177149375243D+00,-.6708322694234459290D+00,
     7  -.5681569459225669555D+00,-.4542322375008505664D+00,
     8  -.3238023990973824142D+00,-.1995013786012745148D+00,
     9  -.1391390562732062486D+00,-.3285512507595560992D-01,
     *  0.6827298171831842140D-01,0.1985789712498753892D+00,
     1  0.3382359344042465329D+00,0.4556213066476002335D+00,
     2  0.5402833028263462196D+00,0.6047670551312823255D+00,
     3  0.6888766628816837367D+00,0.7821570876300537240D+00,
     4  0.8618931089552985711D+00,0.9218529155388627458D+00,
     5  0.9620236528486702717D+00,0.9852845898154829472D+00,
     6  0.9961235141202109686D+00,0.9995698393229103263D+00/
c 
        data ws1/
     1  0.1293861551669072486D-04,0.1060863062136854047D-04,
     2  0.1546232790456628010D-03,0.1089519590236273780D-02,
     3  0.3912928588124452466D-02,0.9913876325051163450D-02,
     4  0.1999952795720435423D-01,0.3428305475131167472D-01,
     5  0.5190015653559924495D-01,0.7091316423320499563D-01,
     6  0.8803978981560874051D-01,0.9910489773632747842D-01,
     7  0.1067745605254463961D+00,0.1224612988641804993D+00,
     8  0.1359456406546648880D+00,0.8919099608936808470D-01,
     9  0.7729377728580819472D-01,0.1057698660278751818D+00,
     *  0.1114453325913364674D+00,0.1418120889311425791D+00,
     1  0.1325151647209786955D+00,0.1001139216251070395D+00,
     2  0.7064916033043500520D-01,0.6884775697950519741D-01,
     3  0.9453803381143852470D-01,0.8841605051947676930D-01,
     4  0.7018632839806681629D-01,0.4977237018586029715D-01,
     5  0.3106747025917853467D-01,0.1621957733754255086D-01,
     6  0.6316596424322844762D-02,0.1328922380555788634D-02/
c 
c       In this run, istart=2, and x0=-.9996124417392792D+00
c 
        data xs2/
     1  -.9999773282033958522D+00,-.9998655853750659076D+00,
     2  -.9997058693514118137D+00,-.9995122919047986270D+00,
     3  -.9982899976860404237D+00,-.9942355421064399170D+00,
     4  -.9844518729198153261D+00,-.9651739045944195348D+00,
     5  -.9324086170462463014D+00,-.8827437237156106001D+00,
     6  -.8143432674583643038D+00,-.7290811301110036820D+00,
     7  -.6410836746884132008D+00,-.5632690653910359007D+00,
     8  -.4542472817293781781D+00,-.3199352400611844613D+00,
     9  -.1807503379090381661D+00,-.5282851417930418388D-01,
     *  0.4944478686825125086D-01,0.1289007507134699264D+00,
     1  0.2172555541082377379D+00,0.3118776667292348299D+00,
     2  0.4351825040002334627D+00,0.5653680080231950215D+00,
     3  0.6849458170304214714D+00,0.7864556539490111728D+00,
     4  0.8667304463929664625D+00,0.9253263990565100488D+00,
     5  0.9639700097723038575D+00,0.9861164021481364246D+00,
     6  0.9963583959909868862D+00,0.9995971175080441337D+00/
c 
        data ws2/
     1  0.6183700329476664857D-04,0.1521646120497065107D-03,
     2  0.1253597767796259899D-03,0.4702474176788674065D-03,
     3  0.2251425816220062248D-02,0.6346543530681097040D-02,
     4  0.1385931404387296541D-01,0.2537574501952319647D-01,
     5  0.4073892566560228949D-01,0.5892317867880716756D-01,
     6  0.7767302112460160019D-01,0.9087836480720680709D-01,
     7  0.7964828844472402408D-01,0.8827611992534590162D-01,
     8  0.1259040516438986926D+00,0.1395503012396125352D+00,
     9  0.1360769523435932704D+00,0.1176302531157546649D+00,
     *  0.8491163488989163113D-01,0.8428534761274498698D-01,
     1  0.8700990544563570420D-01,0.1099021185637406626D+00,
     2  0.1310372621327203298D+00,0.1266485673831039049D+00,
     3  0.1113429653950508572D+00,0.9116858548789565127D-01,
     4  0.6931471772336547353D-01,0.4817179676479204346D-01,
     5  0.2970091904657982130D-01,0.1537121543551382851D-01,
     6  0.5947313359041419216D-02,0.1245556550675120038D-02/
c 
c        In this run, istart=3, and x0=-.9981286463685695D+00
c 
        data xs3/
     1  -.9999995611961697936D+00,-.9998540628284542641D+00,
     2  -.9993742857790042211D+00,-.9986757113392488361D+00,
     3  -.9981663312653065319D+00,-.9967937950377099684D+00,
     4  -.9915204930349573945D+00,-.9789111460911066598D+00,
     5  -.9546223050957622460D+00,-.9143753170284055811D+00,
     6  -.8549122572812931219D+00,-.7748635328958822693D+00,
     7  -.6756971261110818656D+00,-.5630111646050787013D+00,
     8  -.4456309296883692304D+00,-.3239214235057059439D+00,
     9  -.1910948702891348818D+00,-.6010746179139050748D-01,
     *  0.5535955485419704831D-01,0.1839161423601676057D+00,
     1  0.3218475965375455360D+00,0.4440268722178496513D+00,
     2  0.5388280880869608098D+00,0.6102911884152694592D+00,
     3  0.6890754028268101143D+00,0.7793816513964082101D+00,
     4  0.8590169377375837077D+00,0.9197843418869631276D+00,
     5  0.9608518343800309135D+00,0.9847783013241838600D+00,
     6  0.9959792756330179822D+00,0.9995529842015816084D+00/
c 
        data ws3/
     1  0.2710172363899639507D-04,0.3000102957161655325D-03,
     2  0.6406620985879244039D-03,0.7022231589596028216D-03,
     3  0.4169828005843099860D-03,0.2864509508391933312D-02,
     4  0.8262735415316611861D-02,0.1769708235841680747D-01,
     5  0.3160977274445712196D-01,0.4944128317695028472D-01,
     6  0.6972783990564649429D-01,0.9014775719090364990D-01,
     7  0.1072697514428854242D+00,0.1164095833893789395D+00,
     8  0.1181371349315222418D+00,0.1271961515090219517D+00,
     9  0.1360757481844182560D+00,0.1217083746603252972D+00,
     *  0.1171802384534592941D+00,0.1378365653671021968D+00,
     1  0.1335706651933584059D+00,0.1086492228887410573D+00,
     2  0.8154720454650930354D-01,0.6809169781409895527D-01,
     3  0.8916524755178248021D-01,0.8734614740549273387D-01,
     4  0.7074145306519176363D-01,0.5070661378121436890D-01,
     5  0.3187676354419228349D-01,0.1673084922968641930D-01,
     6  0.6542270137711317794D-02,0.1380356526343036681D-02/
c 
c         In this run, istart=4, and x0=-.9939396001308002D+00
c 
        data xs4/
     1  -.9999518260252942597D+00,-.9995823486165553026D+00,
     2  -.9985763320424136363D+00,-.9969430284589902289D+00,
     3  -.9951705791483398649D+00,-.9940247581377933956D+00,
     4  -.9917760340104341425D+00,-.9841659529287217742D+00,
     5  -.9678800748275633120D+00,-.9389124837395019551D+00,
     6  -.8935166857882854987D+00,-.8291599088807323168D+00,
     7  -.7454856774256087065D+00,-.6458959453700616442D+00,
     8  -.5407239524061488631D+00,-.4392434471124403221D+00,
     9  -.3237756692201044043D+00,-.1896145881486838589D+00,
     *  -.5178064829349148467D-01,0.6260916344637568511D-01,
     1  0.1504774340465042191D+00,0.2797335193027957551D+00,
     2  0.4221471072845578805D+00,0.5585788638165317395D+00,
     3  0.6807741784582883355D+00,0.7837767912127285491D+00,
     4  0.8650718008515473055D+00,0.9243895614470932965D+00,
     5  0.9635129127938489915D+00,0.9859385975184402219D+00,
     6  0.9963114603215704353D+00,0.9995919063964139903D+00/
c 
        data ws4/
     1  0.1482951303979139632D-03,0.6509045807777742764D-03,
     2  0.1362615751364411292D-02,0.1819474862344612060D-02,
     3  0.1614608676602235254D-02,0.8409915876182572087D-03,
     4  0.4458047805944593464D-02,0.1131067114809210399D-01,
     5  0.2194693593030545353D-01,0.3662833012039891275D-01,
     6  0.5460031306497258430D-01,0.7419501390358258192D-01,
     7  0.9264465410253860907D-01,0.1048409580361505453D+00,
     8  0.1029442657772958177D+00,0.1044498115458262945D+00,
     9  0.1268972792652605555D+00,0.1387170878033619828D+00,
     *  0.1333115112510702630D+00,0.8845236595250857296D-01,
     1  0.1081192630758840594D+00,0.1410448582477724772D+00,
     2  0.1411901979309974866D+00,0.1303600307909600946D+00,
     3  0.1132126513124913968D+00,0.9237587673848506219D-01,
     4  0.7017212716684915988D-01,0.4876598697739482392D-01,
     5  0.3007281963667329129D-01,0.1556665807672364581D-01,
     6  0.6023738637942412613D-02,0.1261655111411955074D-02/
c 
c          In this run, istart=5, and x0=-.9846319772006164D+00
c 
        data xs5/
     1  -.9998880964533448140D+00,-.9990795666632642042D+00,
     2  -.9969531600562593481D+00,-.9934603161951082106D+00,
     3  -.9893835286357672655D+00,-.9860810938875468307D+00,
     4  -.9839338964118098077D+00,-.9765942336129069578D+00,
     5  -.9593991704279478900D+00,-.9283801280220172386D+00,
     6  -.8794849888060021140D+00,-.8096060103508239962D+00,
     7  -.7172738148025137250D+00,-.6031354901899315987D+00,
     8  -.4705791339726317251D+00,-.3285280484888934910D+00,
     9  -.2115861332279973471D+00,-.1370527331400241856D+00,
     *  -.5155792656817357172D-01,0.6015951387422131238D-01,
     1  0.2016656059544798384D+00,0.3214426163193529887D+00,
     2  0.4310303982865327018D+00,0.5579545784726604034D+00,
     3  0.6798375048163687963D+00,0.7838325388324340161D+00,
     4  0.8657944979285291481D+00,0.9252595313141716944D+00,
     5  0.9641799888917085084D+00,0.9862934807634978918D+00,
     6  0.9964284531008990516D+00,0.9996069283752986919D+00/
c 
        data ws5/
     1  0.3382171667182602330D-03,0.1394669072058553077D-02,
     2  0.2870806760251017210D-02,0.3979712896235321742D-02,
     3  0.3931555220579354793D-02,0.2355088510771782382D-02,
     4  0.3507847642213335527D-02,0.1167538106440296340D-01,
     5  0.2339848786855782793D-01,0.3933214998834108075D-01,
     6  0.5898736988373545389D-01,0.8102283708879746963D-01,
     7  0.1035416036437526673D+00,0.1241977074149953369D+00,
     8  0.1396251542884990522D+00,0.1399844326465464555D+00,
     9  0.8164607843524904320D-01,0.8752247450279872346D-01,
     *  0.8451891830826709951D-01,0.1374248816596136237D+00,
     1  0.1364038319085116959D+00,0.1050266979578899851D+00,
     2  0.1207571447301366072D+00,0.1278126338316818833D+00,
     3  0.1140652175677714969D+00,0.9329591787436481633D-01,
     4  0.7056643898163543457D-01,0.4869840771602232782D-01,
     5  0.2977675476150388616D-01,0.1527008350888181589D-01,
     6  0.5854658834253543965D-02,0.1216838264961800843D-02/
c 
c                In this run, istart=6, and x0=-.9670976881178759D+00
c 
        data xs6/
     1  -.9998760992422212844D+00,-.9989410489337287758D+00,
     2  -.9962940414728256821D+00,-.9914384954596439871D+00,
     3  -.9847152677732668373D+00,-.9773332841472762147D+00,
     4  -.9709991903915781449D+00,-.9673698908047689770D+00,
     5  -.9622353113790208048D+00,-.9469529001193045423D+00,
     6  -.9184072484597393803D+00,-.8734261524509746682D+00,
     7  -.8095657251210775705D+00,-.7263213460812519476D+00,
     8  -.6277062533226002720D+00,-.5312106897482486285D+00,
     9  -.4426101603743825520D+00,-.3203992949215935030D+00,
     *  -.1744845740848946565D+00,-.2030772673875112680D-01,
     1  0.1344273197057819111D+00,0.2853166444497611602D+00,
     2  0.4294662318506000717D+00,0.5635056080087385477D+00,
     3  0.6831668350729235991D+00,0.7844806117525987740D+00,
     4  0.8649372612665845525D+00,0.9240065470570671772D+00,
     5  0.9631901026843597370D+00,0.9857666842378621140D+00,
     6  0.9962559927863751766D+00,0.9995849410006951906D+00/
c 
        data ws6/
     1  0.3780296562463918451D-03,0.1654504939734537773D-02,
     2  0.3729930378803571242D-02,0.5924848793156371891D-02,
     3  0.7313580192676440472D-02,0.7156069836089001940D-02,
     4  0.5265214966413314260D-02,0.2389460697254466982D-02,
     5  0.9645150724222096055D-02,0.2138311540473365440D-01,
     6  0.3625744210712995653D-01,0.5413022245743846691D-01,
     7  0.7370560247855600694D-01,0.9224907163520458221D-01,
     8  0.1022903867078350800D+00,0.8665294598258052515D-01,
     9  0.1024224381750969821D+00,0.1377490003805539191D+00,
     *  0.1517112627676629497D+00,0.1554092570642524679D+00,
     1  0.1533504504690376497D+00,0.1479873615500163231D+00,
     2  0.1397563232217295763D+00,0.1275856422967976874D+00,
     3  0.1110550216517057358D+00,0.9114632820414514031D-01,
     4  0.6967810018278739163D-01,0.4871500778577300073D-01,
     5  0.3020741730608970128D-01,0.1571300547584723895D-01,
     6  0.6105243411650397968D-02,0.1282563098779397411D-02/
c 
c        In this run, istart=7, and x0=-.9378752307252274D+00
c 
        data xs7/
     1  -.9998358374501824680D+00,-.9985899918401509000D+00,
     2  -.9950125217534864597D+00,-.9882499726698025057D+00,
     3  -.9783418463810158910D+00,-.9663335502111788026D+00,
     4  -.9541089646326543389D+00,-.9440125659525870226D+00,
     5  -.9383086294710821140D+00,-.9309022727655344977D+00,
     6  -.9095950036761577989D+00,-.8713660851829913050D+00,
     7  -.8138194555449531548D+00,-.7362296842228479750D+00,
     8  -.6425782892237149884D+00,-.5472787176748581002D+00,
     9  -.4509539791684268448D+00,-.3259658533355079316D+00,
     *  -.1831972288018991167D+00,-.5712743774329706060D-01,
     1  0.5182274372195018960D-01,0.1978763721629104845D+00,
     2  0.3559267091135993880D+00,0.5067396737215971617D+00,
     3  0.6418168043464382988D+00,0.7560939800867180553D+00,
     4  0.8468497544684947066D+00,0.9135952064522908126D+00,
     5  0.9580109473013622088D+00,0.9837077887536866508D+00,
     6  0.9957003934452009366D+00,0.9995221207651164729D+00/
c 
        data ws7/
     1  0.5013528988961932450D-03,0.2213327008852682630D-02,
     2  0.5092339287413238516D-02,0.8427685072789274154D-02,
     3  0.1121132458372952888D-01,0.1248526310064434057D-01,
     4  0.1156360848260765777D-01,0.8288280235163019757D-02,
     5  0.3672544509269183985D-02,0.1369814252053927574D-01,
     6  0.2932105622405707661D-01,0.4756373057822552046D-01,
     7  0.6769113338018584375D-01,0.8694961685967786958D-01,
     8  0.9784156142588736260D-01,0.9117357431340720049D-01,
     9  0.1088371584107809928D+00,0.1381041447562130190D+00,
     *  0.1419490618039433456D+00,0.1060600562619659756D+00,
     1  0.1274250896130018414D+00,0.1570467847162666523D+00,
     2  0.1562859794454257146D+00,0.1440029300479903604D+00,
     3  0.1253177362676318078D+00,0.1027835130069626773D+00,
     4  0.7864719490545090904D-01,0.5512755430770910046D-01,
     5  0.3431430355456293688D-01,0.1792992302263335140D-01,
     6  0.6998306722433790942D-02,0.1475722675682320812D-02/
c 
c            In this run, istart=8, and x0=-.8936088197087079D+00
c 
        data xs8/
     1  -.9998349077262792498D+00,-.9985374242999697415D+00,
     2  -.9945984085317443193D+00,-.9865995991985567417D+00,
     3  -.9738547128894023712D+00,-.9568477220981969793D+00,
     4  -.9373612140109565110D+00,-.9182517932719111363D+00,
     5  -.9028746008549063513D+00,-.8942912853172609206D+00,
     6  -.8853635307104528846D+00,-.8607876786593975995D+00,
     7  -.8177475150101125257D+00,-.7542325337856407087D+00,
     8  -.6689109432888348420D+00,-.5617831991473933825D+00,
     9  -.4345169509249065440D+00,-.2905381956742285855D+00,
     *  -.1350056815147131610D+00,0.2505797705980400142D-01,
     1  0.1797152728241713825D+00,0.3147142574734737765D+00,
     2  0.4302452627611430860D+00,0.5503766170794395689D+00,
     3  0.6694638525882006225D+00,0.7739363166595869601D+00,
     4  0.8578925066195015936D+00,0.9198556877139015681D+00,
     5  0.9610980141271851387D+00,0.9849283159287363594D+00,
     6  0.9960286142921206277D+00,0.9995591552301776364D+00/
c 
        data ws8/
     1  0.5080326120883565890D-03,0.2350966739055602704D-02,
     2  0.5774524267773412400D-02,0.1034503739138198888D-01,
     3  0.1506630243141160232D-01,0.1864683355546726948D-01,
     4  0.1984040118366424373D-01,0.1780474003769358927D-01,
     5  0.1246037205795116761D-01,0.5241897946985599047D-02,
     6  0.1597924111721508091D-01,0.3346068310147962043D-01,
     7  0.5297152641889654173D-01,0.7429123454351759350D-01,
     8  0.9635726711153008767D-01,0.1176215221397241868D+00,
     9  0.1363390221814507287D+00,0.1507593385880129049D+00,
     *  0.1591486782799792178D+00,0.1593758814751191158D+00,
     1  0.1473922752598839616D+00,0.1216893623846639988D+00,
     2  0.1155476970961378191D+00,0.1225566304448899936D+00,
     3  0.1132931866479133643D+00,0.9474139697961808553D-01,
     4  0.7294725309108180866D-01,0.5120127015697637557D-01,
     5  0.3184535627487138120D-01,0.1661043338995972376D-01,
     6  0.6469839435283804348D-02,0.1361795658321789067D-02/
c 
c               In this run, istart=9, and x0=-.8315271563892494D+00
c 
        data xs9/
     1  -.9997709944360571917D+00,-.9979723251500424095D+00,
     2  -.9925114352350112039D+00,-.9813937969716089167D+00,
     3  -.9635464116013281359D+00,-.9393609320574562594D+00,
     4  -.9108778212740556811D+00,-.8816006597205271234D+00,
     5  -.8559436180573416718D+00,-.8383569964319117720D+00,
     6  -.8298006474200089997D+00,-.8110521432165379908D+00,
     7  -.7726784405720947096D+00,-.7132845601128985343D+00,
     8  -.6318798534155719492D+00,-.5285346370026617761D+00,
     9  -.4048020841406228520D+00,-.2638865439772917354D+00,
     *  -.1106997691835749171D+00,0.4755343614627565541D-01,
     1  0.1941156568541422588D+00,0.2836976791449238348D+00,
     2  0.3954688026193675932D+00,0.5337135864587481394D+00,
     3  0.6610626539619686980D+00,0.7692462951205338498D+00,
     4  0.8551947342241217028D+00,0.9183676504699466922D+00,
     5  0.9603658369883745869D+00,0.9846360345271582641D+00,
     6  0.9959489883059117249D+00,0.9995500707993186499D+00/
c 
        data ws9/
     1  0.7045847299554421220D-03,0.3258454604124073978D-02,
     2  0.8009918594877592576D-02,0.1441340780738612672D-01,
     3  0.2121915455276333527D-01,0.2680863070303332818D-01,
     4  0.2956441585214224955D-01,0.2823911021185986675D-01,
     5  0.2231433163471466445D-01,0.1212926786833278002D-01,
     6  0.9753061949782314186D-02,0.2833196945768399991D-01,
     7  0.4866063145130080599D-01,0.7030382123848590487D-01,
     8  0.9250088413133834395D-01,0.1139332083011496836D+00,
     9  0.1329978898681784916D+00,0.1480187285081904580D+00,
     *  0.1572058157643139845D+00,0.1571640423388585067D+00,
     1  0.1248748098770059956D+00,0.7801939708098705628D-01,
     2  0.1367327061021876033D+00,0.1350110180080579368D+00,
     3  0.1185309943101749476D+00,0.9734632417351025743D-01,
     4  0.7446464109090195717D-01,0.5215238765214251654D-01,
     5  0.3242690509805186939D-01,0.1692279154491102336D-01,
     6  0.6597033193846250310D-02,0.1389662299750595687D-02/
c 
c            In this run, istart=10, and x0=-.7498587278001463D+00
c 
        data xs10/
     1  -.9997482743076812306D+00,-.9977582868694432504D+00,
     2  -.9916479053662445706D+00,-.9789976578103545065D+00,
     3  -.9582015409989833513D+00,-.9290584156284426831D+00,
     4  -.8930401959056736776D+00,-.8532564030558507885D+00,
     5  -.8141431812654303172D+00,-.7808702948499468547D+00,
     6  -.7584708978367117192D+00,-.7477382877892638469D+00,
     7  -.7248150687276532122D+00,-.6789234792122278524D+00,
     8  -.6103567206281988420D+00,-.5211908973284752815D+00,
     9  -.4181722548705072484D+00,-.3159808001443904916D+00,
     *  -.2103964907430318886D+00,-.7815650494905596968D-01,
     1  0.7278861534368519747D-01,0.2298572932879338945D+00,
     2  0.3842071554305917266D+00,0.5287392376502982437D+00,
     3  0.6576692159927157031D+00,0.7667190565423534589D+00,
     4  0.8534019924051391706D+00,0.9172286180450999555D+00,
     5  0.9597507574771623688D+00,0.9843752199814909972D+00,
     6  0.9958751667955766379D+00,0.9995414578276628310D+00/
c 
        data ws10/
     1  0.7755406824999389138D-03,0.3618691552041047565D-02,
     2  0.9019931029505426535D-02,0.1655822279090996422D-01,
     3  0.2507192152276275317D-01,0.3296833328083282892D-01,
     4  0.3854469430412849815D-01,0.4027372642756439747D-01,
     5  0.3707614818367742131D-01,0.2861186221566506452D-01,
     6  0.1535905225984147389D-01,0.1203489428935130669D-01,
     7  0.3438003567182247076D-01,0.5737260954343704894D-01,
     8  0.7945554100684543523D-01,0.9786631453184739434D-01,
     9  0.1053868374458596068D+00,0.9897437044783334936D-01,
     *  0.1179475620363505306D+00,0.1441246451674408304D+00,
     1  0.1556956187631841835D+00,0.1570037236791187483D+00,
     2  0.1505166362061773226D+00,0.1375793130988734233D+00,
     3  0.1195625558160842306D+00,0.9811905248964000785D-01,
     4  0.7516436649476952446D-01,0.5275044920318020751D-01,
     5  0.3286960603685586976D-01,0.1718902154259226859D-01,
     6  0.6712773354646412125D-02,0.1415948924660821692D-02/
c 
c              In this run, istart=11, and x0=-.6481288501291113D+00
c 
        data xs11/
     1  -.9997346750359577122D+00,-.9976242444937373907D+00,
     2  -.9910747182132290161D+00,-.9773045619568647451D+00,
     3  -.9541911086427072488D+00,-.9209025230135398428D+00,
     4  -.8782511031687337826D+00,-.8287771578506852229D+00,
     5  -.7765912278261459172D+00,-.7269834877627494453D+00,
     6  -.6857812283415239171D+00,-.6584643687523677150D+00,
     7  -.6459358354666513678D+00,-.6216657485281085780D+00,
     8  -.5732647221513468578D+00,-.5010596447007695831D+00,
     9  -.4061682583253128928D+00,-.2908122108070092433D+00,
     *  -.1584632314224486765D+00,-.1377312499009035958D-01,
     1  0.1376833272414890577D+00,0.2898007564904653967D+00,
     2  0.4364158549993308801D+00,0.5718077418397891629D+00,
     3  0.6911774717957827063D+00,0.7910625959537257057D+00,
     4  0.8696471201052809523D+00,0.9269303255823250859D+00,
     5  0.9647177500910325653D+00,0.9863932723425773001D+00,
     6  0.9964283181495655472D+00,0.9996046289615737365D+00/
c 
        data ws11/
     1  0.8184900702944149204D-03,0.3851429364295016711D-02,
     2  0.9725697388817119945D-02,0.1817906750680321631D-01,
     3  0.2819231829367376590D-01,0.3825078744204345231D-01,
     4  0.4662601182383247989D-01,0.5163033642029724517D-01,
     5  0.5184918370848928020D-01,0.4637984694760218346D-01,
     6  0.3509410204666828855D-01,0.1868423833322444901D-01,
     7  0.1289317272918894971D-01,0.3634508453608744556D-01,
     8  0.6041152311967941467D-01,0.8381805078312024543D-01,
     9  0.1055969545156022034D+00,0.1245352637074865459D+00,
     *  0.1393829797075928822D+00,0.1490614304333705965D+00,
     1  0.1528259095978416523D+00,0.1503732756598519124D+00,
     2  0.1418953443631919090D+00,0.1280816018118827587D+00,
     3  0.1100755316493385729D+00,0.8938915194359487648D-01,
     4  0.6777982976550617288D-01,0.4709119067209038755D-01,
     5  0.2905573616036349808D-01,0.1505360340654082470D-01,
     6  0.5830693311016585454D-02,0.1222162780611257001D-02/
c 
c               In this run, istart=12, and x0=-.5273106402545079D+00
c 
        data xs12/
     1  -.9997309656162974340D+00,-.9975760769450798703D+00,
     2  -.9908130745808852956D+00,-.9763802077804145839D+00,
     3  -.9517019496972759936D+00,-.9153492395704071427D+00,
     4  -.8674537263775063464D+00,-.8098707401156408642D+00,
     5  -.7461224289362385385D+00,-.6811536855931260783D+00,
     6  -.6209109317084894416D+00,-.5717186842218309595D+00,
     7  -.5394387983524397844D+00,-.5248404368888288361D+00,
     8  -.4974301389929197654D+00,-.4433681240540339261D+00,
     9  -.3640503441344544618D+00,-.2619226048079057248D+00,
     *  -.1405776730348100792D+00,-.4725397498031004270D-02,
     1  0.1400063585237026327D+00,0.2874157711474132881D+00,
     2  0.4311697238922658344D+00,0.5652936927797847527D+00,
     3  0.6846587463055912060D+00,0.7854191277850609933D+00,
     4  0.8653578182280043282D+00,0.9240999041379404545D+00,
     5  0.9631519560420232682D+00,0.9857181146242100345D+00,
     6  0.9962349818713379790D+00,0.9995819047631831275D+00/
c 
        data ws12/
     1  0.8311764895867491928D-03,0.3947916558380621223D-02,
     2  0.1010250457972375691D-01,0.1920368144930325970D-01,
     3  0.3039790794455134968D-01,0.4228819834350232140D-01,
     4  0.5319257204121356454D-01,0.6138227733271408312D-01,
     5  0.6528324050475681157D-01,0.6365291181512850766D-01,
     6  0.5576410030457813562D-01,0.4162688206912885737D-01,
     7  0.2201507175239966684D-01,0.1466165804951929547D-01,
     8  0.4089229436537407140D-01,0.6700573848546130567D-01,
     9  0.9122418124794578868D-01,0.1124320729139425083D+00,
     *  0.1294694244725066836D+00,0.1412912991505435381D+00,
     1  0.1471310732813205867D+00,0.1466239477411066868D+00,
     2  0.1398799761375416747D+00,0.1275057457731448440D+00,
     3  0.1105775483902208537D+00,0.9057005846248946160D-01,
     4  0.6924406491882868258D-01,0.4849412278666582391D-01,
     5  0.3015162386431848376D-01,0.1573334404214577812D-01,
     6  0.6131988064120361787D-02,0.1291396667829524816D-02/
c 
c              In this run, istart=13, and x0=-.3898227981812752D+00
c 
        data xs13/
     1  -.9996743663333813201D+00,-.9970760609267320046D+00,
     2  -.9889628252015155401D+00,-.9717322878884907004D+00,
     3  -.9423678706935715530D+00,-.8991440422776586289D+00,
     4  -.8420399650119803740D+00,-.7728884496106738188D+00,
     5  -.6953118728161332065D+00,-.6144929039452541939D+00,
     6  -.5368111946885448058D+00,-.4693426619341480250D+00,
     7  -.4191478009870577313D+00,-.3922073555456071826D+00,
     8  -.3750480918593072856D+00,-.3338794062646430111D+00,
     9  -.2668810963256465319D+00,-.1767264736272751142D+00,
     *  -.6711768073937825164D-01,0.5726092305896516670D-01,
     1  0.1909364933529308899D+00,0.3279498553428699521D+00,
     2  0.4622402392398590684D+00,0.5880782527520530318D+00,
     3  0.7005081477540258884D+00,0.7957585892909187933D+00,
     4  0.8715805268565947785D+00,0.9274697700952386641D+00,
     5  0.9647281851656198115D+00,0.9863093538208400732D+00,
     6  0.9963867517396090201D+00,0.9995984456126154229D+00/
c 
        data ws13/
     1  0.1005159476506734785D-02,0.4750964930760212377D-02,
     2  0.1209203461956900733D-01,0.2288366352946895188D-01,
     3  0.3613956406674540002D-01,0.5031391580066095095D-01,
     4  0.6358758472113901439D-01,0.7410781518379933264D-01,
     5  0.8017549888256400588D-01,0.8039097125108087335D-01,
     6  0.7378171156621664953D-01,0.5996084331074560756D-01,
     7  0.3942028401535790220D-01,0.1498125002099087119D-01,
     8  0.2746458333330075843D-01,0.5444865681429772924D-01,
     9  0.7910601649775081373D-01,0.1005850488542713753D+00,
     *  0.1178485438468243786D+00,0.1299914247566924951D+00,
     1  0.1363600241394352731D+00,0.1366492615774449579D+00,
     2  0.1309691387642801574D+00,0.1198739856019537879D+00,
     3  0.1043515174923653201D+00,0.8577070659868863765D-01,
     4  0.6578798910593158759D-01,0.4620989302458475970D-01,
     5  0.2880625955024519854D-01,0.1506429670180724326D-01,
     6  0.5881317412461502118D-02,0.1240074552928983826D-02/
c 
c          In this run, istart=14, and x0=-.2393822643556076D+00
c 
        data xs14/
     1  -.9996483446197229611D+00,-.9968419849047731070D+00,
     2  -.9880746708852355886D+00,-.9694330773867505890D+00,
     3  -.9375927730530397253D+00,-.8905471732370232284D+00,
     4  -.8280210293277734807D+00,-.7516011922187256671D+00,
     5  -.6646410046012466459D+00,-.5719915029363444641D+00,
     6  -.4796050422955632274D+00,-.3940516011970269292D+00,
     7  -.3219845336256324457D+00,-.2695711284511764571D+00,
     8  -.2418185249083301606D+00,-.2240177578816958987D+00,
     9  -.1813524643620559640D+00,-.1124912724805158392D+00,
     *  -.2092277933298562490D-01,0.8872996497274819526D-01,
     1  0.2108956480379620807D+00,0.3393743897385141768D+00,
     2  0.4677465985899377407D+00,0.5898398901876685144D+00,
     3  0.7002125809969499249D+00,0.7946075488037435573D+00,
     4  0.8703298098762203560D+00,0.9265016248726679838D+00,
     5  0.9641466717717659900D+00,0.9860472543136031925D+00,
     6  0.9963099001501866105D+00,0.9995893047194283139D+00/
c 
        data ws14/
     1  0.1085503482984686437D-02,0.5132035409612069450D-02,
     2  0.1307179611226592979D-01,0.2477770888070323826D-01,
     3  0.3924411676218627585D-01,0.5489499416356494520D-01,
     4  0.6989109783255577184D-01,0.8238509852297265142D-01,
     5  0.9072006985112745781D-01,0.9357872327121396301D-01,
     6  0.9009022642887468066D-01,0.7989843728238347756D-01,
     7  0.6319685942021370750D-01,0.4079409957642599941D-01,
     8  0.1539358971506062435D-01,0.2856110024497893513D-01,
     9  0.5624840788297131721D-01,0.8089386210189439436D-01,
     *  0.1014690970591456904D+00,0.1169024652845943205D+00,
     1  0.1263886986207095094D+00,0.1294909692785292466D+00,
     2  0.1262165476940895023D+00,0.1170554124526219711D+00,
     3  0.1029740469784733222D+00,0.8536235999582521460D-01,
     4  0.6593433774942311156D-01,0.4658212471564336884D-01,
     5  0.2917877754031427311D-01,0.1531978539671338575D-01,
     6  0.5999796126844213738D-02,0.1267854199762016135D-02/
c 
c             In this run, istart=15, and x0=-.8073141183800038D-01
c 
        data xs15/
     1  -.9996253782114384633D+00,-.9966343372888683899D+00,
     2  -.9872798389151078570D+00,-.9673488445610259456D+00,
     3  -.9331915713000219915D+00,-.8824678212239421550D+00,
     4  -.8145694888888045608D+00,-.7307591474244382239D+00,
     5  -.6340856604357428486D+00,-.5291340337214084631D+00,
     6  -.4216613720594536807D+00,-.3181702574683045210D+00,
     7  -.2254813727587647026D+00,-.1504225150141334079D+00,
     8  -.9996900441726395022D-01,-.7852615133575122404D-01,
     9  -.5142326992568013710D-01,0.1258403553683317988D-02,
     *  0.7659410435466550009D-01,0.1704444141541341323D+00,
     1  0.2776484354868234645D+00,0.3923560388525985767D+00,
     2  0.5084265220252387664D+00,0.6198824027166801324D+00,
     3  0.7213849124779990759D+00,0.8086874716071155050D+00,
     4  0.8790195378339900850D+00,0.9313520383770851797D+00,
     5  0.9664932938894065221D+00,0.9869585240275307517D+00,
     6  0.9965514236652007603D+00,0.9996162877908351640D+00/
c 
        data ws15/
     1  0.1156492852085173009D-02,0.5471501064779112251D-02,
     2  0.1395719924832306224D-01,0.2652560110821187983D-01,
     3  0.4218589715715332352D-01,0.5936511098429136823D-01,
     4  0.7621958272388372982D-01,0.9088866562753590146D-01,
     5  0.1016879541017937416D+00,0.1072485362542538069D+00,
     6  0.1066056749905948536D+00,0.9923394632175916017D-01,
     7  0.8500622637910862804D-01,0.6396923867038026820D-01,
     8  0.3552501243799077583D-01,0.1547732383589198180D-01,
     9  0.4033339156551725579D-01,0.6459817843890424655D-01,
     *  0.8537808230760465305D-01,0.1014569010298512688D+00,
     1  0.1119694406373986264D+00,0.1164146687135715873D+00,
     2  0.1147220190960762719D+00,0.1072922751639992201D+00,
     3  0.9499846792567827273D-01,0.7914034957776171434D-01,
     4  0.6135053734364384961D-01,0.4345120490382326552D-01,
     5  0.2725735579112166385D-01,0.1431958918290312659D-01,
     6  0.5607935521669248020D-02,0.1184651405606506351D-02/
c 
c        copy the arrays where tha data are stored to the output arrays
c 
        do 1200 i=1,32
        xs(i,1)=xs1(i)
        xs(i,2)=xs2(i)
        xs(i,3)=xs3(i)
        xs(i,4)=xs4(i)
        xs(i,5)=xs5(i)
        xs(i,6)=xs6(i)
        xs(i,7)=xs7(i)
        xs(i,8)=xs8(i)
        xs(i,9)=xs9(i)
        xs(i,10)=xs10(i)
        xs(i,11)=xs11(i)
        xs(i,12)=xs12(i)
        xs(i,13)=xs13(i)
        xs(i,14)=xs14(i)
        xs(i,15)=xs15(i)
c 
        ws(i,1)=ws1(i)
        ws(i,2)=ws2(i)
        ws(i,3)=ws3(i)
        ws(i,4)=ws4(i)
        ws(i,5)=ws5(i)
        ws(i,6)=ws6(i)
        ws(i,7)=ws7(i)
        ws(i,8)=ws8(i)
        ws(i,9)=ws9(i)
        ws(i,10)=ws10(i)
        ws(i,11)=ws11(i)
        ws(i,12)=ws12(i)
        ws(i,13)=ws13(i)
        ws(i,14)=ws14(i)
        ws(i,15)=ws15(i)
 1200 continue
c 
        do 1600 i=1,15
        do 1400 j=1,32
c 
        ws(32-j+1,30-i+1)=ws(j,i)
        xs(32-j+1,30-i+1)=-xs(j,i)
 1400 continue
 1600 continue
c 
        return
        end
