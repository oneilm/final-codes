        implicit real *8 (a-h,o-z)
        dimension xs(180),whts(180),numspts(100)
        dimension points(34,25),weights(34,25),
     1      coefslrg(20,34,25),
     2      fs(100),fslrg(34,25),fs2(1000),diffs(100)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER n'
        READ *,n
        CALL PRINf('n=*',n,1 )
c 
ccc        iw=21
cccc        call o1rwrt2(iw)
  
cccc        stop
  
cccc        call o1rmatb2(coefslrg,xs,whts,points,weights,numspts)
  
        ir=21
        call o1rrd2(ir,xs,whts,points,weights,
     1      numspts,coefslrg)
  
  
        call prinf('after o1rmatb2, numspts=*',numspts,25)
  
c 
c        construct the test function to be interpolated
c 
        nn=20
        do 1200 i=1,nn
c 
        fs(i)=sigmafun(xs(i))
 1200 continue
c 
        call prin2('test function as created*',fs,nn)
c 
c       interpolate the test function to all of the nodes,
c       and compate the results with the function evaluated directly
c 
        call o1rappb2(coefslrg,numspts,fs,fslrg)
  
c 
        do 2000 k=1,25
c 
        numpts=numspts(k)
        do 1800 i=1,numpts
c 
        fs2(i)=sigmafun(points(i,k))
c 
        diffs(i)=fs2(i)-fslrg(i,k)
  
 1800 continue
c 
        call prinf('k=*',k,1)
        call prin2('and points(1,k)=*',points(1,k),numpts)
        call prin2('fslrg(1,k)=*',fslrg(1,k),numpts)
        call prin2('and fs2 directly is*',fs2,numpts)
        call prin2('and diffs=*',diffs,numpts)
c 
 2000 continue
c 
        stop
        end
  
c 
c 
c 
c 
        function sigmafun(x)
        implicit real *8 (a-h,o-z)
        save
  
        sigmafun=(x+1)*log(x+1)
  
cccc        sigmafun=log(x+1)
  
        sigmafun=sigmafun*sin(4*x)
c 
  
cccc        sigmafun=1
  
  
        return
        end
c 
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This file contains maicinery for the construction of "near"
c        quadratures. It contains 5 user-callable subroutines:
c        o1rnear,o1rrd2,o1rwrt2,o1rmatb2,o1rappb2. Following is
c        a description of these subroutines:
c 
c 
c   o1rmatb2 - returns to the user a set of 25 matrices.
c       Each of these matrices converts the values of a function at
c       the support nodes (20 of them things) into the values of that
c       function at the quadrature nodes corresponding to one zone.
c       Since there are 25 zones, the number of matrices returned by
c       this subroutine is also 25.
c 
c   o1rappb2 - applies to the user-specified vector of length 20
c       (representing the values at the 20 support nodes of an
c       appropriately-behaved function) the 20 interpolation matrices
c       produced by the subroutine o1rmatbg (see). The result is a
c       collection of (approximate) values  of the said function at
c       the nodes points (also produced by the subroutine o1rmatbg).
c       This subroutine is principally a tool for the debugging of the
c       subroutine o1rmatbg, since I am not aware of any other uses for
c       this subroutine.
c 
c   o1rnear - returns to the user the nodes and weights of a quadrature
c        for the "near" integration for integral equations of
c        two-dimensional scattering theory. It assumes that the
c        function is to be integrated on the interval [-1,1] in R^2,
c        and that the points (x,y) at which the integral is to be
c        evaluated is close to the end (1,0) of the interval. The
c        quadratures are valid outside a the cone with the vertex at
c        the point (1,1) and with the angre 5^o (five degrees)
c 
c   o1rwrt2 -  constructs (by calling the subroutine o1rmatbg) and
c        stores on disk on the FORTRAN unit number iw all of the
c        parameters used for the construction of the "near" quadratures.
c        More specifically, it constructs the support nodes and
c        corresponding weights, auxiliary quadrature nodes for each of
c        the 25 zones, interpolation coefficients from the support nodes
c        to the auxiliary quadrature nodes for each of the zones, and
c        the integer array numspts whose i-th entry is the number of
c        auxiliary quadrature nodes for the zone number i
c 
c   o1rrd2 - reads from the FORTRAN unit ir the various types of data
c        to be used in the construction of the "near" quadratures. The
c        data to be read have been (hopefully) stored previously by a
c        call to the subroutine o1rwrt2 (see)
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        subroutine o1rappb2(coefslrg,numspts,fs,fslrg)
        implicit real *8 (a-h,o-z)
        save
        dimension coefslrg(20,34,25),numspts(1),fs(1),
     1      fslrg(34,25)
c 
c       This subroutine applies to the user-specified vector fs of
c       length 20 (representing the values at the 20 support nodes
c       of an appropriately-behaved function f) the 20 interpolation
c       matrices produced by the subroutine o1rmatbg (see). The result
c       is the collection fslrs of (approximate) values  of the function
c       f at the nodes points (also produced by the subroutine o1rmatbg).
c       This subroutine is principally a tool for the debugging of the
c       subroutine o1rmatbg, since I am not aware of any other uses for
c       this subroutine.
c 
c                          Input parameters:
c 
c  coefslrg - a 20 \times 24 \times 20 - array, containing the set of 20
c       matrices described above.
c  numspts - integer array of length 20, containing numbers of quadrature
c       nodes corresponding to the 20 support nodes. Specifically,
c       numspts(i) is the number of quadrature nodes corresponding to
c       the support node number i
c  fs - the values of the function f to be interpolated at the 20 support
c       nodes on the interval [-1,1]
c 
c                          Output parameters:
c 
c  fslrg - the interpolated values of the function f at the quadrature
c       nodes points (the latter are returned by o1rmatbg)
c 
c        . . . given a user-supplied table of function values at the
c              support nodes, evaluate that function at the nodes of all
c              quadratures for all support nodes
c 
        nn=25
        do 2000 k=1,nn
c 
        numpts=numspts(k)
c 
        call prinf('in o1rappbg, k=*',k,1)
        call prinf('and numpts=*',numpts,1)
c 
        do 1800 i=1,numpts
        d=0
        do 1600 j=1,nn
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
c 
c 
c 
        subroutine o1rmatb2(coefslrg,xs,whts,points,weights,numspts)
        implicit real *8 (a-h,o-z)
        save
        dimension points(34,25),weights(34,25),xs(100),
     1      whts(100),dumpts(100),dumwhts(100),coefslrg(20,34,25),
     2      numspts(1),forminte(100)
c 
        dimension thresh(2,26)
c 
c       This subroutine returns to the user a set of 25 matrices.
c       Each of these matrices converts the values of a function at
c       the support nodes (20 of them things) into the values of that
c       function at the quadrature nodes corresponding to one zone.
c       Since there are 25 zones, the number of matrices
c       returned by this subroutine is also 25. On the other hand, there
c       are two anomalies associated with the matrices returned by this
c       subroutine:
c 
c   1. The matrices returned by this subroutine are actually adjoints
c       of the matrices converting the values of a function at the
c       support nodes into its values at the quadrature nodes
c   2. While all of the matrices are declared as 34 \times 25 - matrices,
c       the actual lengths of their rows vary from 22 to 33, and are
c       given by the the elements of the integer array numspts.
c 
c       More specifically, the value f_{k,i} of the function f at the
c       i-th quadrature node corresponding to the k-th support node is
c       given by the formula
c 
c            f_{k,i}=\sum_{j=1}^{20} coefslrg(j,i,k)*f(j),             (1)
c 
c       with k=1,2,...,20, and i=1,2,...,numspts(k).
c 
c       In addition, the subroutine returns a number of other objects
c       that might be useful in the environment for which it has been
c       designed: support nodes xs and corresponding weights whts, the
c       quadrature nodes and weights corresponding to each of the
c       support nodes (all of these are returned in the arrays points,
c       weights, respectively), and the integer array numspts containing
c       numbers of quadrature nodes correponding to all 25 zones.
c 
c 
c                 Input parameters: None
c 
c                 Output parameters:
c 
c  coefslrg - a 20 \times 34 \times 25 - array, containing the set of 25
c       matrices described above.
c  xs - the 20 support nodes on the interval [-1,1]
c  whts - quadrature weights (pretty universal) corresponding to the
c       nodes xs.
c  points - the array containing the quadrature nodes corresponding to
c       all support nodes; the quadrature nodes corresponding to the
c       i-th support nodes are returned in the i-th column of the array
c       points. Please note that the numbers of quadrature nodes are
c       different for different support nodes, so that a part of the
c       array points is returned empty; the number of quadrature nodes
c       actually stored in the i-th column of the array points is
c       returned in the i-th location of the array numspts.
c  weights - the array containing the quadrature weights corresponding to
c       the nodes return in the array points. One weight correponds to
c       each of the nodes in the array points, so that the structures of
c       these two arrays are identical.
c  numspts - integer array of length 25, containing numbers of quadrature
c       nodes corresponding to the 25 zones. Specifically, numspts(i) is
c       the number of quadrature nodes corresponding to the zone number i
c 
        data thresh/
     1  0.10000E-05,0.20000E-05,0.20000E-05,0.40000E-05,
     2  0.40000E-05,0.80000E-05,0.80000E-05,0.16000E-04,
     3  0.16000E-04,0.32000E-04,0.32000E-04,0.64000E-04,
     4  0.64000E-04,0.12800E-03,0.12800E-03,0.25600E-03,
     5  0.25600E-03,0.51200E-03,0.51200E-03,0.10240E-02,
     6  0.10240E-02,0.20480E-02,0.20480E-02,0.40960E-02,
     7  0.40960E-02,0.81920E-02,0.81920E-02,0.16384E-01,
     8  0.16384E-01,0.32768E-01,0.32768E-01,0.49152E-01,
     9  0.49152E-01,0.73728E-01,0.73728E-01,0.11059E+00,
     *  0.11059E+00,0.16589E+00,0.16589E+00,0.24883E+00,
     *  0.24883E+00,0.37325E+00,0.37325E+00,0.55987E+00,
     *  0.55987E+00,0.83981E+00,0.83981E+00,0.12597E+01,
     *  0.12597E+01,0.18896E+01,0.18896E+01,0.00000E+00/
c 
c        . . . retrieve from the subroutine o1rnodes the support nodes
c              and corresponding weights
c 
        call o1rnodes(xs,whts,nn)
c 
c        retrieve from the subroutine o1rnear the nodes and weights of
c        all quadratures
c 
        do 1400 k=1,25
c 
        call prinf('k=*',k,1)
c 
        x=1
        y=(thresh(1,k)+thresh(2,k))/2
c 
        itype=2
        call o1rnear(ier,x,y,itype,dumpts,dumwhts,n)
  
        call prinf('after nearnoe, n=*',n,1)
        call prin2('dumpts=*',dumpts,n)
        call prin2('dumwhts=*',dumwhts,n)
  
        do 1200 i=1,n
c 
        points(i,k)=dumpts(i)
        weights(i,k)=dumwhts(i)
 1200 continue
c 
        numspts(k)=n
 1400 continue
  
  
        call prinf('and numspts=*',numspts,25)
c 
c       for each of the 20 sets of quadrature nodes, construct the matrix
c       of interpolation coefficients connecting the values of a function
c       at the support nodes to its values at the quadrature nodes
c       PLEASE NOTE THAT DIFFERENT quadrature formulae have different
c       sets of nodes
c 
        nnn=25
        do 2000 k=1,nnn
        call prinf('in o1rmatb2, k=*',k,1)
  
c 
        numpts=numspts(k)
  
cccc        call prinf('and numpts=*',numpts,1)
c 
        do 1800 i=1,numpts
  
cccc        call prinf('and i=*',i,1)
c 
        x=points(i,k)
  
        call o1rinter(x,forminte,dumpts,dumwhts,nn)
  
        do 1600 j=1,nn
        coefslrg(j,i,k)=forminte(j)
 1600 continue
  
cccc        call prin2('and forminte =*',forminte,nn)
  
 1800 continue
 2000 continue
c 
  
cccc        stop
        return
        end
  
  
c 
c 
c 
c 
c 
        subroutine o1rnear(ier,x,y,itype,xs,whts,n)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(180),whts(180)
c 
c 
        dimension ns1(25),ns2(25),ns3(25)
        dimension rsmalls(26),thresh(2,30)
c 
c        this subroutine returns to the user the nodes and weights of
c        a quadrature for the "near" integration for integral equations
c        of two-dimensional scattering theory. It assumes that the
c        function is to be integrated on the interval [-1,1] in R^2,
c        and that the points (x,y) at which the integral is to be
c        evaluated is close to the end (1,0) of the interval. The
c        quadratures are valid outside a the cone with the vertex at
c        the point (1,1) and with the angre 5^o (five degrees)
c 
c 
c                  Input parameters:
c 
c  (x,y) - the coordinates of the point where the integral is to be
c        evaluates (assumed to be near the end (1,0).
c  itype - the accuracy requested:
c     itype=1 means four-digit accuracy
c     itype=2 means seven-digit accuracy
c     itype=3 means nine-digit accuracy (much better most of the time)
c 
c 
c                  Output parameters:
c 
c  xs - the nodes of the quadrature
c  xs - the weights corresponding to the nodes xs
c  n - the number of nodes (and weights) returned
c 
         data ifcalled/0/
c 
        data rsmalls/
     1  0.10000000E-05,0.20000000E-05,
     2  0.40000000E-05,0.80000000E-05,
     3  0.16000000E-04,0.32000000E-04,
     4  0.64000000E-04,0.12800000E-03,
     5  0.25600000E-03,0.51200000E-03,
     6  0.10240000E-02,0.20480000E-02,
     7  0.40960000E-02,0.81920000E-02,
     8  0.16384000E-01,0.32768000E-01,
     9  0.49152000E-01,0.73728000E-01,
     *  0.11059200E+00,0.16588800E+00,
     1  0.24883200E+00,0.37324800E+00,
     2  0.55987200E+00,0.83980800E+00,
     3  0.12597120E+01,0.18895680E+01/
c 
        data ns1/13,14,15,15,15,15,17,16,17,17,17,
     1     18,18,19,20,18,18,18,19,19,19,19,20,20,18/
c 
        data ns2/23,23,23,25,25,25,27,28,30,30,30,
     1     31,31,32,33,30,30,30,30,30,30,31,31,29,27/
c 
        data ns3/31,36,34,41,42,41,38,39,37,47,42,
     1     42,50,46,49,41,44,42,43,44,44,45,43,43,39/
c 
c        construct the array of thresholds from the array rsmalls
c 
        if(ifcalled .eq. 1) goto 1300
        thresh(1,1)=rsmalls(1)
        do 1200 i=1,25
c 
        thresh(2,i)=rsmalls(i+1)
        thresh(1,i+1)=rsmalls(i+1)
c 
 1200 continue
        ifcalled=1
 1300 continue
c 
c        find the appropriate storade subroutine to call
c 
        x0=1
        y0=0
c 
        d=sqrt((x-x0)**2+(y-y0)**2)
        do 1400 i=1,25
c 
        ier=0
        ii=-7
        ii=i
        if( (d .ge. thresh(1,i)) .and.
     1      (d .le. thresh(2,i) ) ) goto 1600
 1400 continue
c 
        ier=4
        return
 1600 continue
  
        call prinf('ii as found*',ii,1)
c 
c       return to the user the appropriate weights and nodes
c 
        if(itype .eq. 1) n=ns1(ii)
        if(itype .eq. 2) n=ns2(ii)
        if(itype .eq. 3) n=ns3(ii)
  
  
        call prinf('in o1rnear, n=*',n,1)
c 
        if(ii .eq. 1) call o1rnear1(n,xs,whts)
        if(ii .eq. 2) call o1rnear2(n,xs,whts)
        if(ii .eq. 3) call o1rnear3(n,xs,whts)
        if(ii .eq. 4) call o1rnear4(n,xs,whts)
        if(ii .eq. 5) call o1rnear5(n,xs,whts)
        if(ii .eq. 6) call o1rnear6(n,xs,whts)
        if(ii .eq. 7) call o1rnear7(n,xs,whts)
        if(ii .eq. 8) call o1rnear8(n,xs,whts)
        if(ii .eq. 9) call o1rnear9(n,xs,whts)
        if(ii .eq. 10) call o1rnear10(n,xs,whts)
c 
        if(ii .eq. 11) call o1rnear11(n,xs,whts)
        if(ii .eq. 12) call o1rnear12(n,xs,whts)
        if(ii .eq. 13) call o1rnear13(n,xs,whts)
        if(ii .eq. 14) call o1rnear14(n,xs,whts)
        if(ii .eq. 15) call o1rnear15(n,xs,whts)
        if(ii .eq. 16) call o1rnear16(n,xs,whts)
        if(ii .eq. 17) call o1rnear17(n,xs,whts)
        if(ii .eq. 18) call o1rnear18(n,xs,whts)
        if(ii .eq. 19) call o1rnear19(n,xs,whts)
        if(ii .eq. 20) call o1rnear20(n,xs,whts)
c 
        if(ii .eq. 21) call o1rnear21(n,xs,whts)
        if(ii .eq. 22) call o1rnear22(n,xs,whts)
        if(ii .eq. 23) call o1rnear23(n,xs,whts)
        if(ii .eq. 24) call o1rnear24(n,xs,whts)
        if(ii .eq. 25) call o1rnear25(n,xs,whts)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o1rwrt2(iw)
        implicit real *8 (a-h,o-z)
        save
        real *8 xs(20),whts(20)
        dimension numspts(25),weights(34,25),points(34,25),
     1      coefslrg(20,34,25)
c 
c        This subroutine constructs (by calling the subroutine o1rmatbg)
c        and stores on disk on the FORTRAN unit number iw all of the
c        parameters used for the construction of the matrix totmat. More
c        specifically, it constructs the support nodes and corresponding
c        weights, auxiliary quadrature nodes for each of the support
c        nodes, interpolation coefficients from the support nodes of the
c        auxiliary nodes, and the integer array numspts whose i-th entry
c        is the number of auxiliary quadrature nodes for the support
c        node number i.
c 
        call o1rmatb2(coefslrg,xs,whts,points,weights,numspts)
c 
        write(iw) xs
        write(iw) whts
        write(iw) points
        write(iw) weights
        write(iw) numspts
        write(iw) coefslrg
  
        call prin2('xs=*',xs,20)
        call prin2('whts=*',whts,20)
        call prin2('points=*',points,24*20)
        call prin2('weights=*',weights,24*20)
        call prinf('numspts=*',numspts,20)
        call prin2('coefslrg=*',coefslrg,20*20*24)
        return
        end
c 
c 
c 
c 
c 
        subroutine o1rrd2(ir,xs,whts,points,weights,
     1      numspts,coefslrg)
        implicit real *8 (a-h,o-z)
        save
        real *8 xs(20),whts(20)
        dimension numspts(25),weights(34,25),points(34,25),
     1      coefslrg(20,34,25)
c 
c        This subroutine reads from the FORTRAN unit ir the various
c        types of data to be used in the construction of the "near"
c        quadratures. The data to be read have been (hopefully) stored
c        previously by a call to the subroutine o1rwrt2 (see).
c 
c                        Input parameters:
c 
c  ir - the FORTRAN unit number from which the data are to be read
c 
c                 Output parameters:
c 
c  coefslrg - a 20 \times 24 \times 20 - array, containing the set of 20
c       matrices described above.
c  xs - the 20 support nodes on the interval [-1,1]
c  whts - quadrature weights (pretty universal) corresponding to the
c       nodes xs.
c  points - the array containing the quadrature nodes corresponding to
c       all support nodes; the quadrature nodes corresponding to the
c       i-th support nodes are returned in the i-th column of the array
c       points. Please note that the numbers of quadrature nodes are
c       different for different support nodes, so that a part of the
c       array points is returned empty; the number of quadrature nodes
c       actually stored in the i-th column of the array points is
c       returned in the i-th location of the array numspts.
c  weights - the array containing the quadrature weights corresponding to
c       the nodes return in the array points. One weight correponds to
c       each of the nodes in the array points, so that the structures of
c       these two arrays are identical.
c  numspts - integer array of length 20, containing numbers of quadrature
c       nodes corresponding to the 20 support nodes. Specifically,
c       numspts(i) is the number of quadrature nodes corresponding to
c       the support node number i
  
        read(ir) xs
        read(ir) whts
        read(ir) points
        read(ir) weights
        read(ir) numspts
        read(ir) coefslrg
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o1rnear1(n,xs,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),whts(1),
     1      xs31(31),whts31(31),xs23(23),whts23(23),
     2      xs13(13),whts13(13)
c 
c     Accuracy for the following 31 nodes is about 9.5 digits
c 
        data xs31/
     1  -.9981128616396147E+00,-.9844905261767816E+00,
     2  -.9470206532992628E+00,-.8768076410870395E+00,
     3  -.7698390072017736E+00,-.6268159919761789E+00,
     4  -.4524820039245667E+00,-.2547793479531651E+00,
     5  -.4389163988136613E-01,0.1688033139339491E+00,
     6  0.3718313260561628E+00,0.5547184565768457E+00,
     7  0.7090305529423144E+00,0.8293032497668538E+00,
     8  0.9138116899334188E+00,0.9650948861194200E+00,
     9  0.9899763454156061E+00,0.9983605338448163E+00,
     *  0.9998841099639191E+00,0.9999906258383526E+00,
     *  0.9999968044954566E+00,0.9999978144484076E+00,
     *  0.9999981156270332E+00,0.9999983067293361E+00,
     *  0.9999984775835582E+00,0.9999986343196710E+00,
     *  0.9999987786294895E+00,0.9999989135774598E+00,
     *  0.9999990771536500E+00,0.9999993978463659E+00,
     *  0.9999998398713417E+00/
c 
        data whts31/
     1  0.5686333435430167E-02,0.2368593285186501E-01,
     2  0.5275619155011833E-01,0.8833315367639233E-01,
     3  0.1254755327380922E+00,0.1597658077314174E+00,
     4  0.1875624423212389E+00,0.2061280534934053E+00,
     5  0.2137326372045507E+00,0.2097264920684806E+00,
     6  0.1945664242326698E+00,0.1697865336506889E+00,
     7  0.1379139458262731E+00,0.1023347505107776E+00,
     8  0.6710334718749444E-01,0.3661706970150606E-01,
     9  0.1485419564763227E-01,0.3602792164942678E-02,
     *  0.3479661496927907E-03,0.1596230334522242E-04,
     *  0.1969948612333244E-05,0.4694418938561806E-06,
     *  0.2121168529585860E-06,0.1787602873711198E-06,
     *  0.1634407692332152E-06,0.1502501238767955E-06,
     *  0.1386576492769440E-06,0.1351853649541487E-06,
     *  0.2185766509266733E-06,0.4238798927426579E-06,
     *  0.3753042647288199E-06/
c 
c     Accuracy for the following 23 nodes is 0.31673E-07
c 
        data xs23/
     1  -.9966615621044050E+00,-.9746678957403749E+00,
     2  -.9195553453985501E+00,-.8238614156405282E+00,
     3  -.6867330189551339E+00,-.5125574493978424E+00,
     4  -.3097304974372213E+00,-.8947516282487222E-01,
     5  0.1353761666944833E+00,0.3516450934437583E+00,
     6  0.5470787104100530E+00,0.7115791268478617E+00,
     7  0.8383494945432266E+00,0.9249407205227362E+00,
     8  0.9741831155569070E+00,0.9947210690036457E+00,
     9  0.9995768968820240E+00,0.9999872020307243E+00,
     *  0.9999974587529132E+00,0.9999982297655439E+00,
     *  0.9999986135918029E+00,0.9999989702394969E+00,
     *  0.9999996212663402E+00/
c 
        data whts23/
     1  0.9811461988850771E-02,0.3669138416006996E-01,
     2  0.7475758116414846E-01,0.1167489085373626E+00,
     3  0.1567565253248224E+00,0.1901758288129171E+00,
     4  0.2135872706668617E+00,0.2247634872255227E+00,
     5  0.2227214282131740E+00,0.2077536577973757E+00,
     6  0.1814114493350821E+00,0.1464391551467838E+00,
     7  0.1066775500245265E+00,0.6696483810083773E-01,
     8  0.3297427382841737E-01,0.1038206238971981E-01,
     9  0.1341318640732786E-02,0.3814212889104881E-04,
     *  0.1701986072103584E-05,0.4112131569819004E-06,
     *  0.3662263516989166E-06,0.4239976544005720E-06,
     *  0.7830379404950415E-06/
c 
c     Accuracy for the following 13 nodes is 0.66473E-04
c 
        data xs13/
     1  -.9847619910214469E+00,-.9187399953832465E+00,
     2  -.7988524607832725E+00,-.6289824889296768E+00,
     3  -.4186297973958634E+00,-.1811874899239823E+00,
     4  0.6767919294505775E-01,0.3112810652987894E+00,
     5  0.5336294443485752E+00,0.7206551446419631E+00,
     6  0.8618175384963016E+00,0.9518141009172640E+00,
     7  0.9930479890833670E+00/
c 
        data whts13/
     1  0.3925528214283325E-01,0.9308391338820939E-01,
     2  0.1460379596642505E+00,0.1920750379122539E+00,
     3  0.2263593950793901E+00,0.2458587499495450E+00,
     4  0.2489195130477837E+00,0.2355220338673511E+00,
     5  0.2067911421474080E+00,0.1654687442054645E+00,
     6  0.1159394042506011E+00,0.6434857711041762E-01,
     7  0.2056682299823517E-01/
c 
        if(n .eq. 31) call o1rarrcp(xs31,xs,n)
        if(n .eq. 31) call o1rarrcp(whts31,whts,n)
c 
        if(n .eq. 23) call o1rarrcp(xs23,xs,n)
        if(n .eq. 23) call o1rarrcp(whts23,whts,n)
c 
        if(n .eq. 13) call o1rarrcp(xs13,xs,n)
        if(n .eq. 13) call o1rarrcp(whts13,whts,n)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o1rnear2(n,xs,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),whts(1),
     1      xs36(36),whts36(36),xs23(23),whts23(23),
     2      xs14(14),whts14(14)
c 
c     Accuracy for the following 23 nodes is 0.51239E-07
c 
        data xs23/
     1  -.9965645248482821E+00,-.9740392506674570E+00,
     2  -.9178510322028485E+00,-.8206400580712958E+00,
     3  -.6817409604748312E+00,-.5057526329171300E+00,
     4  -.3012773488675881E+00,-.7972187755967008E-01,
     5  0.1459323376800241E+00,0.3624018105641638E+00,
     6  0.5573810082826585E+00,0.7207801757937395E+00,
     7  0.8458810345474093E+00,0.9303963590159263E+00,
     8  0.9774254710422270E+00,0.9960332277397868E+00,
     9  0.9997880602352909E+00,0.9999894522844417E+00,
     *  0.9999956836630226E+00,0.9999966694871398E+00,
     *  0.9999973436560752E+00,0.9999980221504936E+00,
     *  0.9999992908089117E+00/
c 
        data whts23/
     1  0.1008285726148629E-01,0.3749995856487779E-01,
     2  0.7608107106097164E-01,0.1184273147076740E+00,
     3  0.1585838092117472E+00,0.1919393604949596E+00,
     4  0.2150895802251244E+00,0.2258364289317264E+00,
     5  0.2232365324644566E+00,0.2076307270767622E+00,
     6  0.1806267368235726E+00,0.1450342683374372E+00,
     7  0.1047700308278923E+00,0.6476462335362880E-01,
     8  0.3081728526925711E-01,0.8778718049652176E-02,
     9  0.7773329779728000E-03,0.1805666987411377E-04,
     *  0.1669391461674336E-05,0.7251725486347617E-06,
     *  0.6432677199535479E-06,0.8530295849043926E-06,
     *  0.1483079194304416E-05/
c 
c     Accuracy for the following 14 nodes is  0.21262E-04
c 
        data xs14/
     1  -.9878399644852978E+00,-.9321258231399322E+00,
     2  -.8251137902723260E+00,-.6680422413673242E+00,
     3  -.4687005520060394E+00,-.2392766563994607E+00,
     4  0.5223179080591667E-02,0.2486068759324040E+00,
     5  0.4749236745709073E+00,0.6699559179139825E+00,
     6  0.8227086939644224E+00,0.9268738638283427E+00,
     7  0.9825629018687287E+00,0.9994616130693668E+00/
c 
        data whts14/
     1  0.3179048443051284E-01,0.8087725324803156E-01,
     2  0.1328867494259591E+00,0.1799253340219636E+00,
     3  0.2167024687039262E+00,0.2396172550116982E+00,
     4  0.2467604324860493E+00,0.2374168832559767E+00,
     5  0.2128234240264862E+00,0.1753859926288254E+00,
     6  0.1290398916694222E+00,0.7926950768461478E-01,
     7  0.3360780503145339E-01,0.3745588649704235E-02/
c 
c     Accuracy for the following 36 nodes is 0.99884E-10
c 
        data xs36/
     1  -.9983772278289977E+00,-.9864494116631674E+00,
     2  -.9529551877493310E+00,-.8890147246851811E+00,
     3  -.7900462238948965E+00,-.6558881894930619E+00,
     4  -.4903106327041860E+00,-.3002935764805921E+00,
     5  -.9516331892654206E-01,0.1143781741872285E+00,
     6  0.3172893579536518E+00,0.5032371644873641E+00,
     7  0.6635981243025662E+00,0.7923556060283141E+00,
     8  0.8868390252758898E+00,0.9482331898480540E+00,
     9  0.9816903867635205E+00,0.9956122600933353E+00,
     *  0.9993964806680475E+00,0.9999411482913596E+00,
     *  0.9999873990601245E+00,0.9999937731595161E+00,
     *  0.9999954102532111E+00,0.9999959951994883E+00,
     *  0.9999963208630169E+00,0.9999965915929701E+00,
     *  0.9999968411205213E+00,0.9999970744377382E+00,
     *  0.9999972932299888E+00,0.9999974986837878E+00,
     *  0.9999976923146381E+00,0.9999978796673343E+00,
     *  0.9999980986732990E+00,0.9999984772727477E+00,
     *  0.9999991146157421E+00,0.9999997872371775E+00/
c 
        data whts36/
     1  0.4912315493827864E-02,0.2092025044466446E-01,
     2  0.4757419424625341E-01,0.8106985256139983E-01,
     3  0.1168852407556400E+00,0.1507903288148391E+00,
     4  0.1791869841463837E+00,0.1992764700334052E+00,
     5  0.2091829852633154E+00,0.2080421322785385E+00,
     6  0.1960412351266207E+00,0.1744020557480981E+00,
     7  0.1453046269504940E+00,0.1117564239336102E+00,
     8  0.7740507107098676E-01,0.4625188995082587E-01,
     9  0.2209954295371286E-01,0.7377770434464035E-02,
     *  0.1376127765742079E-02,0.1235743977232746E-03,
     *  0.1300074025897941E-04,0.2818636792863374E-05,
     *  0.8897057842493741E-06,0.3928085016470643E-06,
     *  0.2865228237087405E-06,0.2586203594425509E-06,
     *  0.2410511283984648E-06,0.2258404885261550E-06,
     *  0.2119330057681412E-06,0.1991977453863390E-06,
     *  0.1886766161277666E-06,0.1906741487361456E-06,
     *  0.2706047184182870E-06,0.5091413762232897E-06,
     *  0.7298073370397874E-06,0.5138135612075882E-06/
c 
        if(n .eq. 23) call o1rarrcp(xs23,xs,n)
        if(n .eq. 23) call o1rarrcp(whts23,whts,n)
c 
        if(n .eq. 14) call o1rarrcp(xs14,xs,n)
        if(n .eq. 14) call o1rarrcp(whts14,whts,n)
c 
        if(n .eq. 36) call o1rarrcp(xs36,xs,n)
        if(n .eq. 36) call o1rarrcp(whts36,whts,n)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o1rnear3(n,xs,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),whts(1),
     1      xs34(34),whts34(34),xs23(23),whts23(23),
     2      xs15(15),whts15(15)
c 
c     Accuracy for the following 23 nodes is 0.51239E-07
c 
        data xs23/
     1  -.9961362500918736E+00,-.9713410572792273E+00,
     2  -.9106957025107846E+00,-.8073499537850994E+00,
     3  -.6614665790453957E+00,-.4785596805482994E+00,
     4  -.2681120302268766E+00,-.4228645635908257E-01,
     5  0.1853683769689063E+00,0.4012336307040569E+00,
     6  0.5929319665978940E+00,0.7505993050889477E+00,
     7  0.8680999516153953E+00,0.9441923196881499E+00,
     8  0.9836108808273745E+00,0.9974257542433038E+00,
     9  0.9998334847910908E+00,0.9999843908190218E+00,
     *  0.9999921876029154E+00,0.9999937228674489E+00,
     *  0.9999949040451552E+00,0.9999961778519091E+00,
     *  0.9999986691581317E+00/
c 
        data whts23/
     1  0.1126857975306276E-01,0.4090247745955462E-01,
     2  0.8149949957729889E-01,0.1251370227008908E+00,
     3  0.1656882484322368E+00,0.1985238155917020E+00,
     4  0.2203227906177098E+00,0.2290482587606222E+00,
     5  0.2239696650463769E+00,0.2056767225901795E+00,
     6  0.1760557005273981E+00,0.1382369789059938E+00,
     7  0.9653859977296875E-01,0.5642305690321059E-01,
     8  0.2429833689337336E-01,0.5874349163923770E-02,
     9  0.5053808026674113E-03,0.2129727001058657E-04,
     *  0.2274640806552941E-05,0.1282309758261828E-05,
     *  0.1102746912079082E-05,0.1727749490163716E-05,
     *  0.2819572899405569E-05/
c 
c     Accuracy for the following 15 nodes is 0.11166E-04
c 
        data xs15/
     1  -.9900377201633368E+00,-.9415031481628020E+00,
     2  -.8438954322313958E+00,-.6970111962338310E+00,
     3  -.5073661086798840E+00,-.2859636842214275E+00,
     4  -.4676094604585357E-01,0.1947663087323748E+00,
     5  0.4230273036416003E+00,0.6237864156824986E+00,
     6  0.7855499855355532E+00,0.9010454822059868E+00,
     7  0.9689104490071425E+00,0.9960195936414329E+00,
     8  0.9999690214895467E+00/
c 
        data whts15/
     1  0.2653183008135897E-01,0.7229444037187872E-01,
     2  0.1228736882340916E+00,0.1697777545260704E+00,
     3  0.2076519265465987E+00,0.2328019219283344E+00,
     4  0.2429805160403536E+00,0.2374211226712534E+00,
     5  0.2167016599547381E+00,0.1828755507624078E+00,
     6  0.1394126571562016E+00,0.9131175798771200E-01,
     7  0.4551255049014772E-01,0.1173077841669093E-01,
     8  0.1653926748586335E-03/
c 
c     Accuracy for the following 34 nodes is 0.40276E-09
c 
        data xs34/
     1  -.9981710219118587E+00,-.9849180222299289E+00,
     2  -.9483074460120411E+00,-.8794443294163139E+00,
     3  -.7741971229279815E+00,-.6330884539821686E+00,
     4  -.4606593975401605E+00,-.2646503227463354E+00,
     5  -.5506270749663167E-01,0.1568717715696286E+00,
     6  0.3597761981752989E+00,0.5432156014048830E+00,
     7  0.6987297203036087E+00,0.8207572596508621E+00,
     8  0.9073984276386591E+00,0.9609288965967662E+00,
     9  0.9878164523439925E+00,0.9975831389442149E+00,
     *  0.9997169505554411E+00,0.9999599408880351E+00,
     *  0.9999853648205821E+00,0.9999903934737175E+00,
     *  0.9999919022799062E+00,0.9999926301041433E+00,
     *  0.9999932208404529E+00,0.9999937643575165E+00,
     *  0.9999942703103304E+00,0.9999947422963598E+00,
     *  0.9999951833221332E+00,0.9999955984062336E+00,
     *  0.9999960293952490E+00,0.9999967304591757E+00,
     *  0.9999980244085394E+00,0.9999995102668389E+00/
c 
        data whts34/
     1  0.5516495747927969E-02,0.2308561988599794E-01,
     2  0.5163921555604958E-01,0.8677117836839543E-01,
     3  0.1236249804346497E+00,0.1578217771652442E+00,
     4  0.1857305868637640E+00,0.2046033741364937E+00,
     5  0.2126817000907285E+00,0.2092728465869099E+00,
     6  0.1947795024893955E+00,0.1706737390259073E+00,
     7  0.1394147435005426E+00,0.1043146638918951E+00,
     8  0.6934517883814777E-01,0.3880983687001249E-01,
     9  0.1660540524125302E-01,0.4585341772572744E-02,
     *  0.6447390050428063E-03,0.5864090601512834E-04,
     *  0.9360513135705029E-05,0.2464511738882735E-05,
     *  0.9189334363010587E-06,0.6253451737979969E-06,
     *  0.5642190158311214E-06,0.5240034740411478E-06,
     *  0.4884659239434132E-06,0.4559975149440276E-06,
     *  0.4267466033721926E-06,0.4067550158170924E-06,
     *  0.4975360207365696E-06,0.9753762231747796E-06,
     *  0.1553648088515432E-05,0.1172222809972602E-05/
c 
        if(n .eq. 23) call o1rarrcp(xs23,xs,n)
        if(n .eq. 23) call o1rarrcp(whts23,whts,n)
c 
        if(n .eq. 15) call o1rarrcp(xs15,xs,n)
        if(n .eq. 15) call o1rarrcp(whts15,whts,n)
c 
        if(n .eq. 34) call o1rarrcp(xs34,xs,n)
        if(n .eq. 34) call o1rarrcp(whts34,whts,n)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o1rnear4(n,xs,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),whts(1),
     1      xs25(25),whts25(25),xs41(41),whts41(41),
     2      xs15(15),whts15(15)
c 
c     Accuracy for the following 15 nodes is 0.11166E-04
c 
        data xs15/
     1  -.9890768797788434E+00,-.9374236044884182E+00,
     2  -.8357159123196219E+00,-.6844048349185826E+00,
     3  -.4906174636487516E+00,-.2659140942125569E+00,
     4  -.2469897937549848E-01,0.2172787731358468E+00,
     5  0.4443495698088205E+00,0.6423758388848708E+00,
     6  0.8001602577449412E+00,0.9108893826901951E+00,
     7  0.9738586059413237E+00,0.9971542231456549E+00,
     8  0.9999754049579632E+00/
c 
        data whts15/
     1  0.2882291208120506E-01,0.7602660071333014E-01,
     2  0.1272415756884270E+00,0.1741606349035590E+00,
     3  0.2114608828745896E+00,0.2355186945735555E+00,
     4  0.2442588884418052E+00,0.2370557844468295E+00,
     5  0.2147056622643672E+00,0.1794566625682370E+00,
     6  0.1349496630685366E+00,0.8634677801314870E-01,
     7  0.4088356529226365E-01,0.9026614374963623E-02,
     8  0.9932803274420720E-04/
c 
c     Accuracy for the following 25 nodes is better than 0.39815E-07
c 
        data xs25/
     1  -.9966967455055917E+00,-.9748886172647956E+00,
     2  -.9201293141908157E+00,-.8248955299439776E+00,
     3  -.6882524948965039E+00,-.5145106802751545E+00,
     4  -.3120043994675942E+00,-.9191515974305753E-01,
     5  0.1329444626266816E+00,0.3493936850702681E+00,
     6  0.5451549838680865E+00,0.7100839532534782E+00,
     7  0.8373170129869882E+00,0.9243240415112257E+00,
     8  0.9738571794786985E+00,0.9945304587860095E+00,
     9  0.9994513114010948E+00,0.9999404605854598E+00,
     *  0.9999777496600979E+00,0.9999837552721138E+00,
     *  0.9999863743764827E+00,0.9999888685324001E+00,
     *  0.9999910664527221E+00,0.9999937830988086E+00,
     *  0.9999981659382917E+00/
c 
        data whts25/
     1  0.9713934685982570E-02,0.3641582272477426E-01,
     2  0.7433773803331996E-01,0.1162624515486546E+00,
     3  0.1562852194438519E+00,0.1897899428842759E+00,
     4  0.2133387527984683E+00,0.2246831212004631E+00,
     5  0.2228179719822080E+00,0.2080135957089937E+00,
     6  0.1817991010007123E+00,0.1468974022045207E+00,
     7  0.1071309349037820E+00,0.6732934208518024E-01,
     8  0.3318468340432550E-01,0.1045735043786750E-01,
     9  0.1418977470933769E-02,0.9393024827118409E-04,
     *  0.1179521431908779E-04,0.3072101981333349E-05,
     *  0.2514189376891755E-05,0.2399168116411675E-05,
     *  0.2102637540412266E-05,0.3710707747149037E-05,
     *  0.4150378352293711E-05/
c 
        data xs41/
     1  -.9985878094809333E+00,-.9880800414000348E+00,
     2  -.9581314803409106E+00,-.9001350630770909E+00,
     3  -.8091885557209676E+00,-.6844151976901676E+00,
     4  -.5286408510906125E+00,-.3478053472695580E+00,
     5  -.1502159614065257E+00,0.5431517765165261E-01,
     6  0.2554096676673345E+00,0.4430878329690033E+00,
     7  0.6087045609469190E+00,0.7458000594662464E+00,
     8  0.8508113250231562E+00,0.9235813782759357E+00,
     9  0.9675524340058197E+00,0.9893824515638328E+00,
     *  0.9975189323947853E+00,0.9995616799448241E+00,
     *  0.9999055538619999E+00,0.9999627509490114E+00,
     *  0.9999768676377119E+00,0.9999816443898109E+00,
     *  0.9999836214139108E+00,0.9999847284177392E+00,
     *  0.9999856154422873E+00,0.9999864313735785E+00,
     *  0.9999872016701641E+00,0.9999879324750773E+00,
     *  0.9999886266267322E+00,0.9999892865630109E+00,
     *  0.9999899144837374E+00,0.9999905127997217E+00,
     *  0.9999910862719950E+00,0.9999916531857742E+00,
     *  0.9999923021376753E+00,0.9999933020744822E+00,
     *  0.9999949904705088E+00,0.9999973138829341E+00,
     *  0.9999993928941335E+00/
c 
        data whts41/
     1  0.4287809808909010E-02,0.1854133320467104E-01,
     2  0.4281870064761578E-01,0.7399423157570431E-01,
     3  0.1080388010503877E+00,0.1410269421817988E+00,
     4  0.1695192564658094E+00,0.1907527356105109E+00,
     5  0.2027749182223382E+00,0.2045414079926038E+00,
     6  0.1959685610288552E+00,0.1779326119404314E+00,
     7  0.1522122185812732E+00,0.1213762253995906E+00,
     8  0.8861635911475343E-01,0.5750261827296478E-01,
     9  0.3157177598972227E-01,0.1354349005952113E-01,
     *  0.4043872424759218E-02,0.7669926253905782E-03,
     *  0.1169366713608089E-03,0.2464832223550921E-04,
     *  0.7551601887767750E-05,0.2877611973825730E-05,
     *  0.1363044912265951E-05,0.9489199760546916E-06,
     *  0.8434973185947555E-06,0.7914985689288522E-06,
     *  0.7499650382081211E-06,0.7120776719258857E-06,
     *  0.6766354482352618E-06,0.6435985338972101E-06,
     *  0.6126160839485762E-06,0.5847000572478811E-06,
     *  0.5645588810540337E-06,0.5813134715906090E-06,
     *  0.7628465169868029E-06,0.1301751539072390E-05,
     *  0.2075613211964448E-05,0.2419768569911698E-05,
     *  0.1491137719776313E-05/
c 
        if(n .eq. 25) call o1rarrcp(xs25,xs,n)
        if(n .eq. 25) call o1rarrcp(whts25,whts,n)
c 
        if(n .eq. 15) call o1rarrcp(xs15,xs,n)
        if(n .eq. 15) call o1rarrcp(whts15,whts,n)
c 
        if(n .eq. 41) call o1rarrcp(xs41,xs,n)
        if(n .eq. 41) call o1rarrcp(whts41,whts,n)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o1rnear5(n,xs,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),whts(1),
     1      xs25(25),whts25(25),xs42(42),whts42(42),
     2      xs15(15),whts15(15)
c 
c     Accuracy for the following 15 nodes is 0.20405E-04
c 
        data xs15/
     1  -.9880489827243415E+00,-.9330233659285151E+00,
     2  -.8269025827392392E+00,-.6707911838567631E+00,
     3  -.4723763126718458E+00,-.2437528829066158E+00,
     4  0.1672595244073261E-03,0.2432612192878219E+00,
     5  0.4696339707109045E+00,0.6650937903853229E+00,
     6  0.8186318017674086E+00,0.9238745012516816E+00,
     7  0.9808147202099873E+00,0.9989026776415575E+00,
     8  0.9999743729970705E+00/
c 
        data whts15/
     1  0.3128778953976895E-01,0.8005703442050530E-01,
     2  0.1319431343926215E+00,0.1789655202416013E+00,
     3  0.2158231990136662E+00,0.2389105402334178E+00,
     4  0.2462736098006464E+00,0.2372936702665105E+00,
     5  0.2130633051701266E+00,0.1759985899866526E+00,
     6  0.1299872096450331E+00,0.8045516041594525E-01,
     7  0.3490332722583990E-01,0.4913938127396740E-02,
     8  0.4739698237388969E-04/
c 
c     Accuracy for the following 25 nodes is 0.39815E-07
c 
        data xs25/
     1  -.9967202848700426E+00,-.9750431375327792E+00,
     2  -.9205535894010555E+00,-.8257074378333911E+00,
     3  -.6895268067317051E+00,-.5162714645452304E+00,
     4  -.3142240852302673E+00,-.9451767599530092E-01,
     5  0.1300773995098843E+00,0.3464134446150016E+00,
     6  0.5422341601825968E+00,0.7074018609838697E+00,
     7  0.8350423480160103E+00,0.9225934870061305E+00,
     8  0.9727450827245214E+00,0.9939990244986001E+00,
     9  0.9992972368110185E+00,0.9999036400719004E+00,
     *  0.9999592593475664E+00,0.9999690828064325E+00,
     *  0.9999737874588339E+00,0.9999783084088223E+00,
     *  0.9999823551397342E+00,0.9999876773156913E+00,
     *  0.9999963721135734E+00/
c 
        data whts25/
     1  0.9647808801565419E-02,0.3621505589343760E-01,
     2  0.7400300083436416E-01,0.1158293231034786E+00,
     3  0.1558020916004000E+00,0.1893087981490388E+00,
     4  0.2129102298290880E+00,0.2243530702558781E+00,
     5  0.2226244999709973E+00,0.2079842970066668E+00,
     6  0.1819483994806675E+00,0.1472238562612534E+00,
     7  0.1076141321054303E+00,0.6792430362425141E-01,
     8  0.3380837987808998E-01,0.1096496644032555E-01,
     9  0.1650711246909451E-02,0.1346027313397649E-03,
     *  0.1862689624943742E-04,0.5376228370227973E-05,
     *  0.4553662088812448E-05,0.4358875535936972E-05,
     *  0.3977858559556231E-05,0.7373045285805763E-05,
     *  0.8215531280390175E-05/
c 
c     Accuracy for the following 42 nodes is 0.33129E-10
c 
        data xs42/
     1  -.9984752719123662E+00,-.9872089483222950E+00,
     2  -.9553730158080925E+00,-.8942374804536306E+00,
     3  -.7991087456880392E+00,-.6695337613413506E+00,
     4  -.5088764378348279E+00,-.3236474618830540E+00,
     5  -.1226870905758461E+00,0.8375997620141984E-01,
     6  0.2850361107950476E+00,0.4710762893943901E+00,
     7  0.6333695764840737E+00,0.7658188146412299E+00,
     8  0.8654412659721460E+00,0.9328362866794053E+00,
     9  0.9722722542367914E+00,0.9910670788151386E+00,
     *  0.9977971512780106E+00,0.9995063421502905E+00,
     *  0.9998528741082113E+00,0.9999312637437738E+00,
     *  0.9999546342329384E+00,0.9999632788965436E+00,
     *  0.9999670375127057E+00,0.9999691769905119E+00,
     *  0.9999708753803860E+00,0.9999724295094916E+00,
     *  0.9999738973700063E+00,0.9999752941318438E+00,
     *  0.9999766271749229E+00,0.9999779015724555E+00,
     *  0.9999791199499661E+00,0.9999802833375004E+00,
     *  0.9999813940398283E+00,0.9999824617307549E+00,
     *  0.9999835356729883E+00,0.9999848441589014E+00,
     *  0.9999869449895338E+00,0.9999903894426032E+00,
     *  0.9999949327340329E+00,0.9999988687130312E+00/
c 
        data whts42/
     1  0.4621587803938212E-02,0.1981103376898021E-01,
     2  0.4534291767446982E-01,0.7771137476381454E-01,
     3  0.1126138367144667E+00,0.1459641483627248E+00,
     4  0.1742516885916643E+00,0.1947173704677389E+00,
     5  0.2054792695548559E+00,0.2056216930635331E+00,
     6  0.1952369092145502E+00,0.1754104993433131E+00,
     7  0.1481480489340249E+00,0.1162452060845857E+00,
     8  0.8309740104493285E-01,0.5241299435344087E-01,
     9  0.2770806019630219E-01,0.1136268316802318E-01,
     *  0.3309353288135209E-02,0.6977833330703577E-03,
     *  0.1444047226849663E-03,0.3867078388968266E-04,
     *  0.1328745389923109E-04,0.5378720985603254E-05,
     *  0.2632396865463698E-05,0.1826448787808305E-05,
     *  0.1608613823334316E-05,0.1507068739390921E-05,
     *  0.1430801962615892E-05,0.1363893665187063E-05,
     *  0.1303056578133816E-05,0.1246160180186947E-05,
     *  0.1190714924884866E-05,0.1136382925935161E-05,
     *  0.1086280030354249E-05,0.1055280256626577E-05,
     *  0.1125773931981147E-05,0.1590783006632469E-05,
     *  0.2716338128553401E-05,0.4147159726347520E-05,
     *  0.4643421082210699E-05,0.2787866352149769E-05/
c 
        if(n .eq. 25) call o1rarrcp(xs25,xs,n)
        if(n .eq. 25) call o1rarrcp(whts25,whts,n)
c 
        if(n .eq. 15) call o1rarrcp(xs15,xs,n)
        if(n .eq. 15) call o1rarrcp(whts15,whts,n)
c 
        if(n .eq. 42) call o1rarrcp(xs42,xs,n)
        if(n .eq. 42) call o1rarrcp(whts42,whts,n)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o1rnear6(n,xs,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),whts(1),
     1      xs25(25),whts25(25),xs32(32),whts32(32),
     2      xs15(15),whts15(15),xs41(41),whts41(41)
c 
c     Accuracy for the following 15 nodes is 0.28622E-04
c 
        data xs15/
     1  -.9880198193945755E+00,-.9329173194961299E+00,
     2  -.8267352242897454E+00,-.6705958499363720E+00,
     3  -.4721900835668694E+00,-.2436099255411575E+00,
     4  0.2398279814028736E-03,0.2432439298803519E+00,
     5  0.4695180424179337E+00,0.6648838188317723E+00,
     6  0.8183499146366353E+00,0.9235653300000810E+00,
     7  0.9805527592437113E+00,0.9987790088921580E+00,
     8  0.9999520906850281E+00/
c 
        data whts15/
     1  0.3135467048507572E-01,0.8013171768852304E-01,
     2  0.1319889157729129E+00,0.1789748989099373E+00,
     3  0.2157959381781963E+00,0.2388520381556329E+00,
     4  0.2461868331890938E+00,0.2371970151224504E+00,
     5  0.2129643885266592E+00,0.1759121917436614E+00,
     6  0.1299334585683870E+00,0.8045912025750147E-01,
     7  0.3499891592407035E-01,0.5108299964164832E-02,
     8  0.7520879079617795E-04/
c 
c     Accuracy for the following 25 nodes is 0.71527E-07
c 
        data xs25/
     1  -.9963846553391879E+00,-.9729072857852216E+00,
     2  -.9148771495767487E+00,-.8152051826309615E+00,
     3  -.6736306122424340E+00,-.4951722884916649E+00,
     4  -.2888067019865290E+00,-.6622954021654570E-01,
     5  0.1593983457321493E+00,0.3747385243657300E+00,
     6  0.5675673972287762E+00,0.7280310468527686E+00,
     7  0.8498282152672318E+00,0.9312987742156316E+00,
     8  0.9763450897996706E+00,0.9946496307589938E+00,
     9  0.9992067270851752E+00,0.9998432132165573E+00,
     *  0.9999229457120200E+00,0.9999387288127897E+00,
     *  0.9999476570466082E+00,0.9999566280223113E+00,
     *  0.9999647232967619E+00,0.9999754377108636E+00,
     *  0.9999928190663840E+00/
c 
        data whts25/
     1  0.1058133383793492E-01,0.3892125931898791E-01,
     2  0.7829437260313668E-01,0.1210668063705910E+00,
     3  0.1612243240550845E+00,0.1941758186784442E+00,
     4  0.2165851954728170E+00,0.2263540826444911E+00,
     5  0.2226569285823824E+00,0.2059651551338920E+00,
     6  0.1780305032012778E+00,0.1418283874910818E+00,
     7  0.1014731087836125E+00,0.6211033121506134E-01,
     8  0.2962816314924737E-01,0.9275745065567307E-02,
     9  0.1553380090262816E-02,0.1800247602807994E-03,
     *  0.2881425608729781E-04,0.9528543282450534E-05,
     *  0.8954388817674515E-05,0.8697378901760226E-05,
     *  0.7984716374507983E-05,0.1483100955157301E-04,
     *  0.1631260993969371E-04/
c 
c     Accuracy for the following 32 nodes is 0.14076E-08
c 
        data xs32/
     1  -.9977750129937136E+00,-.9820811392376523E+00,
     2  -.9399962851240464E+00,-.8628536749051005E+00,
     3  -.7474549199410110E+00,-.5955314185296832E+00,
     4  -.4128805717110087E+00,-.2084101866266567E+00,
     5  0.6901573880504275E-02,0.2211207457017101E+00,
     6  0.4225090994163998E+00,0.6006592304507683E+00,
     7  0.7475583691082164E+00,0.8585346108778924E+00,
     8  0.9330412726280064E+00,0.9751678077618147E+00,
     9  0.9934413428267688E+00,0.9987466621556194E+00,
     *  0.9997202614274967E+00,0.9998840640225950E+00,
     *  0.9999232471787824E+00,0.9999361268567009E+00,
     *  0.9999431915236454E+00,0.9999488585286844E+00,
     *  0.9999532003901521E+00,0.9999578649923020E+00,
     *  0.9999619072408896E+00,0.9999656722120758E+00,
     *  0.9999681503857127E+00,0.9999750338849271E+00,
     *  0.9999856136660998E+00,0.9999965454890860E+00/
c 
        data whts32/
     1  0.6664017596078400E-02,0.2699232695615874E-01,
     2  0.5863166356284999E-01,0.9617524688540420E-01,
     3  0.1343190479540185E+00,0.1685458321342675E+00,
     4  0.1952538705837806E+00,0.2118385339035939E+00,
     5  0.2167766646326770E+00,0.2096862513855951E+00,
     6  0.1913380478403114E+00,0.1636108951100376E+00,
     7  0.1293954793067623E+00,0.9245657454538360E-01,
     8  0.5723879573746972E-01,0.2846580858443774E-01,
     9  0.9983059261112908E-02,0.2133321122669258E-02,
     *  0.3362262249505014E-03,0.6958244522100415E-04,
     *  0.2038178591561713E-04,0.8235143083623903E-05,
     *  0.6549617656154890E-05,0.4474626521703055E-05,
     *  0.4850388304768190E-05,0.3992326082464246E-05,
     *  0.4548455783619841E-05,0.1758271112286430E-05,
     *  0.4965277333848369E-05,0.8755239619665425E-05,
     *  0.1190522831028195E-04,0.8341108500711200E-05/
c 
c     Accuracy for the following 41 nodes is 0.74111E-10
c 
        data xs41/
     1  -.9983950501760241E+00,-.9865869015490831E+00,
     2  -.9533917305318137E+00,-.8899578633025609E+00,
     3  -.7916881366936189E+00,-.6583766980105352E+00,
     4  -.4937305787706400E+00,-.3046555915546012E+00,
     5  -.1003988695485474E+00,0.1084158275106764E+00,
     6  0.3108185627782356E+00,0.4965360711724824E+00,
     7  0.6569881803953728E+00,0.7861804447999166E+00,
     8  0.8814371686844644E+00,0.9438989358635859E+00,
     9  0.9786116231642803E+00,0.9937805254797332E+00,
     *  0.9985395313924108E+00,0.9996093272940010E+00,
     *  0.9998376034478849E+00,0.9999006624835334E+00,
     *  0.9999230557669914E+00,0.9999325097322658E+00,
     *  0.9999374937701446E+00,0.9999411278344013E+00,
     *  0.9999443570934439E+00,0.9999473855856480E+00,
     *  0.9999502595991701E+00,0.9999529952172286E+00,
     *  0.9999556005453700E+00,0.9999580806572717E+00,
     *  0.9999604401212265E+00,0.9999626870109097E+00,
     *  0.9999648438735700E+00,0.9999670080936911E+00,
     *  0.9999696239652371E+00,0.9999737957454806E+00,
     *  0.9999806491754566E+00,0.9999897513391216E+00,
     *  0.9999977025361133E+00/
c 
        data whts41/
     1  0.4859541384989906E-02,0.2071996910615604E-01,
     2  0.4717190197538689E-01,0.8046177119735617E-01,
     3  0.1161033327111231E+00,0.1498896607190157E+00,
     4  0.1782372486022373E+00,0.1983553602432649E+00,
     5  0.2083701526480668E+00,0.2074133941119390E+00,
     6  0.1956631209526531E+00,0.1743268631557896E+00,
     7  0.1455659097550883E+00,0.1123640346481451E+00,
     8  0.7833708113555271E-01,0.4743793305608666E-01,
     9  0.2339166781407602E-01,0.8533097484548087E-02,
     *  0.2136162277668805E-02,0.4351990359276113E-03,
     *  0.1069756728162421E-03,0.3484749698571389E-04,
     *  0.1373849485930498E-04,0.6420092451640011E-05,
     *  0.4031230267151538E-05,0.3368774753895754E-05,
     *  0.3115872915536967E-05,0.2947417315968862E-05,
     *  0.2803055280483186E-05,0.2669470558186553E-05,
     *  0.2541993548683422E-05,0.2418930274479307E-05,
     *  0.2301219590424846E-05,0.2195695547706436E-05,
     *  0.2130306464369500E-05,0.2262384416480040E-05,
     *  0.3165770678333008E-05,0.5392586703812168E-05,
     *  0.8274064267394635E-05,0.9343070153654284E-05,
     *  0.5654859094202786E-05/
c 
        if(n .eq. 25) call o1rarrcp(xs25,xs,n)
        if(n .eq. 25) call o1rarrcp(whts25,whts,n)
c 
        if(n .eq. 15) call o1rarrcp(xs15,xs,n)
        if(n .eq. 15) call o1rarrcp(whts15,whts,n)
c 
        if(n .eq. 32) call o1rarrcp(xs32,xs,n)
        if(n .eq. 32) call o1rarrcp(whts32,whts,n)
c 
        if(n .eq. 41) call o1rarrcp(xs41,xs,n)
        if(n .eq. 41) call o1rarrcp(whts41,whts,n)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o1rnear7(n,xs,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),whts(1),
     1      xs27(27),whts27(27),xs38(38),whts38(38),
     2      xs17(17),whts17(17)
c 
c     Accuracy for the following 17 nodes is 0.69154E-05
c 
        data xs17/
     1  -.9900177137363471E+00,-.9414898342290204E+00,
     2  -.8440179950898729E+00,-.6974063033972493E+00,
     3  -.5081497471627521E+00,-.2872153107859523E+00,
     4  -.4851247933417679E-01,0.1925305439122480E+00,
     5  0.4203875552604375E+00,0.6208766728586714E+00,
     6  0.7825656984319447E+00,0.8982506008142679E+00,
     7  0.9666475758087782E+00,0.9946901531638130E+00,
     8  0.9996654443905355E+00,0.9998928829405829E+00,
     9  0.9999541098930145E+00/
c 
        data whts17/
     1  0.2656648690261994E-01,0.7222869155754383E-01,
     2  0.1226676469201295E+00,0.1694422061330057E+00,
     3  0.2072166028603706E+00,0.2323089939005550E+00,
     4  0.2424952440079476E+00,0.2369707771122476E+00,
     5  0.2163547465568427E+00,0.1826926951818074E+00,
     6  0.1394576422960422E+00,0.9165932308772975E-01,
     7  0.4624037921746807E-01,0.1282026972882685E-01,
     8  0.7537609795594991E-03,0.5855565337305199E-04,
     9  0.7787418722268346E-04/
c 
c     Accuracy for the following 27 nodes is 0.41021E-07
c 
        data xs27/
     1  -.9968562208305198E+00,-.9759262692929546E+00,
     2  -.9229464625088394E+00,-.8302210468882064E+00,
     3  -.6965076084939568E+00,-.5257762966821535E+00,
     4  -.3260314201103891E+00,-.1081588316288906E+00,
     5  0.1152724392505742E+00,0.3312569902435896E+00,
     6  0.5276094641538474E+00,0.6941809103011976E+00,
     7  0.8239944773420489E+00,0.9142834921031013E+00,
     8  0.9674068007349204E+00,0.9913547134073479E+00,
     9  0.9984092085241778E+00,0.9996433553683413E+00,
     *  0.9998312329738137E+00,0.9998709030719275E+00,
     *  0.9998872993195731E+00,0.9999034698868858E+00,
     *  0.9999110986934492E+00,0.9999218135024579E+00,
     *  0.9999365136837330E+00,0.9999588933997892E+00,
     *  0.9999889967797062E+00/
c 
        data whts27/
     1  0.9267012160516900E-02,0.3507841148473401E-01,
     2  0.7214870997963756E-01,0.1134875098593003E+00,
     3  0.1532583288081413E+00,0.1868512091643359E+00,
     4  0.2108042438974771E+00,0.2228252000669038E+00,
     5  0.2218487010635486E+00,0.2080691207418691E+00,
     6  0.1829253661245623E+00,0.1490374174309986E+00,
     7  0.1101129940947999E+00,0.7084699832061081E-01,
     8  0.3674474907476611E-01,0.1329653661403178E-01,
     9  0.2781889555011882E-02,0.4009103221370236E-03,
     *  0.7405075359121438E-04,0.2023089256010551E-04,
     *  0.1575168756751966E-04,0.1519027171024972E-04,
     *  0.3897373409098834E-05,0.1456224422814601E-04,
     *  0.1648323245405728E-04,0.2880780467739289E-04,
     *  0.2571793388327816E-04/
c 
c     Accuracy for the following 38 nodes is 0.47932E-09
c 
        data xs38/
     1  -.9981681415934521E+00,-.9849026981200972E+00,
     2  -.9482835363843019E+00,-.8794465835411959E+00,
     3  -.7742929420103868E+00,-.6333697740014463E+00,
     4  -.4612300617384325E+00,-.2656113950405274E+00,
     5  -.5649709163695305E-01,0.1549142719502921E+00,
     6  0.3572923624313253E+00,0.5402599619656736E+00,
     7  0.6954226480378503E+00,0.8172887070772465E+00,
     8  0.9040241076620565E+00,0.9579505183740677E+00,
     9  0.9855277902632755E+00,0.9961547263141863E+00,
     *  0.9990311448066518E+00,0.9996497335665331E+00,
     *  0.9998007376596729E+00,0.9998486683962537E+00,
     *  0.9998673164312402E+00,0.9998771954845926E+00,
     *  0.9998848895939665E+00,0.9998918878762650E+00,
     *  0.9998984554220376E+00,0.9999046579923305E+00,
     *  0.9999105265709785E+00,0.9999160839879546E+00,
     *  0.9999213513096146E+00,0.9999263563254353E+00,
     *  0.9999311820611489E+00,0.9999363195879747E+00,
     *  0.9999438175573525E+00,0.9999569971383682E+00,
     *  0.9999763313264862E+00,0.9999945433315880E+00/
c 
        data whts38/
     1  0.5524275391470685E-02,0.2310037595417527E-01,
     2  0.5163617881903301E-01,0.8671629115259531E-01,
     3  0.1234885806588207E+00,0.1575851581766528E+00,
     4  0.1853889694095641E+00,0.2041672096065964E+00,
     5  0.2121767831290238E+00,0.2087393369579122E+00,
     6  0.1942700161459663E+00,0.1702506593846744E+00,
     7  0.1391466692013736E+00,0.1042709790385823E+00,
     8  0.6958553281825931E-01,0.3936109103079394E-01,
     9  0.1741438865225858E-01,0.5446490576853782E-02,
     *  0.1216142030779512E-02,0.2696964184438221E-03,
     *  0.7777248162618891E-04,0.2786926330972195E-04,
     *  0.1243075595696205E-04,0.8308270261770319E-05,
     *  0.7264093827390408E-05,0.6765746430513971E-05,
     *  0.6378627918677009E-05,0.6031388552946081E-05,
     *  0.5709492322555372E-05,0.5408748963430833E-05,
     *  0.5130057720621550E-05,0.4890358668005580E-05,
     *  0.4819716181227107E-05,0.5789332596413114E-05,
     *  0.9838741029958331E-05,0.1670087305696047E-04,
     *  0.2073740349590718E-04,0.1332936556578782E-04/
c 
        if(n .eq. 27) call o1rarrcp(xs27,xs,n)
        if(n .eq. 27) call o1rarrcp(whts27,whts,n)
c 
        if(n .eq. 17) call o1rarrcp(xs17,xs,n)
        if(n .eq. 17) call o1rarrcp(whts17,whts,n)
c 
        if(n .eq. 38) call o1rarrcp(xs38,xs,n)
        if(n .eq. 38) call o1rarrcp(whts38,whts,n)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o1rnear8(n,xs,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),whts(1),
     1      xs39(39),whts39(39),xs28(28),whts28(28),
     2      xs16(16),whts16(16)
c 
c     Accuracy for the following 16 nodes is 0.36472E-04
c 
        data xs16/
     1  -.9873798444493482E+00,-.9302951959366392E+00,
     2  -.8218478190920653E+00,-.6636054996945422E+00,
     3  -.4635104997925621E+00,-.2338185391984230E+00,
     4  0.1044854218898816E-01,0.2531010136865142E+00,
     5  0.4782447925753603E+00,0.6717639370085680E+00,
     6  0.8228526321019776E+00,0.9255198980894979E+00,
     7  0.9804172367379480E+00,0.9980107540188234E+00,
     8  0.9997607508937723E+00,0.9998969403068600E+00/
c 
        data whts16/
     1  0.3287457657041967E-01,0.8236799622773744E-01,
     2  0.1342215639765009E+00,0.1809040822953232E+00,
     3  0.2172166493171004E+00,0.2396311493896353E+00,
     4  0.2462036424236947E+00,0.2364487324472419E+00,
     5  0.2114567889730753E+00,0.1737555439290988E+00,
     6  0.1273928603090243E+00,0.7800713799975187E-01,
     7  0.3342237722301694E-01,0.5716721725779195E-02,
     8  0.2069788568976650E-03,0.1612286994215925E-03/
c 
c     Accuracy for the following 39 nodes is 0.50024E-09
c 
        data xs39/
     1  -.9980877034943302E+00,-.9843264501166960E+00,
     2  -.9466026437536186E+00,-.8761201026299770E+00,
     3  -.7689981893650086E+00,-.6260570425960182E+00,
     4  -.4521253276778628E+00,-.2551865448940722E+00,
     5  -.4541524668517316E-01,0.1658734784452185E+00,
     6  0.3673227521178664E+00,0.5486295552085565E+00,
     7  0.7015818814462234E+00,0.8209831251138082E+00,
     8  0.9054016426234852E+00,0.9576220080384623E+00,
     9  0.9844585949441422E+00,0.9952227277515491E+00,
     *  0.9985095401483715E+00,0.9993702988443238E+00,
     *  0.9996159420139201E+00,0.9997004019272862E+00,
     *  0.9997348475148500E+00,0.9997535992543962E+00,
     *  0.9997682552871683E+00,0.9997815864119723E+00,
     *  0.9997941161378148E+00,0.9998059763614675E+00,
     *  0.9998172268082206E+00,0.9998279091791955E+00,
     *  0.9998380594712057E+00,0.9998477172602782E+00,
     *  0.9998569551145846E+00,0.9998660366782439E+00,
     *  0.9998762559494845E+00,0.9998918527924427E+00,
     *  0.9999185007547562E+00,0.9999558536849140E+00,
     *  0.9999899286967484E+00/
c 
        data whts39/
     1  0.5757494083510636E-02,0.2389274889339500E-01,
     2  0.5303929572222322E-01,0.8856642112178149E-01,
     3  0.1255286994435814E+00,0.1595345327077015E+00,
     4  0.1869822038572368E+00,0.2051819997504730E+00,
     5  0.2124563102143928E+00,0.2082099320383117E+00,
     6  0.1929540319461140E+00,0.1682753685488598E+00,
     7  0.1367504544702485E+00,0.1018052113473927E+00,
     8  0.6749720407945842E-01,0.3810277628509405E-01,
     9  0.1718085726147279E-01,0.5837771793779370E-02,
     *  0.1565934335252911E-02,0.4184770166682242E-03,
     *  0.1336828531365026E-03,0.5058444948701787E-04,
     *  0.2341770000466219E-04,0.1583019675239588E-04,
     *  0.1383311728195886E-04,0.1289579177517239E-04,
     *  0.1218228211087623E-04,0.1154751688750995E-04,
     *  0.1096017418497792E-04,0.1041044450120727E-04,
     *  0.9896387876687884E-05,0.9429628725563520E-05,
     *  0.9080011914624119E-05,0.9253858600163604E-05,
     *  0.1190595370462656E-04,0.2037229975139319E-04,
     *  0.3302590359453265E-04,0.3929828867347371E-04,
     *  0.2466761784132319E-04/
c 
c     Accuracy for the following 28 nodes is 0.74824E-07
c 
        data xs28/
     1  -.9963281626603377E+00,-.9725386172681839E+00,
     2  -.9138670084204048E+00,-.8132747678382736E+00,
     3  -.6706106576808419E+00,-.4910278917937606E+00,
     4  -.2836420059293417E+00,-.6027710641244244E-01,
     5  0.1657994297919880E+00,0.3811724670305296E+00,
     6  0.5735777879846626E+00,0.7331643243987818E+00,
     7  0.8536913775772343E+00,0.9336467361706713E+00,
     8  0.9772201819948978E+00,0.9945356814795852E+00,
     9  0.9988490165403534E+00,0.9995790853808850E+00,
     *  0.9997194714484159E+00,0.9997593445271709E+00,
     *  0.9997839029160636E+00,0.9998058879434815E+00,
     *  0.9998259159954668E+00,0.9998439261960848E+00,
     *  0.9998604089174878E+00,0.9998806100999196E+00,
     *  0.9999214405307657E+00,0.9999789011168279E+00/
c 
        data whts28/
     1  0.1073959702402644E-01,0.3939894841865883E-01,
     2  0.7909012322474880E-01,0.1220928287281714E+00,
     3  0.1623546709868317E+00,0.1952709262193308E+00,
     4  0.2175090977370362E+00,0.2269877858923984E+00,
     5  0.2229076628806538E+00,0.2057732650432674E+00,
     6  0.1773755470375106E+00,0.1407387451488585E+00,
     7  0.1000462653216561E+00,0.6055370887548284E-01,
     8  0.2831690756719004E-01,0.8668335331668530E-02,
     9  0.1584080208805724E-02,0.2723947327718108E-03,
     *  0.6410492619302724E-04,0.2724609128423976E-04,
     *  0.2300342700567222E-04,0.2101122916832621E-04,
     *  0.1901456038804858E-04,0.1706371312772466E-04,
     *  0.1645530675045620E-04,0.2744169084526995E-04,
     *  0.5450597249184941E-04,0.4928334186417706E-04/
c 
        if(n .eq. 39) call o1rarrcp(xs39,xs,n)
        if(n .eq. 39) call o1rarrcp(whts39,whts,n)
c 
        if(n .eq. 16) call o1rarrcp(xs16,xs,n)
        if(n .eq. 16) call o1rarrcp(whts16,whts,n)
c 
        if(n .eq. 28) call o1rarrcp(xs28,xs,n)
        if(n .eq. 28) call o1rarrcp(whts28,whts,n)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o1rnear9(n,xs,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),whts(1),
     1      xs37(37),whts37(37),xs30(30),whts30(30),
     2      xs17(17),whts17(17)
c 
c     Accuracy for the following 16 nodes is 0.36472E-04
c 
        data xs17/
     1  -.9879414536721423E+00,-.9327104112331326E+00,
     2  -.8266319750144406E+00,-.6708893554092040E+00,
     3  -.4731507078322269E+00,-.2454391215510237E+00,
     4  -.2559723912028078E-02,0.2394709254111413E+00,
     5  0.4648802424624505E+00,0.6596106675981011E+00,
     6  0.8128077699418557E+00,0.9182844598214599E+00,
     7  0.9763041437276696E+00,0.9966004319584969E+00,
     8  0.9994114314994895E+00,0.9996514118113751E+00,
     9  0.9998516627946117E+00/
c 
        data whts17/
     1  0.3152155003481836E-01,0.8016782745459142E-01,
     2  0.1317385392422669E+00,0.1784343470461179E+00,
     3  0.2150121920025710E+00,0.2379123496017204E+00,
     4  0.2451677325151350E+00,0.2362558524509146E+00,
     5  0.2121925373294671E+00,0.1754347982749418E+00,
     6  0.1298992215650774E+00,0.8105291040269537E-01,
     7  0.3652778922143909E-01,0.7736519965125721E-02,
     8  0.5114898124873552E-03,0.1600458492318300E-03,
     9  0.2676447567303698E-03/
c 
c     Accuracy for the following 37 nodes is 0.12644E-09
c 
        data xs37/
     1  -.9983817206528584E+00,-.9864907460041055E+00,
     2  -.9531140793676953E+00,-.8894260466987089E+00,
     3  -.7908902887770732E+00,-.6573731977448763E+00,
     4  -.4926505229004495E+00,-.3036837884472211E+00,
     5  -.9975613210819922E-01,0.1084976372430497E+00,
     6  0.3101260113609903E+00,0.4949080031052364E+00,
     7  0.6543543063567479E+00,0.7826050387082297E+00,
     8  0.8771654803613826E+00,0.9393907082741495E+00,
     9  0.9745060968277157E+00,0.9906798035007554E+00,
     *  0.9966098552848823E+00,0.9984869641175063E+00,
     *  0.9991047014392185E+00,0.9993404677972981E+00,
     *  0.9994436251251668E+00,0.9994967309310028E+00,
     *  0.9995361422401210E+00,0.9995708587273384E+00,
     *  0.9996034753286389E+00,0.9996318485100990E+00,
     *  0.9996570336106132E+00,0.9996825787460009E+00,
     *  0.9997092147926014E+00,0.9997349784271357E+00,
     *  0.9997644384795139E+00,0.9998052349247823E+00,
     *  0.9998624379011878E+00,0.9999300445037895E+00,
     *  0.9999847466057044E+00/
c 
        data whts37/
     1  0.4898341152654785E-02,0.2085216645249485E-01,
     2  0.4739788763999245E-01,0.8073358412312873E-01,
     3  0.1163513029397802E+00,0.1500413904612085E+00,
     4  0.1782292359177023E+00,0.1981407807465310E+00,
     5  0.2079246853443265E+00,0.2067400061892669E+00,
     6  0.1947964213150425E+00,0.1733375846855034E+00,
     7  0.1445659339325843E+00,0.1115112362355294E+00,
     8  0.7783401382833588E-01,0.4750091772604132E-01,
     9  0.2413247219518208E-01,0.9724707760355156E-02,
     *  0.3183120151095869E-02,0.1002128737137244E-02,
     *  0.3582054582768432E-03,0.1487728248049760E-03,
     *  0.6959275426944986E-04,0.4285254814313115E-04,
     *  0.3673263714606619E-04,0.3344185038435590E-04,
     *  0.3103582848977658E-04,0.2610324895872845E-04,
     *  0.2482020953919196E-04,0.2649557477049054E-04,
     *  0.2615810513815243E-04,0.2630343034296744E-04,
     *  0.3391439747095840E-04,0.4870168507422221E-04,
     *  0.6485202397148592E-04,0.6628082185691204E-04,
     *  0.3781927574159915E-04/
c 
c     Accuracy for the following 30 nodes is better than 0.22344E-07
c 
        data xs30/
     1  -.9967418229032425E+00,-.9751971912238029E+00,
     2  -.9210258618576445E+00,-.8267242908379990E+00,
     3  -.6913228525103319E+00,-.5190566827432828E+00,
     4  -.3181517295677918E+00,-.9965599334360231E-01,
     5  0.1237672130869780E+00,0.3390918810432888E+00,
     6  0.5341910208124501E+00,0.6990569923480643E+00,
     7  0.8269383033141506E+00,0.9153713749550143E+00,
     8  0.9670699385198218E+00,0.9903371828243639E+00,
     9  0.9974687158758602E+00,0.9990251433207766E+00,
     *  0.9993816998794087E+00,0.9994865075817809E+00,
     *  0.9995359068982872E+00,0.9995770420610970E+00,
     *  0.9996163822816837E+00,0.9996537588669969E+00,
     *  0.9996871369353666E+00,0.9997152388329824E+00,
     *  0.9997421957517896E+00,0.9997864222907074E+00,
     *  0.9998697795300431E+00,0.9999674931555911E+00/
c 
        data whts30/
     1  0.9585916093982029E-02,0.3599921061286702E-01,
     2  0.7357567501437021E-01,0.1151660930511292E+00,
     3  0.1549110549442312E+00,0.1882310141336129E+00,
     4  0.2117174460859669E+00,0.2231424599801997E+00,
     5  0.2215117094154027E+00,0.2070959877300238E+00,
     6  0.1814153570132682E+00,0.1471741252682770E+00,
     7  0.1081626287106606E+00,0.6914589791749690E-01,
     8  0.3565275767940124E-01,0.1303639841544550E-01,
     9  0.3097178675561558E-02,0.6568618698363948E-03,
     *  0.1739942344003787E-03,0.6236648179073294E-04,
     *  0.4287577092464030E-04,0.4002655094238928E-04,
     *  0.3860166859514679E-04,0.3577115204740642E-04,
     *  0.3070006286538195E-04,0.2602475577139285E-04,
     *  0.3100972243300507E-04,0.6207527558115843E-04,
     *  0.1012058613980530E-03,0.7767858654782477E-04/
  
c 
        if(n .eq. 37) call o1rarrcp(xs37,xs,n)
        if(n .eq. 37) call o1rarrcp(whts37,whts,n)
c 
        if(n .eq. 17) call o1rarrcp(xs17,xs,n)
        if(n .eq. 17) call o1rarrcp(whts17,whts,n)
c 
        if(n .eq. 30) call o1rarrcp(xs30,xs,n)
        if(n .eq. 30) call o1rarrcp(whts30,whts,n)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o1rnear10(n,xs,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),whts(1),
     1      xs47(47),whts47(47),xs30(30),whts30(30),
     2      xs17(17),whts17(17)
c 
c     Accuracy for the following 17 nodes is 0.41687E-04
c 
        data xs17/
     1  -.9866579819098580E+00,-.9271625344676936E+00,
     2  -.8156524229394941E+00,-.6542530369355177E+00,
     3  -.4512794139016348E+00,-.2192680943098000E+00,
     4  0.2652815094257366E-01,0.2697283179803718E+00,
     5  0.4943129231202656E+00,0.6860912169073008E+00,
     6  0.8342780610052585E+00,0.9330659256565144E+00,
     7  0.9835802043419917E+00,0.9975779184912870E+00,
     8  0.9990590056776190E+00,0.9993684809472225E+00,
     9  0.9997676615368145E+00/
c 
        data whts17/
     1  0.3462165406267913E-01,0.8523368724273411E-01,
     2  0.1373995664882517E+00,0.1839753311029881E+00,
     3  0.2198539779212463E+00,0.2415892684271622E+00,
     4  0.2472572193927403E+00,0.2364560335685041E+00,
     5  0.2103134518911102E+00,0.1714183307995363E+00,
     6  0.1239538420555245E+00,0.7376734709074366E-01,
     7  0.2900238009849470E-01,0.3974479994807952E-02,
     8  0.4304633165012101E-03,0.2949802353620657E-03,
     9  0.4823387040390013E-03/
  
c 
c     Accuracy for the following 38 nodes is 0.17546E-09
c 
        data xs47/
     1  -.9986665735405927E+00,-.9886927227589883E+00,
     2  -.9600822983057403E+00,-.9043335576571234E+00,
     3  -.8164299470699833E+00,-.6952471482516047E+00,
     4  -.5432987133363408E+00,-.3661946210310969E+00,
     5  -.1719286215431847E+00,0.2995940108041184E-01,
     6  0.2293046089049755E+00,0.4162721269526408E+00,
     7  0.5822897060387406E+00,0.7209005987999693E+00,
     8  0.8284779938621194E+00,0.9047274679225388E+00,
     9  0.9528429713189490E+00,0.9790388162350232E+00,
     *  0.9910685389364806E+00,0.9958431019219398E+00,
     *  0.9976504543615556E+00,0.9983779735285394E+00,
     *  0.9987015741385601E+00,0.9988590956437622E+00,
     *  0.9989436969032577E+00,0.9989977970392403E+00,
     *  0.9990411212700634E+00,0.9990809041197573E+00,
     *  0.9991194520887743E+00,0.9991577000628791E+00,
     *  0.9991960662640487E+00,0.9992345244118022E+00,
     *  0.9992726394214200E+00,0.9993096995403143E+00,
     *  0.9993448959705100E+00,0.9993775218626078E+00,
     *  0.9994071898110004E+00,0.9994340447430201E+00,
     *  0.9994589734990791E+00,0.9994841400822359E+00,
     *  0.9995143932152688E+00,0.9995578133409708E+00,
     *  0.9996231208276514E+00,0.9997148902367305E+00,
     *  0.9998256511090706E+00,0.9999290844550444E+00,
     *  0.9999890396643043E+00/
c 
        data whts47/
     1  0.4053860562010578E-02,0.1764508962573425E-01,
     2  0.4102213326725003E-01,0.7131607217726220E-01,
     3  0.1046735527808954E+00,0.1372638887210820E+00,
     4  0.1656860168325318E+00,0.1871758198087664E+00,
     5  0.1997507253274121E+00,0.2023151294937307E+00,
     6  0.1947204581906579E+00,0.1777715195366162E+00,
     7  0.1531754482790443E+00,0.1234332125373908E+00,
     8  0.9166927353939007E-01,0.6136867032230332E-01,
     9  0.3592237618601443E-01,0.1780532820536939E-01,
     *  0.7440932559535669E-02,0.2814850577558888E-02,
     *  0.1090645878635309E-02,0.4647657366340785E-03,
     *  0.2179490292813609E-03,0.1113149135236898E-03,
     *  0.6458793633020822E-04,0.4676409853349798E-04,
     *  0.4094867799536366E-04,0.3894980932887033E-04,
     *  0.3829212185462719E-04,0.3827387884354178E-04,
     *  0.3845367053220287E-04,0.3838771772002313E-04,
     *  0.3772062515317964E-04,0.3626152302581010E-04,
     *  0.3400985003402677E-04,0.3117369262126278E-04,
     *  0.2818344944202757E-04,0.2566769926642054E-04,
     *  0.2452419996367963E-04,0.2659548318828946E-04,
     *  0.3531235409917952E-04,0.5305107333107685E-04,
     *  0.7846658569986502E-04,0.1039778262628404E-03,
     *  0.1130600626503551E-03,0.8706186778495462E-04,
     *  0.3124168157049121E-04/
c 
c     Accuracy for the following 30 nodes is 0.22344E-07
c 
        data xs30/
     1  -.9969916639295240E+00,-.9768082736328868E+00,
     2  -.9253386682961912E+00,-.8347378403432630E+00,
     3  -.7035070841979412E+00,-.5353396259603103E+00,
     4  -.3379743667567272E+00,-.1220602287806630E+00,
     5  0.1000286096409917E+00,0.3154263354689280E+00,
     6  0.5120203944234156E+00,0.6796595290618483E+00,
     7  0.8112858601221749E+00,0.9039719573966931E+00,
     8  0.9598281365462769E+00,0.9864770176201352E+00,
     9  0.9956859801574658E+00,0.9981130557180888E+00,
     *  0.9987632176639327E+00,0.9989732967470992E+00,
     *  0.9990756685692255E+00,0.9991466249543235E+00,
     *  0.9992176375079624E+00,0.9992880502988484E+00,
     *  0.9993563010140587E+00,0.9994105544230132E+00,
     *  0.9994831227630896E+00,0.9995926664822059E+00,
     *  0.9997614213220225E+00,0.9999418892803020E+00/
c 
        data whts30/
     1  0.8887177636259062E-02,0.3394192044368893E-01,
     2  0.7029308297588872E-01,0.1111379073623436E+00,
     3  0.1506891824299419E+00,0.1843371455599010E+00,
     4  0.2085987864033999E+00,0.2211464537732629E+00,
     5  0.2208660194851013E+00,0.2078904243407638E+00,
     6  0.1835869080580791E+00,0.1504977434439068E+00,
     7  0.1122514916890019E+00,0.7346175544770192E-01,
     8  0.3953539396814913E-01,0.1581969707938032E-01,
     9  0.4482075018097530E-02,0.1140041323851704E-02,
     *  0.3384569556868720E-03,0.1290388521268538E-03,
     *  0.8476410047991512E-04,0.6414317692781626E-04,
     *  0.7488244448006013E-04,0.6819964063785624E-04,
     *  0.6244635519339861E-04,0.5647976310478788E-04,
     *  0.8813839109490546E-04,0.1371389689465429E-03,
     *  0.1934176379805408E-03,0.1397204953005944E-03/
c 
        if(n .eq. 47) call o1rarrcp(xs47,xs,n)
        if(n .eq. 47) call o1rarrcp(whts47,whts,n)
c 
        if(n .eq. 17) call o1rarrcp(xs17,xs,n)
        if(n .eq. 17) call o1rarrcp(whts17,whts,n)
c 
        if(n .eq. 30) call o1rarrcp(xs30,xs,n)
        if(n .eq. 30) call o1rarrcp(whts30,whts,n)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o1rnear11(n,xs,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),whts(1),
     1      xs42(42),whts42(42),xs30(30),whts30(30),
     2      xs17(17),whts17(17)
c 
c     Accuracy for the following 17 nodes is 0.58787E-04
c 
        data xs17/
     1  -.9861732277641303E+00,-.9251115589637638E+00,
     2  -.8117860128147517E+00,-.6487532453678126E+00,
     3  -.4445632923483561E+00,-.2118767133421081E+00,
     4  0.3401961814011689E-01,0.2767533819231016E+00,
     5  0.5003488390632891E+00,0.6906867308021192E+00,
     6  0.8371070286118363E+00,0.9340340553061168E+00,
     7  0.9830510732397473E+00,0.9966107040112349E+00,
     8  0.9982122247917751E+00,0.9987647961774930E+00,
     9  0.9994933135457043E+00/
c 
        data whts17/
     1  0.3578894557060544E-01,0.8703111358971166E-01,
     2  0.1391709883413029E+00,0.1854311907808789E+00,
     3  0.2208101873822896E+00,0.2419799251537226E+00,
     4  0.2470676243794918E+00,0.2357189342862663E+00,
     5  0.2090843352327998E+00,0.1697864877284461E+00,
     6  0.1220905535074594E+00,0.7197775301325024E-01,
     7  0.2795151075188587E-01,0.4000943707883949E-02,
     8  0.6376955040309462E-03,0.5227489260769685E-03,
     9  0.9669828424786223E-03/
c 
c     Accuracy for the following 42 nodes is 0.25219E-09
c 
        data xs42/
     1  -.9983178371252187E+00,-.9860282902604156E+00,
     2  -.9517583038931584E+00,-.8867473662607669E+00,
     3  -.7866588363707667E+00,-.6516032113241570E+00,
     4  -.4855934026764842E+00,-.2957917292606577E+00,
     5  -.9162684339294003E-01,0.1161890405132626E+00,
     6  0.3167036902044966E+00,0.4997774933824115E+00,
     7  0.6570876313486067E+00,0.7830263537105215E+00,
     8  0.8754315569586888E+00,0.9360457350863420E+00,
     9  0.9704433044979510E+00,0.9868675705020802E+00,
     *  0.9935791790833459E+00,0.9961588243122641E+00,
     *  0.9971993655101383E+00,0.9976589418531587E+00,
     *  0.9978819914315211E+00,0.9980140398702788E+00,
     *  0.9981241274098338E+00,0.9982307571032885E+00,
     *  0.9983261760505613E+00,0.9984079541157225E+00,
     *  0.9984864983347423E+00,0.9985623867229321E+00,
     *  0.9986356996728809E+00,0.9987049428055035E+00,
     *  0.9987696922739624E+00,0.9988300692820160E+00,
     *  0.9988910452474686E+00,0.9989587612361019E+00,
     *  0.9990430391259884E+00,0.9991670365050362E+00,
     *  0.9993516364969539E+00,0.9995906796404954E+00,
     *  0.9998283147879052E+00,0.9999726845191901E+00/
c 
        data whts42/
     1  0.5084272403945542E-02,0.2149190932355263E-01,
     2  0.4853183318310350E-01,0.8221065550365651E-01,
     3  0.1179382898180835E+00,0.1514914788715276E+00,
     4  0.1793196428340043E+00,0.1986957892018996E+00,
     5  0.2078310947671437E+00,0.2059576172510495E+00,
     6  0.1933648741799782E+00,0.1713804887484896E+00,
     7  0.1422927804033902E+00,0.1092161341767479E+00,
     8  0.7588452079485464E-01,0.4629043230038900E-01,
     9  0.2392588398590411E-01,0.1035364909722021E-01,
     *  0.4001129612684838E-02,0.1561478807557253E-02,
     *  0.6627194074690085E-03,0.3078344087396083E-03,
     *  0.1608481673508624E-03,0.1145862980352371E-03,
     *  0.1082403372490350E-03,0.1032008111835435E-03,
     *  0.8682960266798891E-04,0.7938234245667659E-04,
     *  0.7721615656758228E-04,0.7482547033971261E-04,
     *  0.7137263467735063E-04,0.6715183704041952E-04,
     *  0.6221054839699534E-04,0.5947555891246305E-04,
     *  0.6352546648442018E-04,0.7318451562020975E-04,
     *  0.9963148261188941E-04,0.1522268318682163E-03,
     *  0.2164102777093983E-03,0.2523287352667311E-03,
     *  0.2057681223258837E-03,0.7707565368177447E-04/
c 
c     Accuracy for the following 30 nodes is 0.68355E-07
c 
        data xs30/
     1  -.9966942385751845E+00,-.9748894242792693E+00,
     2  -.9202083890429905E+00,-.8252377498460105E+00,
     3  -.6891463979331552E+00,-.5163111455979581E+00,
     4  -.3150798642209822E+00,-.9659012698273373E-01,
     5  0.1264447228875387E+00,0.3409938843634950E+00,
     6  0.5349765021024302E+00,0.6984902599353531E+00,
     7  0.8249591751032208E+00,0.9121777987992437E+00,
     8  0.9632013261287480E+00,0.9866663066693810E+00,
     9  0.9947249542338910E+00,0.9970545313044312E+00,
     *  0.9977537969945296E+00,0.9980065643891700E+00,
     *  0.9981668763801279E+00,0.9983175608740804E+00,
     *  0.9984689554510655E+00,0.9986171126906236E+00,
     *  0.9987524568704231E+00,0.9988694649176624E+00,
     *  0.9989814805696751E+00,0.9991640770184866E+00,
     *  0.9994997272091999E+00,0.9998772589753491E+00/
c 
        data whts30/
     1  0.9719296606387488E-02,0.3639168035012035E-01,
     2  0.7418618353993849E-01,0.1158699383482946E+00,
     3  0.1555630663799083E+00,0.1886966474952290E+00,
     4  0.2118894283101709E+00,0.2229500899171389E+00,
     5  0.2209265964444755E+00,0.2061376683954076E+00,
     6  0.1801580001263742E+00,0.1457560834646163E+00,
     7  0.1067987637377143E+00,0.6813865813221805E-01,
     8  0.3537216490868560E-01,0.1370555769728701E-01,
     9  0.4088567525039503E-02,0.1182716138436121E-02,
     *  0.3807382827436099E-03,0.1786015121816806E-03,
     *  0.1517939458610513E-03,0.1508619807638579E-03,
     *  0.1511742819908609E-03,0.1433738621699895E-03,
     *  0.1262197490150514E-03,0.1091634823915048E-03,
     *  0.1278971595983522E-03,0.2552172470423266E-03,
     *  0.3990639525212012E-03,0.2947287923611244E-03/
c 
        if(n .eq. 42) call o1rarrcp(xs42,xs,n)
        if(n .eq. 42) call o1rarrcp(whts42,whts,n)
c 
        if(n .eq. 17) call o1rarrcp(xs17,xs,n)
        if(n .eq. 17) call o1rarrcp(whts17,whts,n)
c 
        if(n .eq. 30) call o1rarrcp(xs30,xs,n)
        if(n .eq. 30) call o1rarrcp(whts30,whts,n)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o1rnear12(n,xs,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),whts(1),
     1      xs42(42),whts42(42),xs31(31),whts31(31),
     2      xs18(18),whts18(18)
c 
c     Accuracy for the following 18 nodes is 0.46554E-04
c 
        data xs18/
     1  -.9870877523348666E+00,-.9290838751045076E+00,
     2  -.8195909119101820E+00,-.6604179125607113E+00,
     3  -.4596439268578998E+00,-.2296051288697747E+00,
     4  0.1462847975980908E-01,0.2568384272997792E+00,
     5  0.4811274602567563E+00,0.6733980453020409E+00,
     6  0.8229079024493053E+00,0.9237881285853686E+00,
     7  0.9768905457470063E+00,0.9935541270357258E+00,
     8  0.9962234331401463E+00,0.9971067560452385E+00,
     9  0.9978395451998971E+00,0.9991902871790486E+00/
c 
        data whts18/
     1  0.3357220928272799E-01,0.8340927796081961E-01,
     2  0.1352365819060386E+00,0.1817262957614091E+00,
     3  0.2177379256758424E+00,0.2397932456100494E+00,
     4  0.2459592975877400E+00,0.2357977377033408E+00,
     5  0.2103987597613524E+00,0.1723259924789383E+00,
     6  0.1256844444766963E+00,0.7617366311682357E-01,
     7  0.3174207448113480E-01,0.6054459236363827E-02,
     8  0.1105774751329894E-02,0.7658713869627623E-03,
     9  0.8589229862009088E-03,0.1667455336617089E-02/
c 
c     Accuracy for the following 42 nodes is 0.27161E-09
c 
        data xs42/
     1  -.9982527005402485E+00,-.9855512830223195E+00,
     2  -.9503352973650422E+00,-.8838727479538691E+00,
     3  -.7819967038021015E+00,-.6450490079338969E+00,
     4  -.4772877013577660E+00,-.2861037300639935E+00,
     5  -.8111585805505062E-01,0.1268290646742598E+00,
     6  0.3267128149560881E+00,0.5084090518924771E+00,
     7  0.6636973088656528E+00,0.7871722935444025E+00,
     8  0.8769865090553193E+00,0.9353077608300036E+00,
     9  0.9681666280123706E+00,0.9840228459823845E+00,
     *  0.9908422001028885E+00,0.9937075011527119E+00,
     *  0.9949737749950891E+00,0.9955761809593264E+00,
     *  0.9958961314918662E+00,0.9961168828599130E+00,
     *  0.9963106079011510E+00,0.9964919174967066E+00,
     *  0.9966615018633402E+00,0.9968200361296358E+00,
     *  0.9969744255248676E+00,0.9971256290291473E+00,
     *  0.9972717648736302E+00,0.9974104527171985E+00,
     *  0.9975401226257765E+00,0.9976619249569858E+00,
     *  0.9977840803294153E+00,0.9979191201862558E+00,
     *  0.9980884689386897E+00,0.9983396697783542E+00,
     *  0.9987144898339573E+00,0.9991961334289254E+00,
     *  0.9996657784671956E+00,0.9999461629157369E+00/
c 
        data whts42/
     1  0.5274383761023150E-02,0.2215878191096259E-01,
     2  0.4975042311713126E-01,0.8386567211248966E-01,
     3  0.1198189924200746E+00,0.1533535781975158E+00,
     4  0.1809220738571488E+00,0.1998259882848841E+00,
     5  0.2083242119343291E+00,0.2057114062940149E+00,
     6  0.1923513398258586E+00,0.1696557684811457E+00,
     7  0.1400072604610833E+00,0.1066265078849012E+00,
     8  0.7336400726570919E-01,0.4430572465041661E-01,
     9  0.2287626805843534E-01,0.1021143625860855E-01,
     *  0.4265470120528885E-02,0.1831704857101517E-02,
     *  0.8430123277027857E-03,0.4197899107726816E-03,
     *  0.2501374840106424E-03,0.2022876213760088E-03,
     *  0.1867806871329980E-03,0.1759503165184295E-03,
     *  0.1629947834306961E-03,0.1556897457341727E-03,
     *  0.1530447756567300E-03,0.1490541395639488E-03,
     *  0.1427903034064191E-03,0.1343022058925830E-03,
     *  0.1251382686652915E-03,0.1198620868379369E-03,
     *  0.1266368751896641E-03,0.1463312486572503E-03,
     *  0.2010859878662307E-03,0.3091150051125861E-03,
     *  0.4386728478805177E-03,0.5043735743893055E-03,
     *  0.4015883830639555E-03,0.1503617086758107E-03/
c 
c     Accuracy for the following 31 nodes is 0.38201E-07
c 
        data xs31/
     1  -.9966979437193559E+00,-.9749318251650521E+00,
     2  -.9203811460449804E+00,-.8256840108395085E+00,
     3  -.6900411510935957E+00,-.5178371031447744E+00,
     4  -.3174020950059990E+00,-.9983172936964346E-01,
     5  0.1222240790206754E+00,0.3358158719164255E+00,
     6  0.5289591903556270E+00,0.6918616376863866E+00,
     7  0.8180700696933074E+00,0.9055100837231620E+00,
     8  0.9573415684592621E+00,0.9821156609577943E+00,
     9  0.9914458304134657E+00,0.9945365605478098E+00,
     *  0.9956011179451729E+00,0.9960458128289209E+00,
     *  0.9963495080766223E+00,0.9966410751104569E+00,
     *  0.9969395332558666E+00,0.9972337243449985E+00,
     *  0.9975009837986772E+00,0.9977268328351137E+00,
     *  0.9979327159655649E+00,0.9982232195716940E+00,
     *  0.9987257291960874E+00,0.9993932145439282E+00,
     *  0.9998926264201247E+00/
c 
        data whts31/
     1  0.9706535146345072E-02,0.3631696553830411E-01,
     2  0.7399142168838058E-01,0.1155123194547887E+00,
     3  0.1550223935467855E+00,0.1879778698506727E+00,
     4  0.2110227058477561E+00,0.2219886957188950E+00,
     5  0.2199434985797304E+00,0.2052221169458407E+00,
     6  0.1794132746586226E+00,0.1452986986375486E+00,
     7  0.1067573983037211E+00,0.6864244493227589E-01,
     8  0.3647740022370571E-01,0.1511936771259436E-01,
     9  0.5124520684319562E-02,0.1702611472167909E-02,
     *  0.6300846280083208E-03,0.3334983855848152E-03,
     *  0.2904831150429054E-03,0.2950903206443146E-03,
     *  0.2997218996064184E-03,0.2845535240020294E-03,
     *  0.2471550735675485E-03,0.2071045350241164E-03,
     *  0.2228306536207666E-03,0.3830727652151688E-03,
     *  0.6173873099560279E-03,0.6556123968532743E-03,
     *  0.2931417615889141E-03/
c 
        if(n .eq. 42) call o1rarrcp(xs42,xs,n)
        if(n .eq. 42) call o1rarrcp(whts42,whts,n)
c 
        if(n .eq. 18) call o1rarrcp(xs18,xs,n)
        if(n .eq. 18) call o1rarrcp(whts18,whts,n)
c 
        if(n .eq. 31) call o1rarrcp(xs31,xs,n)
        if(n .eq. 31) call o1rarrcp(whts31,whts,n)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o1rnear13(n,xs,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),whts(1),
     1      xs50(50),whts50(50),xs31(31),whts31(31),
     2      xs18(18),whts18(18)
c 
c     Accuracy for the following 18 nodes is  0.64907E-04
c 
        data xs18/
     1  -.9860285324398091E+00,-.9245305429735583E+00,
     2  -.8108008142950511E+00,-.6475771446753340E+00,
     3  -.4435029050081921E+00,-.2112684979457761E+00,
     4  0.3386807816414938E-01,0.2755879936899794E+00,
     5  0.4979704509677876E+00,0.6869531256011041E+00,
     6  0.8319468295840187E+00,0.9275007219317963E+00,
     7  0.9756930450515684E+00,0.9903238322153880E+00,
     8  0.9931731096180161E+00,0.9946661845835074E+00,
     9  0.9961588332944060E+00,0.9986927985372456E+00/
c 
        data whts18/
     1  0.3613330012995496E-01,0.8749517596246079E-01,
     2  0.1394887003296879E+00,0.1854801401541192E+00,
     3  0.2205255971448036E+00,0.2413676481073385E+00,
     4  0.2461717776108003E+00,0.2345973603209431E+00,
     5  0.2077906578888661E+00,0.1683814647172754E+00,
     6  0.1206611263620315E+00,0.7071853041054438E-01,
     7  0.2789978765367582E-01,0.5605727615426861E-02,
     8  0.1642724765234665E-02,0.1402072840591896E-02,
     9  0.1872198728064583E-02,0.2761859995277817E-02/
c 
c 
c     Accuracy for the following 50 nodes is 0.12135E-10
c 
        data xs50/
     1  -.9987582748647483E+00,-.9894097470945746E+00,
     2  -.9623801040261449E+00,-.9093172770864792E+00,
     3  -.8251074963955537E+00,-.7083786790972570E+00,
     4  -.5613182411874900E+00,-.3891775486941846E+00,
     5  -.1995908291312456E+00,-.1761793543295753E-02,
     6  0.1944380438553057E+00,0.3794090346914569E+00,
     7  0.5447450055703078E+00,0.6840785389648856E+00,
     8  0.7937881475377264E+00,0.8734842415584505E+00,
     9  0.9261165529996396E+00,0.9574022019974731E+00,
     *  0.9742845325607215E+00,0.9828663579109151E+00,
     *  0.9871892095426980E+00,0.9894276085132612E+00,
     *  0.9906402636708214E+00,0.9913394420424422E+00,
     *  0.9917867276786334E+00,0.9921232870936709E+00,
     *  0.9924180543147325E+00,0.9926992777840569E+00,
     *  0.9929782657718516E+00,0.9932579480105201E+00,
     *  0.9935342687515263E+00,0.9937990858852421E+00,
     *  0.9940487842076209E+00,0.9942863292015496E+00,
     *  0.9945141353498546E+00,0.9947348099778883E+00,
     *  0.9949511112416272E+00,0.9951635278376358E+00,
     *  0.9953658247591086E+00,0.9955527669903685E+00,
     *  0.9957282078956013E+00,0.9959084175307726E+00,
     *  0.9961268900139355E+00,0.9964333837652219E+00,
     *  0.9968798531074499E+00,0.9974979385365692E+00,
     *  0.9982602683261244E+00,0.9990435189634561E+00,
     *  0.9996526728255870E+00,0.9999500168126687E+00/
c 
        data whts50/
     1  0.3781094933127256E-02,0.1659199828388919E-01,
     2  0.3888985862331132E-01,0.6809448895150340E-01,
     3  0.1005499493314537E+00,0.1325346177029622E+00,
     4  0.1606995730852133E+00,0.1822957836843963E+00,
     5  0.1953285548062040E+00,0.1986684675763820E+00,
     6  0.1921157667110085E+00,0.1764109564543449E+00,
     7  0.1531877848109140E+00,0.1248672666150001E+00,
     8  0.9448450132821515E-01,0.6540781632691536E-01,
     9  0.4083603793128820E-01,0.2293128511687736E-01,
     *  0.1187261909354237E-01,0.5958238883011681E-02,
     *  0.3035539662879050E-02,0.1608270893672365E-02,
     *  0.8978268234776418E-03,0.5420982995207413E-03,
     *  0.3754232487393351E-03,0.3086241925319587E-03,
     *  0.2852032720458576E-03,0.2790139843729888E-03,
     *  0.2794459775681628E-03,0.2792125764681348E-03,
     *  0.2718948752675943E-03,0.2570596292238492E-03,
     *  0.2430676541665311E-03,0.2323127849065295E-03,
     *  0.2237937993554765E-03,0.2179933304869622E-03,
     *  0.2148692926253392E-03,0.2087166059624126E-03,
     *  0.1948694716934430E-03,0.1796274112771004E-03,
     *  0.1738679547522492E-03,0.1921058385167162E-03,
     *  0.2534201760010145E-03,0.3685424664728108E-03,
     *  0.5303562160829928E-03,0.7020271919020871E-03,
     *  0.8018402889754004E-03,0.7298124581703385E-03,
     *  0.4624942412293187E-03,0.1440091293385179E-03/
c 
c     Accuracy for the following 31 nodes is 0.38201E-07
c 
        data xs31/
     1  -.9966074331179129E+00,-.9743524462883377E+00,
     2  -.9188370313580936E+00,-.8228308345478204E+00,
     3  -.6857473724893600E+00,-.5122006083439166E+00,
     4  -.3107337281532927E+00,-.9261830873946719E-01,
     5  0.1293697332220471E+00,0.3422145196667308E+00,
     6  0.5339348472944696E+00,0.6948259023395209E+00,
     7  0.8186306622845954E+00,0.9036182706610636E+00,
     8  0.9534839178893166E+00,0.9773724105705857E+00,
     9  0.9868312031245423E+00,0.9902839154215810E+00,
     *  0.9915892257725500E+00,0.9922317562217866E+00,
     *  0.9927569883010301E+00,0.9933024482112676E+00,
     *  0.9938850488104145E+00,0.9944693265091724E+00,
     *  0.9950028563576900E+00,0.9954545076708514E+00,
     *  0.9958679290362145E+00,0.9964593695078670E+00,
     *  0.9974917917215491E+00,0.9988393657545604E+00,
     *  0.9998027252227544E+00/
c 
        data whts31/
     1  0.9958853738576549E-02,0.3705374603379553E-01,
     2  0.7516002596379130E-01,0.1169254033503688E+00,
     3  0.1564517586422479E+00,0.1891981070709448E+00,
     4  0.2118366750633916E+00,0.2222436768212320E+00,
     5  0.2195417478709255E+00,0.2041301619201255E+00,
     6  0.1776734508027206E+00,0.1430477891233167E+00,
     7  0.1042572683727020E+00,0.6632744263173433E-01,
     8  0.3497377282128245E-01,0.1484885959490648E-01,
     9  0.5489069017832722E-02,0.2009630128998229E-02,
     *  0.8265431101006533E-03,0.5422101404307216E-03,
     *  0.5263978585246255E-03,0.5667276896383815E-03,
     *  0.5919710685490235E-03,0.5672274786145575E-03,
     *  0.4938964779584956E-03,0.4145934516183537E-03,
     *  0.4497214206967682E-03,0.7851868126456974E-03,
     *  0.1265275911639298E-02,0.1296311352343354E-02,
     *  0.5464843711721665E-03/
c 
        if(n .eq. 50) call o1rarrcp(xs50,xs,n)
        if(n .eq. 50) call o1rarrcp(whts50,whts,n)
c 
        if(n .eq. 18) call o1rarrcp(xs18,xs,n)
        if(n .eq. 18) call o1rarrcp(whts18,whts,n)
c 
        if(n .eq. 31) call o1rarrcp(xs31,xs,n)
        if(n .eq. 31) call o1rarrcp(whts31,whts,n)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o1rnear14(n,xs,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),whts(1),
     1      xs46(46),whts46(46),xs32(32),whts32(32),
     2      xs19(19),whts19(19)
c 
c     Accuracy for the following 19 nodes is 0.42910E-04
c 
        data xs19/
     1  -.9870103202852024E+00,-.9288577680584615E+00,
     2  -.8194453466525373E+00,-.6607104563652084E+00,
     3  -.4607658506104918E+00,-.2319154059688000E+00,
     4  0.1084813829882125E-01,0.2513941494983118E+00,
     5  0.4739082788436384E+00,0.6643764328470896E+00,
     6  0.8121519819987928E+00,0.9115035927824444E+00,
     7  0.9637763602906064E+00,0.9814676650572123E+00,
     8  0.9857962636956237E+00,0.9882909327504759E+00,
     9  0.9904548265826771E+00,0.9931459574208845E+00,
     *  0.9979691547262053E+00/
c 
        data whts19/
     1  0.3374265618800822E-01,0.8347821363000787E-01,
     2  0.1349883372092302E+00,0.1810926842460889E+00,
     3  0.2167185800663124E+00,0.2384496967911043E+00,
     4  0.2443726953923462E+00,0.2340653921612032E+00,
     5  0.2085955239557015E+00,0.1705401964426289E+00,
     6  0.1240217372304422E+00,0.7483957490859475E-01,
     7  0.3171042641479768E-01,0.7738532819826905E-02,
     8  0.2746192805133395E-02,0.2330721156025912E-02,
     9  0.2030886024621981E-02,0.3939497416089163E-02,
     *  0.4592092880992467E-02/
c 
c     Accuracy for the following 46 nodes is 0.95466E-10
c 
        data xs46/
     1  -.9984221277377099E+00,-.9868236548579349E+00,
     2  -.9542570248480122E+00,-.8921057135784192E+00,
     3  -.7959526609577885E+00,-.6656895700689063E+00,
     4  -.5050242378505509E+00,-.3207607974501910E+00,
     5  -.1219440895128696E+00,0.8109697488290379E-01,
     6  0.2777678262788976E+00,0.4582379196018616E+00,
     7  0.6144255217119923E+00,0.7408803716121455E+00,
     8  0.8354966622964983E+00,0.8999249741367401E+00,
     9  0.9393557670856505E+00,0.9611440933403357E+00,
     *  0.9723867803508359E+00,0.9780811025160976E+00,
     *  0.9810322678919270E+00,0.9826525969130306E+00,
     *  0.9836490520082082E+00,0.9843932667512305E+00,
     *  0.9850563996747159E+00,0.9856907068306022E+00,
     *  0.9862973785468352E+00,0.9868729267600262E+00,
     *  0.9874181996908205E+00,0.9879320075036118E+00,
     *  0.9884168760277801E+00,0.9888780173124077E+00,
     *  0.9893215065769903E+00,0.9897533741278655E+00,
     *  0.9901789840496571E+00,0.9906006111915240E+00,
     *  0.9910200912741709E+00,0.9914426340699785E+00,
     *  0.9918982771643272E+00,0.9924817964094776E+00,
     *  0.9933281110483823E+00,0.9945509617869686E+00,
     *  0.9961428458959591E+00,0.9978520758625080E+00,
     *  0.9992161043632724E+00,0.9998873254341518E+00/
c 
        data whts46/
     1  0.4776523237463896E-02,0.2034256074598809E-01,
     2  0.4625218632828775E-01,0.7878593206919150E-01,
     3  0.1135293352641394E+00,0.1463667420022203E+00,
     4  0.1738155732017661E+00,0.1931877972601113E+00,
     5  0.2027070367554630E+00,0.2015920328913403E+00,
     6  0.1900946980652543E+00,0.1694836779886065E+00,
     7  0.1419710096540675E+00,0.1105801492013950E+00,
     8  0.7892957812762039E-01,0.5081448895755294E-01,
     9  0.2931834525197834E-01,0.1547570893151718E-01,
     *  0.7840553303124751E-02,0.3999657716531740E-02,
     *  0.2126564422795147E-02,0.1225787733975778E-02,
     *  0.8275880312277462E-03,0.6878639175999714E-03,
     *  0.6460389389061850E-03,0.6219199342477649E-03,
     *  0.5907406951107934E-03,0.5607263052187929E-03,
     *  0.5294325197594654E-03,0.4986413800740513E-03,
     *  0.4720063883487545E-03,0.4513005536194161E-03,
     *  0.4366838378858526E-03,0.4279998150081952E-03,
     *  0.4235681821550365E-03,0.4199191271178129E-03,
     *  0.4196759730576377E-03,0.4293529186394555E-03,
     *  0.4983846354922944E-03,0.6920634384955123E-03,
     *  0.1021168669151439E-02,0.1425836297599754E-02,
     *  0.1715944578126496E-02,0.1619406828259840E-02,
     *  0.1042890131991612E-02,0.3249077634071337E-03/
c 
c     Accuracy for the following 32 nodes is 0.82534E-07
c 
        data xs32/
     1  -.9963805246894953E+00,-.9729168124014321E+00,
     2  -.9150193617237711E+00,-.8157327748721346E+00,
     3  -.6749261171917619E+00,-.4977225841014068E+00,
     4  -.2931614187803758E+00,-.7294894971441390E-01,
     5  0.1497985965211959E+00,0.3618410337227100E+00,
     6  0.5511070019105555E+00,0.7079640979751909E+00,
     7  0.8264588420696311E+00,0.9055523349139916E+00,
     8  0.9502251499055197E+00,0.9710142188550930E+00,
     9  0.9793817488995857E+00,0.9826662959326564E+00,
     *  0.9842177896219917E+00,0.9853234076373951E+00,
     *  0.9863133151594141E+00,0.9872428581006849E+00,
     *  0.9881241707577875E+00,0.9889587256312613E+00,
     *  0.9897427395930891E+00,0.9904743439719168E+00,
     *  0.9911640762747589E+00,0.9918930762897773E+00,
     *  0.9930314370471764E+00,0.9950618127613828E+00,
     *  0.9976998722968224E+00,0.9995919677688692E+00/
c 
        data whts32/
     1  0.1058807244671761E-01,0.3886932722619676E-01,
     2  0.7805813298499946E-01,0.1205102893489760E+00,
     3  0.1602273498711281E+00,0.1926525615091400E+00,
     4  0.2144977731413253E+00,0.2237196253659089E+00,
     5  0.2195471627195384E+00,0.2025046161900928E+00,
     6  0.1744013145319553E+00,0.1383001032818626E+00,
     7  0.9849552890779739E-01,0.6051839234659109E-01,
     8  0.3066949352852654E-01,0.1292750928563102E-01,
     9  0.5010131945786280E-02,0.2079626002806452E-02,
     *  0.1223208095572456E-02,0.1030081240626323E-02,
     *  0.9563326218108897E-03,0.9045150700442578E-03,
     *  0.8583050416775265E-03,0.8100402466376127E-03,
     *  0.7575396454861492E-03,0.7070106901500520E-03,
     *  0.6811912138169024E-03,0.8394095477820972E-03,
     *  0.1535720659109331E-02,0.2488735052459556E-02,
     *  0.2527490925443852E-02,0.1103369702629029E-02/
c 
        if(n .eq. 46) call o1rarrcp(xs46,xs,n)
        if(n .eq. 46) call o1rarrcp(whts46,whts,n)
c 
        if(n .eq. 19) call o1rarrcp(xs19,xs,n)
        if(n .eq. 19) call o1rarrcp(whts19,whts,n)
c 
        if(n .eq. 32) call o1rarrcp(xs32,xs,n)
        if(n .eq. 32) call o1rarrcp(whts32,whts,n)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o1rnear15(n,xs,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),whts(1),
     1      xs49(49),whts49(49),xs33(33),whts33(33),
     2      xs20(20),whts20(20)
c 
c     Accuracy for the following 19 nodes is 0.42910E-04
c 
        data xs20/
     1  -.9880136227046848E+00,-.9332381208580276E+00,
     2  -.8282500145274310E+00,-.6743058276135518E+00,
     3  -.4790343564139045E+00,-.2543384375075934E+00,
     4  -.1485676894201582E-01,0.2235699547436343E+00,
     5  0.4453292669147049E+00,0.6365025429427416E+00,
     6  0.7863296866306306E+00,0.8886794681643668E+00,
     7  0.9442335820066217E+00,0.9647412616734274E+00,
     8  0.9712404078698680E+00,0.9755847858102042E+00,
     9  0.9793820546588288E+00,0.9829460908323309E+00,
     *  0.9888081372900890E+00,0.9970953355098209E+00/
c 
        data whts20/
     1  0.3131288882059081E-01,0.7941992435136289E-01,
     2  0.1302949183267632E+00,0.1762868591583314E+00,
     3  0.2122443503093658E+00,0.2346749305753046E+00,
     4  0.2416223549069021E+00,0.2326023270452124E+00,
     5  0.2085676376397878E+00,0.1719605602064378E+00,
     6  0.1266383374376471E+00,0.7810211645443990E-01,
     7  0.3489844945377392E-01,0.1010847558119056E-01,
     8  0.4781881092208666E-02,0.4032608108752806E-02,
     9  0.3565398073216288E-02,0.4010467042976580E-02,
     *  0.7978686184618873E-02,0.6882346676549376E-02/
c 
c     Accuracy for the following 46 nodes is 0.95466E-10
c 
        data xs49/
     1  -.9985420465181852E+00,-.9877524831815386E+00,
     2  -.9572156517038481E+00,-.8985100731303609E+00,
     3  -.8071105925452144E+00,-.6825990977444534E+00,
     4  -.5282482808382472E+00,-.3503582523696416E+00,
     5  -.1574491315669765E+00,0.4065028946935256E-01,
     6  0.2337775378526268E+00,0.4124252104086101E+00,
     7  0.5686879173927875E+00,0.6971090104428008E+00,
     8  0.7953546436126946E+00,0.8645798216116130E+00,
     9  0.9091965406442268E+00,0.9356650712237669E+00,
     *  0.9504867507894619E+00,0.9586021255669630E+00,
     *  0.9630873846243799E+00,0.9656836346196109E+00,
     *  0.9673736920028797E+00,0.9686856466837142E+00,
     *  0.9698536711542347E+00,0.9709731751136097E+00,
     *  0.9720742159718307E+00,0.9731384609635416E+00,
     *  0.9741471301593763E+00,0.9750981908214449E+00,
     *  0.9759985586704423E+00,0.9768577406621928E+00,
     *  0.9776829290456925E+00,0.9784777107291216E+00,
     *  0.9792439867547653E+00,0.9799862981097029E+00,
     *  0.9807141625162767E+00,0.9814382978039498E+00,
     *  0.9821659189798990E+00,0.9829086544677148E+00,
     *  0.9837063185069359E+00,0.9846799077108093E+00,
     *  0.9860340800758263E+00,0.9879716413087398E+00,
     *  0.9905728908622824E+00,0.9936460247499503E+00,
     *  0.9966427964218787E+00,0.9988398809629191E+00,
     *  0.9998420148627241E+00/
c 
        data whts49/
     1  0.4420971506304451E-02,0.1898601174454147E-01,
     2  0.4351805692483427E-01,0.7464842300044505E-01,
     3  0.1082125708668402E+00,0.1402537416517992E+00,
     4  0.1673843394643252E+00,0.1869604142502019E+00,
     5  0.1972033243213229E+00,0.1972851589629268E+00,
     6  0.1873688518708415E+00,0.1685946604628760E+00,
     7  0.1430098001214617E+00,0.1134366450210872E+00,
     8  0.8324981335953446E-01,0.5595598854188483E-01,
     9  0.3439125624291345E-01,0.1965887406436326E-01,
     *  0.1081104717848545E-01,0.5921468900480820E-02,
     *  0.3328867677888729E-02,0.2022919079218981E-02,
     *  0.1442988204472999E-02,0.1216959010128906E-02,
     *  0.1133829038326517E-02,0.1110359145944752E-02,
     *  0.1087424493928684E-02,0.1037842584606460E-02,
     *  0.9792838033192098E-03,0.9241150644777138E-03,
     *  0.8782874374841000E-03,0.8412948976352410E-03,
     *  0.8096493558200652E-03,0.7801438303667711E-03,
     *  0.7531171117327090E-03,0.7331531396127631E-03,
     *  0.7245188336018507E-03,0.7247403169980735E-03,
     *  0.7322361119624610E-03,0.7586539253903938E-03,
     *  0.8566811640411535E-03,0.1126300794124004E-02,
     *  0.1617736798064893E-02,0.2273481804725220E-02,
     *  0.2899056465208595E-02,0.3153069277626194E-02,
     *  0.2710697636236905E-02,0.1608026563824911E-02,
     *  0.4631479638558236E-03/
c 
c     Accuracy for the following 33 nodes is 0.65551E-07
c 
        data xs33/
     1  -.9965189269728825E+00,-.9738175539185558E+00,
     2  -.9175043938046168E+00,-.8205610996483543E+00,
     3  -.6826768277423354E+00,-.5087383657190856E+00,
     4  -.3075229266475127E+00,-.9047449942702557E-01,
     5  0.1295317195004919E+00,0.3394584154508671E+00,
     6  0.5273873515167808E+00,0.6837839206583003E+00,
     7  0.8027348108229722E+00,0.8831666042639181E+00,
     8  0.9298353127369772E+00,0.9526628457050410E+00,
     9  0.9624997254115235E+00,0.9667470581851196E+00,
     *  0.9691472765532472E+00,0.9710995536257137E+00,
     *  0.9728984738577491E+00,0.9745897710896827E+00,
     *  0.9761855442114595E+00,0.9776928215265523E+00,
     *  0.9791175317517131E+00,0.9804658204650536E+00,
     *  0.9817486011138303E+00,0.9830167886792189E+00,
     *  0.9845691705643382E+00,0.9871851078806402E+00,
     *  0.9913932897476530E+00,0.9962173036376347E+00,
     *  0.9993658965985832E+00/
c 
        data whts33/
     1  0.1020120232826782E-01,0.3770053046962130E-01,
     2  0.7607114907997933E-01,0.1178415390460715E+00,
     3  0.1570911999388290E+00,0.1893026692851739E+00,
     4  0.2112003068743641E+00,0.2207294523499008E+00,
     5  0.2170886785576979E+00,0.2007572259841265E+00,
     6  0.1734909607925052E+00,0.1382962109985931E+00,
     7  0.9940744018642600E-01,0.6225116125749172E-01,
     8  0.3282010010222138E-01,0.1473409180164829E-01,
     9  0.6162038690835221E-02,0.2933529478456595E-02,
     *  0.2077835467662119E-02,0.1861233819633880E-02,
     *  0.1742156599104125E-02,0.1642184686136000E-02,
     *  0.1550491676948883E-02,0.1465036389670489E-02,
     *  0.1385366932009925E-02,0.1312708230606793E-02,
     *  0.1258608286184111E-02,0.1318908634060127E-02,
     *  0.1930372413932418E-02,0.3414969236813076E-02,
     *  0.4840836393754665E-02,0.4370137786771564E-02,
     *  0.1749628926118319E-02/
c 
        if(n .eq. 49) call o1rarrcp(xs49,xs,n)
        if(n .eq. 49) call o1rarrcp(whts49,whts,n)
c 
        if(n .eq. 20) call o1rarrcp(xs20,xs,n)
        if(n .eq. 20) call o1rarrcp(whts20,whts,n)
c 
        if(n .eq. 33) call o1rarrcp(xs33,xs,n)
        if(n .eq. 33) call o1rarrcp(whts33,whts,n)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o1rnear16(n,xs,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),whts(1),
     1      xs41(41),whts41(41),xs30(30),whts30(30),
     2      xs18(18),whts18(18)
c 
c     Accuracy for the following 18 nodes is 0.55303E-04
c 
        data xs18/
     1  -.9864453594276791E+00,-.9264797107596426E+00,
     2  -.8149943975543059E+00,-.6544584198778441E+00,
     3  -.4532994199566620E+00,-.2240127108193909E+00,
     4  0.1831357459253868E-01,0.2574601657150379E+00,
     5  0.4775092271165435E+00,0.6642336167249405E+00,
     6  0.8065220427786908E+00,0.8978716629533068E+00,
     7  0.9406264211185500E+00,0.9541837521533942E+00,
     8  0.9614581851346153E+00,0.9686243105182759E+00,
     9  0.9813559357259296E+00,0.9958300735931185E+00/
c 
        data whts18/
     1  0.3510075109249840E-01,0.8554825638742670E-01,
     2  0.1369864794783037E+00,0.1826432533083985E+00,
     3  0.2175647246183333E+00,0.2384627941174270E+00,
     4  0.2435228611569115E+00,0.2321523826397444E+00,
     5  0.2055420907402908E+00,0.1660322256353989E+00,
     6  0.1174077063973016E+00,0.6547725217700583E-01,
     7  0.2342980525242511E-01,0.8186783440416878E-02,
     8  0.6777113545290807E-02,0.8813027485829366E-02,
     9  0.1598385694531550E-01,0.1041483384976336E-01/
c 
c     Accuracy for the following 41 nodes is 0.92323E-10
c 
        data xs41/
     1  -.9984518832748875E+00,-.9870684859497694E+00,
     2  -.9550990757126782E+00,-.8940907364394988E+00,
     3  -.7997364384516091E+00,-.6719810895723380E+00,
     4  -.5145263719827534E+00,-.3341094543692634E+00,
     5  -.1396487093752201E+00,0.5870848467191703E-01,
     6  0.2505980543982791E+00,0.4264636051453816E+00,
     7  0.5785392658409385E+00,0.7017331747013233E+00,
     8  0.7943268294377219E+00,0.8582935304892487E+00,
     9  0.8988198497161679E+00,0.9226858003802185E+00,
     *  0.9361397447348276E+00,0.9436643160527305E+00,
     *  0.9480098126359273E+00,0.9507637500472228E+00,
     *  0.9528181026385445E+00,0.9545953907107008E+00,
     *  0.9562486195716209E+00,0.9578265190579329E+00,
     *  0.9593454600291069E+00,0.9608124974680545E+00,
     *  0.9622327951392124E+00,0.9636138328060202E+00,
     *  0.9649731539837327E+00,0.9663608922837286E+00,
     *  0.9679164370049591E+00,0.9699326093701812E+00,
     *  0.9728097102788240E+00,0.9768984727305257E+00,
     *  0.9822465452986844E+00,0.9883125475534145E+00,
     *  0.9939606068420603E+00,0.9979458088248330E+00,
     *  0.9997226076762573E+00/
c 
        data whts41/
     1  0.4686881691433312E-02,0.1996777812533774E-01,
     2  0.4540505950202818E-01,0.7732874725927620E-01,
     3  0.1113800962745751E+00,0.1435014853069565E+00,
     4  0.1702707834635562E+00,0.1890606970867717E+00,
     5  0.1981525742856493E+00,0.1968172661130525E+00,
     6  0.1853526282704079E+00,0.1650692657518594E+00,
     7  0.1382213575789085E+00,0.1078736991471098E+00,
     8  0.7764912423838732E-01,0.5117712741434077E-01,
     9  0.3104546986778323E-01,0.1774918537233671E-01,
     *  0.9902767391090210E-02,0.5593786986349183E-02,
     *  0.3353514467476276E-02,0.2300371249337047E-02,
     *  0.1875034692011618E-02,0.1702174954659420E-02,
     *  0.1611370485667605E-02,0.1546843493904640E-02,
     *  0.1492114191939761E-02,0.1442763789736704E-02,
     *  0.1398971749939372E-02,0.1365708954025959E-02,
     *  0.1360361171704336E-02,0.1436805952498844E-02,
     *  0.1723989250263542E-02,0.2377687070206465E-02,
     *  0.3439190758900389E-02,0.4751992721524858E-02,
     *  0.5854689285800602E-02,0.6079340093887482E-02,
     *  0.4994869156919591E-02,0.2871752234119875E-02,
     *  0.8146731383719239E-03/
c 
c     Accuracy for the following 30 nodes is 0.47148E-07
c 
        data xs30/
     1  -.9966314423044634E+00,-.9745513989256621E+00,
     2  -.9195321580760464E+00,-.8245107541426780E+00,
     3  -.6890441622702050E+00,-.5178476793203031E+00,
     4  -.3195101090018845E+00,-.1052877217145572E+00,
     5  0.1121165727000426E+00,0.3198120643972180E+00,
     6  0.5059989958124289E+00,0.6612463933363537E+00,
     7  0.7797662218571697E+00,0.8606592924547376E+00,
     8  0.9087363427760697E+00,0.9334457740749680E+00,
     9  0.9449151256449528E+00,0.9503543509033414E+00,
     *  0.9536743040780151E+00,0.9564155957192240E+00,
     *  0.9589454856513582E+00,0.9613363228929139E+00,
     *  0.9636200631110787E+00,0.9658714278156690E+00,
     *  0.9684208476946770E+00,0.9722560155418367E+00,
     *  0.9785204056853025E+00,0.9868449307440410E+00,
     *  0.9946889572535356E+00,0.9991577343450833E+00/
c 
        data whts30/
     1  0.9886383121060407E-02,0.3674714918415106E-01,
     2  0.7444607473051704E-01,0.1156455382643144E+00,
     3  0.1544807565675923E+00,0.1864608516687242E+00,
     4  0.2083179248281146E+00,0.2179886507631926E+00,
     5  0.2146501488336791E+00,0.1987529047311880E+00,
     6  0.1720279199723024E+00,0.1374809096614142E+00,
     7  0.9939159760827804E-01,0.6320864321027407E-01,
     8  0.3459663453473977E-01,0.1658059477805574E-01,
     9  0.7555989811696561E-02,0.3959107559604990E-02,
     *  0.2914632683530610E-02,0.2614329734104871E-02,
     *  0.2454700840597315E-02,0.2331415640418626E-02,
     *  0.2245830755490532E-02,0.2301885264285725E-02,
     *  0.2968016448948995E-02,0.4927300427824965E-02,
     *  0.7556573798280933E-02,0.8632953403030490E-02,
     *  0.6518184117689359E-02,0.2356380774786553E-02/
c 
        if(n .eq. 41) call o1rarrcp(xs41,xs,n)
        if(n .eq. 41) call o1rarrcp(whts41,whts,n)
c 
        if(n .eq. 18) call o1rarrcp(xs18,xs,n)
        if(n .eq. 18) call o1rarrcp(whts18,whts,n)
c 
        if(n .eq. 30) call o1rarrcp(xs30,xs,n)
        if(n .eq. 30) call o1rarrcp(whts30,whts,n)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o1rnear17(n,xs,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),whts(1),
     1      xs44(44),whts44(44),xs30(30),whts30(30),
     2      xs18(18),whts18(18)
c 
c     Accuracy for the following 18 nodes is 0.57711E-04
c 
        data xs18/
     1  -.9868467740534091E+00,-.9282891950406170E+00,
     2  -.8187204376970794E+00,-.6603082982571096E+00,
     3  -.4612670002720597E+00,-.2339288989702200E+00,
     4  0.6734768098090482E-02,0.2445468656051224E+00,
     5  0.4635821354079790E+00,0.6494438870380099E+00,
     6  0.7905575215306064E+00,0.8796649029216470E+00,
     7  0.9199788076611692E+00,0.9346256583502504E+00,
     8  0.9443536374177530E+00,0.9556809079683510E+00,
     9  0.9755803686783986E+00,0.9947569852486560E+00/
c 
        data whts18/
     1  0.3411714912952722E-01,0.8381795979587045E-01,
     2  0.1349278813886131E+00,0.1804900540378549E+00,
     3  0.2155093318172006E+00,0.2366401488767792E+00,
     4  0.2420138325565007E+00,0.2309925485336837E+00,
     5  0.2046484156527663E+00,0.1651264457913229E+00,
     6  0.1158155730176306E+00,0.6262465708275283E-01,
     7  0.2237415027027990E-01,0.1064639573571556E-01,
     8  0.9343916340506599E-02,0.1520039379788673E-01,
     9  0.2256364119434296E-01,0.1318548170989140E-01/
c 
c     Accuracy for the following 44 nodes is 0.18304E-10
c 
        data xs44/
     1  -.9986752316588204E+00,-.9887997227887785E+00,
     2  -.9606041965008213E+00,-.9059519138341018E+00,
     3  -.8202486330668266E+00,-.7027610334579432E+00,
     4  -.5562911386824991E+00,-.3865826629025033E+00,
     5  -.2015735667349375E+00,-.1054019969908171E-01,
     6  0.1768269348173412E+00,0.3513694677281733E+00,
     7  0.5053793601359298E+00,0.6334354395205690E+00,
     8  0.7330819456505233E+00,0.8051850346754028E+00,
     9  0.8536445584738952E+00,0.8841856072045056E+00,
     *  0.9026154195431411E+00,0.9135588323400563E+00,
     *  0.9201681062699554E+00,0.9244206608658598E+00,
     *  0.9274989668904645E+00,0.9300481316761851E+00,
     *  0.9323601072745218E+00,0.9345464409864941E+00,
     *  0.9366479642786469E+00,0.9386815782632056E+00,
     *  0.9406560761448801E+00,0.9425768094730556E+00,
     *  0.9444485776041720E+00,0.9462811274876255E+00,
     *  0.9481012553002532E+00,0.9499825910837012E+00,
     *  0.9521101693476934E+00,0.9548376777782295E+00,
     *  0.9586263517036531E+00,0.9638933733860153E+00,
     *  0.9707815869949619E+00,0.9788657026904610E+00,
     *  0.9870323519384314E+00,0.9937993135654795E+00,
     *  0.9980499908633349E+00,0.9997533550813731E+00/
c 
        data whts44/
     1  0.4024347901286546E-02,0.1743983764070897E-01,
     2  0.4033522181213042E-01,0.6973816251503946E-01,
     3  0.1017844534029825E+00,0.1327136452007894E+00,
     4  0.1592563122942012E+00,0.1788232885852434E+00,
     5  0.1896333294294869E+00,0.1908024286801919E+00,
     6  0.1823900200009715E+00,0.1653962735216803E+00,
     7  0.1417083312050840E+00,0.1139877827290412E+00,
     8  0.8545339579889200E-01,0.5942371334687570E-01,
     9  0.3848170055029159E-01,0.2359068079911143E-01,
     *  0.1404317953921381E-01,0.8362716620732222E-02,
     *  0.5178640636798962E-02,0.3519456522084853E-02,
     *  0.2741860954834317E-02,0.2401966616557465E-02,
     *  0.2238809341455035E-02,0.2139965962840851E-02,
     *  0.2065680640570693E-02,0.2002954118793208E-02,
     *  0.1946894597430356E-02,0.1895286574335267E-02,
     *  0.1849695711902418E-02,0.1819433830924174E-02,
     *  0.1831983397919642E-02,0.1960098405449950E-02,
     *  0.2355275952881575E-02,0.3179125566551350E-02,
     *  0.4471481113773323E-02,0.6091946398986981E-02,
     *  0.7616382633417922E-02,0.8360961085095939E-02,
     *  0.7711957438371999E-02,0.5622129140202668E-02,
     *  0.2871973287333653E-02,0.7372184952864913E-03/
c 
c     Accuracy for the following 30 nodes is 0.46617E-07
c 
        data xs30/
     1  -.9966558681932611E+00,-.9747251125694724E+00,
     2  -.9200572515398997E+00,-.8256267534675095E+00,
     3  -.6910031477436199E+00,-.5208974415626307E+00,
     4  -.3238821322585812E+00,-.1111911087196527E+00,
     5  0.1044926907386632E+00,0.3102852694179496E+00,
     6  0.4943699327642157E+00,0.6472771375698154E+00,
     7  0.7632103880798176E+00,0.8414532589421106E+00,
     8  0.8873396670186091E+00,0.9108252210662433E+00,
     9  0.9220417664459537E+00,0.9280050251387040E+00,
     *  0.9323154385559151E+00,0.9361456552164746E+00,
     *  0.9397373763193864E+00,0.9431463279685639E+00,
     *  0.9464155986185001E+00,0.9497065336204384E+00,
     *  0.9537683611930838E+00,0.9602974412557399E+00,
     *  0.9704760651138369E+00,0.9828435433568925E+00,
     *  0.9934497987271764E+00,0.9990070793944097E+00/
c 
        data whts30/
     1  0.9816164601256762E-02,0.3650580075600380E-01,
     2  0.7397871101582107E-01,0.1149291127557334E+00,
     3  0.1535119961762497E+00,0.1852508727278372E+00,
     4  0.2068872545946099E+00,0.2163601091697653E+00,
     5  0.2128390250983030E+00,0.1967554297612950E+00,
     6  0.1698137627849856E+00,0.1350105560679975E+00,
     7  0.9671502706260003E-01,0.6067812586322237E-01,
     8  0.3284369654806449E-01,0.1587590395120554E-01,
     9  0.7722377127524687E-02,0.4790342785240247E-02,
     *  0.3994426496343875E-02,0.3696142998526001E-02,
     *  0.3494754427808042E-02,0.3328955021614582E-02,
     *  0.3229472194613484E-02,0.3462365110840222E-02,
     *  0.4974379800680566E-02,0.8323372351189753E-02,
     *  0.1179210507095520E-01,0.1222055784278421E-01,
     *  0.8387206118025289E-02,0.2811987506583751E-02/
c 
        if(n .eq. 44) call o1rarrcp(xs44,xs,n)
        if(n .eq. 44) call o1rarrcp(whts44,whts,n)
c 
        if(n .eq. 18) call o1rarrcp(xs18,xs,n)
        if(n .eq. 18) call o1rarrcp(whts18,whts,n)
c 
        if(n .eq. 30) call o1rarrcp(xs30,xs,n)
        if(n .eq. 30) call o1rarrcp(whts30,whts,n)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o1rnear18(n,xs,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),whts(1),
     1      xs42(42),whts42(42),xs30(30),whts30(30),
     2      xs18(18),whts18(18)
c 
c     Accuracy for the following 18 nodes is 0.66731E-04
c 
        data xs18/
     1  -.9871353616096061E+00,-.9296476274644763E+00,
     2  -.8216720068656471E+00,-.6651804345948726E+00,
     3  -.4682142896619188E+00,-.2429536541138250E+00,
     4  -.4244630195830212E-02,0.2318042667112947E+00,
     5  0.4492272345565596E+00,0.6333524216698774E+00,
     6  0.7718324047737901E+00,0.8561464773174225E+00,
     7  0.8918517488786466E+00,0.9080872258021793E+00,
     8  0.9218494368726249E+00,0.9400603137095496E+00,
     9  0.9690563261465657E+00,0.9935658439883950E+00/
c 
        data whts18/
     1  0.3340138501095948E-01,0.8244777956928043E-01,
     2  0.1331415772285462E+00,0.1784644429214240E+00,
     3  0.2134096420157639E+00,0.2346049644831795E+00,
     4  0.2400859344871837E+00,0.2293259699585359E+00,
     5  0.2030490717559278E+00,0.1631299430535836E+00,
     6  0.1122941262562602E+00,0.5673115707277308E-01,
     7  0.2077472653236653E-01,0.1415168274778144E-01,
     8  0.1406371731586154E-01,0.2439658190577672E-01,
     9  0.3023102720527266E-01,0.1626481605566285E-01/
c 
c     Accuracy for the following 44 nodes is 0.54446E-10
c 
        data xs42/
     1  -.9985467055736298E+00,-.9878124118338006E+00,
     2  -.9575109037556991E+00,-.8994252528607651E+00,
     3  -.8092673607984071E+00,-.6868484607200969E+00,
     4  -.5356327148989311E+00,-.3620520768029110E+00,
     5  -.1746837248049182E+00,0.1667214248546478E-01,
     6  0.2019779498249135E+00,0.3719649353496925E+00,
     7  0.5191146011341834E+00,0.6385624307910005E+00,
     8  0.7288126894470930E+00,0.7919761317539464E+00,
     9  0.8330737978818006E+00,0.8583442023314927E+00,
     *  0.8734452039327061E+00,0.8825643695718024E+00,
     *  0.8884789475193995E+00,0.8929021399570554E+00,
     *  0.8967065345578873E+00,0.9002323630611286E+00,
     *  0.9035948197790975E+00,0.9068335602486506E+00,
     *  0.9099646746819552E+00,0.9129977421440978E+00,
     *  0.9159425985221070E+00,0.9188158468718435E+00,
     *  0.9216563295970484E+00,0.9245681010647047E+00,
     *  0.9278175239105446E+00,0.9319336608452418E+00,
     *  0.9376427880310818E+00,0.9456196687888017E+00,
     *  0.9560694437260828E+00,0.9682728538948835E+00,
     *  0.9805106538379704E+00,0.9906231293904983E+00,
     *  0.9970169186340880E+00,0.9996176069655289E+00/
c 
        data whts42/
     1  0.4404813477702659E-02,0.1886978540698850E-01,
     2  0.4312874215320988E-01,0.7375769753286396E-01,
     3  0.1065813983290898E+00,0.1376674542086582E+00,
     4  0.1636782633185810E+00,0.1820380278987268E+00,
     5  0.1910487772815578E+00,0.1899731501457986E+00,
     6  0.1790772028923415E+00,0.1596284135998307E+00,
     7  0.1338476134561076E+00,0.1047934710012925E+00,
     8  0.7606986542160366E-01,0.5112465382995770E-01,
     9  0.3214346572791276E-01,0.1935495443778440E-01,
     *  0.1154207568303523E-01,0.7151866397605753E-02,
     *  0.4960607813296450E-02,0.4025320466839076E-02,
     *  0.3634671730261994E-02,0.3433962856738796E-02,
     *  0.3296835094119038E-02,0.3183091193824616E-02,
     *  0.3080640708839711E-02,0.2987043478979291E-02,
     *  0.2905217689668704E-02,0.2847124686365088E-02,
     *  0.2849541254867480E-02,0.3016057539446419E-02,
     *  0.3572218709485908E-02,0.4785168486155543E-02,
     *  0.6752709244445959E-02,0.9244072043804489E-02,
     *  0.1153761762445302E-01,0.1257275928769299E-01,
     *  0.1152490759311387E-01,0.8416847432018587E-02,
     *  0.4353475281699349E-02,0.1138417577212781E-02/
c 
c     Accuracy for the following 30 nodes is 0.52658E-07
c 
        data xs30/
     1  -.9967324581000744E+00,-.9752514279816460E+00,
     2  -.9215856085297106E+00,-.8287388052542121E+00,
     3  -.6962276495465892E+00,-.5286647377977138E+00,
     4  -.3344994696976313E+00,-.1248361685884686E+00,
     5  0.8775969161366026E-01,0.2904904979960017E+00,
     6  0.4715794559908917E+00,0.6215533841286649E+00,
     7  0.7346290267514355E+00,0.8102690552809830E+00,
     8  0.8542140986516980E+00,0.8767711283735862E+00,
     9  0.8881951799099304E+00,0.8953612998725675E+00,
     *  0.9012867142153336E+00,0.9067506823246854E+00,
     *  0.9119146714292810E+00,0.9168394561354245E+00,
     *  0.9216331654931744E+00,0.9267575598890252E+00,
     *  0.9337634237668722E+00,0.9449554413037213E+00,
     *  0.9608419846362795E+00,0.9782748587574324E+00,
     *  0.9920570487189349E+00,0.9988397476940007E+00/
c 
        data whts30/
     1  0.9598175454112040E-02,0.3579459691129517E-01,
     2  0.7268216839885201E-01,0.1130676265816752E+00,
     3  0.1511653049637973E+00,0.1825327920327639E+00,
     4  0.2039271587208027E+00,0.2132846578691763E+00,
     5  0.2097532695666283E+00,0.1937243335959939E+00,
     6  0.1668519969866457E+00,0.1321051904857606E+00,
     7  0.9393266245490288E-01,0.5832904928128526E-01,
     8  0.3137233880083117E-01,0.1550724066805889E-01,
     9  0.8494667709689513E-02,0.6305031666890538E-02,
     *  0.5647074108337922E-02,0.5301777913438576E-02,
     *  0.5034711854681108E-02,0.4829216421339751E-02,
     *  0.4819653103641632E-02,0.5686215135347644E-02,
     *  0.8762473108541878E-02,0.1374608119758660E-01,
     *  0.1747109744035879E-01,0.1645429417094116E-01,
     *  0.1049790082013584E-01,0.3321233406836524E-02/
c 
        if(n .eq. 42) call o1rarrcp(xs42,xs,n)
        if(n .eq. 42) call o1rarrcp(whts42,whts,n)
c 
        if(n .eq. 18) call o1rarrcp(xs18,xs,n)
        if(n .eq. 18) call o1rarrcp(whts18,whts,n)
c 
        if(n .eq. 30) call o1rarrcp(xs30,xs,n)
        if(n .eq. 30) call o1rarrcp(whts30,whts,n)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o1rnear19(n,xs,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),whts(1),
     1      xs43(43),whts43(43),xs30(30),whts30(30),
     2      xs19(19),whts19(19)
c 
c     Accuracy for the following 19 nodes is 0.38705E-04
c 
        data xs19/
     1  -.9885146827553295E+00,-.9356416523691863E+00,
     2  -.8336943496300504E+00,-.6837734254787533E+00,
     3  -.4933419959643988E+00,-.2741392859546719E+00,
     4  -.4068469759882876E-01,0.1911299097411616E+00,
     5  0.4053075823165959E+00,0.5868559993712201E+00,
     6  0.7226589190137489E+00,0.8040048000510958E+00,
     7  0.8399365997311896E+00,0.8601323541534251E+00,
     8  0.8774418888382876E+00,0.8968868934850998E+00,
     9  0.9288535426907689E+00,0.9669758112153908E+00,
     *  0.9934621496250802E+00/
c 
        data whts19/
     1  0.3006710784213597E-01,0.7691775055941816E-01,
     2  0.1267371220647617E+00,0.1718264281067604E+00,
     3  0.2070518378119517E+00,0.2289027734480900E+00,
     4  0.2353350994101181E+00,0.2256105058577023E+00,
     5  0.2002287997481147E+00,0.1606605144094307E+00,
     6  0.1092960311912773E+00,0.5446578873108627E-01,
     7  0.2372544951594041E-01,0.1818688474509940E-01,
     8  0.1693711265726067E-01,0.2444084281323050E-01,
     9  0.3809602043982701E-01,0.3481355526211845E-01,
     *  0.1668652891655180E-01/
c 
c     Accuracy for the following 44 nodes is 0.54446E-10
c 
        data xs43/
     1  -.9986359750993765E+00,-.9885172577792105E+00,
     2  -.9598084313656292E+00,-.9045243413137403E+00,
     3  -.8183905310168055E+00,-.7010749824863487E+00,
     4  -.5557919402023199E+00,-.3886550555000944E+00,
     5  -.2078881411487994E+00,-.2293409967655713E-01,
     6  0.1565027390310123E+00,0.3214414182962355E+00,
     7  0.4645946091004631E+00,0.5812669679251670E+00,
     8  0.6700441394674246E+00,0.7329540851414387E+00,
     9  0.7747079801121730E+00,0.8011098313784488E+00,
     *  0.8174851253336874E+00,0.8279211885580120E+00,
     *  0.8352457427992945E+00,0.8411558261250105E+00,
     *  0.8464468435843862E+00,0.8514220024960759E+00,
     *  0.8561887798779870E+00,0.8607862962430528E+00,
     *  0.8652314168374806E+00,0.8695352351892900E+00,
     *  0.8737108895930517E+00,0.8777810166493443E+00,
     *  0.8817909566695936E+00,0.8858504008531456E+00,
     *  0.8902409297631153E+00,0.8955786912569024E+00,
     *  0.9028029144480085E+00,0.9128564743407319E+00,
     *  0.9262473447419204E+00,0.9425430149378604E+00,
     *  0.9601009305796624E+00,0.9764189731214172E+00,
     *  0.9890445997802384E+00,0.9966151283272868E+00,
     *  0.9995752768265327E+00/
c 
        data whts43/
     1  0.4138735243339053E-02,0.1782509340414737E-01,
     2  0.4095005061186219E-01,0.7033212234664993E-01,
     3  0.1019835040117409E+00,0.1321009178512773E+00,
     4  0.1574329924131794E+00,0.1754552078902661E+00,
     5  0.1844887464365902E+00,0.1837854380220086E+00,
     6  0.1735751290315002E+00,0.1550724308139968E+00,
     7  0.1304391508060798E+00,0.1026659418272166E+00,
     8  0.7523901617392915E-01,0.5139507606623314E-01,
     9  0.3310853276757555E-01,0.2059510021894518E-01,
     *  0.1283495760568686E-01,0.8506520667512106E-02,
     *  0.6422194221785145E-02,0.5522452358383964E-02,
     *  0.5105257876806117E-02,0.4861108460722274E-02,
     *  0.4678307114510265E-02,0.4519302851997402E-02,
     *  0.4372678888814090E-02,0.4237021711111099E-02,
     *  0.4117901377514675E-02,0.4029305047798578E-02,
     *  0.4007160329330397E-02,0.4155325479031287E-02,
     *  0.4727889452268546E-02,0.6110629835717449E-02,
     *  0.8501808260267880E-02,0.1169935441568587E-01,
     *  0.1501598107997984E-01,0.1729554804523438E-01,
     *  0.1738687749701429E-01,0.1482209729832222E-01,
     *  0.1018890782265565E-01,0.5026877223195112E-02,
     *  0.1271349136643373E-02/
c 
c     Accuracy for the following 30 nodes is  0.59161E-07
c 
        data xs30/
     1  -.9967985942834480E+00,-.9757307767662382E+00,
     2  -.9230477909040439E+00,-.8318462493433237E+00,
     3  -.7016464813168253E+00,-.5370081990395146E+00,
     4  -.3462935367866341E+00,-.1405069526228966E+00,
     5  0.6787347732414983E-01,0.2660972216052835E+00,
     6  0.4423735668584816E+00,0.5871791606240204E+00,
     7  0.6948106539565526E+00,0.7652884830845084E+00,
     8  0.8054400082321491E+00,0.8264594530752156E+00,
     9  0.8386579551801088E+00,0.8479314753236264E+00,
     *  0.8562497297992062E+00,0.8640613002723030E+00,
     *  0.8714874552982868E+00,0.8786233098400262E+00,
     *  0.8857494179272005E+00,0.8939957114194867E+00,
     *  0.9061617915147147E+00,0.9249058508119505E+00,
     *  0.9490141020824276E+00,0.9729423472512069E+00,
     *  0.9904624916955047E+00,0.9986446404343554E+00/
c 
        data whts30/
     1  0.9406638517722088E-02,0.3512173285632193E-01,
     2  0.7137477263428295E-01,0.1110843899866518E+00,
     3  0.1485349276839223E+00,0.1793288106683525E+00,
     4  0.2002469917655081E+00,0.2092264370231629E+00,
     5  0.2053888873280013E+00,0.1890720957955818E+00,
     6  0.1618682152828599E+00,0.1267587912536841E+00,
     7  0.8848016536578836E-01,0.5365456856568548E-01,
     8  0.2861758469827339E-01,0.1520303089057171E-01,
     9  0.1016346597136840E-01,0.8660709067188455E-02,
     *  0.8034243015803156E-02,0.7605732312517445E-02,
     *  0.7260027087668189E-02,0.7048861570285367E-02,
     *  0.7367110133479281E-02,0.9623641213378586E-02,
     *  0.1523088605641123E-01,0.2206522410868308E-01,
     *  0.2515229667677991E-01,0.2159744668417160E-01,
     *  0.1291318833242226E-01,0.3909116430742726E-02/
c 
        if(n .eq. 43) call o1rarrcp(xs43,xs,n)
        if(n .eq. 43) call o1rarrcp(whts43,whts,n)
c 
        if(n .eq. 19) call o1rarrcp(xs19,xs,n)
        if(n .eq. 19) call o1rarrcp(whts19,whts,n)
c 
        if(n .eq. 30) call o1rarrcp(xs30,xs,n)
        if(n .eq. 30) call o1rarrcp(whts30,whts,n)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o1rnear20(n,xs,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),whts(1),
     1      xs44(44),whts44(44),xs30(30),whts30(30),
     2      xs19(19),whts19(19)
c 
c     Accuracy for the following 19 nodes is 0.43306E-04
c 
        data xs19/
     1  -.9883345542720406E+00,-.9351531335129269E+00,
     2  -.8334540197208966E+00,-.6845949814418042E+00,
     3  -.4961305868724799E+00,-.2798131596992133E+00,
     4  -.5017490256576951E-01,0.1767619951863115E+00,
     5  0.3846251518861712E+00,0.5574356352643761E+00,
     6  0.6802058253726802E+00,0.7459676749111002E+00,
     7  0.7785506014213810E+00,0.8046411716624551E+00,
     8  0.8295903822756552E+00,0.8646913980214672E+00,
     9  0.9156880553074653E+00,0.9639393528801736E+00,
     *  0.9932316223103865E+00/
c 
        data whts19/
     1  0.3045442240081233E-01,0.7701805301184746E-01,
     2  0.1261054993354153E+00,0.1703193370763594E+00,
     3  0.2046242517486370E+00,0.2255617263463072E+00,
     4  0.2310340490800772E+00,0.2201156630399314E+00,
     5  0.1929436490759456E+00,0.1501250846107690E+00,
     6  0.9368999112010613E-01,0.4235712092049752E-01,
     7  0.2775168785191354E-01,0.2475192315393799E-01,
     8  0.2708073264361862E-01,0.4481350508985365E-01,
     9  0.5328540522391802E-01,0.4044953327446994E-01,
     *  0.1752407923279292E-01/
c 
c     Accuracy for the following 44 nodes is 0.15305E-10
c 
        data xs44/
     1  -.9987553554152789E+00,-.9894684186275156E+00,
     2  -.9629370904488509E+00,-.9115223274354191E+00,
     3  -.8309913085213054E+00,-.7208268658013515E+00,
     4  -.5839010775132758E+00,-.4258833862249010E+00,
     5  -.2544996270904573E+00,-.7868637106945331E-01,
     6  0.9232513174882184E-01,0.2499597864929672E+00,
     7  0.3872455774794878E+00,0.4996930078971926E+00,
     8  0.5859481946972438E+00,0.6478953495506203E+00,
     9  0.6898890008921958E+00,0.7172944129093999E+00,
     *  0.7351155411679113E+00,0.7473152861573926E+00,
     *  0.7566419124909289E+00,0.7646516930530260E+00,
     *  0.7720441419543755E+00,0.7790936444300889E+00,
     *  0.7859056175813646E+00,0.7925191761244230E+00,
     *  0.7989443600720882E+00,0.8051800275782470E+00,
     *  0.8112274338271843E+00,0.8171065884499229E+00,
     *  0.8228775828642365E+00,0.8286889150712513E+00,
     *  0.8348961619769355E+00,0.8422344166443514E+00,
     *  0.8518347911465325E+00,0.8649332122851137E+00,
     *  0.8823565187133786E+00,0.9039065607854853E+00,
     *  0.9280070161102052E+00,0.9519447310837106E+00,
     *  0.9726551579885705E+00,0.9877338047985080E+00,
     *  0.9963249022071510E+00,0.9995496416779970E+00/
c 
        data whts44/
     1  0.3781961529452406E-02,0.1640635154636039E-01,
     2  0.3795659284943467E-01,0.6558145531236892E-01,
     3  0.9556059713916012E-01,0.1242795337559090E+00,
     4  0.1486146993140115E+00,0.1661188702996758E+00,
     5  0.1751438658891669E+00,0.1749281709200115E+00,
     6  0.1656497471449701E+00,0.1484425717196822E+00,
     7  0.1253690616318300E+00,0.9929742017415821E-01,
     8  0.7354009240996216E-01,0.5110657939433919E-01,
     9  0.3379959118132995E-01,0.2185578761161119E-01,
     *  0.1444852661534554E-01,0.1040979550239392E-01,
     *  0.8496245937793299E-02,0.7631973129074454E-02,
     *  0.7194592330424264E-02,0.6920505558383585E-02,
     *  0.6709440770893786E-02,0.6519071048961936E-02,
     *  0.6330907512448170E-02,0.6140348638036404E-02,
     *  0.5957497570443658E-02,0.5810141907314115E-02,
     *  0.5754686182725000E-02,0.5923765094626439E-02,
     *  0.6612374888489480E-02,0.8258067670155810E-02,
     *  0.1115716469931943E-01,0.1518887154787908E-01,
     *  0.1963155527415418E-01,0.2320158308543352E-01,
     *  0.2453316541321353E-01,0.2281380678501484E-01,
     *  0.1819076101277430E-01,0.1181569328630641E-01,
     *  0.5560112909407407E-02,0.1356395803463648E-02/
c 
c     Accuracy for the following 30 nodes is 0.56489E-07
c 
        data xs30/
     1  -.9968658077928980E+00,-.9762336538564690E+00,
     2  -.9246312444597489E+00,-.8353114985940511E+00,
     3  -.7078556818388058E+00,-.5468182425401615E+00,
     4  -.3605134470503768E+00,-.1598821606160766E+00,
     5  0.4264212731451811E-01,0.2342937975231325E+00,
     6  0.4031820061934088E+00,0.5396721755727467E+00,
     7  0.6383099880777727E+00,0.7003768377362947E+00,
     8  0.7348015402686319E+00,0.7541559713505426E+00,
     9  0.7680531902293495E+00,0.7803287134480754E+00,
     *  0.7918487304099737E+00,0.8028054566930549E+00,
     *  0.8132957819003213E+00,0.8235187025002095E+00,
     *  0.8342512108587818E+00,0.8480819480237680E+00,
     *  0.8692925020863114E+00,0.8995695226531633E+00,
     *  0.9346547040275897E+00,0.9666038812141811E+00,
     *  0.9885880412387591E+00,0.9984161755106855E+00/
c 
        data whts30/
     1  0.9210068626473613E-02,0.3439942229594591E-01,
     2  0.6991053229308355E-01,0.1087753508963304E+00,
     3  0.1453573123242342E+00,0.1753131556956549E+00,
     4  0.1954540941397354E+00,0.2037168828553824E+00,
     5  0.1991842250608646E+00,0.1821152085883932E+00,
     6  0.1540246554648152E+00,0.1179992233082189E+00,
     7  0.7946141118152804E-01,0.4623609495097078E-01,
     8  0.2483312509328385E-01,0.1553970014655159E-01,
     9  0.1281698242619664E-01,0.1184408011179729E-01,
     *  0.1122153487567408E-01,0.1070584553249380E-01,
     *  0.1030199109763159E-01,0.1025002035080787E-01,
     *  0.1163176931883284E-01,0.1683522638216555E-01,
     *  0.2593734168454342E-01,0.3385707350157659E-01,
     *  0.3489813945693877E-01,0.2780534023584837E-01,
     *  0.1576652859658623E-01,0.4597656221771766E-02/
c 
        if(n .eq. 44) call o1rarrcp(xs44,xs,n)
        if(n .eq. 44) call o1rarrcp(whts44,whts,n)
c 
        if(n .eq. 19) call o1rarrcp(xs19,xs,n)
        if(n .eq. 19) call o1rarrcp(whts19,whts,n)
c 
        if(n .eq. 30) call o1rarrcp(xs30,xs,n)
        if(n .eq. 30) call o1rarrcp(whts30,whts,n)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o1rnear21(n,xs,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),whts(1),
     1      xs44(44),whts44(44),xs30(30),whts30(30),
     2      xs19(19),whts19(19)
c 
c     Accuracy for the following 19 nodes is 0.43306E-04
c 
        data xs19/
     1  -.9886491194165523E+00,-.9367465022733261E+00,
     2  -.8372867074026312E+00,-.6916535867888491E+00,
     3  -.5074080482038008E+00,-.2963248321028831E+00,
     4  -.7306927564278620E-01,0.1459235588398728E+00,
     5  0.3432143510876580E+00,0.5005547773681949E+00,
     6  0.6020375403448768E+00,0.6556969646674000E+00,
     7  0.6938804045190199E+00,0.7284643677130049E+00,
     8  0.7680359807558252E+00,0.8294833455730377E+00,
     9  0.9014663492958592E+00,0.9601162278687159E+00,
     *  0.9928175579195608E+00/
c 
        data whts19/
     1  0.2965922387208696E-01,0.7526292563903619E-01,
     2  0.1233759549111330E+00,0.1665974994955693E+00,
     3  0.1999068931118353E+00,0.2197921015405785E+00,
     4  0.2239631218289788E+00,0.2111186343502881E+00,
     5  0.1804122484086081E+00,0.1312765787801423E+00,
     6  0.7252342676867754E-01,0.4175702764634547E-01,
     7  0.3587139546967119E-01,0.3424901494706096E-01,
     8  0.4908095684200398E-01,0.7111000898291256E-01,
     9  0.6854119264793121E-01,0.4667012679879992E-01,
     *  0.1882105137169266E-01/
c 
c     Accuracy for the following 44 nodes is 0.15305E-10
c 
        data xs44/
     1  -.9988020861813308E+00,-.9898571753636181E+00,
     2  -.9642867829669701E+00,-.9147273703997425E+00,
     3  -.8371419778303326E+00,-.7311375058386664E+00,
     4  -.5996510158709474E+00,-.4483690084364883E+00,
     5  -.2849950430731466E+00,-.1184058721837578E+00,
     6  0.4227972665751256E-01,0.1886990243190942E+00,
     7  0.3142798660350397E+00,0.4152116766003329E+00,
     8  0.4910735518436242E+00,0.5446779332514495E+00,
     9  0.5809412793094231E+00,0.6052958735278480E+00,
     *  0.6224753987152758E+00,0.6359164598242462E+00,
     *  0.6476079324484951E+00,0.6584586514952950E+00,
     *  0.6688311183129476E+00,0.6788620715915439E+00,
     *  0.6885938857861349E+00,0.6980250652431981E+00,
     *  0.7071388070068161E+00,0.7159313510086636E+00,
     *  0.7244391701560211E+00,0.7327591501347979E+00,
     *  0.7410838492796978E+00,0.7498129706361600E+00,
     *  0.7597775506875977E+00,0.7723794229822663E+00,
     *  0.7892828366207720E+00,0.8117765890573818E+00,
     *  0.8400860421842884E+00,0.8728657826336672E+00,
     *  0.9072582283710082E+00,0.9395789740801375E+00,
     *  0.9663218934829174E+00,0.9851427568016000E+00,
     *  0.9956051932693386E+00,0.9994662174362565E+00/
c 
        data whts44/
     1  0.3640696428084237E-02,0.1580713179423444E-01,
     2  0.3658771241186716E-01,0.6320718919134595E-01,
     3  0.9202315050218600E-01,0.1194863943064186E+00,
     4  0.1425265515393822E+00,0.1587394166205011E+00,
     5  0.1665134895472645E+00,0.1651283695775322E+00,
     6  0.1548316142351787E+00,0.1368948004342463E+00,
     7  0.1136200851362542E+00,0.8817525554327417E-01,
     8  0.6403746089522676E-01,0.4401967194709616E-01,
     9  0.2943322368789237E-01,0.2007370942299389E-01,
     *  0.1486172624008938E-01,0.1234431941370883E-01,
     *  0.1118018488474217E-01,0.1057660097728426E-01,
     *  0.1018921425533846E-01,0.9879040801410185E-02,
     *  0.9583909695355134E-02,0.9275170008945812E-02,
     *  0.8951200732432675E-02,0.8639678000786922E-02,
     *  0.8391679583202163E-02,0.8278695663988492E-02,
     *  0.8433216809603397E-02,0.9159783664474504E-02,
     *  0.1100880125366829E-01,0.1448369410269420E-01,
     *  0.1954968809447005E-01,0.2548514693333851E-01,
     *  0.3090864733087710E-01,0.3416224119489661E-01,
     *  0.3399095284193243E-01,0.3004933486489155E-01,
     *  0.2303298170592834E-01,0.1453008610799662E-01,
     *  0.6696858913311461E-02,0.1611222700692569E-02/
c 
c     Accuracy for the following 30 nodes is 0.56489E-07
c 
        data xs30/
     1  -.9969956935609792E+00,-.9771544662016288E+00,
     2  -.9274023982932004E+00,-.8411744188571115E+00,
     3  -.7181106381457538E+00,-.5627585788882246E+00,
     4  -.3833987685413516E+00,-.1909546858561161E+00,
     5  0.2089935900951190E-02,0.1828097038640513E+00,
     6  0.3390917949760024E+00,0.4613867112472509E+00,
     7  0.5456138349107091E+00,0.5964794998480409E+00,
     8  0.6264351511239745E+00,0.6476722579439389E+00,
     9  0.6660228983511731E+00,0.6830962710425181E+00,
     *  0.6992869368368375E+00,0.7147483930956320E+00,
     *  0.7296786527310655E+00,0.7447495835356062E+00,
     *  0.7624495244495859E+00,0.7880054658220872E+00,
     *  0.8253298682420731E+00,0.8719721565429648E+00,
     *  0.9199524673964641E+00,0.9602745319045747E+00,
     *  0.9866961660487654E+00,0.9981761241241478E+00/
c 
        data whts30/
     1  0.8836966965709429E-02,0.3312373438849103E-01,
     2  0.6745698209021021E-01,0.1050364161648781E+00,
     3  0.1403139292604664E+00,0.1689921132309771E+00,
     4  0.1878889488216523E+00,0.1949000963555840E+00,
     5  0.1890147608622727E+00,0.1703788161901379E+00,
     6  0.1405457601827447E+00,0.1032949764618934E+00,
     7  0.6594790458526786E-01,0.3801888467377295E-01,
     8  0.2405492683958094E-01,0.1931824430900169E-01,
     9  0.1760090987101549E-01,0.1659951661524327E-01,
     *  0.1580405331365424E-01,0.1514671746432947E-01,
     *  0.1480433931361098E-01,0.1571465506702384E-01,
     *  0.2064011029103480E-01,0.3123052057559597E-01,
     *  0.4298533018883541E-01,0.4888089246840376E-01,
     *  0.4551819869494300E-01,0.3404707545838514E-01,
     *  0.1859440474083782E-01,0.5309805191216267E-02/
c 
        if(n .eq. 44) call o1rarrcp(xs44,xs,n)
        if(n .eq. 44) call o1rarrcp(whts44,whts,n)
c 
        if(n .eq. 19) call o1rarrcp(xs19,xs,n)
        if(n .eq. 19) call o1rarrcp(whts19,whts,n)
c 
        if(n .eq. 30) call o1rarrcp(xs30,xs,n)
        if(n .eq. 30) call o1rarrcp(whts30,whts,n)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o1rnear22(n,xs,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),whts(1),
     1      xs45(45),whts45(45),xs31(31),whts31(31),
     2      xs19(19),whts19(19)
c 
c     Accuracy for the following 19 nodes is 0.45814E-04
c 
        data xs19/
     1  -.9892324515414290E+00,-.9395612171499370E+00,
     2  -.8437432793489460E+00,-.7031078123581609E+00,
     3  -.5252299838925633E+00,-.3220450948385811E+00,
     4  -.1087196145094177E+00,0.9699527622668636E-01,
     5  0.2745916007798634E+00,0.4021631429628556E+00,
     6  0.4770774334431263E+00,0.5314158382666900E+00,
     7  0.5803061123370192E+00,0.6290389781969077E+00,
     8  0.6991267829626693E+00,0.7933166992620844E+00,
     9  0.8853363574231774E+00,0.9544861686801753E+00,
     *  0.9918466781149764E+00/
c 
        data whts19/
     1  0.2820727624954179E-01,0.7230705810002386E-01,
     2  0.1190529938202133E+00,0.1609251747756683E+00,
     3  0.1928201775173677E+00,0.2110069970781616E+00,
     4  0.2126892802090547E+00,0.1953366198452891E+00,
     5  0.1559647780836748E+00,0.9777631550646480E-01,
     6  0.5940479251282965E-01,0.5104273944805404E-01,
     7  0.4709242582834831E-01,0.5489640434753981E-01,
     8  0.8573101460965485E-01,0.9756495199738591E-01,
     9  0.8302004382520757E-01,0.5378738782845856E-01,
     *  0.2136793520492451E-01/
c 
c     Accuracy for the following 45 nodes is 0.68814E-11
c 
        data xs45/
     1  -.9989625745727619E+00,-.9911576492785744E+00,
     2  -.9686464547695826E+00,-.9246671626942750E+00,
     3  -.8553850442811437E+00,-.7602989432990876E+00,
     4  -.6420205144881472E+00,-.5057689275163770E+00,
     5  -.3587008526219508E+00,-.2091138288228327E+00,
     6  -.6553402265472330E-01,0.6430855079383080E-01,
     7  0.1745797915320115E+00,0.2623353878716486E+00,
     8  0.3279954430588827E+00,0.3748778257886216E+00,
     9  0.4078848682527689E+00,0.4320061234516839E+00,
     *  0.4512041553170076E+00,0.4679476005132304E+00,
     *  0.4834747307776357E+00,0.4983367969273780E+00,
     *  0.5127732946036421E+00,0.5268757105742530E+00,
     *  0.5406384973834215E+00,0.5539927822213995E+00,
     *  0.5668638356152351E+00,0.5792374964515672E+00,
     *  0.5912012814678216E+00,0.6029527606223151E+00,
     *  0.6148334979411780E+00,0.6274833295894950E+00,
     *  0.6421187299456301E+00,0.6606216241961327E+00,
     *  0.6850465950224689E+00,0.7167774831332085E+00,
     *  0.7558066591367032E+00,0.8004032561710779E+00,
     *  0.8473498381522355E+00,0.8926399321526517E+00,
     *  0.9323367368297450E+00,0.9633800017184970E+00,
     *  0.9842538036517137E+00,0.9954404363916625E+00,
     *  0.9994548261762761E+00/
c 
        data whts45/
     1  0.3158794696445065E-02,0.1384329481260523E-01,
     2  0.3233278991141211E-01,0.5626993161050674E-01,
     3  0.8237342957569140E-01,0.1073548455442453E+00,
     4  0.1283209099088607E+00,0.1429761858779525E+00,
     5  0.1497618416354910E+00,0.1479740549648451E+00,
     6  0.1378785248405060E+00,0.1208166514744227E+00,
     7  0.9922590502701004E-01,0.7636324141251029E-01,
     8  0.5553371229766527E-01,0.3907453878591492E-01,
     9  0.2779001989087585E-01,0.2112142675243445E-01,
     *  0.1768010655676475E-01,0.1600290109478750E-01,
     *  0.1513707400742450E-01,0.1462450100562789E-01,
     *  0.1426315003043269E-01,0.1394028916460660E-01,
     *  0.1357279586480262E-01,0.1312211427967135E-01,
     *  0.1261690275178575E-01,0.1214502313582313E-01,
     *  0.1181516831649431E-01,0.1174228983611414E-01,
     *  0.1212042235411597E-01,0.1337901500922706E-01,
     *  0.1621455645420576E-01,0.2114907821179425E-01,
     *  0.2794393651218366E-01,0.3552274896245246E-01,
     *  0.4224283654349604E-01,0.4640162958725748E-01,
     *  0.4680897438549348E-01,0.4310147801079933E-01,
     *  0.3577083597941778E-01,0.2605505357239536E-01,
     *  0.1577124076081942E-01,0.7033620817754581E-02,
     *  0.1652157773896098E-02/
c 
c     Accuracy for the following 31 nodes is 0.30522E-07
c 
        data xs31/
     1  -.9973943878570442E+00,-.9798588668758120E+00,
     2  -.9351388300452994E+00,-.8567220184121042E+00,
     3  -.7439700166179286E+00,-.6010379434045056E+00,
     4  -.4358154962733537E+00,-.2589297693497937E+00,
     5  -.8273025692571989E-01,0.7986040015178153E-01,
     6  0.2169426793916551E+00,0.3203483696034382E+00,
     7  0.3895360855434817E+00,0.4332390647861921E+00,
     8  0.4639871928093120E+00,0.4899139316788479E+00,
     9  0.5138589316088624E+00,0.5365659089830277E+00,
     *  0.5582657072634627E+00,0.5791074080573608E+00,
     *  0.5993827317190608E+00,0.6201198020988843E+00,
     *  0.6447291798260630E+00,0.6795217689568569E+00,
     *  0.7287541830218600E+00,0.7898391494259449E+00,
     *  0.8545489263477521E+00,0.9133321428520280E+00,
     *  0.9586173785334708E+00,0.9865523255959583E+00,
     *  0.9981967738369535E+00/
c 
        data whts31/
     1  0.7706199001203608E-02,0.2950947056900792E-01,
     2  0.6100613941936163E-01,0.9591448085575623E-01,
     3  0.1288794834759976E+00,0.1556532872072299E+00,
     4  0.1730037168127016E+00,0.1786909587147210E+00,
     5  0.1715294231451024E+00,0.1516252983344496E+00,
     6  0.1211051978407105E+00,0.8557001964804584E-01,
     7  0.5436192850482298E-01,0.3531007380495960E-01,
     8  0.2753532247618236E-01,0.2472365932725372E-01,
     9  0.2326848924401927E-01,0.2217867633556779E-01,
     *  0.2124302933526254E-01,0.2048141084724312E-01,
     *  0.2021187637715499E-01,0.2179789060535188E-01,
     *  0.2858053835746929E-01,0.4180986969161327E-01,
     *  0.5619238274341202E-01,0.6452229982603856E-01,
     *  0.6325067807323347E-01,0.5302806376013641E-01,
     *  0.3688152775967361E-01,0.1914745562595466E-01,
     *  0.5281147869522975E-02/
c 
        if(n .eq. 45) call o1rarrcp(xs45,xs,n)
        if(n .eq. 45) call o1rarrcp(whts45,whts,n)
c 
        if(n .eq. 19) call o1rarrcp(xs19,xs,n)
        if(n .eq. 19) call o1rarrcp(whts19,whts,n)
c 
        if(n .eq. 31) call o1rarrcp(xs31,xs,n)
        if(n .eq. 31) call o1rarrcp(whts31,whts,n)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o1rnear23(n,xs,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),whts(1),
     1      xs43(43),whts43(43),xs31(31),whts31(31),
     2      xs20(20),whts20(20)
c 
c     Accuracy for the following 20 nodes is 0.23010E-04
c 
        data xs20/
     1  -.9914969371648314E+00,-.9498785052492272E+00,
     2  -.8661894935045557E+00,-.7409793099811308E+00,
     3  -.5812386848976108E+00,-.3987217568142657E+00,
     4  -.2093205799148439E+00,-.3327727265021796E-01,
     5  0.1058694984640358E+00,0.1988570351738888E+00,
     6  0.2696207742338550E+00,0.3335311447544722E+00,
     7  0.3935611365306859E+00,0.4571338682413826E+00,
     8  0.5455269982986197E+00,0.6601381585847244E+00,
     9  0.7785370912218621E+00,0.8804829805482861E+00,
     *  0.9531847170836978E+00,0.9916486304111969E+00/
c 
        data whts20/
     1  0.2268502204578694E-01,0.6205492066061803E-01,
     2  0.1051614438152380E+00,0.1440703841506077E+00,
     3  0.1734394190366917E+00,0.1889494987169822E+00,
     4  0.1864890450501367E+00,0.1614723871510118E+00,
     5  0.1144182573929132E+00,0.7708763891395419E-01,
     6  0.6657456363149843E-01,0.6154608094037915E-01,
     7  0.5937416709531415E-01,0.7233376805304292E-01,
     8  0.1045143473155380E+00,0.1204623579534492E+00,
     9  0.1129336876457344E+00,0.8883074897818477E-01,
     *  0.5570237311577125E-01,0.2189709940507938E-01/
c 
c     Accuracy for the following 43 nodes is 0.17809E-10
c 
        data xs43/
     1  -.9989924228810084E+00,-.9914516000698774E+00,
     2  -.9698596473412585E+00,-.9280285034022416E+00,
     3  -.8627480943584006E+00,-.7741239541814825E+00,
     4  -.6653264195040884E+00,-.5420738836297294E+00,
     5  -.4119394877095165E+00,-.2834681983571969E+00,
     6  -.1650392900335569E+00,-.6343460957091251E-01,
     7  0.1767769587471335E-01,0.7864948814168771E-01,
     8  0.1233335990224016E+00,0.1571672259508894E+00,
     9  0.1849563321694447E+00,0.2097174282204211E+00,
     *  0.2329066744499093E+00,0.2551115500245793E+00,
     *  0.2765274896712056E+00,0.2972035121127596E+00,
     *  0.3171913679377849E+00,0.3366195767589673E+00,
     *  0.3556420659657904E+00,0.3743487175171503E+00,
     *  0.3927983253550122E+00,0.4111763726382400E+00,
     *  0.4300405175264663E+00,0.4506908861883125E+00,
     *  0.4754755872753678E+00,0.5074168531671073E+00,
     *  0.5490261423323439E+00,0.6011050882713356E+00,
     *  0.6621733691043972E+00,0.7286777603989769E+00,
     *  0.7957721130137372E+00,0.8582858703375205E+00,
     *  0.9116239259982948E+00,0.9525242403856780E+00,
     *  0.9796714473167702E+00,0.9941199420490006E+00,
     *  0.9992960838269744E+00/
c 
        data whts43/
     1  0.3064163731622117E-02,0.1333807457652157E-01,
     2  0.3090202734427960E-01,0.5329915268480300E-01,
     3  0.7723922628317888E-01,0.9947188903413640E-01,
     4  0.1171586026167826E+00,0.1280697373686990E+00,
     5  0.1307561185620206E+00,0.1247657054507977E+00,
     6  0.1109347714792956E+00,0.9164563475859570E-01,
     7  0.7064623470067902E-01,0.5197444798597631E-01,
     8  0.3834067092801446E-01,0.3014803649969441E-01,
     9  0.2592579161591973E-01,0.2382597291832825E-01,
     *  0.2264273015176203E-01,0.2179590382792159E-01,
     *  0.2104170525177809E-01,0.2031746952345040E-01,
     *  0.1968075595847454E-01,0.1920336206088538E-01,
     *  0.1885594726682345E-01,0.1856435607042514E-01,
     *  0.1836191398034234E-01,0.1848072714664014E-01,
     *  0.1945720757664366E-01,0.2224034897335026E-01,
     *  0.2785530724954690E-01,0.3646721097300964E-01,
     *  0.4689487281237557E-01,0.5701223546424012E-01,
     *  0.6451804949228709E-01,0.6766628180370328E-01,
     *  0.6564555094084440E-01,0.5860352232697520E-01,
     *  0.4751862829802698E-01,0.3405525912764510E-01,
     *  0.2041347954128076E-01,0.9068913928794890E-02,
     *  0.2132001711110899E-02/
c 
c     Accuracy for the following 31 nodes is 0.25681E-07
c 
        data xs31/
     1  -.9976850013477988E+00,-.9819617288263989E+00,
     2  -.9415505625162428E+00,-.8703996449344197E+00,
     3  -.7680656951449859E+00,-.6388346175236501E+00,
     4  -.4908164733860206E+00,-.3350825303184389E+00,
     5  -.1846643090999665E+00,-.5289956162605469E-01,
     6  0.4995595054080981E-01,0.1220595699145581E+00,
     7  0.1721227092537763E+00,0.2117779547784717E+00,
     8  0.2473521229280662E+00,0.2808501511608590E+00,
     9  0.3128304707307397E+00,0.3435013449387836E+00,
     *  0.3730457423547122E+00,0.4018845554736854E+00,
     *  0.4314575708205501E+00,0.4661947738997880E+00,
     *  0.5137937217600606E+00,0.5791922665789753E+00,
     *  0.6597024273218304E+00,0.7468264791342076E+00,
     *  0.8303997331823110E+00,0.9014685451713598E+00,
     *  0.9538556977356984E+00,0.9852127893042546E+00,
     *  0.9980332586676575E+00/
c 
        data whts31/
     1  0.6864901760818050E-02,0.2656102895298095E-01,
     2  0.5526519664348753E-01,0.8709424437988358E-01,
     3  0.1168443061270160E+00,0.1402521569508702E+00,
     4  0.1539232738442097E+00,0.1553479358001449E+00,
     5  0.1432120254997078E+00,0.1185078141732978E+00,
     6  0.8683208343633279E-01,0.5898256863785917E-01,
     7  0.4326082307800091E-01,0.3704879068466773E-01,
     8  0.3438286307445142E-01,0.3269172742733316E-01,
     9  0.3129953367013553E-01,0.3006939206360375E-01,
     *  0.2907919735105622E-01,0.2880260480699551E-01,
     *  0.3104986131444356E-01,0.3984168164513571E-01,
     *  0.5627399885832799E-01,0.7402803822407240E-01,
     *  0.8549808502727474E-01,0.8700951680964443E-01,
     *  0.7862076107427637E-01,0.6247092516896409E-01,
     *  0.4189498486082988E-01,0.2121954083503354E-01,
     *  0.5770134209458634E-02/
c 
        if(n .eq. 43) call o1rarrcp(xs43,xs,n)
        if(n .eq. 43) call o1rarrcp(whts43,whts,n)
c 
        if(n .eq. 20) call o1rarrcp(xs20,xs,n)
        if(n .eq. 20) call o1rarrcp(whts20,whts,n)
c 
        if(n .eq. 31) call o1rarrcp(xs31,xs,n)
        if(n .eq. 31) call o1rarrcp(whts31,whts,n)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o1rnear24(n,xs,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),whts(1),
     1      xs43(43),whts43(43),xs29(29),whts29(29),
     2      xs20(20),whts20(20)
c 
c     Accuracy for the following 20 nodes is 0.11337E-04
c 
        data xs20/
     1  -.9931227609653135E+00,-.9581816705065223E+00,
     2  -.8865425440412869E+00,-.7791927984605388E+00,
     3  -.6440685582024894E+00,-.4952210918751745E+00,
     4  -.3524666905543219E+00,-.2342508687335489E+00,
     5  -.1380188697628004E+00,-.5093593221984460E-01,
     6  0.3103990962536796E-01,0.1100538687254136E+00,
     7  0.1976530629770206E+00,0.3156232466229993E+00,
     8  0.4608824235343234E+00,0.6129078781085833E+00,
     9  0.7530370597138361E+00,0.8675359577716115E+00,
     *  0.9476463300565834E+00,0.9904457982896060E+00/
c 
        data whts20/
     1  0.1858096169799480E-01,0.5275915745237006E-01,
     2  0.9029546518876053E-01,0.1230719132002155E+00,
     3  0.1448389975293389E+00,0.1494894204598824E+00,
     4  0.1324532147533561E+00,0.1046567327310067E+00,
     5  0.9036548321003829E-01,0.8428330389442056E-01,
     6  0.7988677598820436E-01,0.7971332659811733E-01,
     7  0.1003051591216060E+00,0.1346342182019034E+00,
     8  0.1521521783628742E+00,0.1487760621818772E+00,
     9  0.1291972918408494E+00,0.9832597107633219E-01,
     *  0.6138358701472247E-01,0.2483080422647526E-01/
c 
c     Accuracy for the following 43 nodes is 0.48038E-11
c 
        data xs43/
     1  -.9991762644044248E+00,-.9930104416979644E+00,
     2  -.9753916285189789E+00,-.9414595322163966E+00,
     3  -.8890857805861597E+00,-.8192057408990638E+00,
     4  -.7355736573363845E+00,-.6441765814134084E+00,
     5  -.5522913979165238E+00,-.4670698108890621E+00,
     6  -.3937209530474515E+00,-.3339731934000060E+00,
     7  -.2859756903243095E+00,-.2459826340115676E+00,
     8  -.2104786683106170E+00,-.1772069850851874E+00,
     9  -.1451017446732662E+00,-.1138329564767681E+00,
     *  -.8338460221368127E-01,-.5380362340716316E-01,
     *  -.2507147692718192E-01,0.2906556860765575E-02,
     *  0.3020629294209942E-01,0.5684025755278110E-01,
     *  0.8287531932929748E-01,0.1086246213929654E+00,
     *  0.1347573486749076E+00,0.1625222896282110E+00,
     *  0.1941899645623675E+00,0.2330205524382063E+00,
     *  0.2822308403142797E+00,0.3437026985027247E+00,
     *  0.4172409551960718E+00,0.5005016524929142E+00,
     *  0.5894494857461967E+00,0.6790647551379518E+00,
     *  0.7640621253565583E+00,0.8395466553583051E+00,
     *  0.9016094685455716E+00,0.9478651659396709E+00,
     *  0.9779143731910300E+00,0.9936609814401445E+00,
     *  0.9992448024785774E+00/
c 
        data whts43/
     1  0.2505514494681699E-02,0.1090233610986151E-01,
     2  0.2516810271770025E-01,0.4305491632468746E-01,
     3  0.6153211064267487E-01,0.7759257979632185E-01,
     4  0.8865840172531551E-01,0.9289940923024194E-01,
     5  0.8964714498751201E-01,0.7989877856277230E-01,
     6  0.6650784263175766E-01,0.5332817099380670E-01,
     7  0.4333136693525287E-01,0.3725697886149845E-01,
     8  0.3412387313313280E-01,0.3258867222014662E-01,
     9  0.3166941579227870E-01,0.3086693234669855E-01,
     *  0.3001978181505174E-01,0.2914498459985405E-01,
     *  0.2833721903166269E-01,0.2763250780469105E-01,
     *  0.2696692216543984E-01,0.2630821592019902E-01,
     *  0.2581297692539384E-01,0.2579419453673887E-01,
     *  0.2666989637392571E-01,0.2923678859947042E-01,
     *  0.3466191852467579E-01,0.4355908050576898E-01,
     *  0.5518744383585182E-01,0.6772467454444249E-01,
     *  0.7895276408225893E-01,0.8688367574797642E-01,
     *  0.9016327960819689E-01,0.8817547131255120E-01,
     *  0.8099352298308712E-01,0.6931507093739443E-01,
     *  0.5440515470316151E-01,0.3803391172809507E-01,
     *  0.2238157049971335E-01,0.9816551260598449E-02,
     *  0.2289874446783353E-02/
c 
c     Accuracy for the following 29 nodes is 0.38451E-07
c 
        data xs29/
     1  -.9976854391759024E+00,-.9823232366007102E+00,
     2  -.9438765457043543E+00,-.8780337242040053E+00,
     3  -.7862931686709251E+00,-.6751130891212024E+00,
     4  -.5550832464007436E+00,-.4395901692224016E+00,
     5  -.3411915250420908E+00,-.2643601544120736E+00,
     6  -.2025664584606489E+00,-.1477170249353966E+00,
     7  -.9618777868183270E-01,-.4689985291449106E-01,
     8  0.4830384831349122E-03,0.4616693966405516E-01,
     9  0.9055585474588601E-01,0.1350190465451493E+00,
     *  0.1841739933000324E+00,0.2478062803980092E+00,
     *  0.3339258112837094E+00,0.4407342722660650E+00,
     *  0.5590864122286330E+00,0.6776986838176289E+00,
     *  0.7860579273414019E+00,0.8757363411585322E+00,
     *  0.9412158557217214E+00,0.9807220437166020E+00,
     *  0.9973389611178993E+00/
c 
        data whts29/
     1  0.6821721813412051E-02,0.2566422479482759E-01,
     2  0.5193847818349775E-01,0.7948228212661382E-01,
     3  0.1029097387088860E+00,0.1176673043600479E+00,
     4  0.1201085195591501E+00,0.1086475464298751E+00,
     5  0.8730670813279295E-01,0.6767813356248281E-01,
     6  0.5735664553285077E-01,0.5288741369479714E-01,
     7  0.5032303736185467E-01,0.4830093583604944E-01,
     8  0.4649502958219283E-01,0.4493016563771814E-01,
     9  0.4403991799238078E-01,0.4558692546813350E-01,
     *  0.5449497392827882E-01,0.7428377918472702E-01,
     *  0.9757649149360294E-01,0.1144206683741766E+00,
     *  0.1203585059163148E+00,0.1150945956479284E+00,
     *  0.1002064772071606E+00,0.7821740954791498E-01,
     *  0.5242749923380588E-01,0.2707474440151640E-01,
     *  0.7700123861867440E-02/
c 
        if(n .eq. 43) call o1rarrcp(xs43,xs,n)
        if(n .eq. 43) call o1rarrcp(whts43,whts,n)
c 
        if(n .eq. 20) call o1rarrcp(xs20,xs,n)
        if(n .eq. 20) call o1rarrcp(whts20,whts,n)
c 
        if(n .eq. 29) call o1rarrcp(xs29,xs,n)
        if(n .eq. 29) call o1rarrcp(whts29,whts,n)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o1rnear25(n,xs,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),whts(1),
     1      xs39(39),whts39(39),xs27(27),whts27(27),
     2      xs18(18),whts18(18)
c 
c     Accuracy for the following 17 nodes is 0.28190E-04
c 
        data xs18/
     1  -.9933969664526930E+00,-.9616885411231008E+00,
     2  -.9004256769292588E+00,-.8144649103402517E+00,
     3  -.7136250446553867E+00,-.6067060472540768E+00,
     4  -.4982709912897290E+00,-.3899690537914176E+00,
     5  -.2780436163129648E+00,-.1452052065586890E+00,
     6  0.2000939421171033E-01,0.2072410619126013E+00,
     7  0.3998321715408963E+00,0.5821154401349111E+00,
     8  0.7406169512111775E+00,0.8648773187843877E+00,
     9  0.9485299478235873E+00,0.9910184203678949E+00/
c 
        data whts18/
     1  0.1754819826625592E-01,0.4662095639585693E-01,
     2  0.7501363555983186E-01,0.9516306217911734E-01,
     3  0.1049792829552219E+00,0.1081579318741594E+00,
     4  0.1084069804298545E+00,0.1085681400037220E+00,
     5  0.1185934660892361E+00,0.1493518856442946E+00,
     6  0.1789351893888686E+00,0.1926724078283414E+00,
     7  0.1898881905047069E+00,0.1724190109485667E+00,
     8  0.1428229428202278E+00,0.1046121626172481E+00,
     9  0.6254376382260343E-01,0.2370277467789050E-01/
c 
c     Accuracy for the following 39 nodes is 0.27180E-11
c 
        data xs39/
     1  -.9994561577994434E+00,-.9954647313742706E+00,
     2  -.9844484590560405E+00,-.9642841448119769E+00,
     3  -.9351702404879451E+00,-.8992154426697946E+00,
     4  -.8592180543580756E+00,-.8174207367511310E+00,
     5  -.7750791565131189E+00,-.7327738950355271E+00,
     6  -.6907733355265636E+00,-.6492150039948181E+00,
     7  -.6081776303686119E+00,-.5677150644606777E+00,
     8  -.5278951393853835E+00,-.4887688614484732E+00,
     9  -.4503293579032636E+00,-.4125029477254274E+00,
     *  -.3751435821480424E+00,-.3380136990937216E+00,
     *  -.3005830964575941E+00,-.2616804912693433E+00,
     *  -.2191328071020189E+00,-.1695715484573959E+00,
     *  -.1088771949413494E+00,-.3383604084265246E-01,
     *  0.5634679013548061E-01,0.1600441578547713E+00,
     *  0.2738645932858692E+00,0.3932348124421783E+00,
     *  0.5128927445965037E+00,0.6273725188703574E+00,
     *  0.7314875704002217E+00,0.8208037165172068E+00,
     *  0.8921075056466639E+00,0.9438695083113322E+00,
     *  0.9766649986655163E+00,0.9934282702279666E+00,
     *  0.9992302528002048E+00/
c 
        data whts39/
     1  0.1647341720305005E-02,0.6974838030267425E-02,
     2  0.1542514012755751E-01,0.2486258334174745E-01,
     3  0.3299156232133211E-01,0.3842823529854743E-01,
     4  0.4118589997174785E-01,0.4220704404585153E-01,
     5  0.4238565791219878E-01,0.4218285833959212E-01,
     6  0.4179645872126193E-01,0.4130735756374849E-01,
     7  0.4075952731114119E-01,0.4015346574844647E-01,
     8  0.3947705568618462E-01,0.3877675112943169E-01,
     9  0.3811415876336881E-01,0.3756352046793975E-01,
     *  0.3719111954878813E-01,0.3715095798776647E-01,
     *  0.3789898769043523E-01,0.4026454609494094E-01,
     *  0.4539637515550262E-01,0.5444495933572292E-01,
     *  0.6750937960005970E-01,0.8270145843404297E-01,
     *  0.9738043302802145E-01,0.1094396506109952E+00,
     *  0.1174310037852353E+00,0.1204243132549715E+00,
     *  0.1179724930263851E+00,0.1101147669293812E+00,
     *  0.9736851122149481E-01,0.8072041041947407E-01,
     *  0.6162499686644399E-01,0.4198997181747528E-01,
     *  0.2409230999334658E-01,0.1029893570401913E-01,
     *  0.2344962994855047E-02/
c 
c     Accuracy for the following 27 nodes is 0.20323E-07
c 
        data xs27/
     1  -.9982772465334109E+00,-.9871026813778605E+00,
     2  -.9603616155764975E+00,-.9176013413528963E+00,
     3  -.8629969020968398E+00,-.8020693631816798E+00,
     4  -.7386757014101658E+00,-.6746768637413109E+00,
     5  -.6109592570939315E+00,-.5480147845255724E+00,
     6  -.4861277950413450E+00,-.4254008050988866E+00,
     7  -.3655949608336955E+00,-.3054231048078454E+00,
     8  -.2404571068062469E+00,-.1612677406143972E+00,
     9  -.5904169565137207E-01,0.6706743737492760E-01,
     *  0.2110790709201720E+00,0.3639286494095073E+00,
     *  0.5157046775414362E+00,0.6568336891606497E+00,
     *  0.7789273639290581E+00,0.8755847318165149E+00,
     *  0.9432525198139541E+00,0.9822105463396082E+00,
     *  0.9976739384350113E+00/
c 
        data whts27/
     1  0.5054052483015096E-02,0.1838708184178852E-01,
     2  0.3514644584881721E-01,0.4961899809820591E-01,
     3  0.5860982246840375E-01,0.6260346917963926E-01,
     4  0.6389547190716333E-01,0.6396297106723305E-01,
     5  0.6339283997279506E-01,0.6244918268139544E-01,
     6  0.6130573131540742E-01,0.6017866254183803E-01,
     7  0.5961398914292023E-01,0.6140690043322124E-01,
     8  0.7020624007629166E-01,0.8980620559211035E-01,
     9  0.1147364008402044E+00,0.1364160123517336E+00,
     *  0.1500619893418370E+00,0.1539671841753172E+00,
     *  0.1479757545916568E+00,0.1328710321752230E+00,
     *  0.1102352995393035E+00,0.8247096074819288E-01,
     *  0.5289186531221959E-01,0.2587946814065509E-01,
     *  0.6855967869996587E-02/
c 
        if(n .eq. 39) call o1rarrcp(xs39,xs,n)
        if(n .eq. 39) call o1rarrcp(whts39,whts,n)
c 
        if(n .eq. 18) call o1rarrcp(xs18,xs,n)
        if(n .eq. 18) call o1rarrcp(whts18,whts,n)
c 
        if(n .eq. 27) call o1rarrcp(xs27,xs,n)
        if(n .eq. 27) call o1rarrcp(whts27,whts,n)
c 
        return
        end
