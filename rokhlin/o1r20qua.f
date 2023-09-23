         implicit real *8 (a-h,o-z)
         dimension points(24,20),weights(24,20),xs(100),
     1      whts(100),coefslrg(20,24,20),
     2      numspts(100),fs(100),fslrg(24,20),fs2(100),diffs(100)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERs
C 
        PRINT *, 'ENTER k'
        READ *,k
        CALL PRINf('k=*',k,1 )
c 
  
        nn=20
  
        call o1rmatbg(coefslrg,xs,whts,points,weights,numspts)
  
        call prin2('after coefslrg, xs=*',xs,nn)
        call prin2('after coefslrg, whts=*',whts,nn)
c 
c       test the obtained interpolation coefficients
c 
c        . . . create the values of the test function at the
c              support nodes
c 
        do 2200 i=1,nn
c 
        call funtest(xs(i),fs(i),k)
  
 2200 continue
  
        call prin2('fs as created*',fs,nn)
c 
c        evaluate the test function at the support nodes via interpolation
c 
        call o1rappbg(coefslrg,numspts,fs,fslrg)
  
        do 2600 i=1,nn
c 
        numpts=numspts(i)
        call prinf('i=*',i,1)
        call prin2('and values of f at corresponding nodes is*',
     1      fslrg(1,i),numpts)
c 
        do 2400 j=1,numpts
c 
        call funtest(points(j,i),fs2(j),k)
c 
        diffs(j)=fs2(j)-fslrg(j,i)
  
  
 2400 continue
  
        call prin2('and fs2=*',fs2,numpts)
        call prin2('and diffs=*',diffs,numpts)
  
  
 2600 continue
  
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine funtest(x,f,k)
        implicit real *8 (a-h,o-z)
        save
  
        f=x
  
        call legepol(x,k,pol,der)
  
        f=(x+1)*log(x+1) * pol
        return
        end
  
c 
c 
c 
c 
c 
        function fun2(point,xgauss,iord)
        implicit real *8 (a-h,o-z)
        save
  
        call legepol(point,iord,pol,der)
  
        fun2=abs(point-xgauss)**2
        fun2=log(abs(point-xgauss)) * pol
        fun2= pol
        return
        end
c 
c 
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code and the beginning of the
c        quadrature code proper
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        The user-callable subroutines in this file are: o1rwrite,
c        o1rread, o1rmatbg, o1rappbg, o1r20qua, o1rnodes, o1rinter,
c        o1runiv, o1rinmat. Following is a description of these
c        subroutines.
c 
c 
c  o1rwrite - constructs (by calling the subroutine o1rmatbg) and stores
c        on disk all of the parameters used for the construction of the
c        matrix totmat. More specifically, it constructs the support
c        nodes and corresponding weights, auxiliary quadrature nodes for
c        each of the support nodes, interpolation coefficients from the
c        support nodes of the auxiliary nodes, and the integer array
c        numspts whose i-th entry is the number of auxiliary quadrature
c        nodes for the support node number i.
c  o1rread - reads from disk the various types of data to be used in
c        the construction of the matrix totmat. The data to be read
c        have been (hopefully) stored previously by a call to the
c        subroutine o1rwrite (see)
c  o1rmatbg - returns to the user a set of 20 matrices. Each of these
c       matrices converts the values of a function at the support nodes
c       (20 of them things) into the values of that function at the
c       quadrature nodes corresponding to one support node. Since there
c       are 20 support nodes, the number of matrices returned by this
c       subroutine is also 20. On the other hand, there
c  o1rappbg - applies to the user-specified vector fs of length 20
c       (representing the values at the 20 support nodes of an
c       appropriately-behaved function f) the 20 interpolation matrices
c       produced by the subroutine o1rmatbg (see). The result is the
c       collection fslrs of (approximate) values  of the function
c       f at the nodes points (also produced by the subroutine o1rmatbg).
c       This subroutine is principally a tool for the debugging of the
c       subroutine o1rmatbg, since I am not aware of any other uses for
c       this subroutine.
c   01r20qua - returns to the user the nodes and weights
c        of the quadrature nodes corresponding to the 20 support nodes
c        returned by the subroutine o1rnodes (see). Different support
c        nodes have different numbers of quadrature nodes corresponding
c        to them. The number of quadratures nodes corresponding to the
c        m-th support node is returned in the array numspts; the
c        maximum number of quadrature nodes in 24, and the minimum
c        number is 20. The resulting quadratures are accurate to a
c        little better than 10 digits.
c 
c   o1rnodes - returns to the user the nodes and weights of
c        the 20-point Generalized Gaussian quadrature on the interval
c        [-1,1] for functions of the forms
c 
c        P_n(x),
c 
c        P_n (x) \cdot (1-x) \cdot log(1-x),
c 
c        P_n (x) \cdot  \cdot(1+x) log(1+x),
c                                                                    (1)
c        P_n (x) \cdot ( (1-x) \cdot log(1-x) )^2,
c 
c        P_n (x) \cdot ( \cdot(1+x) log(1+x)  )^2,
c 
c        with n=0,1,2,...,19. The accuracy of the quadrature is
c        close to 0.2E-9.
c 
c   o1rinmat - returns to the user the matrices a, b of a singular
c        function expansion and its inverse. Specifically, the i-th
c        column of b contains the values at the 20 support nodes of
c        the 20 singular functions of representing on the interval
c        [-1,1] all functions of the forms
c 
c        P_n(x),
c 
c        P_n (x) \cdot (1-x) \cdot log(1-x),
c                                                                    (2)
c        P_n (x) \cdot  \cdot(1+x) log(1+x),
c 
c        with n=0,1,2,...,9.  The accuracy of the obtained expansion
c        is close to 1.0E-7, making it the weakest link among the
c        subroutines in this file.
c   o1rinter - returns to the user the coefficients forminte of the
c        interpolation formula expressing the values of a function
c        at the user-specified point x on the interval [-1,1] via
c       that function's values at the nn nodes points.
c 
c   o1runiv - returns to the user the coefficients of the universal
c        expansions of the first 20 singular functions on the
c        interval representing the functions of the form
c 
c          f(x) = phi(x) + psi(x) * (1-x) * log(1-x) +
c                 eta(x) * (1-x) * log(1+x),                         (3)
c 
c        with phi, psi, eta polynomials of order 9. While each of the
c        functions is approximated to 14 digits or so, the resulting
c        20-function expansion has accuracy of about 7 digits on
c        functions of the form (3). Also, see subroutine o1rinmat.
C        IMPORTANT: THE UNIVERSAL EXPANSION USED IS THE ONE WITH
C        THE PARAMETER NDIGITS SET TO 16!!!!
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        subroutine o1rwrite(iw)
        implicit real *8 (a-h,o-z)
        save
        real *8 xs(20),whts(20)
        dimension numspts(20),weights(24,20),points(24,20),
     1      coefslrg(20,24,20)
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
       call o1rmatbg(coefslrg,xs,whts,points,weights,numspts)
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
        subroutine o1rread(ir,xs,whts,points,weights,
     1      numspts,coefslrg)
        implicit real *8 (a-h,o-z)
        save
        real *8 xs(20),whts(20)
        dimension numspts(20),weights(24,20),points(24,20),
     1      coefslrg(20,24,20)
c 
c        This subroutine reads from the FORTRAN unit ir the various
c        types of data to be used in the construction of the matrix
c        totmat. The data to be read have been (hopefully) stored
c        previously by a call to the subroutine o1rwrite (see)
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
        subroutine o1rappbg(coefslrg,numspts,fs,fslrg)
        implicit real *8 (a-h,o-z)
        save
        dimension coefslrg(20,24,20),numspts(100),fs(1),
     1      fslrg(24,20)
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
        nn=20
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
        subroutine o1rmatbg(coefslrg,xs,whts,points,weights,numspts)
        implicit real *8 (a-h,o-z)
        save
        dimension points(24,20),weights(24,20),xs(100),
     1      whts(100),dumpts(100),dumwhts(100),coefslrg(20,24,20),
     2      numspts(100),forminte(100)
c 
c       This subroutine returns to the user a set of 20 matrices.
c       Each of these matrices converts the values of a function at
c       the support nodes (20 of them things) into the values of that
c       function at the quadrature nodes corresponding to one support
c       node. Since there are 20 support nodes, the number of matrices
c       returned by this subroutine is also 20. On the other hand, there
c       are two anomalies associated with the matrices returned by this
c       subroutine:
c 
c   1. The matrices returned by this subroutine are actually adjoints
c       of the matrices converting the values of a function at the
c       support nodes into its values at the quadrature nodes
c   2. While all of the matrices are declared as 20 \times 24 - matrices,
c       the actual lengths of their rows vary from 20 to 24, and are
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
c       numbers of quadrature nodes correponding to all 20 support nodes.
c 
c 
c                 Input parameters: None
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
c 
c        . . . retrieve from the subroutine o1rnodes the support nodes
c              and corresponding weights
c 
        call o1rnodes(xs,whts,nn)
c 
c        retrieve from the subroutine o1r20qua the nodes and weights of
c        all quadratures
c 
        call o1r20qua(numspts,points,weights)
  
        call prinf('after o1r20qua, numspts=*',numspts,20)
c 
c       for each of the 20 sets of quadrature nodes, construct the matrix
c       of interpolation coefficients connecting the values of a function
c       at the support nodes to its values at the quadrature nodes
c       PLEASE NOTE THAT DIFFERENT quadrature formulae have different
c       of nodes
c 
cccc        nn=20
        do 2000 k=1,nn
        call prinf('in o1rmatbg, k=*',k,1)
  
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
        subroutine o1r20qua(numspts,points,weights)
        implicit real *8 (a-h,o-z)
        save
        dimension numspts(20),points(24,20),weights(24,20)
c 
        dimension pts1(20),wts1(20),pts2(20),wts2(20),
     1      pts3(21),wts3(21),pts4(21),wts4(21),
     2      pts5(22),wts5(22),pts6(23),wts6(23),
     3      pts7(23),wts7(23),pts8(24),wts8(24),
     4      pts9(24),wts9(24),pts10(24),wts10(24)
c 
c        this subroutine returns to the user the nodes and weights
c        of the quadrature nodes corresponding to the 20 support nodes
c        returned by the subroutine o1rnodes (see). Different support
c        nodes have different numbers of quadrature nodes corresponding
c        to them. The number of quadratures nodes corresponding to the
c        m-th support node is returned in the array numspts; the
c        maximum number of quadrature nodes in 24, and the minimum
c        number is 20. PLEASE NOTE THAT THE QUADRATURE NODES FOR THE
C        M-TH SUPPORT NODE ARE RETURNED IN THE M-TH COLUMN OF ARRAY
C        POINTS. THE ARRAY POINTS IS DIMENSIONED (24,20), SO THAT SOME
C        OF ITS COLUMNS HAVE EMPTY ELEMENTS. THE SAME (NATURALLY)
C        APPLIES TO THE ARRAY WEIGHTS, RETURNING THE QUADRATURE WEIGHTS.
c 
c        ALSO, PLEASE NOTE THAT THE ACCURACY OF THE QUADRATURE FORMULAE
C        PROVIDED BY THIS SUBROUTINE IS BETTER THAN 10 DIGITS, WHICH
C        TENDS TO BE AN OVERKILL, GIVEN THE MUCH LOWER ACCURACY PROVIDED
C        BY THE SUBROUTINE O1RINTER (SEE).
C 
c                   Input parameters:  NONE
c 
c                   Output parameters:
c 
c  numspts - a 20-element integer array telling the user how many
c        quadrature nodes correspond to each of the 20 support nodes;
c        numspts(i) is that number for the i-th support node
c  points - the 24 \times 20 matrix containing the quadrature nodes
c        corresponding to the support nodes
c 
c 
        data pts1/
     1  -.9997594587343869E+00,-.9993425970334612E+00,
     2  -.9973353642819323E+00,-.5005362326108683E+00,
     3  -.9856370582941715E+00,-.9537109361866664E+00,
     4  -.8924637131579809E+00,-.7965149427877445E+00,
     5  -.6648362253264097E+00,-.3102450721530261E+00,
     6  -.1032662630553444E+00,0.1094272002330686E+00,
     7  0.3162922954331309E+00,0.5063479572681993E+00,
     8  0.8014192832368719E+00,0.8966489814297456E+00,
     9  0.6702746107233239E+00,0.9878568738632250E+00,
     *  0.9569836828797704E+00,0.9985731906493557E+00/
c 
        data wts1/
     1  0.4773498914216118E-03,0.4043934402147371E-03,
     2  0.5147280037965785E-02,0.1786546785163674E+00,
     3  0.2006757470613778E-01,0.4533781288859885E-01,
     4  0.7806167426711716E-01,0.1140027384713751E+00,
     5  0.1488259976710177E+00,0.2003601847775352E+00,
     6  0.2117501258837568E+00,0.2116939400845513E+00,
     7  0.2001875982329175E+00,0.1783535785430337E+00,
     8  0.1133798056281155E+00,0.7724436118315975E-01,
     9  0.1483765358900122E+00,0.1899467517535772E-01,
     *  0.4433735592056417E-01,0.4342338787068727E-02/
c 
        data pts2/
     1  -.9993309635997905E+00,-.9964422983434834E+00,
     2  -.9853803930507740E+00,-.9936233992228889E+00,
     3  -.9573474632937162E+00,-.9010711962763854E+00,
     4  -.8106634760910542E+00,-.6843644467714681E+00,
     5  -.5246196826700636E+00,-.3375409929054231E+00,
     6  -.1320901518591598E+00,0.8090643136409867E-01,
     7  0.2898491311905617E+00,0.4835010773946603E+00,
     8  0.7884166554263241E+00,0.6521044575454369E+00,
     9  0.8886155696210565E+00,0.9530249891886226E+00,
     *  0.9865470996246114E+00,0.9983977986134084E+00/
c 
        data wts2/
     1  0.1838997486814874E-02,0.3597444772965425E-02,
     2  0.1664309190409963E-01,0.2514304343686053E-02,
     3  0.4087260079387449E-01,0.7268558642965973E-01,
     4  0.1084221661784128E+00,0.1437587113933084E+00,
     5  0.1746983739257376E+00,0.1979488163983982E+00,
     6  0.2111301274429153E+00,0.2129106791735409E+00,
     7  0.2030831565571410E+00,0.1825770450597571E+00,
     8  0.1185676726998664E+00,0.1534074253193022E+00,
     9  0.8187615414236322E-01,0.4776158469819803E-01,
     *  0.2084839941090560E-01,0.4857661874810474E-02/
c 
        data pts3/
     1  -.9993713631873570E+00,-.9952685931821200E+00,
     2  -.9870729485636869E+00,-.9716078834473798E+00,
     3  -.9785356827067643E+00,-.9454157501092378E+00,
     4  -.8885744309030498E+00,-.7967245678780031E+00,
     5  -.6688165185185952E+00,-.5077319854868715E+00,
     6  -.3198846198573493E+00,-.1144287974901833E+00,
     7  0.9772553020966715E-01,0.3050055528124435E+00,
     8  0.4963079221330243E+00,0.7954602002272548E+00,
     9  0.6621046013479924E+00,0.8929061022208991E+00,
     *  0.9551127425445551E+00,0.9872305622873273E+00,
     *  0.9984886284628468E+00/
c 
        data wts3/
     1  0.1880150132738251E-02,0.6462009314903994E-02,
     2  0.9213234624752547E-02,0.1278633570055729E-01,
     3  0.6499229778372345E-02,0.4060924117654783E-01,
     4  0.7384634385246447E-01,0.1100350108189301E+00,
     5  0.1452929247951488E+00,0.1757936522632518E+00,
     6  0.1983587505115858E+00,0.2107161987997118E+00,
     7  0.2116446237879656E+00,0.2010462827166012E+00,
     8  0.1799533660312028E+00,0.1156522518011451E+00,
     9  0.1504694053364112E+00,0.7934480549487519E-01,
     *  0.4593116156862917E-01,0.1987474011888192E-01,
     *  0.4590281371940734E-02/
c 
        data pts4/
     1  -.9989830658265199E+00,-.9925268998896564E+00,
     2  -.9782810272137861E+00,-.9414493999222260E+00,
     3  -.9588914569329073E+00,-.9296348620772918E+00,
     4  -.8926776185725663E+00,-.7032586503077963E+00,
     5  -.8174010751003103E+00,-.5527021299261105E+00,
     6  -.3716195652455425E+00,-.1688868142344255E+00,
     7  0.4450820839127592E-01,0.2565445415130488E+00,
     8  0.4553353192184088E+00,0.7732879624118044E+00,
     9  0.8795694360406875E+00,0.6303016326832730E+00,
     *  0.9487080581935318E+00,0.9851583135163247E+00,
     *  0.9982155356207619E+00/
c 
        data wts4/
     1  0.2982844142145211E-02,0.1036761858887144E-01,
     2  0.1767048648672483E-01,0.1300086068643967E-01,
     3  0.1985579901179103E-01,0.1873000520298717E-01,
     4  0.5584032998436065E-01,0.1330489134945731E+00,
     5  0.9483310507884855E-01,0.1670741941026386E+00,
     6  0.1935965712324753E+00,0.2100193781880644E+00,
     7  0.2147475316747013E+00,0.2073242698833268E+00,
     8  0.1884754259223146E+00,0.1250896782717267E+00,
     9  0.8738945639768830E-01,0.1600810748116556E+00,
     *  0.5162829311727035E-01,0.2284828528436816E-01,
     *  0.5395878420073263E-02/
c 
        data pts5/
     1  -.9993623306184745E+00,-.9950566463079172E+00,
     2  -.9826708900829275E+00,-.9582821481789586E+00,
     3  -.9239162029944230E+00,-.8888261486924922E+00,
     4  -.8662329533309645E+00,-.8411972735449340E+00,
     5  -.5201974834187503E+00,-.7748295460035534E+00,
     6  -.6664000457115757E+00,-.3430702944847769E+00,
     7  -.1442997152597473E+00,0.6503194310234459E-01,
     8  0.2729760750002820E+00,0.4678024635476686E+00,
     9  0.6391217038372761E+00,0.8827568959292468E+00,
     *  0.7789640255546483E+00,0.9501533451566914E+00,
     *  0.9856053057568108E+00,0.9982726445290494E+00/
c 
        data wts5/
     1  0.1867657951580925E-02,0.7483120618374113E-02,
     2  0.1807619066422299E-01,0.3037450920923417E-01,
     3  0.3675960394422788E-01,0.3121633370305890E-01,
     4  0.1455541058438834E-01,0.4440143567481645E-01,
     5  0.1630232864817481E+00,0.8791685432263295E-01,
     6  0.1282376709130577E+00,0.1896732150759352E+00,
     7  0.2060049497437557E+00,0.2106485924965315E+00,
     8  0.2032657365797770E+00,0.1846373584850715E+00,
     9  0.1566564375283680E+00,0.8526547914876201E-01,
     *  0.1222521400143592E+00,0.5026820335990875E-01,
     *  0.2218970517985592E-01,0.5226108363732298E-02/
c 
        data pts6/
     1  -.9992464965780247E+00,-.9933352301551054E+00,
     2  -.9760030155562696E+00,-.8432101062836872E+00,
     3  -.9432848111387815E+00,-.8965650866102012E+00,
     4  -.7953042258028584E+00,-.7665564776929332E+00,
     5  -.7396458265060345E+00,-.6704527773373730E+00,
     6  -.5590402601267311E+00,-.4114615772556809E+00,
     7  -.2359574051419483E+00,-.4257577448506625E-01,
     8  0.1574989855898264E+00,0.3528135229050561E+00,
     9  0.5326181210058640E+00,0.6878546461462419E+00,
     *  0.8120806158547039E+00,0.9022613759250246E+00,
     *  0.9593551426605199E+00,0.9885377751788341E+00,
     *  0.9986544383257266E+00/
c 
        data wts6/
     1  0.2324100819716395E-02,0.1062761609818308E-01,
     2  0.2475055216944430E-01,0.5291803913343549E-01,
     3  0.4046318716131999E-01,0.5172768433138529E-01,
     4  0.4059479325972984E-01,0.1757878266877239E-01,
     5  0.4653671271548307E-01,0.9112421293541946E-01,
     6  0.1306933650439559E+00,0.1630837812030391E+00,
     7  0.1862385234074576E+00,0.1986460043435122E+00,
     8  0.1995781943654980E+00,0.1892426060851714E+00,
     9  0.1688413432364495E+00,0.1405446005217545E+00,
     *  0.1073915175873433E+00,0.7312510388731847E-01,
     *  0.4193204976003831E-01,0.1794112185720176E-01,
     *  0.4096107418252414E-02/
c 
        data pts7/
     1  -.9895527148246533E+00,-.9987324343760389E+00,
     2  -.7811090070425863E+00,-.9645858045603130E+00,
     3  -.9195296639170323E+00,-.8558455116065711E+00,
     4  -.7083820181630259E+00,-.6543129546637614E+00,
     5  -.6299720177071481E+00,-.5926164691150697E+00,
     6  -.5141087139732419E+00,-.3950695488428987E+00,
     7  -.2416649757722363E+00,-.6351743815157737E-01,
     8  0.1278539717868201E+00,0.3202187720663667E+00,
     9  0.5017323932204512E+00,0.6619992009285086E+00,
     *  0.7930785575406252E+00,0.8903953097017818E+00,
     *  0.9535074152065579E+00,0.9866151763023578E+00,
     *  0.9983991227155065E+00/
c 
        data wts7/
     1  0.1594086818366365E-01,0.3826366034727283E-02,
     2  0.7630268404446956E-01,0.3474219073343622E-01,
     3  0.5513028244722615E-01,0.7096552659584025E-01,
     4  0.6630537428373043E-01,0.3896225940169057E-01,
     5  0.1991315091678022E-01,0.5752347057878665E-01,
     6  0.9940927244426295E-01,0.1376008089991970E+00,
     7  0.1675802729819557E+00,0.1867757780002242E+00,
     8  0.1939104065062025E+00,0.1888273088077420E+00,
     9  0.1724504645454698E+00,0.1467444401628223E+00,
     *  0.1146387555951936E+00,0.7991181740075362E-01,
     *  0.4701253478751221E-01,0.2067777678957278E-01,
     *  0.4848189699670056E-02/
c 
        data pts8/
     1  -.9989114181274066E+00,-.9908263966048775E+00,
     2  -.9248645385534834E+00,-.9680116531931305E+00,
     3  -.6840167682475757E+00,-.8600325377489907E+00,
     4  -.7770195928476985E+00,-.5932448040420801E+00,
     5  -.2534205810243058E+00,-.5195991172012450E+00,
     6  -.4778422442420079E+00,-.4466567023701029E+00,
     7  -.3711085922410834E+00,-.1038085618427177E+00,
     8  0.6624108926931165E-01,0.2447884282929578E+00,
     9  0.4202089384557705E+00,0.5819920401729431E+00,
     *  0.7215408668599309E+00,0.9134409255064474E+00,
     *  0.9641674707412207E+00,0.8329127612727895E+00,
     *  0.9899437135491419E+00,0.9988246725687248E+00/
c 
        data wts8/
     1  0.3305353360059100E-02,0.1423309004293877E-01,
     2  0.5416515984696713E-01,0.3236281304748274E-01,
     3  0.9420972164907697E-01,0.7491136151447582E-01,
     4  0.8973875593860681E-01,0.8480629314794851E-01,
     5  0.1354880049097897E+00,0.5998816889482263E-01,
     6  0.2418860668583837E-01,0.5142746899913080E-01,
     7  0.9816931779242091E-01,0.1618084085520072E+00,
     8  0.1762882726742461E+00,0.1788555688602092E+00,
     9  0.1702202045019091E+00,0.1519007081057399E+00,
     *  0.1261959967987976E+00,0.6513172519174726E-01,
     *  0.3713814336179654E-01,0.9609582661491900E-01,
     *  0.1578879020301863E-01,0.3582239310318977E-02/
c 
        data pts9/
     1  -.9989201068540439E+00,-.9908306628467857E+00,
     2  -.9676790427257087E+00,-.8540712929424172E+00,
     3  -.9230170931964715E+00,-.7624421485560813E+00,
     4  -.6541330979461253E+00,-.5392191605625710E+00,
     5  -.4311985591686158E+00,-.3457311526203316E+00,
     6  -.2979586148451097E+00,-.2644151396426910E+00,
     7  -.1851579612693086E+00,-.6422453342164046E-01,
     8  0.8546829348720772E-01,0.2504468077601715E+00,
     9  0.4177249220644361E+00,0.5755736448803112E+00,
     *  0.7142513251753704E+00,0.9092557278046985E+00,
     *  0.8267136498818403E+00,0.9620196005415709E+00,
     *  0.9892255165227969E+00,0.9987290604067465E+00/
c 
        data wts9/
     1  0.3285081819550222E-02,0.1431124160077157E-01,
     2  0.3309868880137767E-01,0.8092727162673100E-01,
     3  0.5669741557269108E-01,0.1013247961248490E+00,
     4  0.1136038602523699E+00,0.1139555100218008E+00,
     5  0.9943290254244755E-01,0.6893774422582047E-01,
     6  0.2731452931929789E-01,0.5437885841616585E-01,
     7  0.1021300809023550E+00,0.1375557491680790E+00,
     8  0.1595680918345117E+00,0.1682122827151778E+00,
     9  0.1643769138749959E+00,0.1496849626552296E+00,
     *  0.1264855049557638E+00,0.6729022734433134E-01,
     *  0.9781293505839811E-01,0.3894612634427702E-01,
     *  0.1680519835865629E-01,0.3864026465840515E-02/
c 
        data pts10/
     1  -.9990364720577143E+00,-.9916405489286161E+00,
     2  -.9697948353744005E+00,-.9262428588598789E+00,
     3  -.8567443830340244E+00,-.7610609532652768E+00,
     4  -.6431652799451363E+00,-.5110268108986621E+00,
     5  -.3761419820725773E+00,-.2528430178006156E+00,
     6  -.1571302562947860E+00,-.1042033041332202E+00,
     7  -.6930207366398028E-01,0.1107303780459912E-01,
     8  0.1313013038346218E+00,0.2761183097514245E+00,
     9  0.4304171787531730E+00,0.5805677267932449E+00,
     *  0.7152362633010458E+00,0.9083511369673156E+00,
     *  0.8260768296084563E+00,0.9614247495044578E+00,
     *  0.9890017128941841E+00,0.9986973827198142E+00/
c 
        data wts10/
     1  0.2948228576574216E-02,0.1324733851987531E-01,
     2  0.3169029455011405E-01,0.5613542493016808E-01,
     3  0.8290393843846391E-01,0.1077968674069402E+00,
     4  0.1266582767694810E+00,0.1356914803876272E+00,
     5  0.1316721412114242E+00,0.1122205261981140E+00,
     6  0.7662686069116662E-01,0.2990699465684278E-01,
     7  0.5553329933196593E-01,0.1027599340275091E+00,
     8  0.1351120192546869E+00,0.1519840800494351E+00,
     9  0.1543291071354313E+00,0.1440715549910145E+00,
     *  0.1238640618914496E+00,0.6741933372242020E-01,
     *  0.9701882705457603E-01,0.3935012621019759E-01,
     *  0.1710320840362771E-01,0.3956075590872511E-02/
c 
        do 1050 i=1,20
        do 1040 j=1,24
c 
        points(j,i)=0
        weights(j,i)=0
 1040 continue
 1050 continue
c 
c       build the matrices points, weights
c 
        call o1rarrcp(pts1,points(1,1),20)
        call o1rarrcp(wts1,weights(1,1),20)
        call o1rarinv(pts1,wts1,20,points(1,20),weights(1,20) )
c 
c 
        call o1rarrcp(pts2,points(1,2),20)
        call o1rarrcp(wts2,weights(1,2),20)
        call o1rarinv(pts2,wts2,20,points(1,19),weights(1,19) )
c 
c 
        call o1rarrcp(pts3,points(1,3),21)
        call o1rarrcp(wts3,weights(1,3),21)
        call o1rarinv(pts3,wts3,21,points(1,18),weights(1,18) )
c 
        call o1rarrcp(pts4,points(1,4),21)
        call o1rarrcp(wts4,weights(1,4),21)
        call o1rarinv(pts4,wts4,21,points(1,17),weights(1,17) )
c 
c 
        call o1rarrcp(pts5,points(1,5),22)
        call o1rarrcp(wts5,weights(1,5),22)
        call o1rarinv(pts5,wts5,22,points(1,16),weights(1,16) )
c 
c 
        call o1rarrcp(pts6,points(1,6),23)
        call o1rarrcp(wts6,weights(1,6),23)
        call o1rarinv(pts6,wts6,23,points(1,15),weights(1,15) )
c 
        call o1rarrcp(pts7,points(1,7),23)
        call o1rarrcp(wts7,weights(1,7),23)
        call o1rarinv(pts7,wts7,23,points(1,14),weights(1,14) )
c 
c 
        call o1rarrcp(pts8,points(1,8),24)
        call o1rarrcp(wts8,weights(1,8),24)
        call o1rarinv(pts8,wts8,24,points(1,13),weights(1,13) )
c 
        call o1rarrcp(pts9,points(1,9),24)
        call o1rarrcp(wts9,weights(1,9),24)
        call o1rarinv(pts9,wts9,24,points(1,12),weights(1,12) )
c 
        call o1rarrcp(pts10,points(1,10),24)
        call o1rarrcp(wts10,weights(1,10),24)
        call o1rarinv(pts10,wts10,24,points(1,11),weights(1,11) )
c 
c       return to the user the numbers of nodes for each column
c 
        numspts(1)=20
        numspts(2)=20
        numspts(20)=20
        numspts(19)=20
c 
        numspts(3)=21
        numspts(4)=21
        numspts(18)=21
        numspts(17)=21
c 
        numspts(5)=22
        numspts(16)=22
c 
        numspts(6)=23
        numspts(7)=23
        numspts(15)=23
        numspts(14)=23
c 
        numspts(8)=24
        numspts(9)=24
        numspts(10)=24
        numspts(13)=24
        numspts(12)=24
        numspts(11)=24
c 
        return
        end
  
c 
c 
c 
c 
        subroutine o1rarinv(ptsin,wtsin,n,ptsout,wtsout)
        implicit real *8 (a-h,o-z)
        save
        dimension ptsin(1),wtsin(1),ptsout(1),wtsout(1)
c 
        do 1200 i=1,n
        ptsout(i)=-ptsin(n+1-i)
        wtsout(i)=wtsin(n+1-i)
 1200 continue
        return
        end
c 
c 
c 
c 
c 
c 
        subroutine o1rarrcp(a,b,n)
        implicit real *8 (a-h,o-z)
        save
        dimension a(1),b(1)
c 
        do 1200 i=1,n
        b(i)=a(i)
 1200 continue
        return
        end
c 
c 
c 
c 
c 
c 
        subroutine o1rnodes(points,weights,n)
        implicit real *8 (a-h,o-z)
        save
        dimension points(20),weights(20),pts(20),whts(20)
c 
c       This subroutine returns to the user the nodes and weights of
c        the 20-point Generalized Gaussian quadrature on the interval
c        [-1,1] for functions of the form
c 
c        P_n(x), P_n (x) \cdot log(1-x)\cdot (1-x),
c        P_n (x) \cdot log(1+x) \cdot(1+x), P_n(x),                 (2)
c        P_n (x) \cdot log(1-x)^2, P_n (x) \cdot log(1+x)^2,
c 
c        with n running from 0 to 19. The accuracy corresponding
c        this quadrature is roughly 0.2E-9.
C 
C           PLEASE NOTE THAT THIS SUBROUTIINE HAS NO INPUT PARAMETERS
c 
        data pts/
     1  -.9993068470120125E+00,-.9933991137301627E+00,
     2  -.9333462961063758E+00,-.9742484922962135E+00,
     3  -.2937620961429727E+00,-.8643432243313858E+00,
     4  -.7641305230113280E+00,-.6329709181566937E+00,
     5  -.4742072478627974E+00,-.9950131107975985E-01,
     6  0.9950131107975985E-01,0.2937620961429727E+00,
     7  0.4742072478627974E+00,0.6329709181566937E+00,
     8  0.8643432243313858E+00,0.9333462961063758E+00,
     9  0.7641305230113280E+00,0.9742484922962135E+00,
     *  0.9933991137301627E+00,0.9993068470120125E+00/
c 
        data whts/
     1  0.2184219177273501E-02,0.1103553719020169E-01,
     2  0.5414947353139745E-01,0.2872900411930890E-01,
     3  0.1887904383174358E+00,0.8437575093430985E-01,
     4  0.1160005877353149E+00,0.1457530064049815E+00,
     5  0.1707765505653168E+00,0.1982054320258153E+00,
     6  0.1982054320258153E+00,0.1887904383174358E+00,
     7  0.1707765505653168E+00,0.1457530064049815E+00,
     8  0.8437575093430985E-01,0.5414947353139745E-01,
     9  0.1160005877353149E+00,0.2872900411930890E-01,
     *  0.1103553719020169E-01,0.2184219177273501E-02/
c 
        n=20
c 
        do 1200 i=1,20
        points(i)=pts(i)
        weights(i)=whts(i)
 1200 continue
  
        call o1rbubbl(n,points,weights)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o1rbubbl(n,xs,ws)
        implicit real *8 (a-h,o-z)
        save
        dimension ws(1),xs(1)
c 
        do 1600 i=1,n
c 
        do 1400 j=1,n-1
c 
        if(xs(j+1) .gt. xs(j)) goto 1400
c 
        d=xs(j+1)
        xs(j+1)=xs(j)
        xs(j)=d
c 
        d=ws(j+1)
        ws(j+1)=ws(j)
        ws(j)=d
c 
 1400 continue
 1600 continue
  
        return
        end
c 
c 
c 
c 
c 
        subroutine o1rinter(x,forminte,points,weights,nn)
        implicit real *8 (a-h,o-z)
        save
        dimension singcoef(120,20),bmatr(20,20),amatr(20,20),
     1      forminte(1),vals(22),ders(22)
c 
        data ifcalled/0/
c 
c       This subroutine returns to the user the coefficients
c       forminte of the interpolation formula expressing the
c       values of a function at the user-specified point x
c       on the interval [-1,1] via that function's values at
c       the nn nodes points.
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
c       The accuracy of the interpolation is about 7 digits.
c       It also returns the nodes points (20 of them things),
c       and the corresponding quadrature weights weights.
c       Finally, just in case the user has forgotten, it returns
c       the number nn of interpolation nodes (which is equal to 20).
c 
c 
c       . . . if needed, retrieve various matrices from the
c             subroutines
c 
        if(ifcalled .eq. 1) goto 1200
  
        call o1runiv(singcoef,nuniv,n7)
        call o1rinmat(bmatr,amatr,n)
        ifcalled=1
c 
 1200 continue
c 
c        evaluate the singular functions at the
c        user-supplied point x
c 
  
  
cccc        do 1600 i=1,n
c 
cccc        ndigits=16
cccc        call protFDER(ndigits,x,VAL,der,singcoef(1,i),Nuniv-1)
cccc        vals(i)=val
 1600 continue
  
        ndigits=16
  
        call protFDEm(ndigits,X,VALs,ders,singcoef,nuniv,n)
  
  
c 
c        multiply the obtained values by the matrix amatr^*
c 
        call o1rmatav(amatr,vals,forminte,n)
c 
c       return to the user the nodes on which the interpolation is
c       based, and the corresponding weights
c 
        call o1rnodes(points,weights,nn)
  
  
        return
        end
  
c 
c 
c 
c 
c 
        subroutine o1runiv(coefs,nuniv,n)
        implicit real *8 (a-h,o-z)
        save
        dimension c1(44),c2(43), c3(46),c4(46),c5(47),c6(47),
     1      c7(49),c8(48),c9(50),c10(48),c11(51),c12(50),
     2      c13(52),c14(58),c15(60),c16(60),c17(60),c18(60),
     3      c19(60),c20(60),coefs(120,20)
c 
c       This subroutine returns to the user the coefficients of the
c       universal expansions of the first 20 singular functions on
c       the interval representing the functions of the form
c 
c          f(x) = phi(x) + psi(x) * (1-x) * log(1-x) +
c                 eta(x) * (1-x) * log(1+x),                          (3)
c 
c       with phi, psi, eta polynomials of order 9. While each of the
c       functions is approximated to 14 digits or so, the resulting
c       20-function expansion has accuracy of about 7 digits on
c       functions of the form (3). Also, see subroutine o1rinmat.
C       IMPORTANT: THE UNIVERSAL EXPANSION USED IS THE ONE WITH
C       THE PARAMETER NDIGITS SET TO 16!!!!
c 
        data c1/
     1   0.1014624426031295E+01,0.2203482759988576E+00,
     2   -.1591009198941792E+00,0.9151667975247263E-01,
     3   -.4097098545323162E-01,0.1337065935992574E-01,
     4   -.2382057928914528E-02,-.5015472825092817E-03,
     5   0.6872456043219344E-03,-.3851047809468416E-03,
     6   0.1631650416487348E-03,-.5870410724312053E-04,
     7   0.1878118151944495E-04,-.5758633362254657E-05,
     8   0.2399449364987992E-05,-.2017236404055450E-05,
     9   0.2197792017795687E-05,-.2187577520361411E-05,
     *   0.1894896544505113E-05,-.1450602755181647E-05,
     1   0.1000737539338477E-05,-.6319979247152753E-06,
     2   0.3697750299891326E-06,-.2022922501335345E-06,
     3   0.1042213889309764E-06,-.5085630280555449E-07,
     4   0.2361233111474219E-07,-.1047058793691034E-07,
     5   0.4448273760412075E-08,-.1815215579825346E-08,
     6   0.7130688833149083E-09,-.2701514948144871E-09,
     7   0.9886515610793474E-10,-.3499674860103651E-10,
     8   0.1199681598438029E-10,-.3986507348173612E-11,
     9   0.1285229454269663E-11,-.4023005737621433E-12,
     *   0.1223420207046610E-12,-.3616502508509836E-13,
     1   0.1039638230700713E-13,-.2907539591025007E-14,
     2   0.7912688227851526E-15,-.2096371852762742E-15/
c 
         data c2/
     1   -.1730851449758185E+01,0.4114890443981512E+00,
     2   -.1105895184200979E+00,0.4704589537206696E-02,
     3   0.2111995339677963E-01,-.1768467517667596E-01,
     4   0.9246010225440549E-02,-.3587091343041205E-02,
     5   0.1011588203328084E-02,-.1257226916357581E-03,
     6   -.8796043399772068E-04,0.9345183649356086E-04,
     7   -.5508197896905778E-04,0.2055804468434017E-04,
     8   0.5441682455882951E-06,-.9900342496229651E-05,
     9   0.1176607333664545E-04,-.1000278262094663E-04,
     *   0.7175477685863621E-05,-.4589612469551027E-05,
     1   0.2687824665861227E-05,-.1463617551319985E-05,
     2   0.7486068767035575E-06,-.3622425027429658E-06,
     3   0.1667241648570250E-06,-.7329330272958642E-07,
     4   0.3087761874912447E-07,-.1250016974006007E-07,
     5   0.4873666493042448E-08,-.1833495347110990E-08,
     6   0.6666186688013322E-09,-.2345498756866206E-09,
     7   0.7995656213161228E-10,-.2643404829503242E-10,
     8   0.8482665835636220E-11,-.2644090010328214E-11,
     9   0.8010601334827089E-12,-.2360096010026406E-12,
     *   0.6764891927916711E-13,-.1887247597774518E-13,
     1   0.5125485691604433E-14,-.1355800350043420E-14,
     2   0.3490282279767343E-15/
c 
         data c3/
     1   0.1432772409833407E+01,0.1159645890375588E+01,
     2   -.7511814924612322E+00,0.3269307686329802E+00,
     3   -.4890174480340014E-01,-.5999892897323880E-01,
     4   0.6548190523401945E-01,-.3808862444599538E-01,
     5   0.1440918929396985E-01,-.2453102275651639E-02,
     6   -.1325664213298014E-02,0.1556497836094835E-02,
     7   -.9086309086301312E-03,0.3634467613400352E-03,
     8   -.8254336817032253E-04,-.1931098835948696E-04,
     9   0.3905547716978073E-04,-.3217467235356567E-04,
     *   0.2130893534062654E-04,-.1286856292583482E-04,
     1   0.7421098130423124E-05,-.4153507790324235E-05,
     2   0.2260572415586706E-05,-.1192045733460103E-05,
     3   0.6063466870596108E-06,-.2965651963261276E-06,
     4   0.1392502150787781E-06,-.6274924907675223E-07,
     5   0.2715098901579164E-07,-.1129184620888009E-07,
     6   0.4519446194682652E-08,-.1743080191363788E-08,
     7   0.6486700085661665E-09,-.2332041107013762E-09,
     8   0.8108569629151228E-10,-.2729552397369431E-10,
     9   0.8903774848129543E-11,-.2816718879498224E-11,
     *   0.8647821873636294E-12,-.2578291966605422E-12,
     1   0.7468727174897921E-13,-.2103063337496637E-13,
     2   0.5758232900581225E-14,-.1533781558921379E-14,
     3   0.3972926765247625E-15,-.1003371440371556E-15/
c 
         data c4/
     1   -.3059821722187122E+01,0.4755481405391927E-02,
     2   0.6346523607997984E+00,-.5664169797336385E+00,
     3   0.2729820240793735E+00,-.3850438243974925E-01,
     4   -.5948279748213747E-01,0.6437736194073687E-01,
     5   -.3779448237019504E-01,0.1410966183686986E-01,
     6   -.1852388777696253E-02,-.1954758973445417E-02,
     7   0.1961600856058159E-02,-.1065585795197142E-02,
     8   0.3763245448353597E-03,-.5548306567792227E-04,
     9   -.3276124796967241E-04,0.2859258004808152E-04,
     *   -.7125472567506281E-05,-.6617171399637703E-05,
     1   0.1072836287045918E-04,-.9486648614032228E-05,
     2   0.6621524833619829E-05,-.4006532507804672E-05,
     3   0.2188409015927523E-05,-.1103001678123626E-05,
     4   0.5200014659397347E-06,-.2314056526304908E-06,
     5   0.9783776474020645E-07,-.3949127546500921E-07,
     6   0.1527459855292282E-07,-.5677813463666238E-08,
     7   0.2033071442190653E-08,-.7026068030900157E-09,
     8   0.2347152704360248E-09,-.7589342758720719E-10,
     9   0.2377776878303350E-10,-.7224943798924909E-11,
     *   0.2130697754619147E-11,-.6102448088868158E-12,
     1   0.1698272599957474E-12,-.4594161155928879E-13,
     2   0.1208563467367206E-13,-.3091869336891848E-14,
     3   0.7699861241969947E-15,-.1860720284458363E-15/
c 
         data c5/
     1   0.1914652304944146E+01,0.1780905128147911E+01,
     2   -.7377758207547431E+00,-.1905366255063196E+00,
     3   0.5726301383223017E+00,-.4704328507539701E+00,
     4   0.1945805562112206E+00,0.1478918602919689E-01,
     5   -.9334766064411922E-01,0.8403506993827883E-01,
     6   -.4617375154834113E-01,0.1427146656763888E-01,
     7   0.2429028053978130E-02,-.7336738389070120E-02,
     8   0.6439699482023649E-02,-.3977392255474604E-02,
     9   0.1889600710821739E-02,-.6367492659773374E-03,
     *   0.5717193038793133E-04,0.1330059389515068E-03,
     1   -.1491544470326029E-03,0.1105609943112641E-03,
     2   -.6836438611726545E-04,0.3773951723049919E-04,
     3   -.1917255697963136E-04,0.9114481594264448E-05,
     4   -.4097069444678026E-05,0.1753748852031948E-05,
     5   -.7184791249109554E-06,0.2827883607411610E-06,
     6   -.1072464483441715E-06,0.3928124393745250E-07,
     7   -.1392125230740185E-07,0.4781029590970716E-08,
     8   -.1593146636807636E-08,0.5156172686173614E-09,
     9   -.1622201927185400E-09,0.4964682746365940E-10,
     *   -.1478889219875492E-10,0.4289810792787025E-11,
     1   -.1212158235096659E-11,0.3337498538831280E-12,
     2   -.8956104601568461E-13,0.2342611360429983E-13,
     3   -.5973853702324757E-14,0.1484485795502937E-14,
     4   -.3602309534918967E-15/
c 
         data c6/
     1   -.3338032858678335E+01,-.3639487274698308E+00,
     2   0.8570914452058929E+00,-.3348876114986857E+00,
     3   -.2836914616119561E+00,0.5185468700493169E+00,
     4   -.3914251562571649E+00,0.1453386877297107E+00,
     5   0.3471385728389443E-01,-.1013149769697706E+00,
     6   0.9041710034501964E-01,-.5228595016701571E-01,
     7   0.1832640181646768E-01,0.1348218603645518E-02,
     8   -.8662286653108970E-02,0.8961270915930708E-02,
     9   -.6639858744289267E-02,0.4101872422019268E-02,
     *   -.2221442799206465E-02,0.1076364041389852E-02,
     1   -.4691263890510002E-03,0.1824439911437178E-03,
     2   -.6140735838243287E-04,0.1625421821382279E-04,
     3   -.1994929365204237E-05,-.1286161978316937E-05,
     4   0.1367864284498641E-05,-.8595332647898792E-06,
     5   0.4449516518470250E-06,-.2056959428451130E-06,
     6   0.8776485677679019E-07,-.3515470959393022E-07,
     7   0.1335320003507464E-07,-.4841274975204363E-08,
     8   0.1682927184510518E-08,-.5627565660943678E-09,
     9   0.1814632876315514E-09,-.5653108109675400E-10,
     *   0.1703946718961004E-10,-.4975096578507079E-11,
     1   0.1408407972674138E-11,-.3868633541224491E-12,
     2   0.1031692302506335E-12,-.2672323033611459E-13,
     3   0.6726597618150572E-14,-.1644762337793412E-14,
     4   0.3918603090905397E-15/
c 
        data c7/
     1   0.1660431894795212E+01,0.1720218854456644E+01,
     2   -.5425827783359173E+00,-.3467187289265464E+00,
     3   0.4785972807421009E+00,-.9199977370891760E-01,
     4   -.3102389309940158E+00,0.4436720444184636E+00,
     5   -.3336698589992716E+00,0.1392833202088531E+00,
     6   0.1437351482067516E-01,-.8661359267774497E-01,
     7   0.9414597621935779E-01,-.7026650614839306E-01,
     8   0.4077552798764419E-01,-.1809772533632844E-01,
     9   0.4775903624365938E-02,0.1252872965760422E-02,
     *   -.2999850760957242E-02,0.2804207339645191E-02,
     1   -.2023292158690905E-02,0.1266647661951969E-02,
     2   -.7190632098052250E-03,0.3784701561127049E-03,
     3   -.1871274366522808E-03,0.8765987496522363E-04,
     4   -.3914042149442893E-04,0.1673128360382293E-04,
     5   -.6870334701250473E-05,0.2717209909272718E-05,
     6   -.1037264510226121E-05,0.3828518486608489E-06,
     7   -.1368261630433233E-06,0.4740486126531153E-07,
     8   -.1593771116287071E-07,0.5204068944156850E-08,
     9   -.1651508155184144E-08,0.5096782416090912E-09,
     *   -.1530395007822310E-09,0.4472805244744742E-10,
     1   -.1272826730275652E-10,0.3527647561856794E-11,
     2   -.9523908746139525E-12,0.2505067481210349E-12,
     3   -.6420129633815424E-13,0.1603164612305047E-13,
     4   -.3901545897249506E-14,0.9241915733312452E-15,
     5   -.2143157619355872E-15/
c 
        data c8/
     1   -.2405037914088231E+01,-.3739367228629750E+00,
     2   0.6296768882409492E+00,-.1287690097402147E+00,
     3   -.2883972329540123E+00,0.2704689723357937E+00,
     4   0.2647061972634297E-01,-.2965770770224379E+00,
     5   0.3801078147905961E+00,-.3003254726765234E+00,
     6   0.1573393583338586E+00,-.3320479313346921E-01,
     7   -.3922946632694415E-01,0.6418673187056088E-01,
     8   -.6028441288459107E-01,0.4499419809748252E-01,
     9   -.2900850530137734E-01,0.1673468272100661E-01,
     *   -.8795677593000201E-02,0.4252682556559478E-02,
     1   -.1899435885173910E-02,0.7832363445577256E-03,
     2   -.2961903643425782E-03,0.1009581107748268E-03,
     3   -.2972393784791080E-04,0.6626442939103431E-05,
     4   -.3807795054105594E-06,-.7240072561261046E-06,
     5   0.5935538277137753E-06,-.3306773233565318E-06,
     6   0.1562828192886403E-06,-.6677388942990924E-07,
     7   0.2651359937970684E-07,-.9927849109156045E-08,
     8   0.3536590205260014E-08,-.1205453072349940E-08,
     9   0.3947104841642986E-09,-.1245141449290538E-09,
     *   0.3792300805767180E-10,-.1116973339474986E-10,
     1   0.3185600670466303E-11,-.8806014251444525E-12,
     2   0.2361288583108840E-12,-.6145451775178403E-13,
     3   0.1553243928000830E-13,-.3812289344177878E-14,
     4   0.9104712648057645E-15,-.2099034144643296E-15/
c 
        data c9/
     1   0.1021870433828945E+01,0.1088256802426908E+01,
     2   -.2903651113598750E+00,-.2677472615743635E+00,
     3   0.2937988512620586E+00,-.1962947320534052E-01,
     4   -.1771380909036348E+00,0.1379879600964094E+00,
     5   0.4943918123838295E-01,-.2270018599215209E+00,
     6   0.3043701513549334E+00,-.2810718668189551E+00,
     7   0.2030029046824277E+00,-.1171314451266044E+00,
     8   0.5051341047322025E-01,-.1014056108138856E-01,
     9   -.8614051625401072E-02,0.1386001890175277E-01,
     *   -.1258324431876837E-01,0.9216568263126583E-02,
     1   -.5936680619428451E-02,0.3486102857107024E-02,
     2   -.1902208288151650E-02,0.9758282644496111E-03,
     3   -.4743408144275735E-03,0.2197068267339273E-03,
     4   -.9737683029980091E-04,0.4143239484378694E-04,
     5   -.1696763929878967E-04,0.6702131065928434E-05,
     6   -.2557805038651916E-05,0.9445244322419493E-06,
     7   -.3378928360639417E-06,0.1172229631209450E-06,
     8   -.3947254857634314E-07,0.1291072794843747E-07,
     9   -.4104453660364809E-08,0.1268945767431188E-08,
     *   -.3816901601117512E-09,0.1117437698103526E-09,
     1   -.3185039421187923E-10,0.8840830758251727E-11,
     2   -.2390237948961727E-11,0.6295326233739610E-12,
     3   -.1615340747411999E-12,0.4038218397700164E-13,
     4   -.9836060258423357E-14,0.2333494815836368E-14,
     5   -.5399362677271628E-15,0.1209431837595344E-15/
c 
        data c10/
     1   -.1754907593448233E+01,-.3080353818572991E+00,
     2   0.4592936757104411E+00,-.5451187419599070E-01,
     3   -.2426148950680476E+00,0.1963032020323486E+00,
     4   0.2008842465174696E-01,-.1547004259533388E+00,
     5   0.1162819133501542E+00,0.2908643123590042E-01,
     6   -.1738949790546872E+00,0.2532137705934164E+00,
     7   -.2593256317858287E+00,0.2174483521875697E+00,
     8   -.1582375846697281E+00,0.1029695968314041E+00,
     9   -.6098997108855168E-01,0.3325387952207986E-01,
     *   -.1681046782905654E-01,0.7912084049171186E-02,
     1   -.3472475174031034E-02,0.1419058043966166E-02,
     2   -.5369545289956644E-03,0.1857097599683943E-03,
     3   -.5702953018216772E-04,0.1441079568149040E-04,
     4   -.2173545906061332E-05,-.5096667538571258E-06,
     5   0.6816393623804375E-06,-.4200023009580141E-06,
     6   0.2066887591110569E-06,-.9008527547650890E-07,
     7   0.3615008338457468E-07,-.1361260827365299E-07,
     8   0.4862504769529535E-08,-.1658937807970825E-08,
     9   0.5430632954483394E-09,-.1711344556148768E-09,
     *   0.5203936922510550E-10,-.1529746311688404E-10,
     1   0.4353148721389528E-11,-.1200473648353511E-11,
     2   0.3210968392529114E-12,-.8335457997058557E-13,
     3   0.2101313685435413E-13,-.5144345234219298E-14,
     4   0.1225325771356653E-14,-.2819677961412056E-15/
c 
        data c11/
     1   0.1318483448292746E+01,0.1424733735077770E+01,
     2   -.3598559754104600E+00,-.3607700324799663E+00,
     3   0.3603772852296439E+00,0.1924050802551145E-01,
     4   -.2666957922395189E+00,0.2013576826050846E+00,
     5   0.2175721397632648E-01,-.1821569237326660E+00,
     6   0.1843136843734384E+00,-.6628241009880048E-01,
     7   -.8243461414457010E-01,0.1922127903330012E+00,
     8   -.2378469044453372E+00,0.2286328338518946E+00,
     9   -.1879737244264955E+00,0.1378793244187425E+00,
     *   -.9236321266870782E-01,0.5734782979434685E-01,
     1   -.3334024215570750E-01,0.1828404996111100E-01,
     2   -.9512037126114067E-02,0.4715091412102689E-02,
     3   -.2234926603969698E-02,0.1015909952531135E-02,
     4   -.4439311848255191E-03,0.1868651229381122E-03,
     5   -.7590042891216936E-04,0.2979251667878218E-04,
     6   -.1131545017095756E-04,0.4163114452849611E-05,
     7   -.1485127050930001E-05,0.5141265415083487E-06,
     8   -.1728433273345039E-06,0.5646569884133910E-07,
     9   -.1793504773398821E-07,0.5541266765644197E-08,
     *   -.1666007850519906E-08,0.4875864691640987E-09,
     1   -.1389480739353443E-09,0.3856342634983132E-10,
     2   -.1042544434408938E-10,0.2745755697773403E-11,
     3   -.7045475346226300E-12,0.1761376157768790E-12,
     4   -.4290281723491354E-13,0.1018047487209511E-13,
     5   -.2353775153109711E-14,0.5294539563552410E-15,
     6   -.1167136592804185E-15/
c 
        data c12/
     1   -.3910265178365073E+01,-.7883317828962073E+00,
     2   0.1004622397314356E+01,-.1036455259832183E-01,
     3   -.5675860104246986E+00,0.3197203164113688E+00,
     4   0.1852664867234964E+00,-.3747498332063957E+00,
     5   0.1873979733738347E+00,0.1055114129667260E+00,
     6   -.2609417196470416E+00,0.2228514514412575E+00,
     7   -.7181731088185798E-01,-.8562090523359557E-01,
     8   0.1853009524293142E+00,-.2153398048587505E+00,
     9   0.1951058311977983E+00,-.1514533971642995E+00,
     *   0.1049778547963525E+00,-.6647749029274782E-01,
     1   0.3902208846509622E-01,-.2144619751256465E-01,
     2   0.1111655336970969E-01,-.5465015093889503E-02,
     3   0.2559283285508628E-02,-.1145737843285124E-02,
     4   0.4917575480702386E-03,-.2028448852026233E-03,
     5   0.8057632063602150E-04,-.3087676223765602E-04,
     6   0.1143087433121312E-04,-.4093578920400741E-05,
     7   0.1419645637697309E-05,-.4772240786979703E-06,
     8   0.1556278059865310E-06,-.4927001332667339E-07,
     9   0.1515214268738738E-07,-.4528853484307677E-08,
     *   0.1316187804827933E-08,-.3720694489186498E-09,
     1   0.1023385805590758E-09,-.2739489698932678E-10,
     2   0.7138348124476685E-11,-.1810824004753341E-11,
     3   0.4472644092542819E-12,-.1075428443905259E-12,
     4   0.2519682065715718E-13,-.5726909664818775E-14,
     5   0.1288648717237312E-14,-.2599541324780394E-15/
c 
        data c13/
     1   0.3252599358928982E+01,0.3739692106038892E+01,
     2   -.6075005194550365E+00,-.1103055204816420E+01,
     3   0.6669719697386277E+00,0.3717735010497032E+00,
     4   -.6232755369364642E+00,0.1344322868183556E+00,
     5   0.3396412171143371E+00,-.3776364638121038E+00,
     6   0.1023402567082807E+00,0.1752537072379733E+00,
     7   -.2724052430647344E+00,0.1947767200179836E+00,
     8   -.4228181230266224E-01,-.9192821238835947E-01,
     9   0.1648645372779238E+00,-.1779469668732115E+00,
     *   0.1533439923644182E+00,-.1143906524589044E+00,
     1   0.7665126105360074E-01,-.4711383724402943E-01,
     2   0.2692241776046681E-01,-.1443694670364089E-01,
     3   0.7315058971319998E-02,-.3520688613730451E-02,
     4   0.1616251392802259E-02,-.7100989258284250E-03,
     5   0.2994026166226441E-03,-.1214278415944417E-03,
     6   0.4746236525378482E-04,-.1790875416002571E-04,
     7   0.6532502242310017E-05,-.2306319026163371E-05,
     8   0.7889322840692030E-06,-.2617179244807965E-06,
     9   0.8426354151030161E-07,-.2634817096135959E-07,
     *   0.8006008771299508E-08,-.2365099983867357E-08,
     1   0.6795633993182915E-09,-.1899793311561509E-09,
     2   0.5168894523788581E-10,-.1368996252631663E-10,
     3   0.3529981816931291E-11,-.8863866572758307E-12,
     4   0.2166161757760794E-12,-.5166750464513320E-13,
     5   0.1187663257163247E-13,-.2784578359024831E-14,
     6   0.5049480078595234E-15,-.2174293782210839E-15/
c 
        data c14/
     1   -.8098576385295602E+01,-.2239439214567230E+01,
     2   0.1865835372009047E+01,0.6043292821872790E+00,
     3   -.1144425524728754E+01,0.2637174844063756E-01,
     4   0.7090845734993222E+00,-.3936749018459129E+00,
     5   -.2218755041902526E+00,0.4305127614909125E+00,
     6   -.2088953409149929E+00,-.1062321927417374E+00,
     7   0.2599961812592745E+00,-.2168263807614357E+00,
     8   0.7247666729234973E-01,0.6796032844600701E-01,
     9   -.1501587089576065E+00,0.1701598250251046E+00,
     *   -.1493410070866749E+00,0.1120741295710143E+00,
     1   -.7504607525251058E-01,0.4589770682343086E-01,
     2   -.2601799606831538E-01,0.1380869005685841E-01,
     3   -.6912131779814782E-02,0.3281493664554221E-02,
     4   -.1483978287086505E-02,0.6415089610007565E-03,
     5   -.2658538081008421E-03,0.1058722367774040E-03,
     6   -.4059648752971395E-04,0.1501402743936643E-04,
     7   -.5363336597315575E-05,0.1852843414596999E-05,
     8   -.6196823113043202E-06,0.2008271330177049E-06,
     9   -.6311590744931833E-07,0.1924907241593046E-07,
     *   -.5700111917717063E-08,0.1639718622361682E-08,
     1   -.4583969598458120E-09,0.1245793665302479E-09,
     2   -.3292211904142351E-10,0.8462294996895608E-11,
     3   -.2115296469508169E-11,0.5149382779786823E-12,
     4   -.1213802996414883E-12,0.2842216279727254E-13,
     5   -.5868727791937126E-14,0.1787903741759284E-14,
     6   0.1426881443909788E-15,0.5018121711718336E-15,
     7   0.4382613591183610E-15,0.4587283061729779E-15,
     8   0.4603065500860079E-15,0.4633628990832915E-15,
     9   0.4640741899622938E-15,0.4632676389688305E-15/
  
        data c15/
     1   0.5726673046411160E+01,0.7169687177370464E+01,
     2   -.2254251183764721E+00,-.2317259890458057E+01,
     3   0.3889868057394278E+00,0.1222795207548875E+01,
     4   -.5697150435381907E+00,-.5600500045086752E+00,
     5   0.6261566900019701E+00,0.1013898735940641E-01,
     6   -.4214679362775739E+00,0.3208985186829021E+00,
     7   0.6170810444057589E-02,-.2301906541814077E+00,
     8   0.2441011274954404E+00,-.1200527398722731E+00,
     9   -.2734298276726888E-01,0.1263336502158476E+00,
     *   -.1612316155398892E+00,0.1493747695350364E+00,
     1   -.1157722794073368E+00,0.7922693804316655E-01,
     2   -.4922495824156513E-01,0.2824122832910955E-01,
     3   -.1513160218000067E-01,0.7633244756318270E-02,
     4   -.3647515907318666E-02,0.1658822852371838E-02,
     5   -.7206974472728137E-03,0.3000451746157989E-03,
     6   -.1200053909774860E-03,0.4620783688694020E-04,
     7   -.1715966927475209E-04,0.6155204033137647E-05,
     8   -.2135416639206045E-05,0.7173246215016878E-06,
     9   -.2335393707607576E-06,0.7375227944506555E-07,
     *   -.2260843655926461E-07,0.6731455395775903E-08,
     1   -.1947667935302690E-08,0.5478663235406307E-09,
     2   -.1498814715742675E-09,0.3988791386728965E-10,
     3   -.1033055112886422E-10,0.2602413441317869E-11,
     4   -.6395149238894020E-12,0.1514948612746738E-12,
     5   -.3644140026561549E-13,0.6973409611298938E-14,
     6   -.2836020732144640E-14,-.7196930566861258E-15,
     7   -.1198106474848615E-14,-.1120446016908565E-14,
     8   -.1152018796733228E-14,-.1156174770767301E-14,
     9   -.1160415758856677E-14,-.1159656464250984E-14,
     *   -.1155174381315158E-14,-.1147083912841387E-14/
c 
        data c16/
     1   -.1145153049283855E+02,-.4017233457776307E+01,
     2   0.2173452147557832E+01,0.1646343424457237E+01,
     3   -.1246341118456992E+01,-.8071284207053313E+00,
     4   0.9478464700555034E+00,0.2301996108929375E+00,
     5   -.7128858080859864E+00,0.2097476417092779E+00,
     6   0.3496237329368914E+00,-.3959932295767214E+00,
     7   0.8676270499147922E-01,0.1980653234898027E+00,
     8   -.2680062093665335E+00,0.1605476960638999E+00,
     9   0.3431078796072744E-03,-.1211911166896989E+00,
     *   0.1720311517263049E+00,-.1665771048551477E+00,
     1   0.1323320231273439E+00,-.9193467434153519E-01,
     2   0.5764683710349547E-01,-.3324100011161576E-01,
     3   0.1784552467091676E-01,-.8997689338258473E-02,
     4   0.4288428584924116E-02,-.1941803478260649E-02,
     5   0.8386402441989569E-03,-.3465806083963781E-03,
     6   0.1374162884144840E-03,-.5238816058148455E-04,
     7   0.1923945450624678E-04,-.6817072097935465E-05,
     8   0.2333610285900675E-05,-.7726438578122213E-06,
     9   0.2476706828274869E-06,-.7692652387041688E-07,
     *   0.2316821218476102E-07,-.6769944497794809E-08,
     1   0.1920323465234400E-08,-.5289806982910325E-09,
     2   0.1415597243395943E-09,-.3680889630426931E-10,
     3   0.9305112917303649E-11,-.2283861854021880E-11,
     4   0.5477289404178580E-12,-.1248376753691083E-12,
     5   0.3058829291583056E-13,-.4239665374675407E-14,
     6   0.3430465358519274E-14,0.1859082221830838E-14,
     7   0.2228923114152798E-14,0.2189642216748260E-14,
     8   0.2221393332870041E-14,0.2229118519414440E-14,
     9   0.2231598512802297E-14,0.2225881126087671E-14,
     *   0.2213148594414772E-14,0.2193880793436374E-14/
c 
        data c17/
     1   0.7015928004094759E+01,0.9363355647154100E+01,
     2   0.6774353637791607E+00,-.3004566446428884E+01,
     3   -.4572447129686147E+00,0.1705235067841963E+01,
     4   0.1459224646944044E+00,-.1117244823660107E+01,
     5   0.1891024623807831E+00,0.6769150812316303E+00,
     6   -.4392709438369876E+00,-.2103111077299733E+00,
     7   0.4448833224118880E+00,-.2034026935565267E+00,
     8   -.1328042004158839E+00,0.2783538957538193E+00,
     9   -.2093537288229894E+00,0.4561686130981522E-01,
     *   0.9721361439755677E-01,-.1692844571010902E+00,
     1   0.1760340631428937E+00,-.1453685707012296E+00,
     2   0.1035130821205743E+00,-.6601340311009098E-01,
     3   0.3852594261656902E-01,-.2086368215017770E-01,
     4   0.1058636960237906E-01,-.5068883752506201E-02,
     5   0.2302775178949637E-02,-.9968456718004277E-03,
     6   0.4126130871124710E-03,-.1637676302976891E-03,
     7   0.6247507050013941E-04,-.2295298437626461E-04,
     8   0.8134931532681211E-05,-.2785299965132564E-05,
     9   0.9224117465727525E-06,-.2957796861187239E-06,
     *   0.9191654670098199E-07,-.2770360012908110E-07,
     1   0.8103667820759825E-08,-.2301832610498519E-08,
     2   0.6352023785123744E-09,-.1703648986722859E-09,
     3   0.4441771496151975E-10,-.1126621511512461E-10,
     4   0.2774384654179121E-11,-.6696559615622344E-12,
     5   0.1520205492233422E-12,-.3887844369081924E-13,
     6   0.4114015314342576E-14,-.5441499904124508E-14,
     7   -.3488794727485867E-14,-.3965264615407952E-14,
     8   -.3920143202930986E-14,-.3961404089712982E-14,
     9   -.3967946405203207E-14,-.3964278366033441E-14,
     *   -.3946857401026535E-14,-.3917534725258092E-14/
c 
        data c18/
     1   0.1359544607506475E+02,0.5541759629876964E+01,
     2   -.1998759990997536E+01,-.2548111593110332E+01,
     3   0.8539930982602127E+00,0.1545200121725729E+01,
     4   -.6117435124360614E+00,-.9368356594577183E+00,
     5   0.6173231543114205E+00,0.4368739311482209E+00,
     6   -.6088842716425164E+00,0.2267293391698684E-01,
     7   0.4111419344966029E+00,-.3241704271859493E+00,
     8   -.1726157377863590E-01,0.2520631690621286E+00,
     9   -.2553129464726864E+00,0.1092243189339845E+00,
     *   0.5480904349995471E-01,-.1578044960850336E+00,
     1   0.1868553580645943E+00,-.1655495986729602E+00,
     2   0.1236012000406463E+00,-.8164113693379900E-01,
     3   0.4896708496722036E-01,-.2710377926134312E-01,
     4   0.1399775690934574E-01,-.6798913885135456E-02,
     5   0.3124477442293590E-02,-.1364904357236918E-02,
     6   0.5688994474629065E-03,-.2269324505821784E-03,
     7   0.8685093178116160E-04,-.3195788081415847E-04,
     8   0.1132587743214933E-04,-.3871716555125241E-05,
     9   0.1278271846464094E-05,-.4080367055543363E-06,
     *   0.1260463083179181E-06,-.3770981411823881E-07,
     1   0.1093345332904487E-07,-.3073828112988151E-08,
     2   0.8383284864388077E-09,-.2218907917161369E-09,
     3   0.5700494716564733E-10,-.1422690487915296E-10,
     4   0.3440453071205411E-11,-.8158298025960488E-12,
     5   0.1798875410701062E-12,-.4662829404250092E-13,
     6   0.3203844201173102E-14,-.7629986996291165E-14,
     7   -.5493756672132915E-14,-.6016106191251181E-14,
     8   -.5975157014498893E-14,-.6016746572945930E-14,
     9   -.6016571392672560E-14,-.5999639673446886E-14,
     *   -.5963125171574929E-14,-.5909390115277737E-14/
c 
         data c19/
     1   -.8156211518456592E+01,-.1135363932100233E+02,
     2   -.1657131762682050E+01,0.3426392886635308E+01,
     3   0.1355763610879835E+01,-.1793151374412115E+01,
     4   -.9585650616547209E+00,0.1178356911871239E+01,
     5   0.5332650555796903E+00,-.9006462047779309E+00,
     6   -.1107887249230735E+00,0.6674273570963126E+00,
     7   -.2519409371133615E+00,-.3197162459209478E+00,
     8   0.4058760166510907E+00,-.1000182793249625E+00,
     9   -.2036635600213701E+00,0.2817569718939009E+00,
     *   -.1650406924530142E+00,-.9879684459124777E-02,
     1   0.1382433684830842E+00,-.1881125752638425E+00,
     2   0.1770659714906852E+00,-.1369051044257003E+00,
     3   0.9248011804423154E-01,-.5630326049257351E-01,
     4   0.3147599157612070E-01,-.1635984039794460E-01,
     5   0.7975829348795704E-02,-.3671498119231175E-02,
     6   0.1604007241400053E-02,-.6677801479302468E-03,
     7   0.2658027373322330E-03,-.1014299640767091E-03,
     8   0.3719118941856447E-04,-.1312838459617479E-04,
     9   0.4468732840542421E-05,-.1468791115059318E-05,
     *   0.4667154149545709E-06,-.1435162188623603E-06,
     1   0.4274485509887953E-07,-.1234020631893638E-07,
     2   0.3455353144936673E-08,-.9388977706255463E-09,
     3   0.2476952264714738E-09,-.6345242526177933E-10,
     4   0.1580316192863073E-10,-.3812187048672816E-11,
     5   0.9063202802946449E-12,-.1965171261332786E-12,
     6   0.5441867853189663E-13,-.7636484833301298E-15,
     7   0.1129787165923233E-13,0.8944592177299925E-14,
     8   0.9540227984956619E-14,0.9496094913918073E-14,
     9   0.9534983716231492E-14,0.9518009066789125E-14,
     *   0.9474181282666398E-14,0.9401044049576369E-14/
c 
        data c20/
     1   -.1593387221591375E+02,-.7165063372625038E+01,
     2   0.1714436996293009E+01,0.3363889413657140E+01,
     3   -.2782434596246360E+00,-.2067162400554984E+01,
     4   0.3803472744950500E-02,0.1377750203443237E+01,
     5   -.1037332793049485E+00,-.9332529580843546E+00,
     6   0.3175601770403527E+00,0.5395498672802699E+00,
     7   -.4774500025139742E+00,-.1119905693442321E+00,
     8   0.4203837151574840E+00,-.2471775880570949E+00,
     9   -.9056935121765201E-01,0.2732100630256711E+00,
     *   -.2293355595876059E+00,0.6524054718771973E-01,
     1   0.9161888388790826E-01,-.1769092654136898E+00,
     2   0.1897346656150921E+00,-.1587786315228070E+00,
     3   0.1135487883565494E+00,-.7228553056565220E-01,
     4   0.4191581517292563E-01,-.2246694728129432E-01,
     5   0.1124544849717342E-01,-.5295680681776912E-02,
     6   0.2359689322392821E-02,-.9993561463904818E-03,
     7   0.4037191217209779E-03,-.1560294814722264E-03,
     8   0.5783046893503979E-04,-.2059729212039715E-04,
     9   0.7061674820014733E-05,-.2333864311553288E-05,
     *   0.7444651936561560E-06,-.2294372177930742E-06,
     1   0.6837774469767886E-07,-.1972047127437895E-07,
     2   0.5507329335046269E-08,-.1490052895144145E-08,
     3   0.3907522598514499E-09,-.9933286448865868E-10,
     4   0.2450506439665299E-10,-.5846307764564711E-11,
     5   0.1370749030197851E-11,-.2936653866189503E-12,
     6   0.7924191081925542E-13,-.1413996305339011E-14,
     7   0.1583575929462441E-13,0.1250403101941887E-13,
     8   0.1329446182769980E-13,0.1321174007033953E-13,
     9   0.1323888199271929E-13,0.1319226737167030E-13,
     *   0.1310965790924587E-13,0.1298839803156286E-13/
c 
c        return the matrix of coefficiets of universal expansions
c        for the 20 basis functions
c 
        do 1400 i=1,20
        do 1200 j=1,120
c 
        coefs(j,i)=0
 1200 continue
 1400 continue
c 
        call o1revenm(c1,coefs(1,1),44)
        call o1revenm(c3,coefs(1,3),46)
        call o1revenm(c5,coefs(1,5),47)
        call o1revenm(c7,coefs(1,7),49)
        call o1revenm(c9,coefs(1,9),50)
        call o1revenm(c11,coefs(1,11),51)
        call o1revenm(c13,coefs(1,13),52)
        call o1revenm(c15,coefs(1,15),60)
        call o1revenm(c17,coefs(1,17),60)
        call o1revenm(c19,coefs(1,19),60)
c 
        call o1roddmv(c2,coefs(1,2),43)
        call o1roddmv(c4,coefs(1,4),46)
        call o1roddmv(c6,coefs(1,6),47)
        call o1roddmv(c8,coefs(1,8),48)
        call o1roddmv(c10,coefs(1,10),48)
  
        call o1roddmv(c12,coefs(1,12),50)
        call o1roddmv(c14,coefs(1,14),58)
        call o1roddmv(c16,coefs(1,16),60)
        call o1roddmv(c18,coefs(1,18),60)
        call o1roddmv(c20,coefs(1,20),60)
c 
        nuniv=120
        n=20
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o1revenm(arrin,arrout,nmove)
        implicit real *8 (a-h,o-z)
        save
        dimension arrin(1),arrout(1)
c 
        do 1200 i=1,nmove
c 
        arrout(i*2-1)=arrin(i)
 1200 continue
        return
c 
c 
c 
c 
        entry o1roddmv(arrin,arrout,nmove)
  
c 
        do 1400 i=1,nmove
c 
        arrout(i*2)=arrin(i)
 1400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine o1rmatav(a,x,y,n)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),x(n),y(n)
c 
        do 1400 i=1,n
        d=0
        do 1200 j=1,n
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
        subroutine o1rinmat(b,a,n)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension a(20,20),part1(100),part2(100),
     1      part1b(10,10),part2b(10,10),b(20,20)
c 
c        This subroutine returns to the user the matrices a, b
c        of a singular function expansion and its inverse. Specifically,
c        the i-th column of b contains the values at the 20 support
c        nodes of the 20 singular functions of representing on the
c        interval [-1,1] all functions of the forms
c 
c 
c        P_n(x),
c 
c        P_n (x) \cdot (1-x) \cdot log(1-x),
c                                                                   (1)
c        P_n (x) \cdot  \cdot(1+x) log(1+x),
c 
c        with n=0,1,2,...,9.
c 
c        The matrix a is the inverse of the matrix b. In other words,
c        the matrix b converts the coefficients of a singular function
c        expansion into the values of that expansion at the support nodes;
c        the matrix a converts the values of a singular function expansion
c        at the support node into the coefficients of that expansion. The
c        accuracy of the obtained expansion is close to 1.0E-7, making it
c        the weakest link among the subroutines in this file.
c 
        data part1b/
     1   0.1136085858789022E+01,0.1125810313747931E+01,
     2   0.1093970859454983E+01,0.1031230125519136E+01,
     3   0.9379948595934643E+00,0.8249502943438656E+00,
     4   0.7089146660363162E+00,0.6070404504875718E+00,
     5   0.5326144413720006E+00,0.4937026909514003E+00,
     6   -.1414445301077430E+01,-.1403421060715730E+01,
     7   -.1361693496177413E+01,-.1269119582206360E+01,
     8   -.1119130411484244E+01,-.9228182697756373E+00,
     9   -.7036173253573554E+00,-.4852640467029058E+00,
     *   -.2810216812044055E+00,-.9160823507012686E-01,
     1   0.2087909309350799E+01,0.2002161136640564E+01,
     2   0.1743173110443594E+01,0.1271238124859312E+01,
     3   0.6707319833804164E+00,0.1027715413020761E+00,
     4   -.3091570540024252E+00,-.5416053836433002E+00,
     5   -.6429191626626658E+00,-.6745590013405117E+00,
     6   -.2750923283343635E+01,-.2563446209398808E+01,
     7   -.2016987317617394E+01,-.1103000684971966E+01,
     8   -.1237105927182900E+00,0.5433038172717423E+00,
     9   0.7620358213568954E+00,0.6542833779261368E+00,
     *   0.4072446498160692E+00,0.1349484386207639E+00,
     1   0.2994814702789162E+01,0.2652141005277476E+01,
     2   0.1688567255164811E+01,0.2560968323751956E+00,
     3   -.8523713305318701E+00,-.1014485980478503E+01,
     4   -.4921126887896030E+00,0.8017562844721874E-01,
     5   0.4270422957467016E+00,0.5762208600946078E+00,
     6   -.3140354239393873E+01,-.2627055663928265E+01,
     7   -.1253586004701550E+01,0.5072425212725544E+00,
     8   0.1295215170446730E+01,0.6255842178331662E+00,
     9   -.3624962969352983E+00,-.6872354376251592E+00,
     *   -.5125136042952411E+00,-.1862573694397564E+00,
     1   0.2726445439086225E+01,0.2164325783925646E+01,
     2   0.7076414394595347E+00,-.9284466081177413E+00,
     3   -.1103068578739196E+01,0.2190420479787183E+00,
     4   0.1003133748643821E+01,0.4401681164097826E+00,
     5   -.3067443950832906E+00,-.5737899941843410E+00,
     6   -.2297271537797024E+01,-.1715837489400672E+01,
     7   -.2963260654183207E+00,0.1007686553761522E+01,
     8   0.5894748401934927E+00,-.8272725711577505E+00,
     9   -.7949934365842408E+00,0.4707678722059481E+00,
     *   0.8608906485622397E+00,0.3167287623437119E+00,
     1   0.1699378209926275E+01,0.1236428653585613E+01,
     2   0.1167325884804135E+00,-.8043529413786893E+00,
     3   -.2203828926755801E+00,0.8316365066278385E+00,
     4   0.1781175430455220E+00,-.9922955831792367E+00,
     5   -.4143638259261813E+00,0.8753599102665315E+00,
     6   -.1683273202767422E+01,-.1163457120507366E+01,
     7   0.8455488359506032E-02,0.7793201340734925E+00,
     8   -.1363444559226335E-02,-.7383317255676822E+00,
     9   0.3322603777310904E+00,0.8248527082640303E+00,
     *   -.7149562927993444E+00,-.8813497338716269E+00/
c 
        data part2b/
     1   0.2204748324528185E+01,0.1533576386325688E+01,
     2   -.2313135831389629E-01,-.1015314981320214E+01,
     3   0.1140108494519059E+00,0.8032765137542496E+00,
     4   -.6634703768820053E+00,-.2759113505520982E+00,
     5   0.9998002789548845E+00,-.5760488530499511E+00,
     6   -.3736741028445106E+01,-.2189276987400755E+01,
     7   0.5846227861122341E+00,0.1259889650392143E+01,
     8   -.7311612524622197E+00,-.4088976351357623E+00,
     9   0.8676648567763588E+00,-.5294039023061462E+00,
     *   -.1653992324593807E+00,0.6871153165299934E+00,
     1   0.5380808470922434E+01,0.2470697225514091E+01,
     2   -.1320223313739641E+01,-.1000212441865773E+01,
     3   0.1113687373369881E+01,-.2699209065745673E+00,
     4   -.4119202130410704E+00,0.6414808076063617E+00,
     5   -.5136379210479249E+00,0.1910428443834634E+00,
     6   -.7171075520389010E+01,-.1703052999767166E+01,
     7   0.2104266128663180E+01,-.1303445973701510E+00,
     8   -.7932332039553352E+00,0.8069941731138706E+00,
     9   -.4459518445129567E+00,0.1861764727301820E-01,
     *   0.3210466965353083E+00,-.5042768421623136E+00,
     1   0.8424876754959929E+01,0.3714157161231233E+00,
     2   -.2116773324712727E+01,0.1073389908317203E+01,
     3   -.9748218964717866E-01,-.4259638193109584E+00,
     4   0.5960646340584487E+00,-.5427993334082303E+00,
     5   0.3649863992953931E+00,-.1274560896965927E+00,
     6   -.8390969425641760E+01,0.1345010992672657E+01,
     7   0.1151012801611762E+01,-.1318003905746022E+01,
     8   0.9213533352210087E+00,-.4536959734456665E+00,
     9   0.5619091613746015E-01,0.2377166551415328E+00,
     *   -.4281140896199655E+00,0.5208452973428596E+00,
     1   0.8018878186615089E+01,-.2723184929566060E+01,
     2   0.2341505920449578E+00,0.6084433612009233E+00,
     3   -.8308658759531968E+00,0.7983695534519644E+00,
     4   -.6617106564839767E+00,0.4844013497110868E+00,
     5   -.2930312636858447E+00,0.9788012934311653E-01,
     6   0.7157773917260531E+01,-.3356242938028188E+01,
     7   0.1364652956080329E+01,-.4586216249341102E+00,
     8   -.4243972020841054E-02,0.2589677333969544E+00,
     9   -.4064021350079253E+00,0.4930966976761602E+00,
     *   -.5417933912143087E+00,0.5637672872461168E+00,
     1   -.6070601775039853E+01,0.3417252171651979E+01,
     2   -.1976952048827501E+01,0.1264452747475735E+01,
     3   -.8559268517728747E+00,0.5938030950649698E+00,
     4   -.4092084078404560E+00,0.2683478218109450E+00,
     5   -.1524240452854454E+00,0.4946512621040585E-01,
     6   -.4849706325102530E+01,0.2865377935518491E+01,
     7   -.1790104392637180E+01,0.1265556438381244E+01,
     8   -.9777177768066870E+00,0.8074280791539377E+00,
     9   -.7023069590239690E+00,0.6371201863815820E+00,
     *   -.5989685298033215E+00,0.5812871847206438E+00/
c 
        data part1/
     1   0.2481424541874686E-02,-.3089398056980311E-02,
     2   0.4560409818448364E-02,-.6008593145291041E-02,
     3   0.6541299855500775E-02,-.6859195193566409E-02,
     4   0.5955098269948300E-02,-.5017634728750107E-02,
     5   0.3711747356647201E-02,-.3676435209360597E-02,
     6   0.4815617878920846E-02,-.8160412682930520E-02,
     7   0.1175072119485147E-01,-.1567360645373594E-01,
     8   0.1844468299642080E-01,-.1848606046727134E-01,
     9   0.1771956523903071E-01,0.1542651013555466E-01,
     *   -.1222192682524448E-01,-.6478745932096129E-02,
     1   0.1242395566315707E-01,-.1548759261170247E-01,
     2   0.2209491671188387E-01,-.2828891013740242E-01,
     3   0.2926770265770080E-01,-.2899083163745936E-01,
     4   0.2388446886499315E-01,-.1893533820428084E-01,
     5   0.1364471868816317E-01,-.1283983809257694E-01,
     6   0.1692377669721253E-01,-.2416380888384870E-01,
     7   0.2727094333556795E-01,-.1876330951052575E-01,
     8   0.3984063191725858E-02,0.1527335213636591E-01,
     9   -.3053870356107755E-01,-.3614881484029446E-01,
     *   0.3441782943437216E-01,0.1934648173893269E-01,
     1   0.3142872513917931E-01,-.3912006677371756E-01,
     2   0.5007977337354214E-01,-.5794635356509160E-01,
     3   0.4851115704121267E-01,-.3601464971194304E-01,
     4   0.2032998907673131E-01,-.8512999797252748E-02,
     5   0.3353576532470814E-02,0.2436122508563179E-03,
     6   -.6643785590032050E-03,0.1680197644734060E-01,
     7   -.3793659108821601E-01,0.6040446537485350E-01,
     8   -.6065429469180030E-01,0.3244498744377431E-01,
     9   0.7290824889665952E-02,0.3723871313902255E-01,
     *   -.5130591134885919E-01,-.3189135248463217E-01,
     1   0.5584056334644521E-01,-.6872224715410856E-01,
     2   0.6883676841887556E-01,-.5972660041307370E-01,
     3   0.1386724715757511E-01,0.2746730717731786E-01,
     4   -.5027501386589216E-01,0.5456544554036038E-01,
     5   -.4355522461208164E-01,0.4219882281673393E-01,
     6   -.5497894745411987E-01,0.6821378074571650E-01,
     7   -.5415173854234991E-01,-.6992628549657585E-02,
     8   0.5794532149498026E-01,-.7061605269262079E-01,
     9   0.3242439011702643E-01,-.2164580320626895E-01,
     *   0.6147662498472979E-01,0.4365524046543624E-01,
     1   0.7914397516722311E-01,-.9442729259891179E-01,
     2   0.5659356695742964E-01,-.1043847601881895E-01,
     3   -.7191920670918011E-01,0.1092842569805610E+00,
     4   -.9307207316455153E-01,0.4973764508208173E-01,
     5   -.1859500653425773E-01,-.1139060634796284E-03,
     6   0.9619971350744342E-02,-.6168175434227817E-01,
     7   0.9395864669254819E-01,-.6700960416144671E-01,
     8   -.8047171496247883E-02,0.7689975702741384E-01,
     9   -.6967578598269444E-01,-.4731676716172609E-02,
     *   -.6463921945324910E-01,-.5424964409556506E-01/
  
c 
        data part2/
     1   0.9569474791254722E-01,-.1070476430686823E+00,
     2   0.1192150541716174E-01,0.6302388542580319E-01,
     3   -.1176811932231306E+00,0.7256866388210932E-01,
     4   0.2540886692641944E-01,-.9596441687623673E-01,
     5   0.9647035372801861E-01,-.8564822174124884E-01,
     6   0.9318032848973431E-01,-.4744457162075896E-01,
     7   -.3130200223794191E-01,0.9370447871611642E-01,
     8   -.4957518626528894E-01,-.5172895339596844E-01,
     9   0.9228418338604099E-01,0.3546333448046599E-01,
     *   0.6156728930267241E-01,0.6350782724806742E-01,
     1   0.1033264154397505E+00,-.1025541359201680E+00,
     2   -.4506052744132011E-01,0.1110686720394982E+00,
     3   -.7172671629578200E-01,-.5283549803459590E-01,
     4   0.1462098804964328E+00,-.1158723351301129E+00,
     5   0.2596114089637135E-01,0.4842940234268293E-01,
     6   -.9670261311639128E-01,0.1264783662554116E+00,
     7   -.6004630744302215E-01,-.6510220521470879E-01,
     8   0.8701558561591962E-01,0.7247869429276138E-02,
     9   -.9621212208083171E-01,-.6552495293807953E-01,
     *   -.5328501181892630E-01,-.7123837912849463E-01,
     1   0.1036683061841981E+00,-.8287194746310948E-01,
     2   -.9249353870540890E-01,0.1117366190225787E+00,
     3   0.1369197878641737E-01,-.1173631004324009E+00,
     4   0.7517030499888842E-01,0.8039574312608407E-01,
     5   -.1694608070569721E+00,0.1408639462640181E+00,
     6   -.4711933839312827E-01,-.9042445622401462E-01,
     7   0.1095557570139962E+00,0.3291215410737904E-02,
     8   -.9280099350035366E-01,0.4156760285590931E-01,
     9   0.8256951337772698E-01,0.9116384436826193E-01,
     *   0.4093874985244449E-01,0.7724720226323167E-01,
     1   0.1005524991471716E+00,-.5305397633509818E-01,
     2   -.1213769712425067E+00,0.7688353864297842E-01,
     3   0.8062158790232987E-01,-.9675829958953501E-01,
     4   -.5791034909846494E-01,0.1625283242545060E+00,
     5   -.7822795200589350E-01,-.1349752763368756E+00,
     6   0.1887528366687529E+00,-.3121039135132966E-01,
     7   -.9697360007698557E-01,0.6049303316466835E-01,
     8   0.6897024346953332E-01,-.8181345738753121E-01,
     9   -.5523385172261883E-01,-.1096880372442434E+00,
     *   -.2570766542577839E-01,-.8136116825412870E-01,
     1   0.9785455792344338E-01,-.1815748420333277E-01,
     2   -.1337012616743720E+00,0.2674787325764243E-01,
     3   0.1142100730089953E+00,-.3691657713757417E-01,
     4   -.1137283112685119E+00,0.6277694562306875E-01,
     5   0.1735010959553567E+00,-.1746899546346138E+00,
     6   -.1141760532330455E+00,0.1361741956432088E+00,
     7   0.3786698403404952E-01,-.9982994430136035E-01,
     8   -.2528428356115602E-01,0.1042329482059855E+00,
     9   0.1937207144813469E-01,0.1193700063945465E+00,
     *   0.8759239924804058E-02,0.8345168860387143E-01/
c 
        n=20
c 
c       return matrix b
c 
        ss=1
c 
        do 1400 i=1,n/2
c 
        do 1200 j=1,n/2
        b(j,i)=part1b(j,i)
        b(n-j+1,i)=ss*b(j,i)
 1200 continue
c 
        ss=-ss
c 
 1400 continue
c 
        ss=1
c 
        do 2400 i=1,n/2
c 
        do 2200 j=1,n/2
        b(j,i+n/2)=part2b(j,i)
        b(n-j+1,i+n/2)=ss*b(j,i+n/2)
 2200 continue
c 
        ss=-ss
c 
 2400 continue
c 
c       return matrix a
c 
        call o1rarrcp(part1,a,n*n/4)
        call o1rarrcp(part2,a(1,n/4+1),n*n/4)
c 
        call o1rarsil(a(1,1),a(1,20),n)
        call o1rarsil(a(1,2),a(1,19),n)
        call o1rarsil(a(1,3),a(1,18),n)
        call o1rarsil(a(1,4),a(1,17),n)
        call o1rarsil(a(1,5),a(1,16),n)
        call o1rarsil(a(1,6),a(1,15),n)
        call o1rarsil(a(1,7),a(1,14),n)
        call o1rarsil(a(1,8),a(1,13),n)
        call o1rarsil(a(1,9),a(1,12),n)
        call o1rarsil(a(1,10),a(1,11),n)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine o1rarsil(a,b,n)
        implicit real *8 (a-h,o-z)
        save
        dimension a(1),b(1)
c 
        c=1
        do 1200 i=1,n
c 
        b(i)=a(i)*c
        c=-c
 1200 continue
        return
        end
