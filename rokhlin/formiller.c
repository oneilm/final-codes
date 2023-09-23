        implicit real *8 (a-h,o-z)
        real *8 w(200 000),khi(10000),rlams(10000),rlams2(10000),
     1      coefs(1000),fs(100 000),ders(100 000)
c
        call prini(6,13)
C
C       SET ALL PARAMETERS
C
        PRINT *, 'ENTER c'
        READ *,c
        CALL PRIN2('c=*',c,1 )
C
        PRINT *, 'ENTER k'
        READ *,k
        CALL PRINf('k=*',k,1 )
c
c       find the eigenvectors
c
        lenw=200 000
c
        call prolcrea(ier,c,w,lenw,nvects,nhigh,
     1    khi,rlams,rlams2,keep,lused)
c
        call prin2('after prolcrea khi are*',khi,nvects)
        call prin2('after prolcrea, rlams2=*',rlams2,nvects)
        call prin2('after prolcrea, rlams=*',rlams,nvects*2)
c
        call prinf('ier=*',ier,1)
        call prinf('nhigh=*',nhigh,1)
        call prinf('keep=*',keep,1)
        call prinf('lused=*',lused,1)
        call prinf('nvects=*',nvects,1)
c
c       test the subroutine for the evaluation of prolate 
c       spheroidal wave functions
c
        delta=1.0d-4
        x=1-delta
cccc        x=0
        
        call proleva(k,x,w,val,der)
c
        call prin2('after proleva, val=*',val,1)        
        call prin2('after proleva, der=*',der,1)
c
c        retrieve from array w the Legendre coefficients 
c        of the -th prolate wave function
c
cccc        k=2
        call prolunpk(k,w,coefs,n)
c
        call prin2('and coefficients as retrieved*',coefs,n)
c
c        evaluate the function at the same point x via the obtained expansion
c
        call legeFDER(X,VAL2,der2,coefs,n)
c
        call prin2('and val2 from legefder is*',val2,1)
        call prin2('and der2 from legefder is*',der2,1)

        call prin2('and difference in val *',val2-val,1)
        call prin2('and difference in der *',der2-der,1)
c
c        plot the obtained spheroidal wave functions, and
c        otherwise test the junk
c
        nn=nhigh+10
ccccc        nn=nhigh*10
        call proltest(w,nvects,c,nhigh,nn,fs,ders)

        stop
        end
c
c
c
c
c
        subroutine proltest(w,nvects,c,nhigh,nn,fs,ders)
        implicit real *8 (a-h,o-z)
        dimension tt(1000),
     1      fs(nn,1),ders(nn,1),ww(1000)
c
c       discretize the interval [-1,1]
c
        itype=1
        call legeexps(itype,nn,tt,u,v,ww)
c
cccc        call prin2('in proltes, tt as created*',tt,nn)
c
        iw=100
c
        do 2000 k=0,nvects-1
c
c       construct the table of values of the k-th function
c
        do 1400 j=1,nn
c
        call proleva(k,tt(j),w,fs(j,k+1),ders(j,k+1) )
 1400 continue
c
c        plot this function
c
        iw=iw+1



        call lotagraph(iw,tt,fs(1,k+1),nn,'the prolate wave function*')

cccc        call RSGRAF(tt,fs(1,k+1),nN,IW,'the prolate wave function*')
c
 2000 continue
c
c        check the inner products between them things
c 
        errmax=0
        do 2600 i=1,nvects
c
cccc        call prinf('i=*',i,1)       
c
        do 2400 j=1,nvects
c
        d=0
        do 2200 k=1,nn
        dd=fs(k,i)*fs(k,j) * ww(k)
        d=d+dd
 2200 continue
c
        ddd=dabs(d)
        if(i .eq. j) ddd=dabs(d-1)
        if(ddd .gt. errmax) errmax=ddd
c
cccc        f(j)=d 
c
 2400 continue
c
cccc        call prin2('and inner products of i-th function*',
cccc     1      f,nvects)
 2600 continue
c
        call prin2('in proltest, maximum orthonormality error is=*',
     1      errmax,1)
c
        return
        end
c
c
c
c
c
        subroutine amatvec(a,x,n,w,y)
        implicit real *8 (a-h,o-z)
        complex *16 a(n,n),x(1),y(1),d
        dimension w(1)
c
        do 2000 i=1,n
        d=0
        do 1800 j=1,n
        d=d+a(i,j)*x(j)*w(j)
 1800 continue
        y(i)=d
 2000 continue
        return
        end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c      this is the end of the debugging code and the beginning of
c      the code for the actual evaluation of the prolate spheroidal
c      wave functions
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
        subroutine prolcrea(ier,c,w,lenw,nvects,nhigh,
     1    khi,rlams,rlams2,keep,lused)
        implicit real *8 (a-h,o-z)
        real *8 khi(1),rlams2(1),w(1)
        complex *16 rlams(1),ima
        data ima/(0.0d0,1.0d0)/
c
        dimension ns(20),nvectss(20)
        data nvectss/25,34,43,50,58,65,72,80,87,94,100,107,114,121,
     1       128,134,141,148,154,161/
        data ns/48,64,80,92,106,120,130,144,156,168,
     1       178,190,202,214,224,236,248,258,268,280/       
c
c        this subroutine constructs the eigenfunctions of the 
c        operator 
c
c         F(\phi) (x) = \int_{-1}^1  \phi (t) * e^{i * c *x * t) dt    (1)
c
c        corresponding to all eigenvalues of (1) whose absolute values
c        are greater than 1.0d-15 (or so). It also produces the eigenvalues
c        of (1) (all of them whose absolute values are greater than 
c        1.0d-15), in decreasing order. In addition, it produces all
c        eigenvalues of the operator
c        
c
c
c         G(\phi) (x) = (1 / \pi) *  \int_{-1}^1  \phi (t) * 
c                                                                      (2)
c         sin(c*(x-t))/(x-t) dt,
c
c         also in decreasing order. The actual evaluation of the 
c         eigenfunctions of (1) (which are also the eigenfunctions 
c         of (2)) is performed by the subroutine prolev0 (see); 
c         prolev0 uses Legendre expansions of these eigenfunctions 
c         created by this subroutine and stored in array w, in a form
c         not readily accessible by the user.
c
c         The subroutine uses the fact that the eigenfunctions of the 
c         operator (1) are also the eigenfunctions of the Sturm-Liouville
c         problem 
c
c      D               D
c      -  ( (1-x**2) * -  \phi (x) ) + (khi-c**2*x**2) * \phi(x) = 0.    (3)
c      Dx              Dx
c 
c        Note that the equation (3) defines the prolate speheroidal 
c        wave functions. The coefficients khi are the so-called 
c        separation coefficients; these are also calculated by this 
c        subroutine and returned to the user. The only other 
c        user-accesible subroutine in this collection are prolunpk
c        and proleva (see both below). Proleva evaluates at a point
c        x \in [-1,1] the user-specified prolate function and its 
c        derivative. Prolunpk returns to the user the coefficients 
c        of the Legendre expansion of the user-specified prolate
c        function.
c
c                          input parameters:
c
c  c - the coefficient (real) in the formula (1)
c  lenw - the amount of storage space provided in the array lw; must
c       be sufficiently large; if it is too small, the error code ier
c       is set to 4, and the execution of the subroutine is terminated.
c       Note that this parameter, normally, does not need to be very large.
c       
c                          output parameters:
c
c  ier - error return code. 
c         ier=0 means successful execution.
c         ier=4 and ier=8 means that the amount of storage space provided in
c                 the array store is insufficient (see input 
c                 parameter lenw above). This is a fatal error.
c         ier=1024 means that the subroutine prolql1 (a standard 
c                 prehistoric subroutine from eispack) has failed
c                 to find the spectrum of a certain tridiagonal matrix;
c                 this has never happened yet, and probably means that memory
c                 has been mismanaged by the user. This is a fatal error.
c  w - the storage area where all the information is stored to be 
c       used by the subroutine prolev0 to evaluate the prolate spheroidal
c       wave functions and their derivatives
c  nvects - the number of solutions of the equation (1) whose Legendre 
c       expansions have been computed
c  nhigh - the highest order of the Legendre expansion used to express 
c       any of the nvects prolate functions.
c  khi - the coefficients in (2) (nvects of them) at which the equation
c       (3) has a non-trivial solution (the "separation coefficients"
c       for the prolate spheroidal wave function).
c  rlams - the eigenvalues (complex) of the operator (1)
c  rlams2 - the eigenvalues of the operator (2)
c  keep - the number of elements of the array w that has to be 
c       unchanged between the call to this subroutine and the subsequent
c       calls to the subroutine prolev0 (see)
c
c  lused - the number of elements of array w that have been used by 
c       this subroutine (it is somewhat greater than keep)
c
c
c        . . . find out how many essentially non-zero eigenvalues 
c              correspond to the user-specified c, and how many
c              Legendre polynomials to choose in the expansions of
c              the eigenvectors
c
         eps=1.0d-15
         n=c*3
         n=n/2
         nvects=c
c
         i=c/10
         if(i .le. 19) n=ns(i+1)
         if(i .le. 19) nvects=nvectss(i+1)

cccc         call prinf('n as determined in prolcrea is*',n,1)                 
cccc         call prinf('nvects as determined in prolcrea is*',nvects,1)
c
c        allocate memory for the subroutine prolvect that 
c        obtains coefficients of Legendre expansions of prolate
c        spheroidal wave functions
c
        ias=1
        las=n/2+6
c
        ibs=ias+las
        lbs=n/2+6
c
        ics=ibs+lbs
        lcs=n/2+6
c
        ixk=ics+lcs
        lxk=n/2+6
c
        iu=ixk+lxk
        lu=n/2+6
c
        iv=iu+lu
        lv=n/2+6
c
        iw=iv+lv
        lw=n/2+6
c
        iladdr=iw+lw
        lladdr=nvects*4
c
        istore=iladdr+lladdr
c
        lleft=lenw-istore
c
        if(lleft .gt. 200) goto 1400
        ier=8
        return
 1400 continue
c
c       obtain the coefficients of Legendre expansions
c       for the first nvects prolate spheroidal wave functions
c
      call prolvect(ier,n,c,w(ias),w(ibs),w(ics),
     1     w(ixk),khi,nvects,w(istore),w(iladdr),nstore,eps,
     2     w(iu),w(iv),w(iw),lleft,nhigh)
c
        lused=istore+nstore
        if(ier .ne. 0) return
c
cccc        call prinf('in prolcrea after prolvect, ier=*',ier,1)
c
c       perform garbage collection
c
        do 2200 i=1,nvects*4
        w(i)=w(iladdr+i-1)
 2200 continue
c
        iladdr=1
c
        istore2=iladdr+nvects*4
c
        do 2400 i=1,nstore+4
        w(istore2+i-1)=w(istore+i-1)
 2400 continue
        istore=istore2
c
        keep=istore+nstore+4
cccc        call prinf('in prolcrea, keep=*',keep,1)
cccc        call prinf('in prolcrea, istore=*',istore,1)
cccc        call prinf('in prolcrea, nstore=*',nstore,1)
c
c        allocate memory for the subroutine prollam0 that will
c        find the spectrum of the integral operator
c
        nn=n+10
        ipexp=keep+1
        lpexp=nn+4
c
        ix=ipexp+lpexp
        lx=nn+2           
c
        iwhts=ix+lx
        lwhts=nn+2
c
c        calculate the largest eigenvalue of the integral operator
c
c         F(\phi) (x) = \int_{-1}^1  \phi (t) * e^{i * c *x * t) dt
c
        call prollam0(nn,c,w(istore),w(iladdr),
     1      w(ipexp),w(ix),w(iwhts),nhigh,rlam0)
c
c       find the absolute values of first nvects eigenvalues 
c       of the integral operator
c
c         F(\phi) (x) = \int_{-1}^1  \phi (t) * e^{i * c *x * t) dt
c
c
        irats77=keep+10
        lrats77=nvects+5
c
        icoefsk=irats77+lrats77
        lcoefsk=nn
c
        icoefskp1=icoefsk+lcoefsk
        lcoefskp1=nn
c
        iderk=icoefskp1+lcoefskp1
        lderk=nn
c
        iderkp1=iderk+lderk
        lderkp1=nn
c
        iscales=iderkp1+lderkp1
        lscales=nhigh+10
c
        ltot=iscales+lscales
        if(ltot .gt. lused) lused=ltot
        if(lused .lt. lenw) goto 3200
        ier=4
        return
 3200 continue
c
        call prolrat(c,w(istore),w(iladdr),nvects,rlam0,
     1      rlams2,nhigh,w(irats77),
     2      w(icoefsk),w(icoefskp1),w(iderk),w(iderkp1),w(iscales) )
c
c       find the first nvects eigenvalues of the integral operator
c
c         F(\phi) (x) = \int_{-1}^1  \phi (t) * e^{i * c *x * t) dt
c
        do 3400 i=1,nvects
c
        rlams(i)=ima**(i-1) * rlams2(i)
 3400 continue
c
cccc        call prin2('actual complex eigenvalues are*',rlams,nvects*2)
c
c       calculate the largest nvects eigenvalues of the integral operator
c
c
c         G(\phi) (x) = \int_{-1}^1  \phi (t) * 
c
c         sin(c*(x-t))/(x-t) dt 
c
        done=1
        pi=datan(done)*4
c
        do 3600 i=1,nvects
        d = rlams2(i) **2
        rlams2(i)=d/pi/2 *c
 3600 continue
cccc        call prin2('and rlams2 are*',rlams2,nvects)        

        return
        end
c
c
c
c
c
        subroutine proleva(k,x,w,val,der)
        implicit real *8 (a-h,o-z)
        dimension w(1),arr(20),iarr(4,4)
        equivalence (arr(1),iarr(1))
c
c        this subroutine evaluates the k-th prolate spheroidal
c        wave function \psi_k (to be referred to simply as prolate 
c        function) at the point x; it also calculates the derivative 
c        \psi_k at the point x. Note that the subroutine is stable 
c        for x on the interval [-1,1]. It works in the small neighborhood
c        of the integral [-1,1], but its accuracy deteriorates rapidly
c        as |x| becomes grater than 1. The precision with which the \psi_k
c        is calculated is abiout 15 digits.
c
c                  input parameters:
c
c  k - the order of the prolate function to be evaluated. Note that k
c        must be between zero and nvects, where nvects has been returned 
c        by a prior call to the subroutine prolcrea (see).
c  x - the point where the prolate function is to be evaluated
c  w - the array that has been created by a prior call to the subroutine 
c        prolcrea (see)
c
c                  output parameters:
c
c  val - the value of the k-th prolate function at the point x
c  der - derivative of the k-th prolate function at the point x
c
c        . . . interpret the first few elements of array w containing its map
c        
        do 1200 i=1,20
        arr(i)=w(i)
 1200 continue
c
        iladdr=1
        istore=iarr(1,4)
        nhigh=iarr(1,1)
        iscale=iarr(1,2) + istore-1
        icoefs=iarr(1,3) +istore-1
c
cccc        call prinf('in proleva, istore=*',istore,1)
cccc        call prinf('in proleva, nhigh=*',nhigh,1)
cccc        call prinf('in proleva, iscale=*',iscale,1)
cccc        call prinf('in proleva, icoefs=*',icoefs,1)
c
cccc        call prin2('in proleva, w(istore)=*',w(istore),nhigh)
cccc        call prin2('in proleva, w(iscale)=*',w(iscale),nhigh)
cccc        call prin2('in proleva, w(icoefs)=*',w(icoefs),nhigh)
c
c       evaluate the k-th prolate function and its derivative at 
c       the point x
c
        call prolev0(k,x,w(istore),w(iladdr),val,der,w(icoefs),
     1      nhigh,w(iscale) )
c
        return
        end
c
c
c
c
c
        subroutine prolunpk(k,w,coefs,n)
        implicit real *8 (a-h,o-z)
        dimension w(1),arr(20),iarr(4,4),coefs(1)
        equivalence (arr(1),iarr(1))
c
c        this subroutine returns to the user the Legendre expansion
c        of the k-th prolate spheroidal wave function \psi_k (to be 
c        referred to simply as prolate function). Note that the 
c        expansion is valid for x on the interval [-1,1]. It works in
c        the small neighborhood  of the integral [-1,1], but its accuracy 
c        deteriorates rapidly  as |x| becomes grater than 1. Also, note
c        that the expansion (with the n terms returned) is accurate to 
c        about 15 digits. 
c
c        IMPORTANT NOTE:
c
c   The coefficients of the expansion returned by this subroutine 
c        correspond to a somewhat unusual normalization of the 
c        Legendre polynomials. Specifically, they are scaled so that 
c        the norm of each polynomial on the interval [-1,1] is equal
c        to 1.
c
c                  input parameters:
c
c  k - the order of the prolate function the coefficients of whose 
c        legendre series are to be evaluated. Note that k must be between 
c        zero and nvects, where nvects has been returned by a prior call
c        to the subroutine prolcrea (see).
c  w - the array that has been created by a prior call to the subroutine 
c        prolcrea (see)
c
c                  output parameters:
c
c  coefs - the first n coefficients in the Legendre expansion of the k-th
c         prolate function.
c  n - the number of terms returned in the array coefs
c
c        . . . interpret the first few wlwmwnts of array w containing its map
c        
        do 1200 i=1,20
        arr(i)=w(i)
 1200 continue
c
        iladdr=1
        istore=iarr(1,4)
        nhigh=iarr(1,1)
        iscale=iarr(1,2) + istore-1
        icoefs=iarr(1,3) +istore-1
c
cccc        call prinf('in prolunpk, istore=*',istore,1)
cccc        call prinf('in prolunpk, nhigh=*',nhigh,1)
cccc        call prinf('in prolunpk, iscale=*',iscale,1)
cccc        call prinf('in prolunpk, icoefs=*',icoefs,1)
c
cccc        call prin2('in prolunpk, w(istore)=*',w(istore),nhigh)
cccc        call prin2('in prolunpk, w(iscale)=*',w(iscale),nhigh)
cccc        call prin2('in prolunpk, w(icoefs)=*',w(icoefs),nhigh)
c
c       retrieve from array w the Legendre coefficients of the 
c       k-th prolate function 
c
        call prolunp0(k,coefs,w(istore),w(iladdr),n,w(icoefs))
c
c       scale the array of Legendre coefficients to make 
c       it properly normalized
c       
        do 1400 i=1,n
        coefs(i)=coefs(i)*w(iscale+i-1)
 1400 continue
        return
        end
c
c
c
c
c
        subroutine prolrat(c,store,laddr,nvects,rlam0,
     1      rlams2,nhigh,rats,coefsk,coefskp1,derk,derkp1,scales)
        implicit real *8 (a-h,o-z)
        dimension store(1),laddr(3,1),rlams2(1),coefsk(1),
     1      coefskp1(1),derk(1),derkp1(1),rats(1),scales(1)
c
c        this subroutine evaluates the absolute values of the first 
c        nvects eigenvalues of the operator
c
c         F(\phi) (x) = \int_{-1}^1  \phi (t) * e^{i * c *x * t) dt   (1)
c
c        Note: this subroutine uses parameters store, laddr,nhigh
c        that must have been produced by the subroutine prolvect (see),
c        and the parameter rlam0 that must have been produced by the 
c        subroutine prollam0. This subroutine has no uses as a stand 
c        alone device.
c
c                         input parameters:
c
c  c - the coefficient (real) in the formula (1)
c  nvects - the number of solutions of the equation (1) whose Legendre 
c       expansions are to be computed
c  store - produced by the subroutine prolvect (see).
c  laddr - produced by the subroutine prolvect
c  rlams0 - the largest eigenvalue of the operator (1) (note that 
c        it is real, as well as all other even-numbered eigenvalues)
c  nhigh - the highest order of the legendre expansion for any of
c       the eigenfunctions (as produced by the subroutine prolvect)
c
c                           output parameters:
c
c  rlams2 - the absolute values of the first nvects eigenvalues 
c       of the operator (1)
c  
c                           work arrays:
c
c  rats - must be at least nvects+1 real *8 elements long
c  coefsk,coefskp1,derk,derkp1 - must be at least nhigh+2 real *8 
c       elements each
c
c
c        . . . one eigenfunction after another, construct the ratios
c              of consecutive pairs of eigenvectors
c
cccc          call prinf('entered prolrat, nvects=*',nvects,1)
        done=1
        do 1100 i=1,nhigh
        scales(i)=dsqrt(2*i-done)
 1100 continue
c
        do 2000 k=1,nvects-1
c
c        retrieve from array store the Legendre expansions for the
c        i-th and i+1-st eigenfunctions
c
        kk=k-1
        call prolunp0(kk,coefsk,store,laddr,nnk,derk)
c
        call prolunp0(k,coefskp1,store,laddr,nnkp1,derk)
c
c
c        rescale both expansion so that the standard differentiation
c        routine would work
c
        do 1200 i=1,nnkp1+4
c
cccc        coefsk(i)=coefsk(i)*dsqrt(2*i-done)
cccc        coefskp1(i)=coefskp1(i)*dsqrt(2*i-done)
        coefsk(i)=coefsk(i)*scales(i)
        coefskp1(i)=coefskp1(i)*scales(i)
 1200 continue
c
c       differentiate both expansions
c
        call legediff(coefsk,nnk,derk)
        call legediff(coefskp1,nnkp1,derkp1)
c
c        rescale back the whole bunch
c
        do 1300 i=1,nnkp1
c
cccc        coefsk(i)=coefsk(i)/dsqrt(2*i-done)
cccc        coefskp1(i)=coefskp1(i)/dsqrt(2*i-done)
c
cccc        derk(i)=derk(i)/dsqrt(2*i-done)
cccc        derkp1(i)=derkp1(i)/dsqrt(2*i-done)
c
        ddd=done/scales(i)
c
        coefsk(i)=coefsk(i)*ddd
        coefskp1(i)=coefskp1(i)*ddd
c
        derk(i)=derk(i)*ddd
        derkp1(i)=derkp1(i)*ddd
 1300 continue
c
c        construct the cross-integrals
c
        nn=nnk
        if(nnkp1 .gt. nn) nn=nnkp1
        nn=nn-1
c
        d1=0
        d2=0
        do 1400 i=1,nn-1
c
        d1=d1+coefsk(i)*derkp1(i)
        d2=d2+coefskp1(i)*derk(i)
 1400 continue
c
        rats(k)=d2/d1
cccc        denoms(k)=d1
        rlams2(k)=d1
 2000 continue
cccc         call prin2('in prolrat, rats=*',rats,nvects)
cccc         call prin2('in prolrat, denominators are=*',rlams2,nvects)
c
c        using the newly obtained ratios and previously obtained 
c        largest eigenvalue of the integral operator, construct the
c        first nvects values of the integral operator
c
        rlams2(1)=rlam0
        do 2200 i=1,nvects-1
        rlams2(i+1)=rlams2(i)*dsqrt(-rats(i))
 2200 continue
cccc        call prin2('and rlams2 in prolrat*',rlams2,nvects)
c
        return
        end
c
c
c
c
c
        subroutine prollam0(n,c,store,laddr,
     1      pexp,x,whts,nhigh,rlam)
        implicit real *8 (a-h,o-z)
        dimension store(1),laddr(4,1),x(1),whts(1),pexp(1)
c
c        starting with the first eigenvector of the integral
c        operator
c
c         F(\phi) (x) = \int_{-1}^1  \phi (t) * e^{i * c *x * t) dt,     (1)
c
c        this subroutine constructs the corresponding eigenvalue
c        of the operator (1). This is done in the dumbest possible 
c        manner; specifically, the subroutine discretizes the interval
c        [-1,1] into Gaussian nodes, and evaluates the integral (1) 
c        for a single point. This point is chosen among the Gaussian 
c        nodes on the interval [-1,1]. The subroutine chooses the node
c        where the eigenvector is the biggest (the obvious choice to 
c        minimize the round-off error).
c        
c                          input parameters:
c
c  n - the number of nodes in the Gaussian discretization of the 
c       interval [-1,1] to be used to evaluate the integral (1)
c  c - the coefficient (real) in the formula (1)
c  store - the array containing the coefficients of the Legendre 
c          expansions of the first nvects non-sero solutions of (1).
c          This array must have been produced by a prior call to 
c          the subroutine prolvect (see)
c  laddr - the map of the array store; this array must have been 
c          produced by a prior call to the subroutine prolvect (see).
c
c                          output parameters:
c
c  rlam - the first eigenvalue of the operator (1)

c                          work arrays:
c
c  x,whts - must be at least n+2 real *8 elements each
c  pexp - must be at least nhigh real *8 elements eash; the parameter
c         nhigh must have been produced by a prior call to the
c         subroutine prolvect (see)
c
c        . . . construct the Gaussian nodes and weights on 
c              the interval [-1,1]
c
        itype=1
        call legeexps(itype,n,x,u,v,whts)
c
c        Construct the first eigenvector at the gaussian nodes,
c        and apply the integral operator to the thing at one point
c
c        . . . construct the eigenfunction
c
        iscale=laddr(1,2)
        imax=n/2
        rat2=0
c
        do 1600 j=1,n
c
cccc        call prolev0(i0,x(j),store,laddr,val,der,pexp,nhigh)
        call prolev0(i0,x(j),store,laddr,val,der,pexp,nhigh,
     1      store(iscale) )
c
        cd=dcos(c*x(imax)*x(j) )
        rat2=rat2+cd*whts(j)*val
c
        if(j .eq. imax) valmax=val
 1600 continue
        rlam=rat2/valmax
cccc        call prin2('in prollam0, rlam=*',rlam,1)

        return
        end
c
c
c
c
c
        subroutine prolunp0(k,coefs,store,laddr,n,w)
        implicit real *8 (a-h,o-z)
        dimension coefs(1),store(1),laddr(4,1),w(1)
c
c       retrieve the legendre coefficients corresponding to the 
c       k-th eigenfunction in their compressed (odd-even) form
c
        call prolret(k,w,store,laddr,n)
c
c        . . . unpack
c
        ii=k/2
        i0=k-ii*2
c
        do 1200 i=1,n*2+2
        coefs(i)=0
 1200 continue
c
        do 1400 i=1,n
        j=(i-1)*2+i0+1
        coefs(j)=w(i)
 1400 continue
c
        n=n*2
        return
        end
c
c
c
c
c
        subroutine prolret(kk,coefs,store,laddr,nn)
        implicit real *8 (a-h,o-z)
        dimension store(1),laddr(4,1),coefs(1)
c
c       retrieve from array store the coefficients of the Legendre
c       expansion of the k-th prolate spherical wave function
c          
        k=kk+1
        i0=laddr(2,k)
        j0=laddr(4,k)
        nj=laddr(3,k)
        nn=j0+nj+1
cccc         call prinf('in prolret, k=*',k,1)
cccc         call prinf('in prolret, i0=*',i0,1)
cccc         call prinf('in prolret, j0=*',j0,1)
cccc         call prinf('in prolret, nj=*',nj,1)
c
        do 1200 i=1,nn
        coefs(i)=0
 1200 continue
cccc          call prin2('coefs after zeroing*',coefs,nn)
c
        do 1400 i=1,nj
        coefs(j0+i-1)=store(i0+i-1)
 1400 continue
        return
        end
c
c
c
c
c
        subroutine prolev0(k,x,w,laddr,val,der,pexp,nhigh,coefs)
        implicit real *8 (a-h,o-z)
        dimension laddr(4,1),w(1),pexp(1),coefs(1)
c
c        this subroutine calculates the values at the point x of the
c        k-th prolate spheroidal function and its derivative. The 
c        subroutine uses the arrays w and laddr that are produced 
c        by the subroutine prolvect (see). Note that the array w in
c        this subroutine is known under the name store in prolvect.
c
c                         input parameters:
c
c  k - the order of the prolate function to be evaluated. note that the 
c       permitted values are 0, 1, ..., nvects-1 (see the subroutine 
c       prolvect for the definition of the parameter nvects).
c  x - the point where the prolate function and its derivative are to
c       be evaluated; should be in the interval [-1,1]
c  w - produced by the subroutine prolvect (see) under the name store
c  laddr - produced by the subroutine prolvect under the same name
c  nhigh - the highest order of the legendre expansion for any of
c       the eigenfunctions (as produced by the subroutine prolvect)
c
c                           output parameters:
c
c  val - the value of the k-th prolate spheroidal wave function at 
c       the point x.
c  der - the value of the derivative of the k-th prolate spheroidal 
c       wave function at the point x.
c
c                           work arrays:
c
c  pexp - must be nhigh real *8 elements long, where lenmax is produced
c       by a preceding call to the subroutine prolvect (see).
c
c
c       . . . expand the chunk of the array w corresponding to 
c             the functions order k, so that it could be fed into
c             the subroutine legefder
c
        k1=k+1
        i=k/2
        i00=k-i*2
c
        j0=laddr(2,k1)-1
        nn=laddr(3,k1)+laddr(4,k1)
c
cccc        call prinf('in prolev0, laddr(k1)=*',laddr(1,k1),3)
c
        do 1200 i=1,nn*2+10
        pexp(i)=0
 1200 continue
c
        i0=(laddr(4,k1)-1)*2
        i0=i0+i00
cccc        call prinf('i0=*',i0,1)
c
        j=0
        do 1400 i=1,nn*2+2,2
        j=j+1
        pexp(i+i0)=w(j+j0)
 1400 continue
c
c        scale the coefficients properly
c
        half=1
        half=half/2
        do 1600 i=1,nn*2+4
        i0=i-1
cccc        pexp(i)=pexp(i)*dsqrt(i0+half)
c  
        pexp(i)=pexp(i)*coefs(i)
 1600 continue
c
        nn2=nn*2-2
c
c       evaluate the value of the k-th function and its derivative
c
        call legeFDER(X,VAL,der,PEXP,Nn2-1)
c
        return
        end
c
c
c
c
c
        subroutine prolvect(ier,n,c,as,bs,cs,
     1     xk,rlamouts,nvects,store,laddr,istore,eps,
     2     u,v,w,lw,nhigh)
        implicit real *8 (a-h,o-z)
        dimension as(1),bs(1),cs(1),u(1),v(1),w(1),xk(1),rlamouts(1),
     1      store(1),laddr(4,1)
c
c        this subroutine evaluates the Legendre coefficients 
c        of the first nvects solutions of the equation
c
c      D               D
c      -  ( (1-x**2) * -  \phi (x) ) + (rlam-c**2*x**2) * \phi(x) = 0,    (1)
c      Dx              Dx
c 
c        and the corresponding coefficients rlam. Note that the equation (1)
c        defines the prolate speheroidal wave functions; the coefficients 
c        rlam are the so-called separation coefficients. The coefficients
c        are stored in array store, and can be used to evaluate the 
c        prolate spheroidfal wave functions and theirr derivatives. The 
c        recommended way to do so is by calling the subroutine prolev0 (see).
c
c                          input parameters:
c
c  n - the maximum permitted number of legendre coefficients in the
c       expansion of any wave function to be computed. Must be
c       sufficiently large (to be supplied by the user!!)
c  c - the coefficient (real) in the formula (1)
c  nvects - the number of solutions of the equation (1) whose Legendre 
c       expansions are to be computed
c  eps - the accuracy to which the calculations are to be performed
c  lw - the amount of storage space provided in the array lw; must
c       be sufficiently large; if it is too small, the error code ier
c       is set to 4, and the execution of the subroutine is terminated.
c       Note that this parameter, normally, does not need to be very large.
c       
c                          output parameters:
c
c  ier - error return code. 
c         ier=0 means successful execution.
c         ier=4 means that the amount of storage space provided in
c                 the array store is insufficient (see input 
c                 parameter lw above). This is a fatal error.
c         ier=1024 means that the subroutine prolql1 (a standard 
c                 prehistoric subroutine from eispack) has failed
c                 to find the spectrum of a certain tridiagonal matrix;
c                 this has never happened yet, and probably means that memory
c                 has been mismanaged by the user. This is a fatal error.
c  rlamouts - the first nvects values of the coefficient rlam in (1) for 
c       which (1) has a non-trivial soulution.
c  store - the array containing the coefficients of the Legendre 
c          expansions of the first nvects non-sero solutions of (1).
c          Note that these solutions are normalized (the L^2 norm of
c          each is equal to 1). Normally, it is to be used by the 
c          subroutine prolev0 (see). Also, note that the first istore
c          elements of the array store should not be altered between the 
c          call to this subroutine and the subsequent calls to prolev0.
c  laddr - the map of the array store; to be used by the subroutine 
c          prolev0 (see). it will be 4*nvects integer *4 elements long.
c          EXPLANATION of the meaning of entries in the array laddr:
c      laddr(1,i) =i
c      laddr(2,i) - the location in array store of the first coefficient
c          of the Legendre expansion for the i-th eigenfunction
c      laddr(3,i) - the number of coefficients of the Legendre expansion 
c          of the i-the eigenfunction whose absolute values are greater
c          than eps; thus, the elements of the array store with numbers 
c          laddr(2,i) through laddr(2,i)+laddr(3,i)-1 contain the 
c          coefficients of the Legendre expansion of the i-th
c          eigenfunction.
c      laddr(4,1) - the number of the first coefficient in the Legendre
c          series of the i-th eigenfunction; this entry exists because 
c          the first several coefficients in the Legendre series of the 
c          higher order eigenfunctions tend to be negligibly small.
c  istore - the total number of elements in array store where something
c          has been stored
c  nhigh - the highest order of a legendre expansion actually used 
c          for any eigenfunction
c
c                             work arrays:
c
c  as,bs,cs,um,v,w - must be at least n/2+6 each
c
c        . . . construct the tridiagonal matrix whose eigenvalues 
c              are the "separation coefficients" for the prolate spheroidal
c              wave functions
c
c         . . . for the even-numbered separation coefficients
c
        ier=0
        delta=1.0d-8
        ifsymm=1
        numit=4
        rlam=0
        ifodd=-1
        call prolmatr(as,bs,cs,n,c,rlam,ifsymm,ifodd)
c   
c       find the spectrum of the tridiagonal matrix
c
        call PROLQL1(N/2,bs,as,IERR)
        if(ierr .ne. 0) ier=1024
        if(ier .ne. 0) return
c
        do 1040 i=1,n/2
        u(n/2-i+1)=-bs(i)
 1040 continue
c
c         . . . for the even-numbered separation coefficients
c
        delta=1.0d-8
        ifsymm=1
        numit=4
        rlam=0
        ifodd=1
        call prolmatr(as,bs,cs,n,c,rlam,ifsymm,ifodd)
c   
c       find the spectrum of the tridiagonal matrix
c
        call PROLQL1(N/2,bs,as,IERR)
        if(ierr .ne. 0) ier=1024
        if(ier .ne. 0) return
c
        do 1080 i=1,n/2
        v(n/2-i+1)=-bs(i)
 1080 continue
c
c
        j=0
        do 1100 i=1,n/2
        j=j+1
        rlamouts(j)=u(i)
        j=j+1
        rlamouts(j)=v(i)
c        
        if(j .eq. nvects) goto 1150
 1100 continue
 1150 continue
c
cccc        call prin2('and rlamouts=*',rlamouts,n/2)         
c
c       use the inverse power method to get the eigenvectors
c
        istore=1
        i=0
        ifodd=1
        do 4000 i7=1,nvects
cccc         call prinf('for vector number*',i7,1)
c
        i=i+1
        ifodd=-ifodd
        if(i .gt. nvects) return
        
c
        do 1200 j=1,n
        xk(j)=1
 1200 continue
c
        done=1
c
c        construct the tridiagonal matrix
c
        rlam=rlamouts(i)+delta
        ifsymm=1
        call prolmatr(as,bs,cs,n,c,rlam,ifsymm,ifodd)
c
cccc        call prin2('after prolmatr, as=*',as,n+1)
cccc        call prin2('after prolmatr, bs=*',bs,n+1)
cccc        call prin2('after prolmatr, cs=*',cs,n+1)
c
c        construct the l-u decomposition of the unsymmtrized tridiagonal
c        matrix
c
         call prolfact(bs,cs,as,n/2,u,v,w)
c
c        repeatedly apply the inverse of the unsymmtrized diagonal matrix 
c        to the seed vector, in the hope of getting 

        do 3000 ijk=1,numit
c
        call prolsolv(u,v,w,n/2,xk)
c
        d=0
        do 2400 j=1,n/2
        d=d+xk(j)**2
c
 2400 continue
c
        d=dsqrt(d)
cccccc        call prin2('d=*',d,1)
        do 2500 j=1,n/2
        xk(j)=xk(j)/d
 2500 continue
c
        err=0
        do 2600 j=1,n/2
        err=err+(as(j)-xk(j))**2
c
        as(j)=xk(j)
 2600 continue
        err=dsqrt(err)
 3000 continue
c
        lneed=istore+n/2+10
cccc          call prinf('in prolvect, lneed=*',lneed,1)
c
        if(lneed .lt. lw) goto 3400
        ier=4
          call prinf('bombing from prolvect, ier=*',ier,1)
        return
 3400 continue
c
c        store the non-zero part of the vector xk in the array store
c
        call prolpack(i,xk,store,laddr,istore,n/2,eps)
 4000 continue
c
c        find the highest order of the Legendre expansion for any 
c        eigenfunction
c
        nhigh=0
c
        do 4200 i=1,nvects
        ihigh=laddr(3,i)+laddr(4,i)
        if(ihigh .gt. nhigh) nhigh=ihigh
 4200 continue
        nhigh=nhigh*2+2
c
c        construct and put in arrays store the scaling coefficients
c        to be used by the subroutine prolev0
c
        iscale=istore+1
        lneed=iscale+nhigh+2
        if(lneed .lt. lw) goto 4400
         call prinf('bombing from prolvect, ier=*',ier,1)
        ier=4
        return
 4400 continue
c
        half=1
        half=half/2
        do 4600 i=1,nhigh
        store(iscale+i-1)=dsqrt(i-half)
 4600 continue
        istore=lneed
        laddr(1,1)=nhigh
        laddr(1,2)=iscale
c
        icoefs=istore+4
        laddr(1,3)=icoefs
        istore=istore+nhigh+8
        laddr(1,4)=nvects*4+1
c
        return
        end
c
c
c
c
c
        subroutine prolpack(k,xk,store,laddr,istore,n,eps)
        implicit real *8 (a-h,o-z)
        dimension store(1),laddr(4,1),xk(1)
c
c       find the first element of the vector of coefficients
c       that is non-zero
c
        do 1200 i=1,n
        i1=i
        if(dabs(xk(i)) .lt. eps) goto 1200
        goto 1400
 1200 continue
 1400 continue
c
c       find the last coefficient of the vector xk that is non-zero
c  
        do 1600 i=n,i1,-1
        i2=i
        if(dabs(xk(i)) .lt. eps) goto 1600
        goto 1800
 1600 continue
 1800 continue
c
c        store in array store the non-zero elements of this vector
c
        nn=i2-i1+1
        do 2000 i=1,nn
        store(istore+i-1)=xk(i1+i-1)
 2000 continue
c
c       enter the appropriate information in the array laddr
c
        laddr(1,k)=k
        laddr(2,k)=istore
        laddr(3,k)=nn
        laddr(4,k)=i1
        istore=istore+nn
c
        return
        end
c
c
c
c
c
        subroutine prolmatr(as,bs,cs,n,c,rlam,ifsymm,ifodd)
        implicit real *8 (a-h,o-z)
        dimension as(1),bs(1),cs(1)
c
c       construct the tridiagonal matrix corresponding 
c       to odd-numbered P_k
c
        if(ifodd .le. 0) goto 1300
        k=0
        do 1200 k0=1,n+2,2
c
        k=k+1
c
        call prolcoef(rlam,k0,c,alpha0,beta0,gamma0,
     1      alpha,beta,gamma)
c
        as(k)=alpha
        bs(k)=beta
        cs(k)=gamma
c 
c        remembering that the norm of P_n is not equal to 1,
c        rescale the matrix to make it symmetric

        if(ifsymm .eq. 0)  goto 1200
c
        if(k0 .gt. 1)
     1    as(k)=as(k)/dsqrt(k0-2+half)*dsqrt(k0+half)
c
        cs(k)=cs(k)*dsqrt(k0+half)/dsqrt(k0+half+2)
c 
 1200 continue
c
        return
 1300 continue
c
c       construct the tridiagonal matrix corresponding 
c       to even-numbered P_k
c
        k=0
        done=1
        half=done/2
        do 1400 k0=0,n+2,2
c
        k=k+1
c
        call prolcoef(rlam,k0,c,alpha0,beta0,gamma0,
     1      alpha,beta,gamma)
c
        as(k)=alpha
        bs(k)=beta
        cs(k)=gamma
c
c        remembering that the norm of P_n is not equal to 1,
c        rescale the matrix to make it symmetric
c
        if(ifsymm .eq. 0) goto 1400
c
        if(k0 .ne. 0) 
     1    as(k)=as(k)/dsqrt(k0-2+half)*dsqrt(k0+half)

c
        cs(k)=cs(k)*dsqrt(k0+half)/dsqrt(k0+half+2)
c 
 1400 continue
c
        return
        end
c
c
c
c
c
        subroutine prolfact(a,b,c,n,u,v,w)
        implicit real *8 (a-h,o-z)
        dimension a(1),b(1),c(1),u(1),v(1),w(1),rhs(1)
c
c        eliminate down
c
        do 1200 i=1,n-1
        d=c(i+1)/a(i)
        a(i+1)=a(i+1)-b(i)*d
        u(i)=d
 1200 continue
c
c        eliminate up
c
        do 1400 i=n,2,-1
        d=b(i-1)/a(i)
        v(i)=d
 1400 continue
c
c       scale the diagonal
c
        done=1
        do 1600 i=1,n
        w(i)=done/a(i)       
 1600 continue
c
        return
c
c
c
c
        entry prolsolv(u,v,w,n,rhs)
c
c        eliminate down
c
        do 2400 i=1,n-1
        rhs(i+1)=rhs(i+1)-u(i)*rhs(i)
 2400 continue
c
c        eliminate up
c
        do 2600 i=n,2,-1
        rhs(i-1)=rhs(i-1)-rhs(i)*v(i)
 2600 continue
c
c       scale
c
        do 2800 i=1,n
        rhs(i)=rhs(i)*w(i)
 2800 continue
        return
        end
c
c
c
c
c
        subroutine prolcoef(rlam,k,c,alpha0,beta0,gamma0,
     1      alpha,beta,gamma)
        implicit real *8 (a-h,o-z)
c
c        this subroutine evaluates the Legendre coefficients 
c        alpha0, beta0, gamma0, alpha, beta, gamma of two functions:
c
c                    
c           (1-x**2)   P_k (x) )  =
c                                                                         (1)
c           alpha0 * P_{k-2} + beta0 * P_{k} + gamma0 * P_{k+2}, 
c
c        and
c
c           D               D
c           -  ( (1-x**2) * -  P_k (x) ) + (rlam-c**2*x**2) * P_k(x) =
c           Dx              Dx
c                                                                         (2)
c           alpha * P_{k-2} + beta * P_{k} + gamma * P_{k+2}.
c
c
c                          input parameters:
c
c  rlam - the coefficient (real) in the formula (2)
c  k - the index in the formulae (1), (2)
c  c - the coefficient (real) in the formula (2)
c  
c                          output parameters:
c
c  alpha0, beta0, gamma0 - coefficients in the expansion (1)
c  alpha, beta, gamma - coefficients in the expansion (2)
c
c
        d=k*(k-1)
        d=d/(2*k+1)/(2*k-1)
        uk=d
c
        d=(k+1)**2
        d=d/(2*k+3)
        d2=k**2
        d2=d2/(2*k-1)
        vk=(d+d2)/(2*k+1)
c
        d=(k+1)*(k+2)
        d=d/(2*k+1)/(2*k+3)
        wk=d
c
        alpha=-c**2*uk
        beta=rlam-k*(k+1)-c**2*vk
        gamma=-c**2*wk
c
        alpha0=uk
        beta0=vk
        gamma0=wk
c
        return
        end
c
c
c
c
c
      SUBROUTINE PROLQL1(N,D,E,IERR)
C
      INTEGER I,J,L,M,N,II,MML,IERR
      real *8  D(N),E(N)
      real *8 B,C,F,G,P,R,S,TST1,TST2
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE IMTQL1,
C     NUM. MATH. 12, 377-383(1968) BY MARTIN AND WILKINSON,
C     AS MODIFIED IN NUM. MATH. 15, 450(1970) BY DUBRULLE.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).
c
c      Very minor changes have been introduced by V. Rokhlin, on
c      5.22.96
c
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
C     TRIDIAGONAL MATRIX BY THE IMPLICIT QL METHOD.
C
C     ON INPUT
C
C        N IS THE ORDER OF THE MATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C      ON OUTPUT
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
C          THE SMALLEST EIGENVALUES.
C
C        E HAS BEEN DESTROYED.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     CALLS PYTHAG FOR  DSQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      E(N) = 0.0D0
C
      DO 290 L = 1, N
         J = 0
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
  105    DO 110 M = L, N
            IF (M .EQ. N) GO TO 120
            TST1 = DABS(D(M)) + DABS(D(M+1))
            TST2 = TST1 + DABS(E(M))
            IF (TST2 .EQ. TST1) GO TO 120
  110    CONTINUE
C
  120    P = D(L)
         IF (M .EQ. L) GO TO 215
         IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         G = (D(L+1) - P) / (2.0D0 * E(L))
cccc         R = PYTHAG(G,1.0D0)
         r=dsqrt(g**2+1)
         G = D(M) - P + E(L) / (G + DSIGN(R,G))
         S = 1.0D0
         C = 1.0D0
         P = 0.0D0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            I = M - II
            F = S * E(I)
            B = C * E(I)
cccc            R = PYTHAG(F,G)
            r=dsqrt(f**2+g**2)
            E(I+1) = R
            IF (R .EQ. 0.0D0) GO TO 210
            S = F / R
            C = G / R
            G = D(I+1) - P
            R = (D(I) - G) * S + 2.0D0 * C * B
            P = S * R
            D(I+1) = G + P
            G = C * R - B
  200    CONTINUE
C
         D(L) = D(L) - P
         E(L) = G
         E(M) = 0.0D0
         GO TO 105
C     .......... RECOVER FROM UNDERFLOW ..........
  210    D(I+1) = D(I+1) - P
         E(M) = 0.0D0
         GO TO 105
C     .......... ORDER EIGENVALUES ..........
  215    IF (L .EQ. 1) GO TO 250
C     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
         DO 230 II = 2, L
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 270
            D(I) = D(I-1)
  230    CONTINUE
C
  250    I = 1
  270    D(I) = P
  290 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END

