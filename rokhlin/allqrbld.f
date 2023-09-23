        implicit real *8 (a-h,o-z)
        dimension par1(100),ab(20 000),
     1      x(10 000),whts(10 000),rlams(10 000),w(100 000),
     2      rints(10 000),xtest(10 000),ftest(10 000)
        external fun2
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER k '
        READ *,k
        CALL PRINF('k=*',k,1 )
c 
        PRINT *, 'ENTER nfuns'
        READ *,nfuns
        CALL PRINF('nfuns=*',nfuns,1 )
c 
        a=-1
        b=1
        eps=1.0d-10
c 
        do 1200 i=1,nfuns
c 
        par1(i)=i**2
ccc        par1(i)=par1(i)/30
 1200 continue
c 
        call prin2('par1 as created*',par1,nfuns)
c 
        lenw=2 000 000
c 
        ipar3=nfuns
c 
c       construct the SVD of the user-specified funcytions
c 
        lab=20 000
        lw=100 000
  
  
        ifab=1
c 
        call allqrbld(ier,a,b,k,fun2,par1,par2,par3,nfuns,eps,
     1    ab,lab,ifab,nn,x,whts,n,rlams,w,lw,ncols,
     2      rints,lused)
  
        call prinf('after allqrbld, ier=*',ier,1)
  
  
        do 1300 i=1,1000
        rlams(i)=0
 1300 continue
  
  
        ifab=0
c 
        call allqrbld(ier,a,b,k,fun2,par1,par2,par3,nfuns,eps,
     1    ab,lab,ifab,nn,x,whts,n,rlams,w,lw,ncols,
     2      rints,lused)
  
        call prinf('after allqrbld, ier=*',ier,1)
  
  
  
  
        if(ier .ne. 0) stop
  
        call prinf('after allqrbld, lused=*',lused,1)
        call prinf('after allqrbld, ncols=*',ncols,1)
        call prin2('after allqrbld, rlams=*',rlams,ncols)
        call prin2('after allqrbld, x=*',x,n)
        call prin2('after allqrbld, ab=*',ab,nn*2)
c 
c       calculate the singular functions at test points and plot them
c 
        ntest=1000
        h=(b-a)/ntest
        do 1400 i=1,ntest
        xtest(i)=a+h/2+(i-1)*h
 1400 continue
c 
cccc        call prin2('and xtest as constructed*',xtest,ntest)
c 
        do 2000 i=1,ncols
c 
        do 1600 j=1,ntest
c 
        call nesteva(ier,w,ab,nn,k,xtest(j),i,ftest(j),der)
 1600 continue
c 
        iw=200+i
        call lotagraph(iw,xtest,ftest,ntest,
     1      'QR-function from coefs*')
c 
 2000 continue
c 
        stop
        end
c 
c 
c 
c 
c 
        subroutine fun2(x,i,par1,par2,par3,f)
        implicit real *8 (a-h,o-z)
        save
        dimension par1(1),par2(1),par3(1)
c 
        done=1
cccc        f=exp(-par2(12-i+1)*x)+( (log((x-1)**2))**2)
cccc     1      **(done/(12-i+1))
c 
cccc        f=1
  
  
cccc        call prin2('in fun2, par1=*',par1,40)
  
cccc        f=exp(-par1(i)*x)
cccc        f=f*sqrt(x)
        f=cos(i*x)
  
cccc        f=x**(i-1)
  
  
        return
        end
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
        subroutine allqrbld(ier,a,b,k,fun,par1,par2,par3,
     1      nfuns,eps,ab,lab,ifab,nn,x,whts,n,rlams,w,lw,
     2      ncols,rints,lused)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),rlams(1),w(1),ab(2,1),par1(1),par2(1),
     1      par3(1),whts(1),rints(1)
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
c      fun(x,i,par1,par2,par3,f),
c 
c       where x is the input argument where the function is to be
c       evaluated, i the sequence number of function, f (output)
c       is the value of the i-th function at the point x, and par
c       is an input parameter, carrying whatever information fun
c       might need.
c 
c  par1,par2,par3 - whatever the user-supplied subroutine fun needs
c  nfuns - the number of functions to be compressed
c  eps - the accuracy to which the legendre expansion must represent
c      the user-supplied function for the subroutine to decide that
c      the representation is accurate
c  ab - the array of ends of subintervals created by a preceding
c      call to this subroutine.  For each i=1,2,...,nn,
c      ab(1,i) is the left end of the i-th subinterval; ab(2,i) is
c      the right end of the i-th subinterval. Note that this is an
c      output parameter only if the input parameter ifab (see above)
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
c      subroutine.  Note that this is an output parameter only
c      if the input parameter ifab (see above) had been set to 0 by
c      the user.
c  x - the discretization of the interval [a,b] consisting of nn
c      subintervals with k Legendre nodes on each,  created by a
c      preceding call to this subroutine.  Note that this is an
c      output parameter only if the input parameter ifab (see above)
c      had been set to 1 by the user.
c  whts - the weights associated with the discretization x, created
c      by a preceding call to this subroutine.   Note that this is
c      an output parameter only if the input parameter ifab (see
c      above) had been set to 1 by the user.
c  lw - the length of the user-provided array w. It has to be at
c      least 17/8*n*nfuns+n+2*k**2+1000. It it is insufficient,
c      ier os set to 64, and the execution of the subroutine
c      terminated
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
c                        not be changed between tie call to this
c                        subroutine and subsequent calls to the
c                        subroutine nesteva.
c 
c  x - the discretization of the interval [a,b] consisting of nn
c        subintervals with k Legendre nodes on each.  Note that this is an
c     output parameter only if the input parameter ifab (see above)
c     had been set to 1 by the user.
c  whts - the weights associated with the discretization x.  Note that
c     this is an output parameter only if the input parameter ifab
c     (see above) had been set to 1 by the user.
c  n - the number of nodes in the discretization x of the interval [a,b];
c        equal to nn*k
c  rlams - the singular values of the matrix of functions given by the
c        function fun. Note that though only ncols elements of array rlams
c        are meaningful, nfuns elements have to be allocated by the
c        user.
c  ncols - the number of QR-functions in the QR-decomposition of the
c      nfuns functions gven by the subroutine fun (see above) whose
c      corresponding weights are greater than eps (see above). It is
c      also the number of meaningful columns of arrays fs, coefs on
c      exit from this subroutine.
c 
c                  Work arrays:
c 
c 
c 
c        . . . construct the nested subdivision of the interval [a,b],
c              such that on each subinterval, the pair-wise products
c              of all user-provided functions are integrated accurately
c              with k nodes
c 
        if(ifab .eq. 0) goto 1200
  
         call manincon(ier,a,b,k,fun,nfuns,par1,par2,par3,
     1    eps,ab,lab,nn,w,lw)
c 
         if(ier .ne. 0) return
c 
 1200 continue
c 
c       . . . construct the all the nested Legendre nodes
c 
        ier=0
        call nodewhts(ab,nn,k,x,whts)
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
        icoefs=1
        lcoefs=n*nfuns+10
c 
        iu=icoefs+lcoefs
        lu=k**2+2
c 
        iv=iu+lu
        lv=k**2+2
c 
        iww=iv+lv
        lww=n+2
c 
        iwwww=iww+lww
        lwwww=n*nfuns*9/8+n+20
c 
        lused=iwwww+lwwww
        if(lused .lt. lw) goto 1400
        ier=64
        return
 1400 continue
c 
        call allffuns(x,whts,n,nfuns,fun,par1,par2,par3,
     1      eps,w(icoefs),ncols,rlams,nn,k,w(iu),w(iv),
     2      w(iww),rints,w(iwwww) )
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine nesteva(ier,coefs,ab,nn,k,x,i,val,der)
        implicit real *8 (a-h,o-z)
        save
        dimension coefs(nn*k,1),ab(2,1)
c 
c        this subroutine evaluates at the point x \in [a,b]
c        one of singular functions specified by their nested Legendre
c        expansions, whose coefficients are stored in the input array
c        coefs. Note that this subroutine uses output of the subroutine
c        allqrbld (see) as its input; it has no function as a
c        stand-alone device.
c 
c                   Input parameters:
c 
c  coefs - the nested Legendre expansions of the functions stored in the
c      array fs.
c  ab - the array of ends of subintervals created. For each i=1,2,...,nn,
c     ab(1,i) is the left end of the i-th subinterval; ab(2,i) is
c     the right end of the i-th subinterval
c  nn - the number of subintervals created on the interval [a,b];
c     also the number of entries in the array ab (see above)
c  k - the number of Gaussian nodes on each of nn subintervals in the
c      nested subdivision of the interval [a,b].
c  x - the point where the QR-function is to be evaluated; must be on
c      the interval [a,b]
c  i - the sequence number of the QR-function to be evaluated; must be
c      between 1 and ncols.
c 
c                    Output parameters:
c 
c  val - the value of the i-th QR-function at the point x
c  der - the derivative of the i-th QR-function at the point x
c 
c       . . . find the subinterval in which the point x lives
c 
        ier=0
        do 1200 j=1,nn
c 
        intnum=j
c 
        if(ab(2,j) .gt. x) goto 1400
 1200 continue
c 
        ier=4
        return
 1400 continue
c 
c       the point x lives in the interval number intnum.
c       evaluate the expansion at this point
c 
        u=2/(ab(2,intnum)-ab(1,intnum))
        v=1-ab(2,intnum)*u
c 
        t=u*x+v
c 
        jj=(intnum-1)*k+1
c 
        call legeFDER(t,VAL,der,coefs(jj,i),k-1)
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
        subroutine allffuns(x,whts,n,nfuns,fun,par1,par2,par3,
     1      eps,fs,ncols,sss,nn,k,u,v,ww,rints,www)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),fs(n,1),whts(1),ww(1),par1(1),par2(1),
     1      par3(1),u(1),v(1),xlege(2000),whtslege(2000),
     2      rints(1),www(1),sss(1)
c 
c        this subroutine constructs a QR-decomposition of a set
c        of user-specified functions on the interval [a,b], and
c        nested Legendre expansions for the first ncols functions
c        in the said QR-decompositions; all elements with higher
c        number than ncols have weights associated with them that
c        are less than eps (user-specified).
c 
c            Input parameters:
c 
c  x - the discretization of the interval [a,b] (normally produced
c        by the subroutine manincon (see)
c  whts - the weights associated with the discretization x; (normally
c        produced by the subroutine manincon (see)
c  n - the number of nodes in the discretization x of the interval [a,b]
c  nfuns - the number of functions to be compressed
c  fun - the subroutine evaluating the functions to be compressed at
c        arbitrary points on the interval [a,b]. The calling sequence
c        of fun must be:
c 
c      fun(x,i,par1,par2,par3,f),
c 
c      where x is the input argument where the function is to be
c      evaluated, i the sequence number of function, f (output)
c      is the value of the i-th function at the point x, and par
c      is an input parameter, carrying whatever information fun
c      might need.
c 
c  par1,par2,par3 - whatever the user-supplied subroutine fun needs
c  eps - the accuracy to which the legendre expansion must represent
c      the user-supplied function for the subroutine to decide that
c      the representation is accurate
c  nn - the number of subintervals in the nested discretization of
c      the interval [a,b] (normally produced by the subroutine
c      manincon (see)
c  k - the number of Gaussian nodes on each of nn subintervals in the
c      nested subdivision of the interval [a,b]; note that n=nn*k
c      - it better be!!!
c 
c           Output parameters:
c 
c  fs - the first ncols of the nfuns QR vectors. Note that the array
c      fs has to be dimensioned n*nfuns, and all n*nfuns elements are
c      used by the subroutine. However, on exit only ncols columns of
c      fs (each n elements long) are meaningful.
c  ncols - the number of QR-functions in the QR-decomposition of the
c      nfuns functions gven by the subroutine fun (see above) whose
c      corresponding weights are greater than eps (see above). It is
c      also the number of meaningful columns of arrays fs, coefs on
c      exit from this subroutine.
c  rnorms - the weights corresponding to the columns of arrays fs,
c      coefs. Note that though only ncols elements of array rnorms
c      are meaningful, nfuns elements have to be allocated by the
c      user.
c  coefs - the nested Legendre expansions of the functions stored in the
c      array fs.
c 
c              Work arrays:
c 
c  u,v - must be at least k**2+2 real *8 elements long each.
c  ww - must be at least n real *8 elements each
c  www - must be at least (n*nfuns)*9/8 +n + 10 real *8 elements each
c 
c        . . . constuct the square roots of weights
c 
        do 1100 i=1,n
        ww(i)=sqrt(whts(i))
 1100 continue
c 
c        construct the values of all functions at all
c        nodes of the discretization
c 
        do 1400 i=1,nfuns
        do 1200 j=1,n
c 
        call fun(x(j),i,par1,par2,par3,d)
        fs(j,i)=d*ww(j)
 1200 continue
 1400 continue
c 
c       conduct the Gram-Schmidt procedure on the constructed functions
c 
        call qrdpivot(fs,n,nfuns,ncols,eps,www,sss)
c 
        call prin2('in alffuns after qrdpivot, sss=*',sss,ncols)
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
c        plot the obtained QR-vectors
c 
        do 2000 i=1,ncols
c 
        iw=100+i
c 
        call lotagraph(iw,x,fs(1,i),n,'QR-function*')
 2000 continue
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
        subroutine manincon(ier,a,b,k,fun,nfuns,par1,par2,par3,
     1    eps,ab,lab,nn,w,lenw)
        save
        dimension par1(1),ab(2,1),w(1),par2(1),par3(1)
        implicit real *8 (a-h,o-z)
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
c      fun(x,i,par1,par2,par3,f),
c 
c      where x is the input argument where the function is to be
c      evaluated, i the sequence number of function, f (output)
c      is the value of the i-th function at the point x, and par
c      is an input parameter, carrying whatever information fun
c      might need.
c 
c  par1,par2,par3 - whatever the user-supplied subroutine fun needs
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
         call allincon(ier,a,b,k,fun,ipar1,par1,par2,par3,eps,
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
         call allincon(ier,a,b,k,fun,ipar1,par1,par2,par3,eps,
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
     1    par2,par3,eps,ab,lab,nn,labused,w,ifinit)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),par1(1),par2(1),ab(2,1),par3(1)
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
c      fun(x,i,par1,par2,par3,f),
c 
c      where x is the input argument where the function is to be
c      evaluated, i is the sequence  number of the function to be
c      evaluated, f  (output) is the function, and par1,par2,apr3
c      are input parameters, carrying whatever information fun
c      might need.
c 
c  par1,par2,par3 - whatever the user-supplied subroutine fun needs
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
         call allinter(ier,a,b,fun,ipar1,par1,par2,par3,k,ab,lab,
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
         subroutine allinter(ier,a,b,fun,ipar1,par1,par2,par3,
     1    k,ab,lab,eps,t,u,v,whts,x,f,coefs,nn,lused,ifinit7)
        implicit real *8 (a-h,o-z)
        save
        dimension t(1),u(2*k,2*k),v(2*k,2*k),whts(1),x(1),
     1      f(1),coefs(1),par1(1),par2(1),ab(2,1),par3(1)
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
     1    par2,par3,eps/d,t,u,v,whts,x,ifinit,f,coefs)
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
     1    par2,par3,eps,t,u,v,whts,x,ifinit,f,coefs)
        implicit real *8 (a-h,o-z)
        save
        dimension t(1),u(2*k,2*k),v(2*k,2*k),whts(1),x(1),
     1      f(1),coefs(1),par1(1),par2(1),par3(1)
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
c      fun(x,i,par1,par2,par3,f),
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
        call fun(x(i),ipar1,par1,par2,par3,f(i) )
 1400 continue
cccc        call prin2('in ifsplit, x=*',x,k*2)
cccc        call prin2('in ifsplit, f=*',f,k*2)
c 
c       construct the legendre expansion of the function on
c       the interval
c 
        call allqrmat(u,f,coefs,k*2)
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
