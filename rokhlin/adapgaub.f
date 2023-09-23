        implicit real *8 (a-h,o-z)
        external fun
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
c 
        PRINT *, 'ENTER m '
        READ *,m
        CALL PRINF('m=*',m,1 )
c 
c        integrate
c 
        a=-1
        b=1
        eps=1.0d-13
  
  
        ifinit=1
cccc        call funask(ncalled,ifinit)
        ddd=funask(ncalled,ifinit)
c 
        ifends=1
        call adapgaub(ier,a,b,fun,jtest,par2,m,eps,
     1      ifends,rint2,maxrec2,numint2)
  
        call prin2('after adapgaub, rint2=*',rint2,1)
        call prinf('after adapgaub, ier=*',ier,1)
        call prinf('after adapgaub, maxrec2=*',maxrec2,1)
        call prinf('after adapgaub, numint2=*',numint2,1)
  
  
        ifinit=1
cccc        call funask(ncalled2,ifinit)
        ddd=funask(ncalled2,ifinit)
  
        ifends=0
        call adapgaub(ier,a,b,fun,jtest,par2,m,eps,
     1      ifends,rint,maxrec,numint)
  
  
cccc        call funask(ncalled1,ifinit)
        ddd=funask(ncalled1,ifinit)
  
  
        call prin2('after adapgaus, rint=*',rint,1)
        call prinf('after adapgaus, ier=*',ier,1)
        call prinf('after adapgaus, maxrec=*',maxrec,1)
        call prinf('after adapgaus, numint=*',numint,1)
  
  
        call prin2('and rint2-rint=*',rint2-rint,1)
  
  
        nfuns=numint*m
  
        call prinf('and nfuns=*',nfuns,1)
        call prinf('and number of function calls for adapgaus=*',
     1      ncalled1,1)
  
        call prinf('and number of function calls for adapgaub=*',
     1      ncalled2,1)
  
c 
        stop
        end
c 
c 
c 
c 
c 
c 
        function fun(x,jtest,dummy)
        implicit real *8 (a-h,o-z)
c 
c 
ccccc        fun=log(x+1+0.00000000001d0)
        save
        fun=log(x+1+0.0001d0)
cccc        fun=log(x+1)
        fun=log(1-x**2+1.0d-12)
  
cccc        fun=1/(1+x)**0.7
  
cccc        fun=cos(300*x)
  
cccc        fun=cos(300*x)+sin(400*x) * fun
  
cccc        fun=1
  
        icount=icount+1
        return
c 
c 
c 
c 
        entry funask(ncalled,ifinit)
  
        ncalled=icount
        if(ifinit .eq. 1) icount=0
        funask=0
        return
        end
c 
c 
c 
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        this is the end of the debugging code and the
c        start of the actual quadrature routines.
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c 
c 
        subroutine adapgaub(ier,a,b,fun,par1,par2,m,eps,
     1      ifends,rint,maxrec,numint)
        implicit real *8 (a-h,o-z)
        save
        dimension t(100),w(100),stack(400),vals(200),
     1      par1(1),par2(1),xs(100),ws(100)
c 
c 
        data m7/-2341034/
c 
c       this subroutine uses the adaptive chebychev quadrature
c       to evaluate the integral of the user-supplied function
c       fun on the user-specified interval [a,b]
c 
c                       input parameters:
c 
c  a,b - the ends of the interval on which the integral is to
c       be evaluated
c  fun - the user-supplied function to be integrated. the calling
c       sequence of fun must be
c 
c        fun(x,par1,par2).                            (1)
c 
c        in (1), x is a point on the interval [a,b] where
c        the function is to be evaluated, and par1, par2
c        are two parameters to be used by fun; they can be
c        variables or arrays, real or integer, as desired.
c        fun is assumed to be real *8.
c  par1, par2 - partameters to be used by the user-supplied
c       function fun (see above)
c  m - the order of the quadrature to me used on each subinterval
c  eps - the accuracy (absolute) to which the integral will be
c       evaluated
c  ifends - an integer parameter telling the subroutine whether it
c       should apply special-purpose quadratures at the ends of the
c       interval, or use adaptive Gaussian quadrature throughout.
c   ifends=1 will cause the subroutine to treat the subintervals at
c       the ends differently
c   ifends=0 will cause the subroutine to use adaptive Gaussian
c       quadratures throughout
c 
c                       output parameters:
c 
c  ier - error return code.
c          ier=0 means normal conclusion
c          ier=8 means that at some point, one subinterval in the
c                subdivision was smaller than (b-a)/2**200. this
c                is a fatal error.
c          ier=16 means that the total number of subintervals in the
c                adaptive subdivision of [a,b] turned out to be greater
c                than 100000.  this is a fatal error.
c 
c  rint - the integral as evaluated
c  maxrec - the maximum depth to which the recursion went at its
c         deepest point. can not be greater than 200, since at that
c         point ier is set to 8 and the execution of the subroutine
c         terminated.
c  numint - the totla number of intervals in the subdivision. can not
c         be greater than 100000,  since at that
c         point ier is set to 16 and the execution of the subroutine
c         terminated.
c 
c 
c          . . . check if the gauss quadarture has been
c                initialized at a preceeding call to this routine
c 
        if(m .eq. m7) goto 1200
        call gauswhts2(m,t,w)
c 
        m7=m
c 
        n=40
        ndigits=16
  
        call adapuniv(ier,ndigits,n,xs,ws)
  
cccc        call prin2('after adapuniv, xs=*',xs,n)
cccc        call prin2('after adapuniv, ws=*',ws,n)
 1200 continue
c 
c        integrate the user-supplied function using the
c        adaptive gaussian quadratures
c 
        nnmax=100000
        maxdepth=200
c 
        call adinrec2(ier,stack,a,b,fun,par1,par2,t,w,m,
     1      vals,nnmax,eps,rint,maxdepth,maxrec,numint,
     2      xs,ws,n,ifends)
        return
        end
c 
c 
c 
c 
c 
        subroutine adinrec2(ier,stack,a,b,fun,
     1      par1,par2,t,w,m,vals,nnmax,eps,
     2      rint,maxdepth,maxrec,numint,
     3      xs,ws,n,ifends)
        implicit real *8 (a-h,o-z)
        save
        dimension stack(2,1),t(1),w(1),vals(1),par1(1),par2(1)
c 
c       start the recursion
c 
        stack(1,1)=a
        stack(2,1)=b
        call oneint2(a,b,fun,par1,par2,t,w,m,vals(1))
c 
c       recursively integrate the thing
c 
        j=1
        rint=0
        ier=0
        maxrec=0
        do 3000 i=1,nnmax
        numint=i
        if(j .gt. maxrec) maxrec=j
c 
c       subdivide the current subinterval
c 
         c=(stack(1,j)+stack(2,j))/2
c 
        call oneint3(stack(1,j),c,fun,
     1      par1,par2,t,w,m,value2,
     2      a,b,xs,ws,n,ifends)
c 
        call oneint3(c,stack(2,j),fun,
     1      par1,par2,t,w,m,value3,
     2      a,b,xs,ws,n,ifends)
c 
        dd=dabs(value2+value3-vals(j))
        ifdone=0
        if(dd .le. eps) ifdone=1
c 
c       if the function on this subinterval has been
c       integrated with sufficient accuracy - add the
c       value to that of the global integral and move up
c       in the stack
c 
        if(ifdone  .eq. 0) goto 2000
c 
        rint=rint+value2+value3
        j=j-1
c 
c        if the whole thing has been integrated - return
c 
        if(j .eq. 0) return
        goto 3000
 2000 continue
c 
c       if the function on this subinterval has not been
c       integrated with sufficient accuracy - move
c       down the stack
c 
        stack(1,j+1)=stack(1,j)
        stack(2,j+1)=(stack(1,j)+stack(2,j))/2
        vals(j+1)=value2
c 
        stack(1,j)=(stack(1,j)+stack(2,j))/2
        vals(j)=value3
c 
        j=j+1
c 
c       if the depth of the recursion has become excessive - bomb
c 
        if(j .le. maxdepth) goto 3000
        ier=8
        return
 3000 continue
        ier=16
        return
        end
c 
c 
c 
c 
c 
        subroutine oneint3(a,b,fun,par1,par2,t,w,m,rint,
     1      aa,bb,xs,ws,n,ifends)
        implicit real *8 (a-h,o-z)
        save
        dimension t(1),w(1),par1(1),par2(1),xs(1),ws(1)
c 
c       if the subinterval being prosessed does not include
c       either the left or the right end of the big interval
c       - use Gaussian quadrature
c 
        if(ifends .eq. 1) goto 1400
        call oneint2(a,b,fun,par1,par2,t,w,m,rint)
        return
c 
 1400 continue
  
        eps=1.0d-14
        ifleft=0
        ifright=0
c 
        d=abs(aa-a)/(bb-aa)
        if(d .lt. eps) ifleft=1
c 
        d=abs(bb-b)/(bb-aa)
        if(d .lt. eps) ifright=1
c 
        if( (ifleft .ne. 0) .or. (ifright .ne. 0) ) goto 2000
c 
        call oneint2(a,b,fun,par1,par2,t,w,m,rint)
c 
        return
c 
 2000 continue
c 
c        if the subinterval contains the left end of the big
c        interval - act accordingly
c 
        if(ifleft .eq. 1)
     1      call oneleft(a,b,fun,par1,par2,xs,ws,n,rint)
c 
        if(ifright .eq. 1)
     1      call oneright(a,b,fun,par1,par2,xs,ws,n,rint)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine oneright(a,b,fun,par1,par2,xs,ws,n,rint)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),ws(1),par1(1),par2(1)
c 
c        map the interval [-1,0] on the interval [a,b]
c 
        alpha=a-b
        beta=a
        rint=0
        do 1200 i=1,n
c 
        ddd=alpha*xs(i)+beta
        rint=rint+fun(ddd,par1,par2)*ws(i)
 1200 continue
c 
        rint=-rint*alpha
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine oneleft(a,b,fun,par1,par2,xs,ws,n,rint)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),ws(1),par1(1),par2(1),xxs(100)
c 
c        map the interval [-1,0] on the interval [a,b]
c 
        alpha=b-a
        beta=b
        rint=0
        do 1200 i=1,n
c 
        xxs(i)=alpha*xs(i)+beta
c 
        rint=rint+fun(xxs(i),par1,par2)*ws(i)
 1200 continue
c 
        rint=rint*alpha
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine oneint2(a,b,fun,par1,par2,t,w,m,rint)
        implicit real *8 (a-h,o-z)
        save
        dimension t(1),w(1),par1(1),par2(1)
c 
c       integrate the function fun on the interval [a,b]
c 
        rint=0
        u=(b-a)/2
        v=(b+a)/2
        do 1200 i=1,m
        tt=u*t(i)+v
        rint=rint+fun(tt,par1,par2)*w(i)
 1200 continue
        rint=rint*u
        return
        end
c 
c 
c 
c 
c 
        subroutine gauswhts2(n,ts,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension ts(1),whts(1)
c 
c        this subroutine constructs the nodes and the
c        weights of the n-point gaussian quadrature on
c        the interval [-1,1]
c 
c                input parameters:
c 
c  n - the number of nodes in the quadrature
c 
c                output parameters:
c 
c  ts - the nodes of the n-point gaussian quadrature
c  w - the weights of the n-point gaussian quadrature
c 
c       . . . construct the array of initial approximations
c             to the roots of the n-th legendre polynomial
c 
        eps=1.0d-14
        ZERO=0
        DONE=1
        pi=datan(done)*4
        h=pi/(2*n)
        do 1200 i=1,n
        t=(2*i-1)*h
        ts(n-i+1)=dcos(t)
1200  CONTINUE
c 
c         use newton to find all roots of the legendre polynomial
c 
        ts(n/2+1)=0
        do 2000 i=1,n/2
c 
        xk=ts(i)
        deltold=1
        do 1400 k=1,10
        call legpol2(xk,n,pol,der)
        delta=-pol/der
        xk=xk+delta
        if(deltold .le. eps) goto 1600
        deltold=dabs(delta)
 1400 continue
 1600 continue
        ts(i)=xk
        ts(n-i+1)=-xk
 2000 continue
c 
c       now, use the explicit integral formulae
c       to obtain the weights
c 
        a=-1
        b=1
        do 2200 i=1,n/2+1
        call prodeng2(a,ts,n,i,fm)
        call prodeng2(b,ts,n,i,fp)
        whts(i)=fp-fm
        whts(n-i+1)=whts(i)
 2200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine legpol2(x,n,pol,der)
        implicit real *8 (a-h,o-z)
c 
        save
        pkm1=1
        pk=x
c 
        pk=1
        pkp1=x
c 
c        if n=0 or n=1 - exit
c 
        if(n .ge. 2) goto 1200
        pol=1
        der=0
        if(n .eq. 0) return
c 
        pol=x
        der=1
        return
 1200 continue
c 
c       n is greater than 1. conduct recursion
c 
        do 2000 k=1,n-1
        pkm1=pk
        pk=pkp1
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
 2000 continue
c 
c        calculate the derivative
c 
        pol=pkp1
        der=n*(x*pkp1-pk)/(x**2-1)
        return
        end
c 
c 
c 
c 
c 
        subroutine prodeng2(x,xs,n,i,f)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1)
c 
c      evaluate the product
c 
        f=1
        dlarge=1.0d20
        dsmall=f/dlarge
c 
        large=0
        do 2000 j=1,n
        dd=dabs(f)
        if( dd .gt. dsmall) goto 1200
         f=f*10000
         large=large-1
 1200 continue
c 
        if( dd .lt. dlarge) goto 1400
        f=f/10000
        large=large+1
 1400 continue
        if(j .eq. i) goto 2000
        f=f*(x-xs(j))/(xs(i)-xs(j))
 2000 continue
        d10000=10000
        f=f*d10000**large
        f=f**2*(x-xs(i))
        return
        end
c 
c 
c 
c 
c 
        subroutine adapuniv(ier,ndigits,n,xs,ws)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),ws(1),xs16by40(40),ws16by40(40)
c 
c        This subroutine returns to the user the nodes and weights of
c        a one-sided "universal" quadrature on the interval [-1,0],
c        with -1 the "universal" end. The quadrature has 40 nodes, and
c        integrates things to 16 digits (in exact arithmetic)
c 
c        Nodes and weights for ndigits= 16, n= 40
c 
         data xs16by40/
     1   -.999999999999999638213D+00,-.999999999999850990523D+00,
     2   -.999999999992186550978D+00,-.999999999827124757114D+00,
     3   -.999999997727871019351D+00,-.999999979263791480363D+00,
     4   -.999999856621124858991D+00,-.999999206545693172313D+00,
     5   -.999996351787170751847D+00,-.999985681269457789468D+00,
     6   -.999951034877759847348D+00,-.999851743339178315595D+00,
     7   -.999597329306835988124D+00,-.999008237700004006109D+00,
     8   -.997764371468272913548D+00,-.995350461097497918423D+00,
     9   -.991015085374625768967D+00,-.983764483687498345634D+00,
     *   -.972409335006481187259D+00,-.955671238715704403634D+00,
     1   -.932338571443200562382D+00,-.901444800414161659578D+00,
     2   -.862432764581415088179D+00,-.815269764207316844432D+00,
     3   -.760489988127791768482D+00,-.699158430372988125010D+00,
     4   -.632767944051531988853D+00,-.563093354764438499422D+00,
     5   -.492031154779559888297D+00,-.421450653742502302812D+00,
     6   -.353074986777930032128D+00,-.288401137975065100288D+00,
     7   -.228659712288983866694D+00,-.174809114334128649934D+00,
     8   -.127555508163444178377D+00,-.873890822901681860001D-01,
     9   -.546279769530087384795D-01,-.294629258638504824256D-01,
     *   -.119974454227582840418D-01,-.227781718570115190455D-02/
c 
         data ws16by40/
     1   0.306323159953996293057D-14,0.688802148757552364188D-12,
     2   0.269033149926221370397D-10,0.484265815622118312601D-09,
     3   0.539651338414544324975D-08,0.427424888191308276147D-07,
     4   0.260216687163476989726D-06,0.127987751960248555570D-05,
     5   0.526289764203711167827D-05,0.185491092719621790396D-04,
     6   0.571078103728377118989D-04,0.155897400539893271980D-03,
     7   0.381977637458643632856D-03,0.848603156477204212050D-03,
     8   0.172427114008955818331D-02,0.322863637784786086469D-02,
     9   0.560850295466324093409D-02,0.909272000271185977379D-02,
     *   0.138333088403933531919D-01,0.198478844161633891520D-01,
     1   0.269813650737590886660D-01,0.349010136783211678237D-01,
     2   0.431292288541331445114D-01,0.511071090134078287838D-01,
     3   0.582731681433076104215D-01,0.641386361403872542852D-01,
     4   0.683437969038799700087D-01,0.706867743339613677743D-01,
     5   0.711240694337088694613D-01,0.697484097820044878554D-01,
     6   0.667527702143166846811D-01,0.623897510768036949789D-01,
     7   0.569337220873735894009D-01,0.506504051417116620242D-01,
     8   0.437758578245055005959D-01,0.365047265616776777117D-01,
     9   0.289863743573970395161D-01,0.213270101039698862548D-01,
     *   0.135963366403290116612D-01,0.584516414634265575243D-02/
c 
        do 1200 i=1,n
c 
        ws(i)=ws16by40(i)
        xs(i)=xs16by40(i)
 1200 continue
c 
        return
        end
