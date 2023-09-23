        implicit real *8 (a-h,o-z)
        complex *16 rint(100),w(10 000),
     1      traps(100),cds(100)
        external funsub
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER m '
        READ *,m
        CALL PRINF('m=*',m,1 )
C 
        PRINT *, 'ENTER n'
        READ *,n
        CALL PRINF('n=*',n,1 )
c 
c        integrate
c 
        jtest=4
c 
        a=0.01
        b=1
        eps=1.0d-10
        lenw=10 000
  
  
c 
        call cadavect(ier,a,b,funsub,n,jtest,par2,m,eps,
     1      rint,maxrec,numint,vals,lenw)
        call prin2('after cadavect, rint=*',rint,2*n)
        call prinf('after cadavect, ier=*',ier,1)
        call prinf('after cadavect, maxrec=*',maxrec,1)
        call prinf('after cadavect, numint=*',numint,1)
 1100 format(' and more accurately, rint= ',e22.16,2x,e22.16)
cccc        write(6,2200) rint
cccc        write(13,2200) rint
cccc          if(2 .ne. 3) stop
c 
c        evaluate the integral via trapezoidal rule
c 
        do 1200 i=1,n
        traps(i)=0
 1200 continue
c 
        nn=100000
        h=(b-a)/(nn+1)
        do 2000 i=1,nn
        t=a+i*h
  
        call funsub(t,n,cds,jtest,dummy)
c 
        do 1400 j=1,n
        traps(j)=traps(j)+cds(j)
 1400 continue
c 
 2000 continue
c 
        do 2200 i=1,n
        traps(i)=traps(i)*h
 2200 continue
c 
        call prin2('and traps=*',traps,2*n)
  
  
  
c 
        stop
        end
c 
c 
c 
c 
c 
        subroutine funsub(x,n,fout,jtest,dummy)
        implicit real *8 (a-h,o-z)
        save
        complex *16 fout(n),ima
        data ima/(0.0d0,1.0d0)/
c 
        do 1200 i=1,n
        fout(i)=x**(i+jtest)+ima*(dlog(x))**i
 1200 continue
        return
        end
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
        subroutine cadavect(ier,a,b,fun,n,par1,par2,m,eps,
cccc     1      rint,maxrec,numint,vals)
     1      rint,maxrec,numint,ww,lenww)
        implicit real *8 (a-h,o-z)
        save
        complex *16 rint(1),ww(1)
        dimension t(100),w(100),stack(400),
     1      par1(1),par2(1)
        data m7/-2341034/
c 
c       this subroutine uses the adaptive Gaussian quadrature
c       to evaluate the integral of the user-supplied vector-valued
c       function fun: [a,b] \to R^n on the user-specified interval [a,b]
c 
c                       input parameters:
c 
c  a,b - the ends of the interval on which the integral is to
c       be evaluated
c  fun - the user-supplied SUBROUTINE evaluating the function to
c       be integrated. the calling sequence of fun must be
c 
c        callfun(x,n,fout,par1,par2).                            (1)
c 
c        in (1), x is a point on the interval [a,b] where
c        the function is to be evaluated, and par1, par2
c        are two parameters to be used by fun; they can be
c        variables or arrays, real or integer, as desired.
c        The output vector fout is assumed to be real *8 vector
c        of length n.
c  n - the dimension of the vector-values function fun to be
c        integrated
c  par1, par2 - parameters to be used by the user-supplied
c        subroutine fun (see above)
c  m - the order of the quadrature to me used on each subinterval
c  eps - the accuracy (absolute) to which the integral will be
c       evaluated
c  lenww - the length of the user-supplied work array ww. While
c        n*400 is alway sufficient, the actual requirements are
c        normally much less. If the array is too short, it will
c        limit the depth of recursion to which the subroutine
c        can go. When it runs out of memory in ww, it sets ier
c        to 4 and bombs.
c 
c                       output parameters:
c 
c  ier - error return code.
c          ier=0 means normal conclusion
c          ier=4 means that the user-supplied work array ww is too
c                short. this is a fatal error.
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
c                        work arrays:
c 
c  ww - must be sufficiently long (see parameter lenww above).
c 
c          . . . check if the chebychev quadarture has been
c                initialized at a preceeding call to this routine
c 
        if(m .eq. m7) goto 1200
        call gauswhts(m,t,w)
  
        m7=m
 1200 continue
c 
c        integrate the user-supplied function using the
c        adaptive gaussian quadratures
c 
        nnmax=100000
        maxdepth=200
c 
        ivalue2=1
        lvalue2=n*2+4
c 
        ivalue3=ivalue2+lvalue2
        lvalue3=n*2+4
c 
        ifout=ivalue3+lvalue3
        lfout=n*2+4
c 
        ival=ifout+lfout
c 
        lleft=lenww-ival
c 
        mmm=(lleft-100)/(n*2)
cccc        call prinf('mmm as calculated*',mmm,1)
        if(maxdepth .lt. mmm) maxdepth=mmm
c 
        call cadvecre(ier,stack,a,b,fun,n,par1,par2,t,w,m,
cccc     1      ww,nnmax,eps,rint,maxdepth,maxrec,numint)
     1      ww(ival),nnmax,eps,rint,maxdepth,maxrec,numint,
     2      ww(ivalue2),ww(ivalue3),ww(ifout) )
        return
        end
c 
c 
c 
c 
c 
        subroutine cadvecre(ier,stack,a,b,fun,n,
     1      par1,par2,t,w,m,vals,nnmax,eps,
     2      rint,maxdepth,maxrec,numint,value2,value3,fout)
        implicit real *8 (a-h,o-z)
        save
        complex *16 rint(1),value2(1000),value3(1000),vals(n,1),
     1      fout(1)
        dimension stack(2,1),t(1),w(1),par1(1),par2(1)
c 
c       start the recursion
c 
        stack(1,1)=a
        stack(2,1)=b
        call convecin(a,b,fun,n,par1,par2,t,w,m,vals(1,1),fout)
c 
c       recursively integrate the thing
c 
        j=1
        do 1200 jjj=1,n
c 
        rint(jjj)=0
 1200 continue
c 
        ier=0
        maxrec=0
        do 3000 i=1,nnmax
ccc        call prinf('i=*',i,1)
        numint=i
        if(j .gt. maxrec) maxrec=j
ccc        call prinf('j=*',j,1)
c 
c       subdivide the current subinterval
c 
         c=(stack(1,j)+stack(2,j))/2
        call convecin(stack(1,j),c,fun,n,
     1      par1,par2,t,w,m,value2,fout)
c 
        call convecin(c,stack(2,j),fun,n,
     1      par1,par2,t,w,m,value3,fout)
c 
        dd=0
        do 1400 jjj=1,n
c 
        ddd=cdabs(value2(jjj)+value3(jjj)-vals(jjj,j))
        if(dd .lt. ddd) dd=ddd
 1400 continue
c 
cccc        dd=cdabs(value2+value3-vals(j))
ccccc         call prin2('in cadvecre, dd=*',dd,1)
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
        do 1600 jjj=1,n
c 
        rint(jjj)=rint(jjj)+value2(jjj)+value3(jjj)
 1600 continue
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
c 
        do 2200 jjj=1,n
        vals(jjj,j+1)=value2(jjj)
 2200 continue
c 
        stack(1,j)=(stack(1,j)+stack(2,j))/2
c 
        do 2400 jjj=1,n
        vals(jjj,j)=value3(jjj)
 2400 continue
c 
cccc        vals(j)=value3
c 
        j=j+1
c 
c       if the depth of the recursion has become excessive - bomb
c 
        if(j .le. maxdepth) goto 3000
c 
        ier=8
        if(maxdepth .lt. 200) ier=4
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
        subroutine convecin(a,b,fun,n,par1,par2,t,w,m,
     1      rint,fout)
        implicit real *8 (a-h,o-z)
        save
        dimension t(1),w(1),par1(1),par2(1)
        complex *16 rint(1),fout(1000)
c 
c       integrate the function fun on the interval [a,b]
c 
        do 1100 i=1,n
c 
        rint(i)=0
 1100 continue
c 
        u=(b-a)/2
        v=(b+a)/2
        do 1200 i=1,m
        tt=u*t(i)+v
c 
        call fun(tt,n,fout,par1,par2)
c 
        do 1150 j=1,n
c 
        rint(j)=rint(j)+fout(j)*w(i)
 1150 continue
c 
ccc       rint=rint+fun(tt,par1,par2)*w(i)
 1200 continue
c 
        do 1400 i=1,n
c 
        rint(i)=rint(i)*u
 1400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine gauswhts(n,ts,whts)
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
        call legpol(xk,n,pol,der)
        delta=-pol/der
cccc         call prin2('delta=*',delta,1)
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
        call prodend(a,ts,n,i,fm)
        call prodend(b,ts,n,i,fp)
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
        subroutine legpol(x,n,pol,der)
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
        subroutine prodend(x,xs,n,i,f)
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
