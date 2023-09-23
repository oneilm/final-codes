        implicit real *8 (a-h,o-z)
        dimension xsout(1000),w(4000 000)
        complex *16 coefs(1000),cint,cint2,funtest
c 
        external funtest,funuser
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER b '
        READ *,b
        CALL PRIN2('b=*',b,1 )
c 
        a=0
        rlmax=20
        rlmin=10
        eps=1.0d-10
        ifinit=1
c 
        lenw=4000 000
        call expoquad(ier,a,b,rlmin,rlmax,funuser,eps,ifinit,
     1      xsout,coefs,npts,w,lenw)
c 
        call prin2('after expoqua0, coefs=*',coefs,npts*2)
        call prin2('while xsout=*',xsout,npts)
c 
c       test the obtained quadratures
c 
        m=20
        eps2=1.0d-10
c 
        call cadapgau(ier,a,b,funtest,par1,par2,m,eps2,
     1      cint,maxrec,numint)
c 
        call prin2('after cadapgau, cint=*',cint,2)
  
c 
        cint2=0
        do 2200 i=1,npts
c 
        cint2=cint2+coefs(i)*funtest(xsout(i),par1,par2)
 2200 continue
  
  
        call prin2('and cint2=*',cint2,2)
        call prin2('and cint2-cint=*',cint2-cint,2)
        call prinf('and npts=*',npts,1)
  
  
        stop
        end
c 
c 
c 
c 
c 
        function funuser(t,k,w)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1)
        complex *16 val,funuser
c 
        funuser=cos(t)
        return
        end
c 
c 
c 
c 
c 
        function funtest(t,k,par2)
        implicit real *8 (a-h,o-z)
        save
        complex *16 funuser,funtest
c 
        funtest=cos(t*2)*funuser(t,k,par)
  
        return
        end
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code and the beginning of the
c        special-purpose quadrature code proper
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine expoquad(ier,a,b,rlmin,rlmax,funuser,eps,
     1      ifinit,xsout,coefs,npts,w,lenw)
        implicit real *8 (a-h,o-z)
        save
        dimension xsout(1),w(1)
c 
c        This subroutine designs a special-purpose quadrature for
c        the integration of a certain type of functions frequently
c        encountered in the numerical complex analysis. The functions
c        are defined on intervals located in C^2 and parallel to the
c        real exis, and have the form
c 
c        f(z)=exp(rlam*x) * phi(x),                                   (1)
c 
c        with the parameter rlam on the (user-supplied) interval
c        [rlmin,rlmax], and x on the (user-supplied) interval [a,b],
c        and phi a reasonably general integrable (complex-valued)
c        function on the interval [a,b].
c 
c                      Input parameters:
c 
c  a,b - the ends of the interval along the real axis in the complex
c        plane on which the integration is to be performed
c  rlmin, rlmax - the interval on which the exponent rlam in (1) has
c        to be for the quadrature to work
c  funuser - the (user-provided) complex-valued function supplying the
c        values of the function phi in (1). The calling sequence of
c        funuser is the simplest possible:
c 
c          funuser(t)                                                 (2)
c 
c  eps - the accuracy to which the function will be integrated
c  ifinit - the (integer) parameter telling the subroutine whether
c        at the present call, the geometry has changed since the preceding
c        one. If the only thing to change since the preceding call is the
c        subroutine funuser, ifinit should be set to zero. In this case,
c        the subroutine will not construct any of the prolate machinery,
c        and will use the first npts**2 +2 * npts +10 elements of the array
c        w as input data. If the geometry (i.e. the parameters a,b,rlmin,
c        rlmax,eps) is also new, the parameter ifinit has to be set to 1.
c  lenw - the length of the user-supplied work array. Should be sufficiently
c        large.
c  w - the first npts**2 +2*npts +10 real *8 locations of this array are
c        viewed as input data if ifinit (see above) had been set to 0.
c 
c                        Output parameters:
c 
c  ier - error return code. ier=0 means successful conclusion
c                           ier=8 means that th adaptive Gaussian
c                                 quadrature inside the code failed
c                                 for some reason
c                           ier > 8 means unsuffcient memory
c  xsout - the nodes of the required quadrature
c  coefs - the (complex) weights of the quadrature
c  npts - the number of nodes in the returned quadrature
c 
c       . . . allocate memory for the construction of a "Chebychev-type"
c             special-purpose quadrature
c 
        ier=0
c 
        c=(b-a)/2 * (rlmax-rlmin)/2
        n=c+100
c 
        iv=1
        lv=n**2
c 
        icints=iv+lv
        lcints=n*2
c 
        ixs=icints+lcints
        lxs=n
c 
        iws=ixs+lxs
        lws=n
c 
        iu=iws+lws
        lu=n**2
c 
        iw=iu+lu
c 
        lleft=lenw-iw
c 
        call prinf('in expoquad, lleft=*',lleft,1)
c 
c       . . . construct the quadrature
c 
        call expoqua0(ier,a,b,rlmin,rlmax,funuser,eps,ifinit,
     1      xsout,coefs,npts,
     2      w(ixs),w(iws),w(iu),w(iv),w(icints),w(iw),lleft)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine expoqua0(ier,a,b,rlmin,rlmax,funuser,eps,
     1      ifinit,xsout,coefs,npts,
     2      xs,ws,u,v,cints,w,lenw)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),ws(1),u(1),v(1),w(1),xsout(1)
c 
        complex *16 cints(1),coefs(1),funuser
c 
        external expaquf1,funuser
c 
c        construct the nodes on the interval [-1,1]
c 
        ier=0
        c=(b-a)/2 * (rlmax-rlmin)/2
c 
        ifrefine=0
        ifeven=0
c 
        call prolexps(jer,c,eps,ifrefine,npts,xs,ws,u,v,
     1      w,lenw,keep,lused)
c 
        if(jer .ne. 0) ier=32
        if(jer .ne. 0) return
c 
c       evaluate the integrals of all prolate functions
c 
        beta2=(b+a)/2
        alpha2=(b-a)/2
        do 2500 i=1,npts
c 
        xsout(i)=alpha2*xs(i)+beta2
 2500 continue
c 
        m=20
        eps2=eps/10
        ii=expaquif(a,b,alpha,beta)
c 
        do 2200 i=1,npts
c 
        ii=expaqf1i(i-1)
c 
        call cadapgau(jer,a,b,expaquf1,funuser,w,m,eps2,
     1      cints(i),maxrec,numint)
 2200 continue
c 
        if(jer .ne. 0) ier=8
        if(jer .ne. 0) return
c 
c       determine coefficients of the quadrature formula
c 
        call expaquma(v,cints,coefs,npts)
c 
        do 2800 i=1,npts
c 
        coefs(i)=coefs(i)/funuser(xsout(i),k,w)
 2800 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine expaquma(a,x,y,n)
        implicit real *8 (a-h,o-z)
        save
        complex *16 x(1),y(1),d
c 
        real *8 a(n,n)
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
        function expaquf1(t,funuser,w)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1)
        complex *16 expaqufp,val,expaquf1,funuser
c 
c 
        x=alpha*t+beta
        val=expaqufp(t,k,w)
        expaquf1=val*funuser(t,k,w)
        return
c 
c 
c 
c 
        entry expaqf1i(k7)
        k=k7
        expaqf1i=0
        return
        end
c 
c 
c 
c 
c 
        function expaqufp(t,k,w)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1)
        complex *16 expaqufp,val
c 
c 
        x=alpha*t+beta
        call prolevv(k,x,w,val)
        expaqufp=val
c 
        return
c 
c 
c 
c 
        entry expaquif(a7,b7,alpha7,beta7)
c 
        a=a7
        b=b7
c 
        alpha=2/(b-a)
        beta=1-alpha*b
c 
        alpha7=alpha
        beta7=beta
c 
        expaquif=0
  
        return
        end
  
  
  
  
  
