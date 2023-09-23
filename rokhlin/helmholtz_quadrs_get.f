        implicit real *8 (a-h,o-z)
        dimension w(20 000 000)
c 
        external funuser,funwht,funuser2,funuser4,funuser3
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
C 
        PRINT *, 'ENTER npols'
        READ *,npols
        CALL PRINf('npols=*',npols,1 )
  
  
  
  
        nxrk=10
        nyrk=10
        nx=10
c 
        nxrk=3
        nyrk=3
        nx=3
  
  
        y1rk=0
        y2rk=1
  
        nxrk=20
        nyrk=20
        nyrk=40
        nx=30
c 
        nxrk=30
        nyrk=30
        nyrk=50
        nx=40
  
  
  
c 
cccc        nxrk=3
cccc        nyrk=3
cccc        nx=3
c 
c 
c 
c 
        epssvd=1.0d-10
        epsquadr=1.0d-6/2 /2
        epsquadr=1.0d-6/2
  
  
        epssvd=1.0d-12
        epsquadr=1.0d-9 *8
  
  
        epssvd=1.0d-6
        epsquadr=1.0d-3*2
  
  
        iijjkk1=28
        iijjkk2=64
  
        iijjkk1=52
        iijjkk2=52
  
  
        iijjkk1=5
        iijjkk2=7
  
cccc        iijjkk1=40
cccc        iijjkk2=64
  
  
  
        lw=20 000 000
  
        call helmholtz_quadrs_get(npols,iijjkk1,iijjkk2,
     1      y1rk,y2rk,nxrk,nyrk,nx,epssvd,epsquadr,
     2      w,lw)
  
        stop
        end
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code, and the beginning of
c        the quadrature generation code proper
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine helmholtz_quadrs_get(npols,iijjkk1,iijjkk2,
     1      y1rk,y2rk,nxrk,nyrk,nx,epssvd,epsquadr,
     2      w,lw)
        implicit real *8 (a-h,o-z)
        save
        dimension w(9 000 000)
c 
        external funuser,funwht,funuser2,funuser4,funuser3
c 
c        This subroutine constructs quadrature formulae for
c        integrals of the form
c 
c 
c         p * \int_0^{\infty} t*exp(-p*x*cd)/cd*J_0(p*y*t) dt        (1)
c 
c        with
c 
c                   cd=sqrt(t^2-(rk/p)^2),                           (2)
c 
c        p=Re(rk), and x, y in the regions
c 
c                   x \in [1,4]                                      (3)
c 
c                   y \in [0,4*sqrt(2)].                             (4)
c 
c        One call to this subroutine generates a collection of
c        iijjkk2-iijjkk1+1 such quadratures; each quadrature
c        is valid for values of the parameter rk in the rectangle
c        in the complex plane defined by the inequalities
c 
c 
c        Re(rk) \in [k,k+1], Im(rk) \in [y1rk,y2rk].                 (5)
c 
c        The quadratures are written in the form of FORTRAN data
c        statements on the files fort.iijjkk1, fort.iijjkk1+1,...,
c        fort.iijjkk2. For example, if iijjkk1=5, iijjkk2=7, the
c        subroutine will generate three sets of quadratures optimal
c        in the rectangular regions
c 
c        Re(rk) \in [5,6], Im(rk) \in [0,1],                         (6a)
c 
c        Re(rk) \in [6,7], Im(rk) \in [0,1],                         (6b)
c 
c        Re(rk) \in [7,8], Im(rk) \in [0,1].                         (6c)
c 
c 
c                          Input parameters:
c 
c  npols - the number of nodes into which the interval is to be
c        discretized on which y*p lives. Must be large: 300 is a
c        good number. PLEASE NOTE THAT THIS PARAMETER AFFECTS THE
C        EXECUTION TIME ONLY WEAKLY; ONE DOES NOT PAY MUCH FOR
C        MAKING IT LARGE.
c  iijjkk1, iijjkk2 - limits on Re(rk) for which the quadratures
c        will be valid (see (6a)-(6c) above)
c  y1rk, y2rk - the limits on Im(rk) between which the quadratures
c        will be valid. ALSO SEE NOTE BELOW.
c  nxrk - the number of nodes into which Re(rk) will be discretized.
c        20 is a good number
c  nyrk - the number of nodes into which Im(rk) will be discretized.
c        30 is a good number
c  nx - the number of nodes into which x will be discretized.
c        30 is a good number
c  epssvd - the accuracy to which the SVDs will be performed inside
c        the algorithm
c  epsquadr - the accuracy to which the quadratures will be
c        constructed
c  lw - the length of the work array w supplied by the user (in
c        real *8 words)
c 
c                          Output parameters: none;
c 
c        All quadratures are written on the FORTRAN files with numbers
c        100+iijjkk1, 100+iijjkk1+1,..., 100+iijjkk2; the report on
c        the execution is written on two files: fort.6 and fort.13.
C        A picky user might feel that the report is poorly designed
c        and excessively verbose; in this case, the picky user is
c        advised to write his own code.
c 
c   IMPORTANT NOTE. Whenever the quadratures are designed to be
c        valid in the region
c 
c        Re(rk) \in [k,k+1], Im(rk) \in [0,1],                       (7)
c 
c        they turn out to be valid in the region
c 
c        Re(rk) \in [k,k+1], Im(rk) \in [0,\infty]                   (8)
c 
c        (for complicated but fairly well-understood reasons).
c        However, in most of the region (8), the resulting
c        quadratures are quite suboptimal.
c 
c 
        do 3000 iijjkk=iijjkk1,iijjkk2
c 
        k=20
c 
        call prinf('originally,  lw/1000=*',lw/1000,1)
c 
        call prinf('before quaevol2, lw=*',lw,1)
c 
        call fileflush(13)
c 
        iw=100+iijjkk
c 
        x1rk=iijjkk
        x2rk=iijjkk+1
c 
        delta=0.01d0
c 
c       evaluate bx
c 
        rlammax=log(epsquadr**2)**2+x2rk**2
        rlammax=sqrt(rlammax)
        bx=rlammax/x1rk
        ax=0
  
  
        call prin2('ax as created*',ax,1)
        call prin2('bx as created*',bx,1)
  
c 
c       construct the parameters of the bump on the
c       integration path
c 
        bexp=30
        if (bx-1 .lt. 1) bexp=-log(1.0d-13)/(bx-1)**2
  
        call prin2('bexp as constructed*',bexp,1)
  
  
  
        call funuser3ini(x1rk,x2rk,nxrk,y1rk,y2rk,nyrk,
     1      nx,nfuns,delta,bexp)
  
        call prinf('and nfuns=*',nfuns,1)
  
        nsings=nfuns
  
        call funuser4ini(x1rk,x2rk,npols,delta,bexp)
  
  
        call fileflush(13)
  
c 
c        store on disk the header for this run
c 
 2100 format('c       ')
 2200 format('c        In this run,  x1rk=',e11.5,
     1    ', x2rk=',e11.5,',',/,
     2       'c                      y1rk=',e11.5,
     3    ', y2rk=',e11.5,',',/,
     6       'c                        bx=',e11.5,',',
     3    ', epsquadr=',e11.5,',',/,
     4       'c                     delta=',e21.15,',',/,
     5       'c                      bexp=',e21.15)
  
  
  
         write(iw,2100)
         write(iw,2200) x1rk,x2rk,y1rk,y2rk,bx,epsquadr,
     1       delta,bexp
  
  
        call fileflush(iw)
  
cccc        stop
  
  
        nnodes=-1
        iwrite=-23
        ifnewton=1
  
        call quaevol4(ier,ax,bx,k,funwht,funuser3,funuser4,
     1      nsings,npols*2,epssvd,epsquadr,w,lw,lused,iw,ncols,
     2      nfintegr,nnodes,iwrite,ifnewton)
  
 2500 format('        ')
 2600 format('        and again, bexp=',e21.15)
        write(iw,2500)
        write(iw,2600) bexp
        write(iw,2500)
c 
 3000 continue
  
        return
        end
c 
c 
c 
c 
c 
        subroutine funwht(x,weight)
        implicit real *8 (a-h,o-z)
c 
c 
        save
  
        weight=1
  
        return
        end
c 
c 
c 
c 
c 
  
        subroutine funuser3(rt,ifun,f)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension xsrk(1000),wxsrk(1000),xs(1000),wxs(1000),
     1      ysrk(1000),wysrk(1000),icontr(4,200 000)
c 
        complex *16 ima,rk,cd,cf,t
c 
c 
        data ima/(0.0d0,1.0d0)/
  
c 
c        unpack the arrays xsrk,ysrk,xs
c 
        ifima=icontr(1,ifun)
        ixrk=icontr(2,ifun)
        iyrk=icontr(3,ifun)
        ix=icontr(4,ifun)
c 
        xrk=xsrk(ixrk)
        yrk=ysrk(iyrk)
        x=xs(ix)
c 
  
        t=rt-delta*exp(-bexp*(rt-1)**2) * ima
  
        goto 1100
  
        call prin2('t=*',t,1)
        call prin2('xrk=*',xrk,1)
        call prin2('yrk=*',yrk,1)
        call prin2('x=*',x,1)
  
 1100 continue
c 
        rk=xrk+ima*yrk
  
cccc        call prin2('rk=*',rk,2)
c 
        p=xrk
c 
c       evaluate the function
c 
        cd=sqrt(t**2-(rk/p)**2)
        cf=exp(-p*x*cd)*t*p/cd
  
  
  
        cf=cf*(1+2*bexp*ima*(rt-1)*delta*
     1      exp(-bexp*(rt-1)**2) )
  
  
c 
        f=cf
        if(ifima .eq. 1) f=-ima*cf
  
  
cccc        f=f*sqrt(wxs(ix))
        f=f*sqrt(wxs(ix) * wxsrk(ixrk) * wysrk(iyrk) )
cccc        f=f*sqrt(wxs(ix) )
  
cccc        call prin2('f=*',f,1)
c 
        return
c 
c 
c 
c 
        entry funuser3ini(x1rk7,x2rk7,nxrk7,y1rk7,y2rk7,nyrk7,
     1      nx7,nfuns,delta7,bexp7)
c 
        x1rk=x1rk7
        x2rk=x2rk7
        nxrk=nxrk7
c 
        y1rk=y1rk7
        y2rk=y2rk7
        nyrk=nyrk7
c 
        nx=nx7
        delta=delta7
        bexp=bexp7
  
  
        call prinf('nxrk=*',nxrk,1)
        call prinf('nyrk=*',nyrk,1)
        call prinf('nx=*',nx,1)
c 
c       discretize the real and imaginary parts of rk
c 
        call interv_discr(x1rk,x2rk,nxrk,xsrk,wxsrk)
c 
        call interv_discr(y1rk,y2rk,nyrk,ysrk,wysrk)
c 
c       discretize x
c 
        x1=1
        x2=4.2
        call interv_discr(x1,x2,nx,xs,wxs)
c 
  
        call prin2('in funuser3ini, xsrk=*',xsrk,nxrk)
        call prin2('in funuser3ini, ysrk=*',ysrk,nyrk)
        call prin2('in funuser3ini, xs=*',xs,nx)
c 
c        construct the control arrays
c 
        ii=0
        do 3000 ifima=0,1
c 
        do 2800 i=1,nxrk
        do 2600 j=1,nyrk
        do 2400 k=1,nx
c 
        ii=ii+1
c 
        icontr(1,ii)=ifima
        icontr(2,ii)=i
        icontr(3,ii)=j
        icontr(4,ii)=k
c 
 2400 continue
 2600 continue
 2800 continue
 3000 continue
c 
        nfuns=ii
c 
        return
        end
  
c 
c 
c 
c 
c 
  
        subroutine interv_discr(a,b,n,xs,ws)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension xs(1),ws(1)
c 
        itype=1
        call legeexps(itype,n,xs,u,v,ws)
c 
        alpha=(b-a)/2
        beta=(b+a)/2
c 
        do 1200 i=1,n
c 
        xs(i)=xs(i)*alpha+beta
        ws(i)=ws(i)*alpha
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
  
        subroutine funuser4(rt,i,f)
        implicit real *8 (a-h,o-z)
        save
        dimension pys(1000),wpys(1000)
c 
        complex *16 t,ima,cd
c 
        data ima/(0.0d0,1.0d0)/
c 
        t=rt-delta*exp(-bexp*(rt-1)**2) * ima
  
cccc        call prin2('in funuser4, t=*',t,2)
  
  
        j=i
        if(j .gt. npols) j=j-npols
  
        cd=cos(pys(j)*t)
        f=cd
        if(i .gt. npols) f=-ima*cd
  
  
        f=f*sqrt(wpys(j))
        return
c 
c 
c 
c 
        entry funuser4ini(x1rk,x2rk,npols7,delta7,bexp7)
c 
        npols=npols7
        delta=delta7
        bexp=bexp7
        ymin=0
cccc        ymax=4*sqrt(2.0d0)*1.05
        ymax=4*sqrt(2.0d0)
  
        pymin=0
        pymax=x2rk*ymax
c 
        call interv_discr(pymin,pymax,npols,pys,wpys)
  
        call prin2('in funuser4ini, pys=*',pys,npols)
  
        return
        end
c 
c 
c 
c 
c 
  
        subroutine quad_copy(bexp,n,xs,ws,cxs,cws)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),ws(1)
        complex *16 cxs(1),cws(1),ima
c 
        data ima/(0.0d0,1.0d0)/
c 
        delta=0.01d0
        do 1200 i=1,n
c 
        rt=xs(i)
        cxs(i)=rt-delta*exp(-bexp*(rt-1)**2) * ima
        cws(i)=ws(i)*
     1      (1+2*bexp*ima*(rt-1)*delta*exp(-bexp*(rt-1)**2) )
c 
 1200 continue
c 
        return
        end
  
  
