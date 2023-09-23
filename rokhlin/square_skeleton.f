        implicit real *8 (a-h,o-z)
c 
        real *8 zsin(2,10 000),
     1      zsout(2,100 000),wsin(10 000),wsout(10000),
     5      rlensin(10000),rlensout(10000),
     6      rlams(10000),
     7      zsin_skel(2,10000),zsout_skel(2,10000)
c 
        complex *16 rk,amatr(1000 000),
     1      w(10 000 000),
     2      v_big(3000 000),u_big(1000 000),
     4      whtsin_skel(10000),whtsout_skel(10000)
c 
        external poteval2
  
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
         PRINT *, 'ENTER rk'
         READ *,rrk
         CALL PRIN2('rrk=*',rrk,1 )
C 
         PRINT *, 'ENTER nts'
         READ *,nts
         CALL PRINf('nts=*',nts,1 )
c 
         coef=3
  
         rk=rrk
  
         gap=0.1d0
         gap=0.2d0
cccc         gap=0.6d0
  
         eps=1.0e-12
  
c 
c       construct the svd
c 
         eps2=eps
        iplot=21
  
  
c 
        lenw=10 000 000
  
  
  
  
  
        call square_skeleton(ier,coef,poteval2,rk,
     1    gap,eps,nts,iplot,zsin_skel,whtsin_skel,
     2      zsout_skel,whtsout_skel,nvects,
     3      zsin,wsin,rlensin,nin,zsout,wsout,
     4      rlensout,nout,w,lenw,lused)
  
        call prinf('after square_skeleton, lused=*',lused,1)
  
  
        iw=28
        call zquaplot2(iw,zsout,nout,2,
     1      zsin,nin,2,'both squares*')
  
  
        call prin2('after square_skeleton, whtsout_skel=*',
     1      whtsout_skel, nvects*2)
  
        call prin2('after square_skeleton, whtsin_skel=*',
     1      whtsin_skel, nvects*2)
c 
        lenw=10 000 000
  
ccccc           goto 3000
  
        call square_basis(ier,poteval2,rk,nts,
     1      coef,gap,eps,iplot,
     2      zsin,wsin,rlensin,nin,zsout,wsout,rlensout,nout,
     3      rlams,nvects2,u_big,v_big,w,lenw)
  
        call prin2('after square_basis, rlams=*',rlams,nvects2*2)
        iw0=100
        call plotem(iw0,v_big,rlensout,nout/2,nvects2)
  
        iw0=200
        call plotem(iw0,u_big,rlensin,nin/2,nvects2)
  
        nn=nvects
        call skeltest(amatr,nn,nvects2,u_big,zsin_skel,zsin,nin)
        call skeltest(amatr,nn,nvects2,v_big,zsout_skel,zsout,nout)
c 
 3000 continue
  
        iw=31
        call zquaplot2(iw,zsout_skel,nvects,2,
     1      zsin_skel,nvects*8,2,'both skeletons expanded*')
  
        iw=32
        call zquaplot2(iw,zsout_skel,nvects/8,2,
     1      zsin_skel,nvects/8,2,'1/8-th of both skeletons*')
  
  
cccc        stop
c 
c        test the obtained quadratures
c 
        call quadtest(zsin,wsin,nin,zsout,wsout,nout,
     1      zsin_skel,whtsin_skel,nvects,zsout_skel,
     2      whtsout_skel,rk)
  
  
        call prinf('and nvects=*',nvects,1)
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine quadtest(zsin,wsin,nin,zsout,wsout,nout,
     1      zsin_skel,wsin_skel,nvects,zsout_skel,wsout_skel,
     2      rk)
        implicit real *8 (a-h,o-z)
        save
        complex *16 wsin_skel(1),wsout_skel(1),cd,cint,cint2,
     1      rk
c 
        dimension zsin(2,1),wsin(1),zsout(2,1),wsout(1),
     1      zsout_skel(2,1),zsin_skel(2,1)
c 
c        for all points on the outer frame, evaluate the
c        integrals of their potentials on the inner frame
c 
        call prin2('in quadtest, rk=*',rk,2)
  
cccc        stop
  
        errmax=0
c 
        do 2000 i=1,nout
c 
        cint=0
        do 1200 j=1,nin/2
c 
        call poteval2(zsout(1,i),zsin(1,j),rk,cd)
c 
        cint=cint+wsin(j)*cd
 1200 continue
c 
cccc        call prin2('cint=*',cint,2)
  
c 
        cint2=0
        do 1400 j=1,nvects
  
        call poteval2(zsout(1,i),zsin_skel(1,j),rk,cd)
c 
        cint2=cint2+wsin_skel(j)*cd
 1400 continue
c 
cccc        call prin2('and cint2=*',cint2,2)
cccc        call prin2('and cint2/cint=*',cint2/cint,2)
c 
        dd=abs(cint2-cint)
  
        if(errmax .lt. dd) errmax=dd
  
 2000 continue
  
        call prin2('and errmax=*',errmax,1)
  
  
cccc        stop
  
c 
c        for all points on the inner frame, evaluate the
c        integrals of their potentials on the outer frame
c 
cccc        call prin2('in quadtest, zsin=*',zsin,nin*2)
cccc        call prin2('in quadtest, zsout=*',zsout,nout*2)
  
  
cccc        call prin2('in quadtest, zsout_skel=*',zsout_skel,nvects*2)
cccc        call prin2('in quadtest, wsout_skel=*',wsout_skel,nvects*2)
  
  
cccc        stop
  
        errmax2=0
  
        do 3000 i=1,nin
c 
        cint=0
  
        do 2200 j=nout/2+1,nout
  
        call poteval2(zsin(1,i),zsout(1,j),rk,cd)
c 
        cint=cint+wsout(j)*cd
 2200 continue
c 
cccc        call prin2('cint=*',cint,2)
  
c 
        cint2=0
        do 2400 j=1,nvects
  
        call poteval2(zsin(1,i),zsout_skel(1,j),rk,cd)
c 
        cint2=cint2+wsout_skel(j)*cd
 2400 continue
c 
cccc        call prin2('and cint2=*',cint2,2)
cccc        call prin2('and cint2/cint=*',cint2/cint,2)
c 
        dd=abs(cint2-cint)
  
        if(errmax2 .lt. dd) errmax2=dd
  
  
cccc        stop
 3000 continue
  
        call prin2('and errmax2=*',errmax2,1)
  
cccc        stop
  
        return
        end
c 
c 
c 
c 
c 
        subroutine plotem(iw0,v,rlens,n,ncols)
        implicit real *8 (a-h,o-z)
        save
        complex *16 v(n,ncols)
        dimension rlens(1)
c 
c       plot all of the functions in the array v
c 
        do 2200 i=1,ncols
c 
        iw=iw0+i
        itype=3
        call comfungraph(iw,rlens,v(1,i),n,itype,
     1      'singular function*')
  
 2200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine skeltest(a,n,nvects,u,zs4,zs,nout)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,1),u(nout,nvects)
        dimension zs(2,1),zs4(2,1)
c 
        complex *16 uu(1000 000),vv(1000 000),ww(1000 000),
     1      ss(10000)
c 
c        construct the matrix
c 
        call prinf('in skeltest, n=*',n,1)
  
        do 1400 i=1,n
  
cccc        call prinf('in skeltest, i=*',i,1)
c 
c        find the true location for this point
c 
        call square_symm_fnd(zs4(1,i),zs,nout,i2)
c 
        do 1200 j=1,nvects
c 
        a(i,j)=u(i2,j)
 1200 continue
c 
 1400 continue
c 
c        skeletonize the matrix bstar
c 
        eps=1.0d-10
        ifpivot=1
c 
        lenww=1000 000
        call csvdpiv(ier,a,n,nvects,uu,vv,ss,ncols,eps,
     1      ww,lenww,ltot)
  
  
  
  
        call prin2('ss=*',ss,ncols*2)
        call prinf('while ncols =*',ncols,1)
        call prinf('and nvects =*',nvects,1)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine comfungraph(iw,ts,zs,n,itype,mes)
        implicit real *8 (a-h,o-z)
        save
        dimension zs(2,1),ts(1),xs(10 000),ys(10 000)
        character *1 mes(1)
c 
        do 1200 i=1,n
c 
        xs(i)=zs(1,i)
        ys(i)=zs(2,i)
 1200 continue
c 
        call quagraph2(iw,ts,xs,n,itype,
     1      ts,ys,n,itype,mes)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine poteval2(z1,z2,rk,pot)
        implicit real *8 (a-h,o-z)
        save
        dimension z1(1),z2(1)
        complex *16 rk,cd,pot,h0,h1
  
c 
        d=(z1(1)-z2(1))**2+(z1(2)-z2(2))**2
        d=sqrt(d)
c 
  
cccc        call prin2('in poteval, d=*',d,1)
  
        cd=d*rk
  
  
cccc        call prin2('in poteval, cd=*',cd,2)
cccc        call prin2('in poteval, rk=*',rk,2)
  
        ifexpon=1
c 
        call hank103(cd,h0,h1,ifexpon)
c 
        pot=h0
c 
        return
        end
c 
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code, and the beginning
c        of the code constructing the skeletonization of the
c        interaction operator between two "concentric" squares
c        in the plane
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c        This file contains only one user-callable subroutine, which
c        can be found below
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        subroutine square_skeleton(ier,coef,poteval,rk,
     1    gap,eps,nts,iplot,zsin_skel,whtsin_skel,
     2      zsout_skel,whtsout_skel,npts,
     3      zsin,wsin,rlensin,nin,zsout,wsout,
     4      rlensout,nout,w2,lenw,lused)
        implicit real *8 (a-h,o-z)
        save
        dimension zsin(2,1),zsout(2,1),wsin(1),wsout(1),
     1      rlensin(1),rlensout(1),w2(1)
c 
        complex *16 rk(1)
c 
c        This subroutine skeletonizes the interactions between
c        two "concentric" squares in R^2. Its expected use is
c        for the Helmholtz interactions, and this is the regime
c        where it has been tested.
c 
c                    Input parameters:
c 
c 
c  coef - the size of the inner layer of the outer square (see
c        general description above)
c  poteval - the user-supplied subroutine defining the interaction
c        between the squares. The calling sequence of poteval is
c 
c                call poteval(z1,z2,rk,f),
c 
c        where the input parameters z21, z2 are the two points
c        between which the interaction is to be calculated,
c        the input parameter rk is whatever parameter poteval
c        might need (supplied by the user), and the output
c        parameter f is the (complex) interaction between the
c        points z1, z2
c 
c 
c  rk - the parameter (real, complex, whatever) used by the
c        subroutine poteval
c  gap - the gap between the inner and outer layers of each
c        square
c  eps - the accuracy to which the SVDs will be performed
c  nts - the number of nodes into which each of the layers
c        comprising the INNER square will be subdivided ON
c        EACH OF THE SIDES OF THE SQUARE. Thus, the total number
c        of nodes discretizing the inner square will be 8*nts.
c  iplot - integer parameter telling the subroutine where to plot
c        the three pictures it wants to plot (the plots are in
c        GNUPLOT-readable format). It will plot the data on gn*
c        files, with *=iplot, iplot+1,iplot+2. Setting iplot=0
c        will suppress the plotting.
c  lenw - the length (in real *8 locations) of the work array
c        w2 supplied by the user
c 
c                Output parameters:
c 
c  ier - error return code
c  zsin_skel - the skeletonized nodes on the inner square
c  whtsin_skel - the weights corresponding to the nodes in
c         the array zsin_skel. PLEASE NOTE THAT THE WEIGHTS IN
C         THIS ARRAY ARE COMPLEX!!!
c  zsOUT_skel - the skeletonized nodes on the OUTER square
c  whtsout_skel - the weights corresponding to the nodes in
c         the array wsout_skel. PLEASE NOTE THAT THE WEIGHTS IN
C         THIS ARRAY ARE COMPLEX!!!
c  npts - the number of skeletonized nodes on each of the square
c  zsin - the nodes on the inner square
c  wsin - the weights corresponding to the nodes in the array wsin
c  rlensin - the arc length on the inner square (distance along the
c        surface from the upper right corner, tabulated at the
c        nodes zsin; expected to be used mostly for plotting
c  nin - the number of nodes on the inner square (number of elements
c        in arrays zsin, wsin, rlensin). Will be set to nts*8
c  zsout - the nodes on the outer square
c  wsout - the weights corresponding to the nodes in the array wsout
c  rlensout - the arc length on the outer square (distance along the
c        surface from the upper right corner, tabulated at the
c        nodes zsout; expected to be used mostly for plotting
c  nout - the number of nodes on the outer square (number of elements
c        in arrays zsout, wsout, rlensout). Will be set to nts*16
c  lused - the amount of space in array w2 (in real *8 words)
c        actually used by the subroutine
c 
c                     Work array:
c 
c  w2 - should be reasonably large
c 
c       . . . allocate memory
c 
        ier=0
c 
        nin=nts*8
        nout=nts*16
c 
        nout2=nts*4
        nin2=nts*2
c 
        nin3=nin2/2
        nout3=nout2/2
  
        iamatr=1
        lamatr=nin3*nout3*2+10
c 
        iu=iamatr+lamatr
        lu=nin3*nout3*2+10
c 
        iv=iu+lu
        lv=nin3*nout3*2+10
c 
        iu2=iv+lv
        lu2=nin3*nout3*2+10
c 
        iwhtsin=iu2+lu2
        lwhtsin=nin3*2+10
c 
        iwhtsout=iwhtsin+lwhtsin
        lwhtsout=nout3*2+10
c 
c 
        izsin2=iwhtsout+lwhtsout
  
        lzsin2=nin2*2+10
c 
        iwsin2=izsin2+lzsin2
        lwsin2=nin2+10
c 
        izsin3=iwsin2+lwsin2
        lzsin3=nin3*2+10
c 
        iwsin3=izsin3+lzsin3
        lwsin3=nin3+10
c 
        izsout2=iwsin3+lwsin3
        lzsout2=nout2*2+10
c 
        iwsout2=izsout2+lzsout2
        lwsout2=nout2+10
c 
        izsout3=iwsout2+lwsout2
        lzsout3=nout2*2+10
c 
        iwsout3=izsout3+lzsout3
        lwsout3=nout3+10
c 
        iipivots=iwsout3+lwsout3
        lpivots=nout3+10
c 
        irnorms=iipivots+lpivots
        lrnorms=nout3+10
c 
        irlams=irnorms+lrnorms
        lrlams=nout3*2+10
c 
        ltot=irlams+lrlams
c 
        iw=ltot+1
  
        lleft=lenw-iw
c 
        call square_skeleton0(coef,poteval,rk,
     1      gap,eps,nts,iplot,zsin_skel,whtsin_skel,
     2      zsout_skel,whtsout_skel,npts,
     3      zsin,wsin,rlensin,nin,zsout,wsout,rlensout,nout,
c 
     4      w2(iamatr),w2(iu),w2(iv),w2(iw),lleft,lused,
     5      w2(iu2),w2(iwhtsin),w2(iwhtsout),
     6      w2(izsin2),w2(izsout2),w2(izsin3),w2(izsout3),
     7      w2(iwsin2),w2(iwsout2),w2(iwsin3),w2(iwsout3),
     8      w2(iipivots),w2(irnorms),w2(irlams) )
c 
        lused=iw+lused*2+10
        return
        end
c 
c 
c 
c 
c 
        subroutine square_skeleton0(coef,poteval,rk,
     1      gap,eps,nts,iplot,zsin_skel,whtsin_skel,
     2      zsout_skel,whtsout_skel,npts,
     3      zsin,wsin,rlensin,nin,zsout,wsout,rlensout,nout,
     4      amatr,u,v,w,lenw,lused,u2,whtsin,whtsout,
     5      zsin2,zsout2,zsin3,zsout3,wsin2,wsout2,
     6      wsin3,wsout3,ipivots,rnorms,rlams)
        save
  
        implicit real *8 (a-h,o-z)
        dimension zsin(2,1),zsout(2,1),wsin(1),wsout(1),
     2      zsin2(2,1),zsout2(2,1),wsin2(1),wsout2(1),
     4      zsin3(2,1),zsout3(2,1),wsin3(1),
     5      wsout3(1),rlensin(1),rlensout(1),
     6      rlams(1),ipivots(1),rnorms(1)
c 
        complex *16 rk,amatr(1),u(1),v(1),w(1),
     1      whtsin(1),u2(1),whtsout(1)
c 
        ntsout=nts*2
c 
        ntsin=nts
        ntsout=nts*2
c 
c        construct the big matrix of interactions and its SVD
c        via the accelerated scheme
c 
        call square_basis04(ier,poteval,rk,nts,
     1      coef,gap,eps,
     2      zsin,wsin,rlensin,nin,zsout,wsout,rlensout,nout,
     3      rlams,nvects,w,lenw,lused,
     4      amatr,u,v,
     5      zsin2,wsin2,zsin3,wsin3,
     6      zsout2,wsout2,zsout3,wsout3,zs2,iplot)
c 
c        . . . skeletonize
c 
        eps2=1.0d-10
c 
        nout3=nout/8
        nin3=nin/8
c 
        ifadj=1
        call square_skelit(v,nout3,nvects,w,rlams,zsout3,zsout2,
     1      ncols1,eps2,wsout3,u2,whtsout,ifadj,
     2      ipivots,rnorms)
c 
        ifadj=0
        call square_skelit(u,nin3,nvects,w,rlams,zsin3,zsin2,
     1      ncols2,eps2,wsin3,u2,whtsin,ifadj,
     2      ipivots,rnorms)
c 
c       expand the obtained discretizations to the whole
c       square
c 
        call square_expand7(zsin2,zsin_skel,whtsin,
     1      whtsin_skel,ncols2)
c 
        call square_expand7(zsout2,zsout_skel,whtsout,
     1      whtsout_skel,ncols1)
c 
c       perform the final reordering of the nodes
c 
        call square_final_reord(zsin_skel,whtsin_skel,
     1      ncols2*8,zsin,nin,ipivots)
c 
        call square_final_reord(zsout_skel,whtsout_skel,
     1      ncols2*8,zsout,nout,ipivots)
c 
        npts=nvects*8
c 
        return
        end
c 
c 
c 
c 
c 
  
        subroutine square_expand7(zs,zs2,ws,ws2,nts)
        implicit real *8 (a-h,o-z)
        save
        dimension zs(2,1),zs2(2,1)
        dimension ws(2,1),ws2(2,1)
c 
c        implement the reflection in the diagonal of the square
c 
        do 1100 i=1,nts
c 
        zs2(1,i)=zs(1,i)
        zs2(2,i)=zs(2,i)
c 
        zs2(1,i+nts)=-zs(2,i)
        zs2(2,i+nts)=-zs(1,i)
c 
c 
        ws2(1,i)=ws(1,i)
        ws2(2,i)=ws(2,i)
c 
        ws2(1,i+nts)=ws(1,i)
        ws2(2,i+nts)=ws(2,i)
c 
 1100 continue
c 
c        implement the x-reflection
c 
        do 1200 i=1,nts*2
        zs2(1,i+nts*2)=-zs2(1,i)
        zs2(2,i+nts*2)=zs2(2,i)
c 
        ws2(1,i+nts*2)=ws2(1,i)
        ws2(2,i+nts*2)=ws2(2,i)
c 
 1200 continue
c 
c        implement the y-reflection
c 
         do 1400 i=1,nts*4
c 
         zs2(1,i+nts*4)=zs2(1,i)
         zs2(2,i+nts*4)=-zs2(2,i)
c 
         ws2(1,i+nts*4)=ws2(1,i)
         ws2(2,i+nts*4)=ws2(2,i)
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
        subroutine square_final_reord(zs,ws,n,zsbig,nbig,inds)
        implicit real *8 (a-h,o-z)
        save
        dimension zs(2,1),ws(2,1),zsbig(2,1),inds(1)
c 
c       find the original locations of all points in array zs
c 
        do 2200 i=1,n
c 
        call square_symm_fnd(zs(1,i),zsbig,nbig,i2)
c 
        inds(i)=i2
 2200 continue
c 
c        sort the original nodes and weights, to position
c        them on the square in the appropriate order
c 
        do 2800 i=1,n
        do 2600 j=1,n-1
c 
        if(inds(j) .lt. inds(j+1)) goto 2600
c 
        d1=zs(1,j)
        d2=zs(2,j)
c 
        zs(1,j)=zs(1,j+1)
        zs(2,j)=zs(2,j+1)
c 
        zs(1,j+1)=d1
        zs(2,j+1)=d2
c 
c 
        d1=ws(1,j)
        d2=ws(2,j)
c 
        ws(1,j)=ws(1,j+1)
        ws(2,j)=ws(2,j+1)
c 
        ws(1,j+1)=d1
        ws(2,j+1)=d2
c 
        ii=inds(j)
        inds(j)=inds(j+1)
        inds(j+1)=ii
c 
 2600 continue
 2800 continue
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine square_skelit(b,n,m,bstar,rlams,zs,zs2,
     1      ncols,eps,ws,u2,whts,ifadj,
     2      ipivots,rnorms)
       implicit real *8 (a-h,o-z)
        save
        complex *16 b(n,m),bstar(m,n),rlams(1),
     1      cd,u2(m,m),whts(1)
        dimension ipivots(1),rnorms(1),
     1      zs(2,1),zs2(2,1),ws(1)
c 
c        transpose b
c 
        do 1400 i=1,m
        do 1200 j=1,n
c 
        bstar(i,j)=b(j,i)*sqrt(ws(j))
 1200 continue
 1400 continue
c 
c        skeletonize the matrix bstar
c 
        ifpivot=1
        call square_grm(bstar,m,n,rnorms,eps,ncols,ifpivot,
     1      ipivots)
c 
c       select the nodes and plot them
c 
        do 2200 i=1,ncols
c 
        j=ipivots(i)
        zs2(1,i)=zs(1,j)
        zs2(2,i)=zs(2,j)
c 
        if(ifadj .eq. 0) goto 2050
c 
        do 2000 jj=1,m
c 
        u2(jj,i)=conjg(b(j,jj))
 2000 continue
        goto 2200
c 
 2050 continue
c 
       do 2100 jj=1,m
c 
        u2(jj,i)=b(j,jj)
 2100 continue
c 
 2200 continue
c 
c        construct the weights corresponding to the
c        constructed nodes
c 
c        . . . integrate the basis functions
c 
        j1=1
        j2=n/2
c 
        if(ifadj .eq. 1) then
            j1=n/2+1
            j2=n
        endif
c 
        do 2600 i=1,m
c 
        cd=0
        do 2400 j=j1,j2
c 
        if(ifadj .eq. 0) cd=cd+ws(j)*b(j,i)
        if(ifadj .eq. 1) cd=cd+ws(j)*conjg(b(j,i))
 2400 continue
c 
        whts(i)=cd
 2600 continue
c 
c       . . . find the weights
c 
        call cqrsolv(u2,m,whts,rcond)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine square_basis04(ier,poteval,rk,nts,
     1      coef,gap,eps,
     2      zsin,wsin,rlensin,nin,zsout,wsout,rlensout,nout,
     3      rlams,nvects,w,lenw,lused,
     4      amatr,u,v,
     5      zsin2,wsin2,zsin3,wsin3,
     6      zsout2,wsout2,zsout3,wsout3,zs2,iplot)
        implicit real *8 (a-h,o-z)
c 
        save
        real *8 zsin(2,1),zsout(2,1),wsin(1),wsout(1),
     1      rlensin(1),rlensout(1),zsin2(2,1),zsout2(2,1),
     3      wsin2(1),wsout2(1),zsin3(2,1),zsout3(2,1),
     4      wsin3(1),wsout3(1),zs2(2,1)
c 
        complex*16 rk,amatr(1),u(1),v(1),w(1),rlams(1)
c 
c       discretize the boundaries of the inner and outer
c       squares
c 
        ier=0
c 
        ntsin=nts
        ntsout=nts*2
c 
        call square_frames_crea(ntsin,ntsout,coef,gap,
     1      zsin,wsin,nin,zsout,wsout,nout,
     2      zsin2,wsin2,nin2,zsout2,wsout2,nout2,
     3      zsin3,wsin3,nin3,zsout3,wsout3,nout3,
     4      rlensin,rlensout,w,iplot)
c 
        nin=ntsin*8
        nout=ntsout*8
c 
c        Construct the totally symmetrized (\times 8, and all "+"
c        signs) singular vectors
c 
        call square_basis4(jer,zsin3,wsin3,nin3,poteval,
     1      zsout3,wsout3,nout3,amatr,rk,eps,nvects,
     2      u,v,rlams,w,lenw,lused)
c 
        if(jer .ne. 0) ier=32
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine square_basis4(ier,zsin,wsin,nin,poteval,
     1      zsout,wsout,nout,a,rk,eps,ncols,
     2      u,v,s,w,lenw,lused)
        implicit real *8 (a-h,o-z)
        save
        dimension zsin(2,1),zsout(2,1),wsin(1),wsout(1),
     1      z1(3),z2(3),zz2(3)
c 
        complex *16 rk(1),w(1),u(nin,1),
     1      v(nout,1),s(1),a(nin,nout),cd00,cd10,cd01,cd11,
     2      cdd00,cdd10,cdd01,cdd11
c 
c        construct the matrix of interactions between the
c        squares
c 
        do 1600 i=1,nin
        do 1400 j=1,nout
c 
        z1(1)=zsin(1,i)
        z1(2)=zsin(2,i)
c 
        z2(1)=zsout(1,j)
        z2(2)=zsout(2,j)
c 
        call poteval(z2,z1,rk,cd00)
c 
        z2(1)=-zsout(1,j)
        z2(2)=zsout(2,j)
  
        call poteval(z2,z1,rk,cd10)
c 
        z2(1)=-zsout(1,j)
        z2(2)=-zsout(2,j)
c 
        call poteval(z2,z1,rk,cd11)
c 
        z2(1)=zsout(1,j)
        z2(2)=-zsout(2,j)
        call poteval(z2,z1,rk,cd01)
c 
        z1(1)=zsin(1,i)
        z1(2)=zsin(2,i)
c 
        zz2(1)=-zsout(2,j)
        zz2(2)=-zsout(1,j)
  
        z2(1)=zz2(1)
        z2(2)=zz2(2)
        call poteval(z2,z1,rk,cdd00)
c 
        z2(1)=-zz2(1)
        z2(2)=zz2(2)
        call poteval(z2,z1,rk,cdd10)
c 
        z2(1)=-zz2(1)
        z2(2)=-zz2(2)
        call poteval(z2,z1,rk,cdd11)
c 
        z2(1)=zz2(1)
        z2(2)=-zz2(2)
        call poteval(z2,z1,rk,cdd01)
c 
        a(i,j)=cd00+cd01+cd10+cd11+
     1         (cdd00+cdd01+cdd10+cdd11)
c 
        a(i,j)=a(i,j)*sqrt(wsin(i)*wsout(j))
c 
 1400 continue
 1600 continue
c 
c      svd the matrix
c 
        call csvdpiv(ier,a,nin,nout,u,v,s,ncols,eps,
     1      w,lenw,lused)
c 
c       Compensate for the square roots of weights. Also,
c       multiply each of the singular vectors by a complex
c       number of absolute value 1, chosen in such a manner
c       as to maximize the real part of each singular vector
c 
        d=0
        do 2600 i=1,ncols
c 
        do 2200 j=1,nin
c 
        u(j,i)=u(j,i)/sqrt(wsin(j))
c 
        d=d+u(j,i)*conjg(u(j,i))*wsin(j)
 2200 continue
c 
        d=1/sqrt(d)
        do 2300 j=1,nin
c 
        u(j,i)=u(j,i)*d
 2300 continue
  
        call square_cfunc_rotate(u(1,i),wsin,u(1,i),nin)
c 
        d=0
        do 2400 j=1,nout
c 
        v(j,i)=v(j,i)/sqrt(wsout(j))
c 
        d=d+v(j,i)*conjg(v(j,i))*wsout(j)
 2400 continue
c 
        d=1/sqrt(d)
        do 2500 j=1,nout
c 
        v(j,i)=v(j,i)*d
 2500 continue
c 
        call square_cfunc_rotate(v(1,i),wsout,v(1,i),nout)
c 
 2600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine square_grm(b,n,m,rnorms,eps,ncols,ifpivot,
     1      ipivots)
        implicit real *8 (a-h,o-z)
        save
        complex *16 b(n,m),cd
        dimension rnorms(1),ipivots(1)
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
        d=d+b(j,i)*dconjg(b(j,i))
 1200 continue
        rnorms(i)=sqrt(d)
cccc        rnorms(i)=d
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
        if(ifpivot .eq. 0) goto 2700
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
        cd=b(j,i)
        b(j,i)=b(j,ipivot)
        b(j,ipivot)=cd
 2600 continue
c 
        iijj=ipivots(i)
        ipivots(i)=ipivots(ipivot)
        ipivots(ipivot)=iijj
c 
        d=rnorms(ipivot)
        rnorms(ipivot)=rnorms(i)
        rnorms(i)=d
 2700 continue
c 
c       orthogonalize the i-th column to all preceeding ones
c 
        if(2 .ne. 3) goto 2790
        if(i .eq. 1) goto 2790
        do 2780 j=1,i-1
c 
        call cleascap(b(1,i),b(1,j),n,cd)
c 
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*cd
 2770 continue
 2780 continue
 2790 continue
c 
c       normalize the i-th column
c 
        call cleascap(b(1,i),b(1,i),n,cd)
c 
c       if the norm of the remaining part of the matrix
c       is sufficiently small - terminate the process
c 
        d=cd
        if( (ifpivot .ne. 0) .and. (d .lt. thresh) ) return
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
        if(rnorms(j) .lt. thresh/100) goto 3200
c 
        call cleascap(b(1,i),b(1,j),n,cd)
        cd=dconjg(cd)
c 
         rrn=0
        do 3000 l=1,n
        b(l,j)=b(l,j)-b(l,i)*cd
        rrn=rrn+b(l,j)*dconjg(b(l,j))
 3000 continue
        rnorms(j)=dsqrt(rrn)
 3200 continue
 3400 continue
 4000 continue
  
 4200 continue
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine square_cfunc_rotate(fin,w,fout,n)
        implicit real *8 (a-h,o-z)
        save
        complex *16 fin(1),fout(1),cd,ima
        dimension w(1)
c 
        data ima/(0.0d0,1.0d0)/
c 
c        determine the rotation that will maximize the real
c        part of the function fin
c 
        call square_cfunc_rotate_fnd(fin,n,w,theta)
c 
        cd=exp(ima*theta)
c 
        do 1200 i=1,n
        fout(i)=fin(i)*cd
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine square_cfunc_rotate_fnd(f,n,w,theta)
        implicit real *8 (a-h,o-z)
        save
        real *8 f(2,1),w(1)
c 
c        determine the rotation that will maximize the real
c        part of the function fin
c 
        alpha=0
        beta=0
        gamma=0
c 
        do 1200 i=1,n
c 
        alpha=alpha+f(1,i)**2*w(i)
        gamma=gamma+f(2,i)**2*w(i)
        beta=beta+f(1,i)*f(2,i)*w(i)
 1200 continue
c 
        theta2=2*beta/(gamma-alpha)
        theta=atan(theta2)/2
c 
        d1=alpha*cos(theta)**2-2*beta*cos(theta)*sin(theta)+
     1      gamma*sin(theta)**2
c 
        done=1
        pi=atan(done)*4
c 
        theta2=theta+pi/2
c 
        d2=alpha*cos(theta2)**2-2*beta*cos(theta2)*sin(theta2)+
     1      gamma*sin(theta2)**2
c 
        theta3=theta-pi/2
c 
        d3=alpha*cos(theta3)**2-2*beta*cos(theta3)*sin(theta3)+
     1      gamma*sin(theta3)**2
c 
        d=d1
        if(d2 .gt. d) then
            d=d2
            theta=theta2
        endif
c 
        if(d3 .gt. d) then
            d=d3
            theta=theta3
        endif
c 
        return
        end
