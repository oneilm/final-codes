        implicit real *8 (a-h,o-z)
        real *8 zsin(2,10 000),zsout(2,100 000),
     1      wsin(10 000),wsout(10000),w(2000 000),
     2      zsout_big(2,10000),zsin_big(2,10000)
c 
        complex *16 rk,amatr(2000 000),ubig(1000 000),
     1      vbig(1000 000),w2(10 000 000)
C 
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
c 
         gap=0.1d0
         gap=0.2d0
c 
c        construct both the inner and the outer squares
c 
        ntsin=nts
        ntsout=nts*2
  
        call rects_cre(gap,ntsin,ntsout,coef,
     1      zsin,wsin,nin,zsout,wsout,nout,w)
  
  
        iw=22
        call zquaplot2(iw,zsin,nin,2,zsout,nout,2,
     1      'both rectangles*')
c 
c        create the svd
c 
        call creasvd(zsin,nin,zsout,nout/2,wsin,wsout,amatr,rk,
     1      ubig,vbig,ncolsbig)
  
  
  
         eps=1.0d-12
  
         lenw=20 000 000
         call rects_skeleton(ier,poteval2,rk,nts,coef,gap,eps,
     1       zsin_big,zsout_big,npts,w2,lenw,lused)
  
  
        call prinf('after rects_skeleton, lused=*',lused,1)
        iw=27
        call zquaplot2(iw,zsin_big,npts,2,
     1      zsout_big,npts,2,
     2      'after rects_skeleton, both skeletons *')
  
        call skeltest(amatr,npts,ncolsbig,ubig,zsin_big,zsin,nin)
        call skeltest(amatr,npts,ncolsbig,vbig,
     1      zsout_big,zsout,nout/2)
  
        stop
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
  
cccc  call prinf('and i2=*',i2,1)
c 
        do 1200 j=1,nvects
c 
        a(i,j)=u(i2,j)
 1200 continue
c 
 1400 continue
  
  
cccc        call prin2('a as created*',a,nvects*n*2)
  
cccc        stop
  
c 
c        skeletonize the matrix bstar
c 
        eps=1.0d-10
        ifpivot=1
cccc        call cleasgrm7(a,n,nvects,rnorms,eps,ncols,ifpivot,
cccc     1      ipivots)
  
        lenww=1000 000
        call csvdpiv(ier,a,n,nvects,uu,vv,ss,ncols,eps,
     1      ww,lenww,ltot)
  
  
  
  
  
  
        call prin2('ss=*',ss,ncols*2)
        call prinf('while ncols =*',ncols,1)
        call prinf('and nvects =*',nvects,1)
        call prinf('and n =*',n,1)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine creasvd(zsin,nin,zsout,nout,wsin,wsout,a,rk,
     1      u,v,ncols)
        implicit real *8 (a-h,o-z)
        save
        dimension zsin(2,1),zsout(2,1),wsin(1),wsout(1)
c 
        complex *16 rk,h0,h1,cd,a(nout,nin),u(nin,1),
     1      v(nout,1),s(10 000),w(3000 000),pot
  
c 
c        construct the matrix of interactions between the
c        squares
c 
        call prinf('entered creasvd, nin=*',nin,1)
        call prinf('entered creasvd, nout=*',nout,1)
  
  
        do 1600 i=1,nin
        do 1400 j=1,nout
  
  
        call poteval2(zsin(1,i),zsout(1,j),rk,pot)
  
        a(j,i)=pot*sqrt(wsin(i)*wsout(j))
  
 1400 continue
 1600 continue
  
cccc        call prin2('a as created*',a,nin*nout*2)
  
  
c 
c      SVD the obtained matrix
c 
        eps=1.0d-12
        lw=3000 000
        call csvdpiv(ier,a,nout,nin,v,u,s,ncols,eps,
     1      w,lw,ltot)
  
  
        call prinf('after svdpivot, ncols=*',ncols,1)
        call prin2('after svdpivot, s=*',s,ncols*2)
        call prin2('and in creasvd, rk=*',rk,2)
c 
        return
  
  
  
c 
c       compensate for the square roots of weights
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
 2600 continue
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
  
cccc        return
  
  
        d=(z1(1)-z2(1)+20)**2+(z1(2)-z2(2))**2
  
c 
  
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
        pot=pot+h0
c 
        return
        end
c 
c 
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code, and the beginning
c        of the code constructing the skeletonization of the
c        interaction operator between two "concentric" rectangles
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
c 
c 
        subroutine rects_skeleton(ier,poteval,rk,
     1       nts,coef,gap,eps,zsin_big,zsout_big,npts,
     2      w2,lenw,lused)
        implicit real *8 (a-h,o-z)
c 
        save
        real *8 zsout_big(2,1),zsin_big(2,1),w2(1)
c 
        complex *16 rk(1)
c 
c        This subroutine skeletonizes the interactions between
c        two "concentric" rectangles in R^2. Its expected use is
c        for the Helmholtz interactions, and this is the regime
c        where it has been tested.
c 
c                    Input parameters:
c 
c 
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
c  rk - the parameter (real, complex, whatever) used by the
c        subroutine poteval
c  nts - the number of nodes into which each of the layers
c        comprising the INNER square will be subdivided ON
c        EACH OF THE HORIZONTAL SIDES OF THE SQUARE. Thus, the
c        total number of nodes discretizing the inner square
C        will be 12*nts.
c  coef - the size of the inner layer of the outer square (see
c        general description above)
c  gap - the gap between the inner and outer layers of each
c        square
c  eps - the accuracy to which the SVDs will be performed
c  lenw - the length (in real *8 locations) of the work array
c        w2 supplied by the user
c 
c                Output parameters:
c 
c  ier - error return code
c  zsin_big - the skeletonized nodes on the inner square
c  zsout_big - the skeletonized nodes on the OUTER square
c  npts - the number of skeletonized nodes on each of the square
c  lused - the amount of space in array w2 (in real *8 words)
c        actually used by this subroutine
c 
c                     Work array:
c 
c  w2 - should be reasonably large
c 
c        allocate memory
c 
        ier=0
c 
        nin=nts*12
        nout=nts*24
c 
        nout2=nout/4
        nin2=nin/4
c 
        iamatr=1
        lamatr=nin2*nout2*2+10
c 
        izsin=iamatr+lamatr
        lzsin=nin*2+10
c 
        iwsin=izsin+lzsin
        lwsin=nin+10
c 
        izsout=iwsin+lwsin
        lzsout=nout*2+10
c 
        iwsout=izsout+lzsout
        lwsout=nout+10
c 
        izsin2=iwsout+lwsout
        lzsin2=lzsin/2+6
c 
        iwsin2=izsin2+lzsin2
        lwsin2=lwsin/2+6
c 
        izsout2=iwsin2+lwsin2
        lzsout2=lzsout/2+6
c 
        iwsout2=izsout2+lzsout2
        lwsout2=lwsout/2+6
c 
        izsin22=iwsout2+lwsout2
        lzsin22=lzsin/2+6
c 
        izsout22=izsin22+lzsin22
        lzsout22=lzsout/2+6
c 
        iicols=izsout22+lzsout22
        licols=nout+10
c 
        iirows=iicols+licols
        lirows=nout+10
c 
        iw=iirows+lirows
c 
        lw=nin2*nout2*2+10
c 
        lused=iw+lw
        if(lused .gt. lenw) then
            ier=16
            call prinf('bombing from rects_skeleton with ier=*',
     1            ier,1)
            return
        endif
c 
        call rects_skeleton0(ier,poteval,rk,nts,coef,gap,eps,
     1       zsin_big,zsout_big,npts,
     2       w2(iamatr),w2(izsin),
     3       w2(izsout),w2(iwsin),w2(iwsout),
     4       w2(iw),w2(izsin2),w2(izsout2),w2(iwsin2),
     5       w2(iwsout2),w2(izsin22),w2(izsout22),
     6       w2(iicols),w2(iirows) )
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine rects_skeleton0(ier,poteval,rk,
     1       nts,coef,gap,eps,zsin_big,zsout_big,npts,
     2       amatr,zsin,zsout,wsin,wsout,
     3       w,zsin2,zsout2,wsin2,wsout2,
     4      zsin22,zsout22,icols,irows)
        implicit real *8 (a-h,o-z)
c 
        save
        real *8 zsin(2,1),zsout(2,1),wsin(1),wsout(1),w(1),
     1      zsin2(2,1),wsin2(1),zsout2(2,1),wsout2(1),
     2      zsout_big(2,1),zsin_big(2,1),zsin22(1),
     3      zsout22(1),icols(1),irows(1)
c 
        complex *16 rk,amatr(1)
c 
c        construct both the inner and the outer squares
c 
        ntsin=nts
        ntsout=nts*2
c 
        call rects_cre(gap,ntsin,ntsout,coef,
     1      zsin,wsin,nin,zsout,wsout,nout,w)
c 
c        extract 1/4-th of the disctretization of
c        each rectangle
c 
cccc        call prin2('zsout=*',zsout,nout*2)
cccc        call prinf('while=*',nout,1)
cccc        call prinf('and nin=*',nin,1)
c 
        iw=23
        call zquaplot2(iw,zsin,nin,2,zsout,nout,2,
     1      'both rectangles*')
c 
        ifouter=0
        call rect_extract(zsin,wsin,zsin2,wsin2,ntsin,ifouter)
c 
        ifouter=0
        call rect_extract(zsout,wsout,zsout2,wsout2,
     1      ntsout,ifouter)
c 
        nin2=nin/4
        nout2=nout/4
c 
        iw=24
        call zquaplot2(iw,zsin2,nin2,2,zsout2,nout2,2,
     1      'both rectangles after extraction*')
c 
c        construct the symmetrized svd
c 
        call rects_skel(ier,zsin2,wsin2,nin2,poteval,
     1      zsout2,wsout2,nout2,amatr,rk,eps,ncols,
     2      w,zsout_big,zsin_big,zsin22,zsout22,zsin,nin,
     3      zsout,nout,icols,irows)
c 
        npts=ncols*4
  
  
        iw=26
        call zquaplot2(iw,zsin_big,npts,2,zsout_big,npts,2,
     1      'both final rectangles*')
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine rects_skel(ier,zsin,wsin,nin,poteval,
     1      zsout,wsout,nout,a,rk,eps,ncols,
     2      w,zsout_big,zsin_big,zsin22,zsout22,
     3      zsin_big_init,nin_big_init,
     4      zsout_big_init,nout_big_init,icols,irows)
c 
        implicit real *8 (a-h,o-z)
        save
        dimension zsin(2,1),zsout(2,1),wsin(1),wsout(1),
     1      z1(3),z2(3),zsin22(2,1),zsout_big(2,1),
     2      zsin_big(2,1),zsout22(2,1),icols(1),irows(1),
     3      zsin_big_init(1)
c 
        complex *16 rk,w(2),
     1      a(nin,nout),cd00,cd10,cd01,cd11
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
        a(i,j)=cd00+cd10+cd01+cd11
        a(i,j)=a(i,j)*sqrt(wsin(i)*wsout(j))
c 
 1400 continue
 1600 continue
c 
        call cskel_basic(a,nin,nout,eps,irows,icols,ncols,w)
  
cccc        call prinf('after cskel_basic, icols=*',icols,ncols)
cccc        call prinf('after cskel_basic, irows=*',irows,ncols)
  
        do 2200 i=1,ncols
c 
        j=irows(i)
        zsin22(1,i)=zsin(1,j)
        zsin22(2,i)=zsin(2,j)
c 
c 
        j=icols(i)
        zsout22(1,i)=zsout(1,j)
        zsout22(2,i)=zsout(2,j)
c 
 2200 continue
  
        call rect_expand(zsin22,zsin_big,ncols)
        call rect_expand(zsout22,zsout_big,ncols)
c 
        npts=ncols*4
cccc        call square_final_reord(zs,ws,n,zsbig,nbig,inds)
        call square_final_reord(zsin_big,w,npts,
     1      zsin_big_init,nin_big_init,icols)
c 
        call square_final_reord(zsout_big,w,npts,
     1      zsout_big_init,nout_big_init,icols)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine rect_expand(zs,zs2,n)
        implicit real *8 (a-h,o-z)
        save
        dimension zs(2,1),zs2(2,1)
c 
c 
        do 1200 i=1,n
c 
        zs2(1,i)=zs(1,i)
        zs2(2,i)=zs(2,i)
 1200 continue
c 
c        implement the x-reflection
c 
        do 1300 i=1,n
c 
        zs2(1,i+n)=-zs2(1,i)
        zs2(2,i+n)=zs2(2,i)
 1300 continue
c 
c        implement the y-reflection
c 
         do 1400 i=1,n*2
c 
         zs2(1,i+n*2)=zs2(1,i)
         zs2(2,i+n*2)=-zs2(2,i)
 1400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine rect_extract(zs,ws,zs3,ws3,nts,ifouter)
        implicit real *8 (a-h,o-z)
        save
        dimension zs(2,1),ws(1),zs3(2,1),ws3(1)
c 
c        extract from the big arrays the geometry pertaining
c        to the upper and left sides of the inner rectangle
c 
        do 1600 i=1,nts*3/2
c 
        zs3(1,i)=zs(1,i+nts/2)
        zs3(2,i)=zs(2,i+nts/2)
        ws3(i)=ws(i+nts/2)
c 
cccc        if(ifouter .eq. 1) goto 1600
  
        zs3(1,i+nts*3/2)=zs(1,i+nts/2+nts*6)
        zs3(2,i+nts*3/2)=zs(2,i+nts/2+nts*6)
        ws3(i+nts*3/2)=ws(i+nts/2+nts*6)
c 
 1600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine rects_cre(gap,ntsin,ntsout,coef,
     1      zsin,wsin,nin,zsout,wsout,nout,w)
        implicit real *8 (a-h,o-z)
c 
        save
        real *8 zsin(2,1),wsin(1),zsout(2,1),wsout(1),w(1)
c 
c        construct the inner rectangle
c 
        its=1
        lts=ntsout+2
c 
        iwhts=its+lts
        lwhts=ntsout*2
c 
        call rectin_cre(gap,ntsin,zsin,wsin,nin,w(its),w(iwhts) )
c 
c        construct the outer rectangle
c 
        call rectout_cre(coef,gap,ntsout,zsout,wsout,nout,
     1      w(its),w(iwhts) )
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine rectout_cre(coef,gap,nts,
     1      zs,weights,nn,ts,whts)
        implicit real *8 (a-h,o-z)
c 
        save
        real *8 zs(2,1),ts(10 000),
     1      whts(10 000),weights(1)
c 
c        construct the rectangle
c 
        itype=1
        call legeexps(itype,nts,ts,u,v,whts)
c 
        x1=-coef
        x2=coef
c 
        y1=x1-1
        y2=x2+1
c 
        y0=0
c 
        call square_lineconstr(x2,y2,x1,y2,ts,nts,zs)
        call square_lineconstr(x1,y2,x1,y0,ts,nts,zs(1,nts+1))
        call square_lineconstr(x1,y0,x1,y1,ts,nts,zs(1,2*nts+1))
        call square_lineconstr(x1,y1,x2,y1,ts,nts,zs(1,3*nts+1))
        call square_lineconstr(x2,y1,x2,y0,ts,nts,zs(1,4*nts+1))
        call square_lineconstr(x2,y0,x2,y2,ts,nts,zs(1,5*nts+1))
c 
ccc        x1=x1+gap
ccc        x2=x2-gap
c 
ccc        y1=y1+gap
ccc        y2=y2-gap
c 
c 
        x1=x1-gap
        x2=x2+gap
c 
        y1=y1-gap
        y2=y2+gap
c 
c 
        call square_lineconstr(x2,y2,x1,y2,ts,nts,zs(1,6*nts+1))
        call square_lineconstr(x1,y2,x1,y0,ts,nts,zs(1,7*nts+1))
        call square_lineconstr(x1,y0,x1,y1,ts,nts,zs(1,8*nts+1))
        call square_lineconstr(x1,y1,x2,y1,ts,nts,zs(1,9*nts+1))
        call square_lineconstr(x2,y1,x2,y0,ts,nts,zs(1,10*nts+1))
        call square_lineconstr(x2,y0,x2,y2,ts,nts,zs(1,11*nts+1))
c 
        nn=nts*12
c 
c        construct the array of weights
c 
        do 1600 i=1,nts
c 
        do 1500 j=0,11
        weights(i+j*nts)=whts(i)
c 
 1500 continue
c 
 1600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine rectin_cre(gap,nts,
     1      zs,weights,nn,ts,whts)
        implicit real *8 (a-h,o-z)
c 
        save
        real *8 zs(2,1),ts(10 000),
     1      whts(10 000),weights(1)
c 
c        construct the rectangle
c 
        itype=1
        call legeexps(itype,nts,ts,u,v,whts)
c 
        x1=-1
        x2=1
c 
        y1=-2
        y2=2
c 
        y0=0
c 
        call square_lineconstr(x2,y2,x1,y2,ts,nts,zs)
        call square_lineconstr(x1,y2,x1,y0,ts,nts,zs(1,nts+1))
        call square_lineconstr(x1,y0,x1,y1,ts,nts,zs(1,2*nts+1))
        call square_lineconstr(x1,y1,x2,y1,ts,nts,zs(1,3*nts+1))
        call square_lineconstr(x2,y1,x2,y0,ts,nts,zs(1,4*nts+1))
        call square_lineconstr(x2,y0,x2,y2,ts,nts,zs(1,5*nts+1))
c 
        x1=x1+gap
        x2=x2-gap
c 
        y1=y1+gap
        y2=y2-gap
c 
c 
        call square_lineconstr(x2,y2,x1,y2,ts,nts,zs(1,6*nts+1))
        call square_lineconstr(x1,y2,x1,y0,ts,nts,zs(1,7*nts+1))
        call square_lineconstr(x1,y0,x1,y1,ts,nts,zs(1,8*nts+1))
        call square_lineconstr(x1,y1,x2,y1,ts,nts,zs(1,9*nts+1))
        call square_lineconstr(x2,y1,x2,y0,ts,nts,zs(1,10*nts+1))
        call square_lineconstr(x2,y0,x2,y2,ts,nts,zs(1,11*nts+1))
c 
        nn=nts*12
c 
c        construct the array of weights
c 
        do 1600 i=1,nts
c 
        do 1500 j=0,11
        weights(i+j*nts)=whts(i)
c 
 1500 continue
c 
 1600 continue
c 
        return
        end
