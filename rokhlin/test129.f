        implicit real *8 (a-h,o-z)
        dimension polsout(1000),dersx(1000),dersy(1000),
     1      nsnonz2(10 000),amatr(1 000 000),xsout(10 000),
     2      wsout(10 000),ysout(10 000),sums(10 000)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER mmax'
        READ *,mmax
        CALL PRINf('mmax=*',mmax,1 )
c 
            iw=41
        ir=18
cccc        ir=151
        call symevafin(ir,mmax,nnout,nsnonz2)
  
  
        call prinf('after symevafin, nnout=*',nnout,1)
        call prinf('after symevafin, nsnonz2=*',nsnonz2,nnout)
  
c 
c       test the evaluation of the symmetrized polynomials
c 
        x=0.1
        y=-0.3
  
        itype=1
        call symeval(itype,mmax,x,y,polsout,dersx,dersy)
  
        call prin2('after symeval, polsout=*',polsout,nnout)
        call prin2('after symeval, dersx=*',dersx,nnout)
        call prin2('after symeval, dersy=*',dersy,nnout)
  
        n=60
        call quadinit(n,mmax,nnout,amatr,xsout,ysout,wsout)
c 
c       test the construction of the matrix and of the RHS
c 
cccc        call quadmatr(xsout,ysout,wsout,nnout,
cccc     1      nnout,amatr,sums,itype,
cccc     1      mmax)
  
  
        imin=5
cccc        call quapoint(ier,xsout,ysout,wsout,nnout,nnout,
cccc     1      amatr,mmax,imin)
  
        call quapoints(ier,iw,xsout,ysout,wsout,nnout,nnout,
     1      amatr,mmax)
  
  
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine quapoints(ier,iw,xs,ys,ws,npts0,nfuns,
     1      amatr,mmax)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),ys(1),ws(1),amatr(1),xs3(10 000),
     1      ys3(10 000),ws3(10 000),signifs(10 000)
c 
        call prinf('in quapoint, nfuns=*',nfuns,1)
c 
c       one point after another, eliminate the nodes
c 
        imin=1
        npts=npts0
c 
c       store on disk the header for this run
c 
 1200 format('        In this run, mmax= ',i6,
     1       '; nfuns= ',i6)
 1400 format(10x)
  
        write(iw,1400)
        write(iw,1200) mmax,nfuns
  
  
  
  
c 
        do 4000 i=1,10000
c 
c        attempt to eliminate this point
c 
        call quapoint(ier,xs,ys,ws,npts,nfuns,
     1      amatr,mmax,imin,xs3,ys3,ws3,istage)
  
        call prinf('after quapoint, ier=*',ier,1)
  
c 
c       if the attempt has been successful - start the elimination of
c       the next point
c 
  
        if(ier .ne. 0) goto 3000
c 
c       the search for this npts has been successful. Start the
c       search for the next smaller npts
c 
        istage=1
c 
        call mvecmove(xs3,xs,npts)
        call mvecmove(ys3,ys,npts)
        call mvecmove(ws3,ws,npts)
  
  
  
        call pntsignif(xs,ys,ws,npts-1,
     1      nfuns,signifs,mmax)
  
        ifanyout=0
        do 2200 jj=1,npts-1
c 
cccc        if(ws(jj) .lt. 0) signifs(jj)=-1
        if(ws(jj) .lt. 0) signifs(jj)=ws(jj)
  
        call ifinside(xs(jj),ys(jj),ifin)
        if(ifin .eq. -1) signifs(jj)=-10
        if(ifin .eq. -1) ifanyout=ifanyout+1
 2200 continue
c 
        call quabubb3(signifs,xs,ys,ws,npts-1)
c 
        call prin2('and signifs after reordering*',signifs,npts-1)
c 
        call prin2('affter quabubb3, ws=*',ws,npts-1)
  
  
            call prin2('success, xs=*',xs,npts-1)
            call prin2('success, ys=*',ys,npts-1)
            call prin2('success, ws=*',ws,npts-1)
  
        npts=npts-1
  
        call prinf('writing on disk data for npts=*',npts,1)
        call prinf('ifanyout=*',ifanyout,1)
  
 2400 format(10x)
c 
 2600 format('    in the following data, some points are ',
     1    'outside the triangle')
  
 2800 format('    in the following data, all points are ',
     1    'inside the triangle')
  
        write(iw,2400)
        if(ifanyout .eq. 0) write(iw,2800)
        if(ifanyout .eq. 1) write(iw,2600)
  
  
        call quaestor(iw,xs,ys,ws,npts)
  
        call prinf('after change, npts=*',npts,1)
        call prinf('and nfuns=*',nfuns,1)
  
        nfuns7=nfuns/3
        jjjjj=nfuns-nfuns*3
        if(jjjjj .ne. 0) nfuns7=nfuns7+1
  
  
cccc        if(npts-1 .lt. nfuns/3) stop
cccccc        if(npts-1 .lt. nfuns7) stop
        if(npts .lt. nfuns7) stop
c 
            imin=1
            goto 4000
  
  
  
  
  
 3000 continue
c 
c       if the attempt has not been successful - increment imin and restart
        if(ier .eq. 0) goto 3600
c 
            imin=imin+1
  
            if(imin .gt. npts) call prinf('failed at imin=*',imin,1)
            if( (imin .gt. npts) .and. (istage .eq. 2) ) stop
  
            if( (imin .gt. npts) .and. (istage .eq. 1) ) then
                imin=1
                istage=2
                goto 4000
            endif
 3600 continue
  
 4000 continue
  
        return
        end
  
  
  
c 
c 
c 
c 
c 
        subroutine pntsignif(xs,ys,ws,npts,nfuns,signifs,mmax)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),ws(1),sums(1),ys(1),
     1      polsout(10 000),signifs(1)
c 
c        evaluate the magnitudes of the contributions of all nodes
c 
        do 1600 j=1,npts
c 
        signifs(j)=0
        itype=0
        call symeval(itype,mmax,xs(j),ys(j),polsout,dersx,dersy)
c 
        do 1400 i=1,nfuns
c 
        signifs(j)=signifs(j)+polsout(i)**2
 1400 continue
c 
        signifs(j)=signifs(j)*ws(j)**2
 1600 continue
  
        call prin2('signifs as evaluated=*',signifs,npts)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine quabubb2(xbase,ybase,xs,ys,ws,npts)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),ys(1),ws(1)
c 
        do 1400 i=1,npts
        do 1200 j=1,npts-1
c 
        dj=(xs(j)-xbase)**2+(ys(j)-ybase)**2
        djp1=(xs(j+1)-xbase)**2+(ys(j+1)-ybase)**2
  
  
  
        if(dj .le. djp1) goto 1200
c 
        d=xs(j)
        xs(j)=xs(j+1)
        xs(j+1)=d
c 
        d=ys(j)
        ys(j)=ys(j+1)
        ys(j+1)=d
c 
        d=ws(j)
        ws(j)=ws(j+1)
        ws(j+1)=d
c 
 1200 continue
 1400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine quapoint(ier,xs,ys,ws,npts0,nfuns,
     1      amatr,mmax,imin,xs3,ys3,ws3,istage)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),ys(1),ws(1),amatr(1),xs3(10 000),
     1      ys3(10 000),ws3(10 000),ws4(10 000)
c 
        call prinf('in quapoint, nfuns=*',nfuns,1)
c 
c       omit the node number imin
c 
        j=0
        do 1200 i=1,npts0
c 
        if(i .eq. imin) goto 1200
c 
        j=j+1
        xs3(j)=xs(i)
        ys3(j)=ys(i)
        ws3(j)=ws(i)
 1200 continue
c 
        epsout=1.0d-20
cccc        do 2400 i=1,60
        do 2400 i=1,60
  
        call prinf('in quapoint, i=*',i,1)
        call prinf('in quapoint, istage=*',istage,1)
  
        npts=npts0-1
        call quaenewt(jer,xs3,ys3,ws3,npts,nfuns,
     1      amatr,mmax,dnew,kcontr)
  
        call prinf('after quaenewt, jer=*',jer,1)
        call prinf('after quaenewt, kcontr=*',kcontr,1)
  
  
        if(jer .ne. 0) goto 2500
  
        if( (i .gt. 4) .and. (dnew/dold .gt. 0.9999) ) goto 2500
       if( (i .gt. 30) .and. (dnew/dold .gt. 0.99d0) ) goto 2500
        if( (i .gt. 4) .and. (kcontr .gt. 3) ) goto 2500
  
  
        if(istage .ne. 1) goto 1600
  
        if( (i .gt. 5) .and. (dnew/dold .gt. 0.3) ) goto 2500
  
        if( (i .gt. 4) .and. (kcontr .gt. 3) ) goto 2500
  
 1600 continue
  
        dold=dnew
  
  
        call prini(0,0)
  
cccc        call quadwhts(ier,xs3,ys3,ws4,npts,nfuns,
        call quadwhts(ier2,xs3,ys3,ws3,npts,nfuns,
     1      amatr,mmax)
  
        call prini(6,13)
  
  
cccc        call prin2('after quadwhts, ws4=*',ws4,npts)
cccc        call prin2('while, ws3=*',ws3,npts)
  
  
c 
        if(jer .ne. 0) then
            ier=8
            return
        endif
c 
        if(dnew .lt. epsout) goto 2600
c 
 2400 continue
c 
 2500 continue
c 
        ier=4
        return
c 
 2600 continue
c 
        ier=0
  
  
  
        return
        end
  
  
c 
c 
c 
c 
c 
        subroutine quadwhts(ier,xs,ys,ws,npts,nfuns,
     1      amatr,mmax)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),ys(1),ws(1),amatr(nfuns,npts),
     1      u(100 000),v(100 000),w(100 000),tt(100 000),
     2      rints(1000),polsout(1000),rnorms(1000)
c 
c        construct the least squares problem defining the weights
c 
        itype=0
        do 1400 i=1,npts
c 
        call symeval(itype,mmax,xs(i),ys(i),polsout,dersx,dersy)
c 
        do 1200 j=1,nfuns
c 
        amatr(j,i)=polsout(j)
 1200 continue
 1400 continue
c 
c       construct the qr-decomposition of the matrix
c 
        eps=1.0d-14
        call leastsq(amatr,u,w,tt,nfuns,npts,ncols,rnorms,eps,v)
c 
        call prin2('in quadwhts, rnorms=*',rnorms,ncols)
c 
c       construct the least squares solution of the Newton equation
c 
        do 1600 i=1,nfuns
c 
        rints(i)=0
 1600 continue
c 
        rints(1)=1
c 
        call leastsq2(u,w,tt,nfuns,npts,ncols,rints,ws,v)
  
        call prin2('after leastsq2, ws=*',ws,ncols)
c 
        return
        end
  
  
  
c 
c 
c 
c 
c 
        subroutine quabubb3(signifs,xs,ys,ws,npts)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),ys(1),ws(1),signifs(1)
c 
        do 1400 i=1,npts
        do 1200 j=1,npts-1
c 
        if(signifs(j) .le. signifs(j+1) ) goto 1200
c 
        d=xs(j)
        xs(j)=xs(j+1)
        xs(j+1)=d
c 
        d=ys(j)
        ys(j)=ys(j+1)
        ys(j+1)=d
c 
        d=ws(j)
        ws(j)=ws(j+1)
        ws(j+1)=d
c 
        d=signifs(j)
        signifs(j)=signifs(j+1)
        signifs(j+1)=d
c 
 1200 continue
 1400 continue
c 
        return
        end
c 
  
c 
c 
c 
c 
c 
        subroutine quabubb(xs,ys,ws,npts)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),ys(1),ws(1)
c 
        do 1400 i=1,npts
        do 1200 j=1,npts-1
c 
        if(ws(j) .le. ws(j+1) ) goto 1200
c 
        d=xs(j)
        xs(j)=xs(j+1)
        xs(j+1)=d
c 
        d=ys(j)
        ys(j)=ys(j+1)
        ys(j+1)=d
c 
        d=ws(j)
        ws(j)=ws(j+1)
        ws(j+1)=d
c 
 1200 continue
 1400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine quaenewt(ier,xs,ys,ws,npts,nfuns,
     1      amatr,mmax,dnew,kcontr)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),ys(1),ws(1),sums(1000),amatr(nfuns,npts*3),
     1      sol(1000),rints(1000),rhs(1000),xs2(1000),ws2(1000),
     2      ys2(1000)
c 
        dimension w(1000 000),u(1000 000),v(1000 000),tt(1000 000),
     1      rnorms(1000)
  
cccc        call prin2('entered quaenewt, xs=*',xs,npts)
cccc        call prin2('entered quaenewt, ys=*',ys,npts)
cccc        call prin2('entered quaenewt, ws=*',ws,npts)
  
c 
        do 1100 i=1,nfuns
c 
        rints(i)=0
 1100 continue
c 
        rints(1)=1
c 
        itype=1
        call quadmatr(xs,ys,ws,npts,nfuns,
     1      amatr,sums,itype,
     1      mmax)
c 
        dold=0
        do 1200 i=1,nfuns
c 
        rhs(i)=sums(i)-rints(i)
        dold=dold+rhs(i)**2
 1200 continue
c 
c       construct the Gram-schmidt decomposition of the matrix a
c 
        eps=1.0d-12
  
        np3=npts*3
        call leastsq(amatr,u,w,tt,nfuns,np3,ncols,rnorms,eps,v)
c 
c       construct the least squares solution of the Newton equation
c 
        call leastsq2(u,w,tt,nfuns,np3,ncols,rhs,sol,v)
  
cccc        call prin2('after leastsq2, sol=*',sol,ncols)
  
c 
c       perform step-length control
c 
        scale=1
        icount=0
        do 4000 k=1,5
cccc        do 4000 k=1,10
c 
        kcontr=k
  
        call prinf('in quaenewt, k=*',k,1)
c 
c       . . . recompute the right-hand sides
c 
        ifout=0
  
        ifrec=0
        do 2400 i=1,npts
c 
        xs2(i)=xs(i)-sol(i)         *scale
        ys2(i)=ys(i)-sol(npts+i)    *scale
        ws2(i)=ws(i)-sol(2*npts+i)    *scale
  
 2400 continue
c 
        ifrec=1
        call prini(0,0)
        if(ifrec .eq. 1)
     1      call quadwhts(ier,xs2,ys2,ws2,npts,nfuns,
     2      amatr,mmax)
  
  
        call prini(6,13)
  
  
        itype=0
        call quadmatr(xs2,ys2,ws2,npts,nfuns,
     1      amatr,sums,itype,
     1      mmax)
c 
        dnew=0
        do 2600 i=1,nfuns
c 
        dnew=dnew+(sums(i)-rints(i))**2
 2600 continue
  
        call prin2('dnew=*',dnew,1)
        call prin2('and dold=*',dold,1)
        call prinf('while nfuns=*',nfuns,1)
        call prinf('and npts=*',npts,1)
  
        if( (dold .gt. dnew) .and. (ifout .eq. 0) )  goto 4200
  
        scale=scale/4
  
  
 4000 continue
  
        ier=4
        return
 4200 continue
  
       ier=0
  
  
cccc        call prin2('and xs2=*',xs2,npts)
cccc        call prin2('and ys2=*',ys2,npts)
cccc        call prin2('and ws2=*',ws2,npts)
  
c 
        call mvecmove(xs2,xs,npts)
        call mvecmove(ys2,ys,npts)
        call mvecmove(ws2,ws,npts)
c 
        return
        end
  
  
c 
c 
c 
c 
c 
        subroutine ifinside(x,y,ifin)
        implicit real *8 (a-h,o-z)
c 
c 
        save
        d12=u12*x+v12-y
  
cccc        call prin2('d12=*',d12,1)
  
        d23=u23*x+v23-y
  
cccc        call prin2('d23=*',d23,1)
  
        d31=u31*x+v31-y
  
cccc        call prin2('d31=*',d31,1)
  
  
  
        ifin=1
        if(d12 .lt. 0) ifin=-1
        if(d23 .gt. 0) ifin=-1
        if(d31 .lt. 0) ifin=-1
  
c 
        return
c 
c 
c 
c 
        entry ifinsin(x1,y1,x2,y2,x3,y3)
c 
        u12=(y2-y1)/(x2-x1)
        v12=y1-u12*x1
c 
        u23=(y3-y2)/(x3-x2)
        v23=y2-u23*x2
c 
        u31=(y1-y3)/(x1-x3)
        v31=y3-u31*x3
c 
  
        return
        end
c 
c 
c 
c 
c 
        subroutine quadmatr(xs,ys,ws,npts,nfuns,a,sums,itype,
     1      mmax)
        implicit real *8 (a-h,o-z)
        save
        dimension a(nfuns,npts*3),xs(1),ws(1),sums(1),ys(1),
     1      polsout(10 000),dersx(10 000),dersy(10 000)
c 
c        construct the matrix of the Jacobian of the mapping
c        connection the points xs and the weights ws with the
c        integrals of the functions
c 
        do 1200 i=1,nfuns
c 
        sums(i)=0
 1200 continue
c 
        do 1600 j=1,npts
c 
        call symeval(itype,mmax,xs(j),ys(j),polsout,dersx,dersy)
c 
        do 1400 i=1,nfuns
c 
        if(itype .eq. 0) goto 1300
c 
        a(i,j)=dersx(i)*ws(j)
        a(i,j+npts)=dersy(i)*ws(j)
        a(i,j+npts*2)=polsout(i)
c 
 1300 continue
        sums(i)=sums(i)+polsout(i)*ws(j)
 1400 continue
c 
 1600 continue
  
cccc        call prin2('and sums=*',sums,nfuns)
  
       return
       end
  
c 
c 
c 
c 
c 
        subroutine quadinit(n,mmax,nn2,amatr,xsout,ysout,wsout)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(10 000),whts(10 000),zs(2,10 000),ys(10 000)
        character *1 mes(100)
c 
        dimension xstria(5),ystria(5),
     1      vert1(2),vert2(2),vert3(2),
     2      vert21(4),vert22(4),vert23(4),w(100 000),rnorms(100 000),
     3      ipivots(100 000),xsout(10 000),wsout(10 000),ysout(10 000)
  
c 
        dimension polsout(100 000),amatr(nn2,n*n)
c 
c        construct the "tensor product" discretization of the triangle
c 
        done=1
        xstria(1)=0
        xstria(2)=-1
        xstria(3)=1
        xstria(4)=0
c 
        ystria(1)=2/sqrt(done*3)
        ystria(2)=-1/sqrt(done*3)
        ystria(3)=-1/sqrt(done*3)
        ystria(4)=2/sqrt(done*3)
c 
c       construct a small number of tensor product nodes
c 
        vert1(1)=xstria(1)
        vert2(1)=xstria(2)
        vert3(1)=xstria(3)
c 
        vert1(2)=ystria(1)
        vert2(2)=ystria(2)
        vert3(2)=ystria(3)
c 
c       construct the 1/3-rd of the triangle
c 
        vert21(1)=0
        vert21(2)=0
c 
        vert22(1)=-1
        vert22(2)=-1/sqrt(done*3)
c 
        vert23(1)=1
        vert23(2)=-1/sqrt(done*3)
c 
  
  
  
cccc        call ifinsin(vert21(1),vert21(2),vert22(1),vert22(2),
cccc     1      vert23(1),vert23(2) )
  
        call ifinsin(vert1(1),vert1(2),vert2(1),vert2(2),
     1      vert3(1),vert3(2) )
  
  
  
  
        ifinit=1
        call triagauc(n,vert21,vert22,vert23,zs,
     1      whts,ifinit,w)
c 
        kk=n**2
  
        do 2400 i=1,kk
c 
        xs(i)=zs(1,i)
        ys(i)=zs(2,i)
 2400 continue
c 
        call msgmerge('lower 1/3 of standard triangle*',
     1     ' with large number of tensor product nodes*',mes)
c 
        iw=22
        call lotaplot2(iw,xstria,ystria,4,xs,ys,kk,mes)
c 
c       construct the matrix ov values of orthogonal polynomials
c       at the nodes of the discretization
c 
  
  
cccc        call prin2('before the loop, xs=*',xs,kk)
cccc        call prin2('before the loop, ys=*',ys,kk)
  
  
cccc            stop
  
        itype=0
        do 3400 i=1,kk
  
  
c 
        call symeval(itype,mmax,xs(i),ys(i),polsout,dersx,dersy)
  
  
cccc        call prin2('and polsout=*',polsout,nn2)
  
c 
        do 3200 j=1,nn2
c 
        amatr(j,i)=polsout(j)*sqrt(whts(i))
cccc        amatr(j,i)=polsout(j)
 3200 continue
 3400 continue
  
  
cccc        call prin2('amatr as constructed*',amatr,nn2*n*n)
  
c 
c        use the pivoted Gram-Schmidt to obtain the initial quadrature
c        nodes
c 
        eps=1.0d-10
        call allsppiv(amatr,nn2,kk,rnorms,eps,ncols,ipivots)
  
        call prin2('after allsppiv, rnorms=*',rnorms,ncols)
        call prinf('after allsppiv, ipivots=*',ipivots,ncols)
  
        call prinf('after allsppiv, ncols=*',ncols,1)
        call prinf('while nn2=*',nn2,1)
        call prinf('and kk=*',kk,1)
c 
c       extract the obtained nodes
c 
        do 3600 i=1,ncols
c 
        j=ipivots(i)
  
        xsout(i)=xs(j)
        ysout(i)=ys(j)
 3600 continue
c 
c       obtain the weights
c 
        itype=0
        do 4400 i=1,ncols
c 
        call symeval(itype,mmax,xsout(i),ysout(i),
     1      polsout,dersx,dersy)
c 
        do 4200 j=1,ncols
c 
cccc        amatr(i,j)=polsout(j)
        amatr(j,i)=polsout(j)
 4200 continue
 4400 continue
c 
        do 4600 i=1,nn2
        wsout(i)=0
 4600 continue
c 
        wsout(1)=1
  
        call qrsolv(amatr,ncols,wsout,rcond)
  
        call prin2('after qrsolv, wsout=*',wsout,ncols)
  
        call prin2('after qrsolv, rcond=*',rcond,1)
  
  
        d=0
        do 4800 i=1,ncols
c 
        d=d+wsout(i)
 4800 continue
  
  
  
        call prin2('and the sum of weights is*',d,1)
  
        dd=d**4-3
        call prin2('and dd=*',dd,1)
  
  
  
  
  
        return
        end
c 
c 
c 
c 
c 
        subroutine symevafin(ir,mmax,nn2out,nsnonz27)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(10 000),whts(10 000),zs(2,10 000),ys(10 000)
        character *1 mes(100)
c 
        dimension xstria(5),ystria(5),
     1      vert1(2),vert2(2),vert3(2),nsnonz27(1),
     2      vert21(4),vert22(4),vert23(4),nsnonz2(100 000)
c 
        dimension pols(1000 000),pols2(500 000),arr(2 000 000)
c 
        dimension z1(2),polsout(1),dersx(1),dersy(1),
     1      polsout1(10 000),dersx1(10 000),dersy1(10 000)
c 
c 
c       construct the standard triangle
c 
        done=1
        xstria(1)=0
        xstria(2)=-1
        xstria(3)=1
        xstria(4)=0
c 
        ystria(1)=2/sqrt(done*3)
        ystria(2)=-1/sqrt(done*3)
        ystria(3)=-1/sqrt(done*3)
        ystria(4)=2/sqrt(done*3)
c 
c       construct a small number of tensor product nodes
c 
        n=4
c 
        vert1(1)=xstria(1)
        vert2(1)=xstria(2)
        vert3(1)=xstria(3)
c 
        vert1(2)=ystria(1)
        vert2(2)=ystria(2)
        vert3(2)=ystria(3)
c 
c       construct the 1/3-rd of the triangle
c 
        vert21(1)=0
        vert21(2)=0
c 
        vert22(1)=-1
        vert22(2)=-1/sqrt(done*3)
c 
        vert23(1)=1
        vert23(2)=-1/sqrt(done*3)
c 
        ifinit=1
        call triagauc(n,vert21,vert22,vert23,zs,
     1      whts,ifinit,w)
c 
        kk=n**2
  
        do 2400 i=1,kk
c 
        xs(i)=zs(1,i)
        ys(i)=zs(2,i)
 2400 continue
c 
c       rotate the nodes
c 
        do 2600 i=1,kk
c 
        call rotate(xs(i),ys(i),xs(i+kk),ys(i+kk))
c 
        whts(kk+i)=whts(i)
 2600 continue
c 
        do 2800 i=1,kk
c 
        call rotate(xs(i+kk),ys(i+kk),xs(i+kk*2),ys(i+kk*2))
c 
        whts(2*kk+i)=whts(i)
 2800 continue
  
        call msgmerge('standard triangle with small number*',
     1     ' of tensor product nodes*',mes)
c 
        iw=21
        call lotaplot2(iw,xstria,ystria,4,xs,ys,kk*3,mes)
c 
c        evaluate all orthogonal polynomials at the nodes of the
c        sparse discretization, and use the btained values to
c        determine which of the polynomials are not killed by
c        the symmetrization
c 
        do 3200 i=1,kk*3
c 
        zs(1,i)=xs(i)
        zs(2,i)=ys(i)
 3200 continue
c 
        call polssym(zs,whts,kk,mmax,ir,pols,
     1      pols2,nsnonz2,nn2,arr)
  
        call prinf('after polssym, nn2=*',nn2,1)
        nn2out=nn2
        mmax7=mmax
        npols=(mmax+1)*(mmax+2)/2
c 
        do 3400 i=1,nn2
c 
        nsnonz27(i)=nsnonz2(i)
 3400 continue
  
  
  
        call prinf('exiting symevafin, nn2out=*',nn2out,1)
c 
        return
c 
c 
c 
c 
        entry symeval(itype,maxord,x,y,polsout,dersx,dersy)
c 
c        evaluate all polynomials at all three points (the original one,
c        plus two rotated images)
c 
        z1(1)=x
        z1(2)=y
c 
        if(itype .eq. 1)
     1    call orthoeva3(maxord,z1,polsout1,dersx1,dersy1,arr,pols)
c 
        if(itype .eq. 0)
     1    call orthoeva(maxord,z1,polsout1,arr,pols)
c 
        if(itype .eq. 0) goto 4400
c 
        do 4200 i=1,nn2
c 
        j=nsnonz2(i)
        polsout(i)=polsout1(j)
        dersx(i)=dersx1(j)
        dersy(i)=dersy1(j)
c 
 4200 continue
c 
        return
c 
 4400 continue
c 
        do 4600 i=1,nn2
c 
        j=nsnonz2(i)
        polsout(i)=polsout1(j)
 4600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine polssym(zs,whts,kk,mmax,ir,pols,
     1      pols2,nsnonz2,nn2,arr)
        implicit real *8 (a-h,o-z)
        save
        dimension zs(2,1),whts(1),pols(kk*3,1),
     1      arr(1),s(200 000),
     2      pols2(kk,1),w(300 000),pp(600 000),
     3      nsnonz2(100 000),dd2s(100 000)
c 
c 
c       initialize the evaluator of orthogonal polynomials
c 
        call ortread(ier,ir,mmaxr,mmm,arr,ns,s)
  
c 
       call prinf('after ortretr, mmaxr=*',mmaxr,1)
       call prinf('after ortretr, mmm=*',mmm,1)
       call prinf('after ortretr, ns=*',ns,1)
       call prinf('after ortretr, mmax=*',mmax,1)
c 
c       evaluate all orthogonal polynomials at all nodes
c 
        npols=(mmax+1)*(mmax+2)/2
  
        do 2400 i=1,kk*3
  
cccc        call prinf(
c 
        call orthoeva(mmax,zs(1,i),pp,arr,w)
c 
        do 2200 j=1,npols
c 
        pols(i,j)=pp(j)
 2200 continue
 2400 continue
c 
c       symmetrize all polynomials with respect to the rotation
c 
        do 2800 i=1,npols
c 
        do 2600 j=1,kk
c 
        pols2(j,i)=pols(j,i)+pols(j+kk,i)+pols(j+kk*2,i)
        pols2(j,i)=pols2(j,i)/3
 2600 continue
c 
 2800 continue
c 
c       now, construct the list of symmetrized functions that are not zero
c 
        nn2=0
        do 3400 i=1,npols
c 
        dd2=0
        do 3200 j=1,kk
c 
        dd2=dd2+pols2(j,i)**2 *whts(j)
 3200 continue
c 
        dd2s(i)=dd2
        if(dd2 .gt. 1.0d-10) then
            nn2=nn2+1
            nsnonz2(nn2)=i
        endif
c 
 3400 continue
  
        call prinf('nn2 as evaluated*',nn2,1)
        call prinf('while npols*',npols,1)
cccc        call prin2('and dd2s=*',dd2s,npols)
        call prinf('and nsnonz2=*',nsnonz2,nn2)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine rotate(xin,yin,xout,yout)
        implicit real *8 (a-h,o-z)
        save
        data ifcalled/0/
c 
        if(ifcalled .eq. 1) goto 2000
c 
c       initialize the matrix of rotation
c 
        done=1
        pi=atan(done)*4
        theta=2*pi/3
        a11=cos(theta)
        a22=cos(theta)
        a12=-sin(theta)
        a21=-a12
        ifcalled=1
 2000 continue
c 
c       apply the rotation matrix to the input vector
c 
        xout=a11*xin+a12*yin
        yout=a21*xin+a22*yin
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine quaestin(iw,ncols,a,b,rlams)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),ws(1),rlams(1),ys(1)
c 
c       store the header of this run
c 
 1200 format('c        ')
 1800 format('c     In this run, the singular values are:')
        write(iw,1200)
        write(iw,1800)
        write(iw,1200)
c 
 1900 format('c',6(2x,e11.5))
c 
        n=ncols/6
        ii=ncols-n*6
c 
        do 2000 i=1,n
c 
        jj=(i-1)*6
        write(iw,1900) (rlams(jj+j),j=1,6)
 2000 continue
c 
        jj=n*6
        if(ii .ne. 0)
     1      write(iw,1900) (rlams(jj+j),j=1,ii)
c 
        write(iw,1200)
c 
        return
c 
c 
c 
c 
        entry quaestor(iw,xs,ys,ws,ns)
c 
 2200 format('c        ')
        write(iw,2200)
c 
 2400 format('c  Data for ',i3,' nodes')
        write(iw,2400) ns
        write(iw,2200)
c 
 2500 format('c        Nodes:')
        write(iw,2500)
        write(iw,2200)
c 
 2600 format(5x,i1,2x,e22.16,',',e22.16,',')
c 
        n=ns/2
        ii=ns-n*2
        do 2800 i=1,n
c 
	i0=mod(i,10)
        if(i0 .eq. 0) i0=11
        write(iw,2600) i0,xs(i*2-1),xs(i*2)
 2800 continue
c 
        if(ii .ne. 0) write(iw,2600) n+1,xs(ns)
c 
        write(iw,2200)
  
  
        n=ns/2
        ii=ns-n*2
        do 2900 i=1,n
c 
	i0=mod(i,10)
        if(i0 .eq. 0) i0=11
        write(iw,2600) i0,ys(i*2-1),ys(i*2)
 2900 continue
c 
        if(ii .ne. 0) write(iw,2600) n+1,ys(ns)
c 
        write(iw,2200)
 3200 format('c        Weights:')
        write(iw,3200)
        write(iw,2200)
c 
        do 3400 i=1,n
c 
	i0=mod(i,10)
        if(i0 .eq. 0) i0=11
        write(iw,2600) i0,ws(i*2-1),ws(i*2)
 3400 continue
c 
        if(ii .ne. 0) write(iw,2600) n+1,ws(ns)
c 
        write(iw,2200)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine quaeretr(ier,ir,n,x,ws)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),ws(1)
        character *1 card(100),Nodes(5),aaa1(100),weights(7)
        character *100 aaa
c 
        equivalence(aaa1(1),aaa)
c 
        data Nodes/'N','o','d','e','s'/,
     1      weights/'W','e','i','g','h','t','s'/
c 
c       form the line to be found in trhe file
c 
 1020 format('c  Data for ',i3,' nodes')
c 
        write(aaa,1020) n
c 
c      keep reading the file till we get to the right line
c 
 1200 format(80a1)
c 
        ier=8
        rewind(ir)
        numrec=0
        do 1350 i=1,10000
c 
        read(ir,1200,end=3000) (card(j),j=1,30)
c 
c       determine if this is the correct line
c 
        ifcorr=1
c 
        do 1300 j=1,20
c 
        if(card(j) .ne. aaa1(j)) ifcorr=0
 1300 continue
c 
        if(ifcorr .eq. 1) ier=0
        if(ifcorr .eq. 1) goto 1360
c 
 1350 continue
c 
 1360 continue
c 
c       find the nodes
c 
        do 1500 i=1,1000
c 
        read(ir,1200,end=3000) (card(j),j=1,30)
c 
        ifline=1
        do 1400 j=1,5
c 
        if(card(j+9) .ne. Nodes(j)) ifline=0
 1400 continue
c 
        if(ifline .eq. 1) goto 1600
 1500 continue
c 
        ier=16
        return
c 
 1600 continue
c 
        read(ir,1200,end=3000) (card(j),j=1,30)
c 
c       read the nodes
c 
        do 2600 i=1,n/2
c 
 2400 format(8x,e22.16,1x,e22.16)
c 
        read(ir,2400) x(2*i-1),x(2*i)
 2600 continue
c 
        n2=n/2
        j=n-n2*2
        if(j .ne. 0) read(ir,2400) x(n)
c 
c       find the weights
c 
        do 3500 i=1,1000
c 
        read(ir,1200,end=3000) (card(j),j=1,30)
c 
cccc        call prina('card as read*',card,20)
  
        ifline=1
        do 3400 j=1,5
c 
        if(card(j+9) .ne. weights(j)) ifline=0
 3400 continue
c 
        if(ifline .eq. 1) goto 3600
 3500 continue
c 
 3600 continue
c 
        read(ir,1200,end=3000) (card(j),j=1,30)
c 
c       read the nodes
c 
        do 4600 i=1,n/2
c 
 4400 format(8x,e22.16,1x,e22.16)
c 
        read(ir,4400) ws(2*i-1),ws(2*i)
 4600 continue
c 
        n2=n/2
        j=n-n2*2
        if(j .ne. 0) read(ir,2400) ws(n)
c 
 3000 continue
        return
        end
