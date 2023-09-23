        implicit real *8 (a-h,o-z)
c 
        real *8 zsin(2,10 000),work(100 000),
     1      zsout(2,100 000),wsin(10 000),wsout(10000),
     2      zsin2(2,10000),zsout2(2,10000),wsin2(10000),
     3      wsout2(10000),
     4      zsin3(2,10000),zsout3(2,10000),wsin3(10000),
     5      wsout3(10000),rlensin(10000),rlensout(10000),
     6      rlams(10000)
c 
        complex *16 rk,amatr(1000 000),u(1000 000),v(1000 000),
     1      s(1000),w(10 000 000),fs5(10000),
     2      v_big(1000 000),u_big(1000 000),
     3      sold(10000)
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
  
         eps=1.0e-14
  
c 
c       construct the skeleton
c 
        eps2=eps
        iplot=21
  
        ntsout=nts*2
  
        ntsin=nts
        ntsout=nts*2
        coef=3
c 
c        construct the big matrix of interactions and its SVD
c        via the accelerated scheme
c 
        lenw=10 000 000
  
  
        call square_basis(ier,poteval2,rk,nts,
     1      coef,gap,eps,iplot,
     2      zsin,wsin,rlensin,nin,zsout,wsout,rlensout,nout,
     3      rlams,nvects,u_big,v_big,w,lenw)
c 
c        construct the big matrix of interactions and its SVD
c        via the slow scheme
c 
        call prin2('after square_basis, rlams=*',rlams,nvects*2)
        call prinf('after square_basis, nvects=*',nvects,1)
c 
        call creaskel(ier,zsin,nin,wsin,
     1      zsout,nout,wsout,amatr,rk,eps,
     2      ncols,u,v,rlensout,fs5,s,w)
  
        call cleascop(s,sold,ncols)
  
        call cases_reord(v,s,ncols,
     1      v_big,rlams,nvects,nout)
c 
        istop=0
        call cases_test2(v,s,ncols,
     1      v_big,rlams,nvects,nout,istop)
c 
c        now, test the left singular vectors
c 
        call cases_reord(u,sold,ncols,
     1      u_big,rlams,nvects,nin)
  
        istop=0
  
        call cases_test2(u,sold,ncols,
     1      u_big,rlams,nvects,nin,istop)
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine cases_test2(v1,rlams1,ncols1,
     1      v2,rlams2,ncols2,n,istop)
        implicit real *8 (a-h,o-z)
        save
        complex *16 v1(n,ncols1),v2(n,ncols2),rlams1(1),
     1      rlams2(1),rat,ww(10000),diffs(10000)
c 
c        Scan the arrays rlams1, rlams2. When the two elements
c        are equal, check if the corresponding columns in
c        arrays v1,v2 are also equal, up to a constant complex
c        multiple of absolute value 1
c 
        ncols=ncols1
        if(ncols .gt. ncols2) ncols=ncols2
c 
        errmax=0
c 
        do 2000 i=1,ncols
c 
        imax=0
        dmax=0
        do 1200 jj=1,n
c 
        if(abs(v1(jj,i)) .gt. dmax) then
            imax=jj
            dmax=abs(v1(jj,i))
        endif
c 
 1200 continue
c 
        rat=v2(imax,i)/v1(imax,i)
c 
        do 1400 jj=1,n
c 
        ww(jj)=v1(jj,i)*rat
c 
        diffs(jj)=(ww(jj)-v2(jj,i))*rlams1(i)
c 
        if(abs(diffs(jj)) .gt. errmax) errmax=abs(diffs(jj))
 1400 continue
c 
cc        call prin2('in cases2_test, v1(1,i)=*',v1(1,i),n*2)
cc        call prin2('in cases2_test, v2(1,i)=*',v2(1,i),n*2)
cc        call prin2('in cases2_test, ww=*',ww,n*2)
cccc        call prin2('in cases2_test, diffs=*',diffs,n*2)
c 
        if( (istop .eq. 1) .and. (i .eq. 2)) stop
c 
 2000 continue
  
        call prin2('in cases2_test, errmax=*',errmax,1)
  
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
c 
c 
        subroutine creaskel(ier,zsin,nin,wsin,
     1      zsout,nout,wsout,a,rk,eps,ncols,u,v,rlensout,
     2      fs5,s,w)
        implicit real *8 (a-h,o-z)
        save
        dimension zsin(2,1),zsout(2,1),wsin(1),wsout(1),
     1      rlensout(1)
c 
        complex *16 rk,w(1),u(nin,1),
     1      v(nout,1),s(1),a(nin,nout),
     2      fs5(1)
c 
c        construct the matrix of interactions between the
c        squares
c 
        do 1600 i=1,nout
        do 1400 j=1,nin
c 
        call poteval2(zsout(1,i),zsin(1,j),rk,a(j,i))
c 
        a(j,i)=a(j,i)*sqrt(wsout(i)*wsin(j))
c 
 1400 continue
 1600 continue
c 
c      svd the matrix
c 
        lw=3000 000
        call csvdpiv(ier,a,nin,nout,u,v,s,ncols,eps,
     1      w,lw,ltot)
c 
        call prinf('in creaskel after csvdpiv, ncols=*',ncols,1)
        call prin2('in creaskel after csvdpiv, s=*',s,ncols*2)
c 
c       scan the eigenvalues. If two consecutive eigenvalues
c       are equal, symmetrize the corresponding singular
c       functions
c 
        ntsin=nin/8
        ntsout=nout/8
c 
        call symmetral(zsout,ntsout,v,ncols,s)
c 
        call symmetral(zsin,ntsin,u,ncols,s)
c 
c 
c        compensate the obtained vectors for the square roots
c        of weights
c 
        do 3800 i=1,ncols
c 
        do 3400 j=1,nin
c 
        u(j,i)=u(j,i)/sqrt(wsin(j))
 3400 continue
c 
        do 3600 j=1,nout
c 
        v(j,i)=v(j,i)/sqrt(wsout(j))
 3600 continue
c 
 3800 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine symmetral(zsout,ntsout,v,ncols,s)
        implicit real *8 (a-h,o-z)
        save
        dimension zsout(2,1)
        complex *16 v(ntsout*8,1),fs1(10000),fs2(10000),s(1),
     1      fs3(10000),fs4(10000),cd1,cd2,cd3,cd4
c 
c       scan the eigenvalues. If two consecutive eigenvalues
c       are equal, symmetrize the corresponding singular
c       functions
c 
        nout=ntsout*8
        do 3000 i=1,ncols-1
c 
        d=s(i+1)/s(i)
c 
        if( abs(d-1) .gt. 0.01) goto 3000
c 
c       symmetrize this singular vector and the
c       following one
c 
        call symmetr1(zsout,v(1,i),ntsout,fs1)
        call symmetr1(zsout,v(1+ntsout*4,i),ntsout,
     1      fs1(ntsout*4+1) )
c 
        call symmetr2(zsout,v(1,i),ntsout,fs2)
        call symmetr2(zsout,v(1+ntsout*4,i),ntsout,
     1      fs2(ntsout*4+1) )
c 
c 
        call symmetr1(zsout,v(1,i+1),ntsout,fs3)
        call symmetr1(zsout,v(1+ntsout*4,i+1),ntsout,
     1      fs3(ntsout*4+1) )
c 
        call symmetr2(zsout,v(1,i+1),ntsout,fs4)
        call symmetr2(zsout,v(1+ntsout*4,i+1),ntsout,
     1      fs4(ntsout*4+1) )
c 
c        depending on the relative sizes of fs1, fs3, replace
c        the i-th column of v with one of them
c 
        call cleascap(fs1,fs1,nout,cd1)
        call cleascap(fs3,fs3,nout,cd3)
c 
        d1=cd1
        d3=cd3
        d1=sqrt(d1)
        d3=sqrt(d3)
  
        do 2200 j=1,nout
c 
        fs1(j)=fs1(j)/d1
        fs3(j)=fs3(j)/d3
 2200 continue
  
        if(d1 .gt. d3) call cleascop(fs1,v(1,i),nout)
        if(d3 .gt. d1) call cleascop(fs3,v(1,i),nout)
c 
c        depending on the relative sizes of fs2, fs4, replace
c        the i+1-st column of v with one of them
c 
        call cleascap(fs2,fs2,nout,cd2)
        call cleascap(fs4,fs4,nout,cd4)
        d2=cd2
        d4=cd4
  
        d2=sqrt(d2)
  
        d4=sqrt(d4)
  
        do 2300 j=1,nout
c 
        fs2(j)=fs2(j)/d2
        fs4(j)=fs4(j)/d4
 2300 continue
c 
        if(d2 .gt. d4) call cleascop(fs2,v(1,i+1),nout)
        if(d4 .gt. d2) call cleascop(fs4,v(1,i+1),nout)
 3000 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine symmetr2(zs,fs,nts,fs2)
        implicit real *8 (a-h,o-z)
        save
        complex *16 fs(nts*4),fs2(nts*4),fswork(10000),cd
        dimension zs(2,1),z(2)
c 
c        perform the x-symmetrization
c 
  
cccc        call prin2('zs=*',zs,nts*8)
  
        do 1200 i=nts/2+1,5*nts/2
c 
        z(1)=-zs(1,i)
        z(2)=zs(2,i)
c 
        call square_symm_fnd(z,zs,nts*4,i2)
  
cccc        call prinf('i=*',i,1)
cccc        call prin2('and z=*',z,2)
ccccc        call prinf('and i2=*',i2,1)
c 
        cd=(fs(i)-fs(i2))/2
        fswork(i)=cd
        fswork(i2)=-cd
 1200 continue
  
cccc        call prin2('in symmetr1, fswork=*',fswork,nts*8)
  
c 
c        perform the y-antisymmetrization
c 
  
cccc        call prin2('zs=*',zs,nts*8)
  
        do 2200 i=nts*3/2+1,7*nts/2
c 
        z(1)=zs(1,i)
        z(2)=-zs(2,i)
c 
        call square_symm_fnd(z,zs,nts*4,i2)
  
cccc        call prinf('i=*',i,1)
cccc        call prin2('and z=*',z,2)
cccc        call prinf('and i2=*',i2,1)
c 
        cd=(fswork(i)+fswork(i2))/2
        fs2(i)=cd
        fs2(i2)=cd
 2200 continue
  
  
  
  
  
cccc        stop
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine symmetr1(zs,fs,nts,fs2)
        implicit real *8 (a-h,o-z)
        save
        complex *16 fs(nts*4),fs2(nts*4),fswork(10000),cd
        dimension zs(2,1),z(2)
c 
c        perform the x-symmetrization
c 
  
cccc        call prin2('zs=*',zs,nts*8)
  
        do 1200 i=nts/2+1,5*nts/2
c 
        z(1)=-zs(1,i)
        z(2)=zs(2,i)
c 
        call square_symm_fnd(z,zs,nts*4,i2)
c 
        cd=(fs(i)+fs(i2))/2
        fswork(i)=cd
        fswork(i2)=cd
 1200 continue
  
cccc        call prin2('in symmetr1, fswork=*',fswork,nts*8)
  
c 
c        perform the y-antisymmetrization
c 
        do 2200 i=nts*3/2+1,7*nts/2
c 
        z(1)=zs(1,i)
        z(2)=-zs(2,i)
c 
        call square_symm_fnd(z,zs,nts*4,i2)
c 
        cd=(fswork(i)-fswork(i2))/2
        fs2(i)=cd
        fs2(i2)=-cd
 2200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cases_reord(v1,rlams1,ncols1,
     1      v2,rlams2,ncols2,n)
        implicit real *8 (a-h,o-z)
        save
        complex *16 v1(n,ncols1),v2(n,ncols2),rlams1(1),
     1      rlams2(1),rat,ww(10000),diffs(10000)
c 
c        Scan the arrays rlams1, rlams2. When the two elements
c        are equal, check if the corresponding columns in
c        arrays v1,v2 are also equal, up to a constant complex
c        multiple of absolute value 1
c 
        ncols=ncols1
        if(ncols .gt. ncols2) ncols=ncols2
c 
        do 2000 i=1,ncols
c 
cccc        call prinf('i=*',i,1)
  
        d=(rlams2(i)-rlams2(i+1))/rlams2(i)
  
  
  
  
cccc        call prin2('d=*',d,1)
  
  
        if(d .gt. 0.0000001) goto 2000
  
  
c 
c       rlams(i) is equal to rlams(i+1). check if the corresponding
c       singular vectors should be reordered
c 
        err=0
  
        imax=0
        dmax=0
        do 1200 jj=1,n
c 
        if(abs(v1(jj,i)) .gt. dmax) then
            imax=jj
            dmax=abs(v1(jj,i))
        endif
c 
 1200 continue
c 
        rat=v2(imax,i)/v1(imax,i)
c 
        do 1400 jj=1,n
c 
        ww(jj)=v1(jj,i)*rat
c 
        diffs(jj)=(ww(jj)-v2(jj,i))
c 
        if(abs(diffs(jj)) .gt. err) err=abs(diffs(jj))
 1400 continue
  
cccc        call prinf('i=*',i,1)
  
cccc        call prin2('err=*',err,1)
c 
c       if the error is large, exchange the order of the
c       appropriate vectors
c 
        if(err .lt. 1.0d-4) goto 2000
c 
  
  
        call cleascop(v1(1,i),ww,n)
        call cleascop(v1(1,i+1),v1(1,i),n)
        call cleascop(ww,v1(1,i+1),n)
  
cccc        call prinf('i=*',i,1)
  
cccc        call prin2('w=*',w,n*2)
  
 2000 continue
c 
        return
        end
c 
c 
c 
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code, and the beginning
c        of the code constructing the SVD of the interaction
c        operator between two "concentric" squares in the plane
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
  
        subroutine square_basis(ier,poteval,rk,nts,
     1      coef,gap,eps,iplot,
     2      zsin,wsin,rlensin,nin,zsout,wsout,rlensout,nout,
     3      rlams,nvects,u_big,v_big,w2,lenw)
        implicit real *8 (a-h,o-z)
c 
        save
        real *8 zsin(2,1),zsout(2,1),wsin(1),wsout(1),
     1      rlensin(1),rlensout(1)
c 
        complex*16 rk,v_big(nts*16,1),rlams(1),u_big(nts*8,1)
c 
        real *8 w2(1)
c 
c       For the user-supplied subroutine poteval, this subroutine
c       constructs the singular values and the coresponding singular
c       vectors of the operator of interactions between two
c       "concentric" squares in the plane (see below for the somewhat
c       unusual definition of the "square" used in this code). The
c       subroutine also returnes the nodes on the two squares, the
c       corresponding quadrature weights, and several other paramters.
c 
c       In order to avoid the issue of "spurious resonances", each of
c       the squares is specified by a "double boundary":
c 
c 
c 
c               * * * * * * * * * * * * * * * * * * * * *
c               * * * * * * * * * * * * * * * * * * * * *
c               * *                                   * *
c               * *                                   * *
c               * *     * * * * * * * * * * * * *     * *
c               * *     * * * * * * * * * * * * *     * *
c               * *     * *                   * *     * *
c               * *     * *                   * *     * *
c               * *     * *                   * *     * *
c               * *     * *                   * *     * *
c               * *     * *                   * *     * *
c               * *     * *                   * *     * *
c               * *     * *                   * *     * *
c               * *     * *                   * *     * *
c               * *     * *                   * *     * *
c               * *     * *                   * *     * *
c               * *     * * * * * * * * * * * * *     * *
c               * *     * * * * * * * * * * * * *     * *
c               * *                                   * *
c               * *                                   * *
c               * * * * * * * * * * * * * * * * * * * * *
c               * * * * * * * * * * * * * * * * * * * * *
c 
c 
c        The outer boundary of the inner square is always the square
c        with the corners (1,1),(-1,1),(-1,-1),(1,-1). The inner
c        boundary of the outer square is the square with the corners
c        (coef,coef),(-coef,coef),(-coef,-coef),(coef,-coef), with
c        coef the user-supplied parameter; the distances between the
c        two layers of each of the double boundaries are equal to the
c        user-supplied parameter gap.
c 
c        The inner square will be discretized into 8*nts nodes, 4*nts
c        nodes on the inner layer, and 4*nts nodes on the outer one.
c        The outer square will be discretized into 16*nts nodes, 8*nts
c        nodes on the inner layer, and 8*nts nodes on the outer one.
c 
c                Input parameters:
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
c 
c  rk - the parameter (real, complex, whatever) used by the
c        subroutine poteval
c  nts - the number of nodes into which each of the layers
c        comprising the INNER square will be subdivided ON
c        EACH OF THE SIDES OF THE SQUARE. Thus, the total number
c        of nodes discretizing the inner square will be 8*nts.
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
c  rlams - singular values of the operator describing the interactions
c        between the squares. PLEASE NOTE THAT THIS IS A COMPLEX ARRAY
C  nvects - the number of elements in array rlams; also the number
c        of columns in each of the arrays u_big, v_big
c  u_big - complex array dimensioned (nin,nvects), containing the
c        left singular vectors of the operator describing the
c        interactions between the squares
c  v_big - complex array dimensioned (nout,nvects), containing the
c        right singular vectors of the operator describing the
c        interactions between the squares
c 
c                     Work array:
c 
c  w2 - should be reasonably large
c 
c       allocate memory
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
        lamatr=nin2*nout2*2+10
c 
        iu=iamatr+lamatr
        lu=nin2*nout2*2+10
c 
        iv=iu+lu
        lv=nin2*nout2*2+10
c 
        is=iv+lv
        ls=nout2*2+10
c 
        irlamsold=is+ls
        lrlamsold=nout2*2+10
c 
        ifs2=irlamsold+lrlamsold
        lfs2=nout*2+10
c 
        ifs6=ifs2+lfs2
        lfs6=nout*2+10
c 
c 
        izsin2=ifs6+lfs6
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
        izs2=iwsout3+lwsout3
        lzs2=nout*2+10
c 
        ltot=izs2+lzs2
        call prinf('ltot=*',ltot,1)
c 
        iww=ltot+1
cccc        stop
c 
        lleft=lenw-ltot
c 
        call prinf('and lleft=*',lleft,1)
c 
        if(lleft .lt. 20 000) then
            ier=16
            call prin2(
     1       'bombing from square_basis with ier=*',
     2          ier,1)
            return
        endif
c 
        call square_basis0(jer,poteval,rk,nts,coef,gap,eps,
     1      zsin,wsin,rlensin,nin,zsout,wsout,rlensout,nout,
     2      rlams,nvects,u_big,v_big,w2(iww),lleft,
     3      w2(iamatr),w2(iu),w2(iv),w2(is),w2(irlamsold),
cccc     4      w2(ifs2),
     4      w2(ifs2),w2(ifs6),
     5      w2(izsin2),w2(iwsin2),w2(izsin3),w2(iwsin3),
     6      w2(izsout2),w2(iwsout2),w2(izsout3),w2(iwsout3),
     7      w2(izs2),iplot)
c 
        if(jer .ne. 0) then
            ier=8
            call prinf(
     1       'bombing from square_basis with ier=*',
     2          ier,1)
        endif
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine square_basis0(ier,poteval,rk,nts,
     1      coef,gap,eps,
     2      zsin,wsin,rlensin,nin,zsout,wsout,rlensout,nout,
     3      rlams,nvects,u_big,v_big,w,lenw,
cccc     4      amatr,u,v,s,rlamsold,fs2,
     4      amatr,u,v,s,rlamsold,fs2,fs6,
     5      zsin2,wsin2,zsin3,wsin3,
     6      zsout2,wsout2,zsout3,wsout3,zs2,iplot)
        implicit real *8 (a-h,o-z)
c 
        save
        real *8 zsin(2,1),zsout(2,1),wsin(1),wsout(1),
     1      rlensin(1),rlensout(1),
     2      zsin2(2,1),zsout2(2,1),wsin2(1),
     3      wsout2(1),
     4      zsin3(2,1),zsout3(2,1),wsin3(1),
     5      wsout3(1),zs2(2,1)
c 
        complex*16 rk,amatr(1),u(1),v(1),w(1),s(1),fs2(1),
     1      v_big(nts*16,1),rlamsold(1),rlams(1),u_big(nts*8,1)
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
c        construct the symmetrized (icase=2) singular vectors
c 
        call square_basis2(jer,zsin2,wsin2,nin2,poteval,
     1      zsout2,wsout2,nout2,amatr,rk,eps,ncols2,
     2      u,v,s,w,lenw,zs2,fs2,zsout,v_big,rlams,
cccc     3      nvects,zsin,u_big)
     3      nvects,zsin,u_big,fs6)
c 
        if(jer .ne. 0) then
            ier=8
            return
        endif
c 
c        Construct the symmetrized (all \times 8 symmetries)
c        singular vectors
c 
        iii=nvects+1
c 
        do 3200 jcase=1,4
c 
        call square_basis3(jer,zsin3,wsin3,nin3,poteval,
     1      zsout3,wsout3,nout3,amatr,rk,eps,ncols3,jcase,
     2      u,v,rlams(iii),w,lenw,zsout,v_big(1,iii),
     3      u_big(1,iii),zsin,zs2)
c 
c 
        if(jer .ne. 0) then
            ier=32
            return
        endif
  
        iii=iii+ncols3
c 
 3200 continue
c 
        nvects=iii-1
        call cleascop(rlams,rlamsold,nvects)
c 
        call square_vects_sort(v_big,nout,rlams,nvects,w,
     1      w(nvects*nout+1))
c 
        call square_vects_sort(u_big,nin,rlamsold,nvects,w,
     1      w(nvects*nout+1))
c 
c        compensate the obtained vectors for the square roots
c        of weights
c 
        do 3800 i=1,nvects
c 
        do 3400 j=1,nin
c 
        u_big(j,i)=u_big(j,i)/sqrt(wsin(j))
 3400 continue
c 
        do 3600 j=1,nout
c 
        v_big(j,i)=v_big(j,i)/sqrt(wsout(j))
 3600 continue
c 
 3800 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine square_vects_sort(v,n,rlams,ncols,w,inds)
        implicit real *8 (a-h,o-z)
        save
        complex *16 v(n,ncols),rlams(1),cd,w(n,1)
        dimension inds(1)
c 
c        sort the singular values, and singular vectors
c        with them
c 
        do 1200 i=1,ncols
c 
        inds(i)=i
 1200 continue
c 
        do 2000 i=1,ncols
c 
        do 1800 j=1,ncols-1
c 
        d1=rlams(j)
        d2=rlams(j+1)
c 
        if(d1 .ge. d2) goto 1800
c 
        cd=rlams(j)
        rlams(j)=rlams(j+1)
        rlams(j+1)=cd
c 
        jj=inds(j)
        inds(j)=inds(j+1)
        inds(j+1)=jj
c 
 1800 continue
 2000 continue
c 
        do 2200 i=1,ncols
c 
        j=inds(i)
        call cleascop(v(1,j),w(1,i),n)
c 
 2200 continue
c 
        call cleascop(w,v,n*ncols)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine square_basis3(ier,zsin,wsin,nin,poteval,
     1      zsout,wsout,nout,a,rk,eps,ncols,icase,
     2      u,v,s,w,lenw,zs_big,v_big,u_big,zsin_big,
     3      zs2)
        implicit real *8 (a-h,o-z)
        save
        dimension zsin(2,1),zsout(2,1),wsin(1),wsout(1),
     1      z1(3),z2(3),zz2(3),zs2(2,1),zs_big(2,1),
     2      zsin_big(2,1)
c 
        complex *16 rk(1),w(1),u(nin,1),
     1      v(nout,1),s(1),a(nin,nout),cd00,cd10,cd01,cd11,
     2      cdd00,cdd10,cdd01,cdd11,cd,
     3      v_big(nout*8,1),u_big(nout*4,1)
c 
c       For the user-supplied subroutine poteval, this subroutine
c       constructs SOME OF the singular values and the corresponding
c       singular vectors of the operator of interactions between two
c       "concentric" squares in the plane. The singular functions
c       returned by this subroutine are those corresponding to simple
c       (not multiple) singular values; one call to this subroutine
c       constructs all singular functions corresponding to one type
c       of symmetry, specified by the user via the parameter icase
c       (see below). The parameter icase can assume values 1,2,3,4,
c       since there are 4 types of symmetry for which the sigular
c       values are simple).
c 
c       In truth, this subroutine works with 1/8-th of each square,
c       instead of working with the whole square (the algorithm
c       requires n^3 operations; replacing n with n/8 has an obvious
c       effect on the CPU time). Thus, it has to be supplied with
c       appropriate parts of the squares; fortunately (and by sheer
c       coincidence) such parts are produced by the subroutine
c       square_frames_crea (see).
c 
c ATTENTION! While this subroutine performs all calculations on a
c       small part of the squares (specifically, 1/8-th), it retuns
c       the results on the whole square! Thus, the columns of the
c       arrays u_big, v_big are of lengths nin*8, nout*8,
c       respectively!
c 
c                Input parameters:
c 
c  zsin - the nodes on the inner square
c  wsin - the weights corresponding to the nodes in the array wsin
c  nin - the number of nodes on the inner square
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
c  zsout - the nodes on the outer square
c  wsout - the weights corresponding to the nodes in the array wsout
c  nout - the number of nodes on the outer square
c 
c  rk - the parameter (real, complex, whatever) used by the
c        subroutine poteval
c  eps - the accuracy to which the SVDs will be performed
c  lenw - the length (in real *8 locations) of the work array
c        w supplied by the user
c  zs_big - discretization of the boundary OF THE OUTER SQUARE,
C        AS PRODUCED BY THE SUBROUTINE SQUARE_FRAMES_CREA (see).
c  zsin_big - discretization of the boundary OF THE INNER SQUARE,
C        AS PRODUCED BY THE SUBROUTINE SQUARE_FRAMES_CREA (see).
c  icase - integer parameter telling the subroutine which class
c        of singular functions to construct. Can assume values
c        1, 2, 3, 4.
c 
c                Output parameters:
c 
c  ier - error return code
C  ncols - the number of singular values of the appropriately
c        symmetrized operator that are greater than eps. Also the
c        number of elements returned in array s, and the number
c        of columns returned in each of the arrays u_big, v_big
c  s - singular values of the operator describing the interactions
c        between the squares. PLEASE NOTE THAT THIS IS A COMPLEX
c        ARRAY
c  u_big - complex array dimensioned (nin*8,ncols), containing the
c        left singular vectors of the operator describing the
c        interactions between the squares
c  v_big - complex array dimensioned (nout*8,ncols), containing the
c        right singular vectors of the operator describing the
c        interactions between the squares
c 
c                     Work arrays:
c 
c  a - should be at least nin*nout*2 real *8 locations
c  u - should be at least nin*nout*2 real *8 locations
c  v - should be at least nin*nout*2 real *8 locations
c  s - should be at least nout*2 real *8 locations
c  w - must be sufficiently long. 3*nin*nout+ 2*nin**2 +
c      2*nout+2*nin+100 complex *16 are always sufficient. However,
c      most of the time it is much less than that
c  zs2 - should be at least nout real *8 locations
c 
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
        if(icase .eq. 1)
     1      a(i,j)=cd00+cd01+cd10+cd11+
     2         (cdd00+cdd01+cdd10+cdd11)
c 
        if(icase .eq. 2)
     1      a(i,j)=cd00+cd01+cd10+cd11-
     2         (cdd00+cdd01+cdd10+cdd11)
c 
        if(icase .eq. 3)
     1      a(i,j)=cd00-cd01-cd10+cd11+
     2         (cdd00-cdd01-cdd10+cdd11)
c 
        if(icase .eq. 4)
     1      a(i,j)=cd00-cd01-cd10+cd11-
     2         (cdd00-cdd01-cdd10+cdd11)
c 
        a(i,j)=a(i,j)*sqrt(wsin(i)*wsout(j))
c 
 1400 continue
 1600 continue
c 
c      svd the matrix
c 
        call csvdpiv(ier,a,nin,nout,u,v,s,ncols,eps,
     1      w,lenw,ltot)
c 
c       one singular function after another, expand them
c       things so that they are defined on the whole
c       square (actually, on two squares)
c 
        do 3000 i=1,ncols
c 
        call square_expand(zsout,v(1,i),zs2,v_big(1,i),
     1      nout/2,zs_big,icase,w)
c 
        call cleascap(v_big(1,i),v_big(1,i),nout*8,cd)
        d=cd
        d=sqrt(d)
        do 2500 j=1,nout*8
c 
        v_big(j,i)=v_big(j,i)/d
 2500 continue
c 
        call square_expand(zsin,u(1,i),zs2,u_big(1,i),
     1      nin/2,zsin_big,icase,w)
c 
        call cleascap(u_big(1,i),u_big(1,i),nin*8,cd)
        d=cd
        d=sqrt(d)
        do 2600 j=1,nin*8
c 
        u_big(j,i)=u_big(j,i)/d
 2600 continue
c 
 3000 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine square_expand(zs,fs,zs2,fs2,nts,
     1      zs_big,icase,w)
        implicit real *8 (a-h,o-z)
        save
        complex *16 fs(1),fs2(1),w(1)
        dimension zs(2,1),zs2(2,1),zs_big(2,1)
c 
c        implement the reflection in the diagonal of the square
c 
        do 1100 i=1,nts*2
c 
        zs2(1,i)=zs(1,i)
        zs2(2,i)=zs(2,i)
c 
        fs2(i)=fs(i)
c 
        zs2(1,i+nts*2)=-zs(2,i)
        zs2(2,i+nts*2)=-zs(1,i)
c 
        fs2(i+nts*2)=fs(i)
c 
        if( (icase .eq. 2) .or. (icase .eq. 4) )
     1      fs2(i+nts*2)=-fs(i)
c 
 1100 continue
c 
c        implement the x-reflection
c 
        do 1200 i=1,nts*4
        zs2(1,i+nts*4)=-zs2(1,i)
        zs2(2,i+nts*4)=zs2(2,i)
c 
        fs2(i+nts*4)=fs2(i)
        if( (icase .eq. 3) .or. (icase .eq. 4) )
     1      fs2(i+nts*4)=-fs2(i)
c 
 1200 continue
c 
c        implement the y-reflection
c 
         do 1400 i=1,nts*8
c 
         zs2(1,i+nts*8)=zs2(1,i)
         zs2(2,i+nts*8)=-zs2(2,i)
c 
         fs2(i+nts*8)=fs2(i)
c 
        if( (icase .eq. 3) .or. (icase .eq. 4) )
     1      fs2(i+nts*8)=-fs2(i)
c 
 1400 continue
c 
c        reorder the obtained values so that they sit in
c        the correct positions on the boundary
c 
        call square_symm_allocate(zs_big,zs2,fs2,nts,w)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine square_symm_allocate(zs,zs2,fs,nts,w)
        implicit real *8 (a-h,o-z)
        save
        complex *16 fs(1),w(1)
        dimension zs(2,1),zs2(2,1)
c 
c        for one value of f after another, find the proper
c        location
c 
        do 1200 i=1,nts*16
c 
        call square_symm_fnd(zs2(1,i),zs,nts*16,i2)
        w(i2)=fs(i)
 1200 continue
c 
        do 1400 i=1,nts*16
        fs(i)=w(i)
 1400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine square_basis2(ier,zsin,wsin,nin,poteval,
     1      zsout,wsout,nout,a,rk,eps,ncols,
     2      u,v,s,w,lenw,zs2,fs2,zsout_big,v_big,
     3      rlams,nvects,zsin_big,u_big,fs6)
        implicit real *8 (a-h,o-z)
        save
        dimension zsin(2,1),zsout(2,1),wsin(1),wsout(1),
     1      z1(3),z2(3),zs2(2,1),zsout_big(2,1),
     2      zsin_big(2,1)
c 
        complex *16 rk(2),w(2),u(nin,1),cd,
     1      v(nout,1),s(1),a(nin,nout),cd00,cd10,
     2      cd01,cd11,fs2(1),v_big(nout*4,1),fs6(1),
     3      rlams(1),u_big(nin*4,1)
c 
c       For the user-supplied subroutine poteval, this subroutine
c       constructs SOME OF the singular values and the corresponding
c       singular vectors of the operator of interactions between two
c       "concentric" squares in the plane. The singular functions
c       returned by this subroutine are those corresponding to double
c       (not simple) singular values; one call to this subroutine
c       constructs all singular functions corresponding to all such
c       singular vectors.
c 
c       In truth, this subroutine works with 1/4-th of each square,
c       instead of working with the whole square (the algorithm
c       requires n^3 operations; replacing n with n/4 has an obvious
c       effect on the CPU time). Thus, it has to be supplied with
c       appropriate parts of the squares; fortunately (and by sheer
c       coincidence) such parts are produced by the subroutine
c       square_frames_crea (see).
c 
c ATTENTION! While this subroutine performs all calculations on a
c       small part of the square (specifically, 1/4-th), it retuns
c       the results on the whole square! Thus, the columns of the
c       arrays u_big, v_big are of lengths nin*4, nout*4,
c       respectively!
c 
c                Input parameters:
c 
c  zsin - the nodes on the inner square
c  wsin - the weights corresponding to the nodes in the array wsin
c  nin - the number of nodes on the inner square
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
c  zsout - the nodes on the outer square
c  wsout - the weights corresponding to the nodes in the array wsout
c  nout - the number of nodes on the outer square
c 
c  rk - the parameter (real, complex, whatever) used by the
c        subroutine poteval
c  eps - the accuracy to which the SVDs will be performed
c  lenw - the length (in real *8 locations) of the work array
c        w supplied by the user
c  zsout_big - discretization of the boundary OF THE OUTER SQUARE,
C        AS PRODUCED BY THE SUBROUTINE SQUARE_FRAMES_CREA (see).
c  zsin_big - discretization of the boundary OF THE INNER SQUARE,
C        AS PRODUCED BY THE SUBROUTINE SQUARE_FRAMES_CREA (see).
c 
c                Output parameters:
c 
c  ier - error return code
c  u_big - complex array dimensioned (nin*8,ncols), containing the
c        left singular vectors of the operator describing the
c        interactions between the squares
c  v_big - complex array dimensioned (nout*8,ncols), containing the
c        right singular vectors of the operator describing the
c        interactions between the squares
c  rlams - singular values of the operator describing the interactions
c        between the squares. PLEASE NOTE THAT THIS IS A COMPLEX
c        ARRAY
C  nvects - the number of singular values of the appropriately
c        symmetrized operator that are greater than eps. Also the
c        number of elements returned in array rlams, and the number
c        of columns returned in each of the arrays u_big, v_big
c 
c                     Work arrays:
c 
c  a - should be at least nin*nout*2 real *8 locations
c  u - should be at least nin*nout*2 real *8 locations
c  v - should be at least nin*nout*2 real *8 locations
c  s - should be at least nout*2 real *8 locations
c  w - must be sufficiently long. 3*nin*nout+ 2*nin**2 +
c      2*nout+2*nin+100 complex *16 are always sufficient. However,
c      most of the time it is much less than that
c  zs2 - should be at least nout real *8 locations
c  fs2 - should be at least nout real *8 locations
c  fs6 - should be at least nout real *8 locations
c 
c 
c        . . . construct the matrix of interactions between the
c              squares
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
        a(i,j)=cd00+cd10-cd01-cd11
        a(i,j)=a(i,j)*sqrt(wsin(i)*wsout(j))
c 
 1400 continue
 1600 continue
c 
c      svd the matrix
c 
        call csvdpiv(ier,a,nin,nout,u,v,s,ncols,eps,
     1      w,lenw,ltot)
c 
c       expand the right singular vectors to make them
c       functions on the whole square
c 
        ntsout=nout/2
c 
        ii=0
        do 2600 i=1,ncols
c 
        call square_expand2(zsout,v(1,i),zs2,fs2,ntsout)
c 
        do 2400 j=1,nout*4
c 
        i2=0
        call square_symm_fnd(zs2(1,j),zsout_big,nout*4,i2)
c 
        fs6(i2)=fs2(j)
c 
 2400 continue
c 
        call cleascap(fs6,fs6,nout*4,cd)
        d=cd
        d=sqrt(d)
        do 2500 j=1,nout*4
c 
        fs6(j)=fs6(j)/d
 2500 continue
c 
        ii=ii+1
        rlams(ii)=s(i)
        call cleascop(fs6,v_big(1,ii),nout*4)
c 
        call square_rot(v_big(1,ii),fs6,nout/2)
c 
        ii=ii+1
        rlams(ii)=s(i)
c 
        call cleascop(fs6,v_big(1,ii),nout*4)
c 
 2600 continue
c 
        nvects=ii
c 
c       expand the left singular vectors to make them
c       functions on the whole square
c 
        ntsin=nin/2
c 
        ii=0
        do 3600 i=1,ncols
c 
        call square_expand2(zsin,u(1,i),zs2,fs2,ntsin)
c 
        do 3400 j=1,nin*4
c 
        i2=0
        call square_symm_fnd(zs2(1,j),zsin_big,nin*4,i2)
c 
        fs6(i2)=fs2(j)
c 
 3400 continue
c 
        call cleascap(fs6,fs6,nin*4,cd)
        d=cd
        d=sqrt(d)
        do 3500 j=1,nin*4
c 
        fs6(j)=fs6(j)/d
 3500 continue
c 
        ii=ii+1
        rlams(ii)=s(i)
        call cleascop(fs6,u_big(1,ii),nin*4)
c 
        call square_rot(u_big(1,ii),fs6,nin/2)
c 
        ii=ii+1
        rlams(ii)=s(i)
c 
        call cleascop(fs6,u_big(1,ii),nin*4)
c 
 3600 continue
c 
        nvects=ii
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine square_expand2(zs,fs,zs2,fs2,nts)
        implicit real *8 (a-h,o-z)
        save
        complex *16 fs(1),fs2(1)
        dimension zs(2,1),zs2(2,1)
c 
c        implement the x-reflection
c 
        do 1200 i=1,nts*2
        zs2(1,i)=zs(1,i)
        zs2(2,i)=zs(2,i)
c 
        fs2(i)=fs(i)
c 
c 
        zs2(1,i+nts*2)=-zs(1,i)
        zs2(2,i+nts*2)=zs(2,i)
        fs2(i+nts*2)=fs(i)
 1200 continue
c 
c        implement the y-reflection
c 
         do 1400 i=1,nts*4
c 
         zs2(1,i+nts*4)=zs2(1,i)
         zs2(2,i+nts*4)=-zs2(2,i)
c 
        fs2(i+nts*4)=-fs2(i)
 1400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine square_rot(fs1,fs2,nts)
        implicit real *8 (a-h,o-z)
        save
        complex *16 fs1(1),fs2(1)
c 
c       rotate the function
c 
        do 1200 j=1,nts
c 
        fs2(j+nts)=fs1(j)
        fs2(j+nts*2)=fs1(j+nts)
        fs2(j+nts*3)=fs1(j+nts*2)
c 
        fs2(j+nts*5)=fs1(j+nts*4)
        fs2(j+nts*6)=fs1(j+nts*5)
        fs2(j+nts*7)=fs1(j+nts*6)
c 
 1200 continue
c 
        do 2200 j=1,nts
c 
        fs2(j)=fs1(j+3*nts)
        fs2(j+4*nts)=fs1(j+7*nts)
 2200 continue
  
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine square_symm_fnd(z,zs,n,i2)
        implicit real *8 (a-h,o-z)
        save
        real *8 z(2),zs(2,1)
c 
c        find in the array zs the element coinciding with z
c 
        eps=1.0d-10
        i2=-1
        do 1200 i=1,n
c 
        d=abs(zs(1,i)-z(1))+abs(zs(2,i)-z(2))
c 
        if(d .gt. eps) goto 1200
c 
        i2=i
        return
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine square_frames_crea(ntsin,ntsout,coef,gap,
     1      zsin,wsin,nin,zsout,wsout,nout,
     2      zsin2,wsin2,nin2,zsout2,wsout2,nout2,
     3      zsin3,wsin3,nin3,zsout3,wsout3,nout3,
     4      rlensin,rlensout,work,iplot)
c 
        implicit real *8 (a-h,o-z)
        save
        dimension zsin(2,1),zsout(2,1),wsin(1),wsout(1),
     1      zsin2(2,1),zsout2(2,1),wsin2(1),wsout2(1),
     2      zsin3(2,1),zsout3(2,1),wsin3(1),wsout3(1),
     3      work(1),rlensin(1),rlensout(1)
c 
c        This subroutine constructs the frames on the whole
c        square, on the square after symmetrization by a factor
c        of 4, and  on the square after symmetrization by the
c        factor of 8. It also returns the cumulative arc-lengths
c        on both squares, to be used predominantly for plotting
c        things
c 
c                Input parameters:
c 
c  ntsin - the number of nodes on one layer of one side of the
c        inner square
c  ntsout - the number of nodes on one layer of one side of the
c        outer square
c  coef - the ratio between the size of the inner layer of the
c        outer square and the size of the outer layer of the
c        inner square; we like to set it to 3
c  gap - the distances between the two layers (outer and inner)
c        on either square
c  iplot - integer parameter telling the subroutine where to plot
c        the three pictures it wants to plot (the plots are in
c        GNUPLOT-readable format). It will plot the data on gn*
c        files, with *=iplot, iplot+1,iplot+2. Setting iplot=0
c        will suppress the plotting.
c 
c                Output parameters:
  
cc        subroutine square_frames_crea(ntsin,ntsout,coef,gap,
cc     1      zsin,wsin,nin,zsout,wsout,nout,
cc     2      zsin2,wsin2,nin2,zsout2,wsout2,nout2,
cc     3      zsin3,wsin3,nin3,zsout3,wsout3,nout3,
cc     4      rlensin,rlensout,work,iplot)
  
  
c 
c  zsin - the nodes on the inner square
c  wsin - the weights corresponding to the nodes in the array wsin
c  nin - the number of nodes on the inner square (number of elements
c        in arrays zsin, wsin, rlensin). Will be set to ntsin*8
c  zsout - the nodes on the outer square
c  wsout - the weights corresponding to the nodes in the array wsout
c  nout - the number of nodes on the outer square (number of elements
c        in arrays zsout, wsout, rlensout). Will be set to ntsout*8
c  zsin2 - the nodes on the one quarter of the inner square, to be
c          used by the subroutine square_basis2 (see)
c  wsin2 - the weights corresponding to the nodes in the array zsin2
c  nin2 - the number of nodes on the one quarter of the inner square
c        (number of elements in arrays zsin2, wsin2). Will be set to
c        ntsin*2
c  zsout2 - the nodes on the one quarter of the outer square
c  wsout2 - the weights corresponding to the nodes in the array zsout2
c  nout2 - the number of nodes on the one quarter of the outer square
c        (number of elements in arrays zsout2, wsout2). Will be set
c        to ntsout*2
c  zsin3 - the nodes on the 18-th of the inner square, to be
c          used by the subroutine quare_basis3 (see)
c  wsin3 - the weights corresponding to the nodes in the array zsin3
c  nin3 - the number of nodes on the 1/-th of the inner square
c        (number of elements in arrays zsin3, wsin3). Will be set to
c        ntsin
c  zsout3 - the nodes on the 1/8-th of the outer square
c  wsout3 - the weights corresponding to the nodes in the array zsout3
c  nout3 - the number of nodes on the 1/8-th of the outer square
c        (number of elements in arrays zsout3, wsout3). Will be set
c        to ntsout
c  rlensin - the arc length on the inner square (distance along the
c        surface from the upper right corner, tabulated at the
c        nodes zsin; expected to be used mostly for plotting
c  rlensout - the arc length on the outer square (distance along the
c        surface from the upper right corner, tabulated at the
c        nodes zsout; expected to be used mostly for plotting
c 
c                       Work array:
c 
c  work - must be at least max(ntsin,ntsout) real *8 elements long
c 
c 
c        . . . construct the frames on the whole square
c 
        call squares_cre(gap,ntsin,ntsout,coef,
     1      zsin,wsin,zsout,wsout,work)
c 
        if(iplot .gt. 0)
     1      call zquaplot2(iplot,zsin,ntsin*8,2,zsout,
     2      ntsout*8,2,'both squares*')
c 
c        Construct the frames on the whole on the square after
c        symmetrization by a factor of 4
c 
  
        call square_extract(zsin,wsin,zsin2,wsin2,ntsin)
        call square_extract(zsout,wsout,zsout2,wsout2,ntsout)
c 
        if(iplot .gt. 0) then
            iw=iplot+1
            call zquaplot2(iw,zsin2,ntsin*2,2,zsout2,ntsout*2,2,
     1      'quarters of both squares*')
        endif
c 
c        Construct the frames on the square after symmetrization
c        by the factor of 8
c 
        call square_extract2(zsin2,wsin2,zsin3,wsin3,ntsin)
        call square_extract2(zsout2,wsout2,zsout3,wsout3,ntsout)
c 
c 
        if(iplot .gt. 0) then
            iw=iplot+2
            call zquaplot2(iw,zsin3,ntsin,2,zsout3,ntsout,2,
     1      '1/8-ths of both squares*')
        endif
c 
        nin=ntsin*8
        nout=ntsout*8
c 
        nin2=ntsin*2
        nout2=ntsout*2
c 
        nin3=ntsin
        nout3=ntsout
c 
c        construct the cumulative arclengths on the squares
c 
        call square_arcs(ntsin,rlensin,work)
        call square_arcs(ntsout,rlensout,work)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine square_arcs(n,rlens,ts)
        implicit real *8 (a-h,o-z)
        save
        real *8 ts(1),rlens(1)
c 
c        construct the cumulative arclengths on all three
c        sets of points, to be used predominantly
c        for plotting
c 
        itype=0
        call legeexps(itype,n,ts,u,v,whts)
c 
  
        do 1200 i=1,n
c 
        rlens(i)=1+ts(i)
        rlens(n+i)=3+ts(i)
        rlens(2*n+i)=5+ts(i)
        rlens(3*n+i)=7+ts(i)
 1200 continue
c 
        do 1400 i=1,n*4
c 
        rlens(n*4+i)=rlens(i)
 1400 continue
  
        return
        end
c 
c 
c 
c 
c 
        subroutine square_extract2(zs,ws,zs3,ws3,nts)
        implicit real *8 (a-h,o-z)
        save
        dimension zs(2,1),ws(1),zs3(2,1),ws3(1)
c 
c        extract from the big arrays the geometry pertaining
c        to the upper and left sides of both squares
c 
        do 1600 i=1,nts/2
c 
        zs3(1,i)=zs(1,i)
        zs3(2,i)=zs(2,i)
        ws3(i)=ws(i)
c 
        zs3(1,i+nts/2)=zs(1,i+nts)
        zs3(2,i+nts/2)=zs(2,i+nts)
        ws3(i+nts/2)=ws(i+nts)
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
        subroutine square_extract(zs,ws,zs3,ws3,nts)
        implicit real *8 (a-h,o-z)
        save
        dimension zs(2,1),ws(1),zs3(2,1),ws3(1)
c 
c        extract from the big arrays the geometry pertaining
c        to the upper and left sides of both squares
c 
        do 1600 i=1,nts
c 
        zs3(1,i)=zs(1,i+nts/2)
        zs3(2,i)=zs(2,i+nts/2)
        ws3(i)=ws(i+nts/2)
c 
        zs3(1,i+nts)=zs(1,i+nts/2*3+nts*3)
        zs3(2,i+nts)=zs(2,i+nts/2*3+nts*3)
        ws3(i+nts)=ws(i+nts/2*3+nts*3)
  
 1600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine squares_cre(gap,ntsin,ntsout,coef,
     1      zsin,wsin,zsout,wsout,w)
        implicit real *8 (a-h,o-z)
c 
        save
        real *8 zsin(2,1),wsin(1),zsout(2,1),wsout(1),w(1)
c 
c        construct the inner square
c 
        its=1
        lts=ntsout+2
c 
        iwhts=its+lts
        lwhts=ntsout*2
c 
        call square_cre(gap,ntsin,zsin,wsin,w(its),w(iwhts) )
c 
c        construct the outer square
c 
        cc=coef+gap
        gapout=gap/cc
c 
        call square_cre(gapout,ntsout,zsout,wsout,w(its),w(iwhts) )
c 
c       scale the outer square appropriately
c 
        do 1600 i=1,ntsout*8
c 
        zsout(1,i)=zsout(1,i)*cc
        zsout(2,i)=zsout(2,i)*cc
        wsout(i)=wsout(i)*coef
 1600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine square_cre(gap,nts,
     1      zs,weights,ts,whts)
        implicit real *8 (a-h,o-z)
c 
        save
        real *8 zs(2,1),ts(10 000),
     1      whts(10 000),weights(1)
c 
c        construct the inner square
c 
        itype=1
        call legeexps(itype,nts,ts,u,v,whts)
c 
        x1=-1
        x2=1
c 
        y1=-1
        y2=1
c 
        call square_lineconstr(x2,y2,x1,y2,ts,nts,zs)
        call square_lineconstr(x1,y2,x1,y1,ts,nts,zs(1,nts+1))
        call square_lineconstr(x1,y1,x2,y1,ts,nts,zs(1,2*nts+1))
        call square_lineconstr(x2,y1,x2,y2,ts,nts,zs(1,3*nts+1))
c 
        x1=x1+gap
        x2=x2-gap
c 
        y1=y1+gap
        y2=y2-gap
c 
        call square_lineconstr(x2,y2,x1,y2,ts,nts,zs(1,4*nts+1))
        call square_lineconstr(x1,y2,x1,y1,ts,nts,zs(1,5*nts+1))
        call square_lineconstr(x1,y1,x2,y1,ts,nts,zs(1,6*nts+1))
        call square_lineconstr(x2,y1,x2,y2,ts,nts,zs(1,7*nts+1))
c 
c        construct the array of weights
c 
        do 1600 i=1,nts
c 
        weights(i)=whts(i)
        weights(i+nts)=whts(i)
        weights(i+2*nts)=whts(i)
        weights(i+3*nts)=whts(i)
        weights(i+4*nts)=whts(i)
        weights(i+5*nts)=whts(i)
        weights(i+6*nts)=whts(i)
        weights(i+7*nts)=whts(i)
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
        subroutine square_lineconstr(x1,y1,x2,y2,ts,nts,zs)
        implicit real *8 (a-h,o-z)
        save
        dimension ts(1),zs(2,1)
c 
c        construct the parametric equation of the line interval
c        connecting the points (x1,y1), (x2,y2), so that the
c        "general" point on that interval has the form
c 
c        (x,y)=(ax*t+bx,ay*t+by),
c 
c        with t \in [t1,t2]
c 
        t1=-1
        t2=1
c 
        ax=(x2-x1)/(t2-t1)
        bx=x1-ax*t1
c 
        ay=(y2-y1)/(t2-t1)
        by=y1-ay*t1
c 
        do 1200 i=1,nts
c 
        zs(1,i)=ax*ts(i)+bx
        zs(2,i)=ay*ts(i)+by
 1200 continue
c 
        return
        end
