      implicit real *8 (a-h,o-z)
        dimension w(2000 000)
  
        dimension zs(2,100 000),zs2(2,100 000),
     1      iz(100 000),izinv(100000),thresh(2,10000),
     2      ichunks(2,10000),bord(2,10 000),thresh2(2,10000),
     3      ichunks2(2,10000)
c 
        dimension expand_inout2(100 000),irows_outin2(100 000),
     1      icols_inout2(100 000),eval_outin2(100 000)
  
c 
        dimension expand_inout1(100 000),irows_outin1(100 000),
     1      icols_inout1(100 000),eval_outin1(100 000)
c 
        dimension expand(200 000),icols(1000),
     1      eval(200 000),
     2      irows(100 000),amatr(1000 000),s1(200 000),
     3      s2(200 000),sdad(500 000)
c 
        dimension t1(200 000),t2(200 000),t1bar(200 000),
     1      t2bar(200 000),t12(200 000),t21(200 000),
     2      sout(200 000),bb(200 000),cc(200 000),
     3      sdadout(200 000)
  
        dimension split1(200 000),split2(200 000),
     1      exchg11(200 000),exchg12(200 000),exchg22(200 000),
     2      exchg21(200 000),cmerge1(200 000),cmerge2(200 000),
     3      spnts_inout_dad(200 000),iwhich(10),
     4      spnts1(200 000),spnts_inout1(200 000),
     5      outin_spnts1(200 000),outin_spnts2(200 000),
     6      spnts2(200 000),spnts_inout2(200 000)
c 
        dimension conds(10)
c 
        complex *16 cd
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
c 
        PRINT *, 'ENTER num'
        READ *,num
        CALL PRINf('num=*',num,1 )
  
c 
c       construct the geometry of the test
c 
        eps=1.0d-11
  
  
c 
        n=num
        rlength=60
        rlength=20
        hight=1
ccc        call creaband(rlength,hight,n,zs,w)
        call creaband3(rlength,hight,n,zs,w)
c 
ccccc        call prin2('after creaband, zs=*',zs,n*2)
c 
  
  
  
        iw=21
        call zquaplot(iw,zs,num,1,'the band*')
c 
c       sort the obtained points
c 
  
  
        nchunks=16
  
        endleft=-rlength/2
        endright=rlength/2
  
  
        call bandgeo(ier,zs,n,nchunks,endleft,endright,
     1      zs2,iz,izinv,thresh,ichunks,w)
        call prinf('after bandgeo, ier=*',ier,1)
c 
c        construct the border of the chunk and plot the whole mess
c 
        ichunk1=3
c 
        i=ichunks(1,ichunk1)
        np=ichunks(2,ichunk1)-ichunks(1,ichunk1)+1
        np1=np
  
        iw=22
        call bandplot(iw,zs2,n,thresh,nchunks,hight,bord,ii)
  
        lbord=ii
  
        iw=23
        call zquaplot2(iw,zs2(1,i),np,1,bord,ii,3,
     1      'one chunk*')
c 
c       construct the skeleton of one chunk
c 
        npts=150
        iplot=33
c 
        call skels_1chunk(n,eps,zs,
     1      nchunks,ichunk1,npts,bord,lbord,iplot,
     2      zs2,iz,izinv,thresh,ichunks,w,
     3      eval_outin1,irows_outin1,ncols_outin1,
     4      expand_inout1,icols_inout1,ncols_inout1)
  
        np1=ichunks(2,ichunk1)-ichunks(1,ichunk1)+1
  
  
 4200 continue
        ichunk2=4
        iplot=34
  
        call prinf('ichunks(1,ichunk2)=*',ichunks(1,ichunk2),1)
  
cccc          stop
  
c 
        call skels_1chunk(n,eps,zs,
     1      nchunks,ichunk2,npts,bord,lbord,iplot,
     2      zs2,iz,izinv,thresh,ichunks,w,
     3      eval_outin2,irows_outin2,ncols_outin2,
     4      expand_inout2,icols_inout2,ncols_inout2)
  
        np2=ichunks(2,ichunk2)-ichunks(1,ichunk2)+1
  
c 
c        construct the threshold and the indices for the
c        next level
c 
 5000 continue
c 
        nchunks2=nchunks/2
        ii=0
        do 5200 i=1,nchunks2
c 
        ii=ii+1
        thresh2(1,i)=thresh(1,ii)
        ichunks2(1,i)=ichunks(1,ii)
        ii=ii+1
        thresh2(2,i)=thresh(2,ii)
        ichunks2(2,i)=ichunks(2,ii)
 5200 continue
c 
        call prin2('thresh2=*',thresh2,nchunks2*2)
        call prin2('while thresh=*',thresh,nchunks*2)
c 
        barrier=thresh(2,1)-thresh(1,1)
  
c 
c        merge the skeletons of the two chunks, obtaining a
c        chunk on the next level
c 
cccc          eps=eps*100
  
        ichunk22=2
  
        call merges_1chunk(zs,zs2,iz,izinv,n,npts,
     1    eps,barrier,w,thresh2,ichunks2,ichunk22,
     2    nchunks/2,irows_outin1,ncols_outin1,
     3    irows_outin2,ncols_outin2,
     4    icols_inout1,ncols_inout1,
     5    icols_inout2,ncols_inout2,
     6    eval,irows,nrows,expand,icols,ncols)
  
  
        ncols_inout=ncols
        ncols_outin=nrows
  
  
        iplot=47
        call plot_inner_skeleton(zs2,irows,nrows,
     1      bord,lbord,iplot,zs,iz)
  
  
  
        iplot=48
        call plot_inner_skeleton(zs2,icols,ncols,
     1      bord,lbord,iplot,zs,iz)
  
c 
c        construct the scattering matrix for the first box
c 
        istart=ichunks(1,ichunk1)
        np=ichunks(2,ichunk1)-ichunks(1,ichunk1)+1
  
        iwhich(1)=1
        iwhich(2)=1
        iwhich(3)=1
c 
        call scatmatr0(ichunks,ichunk1,np,iz,zs,
     1      expand_inout1,eval_outin1,ncols_inout1,
     2      ncols_outin1,spnts1,
     2      s1,spnts_inout1,outin_spnts1,w,iwhich,cond)
  
  
  
  
c 
c        construct the scattering matrix for the second box
c 
        istart=ichunks(1,ichunk2)
  
        np=ichunks(2,ichunk2)-ichunks(1,ichunk2)+1
c 
        call scatmatr0(ichunks,ichunk2,np,iz,zs,
     1      expand_inout2,eval_outin2,ncols_inout2,
     2      ncols_outin2,spnts2,
     2      s2,spnts_inout2,outin_spnts2,w,iwhich,cond)
  
  
c 
c        construct the scattering matrix for the daddy
c 
        istart=ichunks(1,ichunk1)
        np=ichunks(2,ichunk2)-ichunks(1,ichunk1)+1
c 
  
        call prinf('ichunks(1,ichunk1)=*',ichunks(1,ichunk1),1)
        call prinf('ichunks(2,ichunk2)=*',ichunks(2,ichunk2),1)
  
        call prinf('ichunk1=*',ichunk1,1)
        call prinf('ichunk2=*',ichunk2,1)
  
        call prinf('np=*',np,1)
c 
        call scatmatr0dad(istart,np,iz,zs,amatr,
     1      expand_inout1,eval_outin1,ncols_inout,
     2      ncols_outin,sdad)
  
cccc        call prin2('after scatmatr0dad, sdad=*',sdad,
cccc     1      ncols_inout*ncols_outin*2)
  
cccc          stop
  
  
  
c 
c        construct the daddy's scattering matrix directly
c 
  
        np1=ichunks(2,ichunk1)-ichunks(1,ichunk1)+1
        np2=ichunks(2,ichunk2)-ichunks(1,ichunk2)+1
  
  
        call sdad_direct(sdad,ncols_inout,
     1      ncols_outin,eval,ncols_outin1,
     2      ncols_outin2,expand,
     4      expand_inout1,ncols_inout1,np1,expand_inout2,
     5      ncols_inout2,np2,eval_outin1,
     6      eval_outin2,
     7      bb,cc,sdadout,spnts_inout_dad)
c 
c       construct the daddy's scattering matrix by merging the
c       scattering matrices of the kids
c 
  
        lw=2000 000
        iflast=0
        call scattmatr_merge2(ier,expand,eval,
     1      iz,zs,
     2    irows_outin1,ncols_outin1,
     3    irows_outin2,ncols_outin2,
     4    icols_inout1,ncols_inout1,
     5    icols_inout2,ncols_inout2,
     6    irows_outin,ncols_outin,
     7    icols_inout,ncols_inout,s1,s2,sout,split1,split2,
     8      exchg11,exchg12,exchg22,exchg21,cmerge1,
     9      cmerge2,conds,w,lw,ltot,iflast)
  
        call cleasubt(sout,sdadout,w,ncols_outin*ncols_inout)
c 
        call cleascap(w,w,ncols_outin*ncols_inout,cd)
  
        d=cd
        d=sqrt(d)
        call prin2('norm of sout-sdadout is*',d,1)
c 
c        test the merging matrices
c 
        npson1=ichunks(2,ichunk1)-ichunks(1,ichunk1)+1
        npson2=ichunks(2,ichunk2)-ichunks(1,ichunk2)+1
  
        npdad=ichunks(2,ichunk2)-ichunks(1,ichunk1)+1
  
  
        call cmergetest2(spnts_inout_dad,npdad,
     1      ncols_inout,spnts_inout1,npson1,ncols_inout1,
     2      spnts_inout2,npson2,ncols_inout2,
     3      cmerge1,cmerge2)
  
  
        call exchg_test(sdad,npdad,
     1      ncols_inout,spnts_inout1,npson1,ncols_inout1,
     2      spnts_inout2,npson2,ncols_inout2,
     3      exchg11,exchg12,exchg22,exchg21,
     4      ncols_outin1,ncols_outin2,
     5      spnts1,spnts2,outin_spnts1,outin_spnts2)
  
  
  
        call prin2('and conds=*',conds,2)
  
  
  
        call prinf('by the way, ncols_inout=*',
     1      ncols_inout,1)
  
        call prinf('by the way, ncols_outin=*',
     1      ncols_outin,1)
  
  
        call prinf('by the way, ncols_inout1=*',
     1      ncols_inout1,1)
  
        call prinf('by the way, ncols_outin1=*',
     1      ncols_outin1,1)
  
  
  
        call prinf('by the way, ncols_inout2=*',
     1      ncols_inout2,1)
  
        call prinf('by the way, ncols_outin2=*',
     1      ncols_outin2,1)
c 
        stop
        end
c 
c 
c 
c 
c 
        subroutine exchg_test(sdad0,npdad,
     1      ncols_inout,spnts_inout1,npson1,ncols_inout1,
     2      spnts_inout2,npson2,ncols_inout2,
     3      exchg11,exchg12,exchg22,exchg21,
     4      ncols_outin1,ncols_outin2,
     5      spnts1,spnts2,outin_spnts1,outin_spnts2)
c 
        implicit real *8 (a-h,o-z)
        save
        complex *16
     1      spnts_inout1(ncols_inout1,npson1),
     2      spnts_inout2(ncols_inout2,npson2),
     3      sdad0(npdad,npdad),
     4      potdad(1000),chgs(1000),
     5      potson1(1000),potson2(1000),
     6      w1(1000),
     7      w11(1000),w12(1000),
     8      w21(1000),w2(1000),w22(1000)
        complex *16
     1      exchg11(ncols_outin1,ncols_inout1),
     2      exchg12(ncols_outin1,ncols_inout2),
     3      exchg22(ncols_outin2,ncols_inout2),
     4       exchg21(ncols_outin2,ncols_inout1),
     7      spnts1(npson1,npson1),spnts2(npson2,npson2),
     8      w4(1000),w5(1000),w6(1000),w7(1000),w8(1000),
     9      outin_spnts1(npson1,ncols_outin1),
     a      outin_spnts2(npson2,ncols_outin2),
     1      cd1,cd2
c 
c        construct the test set of charges
c 
        do 1200 i=1,npdad
c 
        chgs(i)=i
 1200 continue
c 
        call prin2('in exchg_test, chgs as created*',chgs,npdad*2)
c 
        call prinf('in exchg_test, npdad=*',npdad,1)
c 
c        apply the dad's scattering matrix to the charges
c 
        call cmatvec77(sdad0,npdad,npdad,chgs,potdad)
c 
cccc        call prin2('potdad as created*',potdad,npdad*2)
c 
c        apply the first sonny's scattering matrix to the charges
c 
        call cmatvec77(spnts_inout1,ncols_inout1,npson1,
     1      chgs,potson1)
  
        call cmatvec77(spnts1,npson1,npson1,
     1      chgs,w4)
  
cccc        call prin2('potson1 as created*',potson1,ncols_inout1*2)
cccc        call prin2('and w4 as created*',w4,npson1*2)
c 
c        apply the second sonny's scattering matrix to the charges
c 
        call cmatvec77(spnts_inout2,ncols_inout2,npson2,
     1      chgs(npson1+1),potson2)
c 
        call cmatvec77(spnts2,npson2,npson2,
     1      chgs(npson1+1),w5)
  
cccc        call prin2('potson2 as created*',potson2,ncols_inout2*2)
cccc        call prin2('and w5 as created*',w5,npson2*2)
c 
c       now, apply the exchange matrices to the sonnies' scattered
c       fields
c 
        call cmatvec77(exchg11,ncols_outin1,ncols_inout1,potson1,w11)
        call cmatvec77(exchg12,ncols_outin1,ncols_inout2,potson2,w12)
        call cmatvec77(exchg21,ncols_outin2,ncols_inout1,potson1,w21)
        call cmatvec77(exchg22,ncols_outin2,ncols_inout2,potson2,w22)
c 
        call cleaadd(w11,w12,w12,ncols_outin1)
        call cleaadd(w21,w22,w22,ncols_outin2)
c 
c        at this point, w12 is the incoming field on the first sonny,
c        and w22 is the incoming field on the second sonny
c 
        call cmatvec77(outin_spnts1,npson1,ncols_outin1,w12,w6)
  
cccc        call prin2('w6=*',w6,npson1*2)
  
        call cmatvec77(spnts1,npson1,npson1,chgs,w7)
        call cleaadd(w7,w6,w6,npson1)
c 
        call cmatvec77(outin_spnts2,npson2,ncols_outin2,w22,w7)
        call cmatvec77(spnts2,npson2,npson2,chgs(npson1+1),w8)
        call cleaadd(w8,w7,w7,npson2)
  
cccc        call prin2('and w7=*',w7,npson2*2)
c 
c        concatinate the arrays w6, w7, and subtract them
c        from potdad
c 
        do 2200 i=1,npson1
c 
        w8(i)=w6(i)
 2200 continue
c 
        do 2400 i=1,npson2
c 
        w8(i+npson1)=w7(i)
 2400 continue
c 
        call cleasubt(potdad,w8,w1,npdad)
  
cccc        call prin2('and the difference is*',w1,npdad*2)
  
  
        call cleascap(w1,w1,npdad,cd1)
        call cleascap(potdad,potdad,npdad,cd2)
c 
        d=cd1/cd2
        d=sqrt(d)
        call prin2('relative error in potdad*',d,1)
  
        return
        end
  
c 
c 
c 
c 
c 
        subroutine cmergetest2(spnts_inout_dad,npdad,
     1      ncols_inout,spnts_inout1,npson1,ncols_inout1,
     2      spnts_inout2,npson2,ncols_inout2,
     3      cmerge1,cmerge2)
        implicit real *8 (a-h,o-z)
        save
        complex *16
     1      spnts_inout1(ncols_inout1,npson1),
     2      spnts_inout2(ncols_inout2,npson2),
     3      spnts_inout_dad(ncols_inout,npdad),
     4      potdad(1000),chgs(1000),
     5      potson1(1000),potson2(1000),
     6      cmerge1(ncols_inout,ncols_inout1),
     7      cmerge2(ncols_inout,ncols_inout2),
     8      w1(1000),w2(1000),w3(1000)
c 
c        construct the test set of charges
c 
        do 1200 i=1,npdad
c 
        chgs(i)=i
 1200 continue
  
  
        call prin2('chgs as created*',chgs,npdad*2)
c 
c        apply the dad's scattering matrix to the charges
c 
        call cmatvec77(spnts_inout_dad,ncols_inout,npdad,
     1      chgs,potdad)
  
c 
c        apply the first sonny's scattering matrix to the charges
c 
        call cmatvec77(spnts_inout1,ncols_inout1,npson1,
     1      chgs,potson1)
  
  
        call prin2('potson1 as created*',potson1,ncols_inout1*2)
c 
c        apply the second sonny's scattering matrix to the charges
c 
        call cmatvec77(spnts_inout2,ncols_inout2,npson2,
     1      chgs(npson1+1),potson2)
  
  
        call prin2('potson2 as created*',potson2,ncols_inout2*2)
  
  
c 
c       now, apply the merging matrices to the sonnies' scattered
c       fields
c 
        call cmatvec77(cmerge1,ncols_inout,ncols_inout1,potson1,w1)
  
        call prin2('and w1=*',w1,ncols_inout*2)
  
        call cmatvec77(cmerge2,ncols_inout,ncols_inout2,potson2,w2)
  
        call cleaadd(w1,w2,w3,ncols_inout)
  
        call prin2('and w2=*',w2,ncols_inout*2)
  
        call prin2('and w3=*',w3,ncols_inout*2)
  
        call cleasubt(w3,potdad,w1,ncols_inout)
  
        call prin2('and w3-potdad=*',w1,ncols_inout*2)
  
  
  
        return
        end
c 
c 
c 
c 
c 
        subroutine scatmatr0dad(istart,np,iz,zs,amatr,
     1      expand,eval,ncols_inout,ncols_outin,s)
        implicit real *8 (a-h,o-z)
        save
        complex *16 cd,amatr(np,np),eval(np,ncols_outin),
     1      expand(ncols_inout,np),s(np,np)
        dimension iz(1),zs(2,1),work(500 000)
c 
c        construct the matrix of interactions between nodes
c        in this chunk
c 
        do 2400 i=1,np
        do 2200 j=1,np
c 
        iii=istart+i-1
        jjj=istart+j-1
cccc        call poteval(iz,zs,iii,jjj,cd)
        call poteval(iz,zs,jjj,iii,cd)
c 
        s(i,j)=cd
  
cccc        call prin2('cd=*',cd,2)
 2200 continue
 2400 continue
c 
c        invert the matrix of interactions
c 
        call corthom(s,np,work,cond)
  
  
        call csign_chg(s,s,np*np)
  
cccc        call prin2('and s=*',s,ncols_inout*ncols_outin*2)
cccc        call prin2('and s=*',s,100)
  
  
  
cccc        stop
  
        return
        end
c 
c 
c 
c 
c 
        subroutine poteval(iz,zs,i1,i2,pot)
        implicit real *8 (a-h,o-z)
        save
        complex *16 pot,rk,z,h0,h1
        dimension zs(2,1),iz(1)
c 
        rk=10
c 
        if(i1 .ne. i2) goto 1400
c 
        pot=-17
        pot=1
  
  
cccc        if(i1 .ge. 13) pot=1000000.0
cccc        if(i1 .le. 12) pot=-1000000.0
c 
        return
 1400 continue
  
        j1=iz(i1)
        j2=iz(i2)
  
  
  
cccc        j1=i1
cccc        j2=i2
  
  
  
  
c 
        d=(zs(1,j1)-zs(1,j2))**2+(zs(2,j1)-zs(2,j2))**2
        d=sqrt(d)
        z=d*rk
c 
  
cccc        call prinf('in poteval, j1=*',j1,1)
cccc        call prinf('in poteval, j2=*',j2,1)
cccc        call prin2('in poteval, d=*',d,1)
cccc        call prin2('in poteval, z=*',z,2)
  
        ifexpon=1
        call hank103(z,h0,h1,ifexpon)
c 
cccc        call prin2('in poteval, h0=*',h0,2)
  
  
  
        pot=h0
  
        pot=h0*zs(1,j1)
  
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine creaband3(rlength,hight,n,zs,w)
        implicit real *8 (a-h,o-z)
        save
        dimension zs(2,1),w(1)
c 
c 
c 
        n3=n/3
c 
        c1=0
        c2=0
        c3=0
  
        h=rlength/n3
c 
        do 1200 i=1,n3
c 
        zs(1,i)=-rlength/2+(i-1)*h+h/2
        zs(2,i)=cos(zs(1,i)*c1)*hight/2 *0.9
c 
        zs(1,i+n3)=-rlength/2+(i-1)*h+h/2
        zs(2,i+n3)=sin(zs(1,i)*c2)*hight*0.7 /2
c 
        zs(1,i+n3*2)=-rlength/2+(i-1)*h+h/2
        zs(2,i+n3*2)=-cos(zs(1,i)*c3)*hight *0.8 /2
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine creaband2(rlength,hight,n,zs,w)
        implicit real *8 (a-h,o-z)
        save
        dimension zs(2,1),w(1)
c 
c 
c 
        n3=n/3
c 
        h=rlength/n3
c 
        do 1200 i=1,n3
c 
        zs(1,i)=-rlength/2+(i-1)*h+h/2
        zs(2,i)=cos(zs(1,i)*10)*hight/2
c 
        zs(1,i+n3)=-rlength/2+(i-1)*h+h/2
        zs(2,i+n3)=sin(zs(1,i)*13)*hight*0.7 /2
c 
        zs(1,i+n3*2)=-rlength/2+(i-1)*h+h/2
        zs(2,i+n3*2)=-cos(zs(1,i)*7)*hight *0.8 /2
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine creaband(rlength,hight,n,zs,w)
        implicit real *8 (a-h,o-z)
        save
        dimension zs(2,1),w(1)
c 
c       create 2*n random numbers on the interval [-1,1]
c 
        call newrand(n*2,w)
c 
c       pack them things in the band
c 
        do 1200 i=1,n
c 
        zs(1,i)=w(i)*rlength-rlength/2
        zs(2,i)=w(i+n)*hight-hight/2
 1200 continue
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine sdad_direct(sdad0,ncols_inout,
     1      ncols_outin,evaldad,ncols_outin1,
     2      ncols_outin2,expanddad,
     4      expandson1,ncols_inout1,np1,expandson2,
     5      ncols_inout2,np2,evalson1,
     6      evalson2,
     7      b,c,sdadout,spnts_inout_dad)
c 
        implicit real *8 (a-h,o-z)
        save
        complex *16 sdad0(np1+np2,np1+np2),
     1      evaldad(ncols_outin1+ncols_outin2,ncols_outin),
     2      expanddad(ncols_inout,ncols_inout1+ncols_inout2),
     3      expandson1(ncols_inout1,np1),
     4      expandson2(ncols_inout2,np2),
     5      evalson1(np1,ncols_outin1),
     6      evalson2(np2,ncols_outin2),
     7      b(np1+np2,ncols_outin1+ncols_outin2),
     8      c(ncols_inout1+ncols_inout2,np1+np2),
     9      w1(200 000),w3(200 000),
     a      spnts_inout_dad(ncols_inout,np1+np2)
c 
c        form the matrix b of "diagonally concatinated"
c        matrices evalson1, evalson2
c 
        nbsum=ncols_outin1+ncols_outin2
        do 1600 i=1,nbsum
        do 1400 j=1,np1+np2
c 
        b(j,i)=0
 1400 continue
 1600 continue
c 
        do 2000 i=1,ncols_outin1
c 
        do 1800 j=1,np1
c 
        b(j,i)=evalson1(j,i)
 1800 continue
 2000 continue
  
  
c 
        do 2400 i=1,ncols_outin2
c 
        do 2200 j=1,np2
c 
        b(j+np1,i+ncols_outin1)=evalson2(j,i)
c 
 2200 continue
 2400 continue
  
        call prinf('np1=*',np1,1)
        call prinf('np2=*',np2,1)
  
cccc        call prin2('b as constructed*',b,nbsum*(np1+np2)*2)
c 
c 
c 
c        form the matrix cof "giagonally concatinated"
c        matrices expandson1, expandson2
c 
        ncsum=ncols_inout1+ncols_inout2
        do 3400 i=1,np1+np2
        do 3200 j=1,ncsum
c 
        c(j,i)=0
 3200 continue
 3400 continue
c 
        do 3800 i=1,np1
        do 3600 j=1,ncols_inout1
c 
        c(j,i)=expandson1(j,i)
 3600 continue
 3800 continue
c 
        do 4200 i=1,np2
        do 4000 j=1,ncols_inout2
c 
        c(j+ncols_inout1,i+np1)=expandson2(j,i)
 4000 continue
 4200 continue
c 
c       multiply them things together
c 
        call ncleamult(expanddad,c,w1,ncols_inout,
     1      ncsum,np1+np2)
c 
        call ncleamult(w1,sdad0,spnts_inout_dad,
     1      ncols_inout,np1+np2,np1+np2)
c 
        call ncleamult(spnts_inout_dad,b,w3,ncols_inout,
     1      np1+np2,nbsum)
c 
        call ncleamult(w3,evaldad,sdadout,ncols_inout,
     1      nbsum,ncols_outin)
c 
        return
        end
  
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code, and the beginning of
c        the actual scattering matrix merging code
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c         This file contains 2 user-callable subroutines: scatmatr0,
c         scattmatr_merge2. Following is a brief description of the
c         said two subroutines.
c 
c   scatmatr0 - constructs the scattering matrix (in its four
c         incarnations) for one chunk on level 0. It uses data
c         (hopefully) created by a prior call to the subroutines
c         bandgeo, skels_1chunk (see); it has no known uses as a
c         stand-alone device.
c 
c   scattmatr_merge2 - merges the user-supplied scattering matrices
c         s1, s2 of two brother-chunks, obtaining the scattering
c         matrix sout of the daddy chunk. It also produces the
c         splitting matrices split1, split2, "exchange matrices"
c         exchg11, exchg12, exchg21, exchg22, and the "merging
c         matrices" cmerge1, cmerge2.
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        subroutine scatmatr0(ichunks,ichunk,np,iz,zs,
     1      expand,eval,ncols_inout,ncols_outin,spnts,
     2      s,spnts_inout,outin_spnts,work,iwhich,cond)
        implicit real *8 (a-h,o-z)
        save
        complex *16 cd,spnts(np,np),eval(np,ncols_outin),
     1      expand(ncols_inout,np),
     2      s(ncols_inout,ncols_outin),
     3      spnts_inout(ncols_inout,np),
     4      outin_spnts(np,ncols_outin)
c 
        dimension iz(1),zs(2,1),ichunks(2,1),work(1),iwhich(3)
c 
c        This subroutine constructs the scattering matrix (in its
c        four incarnations) for one chunk on level 0. It uses data
c        (hopefully) created by a prior call to the subroutines
c        bandgeo, skels_1chunk (see); it has no known uses as a
c        stand-alone device.
c 
c                       Input parameters:
c 
c  ichunks - the integer (2,nchunks) array defining the subdivision
c        of the nodes in the structure into chunks. Must have been
c        produced by a prior call to the subroutine bandgeo.
c  ichunk - the sequence number of the chunk for which the scattering
c        matrix is to be created
c  np - the number of nodes inside the chunk number ichunk
c  iz - the integer array describing the transposition that converts
c        the array zs into the array zs2, as produced by a prior call
c        to the subroutine bandgeo.
c  zs - the nodes in the simulation
c  expand - the (ncols,np)-matrix converting a charge distribution
c        at the np nodes in the chunk number ichunk into an equivalent
c        charge distribution on the outgoing skeleton of size ncols;
c        must have been produced by a prior call to the subroutine
c        skels_1chunk (see)
c  eval - the (np,ncols_outin)-matrix converting the potential at
c        the ncols_outin "incoming skeleton" nodes inside the chunk
c        number ichunk into the potential at all np nodes in the
c        said chunk must have been produced by a prior call to the
c        subroutine skels_1chunk (see)
c  ncols_inout - the number of elements in the outgoing (inside to
c        outside) skeleton of the chunk being dealt with
c  ncols_outin - the number of elements in the incoming (outside to
c        inside) skeleton of the chunk being dealt with
c  iwhich - a three-element integer array telling the subroutine which
c        incarnations of the scattering matrix should be created.
c 
c     EXPLANATION: THE MATRIX SPNTS IS ALWAYS CREATED BY THIS SUBROUTINE.
C          THE MATRIX SPNTS_INOUT IS CREATED ONLY IF IWHICH(1) IS NOT
C                  EQUAL TO ZERO
C 
C          THE MATRIX S IS CREATED ONLY IF IWHICH(2) IS NOT EQUAL TO ZERO
C 
C          THE MATRIX OUTIN_SPNTS IS CREATED ONLY IF IWHICH(3) IS NOT
C          EQUAL TO ZERO
C 
c  spnts - the (np,np) matrix converting the incoming potential at the
c        np nodes constituting the chunk into the induced charges at
c        these same np nodes
c  s - the (ncols_inout,ncols_outin) matrix converting the incoming
c        potential tabulated at the ncols_outin "incoming skeleton"
c        nodes into the outgoing scattered potential, tabulated at the
c        ncols_inout "outgoing skeleton" nodes
c  spnts_inout - the (ncols_inout,np) matrix converting the incoming
c        potential tabulated at the np nodes constituting the chunk
c        into the outgoing scattered potential, tabulated at the
c        ncols_inout "outgoing skeleton" nodes
c  outin_spnts - the (np,ncols_outin) matrix converting the incoming
c        potential tabulated at the ncols_outin "incoming skeleton"
c        nodes into the induced charges at the np nodes constituting
c        the chunk
c  cond - the condition number of the linear system that had to be
c        solved to obtain the scattering matrix
c 
c                       Work arrays:
c 
c  work - must be at least 2 * n * (n+1)+2 complex *16 locations long
c 
c        . . . construct the matrix of interactions between nodes
c              in this chunk
c 
        istart=ichunks(1,ichunk)
        do 2400 i=1,np
        do 2200 j=1,np
c 
        iii=istart+i-1
        jjj=istart+j-1
        call poteval(iz,zs,jjj,iii,cd)
c 
        spnts(i,j)=cd
 2200 continue
 2400 continue
c 
c        invert the matrix of interactions
c 
        call corthom(spnts,np,work,cond)
c 
        call prin2('in scatmatr0 after corthom, cond=*',
     1      cond,1)
c 
        call csign_chg(spnts,spnts,np**2)
c 
c       multiply the matrix spnts by the matrices
c       expand, eval, to convert it into the scattering
c       matrix on the skeletons
c 
        if( (iwhich(1) .eq. 0) .and. (iwhich(2) .eq. 0)
     1      .and. (iwhich(3) .eq. 0) ) return
c 
        call ncleamult(expand,spnts,work,
     1      ncols_inout,np,np)
c 
        if(iwhich(1) .ne. 0)
     1      call cleascop(work,spnts_inout,np*ncols_inout)
c 
        if(iwhich(2) .ne. 0)
     1      call ncleamult(work,eval,s,ncols_inout,
     2      np,ncols_outin)
c 
        if(iwhich(3) .ne. 0)
     1      call ncleamult(spnts,eval,outin_spnts,np,
     2      np,ncols_outin)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine scattmatr_merge2(ier,expand,eval,
     1    iz,zs,
     2    irows_outin1,ncols_outin1,
     3    irows_outin2,ncols_outin2,
     4    icols_inout1,ncols_inout1,
     5    icols_inout2,ncols_inout2,
     6    irows_outin,ncols_outin,
     7    icols_inout,ncols_inout,s1,s2,sout,split1,split2,
     8      exchg11,exchg12,exchg22,exchg21,cmerge1,
     9      cmerge2,conds,w,lw,ltot,iflast)
c 
        implicit real *8 (a-h,o-z)
        save
        dimension zs(2,1),iz(1),
     1      irows_outin1(1),irows_outin2(1),icols_inout1(1),
     2      icols_inout(1),irows_outin(1),icols_inout2(1),
     3      conds(2)
c 
        complex *16
     1      expand(ncols_inout,ncols_inout1+ncols_inout2),
     2      eval(ncols_outin1+ncols_outin2,ncols_outin),
     3      s1(ncols_inout1,ncols_outin1),
     4      s2(ncols_inout2,ncols_outin2),sout(1),
     5      split1(ncols_outin1,ncols_outin),
     6      split2(ncols_outin2,ncols_outin),w(1)
c 
        complex *16
     1      exchg11(ncols_outin1,ncols_inout1),
     2      exchg12(ncols_outin1,ncols_inout2),
     3      exchg22(ncols_outin2,ncols_inout2),
     4      exchg21(ncols_outin2,ncols_inout1),
     5      cmerge1(ncols_inout,ncols_inout1),
     6      cmerge2(ncols_inout,ncols_inout2)
c 
c        This subroutine merges the user-supplied scattering
c        matrices s1, s2 of two brother-chunks, obtaining the
c        scattering matrix sout of the daddy chunk. It also
c        produces the splitting matrices split1, split2,
c        "exchange matrices" exchg11, exchg12, exchg21,
c        exchg22, and the "merging matrices" cmerge1,
c        cmerge2.
c 
c        In fact, this subroutine is purely a memory manager;
c        this subroutine calls the subroutine scattmatr_merge
c        (see) that does all the work.
c 
c 
c                      Input parameters:
c 
c  expand - the (ncols_inout,ncols_inout1+ncols_inout2)-matrix
c        converting the vector of concatinated far-field
c        expansions for the sonnies into the far-field expansion
c        for the daddy. Must have been produced by a prior call
c        the subroutine merges_1chunk (see)
c  eval - the (ncols_outin1+ncols_outin2,ncols_outin)-matrix
c        converting the incoming expansion for the daddy into
c        the vector of concatinated incoming expansions for the
c        sonnies. Must have been produced by a prior call
c        the subroutine merges_1chunk (see)
c  iz - the integer array describing the transposition that converts
c        the array zs into the array zs2, as produced by a prior call
c        to the subroutine bandgeo.
c  zs - the user-supplied array of all nodes in the simulation
c  irows_outin1 - coordinates (ncols_outin1 of them things) in
c        array iz of the elements of the inner skeleton for the
c        incoming potential for the first sonny
c  ncols_outin1 - the number of elements in the array irows_outin1
c        (inter alia)
c  irows_outin2 - coordinates (ncols_outin2 of them things) in
c        array iz of the elements of the inner skeleton for the
c        incoming potential for the second sonny
c  ncols_outin2 - the number of elements in the array irows_outin2
c        (inter alia)
c  icols_inout1 - coordinates (ncols_inout1 of them things) in array
c        iz of the elements of the inner skeleton for the outgoing
c        potential of the first sonny
c  ncols_inout1 - the number of elements in the array icols_inout1
c        (inter alia)
c  icols_inout2 - coordinates (ncols_inout1 of them things) in array
c        iz of the elements of the inner skeleton for the outgoing
c        potential of the second sonny
c  ncols_inout1 - the number of elements in the array icols_inout2
c        (inter alia)
c  irows_outin - coordinates (ncols_outin of them things) in
c        array iz of the elements of the inner skeleton for the
c        incoming potential for the daddy
c  ncols_outin1 - the number of elements in the array irows_outin
c        (inter alia)
c  icols_inout - coordinates (ncols_inout1 of them things) in array
c        iz of the elements of the inner skeleton for the outgoing
c        potential of the daddy
c  ncols_inout1 - the number of elements in the array icols_inout
c        (inter alia)
c  s1, s2 - the scattering matrices of the two sonnies (must have
c        been created via preceding calls to one of the subroutines
c        skels_1chunk, merges_1chunk
c  lw - the length of the user-supplied work array w (THE LENGTH
C        ESTIMATES ARE PROBLEM-DEPENDENT AND AT THIS TIME
C        IMPERFECTLY UNDERSTOOD)
c  iflast - integer parameter telling the subroutine whether the
c        daddy chunk is the last (biggest) chunk in the structure.
c        Obviously, the last chunk in the simulation has no
c        scattering matrix,n no incomingpotentials, and no skeletons
c        of any kind. Thus, the matrices sout, split1, split2,
c        cmerge1, cmerge2 should not be generated. The only ones
c        that SHOULD be generated are the exchange matrices for the
c        sonnies: exchg11, exchg12, exchg21, exchg22.
c   iflast = 0 will cause matrices to be generated
c   iflast = 1 will cause only exchg11, exchg12, exchg21, exchg22
c        to be generated
c 
c                      Output parameters:
c 
c  ier - error return code
c  sout - the scattering matrix for the daddy, obtained by merging
c        the scattering matrices for the sonnies
c  split1 - the splitting matrix, converting the incoming potential
c        on the daddy into the corresponding incoming potential
c        (including the son's interactions with the brother) on
c        the first sonny
c  split2 - the splitting matrix, converting the incoming potential
c        on the daddy into the corresponding incoming potential
c        (including the son's interactions with the brother) on
c        the second sonny
c  exchg11 - the exchange matrix, converting the outgoing (khi)
c        potential generated by the first sonny in the absence
c        of the second one into the incoming potential (eta) on
c        the (same) first sonny, produced by the interaction
c        between the brothers
c  exchg12 - the exchange matrix, converting the outgoing (khi)
c        potential generated by the second sonny in the absence
c        of the first one into the incoming potential (eta) on
c        the first sonny, produced by the interaction between the
c        brothers
c  exchg21 - the exchange matrix, converting the outgoing (khi)
c        potential generated by the first sonny in the absence of
c        the second one into the incoming potential (eta) on the
c        second sonny, produced by the interaction between the
c        brothers
c  exchg22 - the exchange matrix, converting the outgoing (khi)
c        potential generated by the second sonny in the absence
c        of the first one into the incoming potential (eta) on
c        the (same) second sonny
c  cmerge1 - the "merging" matrix, converting the outgoing (khi)
c        potential generated by the first sonny in the absence
c        of the second one onto the outgoing (khi) expansion for
c        the daddy
c  cmerge2 - the "merging" matrix, converting the outgoing (khi)
c        potential generated by the second sonny in the absence
c        of the second one onto the outgoing (khi) expansion for
c        the daddy
c  conds - a real array of length 2 containing (rather crude)
c        estimates of the condition numbers of the two linear
c        systems that were solved inside the subroutine
c        scattmatr_merge0, called by this subroutine
c  ltot - the amount of memory (in complex *16 words) in array w
c        actually used by the subroutine.
c 
c                         Work arrays:
c 
c  w - must be longish
c 
c 
c        . . . allocate memory for scattmatr_merge1
c 
        it1=1
        lt1=ncols_outin1*ncols_outin+100
c 
        it2=it1+lt1
        lt2=ncols_outin2*ncols_outin+100
c 
        it1bar=it2+lt2
        lt1bar=ncols_inout*ncols_inout1+100
c 
        it2bar=it1bar+lt1bar
        lt2bar=ncols_inout*ncols_inout2+100
c 
        it12=it2bar+lt2bar
        lt12=ncols_outin1*ncols_inout2+100
c 
        it21=it12+lt12
        lt21=ncols_outin2*ncols_inout1+100
c 
        iw=it21+lt21
c 
        lenw=lw-iw
c 
        if(lenw .lt. 50 000) then
            ier=32
            ltot=iw+50 000
            return
        endif
c 
        call scattmatr_merge(ier,expand,eval,w(it1),w(it2),
     1      w(it1bar),w(it2bar),w(it12),w(it21),iz,zs,
     2    irows_outin1,ncols_outin1,
     3    irows_outin2,ncols_outin2,
     4    icols_inout1,ncols_inout1,
     5    icols_inout2,ncols_inout2,
     6    irows_outin,ncols_outin,
     7    icols_inout,ncols_inout,s1,s2,sout,split1,split2,
     8      exchg11,exchg12,exchg22,exchg21,cmerge1,
     9      cmerge2,conds,w(iw),lw,ltot,iflast)
c 
        ltot=ltot+iw
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine scattmatr_merge(ier,expand,eval,t1,t2,
     1      t1bar,t2bar,t12,t21,iz,zs,
     2    irows_outin1,ncols_outin1,
     3    irows_outin2,ncols_outin2,
     4    icols_inout1,ncols_inout1,
     5    icols_inout2,ncols_inout2,
     6    irows_outin,ncols_outin,
     7    icols_inout,ncols_inout,s1,s2,sout,split1,split2,
     8      exchg11,exchg12,exchg22,exchg21,cmerge1,
     9      cmerge2,conds,w,lw,ltot,iflast)
c 
        implicit real *8 (a-h,o-z)
        save
        dimension zs(2,1),iz(1),
     1      irows_outin1(1),irows_outin2(1),icols_inout1(1),
     2      icols_inout(1),irows_outin(1),icols_inout2(1),
     3      conds(2)
c 
        complex *16
     1      expand(ncols_inout,ncols_inout1+ncols_inout2),
     2      eval(ncols_outin1+ncols_outin2,ncols_outin),
     3      t1(ncols_outin1,ncols_outin),
     4      t2(ncols_outin2,ncols_outin),
     5      t1bar(ncols_inout,ncols_inout1),
     6      t2bar(ncols_inout,ncols_inout2),
     7      t12(ncols_outin1,ncols_inout2),
     8      t21(ncols_outin2,ncols_inout1),cd,
     9      s1(ncols_inout1,ncols_outin1),
     a      s2(ncols_inout2,ncols_outin2),sout(1),
     1      split1(ncols_outin1,ncols_outin),
     2      split2(ncols_outin2,ncols_outin),w(1)
c 
        complex *16
     1      exchg11(ncols_outin1,ncols_inout1),
     2      exchg12(ncols_outin1,ncols_inout2),
     3      exchg22(ncols_outin2,ncols_inout2),
     4      exchg21(ncols_outin2,ncols_inout1),
     5      cmerge1(ncols_inout,ncols_inout1),
     6      cmerge2(ncols_inout,ncols_inout2)
c 
c        This subroutine merges the user-supplied scattering
c        matrices s1, s2 of two brother-chunks, obtaining the
c        scattering matrix sout of the daddy chunk. It also
c        produces the splitting matrices split1, split2,
c        "exchange matrices" exchg11, exchg12, exchg21,
c        exchg22, and the "merging matrices" cmerge1,
c        cmerge2.
c 
c                      Input parameters:
c 
c  expand - the (ncols_inout,ncols_inout1+ncols_inout2)-matrix
c        converting the vector of concatinated far-field
c        expansions for the sonnies into the far-field expansion
c        for the daddy. Must have been produced by a prior call
c        the subroutine merges_1chunk (see)
c  eval - the (ncols_outin1+ncols_outin2,ncols_outin)-matrix
c        converting the incoming expansion for the daddy into
c        the vector of concatinated incoming expansions for the
c        sonnies. Must have been produced by a prior call
c        the subroutine merges_1chunk (see)
c  iz - the integer array describing the transposition that converts
c        the array zs into the array zs2, as produced by a prior call
c        to the subroutine bandgeo.
c  zs - the user-supplied array of all nodes in the simulation
c  irows_outin1 - coordinates (ncols_outin1 of them things) in
c        array iz of the elements of the inner skeleton for the
c        incoming potential for the first sonny
c  ncols_outin1 - the number of elements in the array irows_outin1
c        (inter alia)
c  irows_outin2 - coordinates (ncols_outin2 of them things) in
c        array iz of the elements of the inner skeleton for the
c        incoming potential for the second sonny
c  ncols_outin2 - the number of elements in the array irows_outin2
c        (inter alia)
c  icols_inout1 - coordinates (ncols_inout1 of them things) in array
c        iz of the elements of the inner skeleton for the outgoing
c        potential of the first sonny
c  ncols_inout1 - the number of elements in the array icols_inout1
c        (inter alia)
c  icols_inout2 - coordinates (ncols_inout1 of them things) in array
c        iz of the elements of the inner skeleton for the outgoing
c        potential of the second sonny
c  ncols_inout1 - the number of elements in the array icols_inout2
c        (inter alia)
c  irows_outin - coordinates (ncols_outin of them things) in
c        array iz of the elements of the inner skeleton for the
c        incoming potential for the daddy
c  ncols_outin1 - the number of elements in the array irows_outin
c        (inter alia)
c  icols_inout - coordinates (ncols_inout1 of them things) in array
c        iz of the elements of the inner skeleton for the outgoing
c        potential of the daddy
c  ncols_inout1 - the number of elements in the array icols_inout
c        (inter alia)
c  s1, s2 - the scattering matrices of the two sonnies (must have
c        been created via preceding calls to one of the subroutines
c        skels_1chunk, merges_1chunk
c  lw - the length of the user-supplied work array w (THE LENGTH
C        ESTIMATES ARE PROBLEM-DEPENDENT AND AT THIS TIME
C        IMPERFECTLY UNDERSTOOD)
c 
c                      Output parameters:
c 
c 
c  ier - error return code
c  sout - the scattering matrix for the daddy, obtained by merging
c        the scattering matrices for the sonnies
c  split1 - the splitting matrix, converting the incoming potential
c        on the daddy into the corresponding incoming potential
c        (including the son's interactions with the brother) on
c        the first sonny
c  split2 - the splitting matrix, converting the incoming potential
c        on the daddy into the corresponding incoming potential
c        (including the son's interactions with the brother) on
c        the second sonny
c  exchg11 - the exchange matrix, converting the outgoing (khi)
c        potential generated by the first sonny in the absence
c        of the second one into the incoming potential (eta) on
c        the (same) first sonny, produced by the interaction
c        between the brothers
c  exchg12 - the exchange matrix, converting the outgoing (khi)
c        potential generated by the second sonny in the absence
c        of the first one into the incoming potential (eta) on
c        the first sonny, produced by the interaction between the
c        brothers
c  exchg21 - the exchange matrix, converting the outgoing (khi)
c        potential generated by the first sonny in the absence of
c        the second one into the incoming potential (eta) on the
c        second sonny, produced by the interaction between the
c        brothers
c  exchg22 - the exchange matrix, converting the outgoing (khi)
c        potential generated by the second sonny in the absence
c        of the first one into the incoming potential (eta) on
c        the (same) second sonny
c  cmerge1 - the "merging" matrix, converting the outgoing (khi)
c        potential generated by the first sonny in the absence
c        of the second one onto the outgoing (khi) expansion for
c        the daddy
c  cmerge2 - the "merging" matrix, converting the outgoing (khi)
c        potential generated by the second sonny in the absence
c        of the second one onto the outgoing (khi) expansion for
c        the daddy
c  conds - a real array of length 2 containing (rather crude)
c        estimates of the condition numbers of the two linear
c        systems that were solved inside the subroutine
c        scattmatr_merge0, called by this subroutine
c  ltot - the amount of memory (in complex *16 words) in array w
c        actually used by the subroutine.
c 
c                         Work arrays:
c 
c  t1 - must be at least ncols_outin1*ncols_outin
c        complex *16 elements
c  t2 - must be at least ncols_outin2*ncols_outin
c        complex *16 elements
c  t1bar - must be at least ncols_inout*ncols_inout1
c        complex *16 elements
c  t2bar - must be at least ncols_inout*ncols_inout2
c        complex *16 elements
c  t12 - must be at least ncols_outin1*ncols_inout2
c        complex *16 elements
c  t21 - must be at least ncols_outin2*ncols_inout1
c        complex *16 elements
c  w - must be longish
c 
c        . . . separate the matrices expand, eval into the
c              translation matrices t1, t2, t1bar,t2bar
c 
cccc        call prinf('in scattmatr_merge, iflast=*',iflast,1)
cccc        stop
  
        ier=0
c 
        if(iflast .ne. 0) goto 3100
c 
        do 1600 i=1,ncols_outin
        do 1400 j=1,ncols_outin1
c 
        t1(j,i)=eval(j,i)
 1400 continue
 1600 continue
c 
        do 2000 i=1,ncols_outin
        do 1800 j=1,ncols_outin2
c 
        t2(j,i)=eval(j+ncols_outin1,i)
 1800 continue
 2000 continue
c 
        do 2600 i=1,ncols_inout1
        do 2400 j=1,ncols_inout
c 
        t1bar(j,i)=expand(j,i)
 2400 continue
 2600 continue
c 
        do 3000 i=1,ncols_inout2
        do 2800 j=1,ncols_inout
c 
        t2bar(j,i)=expand(j,i+ncols_inout1)
 2800 continue
 3000 continue
c 
 3100 continue
c 
c        construct the matrices t12, t21
c 
        do 3400 i=1,ncols_inout1
        do 3200 j=1,ncols_outin2
c 
        iii=icols_inout1(i)
        jjj=irows_outin2(j)
        call poteval(iz,zs,iii,jjj,cd)
c 
        t21(j,i)=cd
 3200 continue
 3400 continue
c 
        do 3800 i=1,ncols_inout2
        do 3600 j=1,ncols_outin1
c 
        iii=icols_inout2(i)
        jjj=irows_outin1(j)
        call poteval(iz,zs,iii,jjj,cd)
c 
        t12(j,i)=cd
 3600 continue
 3800 continue
c 
c        construct the scattering and splitting matrices
c 
        iq1=1
        lq1=ncols_outin1**2+100
c 
        iq2=iq1+lq1
        lq2=ncols_outin2**2+100
c 
        it12s2=iq2+lq2
        lt12s2=ncols_outin1*ncols_outin2+100
c 
        it21s1=it12s2+lt12s2
        lt21s1=ncols_outin1*ncols_outin2+100
c 
        max_inout=ncols_inout
        if(ncols_inout1 .gt. max_inout) max_inout=ncols_inout1
        if(ncols_inout2 .gt. max_inout) max_inout=ncols_inout2
c 
        max_outin=ncols_outin
        if(ncols_outin1 .gt. max_outin) max_outin=ncols_outin1
        if(ncols_outin2 .gt. max_outin) max_outin=ncols_inout2
  
  
        iwork3=it21s1+lt21s1
cccc        lwork3=ncols_inout*ncols_outin+100
        lwork3=max_inout*max_outin+100
c 
        iwork4=iwork3+lwork3
cccc        lwork4=ncols_inout*ncols_outin+100
        lwork4=max_inout*max_outin+100
  
c 
        iq1_t12s2=iwork4+lwork4
        lq1_t12s2=ncols_outin1*ncols_outin2+100
c 
        iq2_t21s1=iq1_t12s2+lq1_t12s2
        lq2_t21s1=ncols_outin1*ncols_outin2+100
c 
        ltot=iq2_t21s1+lq2_t21s1
c 
        if(iflast .ne. 0) goto 4000
  
        it2bar_s2_q2=iq2_t21s1+lq2_t21s1
        lt2bar_s2_q2=ncols_inout*ncols_inout2+100
c 
        it1bar_s1_q1=it2bar_s2_q2+lt2bar_s2_q2
        lt1bar_s1_q1=ncols_inout*ncols_inout1+100
c 
        ltot=it1bar_s1_q1+lt1bar_s1_q1
c 
 4000 continue
  
        call prinf('in scattmatr_merge, ltot=*',ltot,1)
        call prinf('while lw=*',lw,1)
c 
        if(ltot .gt. lw) then
            ier=16
            return
        endif
c 
        call scattmatr_merge0(t1,t2,t1bar,t2bar,
     1      t12,t21,s1,s2,w(iq1),w(iq2),
     2      ncols_outin1,ncols_outin,ncols_outin2,
     3      ncols_inout,ncols_inout1,ncols_inout2,sout,
     4      split1,split2,w(it12s2),w(it21s1),
     5      w(iwork3),w(iwork4),conds,iflast)
c 
c        construct the merging and exchange matrices
c 
        call comb_exchg(t1,t2,t1bar,t2bar,
     1      t12,t21,s1,s2,w(iq1),w(iq2),
     2      ncols_outin1,ncols_outin,ncols_outin2,
     3      ncols_inout,ncols_inout1,ncols_inout2,
     4      w(it12s2),w(it21s1),
     5      exchg11,exchg12,exchg22,exchg21,
     6      cmerge1,cmerge2,
     7      w(iq1_t12s2),w(iq2_t21s1),w(it2bar_s2_q2),
     8      w(it1bar_s1_q1),w(iwork3),w(iwork4),iflast)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine comb_exchg(t1,t2,t1bar,t2bar,
     1      t12,t21,s1,s2,q1,q2,
     2      n_outin1,n_outin,n_outin2,
     3      n_inout,n_inout1,n_inout2,t12s2,t21s1,
     4      exchg11,exchg12,exchg22,exchg21,
     5      cmerge1,cmerge2,
     6      q1_t12s2,q2_t21s1,t2bar_s2_q2,
     7      t1bar_s1_q1,work4,work5,iflast)
c 
        implicit real *8 (a-h,o-z)
        save
        complex *16
     1      t1(n_outin1,n_outin),
     2      t2(n_outin2,n_outin),
     3      t1bar(n_inout,n_inout1),
     4      t2bar(n_inout,n_inout2),
     5      t12(n_outin1,n_inout2),
     6      t21(n_outin2,n_inout1),
     7      s1(n_inout1,n_outin1),
     8      s2(n_inout2,n_outin2),
     9      t12s2(n_outin1,n_outin2),
     a      t21s1(n_outin2,n_outin1),
c 
     1      exchg11(n_outin1,n_inout1),
     2      exchg12(n_outin1,n_inout2),
     3      exchg22(n_outin2,n_inout2),
     4      exchg21(n_outin2,n_inout1),
     5      cmerge1(n_inout,n_inout1),
     6      cmerge2(n_inout,n_inout2),
c 
     7      q1_t12s2(n_outin1,n_outin2),
     8      q2_t21s1(n_outin2,n_outin1),
     9      t1bar_s1_q1(n_inout,n_outin1),
     a      t2bar_s2_q2(n_inout,n_outin2),
     1      q1(n_outin1,n_outin1),
     2      q2(n_outin2,n_outin2),
     3      work4(1),work5(1)
c 
c        this subroutine constructs the "exchange matrices"
c        exchg11, exchg12, exchg21, exchg22, and the "merging
c        matrices" cmerge1, cmerge2.
c 
c                      Input parameters:
c 
c  t1,t2,t1bar,t2bar,t12,t21 - translation matrices; this notation
c        in agreement with the notation in the notes.
c  s1,s2 - the scattering matrices of the chunks to be merged
c  q1 - given by the formula q1=(I-t12*s2*t21*s1)^{-1}
c  q2 - given by the formula q2=(I-t21*s1*t12*s2)^{-1}
c  n_outin1 - the number of elements in the "incoming skeleton"
c        of the first of the chunks to be merged (the first sonny)
c  n_outin - the number of elements in the "incoming skeleton"
c        of the daddy
c  n_outin2 - the number of elements in the "incoming skeleton"
c        of the second of the chunks to be merged (the second sonny)
c  n_inout - the number of elements in the "outgoing skeleton"
c        of the daddy
c  n_inout1 - the number of elements in the "outging skeleton"
c        of the first sonny
c  n_inout2 - the number of elements in the "outging skeleton"
c        of the second sonny
c  t12s2 - given by the formula t12s2=t12*s2
c  t21s1 - given by the formula t21s1=t21*s1
c 
c                       Output parameters:
c 
c  exchg11 - the exchange matrix, converting the outgoing (khi)
c        potential generated by the first sonny in the absence
c        of the second one into the incoming potential (eta) on
c        the (same) first sonny, produced by the interaction
c        between the brothers
c  exchg12 - the exchange matrix, converting the outgoing (khi)
c        potential generated by the second sonny in the absence
c        of the first one into the incoming potential (eta) on
c        the first sonny, produced by the interaction between the
c        brothers
c  exchg21 - the exchange matrix, converting the outgoing (khi)
c        potential generated by the first sonny in the absence of
c        the second one into the incoming potential (eta) on the
c        second sonny, produced by the interaction between the
c        brothers
c  exchg22 - the exchange matrix, converting the outgoing (khi)
c        potential generated by the second sonny in the absence
c        of the first one into the incoming potential (eta) on
c        the (same) second sonny
c  cmerge1 - the "merging" matrix, converting the outgoing (khi)
c        potential generated by the first sonny in the absence
c        of the second one onto the outgoing (khi) expansion for
c        the daddy
c  cmerge2 - the "merging" matrix, converting the outgoing (khi)
c        potential generated by the second sonny in the absence
c        of the second one onto the outgoing (khi) expansion for
c        the daddy
c 
c                        Work arrays:
c 
c  q1_t12s2 - must be at least n_outin1*n_outin2 complex *16
c        elements long
c  q2_t21s1 - must be at least n_outin2*n_outin1 complex *16
c        elements long
c  t2bar_s2_q2 - must be at least n_inout*n_outin1 complex *16
c        elements long
c  t1bar_s1_q1 - must be at least n_inout,n_outin2 complex *16
c        elements long
c  work4 - must be at least n_inout*n_outin complex *16 elements
c  work5 - must be at least n_inout*n_outin complex *16 elements
c 
c 
c       . . . construct the exchange matrices
c 
c       q1(n_outin1,n_outin1) * t12s2(n_outin1,n_outin2 =
c       q1_t12s2(n_outin1,n_outin2)
c 
        call ncleamult(q1,t12s2,q1_t12s2,n_outin1,
     1      n_outin1,n_outin2)
c 
c    q1_t12s2(n_outin1,n_outin2) * t21(n_outin2,n_inout1)=
c    exchg11(n_outin1,n_inout1)
c 
        call ncleamult(q1_t12s2,t21,exchg11,n_outin1,
     1      n_outin2,n_inout1)
c 
c       q1(n_outin1,n_outin1) * t12(n_outin1,n_inout2) =
c       exchg12(n_outin1,n_inout2)
c 
        call ncleamult(q1,t12,exchg12,n_outin1,
     1      n_outin1,n_inout2)
c 
c       q2(n_outin2,n_outin2) * t21s1(n_outin2,n_outin1) =
c       q2_t21s1(n_outin2,n_outin1)
c 
        call ncleamult(q2,t21s1,q2_t21s1,n_outin2,
     1      n_outin2,n_outin1)
c 
c 
c       q2_t21s1(n_outin2,n_outin1) * t12(n_outin1,n_inout2) =
c       exchg22(n_outin2,n_inout2)
c 
        call ncleamult(q2_t21s1,t12,exchg22,n_outin2,
     1      n_outin1,n_inout2)
c 
c       q2(n_outin2,n_outin2) * t21(n_outin2,n_inout1)=
c       exchg21(n_outin2,n_inout1)
c 
        call ncleamult(q2,t21,exchg21,n_outin2,
     1      n_outin2,n_inout1)
c 
        if(iflast .ne. 0) return
c 
c        finally, construct the merging matrices
c 
c        . . . cmerge1
c 
c       s1(n_inout1,n_outin1) * q1_t12s2(n_outin1,n_outin2)=
c       work5(n_inout1,n_outin2)
c 
        call ncleamult(s1,q1_t12s2,work5,n_inout1,
     1      n_outin1,n_outin2)
c 
c       t1bar(n_inout,n_inout1) * work5(n_inout1,n_outin2)=
c       work4(n_inout,n_outin2)
c 
        call ncleamult(t1bar,work5,work4,n_inout,
     1      n_inout1,n_outin2)
c 
c        . . . at this point, work4= t1bar * s1 * q1* t12 * s2
c 
c       s2(n_inout2,n_outin2) * q2(n_outin2,n_outin2) =
c       work5(n_inout2,n_outin2)
c 
        call ncleamult(s2,q2,work5,n_inout2,
     1      n_outin2,n_outin2)
c 
c       t2bar(n_inout,n_inout2) * work5(n_inout2,n_outin2)=
c       t2bar_s2_q2(n_inout,n_outin2)
c 
        call ncleamult(t2bar,work5,t2bar_s2_q2,n_inout,
     1      n_inout2,n_outin2)
c 
c       t2bar_s2_q2(n_inout,n_outin2) + work4(n_inout,n_outin2)=
c       work5(n_inout,n_outin2)
c 
        call cleaadd(t2bar_s2_q2,work4,work5,
     1      n_inout*n_outin2)
c 
c       work5(n_inout,n_outin2) * t21(n_outin2,n_inout1)=
c       cmerge1(n_inout,n_inout1)
c 
        call ncleamult(work5,t21,cmerge1,n_inout,
     1      n_outin2,n_inout1)
c 
c       t1bar(n_inout,n_inout1) + cmerge1(n_inout,n_inout1)=
c       cmerge1(n_inout,n_inout1)=
c 
        call cleaadd(cmerge1,t1bar,cmerge1,
     1      n_inout*n_inout1)
c 
c        . . . cmerge2
c 
c 
c       s2(n_inout2,n_outin2) * q2_t21s1(n_outin2,n_outin1)=
c       work5(n_inout2,n_outin1)
c 
        call ncleamult(s2,q2_t21s1,work5,n_inout2,
     1      n_outin2,n_outin1)
c 
c       t2bar(n_inout,n_inout2) * work5(n_inout2,n_outin1)=
c       work4(n_inout,n_outin1)
c 
        call ncleamult(t2bar,work5,work4,n_inout,
     1      n_inout2,n_outin1)
c 
c        . . . at this point, work4= t2bar * s2 * q2* t21 * s1
c 
c       s1(n_inout1,n_outin1) * q1(n_outin1,n_outin1) =
c       work5(n_inout1,n_outin1)
c 
        call ncleamult(s1,q1,work5,n_inout1,
     1      n_outin1,n_outin1)
c 
c       t1bar(n_inout,n_inout1) * work5(n_inout1,n_outin1)=
c       t1bar_s1_q1(n_inout,n_outin1)
c 
        call ncleamult(t1bar,work5,t1bar_s1_q1,n_inout,
     1      n_inout1,n_outin1)
c 
c       t1bar_s1_q1(n_inout,n_outin1) + work4(n_inout,n_outin1)=
c       work5(n_inout,n_outin1)
c 
        call cleaadd(t1bar_s1_q1,work4,work5,
     1      n_inout*n_outin1)
c 
c       work5(n_inout,n_outin1) * t12(n_outin1,n_inout2)=
c       cmerge2(n_inout,n_inout2)
c 
        call ncleamult(work5,t12,cmerge2,n_inout,
     1      n_outin1,n_inout2)
c 
c       t2bar(n_inout,n_inout2) + cmerge2(n_inout,n_inout2)=
c       cmerge2(n_inout,n_inout2)=
c 
        call cleaadd(cmerge2,t2bar,cmerge2,
     1      n_inout*n_inout2)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine scattmatr_merge0(t1,t2,t1bar,t2bar,
     1      t12,t21,s1,s2,q1,q2,
     2      ncols_outin1,ncols_outin,ncols_outin2,
     3      ncols_inout,ncols_inout1,ncols_inout2,sout,
     4      split1,split2,t12s2,t21s1,work3,work4,conds,
     5      iflast)
        implicit real *8 (a-h,o-z)
        save
        complex *16
     1      t1(ncols_outin1,ncols_outin),
     2      t2(ncols_outin2,ncols_outin),
     3      t1bar(ncols_inout,ncols_inout1),
     4      t2bar(ncols_inout,ncols_inout2),
     5      t12(ncols_outin1,ncols_inout2),
     6      t21(ncols_outin2,ncols_inout1),
     7      s1(ncols_inout1,ncols_outin1),
     8      s2(ncols_inout2,ncols_outin2),
     9      q1(ncols_outin1,ncols_outin1),
     a      q2(ncols_outin2,ncols_outin2),
     1      work3(1),work4(1),sout(1),
     2      t12s2(ncols_outin1,ncols_outin2),
     3      t21s1(ncols_outin2,ncols_outin1),
     4      split1(ncols_outin1,ncols_outin),
     5      split2(ncols_outin2,ncols_outin)
c 
        dimension conds(2)
c 
c        This subroutine merges the user-supplied scattering
c        matrices s1, s2 of two brother-chunks, obtaining the
c        scattering matrix sout of the daddy chunk. It also
c        produces the splitting matrices split1, split2. Finally,
c        it returns the matrices q1, q2, to be used by the
c        subroutine comb_exchg, to be used in the construction
c        of the merging and exchange matrices.
c 
c                      Input parameters:
c 
c  t1,t2,t1bar,t2bar,t12,t21 - translation matrices; this notation
c        in agreement with the notation in the notes.
c  s1,s2 - the scattering matrices of the chunks to be merged
c  ncols_outin1 - the number of elements in the "incoming skeleton"
c        of the first of the chunks to be merged (the first sonny)
c  ncols_outin - the number of elements in the "incoming skeleton"
c        of the daddy
c  ncols_outin2 - the number of elements in the "incoming skeleton"
c        of the second of the chunks to be merged (the second sonny)
c  ncols_inout - the number of elements in the "outgoing skeleton"
c        of the daddy
c  ncols_inout1 - the number of elements in the "outging skeleton"
c        of the first sonny
c  ncols_inout2 - the number of elements in the "outging skeleton"
c        of the second sonny
c 
c                      Output parameters:
c 
c  q1 - given by the formula q1=(I-t12*s2*t21*s1)^{-1}
c  q2 - given by the formula q2=(I-t21*s1*t12*s2)^{-1}
c  sout - the scattering matrix for the daddy, obtained by merging
c        the scattering matrices for the sonnies
c  split1 - the splitting matrix, converting the incoming potential
c        on the daddy into the corresponding incoming potential
c        (including the son's interactions with the brother) on
c        the first sonny
c  split2 - the splitting matrix, converting the incoming potential
c        on the daddy into the corresponding incoming potential
c        (including the son's interactions with the brother) on
c        the second sonny
c  t12s2 - given by the formula t12s2=t12*s2
c  t21s1 - given by the formula t21s1=t21*s1
c  conds - a real array of length 2 containing (rather crude)
c        estimates of the condition numbers of the two linear
c        systems that were solved inside the subroutine
c        scattmatr_merge0, called by this subroutine
c 
c 
c                      Work arrays:
c 
c  work3 - must be at least max (ncols-inout**2,ncols_outin**2)
c          complex *16 elements long
c  work4 - must be at least max (ncols-inout**2,ncols_outin**2)
c          complex *16 elements long
c 
c 
c       . . . evaluate the first inverse (q1) given by the formula
c 
c           q1=(I-t12*s2*t21*s1)^{-1}
c 
        call ncleamult(t12,s2,t12s2,ncols_outin1,
     1      ncols_inout2,ncols_outin2)
c 
        call ncleamult(t21,s1,t21s1,ncols_outin2,
     1      ncols_inout1,ncols_outin1)
c 
        call ncleamult(t12s2,t21s1,q1,ncols_outin1,
     1      ncols_outin2,ncols_outin1)
c 
        call csign_chg(q1,q1,ncols_outin1**2)
  
        do 1400 i=1,ncols_outin1
c 
        q1(i,i)=1+q1(i,i)
 1400 continue
c 
        call corthom(q1,ncols_outin1,work3,cond)
c 
        call prin2(
     1  'in scattmatr_merge0 after first corthom, cond-1=*',
     2      cond-1,1)
c 
        conds(1)=cond
c 
c       evaluate the second inverse (q2) given by the formula
c 
c            q2=(I-t21*s1*t12*s2)^{-1}
c 
        call ncleamult(t21s1,t12s2,q2,ncols_outin2,
     1      ncols_outin1,ncols_outin2)
c 
        call csign_chg(q2,q2,ncols_outin2**2)
c 
        do 1800 i=1,ncols_outin2
        q2(i,i)=1+q2(i,i)
 1800 continue
c 
cccc        call prin2('before corthom, q2=*',q2,ncols_outin**2*2)
  
        call corthom(q2,ncols_outin2,work3,cond)
c 
        call prin2(
     1   'in scattmatr_merge0 after second corthom, cond-1=*',
     2      cond-1,1)
c 
        conds(2)=cond
c 
        if(iflast .ne. 0) return
c 
c        calculate work4 = t1 + t12 * s2 * t2
c 
        call ncleamult(t12s2,t2,work3,ncols_outin1,
     1      ncols_outin2,ncols_outin)
c 
        call cleaadd(t1,work3,work4,ncols_outin1*ncols_outin)
c 
c       calculate the product
c 
c       t1bar * s1 * (I-t12*s2*t21*s1)^{-1} * (t1 + t12 * s2 * t2)
c 
        call ncleamult(q1,work4,split1,ncols_outin1,
     1      ncols_outin1,ncols_outin)
c 
        call ncleamult(s1,split1,work4,ncols_inout1,
     1      ncols_outin1,ncols_outin)
c 
        call ncleamult(t1bar,work4,sout,ncols_inout,
     1      ncols_inout1,ncols_outin)
c 
c       . . . now, sout is the first element in the big sum;
c             construct the second element of the said sum
c 
c 
c        calculate work4 = t2 + t21 * s1 * t1
c 
        call ncleamult(t21s1,t1,work3,ncols_outin2,
     1      ncols_outin1,ncols_outin)
c 
        call cleaadd(t2,work3,work4,ncols_outin2*ncols_outin)
c 
c       calculate the product
c 
c       t2bar * s2 * (I-t21*s1*t12*s2)^{-1} * (t2 + t21 * s1 * t1)
c 
        call ncleamult(q2,work4,split2,ncols_outin2,
     1      ncols_outin2,ncols_outin)
c 
        call ncleamult(s2,split2,work4,ncols_inout2,
     1      ncols_outin2,ncols_outin)
c 
        call ncleamult(t2bar,work4,work3,ncols_inout,
     1      ncols_inout2,ncols_outin)
c 
        call cleaadd(sout,work3,sout,ncols_outin*ncols_inout)
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine csign_chg(a,b,n)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(1),b(1)
c 
        do 1200 i=1,n
c 
        b(i)=-a(i)
 1200 continue
c 
        return
        end
  
  
