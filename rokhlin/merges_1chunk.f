        implicit real *8 (a-h,o-z)
        dimension w(2000 000)
  
        dimension zs(2,100 000),zs2(2,100 000),
     1      iz(100 000),izinv(100000),thresh(2,10000),
     2      ichunks(2,10000),bord(2,10 000),thresh2(2,10000),
     3      izskel(100 000),zs3(2,10000),ichunks2(2,10000)
c 
        dimension expand_inout2(100 000),irows_outin2(100 000),
     1      icols_inout2(100 000),eval_outin2(100 000)
  
c 
        dimension expand_inout1(100 000),irows_outin1(100 000),
     1      icols_inout1(100 000),eval_outin1(100 000),
     2      ain(1000 000)
c 
        dimension expand(200 000),icols(1000),work1(10 000),
     1      work2(10 000),work3(10 000),eval(200 000),
     2      irows(100 000)
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
        eps=1.0d-12
  
        lenw=2000 000
  
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
  
        iw=22
        call bandplot(iw,zs2,n,thresh,nchunks,hight,bord,ii)
  
        lbord=ii
  
        iw=23
        call zquaplot2(iw,zs2(1,i),np,1,bord,ii,3,
     1      'one chunk*')
c 
c       construct the skeleton of one chunk
c 
        npts=300
        iplot=33
c 
cccc        goto 5000
  
        call skels_1chunk(ier,n,eps,zs,
     1      nchunks,ichunk1,npts,bord,lbord,iplot,
     2      zs2,iz,izinv,thresh,ichunks,w,lenw,lused,
     3      eval_outin1,irows_outin1,ncols_outin1,
     4      expand_inout1,icols_inout1,ncols_inout1)
  
        np1=ichunks(2,ichunk1)-ichunks(1,ichunk1)+1
c 
        call masstest(expand_inout1,icols_inout1,np1,
     1      ncols_inout1,iz,zs,ichunks,ichunk1,n,
     2      eval_outin1,irows_outin1,ncols_outin1)
  
  
cccc           stop
 4200 continue
        ichunk2=4
        iplot=34
c 
        call skels_1chunk(ier,n,eps,zs,
     1      nchunks,ichunk2,npts,bord,lbord,iplot,
     2      zs2,iz,izinv,thresh,ichunks,w,lenw,lused,
     3      eval_outin2,irows_outin2,ncols_outin2,
     4      expand_inout2,icols_inout2,ncols_inout2)
  
        np2=ichunks(2,ichunk2)-ichunks(1,ichunk2)+1
  
c 
        call masstest(expand_inout2,icols_inout2,np2,
     1      ncols_inout2,iz,zs,ichunks,ichunk2,n,
     2      eval_outin2,irows_outin2,ncols_outin2)
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
        iflev0=0
        ichunk2=2
c 
c        merge the skeletons of the two chunks, obtaining a
c        chunk on the next level
c 
        call merges_1chunk(ier,zs,zs2,iz,izinv,n,npts,
     1    eps,barrier,w,lenw,lused,
     2      thresh2,ichunks2,ichunk2,nchunks/2,
     2    irows_outin1,ncols_outin1,
     3    irows_outin2,ncols_outin2,
     4    icols_inout1,ncols_inout1,
     5    icols_inout2,ncols_inout2,
     6    eval,irows,nrows,expand,icols,ncols)
  
        call prinf('after merges_1chunk, lused=*',lused,1)
  
cccc        stop
  
        ii=598
c 
        call lolo_test2(ii,expand,icols,ncols,iz,zs,
     1      icols_inout1,ncols_inout1,icols_inout2,
     2      ncols_inout2,err)
  
  
  
  
        iplot=47
        call plot_inner_skeleton(zs2,irows,nrows,
     1      bord,lbord,iplot,zs,iz)
  
  
  
        iplot=48
        call plot_inner_skeleton(zs2,icols,ncols,
     1      bord,lbord,iplot,zs,iz)
  
  
cccc        stop
  
  
        ii=2980
        ii=598
  
        call lolo_test3(ii,eval,irows,nrows,iz,zs,
     1      irows_outin1,ncols_outin1,irows_outin2,
     2      ncols_outin2,err)
  
  
        call prinf('and nrows=*',nrows,1)
        call prin2('and barrier=*',barrier,1)
  
  
  
  
        stop
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
c 
c 
c 
c 
        subroutine bandplot(iw,zs,n,thresh,nchunks,
     1      hight,bord,ii)
        implicit real *8 (a-h,o-z)
        save
        dimension zs(2,1),thresh(2,1),bord(2,10000)
c 
c        find the horizontal boundaries of the band
c 
        bordup=-1.0d40
        borddown=1.0d40
c 
        do 1100 i=1,n
c 
        if(bordup .lt. zs(2,i) ) bordup=zs(2,i)
        if(borddown .gt. zs(2,i) ) borddown=zs(2,i)
 1100 continue
c 
        hight=bordup-borddown
cccc        call prin2('bordup=*',bordup,1)
cccc        call prin2('borddown=*',borddown,1)
cccc        call prin2('hight=*',hight,1)
cccc        call prinf('n=*',n,1)
  
        bordup=bordup+hight/2
        borddown=borddown-hight/2
  
cccc        stop
c 
c        construct the border
c 
        h=thresh(1,2)-thresh(1,1)
        bord(2,1)=borddown
        bord(1,1)=thresh(1,1)
        ii=1
        do 1200 i=1,nchunks/2
c 
        ii=ii+1
        bord(1,ii)=bord(1,ii-1)
        bord(2,ii)=bordup
c 
        ii=ii+1
        bord(2,ii)=bordup
        bord(1,ii)=bord(1,ii-1)+h
c 
        ii=ii+1
        bord(2,ii)=borddown
        bord(1,ii)=bord(1,ii-1)
c 
        ii=ii+1
        bord(2,ii)=borddown
        bord(1,ii)=bord(1,ii-1)+h
c 
 1200 continue
c 
        ii=ii+1
        bord(2,ii)=bordup
        bord(1,ii)=bord(1,ii-1)
c 
        ii=ii+1
        bord(2,ii)=bordup
        bord(1,ii)=thresh(1,1)
c 
        ii=ii+1
        bord(2,ii)=borddown
        bord(1,ii)=thresh(1,1)
c 
        ii=ii+1
        bord(2,ii)=borddown
        bord(1,ii)=thresh(1,nchunks)
c 
c       plot the border together with the nodes
c 
        if(iw .eq. 0) return
c 
        call zquaplot2(iw,bord,ii,3,zs,n,1,
     1      'chunks with nodes in them*')
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine lolo_test3(ii,eval,irows,ncols,iz,zs,
     1      irows_inout1,ncols_inout1,irows_inout2,
     2      ncols_inout2,err)
        implicit real *8 (a-h,o-z)
        save
        complex *16 eval(1),vals(10 000),vals2(10 000),
     1      vals3(10 000),val,val2
c 
        dimension iz(1),irows(ncols),zs(2,1),irows0(10000),
     1      irows_inout1(1),irows_inout2(1)
  
  
        call prinf('in lolo_test3, ii=*',ii,1)
  
c 
c       concatinate the two inner skeletons
c 
        do 1100 i=1,ncols_inout1
c 
        irows0(i)=irows_inout1(i)
 1100 continue
c 
        do 1150 i=1,ncols_inout2
c 
        irows0(i+ncols_inout1)=irows_inout2(i)
 1150 continue
  
        ncols0=ncols_inout2+ncols_inout1
  
  
        call prinf('irows0=*',irows0,ncols0)
  
        call chgexpand_in(ii,iz,zs,irows,ncols,vals)
  
  
        call prin2('after chgexpand, vals=*',vals,ncols*2)
  
  
        call cmatvec77(eval,ncols0,ncols,vals,vals2)
  
        call prin2('and vals2=*',vals2,ncols0*2)
  
  
        call chgexpand_in(ii,iz,zs,irows0,ncols0,vals3)
  
        call prin2('and vals3=*',vals3,ncols0*2)
  
        call cleasubt(vals2,vals3,vals3,ncols0)
  
        call prin2('in lolo_test3, the difference=*',
     1      vals3,ncols0*2)
  
  
ccc        stop
  
  
  
  
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine lolo_test2(ii,expand,icols,ncols,iz,zs,
     1      icols_inout1,ncols_inout1,icols_inout2,
     2      ncols_inout2,err)
        implicit real *8 (a-h,o-z)
        save
        complex *16 expand(ncols,1),vals(10 000),
     1      val,val2,charges(10000)
c 
        dimension iz(1),icols(ncols),zs(2,1),icols0(10000),
     1      icols_inout1(1),icols_inout2(1)
  
c 
c       concatinate the two inner skeletons
c 
        do 1100 i=1,ncols_inout1
c 
        icols0(i)=icols_inout1(i)
 1100 continue
c 
        do 1150 i=1,ncols_inout2
c 
        icols0(i+ncols_inout1)=icols_inout2(i)
 1150 continue
  
        ncols0=ncols_inout2+ncols_inout1
c 
c        created the ncols0 test charges
c 
        done=1
  
        do 1400 i=1,ncols0
c 
        charges(i)=i*done/ncols0
 1400 continue
  
cccc        call prin2('charges as created*',charges,ncols0*2)
c 
c       expand the potential of all ncols0 charges in the chunk
c       into the potential of the ncols "significant" charges
c 
        call chunkexpand_out(expand,ncols0,ncols,
     1      charges,vals)
  
cccc        call prin2('vals=*',vals,ncols*2)
c 
c       . . . evaluate the expansion at the node number ii
c 
        call expanout_eval(iz,zs,vals,icols,ncols,ii,val)
c 
c        evaluate the potential directly
c 
        call expanout_eval(iz,zs,charges,
     1      icols0,ncols0,ii,val2)
  
  
        call prin2('in lolo_test2, val2=*',val2,2)
        call prin2('in lolo_test2, val2-val=*',val2-val,2)
  
        err=abs(val2-val)
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine plot_inner_skeleton(zs2,irows,ncols,
     1      bord,lbord,iplot,zs,iz)
        implicit real *8 (a-h,o-z)
        save
        dimension zs2(2,1),zs(2,1),iz(1),bord(2,1),
     1      zsout(2,10000),irows(1)
c 
        do 1400 i=1,ncols
c 
        j=irows(i)
        j=iz(j)
c 
        zsout(1,i)=zs(1,j)
        zsout(2,i)=zs(2,j)
 1400 continue
c 
        call zquaplot2(iplot,bord,lbord,3,zsout,ncols,1,
     1      'one chunk, together with its inner skeleton*')
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine masstest(expand_inout,icols_inout,np,
     1      ncols_inout,iz,zs,ichunks,ichunk,n,
     2      eval_outin,irows_outin,ncols_outin)
        implicit real *8 (a-h,o-z)
        save
        complex *16 expand_inout(1),eval_outin(1)
        dimension ichunks(2,1),iz(1),zs(2,1),icols_inout(1),
     1      errs_inout(10 000),irows_outin(1),errs_outin(10 000)
c 
c       test the expand_inout matrix
c 
        imin=ichunks(1,ichunk)
        imax=ichunks(2,ichunk)
  
        call prinf('imin=*',imin,1)
        call prinf('imax=*',imax,1)
        call prinf('ncols_inout*',ncols_inout,1)
cccc        call prinf('nzskel=*',nzskel,1)
  
cccc        stop
c 
c        test the matrix expand_inout
c 
        do 1600 i=1,n
c 
        if( (i .ge. imin) .and. (i. le. imax) ) goto 1600
c 
        ii=i
  
        call expand_inout_test(ii,expand_inout,icols_inout,np,
     1      ncols_inout,iz,zs,
     2      ichunks,ichunk,errs_inout(i) )
c 
 1600 continue
c 
c        test the matrix eval_outin
c 
        do 2600 i=1,n
c 
cccc        call prinf('in masstest, i=*',i,1)
  
        if( (i .ge. imin) .and. (i. le. imax) ) goto 2600
c 
        ii=i
c 
        call eval_outin_test(ii,eval_outin,irows_outin,np,
     1      ncols_outin,iz,zs,ichunks,ichunk,errs_outin(i) )
c 
 2600 continue
c 
        call prin2('and errs_inout=*',errs_inout,n)
        call prin2('and errs_outin=*',errs_outin,n)
c 
        dinout=0
        doutin=0
        do 2800 i=1,n
c 
        dinout=dinout+errs_inout(i)
        doutin=doutin+errs_outin(i)
 2800 continue
c 
        dinout=dinout/n
        doutin=doutin/n
  
  
        call prin2('and average errs_inout=*',dinout,1)
        call prin2('and average errs_outin=*',doutin,1)
  
  
        return
        end
c 
c 
c 
c 
c 
        subroutine plot_1chunk(zs2,ichunks,ichunk,bord,lbord,
     1      iplot,izskel,nzskel,zs,iz,zsout)
        implicit real *8 (a-h,o-z)
        save
        dimension zs2(2,1),zs(2,1),iz(1),bord(2,1),ichunks(2,1),
     1      izskel(1),zsout(2,10000)
c 
        np=ichunks(2,ichunk)-ichunks(1,ichunk)+1
        izstart=ichunks(1,ichunk)
c 
        do 1400 i=1,nzskel
c 
        j=izskel(i)
        j=iz(j)
c 
        zsout(1,i)=zs(1,j)
        zsout(2,i)=zs(2,j)
 1400 continue
c 
        call zquaplot3(iplot,bord,lbord,3,zsout,nzskel,1,
     1      zs2(1,izstart),np,1,
     1      'one chunk, together with its outer skeleton*')
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine eval_outin_test(ii,eval,irows,np,
     1      ncols,iz,zs,
     2      ichunks,ichunk,err)
        implicit real *8 (a-h,o-z)
        save
        complex *16 eval(np,ncols),vals(10 000),cd,
     1      valsout(10 000),valsout2(10 000)
c 
        dimension iz(1),irows(ncols),zs(2,1),
     1      ichunks(2,1)
c 
c       . . . construct the expansion
c 
  
        izstart=ichunks(1,ichunk)
        call chgexpand_in(ii,iz,zs,irows,ncols,vals)
c 
c       . . . evaluate the expansion
c 
        call cmatvec77(eval,np,ncols,vals,valsout)
c 
c        . . . evaluate the potential directly
c 
        jjj=ii
  
cccc        call prin2('test source is*',zs(1,iz(jjj)),2)
        do 2200 i=1,np
c 
        j=i+izstart-1
cccc        j=i
  
c 
        call poteval(iz,zs,jjj,j,valsout2(i))
 2200 continue
c 
        call cleasubt(valsout,valsout2,valsout,np)
        call cleascap(valsout,valsout,np,cd)
        d=cd
        d=sqrt(d)
c 
  
        call cleascap(valsout2,valsout2,np,cd)
        d2=cd
        d2=sqrt(d2)
c 
        err=d/d2
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine expand_inout_test(ii,expand,icols,np,
     1      ncols,iz,zs,
     2      ichunks,ichunk,err)
        implicit real *8 (a-h,o-z)
        save
        complex *16 expand(np,ncols),vals(10 000),cd,
     1      val,val2,charges(10000)
c 
        dimension iz(1),icols(ncols),zs(2,1),ichunks(2,1)
  
c 
c        created the np test charges
c 
        done=1
        izstart=ichunks(1,ichunk)
  
        do 1400 i=1,np
c 
        charges(i)=i*done/np
 1400 continue
c 
c       expand the potential of all np charges in the chunk
c       into the potential of the ncols "significant" charges
c 
        call chunkexpand_out(expand,np,ncols,
     1      charges,vals)
c 
c       . . . evaluate the expansion at the node number ii
c 
        call expanout_eval(iz,zs,vals,icols,ncols,ii,val)
c 
c        evaluate the potential directly
c 
        val2=0
        do 2200 i=1,np
c 
        j=i+izstart-1
c 
        call poteval(iz,zs,j,ii,cd)
c 
        val2=val2+cd*charges(i)
 2200 continue
  
cccc        call prin2('in expand_inout_test, val2=*',val2,2)
cccc        call prin2('in expand_inout_test, val2-val=*',val2-val,2)
  
        err=abs(val2-val)
c 
        return
        end
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code and the beginning of
c        the skeletonization code proper.
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c        This file contains three user-callable subroutines:
c        bandgeo, skels_1chunk, merges_1chunk. Following are
c        brief descriptions of these three subroutines.
c 
c 
c   bandgeo - performs preliminary geometry processing for
c        the direct solver in a band. Starting with a
c        user-specified collection of n nodes zs, it sorts them
c        in order of increasing abscissae, subdivides the band
c        into chunks on the finest level (the number of chunks
c        to be created is user-supplied, as well as the ends
c        of the band), and allocates the points to appropriate
c        chunks. It also creates the auxiliary control structures
c        iz, izinv, thresh, ichunks (see below). PLEASE NOTE THAT
c        THE BAND IS ASSUMED TO BE HORIZONTAL.
c 
c   skels_1chunk - constructs inner skeletons (for both incoming
c        and outgoing fields) for a single chunk. It uses data
c        constructed via a preceding call to the subroutine
c        bandgeo (see above)
c 
c   merges_1chunk - merges the skeletons of two chunks (both
c        incoming and outgoing skeletons are processed by this
c        subroutine), producing the daddy's skeletons. It also
c        produces the matrix eval, converting incoming expansion
c        on the skeleton of the daddy into equivalent incoming
c        expansions on the skeletons of the two sonnies, and the
c        matrix expand, converting outgoing expansions on the
c        skeletons of the sonnies into an equivalent outgoing
c        expansion on the daddy's skeleton
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        subroutine merges_1chunk(ier,zs,zs2,iz,izinv,
     1    n,npts,eps,barrier,w,lenw,lused,
     2    thresh2,ichunks2,ichunk2,nchunks2,
     3    irows_outin1,ncols_outin1,
     4    irows_outin2,ncols_outin2,
     5    icols_inout1,ncols_inout1,
     6    icols_inout2,ncols_inout2,
     7    eval,irows,nrows,expand,icols,ncols)
c 
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),zs(2,1),zs2(2,1),iz(1),izinv(1),
     1      ichunks2(2,1),thresh2(2,1),
     2      irows_outin1(1),irows_outin2(1),icols_inout1(1),
     3      icols_inout2(1),
     4      eval(1),irows(1),expand(1),icols(1)
c 
c        This subroutine merges the skeletons of two chunks
c        (both incoming and outgoing skeletons are processed
c        by this subroutine), producing the daddy's skeletons.
c        It also produces the matrix eval, converting incoming
c        expansion on the skeleton of the daddy into
c        equivalent incoming expansions on the skeletons of the
c        two sonnies, and the matrix expand, converting outgoing
c        expansions on the skeletons of the sonnies into
c        an equivalent outgoing expansion on the daddy's
c        skeleton
c 
c                   Input parameters:
c 
c 
c  zs - the nodes in the simulation
c  zs2 - the locations of the reordered nodes in the structure, as
c        produced by a prior call to the subroutine bandgeo.
c  iz - the integer array describing the transposition that converts
c        the array zs into the array zs2, as produced by a prior call
c        to the subroutine bandgeo.
c  izinv - array defining the transposition inverse to that defined
c        by iz. Specifically,
c           zs(1,i)=zs2(1,izinv(i)),
c           zs(2,i)=zs2(2,izinv(i)),
c        for all i=1,2,...,n. Must have been produced by a preceding
c        call to bandgeo.
c  n - the number of nodes in the structure
c  nchunks2 - the number of chunks ON THE DADDY'S LEVEL
c  npts - the number of nodes to be contained in the outer skeleton
c        to be created, in addition to the nodes in the two chunks
c        nearest to the chunk ichunk (or in addition to the one node
c        nearest to the chunk ichunk, if ichunk=1 or ichunk=nchunks)
c  eps - the accuracy to which the calculations are to be performed
c  barrier - the distance from the boundary of the chunk of the
c        point where the "far-field region" starts. At this time
c        (7.20.02) it is set to the size of the chunk on level 0.
c  thresh2 - the (2,nchunks2) array of left and right boundaries of
c        all chunks ON THE DADDY'S LEVEL
c  ichunks2 - the integer (2,nchunks2) array defining the subdivision
c        of the nodes in the structure into chunks ON THE DADDY'S
C        LEVEL.
c  ichunk2 - the sequence number OF THE DADDY AMONG CHUNKS ON THE
C        DADDY'S LEVEL
c  irows_outin1 - the integer array containing the addresses in array
c        iz of the skeleton points in the first of the chunks to be
c        merged
c  ncols_outin1 - the number of elements in array icols_inout1
c  irows_outin2 - the integer array containing the addresses in array
c        iz of the skeleton points in the second of the chunks to be
c        merged
c  ncols_outin2 - the number of elements in array icols_inout2
c  icols_inout1 - the integer array containing the addresses in array
c        iz of the skeleton points in the first of the chunks to be
c        merged
c  ncols_inout1 - the number of elements in array icols_inout1
c  icols_inout2 - the integer array containing the addresses in array
c        iz of the skeleton points in the second of the chunks to be
c        merged
c  ncols_inout2 - the number of elements in array icols_inout2
c 
c 
c                         Output parameters:
c 
c  eval - the (ncols_outin1+ncols_outin2,ncols)-matrix converting
c        the incoming expansio for the daddy into the vector of
c        concatinated incoming expansions for the sonnies
c  irows - the array containing the addresses in array iz of the
c        incoming skeleton of the daddy
c  nrows - the number of elements in array irowss
c  expand - the (ncols,ncols_inout1+ncols_inout2)-matrix converting
c        the vector of concatinated far-field expansions for the
c        chunks to be merged into the far-field expansion for the
c        daddy
c  icols - the array containing the addresses in array iz of the
c        outgoing skeleton of the daddy
c  ncols - the number of elements in array icols
c 
c 
c                         Work arrays:
c 
c  w - must be ???? real *8 elements long
c 
c 
c       . . . construct the outer skeleton of the daddy chunk
c 
  
        iizskel=1
        lizskel=n
c 
        iicols=iizskel+lizskel
        licols=ncols_outin1+ncols_outin2
c 
        iirows0=iicols+licols
        lirows0=ncols_outin1+ncols_outin2
c 
        iw=iirows0+lirows0
  
        iflev0=0
c 
        call outer_skeleton(ichunks2,ichunk2,zs2,iz,izinv,n,
     1      thresh2,nchunks2,npts,w(iizskel),nzskel,w(iw),
     2      iflev0,barrier)
c 
c        . . . merge the incoming skeletons of the two chunks,
c              obtaining a chunk on the next level
c 
        ncols0=ncols_outin1+ncols_outin2
        ll=14*nzskel*ncols0+4*(nzskel+ncols0)+1000
        lused=ll+iw
c 
  
cccc        call prinf('in merges_1chunk before merge_outin, lenw=*',
cccc     1      lenw,1)
  
cccc        call prinf('in merges_1chunk before merge_outin, lused=*',
cccc     1      lused,1)
  
        if(lused .gt. lenw) then
            ier=32
            call prinf('bombing from merge_outin with ier=*',
     1            ier,1)
        endif
c 
        call merge_outin(ichunks2,ichunk2,nchunks2,
     1      thresh2,barrier,zs,zs2,iz,izinv,n,npts,w(iw),
     2      irows_outin1,irows_outin2,ncols_outin1,
     2      ncols_outin2,eps,
     4      eval,irows,ncols,
     5      w(iizskel),nzskel,w(iicols),w(iirows0) )
c 
         nrows=ncols
c 
c        merge the outgoing skeletons of the two chunks,
c        obtaining a chunk on the next level
c 
        iizskel=1
        lizskel=n
  
c 
        iirows=iizskel+lizskel
        lirows=ncols_inout1+ncols_inout2
c 
        iicols0=iirows+lirows
        licols0=ncols_inout1+ncols_inout2
c 
        iw=iicols0+licols0
c 
        ncols0=ncols_inout1+ncols_inout2
        ll=14*nzskel*ncols0+4*(nzskel+ncols0)+1000
        lused2=ll+iw
        if(lused .lt. lused2) lused=lused2
c 
        if(lused .gt. lenw) then
            ier=16
            call prinf('bombing from merge_outin with ier=*',
     1            ier,1)
        endif
c 
        call merge_inout(ichunks2,ichunk2,nchunks2,
     1      thresh2,barrier,zs,zs2,iz,izinv,n,npts,w(iw),
     2      icols_inout1,icols_inout2,ncols_inout1,
     3      ncols_inout2,eps,
     4      expand,icols,ncols,
     5      w(iizskel),nzskel,w(iirows),w(iicols0) )
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine merge_outin(ichunks2,ichunk2,nchunks2,
     1      thresh2,barrier,zs,zs2,iz,izinv,n,npts,w,
     2      irows_outin1,irows_outin2,ncols_outin1,
     3      ncols_outin2,eps,
     4      eval,irows,ncols,
     5      izskel,nzskel,icols,irows0)
        implicit real *8 (a-h,o-z)
        save
        dimension ichunks2(2,1),iz(1),zs(2,1),thresh2(2,1),
     1      zs2(2,1),w(1),izskel(1),
     2      irows_outin1(1),irows_outin2(1),irows0(1),
     3      icols(1),irows(1),iwhich(10)
c 
        complex *16 pot,eval(1)
c 
c        This subroutine merges the incoming skeletons of two
c        chunks, producing the daddy's skeleton. It also
c        produces the matrix eval, converting incoming
c        expansion on the skeleton of the daddy into
c        equivalent outgoing expansions on the skeletons
c        of two sonnies
c 
c                   Input parameters:
c 
c 
c  ichunks2 - the integer (2,nchunks2) array defining the subdivision
c        of the nodes in the structure into chunks ON THE DADDY'S
C        LEVEL.
c  ichunk2 - the sequence number OF THE DADDY AMONG CHUNKS ON THE
C        DADDY'S LEVEL
C  nchunks2 - the number of chunks on the daddy's level
c  thresh2 - the (2,nchunks2) array of left and right boundaries of
c        all chunks ON THE DADDY'S LEVEL
c  barrier - the distance from the boundary of the chunk of the
c        point where the "far-field region" starts. At this time
c        (7.20.02) it is set to the size of the chunk on level 0.
c  zs - the nodes in the simulation
c  zs2 - the locations of the reordered nodes in the structure, as
c        produced by a prior call to the subroutine bandgeo.
c  iz - the integer array describing the transposition that converts
c        the array zs into the array zs2, as produced by a prior call
c        to the subroutine bandgeo.
c  izinv - array defining the transposition inverse to that defined
c        by iz. Specifically,
c           zs(1,i)=zs2(1,izinv(i)),
c           zs(2,i)=zs2(2,izinv(i)),
c        for all i=1,2,...,n. Must have been produced by a preceding
c        call to bandgeo.
c  n - the number of nodes in the structure
c  npts - the number of nodes to be contained in the outer skeleton
c        to be created, in addition to the nodes in the two chunks
c        nearest to the chunk ichunk (or in addition to the one node
c        nearest to the chunk ichunk, if ichunk=1 or ichunk=nchunks)
c  irows_outin1 - the integer array containing the addresses in array
c        iz of the skeleton points in the first of the chunks to be
c        merged
c  irows_outin2 - the integer array containing the addresses in array
c        iz of the skeleton points in the second of the chunks to be
c        merged
c  ncols_outin1 - the number of elements in array icols_inout1
c  ncols_outin2 - the number of elements in array icols_inout2
c  eps - the accuracy to which the calculations are to be performed
c 
c                         Output parameters:
c 
c  eval - the (ncols_outin1+ncols_outin2,ncols)-matrix converting
c        the incoming expansio for the daddy into the vector of
c        concatinated incoming expansions for the sonnies
c  irows - the array containing the addresses in array iz of the
c        incoming skeleton of the daddy
c  ncols - the number of elements in array irows
c 
c                         Work arrays:
c 
c  w - must be at least 14*nzskel*ncols0+4*(nzskel+ncols0)+1000
c        real *8 elements long, with
c 
c        ncols0=ncols_outin1+ncols_outin2
c 
c  izskel - should be at least n integer *4 elements long
c  icols - should be at least ncols_inout1+ncols_inout2
c        integer *4 elements long
c  irows0 - should be at least ncols_inout1+ncols_inout2
c        integer *4 elements long
c 
c 
c        . . . merge the inner skeletons of the two sonnies,
c              obtaining the starting point for the daddy's
c              inner skeleton
c 
        do 1200 i=1,ncols_outin1
c 
        irows0(i)=irows_outin1(i)
 1200 continue
c 
        do 1400 i=1,ncols_outin2
c 
        irows0(i+ncols_outin1)=irows_outin2(i)
 1400 continue
c 
        ncols0=ncols_outin1+ncols_outin2
  
        iain=1
        lain=ncols0*nzskel*2+100
c 
        iw=iain+lain
c 
c       construct the incoming interaction matrix
c 
        call ainbld(nzskel,ncols0,irows0,
     1    izskel,iz,zs,w(iain) )
c 
c      skeletonize the incoming matrix for this chunk
c 
        iwhich(1)=0
        iwhich(2)=1
        iwhich(3)=1
        iwhich(4)=0
        iwhich(5)=0
        iwhich(6)=0
c 
        call cskeleton3(w(iain),ncols0,nzskel,eps,iwhich,irows,
     1      icols,ncols,aout,expand,eval,rowsout,colsout,
     2      errout,w(iw))
c 
        call prin2('in merge_out_in, errout=*',errout,1)
c 
        do 2200 i=1,ncols
c 
        irows(i)=irows0(irows(i))
 2200 continue
c 
        nnnp=ncols_outin1+ncols_outin2
        call bubble_rows(irows,ncols,eval,nnnp,w)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine ainbld(nzskel,ncols0,irows0,
     1    izskel,iz,zs,ain)
c 
        implicit real *8 (a-h,o-z)
        save
        dimension zs(2,1),iz(1),izskel(1),irows0(1)
        complex *16 ain(ncols0,1),pot
c 
c       construct the incoming interaction matrix
c 
        do 1800 i=1,nzskel
        do 1600 j=1,ncols0
c 
        iii=izskel(i)
        jjj=irows0(j)
c 
        call poteval(iz,zs,iii,jjj,pot)
        ain(j,i)=pot
cccc        ain(i,j)=pot
 1600 continue
 1800 continue
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine merge_inout(ichunks2,ichunk2,nchunks2,
     1      thresh2,barrier,zs,zs2,iz,izinv,n,npts,w,
     2      icols_inout1,icols_inout2,ncols_inout1,
     3      ncols_inout2,eps,
     4      expand,icols,ncols,
     5      izskel,nzskel,irows,icols0)
        implicit real *8 (a-h,o-z)
        save
        dimension ichunks2(2,1),iz(1),zs(2,1),thresh2(2,1),
     1      zs2(2,1),w(1),izskel(1),
     2      icols_inout1(1),icols_inout2(1),icols0(1),
     3      icols(1),irows(1),iwhich(10)
c 
        complex *16 expand(1)
c 
c        This subroutine merges the outgoing skeletons of two
c        chunks, producing the daddy's skeleton. It also
c        produces the matrix expand, converting outgoing
c        expansions on the skeletons of the sonnies into
c        an equivalent outgoing expansion on the daddy's
c        skeleton
c 
c                   Input parameters:
c 
c 
c  ichunks2 - the integer (2,nchunks2) array defining the subdivision
c        of the nodes in the structure into chunks ON THE DADDY'S
C        LEVEL.
c  ichunk2 - the sequence number OF THE DADDY AMONG CHUNKS ON THE
C        DADDY'S LEVEL
C  nchunks2 - the number of chunks on the daddy's level
c  thresh2 - the (2,nchunks2) array of left and right boundaries of
c        all chunks ON THE DADDY'S LEVEL
c  barrier - the distance from the boundary of the chunk of the
c        point where the "far-field region" starts. At this time
c        (7.20.02) it is set to the size of the chunk on level 0.
c  zs - the nodes in the simulation
c  zs2 - the locations of the reordered nodes in the structure, as
c        produced by a prior call to the subroutine bandgeo.
c  iz - the integer array describing the transposition that converts
c        the array zs into the array zs2, as produced by a prior call
c        to the subroutine bandgeo.
c  izinv - array defining the transposition inverse to that defined
c        by iz. Specifically,
c           zs(1,i)=zs2(1,izinv(i)),
c           zs(2,i)=zs2(2,izinv(i)),
c        for all i=1,2,...,n. Must have been produced by a preceding
c        call to bandgeo.
c  n - the number of nodes in the structure
c  npts - the number of nodes to be contained in the outer skeleton
c        to be created, in addition to the nodes in the two chunks
c        nearest to the chunk ichunk (or in addition to the one node
c        nearest to the chunk ichunk, if ichunk=1 or ichunk=nchunks)
c  icols_inout1 - the integer array containing the addresses in array
c        iz of the skeleton points in the first of the chunks to be
c        merged
c  icols_inout2 - the integer array containing the addresses in array
c        iz of the skeleton points in the second of the chunks to be
c        merged
c  ncols_inout1 - the number of elements in array icols_inout1
c  ncols_inout2 - the number of elements in array icols_inout2
c  eps - the accuracy to which the calculations are to be performed
c 
c                         Output parameters:
c 
c  expand - the (ncols,ncols_inout1+ncols_inout2)-matrix converting
c        the vector of concatinated far-field expansions for the
c        chunks to be merged into the far-field expansion for the
c        daddy
c  icols - the array containing the addresses in array iz of the
c        outgoing skeleton of the daddy
c  ncols - the number of elements in array icols
c 
c                         Work arrays:
c 
c  w - must be at least 14*nzskel*ncols0+4*(nzskel+ncols0)+1000
c        real *8 elements long, with
c 
c        ncols0=ncols_inout1+ncols_inout2
c 
c  izskel - should be at least n integer *4 elements long
c  icols0 - should be at least ncols_inout1+ncols_inout2
c        integer *4 elements long
c 
c        . . . merge the inner skeletons of the two sonnies,
c              obtaining the starting point for the daddy's
c              inner skeleton
c 
        do 1200 i=1,ncols_inout1
c 
        icols0(i)=icols_inout1(i)
 1200 continue
c 
        do 1400 i=1,ncols_inout2
c 
        icols0(i+ncols_inout1)=icols_inout2(i)
 1400 continue
c 
        ncols0=ncols_inout1+ncols_inout2
c 
c       construct the outgoing interaction matrix
c 
        iaaout=1
        laaout=ncols0*nzskel*2+100
c 
        iw=iaaout+laaout
c 
        call aaoutbld(nzskel,ncols0,icols0,
     1    izskel,iz,zs,w(iaaout) )
c 
c      skeletonize the incoming matrix for this chunk
c 
        iwhich(1)=1
        iwhich(2)=0
        iwhich(3)=1
        iwhich(4)=0
        iwhich(5)=0
        iwhich(6)=0
c 
        call cskeleton3(w(iaaout),nzskel,ncols0,eps,iwhich,irows,
     1      icols,ncols,aout,expand,eval,rowsout,colsout,
     2      errout,w(iw) )
c 
cccc        call prinf('in merge_inout after cskeleton, ncols=*',
cccc     1      ncols,1)
cccc        call prinf('in merge_inout after cskeleton, irows=*',
cccc     1      irows,ncols)
cccc        call prin2('in merge_inout, errout=*',errout,1)
  
        do 2200 i=1,ncols
c 
        icols(i)=icols0(icols(i))
 2200 continue
c 
        nnnp=ncols_inout1+ncols_inout2
        call bubble_cols(icols,ncols,expand,nnnp,w)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine aaoutbld(nzskel,ncols0,icols0,
     1    izskel,iz,zs,aaout)
c 
        implicit real *8 (a-h,o-z)
        save
        dimension zs(2,1),iz(1),izskel(1),icols0(1)
        complex *16 aaout(nzskel,ncols0),pot
c 
c       construct the outgoing interaction matrix
c 
        do 1800 i=1,nzskel
        do 1600 j=1,ncols0
c 
        iii=izskel(i)
        jjj=icols0(j)
c 
cccc        call poteval(iz,zs,iii,jjj,pot)
        call poteval(iz,zs,jjj,iii,pot)
        aaout(i,j)=pot
 1600 continue
 1800 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine skels_1chunk(ier,n,eps,zs,
     1      nchunks,ichunk,npts,bord,lbord,iplot,
     2      zs2,iz,izinv,thresh,ichunks,w,lenw,lused,
     3      eval_outin,irows_outin,ncols_outin,
     4      expand_inout,icols_inout,ncols_inout)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension w(1),zs(2,1),zs2(2,1),iz(1),izinv(1),
     1      thresh(2,1),ichunks(2,1),expand_inout(1),
     2      eval_outin(1),irows_outin(1),icols_inout(1),
     3      bord(2,1)
c 
c       This subroutine constructs inner skeletons (for both
c       incoming and outgoing fields) for a single chunk. It
c       uses data constructed via a preceding call to the
c       subroutine bandgeo (see).
c 
c 
c                 Input parameters:
c 
c  n - the number of nodes in the structure
c  eps - the accuracy to which the calculations are to be performed
c  zs - the nodes in the simulation
c  nchunks - the number of chunks in the structure
c  ichunk - the sequence number of the chunk for which the inner
c        skeletons are to be constructed
c  npts - the number of nodes to be contained in the outer skeleton
c        to be created, in addition to the nodes in the two chunks
c        nearest to the chunk ichunk (or in addition to the one node
c        nearest to the chunk ichunk, if ichunk=1 or ichunk=nchunks)
c  bord - the boundary of the region to be simulated; it is used by
c        this subroutine only for plotting, and is not used at all
c        if the parameter iplot (see below) has been set to 0 by
c        the user.
c  lbord - the number of elements in the array lbord; like the
c        latter, it is used by this subroutine only for plotting,
c        and is not used at all if the parameter iplot (see below)
c        has been set to 0 by the user.
c  iplot - the GNUPLOT file number on which the simulation will be
c        plotted. Setting iplot=0 suppresses the plotting.
c  zs2 - the locations of the reordered nodes in the structure, as
c        produced by a prior call to the subroutine bandgeo.
c  iz - the integer array describing the transposition that converts
c        the array zs into the array zs2, as produced by a prior call
c        to the subroutine bandgeo.
c  izinv - array defining the transposition inverse to that defined
c        by iz. Specifically,
c           zs(1,i)=zs2(1,izinv(i)),
c           zs(2,i)=zs2(2,izinv(i)),
c        for all i=1,2,...,n. Must have been produced by a preceding
c        call to bandgeo.
c  thresh - the (2,nchunks) array of left and right boundaries of
c        all chunks in the structure.  Must have been produced by a
c        preceding call to bandgeo.
c  ichunks - the integer (2,nchunks) array defining the subdivision
c        of the nodes in the structure into chunks. Must have been
c        produced by a prior call to the subroutine bandgeo.
c 
c 
c                 Output parameters:
c 
c  eval_outin - the (np,ncols_outin)-matrix converting the potential
c        at the ncols "skeleton" nodes inside the chunk number ichunk
c        (defined by the array irows_outin below) into the potential
c        at all np nodes in that chunk.
c  irows_outin - coordinates (np of them things) in array iz of the
c        elements of the inner skeleton for the incoming potential
c        in the chunk number ichunk.
c  ncols_outin - the number of elements in the array irows_outin
c        (inter alia)
c  expand_inout - the (ncols_inout,np)-matrix converting a charge
c        distribution at the np nodes in the chunk number ichunk
c        into an equivalent charge distribution on the elements
c        of the inner skeleton for the outgoing potential in the
c        chunk number ichunk.
c  icols_inout - coordinates (np of them things) in array iz of the
c        elements of the inner skeleton for the outgoing potential
c        in the chunk number ichunk.
c  ncols_inout - the number of elements in the array icols_out
c        (inter alia)
c 
c                 Work arrays:
c 
c  w - must be at least ?????  elements long
c 
c 
c       . . . construct the outer skeleton of the chunk number
c             ichunk (please note that the same outer skeleton
c             is used for both incoming and outgoing potentials)
c 
        ier=0
c 
        iizskel=1
        lizskel=n+10
c 
        iw=iizskel+lizskel
c 
        iflev0=1
c 
        lused=iw+6*n+2*npts+1000
cccc        call prinf('lused=*',lused,1)
cccc        call prinf('while lenw=*',lenw,1)
c 
        if(lused .gt. lenw) then
            ier=1024
            call prinf('bombing from skels_1chunk with ier=*',
     1              ier,1)
            return
        endif
c 
        call outer_skeleton(ichunks,ichunk,zs2,iz,izinv,n,
     1      thresh,nchunks,npts,w(iizskel),nzskel,w(iw),
     2      iflev0,barrier)
c 
c        if the user so requested, plot this chunk,
c        together with its outer skeleton
c 
        if(iplot .ne. 0)
     1      call plot_1chunk(zs2,ichunks,ichunk,bord,lbord,
     2      iplot,w(iizskel),nzskel,zs,iz,w(iw) )
c 
c        construct the inner skeletons
c 
        np=ichunks(2,ichunk)-ichunks(1,ichunk)+1
        izstart=ichunks(1,ichunk)
c 
        iain=iizskel+lizskel
        lain=np*nzskel*2+100
c 
        iicols=iain+lain
        licols=np+nzskel
c 
        iirows=iicols
c 
        iw=iicols+licols
c 
        lused2=iw+14*nzskel*np+4*(nzskel+np)+1000
  
cccc        call prinf('lused2=*',lused2,1)
cccc        call prinf('while lenw=*',lenw,1)
c 
       if(lused .lt. lused2) lused=lused2
c 
        if(lused2 .gt. lenw) then
            ier=32
            call prinf('bombing from skels_1chunk with ier=*',
     1              ier,1)
            return
        endif
c 
        call skel_out_in(iz,zs,np,nzskel,w(iain),ichunks,
     1        ichunk,w(iizskel),eps,eval_outin,
     2        irows_outin,w(iicols),ncols_outin,w(iw) )
c 
        call skel_in_out(iz,zs,np,nzskel,w(iain),ichunks,
     1        ichunk,w(iizskel),eps,expand_inout,
     2        w(iirows),icols_inout,ncols_inout,w(iw) )
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine skel_out_in(iz,zs,np,nzskel,ain,ichunks,
     1      ichunk,izskel2,eps,eval,irows,icols,ncols,w)
        implicit real *8 (a-h,o-z)
c 
        save
        complex *16 ain(np,nzskel),pot,eval(1),w(1)
        dimension iz(1),zs(2,1),izskel2(1),iwhich(10),
     1      irows(1),icols(1),ichunks(2,1)
c 
c       This subroutine constructs the inner skeleton for
c       the incoming interactions of the chunk number ichunk,
c       i.e. for the interactions between this chunk and the
c       world when the charges outside the chunk create the
c       potential at the nodes inside it.
c 
c                   Input parameters:
c 
c 
c  iz - the integer array describing the transposition that converts
c        the array zs into the array zs2, as produced by a prior call
c        to the subroutine bandgeo.
c  zs - the nodes in the simulation
c  np - the number of nodes in the chunk number ichunk
c  nzskel - the number of elements in the array izskel2 (see below)
c  ichunks - the integer (2,nchunks) array defining the subdivision
c        of the nodes in the structure into chunks. Must have been
c        produced by a prior call to the subroutine bandgeo.
c  ichunk - the sequence number of the chunk for which the outer
c        skeleton is to be constructed
c  izskel2 - the array of coordinates in the array iz of all nodes in
c        the outer skeleton of the chunk number ichunk.
c        Explanation: izskel(i) is the location in the array iz of
c        the coordinate in array zs of the i-th node in the outer
c        skeleton of the chunk number ichunk. PLEASE NOTE THAT THE
C        ARRAY IZSKEL IS OBTAINED FROM THE ARRAY IZSKEL2 VIA THE
C        FORMULA
C 
C        IZSKEL(I)=IZ(IZSKEL2(I))                                   (1)
C 
c  eps - the accuracy to which the calculations are to be performed
c 
c                 Output parameters:
c 
c  eval - the (np,ncols)-matrix converting the potential at the ncols
c        "skeleton" nodes inside the chunk number ichunk into the
c         potential at all np nodes in that chunk.
c  ncols - the number of elements in the skeleton
c  irows - the sequence numbers (np of them things) in array iz of
c          the elements of its inner skeleton. Together with the matrix
c          eval (see above) the principal output of this dubroutine.
c  icols - the sequence numbers (ncols of them things) in array iz
c          of the elements in the OUTER skeleton of this chunk that
c          were used in the construction of the inner skeleton.
c          Probably, will not be useful for anything
c 
c                  Works arrays:
c 
c  w - must be at least 14*nzskel*np+4*(nzskel+np)+1000 real *8
c          elements long
c  ain - must be at least np*ncols*2 real *8 elements long
c 
c       . . . construct the inner skeleton of chunk number ichunk
c 
        izstart=ichunks(1,ichunk)
c 
c        construct the incoming matrix of interactions for this chunk
c 
        do 1600 i=1,nzskel
        do 1400 j=1,np
c 
        iii=izskel2(i)
        jjj=izstart+j-1
c 
        call poteval(iz,zs,iii,jjj,pot)
        ain(j,i)=pot
 1400 continue
 1600 continue
c 
c      skeletonize the incoming matrix for this chunk
c 
        iwhich(1)=0
        iwhich(2)=1
        iwhich(3)=1
        iwhich(4)=0
        iwhich(5)=0
        iwhich(6)=0
  
        call cskeleton3(ain,np,nzskel,eps,iwhich,irows,
     1      icols,ncols,aout,expand,eval,rowsout,colsout,
     2      errout,w)
c 
cccc        call prinf('after cskeleton, ncols=*',ncols,1)
cccc        call prinf('after cskeleton, irows=*',irows,ncols)
cccc        call prin2('in skel_out_in, errout=*',errout,1)
  
        do 2400 i=1,ncols
c 
        irows(i)=irows(i)+izstart-1
 2400 continue
c 
c       reorder the rows in order of increasing number
c 
cccc        call prinf('before bubble_rows, irows=*',irows,ncols)
  
        call bubble_rows(irows,ncols,eval,np,w)
  
  
        return
        end
c 
c 
c 
c 
c 
        subroutine bubble_rows(irows,nrows,eval,np,w)
        implicit real *8 (a-h,o-z)
        save
        complex *16 eval(np,nrows),w(np,nrows)
        dimension irows(nrows),inds(10000)
c 
c       bubble the array irows
c 
        do 1200 i=1,nrows
c 
        inds(i)=i
 1200 continue
c 
cccc        call prinf('in bubble_rows, irows before bubbling*',
cccc     1      irows,nrows)
  
        do 1600 i=1,nrows
c 
        do 1400 j=1,nrows-1
c 
        if(irows(j) .lt. irows(j+1)) goto 1400
c 
        jj=irows(j)
        irows(j)=irows(j+1)
        irows(j+1)=jj
c 
        jj=inds(j)
        inds(j)=inds(j+1)
        inds(j+1)=jj
c 
 1400 continue
 1600 continue
  
  
cccc        call prinf('in bubble_rows, irows after bubbling*',
cccc     1      irows,nrows)
  
cccc        stop
c 
c       reorder the rows of the matrix eval
c 
        do 2400 i=1,nrows
c 
        do 2200 j=1,np
c 
        ii=inds(i)
        w(j,i)=eval(j,ii)
 2200 continue
 2400 continue
c 
        call cleascop(w,eval,np*nrows)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cmatvec77(a,n,m,x,y)
        implicit complex *16 (a-h,o-z)
        save
        dimension a(n,m),x(m),y(n)
c 
        do 1600 i=1,n
c 
        cd=0
        do 1400 j=1,m
c 
        cd=cd+a(i,j)*x(j)
 1400 continue
c 
        y(i)=cd
 1600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine chgexpand_in(ii,iz,zs,irows,ncols,vals)
        implicit real *8 (a-h,o-z)
        save
        complex *16 zs(2,1),vals(ncols)
        dimension irows(ncols),iz(1)
c 
c       construct the values at the "significant" points
c       of the potential located at the point ii
c 
        jj=iz(ii)
        jj=ii
        do 2200 i=1,ncols
c 
        iii=irows(i)
c 
        call poteval(iz,zs,jj,iii,vals(i))
 2200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine skel_in_out(iz,zs,np,nzskel,aaout,ichunks,
     1      ichunk,izskel2,eps,expand,irows,icols,ncols,w)
        implicit real *8 (a-h,o-z)
c 
        save
        complex *16 aaout(nzskel,np),pot,expand(1),w(1)
c 
        dimension iz(1),zs(2,1),izskel2(1),irows(1),
     1      icols(1),iwhich(10),ichunks(2,1)
c 
c       This subroutine constructs the inner skeleton for
c       the outgoing interactions of the chunk number ichunk,
c       i.e. for the interactions between this chunk and the
c       world when the charges inside the chunk create the
c       potential at the nodes outside it.
c 
c                   Input parameters:
c 
c 
c  iz - the integer array describing the transposition that converts
c        the array zs into the array zs2, as produced by a prior call
c        to the subroutine bandgeo.
c  zs - the nodes in the simulation
c  np - the number of nodes in the chunk number ichunk
c  nzskel - the number of elements in the array izskel2 (see below)
c  ichunks - the integer (2,nchunks) array defining the subdivision
c        of the nodes in the structure into chunks. Must have been
c        produced by a prior call to the subroutine bandgeo.
c  ichunk - the sequence number of the chunk for which the outer
c        skeleton is to be constructed
c  izskel2 - the array of coordinates in the array iz of all nodes in
c        the outer skeleton of the chunk number ichunk.
c        Explanation: izskel(i) is the location in the array iz of
c        the coordinate in array zs of the i-th node in the outer
c        skeleton of the chunk number ichunk. PLEASE NOTE THAT THE
C        ARRAY IZSKEL IS OBTAINED FROM THE ARRAY IZSKEL2 VIA THE
C        FORMULA
C 
C        IZSKEL(I)=IZ(IZSKEL2(I))                                   (1)
C 
c  eps - the accuracy to which the calculations are to be performed
c 
c                 Output parameters:
c 
c  expand - the (ncols,np)-matrix converting a charge distribution
c        at the np nodes in the chunk number ichunk into an equivalent
c        charge distribution on the outgoing skeleton of size ncols
c  ncols - the number of elements in the skeleton
c  irows - the sequence numbers (ncols of them things) in array iz
c          of the elements in the OUTER skeleton of this chunk that
c          were used in the construction of the inner skeleton.
c          Probably, will not be useful for anything
c  icols - the sequence numbers (n of them things) in array iz of
c          the elements of its inner skeleton. Together with the matrix
c          expand (see above) the principal output of this dubroutine.
c 
c                  Works arrays:
c 
c  w - must be at least 14*nzskel*np+4*(nzskel+np)+1000 real *8
c          elements long
c  aaout - must be at least np*ncols*2 real *8 elements long
c 
c 
c 
c        . . . construct the matrix of outgoing interactions
c              for this chunk
c 
        izstart=ichunks(1,ichunk)
c 
        do 1600 i=1,nzskel
        do 1400 j=1,np
c 
        iii=izskel2(i)
        jjj=izstart+j-1
c 
        call poteval(iz,zs,jjj,iii,pot)
c 
        aaout(i,j)=pot
 1400 continue
 1600 continue
c 
c      skeletonize the matrix of outgoing interactions
c      for this chunk
c 
        iwhich(1)=1
        iwhich(2)=0
        iwhich(3)=1
        iwhich(4)=0
        iwhich(5)=0
        iwhich(6)=0
c 
        call cskeleton3(aaout,nzskel,np,eps,iwhich,irows,
     1      icols,ncols,aout,expand,eval,rowsout,colsout,
     2      errout,w)
c 
        do 2400 i=1,ncols
c 
        icols(i)=icols(i)+izstart-1
 2400 continue
c 
  
        call bubble_cols(icols,ncols,expand,np,w)
  
  
        return
        end
c 
c 
c 
c 
c 
        subroutine bubble_cols(icols,ncols,expand,np,w)
        implicit real *8 (a-h,o-z)
        save
        complex *16 expand(ncols,np),w(ncols,np)
        dimension icols(ncols),inds(10000)
c 
c       bubble the array icols
c 
        do 1200 i=1,ncols
c 
        inds(i)=i
 1200 continue
c 
cccc        call prinf('in bubble_cols, icols before bubbling*',
cccc     1      icols,ncols)
  
        do 1600 i=1,ncols
c 
        do 1400 j=1,ncols-1
c 
        if(icols(j) .lt. icols(j+1)) goto 1400
c 
        jj=icols(j)
        icols(j)=icols(j+1)
        icols(j+1)=jj
c 
        jj=inds(j)
        inds(j)=inds(j+1)
        inds(j+1)=jj
c 
 1400 continue
 1600 continue
  
  
cccc        call prinf('in bubble_cols, icols after bubbling*',
cccc     1      icols,ncols)
  
cccc        stop
c 
c       reorder the cols of the matrix eval
c 
        do 2400 i=1,ncols
c 
        do 2200 j=1,np
c 
        ii=inds(i)
        w(i,j)=expand(ii,j)
 2200 continue
 2400 continue
c 
        call cleascop(w,expand,np*ncols)
c 
        return
        end
c 
c 
c 
c 
c 
c 
        subroutine expanout_eval(iz,zs,chgs,icols,ncols,ii,val)
        implicit real *8 (a-h,o-z)
        save
        complex *16 chgs(1),val,cd
        dimension icols(1),iz(1),zs(2,1)
c 
c       evaluate the potential at the point number ii
c 
        val=0
        do 1200 i=1,ncols
c 
        j=icols(i)
  
        call poteval(iz,zs,j,ii,cd)
c 
        val=val+cd*chgs(i)
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine chunkexpand_out(expand,np,ncols,
     1      charges,vals)
        implicit real *8 (a-h,o-z)
        save
        complex *16 charges(1),vals(ncols),expand(ncols,np)
c 
c       expand the potential of all np charges in the chunk
c       into the potential of the ncols "significant" charges
c 
        call cmatvec77(expand,ncols,np,charges,vals)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine outer_skeleton(ichunks,ichunk,zs2,iz,izinv,n,
     1      thresh,nchunks,npts,izskel,nzskel,w,
     2      iflev0,barrier)
        save
  
        implicit real *8 (a-h,o-z)
        dimension zs2(2,1),thresh(2,1),izinv(1),
     1      iz(1),izskel(1),w(1),ichunks(2,1)
c 
c        This subroutine constructs the "outer skeleton" for a
c        single chunk.
c 
c                 Input parameters:
c 
c  ichunks - the integer (2,nchunks) array defining the subdivision
c        of the nodes in the structure into chunks. Must have been
c        produced by a prior call to the subroutine bandgeo.
c 
c  ichunk - the sequence number of the chunk for which the outer
c        skeleton is to be constructed
c  zs2 - the locations of the reordered nodes in the structure, as
c        produced by a prior call to the subroutine bandgeo.
c  iz - the integer array describing the transposition that converts
c        the array zs into the array zs2, as produced by a prior call
c        to the subroutine bandgeo.
c  n - the number of nodes in the structure
c  thresh - the (2,nchunks) array of left and right boundaries of
c        all chunks. in the structure
c  nchunks - the number of chunks in the structure
c  npts - the number of nodes to be contained in the outer skeleton
c        to be created, in addition to the nodes in the two chunks
c        nearest to the chunk ichunk (or in addition to the one node
c        nearest to the chunk ichunk, if ichunk=1 or ichunk=nchunks)
c 
c                 Output parameters:
c 
c  izskel - the array of coordinates in the arrays iz, zs2 of all
c        nodes in the outer skeleton of the chunk number ichunk.
c        Explanation: izskel(i) is the location in the arrays
c        iz, zs2 of the i-th node in the outer skeleton of the
c        chunk number ichunk.
c  nzskel - the number of elements in the array izskel.
c 
c                 Work arrays:
c 
c  w - must be at least 6*n+2*npts+1000 real *8 elements long
c 
c        . . . allocate memory for the construction of the outer
c              skeleton of chunk number ichunk
c 
        iiz3=20
        liz3=n+10
c 
        izs3=iiz3+liz3
        lzs3=n*2+10
c 
        iws3=izs3+lzs3
        lws3=n+10
c 
        ithresh3=iws3+lws3
        lthresh3=2*n+10
c 
        ixsrand=ithresh3+lthresh3
        lxsrand=npts+10
c 
        iinums=ixsrand+lxsrand
        linums=npts+10
c 
        call outer_skeleton0(ichunk,zs2,iz,n,thresh,nchunks,
     1      w(iiz3),w(izs3),n3,w(iws3),w(ithresh3),
     2      w(ixsrand),w(iinums),npts,izskel,
     3      iflev0,barrier)
c 
c        to the outer skeleton of the chunk number ichunk,
c        add the elements in the chunks adjoint  to it
c 
        nzskel=npts
c 
c       . . . if this is a level 0 chunk
c 
        if(iflev0 .ne. 1) goto 4000
  
  
        if(ichunk .eq. 1) goto 2400
c 
        j1=ichunks(1,ichunk-1)
        j2=ichunks(2,ichunk-1)
        nj=j2-j1+1
        do 2200 i=1,nj
c 
        izskel(nzskel+i)=iz(j1+i-1)
 2200 continue
c 
        nzskel=nzskel+nj
c 
 2400 continue
c 
        if(ichunk .eq. nchunks) goto 2800
c 
        j1=ichunks(1,ichunk+1)
        j2=ichunks(2,ichunk+1)
        nj=j2-j1+1
        do 2600 i=1,nj
c 
        izskel(nzskel+i)=iz(j1+i-1)
 2600 continue
c 
        nzskel=nzskel+nj
c 
 2800 continue
c 
        do 3200 i=1,nzskel
c 
        izskel(i)=izinv(izskel(i))
 3200 continue
c 
        return
c 
 4000 continue
c 
c       . . . if this is not level 0 chunk
c 
        if(ichunk .eq. 1) goto 4400
c 
        j1=ichunks(1,ichunk-1)
        j2=ichunks(2,ichunk-1)
        nj=j2-j1+1
c 
        njj=0
        do 4200 i=1,nj
c 
        j=j1+i-1
        d=zs2(1,j)
        if(d .lt. thresh(1,ichunk)-barrier) goto 4200
        njj=njj+1
  
        izskel(nzskel+njj)=iz(j1+i-1)
 4200 continue
c 
        nzskel=nzskel+njj
 4400 continue
c 
        if(ichunk .eq. nchunks) goto 5400
c 
        j1=ichunks(1,ichunk+1)
        j2=ichunks(2,ichunk+1)
        nj=j2-j1+1
c 
        njj=0
        do 5200 i=1,nj
c 
        j=j1+i-1
        d=zs2(1,j)
        if(d .gt. thresh(2,ichunk)+barrier) goto 5200
        njj=njj+1
  
        izskel(nzskel+njj)=iz(j1+i-1)
 5200 continue
c 
        nzskel=nzskel+njj
 5400 continue
c 
        do 5600 i=1,nzskel
c 
        izskel(i)=izinv(izskel(i))
 5600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine outer_skeleton0(ichunk,zs2,iz,n,thresh,nchunks,
     1      iz3,zs3,n3,ws3,thresh3,xsrand,inums,npts,izskel,
     2      iflev0,barrier)
        implicit real *8 (a-h,o-z)
        save
        dimension zs2(2,1),thresh(2,1),zs3(2,1),
     1      iz3(1),iz(1),ws3(1),thresh3(2,1),
     2      xsrand(1),inums(1),izskel(1)
c 
c       build the list of all nodes that are outside this chunk
c       and its immediate neighbors
c 
c        . . . if this is a chunk on level 0
c 
        if(iflev0 .eq. 1) goto 2150
  
        i=ichunk
c 
  
        call prinf('ichunk=*',ichunk,1)
        call prin2('barrier=*',barrier,1)
  
cccc         stop
  
        if(i .ge. 2) th1=thresh(1,i)-barrier
        if(i .eq. 1) th1=thresh(1,1)
c 
        if(i .le. nchunks-1) th2=thresh(2,i)+barrier
        if(i .eq. nchunks) th2=thresh(2,nchunks)
c 
        ii=0
        do 2100 i=1,n
c 
        d=zs2(1,i)
        if( (d .gt. th1) .and. (d .lt. th2) ) goto 2100
c 
        ii=ii+1
        zs3(1,ii)=zs2(1,i)
        zs3(2,ii)=zs2(2,i)
c 
        iz3(ii)=iz(i)
 2100 continue
c 
        n3=ii
c 
        goto 2300
c 
 2150 continue
c 
        i=ichunk
c 
        if(i .ge. 3) th1=thresh(2,i-2)
        if(i .lt. 3) th1=thresh(1,1)
c 
        if(i .le. nchunks-2) th2=thresh(2,i+1)
        if(i .gt. nchunks-2) th2=thresh(2,nchunks)
c 
        ii=0
        do 2200 i=1,n
c 
        d=zs2(1,i)
        if( (d .gt. th1) .and. (d .lt. th2) ) goto 2200
c 
        ii=ii+1
        zs3(1,ii)=zs2(1,i)
        zs3(2,ii)=zs2(2,i)
c 
        iz3(ii)=iz(i)
 2200 continue
c 
        n3=ii
c 
 2300 continue
c 
c        find the weights associated with all nodes in zs3
c 
        cent=(thresh(1,ichunk)+thresh(2,ichunk))/2
c 
        dd=0
        do 2400 i=1,n3
c 
        d=1/abs(cent-zs3(1,i))
        ws3(i)=d
        dd=dd+d
 2400 continue
c 
        do 2600 i=1,n3
c 
        ws3(i)=ws3(i)/dd
 2600 continue
c 
c        construct the subdivision of the interval [0,1] such
c        that the length of the i-th subinterval is equal to
c        the weight of the i-th node
c 
        thresh3(1,1)=0
        do 2800 i=1,n3
c 
        thresh3(2,i)=thresh3(1,i)+ws3(i)
        thresh3(1,i+1)=thresh3(2,i)
 2800 continue
c 
c        use the random number generator to select
c        the skeleton of the part from which the chunk
c        number ichunk is separated and with which it
c        interacts
c 
        call newrand(npts,xsrand)
c 
        do 3200 i=1,npts
c 
        call bisect(ier,thresh3,n3,xsrand(i),inums(i))
 3200 continue
c 
        do 3400 i=1,npts
c 
        j=inums(i)
        jj=iz3(j)
        izskel(i)=jj
 3400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine bisect(ier,thr,n3,x,inum)
        implicit real *8 (a-h,o-z)
        save
        dimension thr(2,1)
c 
        ier=0
        i1=1
        i2=n3
c 
        do 2000 i=1,100
c 
        if( i2-i1 .gt. 1) goto 1100
c 
        if( (x .ge. thr(1,i1)) .and. (x .le. thr(2,i1)) ) then
            inum=i1
            return
        endif
c 
        if( (x .ge. thr(1,i2)) .and. (x .le. thr(2,i2)) ) then
            inum=i1
            return
        endif
c 
 1100 continue
c 
        i3=(i1+i2)/2
c 
        if(x .gt. thr(2,i3) ) goto 1200
        if(x .lt. thr(1,i3) ) goto 1200
c 
        inum=i3
        return
 1200 continue
c 
        if(x .gt. thr(1,i3) ) goto 1400
c 
        i2=i3
        goto 2000
c 
 1400 continue
        i1=i3
c 
 2000 continue
c 
        ier=4
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine bandgeo(ier,zs,n,nchunks,endleft,endright,
     1      zs2,iz,izinv,thresh,ichunks,w)
        implicit real *8 (a-h,o-z)
        save
        dimension zs(2,1),zs2(2,1),w(2),iz(1),izinv(1),
     1      thresh(2,1),ichunks(2,1)
c 
        external bandelcomp
c 
c        This subroutine performs the preliminary geometry
c        processing for the direct solver in a band. Starting
c        with the user-specified collection of n nodes zs,
c        it sorts them in order of increasing abscissae,
c        subdivides the band into chunks (the number of chunks
c        to be created is user-supplied, as well as the ends
c        of the band), and allocates the points to appropriate
c        chunks. It also creates the auxiliary control structures
c        iz, izinv, thresh, ichunks (see below). PLEASE NOTE THAT
c        THE BAND IS ASSUMED TO BE HORIZONTAL.
c 
c                     Input parameters:
c 
c  zs - the nodes in the band
c  n - the number of nodes in the band
c  nchunks - the number of equal chunks into which the band is to be
c        subdivided
c  endleft, endright - ends of the band
c 
c                     Output parameters:
c 
c  ier - error return code.
c          ier=0 means successful execution
c          ier=64 means that at least one of the nodes in array zs
c                  was outside the user-specified limits endleft,
c                  endright. This is a fatal error.
c  zs2 - the nodes from zs, sorted in order of increazing zs(1,i)
c  iz - the integer array describing the transposition converting
c        zs into zs2. Specifically,
c           zs2(1,i)=zs(1,iz(i)),
c           zs2(2,i)=zs(2,iz(i)),
c        for all i=1,2,...,n.
c  izinv - array defining the transposition inverse to that defined
c        by iz. Specifically,
c           zs(1,i)=zs2(1,izinv(i)),
c           zs(2,i)=zs2(2,izinv(i)),
c        for all i=1,2,...,n.
c  thresh - the (2,nchunks) array of left and right boundaries of
c        chunks. Specifically,
c           chunks(1,i) is the left boundary of the i-th chunk
c           chunks(2,i) is the right boundary of the i-th chunk
c  ichunks - integer array (2,nchunks) assigning the nodes in
c        array zs2 to the chunks in array thresh. Specifically
c           ichunks(1,i) is the sequence number in the array zs2
c                        of the first node in the chunk number i
c           ichunks(2,i) is the sequence number in the array zs2
c                        of the last node in the chunk number i
c 
c                     Work arrays:
c 
c  w - must be at least n*2+4 real *8 locations long.
c 
c 
c        sort (horizontally) the user-specified nodes
c 
        ier=0
        do 1200 i=1,n
c 
        zs2(1,i)=zs(1,i)
        zs2(2,i)=i
 1200 continue
c 
        call sortrear(zs2,2,n,bandelcomp,w)
c 
        do 1400 i=1,n
c 
        iz(i)=zs2(2,i)+0.1d0
        j=iz(i)
        zs2(1,i)=zs(1,j)
        zs2(2,i)=zs(2,j)
c 
        izinv(j)=i
 1400 continue
c 
c        construct the chunks
c 
        h=(endright-endleft)/nchunks
        thresh(1,1)=endleft
        thresh(2,nchunks)=endright
c 
        do 1600 i=1,nchunks
c 
        thresh(2,i)=thresh(1,i)+h
        thresh(1,i+1)=thresh(2,i)
 1600 continue
c 
c       allocate nodes to the chunks
c 
        if( (zs2(1,1) .lt. thresh(1,1)) .or.
     1      (zs2(1,1) .gt. thresh(2,1)) ) then
c 
            ier=64
            return
        endif
c 
        if( (zs2(1,n) .lt. thresh(1,nchunks)) .or.
     1      (zs2(1,n) .gt. thresh(2,nchunks) ) ) then
c 
            ier=64
            return
        endif
c 
        ichunk=1
c 
        ichunks(1,1)=1
        ichunks(2,nchunks)=n
c 
        do 1800 i=1,n
c 
        if(zs2(1,i) .le. thresh(2,ichunk)) goto 1800
c 
        ichunks(2,ichunk)=i-1
        ichunk=ichunk+1
        ichunks(1,ichunk)=i
 1800 continue
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine bandelcomp(a,b,m,ifagtb)
        implicit real *8 (a-h,o-z)
        save
        dimension a(2),b(2)
c 
        ifagtb =0
        if(a(1) .ge. b(1)) ifagtb=1
        return
        end
c 
c 
c 
c 
c 
        subroutine newrand(num,xs)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1)
        integer  *4 alpha,beta,gamma
c 
c 
        data ifcalled/0/
c 
        done=1
        if(ifcalled .eq. 0) m=37333
        ifcalled=1
c 
        alpha=2**7+9
        beta=1
        gamma=2**22
c 
        do 1800 i=1,100
c 
        m1=m*alpha+beta
        j=m1/gamma
        m=m1-j*gamma
c 
 1800 continue
c 
        do 2200 i=1,num
c 
        m1=m*alpha+beta
        j=m1/gamma
        m=m1-j*gamma
c 
        xs(i)=m*done/gamma
 2200 continue
        return
c 
c 
c 
c 
        entry newrand_init(m7)
        m=m7
        ifcalled=1
        return
        end
  
  
  
