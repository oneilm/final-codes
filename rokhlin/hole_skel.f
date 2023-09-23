        implicit real *8 (a-h,o-z)
        dimension z(2,1000 000),w(10000 000),skelin(2,10000),
     1      skelout(2,10000),int_skelin(10000),
     2      int_skelout(10000),iholes(100000),
     3      skelin2(2,10000)
c 
        complex *16 rk,ima,expand(1000 000),eval(1000 000),
     1      artif(1000 000),uu(1000 000),vv(1000 000),
     2      ss(10000),ww(1000 000),aa(1000 000),
     3      expand_rhs(1000 000),s(1000 000)
c 
        external pot2user10,pot2user01,pot2user11
c 
        data ima/(0.0d0,1.0d0)/
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
         PRINT *, 'ENTER nsmall'
         READ *,nsmall
         CALL PRINF('nsmall=*',nsmall,1 )
c 
c        Construct and plot the test structure
c 
        rk=10
        eps=1.0d-4
        nd_perside=4
        coef=1.2
ccc        coef=1
c 
c 
        call circles_bld(nsmall,nd_perside,coef,z,iholes)
  
        n=nd_perside**2*nsmall
  
        call pot2user0_init(z,n)
c 
        iw=21
        call zquaplot(iw,z,n,2,'test structure*')
c 
        i1=1
        i2=nsmall
  
  
        i1=1+nsmall*3
        i2=nsmall*4
  
        i1=1+nsmall*6
        i2=nsmall*7
  
  
c 
        iw=22
        call zquaplot2(iw,z,n,2,z(1,i1),nsmall,2,
     1      'test structure*')
  
        coedist=200
        coedist=2
cccc        coedist=10
        nframe=100
cccc        nframe=10
c 
        iplot=41
        lenw=10 000 000
  
        call  hole_skel(ier,z,n,i1,i2,eps,
     1      pot2user10,pot2user01,pot2user11,rk,
     2      coedist,nframe,iplot,
     3      skelin,int_skelin,skelout,int_skelout,ncols,
     4      w,lenw,lused,expand,eval,s,s_cond,artif,artif_cond,
     5      expand_rhs)
c 
        call prinf('after hole_skel, ier=*',ier,1)
        call prinf('after hole_skel, lused=*',lused,1)
        call prin2('after hole_skel, s_cond=*',s_cond,1)
        call prin2('after hole_skel, artif_cond=*',artif_cond,1)
c 
c       plot the arrays skelin, skelout
c 
        iw=32
        call zquaplot2(iw,skelin,ncols,2,skelout,ncols,2,
     1      'skelin and skelout*')
c 
        call prinf('and int_skelin=*',int_skelin,ncols)
cccc        call prinf('and int_skelout=*',int_skelout,ncols)
        call prin2('and skelin=*',skelin,ncols*2)
cccc        call prin2('and skelout=*',skelout,ncols*2)
c 
c        reconstruct the arrays skelin, skelout from
c        the maps int_skelin, int_skelout and plot them things
c 
        do 2600 i=1,ncols
c 
        j=int_skelin(i)
        skelin2(1,i)=z(1,j)
        skelin2(2,i)=z(2,j)
c 
 2600 continue
c 
        iw=33
        call zquaplot2(iw,skelin,ncols,2,skelout,ncols,2,
     1      'skelin and skelout reconstructed from int_skelout*')
c 
c        test the matrix expand
c 
        nhole=i2-i1+1
        itest=1
  
        call expand_test(expand,ncols,nhole,rk,
     1      i1,int_skelin,n)
  
  
cccc        stop
  
c 
c        test the matrix eval
c 
        nhole=i2-i1+1
        itest=5
c 
        call eval_test(eval,ncols,nhole,rk,
     1      i1,int_skelin,n,itest)
  
  
cccc        stop
  
c 
c        svd the stupid matrix
c 
        eps2=1.0d-10
        lww=1000 000
        call csvdpiv(ier,artif,ncols,ncols,uu,vv,ss,ncols2,eps2,
     1      ww,lww,ltot)
c 
        call prin2('after svdpiv of artif, ss=*',ss,ncols2*2)
        call prinf('while ncols2=*',ncols2,1)
        call prinf('and ncols=*',ncols,1)
  
        call prin2('and condition number of artif is=*',
     1      ss(1)/ss(ncols),1)
  
        call prin2('and artif_cond=*',artif_cond,1)
        call prin2('and s_cond=*',s_cond,1)
  
  
  
        stop
c 
c       test the "artif" matrix
c 
        call artif_test(i1,nhole,int_skelin,ncols,artif,aa,
     1      rk)
c 
c        test the expand_rhs matrix
c 
        call artif_expand_test(i1,nhole,
     1      ncols,artif,expand,expand_rhs,aa,rk)
  
  
        stop
        end
  
c 
c 
c 
c 
c 
        subroutine artif_expand_test(i1,nhole,
     1      ncols,artif,expand,expand_rhs,a,rk)
        implicit real *8 (a-h,o-z)
c 
        save
        complex *16 sigma(10000),sigma2(10000),
     1      rhs(10000),rhs_skel(10000),a(nhole,nhole),
     2      work(1000 000),b(1000 000),sigma_skel(10 000),
     3      rk,pot,artif(ncols,ncols),expand(1),expand_rhs(1)
c 
c        construct the matrix of interactions of all
c        nodes comprizing the hole
c 
        itest=1
        do 1400 i=1,nhole
        do 1200 j=1,nhole
c 
        ii=i1+i-1
        jj=i1+j-1
c 
        call pot2user11(ii,jj,rk,pot)
        a(j,i)=pot
 1200 continue
 1400 continue
  
  
cccc        call prin2('in artif_expand_test, a=*',a,nhole**2*2)
  
c 
c       convert it into the scattering matrix
c 
        call corthom(a,nhole,work,cond)
c 
c       construct the test rhs
c 
        do 1600 i=1,nhole
c 
        rhs(i)=i
 1600 continue
c 
c       obtain the induced charge
c 
        call cmatvec77(a,nhole,nhole,rhs,sigma)
c 
        call prin2('in artif_expand_test, sigma=*',
     1      sigma,nhole*2)
c 
c        expand the obtained sigma
c 
        call cmatvec77(expand,ncols,nhole,sigma,sigma_skel)
  
        call prin2('in artif_expand_test, sigma_skel=*',
     1      sigma_skel,ncols*2)
c 
c        "artificially expand" rhs
c 
        call cmatvec77(expand_rhs,ncols,nhole,rhs,rhs_skel)
c 
c        invert the user-supplied matrix artif, and apply
c        the inverse to rhs_skel
c 
        call cleascop(artif,b,ncols**2)
c 
        call corthom(b,ncols,work,cond)
c 
        call cmatvec77(b,ncols,ncols,rhs_skel,sigma2)
c 
        call prin2('in artif_expand_test, sigma2=*',
     1      sigma2,ncols*2)
c 
        do 2200 i=1,ncols
c 
        sigma2(i)=sigma2(i)-sigma_skel(i)
 2200 continue
  
        call prin2('and differences=*',sigma2,ncols*2)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine artif_test(i1,nhole,
     1      int_skelin,ncols,artif,a,rk)
        implicit real *8 (a-h,o-z)
c 
        save
        complex *16 sigma(10000),sigma2(10000),
     1      rhs(10000),rhs2(10000),a(nhole,nhole),
     2      work(1000 000),b(1000 000),
     3      cd1,rk,pot,cd2,artif(ncols,ncols)
c 
        dimension int_skelin(1)
c 
c        construct the matrix of interactions of all
c        nodes comprizing the hole
c 
        itest=1
        do 1400 i=1,nhole
        do 1200 j=1,nhole
c 
        ii=i1+i-1
        jj=i1+j-1
c 
        call pot2user11(ii,jj,rk,pot)
        a(j,i)=pot
 1200 continue
 1400 continue
c 
c       convert it into the scattering matrix
c 
        call corthom(a,nhole,work,cond)
c 
c       evaluate the potential of the test charge at the
c       nodes, calculate the induced cgarge, and calculate
c       the induced potential at the test node
c 
        do 1600 i=1,nhole
c 
        ii=i1+i-1
        call pot2user11(itest,ii,rk,pot)
c 
        rhs(i)=pot
 1600 continue
c 
c       obtain the induced charge
c 
        call cmatvec77(a,nhole,nhole,rhs,sigma)
c 
        call prin2('in artif_test, sigma=*',sigma,nhole*2)
c 
        cd1=0
        do 1800 i=1,nhole
c 
        ii=i1+i-1
        call pot2user11(ii,itest,rk,pot)
c 
        cd1=cd1+sigma(i)*pot
 1800 continue
c 
        call prin2('and induced potential=*',cd1,2)
c 
c 
c 
c        invert the user-supplied matrix artif, obtaining
c        the scattering matrix on the skeleton
c 
        call cleascop(artif,b,ncols**2)
c 
        call corthom(b,ncols,work,cond)
c 
c       evaluate the potential of the test charge at the
c       nodes, calculate the induced cgarge, and calculate
c       the induced potential at the test node
c 
        do 2600 i=1,ncols
c 
        ii=int_skelin(i)
        call pot2user11(itest,ii,rk,pot)
c 
        rhs2(i)=pot
 2600 continue
c 
c       obtain the induced charge
c 
        call cmatvec77(b,ncols,ncols,rhs2,sigma2)
c 
        call prin2('in artif_test, sigma2=*',sigma2,ncols*2)
c 
        cd2=0
        do 2800 i=1,nhole
c 
        ii=int_skelin(i)
        call pot2user11(ii,itest,rk,pot)
c 
        cd2=cd2+sigma2(i)*pot
 2800 continue
c 
        call prin2('and second induced potential=*',cd2,2)
        call prin2('and cd2-cd1=*',cd2-cd1,2)
  
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine eval_test(eval,ncols,nhole,
     1      rk,i1,int_skelin,n,itest)
        implicit real *8 (a-h,o-z)
        save
        complex *16 eval(nhole,ncols),sigma(10000),
     1      sigma2(10000),rk,pot,sigma3(10000),
     2      diffs(10000)
c 
        dimension int_skelin(1)
c 
c        construct the array of incoming potentials
c        at the nodes in the skeleton of the hole
c 
        do 1200 i=1,ncols
c 
        call pot2user11(itest,int_skelin(i),rk,pot)
c 
        sigma(i)=pot
 1200 continue
c 
cccc        call prin2('sigma=*',sigma,ncols*2)
c 
c       evaluate the potential at the nodes comprizing
c       the hole
c 
        call cmatvec77(eval,nhole,ncols,sigma,sigma2)
c 
        call prin2('in eval_test, sigma2=*',sigma2,nhole*2)
c 
c       evaluate the potentials at all nodes in the hole
c 
        do 1600 i=1,nhole
c 
        call pot2user11(itest,i1+i-1,rk,pot)
c 
        sigma3(i)=pot
c 
        diffs(i)=sigma3(i)-sigma2(i)
 1600 continue
  
cccc        call prin2('and sigma3=*',sigma3,ncols*2)
        call prin2('in eval_test, diffs=*',diffs,nhole*2)
cccc        call prin2('and eval=*',eval,nhole*2*ncols)
c 
        call prinf('and nhole=*',nhole,1)
        call prinf('and ncols=*',ncols,1)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine expand_test(expand,ncols,nhole,
     1      rk,i1,int_skelin,n)
        implicit real *8 (a-h,o-z)
        save
        complex *16 expand(ncols,nhole),sigma(10000),
     1      sigma2(10000),cd1,rk,pot,cd2
c 
        dimension int_skelin(1),errs(10000)
c 
c        construct the array of charges at the nodes
c        comprizing the hole
c 
        do 1200 i=1,nhole
c 
        sigma(i)=i
 1200 continue
c 
        call prin2('sigma=*',sigma,nhole*2)
c 
c       evaluate the (allegedly) equivalent charges
c       the the skeleton nodes
c 
        call cmatvec77(expand,ncols,nhole,sigma,sigma2)
c 
        call prin2('and sigma2=*',sigma2,ncols*2)
c 
c       evaluate the potential at a test point via both
c       sets of charges and compare them things
c 
        i2=i1+nhole-1
        nitest=0
        dd=0
        ddd=0
        do 2000 itest=1,n
c 
        if( (itest .ge. i1) .and. (itest .le. i2) ) goto 2000
c 
        nitest=nitest+1
c 
        cd1=0
        do 1600 i=1,nhole
c 
        call pot2user11(i1+i-1,itest,rk,pot)
c 
        cd1=cd1+sigma(i)*pot
 1600 continue
c 
cccc        call prin2('cd1=*',cd1,2)
c 
        cd2=0
  
        do 1800 i=1,ncols
c 
        call pot2user11(int_skelin(i),itest,rk,pot)
c 
        cd2=cd2+sigma2(i)*pot
 1800 continue
  
cccc        call prin2('cd2=*',cd2,2)
cccc        call prin2('cd2-cd1=*',cd2-cd1,2)
c 
        errs(nitest)=abs(cd2-cd1)
c 
        dd=dd+abs(cd2-cd1)**2
        ddd=ddd+abs(cd2)**2
  
 2000 continue
  
cccc        call prin2('and errs=*',errs,nitest)
        call prin2('in expand_test, relative error=*',
     1      sqrt(dd/ddd),1)
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
        subroutine pot2user00(z17,z27,rk,pot)
        implicit real *8 (a-h,o-z)
        save
        complex *16 pot,rk
        dimension z1(2),z2(2),z7(2,1),z(2,100000),
     1      z17(2),z27(2)
c 
        z1(1)=z17(1)
        z1(2)=z17(2)
c 
        z2(1)=z27(1)
        z2(2)=z27(2)
c 
        call pot2user_bottom(z1,z2,rk,pot)
c 
        return
c 
c 
c 
c 
        entry pot2user10(i1,z27,rk,pot)
c 
        z2(1)=z27(1)
        z2(2)=z27(2)
c 
        z1(1)=z(1,i1)
        z1(2)=z(2,i1)
c 
        call pot2user_bottom(z1,z2,rk,pot)
c 
        return
c 
c 
c 
c 
        entry pot2user01(z17,i2,rk,pot)
c 
        z1(1)=z17(1)
        z1(2)=z17(2)
c 
        z2(1)=z(1,i2)
        z2(2)=z(2,i2)
  
        call pot2user_bottom(z1,z2,rk,pot)
c 
        return
c 
c 
c 
c 
        entry pot2user11(i1,i2,rk,pot)
c 
        if(i1 .eq. i2) then
            pot=-100
cccc            pot=-1
cccc            pot=-0.1
           return
        endif
c 
        z1(1)=z(1,i1)
        z1(2)=z(2,i1)
c 
        z2(1)=z(1,i2)
        z2(2)=z(2,i2)
c 
        call pot2user_bottom(z1,z2,rk,pot)
c 
        return
c 
c 
c 
c 
        entry pot2user0_init(z7,n7)
        n=n7
c 
        do 1200 i=1,n
c 
        z(1,i)=z7(1,i)
        z(2,i)=z7(2,i)
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine pot2user_bottom(z1,z2,rk,pot)
        implicit real *8 (a-h,o-z)
        save
        complex *16 pot,rk,zz,h0,h1
        dimension z1(2),z2(2)
c 
c 
        d=(z1(1)-z2(1))**2+(z1(2)-z2(2))**2
  
        d=sqrt(d)
        zz=d*rk
c 
        ifexpon=1
        call hank103(zz,h0,h1,ifexpon)
c 
        pot=h0 * z1(1)
cccc        pot=h0 * z2(1)
cccc        pot=h0
  
cccc        pot=2*log(d)
  
cccc        pot=d**2
  
        return
c 
        d=(z1(1)+2000)**2+(z1(2)-z2(2))**2
        d=(z1(1)-z2(1)+20.1)**2+(z1(2)-z2(2))**2
  
cccc        d=(z1(1))**2+(z1(2)-z2(2)-100)**2
  
cccc        d=(z2(1)-500)**2
cccc        d=(z1(1)-500)**2
  
  
cccc        call prin2('d=*',d,1)
  
        d=sqrt(d)
        zz=d*rk
c 
        ifexpon=1
        call hank103(zz,h0,h1,ifexpon)
c 
        pot=pot+h0
  
cccc        d=h0
cccc        pot=d
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine circles_bld(n,nd_perside,coef,zs,iholes)
        implicit real *8 (a-h,o-z)
        save
        integer *4 iholes(1)
        real *8 centers(2,10 000),zs0(2,10 000),zs(2,n,1)
c 
c        construct the centers of disks in the simulation
c 
        done=1
        h=2*done/(nd_perside-1)
        nn=0
        do 1400 i=1,nd_perside
        do 1200 j=1,nd_perside
c 
        nn=nn+1
        centers(1,nn)=-1+(j-1)*h
        centers(2,nn)=-1+(i-1)*h
 1200 continue
 1400 continue
c 
c       construct the basic circle
c 
        pi=atan(done)*4
        h=2*pi/n
        r=done/nd_perside *coef
c 
        do 1600 i=1,n
c 
        zs0(1,i)=r*cos((i-1)*h)
        zs0(2,i)=r*sin((i-1)*h)
 1600 continue
c 
c       create the whole mess of them
c 
        ii=0
        do 2000 i=1,nn
c 
        do 1800 j=1,n
c 
        zs(1,j,i)=zs0(1,j)+centers(1,i)
        zs(2,j,i)=zs0(2,j)+centers(2,i)
c 
        ii=ii+1
c 
        iholes(ii)=i
 1800 continue
 2000 continue
c 
c        . . . rotate them things
c 
        phi=1.1
        a11=cos(phi)
        a22=cos(phi)
        a12=sin(phi)
        a21=-sin(phi)
c 
        do 2400 i=1,nn
        do 2200 j=1,n
c 
        x=a11*zs(1,j,i)+a12*zs(2,j,i)
        y=a21*zs(1,j,i)+a22*zs(2,j,i)
c 
        zs(1,j,i)=x
        zs(2,j,i)=y
 2200 continue
 2400 continue
c 
        return
        end
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code, and the beginning
c        of the single hole skeletonization code proper
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine hole_skel(ier,z,n,i1,i2,eps,
     1      potuser10,potuser01,potuser11,rk,
     2      coedist,nframe,iplot,
     3      skelin,int_skelin,skelout,int_skelout,ncols,
     4      w,lenw,lused,expand,eval,s,s_cond,artif,artif_cond,
     5      expand_rhs)
        implicit real *8 (a-h,o-z)
        save
        dimension z(2,1),w(1),
     1      skelin(2,1),skelout(2,1),int_skelin(1),
     2      int_skelout(1),expand(1),eval(1)
        complex *16 rk,artif(1),expand_rhs(1),s(1)
c 
c 
c        This subroutine skeletonizes one hole. It returns both the
c        inner and the outer skeletons; the skeletons are returned
c        both as nodes in R^2 and as maps of their locations in array
c        z of user-supplied nodes. Ths subroutine also returns the
c        matrices expand, eval, artif, expand_rhs. The matrix expand
c        "expands"  a charge distribution defined on all elements of
c        which the hole consists into an equivalent charge distribu-
c        tions on the skeleton. The matrix eval converts values of an
c        incoming potentials defined on the skeleton into the values
c        of the same incoming potential on all elements of the hole.
c        Finally, the subroutine returns the artificial matrix of
c        interactions artif. EXPLANATION: if the hole is replaced
c        with its skeleton, and the matrix of interactions between
c        the elements of the hole is replaced with the matrix artif,
c        the scattering matrix of the hole (with respect to the rest
c        of the structure) is not changed, to precision eps.
c 
c 
c                           Input parameters:
c 
c 
c  z - the points in R^2 on top of which the tree structure
c        is placed.
c  n - the number of points in array z
c  i1 - the address in z of the first element on the hole to be
c        skeletonized
c  i2 - the address in z of the last element on the hole to be
c        skeletonized
c  eps - the accuracy to which the SVDs will be performed
c  lenw - the length (in real *8 words) of the array w provided
c        by the user
c  potuser10, potuser01, potuser11 - user-supplied subroutines
c        defining the various types of interactions needed by
c        this subroutine. Following are the decriptions of the
c        three
c 
c  potuser10 - the interaction between a (source) node on some apriori
c        discretized boundary, specified by its sequence number among
c        such nodes, and a target that is an point in R^2 that is not
c        a part of any apriori diacretized boundary. The
c        calling sequence is
c 
c                   potuser10(i1,z2,rk,pot)
c 
c        with
c 
c     i1 - the sequence number of the source (input)
c     z2 - the location of the target (input)
c     rk - the Helmholtz coefficient  (input)
c     pot - the potential created by the source on the target (output)
c 
c  potuser01 - the interaction between a (source) point in R^2 that
c        is not a part of any apriori diacretized boundary, and a
c        (target) node on some apriori discretized boundary, specified
c        by its sequence number among such nodes. The calling sequence
c        is
c 
c                   potuser10(z1,i2,rk,pot)
c 
c        with
c 
c     z1 - the location of the source (input)
c     z2 - the sequence number of the target (input)
c     rk - the Helmholtz coefficient  (input)
c     pot - the potential created by the source on the target (output)
c 
c  potuser11 - the interaction between a (source) node on some apriori
c        discretized boundary, specified by its sequence number among
c        such nodes, and a (target) node on some apriori discretized
c        boundary, also specified by its sequence number among such
c        nodes. The calling sequence is
c 
c                   potuser11(i1,i2,rk,pot)
c 
c        with
c 
c     z1 - the sequence number of the source (input)
c     z2 - the sequence number of the target (input)
c     rk - the Helmholtz coefficient  (input)
c     pot - the potential created by the source on the target (output)
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c  rk - the Helmholtz coefficient
c  coedist - the ratio of the radius of the artificial (circular)
c        boundary to the radius of the hole; a good number is 2
c  nframe - the number of equispaced nodes into which the artificial
c        boundary will be discretized
c  iplot - tells the subroutine how to name the files on which the
c        GNUPLOT-readable data are to be stored. For example,
c        specifying iplot=41, will cause the subroutine to produce
c        files with names gn41, gn42 (plus several auxiliary ones).
c        Setting iplot=0 will cause all plotting to be suppressed.
c  lenw - the length (in real *8 words) of the work array provided
c        by the user
c 
c                        Output parameters:
c 
c  ier - error return code.
c     ier=0 means successful conclusion
c     ier=8 means that the amoyunt lenw of space in array w provided
c        by the user is insufficient. This is a fatal error.
c  skelin - the physical locations of points in R^2 comprizing the
c        final inner skeleton
c  int_skelin - the locations in array iz of the sequence numbers of
c        elements in the inner final skeleton (please note that all
c        of these are actual nodes provided by the user)
c  skelout - the physical locations of points in R^2 comprizing the
c        final outer skeleton
c  int_skelout - the locations in array iz of the sequence numbers of
c        elements in the outer final skeleton (please note that NOT
c        all of these are actual nodes provided by the user; the
c        addresses in iz of those elements that are not present there
c        do not make too much sense, and are set to -13)
c  ncols - the number of elements in skelin, skelout, int_skelin,
c        int_skelout
c  lused - the number of elements in the work array w that were
c        actually used by the subroutine
c  expand - the matrix converting an incoming potential on the
c        skeleton of the hole into the incoming potential
c        on all nodes constituting the hole
c  eval - the matrix converting an incoming potential on the
c        skeleton of the hole into the incoming potential
c        on all nodes constituting the hole
c  s - the "original" scattering matrix defined on the original
c        nodes on the hole being skeletonized
c  artif - the "artificial" matrix of interactions between the nodes
c        comprizing the inner skeleton of the hole.
c   EXPLANATION: if the hole is replaced with its skeleton, and the
c        matrix of interactions between the elements of the hole is
c        replaced with the matrix artif, the scattering matrix of
c        the hole (with respect to the rest of the structure) is
c        not changed, to precision eps.
c  expand_rhs - matrix converting a user-supplied right-hand side
c        tabulated at the nhole original nodes into an artificial
c        right-hand side tabulated at the skeleton nodes; the change
c        is invisible (to precision eps) to an observer outside the
c        hole.
c 
c                        Work arrays:
c 
c  w - must be longish
c 
c        . . . allocate memory
c 
        ier=0
c 
        izshole=1
        lzshole=(i2-i1+1)*2+10
c 
        izsothers=izshole+lzshole
        lzsothers=n*2+10
c 
        izsout7=izsothers+lzsothers
        lzsout7=(n+nframe)*2+10
c 
        iframe=izsout7+lzsout7
        lframe=nframe*2+10
c 
        iizsin=iframe+lframe
        lizsin=(i2-i1+1)+10
c 
        iinds=iizsin+lizsin
        linds=n+nframe+10
c 
        iizsout=iinds+linds
        lizsout=n+nframe+10
c 
        iirows=iizsout+lizsout
        lirows=n+nframe+10
c 
        iicols=iirows+lirows
        licols=n+nframe+10
c 
        iamatr=iicols+licols
        lenamatr=lenw-iamatr-10
c 
        call hole_skel0(ier,z,n,i1,i2,eps,
     1      potuser10,potuser01,potuser11,rk,
     2      coedist,nframe,
     3      skelin,int_skelin,skelout,int_skelout,ncols,
     4      w(izshole),w(izsothers),w(izsout7),w(iframe),
     5      w(iizsin),w(iinds),w(iizsout),w(iirows),
     6      w(iicols),w(iamatr),lenamatr,lused,
     7      expand,eval,s,artif,iplot,s_cond,artif_cond)
c 
        lused=lused+iamatr
c 
c        construct the matrix expand_rhs
c 
        nhole=i2-i1+1
        iaa=1
        laa=nhole**2*2+100
c 
        iwww=iaa+laa
c 
        call artif_expand(i1,nhole,potuser11,
     1      ncols,expand,artif,w(iaa),w(iwww),rk,expand_rhs)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine artif_expand(i1,nhole,potuser11,
     1      ncols,expand,artif,a,work,rk,expand_rhs)
        implicit real *8 (a-h,o-z)
c 
        save
        complex *16 a(nhole,nhole),expand(ncols,nhole),
     1      artif(ncols,ncols),work(1),rk,pot
c 
c        This subroutine constructs the matrix expand_rhs,
c        whose purpose in life is to convert a user-supplied
c        right-hand side tabulated at the nhole original
c        nodes into an artificial right-hand side tabulated
c        at the skeleton nodes; the change is invisible (to
c        precision eps) to an observer outside the hole.
c 
c                      Input parameters:
c 
c  i1 - the address in z of the first element on the hole to be
c        skeletonized
c  potuser11 - the interaction between a (source) node on some apriori
c        discretized boundary, specified by its sequence number among
c        such nodes, and a (target) node on some apriori discretized
c        boundary, also specified by its sequence number among such
c        nodes. The calling sequence is
c 
c                   potuser11(i1,i2,rk,pot)
c 
c        with
c 
c     z1 - the sequence number of the source (input)
c     z2 - the sequence number of the target (input)
c     rk - the Helmholtz coefficient  (input)
c     pot - the potential created by the source on the target (output)
c 
c  ncols - rank of interactions of the hole with the world;
c        must have been produced by a preceding call to the
c        subroutine  hole_skel0
c  expand - the matrix converting an incoming potential on the
c        skeleton of the hole into the incoming potential
c        on all nodes constituting the hole; must have been
c        produced by a preceding call to the subroutine  hole_skel0
c  artif - the "artificial" matrix of interactions between the
c        nodes comprizing the inner skeleton of the hole; must have
c        been produced by a preceding call to the subroutine
c        hole_skel0
c  rk - the Helmholtz coefficient
c 
c                        Output parameters:
c 
c  expand_rhs - matric converting a user-supplied right-hand
c        side tabulated at the nhole original nodes into an
c        artificial right-hand side tabulated at the skeleton
c        nodes; the change is invisible (to precision eps) to
c        an observer outside the hole.
c 
c                        Work arrays:
c 
c  a,work - must be sufficiently long
c 
c        . . . construct the matrix of interactions of all
c              nodes comprizing the hole
c 
        itest=1
        do 1400 i=1,nhole
        do 1200 j=1,nhole
c 
        ii=i1+i-1
        jj=i1+j-1
c 
        call potuser11(ii,jj,rk,pot)
        a(j,i)=pot
 1200 continue
 1400 continue
c 
c       convert it into the scattering matrix
c 
        call corthom(a,nhole,work,cond)
c 
c       construct the matrix of "artificial interpolation"
c 
        call cleamat(ier,expand,ncols,nhole,a,nhole,nhole,work)
c 
        call cleamat(ier,artif,ncols,ncols,work,ncols,nhole,
     1      expand_rhs)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine hole_skel0(ier,z,n,i1,i2,eps,
     1      potuser10,potuser01,potuser11,rk,
     2      coedist,nframe,
     3      skelin,int_skelin,skelout,int_skelout,ncols,
     4      zshole,zsothers,zsout,frame,
     5      izsin,inds,izsout,irows,icols,amatr,lenamatr,lused,
     6      expand,eval,s,artif,iplot,s_cond,artif_cond)
c 
        implicit real *8 (a-h,o-z)
        save
        dimension z(2,1),irows(1),icols(1),
     1      zsout(2,1),inds(1),izsin(1),frame(2,1),
     2      skelin(2,1),skelout(2,1),int_skelin(1),
     3      int_skelout(1),amatr(1),
     4      zshole(2,1),zsothers(2,1),izsout(1)
c 
        complex *16 rk,expand(1),eval(1),artif(1),s(1)
c 
        ier=0
c 
        i3=0
        do 1200 i=i1,i2
c 
        i3=i3+1
        zshole(1,i3)=z(1,i)
        zshole(2,i3)=z(2,i)
 1200 continue
c 
        nhole=i3
c 
        i3=0
        do 1400 i=1,n
c 
        if( (i .ge. i1) .and. (i .le. i2) ) goto 1400
c 
        i3=i3+1
        zsothers(1,i3)=z(1,i)
        zsothers(2,i3)=z(2,i)
 1400 continue
c 
        nothers=i3
c 
c       construct the preliminary skeletons
c 
        call prelim_skel(zshole,nhole,zsothers,nothers,
     1      coedist,nframe,zsout,inds,nout,frame,iplot,rk)
c 
c       convert the array inds into the map on the array z
c       (as opposed to the map on the array zsothers)
c 
        i3=0
        do 1800 i=1,nout
c 
        izsout(i)=-13
        if(inds(i) .lt. 0) goto 1800
c 
        if(inds(i) .lt. i1) izsout(i)=inds(i)
        if(inds(i) .gt. i1) izsout(i)=inds(i)+i2-i1+1
c 
 1800 continue
c 
c        construct the final skeletons
c 
        do 2200 i=1,nhole
c 
        izsin(i)=i1+i-1
 2200 continue
c 
        iamatr=1
        lamatr=nhole*nout*8+100
c 
        ia_outin=iamatr+lamatr
        la_outin=nhole*nout*2+100
c 
        iaa2=ia_outin+la_outin
        laa2=nhole**2*2+100
c 
        irnorms=iaa2+laa2
        lrnorms=nout+10
c 
        irnorms3=irnorms+lrnorms
        lrnorms3=nout+10
c 
        iww2=irnorms3+lrnorms3
        lww2=nhole*nout*2+100
c 
        lused=iww2+lww2
c 
        call prinf('in hole_skel0 lenamatr=*',lenamatr,1)
        call prinf('in hole_skel0 lused=*',lused,1)
c 
        if(lused .gt. lenamatr) then
            ier=8
            return
        endif
c 
        call fin_skel(izsin,z(1,i1),nhole,
     1      izsout,zsout,nout,irows,icols,
     2      rk,eps,amatr,skelin,skelout,ncols,
     3      int_skelin,int_skelout,
     4      potuser01,potuser10,potuser11,z,amatr(iww2),
     5      amatr(iaa2),amatr(ia_outin),s_cond,amatr(irnorms),
     6      amatr(ia_outin),amatr(irnorms3))
c 
        if(iplot .eq. 0) goto 2400
c 
        iw=iplot+1
        call zquaplot2(iw,skelin,ncols,2,skelout,ncols,2,
     1      'skelin and skelout*')
c 
 2400 continue
c 
c        construct the matrix expand
c 
        iaa=1
        laa=ncols*nout*2+100
c 
        icc=iaa+laa
        lcc=nout*nhole*2+100
c 
        iww=icc+lcc
        lww=nout*nhole*8 +500
c 
        irnorms=iww+lww
        lrnorms=nhole*2+10
c 
        lused2=irnorms+lrnorms
        if(lused2 .gt. lused) lused=lused2
c 
        if(lused .gt. lenamatr) then
            ier=8
            return
        endif
c 
        call expandmatr(skelin,int_skelin,
     1      zsout,izsout,nout,ncols,izsin,nhole,
     2      potuser10,potuser11,rk,eps/10,amatr(iaa),amatr(icc),
     3      expand,amatr(iww),amatr(irnorms) )
c 
c        construct the matrix eval
c 
        call evalmatr(skelin,int_skelin,
     1      zsout,izsout,nout,ncols,izsin,nhole,
     2      potuser01,potuser11,rk,eps/10,amatr(iaa),amatr(icc),
     3      eval,amatr(iww),amatr(irnorms) )
c 
c        construct the version of the interaction matrix
c        operating on the skeleton of the hole
c 
        iaa=1
        laa=nhole**2*2+100
c 
        iw2=iaa+laa
        lw2=lenamatr-iw2-10
c 
        call artif_matr(ncols,izsin,nhole,
     1      potuser11,rk,s,expand,eval,artif,
     2      amatr(iw2),lw2,ltot,artif_cond)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine artif_matr(nskel,izshole,nhole,
     1      potuser11,rk,a,expand,eval,artif,w2,lw,ltot,
     2      artif_cond)
        implicit real *8 (a-h,o-z)
        save
        complex *16 expand(nskel,nhole),artif(nskel,nskel),
     1      eval(nhole,nskel),w2(1),rk,pot,a(nhole,nhole)
c 
        dimension izshole(1)
c 
c        construct the matrix of interactions of the nodes
c        comprizing the hole with themselves
c 
        ier=0
c 
        do 1400 i=1,nhole
        do 1200 j=1,nhole
c 
        ii=izshole(i)
        jj=izshole(j)
        call potuser11(jj,ii,rk,pot)
c 
        a(i,j)=pot
 1200 continue
 1400 continue
c 
c       invert a
c 
        call corthom(a,nhole,w2,cond)
  
        call prin2('in artif_matr after first corthom, cond=*',
     1      cond,1)
c 
c        multiply the inverse of a by eval
c 
        call ncleamult(a,eval,w2,nhole,nhole,nskel)
c 
c        multiply the result by expand
c 
        call ncleamult(expand,w2,artif,nskel,nhole,nskel)
c 
c        find the maximum singular value of artif^{-1}
c 
        call max_singval(artif,nskel,nskel,svmax,w2,
     1      w2(nskel*2+1))
  
        call prin2('in artif_matr, svmax=*',svmax,1)
c 
c       invert the resulting scattering matrix on the
c       skeleton, getting the matrix artif
c 
        call corthom(artif,nskel,w2,cond)
c 
c        find the maximum singular value of artif
c 
        call max_singval(artif,nskel,nskel,svmin,w2,
     1      w2(nskel*2+1))
  
        call prin2('in artif_matr, svmin=*',svmin,1)
  
        artif_cond=svmin*svmax
  
        call prin2('in artif_matr, artif_cond=*',artif_cond,1)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine evalmatr(skelin,int_skelin,
     1      skelout,int_skelout,nout,nskel,
     2      izshole,nhole,potuser01,potuser11,rk,eps,
     3      a,c,eval,w,rnorms)
        implicit real *8 (a-h,o-z)
        save
        complex *16 eval(nhole,nskel),a(nskel,nout),
     1      c(nhole,nout),w(1),rnorms(1),rk,pot
c 
        dimension skelin(2,1),int_skelin(1),skelout(2,1),
     1      int_skelout(1),izshole(1)
c 
c        This subroutine constructs the "eval" matrix for a
c        single hole
c 
c                         Input parameters:
c 
c  skelin - the physical locations in R^2 of the nodes comprizing
c        the final inner skeleton of the hole being processed
c  int_skelin - the locations in the input array z of the nodes
c        comprizing final the inner skeleton of the hole being
c        processed
c  skelout - the physical locations in R^2 of the nodes comprizing
c        the preliminary outer skeleton of the hole being processed
c  int_skelout - the locations in the input array z of the nodes
c        comprizing the preliminary outer skeleton of the hole
c        being processed
c  nout - the number of elements in the outer skeleton
c  nskel - the number of nodes in the inner skeleton
c  izshole - the locations in array z of all nodes comprizing the
c        hole
c  nhole - the number of elements in array izshole; also the number
c        of nodes of which the hole consists
c  potuser01, potuser11 - user-supplied subroutines defining the
c        various types of interactions; see subroutine r2d_allmatrs
c        for a detailed dscription of these.
c  rk - the Helmholtz coefficient
c  eps - the accuracy to which the SVDs will be performed
C 
C                         OUTPUT PARAMETERS:
C 
c  eval - the matrix converting an incoming potential on the
c        skeleton of the hole into the incoming potential
c        on all nodes constituting the hole
c 
c                         Work arrays:
c 
c  a,c,w,rnorms
c 
c 
c        . . . construct the matrix of interactions of the inner
c              skeleton of the hole with the outer skeleton of
c              the said hole
c 
        do 1400 i=1,nout
        do 1200 j=1,nskel
c 
        jj=int_skelin(j)
        ii=int_skelout(i)
c 
        if(ii .lt. 0) call potuser01(skelout(1,i),jj,rk,pot)
        if(ii .gt. 0) call potuser11(ii,jj,rk,pot)
c 
        a(j,i)=pot
 1200 continue
 1400 continue
c 
c        construct the matrix of potentials of all nodes
c        of which the hole is comprized at the points in
c        the outer skeleton of the said hole
c 
        do 1800 i=1,nout
        do 1600 j=1,nhole
c 
        jj=izshole(j)
        ii=int_skelout(i)
c 
        if(ii .lt. 0) call potuser01(skelout(1,i),jj,rk,pot)
        if(ii .gt. 0) call potuser11(ii,jj,rk,pot)
c 
        c(j,i)=pot
 1600 continue
 1800 continue
c 
c        construct the evaluation matrix
c 
        ifcheck=0
        call ncleamatrr(a,nhole,nskel,nout,c,eval,
     1      eps,ncols,rnorms,w,ifcheck,errl2,errmax)
c 
c        Certain columns of the matrix eval should consist of
c        lots of zeroes and a single one. Enforce this rule
c 
        do 2600 i=1,nhole
        do 2400 j=1,nskel
c 
        jj=int_skelin(j)
        ii=izshole(i)
c 
        if(ii .ne. jj) goto 2400
c 
        do 2200 k=1,nskel
c 
        eval(i,k)=0
 2200 continue
c 
        eval(i,j)=1
 2400 continue
 2600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine expandmatr(skelin,int_skelin,
     1      skelout,int_skelout,nout,nskel,izshole,
     2      nhole,potuser10,potuser11,rk,eps,
     3      a,c,expand,w,rnorms)
        implicit real *8 (a-h,o-z)
        save
        complex *16 expand(nskel,nhole),a(nout,nskel),
     1      c(nout,nhole),w(1),rnorms(1),rk,pot
        dimension skelin(2,1),int_skelin(1),skelout(2,1),
cccc     1      int_skelout(1),zshole(2,1),izshole(1)
     1      int_skelout(1),izshole(1)
c 
c        This subroutine constructs the "expand" matrix for a
c        single hole
c 
c                         Input parameters:
c 
c  skelin - the physical locations in R^2 of the nodes comprizing
c        the final inner skeleton of the hole being processed
c  int_skelin - the locations in the input array z of the nodes
c        comprizing final the inner skeleton of the hole being
c        processed
c  skelout - the physical locations in R^2 of the nodes comprizing
c        the preliminary outer skeleton of the hole being processed
c  int_skelout - the locations in the input array z of the nodes
c        comprizing the preliminary outer skeleton of the hole
c        being processed
c  nout - the number of elements in the outer skeleton
c  nskel - the number of nodes in the inner skeleton
c  izshole - the locations in array z of all nodes comprizing the
c        hole
c  nhole - the number of elements in array izshole; also the number
c        of nodes of which the hole consists
c  potuser10, potuser11 - user-supplied subroutines defining the
c        various types of interactions; see subroutine r2d_allmatrs
c        for a detailed dscription of these.
c  rk - the Helmholtz coefficient
c  eps - the accuracy to which the SVDs will be performed
C 
C                         OUTPUT PARAMETERS:
C 
c  expand - the matrix converting an incoming potential on the
c        skeleton of the hole into the incoming potential
c        on all nodes constituting the hole
c 
c                         Work arrays:
c 
c  a,c,w,rnorms
c 
c 
c        . . . construct the matrix of interactions of the inner
c              skeleton of the hole with the outer skeleton of
c              the said hole
c 
        do 1400 i=1,nout
        do 1200 j=1,nskel
c 
        jj=int_skelin(j)
        ii=int_skelout(i)
c 
        if(ii .lt. 0) call potuser10(jj,skelout(1,i),rk,pot)
        if(ii .gt. 0) call potuser11(jj,ii,rk,pot)
c 
        a(i,j)=pot
 1200 continue
 1400 continue
c 
c        construct the matrix of potentials of all nodes
c        of which the hole is comprized at the points in
c        the outer skeleton of the said hole
c 
        do 1800 i=1,nout
        do 1600 j=1,nhole
c 
        jj=izshole(j)
        ii=int_skelout(i)
c 
        if(ii .lt. 0) call potuser10(jj,skelout(1,i),rk,pot)
        if(ii .gt. 0) call potuser11(jj,ii,rk,pot)
c 
        c(i,j)=pot
 1600 continue
 1800 continue
c 
c        construct the expansion matrix
c 
        ifcheck=0
        call ncleamatll(a,nout,nskel,nhole,c,expand,
     1      eps,ncols,rnorms,w,ifcheck,errl2,errmax)
c 
c        Certain columns of the matrix expand should consist of
c        lots of zeroes and a single one. Enforce this rule
c 
        do 2600 j=1,nskel
        do 2400 i=1,nhole
c 
        jj=int_skelin(j)
        ii=izshole(i)
c 
        if(ii .ne. jj) goto 2400
c 
        do 2200 k=1,nskel
c 
        expand(k,i)=0
 2200 continue
c 
        expand(j,i)=1
 2400 continue
 2600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine fin_skel(izsin,zsin,nin,
     1      izsout,zsout,nout,irows,icols,
     2      rk,eps,amatr,skelin,skelout,ncols,
     3      int_skelin,int_skelout,
     4      potuser01,potuser10,potuser11,z,w,
     5      a,a_outin,s_cond,rnorms,a_inout,rnorms3)
c 
        implicit real *8 (a-h,o-z)
        save
        complex *16 rk,pot,amatr(nin,nout*4),a(nin,nin),cond,
     1      a_outin(nin,nout),w(nin,1),a_inout(nout,nin)
        integer *4 irows(1),icols(1),izsin(1),izsout(1),
     1      int_skelin(1),int_skelout(1)
        dimension zsout(2,1),zsin(2,1),skelin(2,1),
     1      skelout(2,1),z(2,1),rnorms(1),rnorms3(1)
c 
c        This subroutine produces final skeletons (both inner and
c        outer) for one hole. The skeletons are returned both as
c        nodes in R^2 and as maps of their locations in array z of
c        user-supplied nodes.
c 
c                           Input parameters:
c 
c  izsin - the list of (consecutive) sequence numbers in array
c        z of the nodes comprizing the hole
c  zsin - the array containing the physical locations in the plane
c        of the nodes comprizing the hole
c  nin - the number of nodes in the hole (same as nhole)
c  izsout - the locations in array z of the sequence numbers of
c        elements in the outer preliminary skeleton (please note
c        that NOT all of these are actual nodes provided by the
c        user; the addresses in z of those elements that are not
c        present there do not make too much sense, and have been
c        set to -13)
c  zsout - the physical locations of points in R^2 comprizing the
c        preliminary outer skeleton
c  nout - the number of elements in zsout, izsout
c  rk - the Helmholtz coefficient
c  eps - the accuracy to which the SVDs will be performed
c  potuser10, potuser01, potuser11 - user-supplied subroutines
c        defining the various types of interactions needed by
c        this subroutine. Following are the decriptions of the
c        three
c 
c  potuser10 - the interaction between a (source) node on some apriori
c        discretized boundary, specified by its sequence number among
c        such nodes, and a target that is an point in R^2 that is not
c        a part of any apriori diacretized boundary. The
c        calling sequence is
c 
c                   potuser10(i1,z2,rk,pot)
c 
c        with
c 
c     i1 - the sequence number of the source (input)
c     z2 - the location of the target (input)
c     rk - the Helmholtz coefficient  (input)
c     pot - the potential created by the source on the target (output)
c 
c  potuser01 - the interaction between a (source) point in R^2 that
c        is not a part of any apriori diacretized boundary, and a
c        (target) node on some apriori discretized boundary, specified
c        by its sequence number among such nodes. The calling sequence
c        is
c 
c                   potuser10(z1,i2,rk,pot)
c 
c        with
c 
c     z1 - the location of the source (input)
c     z2 - the sequence number of the target (input)
c     rk - the Helmholtz coefficient  (input)
c     pot - the potential created by the source on the target (output)
c 
c  potuser11 - the interaction between a (source) node on some apriori
c        discretized boundary, specified by its sequence number among
c        such nodes, and a (target) node on some apriori discretized
c        boundary, also specified by its sequence number among such
c        nodes. The calling sequence is
c 
c                   potuser11(i1,i2,rk,pot)
c 
c        with
c 
c     z1 - the sequence number of the source (input)
c     z2 - the sequence number of the target (input)
c     rk - the Helmholtz coefficient  (input)
c     pot - the potential created by the source on the target (output)
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c  z - the points in R^2 on top of which the tree structure
c        is placed.
c  n - the number of points in array z
c  i1 - the address in z of the first element on the hole to be
c        skeletonized
c  i2 - the address in z of the last element on the hole to be
c        skeletonized
c 
c                        Output parameters:
c 
c  skelin - the physical locations of points in R^2 comprizing the
c        final inner skeleton
c  skelout - the physical locations of points in R^2 comprizing the
c        final outer skeleton
c  ncols - the number of elements in skelin, skelout, int_skelin,
c  int_skelin - the locations in array z of the sequence numbers of
c        elements in the inner final skeleton (please note that all
c        of these are actual nodes provided by the user)
c  int_skelout - the locations in array z of the sequence numbers of
c        elements in the outer final skeleton (please note that NOT
c        all of these are actual nodes provided by the user; the
c        addresses in z of those elements that are not present there
c        do not make too much sense, and are set to -13)
c        int_skelout
c  s_cond - the condition number of the matrix of (user-supplied)
c        interactions between the elements of the hole
c 
c                        Work arrays:
c 
c  w,a,a_outin,rnorms,rnorms3
c 
c        . . . construct the original scattering matrix
c              (on the uncompressed hole)
c 
        do 1040 i=1,nin
        do 1020 j=1,nin
c 
        ii=izsin(i)
        jj=izsin(j)
        call potuser11(jj,ii,rk,pot)
c 
        a(i,j)=pot
 1020 continue
 1040 continue
c 
c       find the maximim singular value of a
c 
        call max_singval(a,nin,nin,svmax,w,w(1,3) )
c 
        call corthom(a,nin,w,cond)
c 
c       find the maximum singular value of a^{-1}
c 
        call max_singval(a,nin,nin,svmin,w,w(1,3) )
c 
        s_cond=svmax*svmin
c 
c        construct a_outin
c 
        do 1140 i=1,nout
        do 1120 j=1,nin
c 
        jj=izsin(j)
        ii=izsout(i)
c 
        if(ii .lt. 0) call potuser01(zsout(1,i),jj,rk,pot)
c 
        if(ii .gt. 0) call potuser11(ii,jj,rk,pot)
        a_outin(j,i)=pot
c 
 1120 continue
 1140 continue
c 
c        evaluate the product a * a_outin
c 
        call ncleamult(a,a_outin,w,nin,nin,nout)
c 
c        QR the stupid thing
c 
        eps2=1.0d-12
        ifpivot=1
        call cleasgrm(w,nin,nout,rnorms,eps,ncols,ifpivot)
c 
c        construct a_inout
c 
        do 1240 i=1,nout
        do 1220 j=1,nin
c 
        jj=izsin(j)
        ii=izsout(i)
c 
        if(ii .lt. 0) call potuser10(jj,zsout(1,i),rk,pot)
c 
        if(ii .gt. 0) call potuser11(jj,ii,rk,pot)
c 
        a_inout(i,j)=pot
c 
 1220 continue
 1240 continue
c 
c        evaluate the product a_inout * a
c 
        call ncleamult(a_inout,a,amatr,nout,nin,nin)
        call cleastra(amatr,nout,nin,a_outin)
c 
c        QR the stupid thing
c 
        eps2=1.0d-12
        ifpivot=1
        call cleasgrm(a_outin,nin,nout,rnorms3,eps,ncols3,ifpivot)
c 
c        . . . construct the part of the matrix of interactions
c              that does not contain the interactions inside
c 
c 
        do 1800 i=1,nout
        do 1600 j=1,nin
c 
        jj=izsin(j)
        ii=izsout(i)
        if(ii .lt. 0) call potuser10(jj,zsout(1,i),rk,pot)
c 
        if(ii .gt. 0) call potuser11(jj,ii,rk,pot)
c 
        amatr(j,i)=conjg(pot)
c 
        if(ii .lt. 0) call potuser01(zsout(1,i),jj,rk,pot)
c 
        if(ii .gt. 0) call potuser11(ii,jj,rk,pot)
        amatr(j,i+nout)=conjg(pot)
c 
 1600 continue
 1800 continue
c 
c        construct the part of the matrix of interactions
c        involving the self-interactions of the hole
c 
        do 2400 i=1,ncols
        do 2200 j=1,nin
c 
        amatr(j,i+nout*2)=w(j,i)*rnorms(i)/svmin
 2200 continue
 2400 continue
c 
        do 3400 i=1,ncols3
        do 3200 j=1,nin
c 
        amatr(j,i+nout*2+ncols)=a_outin(j,i)*rnorms3(i)/svmin
 3200 continue
 3400 continue
c 
        nnn=nout*2+ncols+ncols3
c 
c        . . . skeletonize
c 
        call cskel_basic(amatr,nin,nnn,eps/100,
     1      irows,icols,ncols2,w)
c 
        ncols=ncols2
c 
        call sortanyn(irows,ncols,w)
        call sortanyn(icols,ncols,w)
c 
c       select the skeletons (outer and inner) out of arrays
c       zsin, zsout
c 
        do 5600 i=1,ncols
c 
        j=irows(i)
        skelin(1,i)=zsin(1,j)
        skelin(2,i)=zsin(2,j)
c 
        int_skelin(i)=izsin(j)
c 
        j=icols(i)
        jj=j
        if(jj .gt. nout) jj=jj-nout
c 
        skelout(1,i)=zsout(1,jj)
        skelout(2,i)=zsout(2,jj)
        int_skelout(i)=izsout(jj)
 5600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine max_singval(a,n,m,svmax,x,y)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,m),x(m),y(n),cd
c 
c         initialize the power method
c 
        numit=100
        done=1
        do 1200 i=1,m
        x(i)=sqrt(done*i)
 1200 continue
c 
        rn=1
        do 3000 k=1,numit
c 
cccc        call prinf('k=*',k,1)
cccc        call prin2('and x=*',x,m*2)
c 
c       normalize x
c 
        d=0
        do 1300 i=1,m
c 
        d=d+x(i)*conjg(x(i))+d
 1300 continue
c 
        rn=sqrt(d)
c 
        if(k .eq. numit) goto 3200
c 
        d=1/rn
        do 1350 i=1,m
c 
        x(i)=x(i)*d
 1350 continue
c 
c        apply a to x
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
c        apply a^* to y
c 
        do 2600 i=1,m
c 
        cd=0
        do 2400 j=1,n
c 
        cd=cd+conjg(a(j,i))*y(j)
 2400 continue
c 
        x(i)=cd
 2600 continue
c 
 3000 continue
c 
 3200 continue
c 
        svmax=sqrt(rn)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine prelim_skel(zshole,nhole,zsothers,nothers,
     1      coedist,nframe,zsaux,inds,nout,frame,iplot,rk)
        implicit real *8 (a-h,o-z)
        save
        dimension zshole(2,1),zsothers(2,1),frame(2,1),
     1      zsaux(2,1),cent(3),inds(1)
c 
c        This subroutine produces a preliminary outer skeleton
c        for one hole, the said skeleton (together with the
c        nodes inside the said hole) to be further compressed
c        by the subroutinef in_skel (see).
c 
c                           Input parameters:
c 
c  zshole - the points discretizing the hole
c  nhole - the number of points in the hole (number of elements in
c        array zshole)
c  zsothers - the points in the simulations, MINUS THE POINTS
C        IN THE HOLE.
C  nothers - the number of elements in array zsothers
c  coedist - the ratio of the radius of the artificial (circular)
c        boundary to the radius of the hole; a good number is 2
c  nframe - the number of equispaced nodes into which the artificial
c        boundary will be discretized
c  zsaux - the physical locations of points in R^2 comprizing the
c        preliminary outer skeleton
c  inds - the locations in array iz of the sequence numbers of
c        elements in the outer preliminary skeleton (please note
c        that NOT all of these are actual nodes provided by the
c        user; the addresses in iz of those elements that are not
c        present there do not make too much sense, and have been
c        set to -13)
c  nout - the number of elements in zsout
c 
c                        Work arrys:
c 
c  frame - must be at least nframe*w2 real *8 elements long
c 
c 
c       . . . find the center of mass of the hole
c 
        cent(1)=0
        cent(2)=0
c 
        do 1200 i=1,nhole
c 
        cent(1)=cent(1)+zshole(1,i)
        cent(2)=cent(2)+zshole(2,i)
c 
 1200 continue
c 
        cent(1)=cent(1)/nhole
        cent(2)=cent(2)/nhole
c 
c        calculate the radius of the hole
c 
        rad22=0
        do 1400 i=1,nhole
c 
        d=(zshole(1,i)-cent(1))**2+(zshole(2,i)-cent(2))**2
c 
        if(d .gt. rad22) rad22=d
 1400 continue
c 
        rad=sqrt(rad22)
c 
c       construct the list of nodes that are within radius*coedist
c       from the center
c 
        nindisk=0
        rc22=(rad*coedist)**2
c 
        do 1600 i=1,nothers
c 
        d=(zsothers(1,i)-cent(1))**2+(zsothers(2,i)-cent(2))**2
c 
        if(d .gt. rc22) goto 1600
c 
        nindisk=nindisk+1
        inds(nindisk)=i
c 
        zsaux(1,nindisk)=zsothers(1,i)
        zsaux(2,nindisk)=zsothers(2,i)
c 
 1600 continue
c 
        call prin2('in prelim_scale, rk=*',rk,2)
c 
c        Construct the circular outer frame
c 
        done=1
        pi=atan(done)*4
        rframe=rad*coedist
c 
        nframe2=nframe/2
        h=2*pi/nframe2
        nout=nindisk
c 
c        . . . the first layer
c 
        do 1800 i=1,nframe2
c 
        d=(i-1)*h
        frame(1,i)=cos(d)*rframe+cent(1)
        frame(2,i)=sin(d)*rframe+cent(2)
c 
        nout=nout+1
        zsaux(1,nout)=frame(1,i)
        zsaux(2,nout)=frame(2,i)
        inds(nout)=-13
 1800 continue
c 
c        . . . the second layer
c 
        rframe=rframe+1/abs(rk)/3
c 
        nframe3=nframe-nframe2
c 
        do 2800 i=1,nframe3
c 
        d=(i-1)*h
        frame(1,i+nframe3)=cos(d)*rframe+cent(1)
        frame(2,i+nframe3)=sin(d)*rframe+cent(2)
c 
        nout=nout+1
        zsaux(1,nout)=frame(1,i+nframe3)
        zsaux(2,nout)=frame(2,i+nframe3)
        inds(nout)=-13
 2800 continue
c 
c 
c       plot them things, if the user so requested
c 
        if(iplot .eq. 0) return
c 
        call zquaplot3(iplot,zshole,nhole,2,zsaux,nindisk,2,
     1      frame,nframe,2,
     2    'preliminary skeletons inside prelim_skel*')
  
        return
        end
c 
c 
c 
c 
c 
        subroutine cmatvec_adj(a,n,x,y)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,n),x(1),y(1),cd
c 
c        apply the matrix a to the vector x obtaining y
c 
        do 1400 i=1,n
        cd=0
        do 1200 j=1,n
cccc        cd=cd+a(i,j)*x(j)
        cd=cd+conjg(a(j,i))*x(j)
 1200 continue
        y(i)=cd
 1400 continue
        return
        end
