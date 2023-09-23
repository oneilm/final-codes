        implicit real *8 (a-h,o-z)
        integer *4 iz(1000 000),box(12)
c 
        real *8 z(2,1000 000),w(100 000 000),corners(2,5),
     2      center(3),zsout(2,10000),zsin(2,10000),
     3      zssons(2,10000),zsdad(2,10000),
     4      errors_expand(10000),errors_eval(10000),
     5      conds(1000)
c 
        complex *16 rk
c 
        external poteval2
        external potuser
  
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
C 
C       SET ALL PARAMETERS
C 
         PRINT *, 'ENTER nsmall'
         READ *,nsmall
         CALL PRINF('nsmall=*',nsmall,1 )
c 
         eps=1.0d-12
         rk=2
         rk=300
  
         rk=20
         rk=100
         rk=600
         rk=20
  
  
         call potuser_init(rk)
  
c 
c       construct and plot the test points
c 
  
        call RANLST(NSMALL,Z,W,N)
cccc          call prinf('after ranlst, n=*',n,1)
        a=10
        b=1
c 
cccc        call creelips(z,n,a,b)
cccc        n=n-1
  
cccc        call prinf('after ranlst, n=*',n,1)
  
        iw=21
        call zquaplot(iw,z,n,2,'random points as created*')
  
  
        lenw=100 000 000
        maxboxes=100000
        nbox=30
  
c 
        call r2d_allmatrs(ier,potuser,poteval2,eps,rk,z,n,
     1      maxboxes,nbox,w,lenw,nall,iz,nlev,
     2      conds,nconds,lused,keep)
c 
        call prin2('after r2d_allmatrs, conds=*',conds,nconds)
        call prinf('after r2d_allmatrs, lused/1000=*',
     1      lused/1000,1)
        call prinf('after r2d_allmatrs, keep/1000=*',
     1      keep/1000,1)
c 
cccc        istore=w(9)
  
        call prinf('nall=*',nall,1)
  
c 
c        plot the whole structure
c 
  
        iw=25
        call allplot3(iw,w,z)
c 
        call irows_icols_test(z,iz,w,nall)
c 
c       one box after another, retrieve their inner and outer
c       skeletons and plot them
c 
        do 3400 ibox=2,nall
c 
c 
        call r2dretr(ibox,w,box,corners,center)
c 
        if(box(4) .eq. 0) goto 3400
c 
        ifinner=0
        call r2d_finskels_retr(w,ibox,ifinner,zsout,lenout)
c 
        ifinner=1
        call r2d_finskels_retr(w,ibox,ifinner,zsin,lenin)
c 
cccc        call prinf('ibox=*',ibox,1)
cccc        call prinf('and lenin=*',lenin,1)
c 
        corners(1,5)=corners(1,1)
        corners(2,5)=corners(2,1)
c 
        iw=700+ibox
        call zquaplot3(iw,zsin,lenin,2,zsout,lenout,2,
     1      corners,5,3,'final skeletons*')
c 
 3400 continue
c 
c        test the scattering matrices
c 
        call scat_matr_test(w,nall,
     1    iz,z,rk,potuser,nlev)
  
        call prin2('and conds=*',conds,nconds)
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine scat_matrs_comp(s,irows_outin,ncols_outin,
     1      icols_inout,ncols_inout,
     2      sdir,np,iz,z,rk,potuser,w,ibox,error)
c 
        implicit real *8 (a-h,o-z)
        save
        integer *4 irows_outin(1),icols_inout(1),iz(1),box(12)
        dimension z(2,1),skelin(2,10000),skelout(2,10000),
     1      targ(3),zs(2,10000),w(1),corners(2,5),center(3)
        complex *16 rk,sdir(np,np),s(ncols_inout,ncols_outin),
     1      pot,potsin(10000),chgs(10000),fld,potsin2(10000),
     2      chgs2(10000),fld2,work(100000)
c 
        targ(1)=100
        targ(2)=100
  
cccc        targ(1)=z(1,1)
cccc        targ(2)=z(2,1)
c 
c        construct the basic scattering matrix (based on actual nodes)
c 
        call r2dretr(ibox,w,box,corners,center)
  
        istart=box(6)
  
        call r2d_scatmatr0(istart,np,potuser,par1,par2,
     1      sdir,conddir,work)
  
  
        call prin2('after r2d_scatmatr0, conddir=*',conddir,1)
  
  
c 
c        construct the physical locations of the skeletons
c        (incoming and outgoing) inside the box
c 
c        . . . incoming
c 
  
        call prinf('in scat_matrs_comp, ncols_outin=*',ncols_outin,1)
  
        do 1200 i=1,ncols_outin
c 
        j=irows_outin(i)
        jj=iz(j)
  
ccc        call
c 
        skelin(1,i)=z(1,jj)
        skelin(2,i)=z(2,jj)
 1200 continue
c 
cccc        call prin2('skelin=*',skelin,ncols_outin*2)
c 
c        . . . outgoing
c 
        do 1400 i=1,ncols_inout
c 
        j=icols_inout(i)
        jj=iz(j)
c 
        skelout(1,i)=z(1,jj)
        skelout(2,i)=z(2,jj)
 1400 continue
c 
cccc        call prin2('skelout=*',skelout,ncols_inout*2)
c 
c        construct the incoming potentials at the incoming skeleton
c 
        do 1600 i=1,ncols_outin
c 
        call poteval2(skelin(1,i),targ,rk,pot)
c 
        potsin(i)=pot
c 
 1600 continue
c 
cccc        call prin2('potsin as created*',potsin,ncols_outin*2)
c 
c       convert them things into a charge distribution at the
c       elements of the outgoing skeleton
c 
        call matvec778(s,ncols_inout,ncols_outin,potsin,chgs)
c 
ccccc        call prin2('and chgs=*',chgs,ncols_inout*2)
c 
c 
c        evaluate the outgoing potential at the target point
c 
        call flds_eval(skelout,ncols_inout,chgs,targ,
     1      rk,fld)
  
        call prin2('and induced potential on target point*',
     1      fld,2)
c 
c       construct the physical locations of nodes inside this box
c 
        do 2200 i=1,np
c 
        j=i+istart-1
        jj=iz(j)
c 
        zs(1,i)=z(1,jj)
        zs(2,i)=z(2,jj)
 2200 continue
c 
cccc        call prin2('and zs=*',zs,np*2)
c 
c        construct the incoming potentials at the nodes
c 
        do 2600 i=1,np
c 
        call poteval2(zs(1,i),targ,rk,pot)
c 
        potsin2(i)=pot
c 
 2600 continue
c 
cccc        call prin2('potsin2 as created*',potsin2,np*2)
c 
c        construct the induced charges at the nodes
c 
c 
        call matvec778(sdir,np,np,potsin2,chgs2)
  
cccc        call prin2('and chgs2=*',chgs2,np*2)
c 
c 
c        finally, evaluate the outgoing potential at the target point
c 
        call flds_eval(zs,np,chgs2,targ,rk,fld2)
  
        call prin2('and fld2=*',fld2,2)
        call prin2('while fld=*',fld,2)
        call prin2('and fld2-fld=*',fld2-fld,2)
c 
        error=abs(fld2-fld)
  
        return
        end
  
  
  
c 
c 
c 
c 
c 
        subroutine flds_eval(sources,n,sigmas,targ,
     1      rk,fld)
        implicit real *8 (a-h,o-z)
        save
        dimension sources(2,1),targ(2)
c 
        complex *16 rk,fld,sigmas(1),pot
c 
c        evaluate the potential of a charge distribution at the
c        point targ
c 
        fld=0
        do 1200 i=1,n
  
cccc        call prinf('in flds_eval, i=*',i,1)
  
        call poteval2(sources(1,i),targ,rk,pot)
  
cccc        call prin2('in flds_eval, pot=*',pot,2)
c 
        fld=fld+sigmas(i)*pot
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine poteval(i1,i2,potuser,par1,par2,pot)
        implicit real *8 (a-h,o-z)
        save
        complex *16 pot
        dimension par1(1),par2(1),iz(100 000),iz7(1),
     1      zs(2,100 000),zs7(2,1)
c 
cccc        call prinf('before potuser, i1=*',i1,1)
cccc        call prinf('before potuser, i2=*',i2,1)
  
        i=iz(i1)
        j=iz(i2)
c 
  
        call potuser(i,j,zs(1,i),zs(1,j),par1,par2,pot)
c 
        return
c 
c 
c 
c 
        entry poteval_init(iz7,zs7,n7)
c 
        n=n7
        call r2d_reacopy(zs7,zs,n*2)
        call r2d_intcopy(iz7,iz,n)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine potuser(j1,j2,z1,z2,par1,par2,pot)
        implicit real *8 (a-h,o-z)
        save
        complex *16 pot,rk,z,h0,h1
        dimension z1(2),z2(2),par1(1),par2(1)
c 
cccc        call prin2('in potuser, z1=*',z1,2)
cccc        call prin2('in potuser, z2=*',z2,2)
  
cccc        call prinf('in potuser, j1=*',j1,1)
cccc        call prinf('in potuser, j2=*',j2,1)
  
c 
        if(j1 .ne. j2) goto 1400
c 
        pot=-1000 000
        pot=-1
c 
        return
 1400 continue
c 
        d=(z1(1)-z2(1))**2+(z1(2)-z2(2))**2
        d=sqrt(d)
        z=d*rk
c 
        ifexpon=1
        call hank103(z,h0,h1,ifexpon)
c 
        pot=h0
cccc        pot=h0*z1(1)
c 
        return
c 
c 
c 
c 
        entry potuser_init(rk7)
        rk=rk7
        return
        end
  
c 
c 
c 
c 
c 
        subroutine irows_icols_test(z,iz,w,nall)
        implicit real *8 (a-h,o-z)
        save
        integer *4 box(12),iz(1)
        dimension w(1),zs(2,10000),center(3),store(1),z(2,1),
     1      corners(2,5),icols(10000),zs2(2,1000)
c 
c        one box after another, find locations in array iz of
c        the elements of skeletons of said boxes
c 
        error=0
        do 2000 ibox=1,nall
c 
c        if this is a childless box - skip
c 
        call r2dretr(ibox,w,box,corners,center)
c 
        if(box(9) .eq. -7) goto 2000
        if(box(4) .eq. 0) goto 2000
c 
c        retrieve from array store the physical locations
c        of nodes in the skeleton
c 
        ifarr=1
        kind=14
        call r2d_finretr(ibox,kind,w,
     1      ncols_inout,ncols_outin,np,zs,len,ifarr)
c 
cccc        call prinf('in irows_icols_test, ibox =*',ibox,1)
cccc        call prin2('zs as retrieved from array store*',zs,len*2)
c 
c        retrieve from array store the locations of the skeleton
c        nodes in array iz
c 
        ifarr=1
        kind=4
        call r2d_finretr(ibox,kind,w,
     1      ncols_inout,ncols_outin,np,icols,len2,ifarr)
c 
cccc        call prinf('and icols as retrieved=*',icols,len2)
c 
c       convert the locations of the nodes in array iz into
c       their physical locations, and compare these with
c       the physical locations extracted from array store
c 
        do 1400 i=1,len2
c 
        j=iz(icols(i))
        zs2(1,i)=z(1,j)
        zs2(2,i)=z(2,j)
 1400 continue
c 
        do 1600 i=1,len2
c 
        zs2(1,i)=zs2(1,i)-zs(1,i)
        zs2(2,i)=zs2(2,i)-zs(2,i)
c 
        error=error+zs2(1,i)**2+zs2(2,i)**2
  
 1600 continue
c 
 2000 continue
c 
        call prin2('total error from irows_icols_test is*',
     1      error,1)
  
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
        subroutine allplot3(iw,w,z)
        implicit real *8 (a-h,o-z)
        save
        integer *4 box(10),labs(100000)
        dimension corners(2,4),center(2),w(1),z(2,1),
     1      w2(100000),centers(2,100000)
c 
c        initialize the plotting routine
c 
         call trplopen(w2)
c 
c        unpack the beginning of array w
c 
        nboxes=w(5)
        n=w(6)
c 
c        plot all boxes
c 
        nlabels=0
        do 2000 i=1,nboxes
c 
        call r2dretr(i,w,box,corners,center)
c 
        if(box(9) .eq. -7) goto 2000
c 
        call plot1box(corners,w2)
  
        call trplpnt(center(1),center(2),w2)
c 
        nlabels=nlabels+1
c 
        centers(1,nlabels)=center(1)
        centers(2,nlabels)=center(2)
c 
        labs(nlabels)=i
 2000 continue
c 
c        create the labels file
c 
        iw7=iw
        call zqua_labels_int(iw7,centers,labs,nlabels)
c 
        call prinf('in allplot3, iw=*',iw,1)
c 
c        plot all nodes
c 
        do 2200 i=1,n
c 
        call trpldot(z(1,i),z(2,i),w2)
 2200 continue
c 
        call trplwrt(iw,w2,'all boxes*')
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine neighb_print(nall,neighbors)
        implicit real *8 (a-h,o-z)
        save
        dimension neighbors(11,1)
c 
        call prinf('and neighbors are*',neighbors,0)
c 
 5600 format(12(1x,i5))
        do 5800 i=1,nall
c 
        write(6,5600) i,(neighbors(j,i),j=1,11)
        write(13,5600) i,(neighbors(j,i),j=1,11)
 5800 continue
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine plot1box(corn,w)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),corn(2,4)
c 
        call trplline(corn(1,1),corn(2,1),
     1      corn(1,2),corn(2,2),w)
c 
        call trplline(corn(1,2),corn(2,2),
     1      corn(1,3),corn(2,3),w)
c 
        call trplline(corn(1,3),corn(2,3),
     1      corn(1,4),corn(2,4),w)
c 
        call trplline(corn(1,4),corn(2,4),
     1      corn(1,1),corn(2,1),w)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine creelips(z,n,a,b)
        implicit real *8 (a-h,o-z)
        save
        dimension z(2,1)
c 
        done=1
        pi=datan(done)*4
        h=2*pi/n
        do 1200 i=1,n
        t=(i-1)*h
        z(1,i)=a*dcos(t)
        z(2,i)=b*dsin(t)
 1200 continue
c 
c       rotate the stupid thing
c 
        theta=1.1
  
cccc        theta=3.14/2
  
cccc         theta=theta/3
ccc        theta=pi/2
        a11=dcos(theta)
        a12=dsin(theta)
        a22=a11
        a21=-a12
        do 1400 i=1,n
        x=a11*z(1,i)+a12*z(2,i)
        y=a21*z(1,i)+a22*z(2,i)
        z(1,i)=x
        z(2,i)=y
 1400 continue
  
        return
        end
c 
C 
C 
C 
C 
        SUBROUTINE RANLST(NSMALL,Z,W,N)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        DIMENSION Z(2,1),CENT(2),W(1)
C 
C       CONSTRUCT FOUR CORNERS TO DEFINE THE BOX
C 
        Z(1,1)=0
        Z(2,1)=1
C 
        Z(1,2)=0
        Z(2,2)=0
C 
        Z(1,3)=1
        Z(2,3)=0
C 
        Z(1,4)=1
        Z(2,4)=1
C 
        NSTART=5
        DONE=1
C 
C       CONSTRUCT BOX 1
C 
       N=NSMALL
       SIZE=0.5
       CENT(1)=0.75
       CENT(2)=0.75
       CALL RIRANB(N,SIZE,CENT,Z(1,NSTART),W)
CCC    CALL PRIN('AFTER RIRANB, Z(1,NSTART)=*',Z(1,NSTART),20)
       NSTART=NSTART+NSMALL
C 
C       CONSTRUCT BOX 2
C 
       SIZE=DONE/4
       CENT(1)=DONE/2+DONE/4+DONE/8
       CENT(2)=DONE/2-DONE/8
       CALL RIRANB(N,SIZE,CENT,Z(1,NSTART),W)
       NSTART=NSTART+NSMALL
C 
C       CONSTRUCT BOX 3
C 
       SIZE=DONE/4
       CENT(1)=DONE/2+DONE/4+DONE/8
       CENT(2)=DONE/8
       CALL RIRANB(N,SIZE,CENT,Z(1,NSTART),W)
       NSTART=NSTART+NSMALL
C 
C       CONSTRUCT BOX 4
C 
       SIZE=DONE/8
       CENT(1)=1-DONE/16
       CENT(2)=DONE/2-DONE/8-DONE/16
       CALL RIRANB(N,SIZE,CENT,Z(1,NSTART),W)
       NSTART=NSTART+NSMALL
C 
C       CONSTRUCT BOX 5
C 
       SIZE=DONE/8
       CENT(1)=DONE/2+DONE/16
       CENT(2)=DONE/2-DONE/16
       CALL RIRANB(N,SIZE,CENT,Z(1,NSTART),W)
       NSTART=NSTART+NSMALL
C 
C       CONSTRUCT BOX 6
C 
       SIZE=DONE/16
       CENT(1)=DONE/2+DONE/16+DONE/32
       CENT(2)=DONE/2-DONE/16-DONE/32
       CALL RIRANB(N,SIZE,CENT,Z(1,NSTART),W)
       NSTART=NSTART+NSMALL
cccc       N=NSTART-1
        n=nstart-5
       RETURN
       END
C 
C 
C 
C 
C 
        SUBROUTINE RIRANB(n,SIZE,CENT,Z,W)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        DIMENSION CENT(2),Z(2,1),W(1)
C 
C       GENERATE PSEUDO-RANDOM NUMBERS
C 
        CALL RIRAND(N*2,W)
        DO 1200 I=1,N
        Z(1,I)=W(I)
        Z(2,I)=W(I+N)
 1200 CONTINUE
C 
C       SHIFT AND STRETCH/COMPRESS THE POINTS IN THE PLANE
C 
        DO 1400 I=1,N
        Z(1,I)=Z(1,I)*SIZE+CENT(1)-0.5*SIZE
        Z(2,I)=Z(2,I)*SIZE+CENT(2)-0.5*SIZE
 1400 CONTINUE
        RETURN
        END
C 
C 
C 
C 
C 
        SUBROUTINE RIRAND(N,Y)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        DIMENSION Y(1)
        DATA IFCALL /0/
C 
C       CONDUCT PRELIMINARY RANDOMIZATION
C 
        IF (IFCALL.EQ.1 ) GOTO  1300
        IFCALL =1
CCCC    LAMBDA=13
        LAMBDA=2**7 +9
        MU=510511
        MU=1
CCC     IP=2**10
        IP=2**20
        D=1
        D=D/IP
        M=17
        DO 1200 I=1,100
        M1=M*LAMBDA  +MU
        J=M1/IP
        M=M1-J*IP
CCCC    CALL PRINF('M=*',M,1)
 1200 CONTINUE
CCC     CALL PRINF('PRELIMINARY RANDOMIZATION DONE, M=*',M,1)
 1300 CONTINUE
C 
C       GENERATE PSEUDO-RANDOM NUMBERS
C 
        DO 1400 I=1,N
        M1=M*LAMBDA +MU
        J=M1/IP
        M=M1-J*IP
        Y(I)=M*D
 1400 CONTINUE
        RETURN
        END
c 
c 
c 
c 
c 
        subroutine plot4(iw,w,ibox,list,nlist,
     1      iz,inband,ninband,zs)
        implicit real *8 (a-h,o-z)
        save
        integer *4 box(12),labs(100000),iz(1),inband(1)
        dimension corners(2,4),center(2),w(1),list(1),
     1      w2(100000),centers(2,10000),zs(2,1)
c 
c        initialize the plotting routine
c 
         call trplopen(w2)
c 
c        plot the big box
c 
        call r2dretr(ibox,w,box,corners,center)
c 
        call plot1box(corners,w2)
c 
c       plot all boxes from the list
c 
        do 2000 i=1,nlist
c 
        j=list(i)
        call r2dretr(j,w,box,corners,centers(1,i))
c 
        call plot1box(corners,w2)
 2000 continue
c 
c        plot the nodes
c 
        do 2200 i=1,ninband
c 
        j=iz(inband(i))
cccc        j=iz(i)
        call trpldot(zs(1,j),zs(2,j),w2)
 2200 continue
c 
c        create the labels file
c 
        nlabels=1
        labs(1)=ibox
        call zqua_labels_int(iw,center,labs,nlabels)
  
c 
c        dump them things on file
c 
        call trplwrt(iw,w2,
     1      'box and separated kids of its neighbors*')
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine evalexpand_test(poteval,rk,ibox,w,
     1      eval,expand,zsdad,zssons,lendad,lensons,eps,
     2      z,n,errmax_expand,errmax_eval)
        implicit real *8 (a-h,o-z)
        save
        integer *4 box(12)
        dimension w(1),zsdad(2,1),zssons(2,1),center(3),
     1      corners(2,5),z(2,1),errs_expand(100000),
     2      errs_eval(100 000)
        complex *16 rk,eval(1),expand(1)
c 
c       scan the array of all nodes. For those that are
c       outside the daddy box, test the matrices eval,expand
c 
        call r2dretr(ibox,w,box,corners,center)
c 
        errmax_eval=0
        errmax_expand=0
c 
        do 2000 i=1,n
c 
        call if_inrectangle(corners,z(1,i),ifinside)
c 
        if(ifinside .eq. 1) then
c 
            errs_eval(i)=0
            errs_expand(i)=0
            goto 2000
        endif
c 
c       this point is outside the box ibox. test the
c       matrices eval, expand
c 
        call evalmat_test(poteval,rk,zsdad,zssons,
     1     lendad,lensons,eval,z(1,i),errmax)
c 
        errs_eval(i)=errmax
c 
        if(errmax_eval .lt. errmax) errmax_eval=errmax
  
        call expanmat_test(poteval,rk,zsdad,zssons,
     1     lendad,lensons,expand,z(1,i),errmax2)
c 
        errs_expand(i)=errmax2
c 
        if(errmax_expand .lt. errmax2) errmax_expand=errmax2
c 
 2000 continue
c 
cccc        call prin2('errs_expand=*',errs_expand,n)
cccc        call prin2('errs_eval=*',errs_eval,n)
c 
  
        return
        end
c 
c 
c 
c 
c 
        subroutine if_inrectangle(corn,z,ifinside)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension corn(2,5),z(2)
c 
        ifinside=1
        if(corn(1,1) .lt. z(1)) ifinside=0
        if(corn(1,2) .gt. z(1)) ifinside=0
c 
        if(corn(2,1) .lt. z(2)) ifinside=0
        if(corn(2,3) .gt. z(2)) ifinside=0
c 
  
        return
        end
c 
c 
c 
c 
c 
c 
        subroutine expanmat_test(poteval,rk,zsdad,zssons,
     1     lendad,lensons,expand,z0,errmax)
        save
  
        implicit real *8 (a-h,o-z)
        complex *16 expand(lensons,lendad),rk,pot,diffs(10000),
     1      valssons(1000),valssons2(1000),valsdad(1000)
c 
        dimension zsdad(2,1),zssons(2,1),z0(2)
c 
c        construct the vector of potentials generated at the
c        point z0 by all elements of the daddy's skeleton
c 
        do 1200 i=1,lendad
c 
        call poteval(z0,zsdad(1,i),rk,pot)
c 
        valsdad(i)=pot
 1200 continue
c 
c        construct the vector of potentials generated at the
c        point z0 by all elements of the kids' skeletons
c 
        do 1400 i=1,lensons
c 
        call poteval(z0,zssons(1,i),rk,pot)
c 
        valssons(i)=pot
 1400 continue
  
cccc        call prin2('and valssons=*',valssons,lensons*2)
c 
c        evaluate the potentials generated by the sonnies
c        via the matrix expand
c 
        call matvec778_tran(expand,lensons,lendad,valsdad,
     1      valssons2)
c 
cccc        call prin2('and valssons2=*',valssons2,lensons*2)
  
        errmax=0
  
        do 1600 i=1,lensons
c 
        diffs(i)=valssons2(i)-valssons(i)
  
        d=abs(diffs(i))
        if(d .gt. errmax) errmax=d
 1600 continue
c 
ccccc        call prin2('and diffs in expanmat_test=*',
cccc     1      diffs,lensons*2)
c 
  
  
cccc        call prin2('and errmax=*',errmax,1)
  
        return
        end
  
c 
c 
c 
c 
c 
        subroutine matvec778_tran(a,m,n,x,y)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,m),x(n),y(m),cd
c 
        do 1400 i=1,m
c 
        cd=0
        do 1200 j=1,n
c 
        cd=cd+a(j,i)*x(j)
 1200 continue
c 
        y(i)=cd
 1400 continue
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine evalmat_test(poteval,rk,zsdad,zssons,
     1     lendad,lensons,eval,z0,errmax)
        implicit real *8 (a-h,o-z)
        save
        complex *16 eval(lensons,lendad),rk,pot,diffs(10000),
     1      valssons(1000),valssons2(1000),valsdad(1000)
c 
        dimension zsdad(2,1),zssons(2,1),z0(2)
c 
c        construct the matrix of potentials generated by the
c        charge z0 at all elements of the daddy's skeleton
c 
        do 1200 i=1,lendad
c 
        call poteval(z0,zsdad(1,i),rk,pot)
c 
        valsdad(i)=pot
 1200 continue
c 
c        construct the matrix of potentials generated by the
c        charge z0 at all elements of the kids' skeletons
c 
        do 1400 i=1,lensons
c 
        call poteval(z0,zssons(1,i),rk,pot)
  
cccc        pot=1
c 
        valssons(i)=pot
 1400 continue
c 
cccc        call prin2('and valssons=*',valssons,lensons*2)
c 
c        evaluate the potential at the sonnies via the matrix eval
c 
        call matvec778(eval,lensons,lendad,valsdad,valssons2)
c 
cccc        call prin2('and valssons2=*',valssons2,lensons*2)
c 
        errmax=0
        do 1600 i=1,lensons
c 
        diffs(i)=valssons2(i)-valssons(i)
        d=abs(diffs(i))
        if (errmax .lt. d) errmax=d
 1600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine matvec778(a,m,n,x,y)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(m,n),x(n),y(m),cd
c 
        do 1400 i=1,m
c 
        cd=0
        do 1200 j=1,n
c 
        cd=cd+a(i,j)*x(j)
 1200 continue
c 
        y(i)=cd
 1400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine scat_matr_test(w,nall,
     1    iz,z,rk,potuser,nlev)
        implicit real *8 (a-h,o-z)
        save
        integer box(12),iz(1)
        real *8 z(2,1),w(1),center(3),
     1           corners(2,5),conds(1000),errors(100000)
        complex *16 rk,sout(100 000),
     1    sdir(100 000)
c 
c 
        integer *4 irows_outin(10000),icols_inout(10000)
  
        iiii=0
        do 4000 lev=nlev,4,-1
c 
        do 3600 ibox=1,nall
c 
  
cccc        call prinf('in scat_matr_test, ibox=*',ibox,1)
  
        call r2dretr(ibox,w,box,corners,center)
  
c 
        if(box(9) .eq. -7) goto 3600
        if(box(2) .ne. lev) goto 3600
        if(box(4) .eq. 0) goto 3600
        if(box(5) .eq. 0) goto 3600
c 
        iiii=iiii+1
c 
c        retrieve the scattering matrix from disk
c 
        ifarr=1
        kind=3
        call r2d_finretr(ibox,kind,w,
     1      ncols_inout7,ncols_outin7,np,sout,
     2      lensout,ifarr)
  
        ifarr=1
        kind=4
        call r2d_finretr(ibox,kind,w,
     1      ncols_inout7,ncols_outin7,np,irows_outin,
     2      ncols_outin,ifarr)
c 
        kind=5
        call r2d_finretr(ibox,kind,w,
     1      ncols_inout7,ncols_outin7,np,icols_inout,
     2      ncols_inout,ifarr)
c 
  
c 
        call r2dretr(ibox,w,box,corners,center)
        np=box(7)
c 
c        comparee the action of the two scattering matrices
c 
        call prinf('before scat_matr_comp, ibox=*',ibox,1)
c 
        call scat_matrs_comp(sout,irows_outin,ncols_outin,
     1      icols_inout,ncols_inout,
     2      sdir,np,iz,z,rk,potuser,w,ibox,error)
  
        call prinf('and by the way, lev=*',lev,1)
        call prin2('and by the way, error=*',error,1)
c 
        errors(iiii)=error
  
 3600 continue
c 
 4000 continue
  
        call prin2('and errors in all s-matrices*',errors,iiii)
  
        return
        end
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code and the beginning
c        of the direct solver code proper.
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine r2d_allmatrs(ier,potuser,poteval,eps,
     1      rk,z,n,maxboxes,nbox,w,lenw,nall,iz,
     2      nlev,conds,nconds,lused,keep)
        implicit real *8 (a-h,o-z)
        save
        integer iz(1)
        real *8 z(2,1),w(1),conds(1)
c 
        complex *16 rk
c 
c        First, this subroutine constructs the "bastardized"
c        binary tree, obtained from the standard quad-tree by
c        the introduction of "bastard" boxes, so that every
c        box is subdivided into two sons (one of which might
c        not exist). Thus, there are two types of boxes: square
c        and rectangular; the rectangular boxes are twice as
c        high as they are wide. Please note that square boxes
c        have even-numbered levels of subdivision: 0 (main
c        box), 2, 4,... . Rectangular boxes have odd-numbered
c        levels of subdivision. The obtained structure is accessed
c        via the subroutine r2dretr (see) returning to the user
c        various types of data associated with the user-specified
c        box in the structure.
c 
c        Second, this subroutine constructs all skeletons on all
c        levels. The said skeletons are also stored in the array
c        w, from which array they are retrieved via the subroutine
c        r2d_finskels_retr (see)
c 
c        Third, this subroutine constructs the sparse inverse of
c        the interaction matrix on the whole structure, in the
c        form of scattering, merging, splitting, and exchange
c        matrices. All of these matrices are stored in the array
c        w, from where they can be retrieved by calling the
c        subroutine r2d_finretr (see).
c 
C        PLEASE NOTE THAT WHILE THE SUBROUTINES R2D_FINSKELS_RETR,
C        R2D_FINRETR ARE (AS OF 10.25.02) A PART OF THIS FILE, THE
C        SUBROUTINE R2DRETR IS A PART OF THE FILE r2dstrcr.f.
C        ATTEMPTS TO FIND A SUBROUTINE IN A FILE WHERE IT IS NOT
C        ARE UNLIKELY TO BE SUCCESSFUL.
c 
c                           Input parameters:
c 
c  poteval - the user-supplied subroutine defining the interaction
c        between the nodes. The calling sequence of poteval is
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
c  eps - the accuracy to which the SVDs will be performed
c  rk - the Helmholtz coefficient
c  z - the points in R^2 on top of which the tree structure
c        is placed.
c  n - the number of points in array z
c  maxboxes - the maximum number of boxes the subroutine is
c        permitted to construct. If this number is insufficient,
c        the error return code ier (see below) is set to 4, and
c        the execution is terminated
c  nbox - the maximum number of points permitted in a childless
c        box
c  lenw - the length (in real *8 words) of the array w provided
c        by the user
c 
c                           Output parameters:
c 
c  ier - error return code:
c    ier=0 means successful execution
c    ier=4 means that the subroutine attempted to create more than
c        maxbox boxes. This is a fatal error.
c    ier=16 means that the subroutine attempted to create more
c        than 199 levels of subdivision. This is a fatal error.
c    ier=1024 means that the amount of space provided by the user
c        in the array w is grossly insufficient (less than n+1000);
c        I would suspect an exceeded array or some other serious
c        nonsense. This is a fatal error.
c    ier=256 means that the amount of space provided by the user
c        in array w is insufficient. This is a fatal error.
c    ier=128 means that the amount of space provided by the user
c        in array w is insufficient. This is a fatal error.
c    ier=32 means that the amount of space provided by the user
c        in array w is insufficient. This is a fatal error.
c  w - array containing the structures created by the subroutine:
c        definitions of boxes, their corners and centers, various
c        types of control information, and inner and outer skeletons
c        for all boxes. The user is expected to access the contents
c        of this array via the subroutines r2dretr, r2d_finskels_retr
c        (see). Please note that keep (see below) first elements of
c        this array must remain unchanged between the call to this
c        subroutine and the subsequent calls to the subroutines
c        r2dretr, r2d_finskels_retr (see).
c  nall - the total number of boxes constructed by the subroutine
c  iz - integer array of length n defining the transposition of
c        the array z (user-specified nodes) such that each of
c        the boxes in the structure contains consecutively numbered
c        points.
c  nlev - the number of levels of subdivision NOT COUNTING THE
C        LEVEL NUMBER 0.
c  conds - the condition numbers (actually, pretty crude estimates)
c        of all of the matrices that had to be inverted in the
c        process of construction the big inverse.
c  nconds - the number of elements returned in the array conds
c  lused - the total number of real *8 elements in array w that
c        have been used by this subroutine
c  keep - the number of real *8 elements of the array w that must
c        remain unchanged between the call to this subroutine and
c        the subsequent calls to the subroutine r2dretr (see).
c 
c 
c        . . . construct and plot the test points
c 
         call r2d_allskels(ier,poteval,rk,eps,z,n,
     1       w,lenw,lused,keep,maxboxes,nbox,nall,iz,nlev)
c 
        call prinf('after r2d_allskels, lused=*',lused,1)
        call prinf('after r2d_allskels, keep=*',keep,1)
  
        istore=keep+10
        lenstore=lenw-istore
c 
        call r2d_evalexpands(ier,poteval,rk,w,eps,
     1      z,iz,nall,w(istore),lenstore,lused4,keep4,maxpnts)
c 
        call prinf('after evalexpands, maxpnts=*',maxpnts,1)
        call prinf('after evalexpands, keep4=*',keep4,1)
c 
c        allocate memory
c 
        mp=maxpnts
        is=istore+keep4
        ls=mp**2*2+10
c 
        is1=is+ls
        ls1=mp**2*2+10
c 
        is2=is1+ls1
        ls2=mp**2*2+10
c 
        ieval=is2+ls2
        leval=mp**2*2+10
        leval=leval*2
c 
        iexpand=ieval+leval
        lexpand=mp**2*2+10
        lexpand=lexpand*2
c 
        isout=iexpand+lexpand
        lsout=mp**2*2+10
c 
        isplit1=isout+lsout
        lsplit1=mp**2*2+10
c 
        isplit2=isplit1+lsplit1
        lsplit2=mp**2*2+10
c 
        iexchg11=isplit2+lsplit2
        lexchg11=mp**2*2+10
c 
        iexchg12=iexchg11+lexchg11
        lexchg12=mp**2*2+10
c 
        iexchg21=iexchg12+lexchg12
        lexchg21=mp**2*2+10
c 
        iexchg22=iexchg21+lexchg21
        lexchg22=mp**2*2+10
c 
        icmerge1=iexchg22+lexchg22
        lcmerge1=mp**2*2+10
c 
        icmerge2=icmerge1+lcmerge1
        lcmerge2=mp**2*2+10
c 
        iinums=icmerge2+lcmerge2
        linums=mp+10
c 
        iirows_outin=iinums+linums
        lirows_outin=mp+10
c 
        iirows_outin1=iirows_outin+lirows_outin
        lirows_outin1=mp+10
c 
        iirows_outin2=iirows_outin1+lirows_outin1
        lirows_outin2=mp+10
c 
        iicols_inout=iirows_outin2+lirows_outin2
        licols_inout=mp+10
c 
        iicols_inout1=iicols_inout+licols_inout
        licols_inout1=mp+10
c 
        iicols_inout2=iicols_inout1+licols_inout1
        licols_inout2=mp+10
c 
        call prinf('before r2d_allmatrs0, istore=*',istore,1)
  
        istore2=iicols_inout2+licols_inout2
        call cleascop(w(istore),w(istore2),keep4/2+100)
c 
        lused0=keep4+100
c 
        call r2d_allmatrs0(ier,poteval,eps,rk,z,n,
     1      w,w(istore2),lenstore,nall,iz,
     2      nlev,lused5,lused0,
     3      w(is),w(is1),w(is2),w(ieval),w(iexpand),
     4      w(isout),w(isplit1),w(isplit2),w(iexchg11),
     5      w(iexchg12),w(iexchg21),w(iexchg22),
     6      w(icmerge1),w(icmerge2),conds,nconds,w(iinums),
     7      w(iirows_outin1),w(iirows_outin2),w(iirows_outin),
     8      w(iicols_inout1),w(iicols_inout2),w(iicols_inout) )
c 
        w(9)=istore2+0.01
c 
c        . . . collect garbage
c 
        call cleascop(w(istore2),w(istore),lused5/2+100)
        w(9)=istore+0.1
c 
        keep=lused5+100
        lused=istore2+lused5
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine r2d_finretr(ibox,kind,w,
     1      ncols_inout,ncols_outin,np,arr,
     2      larr,ifarr)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),arr(1)
c 
c       A call to this subroutine retrieves from the main
c       storage area w one of (many) arrays pertaining to the
c       chunk number ibox; it is fervently hoped that the
c       information to be retrieved had been stores in the
c       array w during a preceding call to the subroutine
c       r2d_allmatrs (see). The chunk number ibox can be a
c       chunk on any level; The type of the array to be
c       retrieved is determined by the parameter kind, in the
c       following manner.
c 
c     kind = 1 - array expand for chunk number ibox
c     kind = 2 - array eval for chunk number ibox
c     kind = 3 - array s for chunk number ibox
c     kind = 4 - integer array irows_outin for chunk number ibox
c     kind = 5 - integer array icols_inout chunk number ibox
c     kind = 6 - array split1 for chunk number ibox
c     kind = 7 - array split2 for chunk number ibox
c     kind = 8 - array exchg11 for chunk number ibox
c 
c     kind = 9 - array exchg12 for chunk number ibox
c     kind = 10 - array exchg21 for chunk number ibox
c     kind = 11 - array exchg22 for chunk number ibox
c     kind = 12 - array cmerge1 for chunk number ibox
c     kind = 13 - array cmerge2 for chunk number ibox
c     kind = 14 - skeleton_inner for chunk number ibox
c 
c                         Input parameters:
c 
c  ibox - the sequence number (address in arrays boxes, index)
c         of the box for which the data are being stored
c  kind - the type of array to be retrieved (see above)
c  w - main storage area for the matrices of all persuasions stored
c        by a preceding call to the subroutine r2d_allmatrs (see)
c 
c 
c                         Output parameters:
c 
c  ncols_inout - the number of nodes in the outgoing inner skeleton
c        of this chunk (also the number of elements in the array
c        icols_inout, etc.)
c  ncols_outin - the number of nodes in the incoming inner skeleton
c        of this chunk (also the number of elements in the array
c        irows_outin, etc.)
c  np - the number of nodes in this chunk
c  arr - the array user has requested (see parameter kind above)
c 
        istore=w(9)
        call r2d_retr(ibox,kind,w(istore),
     1      ncols_inout,ncols_outin,np,arr,
     2      larr,ifarr)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine r2d_allmatrs0(ier,poteval,eps,rk,z,n,
     1      w,store,lenstore,nall,iz,
     2      nlev,lused,lused0,
     3      s,s1,s2,eval,expand,sout,
     4      split1,split2,exchg11,exchg12,exchg21,exchg22,
     5      cmerge1,cmerge2,conds,nconds,inums,
     6      irows_outin1,irows_outin2,irows_outin,
     7      icols_inout1,icols_inout2,icols_inout)
        implicit real *8 (a-h,o-z)
        save
        integer box(12),iz(1),inums(1)
        real *8 z(2,1),w(1),store(1),center(3),
     1           corners(2,5),conds(1)
        complex *16 rk,s(1),s1(1),s2(1),eval(1),expand(1),
     2      sout(1),split1(1),split2(1),
     3    exchg11(1),exchg12(1),exchg22(1),exchg21(1),
     5    cmerge1(1),cmerge2(1)
c 
c 
        integer *4 irows_outin1(1),icols_inout1(1),
     1      irows_outin2(1),icols_inout2(1),
     2      irows_outin(1),icols_inout(1)
c 
        external potuser
c 
c        construct all scattering matrices for all childless
c        boxes, and store them things in array store
c 
        call poteval_init(iz,z,n)
c 
        call r2d_scatmatrs0(potuser,par1,par2,w,nall,
     1      store,icurr,conds,ind,iz,s,inums)
c 
c        progressively merge the scattering matrices, obtaining
c        scattering, splitting, merging, etc. matrices on all
c        levels
c 
        icount=0
        lused=lused0
        do 4000 lev=nlev,1,-1
        do 3600 ibox=1,nall
c 
        call r2dretr(ibox,w,box,corners,center)
c 
        if(box(9) .eq. -7) goto 3600
        if(box(2) .ne. lev) goto 3600
        if(box(4) .eq. 0) goto 3600
c 
c       This box has at least one sonny. If it has ONLY one
c       sonny - copy the sonnie's scattering matrix into daddy's
c       scattering matrix
c 
        if(box(5) .ne. 0) goto 2000
c 
        ison=box(4)
        ifarr=1
        kind=3
        call r2d_retr(ison,kind,store,
     1      ncols_inout,ncols_outin,np,s,lens,ifarr)
c 
        kind=3
        call r2d_store(ier,ibox,kind,
     1      store,lused,lenstore,s,lens)
c 
        goto 3600
c 
 2000 continue
c 
c        This box has 2 sonnies. Merge their scattering matrices
c 
c        . . . retrieve appropriate data from array store
c 
        icount=icount+1
        if(icount .eq. 10) icount=0
c 
        ison1=box(4)
        ifarr=1
        kind=4
        call r2d_retr(ison1,kind,store,
     1      ncols_inout,ncols_outin,np,irows_outin1,
     2      ncols_outin1,ifarr)
c 
        kind=5
        call r2d_retr(ison1,kind,store,
     1      ncols_inout,ncols_outin,np,icols_inout1,
     2      ncols_inout1,ifarr)
c 
        ison2=box(5)
        ifarr=1
        kind=4
        call r2d_retr(ison2,kind,store,
     1      ncols_inout7,ncols_outin7,np,irows_outin2,
     2      ncols_outin2,ifarr)
c 
        kind=5
        call r2d_retr(ison2,kind,store,
     1      ncols_inout7,ncols_outin7,np,icols_inout2,
     2      ncols_inout2,ifarr)
c 
        ifarr=1
        kind=4
        call r2d_retr(ibox,kind,store,
     1      ncols_inout7,ncols_outin7,np,irows_outin,
     2      ncols_outin,ifarr)
c 
        kind=5
        call r2d_retr(ibox,kind,store,
     1      ncols_inout7,ncols_outin7,np,icols_inout,
     2      ncols_inout,ifarr)
c 
        kind=3
        call r2d_retr(ison1,kind,store,
     1      ncols_inout7,ncols_outin7,np,s1,
     2      lens1,ifarr)
c 
        call r2d_retr(ison2,kind,store,
     1      ncols_inout7,ncols_outin7,np,s2,
     2      lens2,ifarr)
c 
        kind=1
        call r2d_retr(ibox,kind,store,
     1      ncols_inout7,ncols_outin7,np,expand,
     2      lenexpand,ifarr)
c 
        kind=2
        call r2d_retr(ibox,kind,store,
     1      ncols_inout7,ncols_outin7,np,eval,
     2      leneval,ifarr)
c 
        iflast=0
        if(ibox .eq. 1) iflast=1
c 
        lenw2=lenstore-lused*2
  
cccc        call prinf('before scattmatr_merge2, lused/1000=*',
cccc     1      lused/1000,1)
cccc        call prinf('before scattmatr_merge2, lenstore/1000=*',
cccc     1      lenstore/1000,1)
c 
        if(icount .eq. 0) then
c 
            call prinf('before scattmatr_merge2, lev=*',lev,1)
            call prinf('while, ibox=*',ibox,1)
            call prinf('and lenw2/1000=*',lenw2/1000,1)
        endif
c 
        iiiwww2=nall*40+lused*2+100
c 
        call scattmatr_merge2(ier,expand,eval,
     1    iz,z,
     2    irows_outin1,ncols_outin1,
     3    irows_outin2,ncols_outin2,
     4    icols_inout1,ncols_inout1,
     5    icols_inout2,ncols_inout2,
     6    irows_outin,ncols_outin,
     7    icols_inout,ncols_inout,s1,s2,sout,split1,split2,
     8    exchg11,exchg12,exchg22,exchg21,cmerge1,
     9    cmerge2,conds(ind+1),store(iiiwww2),lenw2,
     a    ltot,iflast,potuser,par1,par2)
c 
        ind=ind+2
c 
        if(ier .ne. 0) then
            call prinf('after scattmatr_merge2, ier=*',ier,1)
            stop
        endif
c 
c        store the obtained data in array store
c 
        lensout=ncols_inout*ncols_outin
        kind=3
        call r2d_store(ier,ibox,kind,
     1      store,lused,lenstore,sout,lensout)
c 
        kind=6
        lensplit1=ncols_outin1*ncols_outin
        call r2d_store(ier,ibox,kind,
     1      store,lused,lenstore,split1,lensplit1)
c 
        kind=7
        lensplit2=ncols_outin2*ncols_outin
        call r2d_store(ier,ibox,kind,
     1      store,lused,lenstore,split2,lensplit2)
c 
        kind=8
        lenexchg11=ncols_outin1*ncols_inout1
        call r2d_store(ier,ibox,kind,
     1      store,lused,lenstore,exchg11,lenexchg11)
c 
        kind=9
        lenexchg12=ncols_outin1*ncols_inout2
        call r2d_store(ier,ibox,kind,
     1      store,lused,lenstore,exchg12,lenexchg12)
c 
        kind=10
        lenexchg21=ncols_outin2*ncols_inout1
        call r2d_store(ier,ibox,kind,
     1      store,lused,lenstore,exchg21,lenexchg21)
c 
        kind=11
        lenexchg22=ncols_outin2*ncols_inout2
        call r2d_store(ier,ibox,kind,
     1      store,lused,lenstore,exchg22,lenexchg22)
c 
        kind=12
        lencmerge1=ncols_inout1*ncols_inout
        call r2d_store(ier,ibox,kind,
     1      store,lused,lenstore,cmerge1,lencmerge1)
c 
        kind=13
        lencmerge2=ncols_inout2*ncols_inout
        call r2d_store(ier,ibox,kind,
     1      store,lused,lenstore,cmerge2,lencmerge2)
c 
 3600 continue
 4000 continue
c 
        nconds=ind
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine r2d_scatmatrs0(potuser,par1,par2,w,
     1      nboxes,store,icurr,conds,ind,iz,s,inums)
        implicit real *8 (a-h,o-z)
        save
        complex *16 par1(1),par2(1)
        dimension w(1),corners(8),center(2),
     1      store(1),conds(1),s(100000),inums(10000),
     2      work(100000)
c 
        integer box(12),iz(1)
c 
c        construct scattering matrices for all childless boxes
c 
        ind=0
        icount=0
        do 3000 ibox=2,nboxes
c 
        icount=icount+1
        if(icount .eq. 10) icount=0
        if(icount .eq. 0) call prinf('in r2d_scatmatrs0, ibox=*',
     1      ibox,1)
c 
        call r2dretr(ibox,w,box,corners,center)
c 
        if(box(4) .gt. 0) goto 3000
        if(box(9) .eq. -7) goto 3000
c 
c        this box is childless; create the scattering
c        matrix for it
c 
        istart=box(6)
        np=box(7)
c 
        call r2d_scatmatr0(istart,np,potuser,par1,par2,
     1      s,cond,work)
c 
        ind=ind+1
c 
        conds(ind)=cond-1
c 
c        store the obtained scattering matrix
c 
        do 2600 i=1,np
c 
        inums(i)=istart+i-1
 2600 continue
c 
        kind=3
        larr=np**2
c 
        call r2d_store(ier,ibox,kind,
     1      store,lused,lenstore,s,larr)
c 
c 
        kind=4
        larr=np
c 
        call r2d_store(ier,ibox,kind,
     1      store,lused,lenstore,inums,larr)
c 
        kind=5
        larr=np
c 
        call r2d_store(ier,ibox,kind,
     1      store,lused,lenstore,inums,larr)
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
        subroutine r2d_scatmatr0(istart,np,potuser,par1,par2,
     1      s,cond,work)
        implicit real *8 (a-h,o-z)
        save
        complex *16 val,s(np,np),work(1),par1(1),par2(1)
c 
c        construct the scattering matrix in the case when
c        the box contains 1 node
c 
        if(np .ne. 1) goto 2000
c 
        call poteval(istart,istart,potuser,par1,par2,val)
c 
        s(1,1)=-1/val
c 
        cond=1
c 
        return
c 
 2000 continue
c 
c        construct the scattering matrix in the case when
c        the box contains more than 1 node
c 
         do 2400 i=1,np
         do 2200 j=1,np
c 
        iii=istart+i-1
        jjj=istart+j-1
  
cccc        call prinf('before poteval, iii=*',iii,1)
cccc        call prinf('before poteval, jjj=*',jjj,1)
  
        call poteval(iii,jjj,potuser,par1,par2,val)
        s(i,j)=val
c 
 2200 continue
 2400 continue
c 
c        invert the matrix of interactions
c 
        call corthom(s,np,work,cond)
c 
        call csign_chg(s,s,np**2)
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
  
c 
c 
c 
c 
c 
c        Following is a description of a single element of
c        the array index(2,20,nboxes)
c 
c        For any i,j, the element index(1,i,j) is the location
c        in the storade area of the first element of the piece
c        of data specified by the integers (i,j). The element
c        index(2,i,j) is the length (in complex *16 words) of
c        the said piece of data.
c 
c  index(*,1,*) - array expand for this box;  please note
c        that it does not exist for childless boxes
c 
c  index(*,2,*) - array eval for this box; please note
c        that it does not exist for childless boxes
c 
c  index(*,3,*) - array s for this box
c 
c  index(*,4,*) - integer array irows_outin for this box
c 
c  index(*,5,*) - integer array icols_inout for this box
c 
c  index(*,6,*) - array split1 for this box; please note
c        that it does not exist for childless boxes
c 
c  index(*,7,*) - array split2 for this box; please note
c        that it does not exist for childless boxes
c 
c  index(*,8,*) - array exchg11 for this box; please note
c        that it does not exist for childless boxes
c 
c  index(*,9,*) - array exchg12 for this box
c 
c  index(*,10,*) - array exchg21 for this box
c 
c  index(*,11,*) - array exchg22 for this box
c 
c  index(*,12,*) - array cmerge1 for this box; please note
c        that it does not exist for childless boxes
c 
c  index(*,13,*) - array cmerge2 for this box; please note
c        that it does not exist for childless boxes
c 
c  index(1,20,ibox) - the number of points in the box ibox
c 
