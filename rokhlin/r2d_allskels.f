        implicit real *8 (a-h,o-z)
        integer *4 iz(1000 000),box(12)
c 
        real *8 z(2,1000 000),w(1000 000),corners(2,5),
     2      center(3),zsout(2,10000),zsin(2,10000),
     3      zssons(2,10000),zsdad(2,10000),
     4      errors_expand(10000),errors_eval(10000),
     5      store(2000 000)
c 
        complex *16 rk,eval(1000 000),expand(1000 000)
c 
        external poteval2
c 
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
         rk=20
         rk=100
cccc         rk=600
  
  
c 
c       construct and plot the test points
c 
  
        call RANLST(NSMALL,Z,W,N)
          call prinf('after ranlst, n=*',n,1)
        a=10
        b=1
c 
cccc        call creelips(z,n,a,b)
cccc        n=n-1
  
        call prinf('after ranlst, n=*',n,1)
  
  
  
  
  
c 
        iw=21
        call zquaplot(iw,z,n,2,'random points as created*')
  
  
        lenw=1000 000
        maxboxes=100000
        nbox=30
  
        gap=0.2
  
  
         call r2d_allskels(ier,poteval2,rk,eps,gap,z,n,
     1       w,lenw,lused,keep,maxboxes,nbox,nall,iz,nlev)
  
c 
c        plot the whole structure
c 
  
        iw=25
        call allplot3(iw,w,z)
  
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
        call prinf('ibox=*',ibox,1)
        call prinf('and lenin=*',lenin,1)
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
cccc        stop
  
        eps=1.0d-12
  
        lenstore=1000 000
        call r2d_evalexpands(ier,poteval2,rk,w,eps,
     1      z,iz,nall,store,lenstore,lused4,keep4,maxpnts)
  
  
        call irows_icols_test(z,iz,w,store,nall)
  
cccc        stop
  
c 
c        test the constructed skeletons and matrices eval,
c        expand
c 
        do 5400 ibox=2,nall
c 
        call prinf('ibox=*',ibox,1)
c 
        call r2dretr(ibox,w,box,corners,center)
c 
        if(box(4) .eq. 0) goto 5400
c 
  
        ifarr=1
        kind=14
        call r2d_retr(ibox,kind,store,
     1      ncols_inout,ncols_outin,np,zsdad,lendad,ifarr)
c 
        ison1=box(4)
        call r2d_retr(ison1,kind,store,
     1      ncols_inout,ncols_outin,np,zssons,lensons,ifarr)
c 
        if(box(5) .eq. 0) goto 5200
c 
        ison2=box(5)
        call r2d_retr(ison2,kind,store,
     1      ncols_inout,ncols_outin,np,zssons(1,lensons+1),
     2      lenson2,ifarr)
c 
        lensons=lensons+lenson2
c 
 5200 continue
c 
        kind=1
        call r2d_retr(ibox,kind,store,
     1      ncols_inout,ncols_outin,np,expand,lexpand,ifarr)
  
  
        kind=2
        call r2d_retr(ibox,kind,store,
     1      ncols_inout,ncols_outin,np,eval,leval,ifarr)
  
  
c 
c 
c        test the obtained matrices eval, expand
c 
        call evalexpand_test(poteval2,rk,ibox,w,
     1      eval,expand,zsdad,zssons,lendad,lensons,eps,
     2      z,n,errmax_expand,errmax_eval)
c 
        call prin2('errmax_expand=*',errmax_expand,1)
        call prin2('errmax_eval=*',errmax_eval,1)
c 
        errors_expand(ibox)=errmax_expand
        errors_eval(ibox)=errmax_eval
c 
 5400 continue
c 
        call prin2('and errors_expand=*',errors_expand,nall)
        call prin2('and errors_eval=*',errors_eval,nall)
        stop
        end
  
c 
c 
c 
c 
c 
        subroutine irows_icols_test(z,iz,w,store,nall)
        implicit real *8 (a-h,o-z)
        save
        integer *4 box(12),iz(1),nodes(10000)
        dimension w(1),zs(2,10000),center(3),store(1),z(2,1),
     1      corners(2,5),irows(10000),icols(10000),zs2(2,1000)
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
        call r2d_retr(ibox,kind,store,
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
        call r2d_retr(ibox,kind,store,
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
cccc        call prin2('and zs2=*',zs2,len2*2)
  
        do 1600 i=1,len2
c 
        zs2(1,i)=zs2(1,i)-zs(1,i)
        zs2(2,i)=zs2(2,i)-zs(2,i)
c 
        error=error+zs2(1,i)**2+zs2(2,i)**2
  
 1600 continue
c 
cccc        call prin2('and the difference is*',zs2,len2*2)
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end off the debugging code and the beginning
c        of the code for the construction of skeletons and matrices
c        eval, expand
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        There are 4 user-callable subroutines in this file:
c        r2d_allskels, r2d_evalexpands, r2d_store, r2d_retr.
c        Following is a brief description of thse subroutines.
c 
c   r2d_allskels - first, this subroutine constructs the "bastardized"
c        binary tree, obtained from the standard quad-tree by
c        the introduction of "bastard" boxes, so that every
c        box is subdivided into two sons (one of which might
c        not exist).
c 
c        Second, this subroutine constructs all skeletons on all
c        levels. The principal output of this subroutine is the
c        array w, containing the structure created by the subroutine:
c        definitions of boxes, their corners and centers, various
c        types of control information, and inner and outer skeletons
c 
c   r2d_evalexpands - constructs the matrices eval, expand
c        for all boxes on all levels. It also uses the skeletons
c        of all boxes constructed via a preceding call to the
c        subroutine r2d_allskels (see) to construct the arrays
c        irows, icols.
c 
c   r2d_store - a call to this subroutine stores in the storage
c       area store one of (many) possible arrays pertaining to the
c       chunk number ibox.
c 
c   r2d_retr - a call to this subroutine retrieves from the storage
c       area store one of (many) possible arrays pertaining to
c       the chunk number ibox.
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        subroutine r2d_evalexpands(ier,poteval,rk,w,eps,
     1      z,iz,nall,store,lenstore,lused,keep,maxpnts)
        implicit real *8 (a-h,o-z)
        save
        integer *4 box(12)
        dimension w(1),center(3),corners(2,5),store(1)
        complex *16 rk
c 
c        This subroutine constructs the matrices eval, expand
c        for all boxes on all levels. It also uses the skeletons
c        of all boxes constructed via a preceding call to the
c        subroutine r2d_allskels (see) to construct the arrays irows,
c        icols. All of the data generated by this subroutine are
c        returned in the array store (see below) where it can be
c        accessed (and added to) via calls to the subroutines
c        r2d_store, r2d_retr (see). PLEASE NOTE THAT THIS SUBROUTINE
C        USES THE ARRAY W GENERATED BY A PRECEDING CALL TO THE
C        SUBROUTINE R2D_ALLSKELS. THIS SUBROUTINE HAS NO KNOWN USES
C        AS A STAND-ALONE DEVICE.
c 
c 
c                           Input parameters:
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
c  rk - the Helmholtz coefficient
c  w - array containing the structures created by a prior call to
c        the subroutine r2d_allskels (see).
c  eps - the accuracy to which the SVDs will be performed
c  z - the points in R^2 on top of which the tree structure
c        is placed.
c  iz - integer array of length n defining the transposition of
c        the array z (user-specified nodes) such that each of
c        the boxes in the structure contains consecutively numbered
c        points.
c  nall - the total number of boxes constructed by the subroutine
c  lenstore - the amount of space (in real *8 locations) provided
c        by the user in the storage/work array store
  
c                           Output parameters:
c 
c  ier - error return code:
c  store - array containing inner skeletons of all boxes, the
c        arrays irows, icols for all boxes, and matrices eval,
c        expand for all boxes. The data in the array store can
c        (and should) be retrieved and added to vis the subroutines
c        r2d_retr, r2d_store (see)
c  lused - the total number of real *8 elements in array store
c        used by the subroutine
c  keep - the number of real *8 elements in array store that have
c        to be unchanged between the call to this subroutine and
c        subsequent calls to the subroutines r2d_retr, r2d_store
c 
c 
c       . . . scan all boxes, determining the maximum
c             dimensionalities of the arrays eval, expand,
c             zsdadout, zsdad, zssons, c, a, w2, rnorms
c 
        maxpnts=0
        do 1400 i=2,nall
c 
        call r2dretr(i,w,box,corners,center)
c 
        if(box(9) .eq. -7) goto 1400
c 
        ifinner=1
        call r2d_finskels_retr(w,i,ifinner,store,len)
c 
        if(len .gt. maxpnts) maxpnts=len
 1400 continue
c 
        call prinf('maxpnts =*',maxpnts,1)
c 
c       allocate memory
c 
        ia=1
        la=maxpnts**2*2+10
c 
        ic=ia+la
        lc=(maxpnts**2*2+10)*2
c 
        iw2=ic+lc
        lw2=maxpnts**2*2 *8 +1000
c 
        izsdadout=iw2+lw2
        lzsdadout=maxpnts*2+10
c 
        izsdad=izsdadout+lzsdadout
        lzsdad=maxpnts*2+10
c 
        izssons=izsdad+lzsdad
        lzssons=maxpnts*4+10
c 
        ieval=izssons+lzssons
        leval=maxpnts**2*4+10
c 
        iexpand=ieval+leval
        lexpand=maxpnts**2*4+10
c 
        irnorms=iexpand+lexpand
        lrnorms=maxpnts*4+10
c 
        istore=irnorms+lrnorms
c 
        lenstore2=lenstore-istore-2
c 
c        construct all expand and eval matrices, and store them
c        things in array store
c 
        call r2d_evalexpands0(ier,poteval,rk,w,eps,nall,
     1      store(istore),lenstore2,lused,
     2      store(ia),store(ic),store(iw2),store(izsdadout),
     3      store(izsdad),store(izssons),store(ieval),
     4      store(iexpand),store(irnorms) )
c 
c        for each box that is not childless, find the coordinates
c        in array iz of the nodes in its (inner) skeleton
c 
        call r2d_irows_icols(z,iz,w,store(istore),lenstore2,
     1      nall,lused,store(izsdad),store(izssons),
     2      store(irnorms) )
c 
        call cleascop(store(istore),store,lused/2+10)
c 
        keep=lused+10
        lused=lused+istore+10
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine r2d_irows_icols(z,iz,w,store,lenstore,
     1      nall,lused,zs,zs2,nodes)
        implicit real *8 (a-h,o-z)
        save
        integer *4 box(12),iz(1),nodes(1)
        dimension w(1),zs(2,1),center(3),store(1),
     1      z(2,1),corners(2,5),zs2(2,1),zz(3)
c 
c        This subroutine retrieves skeletons of all boxes on
c        all levels, constructs the maps of these skeletons
c        in the array iz, and stores the obtained maps in the
c        array store, from which array they can be retrieved
c        with the help of the subroutine r2d_rest, as the
c        integer arrays irows, icols.
c 
c        PLEASE NOTE THAT IN THIS VERSION, THE ARRAY IROWS
C        IS THE SAME as the array icols for all boxes.
c 
c 
c        . . . one box after another, find locations in array
c              iz of the elements of skeletons of said boxes
c 
        eps=1.0d-60
        icount=0
c 
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
        call r2d_retr(ibox,kind,store,
     1      ncols_inout,ncols_outin,np,zs,len,ifarr)
c 
        j0=box(6)-1
        iii=0
        npts=box(7)
c 
        do 1400 i=1,len
c 
        do 1200 j=1,npts
c 
        jj=j0+j
        zz(1)=z(1,iz(jj))
        zz(2)=z(2,iz(jj))
c 
        d=(zs(1,i)-zz(1))**2+(zs(2,i)-zz(2))**2
c 
        if(d .gt. eps) goto 1200
c 
        nodes(i)=jj
c 
        zs2(1,i)=z(1,iz(jj))
        zs2(2,i)=z(2,iz(jj))
c 
 1200 continue
c 
 1400 continue
c 
c       store the obtained nodes
c 
        kind=4
        call r2d_store(ier,ibox,kind,
     1      store,lused,lenstore,nodes,len)
c 
        kind=5
        call r2d_store(ier,ibox,kind,
     1      store,lused,lenstore,nodes,len)
c 
 2000 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine r2d_evalexpands0(ier,poteval,rk,w,eps,nall,
     1      store,lenstore,lused,
     2      a,c,w2,zsdadout,zsdad,zssons,eval,expand,rnorms)
        implicit real *8 (a-h,o-z)
        save
        integer *4 box(12)
        dimension w(1),zsdad(2,1),zssons(2,1),center(3),
     1      corners(2,5),store(1)
        complex *16 rk,eval(1),expand(1)
  
        complex *16 a(1),zsdadout(1),
     1      c(1),w2(1),rnorms(1)
c 
c        initialize the storage-retrieval facility
c 
        iboxes=2
c 
        call r2d_store_init(store,nall,w)
c 
c       one box after another, construct the matrices eval, expand,
c       and store them in array store
c 
        icount=0
        do 4400 ibox=2,nall
c 
        icount=icount+1
        if(icount .eq. 10) icount=0
        if(icount .eq. 0)
     1      call prinf('in r2d_evalexpands0, ibox=*',ibox,1)
c 
        call r2dretr(ibox,w,box,corners,center)
c 
        if(box(9) .eq. -7) goto 4400
        if(box(4) .eq. 0) goto 4200
c 
        call evalexpand_get(poteval,rk,ibox,w,
     1      eval,expand,zsdad,zssons,lendad,lensons,eps,
     2      zsdadout,a,c,w2,rnorms)
c 
        kind=14
        call r2d_store(ier,ibox,kind,
     1      store,lused,lenstore,zsdad,lendad)
c 
        kind=1
        call r2d_store(ier,ibox,kind,
     1      store,lused,lenstore,expand,lendad*lensons)
c 
        kind=2
        call r2d_store(ier,ibox,kind,
     1      store,lused,lenstore,eval,lendad*lensons)
c 
        goto 4400
c 
 4200 continue
c 
c        this box is childless. store its inner skeleton
c        (consisting of all nodes inside it) in the storage
c        area store
c 
        ifinner=1
        call r2d_finskels_retr(w,ibox,ifinner,zsdad,lendad)
c 
        kind=14
        call r2d_store(ier,ibox,kind,
     1      store,lused,lenstore,zsdad,lendad)
c 
 4400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine r2d_store(ier,ibox,kind,
     1      store,lused,lenstore,arr,larr)
        implicit real *8 (a-h,o-z)
        save
        complex *16 arr(1)
        dimension store(1),w(1)
c 
c 
c       A call to this subroutine stores in the storage area
c       store one of (many) possible arrays pertaining to the
c       chunk number ibox. The chunk number ibox can be a
c       chunk on any level. The type of array to be stored
c       is determined by the parameter kind, in the
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
c  store - storage area, initialized via a prior call to the
c         subroutine (entry, actually) r2d_store_init (see)
c  lenstore - the amount of space (in real *8 words) provided
c         in array store. Appears to be unused at this time
c         (10.25.02)
c  arr - the array to be stored
c  larr - the length of the array to be stored
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
c 
        lenww=(lenstore-istore)/2
        call r2d_store0(ier,ibox,kind,
     1      store(iindex),store(istore),lused,lenww,arr,larr)
c 
        lused=lused*2+istore
  
c 
        return
c 
c 
c 
c 
        entry r2d_store_init(store,nboxes,w)
c 
        iindex=21
        iboxes=w(2)
c 
        istore=iindex+40*nboxes+10
c 
        call r2d_store_init0(store(iindex),nboxes,w(iboxes) )
c 
        return
c 
c 
c 
c 
        entry r2d_retr(ibox,kind,store,
     1      ncols_inout,ncols_outin,np,arr,larr,ifarr)
c 
c       A call to this subroutine retrieves from the storage
c       area store one of (many) possible arrays pertaining to
c       the chunk number ibox. The chunk number ibox can be a
c       chunk on any level. The type of array to be stored
c       is determined by the parameter kind, in the manner
c       described above (see entry r2d_store)
c 
  
        call r2d_retr0(ibox,kind,store(iindex),store(istore),
     1      ncols_inout,ncols_outin,np,arr,larr,ifarr)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine r2d_store0(ier,ibox,kind,
     1      index,w,lused,lenww,arr,larr)
        implicit real *8 (a-h,o-z)
        save
        complex *16 w(1),arr(1)
        integer *4 boxes(10,1),index(2,20,1)
c 
c 
c        store the data about this box in arrays index, w
c 
        if(( kind .eq. 4) .or. (kind .eq. 5) ) goto 2000
c 
c        . . . in the case when the array to be stored is complex
c 
        index(1,kind,ibox)=icurr
        index(2,kind,ibox)=larr
c 
        call cleascop(arr,w(icurr),larr)
c 
        icurr=icurr+larr+5
        lused=icurr
c 
        return
c 
 2000 continue
c 
c        . . . in the case when the array to be stored is integer
c 
        index(1,kind,ibox)=icurr
        index(2,kind,ibox)=larr
c 
        call r2d_intcopy(arr,w(icurr),larr)
c 
        icurr=icurr+larr/2+5
        lused=icurr
c 
        return
c 
c 
c 
c 
        entry r2d_store_init0(index,nboxes,boxes)
c 
        call r2d_intzero(index,40*nboxes)
c 
        do 4200 i=1,nboxes
c 
        index(1,20,i)=boxes(7,i)
 4200 continue
c 
        icurr=1
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine r2d_retr0(ibox,kind,index,w,
     1      ncols_inout,ncols_outin,np,arr,larr,ifarr)
        implicit real *8 (a-h,o-z)
        save
        complex *16 w(1),arr(1)
c 
        integer *4 index(2,20,1)
c 
c       A call to this subroutine retrieves from the main
c       storage area w (with the help of the auxiliary array
c       index)  one of (many) arrays pertaining to the
c       chunk number ibox. The chunk number ibox can be a
c       chunk on any level. The type of matrix to be
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
c  index - integer array dimensioned (2,20,1), containing the
c        addresses in the main storage area of various matrices
c        corresponding to physical interactions between boxes,
c        nodes, etc. For any i,j, the element index(1,i,j) is the
c        location in the storage area of the first element of the
c        piece of data specified by the integers (i,j). The element
c        index(2,i,j) is the length (in complex *16 words) of the
c        said piece of data.
c  w - main storage area for the matrices of all persuasions stored
c        by the subroutine rikers_store_lev0 (see)
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
        np=index(1,20,ibox)
        ncols_outin=index(2,4,ibox)
        ncols_inout=index(2,5,ibox)
c 
        if(kind .eq. 0) return
c 
c       retrieve the user-requested array if it is a complex one
c 
        if(( kind .eq. 4) .or. (kind .eq. 5) ) goto 2000
c 
        i=index(1,kind,ibox)
        larr=index(2,kind,ibox)
c 
        if(ifarr .eq. 1) call cleascop(w(i),arr,larr)
c 
        return
 2000 continue
c 
c       retrieve the user-requested array if it is an integer one
c 
        i=index(1,kind,ibox)
        larr=index(2,kind,ibox)
c 
        if(ifarr .eq. 1) call r2d_intcopy(w(i),arr,larr)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine evalexpand_get(poteval,rk,ibox,w,
     1      eval,expand,zsdad,zssons,lendad,lensons,eps,
     2      zsdadout,a,c,w2,rnorms)
        implicit real *8 (a-h,o-z)
        save
        integer *4 box(12)
        dimension w(1),zsdad(2,1),zssons(2,1),center(3),
     1      corners(2,5),zsdadout(2,1)
        complex *16 rk,eval(1),expand(1),a(1),
     1      c(1),w2(1),rnorms(1)
c 
c        if this is a childless box - exit
c 
        call r2dretr(ibox,w,box,corners,center)
c 
        if(box(9) .eq. -7) then
            ier=8
            return
        endif
c 
        if(box(4) .eq. 0) then
            ier=4
            return
        endif
c 
c       if the box has only one child - both expand and eval
c       are identities. act accordingly
c 
        if(box(5) .ne. 0) goto 2000
c 
        ifinner=1
        call r2d_finskels_retr(w,ibox,ifinner,zsdad,lendad)
c 
        ison=box(4)
        call r2d_finskels_retr(w,ison,ifinner,zssons,lensons)
c 
        neval=lensons
        call r2d_ident_matr(eval,neval)
        call r2d_ident_matr(expand,neval)
c 
        return
 2000 continue
c 
c        box number ibox has two kids. retrieve from array w
c        the inner skeletons for the daddy and the sonnies,
c        and the outer skeleton for the daddy
c 
        ifinner=1
        call r2d_finskels_retr(w,ibox,ifinner,zsdad,lendad)
c 
        ison1=box(4)
        call r2d_finskels_retr(w,ison1,ifinner,zssons,lenson1)
c 
        ison2=box(5)
        call r2d_finskels_retr(w,ison2,ifinner,
     1      zssons(1,lenson1+1),lenson2)
c 
        lensons=lenson1+lenson2
c 
        ifinner=0
        call r2d_finskels_retr(w,ibox,ifinner,zsdadout,lendad)
c 
c        construct the matrix eval
c 
        call evalmatrget(poteval,rk,a,c,eval,lendad,lensons,
     1      zsdad,zsdadout,zssons,eps,w2,rnorms)
c 
c        construct the matrix expand
c 
        call expandmatrget(poteval,rk,a,c,expand,
     1      lendad,lensons,zsdad,zsdadout,zssons,eps,w2,rnorms)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine expandmatrget(poteval,rk,a,c,expand,
     1      nskel,nsons,zsin,zsout,zsons,eps,w,rnorms)
        implicit real *8 (a-h,o-z)
        save
        complex *16 expand(nskel,nsons),a(nskel,nskel),
     1      c(nskel,nsons),w(1),rnorms(1),rk,pot
c 
        dimension zsin(2,1),zsout(2,1),zsons(2,1)
c 
c        construct the matrix of potentials of all points in the
c        outer skeleton on the points in the inner skeleton
c 
        do 1400 i=1,nskel
        do 1200 j=1,nskel
c 
        call poteval(zsin(1,j),zsout(1,i),rk,pot)
c 
        a(i,j)=pot
 1200 continue
 1400 continue
c 
c        construct the matrix of potentials of all points
c        in the skeletons of the sonnies at all points in the
c        outer skeleton
c 
        do 1800 i=1,nskel
        do 1600 j=1,nsons
c 
        call poteval(zsons(1,j),zsout(1,i),rk,pot)
c 
        c(i,j)=pot
 1600 continue
 1800 continue
c 
c        construct the expansion matrix
c 
        ifcheck=1
        call ncleamatll(a,nskel,nskel,nsons,c,expand,
     1      eps,ncols,
     2    rnorms,w,ifcheck,errl2,errmax)
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine evalmatrget(poteval,rk,a,c,eval,
     1      nskel,nsons,zsin,zsout,zsons,eps,w,rnorms)
        implicit real *8 (a-h,o-z)
        save
        complex *16 eval(nsons,nskel),a(nskel,nskel),
     1      c(nsons,nskel),w(1),rnorms(1),rk,pot
c 
        dimension zsin(2,1),zsout(2,1),zsons(2,1)
c 
c        construct the matrix of potentials of all points in the
c        outer skeleton on the points in the inner skeleton
c 
        do 1400 i=1,nskel
        do 1200 j=1,nskel
c 
        call poteval(zsin(1,j),zsout(1,i),rk,pot)
c 
        a(j,i)=pot
 1200 continue
 1400 continue
c 
c        construct the matrix of potentials of all points
c        in the outer skeleton at all points in the skeletons
c        of the sonnies
c 
        do 1800 i=1,nskel
        do 1600 j=1,nsons
c 
        call poteval(zsons(1,j),zsout(1,i),rk,pot)
c 
        c(j,i)=pot
 1600 continue
 1800 continue
c 
c        construct the evaluation matrix
c 
        ifcheck=1
        call ncleamatrr(a,nsons,nskel,nskel,c,eval,
     1      eps,ncols,rnorms,w,ifcheck,errl2,errmax)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine r2d_allskels(ier,poteval,rk,eps,gap,
     1      z,n,w,lenw,lused,keep,maxboxes,nbox,nall,iz,nlev)
c 
        implicit real *8 (a-h,o-z)
        save
        integer *4 iz(1),map(2,1000),ns(100)
        real *8 z(2,1),w(1)
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
c        levels of subdivision.
c 
c        Second, this subroutine constructs all skeletons on all
c        levels. The principal output of this subroutine is the
c        array w, containing the structure created by the subroutine:
c        definitions of boxes, their corners and centers, various
c        types of control information, and inner and outer skeletons
c        for all boxes.  The user is expected to access the contents
c        of array w via the subroutines r2dretr, r2d_finskels_retr
c        (see). Please note that keep (see below) first elements of
c        this array must remain unchanged between the call to this
c        subroutine and the subsequent calls to the subroutines
c        r2dretr, r2d_finskels_retr.
C 
C        PLEASE NOTE THAT WHILE THE SUBROUTINE R2D_FINSKELS_RETR
C        IS (AS OF 10.25.02) A PART OF THIS FILE, THE SUBROUTINE
C        R2DRETR IS A PART OF THE FILE r2dstrcr.f. ATTEMPTS TO
C        FIND A SUBROUTINE IN A FILE WHERE IT IS NOT ARE UNLIKELY
C        TO BE SUCCESSFUL.
c 
c                           Input parameters:
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
c  rk - the Helmholtz coefficient
c  eps - the accuracy to which the SVDs will be performed
c  z - the points in R^2 on top of which the tree structure
c        is placed.
c  n - the number of points in array z
c  lenw - the length (in real *8 words) of the array w provided
c        by the user
c  maxboxes - the maximum number of boxes the subroutine is
c        permitted to construct. If this number is insufficient,
c        the error return code ier (see below) is set to 4, and
c        the execution is terminated
c  nbox - the maximum number of points permitted in a childless
c        box
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
c  nall - the total number of boxes constructed by the subroutine
c  iz - integer array of length n defining the transposition of
c        the array z (user-specified nodes) such that each of
c        the boxes in the structure contains consecutively numbered
c        points.
c  w - array containing the structures created by the subroutine:
c        definitions of boxes, their corners and centers, various
c        types of control information, and inner and outer skeletons
c        for all boxes. The user is expected to access the contents
c        of this array via the subroutines r2dretr, r2d_finskels_retr
c        (see). Please note that keep (see below) first elements of
c        this array must remain unchanged between the call to this
c        subroutine and the subsequent calls to the subroutines
c        r2dretr, r2d_finskels_retr (see).
c  lused - the total number of real *8 elements in array w that
c        have been used by this subroutine
c  keep - the number of real *8 elements of the array w that must
c        remain unchanged between the call to this subroutine and
c        the subsequent calls to the subroutine r2dretr (see).
c  iz - integer array of length n defining the transposition of
c        the array z (user-specified nodes) such that each of
c        the boxes in the structure contains consecutively numbered
c        points.
c  nlev - the number of levels of subdivision NOT COUNTING THE
C        LEVEL NUMBER 0.
c 
c 
c        . . . construct the tree structure (augmenetd with
c              the neighbors)
c 
        call r2dstrcr_neighb(jer,z,n,
     1       maxboxes,nbox,nall,iz,w,lenw,keep,lused,
     2       map,nlev)
c 
        iboxes=w(2)
        call prinf('boxes as created*',w(iboxes),10*nall)
  
  
        call prinf('after r2dstrcr, jer=*',jer,1)
  
        if(ier .ne. 0) stop
        call prinf('after r2dstrcr, nlev=*',nlev,1)
  
        call prinf('and keep=*',keep,1)
        call prinf('and lused=*',lused,1)
c 
c        construct uncompressed outer skeletons for all
c        standard boxes on all levels
c 
        lstore=10 000 000
  
        istore=keep+1
        call r2d_outer_skels(ier,w,nlev,
     1      poteval,rk,gap,eps,w(istore),lstore,lused,keep2,
     2      ns,nptsmax)
  
        call prinf('after outer_skels, lused=*',lused,1)
        call prinf('after outer_skels, ns=*',ns,nlev)
        call prinf('after outer_skels, nptsmax=*',nptsmax,1)
  
cccc        stop
  
c 
c        construct all of the final inner and outer skeletons,
c        and store them things in array w
c 
        izsin=keep2*2+istore+1
        lzsin=n*2+10 +nptsmax*2
c 
        izsout=izsin+lzsin
        lzsout=n*2+10 +nptsmax*2
c 
        iskelin=izsout+lzsout
        lskelin=n*2+10 +nptsmax*2
c 
        iskelout=iskelin+lskelin
        lskelout=n*2+10 +nptsmax*2
c 
        iinband=iskelout+lskelout
        linband=n+10 +nptsmax*2
c 
        ilarge_sepa=iinband+linband
        llarge_sepa=nall+10
c 
        istore2=ilarge_sepa+llarge_sepa
c 
        w(8)=istore2+0.1
c 
        lenstore2=lenw-istore2-1
c 
        call r2d_outer_skel3(ier,w,z,iz,nlev,nall,
     1      poteval,rk,eps,w(ilarge_sepa),w(iinband),
     2      w(istore),w(izsout),w(izsin),
     3      w(iskelout),w(iskelin),w(istore2),lenstore2,
     4      lused)
c 
c        perform garbage collection
c 
        call cleascop(w(istore2),w(izsin),lused)
  
        w(8)=izsin+0.1
c 
        keep=lused*2+izsin+10
        lused=istore2+lused+2+10
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine r2d_finskels_retr(w,ibox,ifinner,zs,len)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),zs(2,1)
C 
C        This subroutine retrieves from the array w the skeleton
c        (inner or outer) of the box number ibox. It is expected
c        that the said skeleton had been stored in the said array
c        w via a preceding call to the subroutine r2d_allskels (see).
c        This subroutine has no known uses as a stand-alone device.
c 
c                Input parameterrs:
c 
c 
c  w - array containing the structures created by the subroutine
c        r2d_allskels.
c  ibox - the box for which the subroutine will return the skeleton
c  ifinner - tells the subroutine which skeleton to return:
c       ifinner=1 will cause the inner skeleton to be returned
c       ifinner=0 will cause the outer skeleton to be returned
c 
c               Output parameters:
c 
c  zs - the skeleton
c  len - the number of nodes in the erturned skeleton
c 
        istore2=w(8)
        call r2d_skels_retr(w(istore2),ibox,ifinner,zs,len)
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine r2d_outer_skel3(ier,w,z,iz,nlev,nboxes,
     1      poteval,rk,eps,large_sepa,inband,
     2      store,zsout,zsin,skelout,skelin,store2,lenstore2,
     3      lused)
        implicit real *8 (a-h,o-z)
        save
        integer *4 iz(1),inband(1),box(12)
        dimension w(1),z(2,1),large_sepa(1),store(1),
     1      zsout(2,1),center(3),
     2      corners(2,5),zsin(2,1),store2(1),
     3      skelin(2,1),skelout(2,1)
c 
        complex *16 rk
c 
c        initialize the storage area store2
c 
        call r2d_skels_store_init(store2,lenstore2,nboxes+1)
        lused=nboxes*5+1000
  
c 
c        one level after another, construct outer and inner
c        skeletons and store them in the storage area store2
c 
        do 4000 lev=nlev,1,-1
c 
        do 3600 ibox=1,nboxes
c 
        call r2dretr(ibox,w,box,corners,center)
c 
        if(box(2) .ne. lev) goto 3600
        if(box(9) .eq. -7) goto 3600
c 
        call prinf('in r2d_outer_skel3, lev=*',lev,1)
        call prinf('in r2d_outer_skel3, ibox=*',ibox,1)
c 
c        if this box is childless, its inner skeleton is
c        composed of all nodes inside it. Act accordingly
c 
        if( (box(4) .ne. 0) .or. (box(5) .ne. 0) ) goto 1600
c 
        i1=box(6)-1
        n1=box(7)
        do 1400 i=1,n1
c 
        j=iz(i+i1)
        skelin(1,i)=z(1,j)
        skelin(2,i)=z(2,j)
 1400 continue
c 
        ifinner=1
c 
        call r2d_skels_store(jer,store2,ibox,ifinner,
     1      skelin,n1,lused)
c 
        goto 3600
c 
 1600 continue
c 
c        if this box has only one son - both skeletons (the
c        inner and the outer) of this box are the same as the
c        corresponding skeletons of the sonny. Copy the latter
c        into the former.
c 
        if(box(5) .ne. 0) goto 2200
c 
        ison=box(4)
c 
        ifinner=1
        call r2d_skels_retr(store2,ison,ifinner,skelin,n1)
c 
        call r2d_skels_store(jer,store2,ibox,ifinner,
     1      skelin,n1,lused)
c 
        ifinner=0
        call r2d_skels_retr(store2,ison,ifinner,skelout,n1)
c 
        call r2d_skels_store(jer,store2,ibox,ifinner,
     1      skelout,n1,lused)
c 
        goto 3600
c 
 2200 continue
c 
c        This box has two sonnies. Construct both (inner and outer)
c        of its preliminary skeletons
c 
        iw2=lused*2+10
c 
        call r2d_prelim_skel(ibox,w,z,iz,
     1      large_sepa,inband,store2(iw2),
     2      store,zsout,nn,zsin,lenin,store2)
c 
c        skeletonize the interactions for this box
c 
        iamatr=lused*2+10
        lamatr=nn*lenin*2+100
c 
        iicols=iamatr+lamatr
        licols=nn
        if(lenin .gt. licols) licols=lenin
        licols=licols+10
c 
        iirows=iicols+licols
        lirows=licols
c 
        iw=iirows+lirows
c 
        call r2d_fin_skelit(zsin,lenin,zsout,nn,
     1      poteval,rk,eps,store2(iamatr),skelin,skelout,ncols,
     2      store2(iicols),store2(iirows),store2(iw) )
c 
c        store both skeletons (inner and outer) in array store2
c 
        ifinner=1
        call r2d_skels_store(jer,store2,ibox,ifinner,
     1      skelin,ncols,lused)
c 
        ifinner=0
        call r2d_skels_store(jer,store2,ibox,ifinner,
     1      skelout,ncols,lused)
c 
 3600 continue
 4000 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine r2d_prelim_skel(ibox,w,z,iz,
     1      large_sepa,inband,w2,
     2      store,zsout,nn,zsin,lenin,store2)
        implicit real *8 (a-h,o-z)
        save
        integer *4 iz(1),inband(1),w2(1),box(12)
        dimension w(1),z(2,1),large_sepa(1),store(1),
     1      zsout(2,1),center(3),
     2      corners(2,5),zsin(2,1),store2(1)
c 
c       retrieve from array store2 the inner skeletons
c       of the sonnies of the box ibox, merging them into
c       the preliminary inner skeleton of the box ibox
c 
        call r2dretr(ibox,w,box,corners,center)
  
        ison1=box(4)
        ifinner=1
        call r2d_skels_retr(store2,ison1,ifinner,zsin,n1)
c 
        ison2=box(5)
        if(ison2 .eq. 0) goto 1080
c 
        ifinner=1
        call r2d_skels_retr(store2,ison2,ifinner,
     1      zsin(1,n1+1),n2)
c 
        n1=n1+n2
c 
 1080 continue
c 
        lenin=n1
c 
c       construct the list of boxes in the preliminary outer
c       skeleton of the box ibox, as well as the "loose" nodes
c       in the said skeleton
c 
        call r2d_outer_skel(ibox,w,z,iz,
     1      large_sepa,nsepa,inband,ninband,w2)
c 
c        retrieve from array store the outer frame
c        of the box ibox, and store it in array zsout
c 
        call r2dretr(ibox,w,box,corners,center)
  
        lev=box(2)
c 
cccc        if(lev .le. 3) then
cccc            nn=0
cccc            goto 1160
cccc        endif
c 
        ifinner=0
        call r2d_skels_retr(store,lev,ifinner,zsout,lenout)
c 
        coef=(corners(1,1)-corners(1,2))/2
c 
        do 1150 i=1,lenout
c 
        zsout(1,i)=zsout(1,i)*coef+center(1)
        zsout(2,i)=zsout(2,i)*coef+center(2)
c 
 1150 continue
c 
        nn=lenout
 1160 continue
c 
c 
c        store in array zsout the "loose" nodes
c 
        do 1200 i=1,ninband
c 
        nn=nn+1
        j=iz(inband(i))
        zsout(1,nn)=z(1,j)
        zsout(2,nn)=z(2,j)
 1200 continue
c 
c       retrieve from array store2 the inner skeletons of the
c       boxes in the preliminary outer skeleton of ibox, and
c       return these to the user, together with the "loose" nodes
c       and the outer frame
c 
        if(nsepa .eq. 0) return
c 
        do 1600 i=1,nsepa
c 
        ii=large_sepa(i)
c 
        ifinner=1
        call r2d_skels_retr(store2,ii,ifinner,
     1      zsout(1,nn+1),n3)
c 
        nn=nn+n3
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
        subroutine r2d_fin_skelit(zsin,nin,zsout,nout,
     1      poteval,rk,eps,amatr,skelin,skelout,ncols,
     2      icols,irows,w)
c 
        implicit real *8 (a-h,o-z)
        save
        complex *16 rk,pot,amatr(nin,nout)
c 
        integer *4 irows(1),icols(1)
c 
        dimension zsout(2,1),zsin(2,1),skelin(2,1),
     1      skelout(2,1),w(1)
c 
c        construct the matrix of interactions
c 
        do 1400 i=1,nout
        do 1200 j=1,nin
c 
        call poteval(zsin(1,j),zsout(1,i),rk,pot)
c 
        amatr(j,i)=pot
 1200 continue
 1400 continue
c 
c        . . . skeletonize
c 
        call cskel_basic(amatr,nin,nout,eps,irows,icols,ncols,w)
c 
c       select the skeletons (outer and inner) out of arrays
c       zsin, zsout
c 
        do 1600 i=1,ncols
c 
        j=irows(i)
        skelin(1,i)=zsin(1,j)
        skelin(2,i)=zsin(2,j)
c 
        j=icols(i)
        skelout(1,i)=zsout(1,j)
        skelout(2,i)=zsout(2,j)
 1600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine r2d_ident_matr(a,n)
        implicit complex *16 (a-h,o-z)
        save
        dimension a(n,n)
c 
        do 1400 i=1,n
        do 1200 j=1,n
c 
        a(j,i)=0
 1200 continue
c 
        a(i,i)=1
 1400 continue
c 
        return
        end
  
