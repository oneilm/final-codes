        implicit real *8 (a-h,o-z)
        integer *4 iz(1000 000),laddr(100),map(2,1000),
     1      states(100 000),izout(1000 000),idiffs(100 000),
     2      izout2(1000 000)
        real *8 z(2,1000 000),w(1000 000),zbord(2,20),
     1      z2(2,1000 000),xys(100)
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
cccc        call zquaplot(iw,z,n,2,'random points as created*')
  
  
        lenw=1000 000
        maxboxes=100000
        nbox=30
c 
  
        call r2dstrcr(ier,z,n,
     1       maxboxes,nbox,nall,iz,w,lenw,keep,lused,
     2      map,nlev)
  
        call prinf('after r2dstrcr, ier=*',ier,1)
  
        if(ier .ne. 0) stop
  
        call prin2('after r2dstrcr, map=*',map,0)
c 
 5200 format(6x,5x,i6,5x,i5,3x,i6)
        write(6,5200) (i,map(1,i),map(2,i),i=1,nlev)
        write(13,5200) (i,map(1,i),map(2,i),i=1,nlev)
  
        call prinf('after r2dstrcr, nlev=*',nlev,1)
  
        call prinf('and keep=*',keep,1)
        call prinf('and lused=*',lused,1)
c 
c        plot the whole structure
c 
        iw=25
cccc        call allplot3(iw,w,z)
c 
        do 5400 ibox=2,nall
        iw=500+ibox
cccc        call standard_plot2(iw,ibox,w,z,iz)
c 
 5400 continue
c 
c       test the sifting procedure
c 
        xleft=0.2
        xright=0.8
c 
        ybottom=0.6
        ytop=0.7
  
        ialarm=0
  
        do 8000 ijk=1,100
  
cccc        goto 5500
  
        call RIRAND(4,xys)
  
        call prin2('xys as constructed*',xys,4)
  
        xleft=xys(1)
        if(xys(2) .lt. xys(1)) xleft=xys(2)
c 
        xright=xys(1)
        if(xys(2) .gt. xys(1)) xright=xys(2)
c 
        ybottom=xys(3)
        if(xys(4) .lt. ybottom) ybottom=xys(4)
c 
        ytop=xys(3)
        if(xys(4) .gt. ytop) ytop=xys(4)
  
 5500 continue
  
c 
c        . . . plot the test box with the points
c 
        zbord(1,1)=xright
        zbord(2,1)=ytop
c 
        zbord(1,2)=xleft
        zbord(2,2)=ytop
c 
        zbord(1,3)=xleft
        zbord(2,3)=ybottom
c 
        zbord(1,4)=xright
        zbord(2,4)=ybottom
c 
        zbord(1,5)=xright
        zbord(2,5)=ytop
c 
  
  
c 
        iw=22
cccc        call zquaplot2(iw,z,n,2,zbord,5,3,
cccc     1      'test rectangle + points*')
  
        t1=clotatim()
  
        do 6100 i=1,100
  
        call r2d_points_in(xleft,xright,ybottom,ytop,
     1      w,z,iz,izout,nout,states)
c 
 6100 continue
  
        t2=clotatim()
  
cccc        call prinf('after r2d_points_in, izout=*',izout,nout)
  
        do 6200 i=1,nout
c 
        j=izout(i)
        z2(1,i)=z(1,iz(j))
        z2(2,i)=z(2,iz(j))
 6200 continue
c 
        iw=23
        call zquaplot3(iw,z,n,1,zbord,5,3,z2,nout,2,
     1    'rectangle + points inside+points outside*')
  
        call prinf('and n=*',n,1)
        call prinf('and nout=*',nout,1)
  
  
c 
c        find all points in the rectangle via direct examination
c 
        t3=clotatim()
  
        do 6400 i=1,100
  
        nout2=0
  
        call r2d_onebox_sift(xleft,xright,ybottom,ytop,
     1      iz,z,1,n,izout2,nout2)
  
 6400 continue
        t4=clotatim()
  
        call prinf('and nout2=*',nout2,1)
  
        call sortanyn(izout,nout,states)
  
        call sortanyn(izout2,nout2,states)
  
        do 7200 i=1,nout
c 
        idiffs(i)=izout(i)-izout2(i)
  
  
        if(idiffs(i) .ne. 0) ialarm=1
 7200 continue
  
  
        call prinf('and izout2=*',izout2,nout2)
        call prinf('and again, izout=*',izout,nout)
        call prinf('and idiffs=*',idiffs,nout)
  
        call prin2('by the way, time for r2d_points_in is*',
     1        t2-t1,1)
  
        call prin2('and time for r2d_onebox_sift is*',
     1        t4-t3,1)
c 
  
 8000 continue
  
        call prinf('and ialarm=*',ialarm,1)
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine standard_plot2(iw,ibox,w2,zs,iz)
        implicit real *8 (a-h,o-z)
        save
        integer *4 box(10),iz(1),dad(10),son1(20),son2(20)
        dimension corners(2,10),w(1000000),center(2),
     1      zs(2,1),corners_dad(2,10),center_dad(2),
     2      corn(20),corners_son1(20),center_son1(20),
     3      corners_son2(20),center_son2(20)
        character *100 dummy
c 
c       description of the "standard" box
c 
c    allboxes(1,i) = i
c    allboxes(2,i) - level of box
c    allboxes(3,i) - address of the box's daddy
c    allboxes(4,i) - address of the box's first sonny
c    allboxes(5,i) - address of the box's second sonny
c    allboxes(6,i) - location in array iz of points inside
c       the box
c    allboxes(7,i) - number of points in the box
c    allboxes(10,i) - number of points in the box
c 
c        . . . initialize the plotting routine
c 
         call trplopen(w)
c 
c        retrieve the data about this box
c 
        call r2dretr(ibox,w2,box,corners,center)
  
cccc        call prin2('corners=*',corners,8)
c 
c        plot this box
c 
        call plot1box(corners,w)
c 
c        find this box's daddy and plot it after
c        blowing it up a bit
c 
        idad=box(3)
c 
        call r2dretr(idad,w2,dad,corners_dad,center_dad)
c 
        scale=1.1
        call boxblow(scale,corners_dad,corn)
c 
        call plot1box(corn,w)
c 
c        find the box's sonnies, and plot them after shrinking
c        them appropriately
c 
        ison1=box(4)
        scale=0.9
c 
        if(ison1 .le. 0) goto 2200
c 
        call r2dretr(ison1,w2,son1,corners_son1,center_son1)
c 
        call boxblow(scale,corners_son1,corn)
        call plot1box(corn,w)
c 
        ison2=box(5)
c 
        scale=0.9
c 
        if(ison2 .le. 0) goto 2200
c 
        call r2dretr(ison2,w2,son2,corners_son2,center_son2)
c 
        call boxblow(scale,corners_son2,corn)
        call plot1box(corn,w)
c 
 2200 continue
c 
c        plot the points inside the box ibox
c 
        i1=box(6)
        ni=box(7)
        do 2400 i=1,ni
c 
        j=iz(i1-1+i)
  
        call trplpnt(zs(1,j),zs(2,j),w)
 2400 continue
c 
c       dump the plot on disk
c 
 2600 format(' box number ', i4,' of type ',i1,
     1    ' with dad and kids *')
c 
        itype=box(10)
        write(dummy,2600) ibox,itype
c 
        call trplwrt(iw,w,dummy)
c 
         return
         end
c 
c 
c 
c 
c 
        subroutine boxblow(scale,cornin,cornout)
        implicit real *8 (a-h,o-z)
        save
        dimension cornin(2,4),cornout(2,4),cent(2)
c 
        cent(1)=cornin(1,1)+cornin(1,2)+cornin(1,3)+
     1      cornin(1,4)
        cent(1)=cent(1)/4
c 
        cent(2)=cornin(2,1)+cornin(2,2)+cornin(2,3)+
     1      cornin(2,4)
        cent(2)=cent(2)/4
c 
        cornout(1,1)=cent(1)+(cornin(1,1)-cent(1))*scale
        cornout(1,2)=cent(1)+(cornin(1,2)-cent(1))*scale
        cornout(1,3)=cent(1)+(cornin(1,3)-cent(1))*scale
        cornout(1,4)=cent(1)+(cornin(1,4)-cent(1))*scale
c 
        cornout(2,1)=cent(2)+(cornin(2,1)-cent(2))*scale
        cornout(2,2)=cent(2)+(cornin(2,2)-cent(2))*scale
        cornout(2,3)=cent(2)+(cornin(2,3)-cent(2))*scale
        cornout(2,4)=cent(2)+(cornin(2,4)-cent(2))*scale
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
        character *8 labels
        character *1 quo
        data labels/'labels '/,quo/'"'/
  
  
  
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
        do 2000 i=1,nboxes
c 
        call r2dretr(i,w,box,corners,center)
c 
        call plot1box(corners,w2)
  
        call trplpnt(center(1),center(2),w2)
c 
        centers(1,i)=center(1)
        centers(2,i)=center(2)
        labs(i)=i
 2000 continue
c 
c        create the labels file
c 
  
cccc        iw7=17
cccc        call zqua_labels_int(iw,centers,labs,nboxes)
  
  
  
        call prinf('in allplot3, iw=*',iw,1)
  
cccc        stop
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code, and the beginning
c        of the code for the construction of the "bastardized"
c        binary tree
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c        This file contains three user-callable subroutines:
c        r2dstrcr, r2dretr, and r2d_points_in. Following
c        is a brief description of these subroutines.
c 
c   r2dstrcr - constructs the "bastardized" binary
c        tree, obtained from the standard quad-tree by the
c        introduction of "bastard" boxes, so that every box
c        is subdivided into two sons (one of which might not
c        exist). Thus, there are two types of boxes: square
c        and rectangular; the rectangular boxes are twice as
c        high as they are wide. Please note that square boxes
c        have even-numbered levels of subdivision: 0 (main
c        box), 2, 4,... . Rectangular boxes have odd-numbered
c        levels of subdivision.
c 
c   r2dretr - retrieves from array w (hopefully)
c        produced via a preceding call to the subroutine
c        r2dstrcr (see) the data corresponding to the
c        user-specified box in the structure.
c 
c   r2d_points_in - for the user-specified conventionally oriented
c        rectangle (defined by the coordinates of its sides),
c        this subroutine finds all particles from the
c        user-specified collection that are inside the soid
c        rectangle; the subroutine uses the "bastardized"
c        tree structure constructed during the preceding call
c        to the subroutine r2dstrcr (see).
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c 
c 
c 
        subroutine r2dstrcr(ier,z,n,
     1      maxboxes,nbox,nall,iz,w,lenw,keep,lused,
     2      map,nlev)
        implicit real *8 (a-h,o-z)
        save
        integer *4 iz(1),laddr(100),map(2,1)
        real *8 z(2,1),center0(2),w(1)
c 
c        This subroutine constructs the "bastardized" binary
c        tree, obtained from the standard quad-tree by the
c        introduction of "bastard" boxes, so that every box
c        is subdivided into two sons (one of which might not
c        exist). Thus, there are two types of boxes: square
c        and rectangular; the rectangular boxes are twice as
c        high as they are wide. Please note that square boxes
c        have even-numbered levels of subdivision: 0 (main
c        box), 2, 4,... . Rectangular boxes have odd-numbered
c        levels of subdivision.
c 
c        The principal output of this subroutine is the array w,
c        containing the structure created by the subroutine:
c        definitions of boxes, their corners and centers, and
c        various types of control information. The user is expected
c        to access the contents of this array via the subroutine
c        r2dretr (see). Please note that keep (see below) first
c        elements of this array must remain unchanged between the
c        call to this subroutine and the subsequent calls to the
c        subroutine r2dretr (see).
c 
c                           Input parameters:
c 
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
c  nall - the total number of boxes constructed by the subroutine
c  iz - integer array of length n defining the transposition of
c        the array z (user-specified nodes) such that each of
c        the boxes in the structure contains consecutively numbered
c        points.
c  w - array containing the structure created by the subroutine:
c        definitions of boxes, their corners and centers, and
c        various types of control information. The user is expected
c        to access the contents of this array via the subroutine
c        r2dretr (see). Please note that keep (see below) first
c        elements of this array must remain unchanged between the
c        call to this subroutine and the subsequent calls to the
c        subroutine r2dretr (see).
c  keep - the number of real *8 elements of the array w that must
c        remain unchanged between the call to this subroutine and
c        the subsequent calls to the subroutine r2dretr (see).
c  lused - the total number of real *8 elements in array w that
c        have been used by this subroutine
c  map - a (2,nlev) integer array telling the user where in the
c        imaginary array allboxes the boxes corresponding to
c        various levels of subdivision are located. Specifivally,
c    map(1,i) is the addresss in the array allboxes of the first
c        box on the level i of subdivision;
c    map(2,i) is the addresss in the array allboxes of the last
c        box on the level i of subdivision. PLEASE NOTE THAT THE
C        MAIN (UNSUBDIVIDED) BOZ HAS THE LEVEL 0 (ZERO) OF SUBDI-
C     VISION. IT IS ASSIGNED THE ADDRESS 1 IN THE ARRAY ALLBOXES.
c 
c  nlev - the number of levels of subdivision NOT COUNTING THE
C        LEVEL NUMBER 0.
c 
c        construct the quad-three
c 
        ier=0
        if(lenw .lt. n+1000) then
            ier=1024
            call prinf('bombing from r2dstrcr, with ier=*',
     1          ier,1)
            return
        endif
c 
        ninire=2
c 
        iboxes=21
        lboxes=maxboxes*10/ninire+10
c 
        iw=iboxes+lboxes
c 
        call r2dallb(ier,z,n,nbox,w(iboxes),maxboxes,
     1    nboxes,iz,laddr,nlev,center0,size,w(iw) )
c 
        if(ier .ne. 0) then
c 
            call prinf('
     1       bombing from r2dstrcr, with ier from r2dallb=*',
     2          ier,1)
c 
            return
        endif
c 
c       . . . collect garbage
c 
        lboxes=nboxes*10/ninire+10
c 
        icenters=iboxes+lboxes
        lcenters=nboxes*2+10
c 
        icorners=icenters+lcenters
        lcorners=nboxes*8+10
c 
        iboxes_bast=icorners+lcorners
        lboxes_bast=nboxes*40/ninire+100
c 
        lused=iboxes_bast+lboxes_bast
c 
        if(lused .gt. lenw) then
            ier=256
            call prinf('bombing from r2dstrcr with ier=*',
     1          ier,1)
            return
        endif
c 
c        construct centers and corners for all boxes in the
c        structure
c 
        call r2d_centcorn(center0,size,w(iboxes),nboxes,
     1      w(icenters),w(icorners))
c 
c        restructure the quad-tree, introducing rectangular
c        (as opposed to square) boxes
c 
        call r2dbast0(w(iboxes),nboxes,w(iboxes_bast),nbast,
     1      w(icenters),w(icorners))
c 
        lboxes_bast=nbast*10/ninire+100
c 
        icenters_bast=iboxes_bast+lboxes_bast
        lcenters_bast=nboxes*4+10
c 
        icorners_bast=icenters_bast+lcenters_bast
        lcorners_bast=nboxes*16+100
c 
c       construct the corners for all bastard boxes
c 
        lused=icorners_bast+lcorners_bast
c 
        if(lused .gt. lenw) then
            ier=128
            call prinf('bombing from r2dstrcr with ier=*',
     1          ier,1)
            return
        endif
c 
        call r2d_bast_corn(w(iboxes_bast),nbast,
     1      w(icorners),w(icorners_bast),w(icenters_bast))
c 
c        construct the "standard" boxes
c 
        nall=nboxes+nbast
c 
        iallboxes=icorners_bast+lcorners_bast
        lallboxes=nall*10/ninire+100
c 
        iallcorners=iallboxes+lallboxes
        lallcorners=nall*8+100
c 
        iallcenters=iallcorners+lallcorners
        lallcenters=nall*2+100
c 
        lused=iallcenters+lallcenters
        if(lused .gt. lenw) then
            ier=32
            call prinf('bombing from r2dstrcr with ier=*',
     1          ier,1)
            return
        endif
c 
        call r2d_boxes_standard(w(iboxes),nboxes,
     1      w(iboxes_bast),nbast,w(icorners),w(icorners_bast),
     2      w(iallboxes),w(iallcorners),w(iallcenters) )
c 
cccc        call prinf('w(iallboxes)=*',w(iallboxes),nall*10)
c 
c        perform garbage collection
c 
        iallboxes2=21
        call r2d_reacopy(w(iallboxes),w(iallboxes2),lallboxes)
        iallboxes=iallboxes2
c 
        iallcorners2=iallboxes+lallboxes
        call r2d_reacopy(w(iallcorners),w(iallcorners2),
     1      lallcorners)
        iallcorners=iallcorners2
c 
        iallcenters2=iallcorners+lallcorners
        call r2d_reacopy(w(iallcenters),w(iallcenters2),
     1      lallcenters)
        iallcenters=iallcenters2
c 
        keep=iallcenters+lallcenters
c 
        w(2)=iallboxes+0.1
        w(3)=iallcorners+0.1
        w(4)=iallcenters+0.1
        w(5)=nall+0.1
        w(6)=n+0.1
c 
c        finally, construct the map of the array allboxes
c 
        call r2dmap_bld(w(iallboxes),nall,map,nlev,nbox)
c 
        w(11)=nlev+0.1
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine r2dretr(ibox,w,box,corners,center)
        implicit real *8 (a-h,o-z)
        save
        integer *4 box(10)
        dimension center(2),corners(2,4),w(1)
c 
c        This subroutine retrieves from array w (hopefully)
c        produced via a preceding call to the subroutine
c        r2dstrcr (see) the data corresponding to the
c        user-specified box in the structure.
c 
c                     Input parameters:
c 
c  ibox - the sequence number in the structure of the box whose
c        description is to be returned
c  w - array containing the structure created by the subroutine
c        r2dstrcr: definitions of boxes, their corners and
c        centers, control information.
c 
c                     Output parameters:
c 
c  box - a 10-element integer array containing the description
c        of the box number ibox, a follows.
c 
c    box(1) = ibox
c    box(2) - level of subdivision of box number ibox
c    box(3) - address of the box's daddy
c    box(4) - address of the box's first sonny
c    box(5) - address of the box's second sonny
c    box(6) - location in array iz of points inside
c       the box
c    box(7) - number of points in the box
c    box(9) - the parameter indicating the status of the box:
c      box(9)=0 means that this box exists
c      box(9)=-7 means that this box has been deleted during
c              post-processing
c    box(10) - type of box (its paternity and its location
c        within its daddy, as follows:
c 
c      box(10)=1 means that this is a legitimate box in the lower
c               position within its daddy (the daddy is, obviously,
c               a bastard)
c      box(10)=2 means that this is a legitimate box in the upper
c               position within its daddy (the daddy is, obviously,
c               a bastard)
c      box(10)=3 means that this is a bastard box in the left
c               position within its daddy (the daddy is, obviously,
c               a legitimate box)
c      box(10)=4 means that this is a bastard box in the right
c               position within its daddy (the daddy is, obviously,
c               a legitimate box)
c  corners - the corners of the box number ibox, starting with the
c      top right corner, and moving counterclockwise
c  center - the center of the box number ibox.
c 
c 
c        . . . retrieve from the beginning of array w the locations
c              in it of the arrays boxes, corners, centers
c 
        iboxes=w(2)
        icorners=w(3)
        icenters=w(4)
c 
c       retrieve all the information
c 
        call r2dretr0(ibox,w(iboxes),w(icorners),
     1      w(icenters),box,corners,center)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine r2dretr0(ibox,boxes,all_corners,centers,
     1      box,corners,center)
        implicit real *8 (a-h,o-z)
        save
        integer *4 box(10),boxes(10,1)
        dimension center(2),corners(2,4),centers(2,1),
     1      all_corners(2,4,1)
c 
        call r2d_intcopy(boxes(1,ibox),box,10)
        call r2d_reacopy(all_corners(1,1,ibox),corners,8)
        call r2d_reacopy(centers(1,ibox),center,2)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine r2d_points_in(xleft,xright,ybottom,ytop,
     1      w,z,iz,izout,nout,states)
        implicit real *8 (a-h,o-z)
        save
        integer *4 iz(1),box(12),izout(1)
        dimension corners(2,4),w(1),z(2,1),center(3),
     1      zsout(2,1),states(1)
c 
c        for the user-specified conventionally oriented
c        rectangle (defined by the coordinates of its sides),
c        this subroutine finds all particles from the
c        user-specified collection; that are insied the said
c        rectangle. The subroutine uses the "bastardized"
c        tree structure constructed during the preceding
c        call to the subroutine r2dstrcr (see).
c 
c                     Input parameters:
c 
c  xleft - the x-coordinate of the left (vertical) side of
c        the rectangle
c  xright - the x-coordinate of the right (vertical) side of
c        the rectangle
c  ybottom - the y-coordinate of the bottom (horizontal) side
c        of the rectangle
c  ytop - the y-coordinate of the top (horizontal) side
c        of the rectangle
c  w - the array containing the tree structure corresponding
c        to the collection of points z (constructed via a
c        preceding call to the subroutine r2dstrcr (see)
c  z - the collection of points to be scanned
c  iz - integer array of length n defining the transposition of
c        the array z (user-specified nodes) such that each of
c        the boxes in the structure contains consecutively numbered
c        points (constructed via a preceding call to the subroutine
c        r2dstrcr (see)
c 
c                      Output parameters:
c 
c  izout - the subset of elements of the array iz indexing the
c        points in array z that are inside the user-specified
c        rectangle
c  nout - the number of elements in array izout
c 
c                      Work arrays:
c 
c  states - must contain at least nboxes integer *4 elements, with
c        nboxes the total number of boxes in the structure (as
c        returned by the subroutine r2dstrcr (see)
c 
c 
c        . . . if the whole big box is inside the user-specified
c              rectangle - copy all the points and exit
c 
        nboxes=w(5)
c 
        call r2dretr(1,w,box,corners,center)
c 
        call r2d_relloc(xleft,xright,ybottom,ytop,
     1      corners,ifinside)
c 
        nout=0
        if(ifinside .eq. -1) return
        if(ifinside .eq. 1) then
            nout=box(7)
c 
cccc            call r2d_intcopy(iz,izout,nout)
c 
        do 1100 i=1,nout
        izout(i)=i
 1100 continue
c 
            return
        endif
c 
c        the big box is neither inside nor outside the
c        user-supplied rectangle - perform recursive
c        search
c 
c        . . . initialize the state array
c 
        call r2d_reazero(states,nboxes)
c 
        ibox=1
c 
c       . . . perform the recursive search
c 
        do 5000 ijk=1,10 000 000
c 
c        if the box number is equal to zero - terminate
c        the search
c 
         if(ibox .le. 0) return
c 
c        if this box has not been analyzed - analyze it
c 
        call r2dretr(ibox,w,box,corners,center)
        ist=states(ibox)
c 
        if(ist .ne. 0) goto 3000
c 
        call r2d_relloc(xleft,xright,ybottom,ytop,
     1      corners,ifinside)
c 
c        if this box is completely inside the rectangle
c        - copy its contents into the output array and
c        transfer control to its daddy
c 
        if(ifinside .ne. 1) goto 1400
c 
        ic=box(6)
        nc=box(7)
c 
cccc        call r2d_intcopy(iz(ic),izout(nout+1),nc)
  
        do 1200 j=1,nc
c 
        izout(nout+j)=ic+j-1
 1200 continue
  
  
        nout=nout+nc
        states(ibox)=999
        ibox=box(3)
c 
        goto 5000
c 
 1400 continue
c 
c        if this box is completely outside the rectangle
c        -transfer control to its daddy
c 
        if(ifinside .ne. -1) goto 1600
c 
        states(ibox)=999
        ibox=box(3)
        goto 5000
c 
 1600 continue
c 
c        This box is neither fully inside the rectangle
c        nor fully outside it. If this box is childless
c        - process it individually
c 
        if(box(4) .ne. 0) goto 1800
c 
        ic=box(6)
        nc=box(7)
  
        call r2d_onebox_sift(xleft,xright,ybottom,ytop,
     1      iz,z,ic,nc,izout,nout)
c 
        states(ibox)=999
        ibox=box(3)
        goto 5000
c 
 1800 continue
c 
c        This box is neither fully inside nor fully outside,
c        and it has children. transfer control to the first
c        son
c 
        states(ibox)=1
        ibox=box(4)
        goto 5000
c 
 3000 continue
c 
c        Control has been transfered to this box from one
c        of its sons. If it was the second son, or the
c        the second son does not exist, transfer control
c        to the daddy
c 
        if( (ist .eq. 2) .or. (box(5) .eq. 0) ) then
c 
            states(ibox)=999
            ibox=box(3)
            goto 5000
        endif
c 
c        Control has been transfered to this box from its
c        of first son, and the second son exists. Transfer
c        control the second sonny
c 
            states(ibox)=2
            ibox=box(5)
            goto 5000
c 
 5000 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine r2d_relloc(xleft,xright,ybottom,ytop,
     1      corners,ifinside)
        implicit real *8 (a-h,o-z)
        save
        dimension corners(2,4)
c 
c       determine whether the box specified by its corners
c       is all inside the rectangle specified by its sides,
c       all outside the said rectangle, or neither
c 
        x1=corners(1,2)
        x2=corners(1,1)
c 
        y1=corners(2,3)
        y2=corners(2,1)
c 
        ifinside=0
c 
c       . . . check if it is all inside
c 
        if( (x1 .ge. xleft) .and. (x2 .le. xright) .and.
     1      (y1 .ge. ybottom) .and. (y2 .le. ytop ) ) then
c 
            ifinside=1
            return
        endif
c 
c       check if it is all outside
c 
        if( (x2 .le. xleft) .or. (x1 .ge. xright) .or.
     1      (y2 .le. ybottom) .or. (y1 .ge. ytop ) )
     2      ifinside=-1
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine r2d_onebox_sift(xleft,xright,ybottom,ytop,
     1      iz,z,ic,nc,izout,nout)
        implicit real *8 (a-h,o-z)
        save
        integer *4 iz(1),izout(1)
        dimension z(2,1)
c 
c        Sift through nc elements in array iz, starting
c        with element number ic. Store in array izout
c        those whose corresponding elements in array z
c        are inside the rectangle
c 
        do 1200 i1=1,nc
c 
        i=ic+i1-1
        j=iz(i)
        x=z(1,j)
        y=z(2,j)
c 
        if( (xleft .gt. x) .or. (xright .lt. x) .or.
     1      (ybottom .gt. y) .or. (ytop .lt. y) ) goto 1200
c 
        nout=nout+1
cccc        izout(nout)=iz(i)
        izout(nout)=i
 1200 continue
  
         return
         end
c 
c 
c 
c 
c 
        subroutine r2dmap_bld(boxes,nboxes,map,nlev,nbox)
        implicit real *8 (a-h,o-z)
        save
        integer *4 boxes(10,1),map(2,1)
c 
c       scan the array boxes, entering the beginning and
c       the end of each level
c 
        nlev=0
c 
        do 1400 i=2,nboxes
c 
        k=boxes(2,i)
        if(k .ne. nlev) nlev=k
c 
        if(boxes(2,i-1) .ne. k) map(1,k)=i
        if(i .eq. nboxes) map(2,k)=i
c 
        if( (i .ne. nboxes) .and. (boxes(2,i+1) .ne. k) )
     1      map(2,k)=i
  
cccc        goto 1400
  
c 
c        if this box has children and fewer than nbox particles
c        - kill its kids
c 
        if( (boxes(4,i) .gt. 0) .and. (boxes(7,i) .lt. nbox) ) then
c 
            ison1=boxes(4,i)
            ison2=boxes(5,i)
            boxes(4,i)=0
            boxes(5,i)=0
c 
            if(ison1 .gt. 0) boxes(9,ison1)=-7
            if(ison2 .gt. 0) boxes(9,ison2)=-7
        endif
c 
 1400 continue
c 
        nlev=nlev+1
  
cccc        call prinf('exiting r2dmap_bld, boxes=*',boxes,10*nboxes)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine r2d_boxes_standard(boxes,nboxes,
     1      boxes_bast,nbast,corners,corners_bast,
     2      allboxes,allcorners,allcenters)
        implicit real *8 (a-h,o-z)
        save
        integer *4 boxes(10,1),boxes_bast(10,1),
     1      allboxes(10,1)
        dimension corners(2,4,1),corners_bast(2,4,1),
     1      allcorners(2,4,1),allcenters(2,1)
c 
c        This subroutine converts the arrays boxes, boxes_bast
c        describing respectively the legitimate boxes and the
c        bastard boxes into the array allboxes, describing both
c        types of boxes in a single uniform format. It also
c        merges the arrays of corners of the two types of boxes
c        into a single array corners, and constructs the array
c        centers of all boxes in array allboxes.
c 
c                     Input parameters:
c 
c  boxes - an integer array dimensioned (10,nboxes). Each 10-element
c        column describes one box, as follows:
c 
c       1. level - the level of subdivision on which this box
c             was constructed;
c       2, 3.  - the coordinates of this box among  all
c             boxes on this level
c       4 - the daddy of this box, identified by it address
c             in array boxes
c       5,6,7,8 - the list of children of this box (four of them,
c             and the child is identified by its address in the
c             array boxes; if a box has only one child, only the
c             first of the four child entries is non-zero, etc.)
c       9 - the location in the array iz of the particles
c             living in this box
c       10 - the number of particles living in this box
c 
c       As as rule, this array will have been produced by a prior
c       call to the subroutine r2dallb (see).
c 
c   PLEASE NOTE THAT THE ARRAY BOXES IS DESTROYED BY THIS SUBROUTINE!!!
c 
c  nboxes - the number of columns in array boxes
c  boxes_bast - integer array dimensioned (10,nbast). Each 10-element
c        column describes one bastard box, as follows:
c 
c  (boxes_bast(2,i) - the lavel of subdivision of the i-th
c        bastard's daddy
c  (boxes_bast(3,i) - the addresse in the array boxes of the i-th
c        bastard's daddy
c  (boxes_bast(4,i) - the addresse in the array boxes of the i-th
c        bastard's first sonny
c  (boxes_bast(5,i) - the addresse in the array boxes of the i-th
c        bastard's second sonny
c  (boxes_bast(6,i) - the location in array is of the points living
c        in the i-th bastard
c  (boxes_bast(7,i) - the number of points living in the i-th bastard
c 
c  nbast - the number of bastard boxes described in array boxes_bast
c  corners - the array of corners of all legitimate boxes; for each
c      box, the corners start with the top right corner, and
c      move counterclockwise
c  corners_bast - the array of corners of all bastard boxes; for each
c      box, the corners start with the top right corner, and
c      move counterclockwise
c 
c                        Output parameters:
c 
c  allboxes - integer array dimensioned (10,nboxes+nbast). The
c      10-element column number ibox describes the box, number
c      ibox, as follows:
c 
c    allboxes(1,ibox) = ibox
c    allboxes(2,ibox - level of subdivision of box number ibox
c    allboxes(3,ibox) - address of the box's daddy
c    allboxes(4,ibox) - address of the box's first sonny
c    allboxes(5,ibox) - address of the box's second sonny
c    allboxes(6,ibox) - location in array iz of points inside
c       the box
c    allboxes(7,ibox) - number of points in the box
c      box(9)=0 means that this box exists
c      box(9)=-7 means that this box has been deleted during
c              post-processing
c    allboxes(10,ibox) - type of box (its paternity and its location
c        within its daddy, as follows:
c 
c      allboxes(10,ibox)=1 means that this is a legitimate box in the lower
c               position within its daddy (the daddy is, obviously,
c               a bastard)
c      allboxes(10,ibox)=2 means that this is a legitimate box in the upper
c               position within its daddy (the daddy is, obviously,
c               a bastard)
c      allboxes(10,ibox)=3 means that this is a bastard box in the left
c               position within its daddy (the daddy is, obviously,
c               a legitimate box)
c      allboxes(10,ibox)=4 means that this is a bastard box in the right
c               position within its daddy (the daddy is, obviously,
c               a legitimate box)
c 
c  allcorners - the array of corners of all boxes; for each
c      box, the corners start with the top right corner, and
c      move counterclockwise
c  allcenters - the centers of all boxes
c 
c 
c        . . . for each of the bastards, find its nonbastard daddy
c              and his legitimate sons, and enter the information
c              in their entries in array boxes
c 
        nz=(nboxes+nbast)*10
        call r2d_intzero(allboxes,nz)
c 
        do 1200 i=1,nboxes
c 
        call r2d_intzero(boxes(4,i),5)
 1200 continue
c 
        do 1400 i=1,nbast
c 
        idad=boxes_bast(3,i)
        if(boxes(5,idad) .ne. 0) boxes(6,idad)=i
        if(boxes(5,idad) .eq. 0) boxes(5,idad)=i
c 
        ison1=boxes_bast(4,i)
        boxes(4,ison1)=i
c 
        ison2=boxes_bast(5,i)
        if(ison2 .gt. 0) boxes(4,ison2)=i
c 
 1400 continue
c 
c       copy the "boxes" information to the array allboxes,
c       in the process converting it into "standard box"
c       format
c 
        do 2000 i=1,nboxes
c 
        call r2d_intzero(allboxes(1,i),10)
c 
        allboxes(1,i)=i
        allboxes(2,i)=boxes(1,i)*2
c 
        allboxes(3,i)=boxes(4,i)+nboxes
        if(i .eq. 1) allboxes(3,i)=0
c 
        if(boxes(5,i) .gt. 0)
     1      allboxes(4,i)=boxes(5,i)+nboxes
c 
        if(boxes(6,i) .ne. 0)
     1      allboxes(5,i)=boxes(6,i)+nboxes
c 
        allboxes(6,i)=boxes(9,i)
        allboxes(7,i)=boxes(10,i)
c 
 2000 continue
c 
c       copy the "boxes_bast" information to the array
c       allboxes, in the process converting it into
c       "standard box" format
c 
        do 2400 i=1,nbast
c 
        j=i+nboxes
        allboxes(1,j)=j
c 
        idad=boxes_bast(3,i)
        allboxes(2,j)=boxes(1,idad)*2+1
c 
        allboxes(3,j)=idad
c 
        allboxes(4,j)=boxes_bast(4,i)
        allboxes(5,j)=boxes_bast(5,i)
c 
        allboxes(6,j)=boxes_bast(6,i)
        allboxes(7,j)=boxes_bast(7,i)
c 
 2400 continue
c 
c       merge the arrays of corners of normal boxes
c       and bastards, obtaining the corners for all
c 
        nc=nboxes*8
        call r2d_reacopy(corners,allcorners,nc)
c 
        ncb=nbast*8
        call r2d_reacopy(corners_bast,
     1      allcorners(1,1,nboxes+1),ncb)
c 
c        construct the centers of all boxes
c 
        nall=nboxes+nbast
c 
        do 2600 i=1,nall
c 
        d=allcorners(1,1,i)+allcorners(1,2,i)+
     1      allcorners(1,3,i)+allcorners(1,4,i)
c 
        allcenters(1,i)=d/4
c 
        d=allcorners(2,1,i)+allcorners(2,2,i)+
     1      allcorners(2,3,i)+allcorners(2,4,i)
c 
        allcenters(2,i)=d/4
 2600 continue
c 
c        determine the types of all boxes with respect to
c        parents, and store the obtained information in
c        array allboxes
c 
        do 2800 i=2,nall
c 
        call r2d_box_type(i,allboxes,allcenters,itype)
c 
        allboxes(10,i)= itype
 2800 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine r2d_box_type(ibox,boxes,centers,itype)
        implicit real *8 (a-h,o-z)
        save
        integer *4 boxes(10,1)
        dimension centers(2,1)
c 
c        dtermine whether it is a "legitimate" box or a "bastard"
c 
        i=boxes(2,ibox)
        j=i/2
        j=i-j*2
c 
        ifbast=0
        if(j .ne. 0) ifbast=1
c 
c       determine the kind of this box if it is "legitimate"
c 
        if(ifbast .eq. 1) goto 2000
c 
        idad=boxes(3,ibox)
        itype=2
        if(centers(2,idad) .ge. centers(2,ibox)) itype=1
        return
c 
 2000 continue
c 
c       this box is a bastard. act accordingly
c 
        itype=4
        idad=boxes(3,ibox)
        if(centers(1,idad) .ge. centers(1,ibox)) itype=3
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine r2d_bast_corn(boxes_bast,nbast,
     1      corners,corners_bast,centers_bast)
        implicit real *8 (a-h,o-z)
        save
        integer *4 boxes_bast(10,1)
        dimension corners(2,4,1),corners_bast(2,4,1),
     1      xs(10),ys(10),centers_bast(2,1)
c 
c        scan the array boxes_bast, constructing the corners
c 
        do 3000 i=1,nbast
c 
c       construct the corners for this bastard if this
c       bastard has two children
c 
        if(boxes_bast(5,i) .le. 0) goto 2000
c 
        is1=boxes_bast(4,i)
        is2=boxes_bast(5,i)
c 
        do 1200 j=1,4
c 
        xs(j)=corners(1,j,is1)
        xs(j+4)=corners(1,j,is2)
c 
        ys(j)=corners(2,j,is1)
        ys(j+4)=corners(2,j,is2)
 1200 continue
c 
        call r2dminmax(xs,8,xmin,xmax)
        call r2dminmax(ys,8,ymin,ymax)
c 
        corners_bast(1,1,i)=xmax
        corners_bast(2,1,i)=ymax
c 
        corners_bast(1,2,i)=xmin
        corners_bast(2,2,i)=ymax
c 
        corners_bast(1,3,i)=xmin
        corners_bast(2,3,i)=ymin
c 
        corners_bast(1,4,i)=xmax
        corners_bast(2,4,i)=ymin
c 
        goto 2600
  
 2000 continue
c 
c        this bastard has only one child; find its corners
c 
        is1=boxes_bast(4,i)
        do 2200 j=1,4
c 
        xs(j)=corners(1,j,is1)
        ys(j)=corners(2,j,is1)
 2200 continue
c 
        call r2dminmax(xs,4,xmin,xmax)
        call r2dminmax(ys,4,ymin,ymax)
c 
        j=boxes_bast(10,i)
c 
        if( (j .eq. 1) .or. (j .eq. 3) ) ymax=2*ymax-ymin
        if( (j .eq. 2) .or. (j .eq. 4) ) ymin=2*ymin-ymax
c 
        corners_bast(1,1,i)=xmax
        corners_bast(2,1,i)=ymax
c 
        corners_bast(1,2,i)=xmin
        corners_bast(2,2,i)=ymax
c 
        corners_bast(1,3,i)=xmin
        corners_bast(2,3,i)=ymin
c 
        corners_bast(1,4,i)=xmax
        corners_bast(2,4,i)=ymin
c 
 2600 continue
c 
c       construct the bastard's center
c 
        centers_bast(1,i)=(xmin+xmax)/2
        centers_bast(2,i)=(ymin+ymax)/2
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
        subroutine r2dminmax(arr,n,arrmin,arrmax)
        implicit real *8 (a-h,o-z)
        save
        dimension arr(1)
c 
        arrmin=arr(1)
        arrmax=arr(1)
c 
        do 1200 i=2,n
c 
        if(arr(i) .gt. arrmax) arrmax=arr(i)
        if(arr(i) .lt. arrmin) arrmin=arr(i)
c 
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine r2dbast0(boxes,nboxes,boxes_bast,nbast,
     1      centers,corners)
       implicit real *8 (a-h,o-z)
        save
        dimension centers(2,1),corners(2,4,1)
        integer *4 bast1(12),bast2(12),boxes(10,1),
     1      boxes_bast(10,1),isons(10),iwhich(10)
c 
c        note that we are using a strange numbering convention
c        for the children:
c 
c             2     4
c 
c             1     3
c 
c        The entries in the array boxes_bast are as follows:
c 
c  (boxes_bast(2,i) - the lavel of subdivision of the i-th
c        bastard's daddy
c  (boxes_bast(3,i) - the addresse in the array boxes of the i-th
c        bastard's daddy
c  (boxes_bast(4,i) - the addresse in the array boxes of the i-th
c        bastard's first sonny
c  (boxes_bast(5,i) - the addresse in the array boxes of the i-th
c        bastard's second sonny
c  (boxes_bast(6,i) - the location in array is of the points living
c        in the i-th bastard
c  (boxes_bast(7,i) - the number of points living in the i-th bastard
c 
c 
c        scan the array boxes created by r2dallb, introducing
c        intermediate levels of subdivision (consisting of
c        rectangular boxes, otherwise known as bastards)
c 
        ibast=0
c 
        do 3000 k=1,nboxes
c 
c       if this box is childless - do nothing
c 
        if(boxes(5,k) .eq. 0) goto 3000
c 
c        This box has children. Generate the intermediate
c        (rectangular) children, otherwise known as bastards,
c        converting the former into grandchildren
c 
        call r2d_intzero(bast1,10)
        call r2d_intzero(bast2,10)
c 
c        construct the table of children that the k-th box has,
c        and the table of their locations within the k-th box
c 
        isons(1)=boxes(5,k)
        isons(2)=boxes(6,k)
        isons(3)=boxes(7,k)
        isons(4)=boxes(8,k)
c 
        do 1400 i=1,4
        if(isons(i) .eq. 0) goto 1600
c 
        nsons=i
        call r2d_whichkid(isons(i),k,centers,iwhich(i))
 1400 continue
 1600 continue
c 
c        construct the first bastard, if it exists
c 
        ifbast1=0
        nkids1=0
c 
        do 2000 i=1,nsons
c 
        if(iwhich(i) .eq. 1) then
c 
            nkids1=1
c 
            bast1(2)=boxes(1,k)
            bast1(3)=k
            bast1(4)=isons(i)
            bast1(6)=boxes(9,isons(i))
            bast1(7)=boxes(10,isons(i))
            bast1(10)=1
c 
        endif
c 
 2000 continue
c 
        do 2200 i=1,nsons
c 
        if(iwhich(i) .eq. 2) then
c 
            nkids1=nkids1+1
c 
            bast1(2)=boxes(1,k)
            bast1(3)=k
c 
            if(nkids1 .eq. 1) bast1(4)=isons(i)
            if(nkids1 .eq. 2) bast1(5)=isons(i)
c 
            if(nkids1 .eq. 1) bast1(6)=boxes(9,isons(i))
  
            bast1(7)=bast1(7)+boxes(10,isons(i))
            bast1(10)=2
c 
        endif
c 
 2200 continue
c 
c        construct the second bastard, if it exists
c 
        ifbast2=0
        nkids2=0
c 
        do 2600 i=1,nsons
c 
        if(iwhich(i) .eq. 3) then
c 
            nkids2=1
c 
            bast2(2)=boxes(1,k)
            bast2(3)=k
            bast2(4)=isons(i)
            bast2(6)=boxes(9,isons(i))
            bast2(7)=boxes(10,isons(i))
            bast2(10)=3
c 
        endif
c 
 2600 continue
c 
        do 2800 i=1,nsons
c 
        if(iwhich(i) .eq. 4) then
c 
            nkids2=nkids2+1
c 
            bast2(2)=boxes(1,k)
            bast2(3)=k
c 
            if(nkids2 .eq. 1) bast2(4)=isons(i)
            if(nkids2 .eq. 2) bast2(5)=isons(i)
c 
            if(nkids2 .eq. 1) bast2(6)=boxes(9,isons(i))
            bast2(7)=bast2(7)+boxes(10,isons(i))
            bast2(10)=4
c 
        endif
c 
 2800 continue
c 
c        store the information about the generated bastards
c        in the array boxes_bast
c 
        if(nkids1 .ne. 0) then
            ibast=ibast+1
            bast1(1)=ibast
            call r2d_intcopy(bast1,boxes_bast(1,ibast),10)
        endif
c 
        if(nkids2 .ne. 0) then
            ibast=ibast+1
            bast2(1)=ibast
            call r2d_intcopy(bast2,boxes_bast(1,ibast),10)
        endif
 3000 continue
c 
        nbast=ibast
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine r2d_whichkid(ison,idad,centers,iwhich)
        implicit real *8 (a-h,o-z)
        save
        dimension centers(2,1)
  
        if( (centers(1,ison) .le. centers(1,idad)) .and.
     1      (centers(2,ison) .le. centers(2,idad)) ) then
c 
            iwhich=1
            return
        endif
c 
c 
        if( (centers(1,ison) .le. centers(1,idad)) .and.
     1      (centers(2,ison) .ge. centers(2,idad)) ) then
c 
            iwhich=2
            return
        endif
c 
c 
        if( (centers(1,ison) .ge. centers(1,idad)) .and.
     1      (centers(2,ison) .ge. centers(2,idad)) ) then
c 
            iwhich=4
            return
        endif
c 
c 
        if( (centers(1,ison) .ge. centers(1,idad)) .and.
     1      (centers(2,ison) .le. centers(2,idad)) ) then
c 
            iwhich=3
            return
        endif
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine r2d_centcorn(center0,size,boxes,nboxes,
     1      centers,corners)
        implicit real *8 (a-h,o-z)
        save
        integer *4 boxes(10,1)
        dimension centers(2,1),corners(2,4,1),center(2),
     1      center0(2)
c 
c       this subroutine produces arrays of centers and
c       corners for all boxes in the quad-tree structure.
c 
c              input parameters:
c 
c  center0 - the center of the box on the level 0, containing
c         the whole simulation
c  size - the side of the box on the level 0
c  boxes - an integer array dimensioned (10,nboxes), as produced
c        by the subroutine r2dallb (see)
c        column describes one box, as follows:
c  nboxes - the total number of boxes created
c 
c 
c              output parameters:
c 
c  centers - the centers of all boxes in the array boxes
c  corners - the corners of all boxes in the array boxes
c 
c       . . . construct the corners for all boxes
c 
        x00=center0(1)-size/2
        y00=center0(2)-size/2
        do 1400 i=1,nboxes
        level=boxes(1,i)
        side=size/2**level
        side2=side/2
        ii=boxes(2,i)
        jj=boxes(3,i)
        center(1)=x00+(ii-1)*side+side2
        center(2)=y00+(jj-1)*side+side2
c 
        centers(1,i)=center(1)
        centers(2,i)=center(2)
c 
        corners(1,1,i)=center(1)+side/2
        corners(1,2,i)=corners(1,1,i)-side
        corners(1,3,i)=corners(1,2,i)
        corners(1,4,i)=corners(1,3,i)+side
c 
        corners(2,1,i)=center(2)+side/2
        corners(2,2,i)=corners(2,1,i)
        corners(2,3,i)=corners(2,2,i)-side
        corners(2,4,i)=corners(2,3,i)
 1400 continue
         return
         end
c 
c 
c 
c 
c 
        subroutine r2dallb(ier,z,n,nbox,boxes,maxboxes,
     1    nboxes,iz,laddr,nlev,center0,size,iwork)
        implicit real *8 (a-h,o-z)
        save
        integer *4 boxes(10,1),iz(1),laddr(2,1),iwork(1),
     1      is(4),ns(4),iisons(4),jjsons(4)
        real *8 z(2,1),center0(1),center(2)
        data iisons/1,1,2,2/,jjsons/1,2,1,2/
c 
c        this subroutine constructs a quad-tree corresponding
c        to the user-specified collection of points in the plane
c 
c              input parameters:
c 
c  z - the set of points in the plane
c  n - the number of elements in z
c  nbox - the maximum number of points permitted in a box on
c        the finest level. in other words, a box will be further
c        subdivided if it contains more than nbox points.
c  maxboxes - the maximum total number of boxes the subroutine
c        is permitted to create. if the points z are such that
c        more boxes are needed, the error return code ier is
c        set to 4, and the execution of thye subroutine is
c        terminated.
c 
c              output parameters:
c 
c  ier - the error return code.
c    ier=0 means successful execution
c    ier=4 means that the subroutine attempted to create more
c        than maxboxes boxes
c    ier=16 means that the subroutine attempted to construct more
c        than 199 levels of subdivision.
c  boxes - an integer array dimensioned (10,nboxes). each 10-element
c        column describes one box, as follows:
c 
c       1. level - the level of subdivision on which this box
c             was constructed;
c       2, 3.  - the coordinates of this box among  all
c             boxes on this level
c       4 - the daddy of this box, identified by it address
c             in array boxes
c       5,6,7,8 - the list of children of this box (four of them,
c             and the child is identified by its address in the
c             array boxes; if a box has only one child, only the
c             first of the four child entries is non-zero, etc.)
c       9 - the location in the array iz of the particles
c             living in this box
c       10 - the number of particles living in this box
c 
c    important warning: the array boxes has to be dimensioned
c                       at least (10,maxboxes)!! otherwise,
c                       the subroutine is likely to bomb, since
c                       it assumes that it has that much space!!!!
c  nboxes - the total number of boxes created
c  iz - the integer array addressing the particles in
c         all boxes.
c       explanation: for a box ibox, the particles living in
c         it are:
c 
c         (z(1,j),z(2,j)),(z(1,j+1),z(2,j+1)),(z(1,j+2),z(2,j+2)),
c         . . . (z(1,j+nj-1),z(2,j+nj-1)),(z(1,j+nj),z(2,j+nj)),
c 
c         with j=boxes(9,ibox), and nj=boxes(10,ibox)
c 
c  laddr - an integer array dimensioned (2,numlev<), containing
c         the map of array boxes, as follows:
c       laddr(1,i) is the location in array boxes of the information
c         pertaining to level=i-1
c       laddr(2,i) is the number of boxes created on the level i-1
c 
c  nlev - the maximum level number on which any boxes have
c         been created. the maximim number possible is 200.
c         it is recommended that the array laddr above be
c         dimensioned at least (2,200), in case the user underestimates
c         the number of levels required.
c  center0 - the center of the box on the level 0, containing
c         the whole simulation
c  size - the side of the box on the level 0
c 
c                  work arrays:
c 
c  iwork - must be at least n+2 integer*4 elements long.
c 
c        . . . find the main box containing the whole picture
c 
        ier=0
        xmin=1.0d50
        xmax=-xmin
        ymin=1.0d50
        ymax=-ymin
c 
        do 1100 i=1,n
        if(z(1,i) .lt. xmin) xmin=z(1,i)
        if(z(1,i) .gt. xmax) xmax=z(1,i)
        if(z(2,i) .lt. ymin) ymin=z(2,i)
        if(z(2,i) .gt. ymax) ymax=z(2,i)
 1100 continue
        size=xmax-xmin
        sizey=ymax-ymin
        if(sizey .gt. size) size=sizey
c 
        center0(1)=(xmin+xmax)/2
        center0(2)=(ymin+ymax)/2
c 
        boxes(1,1)=0
        boxes(2,1)=1
        boxes(3,1)=1
        boxes(4,1)=0
        boxes(5,1)=0
        boxes(6,1)=0
        boxes(7,1)=0
        boxes(8,1)=0
        boxes(9,1)=1
        boxes(10,1)=n
c 
        laddr(1,1)=1
        laddr(2,1)=1
c 
        do 1200 i=1,n
        iz(i)=i
 1200 continue
c 
c       recursively (one level after another) subdivide all
c       boxes till none are left with more than nbox particles
c 
        maxson=maxboxes
        maxlev=200
        ison=1
        nlev=0
c 
        do 3000 level=0,maxlev
c 
        laddr(1,level+2)=laddr(1,level+1)+laddr(2,level+1)
        nlevson=0
        idad0=laddr(1,level+1)
        idad1=idad0+laddr(2,level+1)-1
c 
        do 2000 idad=idad0,idad1
c 
c       subdivide the box number idad (if needed)
c 
        numpdad=boxes(10,idad)
        if(numpdad .le. nbox) goto 2000
c 
        ii=boxes(2,idad)
        jj=boxes(3,idad)
        call r2dcentf(center0,size,level,ii,jj,center)
c 
        iiz=boxes(9,idad)
        nz=boxes(10,idad)
        call r2dsepa1(center,z,iz(iiz),nz,iwork,
     1    is,ns)
c 
c       store in array boxes the sons obtained by the routine
c       r2dsepa1
c 
         idadson=5
        do 1600 i=1,4
        if(ns(i) .eq. 0) goto 1600
        nlevson=nlevson+1
        ison=ison+1
        nlev=level+1
c 
c       . . . if the user-provided array boxes is too
c             short - bomb out
c 
        if(ison .le. maxson) goto 1400
        ier=4
        return
 1400 continue
c 
c        store in array boxes all information about this son
c 
        boxes(5,ison)=0
        boxes(6,ison)=0
        boxes(7,ison)=0
        boxes(8,ison)=0
c 
        boxes(1,ison)=level+1
        iison=(ii-1)*2+iisons(i)
        jjson=(jj-1)*2+jjsons(i)
        boxes(2,ison)=iison
        boxes(3,ison)=jjson
        boxes(4,ison)=idad
        boxes(9,ison)=is(i)+iiz-1
        boxes(10,ison)=ns(i)
c 
        boxes(idadson,idad)=ison
        idadson=idadson+1
        nboxes=ison
 1600 continue
 2000 continue
        laddr(2,level+2)=nlevson
         if(nlevson .eq. 0) goto 4000
 3000 continue
        ier=16
 4000 continue
        nboxes=ison
        return
        end
c 
c 
c 
c 
c 
        subroutine r2dsepa1(cent,z,iz,n,iwork,
     1    is,ns)
        implicit real *8 (a-h,o-z)
        save
        dimension cent(1),z(2,1),iz(1),iwork(1),is(1),ns(1)
c 
c        this subroutine reorders the particles in a box,
c        so that each of the children occupies a contigious
c        chunk of array iz
c 
c        note that we are using a strange numbering convention
c        for the children:
c 
c             2     4
c 
c             1     3
c 
c 
c        . . . separate all particles in this box horizontally
c 
        n1=0
        n2=0
        n3=0
        n4=0
        itype=1
        thresh=cent(1)
        call r2dsepa0(z,iz,n,itype,thresh,iwork,n12)
        n34=n-n12
c 
c       at this point, the contents of sons number 1,2 are in
c       the part of array iz with numbers 1,2,...n12
c       the contents of sons number 3,4  are in
c       the part of array iz with numbers n12+1,n12+2,...,n
c 
c        . . . seperate the boxes 1, 2, and boxes 3, 4
c 
        itype=2
        thresh=cent(2)
        if(n12 .ne. 0)
     1    call r2dsepa0(z,iz,n12,itype,thresh,iwork,n1)
        n2=n12-n1
c 
        if(n34 .ne. 0)
     1    call r2dsepa0(z,iz(n12+1),n34,itype,thresh,iwork,n3)
        n4=n34-n3
c 
c       store the information about the sonnies in appropriate arrays
c 
        is(1)=1
        ns(1)=n1
c 
        is(2)=is(1)+ns(1)
        ns(2)=n2
c 
        is(3)=is(2)+ns(2)
        ns(3)=n3
c 
        is(4)=is(3)+ns(3)
        ns(4)=n4
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine r2dsepa0(z,iz,n,itype,thresh,iwork,n1)
        implicit real *8 (a-h,o-z)
        save
        dimension z(2,1),iz(1),iwork(1)
c 
c       subdivide the points in this box, horizontally
c       or vertically, depending on the parameter itype
c 
        i1=0
        i2=0
        do 1400 i=1,n
        j=iz(i)
        if(z(itype,j) .gt. thresh) goto 1200
        i1=i1+1
        iz(i1)=j
        goto 1400
c 
 1200 continue
        i2=i2+1
        iwork(i2)=j
 1400 continue
c 
        do 1600 i=1,i2
        iz(i1+i)=iwork(i)
 1600 continue
        n1=i1
        return
        end
c 
c 
c 
c 
c 
        subroutine r2dcentf(center0,size,level,i,j,center)
        implicit real *8 (a-h,o-z)
        save
        dimension center(1),center0(1)
        data level0/-1/
c 
c       this subroutine finds the center of the box
c       number (i,j) on the level level. note that the
c       box on level 0 is assumed to have the center
c       center0, and the side size
c 
        if(level .eq. level0) goto 1200
        side=size/2**level
        side2=side/2
        x0=center0(1)-size/2
        y0=center0(2)-size/2
        level0=level
 1200 continue
        center(1)=x0+(i-1)*side+side2
        center(2)=y0+(j-1)*side+side2
        return
        end
c 
c 
c 
c 
c 
        subroutine r2d_intcopy(ia,ib,n)
        implicit real *8 (a-h,o-z)
        save
        dimension ia(1),ib(1),a(1),b(1)
c 
        do 1200 i=1,n
c 
        ib(i)=ia(i)
 1200 continue
  
        return
c 
c 
c 
c 
        entry r2d_intzero(ia,n)
c 
        do 2200 i=1,n
c 
        ia(i)=0
 2200 continue
c 
        return
c 
c 
c 
c 
        entry r2d_reacopy(a,b,n)
c 
        do 3200 i=1,n
c 
        b(i)=a(i)
 3200 continue
c 
        return
c 
c 
c 
c 
        entry r2d_reazero(a,n)
c 
        do 4200 i=1,n
c 
        a(i)=0
 4200 continue
        return
        end
  
