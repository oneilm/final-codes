        implicit real *8 (a-h,o-z)
        integer *4 iz(1000 000),map(2,1000),
     1      large_sepa(100 000),inband(100000),
     2      w2(100000),box(12)
        real *8 z(2,1000 000),w(1000 000),
     1      corners(2,5),center(3)
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
        call zquaplot(iw,z,n,2,'random points as created*')
  
  
        lenw=1000 000
        maxboxes=100000
        nbox=30
c 
  
        call r2dstrcr_neighb(ier,z,n,
     1       maxboxes,nbox,nall,iz,w,lenw,keep,lused,
     2       map,nlev)
  
  
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
c        print the neighbors
c 
        iboxes=21
        call prinf('boxes are*',w(iboxes),nall*10)
        ineighbors=w(7)
        call neighb_print(nall,w(ineighbors) )
  
  
  
c 
c        plot the whole structure
c 
  
        iw=25
        call allplot3(iw,w,z)
c 
c        construct the outer skeleton for one box
c 
  
        ibox=17
cccc        ibox=125
  
cccc        ibox=1908
  
cccc        stop
  
cccc        do 6400 i=2,nall
        do 6400 i=2,20
c 
        ibox=i
  
        call r2dretr(ibox,w,box,corners,center)
  
        if(box(9) .eq. -7) goto 6400
  
        call r2d_outer_skel(ibox,w,z,iz,
     1      large_sepa,nsepa,inband,ninband,w2)
c 
  
  
  
        iw=100+i
  
        call plot4(iw,w,ibox,large_sepa,nsepa,
     1      iz,inband,ninband,z)
  
 6400 continue
  
        stop
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
  
        theta=3.14/2
  
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
c 
c 
c 
c 
        subroutine r2d_all_neighbors_dumb(w,nall,nlev,neighbors)
        implicit real *8 (a-h,o-z)
        save
        integer *4 neighbors(11,1),boxi(12),boxj(12)
        dimension w(1),centi(3),corni(10),centj(3),cornj(10)
c 
        do 1400 i=1,nall
        do 1200 j=1,11
c 
        neighbors(j,i)=0
 1200 continue
 1400 continue
c 
c        construct the lists of neighbors for all
c        legitimate boxes
c 
        do 2000 i=2,nall
c 
c       if this box is on a wrong level - bypass it
c 
        call r2dretr(i,w,boxi,corni,centi)
c 
        if(boxi(9) .eq. -7) goto 2000
  
  
        ii=boxi(2)
        iii=ii/2
        iii=ii-iii*2
        if(iii .ne. 0) goto 2000
c 
        do 1800 j=1,nall
c 
        if(j .eq. i) goto 1800
c 
        call r2dretr(j,w,boxj,cornj,centj)
c 
        if(boxj(9) .eq. -7) goto 1800
c 
        jj=boxj(2)
        if(jj .ne. ii) goto 1800
c 
c        these boxes are on the same level. Check if they
c        are neighbors
c 
        call r2d_ifneighb(corni,cornj,ifclose)
c 
        if(ifclose .eq. 0) goto 1800
c 
        jjj=neighbors(1,i)+1
        neighbors(1,i)=jjj
        neighbors(jjj+1,i)=j
c 
 1800 continue
 2000 continue
  
cccc          return
c 
c        construct the lists of neighbors for all
c        bastard boxes
c 
        do 3000 i=2,nall
c 
c       if this box is on a wrong level - bypass it
c 
        call r2dretr(i,w,boxi,corni,centi)
c 
        if(boxi(9) .eq. -7) goto 3000
c 
        ii=boxi(2)
        iii=ii/2
        iii=ii-iii*2
        if(iii .eq. 0) goto 3000
c 
        do 2800 j=1,nall
c 
c 
        call r2dretr(j,w,boxj,cornj,centj)
c 
        if(boxj(9) .eq. -7) goto 2800
c 
        jj=boxj(2)
        if(jj .ne. ii+1) goto 2800
        if(boxj(3) .eq. i) goto 2800
c 
c        the j-th box is on the same level as the children of the
c        i-th one. Check if they are neighbors
c 
        call r2d_ifneighb(corni,cornj,ifclose)
c 
        if(ifclose .eq. 0) goto 2800
c 
        jjj=neighbors(1,i)+1
        neighbors(1,i)=jjj
        neighbors(jjj+1,i)=j
c 
 2800 continue
 3000 continue
c 
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code and the beginning
c        of neighborhood construction code proper
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This file contains two user-callable subroutines:
c        r2dstrcr_neighb and r2d_outer_skel. Following is a
c        brief description of these subroutines.
c 
c   r2dstrcr_neighb - an augmentation of the subroutine r2dstrcr,
c        constructing (in addition to the machinery built by r2dstrcr)
c        lists of neighbors for all boxes in the structure.
c 
c   r2d_outer_skel - for the box number ibox, this subroutine uses the
c        various data constructed by a prior call to the
c        subroutine r2dstrcr_neighb (see) to obtain the outer
c        skeleton for the box number ibox. The outer skeleton
c        consists of a bunch of boxes on various levels (higher
c        than the level of ibox), plus a collection of "free"
c        points.
c 
c 
        subroutine r2dstrcr_neighb(ier,z,n,
     1      maxboxes,nbox,nall,iz,w,lenw,keep,lused,
     2      map,nlev)
        implicit real *8 (a-h,o-z)
        save
        integer *4 iz(1),map(2,1)
        real *8 z(2,1),w(1)
c 
c       PLEASE NOTE THAT THIS SUBROUTINE IS A SIMPLE AUGMENTATION
C       OF THE SUBROUTINE R2DSTRCR, WHICH THIS SUBROUTINE CALLS.
C       THE AUGMENTATION CONSISTS OF CONSTRUCTING THE NEIGHBORS
C       FOR EACH OF THE BOXES. NEIGHBORS FOR A LEGITIMATE BOX B ARE
C       (REASONABLY ENOUGH) BOXES ON THE SAME LEVEL THAT TOUCH B.
C       FOR A BASTARD B, NEIGHBORS ARE (LESS REASONABLY) SQUARES
C       ON THE NEXT LEVEL OF SUBDIVISION (WHO THEMSELVES ARE
C       LEGITIMATE BOXES) TOUCHING B. The neighbors are stored in
c       the array w, starting with the element number w(7),and
c       are formatted as an (11,nall) integer *4 array.
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
C        MAIN (UNSUBDIVIDED) BOX HAS THE LEVEL 0 (ZERO) OF SUBDI-
C     VISION. IT IS ASSIGNED THE ADDRESS 1 IN THE ARRAY ALLBOXES.
c 
c  nlev - the number of levels of subdivision NOT COUNTING THE
C        LEVEL NUMBER 0.
c 
c       construct the "bastardized" structure
c 
        call r2dstrcr(jer,z,n,
     1      maxboxes,nbox,nall,iz,w,lenw,keep,lused,
     2      map,nlev)
c 
        if(jer .ne. 0) then
            ier=64
            return
        endif
c 
c        construct the neighbors
c 
        ineighbors=keep+4
c 
        lneighbors=nall*11+100
c 
        iw2=ineighbors+lneighbors
c 
        lw2=3100
c 
        lused2=iw2+lw2
c 
        if(lused .lt. lused2) lused=lused2
c 
        if(lused .gt. lenw) then
            ier=32
            return
        endif
c 
        call r2d_all_neighbors(w,nall,nlev,w(ineighbors),w(iw2) )
c 
        w(7)=ineighbors+0.1
        keep=keep+nall*11+100
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine r2d_outer_skel(ibox,w,z,iz,
     1      large_sepa,nsepa,inband,ninband,w2)
        implicit real *8 (a-h,o-z)
        save
        integer *4 iz(1),inband(1),w2(1)
        dimension w(1),z(2,1),large_sepa(1)
c 
c        For the box number ibox, this subroutine uses the
c        various data constructed by a prior call to the
c        subroutine r2dstrcr_neighb (see) to obtain the outer
c        skeleton for the box number ibox. The outer skeleton
c        consists of a bunch of boxes on various levels (higher
c        than the level of ibox), plus a collection of "free"
c        points.
c 
c                      Input parameters:
c 
c  ibox - the box whose outer skeleton is to be constructed
c  w - array returned by a prior call to the subroutine
c        r2dstrcr_neighb (see)
c  z - the points in R^2 on top of which the tree structure
c        is placed.
c  iz - integer array of length n defining the transposition of
c        the array z (user-specified nodes) such that each of
c        the boxes in the structure contains consecutively numbered
c        points (returned by a prior call to the subroutine
c        r2dstrcr_neighb)
c 
c                       Output parameters:
c 
c  large_sepa - the list of sequence numbers  of boxes
c        constituting (together with the nodes in array inband) the
c        outer skeleton of the box ibox
c  nsepa - the number of elements in array large_sepa
c  inband - the list of sequence numbers (in array iz) of the "free"
c        nodes in the outer skeleton of the box ibox
c  ninband - the number of elements in array inband
c 
c                       Work arrays:
c 
c  w2 - must be at least 4*n+5*nall+1000 integer *4 elements long
c 
c 
        ineighb=w(7)
        call r2d_outer_skel0(ibox,w,z,iz,
     1      large_sepa,nsepa,inband,ninband,
     2      w(ineighb),w2)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine r2d_outer_skel0(ibox,w,z,iz,
     1      large_sepa,nsepa,inband,ninband,
     2      neighbors,w2)
        implicit real *8 (a-h,o-z)
        save
        integer *4 iz(1),box(12),inband(1),w2(1)
        dimension w(1),z(2,1),center(3),neighbors(11,1),
     1      corners(2,5),large_sepa(1)
c 
c       construct the list of all points in the potential
c       interaction list of this box
c 
        nall=w(5)
        n=w(6)
        iizz=1
        lizz=n*2+100
c 
        iww=iizz+lizz
        lww=n*2+100
c 
        istates=iww+lww
        lstates=nall+100
c 
        iiarr=istates+lstates
        liarr=nall*2+100
c 
        call r2dretr(ibox,w,box,corners,center)
c 
        corners(1,5)=corners(1,1)
        corners(2,5)=corners(2,1)
c 
        xleft=2*corners(1,2)-corners(1,1)
        xright=2*corners(1,1)-corners(1,2)
c 
        ybottom=2*corners(2,3)-corners(2,2)
        ytop=2*corners(2,2)-corners(2,3)
c 
        call r2d_points_in(xleft,xright,ybottom,ytop,
     1      w,z,iz,w2(iizz),nband,w2(istates) )
c 
c        throw out the nodes inside the box ibox
c 
        call r2d_sift(corners,
     1      iz,z,w2(iizz),nband,inband,nband2)
c 
        nband=nband2
c 
c       Scan the descendants of the neighbors of the box ibox.
c       Select those that are separated from ibox.
c 
        nneighb=neighbors(1,ibox)
        nsepa=0
c 
        do 3000 i=1,nneighb
  
        j=neighbors(i+1,ibox)
c 
        call r2d_one_neighb(ibox,j,w,
     2      large_sepa(1+nsepa),nsepa0,w2(iiarr) )
c 
        nsepa=nsepa+nsepa0
 3000 continue
c 
c       remove from the list of nodes in the band all
c       of the said nodes that are also inside one of
c       separated boxes
c 
        call r2d_nodes_subtr(w,large_sepa,nsepa,
     1      iz,inband,nband2,ninband2,w2(iizz),w2(iww) )
c 
        ninband=ninband2
        return
        end
c 
c 
c 
c 
c 
        subroutine r2d_nodes_subtr(w,large_sepa,nsepa,
     1      iz,inband,nbandin,nbandout,izz,ww)
        implicit real *8 (a-h,o-z)
        save
        integer *4 box(11),large_sepa(1),izz(1),iz(1),
     1      inband(1),ww(1)
        dimension corners(2,5),center(3),w(1)
c 
c        construct the list of all nodes in the separated
c        rectangles
c 
        nn=0
        do 1400 i=1,nsepa
c 
        ibox=large_sepa(i)
        call r2dretr(ibox,w,box,corners,center)
c 
        i1=box(6)-1
        n1=box(7)
        do 1200 j=1,n1
c 
        nn=nn+1
        izz(nn)=i1+j
c 
 1200 continue
c 
 1400 continue
c 
c        sort both the list of all nodes in the band, and
c        the list of nodes in the band that are inside
c        separated boxes
c 
        do 1600 i=1,nbandin
c 
        izz(nn+i)=inband(i)
 1600 continue
c 
        nnn=nbandin+nn
c 
        call sortanyn(izz,nnn,ww)
c 
c       eliminate the duplicate nodes from array izz
c 
        j=0
        i1=1
        do 1800 i=1,nnn
c 
        if(i1 .eq. nnn) then
            j=j+1
            inband(j)=izz(i1)
            goto 2200
        endif
c 
        if(i1 .gt. nnn) goto 2200
c 
        if(izz(i1+1) .ne. izz(i1)) then
            j=j+1
            inband(j)=izz(i1)
            i1=i1+1
            goto 1800
        endif
c 
c 
        if(izz(i1+1) .eq. izz(i1)) then
            i1=i1+2
            goto 1800
        endif
c 
 1800 continue
c 
 2200 continue
c 
        nbandout=j
        return
        end
c 
c 
c 
c 
c 
        subroutine r2d_one_neighb(ibox,ineighb,w,
     1      large_sepa,nsepa,iarr)
        implicit real *8 (a-h,o-z)
        save
        integer *4 iarr(2,1),box(12),boxj(12),large_sepa(1)
        dimension w(1),center(3),
     1      cornersj(2,5),centerj(3),corners(2,5)
c 
c        For the user-supplied box number ibox, this subroutine
c        finds all descendants of the neighbors of ibox that
c        are separeted from ibox and whose parents are not
c        separated from ibox
c 
c        . . . initialize the refinement process
c 
        call r2dretr(ibox,w,box,corners,center)
c 
        iarr(1,1)=ineighb
        iarr(2,1)=0
c 
c        Scan the array iarr. Whenever we find a box that
c        has not been processed, process it
c 
        narr=1
c 
c       ... The element iarr(1,i) of the array iarr is the sequence
c           number of the box. The interpretation of the element
c           iarr(2,i) is as follows.
c 
c       iarr(2,i)=0 means that the box number iarr(1,i) has not been
c                   examined
c       iarr(2,i)=7 means that the box number iarr(1,i) is separated
c                   from the box number ibox
c       iarr(2,i)=-7 means that the box number iarr(1,i) is childless
c                   and touches the box number ibox
c       iarr(2,i)=11 means that the box number iarr(1,i) has kids
c                   and touches the box number ibox
c 
        do 3200 ijk=1,100
c 
        ifdid=0
        narr2=narr
c 
        do 3000 i=1,narr2
c 
        if(iarr(2,i) .ne. 0) goto 3000
c 
        ifdid=1
c 
c        if this box does not touch the box number ibox
c        - mark it accordingly
c 
        jbox=iarr(1,i)
        call r2dretr(jbox,w,boxj,cornersj,centerj)
        call r2d_ifinters(cornersj,corners,ifclose)
c 
        if(ifclose .eq. 1) goto 1800
c 
        iarr(2,i)=7
        goto 3000
c 
 1800 continue
c 
c        This box touches the box number ibox. If it has kids,
c        store their numbers in array iarr, and mark it
c        accordingly. If it is childless, mark it as such
c 
        if(boxj(4) .ne. 0) goto 2400
c 
        iarr(2,i)=-7
        goto 3000
c 
 2400 continue
c 
c      This box has kids. Enter the kids in the array iarr
c 
        iarr(2,i)=11
c 
        narr=narr+1
        iarr(1,narr)=boxj(4)
        iarr(2,narr)=0
c 
        if(boxj(5) .eq. 0) goto 3000
c 
        narr=narr+1
        iarr(1,narr)=boxj(5)
        iarr(2,narr)=0
c 
 3000 continue
c 
        if(ifdid .eq. 0) goto 3400
c 
 3200 continue
c 
 3400 continue
c 
c        scan array iarr, storing in array large_sepa
c        all childless boxes separated from the box ibox
c 
        nsepa=0
c 
        do 3600 i=1,narr
c 
        if( iarr(2,i) .ne. 7 )  goto 3600
c 
c        this box is separated from box ibox and has no kids.
c        Record its name
c 
        nsepa=nsepa+1
        large_sepa(nsepa)=iarr(1,i)
c 
 3600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine r2d_ifinters(c1,c2,ifinter)
        implicit real *8 (a-h,o-z)
        save
        dimension c1(2,1),c2(2,1)
c 
c        this subroutine determines if two boxes in the square
c        intersect or touch.
c 
c                input parameters:
c 
c  c1 - the four corners of the first box
c  c2 - the four corners of the second box
c 
c                output parametes:
c 
c  ifinter - the indicator.
c      ifinter=1 means that the boxes intersect
c      ifinter=0 means that the boxes do not intersect
c 
c       . . . find the maximum and minimum coordinates
c             for both boxes
c 
        xmin1=dmin1(c1(1,1),c1(1,2),c1(1,3),c1(1,4) )
        xmin2=dmin1(c2(1,1),c2(1,2),c2(1,3),c2(1,4) )
c 
        xmax1=dmax1(c1(1,1),c1(1,2),c1(1,3),c1(1,4) )
        xmax2=dmax1(c2(1,1),c2(1,2),c2(1,3),c2(1,4) )
c 
        ymin1=dmin1(c1(2,1),c1(2,2),c1(2,3),c1(2,4) )
        ymin2=dmin1(c2(2,1),c2(2,2),c2(2,3),c2(2,4) )
c 
        ymax1=dmax1(c1(2,1),c1(2,2),c1(2,3),c1(2,4) )
        ymax2=dmax1(c2(2,1),c2(2,2),c2(2,3),c2(2,4) )
c 
c        decide if the boxes intersect
c 
        eps=xmax1-xmin1
        if(eps .gt. xmax2-xmin2) eps=xmax2-xmin2
        if(eps .gt. ymax2-ymin2) eps=ymax2-ymin2
        if(eps .gt. ymax1-ymin1) eps=ymax1-ymin1
c 
        eps=eps/100
cccc         eps=1.0d-5
c 
        ifinter=1
        if(xmin1 .gt. xmax2+eps) ifinter=0
        if(xmin2 .gt. xmax1+eps) ifinter=0
c 
        if(ymin1 .gt. ymax2+eps) ifinter=0
        if(ymin2 .gt. ymax1+eps) ifinter=0
        return
        end
c 
c 
c 
c 
c 
        subroutine r2d_sift(corners,
     1      iz,z,izs,nzs,izout,nout)
        implicit real *8 (a-h,o-z)
        save
        integer *4 iz(1),izout(1),izs(1)
        dimension z(2,1),corners(2,4)
c 
c        Sift through nc elements in array iz, indexed by
c        the elements of the array izs. Store in array izout
c        those whose corresponding elements in array z
c        are outside the rectangle with the corners corners
c 
        xright=corners(1,1)
        xleft=corners(1,2)
        ybottom=corners(2,3)
        ytop=corners(2,1)
c 
        nout=0
        do 1200 i=1,nzs
c 
        j=iz(izs(i))
        x=z(1,j)
        y=z(2,j)
c 
        if( (xleft .lt. x) .and. (xright .ge. x) .and.
     1      (ybottom .le. y) .and. (ytop .ge. y) ) goto 1200
c 
        nout=nout+1
        izout(nout)=izs(i)
 1200 continue
c 
         return
         end
c 
c 
c 
c 
c 
        subroutine r2d_all_neighbors(w,nall,nlev,neighbors,w2)
        implicit real *8 (a-h,o-z)
        save
        integer *4 neighbors(11,1),box(12),w2(1)
        dimension w(1),corners(2,5),center(3)
c 
        do 1400 i=1,nall
        do 1200 j=1,11
c 
        neighbors(j,i)=0
 1200 continue
 1400 continue
c 
        do 3000 lev=1,nlev
c 
        do 2800 i=1,nall
c 
c       if this box is on a wrong level - bypass it
c 
        call r2dretr(i,w,box,corners,center)
c 
        if(box(9) .eq. -7) goto 2800
c 
        if(box(2) .ne. lev) goto 2800
c 
c        find the neighbors of this box
c 
        ilist=1
        llist=1000
c 
        ilist2=ilist+llist
        llist2=1000
c 
        ilist3=ilist2+llist2
        llist3=1000
c 
        call r2d_neighbors_find(i,w,neighbors,
     1      w2(ilist),w2(ilist2),w2(ilist3) )
 2800 continue
 3000 continue
c 
c       sort the neighbors for each of the boxes
c 
        do 3200 i=1,nall
c 
        nj=neighbors(1,i)
c 
        if(nj .lt. 2) goto 3200
        call r2d_bubble77(neighbors(2,i),nj)
 3200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine r2d_bubble77(ia,n)
        implicit real *8 (a-h,o-z)
        save
        integer *4 ia(n)
  
        do 1400 i=1,n-1
        do 1200 j=1,n-1
c 
        if(ia(j) .le. ia(j+1)) goto 1200
c 
        k=ia(j)
        ia(j)=ia(j+1)
        ia(j+1)=k
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
        subroutine r2d_neighbors_find(ibox,w,neighbors,
     1      list,list2,list3)
        implicit real *8 (a-h,o-z)
        save
        integer *4 box(11),list(1),list2(1),
     1       neighbors(11,1),box3(11),list3(1)
        dimension center(3),corners(2,5),w(1),centerbox(3),
     1      corners3(2,5)
c 
        nn=0
        do 1200 i=1,11
c 
        neighbors(i,ibox)=0
 1200 continue
c 
c        . . . retrieve from the array w the the various types
c              of information pertaining to the box number ibox
c 
        call r2dretr(ibox,w,box,corners,centerbox)
c 
c        if this is a bastard box (as opposed to a legitimate one),
c        call the appropriate subroutine and exit
c 
        if(box(10) .gt. 2) then
            ilist=1
            llist=1000
c 
            ilist2=ilist+llist
            llist2=1000
c 
            ilist3=ilist2+llist2
            llist3=1000
c 
            call r2d_neighbors_find2(ibox,w,neighbors,
     1      list,list2,list3)
            return
        endif
c 
 1300 continue
c 
        if(box(9) .eq. -7) return
c 
c       This box actually exists (it has not been deleted);
c       find this box's daddy and daddy's neighbors; construct
c       the list of their kids; the resulting list will contain
c       the box ibox, its brother and its cousins
c 
c 
        idad=box(3)
c 
        ii=0
c 
        nnn=neighbors(1,idad)
        if(nnn .eq. 0) goto 1600
c 
        do 1400 i=1,nnn
c 
        ii=ii+1
        list(ii)=neighbors(i+1,idad)
 1400 continue
 1600 continue
c 
        call r2dretr(idad,w,box3,corners3,center)
c 
        i0=0
c 
        if( (box3(4) .ne. ibox) .and. (box3(4) .ne. 0) ) then
            i0=1
            list2(1)=box3(4)
        endif
c 
        if( (box3(5) .ne. ibox) .and. (box3(5) .ne. 0) ) then
            i0=1
            list2(1)=box3(5)
        endif
c 
        if(ii .eq. 0) goto 1800
        do 1700 i=1,ii
c 
        list2(i+i0)=list(i)
 1700 continue
c 
 1800 continue
c 
        iii=ii+i0
c 
c       now, look at the list list2. Find in it the neighbors
c       of ibox, and copy them into the array neighbors
c 
        nn=0
        do 2200 i=1,iii
c 
        j=list2(i)
        call r2dretr(j,w,box3,corners3,center)
c 
        call r2d_ifneighb(corners,corners3,ifclose)
c 
        if(ifclose .eq. 0) goto 2200
c 
        nn=nn+1
        neighbors(nn+1,ibox)=j
 2200 continue
c 
        neighbors(1,ibox)=nn
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine r2d_neighbors_find2(ibox,w,neighbors,
     1      list,list2,list3)
        implicit real *8 (a-h,o-z)
        save
        integer *4 box(11),list(1),list2(1),
     1       neighbors(11,1),box3(11),list3(1)
        dimension center(3),corners(2,5),w(1),centerbox(3),
     1      corners3(2,5)
c 
        nn=0
        do 1200 i=1,11
c 
        neighbors(i,ibox)=0
 1200 continue
c 
c        . . . retrieve from the array w the the various types
c              of information pertaining to the box number ibox
c 
        call r2dretr(ibox,w,box,corners,centerbox)
c 
        if(box(9) .eq. -7) return
c 
c       This box actually exists (it has not been deleted);
c       find this box's daddy and daddy's neighbors; construct
c       the list of their kids; the resulting list will contain
c       the box ibox, its brother and its cousins
c 
        idad=box(3)
c 
        ii=1
        list(ii)=idad
c 
        nnn=neighbors(1,idad)
c 
        if(nnn .eq. 0) goto 1600
c 
        if(nnn .eq. 0) goto 1700
        do 1400 i=1,nnn
c 
        ii=ii+1
        list(ii)=neighbors(i+1,idad)
 1400 continue
 1600 continue
c 
 1700 continue
c 
c       Now, the array list contains the numbers of daddy and
c       daddy's neighbors; construct the list of their children
c 
        iii=0
        do 1800 i=1,ii
c 
        j=list(i)
c 
        call r2dretr(j,w,box3,corners3,center)
c 
        if( (box3(4) .ne. 0) .and. (box3(4) .ne. ibox) ) then
            iii=iii+1
            list2(iii)=box3(4)
        endif
c 
        if( (box3(5) .ne. 0) .and. (box3(5) .ne. ibox) ) then
            iii=iii+1
            list2(iii)=box3(5)
        endif
 1800 continue
c 
c       Now, the array list2 contains the numbers of daddy's
c       children and those of the children of daddy's neighbors.
c       Find the daddy's grandchildren, and the grandchildren
c       of daddy's neighbors
c 
        jj=0
        do 2000 i=1,iii
c 
        j=list2(i)
c 
        call r2dretr(j,w,box3,corners3,center)
c 
        if( (box3(4) .ne. 0) .and. (box3(4) .ne. ibox) ) then
            jj=jj+1
            list3(jj)=box3(4)
        endif
c 
        if( (box3(5) .ne. 0) .and. (box3(5) .ne. ibox) ) then
            jj=jj+1
            list3(jj)=box3(5)
        endif
 2000 continue
c 
c       now, look at the list list3. Find in it the neighbors
c       of ibox, and copy them into the array neighbors
c 
        nn=0
        do 2200 i=1,jj
c 
        j=list3(i)
        call r2dretr(j,w,box3,corners3,center)
c 
        call r2d_ifneighb(corners,corners3,ifclose)
c 
        if(ifclose .eq. 0) goto 2200
        nn=nn+1
        neighbors(nn+1,ibox)=j
 2200 continue
c 
        neighbors(1,ibox)=nn
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine r2d_ifneighb(corn1,corn2,ifclose)
        implicit real *8 (a-h,o-z)
        save
        dimension corn1(2,1),corn2(2,1)
c 
        dmin=1.0d20
        dmax=0
c 
        do 1400 i=1,4
        do 1200 j=1,4
c 
        d=(corn1(1,i)-corn2(1,j))**2+
     1      (corn1(2,i)-corn2(2,j))**2
c 
        if(d .lt. dmin) dmin=d
        if(d .gt. dmax) dmax=d
c 
 1200 continue
 1400 continue
c 
        rat=sqrt(dmin/dmax)
        ifclose=0
        if(rat .lt. 1.0d-10) ifclose=1
c 
        return
        end
c 
c 
c 
c 
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
c 
