        implicit real *8 (a-h,o-z)
        integer *4 wint(10000),
     1      iz(10000),laddr(100)
        real *8 z(2,10000),center0(2),x(10000),y(10000),
     1      w(200000)
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
         PRINT *, 'ENTER nbox'
         READ *,nbox
         CALL PRINF('nbox=*',nbox,1 )
c 
         PRINT *, 'ENTER nbtest'
         READ *,nbtest
         CALL PRINF('nbtest=*',nbtest,1 )
c 
c       construct and plot the test points
c 
        call RANLST(NSMALL,Z,W,N)
          call prinf('after ranlst, n=*',n,1)
        a=10
        b=1
c 
cccc        call creelips(z,n,a,b)
ccc          n=n-1
         call prinf('after ranlst, n=*',n,1)
c 
        do 1200 i=1,n
        x(i)=z(1,i)
        y(i)=z(2,i)
 1200 continue
        iw=21
        call RSPLOT(X,Y,N,IW,'random points as created*')
c 
c       now, construct the lists
c 
        lw=200000
        call d2mstrcr(ier,z,n,nbox,
     1    nboxes,iz,laddr,nlev,center0,size,w,lw,lused)
         call prinf('after d2mstrcr, ier=*',ier,1)
         call prin2('after d2mstrcr, center0=*',center0,2)
         call prinf('after d2mstrcr, lused=*',lused,1)
c 
c       for a chosen box, display the lists 2 and 5
c 
         iw=33
         ibox=nbtest
c 
        call allplt12(iw,z,n,nboxes,
     1    center0,size,
     2      'boxes 1, 2, 3, 4 for the marked box*',w,ibox)
c 
c       check the numbers of lists created and the amount of storage used
c 
         if(2 .ne. 3) goto 2400
c 
c       now, display all lists for all boxes
c 
        iw0=40
        do 2200 i=2,nboxes
        iw=iw0+i
        call allplt12(iw,z,n,nboxes,
     1    center0,size,
     2      'boxes 1, 2, 3, 4 for the marked box*',w,i)
c 
 2200 continue
 2400 continue
c 
  
        call d2mlinfo(w,lused,wint)
        call prinf('finally, lused=*',lused,1)
        call prinf('and lists are this big:*',wint,5)
  
  
        stop
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
c 
        subroutine allplt12(iw,z,n,nboxes,
     1    center0,size,mes,w,ibox)
        implicit real *8 (a-h,o-z)
        save
        integer *4 box(11),list(1000),w(1)
        real *8 z(2,1),xx(5),yy(5),center0(2),center(2),
     1      corners(2,5)
        character *1 mes(1)
cccc         call prinf('entered allplt12, iw=*',iw,1)
c 
c        plot the title
c 
        IF(IW.EQ.0) RETURN
        call titlpr(MES,Iw)
C 
C 
C . . . FIND THE LIMITS TO BE SET
C 
        XMIN=1.0D20
        YMIN=1.0D20
        XMAX=-1.0D20
        YMAX=-1.0D20
        DO 1000 I=1,N
        IF(XMIN.GT.Z(1,I))XMIN=Z(1,I)
        IF(XMAX.LT.Z(1,I))XMAX=Z(1,I)
        IF(YMIN.GT.Z(2,I))YMIN=Z(2,I)
        IF(YMAX.LT.Z(2,I))YMAX=Z(2,I)
 1000 CONTINUE
        DX=(XMAX-XMIN)
        DY=YMAX-YMIN
        DD=DX
        IF(DY.GT.DX) DD=DY
        DMIN=XMIN
        IF(YMIN.LT.XMIN) DMIN=YMIN
        DMAX=XMAX
        IF(YMAX.GT.XMAX) DMAX=YMAX
        DMIN=DMIN-DD/10
        DMAX=DMAX+DD/10
 1200 FORMAT(1X,37H#SET DEVICE postscript FILE 'plot.ps',/,
     2  1X,'SET AXIS TOP OFF',/,
     3  1X,'SET AXIS RIGHT OFF',/,
     4  1X,'SET LIMITS X ', 2(1X,E11.5),' Y ',2(1X,E11.5),/,
     5  1X,'SET SCALE EQUAL')
c 
          WRITE(IW,1200) DMIN,DMAX,DMIN,DMAX
 1400 FORMAT(2(2X,E15.7))
ccccc         if(2 .ne. 3) goto 1510
c 
c       plot all particles in the simulation
c 
       WRITE(IW,1400) (Z(1,I),Z(2,I),I=1,N)
 1500 FORMAT(1X,'PLOT dot')
       write(iw,1500)
 1510 continue
c 
c       plot all boxes
c 
        do 3000 i=1,nboxes
c 
c       get the information for this box
c 
        call d2mgetb(ier,i,box,center,corners,w)
  
        level=box(1)
        ii=box(2)
        jj=box(3)
        side=size/2**level
c 
cccc        call d2mcentf(center0,size,level,ii,jj,center)
        xx(1)=center(1)+side/2
        xx(2)=xx(1)-side
        xx(3)=xx(2)
        xx(4)=xx(3)+side
        xx(5)=xx(1)
c 
        yy(1)=center(2)+side/2
        yy(2)=yy(1)
        yy(3)=yy(2)-side
        yy(4)=yy(3)
        yy(5)=yy(1)
 1410 FORMAT(2(2X,E15.7))
        write(iw,1410) (xx(j),yy(j),j=1,5)
        write(iw,1450)
 3000 continue
C1600 FORMAT(1X,'FRAME')
C1800 FORMAT(1X,'EXIT')
 1450 FORMAT(1X,'JOIN')
c1500 FORMAT(1X,'PLOT')
 1600 FORMAT(1X,'#FRAME')
 1800 FORMAT(1X,'#EXIT')
cccc        WRITE(IW,1500)
c 
c        mark the center of the box ibox
c 
        call d2mgetb(ier,ibox,box,center,corners,w)
        write(iw,1410) center(1),center(2)
 3200 format(1x,'plot ''*'' ')
        write(iw,3200)
c 
c       now, mark list 2 of the box ibox
c 
        itype=2
        call d2mgetl(ier,ibox,itype,list,nlist,w)
cccc         call prinf('in allplt12, after d2mgetl, ier=*',ier,1)
cccc         call prinf('in allplt12, list2 is*',list,nlist)
        do 3600 i=1,nlist
        jbox=list(i)
        call d2mgetb(ier,jbox,box,center,corners,w)
 3400 format(1x,' text  ','x ',e11.5, ' y ',e11.5,,
     1    ' center', ' '' ', i1, ' '' ')
cccc         write(iw,3400) centers(1,jbox), centers(2,jbox),
cccc     1       itype
         write(iw,3400) center(1),center(2),itype
 3600 continue
c 
c        plot list 1
c 
  
        itype=1
        call d2mgetl(ier,ibox,itype,list,nlist,w)
cccc         call prinf('in allplt12, after d2mgetl, ier=*',ier,1)
cccc         call prinf('in allplt12, list1 is*',list,nlist)
        do 4200 i=1,nlist
        jbox=list(i)
        call d2mgetb(ier,jbox,box,center,corners,w)
cccc         write(iw,3400) centers(1,jbox), centers(2,jbox),
cccc     1       itype
         write(iw,3400) center(1),center(2),itype
 4200 continue
c 
        itype=3
        call d2mgetl(ier,ibox,itype,list,nlist,w)
cccc        call prinf('in allplt12, after d2mgetl, ier=*',ier,1)
cccc        call prinf('in allplt12, list3 is*',list,nlist)
        do 4400 i=1,nlist
        jbox=list(i)
        call d2mgetb(ier,jbox,box,center,corners,w)
cccc         write(iw,3400) centers(1,jbox), centers(2,jbox),
cccc     1       itype
         write(iw,3400) center(1),center(2),itype
 4400 continue
c 
        itype=4
        call d2mgetl(ier,ibox,itype,list,nlist,w)
cccc         call prinf('in allplt12, after d2mgetl, ier=*',ier,1)
cccc         call prinf('in allplt12, list4 is*',list,nlist)
        do 4600 i=1,nlist
        jbox=list(i)
        call d2mgetb(ier,jbox,box,center,corners,w)
cccc         write(iw,3400) centers(1,jbox), centers(2,jbox),
cccc     1       itype
         write(iw,3400) center(1),center(2),itype
 4600 continue
c 
        WRITE(IW,1600)
        WRITE(IW,1800)
        close(iw)
        return
        end
  
c 
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c         this is the end of the debugging code and the beginning
c         of the actual logic subroutines for the FMM in the plane
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        subroutine d2mstrcr(ier,z,n,nbox,
     1    nboxes,iz,laddr,nlev,center,size,w,lw,lused777)
        implicit real *8 (a-h,o-z)
        save
        integer *4 iz(1),w(1),laddr(1),box(1),nums(1)
        real *8 z(2,1),center(2)
c 
c        this subroutine constructs the logical structure for the
c        fully adaptive FMM in two dimensions and stores it in the
c        array w in the form of a linked list. after that, the user
c        can obtain the information about various boxes and lists
c        in it by calling the entries d2mgetb, d2mgetl, d2mlinfo
c        of this subroutine (see).
c 
c              note on the list conventions.
c 
c    list 1 of the box ibox - the list of all boxes with which the
c           box ibox interacts directly, including the boxes on the
c           same level as ibox, boxes on the finer levels, and boxes
c           on the coarser levels. obviously, list 1 is empty for any
c           box that is not childless.
c 
c    list 2 of the box ibox - the list of all boxes with which the
c           box ibox interacts in the regular multipole fashion, i.e.
c           boxes on the same level as ibox that are separated from it
c           but whose daddies are not separated from the daddy of ibox.
c 
c    list 3 of the box ibox - for a childless ibox, the list of all
c           boxes on the levels finer than that of ibox, which are
c           separated from ibox, and whose daddys are not separated
c           from ibox. for a box with children, list 3 is empty.
c 
c    list 4 is dual to list 3, i.e. jbox is on the list 4 of ibox if
c           and only if ibox is on the list 3 of jbox.
c 
c                            input parameters:
c 
c  z - the user-specified points in the plane
c  n - the number of elements in array z
c  nbox - the maximum number of points in a box on the finest level
c  lw - the amount of memory in the array w (in integer *4 elements)
c 
c                            output parameters:
c 
c  ier - error return code
c    ier=0   means successful execution
c    ier=32  means that the amount lw of space in array w
c                 is insufficient
c    ier=16 means that the subroutine attempted to construct more
c        than 199 levels of subdivision; indicates bad trouble.
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
c  laddr - an integer array dimensioned (2,nlev), describing the
c         numbers of boxes on various levels of sybdivision, so that
c         the first box on level (i-1) has sequence number laddr(1,i),
c         and there are laddr(2,i) boxes on level i-1
c  nlev - the maximum level number on which any boxes have
c         been created. the maximim number possible is 200.
c         it is recommended that the array laddr above be
c         dimensioned at least (2,200), in case the user underestimates
c         the number of levels required.
c  center - the center of the box on the level 0, containing
c         the whole simulation
c  size - the side of the box on the level 0
c  w - the array containing all tables describing boxes, lists, etc.
c         it is a linked list (for the most part), and can only be accessed
c         via the entries d2mgetb, d2mgetl, d2mlinfo, of this  subroutine
c         (see below). the first lused 777 integer *4 locations of
c         this array should not be altered between the call to this
c         entry and subsequent calls to the entries d2mgetb, d2mgetl,
c         d2mlinfo, of this  subroutine
c 
c  lused777 - the amount of space in the array w (in integer *4 words)
c        that is occupied by various tables on exit from this
c        subroutine. this space should not be altered between the
c        call to this entry and subsequent calls to entries d2mgetb,
c        d2mgetl, d2mlinfo, of this  subroutine (see below).
c 
c        . . . construct the quad-tree structure for the user-specified
c              set of points
c 
        nboxes7=nboxes
        ninire=2
        iiwork=1
        liwork=n+4
c 
        iboxes=iiwork+liwork
        lboxes=lw-n-5
        maxboxes=lboxes/10-1
c 
        call d2mallb(ier,z,n,nbox,w(iboxes),maxboxes,
     1    nboxes,iz,laddr,nlev,center,size,w(iiwork) )
c 
c        if the memory is insufficient - bomb
c 
        if(ier .eq. 0) goto 1100
           call prinf('in d2mstrcr after d2mallb, ier=*',ier,1)
        if(ier .eq. 4) ier=32
        return
 1100 continue
c 
c       compress the array w
c 
cccc          if(2 .ne. 3) goto 1300
        nn=nboxes*10
        do 1200 i=1,nn
        w(i)=w(iboxes+i-1)
 1200 continue
        iboxes=1
 1300 continue
        lboxes=nboxes*10+11
c 
c       construct the centers and the corners for all boxes
c       in the quad-tree
c 
        icenters=iboxes+lboxes
        lcenters=(nboxes*2+2)*ninire
c 
        icorners=icenters+lcenters
        lcorners=(nboxes*8+2)*ninire
c 
        iwlists=icorners+lcorners
        lwlists=lw-iwlists-6
c 
        call centcorn(center,size,w(iboxes),nboxes,
     1      w(icenters),w(icorners) )
c 
c       now, construct all lists for all boxes
c 
  
           call prinf('before crelists, lwlists=*',lwlists,1)
        call crelists(ier,w(iboxes),nboxes,w(icorners),
     1        w(iwlists),lw,lused)
c 
        lused777=lused+iwlists
        call prinf('after crelists, ier=*',ier,1)
c 
        return
c 
c 
c 
c 
         entry d2mgetb(ier,ibox,box,center,corners,w)
c 
c        this entry returns to the user the characteristics of
c        user-specified box  ibox.
c 
c                     input parameters:
c 
c  ibox - the box number for which the information is desired
c  w - storage area as created by the entry d2mstrcr (see above)
c 
c                     output parameters:
c 
c  ier - the error return code.
c    ier=0 means successful execution
c    ier=4 means that ibox is either greater than the number of boxes
c           in the structure or less than 1.
c  boxes - an integer array of length 10 containing the description
c           of the box number ibox in the quad-tree. following is a
c           description of the entries in array box.
c 
c       box(1) - the level of subdivision on which this box
c             was constructed;
c       box(2), box(3) - the coordinates of this box among  all
c             boxes on this level
c       box(4) - the daddy of this box, identified by it address
c             in array boxes
c       box(5), box(6), box(7), box(8) - the list of children of
c              this box (four of them, and the child is identified by
c              its address in the array boxes; if a box has only one
c              child, only the first of the four child entries is
c              non-zero, etc.)
c 
c       box(9) - the location in the array iz of the particles
c             living in this box
c       box(10) - the number of particles living in this box
c       center - the center of the box number ibox
c       corners - the corners of the box number ibox
c 
c       . . . return to the user all information about the box ibox
c 
        ier=0
ccc        if( (ibox .ge. 1)  .and. (ibox .le. nboxes7) ) goto 2100
ccc        ier=4
ccc        return
ccc 2100 continue
c 
        ibox0=iboxes+(ibox-1)*10-1
        do 2200 i=1,10
        box(i)=w(ibox0+i)
 2200 continue
c 
c      return to the user the center and the corners of the box ibox
c 
        call d2mcpcc(w(icenters),w(icorners),ibox,center,corners)
c 
        return
c 
c 
c 
c 
         entry d2mgetl(ier,ibox,itype,list,nlist,w)
c 
c  ibox - the box number for which the information is desired
c  itype - the type of the desired list for the box ibox
c  w - storage area as created by the entry d2mstrcr (see above)
c 
c                     output parameters:
c 
c  ier - the error return code.
c    ier=0 means successful execution
c    ier=4 means that the list  itype  for the box  ibox  is empty
c  list - the list  itype  for the box  ibox
c  nlist - the number of elements in array  list
c 
c       return to the user the list number itype for the box ibox
c 
        call linkretr(ier,itype,ibox,list,nlist,w(iwlists),lused)
        return
c 
c 
c 
c 
        entry d2mlinfo(w,lused77,nums)
c 
c       this entry returns to the user some of the information
c       about the storage area. namely, it returns the total
c       amount lused777 of memory utilized in the array w (in integer *4
c       locations), and the integer array nums containing the numbers
c       of elements in each of the lists
  
c 
        call linkinfo(w(iwlists),lused77,nums)
        return
        end
c 
c 
c 
c 
c 
        subroutine d2mcpcc(centers,corners,ibox,center,corner)
        implicit real *8 (a-h,o-z)
        save
        dimension centers(2,1),corners(2,4,1),center(2),corner(2,4)
c 
        center(1)=centers(1,ibox)
        center(2)=centers(2,ibox)
c 
        corner(1,1)=corners(1,1,ibox)
        corner(2,1)=corners(2,1,ibox)
c 
        corner(1,2)=corners(1,2,ibox)
        corner(2,2)=corners(2,2,ibox)
c 
        corner(1,3)=corners(1,3,ibox)
        corner(2,3)=corners(2,3,ibox)
c 
        corner(1,4)=corners(1,4,ibox)
        corner(2,4)=corners(2,4,ibox)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine crelists(ier,boxes,nboxes,corners,w,lw,lused)
        implicit real *8 (a-h,o-z)
        save
        integer *4 boxes(10,1),collkids(50),w(1),
     1      dadcolls(20),list5(500),stack(600)
        real *8 corners(2,4,1)
c 
c        this subroutine constructs all lists for all boxes
c        and stores them in the storage area w in the form
c        of a linked list. the resulting data can be accessed
c        by calls to various entries of the subroutine linkinit (see).
c 
c                          input parameters:
c 
c  boxes - an integer array dimensioned (10,nboxes), as created by
c        the subroutine d2mallb (see).  each 10-element column
c         describes one box, as follows:
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
c  nboxes - the total number of boxes created
c  corners - the array of corners of all the boxes to in array boxes
c  lw - the total amount of storage in array w (in integer *4 words)
c 
c              output parameters:
c 
c  ier - the error return code.
c    ier=0 means successful execution
c    ier=32 means that the amount lw of space in array w is
c           insufficient. it is a fatal error.
c  w - storage area containing all lists for all boxes in
c        the form of linked lists, accessible by the subroutine
c        linkretr (see).
c  lused - the amount of space in the array w (in integer *4 words)
c        that is occupied by various tables on exit from this
c        subroutine. this space should not be altered between the
c        call to this subroutine and subsequent calls to the
c        entries linkretr, etc. of the subroutine linkinit (see).
c 
c       . . . initialize the storage-retrieval routine for all
c             boxes
c 
        ier=0
        ntypes=5
        call linkinit(ier,nboxes,ntypes,w,lw)
         call prinf('in crelists after linkinit, ier=*',ier,1)
c 
c        construct lists 5, 2 for all boxes
c 
        do 2000 ibox=2,nboxes
c 
c       find this guy's daddy
c 
        idad=boxes(4,ibox)
c 
c       find daddy's collegues, including daddy himself
c 
        dadcolls(1)=idad
        itype5=5
        itype2=2
        call linkretr(ier,itype5,idad,dadcolls(2),ncolls,w,lused)
        ncolls=ncolls+1
c 
c        find the children of the daddy's collegues
c 
        nkids=0
        do 1600 i=1,ncolls
        icoll=dadcolls(i)
        do 1400 j=1,4
        kid=boxes(4+j,icoll)
        if(kid .le. 0) goto 1600
        if(kid .eq. ibox) goto 1400
        nkids=nkids+1
        collkids(nkids)=kid
 1400 continue
 1600 continue
c 
c       sort the kids of the daddy's collegues into the
c       lists 2, 5 of the box ibox
c 
        nlist1=1
        do 1800 i=1,nkids
c 
c       check if this kid is touching the box ibox
c 
        kid=collkids(i)
        call ifinters(corners(1,1,kid),corners(1,1,ibox),ifinter)
c 
        if(ifinter .eq. 1)
     1    call linkstor(ier,itype5,ibox,kid,nlist1,w,lused)
c 
        if(ifinter .eq. 0)
     1    call linkstor(ier,itype2,ibox,kid,nlist1,w,lused)
 1800 continue
c 
c        if storage capacity has been exceeed - bomb
c 
        if(ier .eq. 32) return
c 
 2000 continue
c 
c       now, construct lists 1, 3
c 
        do 3000 i=1,nboxes
c 
c       if this box has kids - its lists 1, 3 are empty;
c       do not construct them
c 
        if(boxes(5,i) .gt. 0) goto 3000
  
  
        call linkretr(jer,itype5,i,list5,nlist,w,lused)
  
        if(jer .eq. 4) goto 3000
c 
        do 2200 j=1,nlist
        jbox=list5(j)
        call lists31(ier,i,jbox,boxes,nboxes,
     1    corners,w,stack,lused)
c 
c        if storage capacity has been exceeed - bomb
c 
        if(ier .eq. 32) return
 2200 continue
 3000 continue
c 
c       copy all elements of lists 1, 2, 3, while skipping list 4
c 
        ntypes4=4
        call linkinit(ier,nboxes,ntypes4,w(lused+1),lw-(lused+5))
        do 3600 ibox=1,nboxes
        do 2400 itype=1,3
c 
        call linkretr(jer,itype,ibox,list5,nlist,w,lused)
        if(jer .eq. 4) goto 2400
        call linkstor(jer,itype,ibox,list5,nlist,w(lused+1),lused2)
 2400 continue
 3600 continue
c 
c       compress array w
c 
        do 4000 i=1,lused2
        w(i)=w(lused+i)
 4000 continue
        lused=lused2
c 
c        finally, construct the lists 4 for all boxes
c        that need them
c 
        itype3=3
        itype4=4
        nlist1=1
        do 4400 ibox=1,nboxes
c 
        call linkretr(jer,itype3,ibox,list5,nlist,w,lused)
        if(jer .eq. 4) goto 4400
        do 4200 j=1,nlist
        call linkstor(ier,itype4,list5(j),ibox,nlist1,w,lused2)
c 
c        if storage capacity has been exceeed - bomb
c 
        if(ier .eq. 32) return
 4200 continue
 4400 continue
        lused=lused2
  
         call prinf('exiting crelists, lused=*',lused,1)
        return
        end
c 
c 
c 
c 
c 
        subroutine lists31(ier,ibox,jbox0,boxes,nboxes,
     1    corners,w,stack,lused)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),corners(2,4,1)
        integer *4 boxes(10,1),stack(3,200)
        data itype1/1/,itype2/2/,itype3/3/,itype4/4/,itype5/5/,
     1      nlist1/1/
c 
c       this subroutine constructs all elements of lists 1 and 3
c       resulting from the subdivision of one element of list 5
c       of the box ibox. all these elements of lists 1, 3 are
c       stored in the linked lists by the subroutine linstro (see)
c 
c        input parameters:
c 
c  ibox - the box whose lists are being constructed
c 
c  jbox0 - the element of list 5 of the box ibox that is being
c          subdfivided
c  boxes - the array boxes as created by the subroutine d2mallb (see)
c  nboxes - the number of boxes in array boxes
c  corners - the array of corners of all the boxes to in array boxes
c  w - the storage area formatted by the subroutine linkinit (see)
c          to be used to store the elements of lists 1, 3 constructed
c          by this subroutine. obviously, by this time, it contains
c          plenty of other lists.
c 
c                      output parameters:
c 
c  ier - the error return code.
c    ier=0 means successful execution
c    ier=32 means that the amount lw of space in array w is
c           insufficient. it is a fatal error.
c  w - the augmented storage area, containing all the boxes just
c          created, in addition to whatever had been stored previously
c  lused - the total length of array w (in integer *4 words) used
c          on exit from this subroutine
c 
c                      work arrays:
c 
c  stack - must be at least 600 integer *4 locations long
c 
c       . . . starting with the initial element of list 5,
c             subdivide the boxes recursively and store
c             the pieces where they belong
c 
c        . . . initialize the process
c 
        jbox=jbox0
        istack=1
        stack(1,1)=1
        stack(2,1)=jbox
        nsons=0
  
        if(boxes(5,jbox) .gt. 0) nsons=nsons+1
        if(boxes(6,jbox) .gt. 0) nsons=nsons+1
        if(boxes(7,jbox) .gt. 0) nsons=nsons+1
        if(boxes(8,jbox) .gt. 0) nsons=nsons+1
        stack(3,1)=nsons
c 
c       . . . move up and down the stack, generating the elements
c             of lists 1, 3 for the box jbox, as appropriate
c 
        do 5000 ijk=1,10000000
c 
c       if this box is separated from ibox - store it in list 3;
c       enter this fact in the daddy's table; pass control
c       to the daddy
c 
        call ifinters(corners(1,1,ibox),corners(1,1,jbox),ifinter)
c 
        if(ifinter .eq. 1) goto 2000
        call linkstor(ier,itype3,ibox,jbox,nlist1,w,lused)
c 
c        if storage capacity has been exceeed - bomb
c 
        if(ier .eq. 32) return
        istack=istack-1
        stack(3,istack)=stack(3,istack)-1
        jbox=stack(2,istack)
        goto 5000
 2000 continue
c 
c       this box is not separated from ibox. if it is childless
c       - enter it in list 1; enter this fact in the daddy's table;
c       pass control to the daddy
c 
        if(boxes(5,jbox) .ne. 0) goto 3000
        call linkstor(ier,itype1,ibox,jbox,nlist1,w,lused)
c 
c        if storage capacity has been exceeed - bomb
c 
        if(ier .eq. 32) return
c 
c       . . . entered jbox in the list1 of ibox; if jbox
c             is on the finer level than ibox - enter ibox
c             in the list 1 of jbox
c 
        if(boxes(1,jbox) .eq. boxes(1,ibox)) goto 2400
        call linkstor(ier,itype1,jbox,ibox,nlist1,w,lused)
c 
c        if storage capacity has been exceeed - bomb
c 
        if(ier .eq. 32) return
 2400 continue
c 
c       if we have processed the whole box jbox0, get out
c       of the subroutine
c 
        if(jbox .eq. jbox0) return
c 
        istack=istack-1
        stack(3,istack)=stack(3,istack)-1
        jbox=stack(2,istack)
        goto 5000
 3000 continue
c 
c       this box is not separated from ibox, and has children. if
c       the number of unprocessed sons of this box is zero
c       - pass control to his daddy
c 
        nsons=stack(3,istack)
        if(nsons .ne. 0) goto 4000
c 
        if(jbox .eq. jbox0) return
c 
        istack=istack-1
        stack(3,istack)=stack(3,istack)-1
        jbox=stack(2,istack)
        goto 5000
 4000 continue
c 
c       this box is not separated from ibox; it has sons, and
c       not all of them have been processed. construct the stack
c       element for the appropriate son, and pass the control
c       to him.
c 
        jbox=boxes(4+nsons,jbox)
        istack=istack+1
        nsons=0
        if(boxes(5,jbox) .gt. 0) nsons=nsons+1
        if(boxes(6,jbox) .gt. 0) nsons=nsons+1
        if(boxes(7,jbox) .gt. 0) nsons=nsons+1
        if(boxes(8,jbox) .gt. 0) nsons=nsons+1
        stack(1,istack)=istack
        stack(2,istack)=jbox
        stack(3,istack)=nsons
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
        subroutine ifinters(c1,c2,ifinter)
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
         eps=1.0d-5
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
        subroutine linkinit(ier,nboxes,ntypes,w,lw)
        save
        integer *4 w(1),nums(1),inums(20)
        data iilistad/1/,iilists/2/,inumele/3/,inboxes/4/,
     1      intypes/5/,ilw/6/,iltot/7/,
     2      inums/11,12,13,14,15,16,17,18,19,20,21,22,23,24,
     3            25,26,27,28,29,30/
c 
c        this is the initialization entry point for the linked list
c        storage-retrieval facility. it formats the array w, to
c        be used by the entries linkstor, linkretr, linkrem,
c        linkinfo below.
c 
c                     input parameters:
c 
c  nboxes - the number of boxes for which the various types of lists
c        will be stored
c  ntypes - the number of types of lists that will be stored
c  lw - the total amount of space in the area w to be used for storage
c        (in integer *4 locations)
c 
c                     output parameters:
c 
c  ier - error return code;
c    ier=0 means successful execution
c    ier=1024 means that the amount of space in array w is grossly
c        insufficient
c  w - the formatted area for storage
c 
c 
c       . . . allocate memory for the storage facility
c 
        ier=0
c 
        ilistadd=32
        nlistadd=nboxes*ntypes+10
c 
        ilists=ilistadd+nlistadd
        numele=0
        ltot=ilists+numele*2+10
c 
         if(ltot+100 .lt. lw) goto 1200
         ier=1024
         return
 1200 continue
c 
         do 1400 i=1,20
         w(inums(i))=0
 1400 continue
c 
        w(iilistad)=ilistadd
        w(iilists)=ilists
        w(inumele)=numele
        w(inboxes)=nboxes
        w(intypes)=ntypes
        w(ilw)=lw
        w(iltot)=ltot
c 
        call linkini0(w(ilistadd),nboxes,ntypes)
        return
c 
c 
c 
c 
        entry linkstor(ier,itype,ibox,list,nlist,w,lused)
c 
c       this entry stores dynamically a list of positive numbers
c       in the storage array w.
c 
c                      input parameters:
c 
c  itype - the type of the elements being stored
c  ibox - the box to which these elements corresponds
c  list - the list of positive integer elements to be stored
c  nlist - the number of elements in the array list
c  w - the storage area used by this subroutine; must be first
c          formatted by the entry linkinit of this subroutine
c          (see above)
c 
c                      output parameters:
c 
c  ier - error return code;
c         ier=0 means successful execution
c         ier=32 means that the storage area w would be exceeded
c                 by this storage request
c  w - the storage area used by this subroutine
c  lused - the number of integer *4 elements used in array
c          w after this call.
c 
c       . . . if this storage request exceeds the available memory - bomb
c 
        ier=0
        if(w(iilists)+w(inumele)*2+nlist*2 .lt. w(ilw) ) goto 2200
        ier=32
        return
 2200 continue
c 
c       store the user-specified list in array w
c 
        call linksto0(itype,ibox,list,nlist,w(w(iilistad)),
     1      w(inboxes),w(w(iilists)),w(inumele),w(inums(itype)) )
c 
c       augment the amount of storege used
c 
  
        lused=w(iilists)+w(inumele)*2+10
        return
c 
c 
c 
c 
        entry linkretr(ier,itype,ibox,list,nlist,w,lused)
c 
c       this entry retrieves from the storage area  w
c       a list of positive numbers that has been stored there
c       by the entry linkstor (see above).
c 
c                      input parameters:
c 
c  itype - the type of the elements to be retrieved
c  ibox - the box to which these elements correspond
c  w - the storage area from which the information is to be
c          retrieved
c 
c                      output parameters:
c 
c  ier - error return code;
c          ier=0 means successful execution, nothing to report
c          ier=4 means that no elements are present of the type
c                  itype and the box ibox
c  list - the list of positive integer elements retrieved
c  nlist - the number of elements in the array list
c  lused - the number of integer *4 elements used in array
c          w after this call.
c 
        call linkret0(ier,itype,ibox,w(w(iilistad)),
     1      w(w(iilists)),list,w(inboxes),nlist)
c 
        lused=w(iilists)+w(inumele)*2+10
c 
        return
c 
c 
c 
c 
        entry linkrem(ier,itype,ibox,list,nlist,w,lused)
c 
c       this entry deletes elements in array lists corresponding
c       the user-specified list  list. actually, it does not
c       delete anything, but rather marks these elements for
c       future destruction by making them negative.
c 
c                      input parameters:
c 
c  itype - the type of the elements to be destroyed
c  ibox - the box to which these elements correspond
c  list - the list of positive integer elements to be destroyed
c  nlist - the number of elements in the array list
c  w - the storage area from which the information is to be
c          retrieved
c 
c                      output parameters:
c 
c  ier - error return code;
c          ier=0 means successful execution, nothing to report
c          ier=4*k means that k of the elements on the user-specified
c                   list  list were not present
c          ier=22 means that  no elements whatsoever were present
c                   for the type  itype   and the box  ibox
c  w - the storage area from which the information is to be
c          retrieved
c  lused - the number of integer *4 elements used in array
c          w both before and after this call.
c 
c       mark for destruction the user-specified elements
c 
        call linkrem0(ier,itype,ibox,list,nlist,w(w(iilistad)),
     1      w(inboxes),w(w(iilists)),w(inums(itype)) )
c 
        lused=w(iilists)+w(inumele)*2+10
        return
c 
c 
c 
c 
        entry linkinfo(w,lused,nums)
c 
c       this entry returns to the user some of the information
c       about the storage area. namely, it returns the total
c       amount lused of memory utilized in the array w (in integer *4
c       locations), and the integer array nums containing the numbers
c       of elements in each of the lists
c 
        lused=w(iilists)+w(inumele)*2+10
        ntypes7=w(intypes)
         call prinf('in linkinfo, lused=*',lused,1)
         call prinf('in linkinfo, ntypes7=*',ntypes7,1)
         call prinf('in w(inumele)=*',w(inumele),1)
        do 6200 i=1,ntypes7
        nums(i)=w(inums(i))
 6200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine linksto0(itype,ibox,list,nlist,listaddr,
     1      nboxes,lists,numele,numtype)
        save
        integer *4 listaddr(nboxes,1),lists(2,1),list(1)
c 
c       this entry stores dynamically a list of positive numbers
c       in the storage array lists, while entering the information
c       about this event in the array listaddr.
c 
c                      input parameters:
c 
c  itype - the type of the elements being stored
c  ibox - the box to which these elements corresponds
c  list - the list of positive integer elements to be stored
c  nlist - the number of elements in the array list
c  listaddr - the addressing array for the main storage array lists;
c             it is assumed that it has been formatted by a call
c             to the entry linkini0 of this subroutine (see below).
c  nboxes - the total number of boxes indexing the elements
c           in array lists
c  lists - the main storage area used by this subroutine
c  numele - the number of elements stored in array lists on entry
c             to this subroutinec
  
c                      output parameters:
c 
c  listaddr - the addressing array for the main storage array lists;
c             it is assumed that it has been formatted by a call
c             to the entry linkini0 of this subroutine (see below).
c  lists - the main storage area used by this subroutine
c  numele - the number of elements stored in array lists on exit
c             from this subroutine
c  numtype - the total number of elements in all lists of the
c           type itype AFTER this call
c 
c .........................................................................
c 
c        interpretation of the entries in arrays lists, listaddr:
c 
c  lists(1,i) - the location in array lists of the preceding
c        element in the list with (box,type) combination as the
c        user-supplied (ibox,itype)
c     lists(1,i) .leq. 0 means that this is the first element
c        of its type,
c  lists(2,i) - the box number on the interaction list.
c 
c 
c  listaddr(ibox,itype) - the address in the array lists of the last element
c        of this type for this box;
c     listaddr(ibox,itype) .leq. 0 means that there are no elements
c        in the list of this type for this box.
c 
c       . . . store the user-supplied list elements in the array lists,
c             and enter information about this change in array listaddr
c 
        ilast=listaddr(ibox,itype)
        do 1200 i=1,nlist
        numele=numele+1
        numtype=numtype+1
        lists(1,numele)=ilast
        lists(2,numele)=list(i)
c 
        ilast=numele
 1200 continue
        listaddr(ibox,itype)=ilast
c 
        return
c 
c 
c 
c 
        entry linkret0(ier,itype,ibox,listaddr,lists,list,
     1      nboxes,nlist)
c 
c       this entry retrieves from the main storage array lists
c       a list of positive numbers that has been stored there
c       by the entry linksto0 (see above).
c 
c                      input parameters:
c 
c  itype - the type of the elements being to be retrieved
c  ibox - the box to which these element corresponds
c  listaddr - the addressing array for the main storage array lists;
c             it is assumed that it has been formatted by a call
c             to the entry linkini0 of this subroutine (see below).
c  nboxes - the total number of boxes indexing the elements
c           in array lists
c  lists - the main storage area used by this subroutine
c 
c                      output parameters:
c 
c  ier - error return code;
c          ier=0 means successful execution, nothing to report
c          ier=4 means that no elements are present of the type
c                  itype and the box ibox
c  list - the list of positive integer elements to be stored
c  nlist - the number of elements in the array list
c 
c 
c       . . . retrieve and store in array list the list of the
c             type itype for the box ibox
c 
        ier=0
        ilast=listaddr(ibox,itype)
        if(ilast .gt. 0) goto 2200
        nlist=0
        ier=4
        return
 2200 continue
c 
        nlist=0
        do 2400 i=1,10000000
c 
        if(lists(2,ilast) .le. 0) goto 2300
        nlist=nlist+1
        list(nlist)=lists(2,ilast)
 2300 continue
        ilast=lists(1,ilast)
        if(ilast .le. 0) goto 2600
 2400 continue
 2600 continue
c 
        if(nlist .gt. 0) goto 2650
        ier=4
        return
 2650 continue
c 
c        flip the retrieved array
c 
        if(nlist .eq. 1) return
        do 2700 i=1,nlist/2
        j=list(i)
        list(i)=list(nlist-i+1)
        list(nlist-i+1)=j
 2700 continue
        return
c 
c 
c 
c 
        entry linkini0(listaddr,nboxes,ntypes)
c 
c       this subroutine initializes the array listaddr to be used
c       later by other entries of this subroutine
c 
c                        input parameters:
c 
c  nboxes - the number of boxes for which the various types of lists
c        will be stored
c  ntypes - the number of types of lists that will be stored
c 
c                        output parameters:
c 
c  listaddr - the addressing array for the main storage array lists;
c             it is assumed that it has been formatted by a call
c             to the entry linkini0 of this subroutine (see below).
c 
c       . . . initialize the array listaddr
c 
        do 3000 k=1,ntypes
        do 2800 i=1,nboxes
        listaddr(i,k)=-1
 2800 continue
 3000 continue
        return
c 
c 
c 
c 
        entry linkrem0(ier,itype,ibox,list,nlist,listaddr,
     1      nboxes,lists,numtype)
c 
c       this entry deletes elements in array lists corresponding
c       the user-specified list  list. actually, it does not
c       delete anything, but rather marks these elements for
c       future destruction by making them negative.
c 
c                      input parameters:
c 
c  itype - the type of the elements to be destroyed
c  ibox - the box to which these elements correspond
c  list - the list of positive integer elements to be destroyed
c  nlist - the number of elements in the array list
c  listaddr - the addressing array for the main storage array lists;
c             it is assumed that it has been formatted by a call
c             to the entry linkini0 of this subroutine (see below).
c  nboxes - the total number of boxes indexing the elements
c           in array lists
c  lists - the main storage area used by this subroutine
c 
c                      output parameters:
c 
c  ier - error return code;
c          ier=0 means successful execution, nothing to report
c          ier=4*k means that k of the elements on the user-specified
c                   list  list were not present
c          ier=22 means that  no elements whatsoever were present
c                   for the type  itype   and the box  ibox
c  listaddr - the addressing array for the main storage array lists;
c             it is assumed that it has been formatted by a call
c             to the entry linkini0 of this subroutine (see below).
c  lists - the main storage area used by this subroutine
c 
c       . . . mark for destruction the elements of type (itype,ibox)
c             that are in the list  list.
c 
        ier=0
        do 4000 i=1,nlist
        ilast=listaddr(ibox,itype)
        if(ilast. gt. 0) goto 3200
        ier=22
        return
 3200 continue
c 
        iffound=0
        do 3600 j=1,10000000
c 
        if(ilast .le. 0) goto 3800
        if(lists(2,ilast) .ne. list(i)) goto 3400
        lists(2,ilast)=-lists(2,ilast)
        numtype=numtype-1
        iffound=1
 3400 continue
        ilast=lists(1,ilast)
 3600 continue
 3800 continue
         if(iffound .eq. 0) ier=ier+4
c 
 4000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine d2mallb(ier,z,n,nbox,boxes,maxboxes,
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
c         call prin2('in d2mallb, center0=*',center0,2)
         call prin2('in d2mallb, size=*',size,1)
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
ccccc         call prinf('in d2mallb, nbox=*',nbox,1)
         call prinf('in d2mallb, n=*',n,1)
        do 3000 level=0,maxlev
ccccc         call prinf('in d2mallb, level=*',level,1)
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
        call d2mcentf(center0,size,level,ii,jj,center)
c 
        iiz=boxes(9,idad)
        nz=boxes(10,idad)
        call d2msepa1(center,z,iz(iiz),nz,iwork,
     1    is,ns)
c 
c       store in array boxes the sons obtained by the routine
c       d2msepa1
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
        subroutine d2msepa1(cent,z,iz,n,iwork,
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
        call d2msepa0(z,iz,n,itype,thresh,iwork,n12)
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
     1    call d2msepa0(z,iz,n12,itype,thresh,iwork,n1)
        n2=n12-n1
c 
        if(n34 .ne. 0)
     1    call d2msepa0(z,iz(n12+1),n34,itype,thresh,iwork,n3)
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
        subroutine d2msepa0(z,iz,n,itype,thresh,iwork,n1)
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
        subroutine d2mcentf(center0,size,level,i,j,center)
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
        subroutine centcorn(center0,size,boxes,nboxes,
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
c        by the subroutine d2mallb (see)
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
  
  
  
  
  
  
  
  
