        implicit real *8 (a-h,o-z)
        dimension z(30000),
     1      w(30000),x(30000),y(30000),z(2,30000),
     2      iz(30000),cent0(2),
     3      boxes(10,30000),laddr(2,200),
     4      iwork(30000)
c 
        call prini(6,13)
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
c       construct and plot the test points
c 
        call RANLST(NSMALL,Z,W,N)
        a=10
        b=1
c 
        call creelips(z,n,a,b)
         call prinf('after ranlst, n=*',n,1)
c 
        do 1200 i=1,n
        x(i)=z(1,i)
        y(i)=z(2,i)
 1200 continue
        iw=21
        call RSPLOT(X,Y,N,IW,'random points as created*')
c 
c      construct all boxes on all levels
c 
        call prin2('before d2mallb, z=*',z,n*2)
ccccc        call pairtest(z,n)
         maxboxes=2999
cccc         maxboxes=299
        call d2mallb(ier,z,n,nbox,
     1    boxes,maxboxes,nboxes,iz,laddr,nlev,
     2      cent0,size,iwork)
c 
        call prinf('after d2mallb, ier=*',ier,1)
        call prinf('after d2mallb, nlev=*',nlev,1)
        call prinf('after d2mallb, nboxes=*',nboxes,1)
cccc        call prinf('after d2mallb, boxes=*',boxes,nboxes*10)
c 
c        now, plot the whole structure
c 
        iw=25
        call allplt(iw,z,n,boxes,nboxes,
     1    cent0,size,'structure as created*')
  
  
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
C 
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
       N=NSTART-1
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
        subroutine boxprt(boxes,ibox,z,iz,work,
     1      center0,size)
        implicit real *8 (a-h,o-z)
        save
        integer *4 iz(1),boxes(10,1)
        real *8 z(2,1),work(2,1),center(2),center0(2),
     1      xx(4),yy(4)
c 
c       find corners of the box
c 
         call prinf('printing a box on level*',level,1)
        level=boxes(1,ibox)
        ii=boxes(2,ibox)
        jj=boxes(3,ibox)
         call prinf('ii=*',ii,1)
         call prinf('jj=*',jj,1)
        call d2mcentf(center0,size,level,ii,jj,center)
c 
        side=size/2**level
        side2=side/2
        xx(1)=center(1)+side2
        xx(2)=xx(1)-side
        xx(3)=xx(2)
        xx(4)=xx(1)
c 
        yy(1)=center(2)+side2
        yy(2)=yy(1)
        yy(3)=yy(2)-side
        yy(4)=yy(1)
         call prin2('xx are*',xx,4)
         call prin2('yy are*',yy,4)
c 
c        and print the particles in it
c 
        iiz=boxes(9,ibox)
        nz=boxes(10,ibox)
        do 1200 i=1,nz
        j=iz(iiz+i-1)
        work(1,i)=z(1,j)
        work(2,i)=z(2,j)
 1200 continue
        call prin2('in boxprt particles are*',work,nz*2)
        return
        end
c 
c 
c 
c 
c 
        subroutine pairtest(z,n)
        implicit real *8 (a-h,o-z)
        save
        dimension z(2,1)
c 
c       test all pairs of points in z for coincidences
c 
        distmin=1.0d20
        do 1400 i=2,n
        do 1200 j=1,i-1
        d=(z(1,i)-z(1,j))**2+(z(2,i)-z(2,j))**2
        d=dsqrt(d)
        if(d .lt. distmin) distmin=d
 1200 continue
 1400 continue
         call prin2('in pairtest, distmin=*',distmin,1)
        return
        end
c 
c 
c 
c 
c 
        subroutine allplt(iw,z,n,boxes,nboxes,
     1    center0,size,mes)
        implicit real *8 (a-h,o-z)
        save
        integer *4 boxes(10,1)
        real *8 z(2,1),xx(5),yy(5),center0(2),center(2)
        character *1 mes(1)
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
c 
c       plot all particles in the simulation
c 
       WRITE(IW,1400) (Z(1,I),Z(2,I),I=1,N)
 1500 FORMAT(1X,'PLOT')
       write(iw,1500)
c 
c       plot all boxes
c 
        do 3000 i=1,nboxes
        level=boxes(1,i)
        ii=boxes(2,i)
        jj=boxes(3,i)
        side=size/2**level
c 
        call d2mcentf(center0,size,level,ii,jj,center)
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
        WRITE(IW,1600)
        WRITE(IW,1800)
        return
        end
c 
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c  this is the end of the debugging code and the beginning of the
c  actual quad-tree subroutine
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
c  numlev - the maximum level number on which any boxes have
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
c       boxes till none are left more than nbox particles
c 
        maxson=maxboxes
        maxlev=200
cccc        maxlev=4
        ison=1
        nlev=0
        do 3000 level=0,maxlev
         call prinf('in d2mallb, level=*',level,1)
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
        boxes(1,ison)=level+1
        iison=(ii-1)*2+iisons(i)
        jjson=(jj-1)*2+jjsons(i)
        boxes(2,ison)=iison
        boxes(3,ison)=jjson
        boxes(4,ison)=idad
        boxes(5,ison)=0
        boxes(6,ison)=0
        boxes(7,ison)=0
        boxes(8,ison)=0
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
c        so that each of the children accupies a contigious
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
c       or vertically depending on the parameter itype
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
