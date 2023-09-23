        implicit real *8 (a-h,o-z)
        real *8 z(3,10000)
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
        itype=2
        shift=0.1
        slice=0.1
c 
        call sepatest(nsmall,z,itype,shift,slice)
c 
        stop
        end
c 
c 
c 
c 
c 
        subroutine sepatest(nsmall,z,itype,shift,slice)
        implicit real *8 (a-h,o-z)
        save
        dimension z(3,1),w(500 000),cent(3),z2(3,10 000),center0(3)
        integer *4 iz(10 000)
c 
        dimension laddr(1000),centers(10 000),corners(40 000)
        integer *4 boxes(15,10 000)
c 
c       construct the test box to be subdivided
c 
        do 1400 i=1,10 000
        do 1200 j=1,15
        boxes(j,i)=0
c 
 1200 continue
 1400 continue
c 
        size=1
        n=nsmall
        cent(1)=0
        cent(2)=0
        cent(3)=0
c 
        call cornplot(cent,size)
c 
        call RIRANB3(n,SIZE,CENT,Z,W)
  
  
       if(2 .ne. 3) goto 2550
c 
c       set the coordinate number itype to a constant
c 
        call ran2in3(nsmall,z,w,n,itype,shift,z2)
c 
        call prinf('after ran2in3, n=*',n,1)
  
        do 2500 i=1,n
        z(itype,i)=shift
 2500 continue
c 
 2550 continue
        a=1
        b=3
        call creelips(z,n,a,b)
  
c 
        call prin2('z as constructed*',z,100)
c 
c       plot all three projections of the constructed box
c 
        iw=20
c 
        call boxplt1(iw,z,n,'box as created*')
c 
c        test the subroutine separating the whole mess into
c        an oct-tree
c 
        nbox=3
        maxboxes=10000
        lw=500 000
c 
        call d3mstrcr(ier,z,n,nbox,
     1    nboxes,iz,laddr,nlev,center0,size,w,lw,lused777)
c 
        call prinf('after d3mstrcr, ier=*',ier,1)
        call prinf('after d3mstrcr, nboxes=*',nboxes,1)
        call prinf('after d3mstrcr, nlev=*',nlev,1)
c 
c        plot a slice of the structure
c 
        iw=301
        t=slice
        call sliceplt(iw,corners,centers,nboxes,
     1    itype,t,w)
c 
c        test the subroutine constructing all lists for all
c        boxes
c 
        call lststst0(boxes,nboxes,z,iz,n,w)
  
cccc        call lststest(boxes,nboxes,corners,centers,z,iz,n,w)
  
  
  
        return
        end
  
c 
c 
c 
c 
c 
  
        subroutine lststst0(boxes,nboxes,z,iz,n,w)
        implicit real *8 (a-h,o-z)
        save
        dimension z(3,1),f(10 000),f2(10 000),q(10 000),z0(3),w(1),
     1      f3(10 000),rats(10 000),f4(10 000)
        integer *4 ind(10 000),boxes(15,1),iz(1)
c 
c       construct the array ind assigning points to childless boxes
c 
        call wherepts(w,nboxes,iz,ind,n,boxes)
c 
c        construct the charges whose potential is to be evaluated at
c        the nodes
c 
        h=0.7
        call oneinfini(n)
  
        do 1200 i=1,n
        q(i)=dabs(dsin(h*i**2))
 1200 continue
c 
        q(n-1)=1
c 
c       evaluate their potential on the whole ensemble directly
c 
        call poteva(z,q,f,n)
c 
        call prin2('potential evaluated directly*',f,n)
c 
c      . . . and via the structure
c 
        do 2200 i=1,n
c 
        ibox=ind(i)
        z0(1)=z(1,i)
        z0(2)=z(2,i)
        z0(3)=z(3,i)
c 
        call poteva2(z,q,f2(i),ibox,z0,w,boxes,iz)
c 
 2200 continue
        call prin2('potential from poteva2*',f2 ,n)
c 
        call potlist4(boxes,nboxes,w,
     1      z,q,iz,f4,n)
c 
        call prin2('potential from lists 4*',f4,n)
  
        do 2300 i=1,n
        f2(i)=f2(i)+f4(i)
 2300 continue
c 
        call prin2('potential evaluated via the structure*',f2,n)
c 
        errmax=0
        do 2400 i=1,n
        f3(i)=f2(i)-f(i)
c 
        rats(i)=f(i)/f2(i)
c 
        d=dabs(f2(i)-f(i))
        if(errmax .lt. d) errmax=d
 2400 continue
c 
        call prin2('and f3=*',f3,n)
cccc        call prin2('and rats=*',rats,n)
c 
        call prin2('and errmax=*',errmax,1)
        call oneinfo(ier)
        call prinf('and n^2=*',n**2,1)
        return
        end
  
c 
c 
c 
c 
c 
        subroutine potlist4(boxes,nboxes,w,
     1      z,q,iz,f,n)
        implicit real *8 (a-h,o-z)
        save
        integer *4 boxes(15,1),w(1),list(5000)
        dimension z(3,1),q(1),iz(1),z0(3),f(1)
c 
c       set the array f to zero
c 
        do 1200 i=1,n
        f(i)=0
 1200 continue
c 
c        scan all boxes; for those that have lists 4.
c        calculate potentials due to those
c 
        do 2000 ibox=1,nboxes
c 
        itype4=4
        call d3mgetl(ier,ibox,itype4,list,nlist,w)
c 
        if(ier .eq. 4) goto 2000
c 
c       this box has list 4; account for all its elements,
c       with respect to all points in the box ibox
c 
        j1=boxes(14,ibox)
        nj=boxes(15,ibox)
c 
        do 1400 i=j1,j1+nj-1
c 
        j=iz(i)
        z0(1)=z(1,j)
        z0(2)=z(2,j)
        z0(3)=z(3,j)
c 
        call onelist(z,q,iz,z0,f2,list,nlist,boxes)
c 
        j=iz(i)
        f(j)=f(j)+f2
 1400 continue
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
        subroutine wherepts(w,nboxes,iz,ind,n,boxes)
        implicit real *8 (a-h,o-z)
        save
        integer *4 boxes(15,1),ind(1),iz(1),w(1)
        dimension center(3),corners(24)
c 
c        cosntruct the array boxes
c 
        do 1100 ibox=1,nboxes
c 
         call d3mgetb(ier,ibox,boxes(1,ibox),
     1       center,corners,w)
 1100  continue
c 
c       for each point, find the sequence number of the
c       childless box it lives in
c 
        do 1200 i=1,n
c 
        ind(i)=-77
 1200 continue
c 
        do 2000 i=1,nboxes
c 
        if(boxes(6,i) .ne. 0) goto 2000
c 
c       this is a childless box. mark all points that are in it
c 
        j0=boxes(14,i)
        nj=boxes(15,i)
c 
        do 1400 j=j0,j0+nj-1
c 
        k=iz(j)
        ind(k)=i
 1400 continue
c 
 2000 continue
c 
        call prinf('exiting wherepts, ind=*',ind,n)
        return
        end
c 
c 
c 
c 
c 
  
        subroutine poteva(z,q,f,n)
        implicit real *8 (a-h,o-z)
        save
        dimension z(3,1),q(1),f(1)
c 
c       evaluate the potential of charges q at the points z at
c       these same points
c 
        do 1600 i=1,n
        d=0
        do 1200 j=1,n
        if(j .eq. i) goto 1200
c 
        dd=(z(1,i)-z(1,j))**2+(z(2,i)-z(2,j))**2+
     1      (z(3,i)-z(3,j))**2
        d=d+q(j)/dsqrt(dd)
c 
 1200 continue
c 
        f(i)=d
 1600 continue
c 
        return
        end
c 
C 
C 
C 
C 
        SUBROUTINE poteva2(z,q,f,ibox,z0,w,boxes,iz)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        DIMENSION Z(3,1),q(1),z0(3),list(1000)
        integer *4 boxes(15,1),iz(1),w(1)
c 
        f=0
c 
c       account for the interactions of this point with all
c       the points in the list 1 of the childless box containing it
c 
c       . . . find the list 1
c 
        itype1=1
        list(1)=ibox
c 
         call d3mgetl(ier,ibox,itype1,list(2),nlist,w)
c 
        nlist=nlist+1
        if(ier .eq. 4) nlist=1
c 
c 
c        account for the interactions of this point with those in
c        elements of list1 of the childless box containing this point
c 
        call onelist(z,q,iz,z0,f2,list,nlist,boxes)
c 
        f=f+f2
c 
 1400 continue
c 
c        account for the interactions due to the points
c        in elements of the lists 2 of this box and its ancestors
c 
        ibb=ibox
c 
        do 2800 i=1,100
c 
c      get the list 2
c 
         itype2=2
         call d3mgetl(ier,ibb,itype2,list,nlist,w)
c 
        if(ier .eq. 4) goto 2600
c 
c       account for the interactions of this point with all
c       the points in this list2
c 
        call onelist(z,q,iz,z0,f2,list,nlist,boxes)
c 
        f=f+f2
c 
 2600 continue
c 
c       . . . go to the daddy
c 
        ibb2=boxes(5,ibb)
        ibb=ibb2
        if(ibb .eq. 1) goto 3000
 2800 continue
 3000 continue
c 
c        account for the interactions with the list 3 of this box
c 
         itype3=3
         call d3mgetl(ier,ibox,itype3,list,nlist,w)
c 
        if(ier .eq. 4) goto 3200
c 
        call onelist(z,q,iz,z0,f2,list,nlist,boxes)
c 
        f=f+f2
c 
 3200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine onelist(z,q,iz,z0,f,list,nlist,boxes)
        implicit real *8 (a-h,o-z)
        save
        integer *4 boxes(15,1)
        dimension z(3,1),z0(3),iz(1),q(1),list(1)
c 
c       . . . account for the interactions in elements of this
c             list
c 
        f=0
        do 1400 i=1,nlist
c 
        k=list(i)
c 
        j1=boxes(14,k)
        nj=boxes(15,k)
c 
c       account for the interactions
c 
        call onebox(z,q,iz(j1),nj,z0,d)
c 
        f=f+d
c 
 1400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine onebox(z,q,iz,nz,z0,f)
        implicit real *8 (a-h,o-z)
        save
        dimension z(3,1),z0(3),iz(1),q(1)
        data inf0/0/
c 
c       construct the potential at the point z0
c       created by the points in this box
c 
        f=0
c 
        do 1200 i=1,nz
  
        info=info+1
        j=iz(i)
c 
        d=(z(1,j)-z0(1))**2+(z(2,j)-z0(2))**2+
     1      (z(3,j)-z0(3))**2
c 
        if(d .lt. 1.0d-24) goto 1200
c 
        f=f+q(j)/dsqrt(d)
 1200 continue
        return
c 
c 
c 
c 
        entry oneinfo(ier)
        call prinf('and info=*',info,1)
        return
c 
c 
c 
c 
        entry oneinfini(n7)
        n=n7
        return
        end
c 
c 
c 
c 
c 
        subroutine lststest(boxes,nboxes,corners,centers,
     1      z,iz,n,w)
        implicit real *8 (a-h,o-z)
        save
        integer *4 boxes(15,1),iz(1)
        dimension centers(3,1),corners(3,8,1),w(1 000 000),
     1      nums(10),z(3,1)
c 
c        construct all lists for all boxes
c 
        lw=2000 000
        call d3mlsts(ier,boxes,nboxes,corners,w,lw,lused)
c 
        call prinf('after d3mlsts, ier=*',ier,1)
        call prinf('after d3mlsts, lused=*',lused,1)
  
  
        call linkinfo(w,lused2,nums)
  
        call prinf('after linkinfo, lused2=*',lused2,1)
        call prinf('after linkinfo, nums=*',nums,5)
c 
c        test the constructed lists
c 
        call lststst0(boxes,nboxes,z,iz,n,w)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine sliceplt(iw,corners,centers,nboxes,
     1    itype,t,w)
        implicit real *8 (a-h,o-z)
        save
        dimension corners(3,8,1),centers(3,1),icuts(10 000),
     1      xx(5,10 000),yy(5,10 000),cent(3,1000),box(15),w(1)
c 
c       construct the centers and corners of all boxes
c       on all levels
c 
        do 1100 ibox=1,nboxes
c 
         call d3mgetb(ier,ibox,box,centers(1,ibox),
     1       corners(1,1,ibox),w)
 1100 continue
c 
c        find all boxes on all levels cut by this plane
c 
        do 1400 i=1,nboxes
c 
c        . . . find the minimum and maximum extents of this
c              box along the user-specified coordinate
c 
        tmin=1.0d20
        tmax=-1.0d20
c 
        do 1200 j=1,8
c 
        if(corners(itype,j,i) .lt. tmin) tmin=corners(itype,j,i)
        if(corners(itype,j,i) .gt. tmax) tmax=corners(itype,j,i)
 1200 continue
c 
c       check if this box is cut by this plane
c 
        ifcut=1
        if(t .lt. tmin) ifcut=0
        if(t .gt. tmax) ifcut=0
c 
        if(ifcut .eq. 0) goto 1400
c 
        numcut=numcut+1
        icuts(numcut)=i
c 
 1400 continue
c 
        call prinf('sequence numbers of boxes cut by this plane*',
     1      icuts,numcut)
c 
c        now, plot them things
c 
        if(itype .ne. 3) goto 1600
        ix=1
        iy=2
c 
 1600 continue
c 
        if(itype .ne. 2) goto 1800
        ix=1
        iy=3
c 
 1800 continue
c 
        if(itype .ne. 1) goto 2000
        ix=2
        iy=3
c 
 2000 continue
c 
        call prinf('ix=*',ix,1)
        call prinf('iy=*',iy,1)
  
        do 3000 i=1,numcut
c 
        k=icuts(i)
c 
c       construct the square to be plotted
c 
        xmax=-1.0d20
        xmin=1.0d20
c 
        ymax=-1.0d20
        ymin=1.0d20
c 
        do 2200 j=1,8
        if(corners(ix,j,k) .lt. xmin) xmin=corners(ix,j,k)
        if(corners(ix,j,k) .gt. xmax) xmax=corners(ix,j,k)
c 
        if(corners(iy,j,k) .lt. ymin) ymin=corners(iy,j,k)
        if(corners(iy,j,k) .gt. ymax) ymax=corners(iy,j,k)
 2200 continue
c 
c        . . . plot
c 
        xx(1,i)=xmin
        yy(1,i)=ymin
c 
        xx(2,i)=xmax
        yy(2,i)=yy(1,i)
c 
        xx(3,i)=xx(2,i)
        yy(3,i)=ymax
c 
        xx(4,i)=xx(1,i)
        yy(4,i)=yy(3,i)
c 
        xx(5,i)=xx(1,i)
        yy(5,i)=yy(1,i)
c 
        cent(1,i)=centers(1,k)
        cent(2,i)=centers(2,k)
c 
 3000 continue
c 
c        now, plot all this junk
c 
cccc        call prin2('before cutplt, xx=*',xx,numcut*5)
cccc        call prin2('before cutplt, yy=*',yy,numcut*5)
c 
        call cutplt(iw,xx,yy,numcut,'slice*')
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine cutplt(iw,xx,yy,numcut,mes)
        implicit real *8 (a-h,o-z)
        save
        real *8 xx(5,1),yy(5,1)
        character *1 mes(1)
c 
         call prinf('entered allplt12, iw=*',iw,1)
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
        DO 1100 I=1,Numcut
        do 1000 j=1,5
        IF(XMIN. GT. xx(j,i) ) xmin=xx(j,i)
        IF(XMax. lT. xx(j,i) ) xmax=xx(j,i)
c 
        IF(yMIN. GT. yy(j,i) ) ymin=yy(j,i)
        IF(yMax. lT. yy(j,i) ) ymax=yy(j,i)
c 
 1000 CONTINUE
 1100 CONTINUE
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
c       plot all boxes intersecting the user-specified plane
c 
        do 3000 i=1,numcut
c 
 1410 FORMAT(2(2X,E15.7))
        write(iw,1410) (xx(j,i),yy(j,i),j=1,5)
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
        WRITE(IW,1600)
        WRITE(IW,1800)
        close(iw)
        return
        end
c 
c 
c 
c 
c 
        subroutine bwhere(z,iz,n,numb)
        implicit real *8 (a-h,o-z)
        save
        dimension z(3,1),iz(1)
c 
c       find the maximum and minimum x, y, z for all
c       points in this box
c 
        call prinf('looking at box number*',numb,1)
        call prinf('n=*',n,1)
        if(n .eq. 0) return
  
  
        xmin=1.0d20
        ymin=1.0d20
        zmin=1.0d20
c 
        xmax=-1.0d20
        ymax=-1.0d20
        zmax=-1.0d20
c 
        do 1100 i=1,n
        j=iz(i)
c 
        if(z(1,j) .gt. xmax) xmax=z(1,j)
        if(z(2,j) .gt. ymax) ymax=z(2,j)
        if(z(3,j) .gt. zmax) zmax=z(3,j)
c 
        if(z(1,j) .lt. xmin) xmin=z(1,j)
        if(z(2,j) .lt. ymin) ymin=z(2,j)
        if(z(3,j) .lt. zmin) zmin=z(3,j)
c 
 1100 continue
c 
ccc 1200 format(2x, 'xmin=',f6.3,'  xmax= ',f6.3, /,2x,
ccc     1  'ymin= ',f6.3,' ymax= ',f6.3,/,2x,
ccc     2    'zmin= ',f6.3,'  zmax= ',f6.3 )
c 
 1200 format(2x, 'xmin=',e11.5,'  xmax= ',e11.5, /,2x,
     1  'ymin= ',e11.5,' ymax= ',e11.5,/,2x,
     2    'zmin= ',e11.5,'  zmax= ',e11.5 )
        write(6,1200) xmin,xmax,ymin,ymax,zmin,zmax
        write(13,1200) xmin,xmax,ymin,ymax,zmin,zmax
c 
        return
        end
c 
c 
c 
        subroutine cornplot(cent,size)
        implicit real *8 (a-h,o-z)
        save
        dimension cent(1),x(4),y(4)
c 
c       plot the corners for the xy-section
c 
        x(1)=cent(1)-size/2
        y(1)=cent(2)-size/2
c 
        x(2)=cent(1)+size/2
        y(2)=cent(2)-size/2
c 
        x(3)=cent(1)-size/2
        y(3)=cent(2)+size/2
c 
        x(4)=cent(1)+size/2
        y(4)=cent(2)+size/2
c 
        iw=201
        n=4
        call RSPLOT(X,Y,N,Iw,'  *')
        return
        end
c 
c 
c 
c 
c 
        subroutine boxplt1(iw,z,n,mes)
        implicit real *8 (a-h,o-z)
        save
        real *8 z(3,1),x(10 000),y(10 000)
        character *1 mes(1),line(1000)
c 
        iw1=iw+1
        iw2=iw+2
        iw3=iw+3
c 
c        plot the xy-projection (the z-direction collapsed)
c 
        do 1400 i=1,n
        x(i)=z(1,i)
        y(i)=z(2,i)
 1400 continue
c 
        call msgmerge(mes,', xy-projection*',line)
c 
        call RSPLOT(X,Y,N,IW1,line)
c 
c        plot the xz-projection (the z-direction collapsed)
c 
        do 1600 i=1,n
        x(i)=z(1,i)
        y(i)=z(3,i)
 1600 continue
c 
        call msgmerge(mes,', xz-projection*',line)
c 
        call RSPLOT(X,Y,N,IW2,line)
c 
c        plot the yz-projection (the z-direction collapsed)
c 
        do 1800 i=1,n
        x(i)=z(2,i)
        y(i)=z(3,i)
 1800 continue
c 
        call msgmerge(mes,', yz-projection*',line)
c 
        call RSPLOT(X,Y,N,IW3,line)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine ran2in3(nsmall,z,w,n,itype,shift,z2)
        implicit real *8 (a-h,o-z)
        save
        dimension z(3,1),z2(2,1)
c 
c       construct the two-dimensional structure
c 
        call RANLST(NSMALL,Z2,W,N)
c 
c       pack it into the three-dimensional array, shifting
c       it by shift in the direction itype
c 
        if(itype .ne. 3) goto 1600
        ix=1
        iy=2
c 
 1600 continue
c 
        if(itype .ne. 2) goto 1800
        ix=1
        iy=3
c 
 1800 continue
c 
        if(itype .ne. 1) goto 2000
        ix=2
        iy=3
c 
 2000 continue
c 
        call prinf('in ran2in3, ix=*',ix,1)
        call prinf('in ran2in3, iy=*',iy,1)
c 
        do 2200 i=1,n
c 
        z(ix,i)=z2(1,i)
        z(iy,i)=z2(2,i)
        z(itype,i)=shift
 2200 continue
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
c 
c 
c 
c 
c 
        SUBROUTINE RIRANB3(n,SIZE,CENT,Z,W)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        DIMENSION CENT(3),Z(3,1),W(1)
C 
C       GENERATE PSEUDO-RANDOM NUMBERS
C 
        CALL RIRAND(N*3,W)
        DO 1200 I=1,N
        Z(1,I)=W(I)
        Z(2,I)=W(I+N)
        Z(3,I)=W(I+N*2)
 1200 CONTINUE
C 
C       SHIFT AND STRETCH/COMPRESS THE POINTS IN THE PLANE
C 
        DO 1400 I=1,N
        Z(1,I)=Z(1,I)*SIZE+CENT(1)-0.5*SIZE
        Z(2,I)=Z(2,I)*SIZE+CENT(2)-0.5*SIZE
        Z(3,I)=Z(3,I)*SIZE+CENT(3)-0.5*SIZE
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
        subroutine creelips(z,n,a,b)
        implicit real *8 (a-h,o-z)
        save
        dimension z(3,1)
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
        z(3,i)=z(2,i)*1.2
c 
cccc        z(2,i)=0.1+0.01*z(3,i)
 1400 continue
  
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
        subroutine d3mstrcr(ier,z,n,nbox,
     1    nboxes,iz,laddr,nlev,center,size,w,lw,lused777)
        implicit real *8 (a-h,o-z)
        save
        integer *4 iz(1),w(1),laddr(1),box(1),nums(1)
        real *8 z(3,1),center(3),corners(3,8)
c 
c        this subroutine constructs the logical structure for the
c        fully adaptive FMM in two dimensions and stores it in the
c        array w in the form of a link-list. after that, the user
c        can obtain the information about various boxes and lists
c        in it by calling the entries d3mgetb, d3mgetl, d3mlinfo
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
c         (z(1,j),z(2,j),z(3,j)),(z(1,j+1),z(2,j+1),z(3,j+1)),
c         (z(1,j+2),z(2,j+2),z(3,j+3)), . . .
c         (z(1,j+nj-1),z(2,j+nj-1),z(3,j+nj-1)),
c         (z(1,j+nj),z(2,j+nj),z(3,j+nj)),
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
c         it is a link-list (for the most part), and can only be accessed
c         via the entries d3mgetb, d3mgetl, d3mlinfo, of this  subroutine
c         (see below). the first lused 777 integer *4 locations of
c         this array should not be altered between the call to this
c         entry and subsequent calls to the entries d3mgetb, d3mgetl,
c         d3mlinfo, of this  subroutine
c 
c  lused777 - the amount of space in the array w (in integer *4 words)
c        that is occupied by various tables on exit from this
c        subroutine. this space should not be altered between the
c        call to this entry and subsequent calls to entries d3mgetb,
c        d3mgetl, d3mlinfo, of this  subroutine (see below).
c 
c        . . . construct the quad-tree structure for the user-specified
c              set of points
c 
        ninire=2
        iiwork=1
        liwork=n+4
c 
        iboxes=iiwork+liwork
        lboxes=lw-n-5
        maxboxes=lboxes/15-1
c 
        call d3mallb(ier,z,n,nbox,w(iboxes),maxboxes,
     1    nboxes,iz,laddr,nlev,center,size,w(iiwork) )
c 
c        if the memory is insufficient - bomb
c 
        if(ier .eq. 0) goto 1100
           call prinf('in d3mstrcr after d3mallb, ier=*',ier,1)
        if(ier .eq. 4) ier=32
        return
 1100 continue
c 
c       compress the array w
c 
        nn=nboxes*15
        do 1200 i=1,nn
        w(i)=w(iboxes+i-1)
 1200 continue
        iboxes=1
 1300 continue
        lboxes=nboxes*15+16
c 
c       construct the centers and the corners for all boxes
c       in the quad-tree
c 
        icenters=iboxes+lboxes
        lcenters=(nboxes*3+2)*ninire
c 
        icorners=icenters+lcenters
        lcorners=(nboxes*24+2)*ninire
c 
        iwlists=icorners+lcorners
        lwlists=lw-iwlists-6
c 
        call d3mcentc(center,size,w(iboxes),nboxes,
     1      w(icenters),w(icorners) )
c 
c       now, construct all lists for all boxes
c 
  
           call prinf('before d3mlsts, lwlists=*',lwlists,1)
        call d3mlsts(ier,w(iboxes),nboxes,w(icorners),
     1        w(iwlists),lwlists,lused)
c 
        lused777=lused+iwlists
        call prinf('after d3mlsts, ier=*',ier,1)
c 
        nboxes7=nboxes
c 
        return
c 
c 
c 
c 
         entry d3mgetb(ier,ibox,box,center,corners,w)
c 
c        this entry returns to the user the characteristics of
c        user-specified box  ibox.
c 
c                     input parameters:
c 
c  ibox - the box number for which the information is desired
c  w - storage area as created by the entry d3mstrcr (see above)
c 
c                     output parameters:
c 
c  ier - the error return code.
c    ier=0 means successful execution
c    ier=4 means that ibox is either greater than the number of boxes
c           in the structure or less than 1.
c  box - an integer array dimensioned box(15). its elements describe
c        the box number ibox, as follows:
c 
c       1. level - the level of subdivision on which this box
c             was constructed;
c       2, 3, 4  - the coordinates of this box among  all
c             boxes on this level
c       5 - the daddy of this box, identified by it address
c             in array boxes
c       6,7,8,9,10,11,12,13 - the  list of children of this box
c             (eight of them, and the child is identified by its address
c             in the array boxes; if a box has only one child, only the
c             first of the four child entries is non-zero, etc.)
c       14 - the location in the array iz of the particles
c             living in this box
c       15 - the number of particles living in this box
c  center - the center of the box number ibox
c  corners - the corners of the box number ibox
c 
c       . . . return to the user all information about the box ibox
c 
        ier=0
        if( (ibox .ge. 1)  .and. (ibox .le. nboxes7) ) goto 2100
        ier=4
        return
 2100 continue
c 
        ibox0=iboxes+(ibox-1)*15-1
        do 2200 i=1,15
        box(i)=w(ibox0+i)
 2200 continue
c 
c      return to the user the center and the corners of the box ibox
c 
        call d3mcpcc(w(icenters),w(icorners),ibox,center,corners)
c 
        return
c 
c 
c 
c 
         entry d3mgetl(ier,ibox,itype,list,nlist,w)
c 
c  ibox - the box number for which the information is desired
c  itype - the type of the desired list for the box ibox
c  w - storage area as created by the entry d3mstrcr (see above)
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
        entry d3mlinfo(w,lused77,nums)
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
        subroutine d3mlsts(ier,boxes,nboxes,corners,w,lw,lused)
        implicit real *8 (a-h,o-z)
        save
        integer *4 boxes(15,1),collkids(5000),w(1),
     1      dadcolls(200),list5(2000),stack(6000)
        real *8 corners(3,8,1)
c 
c        this subroutine constructs all lists for all boxes
c        and stores them in the storage area w in the form
c        of a link list. the resulting data can be accessed
c        by calls to various entries of the subroutine linkinit (see).
c 
c                          input parameters:
c 
c  boxes - an integer array dimensioned (15,nboxes), as created by
c        the subroutine d3mallb (see).  each 11-element column
c         describes one box, as follows:
c 
c 
c       1. level - the level of subdivision on which this box
c             was constructed;
c       2, 3, 4  - the coordinates of this box among  all
c             boxes on this level
c       5 - the daddy of this box, identified by it address
c             in array boxes
c       6,7,8,9,10,11,12,13 - the  list of children of this box
c             (eight of them, and the child is identified by its address
c             in the array boxes; if a box has only one child, only the
c             first of the four child entries is non-zero, etc.)
c       14 - the location in the array iz of the particles
c             living in this box
c       15 - the number of particles living in this box
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
c        the form of link-lists, accessible by the subroutine
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
         call prinf('in d3mlsts, nboxes=*',nboxes,1)
c 
        ier=0
        ntypes=5
        call linkinit(ier,nboxes,ntypes,w,lw)
cccc         call prinf('in d3mlsts after linkinit, ier=*',ier,1)
c 
c        construct lists 5, 2 for all boxes
c 
        do 2000 ibox=2,nboxes
c 
c       find this guy's daddy
c 
        idad=boxes(5,ibox)
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
        do 1400 j=1,8
        kid=boxes(5+j,icoll)
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
        call d3mifint(corners(1,1,kid),corners(1,1,ibox),ifinter)
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
        call prinf('constructed lists 5, 2; lused=*',lused,1)
c 
c       now, construct lists 1, 3
c 
        do 3000 i=1,nboxes
c 
c       if this box has kids - its lists 1, 3 are empty;
c       do not construct them
c 
        if(boxes(6,i) .gt. 0) goto 3000
c 
        call linkretr(jer,itype5,i,list5,nlist,w,lused)
c 
        if(jer .eq. 4) goto 3000
c 
        do 2200 j=1,nlist
        jbox=list5(j)
        call d3mlst31(ier,i,jbox,boxes,nboxes,
     1    corners,w,stack,lused)
c 
c        if storage capacity has been exceeded - bomb
c 
        if(ier .eq. 32) return
 2200 continue
 3000 continue
c 
        call prinf('constructed lists 1,3, lused=*',lused,1)
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
  
         call prinf('exiting d3mlsts, lused=*',lused,1)
        return
        end
c 
c 
c 
c 
c 
        subroutine d3mlst31(ier,ibox,jbox0,boxes,nboxes,
     1    corners,w,stack,lused)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),corners(3,8,1)
        integer *4 boxes(15,1),stack(3,1)
        data itype1/1/,itype2/2/,itype3/3/,itype4/4/,itype5/5/,
     1      nlist1/1/
c 
c       this subroutine constructs all elements of lists 1 and 3
c       resulting from the subdivision of one element of list 5
c       of the box ibox. all these elements of lists 1, 3 are
c       stored in the link-lists by the subroutine linstro (see)
c 
c        input parameters:
c 
c  ibox - the box whose lists are being constructed
c 
c  jbox0 - the element of list 5 of the box ibox that is being
c          subdfivided
c  boxes - the array boxes as created by the subroutine d3mallb (see)
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
c 
        nsons=0
        do 1200 j=6,13
c 
        if(boxes(j,jbox) .gt. 0) nsons=nsons+1
 1200 continue
c 
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
        call d3mifint(corners(1,1,ibox),corners(1,1,jbox),ifinter)
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
        if(boxes(6,jbox) .ne. 0) goto 3000
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
        jbox=boxes(5+nsons,jbox)
        istack=istack+1
c 
        nsons=0
        do 4600 j=6,13
c 
        if(boxes(j,jbox) .gt. 0) nsons=nsons+1
 4600 continue
c 
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
        subroutine d3mifint(c1,c2,ifinter)
        implicit real *8 (a-h,o-z)
        save
        dimension c1(3,8),c2(3,8)
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
        xmin1=1.0d20
        ymin1=1.0d20
        zmin1=1.0d20
c 
        xmax1=-1.0d20
        ymax1=-1.0d20
        zmax1=-1.0d20
c 
        xmin2=1.0d20
        ymin2=1.0d20
        zmin2=1.0d20
c 
        xmax2=-1.0d20
        ymax2=-1.0d20
        zmax2=-1.0d20
c 
        do 1200 i=1,8
c 
        if(xmin1 .gt. c1(1,i)) xmin1=c1(1,i)
        if(ymin1 .gt. c1(2,i)) ymin1=c1(2,i)
        if(zmin1 .gt. c1(3,i)) zmin1=c1(3,i)
c 
        if(xmax1 .lt. c1(1,i)) xmax1=c1(1,i)
        if(ymax1 .lt. c1(2,i)) ymax1=c1(2,i)
        if(zmax1 .lt. c1(3,i)) zmax1=c1(3,i)
c 
c 
        if(xmin2 .gt. c2(1,i)) xmin2=c2(1,i)
        if(ymin2 .gt. c2(2,i)) ymin2=c2(2,i)
        if(zmin2 .gt. c2(3,i)) zmin2=c2(3,i)
c 
        if(xmax2 .lt. c2(1,i)) xmax2=c2(1,i)
        if(ymax2 .lt. c2(2,i)) ymax2=c2(2,i)
        if(zmax2 .lt. c2(3,i)) zmax2=c2(3,i)
c 
 1200 continue
c 
c        decide if the boxes intersect
c 
        eps=xmax1-xmin1
        if(eps .gt. xmax2-xmin2) eps=xmax2-xmin2
        if(eps .gt. ymax2-ymin2) eps=ymax2-ymin2
        if(eps .gt. zmax2-zmin2) eps=zmax2-zmin2
c 
        if(eps .gt. ymax1-ymin1) eps=ymax1-ymin1
        if(eps .gt. zmax1-zmin1) eps=zmax1-zmin1
c 
        eps=eps/100
c 
        ifinter=1
        if(xmin1 .gt. xmax2+eps) ifinter=0
        if(xmin2 .gt. xmax1+eps) ifinter=0
c 
        if(ymin1 .gt. ymax2+eps) ifinter=0
        if(ymin2 .gt. ymax1+eps) ifinter=0
c 
        if(zmin1 .gt. zmax2+eps) ifinter=0
        if(zmin2 .gt. zmax1+eps) ifinter=0
        return
        end
c 
c 
c 
c 
c 
        subroutine d3mcpcc(centers,corners,ibox,center,corner)
        implicit real *8 (a-h,o-z)
        save
        dimension centers(3,1),corners(3,8,1),center(3),
     1      corner(3,8)
c 
        center(1)=centers(1,ibox)
        center(2)=centers(2,ibox)
        center(3)=centers(3,ibox)
c 
        do 1200 i=1,8
        corner(1,i)=corners(1,i,ibox)
        corner(2,i)=corners(2,i,ibox)
        corner(3,i)=corners(3,i,ibox)
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
        subroutine d3mallb(ier,z,n,nbox,boxes,maxboxes,
     1    nboxes,iz,laddr,nlev,center0,size,iwork)
        implicit real *8 (a-h,o-z)
        save
        integer *4 boxes(15,1),iz(1),laddr(2,1),iwork(1),
     1      is(8),ns(8),iisons(8),jjsons(8),kksons(8)
        real *8 z(3,1),center0(1),center(3)
        data iisons/1,1,1,1,2,2,2,2/,jjsons/1,1,2,2,1,1,2,2/,
     1      kksons/1,2,1,2,1,2,1,2/
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
c        set to 4, and the execution of the subroutine is
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
c  boxes - an integer array dimensioned (15,nboxes). each 15-element
c        column describes one box, as follows:
c 
c       1. level - the level of subdivision on which this box
c             was constructed;
c       2, 3, 4  - the coordinates of this box among  all
c             boxes on this level
c       5 - the daddy of this box, identified by it address
c             in array boxes
c       6,7,8,9,10,11,12,13 - the  list of children of this box
c             (eight of them, and the child is identified by its address
c             in the array boxes; if a box has only one child, only the
c             first of the four child entries is non-zero, etc.)
c       14 - the location in the array iz of the particles
c             living in this box
c       15 - the number of particles living in this box
c 
c    important warning: the array boxes has to be dimensioned
c                       at least (15,maxboxes)!! otherwise,
c                       the subroutine is likely to bomb, since
c                       it assumes that it has that much space!!!!
c  nboxes - the total number of boxes created
c  iz - the integer array addressing the particles in
c         all boxes.
c       explanation: for a box ibox, the particles living in
c         it are:
c 
c         (z(1,j),z(2,j),z(3,j)),(z(1,j+1),z(2,j+1),z(3,j+1)),
c         (z(1,j+2),z(2,j+2),z(3,j+3)), . . .
c         (z(1,j+nj-1),z(2,j+nj-1),z(3,j+nj-1)),
c         (z(1,j+nj),z(2,j+nj),z(3,j+nj)),
c 
c         with j=boxes(9,ibox), and nj=boxes(10,ibox)
c 
c  laddr - an integer array dimensioned (2,numlev), containing
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
cccc          call prinf('in d3mallb, maxboxes=*',maxboxes,1)
c 
        ier=0
        xmin=1.0d50
        xmax=-xmin
        ymin=1.0d50
        ymax=-ymin
        zmin=1.0d50
        zmax=-zmin
c 
        do 1100 i=1,n
        if(z(1,i) .lt. xmin) xmin=z(1,i)
        if(z(1,i) .gt. xmax) xmax=z(1,i)
        if(z(2,i) .lt. ymin) ymin=z(2,i)
        if(z(2,i) .gt. ymax) ymax=z(2,i)
        if(z(3,i) .lt. zmin) zmin=z(3,i)
        if(z(3,i) .gt. zmax) zmax=z(3,i)
 1100 continue
        size=xmax-xmin
        sizey=ymax-ymin
        sizez=zmax-zmin
        if(sizey .gt. size) size=sizey
        if(sizez .gt. size) size=sizez
c 
        center0(1)=(xmin+xmax)/2
        center0(2)=(ymin+ymax)/2
        center0(3)=(zmin+zmax)/2
c 
cccc         call prin2('in d3mallb, center0=*',center0,3)
cccc         call prin2('in d3mallb, size=*',size,1)
c 
        boxes(1,1)=0
        boxes(2,1)=1
        boxes(3,1)=1
        boxes(4,1)=1
        boxes(5,1)=0
        boxes(6,1)=0
        boxes(7,1)=0
        boxes(8,1)=0
        boxes(9,1)=0
        boxes(10,1)=0
        boxes(11,1)=0
        boxes(12,1)=0
        boxes(13,1)=0
        boxes(14,1)=1
        boxes(15,1)=n
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
c 
cccc         call prinf('in d3mallb, maxson=*',maxson,1)
c 
        maxlev=200
        ison=1
        nlev=0
cccc         call prinf('in d3mallb, nbox=*',nbox,1)
cccc         call prinf('in d3mallb, n=*',n,1)
        do 3000 level=0,maxlev
cccc          call prinf('in d3mallb, level=*',level,1)
        laddr(1,level+2)=laddr(1,level+1)+laddr(2,level+1)
        nlevson=0
        idad0=laddr(1,level+1)
        idad1=idad0+laddr(2,level+1)-1
c 
        do 2000 idad=idad0,idad1
c 
c       subdivide the box number idad (if needed)
c 
        numpdad=boxes(15,idad)
        if(numpdad .le. nbox) goto 2000
c 
        ii=boxes(2,idad)
        jj=boxes(3,idad)
        kk=boxes(4,idad)
        call d3mcentf(center0,size,level,ii,jj,kk,center)
c 
        iiz=boxes(14,idad)
        nz=boxes(15,idad)
c 
cccc        call prinf('before d2msepa1, nz=*',nz,1)
c 
        call d3msepa1(center,z,iz(iiz),nz,iwork,
     1    is,ns)
c 
cccc        call prinf('after d3msepa1, is=*',is,8)
cccc        call prinf('after d3msepa1, ns=*',ns,8)
c 
c 
c       store in array boxes the sons obtained by the routine
c       d3msepa1
c 
         idadson=6
        do 1600 i=1,8
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
        do 1500 lll=6,13
        boxes(lll,ison)=0
 1500 continue
c 
        boxes(1,ison)=level+1
        iison=(ii-1)*2+iisons(i)
        jjson=(jj-1)*2+jjsons(i)
        kkson=(kk-1)*2+kksons(i)
        boxes(2,ison)=iison
        boxes(3,ison)=jjson
        boxes(4,ison)=kkson
        boxes(5,ison)=idad
c 
        boxes(14,ison)=is(i)+iiz-1
        boxes(15,ison)=ns(i)
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
        subroutine d3mcentc(center0,size,boxes,nboxes,
     1      centers,corners)
        implicit real *8 (a-h,o-z)
        save
        integer *4 boxes(15,1)
        dimension centers(3,1),corners(3,8,1),center(3),
     1      center0(3)
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
c        by the subroutine d3mallb (see)
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
        z00=center0(3)-size/2
        do 1400 i=1,nboxes
        level=boxes(1,i)
        side=size/2**level
        side2=side/2
        ii=boxes(2,i)
        jj=boxes(3,i)
        kk=boxes(4,i)
        center(1)=x00+(ii-1)*side+side2
        center(2)=y00+(jj-1)*side+side2
        center(3)=z00+(kk-1)*side+side2
c 
        centers(1,i)=center(1)
        centers(2,i)=center(2)
        centers(3,i)=center(3)
c 
        corners(1,1,i)=center(1)-side/2
        corners(1,2,i)=corners(1,1,i)
        corners(1,3,i)=corners(1,1,i)
        corners(1,4,i)=corners(1,1,i)
c 
        corners(1,5,i)=corners(1,1,i)+side
        corners(1,6,i)=corners(1,5,i)
        corners(1,7,i)=corners(1,5,i)
        corners(1,8,i)=corners(1,5,i)
c 
c 
        corners(2,1,i)=center(2)-side/2
        corners(2,2,i)=corners(2,1,i)
        corners(2,5,i)=corners(2,1,i)
        corners(2,6,i)=corners(2,1,i)
c 
        corners(2,3,i)=corners(2,1,i)+side
        corners(2,4,i)=corners(2,3,i)
        corners(2,7,i)=corners(2,3,i)
        corners(2,8,i)=corners(2,3,i)
c 
c 
        corners(3,1,i)=center(3)-side/2
        corners(3,3,i)=corners(3,1,i)
        corners(3,5,i)=corners(3,1,i)
        corners(3,7,i)=corners(3,1,i)
c 
        corners(3,2,i)=corners(3,1,i)+side
        corners(3,4,i)=corners(3,2,i)
        corners(3,6,i)=corners(3,2,i)
        corners(3,8,i)=corners(3,2,i)
c 
 1400 continue
         return
         end
c 
c 
c 
c 
c 
        subroutine d3mcentf(center0,size,level,i,j,k,center)
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
        z0=center0(3)-size/2
        level0=level
 1200 continue
        center(1)=x0+(i-1)*side+side2
        center(2)=y0+(j-1)*side+side2
        center(3)=z0+(k-1)*side+side2
        return
        end
  
c 
c 
c 
c 
c 
        subroutine d3msepa1(cent,z,iz,n,iwork,
     1    is,ns)
        implicit real *8 (a-h,o-z)
        save
        dimension cent(1),z(3,1),iz(1),iwork(1),is(1),ns(1)
c 
c        this subroutine reorders the particles in a box,
c        so that each of the children occupies a contigious
c        chunk of array iz
c 
c        note that we are using a strange numbering convention
c        for the children:
c 
c             3,4     7,8
c 
c                                     <- looking down the z-axis
c             1,2     5,6
c 
c 
c                        input parameters:
c 
c  cent - the center of the box to be subdivided
c  z - the list of all points in the box to be subdivided
c  iz - the integer array specifying the transposition already
c       applied to the points z, before the subdivision of
c       the box into children
c  n - the total number of points in array z
c 
c                        output parameters:
c 
c  iz - the integer array specifying the transposition already
c       applied to the points z, after the subdivision of
c       the box into children
c  is - an integer array of length 8 containing the locations
c       of the sons in array iz
c  ns - an integer array of length 8 containig the numbers of
c       elements in the sons
c 
c                        work arrays:
c 
c  iwork - must be n integer *4 elements long
c 
c        . . . separate all particles in this box in x
c 
        n1=0
        n2=0
        n3=0
        n4=0
        n5=0
        n6=0
        n7=0
        n8=0
c 
        n12=0
        n34=0
        n56=0
        n78=0
c 
        n1234=0
        n5678=0
c 
        itype=1
        thresh=cent(1)
        call d3msepa0(z,iz,n,itype,thresh,iwork,n1234)
        n5678=n-n1234
c 
c       at this point, the contents of sons number 1,2,3,4 are in
c       the part of array iz with numbers 1,2,...n1234
c       the contents of sons number 5,6,7,8  are in
c       the part of array iz with numbers n1234+1,n12+2,...,n
c 
c        . . . separate the boxes 1, 2, 3, 4 and boxes 5, 6, 7, 8
c 
        itype=2
        thresh=cent(2)
        if(n1234 .ne. 0)
     1    call d3msepa0(z,iz,n1234,itype,thresh,iwork,n12)
        n34=n1234-n12
c 
        if(n5678 .ne. 0)
     1    call d3msepa0(z,iz(n1234+1),n5678,itype,thresh,iwork,n56)
        n78=n5678-n56
c 
c       perform the final separation of pairs of sonnies into
c       individual ones
c 
        itype=3
        thresh=cent(3)
        if(n12 .ne. 0)
     1    call d3msepa0(z,iz,n12,itype,thresh,iwork,n1)
        n2=n12-n1
c 
        if(n34 .ne. 0)
     1    call d3msepa0(z,iz(n12+1),n34,itype,thresh,iwork,n3)
        n4=n34-n3
c 
        if(n56 .ne. 0)
     1    call d3msepa0(z,iz(n1234+1),n56,itype,thresh,iwork,n5)
        n6=n56-n5
c 
        if(n78 .ne. 0)
     1    call d3msepa0(z,iz(n1234+n56+1),n78,itype,thresh,iwork,n7)
        n8=n78-n7
c 
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
        is(5)=is(4)+ns(4)
        ns(5)=n5
c 
        is(6)=is(5)+ns(5)
        ns(6)=n6
c 
        is(7)=is(6)+ns(6)
        ns(7)=n7
c 
        is(8)=is(7)+ns(7)
        ns(8)=n8
c 
cccc        call prinf('is as created*',is,8)
cccc        call prinf('ns as created*',ns,8)
        return
        end
c 
c 
c 
c 
c 
        subroutine d3msepa0(z,iz,n,itype,thresh,iwork,n1)
        implicit real *8 (a-h,o-z)
        save
        dimension z(3,1),iz(1),iwork(1)
c 
c       subdivide the points in this box, horizontally
c       or vertically, depending on the parameter itype
c 
cccc        call prin2('in d3msepa0,thresh=*',thresh,1)
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
        subroutine linkinit(ier,nboxes,ntypes,w,lw)
        save
        integer *4 w(1),nums(1),inums(20)
        data iilistad/1/,iilists/2/,inumele/3/,inboxes/4/,
     1      intypes/5/,ilw/6/,iltot/7/,
     2      inums/11,12,13,14,15,16,17,18,19,20,21,22,23,24,
     3            25,26,27,28,29,30/
c 
c        this is the initialization entry point for the link-list
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
