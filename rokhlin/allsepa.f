        implicit real *8 (a-h,o-z)
        dimension
     1      w(300000),x(30000),y(30000),z(2,30000),
     2      iz(30000),
     6      w1(30000),w2(30000),
     7      z2(2,30000),center0(2)
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
         PRINT *, 'ENTER ntheta'
         READ *,ntheta
         CALL PRINF('ntheta=*',ntheta,1 )
c 
c       construct and plot the test points
c 
        call RANLST(NSMALL,Z,W,N)
        a=10
        b=1
c 
        call creelips(z,n,a,b)
          n=n-1
         call prinf('after ranlst, n=*',n,1)
c 
        do 1200 i=1,n
        x(i)=z(1,i)
        y(i)=z(2,i)
 1200 continue
        iw=21
        call RSPLOT(X,Y,N,IW,'random points as created*')
c 
c       construct the optimal separator
c 
       eps=1.0d-12
       ifboxes=1
        maxboxes=30000
        iw1=38
        iw2=0
        iw3=42
c 
        lw=300000
        call allsepa(ier,n,z,abest,bbest,cbest,iz,n1,
cccc     1    w,lw,ntheta,eps,iw1,iw2,iw3,nbox)
     1    w,lw,ntheta,eps,nbox)
            call prinf('after allsepa, ier=*',ier,1)
c 
c       now, plot the separator together with all points
c 
         iw3=42
           call allplt3(iw3,z,n,
     2    center0,size,abest,bbest,cbest,
     3      ' the best separator*')
c 
c       now, plot the separator together with the points
c       on one side
c 
        do 2200 i=1,n1
        j=iz(i)
        z2(1,i)=z(1,j)
        z2(2,i)=z(2,j)
         w1(i)=abest*z2(1,i)+bbest*z2(2,i)+cbest
 2200 continue
cccc          call prin2('w1=*',w1,n1)
c 
         iw3=43
           call allplt3(iw3,z2,n1,
     2    center0,size,abest,bbest,cbest,
     3      ' the best separator*')
c 
c       now, plot the separator together with the points
c       on the other side
c 
        do 2400 i=1,n-n1
        ii=n1+i
        j=iz(ii)
        z2(1,i)=z(1,j)
        z2(2,i)=z(2,j)
         w2(i)=abest*z2(1,i)+bbest*z2(2,i)+cbest
 2400 continue
cccc          call prin2('w2=*',w2,n-n1)
c 
         iw3=44
           call allplt3(iw3,z2,n-n1,
     2    center0,size,abest,bbest,cbest,
     3      ' the best separator*')
c 
  
  
  
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine allsepa(ier,n,z,a,b,c,iz,n1,
     1    w,lw,ntheta,eps,nbox)
        implicit real *8 (a-h,o-z)
        save
        dimension z(2,1),w(1),iz(1),center0(2)
c 
c       given a set of n user-specified points z(i) in the
c       plane, this subroutine constructs a line that will
c       separate  them in an optimal fashion into two
c       subsets, each of which containing  n/2 points; if n
c       is odd, one subset will contain (n-1)/2 points
c       and the other will contain (n+1)/2 points.
c 
c 
c              input parameters:
c 
c  n - the number of elements in z
c  z - the set of points in the plane
c  lw - the length of the work array w in real *8 locations
c  ntheta - the number of points into which the angle [0, pi]
c         is to be subdivided while searching for the optimal
c         separator. recommended value: 100.
c  eps - the machine zero; it is assumed that anything that is
c        closer than eps to the separator line might be on either
c        side of it. recommended value: 1.0d-12, provided that
c        the size of tyhe main box is 1.
c  nbox - the localization parameter in the estimation of the
c        interaction rank. recommended value: 3
c 
c              output parameters:
c 
c  ier - the error return code.
c    ier=0 means successful execution
c    ier=2 means that the subroutine has run out of memory, after
c          having created the quad-tree structure succcessfully.
c          doubling the available space in w PROBABLY will be
c          sufficient.
c    ier=4 means that the subroutine has run out of memory during
c          the construction of the quad-tree structure .
c          doubling the available space in w PROBABLY will not be
c          sufficient.
c    ier=16 means that the subroutine attempted to construct more
c        than 199 levels of subdivision. serious trouble.
c  abest, bbest, cbest - the normal coordinates of the separator
c        line
c  iz - the integer array addressing the particles in
c         all boxes as produced by the subroutine d2mallb (see)
c  n1 - the number of points in the array z on the negative
c        side of the separator line
c  iz - an integer array of length n containing the coordinates
c        in the array  z  of points on both sides of the separator.
c 
c      explanation: for all  1 \leq   i \leq n1,
c            all points z(iz(i)) are on the negative side
c            of the separator line.
c            for all  n1+1 \leq   i \leq n,
c            all points z(iz(i)) are on the positive side
c            of the separator line.
c 
c                  work arrays:
c 
c  nums - must be at least ntheta+2 integer*4 elements long.
c  w1, w2 - must be at least ntheta+2 integer*4 elements long each
c 
c        allocate memory for the construction of the
c        quad-tree structure
c 
        ninire=2
        iw1=0
c 
        ier=0
c 
        iladdr=1
        lladdr=410
c 
        iiwork=iladdr+lladdr
        liwork=n/ninire+5
c 
        iboxes=iiwork+liwork
        maxboxes=(lw-iboxes)/10-3
c 
c       construct the quad-tree structure
c 
         call prinf('before d2mallb, n=*',n,1)
        call d2mallb(ier,z,n,nbox,w(iboxes),maxboxes,
     1    nboxes,iz,w(iladdr),nlev,center0,size,w(iiwork) )
        call prinf('after d2mallb, ier=*',ier,1)
        call prinf('after d2mallb, nlev=*',nlev,1)
        call prinf('after d2mallb, nboxes=*',nboxes,1)
         if(ier .ne. 0) return
c 
c        now, plot the whole structure
c 
        if(iw1 .ne. 0)
     1      call allplt(iw1,z,n,w(iboxes),nboxes,
     2    center0,size,'structure as created*')
c 
c        find the optimal separator
c 
c        . . . allocate memory
c 
  
        lboxes=nboxes*10/ninire+10
c 
        icorners=iboxes+lboxes
        lcorners=nboxes*4*2+10
c 
        ithetas=icorners+lcorners
        lthetas=ntheta+10
         lthetas=1
c 
        inums=ithetas+lthetas
        lnums=ntheta/ninire+10
c 
        iww1=inums+lnums
        lww1=n+10
c 
        iww2=iww1+lww1
        lww2=n+10
c 
        icenters=iww2+lww2
        lcenters=nboxes*2
c 
        ltot=icenters+lcenters+1
         call prinf('ltot as computed is*',ltot,1)
        if(ltot .le. lw) goto 2000
        ier=2
        return
 2000 continue
c 
c         . . . construct the corners and centers of all boxes
c 
        call centcorn(center0,size,w(iboxes),nboxes,
     1      w(icenters),w(icorners) )
c 
c        . . .  perform the actual search for the separator
c 
        call sepaopt(ier,w(iboxes),nlev,w(icorners),n,z,
     1      ntheta,w(inums),eps,w(iww1),w(iww2),
ccc     2      w(icenters),maxboxes,a,b,c,n1,center0,size)
     2      w(icenters),maxboxes,center0,size,a,b,c,n1,iz)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine sepaopt(ier,boxes,nlev,corners,n,z,
     1      ntheta,nums,eps,w1,w2,
     2      centers,maxboxes,center0,size,
     3      abest,bbest,cbest,n1,iz)
        implicit real *8 (a-h,o-z)
        save
        integer *4 boxes(10,1),nums(1)
        dimension z(2,1),corners(2,4,1),center0(2),
     1   w1(1),w2(1),centers(2,1),iz(1)
c 
c       given a set of n user-specified points in the plane,
c       this subroutine constructs a line that will
c       separate  them in an optimal fashion into two
c       subsets, each of which contain n/2 points; if n
c       is odd, one subset will contain (n-1)/2 points
c       and the other will contain (n+1)/2 points
c 
c 
c              input parameters:
c 
c  boxes - an integer array dimensioned (10,nboxes), as produced
c        by the subroutine d2mallb (see)
c        column describes one box, as follows:
c  nlev - the maximum level number on which any boxes have
c         been created by the subroutine d2mallb (see)
c  corners - the corners of all boxes in the array boxes
c  n - the number of elements in z
c  z - the set of points in the plane
c  iz - the integer array addressing the particles in
c         all boxes as produced by the subroutine d2mallb (see)
c  ntheta - the number of points into which the angle [0, pi]
c         is to be subdivided while searching for the optimal
c         separator
c  eps - the machine zero; it is assumed that anything that is
c        closer than eps to the separator line might be on either
c        side of it
c  centers - the centers of all boxes in the array boxes
c  maxboxes - the maximum total number of boxes the subroutine
c        is permitted to create. if the points z are such that
c        more boxes are needed, the error return code ier is
c        set to 4, and the execution of thye subroutine is
c        terminated.
c  center0 - the center of the box on the level 0, containing
c         the whole simulation
c  size - the side of the box on the level 0
c 
c 
c              output parameters:
c 
c  ier - the error return code.
c    ier=0 means successful execution
c    ier=4 means that the subroutine attempted to create more
c        than maxboxes boxes
c    ier=16 means that the subroutine attempted to construct more
c        than 199 levels of subdivision.
c  abest, bbest, cbest - the normal coordinates of the separator
c        line
c  n1 - the number of points in the array z on the negative
c        side of the separator line
c  iz - an integer array of length n containing the coordinates
c        in the array  z  of points on both sides of the separator.
c 
c      explanation: for all  1 \leq   i \leq n1,
c            all points z(iz(i)) are on the negative side
c            of the separator line.
c            for all  n1+1 \leq   i \leq n,
c            all points z(iz(i)) are on the positive side
c            of the separator line.
c 
c                  work arrays:
c 
c  nums - must be at least ntheta+2 integer*4 elements long.
c  w1, w2 - must be at least ntheta+2 integer*4 elements long each
c 
c       . . . one angle after another, find the separator and
c             determine its cost
c 
        done=1
        pi=datan(done)*4
        htheta=pi/ntheta
        numbest=10000000
c 
        do 3000 i=1,ntheta
        theta=(i-1)*htheta+htheta/2
c 
c       use bisection to construct the separator in the direction
c       theta
c 
        ifcheck=0
        iw=0
        call fndsep(boxes,nlev,corners,n,z,iz,
     1    theta,w1,w2,eps,a,b,c,iw,ifcheck)
c 
c       estimate the interaction rank between the parts
c 
        call rankest(boxes,nlev,iz,centers,n,z,
     1    a,b,c,num11,num12,w1,w2,size)
         nums(i)=num11
         if(num12 .lt. num11) nums(i)=num12
c 
cccc           call prinf('i=*',i,1)
cccc           call prinf('nums(i)=*',nums(i),1)
         if(nums(i) .ge. numbest) goto 3000
         numbest=nums(i)
         ibest=i
         abest=a
         bbest=b
         cbest=c
 3000 continue
         call prinf('and nums as obtained are*',nums,ntheta)
         call prinf('and numbest=*',numbest,1)
  
        call prinf('ibest=*',ibest,1)
c 
c       finally, subdivide the points in z into two
c       subsets, each in its own part of the plane
c 
        call d2msep(z,iz,n,abest,bbest,cbest,w1,n1)
         call prinf('after d2msep, n1=*',n1,1)
c 
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
c 
c 
c 
c 
c 
        subroutine d2msep(z,iz,n,a,b,c,iwork,n1)
        implicit real *8 (a-h,o-z)
        save
        dimension z(2,1),iz(1),iwork(1)
c 
c       subdivide the points in this box, putting those on
c       the negative side of the separator in the first part
c       of the array iz, and the rest in the second part
c 
         call prin2('in d2msep, a=*',a,1)
         call prin2('in d2msep, b=*',b,1)
         call prin2('in d2msep, c=*',c,1)
        i1=0
        i2=0
        do 1400 i=1,n
        d=a*z(1,i)+b*z(2,i)+c
ccc        if(d .gt. 0) goto 1200
       if(d .ge. 0) goto 1200
        i1=i1+1
        iz(i1)=i
        goto 1400
c 
 1200 continue
        i2=i2+1
        iwork(i2)=i
 1400 continue
c 
        do 1600 i=1,i2
        iz(i1+i)=iwork(i)
 1600 continue
        n1=i1
          call prinf('in d2msep, n1=*',n1,1)
        return
        end
c 
c 
c 
c 
c 
        subroutine rankest(boxes,nlev,iz,centers,n,z,
     1    a,b,c,num1,num2,kids,dads,size)
        implicit real *8 (a-h,o-z)
        save
        integer *4 boxes(10,1),kids(1),dads(1)
        dimension z(2,1),centers(2,1),iz(1)
c 
c       given the user-specified collection of points z,
c       this subroutine estimates the rank of interactions
c       of each of the subset with the separator. it is assumed that
c       a quad-tree structure has been imposed on the nodes
c       z by a prior call to the subroutine d2mallb (see).
c       note also that the separator has the form
c 
c          a*x+b*y+c=0,                                       (1)
c 
c          with a**2+b**2=0.                                  (2)
c 
c               input parameters:
c 
c  boxes, nlev, iz - elements of the quad-tree structure as
c        created by the subroutine d2mallb
c  centers - the centers of all boxes in the array boxes
c  n - the number of points in the array z
c  a,b,c - the coefficients in the normal equation (1) of the
c        separator line
c  size - the side of the box on the level 0, as produced by the
c        subroutine d2mallb
c 
c               output parameters:
c 
c  num1 - the estimated rank on the negative side of the separator
c  num2 - the estimated rank on the positive side of the separator
c 
c               work arrays:
c 
c  kids, dads - each must be at least n integer *4 elements long
c 
c      . . . initialize the search of the boxes on the first level
c 
        icost=1
        kids(1)=1
        numkids=1
        two=2
        num1=0
        num2=0
c 
c       process levels one after another
c 
cccccc         call prin2('inside rankest, a=*',a,1)
cccccc         call prin2('inside rankest, b=*',b,1)
cccccc         call prin2('inside rankest, c=*',c,1)
        do 4000 level=1,nlev
         thresh=size/2**level*dsqrt(two)*2
c 
        do 1100 i=1,numkids
        dads(i)=kids(i)
 1100 continue
        numdads=numkids
        numkids=0
c 
        do 3000 j=1,numdads
        i=dads(j)
c 
c       determine the distance from the center of this
c       box to the separator
c 
        d=a*centers(1,i)+b*centers(2,i)+c
c 
c       determine if this box is separated from the separator
c 
        if(dabs(d) .lt. thresh) goto 1600
c 
c       this box is separated from the separator.
c       act accordingly
c 
        if(d .gt. 0) goto 1200
        num1=num1+icost
        goto 3000
 1200 continue
        num2=num2+icost
  
        goto 3000
 1600 continue
c 
c       this box is not separated from the separator.
c       act accordingly
c 
c       . . . if this box has kids - store this
c             fact in array kids
c 
        if(boxes(5,i) .eq. 0) goto 1800
        numkids=numkids+1
        kids(numkids)=boxes(5,i)
c 
        if(boxes(6,i) .eq. 0) goto 3000
        numkids=numkids+1
        kids(numkids)=boxes(6,i)
c 
        if(boxes(7,i) .eq. 0) goto 3000
        numkids=numkids+1
        kids(numkids)=boxes(7,i)
c 
        if(boxes(8,i) .eq. 0) goto 3000
        numkids=numkids+1
        kids(numkids)=boxes(8,i)
        goto 3000
 1800 continue
c 
c        if this box has no kids - look at the particles in it
c 
        j1=boxes(9,i)
        nj=boxes(10,i)
        do 2200 l=j1,j1+nj-1
c 
        jj=iz(l)
        x=z(1,jj)
        y=z(2,jj)
        d1=a*x+b*y+c
c 
c       this particle is not close to the separator.
c       act accordingly
c 
        if(d1 .gt. 0) num2=num2+1
        if(d1 .le. 0) num1=num1+1
 2200 continue
 3000 continue
        if(numkids .eq. 0) goto 4100
 4000 continue
 4100 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine fndsep(boxes,nlev,corners,n,z,iz,
     1    theta,w1,w2,eps,a,b,c,iw0,ifcheck)
        implicit real *8 (a-h,o-z)
        save
        integer *4 boxes(10,1)
        dimension z(2,1),corners(2,4,1),cc(4),w1(1),w2(1)
c 
c       find the starting point for the bisection process
c 
        a=dsin(theta)
        b=dcos(theta)
c 
        cc(1)=-(a*corners(1,1,1)+b*corners(2,1,1))
        cc(2)=-(a*corners(1,2,1)+b*corners(2,2,1))
        cc(3)=-(a*corners(1,3,1)+b*corners(2,3,1))
        cc(4)=-(a*corners(1,4,1)+b*corners(2,4,1))
c 
        c1=cc(1)
        c2=cc(1)
        do 1200 i=2,4
        if(c1 .gt. cc(i)) c1=cc(i)
        if(c2 .lt. cc(i)) c2=cc(i)
 1200 continue
c 
c       check the values of the two obtained separator
c       functions at all points z
c 
         if(ifcheck .eq. 0) goto 1500
        do 1400 i=1,n
        w1(i)=a*z(1,i)+b*z(2,i)+c1
        w2(i)=a*z(1,i)+b*z(2,i)+c2
 1400 continue
         call prin2('first separator function at all z*',
     1       w1,n)
         call prin2('second separator function at all z*',
     1       w2,n)
 1500 continue
c 
c       perform the bisection
c 
          if(iw0 .eq. 0) goto 1600
         iw=iw0+1
        call allplt3(iw,z,n,
     1      center0,size,a,b,c1,'structure with separator*')
         iw=iw+1
        call allplt3(iw,z,n,
     1      center0,size,a,b,c2,'structure with separator*')
 1600 continue
        do 3000 i=1,100
        niter=i
        c3=(c1+c2)/2
          if(iw0 .eq. 0) goto 1800
         iw=iw+1
        call allplt3(iw,z,n,
     1      center0,size,a,b,c3,'structure with separator*')
 1800 continue
c 
        call countpts(boxes,nlev,iz,corners,n,z,
     1    a,b,c3,num1,num2,numclose,w1,w2,eps)
c 
        idiff=num2-num1
        if((idiff .le. 1) .and. (idiff .ge. -1) ) goto 3200
        if(idiff .gt. 0) c2=c3
        if(idiff .lt. 0) c1=c3
 3000 continue
 3200 continue
cccc         call prinf('in fndsep, niter=*',niter,1)
        c=c3
c 
         return
         end
c 
c 
c 
c 
c 
        subroutine countpts(boxes,nlev,iz,corners,n,z,
     1    a,b,c,num1,num2,numclose,kids,dads,eps)
        implicit real *8 (a-h,o-z)
        save
        integer *4 boxes(10,1),kids(1),dads(1)
        dimension z(2,1),corners(2,4,1),iz(1),coecorn(4)
c 
c       given the user-apecified collection of points z,
c       this subroutine evaluates the numbers of nodes z on
c       each sides of the user-specified line. it is assumed that
c       a quad-tree structure has been imposed on the nodes
c       z by a prior call to the subroutine d2mallb (see).
c       note also that the separator has the form
c 
c          a*x+b*y+c=0,                                       (1)
c 
c          with a**2+b**2=0.                                  (2)
c 
c               input parameters:
c 
c  boxes, nlev, iz - elements of the quad-tree structure as
c        created by the subroutine d2mallb
c  corners - the corners of all boxes in the array boxes
c  n - the number of points in the array z
c  a,b,c - the coefficients in the normal equation (1) of the
c        separator line
c  eps - the machine zero; it is assumed that anything that is
c        closer than eps to the separator line might be on either
c        side of it
c 
c               output parameters:
c 
c  num1 - the number of points in the array z such that
c           a*z(1,i)+b*z(2,i)+c < 0
c  num2 - the number of points in the array z such that
c           a*z(1,i)+b*z(2,i)+c > 0
c  numclose - the number of nodes in array z that are closer to
c        the separator than eps; these have been allocated between
c        num1, num2 in such a manner as to minimize the difference
c        between num1, num2
c 
c               work arrays:
c 
c  kids, dads - each must be at least n integer *4 elements long
c 
c        initialize the search of the boxes on the first level
c 
        kids(1)=1
        numkids=1
c 
c       . . . process levels one after another
c 
        numclose=0
        num1=0
        num2=0
        do 4000 level=1,nlev
cccc        call prinf('in countpts, level=*',level,1)
cccc        call prinf('numkids=*',numkids,1)
c 
        do 1100 i=1,numkids
        dads(i)=kids(i)
 1100 continue
        numdads=numkids
        numkids=0
c 
        do 3000 j=1,numdads
cccc        call prinf(' j=*',j,1)
        i=dads(j)
cccc        call prinf('and i=*',i,1)
cccc        call prinf('while boxes(1,i)=*',boxes(1,i),10)
c 
c       determine the locations of all four corners
c       of this box with respect to the separator
c 
        c1=a*corners(1,1,i)+b*corners(2,1,i)+c
        c2=a*corners(1,2,i)+b*corners(2,2,i)+c
        c3=a*corners(1,3,i)+b*corners(2,3,i)+c
        c4=a*corners(1,4,i)+b*corners(2,4,i)+c
         coecorn(1)=c1
         coecorn(2)=c2
         coecorn(3)=c3
         coecorn(4)=c4
cccc         call prin2('coecorn are*',coecorn,4)
c 
c       find the minimum and maximum of the coefficients c
c 
        cmin=c1
        cmax=c1
        if(c2 .lt. cmin) cmin=c2
        if(c2 .gt. cmax) cmax=c2
c 
        if(c3 .lt. cmin) cmin=c3
        if(c3 .gt. cmax) cmax=c3
c 
        if(c4 .lt. cmin) cmin=c4
        if(c4 .gt. cmax) cmax=c4
c 
c       determine if this box is on one side of the separator
c 
cccc         call prin2('cmin=*',cmin,1)
cccc         call prin2('cmax=*',cmax,1)
        d=cmin*cmax
cccc         call prin2('d=*',d,1)
cccc         call prin2('eps=*',eps,1)
        if(dabs(d) .lt. eps) goto 1600
        if(d .lt. 0) goto 1600
c 
c       this box is on one side of the separator.
c       act accordingly
c 
        if(cmin .gt. 0) goto 1200
        num1=num1+boxes(10,i)
cccc         call prinf('i=*',i,1)
cccc         call prinf('num1 as created*',num1,1)
        goto 3000
 1200 continue
        num2=num2+boxes(10,i)
cccc         call prinf('num2 as created*',num2,1)
        goto 3000
 1600 continue
c 
c       this box is crossed by the separator.
c       act accordingly
c 
c       . . . if this box has kids - store this
c             fact in array kids
c 
cccc          call prinf('box crossed, boxes(5,i)=*',
cccc     1        boxes(5,i),4)
        if(boxes(5,i) .eq. 0) goto 1800
        numkids=numkids+1
        kids(numkids)=boxes(5,i)
c 
        if(boxes(6,i) .eq. 0) goto 3000
        numkids=numkids+1
        kids(numkids)=boxes(6,i)
c 
        if(boxes(7,i) .eq. 0) goto 3000
        numkids=numkids+1
        kids(numkids)=boxes(7,i)
c 
        if(boxes(8,i) .eq. 0) goto 3000
        numkids=numkids+1
        kids(numkids)=boxes(8,i)
        goto 3000
 1800 continue
c 
c        if this box has no kids - look at the particles in it
c 
        j1=boxes(9,i)
        nj=boxes(10,i)
        do 2200 l=j1,j1+nj-1
c 
        jj=iz(l)
        x=z(1,jj)
        y=z(2,jj)
        d=a*x+b*y+c
c 
c       if this particle is close to the separator - act accordingly
c 
        if(dabs(d) .gt. eps) goto 2000
        numclose=numclose+1
        goto 2200
 2000 continue
c 
c       this particle is not close to the separator.
c       act accordingly
c 
        if(d .gt. 0) num2=num2+1
        if(d .le. 0) num1=num1+1
  
cccc        if(d .gt. 0) num1=num1+1
cccc        if(d .le. 0) num2=num2+1
 2200 continue
 3000 continue
        if(numkids .eq. 0) goto 4100
 4000 continue
 4100 continue
c 
c       allocate points in doubt ( those within eps
c       of the separator) between the two classes
c 
         if(numclose .eq. 0) return
        ndif=num2-num1
        if(ndif .lt. 0) ndif=-ndif
c 
        if(ndif .lt. numclose) goto 4200
        if(num1 .lt. num2) num1=num1+numclose
        if(num2 .le. num1) num2=num1+numclose
        return
 4200 continue
c 
        num1=n/2
        num2=n-num1
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
cccc        theta=pi/2
         theta=theta/3
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
         call prinf('in d2mallb, nbox=*',nbox,1)
         call prinf('in d2mallb, n=*',n,1)
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
c 
c 
        subroutine allplt3(iw,z,n,
     1    center0,size,a,b,c,mes)
        implicit real *8 (a-h,o-z)
        save
        real *8 z(2,1),xx(5),yy(5),center0(2)
        character *1 mes(1)
c 
c        plot the title
c 
         call prin2('in allplt3, a=*',a,1)
         call prin2('in allplt3, b=*',b,1)
         call prin2('in allplt3, c=*',c,1)
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
c       now, plot the separator
c 
        if(b .lt. a) goto 3400
        xx(1)=-100000
        xx(2)= 100000
        xx(1)=-10
        xx(2)= 10
        yy(1)=-(c+a*xx(1))/b
        yy(2)=-(c+a*xx(2))/b
 1410 FORMAT(2(2X,E15.7))
        write(iw,1410) xx(1),yy(1),xx(2),yy(2)
 3200 format(1x,'join ')
        write(iw,3200)
        goto 5000
 3400 continue
        yy(1)=-100000
        yy(2)= 100000
        yy(1)=-10
        yy(2)= 10
        xx(1)=-(c+b*yy(1))/a
        xx(2)=-(c+b*yy(2))/a
        write(iw,1410) xx(1),yy(1),xx(2),yy(2)
        write(iw,3200)
C1600 FORMAT(1X,'FRAME')
C1800 FORMAT(1X,'EXIT')
 1450 FORMAT(1X,'JOIN ')
c1500 FORMAT(1X,'PLOT')
 1600 FORMAT(1X,'#FRAME')
 1800 FORMAT(1X,'#EXIT')
cccc        WRITE(IW,1500)
        WRITE(IW,1600)
        WRITE(IW,1800)
 5000 continue
        return
        end
