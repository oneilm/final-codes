        implicit real *8 (a-h,o-z)
c
        call prini(6,13)
C
C       SET ALL PARAMETERS
C
         PRINT *, 'ENTER n'
         READ *,n
         CALL PRINF('n=*',n,1 )
c
        call testit(n)
       
c
        stop
        end
c
c
c
c
c
        subroutine testit(n)
        implicit real *8 (a-h,o-z)
        save
c
        real *8 xs(10000)
c
        integer *4 ifbombs(100000),w(1000 000)
c
c        create the test nodes
c
        done=1
        h=done/n

        xs(1)=-0.1
        do 1200 i=1,n
c
        xs(i+1)=(i-1)*h+h/2
 1200 continue
c
        call bidiag_rand(n,xs(2))
        call bidiag_rand(n,xs(2))
        call bidiag_rand(n,xs(2))
c
        call rsortanyn(xs,n,w)
c
        call prin2('xs as created *',xs,n)
c
c       test the subdivider routine
c
        a=-0.3
        b=1
        nmax=3
        maxchunks=10000

        do 1400 i=1,1000 000
c
        w(i)=777777
 1400 continue

        iplot=21
        lenw=1000 000
        call d1mstrcr(ier,n,xs,a,b,nmax,iplot,nn,
     1      w,lenw,keep)
c
c        conduct comprehensive testing of the obtained structure
c
        ifbombed=0

        do 4200 i=1,n
c
        call prinf('i=*',i,1)

        call prini(0,0)
        ipnt=i
c

        call all_lookup(ipnt,nn,n,ifbombs(i),w)

        if(ifbombs(i) .ne. 0) ifbombed=1
        call prini(6,13)

 4200 continue

        call prinf('ifbombs=*',ifbombs,n)
        call prinf('and ifbombed=*',ifbombed,1)
        call prinf('and keep=*',keep,1)


        return
        end
c
c
c
c
c
        subroutine all_lookup(ipnt,nboxes,n,ifbomb,w)
        implicit real *8 (a-h,o-z)
        integer *4 allnodes(100000),list(1000),ww(100 000),
     1       w(1),box(10)
        real *8 ab(3)
        save
c
c        find the childless box where this point lives
c
        do 1200 i=1,nboxes
c
        ibox=i

        call d1m_onebox(ibox,w,box,ab)

        if( (box(4) .gt. 0) .or. (box(5) .gt. 0) ) goto 1200
c
        i1=box(6)
        i2=i1+box(7)-1
c
        if( (ipnt .ge. i1) .and. (ipnt .le. i2) ) goto 1400
c
 1200 continue
 1400 continue
c
        call prinf('ipnt lives in childless chunk number*',ibox,1)
c
c        account for the nodes in the list 1 of the box ichnk
c
        ii=0
c
        itype=1
c
        call d1m_onelist(ibox,itype,w,list,nlist)

        call prinf('list=*',list,nlist)
c
        if(nlist .eq. 0) goto 1600
c
        do 1500 i=1,nlist
c
        jbox=list(i)
        call pnts_inbox2(jbox,w,allnodes(ii+1),np)
c
        ii=ii+np      
 1500 continue
c
 1600 continue
c
        call prinf('after processing list1, allnodes=*',
     1      allnodes,ii)
c
c        account for the nodes in the list 11 of the box ichnk
c
        itype=11
c
        call d1m_onelist(ibox,itype,w,list,nlist)
c
        call prinf('processing list11, list=*',list,nlist)
c
        if(nlist .eq. 0) goto 1800
c
        do 1700 i=1,nlist
c
        jbox=list(i)
        call pnts_inbox2(jbox,w,allnodes(ii+1),np)
c
        ii=ii+np      
 1700 continue

 1800 continue

        call prinf('after processing list11, allnodes=*',
     1      allnodes,ii)
c
c        account for the nodes in the list 3 of the box ichnk
c
        itype=3
c
        call d1m_onelist(ibox,itype,w,list,nlist)
c
        call prinf('processing list3, list=*',list,nlist)
c
        if(nlist .eq. 0) goto 2600
c
        do 2500 i=1,nlist
c
        jbox=list(i)
        call pnts_inbox2(jbox,w,allnodes(ii+1),np)
c
        ii=ii+np      
 2500 continue

 2600 continue

        call prinf('after processing list3, allnodes=*',
     1      allnodes,ii)
c
c        account for the nodes in all lists4 of all ancestors
c        of the box ibox
c
        icurr=ibox

        do 3600 ijk=1,100
c
c        obtain the list 2 of the box icurr
c
        itype=4

        call d1m_onelist(icurr,itype,w,list,nlist)
c
        call prinf('list=*',list,nlist)
c
        if(nlist .eq. 0) goto 3400
c
        do 3200 i=1,nlist
c
        jbox=list(i)
        call pnts_inbox2(jbox,w,allnodes(ii+1),np)
c
        ii=ii+np      
 3200 continue

 3400 continue

        call prinf('after processing list4, allnodes=*',
     1      allnodes,ii)


        call d1m_onebox(icurr,w,box,ab)
        icurr=box(3)
c
        if(icurr .eq. 1) goto 3800
c
 3600 continue
 3800 continue

c
c        account for the nodes in all lists2 of all ancestors
c        of the box ibox
c
        icurr=ibox

        do 5000 ijk=1,100
c
c        obtain the list 2 of the box icurr
c
        itype=2                        

        call d1m_onelist(icurr,itype,w,list,nlist)
c
        call prinf('list=*',list,nlist)
c
        if(nlist .eq. 0) goto 4600
c
        do 4500 i=1,nlist
c
        jbox=list(i)
        call pnts_inbox2(jbox,w,allnodes(ii+1),np)
c
        ii=ii+np      
 4500 continue

 4600 continue

        call prinf('after processing list4, allnodes=*',
     1      allnodes,ii)


        call d1m_onebox(icurr,w,box,ab)
        icurr=box(3)
c
        if(icurr .eq. 1) goto 5200
c
 5000 continue
 5200 continue

c

c        add to the list the nodes in the box ibox itself
c
        call d1m_onebox(ibox,w,box,ab)
c
        i1=box(6)
        n1=box(7)
c
        do 5400 i=1,n1
c
        ii=ii+1
c
        allnodes(ii)=i1+i-1
 5400 continue

        call prinf('after processing ibox itself, allnodes=*',
     1      allnodes,ii)

        call sortanyn(allnodes,ii,ww)

        call prinf('after sorting, allnodes=*',
     1      allnodes,ii)
c
c       check that allnodes consists of increasing numbers 
c       from 1 to n
c
        ifbomb=0
        do 6200 i=1,n
c
        j=i-allnodes(i)
        if(j .ne. 0) ifbomb=1
 6200 continue

        call prinf('and ifbomb=*',ifbomb,1)
c
        return
        end 
c
c
c
c
c
        subroutine pnts_inbox2(ibox,w,pnts,npnts)
        implicit real *8 (a-h,o-z)
        integer *4 pnts(1),w(1),box(10)
        real *8 ab(2)
        save
c
c        extract points from the box number ibox and return 
c        them to the user
c
        call d1m_onebox(ibox,w,box,ab)
        i1=box(6)
        npnts=box(7)
c
        do 1200 i=1,npnts
c
        pnts(i)=i1+i-1
 1200 continue
c
        return
        end
c 
c 
c 
c 
c 
        subroutine bidiag_rand(n,y)
        implicit real *8 (a-h,o-z)
        save
        dimension y(1)
        data ifcalled/0/
c 
c       generate pseudo-random numbers
c 
        if(ifcalled .ne. 0) goto 1100
c 
        done=1
        pi=atan(done)*4
        add=sqrt(7*done)
        x=1.1010101010101
        ifcalled=1
 1100 continue
c 
        do 1200 i=1,n
c 
        phi=x*100000*pi+add
        j=phi
        phi=(phi-j)
        x=phi
c 
        y(i)=x
 1200 continue
c 
        return
c 
c 
c 
c 
        entry bidiag_rand_reset(x7)
c 
        done=1
        pi=atan(done)*4
        add=sqrt(7*done)
        x=x7
        ifcalled=1
        return
        end
c
c
c 
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c         this is the end of the debugging code and the beginning
c         of the actual logic subroutines for the FMM on the line
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c         There are three user-callable subroutines in this file:
c         d1mstrcr, d1m_onebox, d1m_onelist. Following is a brief 
c         description of these three subroutines.
c
c   d1mstrcr - constructs the logical structure for the fully 
c        adaptive FMM in one dimension and stores it in the
c        array w. After that, the user can obtain the information 
c        about various boxes and lists in it by calling the 
c        subroutines d1m_onebox, d1m_onelist (see).
c
c   d1m_onebox - the standard data for the user-specified box number 
c        ibox (see the subroutine d1m_onebox for a complete description)
c
c   d1m_onelist - returns to the user a single user-specified list of 
c        the type itype (the possible types are 1,11,2,3,4) for the 
c        user-specified box number ibox.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
        subroutine d1mstrcr(ier,n,xs,a,b,nmax,iplot,nn,w,lenw,keep)
        implicit real *8 (a-h,o-z)
        real *8 xs(1)
        integer *4 w(1)
        save
c
c 
c        this subroutine constructs the logical structure for the
c        fully adaptive FMM in one dimension and stores it in the
c        array w. After that, the user can obtain the information 
c        about various boxes and lists in it by calling the 
c        subroutines d1m_onebox, d1m_onelist (see).
c
c     IMPORTANT!!!  PLEASE NOTE THAT THE NODES IN THE SIMULATION 
C        MUST BE IN INCREASING ORDER!!!
c 
c              Note on the list conventions.
c 
c    list 1 of the box ibox - the list of all boxes with which the
c           box ibox interacts directly that are either on the
c           same level as ibox or on the finer levels; the boxes 
c           on the coarser levels are not a part of list 1 of the
c           box ibox (see list 11 below). Obviously, list 1 is empty 
c           for any box that is not childless.
c 
c    list 11 of the box ibox - the list of all boxes with which the
c           box ibox interacts directly that are on coarser levels 
c           than ibox.
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
c  n - the number of elements (points) in the simulation
c  xs - the user-specified points on the line
c  a,b - the left and right ends of the structure respectively
c  nmax - the maximum number of points in a box on the finest level
c  lenw - the amount of memory in the array w (in integer *4 elements)
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
c  w - the array containing all tables describing boxes, lists, etc.
c         it contains a linked list (inter alia), and can only be 
c         accessed via the subroutines d1m_onebox, d1m_onelist (see).
c         The first keep integer *4 locations of this array should 
c         not be altered between the call to this subroutine and the
c         subsequent calls to the subroutines d1m_onebox, d1m_onelist.
c  keep - the amount of space in the array w (in integer *4 words)
c        that is occupied by various tables on exit from this
c        subroutine. This space should not be altered between the
c        call to this entry and subsequent calls to the subroutines 
c        d1m_onebox, d1m_onelist.
c
c       . . . construct all chunks on all levels
c
        ier=0
        ninire=2
c
        iaball=31
        laball=2*n*ninire*2+100

        laball=laball*4
c
        iifdone=iaball+laball
        lifdone=n*2+100
c
        ich=iifdone+lifdone
        lench=lenw-ich
c
        maxchunks=lench/10
c
        if(maxchunks .lt. 20) then
            ier=2048
            return
        endif
c
        call d1m_allsubdiv(jer,a,b,xs,n,nmax,maxchunks,
     1      w(iaball),w(ich),nn,w(iifdone) )
c
        if(jer .ne. 0) then
            ier=1024
            return
        endif
c
c        . . . compress memory
c
        laball=2*nn*ninire+2
c
        ich2=iaball+laball
        lch=nn*10+2
c
        call d1marrmove(w(ich),w(ich2),lch)
c
        ich=ich2
        w(1)=iaball
        w(2)=ich

        call prinf('after d1m_allsubdiv, w(ich)=*',w(ich),nn*10)
        call prin2('after d1m_allsubdiv, w(iaball)=*',w(iaball),nn*2)
c
        ilab1=ich+lch
        llab1=(n+nn)*2
c
        ilab2=ilab1+llab1
        llab2=(n+nn)*4
c
        ilab3=ilab2+llab2
        lenlab3=lenw-ilab3-10
c
        if(lenlab3 .lt. 1000) then
            ier=512
            return
        endif
c
        if(iplot .gt. 0) 
     1      call d1m_allplot(jer,iplot,a,b,xs,n,w(iaball),
     2      w(ich),nn,w(ilab1),w(ilab2),w(ilab3),lenlab3)
c
        if(jer .ne. 0) then
            ier=256
            return
        endif
c
c        find colleagues of all boxes
c
        iicolls=ich+lch
        licolls=nn*2+2
c
        call d1m_colleagues(w(ich),nn,w(iicolls),w(iaball) )
c
c        find lists 2 for all boxes
c
        ilists2=iicolls+licolls
        llists2=nn*3+10
c
        ltot=ilists2+llists2
c
        if(ltot .gt. lenw) then
            ier=128
            return
        endif
c
        w(3)=ilists2
c
        call d1m_lists2(w(ich),nn,w(iicolls),w(iaball),
     1      w(ilists2))
c
c        find the lists 1, 3 for all boxes
c
        ilists1=ilists2+llists2
        llists1=nn*2+2
c
        imaplist3=ilists1+llists1
        lmaplist3=nn*2+4
c
        iwork1=imaplist3+lmaplist3
        lwork1=300
c
        iwork2=iwork1+lwork1
        lwork2=300
c
        iwork=iwork2+lwork2
        lwork=nn*2+10
c
        iarr=iwork+lwork
c
        lenarr=lenw-iarr-10
c
        if(lenarr .lt. 20) then
            ier=64
            return
        endif
c
        w(4)=ilists1
        w(5)=imaplist3
c
        call d1m_lists31(jer,nn,w(ich),w(iicolls),w(iaball),
     1      w(imaplist3),w(iarr),lenarr,nstore,w(ilists1),
     2      w(iwork1),w(iwork2),w(iwork) )
c
        if(jer .ne. 0) then
            ier=64
            return
        endif
c
        larr=nstore+2
c
c        collect garbage
c
        iarr2=imaplist3+lmaplist3
        call d1marrmove(w(iarr),w(iarr2),larr)
        iarr=iarr2
        w(6)=iarr
c
c       construct the lists 11
c
        ilists11=iarr+larr
        llists11=nn*2+2
c
        ltot=ilists11+llists11
c
        if(ltot .gt. lenw) then
            ier=32
            return
        endif
c
        call d1m_lists11(w(ich),nn,w(ilists1),w(ilists11),
     1      w(iicolls) )
c
        w(7)=ilists11
c
        iaddr40=ilists11+llists11
        laddr40=nn+2
c
        w(8)=iaddr40
c
        iarr40=iaddr40+laddr40
        w(9)=iarr40
c
c        construct lists4
c
        lspace=lenw-iarr40-2

        call d1m_lists4(jer,nn,w(imaplist3),w(iarr),
     1      w(iarr40),lspace,w(iaddr40),nstored)
c
        if(jer .ne. 0) then 
            ier=8
            return
        endif
c
        keep=nstored+2+iarr40
c
        return
        end
c
c
c
c
c
        subroutine d1marrmove(a,b,n)
        integer *4 a(1),b(1)
        save
c
        do 1200 i=1,n
        b(i)=a(i)
 1200 continue
        return
        end
c
c
c
c
c
        subroutine d1m_onelist(ibox,itype,w,list,nlist)
        implicit real *8 (a-h,o-z)
        integer *4 w(1)
        save
c
c       this subroutine returns to the user a single user-specified 
c       list of the type itype (the possible types are 1,11,2,3,4) 
c       for the user-specified box number ibox.
c
c                Input parameters:
c
c  ibox - the box number for which the list is to be returned
c  itype - the type of the list to be returned
c  w - the array contraining all kinds of information about the
c       structure; hopefully, it has been created by a preceding
c       call to the subroutine d1mstrcr (see)
c
c                Output parameters:
c
c  list - the list of boxes making up the list of the type itype
c       of the box ibox
c  nlist - the number ofe elements in the array list. Please note
c       that nlist=0 is a possibel and frequently encountered 
c       result
c
c
c        . . . unpack the beginning of the array w
c
        iaball=w(1)
        ich=w(2)
        ilists2=w(3)
        ilists1=w(4)
        imaplist3=w(5)
        iarr=w(6)
        ilists11=w(7)
        iarr40=w(9)
        iaddr40=w(8)
c
c       retrieve the requested list
c
        call d1m_onelist0(ibox,itype,
     1      w(ilists1),w(ilists11),w(ilists2),w(imaplist3),
     2      w(iarr),w(iaddr40),w(iarr40),list,nlist)
c
        return
        end          
c
c
c
c
c
        subroutine d1m_onebox(ibox,w,box,ab)
        implicit real *8 (a-h,o-z)
        integer *4 w(1),box(7)
        real *8 ab(2)
        save
c
c       this subroutine returns to the user the standard data for
c       the user-specified box number ibox.
c
c                Input parameters:
c
c  ibox - the box number for which the data are to be returned
c  w - the array contraining all kinds of information about the
c       structure; hopefully, it has been created by a preceding
c       call to the subroutine d1mstrcr (see)
c
c                Output parameters:
c
c  box - integer array containing 7 elements, describing the box
c       number ibox, as follows:
c
c       box(1) - sequence number of the box number ibox
c       box(2) - level of subdivision of the box number ibox
c       box(3) - daddy of the box number ibox
c       box(4) - first sonny of the box number ibox
c       box(5) - second sonny of the box number ibox
c       box(6) - location of the first node in the box number ibox
c       box(7) - number of nodes in the box number ibox
c  ab - 2-element array contaning the ends of the box number ibox
c
c        . . . unpack the beginning of the array w
c
        iaball=w(1)
        ich=w(2)
c
c       retrieve the data about the box ibox
c
        call d1m_onebox0(ibox,w(ich),w(iaball),box,ab)
c
        return
        end          
c
c
c
c
c
        subroutine d1m_onebox0(ibox,ch,aball,box,ab)
        implicit real *8 (a-h,o-z)
        real *8 aball(2,1),ab(2)
        integer *4 ch(10,1),box(7)
        save
c
        do 1200 i=1,7
c
        box(i)=ch(i,ibox)
 1200 continue
c
        ab(1)=aball(1,ibox)
        ab(2)=aball(2,ibox)
c
        return
        end          
c
c
c
c
c
        subroutine d1m_onelist0(ibox,itype,
     1      lists1,lists11,lists2,maplist3,arr3,addr40,arr40,
     2      list,nlist)
        implicit real *8 (a-h,o-z)
        integer *4 lists1(2,1),lists2(3,1),maplist3(2,1),
     1      arr3(1),addr40(1),arr40(1),list(1),lists11(2,1)
        save
c
c        return list 1 of the box ibox
c
        if(itype .ne. 1) goto 1600
c
        nlist=0
        do 1200 i=1,2
c
        i1=lists1(i,ibox)
c
        if(i1 .gt. 0) then
            nlist=nlist+1
            list(nlist)=i1
        endif
c
 1200 continue
c
        return
c
 1600 continue
c
c        return list 11 of the box ibox
c
        if(itype .ne. 11) goto 2000
c
        nlist=0
        do 1800 i=1,2
c
        i1=lists11(i,ibox)
c
        if(i1 .gt. 0) then
            nlist=nlist+1
            list(nlist)=i1
        endif
c
 1800 continue
c
        return
c
 2000 continue
c
c        return list 2 of the box ibox
c
        if(itype .ne. 2) goto 3000
c
        nlist=0
        do 2200 i=1,3
c
        i1=lists2(i,ibox)
c
        if(i1 .gt. 0) then
            nlist=nlist+1
            list(nlist)=i1
        endif
c
 2200 continue
c
        return
c
 3000 continue
c
c        return list 3 of box ibox
c
        if(itype .ne. 3) goto 4000
c
        i1=maplist3(1,ibox)
        n1=maplist3(2,ibox)
c
        do 3200 i=1,n1
c
        list(i)=arr3(i1+i-1)
 3200 continue
c
        nlist=n1
        return
c
 4000 continue
c
c        return list 4 of this box
c
        if(itype .ne. 4) goto 5000
c
        call d1m_linkret(ibox,addr40,arr40,list,nlist)
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
        subroutine d1m_lists11(ch,nboxes,lists1,
     1      lists11,nums)
        implicit real *8 (a-h,o-z)
        integer *4 ch(10,1),lists1(2,1),lists11(2,1),nums(1)
        save
c
c        This subroutine returns to the user lists 11 for all 
c        boxes on all levels. Please note that list 11 of the
c        childless box ibox consists of boxes jbox that are
c        LARGER than ibox and are touching it; for a box with
c        kids, list 11 is empty. I STRONGLY SUSPECT THAT THE 
C        ONLY POSSIBLE SIZES FOR THE LIST 11 ARE 0 AND 1; 
C        AT THIS TIME, I SEE NOT QUITE SURE THAT IT CAN NOT 
C        CONTAIN TWO ELEMENTS (02.04.04).
C
C                     Input parameters:
c
c  ch - an integer array dimensioned (10,nn). each 10-element
c        column describes one chunk, as follows:
c 
c       ch(1,i) - sequence number of the i-th chunk
c       ch(2,i) - level of subdivision of the i-th chunk
c       ch(3,i) - daddy of the i-th chunk
c       ch(4,i) - first sonny of the i-th chunk
c       ch(5,i) - second sonny of the i-th chunk
c       ch(6,i) - location of the first node in the i-th chunk
c       ch(7,i) - number of nodes in the i-th chunk
c       ch(8,i)-ch(10,i) - not used at this time (3.02.04)
c  nboxes - the total number of boxes in the structure
c  lists1 - integer array dimensioned (2,nboxes) containing lists
c       1 for all boxes on all levels
c
c                     Output parameters:
c
c  lists11 - array containing lists 11 for all boxes on all levels
c
c                     Work arrays:
c
c  nums - must be at least nboxes*2 integer *4 elements long
c
        do 1200 i=1,nboxes
c
        lists11(1,i)=-7
        lists11(2,i)=-7
c
        nums(i)=0
 1200 continue
c
c        scan lists 1 for all boxes. Whenever a smaller box is 
c        found in such a list of the box ibox, enter ibox in the
c        appropriate location in array lists11
c
        do 2000 ibox=1,nboxes
c
        do 1600 i=1,2
c
        i1=lists1(i,ibox)
c
        if(i1 .le. 0) goto 1600
c
        levbox=ch(2,ibox)
        levi1=ch(2,i1)
c
        if(levi1 .le. levbox) goto 1600
c
        jj=nums(i1)+1
        nums(i1)=jj
        lists11(jj,i1)=ibox
c
 1600 continue
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
        subroutine d1m_lists4(ier,nboxes,maplist3,arr,
     1      arr40,lspace,addr40,nstored)
        implicit real *8 (a-h,o-z)
        integer *4 maplist3(2,1),arr(1),arr40(2,1),addr40(1)
        save
c
c        This subroutine constructs lists 4 for all boxes on all 
c        levels. The obtained boxes are stored in a linked list
c        managed by the subroutine d1m_linkinit (see).
c
c                      Input parameters:
c
c  
c  nboxes - the total number of boxes in the structure
c  maplist3 - the map of the array arr containing lists 3 for all
c       boxes on all levels (see the "please note" explanation 
c       above
c  arr - the storage area in which the lists 3 for all boxes on
c       all levels are stored, indexed by the map maplist3 (see
c       explanation above)
c  lspace - the amount of space provided (in integer *4 words) in 
c       the array arr40 to be used for storing the linked list.
c
c
c                     Output parameters:
c
c  ier - error return code
c  arr40 - the linked list containing lists 4 for all boxes on 
c       all levels
c  addr40 - integer array of length nboxes indexing the linked
c       list in array arr40
c  nstored - the amount of space (in integer *4 elements) used 
c       in the array arr40
c
c        . . . initialize the linked list storage-retrieval facility
c              for the lists 4
c
        ier=0
        call d1m_linkinit(addr40,nboxes,lspace)        
c
        do 2000 i=1,nboxes
c
        i3=maplist3(1,i)
        n3=maplist3(2,i)
c
        if(i3 .le. 0) goto 2000
c
        do 1600 j=1,n3
c
        jj=arr(i3+j-1)
        call d1m_linksto(ier,jj,i,addr40,arr40,larr)

        nstored=larr
c
 1600 continue
 2000 continue
c
        return
        end
c
c
c
c
c
        subroutine d1m_lists31(ier,nboxes,ch,icolls,aball,
     1      maplist3,arr,lenarr,nstore,lists1,l3,list3,work)
        implicit real *8 (a-h,o-z)
        real *8 aball(2,1)
        integer *4 icolls(2,1),ch(10,1),
     1      maplist3(2,1),arr(1),list1(10),list3(1),
     2      l1(10),l3(1),lists1(2,1),work(1)
        save
c
c        This subroutine constructs lists 3,1 for all boxes 
c        on all levels. Please note that list1 for the box ibox
c        consists of all boxes with which the box number ibox 
c        interacts directly, except for those boxes that are
c        bigger than box number ibox. List 3 of the box ibox 
c        consists of boxes jbox that are smaller than ibox, 
c        from which ibox is not separated, but that are separated 
c        from ibox.
c
c        PLEASE NOTE THAT THE LENGTHS OF LISTS 3 CAN NOT BE
C        KNOWN APRIORI, AND SOME OF THEM THINGS CAN BE FAIRLY
C        LONG. THUS, THEY ARE STORED IN AN INDEXED STORAGE ARRAY,
C        AS FOLLOWS:
C
C        For any i \in [1,nboxes], maplist(1,i) is the location
c        in array arr of the first element of the list 3 for the
c        box number i; maplist(2,i) is the number of elements 
c        (boxes) in the list 3 of the box i.
c
c                     Input parameters:
c
c  nboxes - the total number of boxes in the structure
c  ch - an integer array dimensioned (10,nn). each 10-element
c        column describes one chunk, as follows:
c 
c       ch(1,i) - sequence number of the i-th chunk
c       ch(2,i) - level of subdivision of the i-th chunk
c       ch(3,i) - daddy of the i-th chunk
c       ch(4,i) - first sonny of the i-th chunk
c       ch(5,i) - second sonny of the i-th chunk
c       ch(6,i) - location of the first node in the i-th chunk
c       ch(7,i) - number of nodes in the i-th chunk
c       ch(8,i)-ch(10,i) - not used at this time (3.02.04)
c  icolls - the array dimensioned (2,nboxes) containing colleagues
c       of all boxes in the structure (please note that the colleagues
c       of ibox are boxes of the same size as ibox, touching ibox).
c  aball - the array dimensioned (2,nn) containing ends of all
c        chunks in the structure (see array ch below)
c
c                     Output parameters:
c
c  maplist3 - the map of the array arr containing lists 3 for all
c       boxes on all levels (see the "please note" explanation 
c       above
c  arr - the storage area in which the lists 3 for all boxes on
c       all levels are stored, indexed by the map maplist3 (see
c       explanation above)
c  nstored - the number of integer *4 elements in array arr that
c       has been used for storage
c  lists1 - integer array dimensioned (2,nboxes) containing lists
c       1 for all boxes on all levels
c
c                     Work arrays:
c
c  l3,list3 - must be 200 integer *4 elements each
c  work - must be at least nboxes*2 integer *4 elements
c
c
        ier=0
        istore=1
c
        do 3000 ibox=1,nboxes
c
c        if this box has kids - its lists 1, 3 are empty.
c        act accordingly
c
        if( (ch(4,ibox) .gt. 0) .or. (ch(5,ibox) .gt. 0) ) then
            maplist3(1,ibox)=-1        
            maplist3(2,ibox)=0        
c
            lists1(1,ibox)=-7
            lists1(2,ibox)=-7
c
            goto 3000
        endif
c
c        this box is childless. find its colleagues and 
c        lists 1, 3
c
        nlist1=0
        nlist3=0
c
        do 1800 i=1,2
c
        icoll=icolls(i,ibox)
        if(icoll .le. 0) goto 1800
c
c        find the elements of lists 1, 3 contributed by this 
c        colleague
c
        nl1=0
        nl3=0
        call d1m_list31_0(ibox,icoll,ch,aball,
     1      l1,l3,nl1,nl3,work)
c
        do 1200 j=1,nl1
c
        list1(nlist1+j)=l1(j)
 1200 continue
c
        nlist1=nlist1+nl1
c
        do 1400 j=1,nl3
c
        list3(nlist3+j)=l3(j)
 1400 continue
c
        nlist3=nlist3+nl3
c
 1800 continue
c
c        store the obtained lists 1, 3 in the array lists1
c
        lists1(1,ibox)=-7
        lists1(2,ibox)=-7
c
        if(nlist1 .ge. 1) lists1(1,ibox)=list1(1)
        if(nlist1 .eq. 2) lists1(2,ibox)=list1(2)
c
c        store the obtained lists 3 in the storage area,
c        and enter their addresses in the arrays maplist3
c
        maplist3(1,ibox)=-1
        maplist3(2,ibox)=-1
c
        if(nlist3 .le. 0) goto 2600
c
        maplist3(1,ibox)=istore
        maplist3(2,ibox)=nlist3
c
        if(istore+nlist3 .gt. lenarr) then
            ier=16
            return
        endif
c
        do 2400 j=1,nlist3
c
        arr(istore+j-1)=list3(j)
 2400 continue
c
        istore=istore+nlist3
c
 2600 continue
c
 3000 continue
c
        nstore=istore+2
c
        return
        end
c
c
c
c
c
        subroutine d1m_list31_0(ibox,icoll,ch,aball,
     1      list1,list3,nlist1,nlist3,index_arr)
        implicit real *8 (a-h,o-z)
        real *8 aball(2,1)
        integer *4 ch(10,1),index_arr(2,1),list1(1),
     1      list3(1)
        save
c
c        This subroutine recursively subdivides the colleague icoll
c        of the childless box ibox, producing the part of lists 1 
c        and 3 of the box ibox resulting from the colleague icoll
c
c                Input parameters:
c
c  ibox - the box for which the lists 1, 3 are to be constructed;
c        must be childless
c  icoll - the colleague of the box ibox (one of possible two) 
c        which is to be subdivided into lists 1, 3
c  ch - an integer array dimensioned (10,nn). each 10-element
c        column describes one chunk, as follows:
c 
c       ch(1,i) - sequence number of the i-th chunk
c       ch(2,i) - level of subdivision of the i-th chunk
c       ch(3,i) - daddy of the i-th chunk
c       ch(4,i) - first sonny of the i-th chunk
c       ch(5,i) - second sonny of the i-th chunk
c       ch(6,i) - location of the first node in the i-th chunk
c       ch(7,i) - number of nodes in the i-th chunk
c       ch(8,i)-ch(10,i) - not used at this time (3.02.04)
c  icolls - the array dimensioned (2,nboxes) containing colleagues
c       of all boxes in the structure (please note that the colleagues
c       of ibox are boxes of the same size as ibox, touching ibox).
c  aball - the array dimensioned (2,nn) containing ends of all
c        chunks in the structure 
c
c                Output parameters:
c
c  list1 - the elements (boxes) of the list 1 of the box ibox
c        resulting from his colleague icoll
c  list3 - the elements (boxes) of the list 3 of the box ibox
c        resulting from his colleague icoll
c  nlist1 - the number of elements in list 1 (can be 0 or 1)
c  nlist3 - the number of elements in list 3 (can be 0 or more)
c
c                Work arrays:
c
c  index_arr - must be at least 2*nboxes integer *4 elements long
c
c        . . . initialize the process of subdividing the colleague
c
        ii=1
        nn=ii
        nlist1=0
        nlist3=0
        index_arr(1,1)=icoll
        index_arr(2,1)=-1
        
        do 2000 ijk=1,100
c
        ifacted=0
        do 1800 i=1,nn
c
c        if this box already has been processed - skip
c
        if(index_arr(2,i) .eq. 1) goto 1800
c
        ifacted=1
        ib=index_arr(1,i)
c
c        if this box does not touch ibox 
c        - enter it in list 3 of ibox
c
        eps=1.0d-10
        d1=abs(aball(1,ibox)-aball(2,ib))
        d2=abs(aball(2,ibox)-aball(1,ib))
        d3=aball(2,ibox)-aball(1,ibox)
c
        if( (d1 .gt. d3*eps) .and. (d2 .gt. d3*eps) ) then
c
            nlist3=nlist3+1
            list3(nlist3)=ib
            index_arr(2,i)=1
            goto 1800
        endif 
c
        i1=ch(4,ib)
        i2=ch(5,ib)
c
        if( (i1 .le. 0) .and. (i2 .le. 0) ) goto 1400
c
        if(i1 .gt. 0) then
            ii=ii+1
            index_arr(1,ii)=i1
            index_arr(2,ii)=-1
        endif
c
        if(i2 .gt. 0) then
            ii=ii+1
            index_arr(1,ii)=i2
            index_arr(2,ii)=-1
        endif
c
        index_arr(2,i)=1                    
        goto 1800
c
 1400 continue
c
c       this box is childless, and touches the box ibox.
c       assign it to list 1
c
        nlist1=nlist1+1
        list1(nlist1)=ib
        index_arr(2,i)=1
        goto 1800
c
 1800 continue
c
        if(ifacted .eq. 0) goto 2200
c
        nn=ii
 2000 continue
 2200 continue
        return
        end          
c
c
c
c
c
        subroutine d1m_lists2(ch,nn,icolls,aball,lists2)
        implicit real *8 (a-h,o-z)
        real *8 aball(2,1)
        integer *4 ch(10,1),icolls(2,1),kids(100),
     1      lists2(3,1)
        save
c
c        construct lists 2 for all boxes
c
        lists2(1,1)=-7
        lists2(2,1)=-7
        lists2(3,1)=-7
c
        do 2000 i=2,nn
c
        lists2(1,i)=-7
        lists2(2,i)=-7
        lists2(3,i)=-7
c
c        find the daddy of this box and the daddy's colleagues
c
        idad=ch(3,i)
        idadcol1=icolls(1,idad)
        idadcol2=icolls(2,idad)
c
c       find the kids of the daddy and daddy's colleagues
c
        nkids=0

        kid1=ch(4,idad)                
        kid2=ch(5,idad)                
c
        if( (kid1 .ne. i) .and. (kid1 .gt. 0) ) then
            nkids=nkids+1
            kids(nkids)=kid1
        endif
c
        if( (kid2 .ne. i) .and. (kid2 .gt. 0) ) then
            nkids=nkids+1
            kids(nkids)=kid2
        endif
c
c       . . . look at the first colleague
c
        if(idadcol1 .le. 0) goto 1400
c
        kid1=ch(4,idadcol1)        
        kid2=ch(5,idadcol1)        
c
        if( (kid1 .ne. i) .and. (kid1 .gt. 0) ) then
            nkids=nkids+1
            kids(nkids)=kid1
        endif
c
        if( (kid2 .ne. i) .and. (kid2 .gt. 0) ) then
            nkids=nkids+1
            kids(nkids)=kid2
        endif
 1400 continue
c
c       . . . look at the second colleague
c
        if(idadcol2 .le. 0) goto 1600
c
        kid1=ch(4,idadcol2)        
        kid2=ch(5,idadcol2)        
c
        if( (kid1 .ne. i) .and. (kid1 .gt. 0) ) then
            nkids=nkids+1
            kids(nkids)=kid1
        endif
c
        if( (kid2 .ne. i) .and. (kid2 .gt. 0) ) then
            nkids=nkids+1
            kids(nkids)=kid2
        endif
 1600 continue
c
c       look at the kids of daddy's and his coleagues, and 
c       sort them into lists 2, 5 for the box number i 
c
        d3=aball(2,i)-aball(1,i)
        eps=1.0d-10
        jj=0
        jj5=0
        do 1800 j=1,nkids
c
        kid=kids(j)
        d1=abs(aball(1,kid)-aball(2,i))
        d2=abs(aball(2,kid)-aball(1,i))
c
        if( (d1 .gt. eps*d3) .and. (d2 .gt. eps*d3) ) then
            jj=jj+1
            lists2(jj,i)=kid
        endif
c
 1800 continue
c
 2000 continue

        return
        end          
c
c
c
c
c
        subroutine d1m_colleagues(ch,nn,icolls,aball)
        implicit real *8 (a-h,o-z)
        real *8 aball(2,1)
        integer *4 ch(10,1),icolls(2,1),kids(100)
        save
c
c        construct lists of colleagues for all boxes on all
c        levels
c
        icolls(1,1)=-7
        icolls(2,1)=-7
c
        do 2000 i=2,nn
c
        icolls(1,i)=-7
        icolls(2,i)=-7
c
c        find the daddy of this box and the daddy's colleagues
c
        idad=ch(3,i)
        idadcol1=icolls(1,idad)
        idadcol2=icolls(2,idad)
c
c       find the kids of the daddy and his colleagues
c
c       . . . look at the daddy
c
        nkids=0
        kid1=ch(4,idad)                
        kid2=ch(5,idad)                
c
        if( (kid1 .ne. i) .and. (kid1 .gt. 0) ) then
            nkids=nkids+1
            kids(nkids)=kid1
        endif
c
        if( (kid2 .ne. i) .and. (kid2 .gt. 0) ) then
            nkids=nkids+1
            kids(nkids)=kid2
        endif
c
c       . . . look at the first colleague
c
        if(idadcol1 .le. 0) goto 1400
c
        kid1=ch(4,idadcol1)        
        kid2=ch(5,idadcol1)        
c
        if( (kid1 .ne. i) .and. (kid1 .gt. 0) ) then
            nkids=nkids+1
            kids(nkids)=kid1
        endif
c
        if( (kid2 .ne. i) .and. (kid2 .gt. 0) ) then
            nkids=nkids+1
            kids(nkids)=kid2
        endif
 1400 continue
c
c       . . . look at the second colleague
c
        if(idadcol2 .le. 0) goto 1600
c
        kid1=ch(4,idadcol2)        
        kid2=ch(5,idadcol2)        
c
        if( (kid1 .ne. i) .and. (kid1 .gt. 0) ) then
            nkids=nkids+1
            kids(nkids)=kid1
        endif
c
        if( (kid2 .ne. i) .and. (kid2 .gt. 0) ) then
            nkids=nkids+1
            kids(nkids)=kid2
        endif
 1600 continue
c
c       look at the kids of daddy and his coleagues, and 
c       extract the colleagues of the box number i from
c       among them
c
        d3=aball(2,i)-aball(1,i)
        eps=1.0d-10
        jj=0
        do 1800 j=1,nkids
c
        kid=kids(j)
        d1=abs(aball(1,kid)-aball(2,i))
        d2=abs(aball(2,kid)-aball(1,i))
c
        if( (d1 .lt. eps*d3) .or. (d2 .lt. eps*d3) ) then
            jj=jj+1
            icolls(jj,i)=kid
        endif
c
 1800 continue
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
        subroutine d1m_allplot(ier,iw,a,b,xs,n,aball,ch,nn,
     1     labs,zslabs,w,lenw)
        implicit real *8 (a-h,o-z)
        real *8 xs(1),aball(2,1),w(1),zslabs(2,1)
        integer *4 ch(10,1),labs(1)
        save
c
        ier=0
        ll=0
c
        call trplopen(w)
c
c       plot the ends of the big interval
c
        rlen=b-a

        x1=a
        y1=-0.5*rlen
c
        x2=a
        y2=0.5*rlen
c
        ll=ll+1
        call trplline(x1,y1,x2,y2,w)
c
        x1=b
        y1=-0.5*rlen
c
        x2=b
        y2=0.5*rlen
c
        ll=ll+1
        call trplline(x1,y1,x2,y2,w)
c
c        plot the points
c
        do 1400 i=1,n
c
        zero=0
        ll=ll+1
        if(ll*5 .gt. lenw) then
            ier=16
            return
        endif
c
cccc        call trplpnt(xs(i),zero,w)
        call trpldot(xs(i),zero,w)
 1400 continue
c
c        plot the boundaries of all chunks
c
        do 2000 i=1,nn
c
        lev=ch(2,i)
        x1=aball(1,i)
        x2=aball(1,i)
c
        y1=0.5/2**lev *rlen
        y2=-y1
c
        ll=ll+1
        if(ll*5 .gt. lenw) then
            ier=16
            return
        endif
c
        call trplline(x1,y1,x2,y2,w)
c
        lev=ch(2,i)
        x1=aball(2,i)
        x2=aball(2,i)
c
        y1=0.5/2**lev *rlen
        y2=-y1
c
        ll=ll+1
        if(ll*5 .gt. lenw) then
            ier=16
            return
        endif
c
        call trplline(x1,y1,x2,y2,w)
c
        zslabs(1,i)=(aball(1,i)+aball(2,i))/2
        zslabs(2,i)=y1
c
        labs(i)=i
 2000 continue
c
        call zqua_labels_int(iw,zslabs,labs,nn)
        call trplwrt(iw,w,'structure*')
c
        return
        end
c
c
c
c
c
        subroutine d1m_allsubdiv(ier,a,b,xs,n,nmax,maxchunks,
     1      aball,ch,nn,ifdone)
        implicit real *8 (a-h,o-z)
        real *8 xs(1),aball(2,1),abson1(2),abson2(2)
        integer *4 ch(10,1),ifdone(1)
        save
c
c        For a collection xs of n points on the interval [a,b],
c        this subroutine constructs the binary tree.
c
c              input parameters:
c 
c  
c  a,b - the ends of the interval on which the nodes live
c  xs - the set of points in on the interval [a,b]
c  n - the number of elements in xs
c  nmax - the maximum number of points permitted in a chunk on
c        the finest level. In other words, a chunk will be further
c        subdivided if it contains more than nmax points.
c  maxchunks - the maximum total number of boxes the subroutine
c        is permitted to create. If the points xs are such that
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
c  aball - the array dimensioned (2,nn) containing ends of all
c        chunks in the structure (see array ch below)
c  ch - an integer array dimensioned (10,nn). each 10-element
c        column describes one chunk, as follows:
c 
c       ch(1,i) - sequence number of the i-th chunk
c       ch(2,i) - level of subdivision of the i-th chunk
c       ch(3,i) - daddy of the i-th chunk
c       ch(4,i) - first sonny of the i-th chunk
c       ch(5,i) - second sonny of the i-th chunk
c       ch(6,i) - location of the first node in the i-th chunk
c       ch(7,i) - number of nodes in the i-th chunk
c       ch(8,i)-ch(10,i) - not used at this time (2.28.04)
c
c    important warning: the array boxes has to be dimensioned
c                       at least (10,maxchunks)!! otherwise,
c                       the subroutine is likely to bomb, since
c                       it assumes that it has that much space!!!!
c  nn - the total number of boxes created
c
c                  work arrays:
c
c  ifdone - must be at least maxchunks integer *4 elements long
c
c
c        . . . enter in array chunks the data for the 
c              whole big interval
c
        ier=0
c
        ch(1,1)=1
        ch(2,1)=0
        ch(3,1)=-7
        ch(4,1)=-1
        ch(5,1)=-1
        ch(6,1)=1
        ch(7,1)=n
c
        aball(1,1)=a
        aball(2,1)=b
c
        do 1200 i=1,n
c
        ifdone(i)=0
 1200 continue
c
c        recursively subdivide the chunks
c
        nn=1
        ii=1
c
        do 3000 ijk=1,100
c
        if_enough=1
        do 2000 i=1,nn
c
c       if this chunk has already been processed (subdivided 
c       or not) - skip it
c
        if(ifdone(i) .eq. 1) goto 2000
c
c       if this chunk contains fewer than nmax nodes 
c       declare it processed and don't bother with it
c       anymore
c
        kk=ch(7,i)
        if(kk .le. nmax) then
            ifdone(i)=1
            goto 2000
        endif
c
c        this chunk needs to be subdivided. act accordingly
c
        idad=ch(6,i)
        ndad=ch(7,i)
c
        call d1m_onesubdiv(aball(1,i),idad,ndad,xs,
     1      abson1,ison1,nson1,abson2,ison2,nson2)
c
c        create the sonnies as needed
c
        if(ison1 .gt. 0) then
            ii=ii+1
            if(ii .gt. maxchunks) then
                ier=4
                return
            endif
c
            aball(1,ii)=abson1(1)                    
            aball(2,ii)=abson1(2)                    
c
            ch(1,ii)=ii
            ch(2,ii)=ch(2,i)+1
            ch(3,ii)=i
            ch(4,ii)=-1
            ch(5,ii)=-1
            ch(6,ii)=ison1
            ch(7,ii)=nson1
c
            ch(4,i)=ii
        if_enough=0
        endif
c
        if(ison2 .gt. 0) then
            ii=ii+1
            if(ii .gt. maxchunks) then
                ier=4
                return
            endif
c
            aball(1,ii)=abson2(1)                    
            aball(2,ii)=abson2(2)                    
c
            ch(1,ii)=ii
            ch(2,ii)=ch(2,i)+1
            ch(3,ii)=i
            ch(4,ii)=-1
            ch(5,ii)=-1
            ch(6,ii)=ison2
            ch(7,ii)=nson2
c
            ch(5,i)=ii
        if_enough=0
        endif
c
        ifdone(i)=1
c
 2000 continue
c
        nn=ii
c
        if(if_enough .eq. 1) goto 3200
 3000 continue
 3200 continue
c
        return
        end
c
c
c
c
c
        subroutine d1m_onesubdiv(abdad,idad,ndad,xs,
     1      abson1,ison1,nson1,abson2,ison2,nson2)
c
        implicit real *8 (a-h,o-z)
        real *8 xs(1),abdad(2),abson1(2),abson2(2)
        save
c
c        subdivide the daddy chunk into two sonnies
c
        dd=(abdad(1)+abdad(2))/2
        abson1(1)=abdad(1)
        abson1(2)=dd
c
        abson2(1)=dd
        abson2(2)=abdad(2)
c
        ii=0
        do 1200 i=1,ndad
c
        j=idad+i-1
        ii=i
        if(xs(j) .gt. dd) goto 1400
c
 1200 continue
c
        ii=ndad+1
c
 1400 continue
c
c       if one of the two sonnies does not exist - act accordingly
c     
        if(ii .eq. ndad+1) then
            ison2=-7    
            nson2=-7    
c 
            ison1=idad
            nson1=ndad
c 
            return
        endif
c
c
        if(ii .eq. 1) then
            ison1=-7    
            nson1=-7    
c 
            ison2=idad
            nson2=ndad
c 
            return
        endif
c
c        both sonnies DO exist. construct them things
c
        ison1=idad
        nson1=ii-1
c
        ison2=idad+ii-1
        nson2=ndad-ii+1
c
        return
        end

c
c            Interpretation of various entries in array chunks
c
c
c       ch(1,i) - sequence number of the i-th chunk
c       ch(2,i) - level of subdivision of the i-th chunk
c       ch(3,i) - daddy of the i-th chunk
c       ch(4,i) - first sonny of the i-th chunk
c       ch(5,i) - second sonny of the i-th chunk
c       ch(6,i) - location of the first node in the i-th chunk
c       ch(7,i) - number of nodes in the i-th chunk
c
c
c
c
c
c
        subroutine d1m_linkinit(addr,n7,lspace7)        
        implicit real *8 (a-h,o-z)
        integer *4 addr(1),arr(2,1),vals(1)
        save
c
c        This entry initializes a linked list storage-retrieval 
c        facility. For n7 boxes, the facility stores one integer *4
c        number per storage call. A retrieval call retrieves all 
c        such numbers stored for the specified box by all preceding
c        storage calls.
c
c                   Input parameters:
c
c  n7 - the number of types of boxes
c  lspace7 - the amount of space (in integer *4 elements) supplied
c        by the user in the storage array arr (see entries d1m_linksto,
c        d1m_linkret of this subroutine for the description of the 
c        parameter arr)
c
c                   Output parameters:
c
c  addr - the initialized integer *4 array of length n7
c
        n=n7
c
        do 1200 i=1,n
c
        addr(i)=-7
 1200 continue
c
        nstored=0
        lspace=lspace7
c
        return
c
c
c
c
        entry d1m_linksto(ier,ibox,ival,addr,arr,larr)
c
c        This entry stores in the arrays addr, arr the user-supplied
c        integer *4 number associated with the box ibox
c
c                   Input parameters:
c
c  ibox - the box (one of n) with which the number to be stored
c        is associated
c  ival - the value to be stored
c  addr - the pointer array: the element addr(i) contains the address 
c        in array arr of the last location associated with the box
c        number i
c  arr - the main storage area. It is formatted as arr(2,???), with
c        the element arr(1,i) containing the value stored, and 
c        the element arr(2,i) containing the address in array arr
c        of the preceding element associated with the same box as 
c        the one with which the current element is associated
c
c                   Output parameters:
c
c  ier - error return code. 
c     ier=0 means succesful conclusion
c     ier=4 means that the amount of space allocated by the user in
c        array arr has been exceeded. Nothing has been stored
c  
c  addr - see above
c  arr - see above
c  larr - the amount of space in array arr that has been consumed 
c        by the subroutine
c
c       store the value ival assiciated with the box ibox
c        
        ier=0
        nstored=nstored+1
c
        if(nstored*2+2 .gt. lspace) then
            ier=4
            return
        endif
c
        arr(1,nstored)=ival
        arr(2,nstored)=addr(ibox)
        addr(ibox)=nstored
c
        larr=nstored*2
c
        return
c
c
c
c
        entry d1m_linkret(ibox,addr,arr,vals,nvals)
c
c        This entry retrieves from the arrays addr, arr the array
c        of integer *4 numbers associated with the box ibox;
c        hopefully, these numbers have been stored there by preceding
c        calls to the entry d1m_linksto
c
c                   Input parameters:
c
c  ibox - the box (one of n) with which the numbers to be retrieved
c        are associated
c  addr - the pointer array: the element addr(i) contains the address 
c        in array arr of the last location associated with the box
c        number i
c  arr - the main storage area. It is formatted as arr(2,???), with
c        the element arr(1,i) containing the value stored, and 
c        the element arr(2,i) containing the address in array arr
c        of the preceding element associated with the same box as 
c        the one with which the current element is associated
c
c                   Output parameters:
c
c  vals - the list of all values associated with the box ibox that
c        have been stored in the facility
c  nvals - the number of elements returned in the array vals
c
c       . . . retrieve all values associated with the box ibox
c
        nvals=0
        ii=addr(ibox)
        if(ii .lt. 0) return
c
        do 1400 i=1,1000000
c
        nvals=i
        vals(nvals)=arr(1,ii)
c
        ii=arr(2,ii)
        if(ii .lt. 0) return
 1400 continue
c
        return
        end
