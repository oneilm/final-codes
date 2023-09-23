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
        external matfun2
c
        real *8 xs(100 000),charge(100 000),fld(500 000),w2(3000 000),
cccc     1      vals(100 000),diffs(100 000),ww(3000 000),wsave(1500 000)
     1      vals(100 000),diffs(100 000),wsave(1500 000)
c
cccc        integer *4 w(20 000 000)
        integer *4 w(60 000 000)
ccc        integer *4 w(7 960 000)
ccc        integer *4 w(7 940 000)



cccc        integer *4 w(5 000 000)
c
c        create the test nodes
c
        done=1
        h=done/n

        xs(1)=-0.4
        xs(2)=-0.1
        do 1200 i=2,n
c
        xs(i+1)=(i-1)*h+h/2
 1200 continue
c
        call rirand(n,xs(2))
        call rirand(n,xs(2))
        call rirand(n,xs(2))
c
        call rsortanyn(xs,n,w)
c

        do 1300 j=1,n
        if(j .gt. n/2) xs(j)=xs(j)+2
cc        if(j .gt. 10) xs(j)=xs(j)+2
 1300 continue


cccc        call prin2('xs as created *',xs,n)
        call prin2('xs as created *',xs,100)


        call matfuni2(xs,n)


        call matfun2_set(1)

c
c       test the subdivider routine
c
        a=-0.3
        b=1
        nmax=35
        nmax=45
        nmax=50
        nmax=60
cccc        nmax=40

        do 1400 i=1,1000 000
c
        w(i)=777777
 1400 continue

        iplot=21
        nts=20
cccc        nts=30
        eps=1.0d-13
cccc        eps=1.0d-10
cccc        eps=1.0d-6

        ifrelskel=1

c

        lenw=60 000 000
        lenw=30 000 000


cccc        call prini(0,0)

        call hudson_allchunks(ier,xs,n,matfun2,nmax,nts,
     1      eps,ifrelskel,iplot,nn,w,lenw,keep,lused)
c
       
        call prinf('after hudson_allchunks, ier=*',ier,1)


        do 5200 i=1,1000
        charge(i)=0
cccc        charge(i)=1
 5200 continue
c
        charge(50)=1
cc        charge(31)=1
c
        call bidiag_rand(99000,charge)

        charge(1)=0
        charge(2)=0

        do 5300 j=1,22
cc        charge(i)=0
 5300 continue

c
c       evaluate the potential via the subroutine
c
        t1=clotatim()

        lenw2=3000 000
        do 5400 ijk=1,100
c
        call hudson_eval(w,charge,vals,w2,lenw2,lusedw2)

 5400 continue
        t2=clotatim()

cccc        call prini(6,13)


        call prin2('after hudson_eval, vals=*',vals,n)

cccc        call comp_fld(xs,charge,n,amatr,fld)

        t3=clotatim()


        call comp_fld(xs,charge,n,fld)

        iw31=33
cc        call res_rd(iw31,fld,n)
cccc        call res_wrt(iw31,fld,n)



        t4=clotatim()

cccc        call prin2('after comp_fld, fld=*',fld,n)

        errmax=0
        d2=0
        d3=0
        imax=0
        do 6400 i=1,n
c
        diffs(i)=vals(i)-fld(i)
        d=abs(diffs(i))
        if(d .gt. errmax) imax=i
        if(d .gt. errmax) errmax=d
c
        d2=d2+vals(i)**2
        d3=d3+d**2
 6400 continue

        call prin2('and diffs=*',diffs,n)


        call prin2('and errmax=*',errmax,1)

        errrel=sqrt(d3/d2)
        call prin2('and errrel=*',errrel,1)
        call prinf('and imax=*',imax,1)
        call prin2('and vals(imax)=*',vals(imax),1)

cccc        stop


        call prin2('time for hudson_eval=*',t2-t1,1)
        call prin2('time for comp_fld=*',t4-t3,1)


        rat=100*(t4-t3)/(t2-t1)
        
        call prin2('and rat=*',rat,1)

        breakeven=n/rat

        call prin2('and breakeven=*',breakeven,1)


        call prinf('and by the way, keep=*',keep,1)
        call prinf('and keep/1000=*',keep/1000,1)

        call prinf('and lused=*',lused,1)
        call prinf('and lused/1000=*',lused/1000,1)
c
c        compare the speed of this thing with the speed of FFT
c
cccc        stop


        call DCFFTI(N,WSAVE)

        t6=clotatim()

        do 6600 ijk=1,100

        call DCFFTf(N,fld,WSAVE)


        do 6500 i=1,n*2
c
cccc        fld(i)=vals(i)
        fld(i)=i
 6500 continue

 6600 continue

        t7=clotatim()

        call prin2('and the time for fft=*',t7-t6,1)

cccc        call prin2('and fld=*',fld,n*2)

cc        stop       
c
c        compare the speed of this thing with the speed of hudson_matvec
c
        do 7200 i=1,1000 000
c
        wsave(i)=i
 7200 continue


        t8=clotatim()

        do 7400 i=1,100

        call hudson_matvec(wsave,1000,1000,vals,fld)


 7400 continue

        t9=clotatim()

 
        call prin2('and the time for 1000 by 1000 hudson_matvec=*',
     1      t9-t8,1)





        return
        end
c
c
c
c
c
        subroutine res_wrt(iw,a,n)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n)
c
        write (iw,1200) (a(i),i=1,n)
 1200 format(10x,e23.17)
c
        return
c
c
c
c
        entry res_rd(iw,a,n)
        read (iw,1200) (a(i),i=1,n)
        return
        end
c
c
c
c
c
        subroutine matfun2(x8,ix,y8,iy,f)
        implicit real *8 (a-h,o-z)
        save
        dimension xs7(1),xs(100 000)
c
        x=x8
        y=y8
c
        if(ifuse .eq. 1) then
            x=xs(ix)
            y=xs(iy)


        goto 1100
            call prin2('x8=*',x8,1)
            call prin2('y8=*',y8,1)
c
         
            call prin2('while x=*',x,1)
            call prin2('and y=*',y,1)

            call prinf('and ix=*',ix,1)
            call prinf('and iy=*',iy,1)

            call prin2(' *',y,0)
 1100 continue

        endif

cccc        call prinf('in  matfun2, ifuse=*',ifuse,1)

c
c        construct lota for this pair of chunks
c
ccc        if( abs(x-y) .lt. 1.0d-9) then
        if( abs(x-y) .lt. 1.0d-7) then
            f=0
            return
        endif
c
cccc        f=1/(x-y)*x**2*y
cccc        f=1/(x-y)*x*y
cccc        f=1/(x-y)*(x+y)

cc        f=1/(x-y)*exp(-y*10)


cccc        f=1/(x-y)*x**2
        f=1/(x-y)
cc        if(x .gt. y) f=(x-y)+x**2
ccc        if(x .gt. y) f=(x-y)

        if(x .gt. y) f=1	
cc        if(x .gt. y) f=-f

cccc        if(x .gt. y) f=-f

cccc        f=x


c
        return
c
c
c
c    
        entry matfuni2(xs7,n7)
c
        do 2200 i=1,n7
c
        xs(i)=xs7(i)
 2200 continue
c
        n=n7
        ifuse=0
        return
c
c
c
c
        entry matfun2_set(ifuse7)
c
        ifuse=ifuse7
        return
        end
c
c
c
c
c
        subroutine comp_fld(xs,charge,n,fld)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),fld(1),charge(1)
c
c        construct the box matrix of interactions
c
        do 1400 i=1,n
c
        dd=0

        do 1200 j=1,n
c


ccc        call matfun1(xs(i),i,xs(j),i,d)
ccccc        call matfun1(xs(j),j,xs(i),i,d)
        call matfun2(xs(j),j,xs(i),i,d)
c
        dd=dd+charge(j)*d
 1200 continue
c
        fld(i)=dd
 1400 continue
c
c        apply the matrix to the vector
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
        x=1.1010101010101d0
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c         This is the end of the debugging code, and the beginning of
c         the one-dimensional FMM code proper
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c
c
c         This file contains (on 06.04.04) two user-callable subroutines:
c         hudson_allchunks and hudson_eval. Following is a brief 
c         description of these two subroutines.
c
c   hudson_allchunks - for the user-specified nodes xs and several 
c        other parameters, constructs the box structure, incoming and 
c        outgoing skeletons on all levels, and all translation matrices
c        (LOTA, TATA, LOTA, etc.)for all boxes on all levels. The 
c        information created by this is expected to be used by the 
c        subroutine hudson_eval (see)
c
c  hudson_eval - evaluates the potential at the nodes xs of a "charge"
c        distribution at these same nodes. It is fervently hoped that 
c        the necessary precomputation has been performed via a preceding 
c        call to the subroutine hudson_allchunks (see)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
        subroutine hudson_allchunks(ier,xs,n,matfun,nmax,
     1      nts,eps,ifrelskel,iplot,nn,w,lenw,keep,lused)
        implicit real *8 (a-h,o-z)
        save
        real *8 xs(1)
        integer *4 w(1)
c
c
c        For the user-specified nodes xs and several other parameters,
c        this subroutine constructs the box structure, incoming and 
c        outgoing skeletons on all levels, and all translation matrices
c        (LOTA, TATA, LOTA, etc.)for all boxes on all levels. The 
c        information created by this is expected to be used by the 
c        subroutine hudson_eval (see)
c
c                  Input parameters:
c
c  xs - the used-specified nodes in the simulation; must be in increasing
c        order
c  n - the number of elements in array xs
c  matfun - the user-supplied subroutine evaluating the interactions
c        to be "accelerated"
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        The calling sequence of matfun is
c
c        subroutine matfun(x,ix,y,iy,f),
c
c        where 
c
c    x - the SOURCE node (its location in R^1)        (input)
c    ix - the sequence number of the source node      (input)
c    y - the target node in R^1                       (input)
c    iy - the sequence number of the target node      (input)
c
c    f - the value of potential (output)
c        
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  nmax - the maximum number of points in a box on the finest level
c  nts - the number of chebychev nodes in the auxiliary discretization
c        of each box in the list 2 of the box being skeletonized
c  eps - the accuracy of the skeletonization
c  ifrelskel - set to 0 for skeletonization to absolute precision eps;
c              set to 1 for skeletonization to relative precision eps,
c              relative to Frobenius norms
c  iplot - the number of the GNUPLOT file to be produced by the
c        subroutine. Setting iplot .leq. 0 will suppress the plotting
c  lenw - the total length of the user-supplied array w (to be used 
c        for all purposes), in integer *4 elements
c
c
c                    Output parameters:
c
c  ier - error return code
c  nn - the total number of boxes created
c  w - the array containing all LOLO and TATA matrices, in addition 
c        to the structure created previously via preceding calls to 
c        the subroutine d1mstrcr, hudson_allskels (see)
c  keep - the number of elements in array w (in integer *4 words)
c        that should not be changed between the call to this subroutine
c        and subsequent use of the contents of array w
c  lused - the total length of the part of array that was used by
c        the subroutine
c
c        find the ends of the interval on which the structure is
c        to be built
c
        a=xs(1)-(xs(n)-xs(1))/100000
        b=xs(n)+(xs(n)-xs(1))/100000
c
c        build the boxes and the corresponding interaction lists
c
cccc        call prini(6,13)

        call d1mstrcr(ier,n,xs,a,b,nmax,iplot,nn,
     1      w,lenw,keep)


cccc        call prini(0,0)
c
c        determine the maximum sizes of various arrays to be used
c        by hudson_allskels
c
        iwww7=keep+10
        iskelout7=iwww7+n+1
c
        ilist7=iskelout7+n+1
        ilists7=ilist7+n+1  
        call hudson_stats234(w,xs,nn,
     1      nlists2_max,nlists3_max,nlists4_max,npts4_max,
     2      w(iwww7),w(iskelout7),w(ilist7),w(ilists7) )
c
c        allocate memory for hudson_allskels
c
        ninire=2
        nskelmax=nts*4
        if(nskelmax .gt. n) nskelmax=n
c
        nskeloutmax=nmax+(nlists2_max+nlists3_max)*nts+npts4_max
c
        iskel=keep+10
        lskel=n
c
        iipts_inbox=iskel+lskel
        lipts_inbox=nmax+2+nts
c
        iipivots=iipts_inbox+lipts_inbox
        lipivots=n
        if(nskeloutmax .gt. lipivots) lipivots=nskeloutmax
c
        ilist=iipivots+lipivots
        llist=0
c
        ilists=ilist+llist
        llists=0
c
        iamatr=ilists+llists
        lamatr=nskeloutmax*nskelmax*ninire+2*n
c
        iapts_inbox=iamatr+lamatr
        lapts_inbox=(nmax+2+nts)*ninire
c
        iskelout=iapts_inbox+lapts_inbox
        lskelout=nskeloutmax*ninire
c
        irnorms=iskelout+lskelout
        lrnorms=lipivots*ninire
c
        iinds=irnorms+lrnorms
        linds=nn+4
c
        iiiiskelout=iinds+linds
        liiiskelout=nskeloutmax+n
c
        iaddr=iiiiskelout+liiiskelout
        laddr=nn*8+10
c
        iarr=iaddr+laddr
c
c
        lused=0
c
c
c        skeletonize all boxes on all levels, and store the results
c
        w(21)=iaddr
        w(22)=iarr
c
        call hudson_allskels(ier,w,xs,nn,nts,eps,ifrelskel,
     a      w(iaddr),w(iarr),lused,
     1      w(iskel),w(iipts_inbox),w(iipivots),w(ilist),w(ilists),
     2      w(iamatr),lamatr,w(iapts_inbox),w(iskelout),w(irnorms),
     3      ncols_max)
c
        call hudson_allskels2(ier,w,xs,nn,nts,eps,ifrelskel,
     a      w(iaddr),w(iarr),lused,
     1      w(iskel),w(iipts_inbox),w(iipivots),w(ilist),w(ilists),
     2      w(iamatr),w(iapts_inbox),w(iskelout),w(irnorms),
     3      ncols_max,matfun,w(iinds),w(iiiiskelout) )

 2200 continue
c
        larr=lused+10
c
c        perform garbage collection
c
        iaddr2=keep+10
        call d1marrmove(w(iaddr),w(iaddr2),laddr)
        iaddr=iaddr2
c
        iarr2=iaddr+laddr+2
        call d1marrmove(w(iarr),w(iarr2),larr)
c
        iarr=iarr2
        w(21)=iaddr                
        w(22)=iarr                
c
        keep=iarr+larr
c
c        construct all lolos and tatas, and store them things
c
        iww=keep+10
c
        lw_used=keep+10
        call hudson_all_matrices(ier,nn,w,lenw,lw_used,xs,eps/10,
     1      ncols_max,nlists2_max,nlists3_max,nlists4_max,
     2      keep_out,lused,n,matfun,nmax)
c
        keep=keep_out
c
        w(24)=n
        w(25)=nn

c
        return
        end
c
c
c
c
c
        subroutine hudson_eval(w,charge,vals,w2,lenw2,lusedw2)
        implicit real *8 (a-h,o-z)
        save
        dimension charge(1),vals(1),w2(1)
        integer *4 w(1)
c
c        This subroutine evaluates the potential at the nodes
c        xs of a charge distribution at these same nodes.
c        It is fervently hoped that the necessary precomputation
c        has been performed via a preceding call to the subroutine
c        hudson_allchunks (see)
c
c        Actually, this is merely a memory-management routine;
c        all of the work is done by the subroutine hudson_eval0
c        (see)
c
c                     Input paramaters:
c  
c  w - the array containing all LOLO, TATA, LOTA, etc. matrices,
c         in addition to the structure created via calls to the
c         subroutines d1mstrcr, hudson_allskels (see). Must have been
c         created via a preceding call to the subroutine hudson_allchunks
c         (see)
c  charge - the charge distribution whose potential is to be
c         evaluated
c  lenw2 - the length of the work array w2, in real *8 words
c
c                     Output parameters:
c
c  vals - the potential at the nodes xs of the charges at these
c         same nodes specified via array charges 
c  lusedw2 - the part of the array w2 actually used by the subroutine
c         (in real *8 words)
c
c                     Work arrays:
c
c  w2 - must be suffciently long
c
c        . . . allocate memory for work arrays
c
        ier=0
        n=w(24)
        nboxes=w(25)

        lenfld=10 000
        ifld=1
        lfld=lenfld
c
        ifld2=ifld+lfld
        lfld2=lenfld
c
        ifld3=ifld2+lfld2
        lfld3=lenfld
c
        iwrea=ifld3+lfld3
        lwrea=lenfld
c
        iwint1=iwrea+lwrea
        lwint1=lenfld
c
        iwint2=iwint1+lwint1
        lwint2=lenfld
c
        iinds=iwint2+lwint2
        linds=nboxes+2
c
        ilist=iinds+linds
        llist=lenfld
c
        iw2=ilist+llist
c
        lenw2_2=lenw2-iw2-2
c
        call hudson_eval0(ier,n,w,nboxes,charge,vals,
     1      w2(ifld),w2(ifld2),w2(ifld3),w2(ilist),
     2      w2(iinds),w2(iwrea),w2(iwint1),w2(iwint2),
     3      w2(iw2),lenw2_2,lusedw2)
c
        return
        end
c
c
c
c
c
        subroutine hudson_eval0(ier,n,w,nboxes,charge,vals,
     1      fld,fld2,fld3,list,inds,wrea,wint1,wint2,w2,
     2      lenw2,lusedw2)
        implicit real *8 (a-h,o-z)
        save
        dimension fld(1),fld2(1),fld3(1),list(1),charge(1),
     1      vals(1),inds(1),wrea(1),wint1(1),wint2(1),w2(1)
c
        integer *4 w(1)
c
        ier=0
        iww=w(23)
c
c       initialize the storage-retrieval subroutine for 
c       the fields
c
        call hudson_init_flds(nboxes,w2,lenw2)
c
c       construct all outgoing expansions
c
        call hudson_all_lolos(ier,w,nboxes,charge,
     1      w2,icurr,w(iww),
     2      fld,fld2,fld3,inds,wrea,wint1,wint2)
c
c       apply all lotas
c
        call hudson_all_lotas(ier,w,nboxes,w2,icurr,
     1      fld,fld2,fld3,list,w(iww))
c
c        account for lists 3
c
        call hudson_all_lists3(ier,w,nboxes,w2,icurr,w(iww),
     1      fld,fld3,list)
c
c        account for lists 4
c
        call hudson_all_lists4(ier,w,nboxes,w2,icurr,w(iww),
     1      fld,fld2,fld3,list)
c
c        apply all tatas
c
        call hudson_all_tatas2(ier,w,nboxes,w2,icurr,
     1      vals,w(iww),
     2      fld,fld2,fld3,wrea,inds,wint1,wint2)
c
c
cccc        call prini(0,0)
        call hudson_all_tatas3(ier,w,nboxes,w2,icurr,
     1      vals,w(iww),
     2      fld,fld2,fld3,wrea,inds,wint1,wint2)

ccc        call prini(6,13)
c
c        account for all lists 1
c
        call hudson_selves_appl(ier,w,nboxes,charge,vals,w(iww),
     1      fld3,list)
c
        return
        end

c
c
c
c
c
        subroutine hudson_all_tatas2(ier,w,nboxes,w3,icurr,
     1      vals,ww,
     2      fld,fld2,fld4,fld5,inds,ito_infrom,ito_notinfrom)
        implicit real *8 (a-h,o-z)
        save
        real *8 fld(1),fld2(1),ab(2),vals(1),
     1      fld4(1),fld5(1)
        integer *4 box(10),inds(1),ww(1),w3(1),w(1),
     1      ito_infrom(2,1),ito_notinfrom(1)
c
c        . . . set the array of indices to -1
c
        n=w(24)
c
        do 1200 i=1,nboxes
c
        inds(i)=-1
c
cccc        if(i .le. 7) inds(i)=2
cccc        if(i .le. 5) inds(i)=2
 1200 continue

        inds(1)=2
c
c        apply tatas for all boxes
c
        do 3000 ijk=1,100
c
        ifacted=0
c
        do 2000 ibox=2,nboxes
c

        if(inds(ibox) .gt. 0) goto 2000
c
c        find this box's daddy; if the daddy has not been 
c        processed - skip
c
        call d1m_onebox(ibox,w,box,ab)
c
        idad=box(3)

        if( (inds(ibox) .le. 0) .and. (inds(idad) .gt. 0) ) then
            ifacted=1 
            inds(ibox)=1
        endif


c
        if(inds(idad) .lt. 0) goto 2000
        if(idad .lt. 2) goto 2000
c
c        retrieve the ibox's and dad's incoming expansions
c
        kind=2
        call hudson_retr_flds(jer,idad,kind,fld,
     1      lenarr,w3)
c
        if(jer .ne. 0) goto 2000
c
        kind=2
        call hudson_retr_flds(jer,ibox,kind,fld2,
     1      lenarr,w3)
c
        if(jer .ne. 0) then
            itype=2
            call hudson_skel_size(ier,w,ibox,itype,nskelbox)
            call hudson_setzero(fld2,nskelbox)
        endif
c
c        apply the tata matrix, and store the result in arrays
c        index2, w2
c
        iboxfrom=idad
        itypefrom=6
        iboxto=ibox
c
        call hudson_tata_funny(ier,iboxfrom,itypefrom,iboxto,
     1      w,nskelfrom,nskelto,fld,fld4,ww,
     2      ito_infrom,ito_notinfrom,fld5)
c
        do 1400 j=1,nskelto
c
        fld2(j)=fld2(j)+fld4(j)
 1400 continue
c
        kind=2
        call hudson_store_flds(ier,ibox,kind,fld2,
     1      nskelbox,w3,icurr)
c 
        inds(ibox)=1
        ifacted=1
c
 2000 continue
c
        if(ifacted .eq. 0) goto 3200
c
 3000 continue
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
        subroutine hudson_all_tatas3(ier,w,nboxes,w3,icurr,
     1      vals,ww,
     2      fld,fld2,fld4,fld5,inds,ito_infrom,ito_notinfrom)
        implicit real *8 (a-h,o-z)
        save
        real *8 fld(1),fld2(1),ab(2),vals(1),
     1      fld4(1),fld5(1)
        integer *4 box(10),inds(1),ww(1),w3(1),w(1),
     1      ito_infrom(2,1),ito_notinfrom(1)
c
        n=w(24)
c
c        apply tatas for all boxes
c
        do 2000 ibox=2,nboxes
c
cccc        if( (inds(ibox) .eq. 2) .and. (ijk .gt. 1) ) goto 2000
        if(inds(ibox) .eq. 2) goto 2000
c
c        if this is a childless box - convert the "regular"
c        incoming potential into the "big" one
c
        call d1m_onebox(ibox,w,box,ab)
c
        if(box(4) .gt. 0) goto 2000
        if(box(5) .gt. 0) goto 2000
c
        iboxfrom=ibox
        itypefrom=18
        iboxto=ibox
c
        kind=2
        call hudson_retr_flds(jer,ibox,kind,fld2,
     1      lenarr,w3)
c
        if(jer .ne. 0) goto 2000
c
        call hudson_tata_funny(ier,iboxfrom,itypefrom,iboxto,
     1      w,nskelfrom,nskelto,fld2,fld4,ww,
     2      ito_infrom,ito_notinfrom,fld5)
c
        kind=4
        call hudson_retr_flds(jer,ibox,kind,fld,
     1      lenarr,w3)
c
        if(jer .ne. 0) call hudson_setzero(fld,nskelto)

        do 1900 j=1,nskelto
c
        fld4(j)=fld4(j)+fld(j)
 1900 continue
c
        kind=4
        call hudson_store_flds(ier,ibox,kind,fld4,
     1      nskelto,w3,icurr)
c 
 2000 continue
c
c        account for all list 1's
c
        call hudson_all_lists1(ier,w,nboxes,
     1      w3,icurr,ww,
     2      fld,fld2,fld4,ito_infrom)
c
        do 3400 i=1,n
c
        vals(i)=0
 3400 continue
c
c        evaluate the potentials at all points
c
        do 4000 ibox=1,nboxes
c
c       if this box has kids - skip
c
        call d1m_onebox(ibox,w,box,ab)
c
        i1=box(4)
        i2=box(5)
c
        if(i1 .gt. 0) goto 4000
        if(i2 .gt. 0) goto 4000
c
c        retrieve from appropriate arrays the "big" eval matrix for 
c        this box, and the incoming potential on its skeleton
c
        kind=4
        call hudson_retr_flds(jer,ibox,kind,fld,
     1      lenarr,w3)
c
c       evaluate the incoming potentials at the appropriate nodes
c
        ip=box(6)
        nptsbox=box(7)
c
        iboxfrom=ibox
        itypefrom=16
        iboxto=ibox
        nskelto=nptsbox
c
        call hudson_tata_funny(ier,iboxfrom,itypefrom,iboxto,
     1      w,nskelfrom,nskelto,fld,vals(ip),ww,
     2      ito_infrom,ito_notinfrom,fld5)
c
 3800 continue
c
 4000 continue
c
        return
        end
c
c
c
c
c
        subroutine hudson_selves_appl(ier,w,nboxes,charge,vals,ww,
     1      fld3,list)
        implicit real *8 (a-h,o-z)
        save
        real *8 w(1),ab(2),fld3(1),vals(1),charge(1)
        integer *4 list(1),box(10)
c
c        This subrouutine evaluate self-interactions for 
c        all childless boxes
c
c                  Input parameters:
c
c  w - the array containing the structure
c  nboxes - the number of boxes created (one fervently hopes) via 
c        a preceding call to the subroutine d1mstrcr (see)
c  charge - the user-supplied charge distribution at the nodes xs
c  vals - the output array of values to which the self-interactions
c        will be added by this subroutine
c  ww - the array containing operators LOLO, TATA, LOTA, as
c       (hopefully) created by a preceding call to the subroutine
c       hudson_allchunks (see)
c
c                  Output parameters:
c
c  ier - error return code
c  vals - the output array of values to which the self-interactions
c       of all boxes have been added
c
c
c                  Work arrays:
c
c  fld3,list
c
        do 2000 ibox=1,nboxes
c
        call d1m_onebox(ibox,w,box,ab)
c
        if(box(4) .gt. 0) goto 2000
        if(box(5) .gt. 0) goto 2000
c
c        retrieve the ibox's nodes and charges
c
        i1_ibox=box(6)
        npts_ibox=box(7)
c
c        account for the interaction of ibox with itself
c
        iboxfrom=ibox
        itypefrom=23
        iboxto=ibox
c
        call hudson_lota_funny(iboxfrom,itypefrom,iboxto,
     1      w,nskelfrom,nskelto,charge(i1_ibox),fld3,ww)
c
        do 1400 j=1,npts_ibox
        vals(i1_ibox+j-1)=vals(i1_ibox+j-1)+fld3(j)
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
        subroutine hudson_all_lists4(ier,w,nboxes,w3,icurr,ww,
     1      fld,fld2,fld3,list)
        implicit real *8 (a-h,o-z)
        save
        real *8 w(1),fld(1),fld3(1),fld2(1),ww(1)
        integer *4 list(1),w3(1)
c
c        This subroutine applies the list 4 operators
c        to all chunks on all levels
c
c        
c                  Input parameters:
c
c  w - the array containing the structure; must have been created via 
c        a preceding call to the subroutine d1mstrcr (see)
c  nboxes - the total number of boxes in the structure
c  w3 - storage area containing all fields; please note that this
c       is also an output parameter
c  ww - storage area containing all sorts of matrices
c
c                  Output parameters:
c
c  ier - error return code
c  w3 - storage area containing all fields; please note that this
c       is also an input parameter
c  icurr - the number of elements in array w3 (in integer *4 elements)
c       used in array w3 at this time
c
c                  Work arrays:
c
c  fld,fld2,fld3,list
c
c
c        . . . evaluate list potentials for all boxes
c
        do 2000 ibox=1,nboxes
c
c        retrieve list 4 for this box
c
        itype=4
        call d1m_onelist(ibox,itype,w,list,nlist)
c
        if(nlist .le. 0) goto 2000
c
c        retrieve the ibox's outgoing expansion
c
        kind=1
        call hudson_retr_flds(jer,ibox,kind,fld,
     1      lenarr,w3)
c
        if(jer .gt. 0) goto 2000
c
c        process the elements of ibox's list 4 one after another
c
        do 1800 i=1,nlist
c
        jbox=list(i)
c     
c        evaluate the potential of the regular outgoing skeleton 
c        of ibox at the nodes in the big incoming skeleton of jbox
c
        iboxfrom=ibox
        itypefrom=25
        iboxto=jbox
c
        call hudson_lota_funny(iboxfrom,itypefrom,iboxto,
     1      w,nskelfrom,nskelto,fld,fld3,ww)
c
        nskeljbox=nskelto
c
c        store the obtained incoming potential
c
        kind=4
        call hudson_retr_flds(jer,jbox,kind,fld2,
     1      lenarr,w3)
c
        if(jer .ne. 0) call hudson_setzero(fld2,nskeljbox)
c
        do 1400 j=1,nskeljbox
        fld2(j)=fld2(j)+fld3(j)
 1400 continue
c
c        store the result
c
        kind=4
        call hudson_store_flds(ier,jbox,kind,fld2,
     1      nskeljbox,w3,icurr)
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
        subroutine hudson_all_lists3(ier,w,nboxes,w3,icurr,ww,
     1      fld1,fld3,list)
        implicit real *8 (a-h,o-z)
        save
        real *8 w(1),ab(2),fld1(1),fld3(1)
        integer *4 list(1),box(10),w3(1)
c
c        This subroutine applies the list 3 operators
c        to all chunks on all levels
c
c        
c                  Input parameters:
c
c  w - the array containing the structure; must have been created via 
c        a preceding call to the subroutine d1mstrcr (see)
c  nboxes - the total number of boxes in the structure
c  w3 - storage area containing all fields; please note that this
c       is also an output parameter
c  ww - storage area containing all sorts of matrices
c
c                  Output parameters:
c
c  ier - error return code
c  w3 - storage area containing all fields; please note that this
c       is also an input parameter
c  icurr - the number of elements in array w3 (in integer *4 elements)
c       used in array w3 at this time
c
c                  Work arrays:
c
c  fld1,fld3,list
c
c        . . . apply list3 for all boxes
c
        do 2000 ibox=1,nboxes
c
c        retrieve list 3 for this box
c
        itype=3
        call d1m_onelist(ibox,itype,w,list,nlist)
c
        if(nlist .le. 0) goto 2000
c
c        retrieve the ibox's points and charges
c
        call d1m_onebox(ibox,w,box,ab)
c
        i1=box(6)
        npts=box(7)
c
c        process the elements of ibox's list 3 one after another
c
        do 1800 i=1,nlist

        jbox=list(i)
c
c        evaluate the potential of charges in the "big" skeleton 
c        in ibox at the incoming skeleton of jbox
c
c        . . . retrieve the "big" outgoing potential on ibox
c
        kind=3
        call hudson_retr_flds(jer,ibox,kind,fld1,
     1      lenarr,w3)
c
c        evaluate the potential of the "big" outgoing skeleton 
c        of ibox at the nodes in the "regular" incoming skeleton 
c        of jbox
c
        iboxfrom=ibox
        itypefrom=27
        iboxto=jbox
c
        call hudson_lota_funny(iboxfrom,itypefrom,iboxto,
     1      w,nskelfrom,nskelto,fld1,fld3,ww)
c
        nskeljbox=nskelto
c
c       retrieve the existing incoming potential at the 
c       skeleton of jbox
c
        kind=2
        call hudson_retr_flds(jer,jbox,kind,fld1,
     1      lenarr,w3)
c
        if(jer .ne. 0) call hudson_setzero(fld1,nskeljbox)
c
        do 1400 j=1,nskeljbox
        fld1(j)=fld1(j)+fld3(j)
 1400 continue
c
c        store the result
c
        kind=2
        call hudson_store_flds(ier,jbox,kind,fld1,
     1      nskeljbox,w3,icurr)
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
        subroutine hudson_all_lists1(ier,w,nboxes,
     1      w3,icurr,ww,
     2      fld,fld2,fld4,list)
        implicit real *8 (a-h,o-z)
        save
        real *8 w(1),ab(2),fld2(1),fld(1),fld4(1)
        integer *4 list(1),box(10),w3(1)
c
c        This subroutine applies to the appropriate charge 
c        distributions the matrices accounting for the 
c        interactions of all childless boxes with their
c        list 1's. Please note that this subroutine DOES NOT
c        deal with self-interactions of childless boxes; for
c        that, see the subroutine hudson_selves_appl.
c
c                  Input parameters:
c
c  w - the array containing the structure; must have been created via 
c        a preceding call to the subroutine d1mstrcr (see)
c  nboxes - the number of boxes created (one fervently hopes) via 
c        a preceding call to the subroutine d1mstrcr (see)
c  w3 - storage area containing all fields; please note that this
c       is also an output parameter
c  ww - storage area containing all sorts of matrices
c  
c                  Output parameters:
c
c  ier - error return code
c  w3 - storage area containing all fields; please note that this
c       is also an output parameter
c  icurr - the number of elements in array w3 (in integer *4 elements)
c       used in array w3 at this time
c
c                  Work arrays:
c
c  fld,fld2,fld4,list
c
c
c        . . . apply matrices in lists 1 for all boxes
c
        do 3000 ijk=1,2
c
        do 2000 ibox=1,nboxes
c
        call d1m_onebox(ibox,w,box,ab)
c
        if(box(4) .gt. 0) goto 2000
        if(box(5) .gt. 0) goto 2000
c
c        retrieve list 1 for this box
c
        itype=1
        if(ijk .eq. 2) itype=11
        call d1m_onelist(ibox,itype,w,list,nlist)
c
        if(nlist .le. 0) goto 2000
c
c        retrieve the "big" outgoing field on ibox
c
        kind=3
        call hudson_retr_flds(jer,ibox,kind,fld,
     1      lenarr,w3)
c
c        process the elements of ibox's list 1 one after another
c
        do 1800 i=1,nlist
c
        jbox=list(i)
c
c        retrieve the incoming "big" skeleton inside jbox
c
        iboxto=jbox
        itypeto=4
c
c        apply the "lota" matrix to the "big" outgoing potential
c        on ibox
c
        itypefrom=21
        iboxfrom=ibox
        call hudson_lota_funny(iboxfrom,itypefrom,iboxto,
     1      w,nskelfrom,nskelto,fld,fld2,ww)
c
c        add the obtained incoming potential on jbox to
c        the potential already there
c
        kind=4
        call hudson_retr_flds(jer,jbox,kind,fld4,
     1      lenarr,w3)
c
        if(jer .ne. 0) call hudson_setzero(fld4,nskelto)
c
        do 1400 j=1,nskelto
c
        fld4(j)=fld2(j)+fld4(j)
 1400 continue
c
c       . . . store the result
c
        kind=4
        call hudson_store_flds(ier,jbox,kind,fld4,
     1      nskelto,w3,icurr)
c
 1800 continue
c
 2000 continue
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
        subroutine hudson_lota_funny(iboxfrom,itypefrom,iboxto,
     1      w,nskelfrom,nskelto,fld,fld2,ww)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),fld(1),fld2(1)
        integer *4 ww(1)
c
c       retrieve alota from array ww
c
        key=0
        call hudson_linret(ier,iboxfrom,itypefrom,iboxto,key,ww,
     1      ip_intarr,lenint,ip_reaarr,lenrea)
c
c       unpack integer data from array ww
c
        i=ip_intarr
        nskelfrom=ww(i)
        nskelto=ww(i+1)
c
c       apply the lota operator to the "from" potential
c
        call hudson_matvec(ww(ip_reaarr),nskelto,nskelfrom,fld,fld2)
c
        return
        end
c
c
c
c
c
        subroutine hudson_tata_funny(ier,iboxfrom,itypefrom,iboxto,
     1      w,nskelfrom,nskelto,fld,fld4,ww,
     2      ito_infrom,ito_notinfrom,fld5)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),ito_infrom(2,1),
     1      fld(1),fld4(1),fld5(1),ito_notinfrom(1)
        integer *4 ww(1)
c
c       retrieve the tata operator from array ww
c

        lenint=0
        lenrea=0
        key=0
        call hudson_linret(ier,iboxfrom,itypefrom,iboxto,key,ww,
     1      ip_intarr,lenint,ip_reaarr,lenrea)
c 
        if( (lenint .eq. 0) .and. (lenrea .eq. 0) ) then
           ier=4
            return
        endif
c
c
c       unpack integer data from array ww
c
        i=ip_intarr
        nskelfrom=ww(i)
        nskelto=ww(i+1)
c
        nto_infrom=ww(i+2)
        nto_notinfrom=ww(i+3)
c
        iito_infrom=ww(i+4)
        iito_notinfrom=ww(i+5)
c
        call hudson_iarrmove(ww(i+iito_infrom-1),
     1      ito_infrom,nto_infrom*2)
c
        call hudson_iarrmove(ww(i+iito_notinfrom-1),
     1      ito_notinfrom,nto_notinfrom)
c
c       apply the TATA operator to the "from" potential
c
        call hudson_matvec_adj(ww(ip_reaarr),nskelfrom,nto_notinfrom,
     1      fld,fld5)
c
c        finished applying the compressed matrix; shuffle
c
        do 2400 i=1,nto_notinfrom
c
        j=ito_notinfrom(i)
        fld4(j)=fld5(i)
 2400 continue
c
c
        do 2600 i=1,nto_infrom
c
        j1=ito_infrom(1,i)
        j2=ito_infrom(2,i)
c
        fld4(j1)=fld(j2)
 2600 continue
c
        return
        end
c
c
c
c
c
        subroutine hudson_all_lotas(ier,w,nboxes,w3,icurr,
     1      fld,fld2,fld3,list,ww)
        implicit real *8 (a-h,o-z)
        save
        real *8 fld(1),fld2(1),fld3(1)
        integer *4 list(1),w(1),w3(1),ww(1)
c
c        apply lotas for all boxes
c
        do 2000 ibox=1,nboxes
c
c        retrieve list 2 for this box
c
        itype=2
        call d1m_onelist(ibox,itype,w,list,nlist)
c
        if(nlist .le. 0) goto 2000
c
c        retrieve the ibox's outgoing expansion
c
        kind=1
        call hudson_retr_flds(jer,ibox,kind,fld,
     1      lenarr,w3)
c
        if(jer .ne. 0) goto 2000
c
        itype=1
        call hudson_skel_size(ier,w,ibox,itype,nskelibox)
c
        do 1600 i=1,nlist
c
        jbox=list(i)
        itype=2
        call hudson_skel_size(jer,w,jbox,itype,nskeljbox)
c
        kind=2
        call hudson_retr_flds(jer,jbox,kind,fld2,
     1      lenarr,w3)
c
        if(jer .ne. 0) call hudson_setzero(fld2,nskeljbox)
c
c
c        evaluate the potential of the charges 
c        at the skeleton nodes in jbox
c
        iboxfrom=ibox
        itypefrom=33
        iboxto=jbox
c
        call hudson_lota_funny(iboxfrom,itypefrom,iboxto,
     1      w,nskelibox,nskeljbox,fld,fld3,ww)
c
        do 1400 j=1,nskeljbox
c
        fld2(j)=fld2(j)+fld3(j)
 1400 continue
c
        kind=2
        call hudson_store_flds(ier,jbox,kind,fld2,
     1      nskeljbox,w3,icurr)
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
        subroutine hudson_all_lolos(ier,w,nboxes,charges,
     1      w3,icurr,ww,
     2      fld,fld2,fld4,inds,wrea,wint1,wint2)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),ab(2),charges(1),fld(1),fld2(1),
     1      fld4(1),wrea(1),wint1(1),wint2(1)
        integer *4 box(10),inds(1),ww(1),w3(1)
c
c        initialize the array inds storing the information about
c        the condition of boxes (whether they have been processed 
c        or not)
c
        do 1200 i=1,nboxes
c
        inds(i)=0
 1200 continue

c
c        construct outgoing expansions for all childless boxes
c
        do 2000 ibox=1,nboxes
c
c        if this box has kids - skip it
c
        call d1m_onebox(ibox,w,box,ab)
c
        if(box(4) .gt. 0) goto 2000
        if(box(5) .gt. 0) goto 2000
c
c        this box is childless. construct its "big"
c        outgoing expansion
c
        i1box=box(6)
        itypefrom=17
c
        call hudson_lolo_funny(ier,ibox,itypefrom,ibox,
     1      w,nskelbox,ndummy,charges(i1box),fld4,ww,
     2      wrea,wint1,wint2)
c
        kind=3
        call hudson_store_flds(ier,ibox,kind,fld4,
     1      nskelbox,w3,icurr)
c
c        convert the "big" outgoing expansion into the
c        "standard" one
c
        iboxfrom=ibox
        itypefrom=19
        iboxto=ibox
c
        call hudson_lolo_funny(ier,iboxfrom,itypefrom,ibox,
     1      w,nskelbox2,nskelbox,fld4,fld,ww,
     2      wrea,wint1,wint2)
c
        kind=1
        call hudson_store_flds(ier,ibox,kind,fld,
     1      nskelbox,w3,icurr)
c
        inds(ibox)=1
 2000 continue
c
c       construct all outgoing expansions on all levels
c
        do 4000 ijk=1,100
c
        ifacted=0
c
        do 3600 ibox=2,nboxes
c
c        if this box has been processed - skip
c
        if(inds(ibox) .eq. 1) goto 3600
c
c       if at least one of this box's kids has not been processed
c       - skip
c
        call d1m_onebox(ibox,w,box,ab)
c
        i1=box(4)
        i2=box(5)
c
        if(i1 .gt. 0) then
            if(inds(i1) .ne. 1) goto 3600
        endif
c
        if(i2 .gt. 0) then
            if(inds(i2) .ne. 1) goto 3600
        endif
c
c       this box needs to be processed - act accordingly
c
c        . . . if this box has both kids - merge them things
c
        if( (i1 .gt. 0) .and. (i2. gt. 0) ) goto 3000
c
c        . . . if this box has only the first son - copy the 
c              sonnie's outgoing expansion into daddy's
c
        if(i1 .le. 0) goto 2400
c
        kind=1
        call hudson_retr_flds(ier,i1,kind,fld,
     1      lenarr,w3)
c
        call hudson_store_flds(ier,ibox,kind,fld,
     1      lenarr,w3,icurr)
c
        inds(ibox)=1
        ifacted=1
c
        goto 3600
c
 2400 continue
c
c        . . . if this box has only the second son - copy the 
c              sonnie's outgoing expansion into daddy's
c
        if(i2 .le. 0) then
            call prin2('nonsense in hudson_all_lolos; bombing*',x,0)
            stop
        endif
c
        kind=1
        call hudson_retr_flds(ier,i2,kind,fld,
     1      lenarr,w3)
c
        call hudson_store_flds(ier,ibox,kind,fld,
     1      lenarr,w3,icurr)
c
        inds(ibox)=1
        ifacted=1
c
        goto 3600
c
 3000 continue
c
c        this box has both sons. Act accordingly
c
        ison1=box(4)
        ison2=box(5)
c
c       . . . .translate the first sonny
c
        kind=1
        call hudson_retr_flds(ier,ison1,kind,fld,
     1      lenarr,w3)
c
        itypefrom=5
        call hudson_lolo_funny(ier,ison1,itypefrom,ibox,
     1      w,nskelson1,nskelbox,fld,fld2,ww,
     2      wrea,wint1,wint2)
c
c       . . . .translate the second sonny
c
        call hudson_retr_flds(ier,ison2,kind,fld,
     1      lenarr,w3)
c
        itypefrom=5
        call hudson_lolo_funny(ier,ison2,itypefrom,ibox,
     1      w,nskelson2,nskelbox,fld,fld4,ww,
     2      wrea,wint1,wint2)

c
        do 3200 j=1,nskelbox
c
        fld2(j)=fld2(j)+fld4(j)
 3200 continue
c
c       . . . store the result
c
        kind=1
        call hudson_store_flds(ier,ibox,kind,fld2,
     1      nskelbox,w3,icurr)
c
        inds(ibox)=1
        ifacted=1
 3600 continue
c
        if(ifacted .eq. 0) goto 4200
c
 4000 continue
 4200 continue
c
        return
        end
c
c
c
c
c
        subroutine hudson_lolo_funny(ier,iboxfrom,itypefrom,iboxto,
     1      w,nskelfrom,nskelto,fld_in,fld_out,ww,
     2      fld2,from_into,from_notinto)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),fld_in(1),fld_out(1),fld2(1)
        integer *4 ww(1),from_into(2,1),from_notinto(1)
c
c       retrieve alolo_out from array ww
c
        ier=0
c
        key=0
        call hudson_linret(jer,iboxfrom,itypefrom,iboxto,key,ww,
     1      ip_intarr,lenint,ip_reaarr,lenrea)
c
        if( (lenint .eq. 0) .and. (lenrea .eq. 0) ) then
            ier=4
            return
        endif
c
        if(itypefrom .ne. 19) goto 1100
c
 1100 continue
c
c       unpack integer data from array ww
c
        i=ip_intarr
        nskelfrom=ww(i)
        nskelto=ww(i+1)
c
        nfrom_into=ww(i+2)
        nfrom_notinto=ww(i+3)
c
        ifrom_into=ww(i+4)
        ifrom_notinto=ww(i+5)
c
        call hudson_iarrmove(ww(i+ifrom_into-1),
     1      from_into,nfrom_into*2)

        call hudson_iarrmove(ww(i+ifrom_notinto-1),
     1      from_notinto,nfrom_notinto)
c
c        apply the LOLO operator to the sonny's potential, and 
c        add the result to daddy's potential
c
        do 2000 i=1,nfrom_notinto
        j=from_notinto(i)
        fld2(i)=fld_in(j)
 2000 continue
c
        call hudson_matvec(ww(ip_reaarr),nskelto,nfrom_notinto,
     1      fld2,fld_out)
c
c       add in the diagonal part of the matrix
c
        do 2200 i=1,nfrom_into
c
        j1=from_into(1,i)
        j2=from_into(2,i)
c
        fld_out(j2)=fld_out(j2)+fld_in(j1)
 2200 continue
c
        return
        end
c
c
c        Glossary of elements in array addr
c
c        addr(1,i) = i - sequence number of the column of addr
c        addr(2,i) = ilast   - link to the preceding column of addr 
c                       related to  the box number ibox
c        addr(3,i) = ibox - the box for which the information is being 
c                       stored/retrieved (the key on which the linking
c                       is done)
c        addr(4,i) = jbox - the first auxiliary key
c        addr(5,i) = itype - the second auxiliary key
c        addr(6,i) = key - the third auxiliary key
c        addr(7,i) = i_intarr - the address in array arr of the integer 
c                       data being stored/retrieved
c        addr(8,i) = lenint - the length (in integer *4 elements) of the 
c                       integer array being stored/retrieved
c        addr(9,i) = i_reaarr - the address in array arr of the real
c                       data being stored/retrieved
c        addr(10,i) = lenrea - the length  (in integer *4 elements) of 
c                       the real array being stored/retrieved
c
c
c
c
c


c
c
c
c
c
        subroutine hudson_store_flds(ier,ibox,kind,arr,
     1      lenarr,w,icurr_out)
        implicit real *8 (a-h,o-z)
        save
        integer *4 w(1),arr(1)
c
        ier=0
        ninire=2
c
        iindex2=w(1)
        iw=w(2)
        icurr=w(3)
        lenw=w(4)
c
        if(lenw .lt. (icurr+lenarr)*ninire+2) then
            ier=16
            call prinf('bombing from hudson_store_flds with ier=*',
     1                 ier,1)
        endif
c
        call hudson_store_flds0(ier,ibox,kind,arr,
     1      lenarr,w(iindex2),w(iw),icurr)
c
        w(3)=icurr
c
        icurr_out=icurr+iw+2
c
        return
        end
c
c
c
c
c
        subroutine hudson_init_flds(nboxes,w2,lenw2)
        implicit real *8 (a-h,o-z)
        save
        integer *4 w2(1)
c
c       initialize the storage-retrieval subroutine for 
c       the fields
c
        iin2=11
        lin2=nboxes*12+10
c
        iw2=iin2+lin2
c
        w2(1)=iin2
        w2(2)=iw2
        w2(3)=1
        w2(4)=lenw2
c
        call hudson_init_flds0(w2(iin2),nboxes)
c
        return
        end
c
c
c
c
c
        subroutine hudson_retr_flds(ier,ibox,kind,arr,
     1      lenarr,w)
        implicit real *8 (a-h,o-z)
        save
        integer *4 w(1),arr(1)
c
        iindex2=w(1)
        iw=w(2)
c
        call hudson_retr_flds0(ier,ibox,kind,arr,lenarr,
     1      w(iindex2),w(iw))
c
        return
        end


c
c
c
c
c
        subroutine hudson_store_flds0(ier,ibox,kind,arr,
     1      lenarr,index2,w,icurr)
        implicit real *8 (a-h,o-z)
        save
        real *8 w(1),arr(1)
        integer *4 index2(2,6,1)
c
c        One call to this subroutine stores in the array w one
c        of field representations, to be later retrieved by the
c        entry hudson_retr_flds0 of this subroutine (see). The type of 
c        the field to be stored is specified by the parameter
c        kind, as follows:
c
c     kind = 1 - regular outgoing expansion 
c     kind = 2 - regular incoming expansion
c     kind = 3 - "big" outgoing expansion; exists for childless 
c        chunks only
c     kind = 4 "big" incoming expansion; exists for childless 
c        chunks only
c    
c  w - main storage area for the potentials of all persuasions 
c
c        . . . store the data about this array in arrays 
c              index2, w
c
        if(index2(1,kind,ibox) .gt. 0) goto 1400
c
        index2(1,kind,ibox)=icurr
        index2(2,kind,ibox)=lenarr
c
        call hudson_rarrcopy(arr,w(icurr),lenarr)
        icurr=icurr+lenarr+5       
c
        return
c
 1400 continue
c
        ii=index2(1,kind,ibox) 
        ll=index2(2,kind,ibox) 
        call hudson_rarrcopy(arr,w(ii),ll)
c
        return
c
c
c
c
        entry hudson_retr_flds0(ier,ibox,kind,arr,lenarr,
     1      index2,w)
c
        ier=0
c
        ii=index2(1,kind,ibox)
        lenarr=index2(2,kind,ibox)
c
        if(ii .le. 0) then
            ier=4
            return
        endif
c
        call hudson_rarrcopy(w(ii),arr,lenarr)
c
        return    
c
c
c
c        
        entry hudson_init_flds0(index2,nboxes)
c
        do 4400 i=1,nboxes
        do 4200 j=1,6
c
        index2(1,j,i)=-1
        index2(2,j,i)=-1
 4200 continue
 4400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine hudson_rarrcopy(x,y,n)
        implicit real *8 (a-h,o-z)
        save
        real *8 x(1),y(1)
c 
        do 2400 i=1,n
        y(i)=x(i)
 2400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine hudson_setzero(x,n)
        implicit real *8 (a-h,o-z)
        save
        real *8 x(1)
c 
        do 2400 i=1,n
        x(i)=0
 2400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine hudson_matvec_adj(a,n,m,x,y)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(n,m),x(n),y(m),cd
c 
c        apply the matrix a to the vector x obtaining y
c 
        do 1400 i=1,m
        cd=0
        do 1200 j=1,n
        cd=cd+a(j,i)*x(j)
 1200 continue
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
        subroutine hudson_matvec(a,n,m,x,y)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(n,m),x(m),y(n),cd
c 
c        apply the matrix a to the vector x obtaining y
c 
        do 1400 i=1,n
        cd=0
        do 1200 j=1,m
        cd=cd+a(i,j)*x(j)
 1200 continue
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
        subroutine hudson_all_matrices(ier,nboxes,w,lenw,
     1      lw_used,xs,eps,ncols_max7,
     2      nlists2_max,nlists3_max,nlists4_max,
     3      keep_out,lused,n,matfun,nmax)
        implicit real *8 (a-h,o-z)
        save
        dimension ab(2),xs(1)
        integer *4 box(10),w(1)
c
c        For the user-specified nodes xs and box structure (hopefully, 
c        constructed via a preceding call to subroutines d1mstrcr, 
c        hudson_allskels, hudson_allskels2) this subroutine builds 
c        all interaction matrices (LOLO, TATA, EVAL, EXPAND, etc. ) 
c        for all boxes on all levels. The matrices are stored in the 
c        same array w in a linked list, from which sad condition they 
c        can be extracted via the subroutine linret (see).
c
c                  Input parameters:
c
c  nboxes - the total number of boxes in the structure
c  w - the array containing the structure; must have been created via 
c        a preceding call to the subroutine d1mstrcr (see)
c  lenw - the total length of the user-supplied array w (to be used 
c        for all purposes), in integer *4 elements
c  lw_used - the part of the array w that is already occupied by the 
c        data created by the preceding calls to subroutines d1mstrcr,
c        hudson_allskels (see)
c  xs - the used-specified nodes in the simulation
c  eps - the accuracy of the skeletonization
c  ncols_max - the maximum number of elements in any skeleton for any
c        box on any level
c  nlists2_max - the maximum total number of boxes in the lists 2 
c        of the box ibox and all of his ancestors, among all boxes
c        ibox n the structure
c  nlists4_max - the maximum total number of boxes in the lists 4 
c        of the box ibox and all of his ancestors, among all boxes
c        ibox n the structure
c
c                    Output parameters:
c
c  ier - error return code
c  w - the array containing all interaction matrices, in addition 
c        to the structure created previously via preceding calls to 
c        the subroutine d1mstrcr, hudson_allskels hudson_allskels2
c  keep_out - the number of elements in array w (in integer *4 words)
c        that should not be changed between the call to this subroutine
c        and subsequent use of the contents of array w
c  lused - the total length of the part of array that was used by
c        the subroutine (in integer *4 words)
c
c
c
c        . . . allocate memory for the construction of LOLOs,
c              TATAs, LOTAs, etc.
c
        ncols_max=ncols_max7
        if(nmax+2 .gt. ncols_max) ncols_max=nmax+2
c
        nlists234_max=nlists2_max+nlists3_max+nlists4_max
c
        ninire=2
        lskelbox=(ncols_max+2)*ninire
        lskelout=(nlists234_max*ncols_max+100)*ninire
c
cc        lskelbox=lskelbox*2
cc        lskelout=lskelout*2
c
        ialolo=lw_used+1
        lalolo=lskelbox*lskelout/ninire
c
        ievalfrom=ialolo+lalolo
        levalfrom=(lskelbox*lskelout+100)*ninire 
c
        ievalto=ievalfrom+levalfrom
        levalto=(lskelbox*lskelout+100)*ninire 
c
        ialolo_out=ievalto+levalto
        lalolo_out=lalolo
c
        iwork=ialolo_out+lalolo_out
        lwork=lalolo*16
c
        iskelfrom=iwork+lwork
        lskelfrom=(ncols_max+2)*ninire +1000
c
        iskelto=iskelfrom+lskelfrom
        lskelto=(ncols_max+2)*ninire +1000
c
        iskelout=iskelto+lskelto
        lskelout=(nlists234_max*ncols_max+100)*ninire
        lskelout=lskelout+4000
c
        irnorms=iskelout+lskelout
        lrnorms=nlists234_max*ncols_max+100
        lrnorms=lrnorms*ninire  +4000
cccc        lrnorms=lrnorms*ninire  
c
        iipnts_from=irnorms+lrnorms
        lipnts_from=ncols_max+2  +1000
cccc        lipnts_from=ncols_max+2 
c
        iipnts_to=iipnts_from+lipnts_from
        lipnts_to=ncols_max+2 +1000
cccc        lipnts_to=ncols_max+2 
c
        iipntsout=iipnts_to+lipnts_to        
        lipntsout=(nlists234_max*ncols_max+100)+1000
c
        iipnts_jbox=iipntsout+lipntsout
        lipnts_jbox=ncols_max+2 +1000
cccc        lipnts_jbox=ncols_max+2 
c
        iiarr_arr=iipnts_jbox+lipnts_jbox
        liarr_arr=ncols_max*3+100+1000
cccc        liarr_arr=ncols_max*3+100
c
        iifrom_into=iiarr_arr+liarr_arr
        lifrom_into=2*ncols_max+60 +1000
cccc        lifrom_into=2*ncols_max+60 
c
        iifrom_notinto=iifrom_into+lifrom_into
        lifrom_notinto=ncols_max+30 +1000
cccc        lifrom_notinto=ncols_max+30 
c
        ltot=iifrom_notinto+lifrom_notinto
c
c        one after another, construct the LOLO matrices
c        for all boxes in the structure
c
        iww=ltot+1+200 000
        iww=ltot+1

        lenww=lenw-iww-2
        call hudson_linini(nboxes,w(iww),lenww)
c
        w(23)=iww
c
        do 2000 ibox=2,nboxes
c
c        construct the LOLO matrix for this box
c
        call d1m_onebox(ibox,w,box,ab)

        idad=box(3)
c
        iboxfrom=ibox
        itypefrom=1
c
        iboxto=idad
        itypeto=1
c
        itypeout=1
c
        call hudson_lolo_bld(jer,matfun,iboxfrom,itypefrom,
     1      iboxto,itypeto,itypeout,w,xs,w(ialolo),
     2      nskelbox,nskeldad,eps,
     3      w(iskelfrom),w(iskelto),w(iskelout),w(ievalfrom),w(ievalto),
     4      w(irnorms),w(iwork),w(iipnts_from),w(iipnts_to),
     5      w(iipntsout))
c
c        . . . find this guy's daddy
c
        call d1m_onebox(ibox,w,box,ab)
        idad=box(3)
c
c        compress the obtained lolo matrix and store it
c
        itypeout=5
        itypefrom=1
        itypeto=1
        iboxfrom=ibox
        iboxto=idad
        call hudson_lolo_funny0(jer,itypeout,iboxfrom,itypefrom,
     1      iboxto,itypeto,w,xs,w(ialolo),nskelbox,nskeldad,
     2      w(ialolo_out),w(iww),keep,
     3      w(iskelfrom),w(iskelto),w(iipnts_from),w(iipnts_to),
     4      w(iifrom_into),w(iifrom_notinto),w(iiarr_arr))
c
 2000 continue
c
c        one after another, construct the TATA matrices
c        for all boxes in the structure
c
        do 3000 ibox=2,nboxes
c
c        construct the TATA matrix for this box
c
        call d1m_onebox(ibox,w,box,ab)
c
        idad=box(3)
c
        iboxfrom=idad
        itypefrom=2
c
        iboxto=ibox
        itypeto=2
c
        itypeout=1
c
        call hudson_tata_bld(ier,matfun,iboxfrom,itypefrom,
     1      iboxto,itypeto,itypeout,w,xs,w(ialolo),
     2      nskeldad,nskelbox,eps,
     3      w(iskelfrom),w(iskelto),w(iskelout),w(ievalfrom),w(ievalto),
     4      w(irnorms),w(iwork),w(iipnts_from),w(iipnts_to),
     5      w(iipntsout))
c
c        compress the obtained TATA matrix and store it
c
        itypeout=6
        iboxfrom=idad
        itypefrom=2
        iboxto=ibox
        itypeto=2
c
        call hudson_tata_funny0(ier,itypeout,iboxfrom,itypefrom,
     1      iboxto,itypeto,w,xs,w(ialolo),
     2      nskeldad,nskelbox,w(ialolo_out),w(iww),keep,
     3      w(iskelfrom),w(iskelto),w(iipnts_from),w(iipnts_to),
     4      w(iifrom_into),w(iifrom_notinto),w(iiarr_arr))
c
 3000 continue
c
c        one after another, construct the expand matrices
c        for all childless boxes in the structure
c
        do 4000 ibox=2,nboxes
c
c        construct the expand matrix for this box
c
        call d1m_onebox(ibox,w,box,ab)
        if( (box(4) .gt. 0) .or.(box(5) .gt. 0) ) goto 4000
c
        iboxfrom=ibox
        itypefrom=5
c
        iboxto=ibox
        itypeto=3
c
        itypeout=2
c
        call hudson_lolo_bld(jer,matfun,iboxfrom,itypefrom,
     1      iboxto,itypeto,itypeout,w,xs,w(ialolo),
     2      nskelbox2,nskelbox,eps,
     3      w(iskelfrom),w(iskelto),w(iskelout),w(ievalfrom),w(ievalto),
     4      w(irnorms),w(iwork),w(iipnts_from),w(iipnts_to),
     5      w(iipntsout))
c
c       store the obtained expand matrix
c
        itypeout=17
        itypefrom=5
        itypeto=3
        iboxfrom=ibox
        iboxto=ibox
        call hudson_lolo_funny0(ier,itypeout,iboxfrom,itypefrom,
     1      iboxto,itypeto,w,xs,w(ialolo),nskelbox2,nskelbox,
     2      w(ialolo_out),w(iww),keep,
     3      w(iskelfrom),w(iskelto),w(iipnts_from),w(iipnts_to),
     4      w(iifrom_into),w(iifrom_notinto),w(iiarr_arr))
c
c        construct the LOLO matrix converting the 
c        large expansion into the "standard" one
c
        iboxfrom=ibox
        itypefrom=3
c
        iboxto=ibox
        itypeto=1
c
        itypeout=1
c
        call hudson_lolo_bld(jer,matfun,iboxfrom,itypefrom,
     1      iboxto,itypeto,itypeout,w,xs,w(ialolo),
     2      nskelbox1,nskelbox,eps,
     3      w(iskelfrom),w(iskelto),w(iskelout),w(ievalfrom),w(ievalto),
     4      w(irnorms),w(iwork),w(iipnts_from),w(iipnts_to),
     5      w(iipntsout))
c
c        handle the case of an empty skeleton
c
        if(jer .ne. 0) then
            nskeldad=0
        endif
c
        itypeout=19
        itypefrom=3
        itypeto=1
        iboxfrom=ibox
        iboxto=ibox
c
        call hudson_lolo_funny0(ier,itypeout,iboxfrom,itypefrom,
     1      iboxto,itypeto,w,xs,w(ialolo),nskelbox1,nskelbox,
     2      w(ialolo_out),w(iww),keep,
     3      w(iskelfrom),w(iskelto),w(iipnts_from),w(iipnts_to),
     4      w(iifrom_into),w(iifrom_notinto),w(iiarr_arr))
c
 4000 continue
c
c        one after another, construct the eval matrices
c        for all childless boxes in the structure
c
        do 5000 ibox=2,nboxes
c
c        construct the "big" eval matrix for this box
c
        call d1m_onebox(ibox,w,box,ab)
        if( (box(4) .gt. 0) .or.(box(5) .gt. 0) ) goto 5000
c
        iboxfrom=ibox
        itypefrom=4
c
        iboxto=ibox
        itypeto=5
c
        itypeout=2
c
        call hudson_tata_bld(ier,matfun,iboxfrom,itypefrom,
     1      iboxto,itypeto,itypeout,w,xs,w(ialolo),
     2      nskelbox,nptsbox,eps,
     3      w(iskelfrom),w(iskelto),w(iskelout),w(ievalfrom),w(ievalto),
     4      w(irnorms),w(iwork),w(iipnts_from),w(iipnts_to),
     5      w(iipntsout))
c
c        handle the case of an empty skeleton
c
        if(ier .ne. 0) then
            nskeldad=0
        endif
c
        itypeout=16
        iboxfrom=ibox
        itypefrom=4
        iboxto=ibox
        itypeto=5
c
        call hudson_tata_funny0(ier,itypeout,iboxfrom,itypefrom,
     1      iboxto,itypeto,w,xs,w(ialolo),
     2      nskelbox,nptsbox,w(ialolo_out),w(iww),keep,
     3      w(iskelfrom),w(iskelto),w(iipnts_from),w(iipnts_to),
     4      w(iifrom_into),w(iifrom_notinto),w(iiarr_arr))
c
c        construct the TATA matrix converting the incoming potential 
c        from the regular incoming skeleton to the "big" one
c
        iboxfrom=ibox
        itypefrom=2
c
        iboxto=ibox
        itypeto=4
c
        itypeout=1
c
        call hudson_tata_bld(ier,matfun,iboxfrom,itypefrom,
     1      iboxto,itypeto,itypeout,w,xs,w(ialolo),
     2      nskelbox,nskelbox2,eps,
     3      w(iskelfrom),w(iskelto),w(iskelout),w(ievalfrom),w(ievalto),
     4      w(irnorms),w(iwork),w(iipnts_from),w(iipnts_to),
     5      w(iipntsout))
c
        itypeout=18
        iboxfrom=ibox
        itypefrom=2
        iboxto=ibox
c
        call hudson_tata_funny0(ier,itypeout,iboxfrom,itypefrom,
     1      iboxto,itypeto,w,xs,w(ialolo),
     2      nskelbox,nskelbox2,w(ialolo_out),w(iww),keep,
     3      w(iskelfrom),w(iskelto),w(iipnts_from),w(iipnts_to),
     4      w(iifrom_into),w(iifrom_notinto),w(iiarr_arr))
c
 5000 continue
c
c        construct all list 1 matrices
c
        call hudson_selves_bld(ier,w,nboxes,matfun,
     1      xs,w(iww),w(ialolo),w(iipntsout))
c
        call hudson_lists1_bld(ier,w,nboxes,xs,w(iww),
     1      w(ialolo),w(iskelfrom),w(iskelto),w(iiarr_arr),
     2      w(iipnts_from),w(iipnts_to),matfun)
c
c        construct the list 4 operators
c
        call hudson_lists4_bld(ier,w,nboxes,xs,w(iww),w(ialolo),
     1      w(iskelfrom),w(iskelto),w(iiarr_arr),
     2      w(iipnts_from),w(iipnts_to),matfun)
c
c        construct the list 3 operators
c
        call hudson_lists3_bld(ier,w,nboxes,xs,w(iww),w(ialolo),
     1      w(iskelfrom),w(iskelto),w(iiarr_arr),
     2      w(iipnts_from),w(iipnts_to),matfun)
c
c        construct all lota matrices
c
        call hudson_lotas_bld(ier,w,nboxes,xs,keep7,
     1      w(ialolo),w(iskelfrom),w(iskelto),w(iipnts_from),
     2      w(iipnts_to),matfun,w(iipnts_jbox),w(iww))
c
c        perform garbage collection
c
        lused=iww+keep7
c
        iww2=ialolo
        call d1marrmove(w(iww),w(iww2),keep7)
c
        w(23)=iww2
        keep_out=iww2+keep7
c
        return
        end

c
c
c
c
c
        subroutine hudson_lists4_bld(ier,w,nboxes,xs,ww,
     1      amatr,skelibox,skeljbox,list,ipnts_box,ipnts_jbox,
     2      matfun)
        implicit real *8 (a-h,o-z)
        save
        real *8 w(1),amatr(1),skeljbox(1),skelibox(1),xs(1),ww(1)
        integer *4 list(1),ipnts_box(1),ipnts_jbox(1)
c
c        This subroutine constructs the matrices accounting
c        for the interactions of all boxes with their list 4's. 
c
c                  Input parameters:
c
c  w - the array containing the structure; must have been created via 
c        a preceding call to the subroutine d1mstrcr (see)
c  nboxes - the number of boxes created (one fervently hopes) via 
c        a preceding call to the subroutine d1mstrcr (see)
c  xs - the used-specified nodes in the simulation
c  
c                  Output parameters:
c
c  ww - the array where the operators are returned, to be retrieved
c        from by the subroutine hudson_lota_funny (see)
c
c                  Work arrays:
c
c  amatr,skelibox,skejbox,list,ipnts_box,ipnts_ibox
c
c        apply list4 for all boxes
c
        do 2000 ibox=1,nboxes
c
c        retrieve list 4 for this box
c
        itype=4
        call d1m_onelist(ibox,itype,w,list,nlist)
c
        if(nlist .le. 0) goto 2000
c
c        . . . retrieve the ibox's "regular" outgoing skeleton
c
        itype=1
        call hudson_skel_retr2(jer,w,ibox,itype,xs,
     1      ipnts_box,skelibox,nskelibox)
c
        if(jer .ne. 0) goto 2000
c
c        process the elements of ibox's list 4 one after another
c
        do 1800 i=1,nlist
c
        jbox=list(i)
c
c        retrieve the jbox's "big" incoming skeleton
c
        itype=4
        call hudson_skel_retr2(jer,w,jbox,itype,xs,
     1      ipnts_jbox,skeljbox,nskeljbox)
c
c        evaluate the potential of the regular outgoing skeleton 
c        of ibox at the nodes in the big incoming skeleton of jbox
c
        call hudson_crea_lota(matfun,ipnts_box,skelibox,nskelibox,
     1      ipnts_jbox,skeljbox,nskeljbox,amatr)
c
c        store the constructed matrix
c
        length=nskelibox*nskeljbox
        itypeout=25
        iboxfrom=ibox
        iboxto=jbox
        nskelfrom=nskelibox
        nskelto=nskeljbox
        call hudson_lota_funny0(ier,itypeout,iboxfrom,
     1      iboxto,nskelto,nskelfrom,amatr,ww,keep)
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
        subroutine hudson_lists3_bld(ier,w,nboxes,xs,ww,
     1      amatr,skelibox,skeljbox,list,ipnts_box,
     2      ipnts_ibox,matfun)
        implicit real *8 (a-h,o-z)
        save
        real *8 w(1),amatr(1),skelibox(1),skeljbox(1),xs(1)
        integer *4 list(1),ipnts_box(1),ipnts_ibox(1)
c
c        This subroutine constructs the matrices accounting
c        for the interactions of all boxes with their list 3's. 
c
c                  Input parameters:
c
c  w - the array containing the structure; must have been created via 
c        a preceding call to the subroutine d1mstrcr (see)
c  nboxes - the number of boxes created (one fervently hopes) via 
c        a preceding call to the subroutine d1mstrcr (see)
c  xs - the used-specified nodes in the simulation
c  
c                  Output parameters:
c
c  ww - the array where the operators are returned, to be retrieved
c        from by the subroutine hudson_lota_funny (see)
c
c                  Work arrays:
c
c  amatr,skelibox,skejbox,list,ipnts_box,ipnts_ibox
c
c        . . . apply list1 for all boxes
c
c        apply list3 for all boxes
c
        do 2000 ibox=1,nboxes
c
c        retrieve list 3 for this box
c
        nlist=0
        itype=3
        call d1m_onelist(ibox,itype,w,list,nlist)
c
        if(nlist .le. 0) goto 2000
c
c        . . . retrieve the ibox's "big" skeleton
c
        itype=3
        call hudson_skel_retr2(jer,w,ibox,itype,xs,
     1      ipnts_ibox,skelibox,nskelibox)
c
c        process the elements of ibox's list 3 one after another
c
        do 1800 i=1,nlist

        jbox=list(i)
c
c        . . . retrieve the jbox's "regular" skeleton 
c
        itype=2
        call hudson_skel_retr2(jer,w,jbox,itype,xs,
     1      ipnts_box,skeljbox,nskeljbox)
c
c        construct the interaction matrix
c
        call hudson_crea_lota(matfun,ipnts_ibox,skelibox,nskelibox,
     1      ipnts_box,skeljbox,nskeljbox,amatr)
c
c        store the obtained matrix
c
        length=nskelibox*nskeljbox
        itypeout=27
        iboxfrom=ibox
        iboxto=jbox
        nskelfrom=nskelibox
        nskelto=nskeljbox
        call hudson_lota_funny0(ier,itypeout,iboxfrom,
     1      iboxto,nskelto,nskelfrom,amatr,ww,keep)
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
        subroutine hudson_lists1_bld(ier,w,nboxes,xs,ww,
     1      amatr,skelfrom,skelto,list,ipntsfrom,ipntsto,
     2      matfun)
        implicit real *8 (a-h,o-z)
        save
        real *8 w(1),amatr(1),ab(2),skelfrom(1),
     1      skelto(1),xs(1)
        integer *4 list(1),box(10),ipntsfrom(1),ipntsto(1)
c
c        This subroutine constructs the matrices accounting
c        for the interactions of all childless boxes with their
c        list 1's. Please note that this subroutine DOES NOT
c        deal with self-interactions of childless boxes; for
c        that, see the subroutine hudson_selves_bld.
c
c                  Input parameters:
c
c  w - the array containing the structure; must have been created via 
c        a preceding call to the subroutine d1mstrcr (see)
c  nboxes - the number of boxes created (one fervently hopes) via 
c        a preceding call to the subroutine d1mstrcr (see)
c  xs - the used-specified nodes in the simulation
c  
c                  Output parameters:
c
c  ww - the array where the operators are returned, to be retrieved
c        from by the subroutine hudson_lota_funny (see)
c
c                  Work arrays:
c
c  amatr,skelfrom,skelto,list,ipntsfrom,ipntsto
c
c        . . . apply list1 for all boxes
c
        do 3000 ijk=1,2
c
        do 2000 ibox=1,nboxes
c
        call d1m_onebox(ibox,w,box,ab)
c
        if(box(4) .gt. 0) goto 2000
        if(box(5) .gt. 0) goto 2000
c
c        retrieve list 1 for this box
c
        itype=1
        if(ijk .eq. 2) itype=11
        call d1m_onelist(ibox,itype,w,list,nlist)
c
        if(nlist .le. 0) goto 2000
c
c        retrieve the outgoing "big" skeleton inside ibox
c
        iboxfrom=ibox
        itypefrom=3        
c
        call hudson_skel_retr2(ier,w,iboxfrom,itypefrom,xs,
     1      ipntsfrom,skelfrom,nskelfrom)
c
c        process the elements of ibox's list 1 one after another
c
        do 1800 i=1,nlist
c
        jbox=list(i)
c
c        retrieve the incoming "big" skeleton inside jbox
c
        iboxto=jbox
        itypeto=4
        call hudson_skel_retr2(ier,w,iboxto,itypeto,xs,
     1      ipntsto,skelto,nskelto)
c
c       construct the corresponding "lota" matrix
c
        call hudson_crea_lota(matfun,ipntsfrom,skelfrom,nskelfrom,
     1      ipntsto,skelto,nskelto,amatr)
c
c        store the lota matrix
c
        itypeout=21
c        
        call hudson_lota_funny0(ier,itypeout,iboxfrom,
     1      iboxto,nskelto,nskelfrom,amatr,ww,keep)
c
 1800 continue
c
 2000 continue
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
        subroutine hudson_lolo_bld(ier,matfun,iboxfrom,itypefrom,
     1      iboxto,itypeto,itypeout,w,xs,alolo,
     2      nskelfrom,nskelto,eps,
     3      skelfrom,skelto,skelout,evalfrom,evalto,
     4      rnorms,ww,ipntsfrom,ipntsto,ipntsout)
        implicit real *8 (a-h,o-z)
        save
        dimension skelfrom(1),skelto(1),xs(1),
     1      skelout(1),evalfrom(1),evalto(1),
     2       alolo(1),rnorms(1),ww(1)
        integer *4 w(1),ipntsfrom(1),ipntsto(1),ipntsout(1)
c
c        For the pair of user-specified boxes iboxfrom, iboxto,
c        and the types itypefrom, itypeto of their respective
c        skeletons, this subroutine constructs the LOLO-type 
c        matrix alolo. Actually, this is merely a bookkeeping 
c        interface; all of the hard work is performed by the 
c        subroutine hudson_lolo_bld0 (see) called by this subroutine
c
c                  Input parameters:
c
c  iboxfrom - the box from which the LOLO matrix will map
c  itypefrom - the type of the "from" skeleton from which the
c        matrix will map. The possible values for itypefrom are:
c      itypefrom=1 means the regular outgoing skeleton
c      itypefrom=3 means the "big" outgoing skeleton (only exists 
c        for childless boxes)
c  iboxto - the box to whose skeleton the LOLO matrix will map
c  itypeto - the type of the skeleton on the box iboxto to which the
c        tata matrix will map. The possible values of itypeto are:
c      itypeto=1 means the regular outgoing skeleton
c      itypeto=3 means the "big" outgoing skeleton (only exists 
c        for childless boxes)
c  itypeout - the type of the outer skeleton with respect to which 
c        operator LOLO will be constructed:
c      itypeout=1 means that the outer skeleton will EXCLUDE the
c        list 1 of the relevant box (the "classical" outer skeleton)
c      itypeout .neq. 1 means that the outer skeleton will INCLUDE the
c        list 1 of the relevant box (the "new" outer skeleton, 
c        corresponding to the "big" inner skeleton)
c  w - the array containing the structure; must have been created via 
c        a preceding call to the subroutine d1mstrcr (see)
c  xs - the used-specified nodes in the simulation
c  eps - the accuracy of the skeletonization
c
c                    Output parameters:
c
c  alolo - the LOLO matrix leading from the skeleton of the type 
c        itypefrom in the box iboxfrom to the skeleton of the type
c        itypeto in the box iboxto
c  nskelfrom - the number of elements in the "from" skeleton
c  nskelto - the number of elements in the "to"
c
c                    Work arrays:
c
c  skelfrom,skelto,skelout,evalfrom,evalto,rnorms,ww,ipntsfrom,
c        ipntsto
c
c        retrieve the "from" skeleton
c
        ier=0
c
        call hudson_skel_retr2(jer,w,iboxfrom,itypefrom,xs,
     1      ipntsfrom,skelfrom,nskelfrom)
c
        if(nskelfrom .eq. 0) then
            ier=8
            return
        endif
c
c        retrieve the "to" skeleton
c
        call hudson_skel_retr2(jer2,w,iboxto,itypeto,xs,
     1      ipntsto,skelto,nskelto)
c
        if(nskelto .eq. 0) then
            ier=4
            return
        endif
c
c       depending n the parameter itypeout, construct
c       the outer skeleton
c
        list_type=2
        ipts_type=2
        itype_pts=2
        call hudson_pnts_inlists(iboxto,w,list_type,ipts_type,
     1      nts,skelout,nskelout,xs,ipntsout,itype_pts,
     2      rnorms,alolo)
c
        list_type=4
        ipts_type=1
c
        call hudson_pnts_inlists(iboxto,w,list_type,ipts_type,
     1      nts,skelout(nskelout+1),nskel4,xs,
     2      ipntsout(nskelout+1),itype_pts,rnorms,alolo)
c
        nskelout=nskelout+nskel4
c
        if(itypeout .eq. 1) goto 2000
c
        list_type=3
        ipts_type=2
c
        call hudson_pnts_inlists(iboxto,w,list_type,ipts_type,
     1      nts,skelout(nskelout+1),nskel4,xs,
     2      ipntsout(nskelout+1),itype_pts,rnorms,alolo)
c
        nskelout=nskelout+nskel4
c
        call hudson_list1_pts(iboxto,w,xs,ipntsout(nskelout+1),
     1      skelout(nskelout+1),nptsout,rnorms)
c
        nskelout=nptsout+nskelout
c
 2000 continue
c
        if(nskelout .eq. 0) then
            ier=4
            return
        endif  
c
c       construct the matrix lolo connecting iboxfrom with
c       iboxto
c
        call hudson_lolo_bld0(matfun,skelfrom,nskelfrom,
     1     skelto,nskelto,skelout,nskelout,
     2      evalfrom,evalto,alolo,eps,rnorms,ww,
     3      ipntsfrom,ipntsto,ipntsout)
c
        return
        end
c     
c
c
c
c
        subroutine hudson_lolo_bld0(matfun,skelfrom,nskelfrom,
     1     skelto,nskelto,skelout,nskelout,
     2      evalfrom,evalto,alolo,eps,rnorms,ww,
     3      ipntsfrom,ipntsto,ipntsout)
        implicit real *8 (a-h,o-z)
        save
        dimension skelfrom(1),skelto(1),skelout(1),
     1      evalfrom(nskelout,nskelfrom),rnorms(1),ww(1),
     2      evalto(nskelout,nskelto),
     3      alolo(nskelto,nskelfrom),
     4      ipntsfrom(1),ipntsto(1),ipntsout(1)
c
c        For the user-specified box iboxfrom, this subroutine 
c        constructs the LOLO matrix alolo leading FROM the box 
c        iboxfrom to the box iboxto
c
c                  Input parameters:
c
c  skelfrom - the skeleton of the box iboxfrom (actual lications in R^1)
c  nskelfrom - the number of elements in skelfrom
c  skelto - the skeleton of iboxto (actual lications in R^1)
c  nskelto - the number of elements in skelto
c  skelout - outer skeleton of the box iboxto (actual lications in R^1)
c  nskelout - number of elements in skelout
c  eps - the accuracy of the skeletonization
c  ipntsfrom - the skeleton of the box iboxfrom (locations of 
c        nodes in array xs)
c  ipntsto - the skeleton of the box iboxto (locations of 
c        nodes in array xs)
c  ipntsout - the outer skeleton of the box iboxto (locations of 
c        nodes in array xs)
c
c                    Output parameters:
c
c  alolo - the LOLO matrix leading from the box ibox to his dad
c
c                    Work arrays:
c
c  evalfrom,evalto,rnorms,ww
c
c
c        . . . construct the matrix evalfrom
c
        do 1400 i=1,nskelfrom
        do 1200 j=1,nskelout
c
        ix=ipntsfrom(i)
        iy=ipntsout(j)

        call matfun(skelfrom(i),ix,skelout(j),iy,d)
        evalfrom(j,i)=d
 1200 continue
 1400 continue
c
c        construct the matrix evaldad
c
        do 2400 i=1,nskelto
        do 2200 j=1,nskelout
c
        ix=ipntsto(i)
        iy=ipntsout(j)

        call matfun(skelto(i),ix,skelout(j),iy,d)
c
        evalto(j,i)=d        
 2200 continue
 2400 continue
c
c       use least squares to construct the LOLO matrix
c
        ifcheck=0
        call nrleamatll(evalto,nskelout,nskelto,nskelfrom,
     1      evalfrom,alolo,eps,ncols,
     2    rnorms,ww,ifcheck,errl2,errmax)
c
        return
        end
c
c
c
c
c
        subroutine hudson_tata_funny0(ier,itypeout,iboxfrom,itypefrom,
     1      iboxto,itypeto,w,xs,atata,
     2      nskelfrom,nskelto,atata_out,ww,keep,
     3      skel_from,skel_to,iskel_from,iskel_to,
     4      ito_infrom,ito_notinfrom,iarr)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),iskel_from(1),iskel_to(1),
     1      ito_infrom(2,1),iarr(1),xs(1),
     2      atata(nskelto,nskelfrom),atata_out(nskelfrom,1),
     3      ito_notinfrom(1),skel_from(1),skel_to(1)
c
        integer *4 ww(1)
c
c        This subroutine examines the user-specified TATA-type
c        matrix, determines which of its columns can be replaced
c        with columns of the identity matrix, constructs tables
c        to be used to control the accelerated application of
c        the said matrix to vectors, and stores the whole mess
c        of them in array ww, to be used by the subroutine 
c        hudson_tata_funny (see).
c
c                  Input parameters:
c
c  itypeout - the type of  the matrix to be stored. The possible 
c        types are as follows:
c     itype=6 denotes a TATA matrix where both skeletons are "standard"
c     itype=18 denotes a TATA matrix converting the incoming potential
c        on a regular skeleton into incoming expansion on the "big"
c        skeleton    
c     itype=16 denotes a TATA-type matrix, converting a "big" incoming
c        expansion into the potential at all nodes in the box
c  iboxfrom - the box from which the LOLO-type  matrix will map
c  itypefrom - the type of the "from" skeleton from which the
c        matrix will map. The possible values for itypefrom are:
c     itypefrom=1 means the regular outgoing skeleton
c     itypefrom=3 means the "big" outgoing skeleton (only exists 
c        for childless boxes)
c     itypefrom=5 means the actual nodes on the interval (as opposed
c        to a skeleton, which is a subset of the said nodes)
c  iboxto - the box to whose skeleton the LOLO-type matrix will map
c  itypeto - the type of the skeleton on the box iboxto to which the
c        tata matrix will map. The possible values of itypeto are:
c     itypeto=1 means the regular outgoing skeleton
c     itypeto=3 means the "big" outgoing skeleton (only exists 
c        for childless boxes)
c  w - the array containing the structure; must have been created via 
c        a preceding call to the subroutine d1mstrcr (see)
c  xs - the used-specified nodes in the simulation
c  alolo - the LOLO matrix leading from the skeleton of the type 
c        itypefrom in the box iboxfrom to the skeleton of the type
c        itypeto in the box iboxto
c
c                    Output parameters:
c
c  ier - error return code
c  ww - the array where the compressed matrix has been stored by the
c        subroutine (among other matrices)
c
c                    work arrays:
c
c  alolo_out,skel_from,skel_to,iskel_from,iskel_to,ifrom_into,
c        ifrom_notinto,iarr
c
c
c        . . . retrive both "from" and "to" skeletons
c
        call hudson_skel_retr2(jer,w,iboxfrom,itypefrom,xs,
     1      iskel_from,skel_from,nskelfrom)
c
        call hudson_skel_retr2(jer2,w,iboxto,itypeto,xs,
     1      iskel_to,skel_to,nskelto)
c
        lenint=0
        lenrea=0
c
        if( (jer .ne. 0) .or. (jer2 .ne. 0) ) goto 2800        
c
c        construct the list of nodes inside the "to"
c        skeleton that are also elements in the "from"
c        skeleton
c
        ii=0
        iii=0
c
        do 1400 i=1,nskelto
c
        do 1200 j=1,nskelfrom
c
        if(iskel_to(i) .ne. iskel_from(j)) goto 1200
c
        ii=ii+1
        ito_infrom(1,ii)=i
        ito_infrom(2,ii)=j
        goto 1400
 1200 continue
c
        iii=iii+1
 1400 continue
c
        nn=ii
        nto_infrom=ii
        nnn=iii
c
        ito_infrom(1,nn+1)=-7
        ito_infrom(2,nn+1)=-7
c
c        separate the part of the tata matrix that will 
c        be used (as opposed to the part that will be
c        replaced with identity)
c
        ii=1
        i1=0
        do 1800 i=1,nskelto
c
        if (i .eq. ito_infrom(1,ii)) then
            ii=ii+1
            goto 1800
        endif
c
        i1=i1+1
        do 1600 j=1,nskelfrom
c
        atata_out(j,i1)=atata(i,j)
 1600 continue
c
 1800 continue
c
        nto_notinfrom=i1
c
c        repack the potentials on the sonnie's skeleton
c
        ii=1
        iii=1
        do 2200 i=1,nskelto
c
        if(ito_infrom(1,ii) .eq. i) then
            jj=ito_infrom(2,ii) 
            ii=ii+1
            goto 2200
        endif
c
        ito_notinfrom(iii)=i
        iii=iii+1
c
 2200 continue
c
c        prepare the integer array to be stored by the 
c        subroutine hudson_linsto
c
        iarr(1)=nskelfrom
        iarr(2)=nskelto
        iarr(3)=nto_infrom
        iarr(4)=nto_notinfrom
c
        iito_infrom=11
        lito_infrom=(nto_infrom+2)*2
c
        iito_notinfrom=iito_infrom+lito_infrom
        lito_notinfrom=nto_notinfrom
c
        iarr(5)=iito_infrom
        iarr(6)=iito_notinfrom
c
        jj=0
        do 2400 i=1,nto_infrom
c
        iarr(iito_infrom+jj)=ito_infrom(1,i)
        jj=jj+1
        iarr(iito_infrom+jj)=ito_infrom(2,i)
        jj=jj+1
c
 2400 continue
c
        do 2600 i=1,nto_notinfrom
c
        iarr(iito_notinfrom+i-1)=ito_notinfrom(i)
 2600 continue
c
        lenint=iito_notinfrom+lito_notinfrom
        lenrea=nskelfrom*nto_notinfrom
c
c        store the whole thing in array ww
c
 2800 continue
c
        key=0
        call hudson_linsto(ier,iboxfrom,itypeout,iboxto,key,
     1      iarr,lenint,atata_out,lenrea,
     2      ww,keep)
c
        return
        end
c
c
c
c
c
        subroutine hudson_lolo_funny0(ier,itypeout,iboxfrom,
     1      itypefrom,iboxto,itypeto,w,xs,alolo,nskelfrom,nskelto,
     2      alolo_out,ww,keep,
     3      skel_from,skel_to,iskel_from,iskel_to,
     4      ifrom_into,ifrom_notinto,iarr)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),iskel_from(1),iskel_to(1),
     1      ifrom_into(2,1),ifrom_notinto(1),
     2      alolo(nskelto,nskelfrom),alolo_out(nskelto,1),
     3      iarr(1),xs(1),skel_from(1),skel_to(1)
        integer *4 ww(1)
c
c        This subroutine examines the user-specified LOLO-type
c        matrix, determines which of its rows can be replaced
c        with rows of the identity matrix, constructs tables
c        to be used to control the accelerated application of
c        the said matrix to vectors, and store the whole mess
c        of them in array ww, to be used by the subroutine 
c        hudson_lolo_funny (see).
c
c                  Input parameters:
c
c  itypeout - the type of  the matrix to be stored. The possible 
c        types are as follows:
c     itype=5 denotes a LOLO matrix where both skeletons are "standard"
c     itype=17 denotes an EXPAND matrix converting nodes on the 
c        interval "ibox" into the "big" outgoing expansion
c     itype=19 denotes a "LOLO" matrix, converting a "big" outgoing
c        expansion into the "standard" outgoing expansion
c  iboxfrom - the box from which the LOLO-type  matrix will map
c  itypefrom - the type of the "from" skeleton from which the
c        matrix will map. The possible values for itypefrom are:
c     itypefrom=1 means the regular outgoing skeleton
c     itypefrom=3 means the "big" outgoing skeleton (only exists 
c        for childless boxes)
c     itypefrom=5 means the actual nodes on the interval (as opposed
c        to a skeleton, which is a subset of the said nodes)
c  iboxto - the box to whose skeleton the LOLO-type matrix will map
c  itypeto - the type of the skeleton on the box iboxto to which the
c        tata matrix will map. The possible values of itypeto are:
c     itypeto=1 means the regular outgoing skeleton
c     itypeto=3 means the "big" outgoing skeleton (only exists 
c        for childless boxes)
c  w - the array containing the structure; must have been created via 
c        a preceding call to the subroutine d1mstrcr (see)
c  xs - the used-specified nodes in the simulation
c  alolo - the LOLO matrix leading from the skeleton of the type 
c        itypefrom in the box iboxfrom to the skeleton of the type
c        itypeto in the box iboxto
c
c                    Output parameters:
c
c  ier - error return code
c  ww - the array where the compressed matrix has been stored by the
c        subroutine (among other matrices)
c
c                    work arrays:
c
c  alolo_out,skel_from,skel_to,iskel_from,iskel_to,ifrom_into,
c        ifrom_notinto,iarr)
c
c
c        . . . retrieve skeletons for both source and target
c
        call hudson_skel_retr2(ier,w,iboxfrom,itypefrom,xs,
     1      iskel_from,skel_from,nskelfrom)
c
        call hudson_skel_retr2(ier,w,iboxto,itypeto,xs,
     1      iskel_to,skel_to,nskelto)
c
        lenrea=0
        lenint=0
        if( (nskelfrom .eq. 0) .or. (nskelto .eq. 0) ) goto 2400
c
c        construct the list of nodes inside the source's
c        skeleton that are also elements in the target's
c        skeleton
c
        ii=0
        do 1400 i=1,nskelfrom
c
        do 1200 j=1,nskelto
c
        if(iskel_from(i) .eq. iskel_to(j)) goto 1300
c
 1200 continue
c
        goto 1400
 1300 continue
c
        ii=ii+1
        ifrom_into(1,ii)=i
        ifrom_into(2,ii)=j
        goto 1400
 1400 continue
c
        nn=ii
        nfrom_into=ii
c
        ifrom_into(1,nn+1)=-7
        ifrom_into(2,nn+1)=-7
c
c        separate the part of the lolo matrix that will 
c        be used (as opposed to the part that will be
c        replaced with identity)
c
        ii=1
        i1=0
        do 1800 i=1,nskelfrom
c
        if (i .eq. ifrom_into(1,ii)) then
            ii=ii+1
            goto 1800
        endif
c
        i1=i1+1
        do 1600 j=1,nskelto
c
        alolo_out(j,i1)=alolo(j,i)
 1600 continue
c
        ifrom_notinto(i1)=i
 1800 continue
c
        nfrom_notinto=i1
c
c        prepare the integer array to be stored by the 
c        subroutine hudson_linsto
c
        iifrom_into=11
        lifrom_into=(nfrom_into+2)*2
c
        iifrom_notinto=iifrom_into+lifrom_into
        lifrom_notinto=nfrom_notinto
c
        iarr(1)=nskelfrom
        iarr(2)=nskelto
        iarr(3)=nfrom_into
        iarr(4)=nfrom_notinto
c
        iarr(5)=iifrom_into
        iarr(6)=iifrom_notinto
c
        jj=0
        do 2000 i=1,nfrom_into
c
        iarr(iifrom_into+jj)=ifrom_into(1,i)
        jj=jj+1
        iarr(iifrom_into+jj)=ifrom_into(2,i)
        jj=jj+1
c
 2000 continue
c
        do 2200 i=1,nfrom_notinto
c
        iarr(iifrom_notinto+i-1)=ifrom_notinto(i)
 2200 continue
c
        lenint=iifrom_notinto+nfrom_notinto+1
        lenrea=nskelto*nfrom_notinto+2
c
c        store the whole thing in array ww
c
 2400 continue
c
        key=0
        call hudson_linsto(ier,iboxfrom,itypeout,iboxto,key,
     1      iarr,lenint,alolo_out,lenrea,
     2      ww,keep)
c
        return
        end
c
c
c
c
c
        subroutine hudson_selves_bld(ier,w,nboxes,matfun,
     1      xs,ww,amatr,iskelibox)
        implicit real *8 (a-h,o-z)
        save
        real *8 w(1),amatr(1),ab(2),xs(1)
        integer *4 box(10),iskelibox(1)
c
c        This subroutine constructs the matrices accounting
c        for the self-interactions of all childless boxes
c
c                  Input parameters:
c
c  w - the array containing the structure; must have been created via 
c        a preceding call to the subroutine d1mstrcr (see)
c  nboxes - the number of boxes created (one fervently hopes) via 
c        a preceding call to the subroutine d1mstrcr (see)
c  xs - the used-specified nodes in the simulation
c  
c                  Output parameters:
c
c  ww - the array where the operators are returned, to be retrieved
c        from by the subroutine hudson_lota_funny (see)
c
c                  Work arrays:
c
c  amatr - 
c
c        construct list1 interactions for all boxes
c
        do 2000 ibox=1,nboxes
c
        call d1m_onebox(ibox,w,box,ab)
c
        if(box(4) .gt. 0) goto 2000
        if(box(5) .gt. 0) goto 2000
c
c        retrieve the ibox's nodes 
c
        i1_ibox=box(6)
        npts_ibox=box(7)
c
        do 2200 i=1,npts_ibox
c
        iskelibox(i)=i1_ibox-1+i
 2200 continue
c
c        process the elements of ibox's list 1 one after another
c
        jbox=ibox
c
c        retrieve the nodes and potentials inside jbox
c
        call hudson_crea_lota(matfun,iskelibox,xs(i1_ibox),npts_ibox,
     1      iskelibox,xs(i1_ibox),npts_ibox,amatr)
c
c       store the obtained list 1 matrix
c
        length=npts_ibox**2
        itypeout=23
        iboxfrom=ibox
        iboxto=ibox
        nskelto=npts_ibox
        nskelfrom=npts_ibox
        call hudson_lota_funny0(ier,itypeout,iboxfrom,
     1      iboxto,nskelto,nskelfrom,amatr,ww,keep)
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
        subroutine hudson_lota_funny0(ier,itypeout,iboxfrom,
     1      iboxto,nskelto,nskelfrom,alota,ww,keep)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),alota(nskelto,nskelfrom)
        integer *4 ww(1),iarr(4)
c
c        This subroutine stores in array ww a single matrix of 
c        type LOTA; the said matrix can be retrieved from the
c        said state via the subroutine hudson_lota_funny (see)
c
c
c                  Input parameters:
c
c  itypeout - the type of  the matrix to be stored. The possible 
c        types are as follows:
c  iboxfrom - the box from which the LOLO-type  matrix will map
c  iboxto - the box to whose skeleton the LOLO-type matrix will map
c  nskelto - the number of elements in the incoming skeleton 
c       of the box iboxto
c  nskelfrom - the number of elements in the outgoing skeleton 
c       of the box iboxfrom
c  alota - the LOTA matrix leading from the skeleton of the type 
c        itypefrom in the box iboxfrom to the skeleton of the type
c        itypeto in the box iboxto
c
c                    Output parameters:
c
c  ier - error return code
c  ww - the array where the compressed matrix has been stored by the
c        subroutine (among other matrices)
c  keep - the number of integer *4 elements in array ww that have to
c        be kept unchanged between this call to this subroutine and
c        the subsequent calls to the subroutine hudson_lota_funny
c
c
c        . . . prepare the integer array to be stored by the 
c              subroutine hudson_linsto
c
        iarr(1)=nskelfrom
        iarr(2)=nskelto
c
        lenint=3
        lenrea=nskelfrom*nskelto
c
c        store the whole thing in array ww
c
        key=0
        call hudson_linsto(ier,iboxfrom,itypeout,iboxto,key,
     1      iarr,lenint,alota,lenrea,
     2      ww,keep)
c
        return
        end
c
c
c
c
c
        subroutine hudson_lotas_bld(ier,w,nboxes,xs,keep,
     1      alota,skelibox,skeljbox,list,ipnts_box,matfun,
     2      ipnts_jbox,ww)
        implicit real *8 (a-h,o-z)
        save
        real *8 w(1),alota(1),skelibox(1),skeljbox(1),xs(1)
        integer *4 list(1),ipnts_box(1),ipnts_jbox(1),ww(1)
c
c        This subroutine constructs all of the "classical"
c        LOTA matrices, and stores them in the array w3,
c        from which said array it can be retrieved via the 
c        subroutine linkret (see)
c
c                 Input parameters:
c
c  w - the array containing the structure; must have been created via 
c        a preceding call to the subroutine d1mstrcr (see)
c  nboxes - the number of boxes created (one fervently hopes) via 
c        a preceding call to the subroutine d1mstrcr (see)
c  xs - the used-specified nodes in the simulation
c
c                    Output parameters:
c  
c  ier - error return code
c  keep - the length of the piece of the array w3 to be kept between
c        the call to this subroutine and subsequent calls to the 
c        subroutine hudson_all_lotas (see) that uses the data stored in w3
c  w3 - storage area where all of the LOTA matrices have been stored;
c        must have been initialized via a preceding call to linkini 
c        (see)
c
c                    work arrays:
c
c  alota,skelibox,skeljbox,list,ipnts_box
c  
c
c        . . . construct all lota matrices, and store them things
c
        do 2000 ibox=1,nboxes
c
c        retrieve list 2 for this box
c
        itype=2
        call d1m_onelist(ibox,itype,w,list,nlist)
c
        if(nlist .le. 0) goto 2000
c
c        retrieve ibox's outgoing skeleton
c
        itype=1
        call hudson_skel_retr2(ier,w,ibox,itype,xs,
     1      ipnts_box,skelibox,nskelibox)
c
c        process ibox's list 2
c
        do 1600 i=1,nlist
c
        jbox=list(i)
        itype=2
        call hudson_skel_retr2(jer,w,jbox,itype,xs,
     1      ipnts_jbox,skeljbox,nskeljbox)
c
        call hudson_crea_lota(matfun,ipnts_box,skelibox,nskelibox,
     1      ipnts_jbox,skeljbox,nskeljbox,alota)
c
c       store the obtained lota matrix
c
        itypeout=33
        iboxfrom=ibox
        iboxto=jbox
        nskelfrom=nskelibox
        nskelto=nskeljbox
        call hudson_lota_funny0(ier,itypeout,iboxfrom,
     1      iboxto,nskelto,nskelfrom,alota,ww,keep)
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
        subroutine hudson_tata_bld(ier,matfun,iboxfrom,itypefrom,
     1      iboxto,itypeto,itypeout,w,xs,atata,
     2      nskelfrom,nskelto,eps,
     3      skelfrom,skelto,skelout,evalfrom,evalto,
     4      rnorms,ww,ipnts_from,ipnts_to,ipntsout)
        implicit real *8 (a-h,o-z)
        save
        dimension skelfrom(1),skelto(1),xs(1),skelout(1),
     1      evalfrom(1),evalto(1),atata(1),rnorms(1),ww(1)
        integer *4 w(1),ipnts_from(1),ipnts_to(1),ipntsout(1)
c
c        For the pair of user-specified boxes iboxfrom, iboxto,
c        and the types itypefrom, itypeto of their respective
c        skeletons, this subroutine constructs the TATA-type 
c        matrix atata. Actually, this is merely a bookkeeping 
c        interface; all of the hard work is performed by the 
c        subroutine hudson_tata_bld0 (see) called by this subroutine
c
c                  Input parameters:
c
c  iboxfrom - the box from which the tata matrix will map
c  itypefrom - the type of the "from" skeleton from which the
c        matrix will map. The possible values for itypefrom are:
c      itypefrom=2 means the regular incoming skeleton
c      itypefrom=4 means the "big" incoming skeleton (only exists 
c        for childless boxes)
c  iboxto - the box to whose skeleton the tata matrix will map
c  itypeto - the type of the skeleton on the box iboxto to which the
c        tata matrix will map. The possible values of itypeto are:
c      itypeto=2 means the regular incoming skeleton
c      itypeto=4 means the "big" incoming skeleton (only exists 
c        for childless boxes)
c  itypeout - the type of the outer skeleton with respect to which 
c        operator TATA will be constructed:
c      itypeout=1 means that the outer skeleton will EXCLUDE the
c        list 1 of the relevant box (the "classical" outer skeleton)
c      itypeout .neq. 1 means that the outer skeleton will INCLUDE the
c        list 1 of the relevant box (the "new" outer skeleton, 
c        corresponding to the "big" inner skeleton)
c  w - the array containing the structure; must have been created via 
c        a preceding call to the subroutine d1mstrcr (see)
c  xs - the used-specified nodes in the simulation
c  eps - the accuracy of the skeletonization
c
c                    Output parameters:
c
c  atata - the TATA matrix leading from the skeleton of the type 
c        itypefrom in the box iboxfrom to the skeleton of the type
c        itypeto in the box iboxto
c  nskelfrom - the number of elements in the "from" skeleton
c  nskelto - the number of elements in the "to"
c
c                    Work arrays:
c
c  skelfrom,skelto,skelout,evalfrom,evalto,rnorms,ww,ipnts_from,
c        ipnts_to
c
c       retrieve the "to" skeleton
c
        ier=0
c
        call hudson_skel_retr2(jer,w,iboxto,itypeto,xs,
     1      ipnts_to,skelto,nskelto)
c
        if(jer .ne. 0) then
            ier=4
            return
        endif
c
        if(nskelto .eq. 0) then
            ier=4
            return
        endif
c
c        retrieve the "from" skeleton
c
        call hudson_skel_retr2(jer,w,iboxfrom,itypefrom,xs,
     1      ipnts_from,skelfrom,nskelfrom)
c
        if(jer .ne. 0) then
            ier=4
            return
        endif
c
        if(nskelfrom .eq. 0) then
            ier=4
            return
        endif
c
c        retrieve the outer skeleton
c
        list_type=2
        ipts_type=2
        itype_pts=1
        call hudson_pnts_inlists(iboxfrom,w,list_type,ipts_type,
     1      nts,skelout,nskelout,xs,ipntsout,itype_pts,
     3      rnorms,atata)
c
        list_type=4
        ipts_type=1
c
        call hudson_pnts_inlists(iboxfrom,w,list_type,ipts_type,
     1      nts,skelout(nskelout+1),nskel4,xs,
     2      ipntsout(nskelout+1),itype_pts,
     3      rnorms,atata)
c
        nskelout=nskelout+nskel4
c
        if(itypeout .eq. 1) goto 2000
c
        list_type=3
        ipts_type=2
c
        call hudson_pnts_inlists(iboxfrom,w,list_type,ipts_type,
     1      nts,skelout(nskelout+1),nskel3,xs,
     2      ipntsout(nskelout+1),itype_pts,
     3      rnorms,atata)
c
        nskelout=nskelout+nskel3
c
        call hudson_list1_pts(iboxto,w,xs,ipntsout(nskelout+1),
     1      skelout(nskelout+1),nptsout,rnorms)
        nskelout=nptsout+nskelout
c
 2000 continue
c
c       construct the matrix tata for this chunk
c

        call hudson_tata_bld0(matfun,skelto,nskelto,
     1     skelfrom,nskelfrom,skelout,nskelout,
     2      evalfrom,evalto,atata,rnorms,ww,eps,
     3      ipnts_from,ipnts_to,ipntsout)
c
        return
        end
c     
c
c
c
c
        subroutine hudson_tata_bld0(matfun,skelto,nskelto,
     1     skelfrom,nskelfrom,skel_out,nskel_out,
     2      evalto,evalfrom,atata,rnorms,ww,eps,
     3      ipntsfrom,ipntsto,ipntsout)
        implicit real *8 (a-h,o-z)
        save
        dimension skelto(1),skelfrom(1),
     1      skel_out(1),evalto(nskelto,nskel_out),
     2      evalfrom(nskelfrom,nskel_out),
     3      atata(nskelto,nskelfrom),rnorms(1),
     4      ww(1),ipntsfrom(1),ipntsto(1),ipntsout(1)
c
c        For the user-specified box ibox, this subroutine 
c        constructs the TATA matrix atata leading TO the box 
c        iboxto from the box iboxfrom
c
c                  Input parameters:
c
c  skelto - the skeleton of iboxto
c  nskelto - the number of elements in skelto
c  skelfrom - the skeleton of the box ibox
c  nskelfrom - the number of elements in skelfrom
c  skel_out - outer skeleton of the box iboxfrom
c  nskel_out - number of elements in skel_out
c  eps - the accuracy of the skeletonization
c
c                    Output parameters:
c
c  lolo - the TATA matrix leading TO the box ibox from his dad
c
c                    Work arrays:
c
c  evalto,evalfrom,rnorms,ww
c
c
c        . . . construct the matrix evalto
c
        do 1400 i=1,nskelto
        do 1200 j=1,nskel_out
c
        iy=ipntsto(i)
        ix=ipntsout(j)
c
        call matfun(skel_out(j),ix,skelto(i),iy,d)
        evalto(i,j)=d
 1200 continue
 1400 continue
c
c        construct the matrix evalfrom
c
        do 2400 i=1,nskelfrom
        do 2200 j=1,nskel_out
c
        iy=ipntsfrom(i)
        ix=ipntsout(j)
c
        call matfun(skel_out(j),ix,skelfrom(i),iy,d)
        evalfrom(i,j)=d
 2200 continue
 2400 continue
c
c       use least squares to construct the TATA matrix
c
        ifcheck=0
c
        call nrleamatrr(evalfrom,nskelto,nskelfrom,nskel_out,
     1      evalto,atata,eps,ncols,
     2    rnorms,ww,ifcheck,errl2,errmax)
c
        return
        end
c
c
c
c
c
        subroutine hudson_crea_lota(matfun,iskelibox,skelibox,
     1      nskelibox,iskeljbox,skeljbox,nskeljbox,alota)
        implicit real *8 (a-h,o-z)
        save
        dimension alota(nskeljbox,nskelibox),
     1      skelibox(1),skeljbox(1),iskelibox(1),iskeljbox(1)
c
c        This subroutine constructs the matrix of interactions
c        between the elements of two user-specified sets of nodes
c
c                  Input parameters:
c
c  matfun - the user-supplied subroutine evaluating the interactions
c        to be "accelerated"
c  iskelibox - the addresses in array xs of the sources 
c  skelibox - the location in R^1 of the sources 
c  nskelibox - the number of sources
c  iskeljbox - the addresses in array xs of the targets
c  skeljbox - the location in R^1 of the targets
c  nskeljbox - the number of targets
c
c                  Output parameters:
c
c  alota - the matrix of interactions
c
c        . . . construct lota for this pair of chunks
c
        do 1400 i=1,nskelibox
        do 1200 j=1,nskeljbox
c
        x=skelibox(i)
        ix=iskelibox(i)
c
        y=skeljbox(j)
        iy=iskeljbox(j)
c
        call matfun(x,ix,y,iy,alota(j,i))
c
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
        subroutine hudson_allskels(ier,w,xs,nboxes,nts,eps,
     1      ifrelskel,addr,arr,lused,
     2      skel,ipts_inbox,ipivots,list,lists,
     3      amatr,lamatr,pts_inbox,skelout,rnorms,ncols_max)
        implicit real *8 (a-h,o-z)
        save
        integer *4 w(1),box(10),skel(1),ipts_inbox(1),
     1      addr(8,1),arr(1),ipivots(1),list(1),lists(1)
        real *8 ab(2),amatr(1),pts_inbox(1),
     1      skelout(1),rnorms(1)
c
c        This subroutine constructs the preliminary skeletons 
c        (both incoming and outgoing) for all boxes on all levels.
c        The skeletons are preliminary in the sense that both 
c        the OUTER skeletons and the kernel functions used to 
c        construct them are artificial. The output of this 
c        subroutine is used by the subroutine hudson_allskels2 (see) to
c        construct the final ("true") skeletons for all boxes on 
c        all levels.
c
c                    Input parameters:
c
c  w - the array containing the structure; must have been created via 
c        a preceding call to the subroutine d1mstrcr (see)
c  xs - the used-specified nodes in the simulation
c  nboxes - the total number of boxes in the structure
c  nts - the number of chebychev nodes in the auxiliary discretization
c        of each box in the list 2 of the box being skeletonized
c  eps - the accuracy of the skeletonization
c
c                    Output parameters:
c
c  ier - error return code
c  addr, arr, lused - parameters used by the subroutine hudson_skel_store0 (see)
c  ncols_max - the maximum number of elements in any skeleton for any
c        box on any level
c
c                    Work arrays:
c
c  skel - must be at least n (number of nodes in the simulation)
c         real *8 locations long
c  ipts_inbox - must be at least max(nmax,nts) integer *4 elements, where
c        nmax is the maximum number of points in a childless box
c  ipivots - must be at least n (the number of points in the simulation)
c        integer *4 locations long
c  list - must be at least nlists2_max_nlists4_max integer *4 locations 
c        long, where nlists2_max is the maximum (among all boxes ibox 
c        in the structure) total number of elements in the lists2 of ibox 
c        and all its ancestors. Similarly, nlists4_max is the maximum 
c        (among all boxes ibox in the structure) total number of elements 
c        in the lists4 of ibox and all its ancestors
c  lists - must be at least nlists2_max_nlists4_max integer *4 locations 
c        long, where nlists2_max is the maximum (among all boxes ibox 
c        in the structure) total number of elements in the lists2 of ibox 
c        and all its ancestors. Similarly, nlists4_max is the maximum 
c        (among all boxes ibox in the structure) total number of elements 
c        in the lists4 of ibox and all its ancestors
c  amatr - must be at least (nlists2_max*nts+npts4_max)*nskelmax real *8 
c        elements long
c  pts_inbox - must be at least max(nmax,nts) real *8 elements long, 
c        where nmax is the maximum number of points in a childless box
c  skelout - must be at least n (the number of points in the simulation)
c        real *8 locations long
c  rnorms - must be at least n (the number of points in the simulation)
c        real *8 locations long
c
c        . . . skeletonize all boxes without kids
c
        nn=nboxes
c
        call hudson_skel_store0init(addr,nn)
c        
        ncols_max=0
        do 1400 i=1,nn
c
        ibox=i
c
        call hudson_nokids(jer,ibox,w,xs,nts,eps,ifrelskel,
     1      addr,arr,lused,
     2      ipts_inbox,skel,ipivots,list,lists,pts_inbox,amatr,lamatr,
     3      skelout,rnorms,ncols)
c
        if(ncols_max .lt. ncols) ncols_max=ncols
 1400 continue
c
c        Scan all boxes. Whenever a box is encountered whose
c        both kids have been skeletonized, skeletonize the
c        said box
c
        do 4000 ijk=1,100
c
        ifacted=0
c
        do 3000 i=1,nboxes
c
        if(addr(1,i) .gt. 0) goto 3000
c
        ibox=i
c
c        if this box is childless - skip it
c
        call d1m_onebox(ibox,w,box,ab)
c
        ison1=box(4)
        ison2=box(5)
c
        if( (ison1 .le. 0) .and. (ison2 .le. 0) ) goto 3000
c
c        this box has kids. check if both have been skeletonized
c
        ifacted=1
        nsons=0
        ifbadsons=0
c        
        if(ison1 .gt. 0) then
            nsons=nsons+1
            if(addr(1,ison1) .le. 0) ifbadsons=1
        endif
c
c        
        if(ison2 .gt. 0) then
            nsons=nsons+1
            if(addr(1,ison2) .le. 0) ifbadsons=1
        endif
c
        if(ifbadsons .eq. 1) goto 3000
c
c        This box has at least one kid, and no kids that have not
c        been skeletonized. If this box has only one kid, copy
c        the kid's skeleton into daddy, and be done
c
        if(nsons. ne. 1) goto 2000
c
        if(ison1 .gt. 0) then
            itype=1
            call hudson_skel_retr0(ier,ison1,itype,skel,nskel,
     1          addr,arr)
            call hudson_skel_store0(ibox,itype,skel,nskel,addr,
     1          arr,lused)
            itype=2
            call hudson_skel_retr0(ier,ison1,itype,skel,nskel,
     1          addr,arr)
            call hudson_skel_store0(ibox,itype,skel,nskel,addr,
     1          arr,lused)
            goto 3000
        endif
c
        if(ison2 .gt. 0) then
            itype=1
            call hudson_skel_retr0(ier,ison2,itype,skel,nskel,
     1          addr,arr)
            call hudson_skel_store0(ibox,itype,skel,nskel,addr,
     1          arr,lused)
            itype=2
            call hudson_skel_retr0(ier,ison2,itype,skel,nskel,
     1          addr,arr)
            call hudson_skel_store0(ibox,itype,skel,nskel,addr,
     1          arr,lused)
            goto 3000
        endif
c
 2000 continue
c
c        the box number ibox has two sonnies, and both have been
c        skeletonized. Merge them things

        if(ibox .eq. 1) goto 4200
c
        call hudson_withkids(ier,ibox,w,xs,addr,arr,
     1      lused,nts,eps,ifrelskel,
     2      skel,ipivots,list,lists,amatr,lamatr,
     3      skelout,rnorms,ipts_inbox,pts_inbox,ncols)
c        
        if(ncols_max .lt. ncols) ncols_max=ncols
 3000 continue
c
        if(ifacted .eq. 0) goto 4200
c        
 4000 continue
 4200 continue
c
        return
        end
c
c
c
c
c
        subroutine hudson_list1_pts(ibox,w,xs,
     1      iptsout,ptsout,nptsout,list)
        implicit real *8 (a-h,o-z)
        save
        integer *4 w(1),box(10),list(1),iptsout(1)
        real *8 ab(2),ptsout(1),xs(1)
c
c        construct list 1 for this box
c
        ii=0

        do 2200 ijk=1,2
c
        itype=1
        if(ijk .eq. 2) itype=11
        call d1m_onelist(ibox,itype,w,list,nlist)
c
        do 2000 i=1,nlist
c
c       extract data for this element of list 1
c
        jbox=list(i)
c
        call d1m_onebox(jbox,w,box,ab)
c
        i1=box(6)
        n1=box(7)
        do 1400 j=1,n1
c
        ii=ii+1
        ptsout(ii)=xs(i1+j-1)
        iptsout(ii)=i1+j-1
 1400 continue
 2000 continue
c
 2200 continue
c
        nptsout=ii
c
        return
        end
c
c
c
c
c
        subroutine hudson_withkids(ier,ibox,w,xs,
     1      addr,arr,lused,nts,eps,ifrelskel,
     2      skel,ipivots,list,lists,amatr,lamatr,
     3      skelout,rnorms,map_xsskel,xsskel,ncols)
        implicit real *8 (a-h,o-z)
        save
        integer *4 w(1),box(10),skel(1),addr(8,1),arr(1),
     1      ipivots(1),map_xsskel(1),list(1),lists(1)
        real *8 ab(2),amatr(1),xs(1),rnorms(1),xsskel(1),skelout(1)
c
c        This subroutine skeletonizes one box having two kids, and 
c        stores the skeletons in arrays addr, arr. PLEASE NOTE THAT 
c        AT THIS TIME (3.21.04), A SINGLE SKELETON IS CREATED, AND 
C        STORED BOTH AS THE INNER AND OUTER ONES. OBVIOUSLY, THIS 
C        WILL HAVE TO BE FIXED.
c
c                    Input parameters:
c
c  ibox - the box to be skeletonized
c  w - the array containing the structure; must have been created via 
c        a preceding call to the subroutine d1mstrcr (see)
c  xs - the used-specified nodes in the simulation
c     addr, arr, lused - parameters used by the subroutine hudson_skel_store0 (see)
c  nts - the number of chebychev nodes in the auxiliary discretization
c        of each box in the list 2 of the box being skeletonized
c  eps - the accuracy of the skeletonization
c
c                    Output parameters:
c
c  ier - error return code
c  addr, arr, lused - parameters used by the subroutine hudson_skel_store0 (see)
c  ncols - the number of elements in (both) inner skeletons created
c
c                    Work arrays:
c
c  skel - must be at least n (the number of points in the simulation)
c        real *8 locations long
c  ipivots - must be at least n (the number of points in the simulation)
c        integer *4 locations long
c  list - must be at least n (the number of points in the simulation)
c        integer *4 locations long
c  lists - must be at least n (the number of points in the simulation)
c        integer *4 locations long
c  amatr - must be at least ????? real *8 locations long
c  skelout - must be at least n (the number of points in the simulation)
c        real *8 locations long
c  rnorms - must be at least n (the number of points in the simulation)
c        real *8 locations long
c  map_xsskel - must be at least nmax integer *4 elements, where
c        nmax is the maximum number of points in a childless box
c  xsskel - must be at least nmax real *8 elements, where
c        nmax is the maximum number of points in a childless box
c
c        . . . if this box fewer than 2 kids - get out
c
        ier=0
c
        call d1m_onebox(ibox,w,box,ab)
c       
        if( (box(4) .le. 0) .or. (box(5) .le. 0) ) then
c
            ier=2
            return
        endif
c
c        retrieve the kiddies' skeletons
c
        ison1=box(4)
        ison2=box(5)
c
        itype=1
        call hudson_skel_retr0(ier,ison1,itype,skel,npts1,addr,arr)
c
        itype=1
        call hudson_skel_retr0(ier,ison2,itype,skel(npts1+1),npts2,
     1      addr,arr)
c
        npts=npts1+npts2
c
        do 1400 i=1,npts
c
        j=skel(i)
        xsskel(i)=xs(j)
        map_xsskel(i)=j
 1400 continue
c
c        construct the outer skeleton
c
        list_type=2
        ipts_type=3
cccc        itype_pts=2
c
        call hudson_pnts_inlists(ibox,w,list_type,ipts_type,
     1      nts,skelout,nskel,xs,amatr,itype_pts,
     2      rnorms,amatr(lamatr/2+1) )
c
        list_type=4
        ipts_type=1
c
        call hudson_pnts_inlists(ibox,w,list_type,ipts_type,
     1      nts,skelout(nskel+1),nskel4,xs,amatr,itype_pts,
     2      rnorms,amatr(lamatr/2+1) )
c
        nskel=nskel+nskel4
c
        ncols=0
        if(nskel .eq. 0) goto 1800
c
c        . . . skeletonize
c
        call hudson_skel(xsskel,npts,skelout,nskel,amatr,
     1      eps,ifrelskel,ipivots,ncols,rnorms)
c
c       find the locations in array xs of the points 
c       in ipivots
c
        call d1m_onebox(ibox,w,box,ab)
        j0=box(6)
        do 1600 i=1,ncols
c
        j=ipivots(i)
        skel(i)=map_xsskel(j)
 1600 continue
c
 1800 continue
c
        itype=1
        call hudson_skel_store0(ibox,itype,skel,ncols,addr,arr,lused)
        itype=2
c
        call hudson_skel_store0(ibox,itype,skel,ncols,addr,arr,lused)
        return
        end
c
c
c
c
c
        subroutine hudson_nokids(jer,ibox,w,xs,nts,eps,ifrelskel,
     1      addr,arr,lused,
     2      ipts_inbox,skel,ipivots,list,lists,pts_inbox,amatr,lamatr,
     3      skelout,rnorms,ncols)
        implicit real *8 (a-h,o-z)
        save
        integer *4 w(1),box(10),ipts_inbox(1),skel(1),
     1      addr(8,1),arr(1),ipivots(1),list(1),lists(1)
        real *8 ab(2),pts_inbox(1),amatr(1),xs(1),
     1      skelout(1),rnorms(1)
c
c        This subroutine skeletonizes one childless box, and stores
c        the skeletons in arraya addr, arr. PLEASE NOTE THAT this 
c        subroutine creates preliminary skeletons; as such, the 
C        INCOMING AND OUTGOING SKELETON COINCIDE. 
c
c                    Input parameters:
c
c  ibox - the box to be skeletonized
c  w - the array containing the structure; must have been created via 
c        a preceding call to the subroutine d1mstrcr (see)
c  xs - the used-specified nodes in the simulation
c  nts - the number of chebychev nodes in the auxiliary discretization
c        of each box in the list 2 of the box being skeletonized
c  eps - the accuracy of the skeletonization
c  addr, arr, lused - parameters used by the subroutine hudson_skel_store0 (see)
c
c                    Output parameters:
c
c  ier - error return code
c  addr, arr, lused - parameters used by the subroutine hudson_skel_store0 (see)
c  ncols - the number of elements in (both) inner skeletons created
c
c                    Work arrays:
c
c  ipts_inbox - must be at least nmax integer *4 elements, where
c        nmax is the maximum number of points in a childless box
c  skel - must be at least n (the number of points in the simulation)
c        real *8 locations long
c  ipivots - must be at least n (the number of points in the simulation)
c        integer *4 locations long
c  list - must be at least n (the number of points in the simulation)
c        integer *4 locations long
c  lists - must be at least n (the number of points in the simulation)
c        integer *4 locations long
c  pts_inbox - must be at least nmax real *8 elements, where
c        nmax is the maximum number of points in a childless box
c  amatr - must be at least ????? real *8 locations long
c  skelout - must be at least n (the number of points in the simulation)
c        real *8 locations long
c  rnorms - must be at least n (the number of points in the simulation)
c        real *8 locations long
c
c        . . . if this box has kids - get out
c
        ier=0
        ncols=0
c
        call d1m_onebox(ibox,w,box,ab)
c       
        if( (box(4) .gt. 0) .or. (box(5) .gt. 0) ) then
c
            ier=2
            return
        endif
c
c        construct the list of all points in this box
c
        call hudson_pnts_inbox3(ibox,w,xs,ipts_inbox,pts_inbox,np)
c
c        construct the outer skeleton
c
        list_type=2
        ipts_type=3
c
        itype_pts=1
        call hudson_pnts_inlists(ibox,w,list_type,ipts_type,
     1      nts,skelout,nskel,xs,amatr,itype_pts,
     2      rnorms,amatr(lamatr/2+1) )
c
        list_type=3
        ipts_type=3
c
        call hudson_pnts_inlists(ibox,w,list_type,ipts_type,
     1      nts,skelout(nskel+1),nskel4,xs,amatr,itype_pts,
     2      rnorms,amatr(lamatr/2+1) )
c
        nskel=nskel+nskel4
        ncols=0
        if(nskel .eq. 0) goto 1800
c
c        . . . skeletonize
c
        call hudson_skel(pts_inbox,np,skelout,nskel,amatr,
     1      eps,ifrelskel,ipivots,ncols,rnorms)
c
c       find the locations in array xs of the points 
c       in ipivots
c
        j0=box(6)
        do 1600 i=1,ncols
c
        j=ipivots(i)
        skel(i)=j0+j-1
c
 1600 continue
c
 1800 continue
c
c        store the obtained skeleton
c
        itype=1
        call hudson_skel_store0(ibox,itype,skel,ncols,addr,arr,lused)
        itype=2
c
        call hudson_skel_store0(ibox,itype,skel,ncols,addr,arr,lused)
c
        return
        end
c
c
c
c
c
        subroutine hudson_allskels2(ier,w,xs,nboxes,nts,eps,
     1      ifrelskel,addr,arr,lused,
     2      skel,ipts_inbox,ipivots,list,lists,
     3      amatr,pts_inbox,skelout,rnorms,ncols_max,matfun,
     4      inds,iskelout)
        implicit real *8 (a-h,o-z)
        save
        integer *4 w(1),box(10),skel(1),ipts_inbox(1),
     1      addr(8,1),arr(1),ipivots(1),list(1),lists(1),
     2      inds(1),iskelout(1)
        real *8 ab(2),amatr(1),pts_inbox(1),
     1      skelout(1),rnorms(1)
c
c        . . . skeletonize all boxes without kids
c
        do 1050 i=1,nboxes
c
        inds(i)=0
 1050 continue
c
        nn=nboxes
c
        ncols_max=0
        do 1400 i=1,nn
c
        ibox=i
c
        itype=1
        call hudson_nokids2(jer,ibox,w,xs,nts,eps,ifrelskel,
     1      addr,arr,lused,
     2      ipts_inbox,skel,ipivots,list,lists,pts_inbox,amatr,
     3      skelout,rnorms,ncols,itype,matfun,iskelout)
c
        if(jer .eq. 2) goto 1400
c
        if(ncols_max .lt. ncols) ncols_max=ncols
c
        inds(i)=1
c
        itype=2
        call hudson_nokids2(jer,ibox,w,xs,nts,eps,ifrelskel,
     1      addr,arr,lused,
     2      ipts_inbox,skel,ipivots,list,lists,pts_inbox,amatr,
     3      skelout,rnorms,ncols,itype,matfun,iskelout)
c
        if(ncols_max .lt. ncols) ncols_max=ncols
 1400 continue
c
c        Scan all boxes. Whenever a box is encountered whose
c        both kids have been skeletonized, skeletonize the
c        said box
c
        do 4000 ijk=1,100
c
        ifacted=0
c
        do 3000 i=1,nboxes
c
        if(inds(i) .gt. 0) goto 3000
c
        ibox=i
c
c        if this box is childless - skip it
c
        call d1m_onebox(ibox,w,box,ab)
c
        ison1=box(4)
        ison2=box(5)
c
        if( (ison1 .le. 0) .and. (ison2 .le. 0) ) goto 3000
c
c        this box has kids. check if both have been skeletonized
c
        ifacted=1
        nsons=0
        ifbadsons=0
c        
        if(ison1 .gt. 0) then
            nsons=nsons+1
            if(inds(ison1) .le. 0) ifbadsons=1
        endif
c
c        
        if(ison2 .gt. 0) then
            nsons=nsons+1
            if(inds(ison2) .le. 0) ifbadsons=1
        endif
c
        if(ifbadsons .eq. 1) goto 3000
c
        inds(i)=1
c
c        This box has at least one kid, and no kids that have not
c        been skeletonized. If this box has only one kid, copy
c        the kid's skeleton into daddy, and be done
c
        if(nsons. ne. 1) goto 2000
c
        if(ison1 .gt. 0) then
            itype=1
            call hudson_skel_retr0(ier,ison1,itype,skel,nskel,
     1          addr,arr)
            call hudson_skel_store0(ibox,itype,skel,nskel,addr,
     1          arr,lused)
            itype=2
            call hudson_skel_retr0(ier,ison1,itype,skel,nskel,
     1          addr,arr)
            call hudson_skel_store0(ibox,itype,skel,nskel,addr,
     1          arr,lused)
            goto 3000
        endif
c
        if(ison2 .gt. 0) then
            itype=1
            call hudson_skel_retr0(ier,ison2,itype,skel,nskel,
     1          addr,arr)
            call hudson_skel_store0(ibox,itype,skel,nskel,addr,
     1          arr,lused)
            itype=2
            call hudson_skel_retr0(ier,ison2,itype,skel,nskel,
     1          addr,arr)
            call hudson_skel_store0(ibox,itype,skel,nskel,addr,
     1          arr,lused)
            goto 3000
        endif
c
 2000 continue
c
c        the box number ibox has two sonnies, and both have been
c        skeletonized. Merge them things

        if(ibox .eq. 1) goto 4200
c
        itype=1
        call hudson_withkids2(ier,ibox,w,xs,addr,arr,
     1      lused,nts,eps,ifrelskel,
     2      skel,ipivots,list,lists,amatr,
     3      skelout,rnorms,ipts_inbox,pts_inbox,ncols,
     4      itype,matfun,iskelout)
c        
c
        itype=2
        call hudson_withkids2(ier,ibox,w,xs,addr,arr,
     1      lused,nts,eps,ifrelskel,
     2      skel,ipivots,list,lists,amatr,
     3      skelout,rnorms,ipts_inbox,pts_inbox,ncols,
     4      itype,matfun,iskelout)
c        
        if(ncols_max .lt. ncols) ncols_max=ncols
 3000 continue
c
        if(ifacted .eq. 0) goto 4200
c        
 4000 continue
 4200 continue
c
c        construct big skeletons for all childless boxes
c
        do 4400 ibox=1,nn
c
        itype=3
c
        call hudson_bigskel(ier,ibox,w,xs,nts,eps,ifrelskel,
     1      addr,arr,lused,
     2      ipts_inbox,skel,ipivots,list,lists,pts_inbox,amatr,
     3      skelout,rnorms,ncols,itype,matfun,iskelout)
c
        if(ncols_max .lt. ncols) ncols_max=ncols
c
        itype=4
c
        call hudson_bigskel(ier,ibox,w,xs,nts,eps,ifrelskel,
cccc        call one_bigskel2(ier,ibox,w,xs,nts,eps/10,
     1      addr,arr,lused,
     2      ipts_inbox,skel,ipivots,list,lists,pts_inbox,amatr,
     3      skelout,rnorms,ncols,itype,matfun,iskelout)
c
        if(ncols_max .lt. ncols) ncols_max=ncols
 4400 continue
c        
        return
        end
c
c
c
c
c
        subroutine hudson_bigskel(ier,ibox,w,xs,nts,eps,ifrelskel,
     1      addr,arr,lused,
     2      ipts_inbox,skel,ipivots,list,lists,pts_inbox,amatr,
     3      skelout,rnorms,ncols,itype,matfun,ipntsout)
        implicit real *8 (a-h,o-z)
        save
        integer *4 w(1),box(10),ipts_inbox(1),skel(1),
     1      addr(8,1),arr(1),ipivots(1),list(1),lists(1),ipntsout(1)
        real *8 ab(2),pts_inbox(1),amatr(1),xs(1),
     1      skelout(1),rnorms(1)
c
c        . . . if this box has kids - get out
c
        ier=0
        ncols=0
c
        call d1m_onebox(ibox,w,box,ab)
c       
        if( (box(4) .gt. 0) .or. (box(5) .gt. 0) ) then
c
            ier=2
            return
        endif
c
c        construct the list of all points in this box
c
        call hudson_pnts_inbox3(ibox,w,xs,ipts_inbox,pts_inbox,np)
c
c        construct the outer skeleton
c
c        . . . take into account the lists 2
c
        list_type=2
        ipts_type=2
        itype_pts=2
        if(itype .eq. 4) itype_pts=1
c
        call hudson_pnts_inlists(ibox,w,list_type,ipts_type,
     1      nts,skelout,nskel,xs,ipntsout,itype_pts,
     2      rnorms,amatr)
c
c        . . . take into account the lists 4
c
        list_type=4
        ipts_type=1
c
        call hudson_pnts_inlists(ibox,w,list_type,ipts_type,
     1      nts,skelout(nskel+1),nskel4,xs,ipntsout(nskel+1),
     2      itype_pts,rnorms,amatr)

c
        nskel=nskel+nskel4


cccc        itype_pts=1



c
c        . . . take into account the lists 3
c
        list_type=3
        ipts_type=2
c
        call hudson_pnts_inlists(ibox,w,list_type,ipts_type,
     1      nts,skelout(nskel+1),nskel3,xs,ipntsout(nskel+1),
     2      itype_pts,rnorms,amatr)

c
        nskel=nskel+nskel3
c
c        . . . take into account the list 1
c
        call hudson_list1_pts(ibox,w,xs,ipntsout(nskel+1),
     1      skelout(nskel+1),nskel1,rnorms)
c
        nskel=nskel+nskel1
        ncols=0
        if(nskel .eq. 0) goto 1800
c
c        . . . skeletonize
c
        itype_skel=itype-2 
c
        call hudson_skel2(ipts_inbox,pts_inbox,np,
     1      ipntsout,skelout,nskel,amatr,eps,ifrelskel,
     2      ipivots,ncols,rnorms,matfun,itype_skel)
c
c       find the locations in array xs of the points 
c       in ipivots
c
        j0=box(6)
        do 1600 i=1,ncols
c
        j=ipivots(i)
        skel(i)=j0+j-1
c
 1600 continue
c
 1800 continue
c
c        store the obtained skeleton
c
        call hudson_skel_store0(ibox,itype,skel,ncols,addr,arr,lused)
c
        return
        end
c
c
c
c
c
        subroutine hudson_skel2(ipts_inbox,pts_inbox,np,
     1      iskelout,skelout,nskel,amatr,eps,ifrelskel,
     2      ipivots,ncols,rnorms,matfun,itype)
        implicit real *8 (a-h,o-z)
        save
        integer *4 ipivots(1),ipts_inbox(1),iskelout(1)
        real *8 pts_inbox(1),skelout(1),amatr(np,nskel),rnorms(1)
c
c        This subroutine skeletonizes a single box. 
c
c
c        . . . construct the matrix of interactions
c
        if(itype .eq. 2) goto 1600
        do 1400 i=1,nskel
        do 1200 j=1,np
c
        x=skelout(i)
        ix=iskelout(i)
c
        y=pts_inbox(j)
        iy=ipts_inbox(j)
        call matfun(y,iy,x,ix,amatr(j,i))
 1200 continue
 1400 continue
c
 1600 continue
c
        if(itype .eq. 1) goto 2600
c
        do 2400 i=1,nskel
        do 2200 j=1,np
c
        x=skelout(i)
        ix=iskelout(i)
        y=pts_inbox(j)
        iy=ipts_inbox(j)
c
        call matfun(x,ix,y,iy,amatr(j,i))
 2200 continue
 2400 continue
c
 2600 continue
c
c        . . . skeletonize
c
        call jacskel(amatr,np,nskel,eps,ifrelskel,
     1      ipivots,ncols,rnorms)
c
c        handle empty skeletons
c
        if(ncols .eq. 0) then
            ncols=1
            ipivots(1)=1
c
c           the following line is included for completeness,
c           even though hudson doesn't use the values that
c           jacskel puts into rnorms
c
            rnorms(1)=0
        endif
c
        return
        end
c
c
c
c
c
        subroutine hudson_withkids2(ier,ibox,w,xs,
     1      addr,arr,lused,nts,eps,ifrelskel,
     2      skel,ipivots,list,lists,amatr,
     3      skelout,rnorms,map_xsskel,xsskel,ncols,itype,matfun,
     4      iskelout)
        implicit real *8 (a-h,o-z)
        save
        integer *4 w(1),box(10),skel(1),addr(8,1),arr(1),
     1      ipivots(1),map_xsskel(1),list(1),lists(1),iskelout(1)
        real *8 ab(2),amatr(1),xs(1),rnorms(1),xsskel(1),skelout(1)
c
c        . . . if this box fewer than 2 kids - get out
c
        ier=0
c
        call d1m_onebox(ibox,w,box,ab)
c       
        if( (box(4) .le. 0) .or. (box(5) .le. 0) ) then
c
            ier=2
            return
        endif
c
c        retrieve the kiddies' skeletons
c
        ison1=box(4)
        ison2=box(5)
c
        call hudson_skel_retr2(ier,w,ison1,itype,xs,
     1      map_xsskel,xsskel,npts1)
c
        call hudson_skel_retr2(ier,w,ison2,itype,xs,
     1      map_xsskel(npts1+1),xsskel(npts1+1),npts2)
c
        npts=npts1+npts2
c
c        construct the outer skeleton
c
        list_type=2
        ipts_type=2
c
        itype_pts=2
        if(itype .eq. 2) itype_pts=1
c
        call hudson_pnts_inlists(ibox,w,list_type,ipts_type,
     1      nts,skelout,nskel,xs,iskelout,itype_pts,
     2      rnorms,amatr)
c
        list_type=4
        ipts_type=1

c
        call hudson_pnts_inlists(ibox,w,list_type,ipts_type,
     1      nts,skelout(nskel+1),nskel4,xs,iskelout(nskel+1),
     2      itype_pts,rnorms,amatr)

c
        nskel=nskel+nskel4
c
        ncols=0
        if(nskel .eq. 0) goto 1800
c
c        . . . skeletonize
c
        call hudson_skel2(map_xsskel,xsskel,npts,
     1      iskelout,skelout,nskel,amatr,eps,ifrelskel,
     2      ipivots,ncols,rnorms,matfun,itype)
c
c       find the locations in array xs of the points 
c       in ipivots
c
        call d1m_onebox(ibox,w,box,ab)
        j0=box(6)
        do 1600 i=1,ncols
c
        j=ipivots(i)
        skel(i)=map_xsskel(j)
 1600 continue
c
 1800 continue
c
c        store the obtained skeleton
c
        call hudson_skel_store0(ibox,itype,skel,ncols,addr,arr,lused)
c
        return
        end
c
c
c
c
c
        subroutine hudson_nokids2(ier,ibox,w,xs,nts,eps,ifrelskel,
     1      addr,arr,lused,
     2      ipts_inbox,skel,ipivots,list,lists,pts_inbox,amatr,
     3      skelout,rnorms,ncols,itype,matfun,iskelout)
        implicit real *8 (a-h,o-z)
        save
        integer *4 w(1),box(10),ipts_inbox(1),skel(1),
     1      addr(8,1),arr(1),ipivots(1),list(1),lists(1),
     2      iskelout(1)
        real *8 ab(2),pts_inbox(1),amatr(1),xs(1),
     1      skelout(1),rnorms(1)
c
c        . . . if this box has kids - get out
c
        ier=0
        ncols=0
c
        call d1m_onebox(ibox,w,box,ab)
c       
        if( (box(4) .gt. 0) .or. (box(5) .gt. 0) ) then
c
            ier=2
            return
        endif
c
c        construct the list of all points in this box
c
        call hudson_pnts_inbox3(ibox,w,xs,ipts_inbox,pts_inbox,np)
c
c        construct the outer skeleton
c
        list_type=2
        ipts_type=2
c
        itype_pts=2
        if(itype .eq. 2) itype_pts=1
        call hudson_pnts_inlists(ibox,w,list_type,ipts_type,
     1      nts,skelout,nskel,xs,iskelout,itype_pts,
     2      rnorms,amatr)
c
        list_type=4
        ipts_type=1
c
        call hudson_pnts_inlists(ibox,w,list_type,ipts_type,
     1      nts,skelout(nskel+1),nskel4,xs,iskelout(nskel+1),
     2      itype_pts,rnorms,amatr)
c
        nskel=nskel+nskel4
        ncols=0
        if(nskel .eq. 0) goto 1800

c
c        . . . skeletonize
c
        call hudson_skel2(ipts_inbox,pts_inbox,np,
     1      iskelout,skelout,nskel,amatr,eps,ifrelskel,
     2      ipivots,ncols,rnorms,matfun,itype)
c
c       find the locations in array xs of the points 
c       in ipivots
c
        j0=box(6)
        do 1600 i=1,ncols
c
        j=ipivots(i)
        skel(i)=j0+j-1
c
 1600 continue
c
 1800 continue

c
c        store the obtained skeleton
c
        call hudson_skel_store0(ibox,itype,skel,ncols,addr,arr,lused)
c
        return
        end
c
c
c
c
c
        subroutine hudson_pnts_inlists(ibox,w,list_type,ipts_type,
     1      nts,skelout,nskel,xs,iskel_out,itype,list,lists)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),ts(100),ab(2),skelout(1),xs(1)
        integer *4 list(1),box(10),lists(1),iskel_out(1)
c
        data ntsold/-1/
c
c        This subroutine returns to the user a list of nodes in
c        a collection of boxes. The boxes to be searched are 
c        determined by the parameters (ibox,list_type), as follows
c
c        list_type=2 means that the subroutine will search lists 2
c                    of the box ibox and all its ancestors
c
c        list_type=4 means that the subroutine will search lists 4
c                    of the box ibox and all its ancestors
c
c        PLEASE NOTE THAT 2 AND 4 ARE THE ONLY PERMITTED VALUES OF
C        THE PARAMETER LIST_TYPE
C
c        The type of nodes to be returned is determined by the 
c        parameter ipts_type, as follows.
c
c        ipts_type=1 means that all nodes in each of the boxes will
c                    be returned
c        ipts_type=2 means that the skeletons will be returned
c
c        ipts_type=3 means that nts chebychev nodes will be constructed
c                    on the box number ibox, and returned to the user
c
c     
c                      Input parameters:
c
c  ibox - the box whose ancestors will be scanned
c  w - array containing the data constructed vis a preceding call to
c        d1mstrcr (see)
c  list_type - determines the type of lists to be searched (see above)
c  ipts_type - determines the type of nodes to be returned (see above)
c  nts - the number of chebychev nodes to be constructed on the box,
c        if ipts_type=3; any other value of ipts_type will cause this
c        parameter to be ignored
c  xs - the nodes in the simulation; not used if ipts_type=3
c  itype - the type of skeletons to be returned:
c
c    itype=1 means that the outgoing skeletons will be returned
c    itype=2 means that the incoming skeletons will be returned
c
c
c                      Output parameters:
c
c  skelout - the actual locations in R^1 of the points returned
c  nskel - the number of elements in skelout
c  iskel_out - the locations in array xs of the nodes in skelout        
c
c
c        . . . if needed, construct the Chebychev nodes
c
        if(ipts_type .ne. 3) goto 1400
        if(nts  .le. 3) goto 1400
        if(ntsold .eq. nts) goto 1400
c
        itypecheb=0
        call chebexps(itypecheb,nts,ts,u,v,whts)
        ntsold=nts
c
 1400 continue
c       
c       construct the list of all boxes in the lists list_type 
c       of the box ibox and its ancestors
c       
        call hudson_ancestors(ibox,w,list_type,
     1      lists,ii,list)
c
c        if the points to be constructed in each of boxes are
c        Chebychev - act accordingly
c
        if(ipts_type .ne. 3) goto 3200
c
c        construct preliminary (Chebychev) skeletons
c        for all boxes in the obtained list
c
        jj=0
        do 3000 i=1,ii
c
        jbox=lists(i)
        call d1m_onebox(jbox,w,box,ab)
c
        u=(ab(2)-ab(1))/2        
        v=(ab(2)+ab(1))/2        
c
        do 2600 j=1,nts
c
        jj=jj+1
        skelout(jj)=u*ts(j)+v
c
 2600 continue
c
 3000 continue
c
        nskel=jj
c
        return
c
 3200 continue
c
c        if the points to be constructed in each of boxes are
c        skeletons - act accordingly
c
        if(ipts_type .ne. 2) goto 4200
c
c        extract the skeletons for all boxes in the obtained list
c
        itype_pts=1
        jj=0
        do 3600 i=1,ii
c
        jbox=lists(i)
c
        call hudson_skel_retr(ier,w,jbox,itype,iskel_out(jj+1),np)
c
        jj=jj+np      
c
 3600 continue
c
        nskel_out=jj
c
        do 3800 i=1,nskel_out
c
        skelout(i)=xs(iskel_out(i))
 3800 continue
c
        nskel=nskel_out
        return
c
 4200 continue
c
c        the points to be constructed are actual points 
c        - act accordingly
c
        jj=0
        do 4600 i=1,ii
c
        jbox=lists(i)
c
        call hudson_pnts_inbox3(jbox,w,xs,iskel_out(jj+1),
     1      skelout(jj+1),np)
c
        jj=jj+np      
c
 4600 continue
c
        nskel=jj
c
        return
c
 5200 continue
c
        return
        end
c
c
c
c
c
        subroutine hudson_ancestors(ibox,w,list_type,
     1      lists,nlist,list)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),ab(2)
        integer *4 list(1),box(10),lists(1)
c
c        This subroutine returns to the user the list of boxes 
c        in lists 2 or 4 of box ibox and all of its ancestors
c        The boxes to be returned are determined by the parameters 
c        (ibox,list_type), as follows
c
c        list_type=2 means that the subroutine will search lists 2
c                    of the box ibox and all its ancestors
c
c        list_type=4 means that the subroutine will search lists 4
c                    of the box ibox and all its ancestors
c
c        PLEASE NOTE THAT 2 AND 4 ARE THE ONLY PERMITTED VALUES OF
C        THE PARAMETER LIST_TYPE
C
c     
c                      Input parameters:
c
c  ibox - the box the lists 2 or 4 of whose ancestors will be scanned
c  w - array containing the data constructed vis a preceding call to
c        d1mstrcr (see)
c  list_type - determines the type of lists to be searched (see above)
c
c
c                      Output parameters:
c
c  lists - the list of boxes in lists 2 or 4 (depending on the parameter
c        list_type) of ibox and all its ancestors
c
c                      Work arrays:
c
c  list - must contain at least as many integer *4 elements as the 
c         length of the longest list 2 or 4
c       
c       
c       . . . conduct the recursion
c       
        icurr=ibox
c
        ii=0
        do 3000 ijk=1,100
c
c        obtain the list 2 of the box icurr
c
        itype=list_type
        call d1m_onelist(icurr,itype,w,list,nlist)
c
        if(nlist .eq. 0) goto 2400
c
        do 2200 i=1,nlist
c
        ii=ii+1
        lists(ii)=list(i)
 2200 continue
 2400 continue
c
        call d1m_onebox(icurr,w,box,ab)
        icurr=box(3)
c
        if(icurr .eq. 1) goto 3200
c
 3000 continue
 3200 continue
c
       nlist=ii
c
        return
        end
c
c
c
c
c
        subroutine hudson_skel(pts_inbox,np,skelout,nskel,amatr,
     1      eps,ifrelskel,ipivots,ncols,rnorms)
        implicit real *8 (a-h,o-z)
        save
        integer *4 ipivots(1)
        real *8 pts_inbox(1),skelout(1),amatr(np,nskel),rnorms(1)
c
c        This subroutine skeletonizes a single box. 
c
c
c        . . . construct the matrix of fake interactions
c
        do 1400 i=1,nskel
        do 1200 j=1,np
c
        amatr(j,i)=1/(skelout(i)-pts_inbox(j))
 1200 continue
 1400 continue
c
c        . . . skeletonize
c
        call jacskel(amatr,np,nskel,eps,ifrelskel,
     1      ipivots,ncols,rnorms)
c
c        handle empty skeletons
c
        if(ncols .eq. 0) then
            ncols=1
            ipivots(1)=1
c
c           the following line is included for completeness,
c           even though hudson doesn't use the values that
c           jacskel puts into rnorms
c
            rnorms(1)=0
        endif
c
        return
        end
c
c
c
c
c
        subroutine hudson_skel_retr2(ier,w,ibox,itype,xs,
     1      ipnts_box,pnts_box,nskelbox)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),xs(1),ipnts_box(1),pnts_box(1),ab(2)
        integer *4 box(10)
c
c        This subroutine returns to the user the skeleton of
c        the user-specified type "itype" for the user-specified 
c        box "ibox". Following are interpretations of the 
c        permitted values of itype:
c
c     itype=1 - outgoing skeleton of the "standard" size (valid
c        outside ibox itself and its neighbors)
c     itype=2 - incoming skeleton of the "standard" size (valid
c        outside ibox itself and its neighbors)
c     itype=3 - outgoing skeleton of the "big" size (valid
c        everywhere outside ibox)
c     itype=4 - incoming skeleton of the "big" size (valid
c        everywhere outside ibox)
c     itype=5 - all of the nodes in the box (arguably, this is
c        not a skeleton at all)
c
c                    Input parameters:
c
c  w - the array containing all of the information about thye 
c        structure; it is fervently hoped that it has been created 
c        during a preceding call to the subroutine hudson_allskels (see)
c  ibox - the sequence number of the box for which the skeleton is
c        to be returned
c  itype - the type of the skeleton to be returned (see above)
c  xs - the nodes in the simulation
c  
c                    Output parameters:
c
c  ier - error return code; 
c      ier = 0 means successful execution
c      ier > 0 means that the requested skeleton has not been found
c  ipnts_box - the list of sequence numbers of the nodes (in array 
c        xs) in the requested skeleton
c  pnts_box - the list of locations in R^1 of the nodes in the 
c        requested skeleton
c  nskelbox - the number of nodes returned in arrays ipnts_box, pnts_box
c
c
        if(itype .ge. 5) goto 1400
c
        call hudson_skel_retr(ier,w,ibox,itype,ipnts_box,nskelbox)
c
        if(ier .ne. 0) return
c
        do 1200 i=1,nskelbox
c
        pnts_box(i)=xs(ipnts_box(i))
 1200 continue
c
        return
c
 1400 continue
c
c       retrieve the skeleton if it is in fact all nodes 
c       on the interval ibox
c
        if(itype .ne. 5) then
            ier=128
            return
        endif
c
        call d1m_onebox(ibox,w,box,ab)
c
        i1=box(6)
        n1=box(7)
        do 1600 i=1,n1
c
        pnts_box(i)=xs(i1+i-1)
        ipnts_box(i)=i1+i-1
 1600 continue
c
        nskelbox=n1
c
        return
        end
c
c
c
c
c
        subroutine hudson_skel_retr(ier,w,ibox,itype,skel,nskel)
        implicit real *8 (a-h,o-z)
        integer *4 skel(1),w(1)
        save
c
c        retrieve from array arr the skeleton of type itype
c        of the box ibox
c        
        iaddr=w(21)
        iarr=w(22)
c
        call hudson_skel_retr0(ier,ibox,itype,skel,nskel,
     1      w(iaddr),w(iarr))
c
        return
c
c
c
c
        entry hudson_skel_size(ier,w,ibox,itype,nskel)
c
c        retrieve from array addr the size of the skeleton of 
c        type itype of the box ibox
c        
        iaddr=w(21)
        call hudson_skel_size0(ier,ibox,itype,nskel,w(iaddr) )
c
        return
        end
c
c
c
c
c
        subroutine hudson_skel_store0(ibox,itype,skel,nskel,
     1      addr,arr,lused)
        implicit real *8 (a-h,o-z)
        integer *4 skel(1),addr(8,1),arr(1)
        save
c
c        store in array arr the skeleton of type itype of 
c        the box ibox; enter the information about this event 
c        in array addr
c
        if(itype .eq. 1) then
            addr(1,ibox)=lused+1
            addr(2,ibox)=nskel
        endif
c
        if(itype .eq. 2) then
            addr(3,ibox)=lused+1
            addr(4,ibox)=nskel
        endif
c
c
        if(itype .eq. 3) then
            addr(5,ibox)=lused+1
            addr(6,ibox)=nskel
        endif
c
        if(itype .eq. 4) then
            addr(7,ibox)=lused+1
            addr(8,ibox)=nskel
        endif
c
        do 1200 i=1,nskel
c
        arr(lused+i)=skel(i)
 1200 continue
c
        lused=lused+nskel
c
        return
c
c
c
c
        entry hudson_skel_size0(ier,ibox,itype,nskel,addr)
c
c        retrieve from array addr the size of the skeleton 
c        of type itypeof the box ibox
c        
        ier=0
c
        if(itype .eq. 1) then
            ii=addr(1,ibox)
            nn=addr(2,ibox)
        endif
c
        if(itype .eq. 2) then
            ii=addr(3,ibox)
            nn=addr(4,ibox)
        endif
c
        if(itype .eq. 3) then
            ii=addr(5,ibox)
            nn=addr(6,ibox)
        endif
c
        if(itype .eq. 4) then
            ii=addr(7,ibox)
            nn=addr(8,ibox)
        endif
c
        if(ii .le. 0) then
            ier=8
            return
        endif
c
        nskel=nn
        return
c
c
c
c
        entry hudson_skel_retr0(ier,ibox,itype,skel,nskel,addr,arr)
c
c        retrieve from array arr the skeleton of type itype
c        of the box ibox
c        
        ier=0
c
        if(itype .eq. 1) then
            ii=addr(1,ibox)
            nn=addr(2,ibox)
        endif
c
        if(itype .eq. 2) then
            ii=addr(3,ibox)
            nn=addr(4,ibox)
        endif
c
        if(itype .eq. 3) then
            ii=addr(5,ibox)
            nn=addr(6,ibox)
        endif
c
        if(itype .eq. 4) then
            ii=addr(7,ibox)
            nn=addr(8,ibox)
        endif
c
        if(ii .le. 0) then
            ier=8
            nskel=0
            return
        endif
c
        do 2200 i=1,nn
c
        skel(i)=arr(ii+i-1)
 2200 continue
c
        nskel=nn
        return
c
c
c
c
        entry hudson_skel_store0init(addr,nboxes)
c
        do 3200 i=1,nboxes
c
        addr(1,i)=-1
        addr(2,i)=-1
        addr(3,i)=-1
        addr(4,i)=-1
c
        addr(5,i)=-1
        addr(6,i)=-1
        addr(7,i)=-1
        addr(8,i)=-1
 3200 continue
c
        return
        end
c
c
c
c
c
        subroutine hudson_pnts_inbox2(ibox,w,pnts,npnts)
        implicit real *8 (a-h,o-z)
        save
        integer *4 pnts(1),w(1),box(10)
        real *8 ab(2)
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
        subroutine hudson_pnts_inbox3(ibox,w,xs,pnts,rpnts,npnts)
        implicit real *8 (a-h,o-z)
        save
        integer *4 pnts(1),w(1),box(10)
        real *8 ab(2),xs(1),rpnts(1)
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
        rpnts(i)=xs(pnts(i))
 1200 continue
c
        return
        end
c
c
c
c
c
        subroutine hudson_stats234(w,xs,nboxes,
     1      nlists2_max,nlists3_max,nlists4_max,npts4_max,
     2      www,skelout,list,lists)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),www(1),xs(1),skelout(1)
        integer *4 list(1),lists(1)
c
c        The data produced by this (admittedly, fairly idiosyncratic)
c        routine are to be used for the purposes of memory allocation.
c
c        For the structure built during a preceding call to the 
c        subroutine d1mstrcr (see), this subroutine returns to 
c        the user three integer parameters of the said structure:
c
c  nlists2_max - the maximum total number of boxes in the lists 2 
c        of the box ibox and all of his ancestors, among all boxes
c        ibox in the structure
c  nlists3_max - the maximum total number of boxes in the lists 3 
c        of the box ibox and all of his ancestors, among all boxes
c        ibox in the structure
c  nlists4_max - the maximum total number of boxes in the lists 4 
c        of the box ibox and all of his ancestors, among all boxes
c        ibox in the structure
c  npts4_max - the maximum total number of nodes in all boxes in the 
c        lists 4 of the box ibox and all of his ancestors, among all 
c        boxes ibox in the structure
c
c                      Input parameters:
c
c  
c
c       . . . for each of the boxes ibox in the structrure, find the 
c             total number of boxes in the lists 2 of the box ibox 
c             and all of its ancestors
c
        nlists2_max=0
        nlists3_max=0
        nlists4_max=0
        npts4_max=0
        do 2000 ibox=2,nboxes
c
        list_type=2
        call hudson_ancestors(ibox,w,list_type,
     1      lists,nlist,list)
        if(nlist .gt. nlists2_max) nlists2_max=nlist
c
        list_type=3
        call hudson_ancestors(ibox,w,list_type,
     1      lists,nlist,list)
        if(nlist .gt. nlists3_max) nlists3_max=nlist
c
        list_type=4
        call hudson_ancestors(ibox,w,list_type,
     1      lists,nlist,list)
        if(nlist .gt. nlists4_max) nlists4_max=nlist
c
c
        list_type=4
        ipts_type=1
c
        itype_pts=1
        call hudson_pnts_inlists(ibox,w,list_type,ipts_type,
     1      nts,skelout,nskel,xs,www,itype_pts,
     2      list,lists)
c
        itype_pts=2
        call hudson_pnts_inlists(ibox,w,list_type,ipts_type,
     1      nts,skelout,nskel,xs,www,itype_pts,
     2      list,lists)
c
        if(nskel .gt. npts4_max) npts4_max=nskel
c
c
 2000 continue


        return
        end
c  
c  
c 
c 
c 
        subroutine hudson_linsto(ier,ibox,itype,jbox,key,
     1      intarr,lenint,reaarr,lenrea,
     2      w,keep)
        implicit real *8 (a-h,o-z)
        save
        integer *4 intarr(1),reaarr(1),w(1)
c
c        This subroutine stores in array w a pair of arrays:
c        the integer array intarr, and the real array reaarr;
c        the arrays are indexed by the integer aparmeters ibox,
c        itype, jbox, key.
c       
c        Actually, the indexing parameters are stored in a linked 
c        list (located in the beginning of the array w), together 
c        with the addresses (in the remainder of array w) of the 
c        corresponding integer (intarr) and real (reaarr) arrays.
c
c        The meanings of the parameter itype are as follows.
c
c     itype=5 denotes a LOLO matrix where both skeletons are "standard"
c
c     itype=17 denotes an EXPAND matrix converting nodes on the 
c        interval "ibox" into the "big" outgoing expansion
c
c     itype=19 denotes a "LOLO" matrix, converting a "big" outgoing
c        expansion into the "standard" outgoing expansion
c
c     itype=6 denotes a TATA matrix where both skeletons are "standard"
c
c     itype=8 denotes a EVAL (technically, TATA) matrix converting 
c        potential on a "regular" incoming skeleton on a childless
c        box into the values of potential at all nodes on the said box
c     itype=18 denotes a TATA matrix converting incoming expansion 
c        on a regular skeleton into incoming expansion on the "big"
c        skeleton    
c     itype=16 denotes a TATA matrix converting incoming expansion 
c        on a "big" skeleton into the potential at all nodes in the
c        box
c     itype=21 denotes a LOTA matrix converting outgoing expansion on
c        a "big" skeleton of ibox into an incoming expansion on a big 
c        skeleton of a box in list 1 of ibox.
c     itype=23 denotes a LOTA matrix converting outgoing expansion on
c        the NODES of a chunk with the same nodes on the same said 
c        chunk (self-interaction of a childless chunk)
c     itype=25 denotes a LOTA-type matrix converting outgoing expansion 
c        on the regular outgoing skeleton of a box on the "big" incoming
c        skeleton of another chunk
c
c     itype=27 denotes a LOTA-type matrix converting outgoing expansion 
c        on the "big" outgoing skeleton of a box on the regular incoming
c        skeleton of another chunk
c
c     itype=33 denotes a "regular" LOTA matrix, converting an outgoing
c        expansion on a "regular" skeleton into an incoming expansion
c        on another "regular" skeleton - encountered when the standard
c        list 2's are processed.
c
c
c        . . . retrieve from array w the miscellaneous data
c
        nstored2=w(1)
        keep2=w(2)
c
        ilistaddr2=w(3)
        iaddr2=w(4)
        iarr2=w(5)
        lenw=w(7)
c
c       if we are about to run out of space in array w - bomb
c
        ier=0
        ninire=2
        nm=iarr2+keep2+lenrea*ninire+lenint
c
        if(nm .gt. lenw) then
            ier=16
            call prinf('bombing from hudson_linsto with ier=*',ier,1)
            stop
        endif
c 
c     store the user-supplied arrays in the array w
c 
        call hudson_linsto0(ibox,itype,jbox,key,
     1      intarr,lenint,reaarr,lenrea,
     2      w(iaddr2),w(ilistaddr2),w(iarr2),nstored2,keep2)
c
        w(1)=nstored2
        w(2)=keep2
c
         keep=keep2+iarr2+2
c
c        if we are about to run out of storage for the array
c        addr, allocate more
c
        laddr=w(6)
        ii=laddr/10
c
        if(ii .gt. nstored2+10) return
        ninire=2
c
cccc        call prinf('before shifting, w(iarr2)=*',w(iarr2),12)
c
c        . . . if this shift will cause us to run out oof space 
c              - bomb
c
        lshift=nstored2*10+110*ninire

        nm=iarr2+keep2+lshift+2
c
        if(nm .gt. lenw) then
            ier=8
            call prinf('bombing from hudson_linsto with ier=*',ier,1)
            stop
        endif
c
        call hudson_lin_shuffle(w(iarr2),keep,nstored2,
     1     ninire,ishift)
c            
        iarr2=iarr2+ishift*ninire 
        w(5)=iarr2
        laddr=laddr+ishift*ninire 
        w(6)=laddr
c
        return
        end
c
c
c
c
c
c
        subroutine hudson_lin_shuffle(arr,keep,nstored,
     1     ninire,ishift)
        implicit real *8 (a-h,o-z)
        save
        dimension arr(1)
c
c        shift the array arr
c
        ishift=nstored*10/ninire+10
        ishift=nstored*10/ninire*2+110
        ishift=nstored*10/ninire+110
c
        do 1200 i=0,keep+1+1
c
        arr(keep+ishift+1-i)=arr(keep-i+1)
 1200 continue
c
        return
        end
c
c  
c 
c 
c 
        subroutine hudson_linret(ier,ibox,itype,jbox,key,w,
     1      ip_intarr,lenint,ip_reaarr,lenrea)
        implicit real *8 (a-h,o-z)
        save
        integer *4 w(1)
c
c        retrieve from array w the various data
c
        nstored2=w(1)
        keep2=w(2)
c
        ilistaddr2=w(3)
        iaddr2=w(4)
        iarr2=w(5)
c 
c       . . . find in the array addr the column corresponding 
c             to the used-specified 4 keys
c
        call hudson_linret0(ier,ibox,itype,jbox,key,
     1      ip_intarr,lenint,ip_reaarr,lenrea,
     2      w(iaddr2),w(ilistaddr2) )

cccc        call prinf('after linret0, lenrea=*',lenrea,1)
c
        ip_reaarr=ip_reaarr+iarr2-1
        ip_intarr=ip_intarr+iarr2-1
c
        return
        end
c 
c 
c 
c 
        subroutine hudson_linini(nboxes,w,lenw)
        implicit real *8 (a-h,o-z)
        save
        integer *4 w(1)
c
c        allocate memory in array w for the arrays listaddr, addr, arr
c
        ilistaddr=21
        llistaddr=nboxes+2
c
        iaddr=ilistaddr+llistaddr
        laddr=200 000    
        laddr=4000
c
        iarr=iaddr+laddr
c
c        initialize the storage-retrieval procedure
c
        nstored=0
        keep=0
c
        w(1)=nstored
        w(2)=keep
        w(3)=ilistaddr
        w(4)=iaddr
        w(5)=iarr
c
        w(6)=laddr
        w(7)=lenw
c
        call hudson_linini0(w(ilistaddr),nboxes)
c
        return
        end
c  
c  
c 
c 
c 
        subroutine hudson_linsto0(ibox,itype,jbox,key,
     1      intarr,lenint,reaarr,lenrea,
     2      addr,listaddr,arr,nstored,keep)
        implicit real *8 (a-h,o-z)
        save
        integer *4 listaddr(1),addr(10,1),intarr(1),arr(1),
     1      reaarr(1)
c 
c       store the user-supplied arrays in the storage area arr,
c       and enter information about this change in arrays listaddr
c 
c        . . . take care of the address array
c
        ninire=2
c
        ilast=listaddr(ibox)
        nstored=nstored+1
c
        i=nstored
c
        addr(1,i)=i
        addr(2,i)=ilast
        addr(3,i)=ibox
        addr(4,i)=jbox
        addr(5,i)=itype
        addr(6,i)=key
        addr(7,i)=keep+1
        addr(8,i)=lenint
        addr(9,i)=keep+1+lenint+2
        addr(10,i)=lenrea
c 
        listaddr(ibox)=i
c
c        store the integer data 
c
        iarr=addr(7,i)
        call hudson_iarrmove(intarr,arr(iarr),lenint)
        keep=keep+lenint
c
c       store the real data 
c
        iarr=addr(9,i)
        call hudson_iarrmove(reaarr,arr(iarr),lenrea*ninire)
        keep=keep+lenrea*ninire+2
c
        return
        end
c  
c 
c 
c 
c 
        subroutine hudson_linret0(ier,ibox,itype,jbox,key,
     1      ip_intarr,lenint,ip_reaarr,lenrea,
     2      addr,listaddr)
        implicit real *8 (a-h,o-z)
        save
        integer *4 listaddr(1),addr(10,1)
c 
c       . . . find in the array addr the column corresponding 
c             to the used-specified 4 keys
c
        ier=0
        ilast=listaddr(ibox)
        if(ilast .gt. 0) goto 2200
        nlist=0
        ier=4
        return
 2200 continue
c 
        ii=0
        do 2800 i=1,10000000
c 
        if( addr(3,ilast) .ne. ibox) goto 2400
        if( addr(4,ilast) .ne. jbox) goto 2400
        if( addr(5,ilast) .ne. itype) goto 2400
c
        ii=ilast
        goto 3400
c
 2400 continue
c
        ilast=addr(2,ilast)
        if(ilast .le. 0) goto 3400
 2800 continue
c
 3400 continue
c
        if(ii .eq. 0) then
            ier=4
            return
        endif
c
c       The requested record has been found. return it
c       to the user
c
        i=ii
c
        ip_intarr=addr(7,i)
        lenint=addr(8,i)
c
        ip_reaarr=addr(9,i)
        lenrea=addr(10,i)

cccc        call prinf('addr(1,i)=*',addr(1,i),10)
c
        return
        end
c 
c 
c 
c 
c 
        subroutine hudson_linini0(listaddr,nboxes)
        implicit real *8 (a-h,o-z)
        save
        integer *4 listaddr(1)
c 
c       . . . initialize the array listaddr
c 
        do 2800 i=1,nboxes
        listaddr(i)=-1
c
 2800 continue
c
        return
        end

c
c
c
c
c
        subroutine hudson_iarrmove(a,b,n)
        integer *4 a(1),b(1)
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
cc        subroutine hudson_skel_retr2(ier,w,ibox,itype,xs,
cc     1      ipnts_box,pnts_box,nskelbox)
c
c        This subroutine returns to the user the skeleton of
c        the user-specified type "itype" for the user-specified 
c        box "ibox". Following are interpretations of the 
c        permitted values of itype:
c
c     itype=1 - outgoing skeleton of the "standard" size (valid
c        outside ibox itself and its neighbors)
c     itype=2 - incoming skeleton of the "standard" size (valid
c        outside ibox itself and its neighbors)
c     itype=3 - outgoing skeleton of the "big" size (valid
c        everywhere outside ibox)
c     itype=4 - incoming skeleton of the "big" size (valid
c        everywhere outside ibox)
c     itype=5 - all of the nodes in the box (arguably, this is
c        not a skeleton at all)
c
cc              . . . . . .
c  
cc         return
cc         end
c  
c 
c 
c 
cc        subroutine hudson_linsto(ier,ibox,itype,jbox,key,
cc     1      intarr,lenint,reaarr,lenrea,
cc     2      w,keep)
c
c        This subroutine stores in array w a pair of arrays:
c        the integer array intarr, and the real array reaarr;
c        the arrays are indexed by the integer aparmeters ibox,
c        itype, jbox, key.
c       
c        Actually, the indexing parameters are stored in a linked 
c        list (located in the beginning of the array w), together 
c        with the addresses (in the remainder of array w) of the 
c        corresponding integer (intarr) and real (reaarr) arrays.
c
c        The meanings of the parameter itype are as follows.
c
c     itype=5 denotes a LOLO matrix where both skeletons are "standard"
c
c     itype=17 denotes an EXPAND matrix converting nodes on the 
c        interval "ibox" into the "big" outgoing expansion
c
c     itype=19 denotes a "LOLO" matrix, converting a "big" outgoing
c        expansion into the "standard" outgoing expansion
c
c     itype=6 denotes a TATA matrix where both skeletons are "standard"
c
c     itype=8 denotes a EVAL (technically, TATA) matrix converting 
c        potential on a "regular" incoming skeleton on a childless
c        box into the values of potential at all nodes on the said box
c     itype=18 denotes a TATA matrix converting incoming expansion 
c        on a regular skeleton into incoming expansion on the "big"
c        skeleton    
c     itype=16 denotes a TATA matrix converting incoming expansion 
c        on a "big" skeleton into the potential at all nodes in the
c        box
c     itype=21 denotes a LOTA matrix converting outgoing expansion on
c        a "big" skeleton of ibox into an incoming expansion on a big 
c        skeleton of a box in list 1 of ibox.
c     itype=23 denotes a LOTA matrix converting outgoing expansion on
c        the NODES of a chunk with the same nodes on the same said 
c        chunk (self-interaction of a childless chunk)
c     itype=25 denotes a LOTA-type matrix converting outgoing expansion 
c        on the regular outgoing skeleton of a box on the "big" incoming
c        skeleton of another chunk
c
c     itype=27 denotes a LOTA-type matrix converting outgoing expansion 
c        on the "big" outgoing skeleton of a box on the regular incoming
c        skeleton of another chunk
c
c     itype=33 denotes a "regular" LOTA matrix, converting an outgoing
c        expansion on a "regular" skeleton into an incoming expansion
c        on another "regular" skeleton - encountered when the standard
c        list 2's are processed.
c
c   
c           . . . . . . .	
c
cc         return
cc         end
c
c
c
c
c
cc        subroutine hudson_store_flds0(ier,ibox,kind,arr,
cc     1      lenarr,index2,w,icurr)
c
c        One call to this subroutine stores in the array w one
c        of field representations, to be later retrieved by the
c        entry hudson_retr_flds0 of this subroutine (see). The type of 
c        the field to be stored is specified by the parameter
c        kind, as follows:
c
c     kind = 1 - regular outgoing expansion 
c     kind = 2 - regular incoming expansion
c     kind = 3 - "big" outgoing expansion; exists for childless 
c        chunks only
c     kind = 4 "big" incoming expansion; exists for childless 
c        chunks only
c    
c            . . . . . .
c
c        return
c        end
c
 
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
