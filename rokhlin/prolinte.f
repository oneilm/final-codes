        implicit real *8 (a-h,o-z)
  
        real *8 w(10 000 000),coefs(10 000),
     1      ts(10 000),whts(10 000)
c 
        real *8 rlams32(10 000)
c 
        real *8 ws(10 000),rlams3(10 000)
  
  
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER c'
        READ *,c
        CALL PRIN2('c=*',c,1 )
C 
        PRINT *, 'ENTER n'
        READ *,n
        CALL PRINf('n=*',n,1 )
  
  
        pointint=0
cccc        pointint=0.9
  
c 
c       discretize the interval [-1,1] into n equispaced nodes
c 
        done=1
        h=2*done/(n-1)
        do 1200 i=1,n
        ts(i)=-1+(i-1)*h
 1200 continue
  
  
  
        lenw=10 000 000
  
        eps=1.0d-14
  
        if(2 .ne. 3) goto 1400
        c2=c*2
        call proquadr(ier,c2,eps,npts,ts,ws,rlams3,nvects,
     1      rlams32,nvects2,w,lenw,lused)
  
  
        n=npts
  
 1400 continue
        ifinit=1
        ifquadr=1
cccc        call prolinte(ier,c,n,ts,coefs,pointint,whts,w,ifinit,keep)
        call prolinte(ier,c,n,ts,coefs,pointint,whts,w,lenw,
     1      ifinit,keep)
  
  
        call prinf('after first prolinte, keep=*',keep,1)
  
        iw=22
        call lotagraph(iw,ts,coefs,n,
     1       'coefficients of the filter*')
  
        ifmasstest=0
  
c 
c       test the obtained interpolation formula on the exponentials
c 
        call coefstest(n,c,ts,coefs,pointint,whts,w,ifmasstest)
  
        call prin2('and whts=*',whts,n)
  
        call prinf('and remember, n=*',n,1)
  
        call prin2('and coefs=*',coefs,n)
  
  
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine coefstest(n,c,ts,coefs,pointint,whts,w,ifmasstest)
        implicit real *8 (a-h,o-z)
        save
        dimension ts(1),coefs(1),freqs(1000),vals(1000),
     1      errs(1000),coeflag1(1000),coeflag2(1000),
     2      vals3(1000),diff3(1000),whts(1000),rints(1000),
     3      errsinte(1000),ptstest(10 000),w(1),errsfreq(10 000)
c 
c       construct the values of the frequency for which the
c       quadrature will be tested
c 
        ntest=100
        itype=0
        call legeexps(itype,ntest,freqs,u,v,wtest)
c 
        do 1200 i=1,ntest
cccc        freqs(i)=(freqs(i)+1)*c/2*0.99
        freqs(i)=(freqs(i)+1)*c/2
 1200 continue
c 
        call prin2('and freqs are*',freqs,ntest)
c 
c       one frequency after another, test the interpolation
c 
        do 2000 i=1,ntest
c 
        d=0
        do 1600 j=1,n
c 
        d=d+cos(freqs(i)*ts(j))*coefs(j)
 1600 continue
c 
        vals(i)=d
        errs(i)=vals(i)-cos(freqs(i)*pointint)
 2000 continue
c 
        call prin2('and vals=*',vals,ntest)
        call prin2('and errs=*',errs,ntest)
c 
c       construct n-point Lagrange interpolation, apply it
c       to each of the cosines, and evaluate the errors
c 
        call lagrini(ts,n,coeflag1)
c 
        x=pointint
        call lagrcoef(ts,coeflag1,n,x,coeflag2)
  
        call prin2('and Lagrange coefficients are*',coeflag2,n)
  
        do 3400 i=1,ntest
c 
        d=0
        do 3200 j=1,n
c 
        dd=cos(ts(j)*freqs(i))
c 
        d=d+coeflag2(j)*dd
 3200 continue
c 
        vals3(i)=d
c 
cccc        diff3(i)=vals3(i)-vals(i)
        diff3(i)=vals3(i)-cos(freqs(i)*pointint)
 3400 continue
  
        call prin2('and vals3=*',vals3,ntest)
        call prin2('and diff3=*',diff3,ntest)
c 
c        plot the logarithms of the two errors
c 
        do 3600 i=1,ntest
c 
        errs(i)=log10(abs(errs(i))+1.0d-15)
        diff3(i)=log10(abs(diff3(i))+1.0d-15)
 3600 continue
  
  
        iw=21
        call lotagraph2(iw,freqs,errs,ntest,freqs,diff3,ntest,
     1       'both errors*')
c 
c       test the weights of the obtained quadrature
c 
        do 4400 i=1,ntest
c 
        d=0
        do 4200 j=1,n
c 
        d=d+whts(j)*cos(freqs(i)*ts(j))
 4200 continue
c 
        rints(i)=d
        errsinte(i)=rints(i)-2/freqs(i)*sin(freqs(i))
 4400 continue
c 
        call prin2('and rints=*',rints,ntest)
        call prin2('and errsinte=*',errsinte,ntest)
c 
        if(ifmasstest .eq. 0) return
c 
c       conduct massive testing of the interpolation
c 
        h=2
        h=h/ntest
        do 4600 i=1,ntest
c 
        ptstest(i)=-1+h*i -h/2
 4600 continue
  
        call prin2('and ptstest are*',ptstest,ntest)
  
        ifinit=0
        ifquadr=0
  
  
        do 5400 ijk=1,ntest
  
        point=ptstest(ijk)
  
        call prolinte(ier,c,n,ts,coefs,point,whts,w,lenw,
     1      ifinit,keep)
  
  
c 
c       one frequency after another, test the interpolation
c 
        do 5000 i=1,ntest
  
  
  
c 
        d=0
        do 4800 j=1,n
c 
        d=d+cos(freqs(i)*ts(j))*coefs(j)
 4800 continue
c 
        vals(i)=d
        errs(i)=vals(i)-cos(freqs(i)*point)
 5000 continue
c 
        call prin2('and vals=*',vals,ntest)
        call prin2('and errs=*',errs,ntest)
        call prin2('and point=*',point,1)
  
  
        errmax=0
        do 5200 k=1,ntest
        if(errmax .lt. errs(k)) errmax=errs(k)
 5200 continue
c 
        errsfreq(ijk)=errmax
  
 5400 continue
  
  
        call prin2('and errsfreq=*',errsfreq,ntest)
  
  
  
        return
        end
c 
c 
c 
c 
c 
        subroutine prolinte(ier,c,n,ts,coefs,pointint,whts,w,lenw,
     1      ifinit,keep)
        implicit real *8 (a-h,o-z)
        save
        real *8 w(1),ts(1),whts(1),coefs(1)
c 
c 
c          This subroutine constructs prolate function-based
c          interpolation formulae with user-supplied nodes.
c          It also produces quadrature formulae (also prolate
c          function - based), with the same user-supplied nodes.
c          Given the user-supplied nodes ts and the point pointint,
c          the subroutine constructs the coefficients coefs of the
c          optimal interpolation formula, interpolating a band-limited
c          function from the nodes ts to the point pointint; the
c          band-limit of the function to be interpolated is c (also
c          user-supplied). The nodes have to be on the interval [-1,1],
c          so that at the band-limit c, the length of the interval
c          on which the interpolation is to be performed will be
c 
c          c*(1-(-1))/2/pi=c/pi
c 
c          wavelengths. PLEASE NOTE THE IMPORTANCE OF THE PARAMETER
c          IFINIT (SEE BELOW).
c 
c                Input parameters:
c 
c  c - the band-limit of the functions to be interpolated
c  n - the number of interpolation nodes in the array ts (also, the
c          number of coefficients returned in the array coefs, and the
c          number of quadrature weights reurned in the array weights)
c  the n nodes on which the interpolation formula will be based
c  pointint - the point at which the function will be interpolated
c  ifinit - the paramater telling the subroutine whether this is the
c          first call with these points ts and bandlimit c.
c       ifinit=1 means that the procedure should be initialized
c       ifinit=0 means that the initialization should be skipped
c 
c      EXPLANATION: the procedure performed by this subroutine consists
c          of two fairly distinct parts. The first part can be viewed as
c          precomputation, and uses as its input the points ts and the
c          bandlimit c; the second part uses as its input certain arrays
c          generated by the first (and safely hidden somewhere in the
c          array w), and the user-supplied point pointint. On the first
c          call with the parameters (ts, c), ifinit should be set to 1;
c          on all subsequent calls, ifinit should be set to 0.
c  w - contains various types of information produced during a preceding
c          call to this subroutine with ifinit=1; please note that this
c          is only an input parameter if ifinit=0; otherwise, this is
c          an output parameter
c  lenw - the amount of space in array w, in real *8 words
c 
c                Output parameters:
c 
c  ier - error return code
c  coefs - the n interpolating coefficients
c  whts - the weights of the quadrature formula
c  w - contains various types of information to be used during subsequent
c          calls to this subroutine; pleasse note that this is only an
c          output parameter if ifinit=1; otherwise, this is an input
c          parameter
c  keep - the number of elements in the array w that shiuld not be
c 
c 
c       . . . allocate space for arrays khi,rlams,rlams2
c 
        imap=10
c 
        if(ifinit .eq. 0) goto 1200
c 
        lmap=20
c 
        ikhi=imap+lmap
        lkhi=c+200
c 
        irlams=ikhi+lkhi
        lrlams=2*c+400
c 
        irlams2=irlams+lrlams
        lrlams2=c+200
c 
        iw=irlams2+lrlams2
  
        w(1)=ikhi+0.1
        w(2)=irlams+0.1
        w(3)=irlams2+0.1
c 
 1200 continue
c 
        if(ifinit .ne. 0) goto 1400
c 
        ikhi=w(1)
        irlams=w(2)
        irlams2=w(3)
c 
 1400 continue
c 
        call prolint1(ier,c,n,ts,coefs,pointint,whts,w(iw),ifinit,
     1      w(ikhi),w(irlams),w(irlams2),w(imap),ltot,lenw)
c 
        if(ifinit .ne. 0) keep=iw+ltot+10
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine prolint1(ier,c,n,ts,coefs,pointint,whts,w,
     1      ifinit,khi,rlams,rlams2,map,ltot,lenw)
        implicit real *8 (a-h,o-z)
        save
        real *8 w(1),ts(1),whts(1),coefs(1)
c 
        real *8 khi(10 000),rlams(20 000),
     1      rlams2(10 000)
c 
        integer *4 map(1)
c 
c 
c       construct prolate functions using the new version
c       of the subroutine
c 
c 
c       if the initialization is not necessary - bypass it
c 
        if(ifinit .eq. 0) goto 3000
  
cccc        lenw=10 000 000
c 
        iw=1
c 
        call prolcrea(ier,c,w(iw),lenw,nvects,nhigh,
     1    khi,rlams,rlams2,keep,lused)
  
        call prin2('after prolcrea khi are*',khi,nvects)
        call prin2('after prolcrea, rlams2=*',rlams2,nvects)
        call prin2('after prolcrea, rlams=*',rlams,nvects*2)
c 
        call prinf('ier=*',ier,1)
        call prinf('nhigh=*',nhigh,1)
        call prinf('keep=*',keep,1)
        call prinf('lused=*',lused,1)
        call prinf('nvects=*',nvects,1)
c 
        if(ier .ne. 0)  return
c 
c        determine the number of prolate functions to use
c 
        eps=1.0d-15
c 
        do 2800 i=1,nvects
        if(rlams2(i) .gt. eps**2) num=i-1
 2800 continue
c 
        call prinf('and num=*',num,1)
        call prin2('and c/num=*',c/num,1)
c 
c       construct the interpolation scheme based on the prolate functions
c 
        iu=iw+keep
        lu=n*num+10
c 
        iv=iu+lu
        lv=n*num+10
c 
        itt=iv+lv
        ltt=n*num+10
c 
        itt=iv+lv
        ltt=n*num+10
c 
        iww=itt+ltt
        lww=n*num+10
c 
        iamatr=iww+lww
        lamatr=n*num+10
c 
        irints=iamatr+lamatr
        lrints=num+2
c 
        iwork=irints+lrints
        lwork=nhigh+2 +1000
c 
        ivals=iwork+lwork
        lvals=num+2
c 
        iscales=ivals+lvals
        lscales=num+2
c 
        irnorms=iscales+lscales
        lrnorms=num+n
c 
        ltot=irnorms+lrnorms
c 
 3000 continue
c 
        if(ifinit .eq. 1) goto 3200
c 
        iw=map(5)
        iu=map(6)
        iv=map(7)
        itt=map(8)
        iww=map(9)
        iamatr=map(10)
        irints=map(11)
        iwork=map(12)
        ivals=map(13)
        iscales=map(14)
        irnorms=map(15)
c 
 3200 continue
c 
        call prolint0(ier,n,num,w(iamatr),w(iw),rlams2,pointint,
     1      ts,coefs,whts,ifinit,w(iu),w(iv),w(itt),w(iww),
     2      w(irints),w(iwork),w(ivals),w(iscales),
     3      w(irnorms) )
c 
c       if the algorithm has been initialized during this call,
c       store the map of the array w in the array map
c 
        if(ifinit .eq. 0) return
c 
        map(5)=iw
        map(6)=iu
        map(7)=iv
        map(8)=itt
        map(9)=iww
        map(10)=iamatr
        map(11)=irints
        map(12)=iwork
        map(13)=ivals
        map(14)=iscales
        map(15)=irnorms
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine prolint0(ier,n,num,amatr,w,rlams2,pointint,
     1     ts,sol,whts,ifinit,u,v,tt,ww,rints,work,
     2      vals,scales,rnorms)
        implicit real *8 (a-h,o-z)
        save
        dimension ts(1),amatr(num,n),w(1),rints(1),
     1      work(1),vals(1),rlams2(1)
c 
c 
        dimension u(1),v(1),ww(1),tt(1),rnorms(1),sol(1),
     1      whts(1),scales(1)
c 
c       if the initialization is not necessary - bypass it
c 
        if(ifinit .eq. 0) goto 2200
c 
        do 1300 i=1,num
c 
        scales(i)=sqrt(rlams2(i))
 1300 continue
c 
c       construct the matrix of values of prolate functions at the
c       equispaced nodes
c 
        do 1600 i=1,n
        do 1400 j=1,num
c 
        call proleva(j-1,ts(i),w,d,der )
c 
        amatr(j,i)=d * scales(j)
 1400 continue
 1600 continue
c 
c       evaluate the integrals of the requisite prolate functions
c 
        do 1800 i=1,num
c 
        call prolunpk(i-1,w,work,nterms)
c 
        rints(i)=work(1)*2 *scales(i)
c 
 1800 continue
c 
 2200 continue
c 
c       evaluate the requisite prolate functions at the
c       test point
c 
        d2=pointint
        do 2400 i=1,num
c 
        call proleva(i-1,d2,w,vals(i),der )
        vals(i)=vals(i)*scales(i)
 2400 continue
c 
c        if the user so requested, construct the decomposition of
c        the appropriate matrices
c 
        eps=1.0d-14
c 
        if(ifinit .ne. 0)
     1      call leastsq(amatr,u,ww,tt,num,n,ncols,rnorms,eps,v)
c 
        call prin2('after leastsq, rnorms=*',rnorms,n)
c 
c       evaluate the coefficients of the interpolation formula
c 
        call leastsq2(u,ww,tt,num,n,ncols,vals,sol,v)
c 
c       evaluate the coefficients of the quadrature formula
c 
        if(ifinit .ne. 0)
     1      call leastsq2(u,ww,tt,num,n,ncols,rints,whts,v)
        return
        end
  
  
  
  
  
  
  
  
  
  
  
  
  
