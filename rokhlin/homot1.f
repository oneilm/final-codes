        implicit real *8 (a-h,o-z)
        dimension w(9 000 000),par1(100),par2(100),
     1      iws(2)
c 
        external fun1,fun3,fun2
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER nord'
        READ *,nord
        CALL PRINF('nord=*',nord,1 )
c 
c 
        PRINT *, 'ENTER nsings'
        READ *,nsings
        CALL PRINF('nsings=*',nsings,1 )
c 
        PRINT *, 'ENTER nshifts'
        READ *,nshifts
        CALL PRINF('nshifts=*',nshifts,1 )
c 
c        initialize the function computer
c 
c 
        eps=1.0d-12
c 
        singmin=-0.65
        singmin=-0.01
cccc        singmin=-0.85
        singmax=1
  
c 
        call fun3ini(singmin,singmax,nsings,nshifts,nord)
c 
        call fun1ini(singmin,singmax,nsings,nshifts)
  
        nfuns=nsings*nshifts*nord
c 
        nsings1=nsings*nshifts
        npols=nord
  
  
c 
  
        ax=0
        bx=70
c 
        k=20
        lab=20 000
        lw=9 000 000
  
        call prinf('originally,  lw=*',lw,1)
c 
        iws(1)=45
        iws(2)=0
  
        call prinf('before homot0, lw=*',lw,1)
  
        time1=clotatim()
  
cccc        call homot0(ier,ax,bx,k,fun3,par1,par2,par3,
cccc     1      nfuns,eps,w,lw,lused,iws,iterstot)
  
  
        call homot1(ier,ax,bx,k,fun1,fun2,par1,par2,
     1      nsings1,npols,eps,w,lw,lused,iws,iterstot)
  
        call prinf('and after homot0, ier=*',ier,1)
        call prinf('and after homot0, lused=*',lused,1)
        call prinf('and after homot0, iterstot=*',iterstot,1)
  
        time2=clotatim()
  
        call prin2('and time for homot0 is*',time2-time1,1)
        stop
        end
c 
c 
c 
c 
c 
        subroutine fun1(x,i,par1,par2,f)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension sings(10 000),shifts(10 000)
        integer *4 ishifts(200 000),isings(200 000),
     1      whtsings(100 000)
c 
c        decode the function number i to obtain the number
c        of the singularity, the number of the power, and the
c        number of the shift
c 
cccc        call prin2('in fun1, x=*',x,1)
cccc        call prinf('in fun1, i=*',i,1)
  
         t=exp(-x)
  
        ising=isings(i)
        ishift=ishifts(i)
c 
        f=(t+shifts(ishift))**sings(ising)
  
        if(ising .eq. nsings) f=log(t+shifts(ishift))
  
  
cccc        call prin2('in fun1, f=*',f,1)
  
  
c 
        f=f*t
c 
        return
c 
c 
c 
c 
        entry fun1ini(singmin,singmax,nsings7,nshifts)
c 
c        create the array of singularities
c 
        nsings=nsings7
        h=(singmax-singmin)/(nsings-1)
c 
        itype=1
        call legeexps(itype,nsings,sings,u,v,whtsings)
  
        alpha=(singmax-singmin)/2
        beta=(singmax+singmin)/2
        do 2200 j=1,nsings
        sings(j)=alpha*sings(j)+beta
 2200 continue
c 
c       create the array of shifts
c 
        done=1
        h=done/nshifts
        do 2400 j=1,nshifts
        shifts(j)=(j-1)*h
  
cccc        shifts(j)=exp(-20*shifts(j))
        shifts(j)=exp(-40*shifts(j))
 2400 continue
c 
        do 2500 j=1,nshifts
c 
        shifts(j)=shifts(j)-shifts(nshifts)
 2500 continue
c 
        call prin2('in fun1ini, sings=*',sings,nsings)
        call prin2('in fun1ini, shifts=*',shifts,nshifts)
c 
c       create the index tables for the singularities, powers, and
c       shift
c 
        iijjkk=0
c 
        do 2800 ii=1,nsings
        do 2600 jj=1,nshifts
c 
        iijjkk=iijjkk+1
c 
        isings(iijjkk)=ii
        ishifts(iijjkk)=jj
 2600 continue
 2800 continue
  
cccc        call prinf('isings as created*',isings,iijjkk)
cccc        call prinf('ishifts as created*',ishifts,iijjkk)
c 
  
cccc        if(2 .ne. 3) stop
        return
        end
c 
c 
c 
c 
c 
        subroutine fun2(x,i,f)
        implicit real *8 (a-h,o-z)
c 
        save
         t=exp(-x)
  
cccc        call prinf('in fun2, i=*',i,1)
cccc        call prin2('in fun2, t=*',t,1)
  
        d=t*2-1
        call legepol(d,i-1,f,der)
c 
  
cccc        call prin2('in fun2, f=*',f,1)
  
cccc        f=t**i
  
        return
        end
c 
c 
c 
c 
c 
        subroutine fun3(x,i,par1,par2,f)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension sings(10 000),shifts(10 000)
        integer *4 ishifts(2000 000),isings(2000 000),
     1      iords(2000 000),whtsings(100 000)
c 
c        decode the function number i to obtain the number
c        of the singularity, the number of the power, and the
c        number of the shift
c 
cccc        call prin2('in fun3, x=*',x,1)
cccc        call prinf('in fun3, i=*',i,1)
  
         t=exp(-x)
  
        ising=isings(i)
        ishift=ishifts(i)
        iord=iords(i)
c 
        tt=2*t-1
        call legepol(tt,iord,pol,der)
  
cccc        call prin2('in fun3, pol=*',pol,1)
  
  
c 
        f=(t+shifts(ishift))**sings(ising)
  
        if(ising .eq. nsings) f=log(t+shifts(ishift))
  
  
cccc        call prin2('in fun3, f=*',f,1)
  
  
c 
        f=f*pol *t
  
cccc        f=pol
c 
        return
c 
c 
c 
c 
        entry fun3ini(singmin,singmax,nsings7,nshifts,nord)
c 
c        create the array of singularities
c 
        nsings=nsings7
        h=(singmax-singmin)/(nsings-1)
c 
        itype=1
        call legeexps(itype,nsings,sings,u,v,whtsings)
  
        alpha=(singmax-singmin)/2
        beta=(singmax+singmin)/2
        do 2200 j=1,nsings
        sings(j)=alpha*sings(j)+beta
 2200 continue
c 
c       create the array of shifts
c 
        done=1
        h=done/nshifts
        do 2400 j=1,nshifts
        shifts(j)=(j-1)*h
  
cccc        shifts(j)=exp(-20*shifts(j))
        shifts(j)=exp(-40*shifts(j))
 2400 continue
c 
        do 2500 j=1,nshifts
c 
        shifts(j)=shifts(j)-shifts(nshifts)
 2500 continue
c 
        call prin2('in fun3ini, sings=*',sings,nsings)
        call prin2('in fun3ini, shifts=*',shifts,nshifts)
c 
c       create the index tables for the singularities, powers, and
c       shift
c 
        iijjkk=0
c 
        do 3000 kk=1,nord
        do 2800 ii=1,nsings
        do 2600 jj=1,nshifts
c 
        iijjkk=iijjkk+1
c 
        iords(iijjkk)=kk-1
        isings(iijjkk)=ii
        ishifts(iijjkk)=jj
 2600 continue
 2800 continue
 3000 continue
  
cccc        call prinf('iords as created*',iords,iijjkk)
cccc        call prinf('isings as created*',isings,iijjkk)
cccc        call prinf('ishifts as created*',ishifts,iijjkk)
c 
  
cccc        if(2 .ne. 3) stop
        return
        end
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c       This is the end of the debugging code and the beginning of the
c       code for the design of quadratures
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
  
        subroutine homot1(ier,a,b,k,fun1,fun2,par1,par2,
     1      nsings,npols,eps,www,lw,lused2,iws,iterstot)
        implicit real *8 (a-h,o-z)
        save
        dimension par1(1),par2(1),ab(2,100 000),
     1      rlams(1000),rints(1000),xs0(1000),xsout(1000),
     2      wsout(1000),xs(1000),ws(1000),zeroes(1000),www(1),
     3      iws(2)
c 
c       ATTENTION: This subroutine is a somewhat specialized
c       version of the subroutine homot0; in many cases when
c       one would thing of using this subroutine, homot0 should
c       used instead. Basically, the difference between homot1 and
c       homot0 is that homot0 uses the subroutine allsvbld to
c       construct the singular value decomposition, and homot1
c       uses the more specialized allsvcmb.
c 
c       This subroutine designs Generalized Gaussian quadratures
c       for a user-specified set of basis functions. The user
c       supplies the functions to be integrated via the subroutines
c       fun1, fun2 (see below). The exact conditions under which this
c       subroutine is guaranteed to work are complicated and given
c       by Markov's theorems (totally positive kernels, Chebychev
c       systems, and other nonsense). On the other hand, this
c       subroutine often works when there appear to be no known
c       theorems to guarantee its successful operation. In other
c       terms, THIS IS NOT A BLACK BOX; HUMAN INTERVENTION IS
c       OFTEN NEEDED.
c                          Input parameters:
c 
c  a - the left end of the interval on which the quadratures are
c       to be designed
c  b - the right end of the interval on which the quadratures are
c       to be designed
c  k - the number of Gaussian nodes on each of nn subintervals in the
c       nested subdivision of the interval [a,b].
c  fun3 - the user-supplied function to be integrated
c      The calling sequence of fun must be:
c 
c      fun3(x,i,par1,par2,par3,f),
c 
c       where x is the input argument where the function is to be
c       evaluated, i the sequence number of function, f (output)
c       is the value of the i-th function at the point x, and par
c       is an input parameter, carrying whatever information fun
c       might need.
c 
c  par1,par2,par3 - whatever the user-supplied subroutine fun3 needs
c  nfuns - the number of functions to be compres
c  eps - the accuracy to which the legendre expansion must represent
c      the user-supplied function for the subroutine to decide that
c      the representation is accurate
c  lw - the length of the user-provided array w. It has to be at
c      least 17/8*n*nfuns+n+2*k**2+1000. It it is insufficient,
c      ier is set to 512, and the execution of the subroutine
c      terminated
c  iw - an integer 2-element array specifying the FORTRAN file numbers
c      on which the obtained nodes, weights and (possibly) the
c      singular functions will be stored. EXPLANATION: The file
c      iws(1) will contain the quadrature nodes and weights, as well as
c      a small amount of relavant information. The file iws(2) will
c      contain the same information, and, in addition, the nested
c      Legendre expansions of all the obtained singular functions.
c      All these data are expected to be retrieved by means of the
c      subroutines xssretr, xsscoeret (see). The information on the
c      units iws(1), iws(2) is also readable by humans (some of the
c      humans, anyway)
c 
c                    Output parameters:
c 
c  ier - error return code.
c            ier = 0 means that at least one quadrature believed to
c                   be valid has been created
c            ier = 16 means that no quadratures  believed to
c                   be valid have been created
c            ier = 1024 means that the subroutine newton subroutine
c                   has run out of memory; an unlikely and perhaps
c                   impossible event
c            ier = 2048 means that the subroutine nodes0 preparing the
c                   starting point for the newton subroutine
c                   has run out of memory.
c            ier = 4096 means that the subroutine allsvbld at the
c                   heart of this code has run out of memory
c  lused2 - the maximum amount of space in the array www used by the
c      subroutine at any time
c  iterstot - the total number of Newton iterations performed by
c      the subroutine
c 
c                    Work arrays:
c 
c  w - must be LOOOONG. Also, see the description of the parameter
c       lw above.
c 
        external fun1,fun2
c 
c        construct the basis functions
c 
        ix=1
        lx=100 000
c 
        iwhts=ix+lx
        lwhts=100 000
c 
        iwallq=iwhts+lwhts
c 
        ier=0
        lab=100 000
        ifab=1
  
        call prinf('before allsvbld, lw=*',lw,1)
  
        lrlams=1000
        igraph=0
c 
        call allsvcmb(ier,a,b,k,fun1,par1,par2,fun2,
     1      nsings,npols,eps,ab,lab,ifab,nn,www(ix),www(iwhts),
     2      n,rlams,lrlams,
     3      www(iwallq),lw,ncols,rints,lused,igraph)
c 
        call prinf('in homot1 after allsvcmb, jer=*',jer,1)
        call prinf('in homot1 after allsvcmb, nn=*',nn,1)
        call prinf('in homot1 after allsvcmb, k=*',k,1)
        call prinf('in homot1 after allsvcmb, ncols=*',ncols,1)
c 
        if(jer .eq. 0) goto 1200
        ier=4098
        return
c 
 1200 continue
c 
         if(iws(1) .gt. 0) call xsstorei(iws(1),ncols,a,b,rlams)
         if(iws(2) .gt. 0) call xsstorei(iws(2),ncols,a,b,rlams)
c 
c        if the user so requested - store on disk all of the
c        obtained basis functions, together with the related
c        information
c 
        if(iws(2) .gt. 0) call xssclose
     1    (iws(2),nn,ab,www(iwallq),k,ncols,www(ix),www(iwhts) )
c 
c        perform garbage collection
c 
        do 1300 i=1,lused
        www(i)=www(iwallq+i-1)
 1300 continue
c 
        lused2=lused
        lused=n*ncols+10
c 
        iwallq=1
c 
        lwleft=lw-lused
        call prinf('after allsvcmb, lw=*',lw,1)
        call prinf('after allsvcmb, lused=*',lused,1)
        call prinf('after allsvcmb, lwleft=*',lwleft,1)
c 
c        write on disk the header for the file to contain
c        nodes and weights of the quadratures obtained by
c        this subroutine
c 
        do 1400 i=1,1000
        zeroes(i)=0
 1400 continue
c 
c        for all possible numbers of nodes starting with one,
c        find the Gaussian quadratures and store them on disk
c 
        nnnmax=ncols/2
        iterstot=0
        do 3000 ijk=1,nnnmax
c 
        call prinf('in homot1, ijk=*',ijk,1)
c 
        nfuse=ijk*2
c 
c       allocate memory for the determination of the initial
c       approximation to the nodes
c 
        ninire=2
        nmult=8
        numpts=nmult*nn*k
c 
        iamatr=lused+10
        lamatr=numpts*nfuse+10
c 
        ibmatr=iamatr+lamatr
        lbmatr=numpts*nfuse+10
c 
        ix2=ibmatr+lbmatr
        lx2=numpts+10
c 
        iwhts2=ix2+lx2
        lwhts2=numpts+10
c 
        irnorms=iwhts2+lwhts2
        lrnorms=numpts*2+10
c 
        iipivots=irnorms+lrnorms
        lipivots=(numpts+10)/ninire
c 
        iiwhich=iipivots+lipivots
        liwhich=(numpts+10)/ninire
c 
        ltot=iiwhich+liwhich
c 
        call prinf('ltot before nodes0 is*',ltot,1)
c 
        if(lused2 .lt. ltot) lused2=ltot
        if(ltot .lt. lw) goto 1600
        ier=2048
        return
 1600 continue
c 
c       find initial approximations to Gaussian nodes
c 
        eps3=eps*10
c 
        call nodes0(nfuse,www(iwallq),ab,nn,k,www(ix2),www(iwhts2),
     1      www(iamatr),www(ibmatr),eps3,xs0,nmult,
     2      www(irnorms),www(iipivots),www(iiwhich) )
c 
c        attempt to construct the Gaussian nodes and weights
c 
        n=ijk
c 
        dirmass0=1.0d7
        numit=15
        eps1=eps*100
c 
        nnn=nfuse/2
c 
        call prin2('before newtons, xs0=*',xs0,nnn)
c 
c       allocate memory for the Newton refinement of the obtained
c       approximate quadrature
c 
        iamatr=lused+10
        lamatr=n*n*4+10
c 
        ibmatr=iamatr+lamatr
        lbmatr=n*n+10
c 
        ibmatr2=ibmatr+lbmatr
        lbmatr2=n*n*2+10
c 
        iuu=ibmatr2+lbmatr2
        luu=n**2*4+10
c 
        ivv=iuu+luu
        lvv=n**2*4+10
c 
        itt=ivv+lvv
        ltt=n**2*4+10
c 
        iww=itt+ltt
        lww=n**2*4+10
c 
        irnorms=iww+lww
        lrnorms=n*2+10
c 
        ltot=irnorms+lrnorms
c 
        call prinf('iamatr before newtons is*',iamatr,1)
        call prinf('ltot before newtons is*',ltot,1)
c 
        if(ltot .gt. lused2) lused2=ltot
        if(ltot .lt. lw) goto 1800
        ier=1024
        return
 1800 continue
c 
        call newtons(jer,nnn,www(iwallq),ab,nn,k,www(iamatr),
     1      xs0,numit,www(ibmatr),www(ibmatr2),xs,ws,eps1,xsout,wsout,
     2      dirmass0,rints,nitertot,
     3      www(iuu),www(ivv),www(itt),www(irnorms),www(iww) )
c 
        call prinf('after newtons, jer=*',jer,1)
        call prinf('after newtons, nitertot=*',nitertot,1)
c 
        iterstot=iterstot+nitertot
c 
 2200 format('c   Please note that for the number ',
     1    'of nodes equal to ',
     2    i4,/,'c   the non-linear search has failed. Both the ',
     3    'nodes and',/,'c   the weights are garbage')
c 
 2400 format('c            ')
c 
        if(iws(1) .le. 0) goto 2600
c 
        if(jer .ne. 0) write(iws(1),2400)
        if(jer .ne. 0) write(iws(1),2200) nnn
        if(jer .ne. 0) write(iws(1),2400)
c 
        if(jer .eq. 0) call xsstore(iws(1),xsout,wsout,nnn)
        if(jer .ne. 0) call xsstore(iws(1),zeroes,zeroes,nnn)
c 
 2600 continue
c 
        if(iws(2) .le. 0) goto 2800
c 
        if(jer .ne. 0) write(iws(2),2400)
        if(jer .ne. 0) write(iws(2),2200) nnn
        if(jer .ne. 0) write(iws(2),2400)
c 
        if(jer .eq. 0) call xsstore(iws(2),xsout,wsout,nnn)
        if(jer .ne. 0) call xsstore(iws(2),zeroes,zeroes,nnn)
c 
 2800 continue
c 
        call prin2('after newtons, xsout=*',xsout,nnn)
        call prin2('after newtons, wsout=*',wsout,nnn)
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
  
        subroutine homot0(ier,a,b,k,fun3,par1,par2,par3,
     1      nfuns,eps,www,lw,lused2,iws,iterstot)
        implicit real *8 (a-h,o-z)
        save
        dimension par1(1),par2(1),par3(1),ab(2,100 000),
     1      rlams(1000),rints(1000),xs0(1000),xsout(1000),
     2      wsout(1000),xs(1000),ws(1000),zeroes(1000),www(1),
     3      iws(2)
c 
c       This subroutine designs Generalized Gaussian quadratures
c       for a user-specified set of basis functions. The user
c       supplies the functions to be integrated via the subroutine
c       fun3 (see below). The exact conditions under which this
c       subroutine is guaranteed to work are complicated and given
c       by Markov's theorems (totally positive kernels, Chebychev
c       systems, and other nonsense). On the other hand, this
c       subroutine often works when there appear to be no known
c       theorems to guarantee its successful operation. In other
c       terms, THIS IS NOT A BLACK BOX; HUMAN INTERVENTION IS
c       OFTEN NEEDED.
c                          Input parameters:
c 
c  a - the left end of the interval on which the quadratures are
c       to be designed
c  b - the right end of the interval on which the quadratures are
c       to be designed
c  k - the number of Gaussian nodes on each of nn subintervals in the
c       nested subdivision of the interval [a,b].
c  fun3 - the user-supplied function to be integrated
c      The calling sequence of fun must be:
c 
c      fun3(x,i,par1,par2,f),
c 
c       where x is the input argument where the function is to be
c       evaluated, i the sequence number of function, f (output)
c       is the value of the i-th function at the point x, and par1,
c       par2 are input parameters, carrying whatever information fun3
c       might need.
c 
c  par1,par2 - whatever the user-supplied subroutine fun3 needs
c  nfuns - the number of functions to be compresed and integrated
c  eps - the accuracy to which the legendre calculations will
c       be performed
c  lw - the length of the user-provided array w. It should be LARGE;
c      if it is insufficient, the error return code ier (see below)
c      is set appropriately, and the execution of the subroutine is
c      terminated.
c  iws - an integer 2-element array specifying the FORTRAN file numbers
c      on which the obtained nodes, weights and (possibly) the
c      singular functions will be stored. EXPLANATION: The file
c      iws(1) will contain the quadrature nodes and weights, as well as
c      a small amount of relavant information. The file iws(2) will
c      contain the same information, and, in addition, the nested
c      Legendre expansions of all the obtained singular functions.
c      All these data are expected to be retrieved by means of the
c      subroutines xssretr, xsscoeret (see). The information on the
c      units iws(1), iws(2) is also readable by humans (some of the
c      humans, anyway)
c 
c                    Output parameters:
c 
c  ier - error return code.
c            ier = 0 means that at least one quadrature believed to
c                   be valid has been created
c            ier = 16 means that no quadratures  believed to
c                   be valid have been created
c            ier = 1024 means that the subroutine newton subroutine
c                   has run out of memory; an unlikely and perhaps
c                   impossible event
c            ier = 2048 means that the subroutine nodes0 preparing the
c                   starting point for the newton subroutine
c                   has run out of memory.
c            ier = 4096 means that the subroutine allsvbld at the
c                   heart of this code has run out of memory
c  lused2 - the maximum amount of space in the array www used by the
c      subroutine at any time
c  iterstot - the total number of Newton iterations performed by
c      the subroutine
c 
c                    Work arrays:
c 
c  w - must be LOOOONG. Also, see the description of the parameter
c       lw above.
c 
        external fun3
c 
c        construct the basis functions
c 
        ix=1
        lx=100 000
c 
        iwhts=ix+lx
        lwhts=100 000
c 
        iwallq=iwhts+lwhts
c 
        ier=0
        lab=100 000
        ifab=1
  
        call prinf('before allsvbld, lw=*',lw,1)
  
        lrlams=1000
        igraph=0
c 
        call allsvbld(jer,a,b,k,fun3,par1,par2,nfuns,eps,
     1    ab,lab,ifab,nn,www(ix),www(iwhts),n,rlams,lrlams,
     2      www(iwallq),lw,ncols,rints,lused,igraph)
c 
        call prinf('in homot0 after allsvbld, jer=*',jer,1)
        call prinf('in homot0 after allsvbld, nn=*',nn,1)
        call prinf('in homot0 after allsvbld, k=*',k,1)
        call prinf('in homot0 after allsvbld, ncols=*',ncols,1)
c 
        if(jer .eq. 0) goto 1200
        ier=4098
        return
c 
 1200 continue
c 
         if(iws(1) .gt. 0) call xsstorei(iws(1),ncols,a,b,rlams)
         if(iws(2) .gt. 0) call xsstorei(iws(2),ncols,a,b,rlams)
c 
c        if the user so requested - store on disk all of the
c        obtained basis functions, together with the related
c        information
c 
        if(iws(2) .gt. 0) call xssclose
     1    (iws(2),nn,ab,www(iwallq),k,ncols,www(ix),www(iwhts) )
c 
c        perform garbage collection
c 
        do 1300 i=1,lused
        www(i)=www(iwallq+i-1)
 1300 continue
c 
        lused2=lused
        lused=n*ncols+10
c 
        iwallq=1
c 
        lwleft=lw-lused
        call prinf('after allsvbld, lw=*',lw,1)
        call prinf('after allsvbld, lused=*',lused,1)
        call prinf('after allsvbld, lwleft=*',lwleft,1)
c 
c        write on disk the header for the file to contain
c        nodes and weights of the quadratures obtained by
c        this subroutine
c 
        do 1400 i=1,1000
        zeroes(i)=0
 1400 continue
c 
c        for all possible numbers of nodes starting with one,
c        find the Gaussian quadratures and store them on disk
c 
        nnnmax=ncols/2
        iterstot=0
        do 3000 ijk=1,nnnmax
c 
        call prinf('in homot0, ijk=*',ijk,1)
c 
        nfuse=ijk*2
c 
c       allocate memory for the determination of the initial
c       approximation to the nodes
c 
        ninire=2
        nmult=8
        numpts=nmult*nn*k
c 
        iamatr=lused+10
        lamatr=numpts*nfuse+10
c 
        ibmatr=iamatr+lamatr
        lbmatr=numpts*nfuse+10
c 
        ix2=ibmatr+lbmatr
        lx2=numpts+10
c 
        iwhts2=ix2+lx2
        lwhts2=numpts+10
c 
        irnorms=iwhts2+lwhts2
        lrnorms=numpts*2+10
c 
        iipivots=irnorms+lrnorms
        lipivots=(numpts+10)/ninire
c 
        iiwhich=iipivots+lipivots
        liwhich=(numpts+10)/ninire
c 
        ltot=iiwhich+liwhich
c 
        call prinf('ltot before nodes0 is*',ltot,1)
c 
        if(lused2 .lt. ltot) lused2=ltot
        if(ltot .lt. lw) goto 1600
        ier=2048
        return
 1600 continue
c 
c       find initial approximations to Gaussian nodes
c 
        eps3=eps*10
c 
        call nodes0(nfuse,www(iwallq),ab,nn,k,www(ix2),www(iwhts2),
     1      www(iamatr),www(ibmatr),eps3,xs0,nmult,
     2      www(irnorms),www(iipivots),www(iiwhich) )
c 
c        attempt to construct the Gaussian nodes and weights
c 
        n=ijk
c 
        dirmass0=1.0d7
        numit=15
        eps1=eps*100
c 
        nnn=nfuse/2
c 
        call prin2('before newtons, xs0=*',xs0,nnn)
c 
c       allocate memory for the Newton refinement of the obtained
c       approximate quadrature
c 
        iamatr=lused+10
        lamatr=n*n*4+10
c 
        ibmatr=iamatr+lamatr
        lbmatr=n*n+10
c 
        ibmatr2=ibmatr+lbmatr
        lbmatr2=n*n*2+10
c 
        iuu=ibmatr2+lbmatr2
        luu=n**2*4+10
c 
        ivv=iuu+luu
        lvv=n**2*4+10
c 
        itt=ivv+lvv
        ltt=n**2*4+10
c 
        iww=itt+ltt
        lww=n**2*4+10
c 
        irnorms=iww+lww
        lrnorms=n*2+10
c 
        ltot=irnorms+lrnorms
c 
        call prinf('iamatr before newtons is*',iamatr,1)
        call prinf('ltot before newtons is*',ltot,1)
c 
        if(ltot .gt. lused2) lused2=ltot
        if(ltot .lt. lw) goto 1800
        ier=1024
        return
 1800 continue
c 
        call newtons(jer,nnn,www(iwallq),ab,nn,k,www(iamatr),
     1      xs0,numit,www(ibmatr),www(ibmatr2),xs,ws,eps1,xsout,wsout,
     2      dirmass0,rints,nitertot,
     3      www(iuu),www(ivv),www(itt),www(irnorms),www(iww) )
c 
        call prinf('after newtons, jer=*',jer,1)
        call prinf('after newtons, nitertot=*',nitertot,1)
c 
        iterstot=iterstot+nitertot
c 
 2200 format('c   Please note that for the number ',
     1    'of nodes equal to ',
     2    i4,/,'c   the non-linear search has failed. Both the ',
     3    'nodes and',/,'c   the weights are garbage')
c 
 2400 format('c            ')
c 
        if(iws(1) .le. 0) goto 2600
c 
        if(jer .ne. 0) write(iws(1),2400)
        if(jer .ne. 0) write(iws(1),2200) nnn
        if(jer .ne. 0) write(iws(1),2400)
  
  
        caccur=rlams(nnn*2)/rlams(1)
c 
        if(jer .eq. 0) call xsstore(iws(1),xsout,wsout,nnn)
        if(jer .ne. 0) call xsstore(iws(1),zeroes,zeroes,nnn)
c 
 2600 continue
c 
        if(iws(2) .le. 0) goto 2800
c 
        if(jer .ne. 0) write(iws(2),2400)
        if(jer .ne. 0) write(iws(2),2200) nnn
        if(jer .ne. 0) write(iws(2),2400)
c 
        if(jer .eq. 0) call xsstore(iws(2),xsout,wsout,nnn)
        if(jer .ne. 0) call xsstore(iws(2),zeroes,zeroes,nnn)
c 
 2800 continue
c 
        call prin2('after newtons, xsout=*',xsout,nnn)
        call prin2('after newtons, wsout=*',wsout,nnn)
        call prin2('and corresponding accuracy is*',caccur,1)
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
        subroutine nodes0(n,w,ab,nn,k,x,whts,amatr,
     1      bmatr,eps,xsout,nmult,
     2      rnorms,ipivots,iwhich)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),ab(2,1),x(1),whts(1),amatr(n,1),
     1       bmatr(n,1),rnorms(1),ipivots(1),
     2       iwhich(1),xsout(1)
c 
c       This subroutine uses a perverted vesrion of the
c       Gram-Schmidt process to construct a first approximation
c       to the quadrature nodes for the user-provides set of functions.
c       The functions and the interval on which they are to be
c       integrated are specified by the parameters w,ab,nn,k;
c       It is assumed that these parameters have been created
c       by a prior call to the subroutine allsvbld (see). It is
c       also expected that the output of this subroutine (nodes
c       xsout) will be used by the subroutine homot0 (see). To the
c       best of the author's knowledge, this subroutine has no
c       function as a stand-alone device.
c 
c                     Input parameters:
c 
c 
c  n - the number of functions for which the quadratures are to be
c       designed. Must be even
c  w - the nested Legendre expansions of the basis functions;
c      normally, produced by a prior call to the subroutine allsvbld
c  ab - the array of ends of subintervals created. For each i=1,2,...,nn,
c      ab(1,i) is the left end of the i-th subinterval; ab(2,i) is
c      the right end of the i-th subinterval; normally, produced by a
c      prior call to the subroutine allsvbld
c  nn - the number of subintervals created on the interval [a,b];
c     also the number of entries in the array ab (see above); normally,
c      produced by a prior call to the subroutine allsvbld
c  k - the number of Gaussian nodes on each of nn subintervals in the
c      nested subdivision of the interval [a,b].
c  eps - the precision to which the Gram-Schmidt process is to be
c      conducted
c  nmult - the oversampling coefficient to be used by the subroutine
c      in the construction of array x. Explanation: The subroutine
c      allsvbld (see) discretizes adaptively the interval [a,b]; the
c      array ab and the parameters nn,k define the discretization.
c      This subrouytine needs a tabulation of the functions to be
c      integrated that resolves them "reasonably well". This is
c      accomplished by taking the discretization generated by allsvbld
c      and oversamplinng  it by the factor nmult, on each subinterval.
c      Normally, we use nmult=8 or so.
c 
c                     Output parameters:
c 
c  xsout - the approximate location of the quadrature nodes (n/2 of them)
c      for the integration of the n input funcctions.
c 
c                     Work arrays:
c 
c  x, whts - must be at least nn*k*nmult real *8 elements each
c  amatr,bmatr - must be at least nn*k*nmult*n real *8 elements each
c  irnorms,ipivots, iwhich - must be at least nn*k*nmult integer *4
c      elements each
c 
c        . . . oversample the discretization of the interval by an
c              appropriate factor
c 
        call nodewhts(ab,nn,k*nmult,x,whts)
c 
        numpts=nn*k*nmult
c 
c       construct the matrix of values and derivatives of the basis
c       functions at the nodes of the discretization
c 
        do 1400 i=1,numpts
        do 1200 j=1,n
c 
        call nesteva(ier,w,ab,nn,k,x(i),j,f,der)
c 
        amatr(j,i)=f
        bmatr(j,i)=der
c 
 1200 continue
 1400 continue
c 
c       use the pivoted Chebychev to find the approximate
c       quarature nodes
c 
        call chebpiv(amatr,bmatr,n,numpts,rnorms,eps,
     1      ncols,ipivots,iwhich)
c 
c        reorder the obtained indices and return to the user the
c        corresponding nodes
c 
        call prinf('after chebpiv, ncols=*',ncols,1)
        call prinf('while n=*',n,1)
        call bubblint(ipivots,ncols)
c 
        do 2200 i=1,ncols
        xsout(i)=x(ipivots(i))
 2200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine bubblint(ia,n)
        implicit real *8 (a-h,o-z)
        save
        dimension ia(1)
c 
        do 1400 i=1,n
        do 1200 j=1,n-1
        if(ia(j) .lt. ia(j+1)) goto 1200
c 
        k=ia(j)
        ia(j)=ia(j+1)
        ia(j+1)=k
 1200 continue
 1400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine chebpiv(a,b,n,m,rnorms,eps,ncols,ipivots,iwhich)
        implicit real *8 (a-h,o-z)
        save
        dimension b(n,m),a(n,m),rnorms(2,1),ipivots(1),
     1      iwhich(1)
c 
c        For two user-supplied matrices a(n,m), b(n,m), this
c        subroutine conducts a weird version of the pivoted
  
  
c 
c        . . . initialize the array of pivots
c 
        do 1200 i=1,m
        ipivots(i)=i
 1200 continue
c 
c        orthogonalize the pairs of vectors a(i),b(i) to each other
c 
        do 1600 i=1,m
c 
        call pivscapr(a(1,i),a(1,i),n,da)
        call pivscapr(b(1,i),b(1,i),n,db)
        call pivscapr(a(1,i),b(1,i),n,dab)
c 
        if(db .gt. da) goto 1400
c 
        d=dab/da
        call  svdsubt(a(1,i),b(1,i),d,n)
        iwhich(i)=1
        goto 1600
c 
 1400 continue
c 
        d=dab/db
        call  svdsubt(b(1,i),a(1,i),d,n)
        iwhich(i)=2
c 
 1600 continue
c 
c       . . . prepare the array of values of norms
c             of columns
c 
        done=1
        dtot=0
        do 1800 i=1,m
        da=0
        db=0
        do 1700 j=1,n
        da=da+a(j,i)**2
        db=db+b(j,i)**2
 1700 continue
        rnorms(1,i)=da
        rnorms(2,i)=db
        dtot=dtot+da+db
 1800 continue
c 
c       . . . conduct the weird gram-schmidt process
c 
        thresh=dtot*eps**2
cccc        do 5000 i=1,m
        do 5000 i=1,n/2
c 
c       find the pivot
c 
        ipivot=i
        rn=rnorms(1,i)+rnorms(2,i)
c 
        do 2200 j=i+1,m
        if(rnorms(1,j)+rnorms(2,j) .le. rn) goto 2200
        rn=rnorms(1,j)+rnorms(2,j)
        ipivot=j
 2200 continue
 2400 continue
c 
c       put the columns number ipivot in matrices a,b
c       in the i-th place
c 
  
cccc        call prinf('i=*',i,1)
cccc        call prinf('ipivot=*',ipivot,1)
  
        do 2600 j=1,n
        d=a(j,i)
        a(j,i)=a(j,ipivot)
        a(j,ipivot)=d
c 
        d=b(j,i)
        b(j,i)=b(j,ipivot)
        b(j,ipivot)=d
 2600 continue
c 
        iijj=ipivots(i)
        ipivots(i)=ipivots(ipivot)
        ipivots(ipivot)=iijj
c 
        d=rnorms(1,ipivot)
        rnorms(1,ipivot)=rnorms(1,i)
        rnorms(1,i)=d
c 
        d=rnorms(2,ipivot)
        rnorms(2,ipivot)=rnorms(2,i)
        rnorms(2,i)=d
c 
c       orthogonalize the i-th column to all preceeding ones
c 
        if(i .eq. 1) goto 3000
        do 2800 j=1,i-1
c 
        call pivscapr(a(1,j),a(1,i),n,d)
        call svdsubt(a(1,j),a(1,i),d,n)
c 
        call pivscapr(a(1,j),b(1,i),n,d)
        call svdsubt(a(1,j),b(1,i),d,n)
c 
c 
        call pivscapr(b(1,j),a(1,i),n,d)
        call svdsubt(b(1,j),a(1,i),d,n)
c 
        call pivscapr(b(1,j),b(1,i),n,d)
        call svdsubt(b(1,j),b(1,i),d,n)
c 
 2800 continue
 3000 continue
c 
c       orthonormalize the i-th columns of a,b
c 
        call pivscapr(a(1,i),a(1,i),n,da)
        call pivscapr(b(1,i),b(1,i),n,db)
        call pivscapr(a(1,i),b(1,i),n,dab)
c 
        if(db .gt. da) goto 3400
c 
        d=dab/da
        call  svdsubt(a(1,i),b(1,i),d,n)
        iwhich(i)=1
        goto 3600
c 
 3400 continue
c 
        d=dab/db
        call  svdsubt(b(1,i),a(1,i),d,n)
        iwhich(i)=2
c 
 3600 continue
c 
c       if the norm of the remaining part of the matrix
c       is sufficiently small - terminate the process
c 
        if(da+db .lt. thresh) return
  
        ncols=i
c 
        call pivscapr(a(1,i),a(1,i),n,da)
        call pivscapr(b(1,i),b(1,i),n,db)
c 
        da=1/sqrt(da)
        db=1/sqrt(db)
c 
        do 3800 j=1,n
        a(j,i)=a(j,i)*da
        b(j,i)=b(j,i)*db
 3800 continue
c 
c        orthogonalize everything else to it
c 
        do 4200 j=i+1,m
c 
        if(rnorms(1,j) + rnorms(2,j) .lt. thresh) goto 4200
c 
        call pivscapr(a(1,i),a(1,j),n,d)
        call  svdsubt(a(1,i),a(1,j),d,n)
c 
        call pivscapr(a(1,i),b(1,j),n,d)
        call  svdsubt(a(1,i),b(1,j),d,n)
c 
        call pivscapr(b(1,i),a(1,j),n,d)
        call  svdsubt(b(1,i),a(1,j),d,n)
c 
        call pivscapr(b(1,i),b(1,j),n,d)
        call  svdsubt(b(1,i),b(1,j),d,n)
c 
        call pivscapr(a(1,j),a(1,j),n,da)
        call pivscapr(b(1,j),b(1,j),n,db)
c 
        rnorms(1,j)=da
        rnorms(2,j)=db
 4200 continue
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
        subroutine svdsubt(x,y,d,n)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1)
c 
        do 1200 i=1,n
        y(i)=y(i)-d*x(i)
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine svdadd(x,y,d,n)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1)
c 
        do 1200 i=1,n
        y(i)=y(i)+d*x(i)
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine pivscapr(x,y,n,prod)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1)
c 
        prod=0
        do 1200 i=1,n
        prod=prod+x(i)*y(i)
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine newtons(ier,n,w,ab,nn,k,amatr,xs0,
     1      numit,bmatr,bmatr2,xs,ws,eps,xsout,wsout,
     2      dirmass0,rints,nitertot,
     3      uu,vv,tt,rnorms,ww)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),ab(2,1),amatr(n*2,n*2),
     1      xs0(1),bmatr(n,n),sol(1000),
     2      ws(1000),xs(1000),bmatr2(n,n*2),rhs(1000),
     3      xsout(1000),wsout(1000),rints(1)
c 
        dimension uu(1),vv(1),tt(1),rnorms(1),ww(1)
c 
c       This subroutine evaluates (or rather, attempts to evaluate)
c       n Generalized Gaussian nodes on the interval [a,b]. Please
c       note that this subroutine uses the parameters w, ab, nn, etc.
c       produced by a previous call to the subroutine allsvbld (see);
c       this subroutine has no function as a stand-alone device.
c 
c       This subroutine conducts a homotopy IN THE WEIGHT FUNCTION
c       of the quadrature; the functions to be integrated are not
c       changed. In the beginning, the subroutine
c 
c                          Input parameters:
c 
c  n - the number of nodes to be otained
c  w - the nested Legendre expansions of the basis functions;
c      normally, produced by a prior call to the subroutine allsvbld
c  ab - the array of ends of subintervals created. For each i=1,2,...,nn,
c      ab(1,i) is the left end of the i-th subinterval; ab(2,i) is
c      the right end of the i-th subinterval; normally, produced by a
c      prior call to the subroutine allsvbld
c  nn - the number of subintervals created on the interval [a,b];
c     also the number of entries in the array ab (see above); normally,
c      produced by a prior call to the subroutine allsvbld
c  k - the number of Gaussian nodes on each of nn subintervals in the
c      nested subdivision of the interval [a,b].
c  x - the discretization of the interval [a,b] consisting of nn
c        subintervals with k Legendre nodes on each
c  whts - the weights associated with the discretization x; (normally
c        produced by the subroutine manincon (see)
c  xs0 - the initial positions of the nodes to be used by the
c       continuation process
c  numit - the number of iterations the Newton process is permitted
c      to perform for every value of the continuation parameter before
c      it surrenders
c  eps - the accuracy to which the SVD of the input functions will be
c      constructed; the best accuracy obtainable with the resulting
c      quadratures is usually at least 100 times worse
c  dirmass0 - the Dirac masses located at the original (starting)
c      starting point nodes at the beginning of the homotopy. A good
c      value: 1.0e-6
c 
c                      Output prameters:
c 
c  ier - error return code;
c              ier=0 means successful conclusion
c              ier=16 means that at some point, the subroutine had to
c                     reduce the step in the homotopy to under 1.0d-5;
c                     the homotopy has failed
c              ier=64 means that the subroutine failed to get to
c                     the other end of the homotopy after 500 steps of
c                     continuation. This is an extremely unlikely
c                     situation
c              ier=128 means that the first attempt to use Newton failed.
c                     The initial approximation to the nodes is terribly
c                     out of wack with the functions to be integrated
c  xs - the nodes of the obtained quadrature
c  ws - the weghts of the obtained quadrature
c 
c                      Work arrays:
c 
c  uu,vv,tt,ww - must be at least n**2 * 4 real *8 elements each
c  rnorms - must be at least n*2 real *8 elements
c 
        ier=0
        dirmass=dirmass0
        nitertot=0
c 
cccc        call prin2('in newtons, xs0=*',xs0,n)
        call prin2('in newtons, dirmass0=*',dirmass0,1)
c 
c       construct the matrix of values of all functions at the
c       original (initialization) nodes, to be used to
c       stabilize the Newton process
c 
  
        do 1140 i=1,n*2
        do 1120 j=1,n
c 
        call nesteva(ier,w,ab,nn,k,xs0(j),i,f,der)
c 
        bmatr2(j,i)=f
 1120 continue
 1140 continue
c 
c       find the weights for the initialization of the Newton
c 
        do 1400 i=1,n
        do 1200 j=1,n
c 
        call nesteva(ier,w,ab,nn,k,xs0(j),i,f,der)
c 
        bmatr(i,j)=f
 1200 continue
 1400 continue
c 
        numpts=nn*k
c 
        do 2400 i=1,n
        d=0
        do 2300 j=1,n
        d=d+bmatr2(j,i)
 2300 continue
c 
        rhs(i)=rints(i)+d*dirmass
c 
 2400 continue
c 
  
        call grmsol(bmatr,n,rhs,ws,rcond)
  
  
  
cccc        call leastsq(bmatr,uu,ww,tt,n,n,ncols,rnorms,eps,vv)
c 
cccc        call leastsq2(uu,ww,tt,n,n,ncols,rhs,ws,vv)
c 
        do 2600 i=1,n
        xs(i)=xs0(i)
 2600 continue
c 
cccc        call prin2('initializing Newton, rnorms=*',rnorms,ncols)
cccc        call prin2('initializing Newton, ws=*',ws,n)
cccc        call prin2('while xs=*',xs,n)
  
        coef=0.1
c 
        do 4000 ijk=1,500
c 
        dirmass2=dirmass/(1+coef)
c 
c       conduct the Newton process
c 
        dirmass2=dirmass/(1+coef)
        eps2=eps*10
c 
        call newton(jer,n,w,ab,nn,k,amatr,xs,ws,sol,
     1      numit,niter,bmatr2,xsout,wsout,dirmass2,eps2,
     2      rints,uu,vv,tt,rnorms,ww)
c 
        nitertot=nitertot+niter
c 
       if(jer .ne. 0) goto 3000
c 
       do 2800 i=1,n
       xs(i)=xsout(i)
       ws(i)=wsout(i)
 2800 continue
c 
       dirmass=dirmass2
       coef=coef*2
       if(dirmass .lt. 1) coef=coef*2
c 
       goto 3600
c 
 3000 continue
c 
c       if this is the first time we are trying to conduct Newton,
c       and it failed - bomb immediately
c 
        if(ijk .ne. 1) goto 3200
        ier=128
        return
 3200 continue
c 
       coef=coef/2
c 
       eps3=1.0d-5
  
       if(coef .gt. eps3) goto 4000
       ier=16
       return
c 
 3600 continue
c 
       if(dirmass .lt. 1.0d-20) return
 4000 continue
c 
       ier=64
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine newton(ier,n,w,ab,nn,k,amatr,xs0,ws0,sol,
     1      numit,niter,bmatr2,xs,ws,dirmass,eps,
     2      rints,uu,vv,tt,rnorms,ww)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),ab(2,1),amatr(n*2,n*2),
     1      ws0(1),xs0(1),sol(1),
     2      ws(1),xs(1),bmatr2(n,n*2),rhs(1000),rints(1)
c 
        dimension uu(1),vv(1),tt(1),rnorms(1),ww(1)
c 
        ier=0
c 
        do 1200 i=1,n
        xs(i)=xs0(i)
        ws(i)=ws0(i)
 1200 continue
c 
c       conduct Newton iterations
c 
        stepold=1.0d20
  
        do 3400 ijk=1,numit
  
        niter=ijk
  
        call nodegau(n,w,ab,nn,k,amatr,xs,ws,sol,
     1      bmatr2,dirmass,rhs,rints,
     2      uu,vv,tt,rnorms,ww)
c 
        do 2600 i=1,n
        xs(i)=xs(i)-sol(i+n)
  
        if(xs(i) .lt. ab(1,1) ) xs(i)=ab(1,1)
  
        ws(i)=ws(i)-sol(i)
 2600 continue
c 
c       if the solution is sufficiently small - exit
c 
        d=0
        do 2800 i=1,n
        d=d+sol(i+n)**2
 2800 continue
c 
        step=sqrt(d)
        if(step .lt. stepold*10) goto 3000
cccc        if(step .lt. stepold*2) goto 3000
c 
        ier=8
        return
 3000 continue
c 
        stepold=step
c 
        if(step .lt. eps) return
 3400 continue
c 
        ier=4
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine nodegau(n,w,ab,nn,k,amatr,xs,ws,sol,
     1      bmatr2,dirmass,rhs,rints,
     2      uu,vv,tt,rnorms,ww)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),ab(2,1),amatr(n*2,n*2),
     1      ws(1),xs(1),rhs(1),sol(1),bmatr2(n,n*2),
     2      rints(1)
c 
        dimension uu(1),vv(1),tt(1),rnorms(1),ww(1)
c 
c 
c       construct the Jacobian matrix
c 
        do 1400 i=1,n
        do 1200 j=1,n*2
c 
        call nesteva(ier,w,ab,nn,k,xs(i),j,f,der)
c 
        amatr(j,i)=f
        amatr(j,n+i)=der*ws(i)
 1200 continue
 1400 continue
c 
c        construct the right-hand sides
c 
        numpts=nn*k
c 
        do 1800 i=1,n*2
c 
        d=0
        do 1700 j=1,n
c 
        d=d+bmatr2(j,i)
 1700 continue
c 
        rhs(i)=rints(i)+d*dirmass
c 
 1800 continue
c 
        do 2400 i=1,n*2
        d=0
        do 2200 j=1,n
c 
        d=d+amatr(i,j)*ws(j)
c 
 2200 continue
c 
        rhs(i)=d-rhs(i)
 2400 continue
c 
c       perform the Newton step
c 
        call grmsol(amatr,n*2,rhs,sol,rcond)
c 
        eps7=1.0d-14
cccc        call leastsq(amatr,uu,ww,tt,n*2,n*2,ncols,rnorms,eps7,vv)
c 
cccc        call leastsq2(uu,ww,tt,n*2,n*2,ncols,rhs,sol,vv)
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine svdscapr(x,y,n,prod)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1)
c 
        prod=0
        do 1200 i=1,n
        prod=prod+x(i)*y(i)
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine xscoeret(ir,nn,ab,coefs,k,ncols,x,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension ab(2,1),coefs(nn*k*ncols),x(k*nn),
     1      whts(k*nn)
c 
        character *1 card(100)
c 
c       this subroutine retrieves from disk the coefficients
c       of nested Legendre expansions of a bunch of functions;
c       hopefully, these have been stored previously by the
c       subroutines xsstorei, xsstore, xssclose (see). The
c       subroutine also retrieves several other pieces of data
c       pertaining to the expansions stored on disk.
c 
c                  Input parameters:
c 
c  ir - the FORTRAN unit number from which the data are to be read
c 
c                  Output prameters:
c 
c  nn - the number of subintervals created on the interval [a,b];
c     also the number of entries in the array ab (see above)
c  ab - the array of ends of subintervals created. For each i=1,2,...,nn,
c     ab(1,i) is the left end of the i-th subinterval; ab(2,i) is
c     the right end of the i-th subinterval
c  coefs - the nested Legendre expansions of the functions stored in the
c      array fs.
c  k - the number of Gaussian nodes on each of nn subintervals in the
c  ncols - the number of singular functions whose expansions are
c      stored on the unit ir
c  x - the discretization of the interval [a,b] consisting of nn
c        subintervals with k Legendre nodes on each
c  whts - the weights associated with the discretization x; (normally
c        produced by the subroutine manincon (see)
c 
c       . . . find in the file ir the line where the integer data
c             are stored
c 
 1200 format(80a1)
c 
        ifrewind=1
        call txtfnd(ier,ir,'As returned',6,3,
     1      ifrewind,lineout)
c 
        read(ir,1200,end=6000) (card(j),j=1,30)
c 
 2400 format(11x,i6,5x,i3,8x,i3)
c 
        read(ir,2400) nn,k,ncols
c 
c       find in the file the line where the nodes of the large
c       discretization are stored
c 
        ifrewind=1
        call txtfnd(ier,ir,'nodes of the large',14,7,
     1      ifrewind,lineout)
c 
        read(ir,1200,end=6000) (card(j),j=1,30)
c 
        numpts=nn*k
c 
 2600 format(8x,e22.16,1x,e22.16)
c 
        kk=numpts/2
        ii=numpts-kk*2
        do 2800 i=1,kk
c 
        read(ir,2600) x(2*i-1),x(2*i)
 2800 continue
c 
        if(ii .ne. 0) read(ir,2600) x(numpts)
c 
c       find in the file the line where the weights of the large
c       discretization are stored
c 
        ifrewind=1
        call txtfnd(ier,ir,'weights of the large',15,7,
     1      ifrewind,lineout)
c 
        read(ir,1200,end=6000) (card(j),j=1,30)
c 
        numpts=nn*k
c 
 3000 format(8x,e22.16,1x,e22.16)
c 
        kk=numpts/2
        ii=numpts-kk*2
        do 3200 i=1,kk
c 
        read(ir,3000) whts(2*i-1),whts(2*i)
 3200 continue
c 
        if(ii .ne. 0) read(ir,3000) whts(numpts)
c 
c       find in the file the line where the subintervals of the
c       discretization are stored
c 
        ifrewind=1
        call txtfnd(ier,ir,'ends of the subintervals',18,7,
     1      ifrewind,lineout)
c 
        read(ir,1200,end=6000) (card(j),j=1,30)
c 
        do 3400 i=1,nn
c 
        read(ir,3000) ab(1,i),ab(2,i)
 3400 continue
c 
c       find in the file the line where the coefficients of all
c       expansions are stored
c 
        ifrewind=1
        call txtfnd(ier,ir,'The coefficients',15,3,
     1      ifrewind,lineout)
c 
        read(ir,1200,end=6000) (card(j),j=1,30)
c 
        numcoefs=nn*k*ncols
c 
 3600 format(8x,e22.16,1x,e22.16)
c 
        kk=numcoefs/2
        ii=numcoefs-kk*2
        do 3800 i=1,kk
c 
        read(ir,3600) coefs(2*i-1),coefs(2*i)
 3800 continue
c 
        if(ii .ne. 0) read(ir,3600) coefs(numcoefs)
c 
 6000 continue
  
        return
        end
c 
c 
c 
c 
c 
        subroutine xssretr(ier,ir,n,x,ws)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),ws(1)
        character *1 card(100)
c 
c       This subroutine retrieves from disk the nodes and weights
c       of a single quadrature formula
c 
c               Input parameters:
c 
c  ir - the FORTRAN unit number from which the data are to be read
c  n - the number of nodes in the quadrature to be retrieved
c 
c               Output parameters:
c 
c  x - the n nodes
c  ws - the n weights corresponding to the nodes x
c 
c       . . . find in the file the line where the requested
c             quadrature nodes are stored
c 
        ier=0
c 
        ifrewind=1
        do 1200 i=1,n
c 
        call txtfnd(ier,ir,'Nodes',5,9,ifrewind,lineout)
        ifrewind=0
 1200 continue
c 
 1400 format(80a1)
c 
c        found the record number n. read it
c 
        read(ir,1400,end=6000) (card(j),j=1,30)
c 
        do 2600 i=1,n/2
c 
 2400 format(8x,e22.16,1x,e22.16)
c 
        read(ir,2400) x(2*i-1),x(2*i)
 2600 continue
c 
        n2=n/2
        j=n-n2*2
        if(j .ne. 0) read(ir,2400) x(n)
c 
c       find in the file the line where the requested quadrature
c       weights are stored
c 
        ifrewind=1
        do 3200 i=1,n
c 
        call txtfnd(ier,ir,'Weights',5,9,ifrewind,lineout)
c 
cccc        call prinf('after txtfnd, ier=*',ier,1)
        ifrewind=0
 3200 continue
c 
cccc 3400 format(80a1)
c 
c        found the record number n. read it
c 
        read(ir,1400,end=6000) (card(j),j=1,30)
c 
        do 3600 i=1,n/2
c 
 3400 format(8x,e22.16,1x,e22.16)
c 
        read(ir,3400) ws(2*i-1),ws(2*i)
 3600 continue
c 
        n2=n/2
        j=n-n2*2
        if(j .ne. 0) read(ir,3400) ws(n)
c 
        ier=0
        return
c 
 6000 continue
        ier=1024
        return
        end
c 
c 
c 
c 
c 
        subroutine xssclose(iw,nn,ab,coefs,k,ncols,x,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension ab(2,1),coefs(nn*k*ncols),x(k*nn),
     1      whts(k*nn)
c 
c       this subroutine stores on disk the coefficients
c       of nested Legendre expansions of a bunch of functions;
c       hopefully, these will be read later by the
c       subroutines xssretr, xsscoeret (see). The
c       subroutine also stores several other pieces of data
c       pertaining to the expansions stored on disk.
c 
c                  Input parameters:
c 
c  ir - the FORTRAN unit number on which the data are to be stored
c  nn - the number of subintervals created on the interval [a,b];
c     also the number of entries in the array ab (see above)
c  ab - the array of ends of subintervals created. For each i=1,2,...,nn,
c     ab(1,i) is the left end of the i-th subinterval; ab(2,i) is
c     the right end of the i-th subinterval
c  coefs - the nested Legendre expansions of the functions stored in the
c      array fs.
c  k - the number of Gaussian nodes on each of nn subintervals in the
c  ncols - the number of singular functions whose expansions are
c      stored on the unit ir
c  x - the discretization of the interval [a,b] consisting of nn
c        subintervals with k Legendre nodes on each
c  whts - the weights associated with the discretization x; (normally
c        produced by the subroutine manincon (see)
c 
c                  Output prameters: none
c 
c       . . . store on disk the information to be used in the
c             evaluation of obtained singular functions
c 
 1200 format('c        ')
 1400 format('c  As returned from allqrbld, ',/,'c   ',/,
     1     6x,' nn= ',i6,',  k=',i3,', ncols=',i3)
c 
        write(iw,1200)
        write(iw,1400) nn,k,ncols
c 
        write(iw,1200)
c 
 1600 format('c  The nodes of the large discretization are')
c 
        write(iw,1600)
        write(iw,1200)
c 
        numpts=nn*k
c 
 2600 format(5x,i1,2x,e22.16,',',e22.16,',')
c 
        kk=numpts/2
        ii=numpts-kk*2
        do 2800 i=1,kk
c 
        jj=i/10
        j=i-jj*10
        if(j .eq. 0) j=10
c 
        write(iw,2600) j,x(2*i-1),x(2*i)
 2800 continue
c 
        if(ii .ne. 0) write(iw,2600) j+1,x(nn)
c 
 3200 format('c  The weights of the large discretization are')
c 
        write(iw,1200)
        write(iw,1200)
        write(iw,3200)
        write(iw,1200)
c 
        numpts=nn*k
c 
 3400 format(5x,i1,2x,e22.16,',',e22.16,',')
c       kk=numpts/2
        ii=numpts-kk*2
        do 3800 i=1,kk
c 
        jj=i/10
        j=i-jj*10
        if(j .eq. 0) j=10
c 
        write(iw,3400) j,whts(2*i-1),whts(2*i)
 3800 continue
c 
        if(ii .ne. 0) write(iw,3400) j+1,whts(nn)
c 
        write(iw,1200)
        write(iw,1200)
c 
c 
 4200 format('c  The ends of the subintervals are')
c 
        write(iw,4200)
        write(iw,1200)
c 
        do 4400 i=1,nn
c 
        jj=i/10
        j=i-jj*10
        if(j .eq. 0) j=10
c 
c 
        write(iw,3400) j,ab(1,i),ab(2,i)
c 
 4400 continue
c 
c 
 5200 format('c  The coefficients are')
c 
        write(iw,1200)
        write(iw,5200)
        write(iw,1200)
c 
        ncoefs=nn*k*ncols
        kk=ncoefs/2
        ii=ncoefs-kk*2
c 
        do 5400 i=1,kk
c 
        jj=i/10
        j=i-jj*10
        if(j .eq. 0) j=10
c 
        write(iw,3400) j,coefs(i*2-1),coefs(i*2)
c 
 5400 continue
c 
        if(jj .eq. 1) write(iw,3400) j+1,coefs(ncoefs)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine xsstorei(iw,ncols,a,b,rlams)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),ws(1),rlams(1)
  
c 
c       store the header of this run
c 
 1200 format('c        ')
 1800 format('c     In this run, the singular values are:')
        write(iw,1200)
        write(iw,1800)
        write(iw,1200)
c 
 1900 format('c',6(2x,e11.5))
c 
        n=ncols/6
        ii=ncols-n*6
c 
        do 2000 i=1,n
c 
        jj=(i-1)*6
        write(iw,1900) (rlams(jj+j),j=1,6)
 2000 continue
c 
        jj=n*6
        if(ii .ne. 0)
     1      write(iw,1900) (rlams(jj+j),j=1,ii)
c 
        write(iw,1200)
c 
        return
c 
c 
c 
c 
        entry xsstore(iw,xs,ws,ns)
c 
 2200 format('c        ')
        write(iw,2200)
c 
 2400 format('c  Data for ',i3,' nodes')
        write(iw,2400) ns
        write(iw,2200)
c 
 2500 format('c        Nodes:')
        write(iw,2500)
        write(iw,2200)
c 
 2600 format(5x,i1,2x,e22.16,',',e22.16,',')
c 
        n=ns/2
        ii=ns-n*2
        do 2800 i=1,n
c 
        write(iw,2600) i,xs(i*2-1),xs(i*2)
 2800 continue
c 
        if(ii .ne. 0) write(iw,2600) n+1,xs(ns)
c 
        write(iw,2200)
 3200 format('c        Weights:')
        write(iw,3200)
        write(iw,2200)
c 
        do 3400 i=1,n
c 
        write(iw,2600) i,ws(i*2-1),ws(i*2)
 3400 continue
c 
        if(ii .ne. 0) write(iw,2600) n+1,ws(ns)
c 
        write(iw,2200)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine txtfnd(ier,ir,text,ltext,istart,
     1      ifrewind,lineout)
        implicit real *8 (a-h,o-z)
        save
        character *1 text(1),card(800),blank
        data blank/' '/
c 
c        This subroutine finds in the file ir the nearest line on
c        which the string text of length ltext can be found in
c        the position specified by the parameter istart. The
c        subroutine either rewinds the unit ir before the search
c        or not, depending
c        on the value of the user-specified parameter ifrewind.
c        The subroutine does two things:
c 
c    1. It positions the unit ir to read the line {\it following}
c        the line where the text has been found
c    2. It returns to the user the parameter lineout telling him
c        how many lines it has read to get to the text to be found
c        (inclusive)
c 
c                 Input parameters:
c 
c  ir - the FORTRAN unit number on which the search is to be
c        performed
c  text - the string to be found
c  ltext - the length of the string to be found
c  istart - the location of the string on the line: the first character
c        of the istring will be looked for in position istart+1
c  ifrewind - tells the subroutine whether the unit ir is to be rewound
c        the search: ifrewind=1 will cause the unit to be rewound
c                    ifrewind=0 will prevent the rewinding
c 
c                 Output parameters:
c 
c  ier - tells whether the subroutine has found the text:
c           ier=1 means that
c  lineout - the number of lines the subroutine has read to find
c           the requested text
c 
c       . . . find in the file ir the line on which the text
c             "text" is found, starting with the position
c             istart+1
c 
        ier=8
        if(ifrewind .eq. 1) rewind(ir)
        do 2000 i=1,1000000
c 
 1200 format(80a1)
c 
        do 1400 j=1,130
        card(j)=blank
 1400 continue
c 
        read(ir,1200,end=3000) (card(j),j=1,30)
c 
        iffound=0
        do 1600 j=1,ltext
        if(card(j+istart) .ne. text(j)) goto 2000
 1600 continue
c 
        lineout=i
        ier=0
        return
c 
 2000 continue
c 
        ier=1024
 3000 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine txtwrite(iw,nmess,t1,t2,t3,t4,t5,t6,t7,t8,t9,
     1      t10,t11,t12,t13,t14,t15,t16,t17,t18,t19)
        implicit real *8 (a-h,o-z)
        save
        character *1 t1(1),t2(1),t3(1),t4(1),t5(1),t6(1),t7(1),
     1      string(100 000),blank
        data blank/' '/
c 
c        This subroutine writes out on disk a text provided to it by
c        the user in the
c 
c        concatenate all the message lines supplied by the
c        user into one long string
c 
        j=1
        call textlen(t1,nchar,string(j))
        j=j+nchar
        if(nmess .eq. 1) goto 3000
c 
        call textlen(t2,nchar,string(j))
        j=j+nchar
        if(nmess .eq. 2) goto 3000
c 
        call textlen(t3,nchar,string(j))
        j=j+nchar
        if(nmess .eq. 3) goto 3000
c 
        call textlen(t4,nchar,string(j))
        j=j+nchar
        if(nmess .eq. 4) goto 3000
c 
        call textlen(t5,nchar,string(j))
        j=j+nchar
        if(nmess .eq. 5) goto 3000
c 
        call textlen(t6,nchar,string(j))
        j=j+nchar
        if(nmess .eq. 6) goto 3000
c 
        call textlen(t7,nchar,string(j))
        j=j+nchar
        if(nmess .eq. 7) goto 3000
c 
        call textlen(t8,nchar,string(j))
        j=j+nchar
        if(nmess .eq. 8) goto 3000
c 
        call textlen(t9,nchar,string(j))
        j=j+nchar
        if(nmess .eq. 9) goto 3000
c 
 3000 continue
c 
 3100 format('c              ')
        write(iw,3100)
c 
c       format the obtained long line into a sequence of comment
c       lines of FORTRAN
c 
        nchartot=j-1
c 
cccc        call prina('in txtwrite, string=*',string,j-1)
c 
        nlines=0
        j=0
        linelen=0
        linestrt=1
        do 3400 i=1,10000
c 
        j=j+1
c 
        if(j .gt. nchartot) goto 3600
c 
        if(string(j) .eq. blank) lastblnk=j
c 
        linelen =linelen+1
        if(linelen .ne. 60) goto 3400
c 
c       we have come to the end of the line. Print this
c       line and start another one
c 
 3200 format('c',5x,4x,80a1)
        write(iw,3200) (string(jj),jj=linestrt,lastblnk)
c 
        linestrt=lastblnk+1
        nlines=nlines+1
        linelen=0
 3400 continue
c 
 3600 continue
c 
c       If there is an incomplete line left - print it
c 
        write(iw,3200) (string(jj),jj=linestrt,nchartot)
        write(iw,3100)
c 
        return
        end
c 
c 
c 
c 
c 
        SUBROUTINE textlen(MES,nchar,line)
        save
        CHARACTER *1 MES(1),AST,line(1)
        DATA AST/'*'/
C 
C         DETERMINE THE LENGTH OF THE MESSAGE
C 
        I=0
        DO 1400 I=1,10000
        IF(MES(I).EQ.AST) GOTO 1600
        I1=I
 1400 CONTINUE
 1600 CONTINUE
c 
        nchar=i1
        do 1800 i=1,nchar
        line(i)=mes(i)
 1800 continue
         RETURN
         END
