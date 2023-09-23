        implicit real *8 (a-h,o-z)
        dimension ts(1000),whts(1000),w(1000 000)
  
        real *8 xs(10 000),xsaux(1000),whtsaux(2000),
     1      totmat(100 000)
c 
        external interact
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER ndigits'
        READ *,ndigits
        CALL PRINf('ndigits=*',ndigits,1 )
c 
        PRINT *, 'ENTER n'
        READ *,n
        CALL PRINf('n=*',n,1 )
  
        naux=n*6
  
  
        ir=17
        call uniptret(ier,ir,ndigits,n,naux,xs,xsaux,whts,whtsaux)
  
        call prinf('after uniptret, ier=*',ier,1)
        call prin2('after uniptret, xs=*',xs,n)
        call prin2('after uniptret, whts=*',whts,n)
        call prin2('after uniptret, xsaux=*',xsaux,naux)
        call prin2('after uniptret, whtsaux=*',whtsaux,naux)
  
  
cccc          stop
  
        ifinit=1
        call unisegm(ier,ir,ifinit,ndigits,n,naux,interact,
     1      totmat,w)
  
        if(2 .ne. 3) goto 3200
  
        ifinit=0
        call unisegm(ier,ir,ifinit,ndigits,n,naux,interact,
     1      totmat,w)
  
  
 3200 continue
  
  
cccc         stop
  
c 
c        test the obtained matrix
c 
  
        call mattest(totmat,n,ts,xs,whts)
  
  
  
        stop
        end
  
c 
c 
c 
c 
c 
        subroutine mattest(totmat,n,
     1      ts,xs,whts2)
        implicit real *8 (a-h,o-z)
        save
        dimension ts(1),sigma(1000),vals(1000),
     1      xs(1),whts2(1),
     2      vals2(1000),diffs(1000)
c 
        external funeval
  
        real *8 totmat(n,n)
c 
c        construct the values of the charge at the support nodes
c 
        do 1200 i=1,n
c 
        sigma(i)=sigmafun(xs(i))
 1200 continue
c 
c        evaluate the integral at the support nodes
c 
        call matvec(totmat,sigma,n,vals)
  
        call prin2('in mattest, vals are*',vals,n)
  
        iw=22
        call lotagraph(iw,xs,vals,n,'potential*')
  
c 
c       calculate the first several moments of the function vals
c 
        d0=0
        d1=0
        d2=0
        do 2200 i=1,n
c 
        d0=d0+vals(i)*whts2(i)
        d2=d2+vals(i)*whts2(i)*xs(i)**2
 2200 continue
  
        call prin2('in totmat, d0=*',d0,1)
        call prin2('in totmat, d2=*',d2,1)
c 
c        now, calculate the same values via adapgaus
c 
        eps=1.0d-12
        a=-1
        b=1
        mm=12
  
        do 2400 i=1,n
  
        call adapgaus(ier,a,b,funeval,xs(i),par2,mm,eps,
     1      rint,maxrec,numint)
  
  
        vals2(i)=rint
        diffs(i)=vals(i)-vals2(i)
 2400 continue
  
  
        call prin2('and after adapgaus, vals2=*',vals2,n)
        call prin2('and diffs=*',diffs,n)
  
  
        return
        end
c 
c 
c 
c 
c 
        function funeval(x,x0,dummy)
        implicit real *8 (a-h,o-z)
c 
cccc        call interact(x0,x,f)
        save
        call interact(j,i,x,x0,f)
c 
        funeval=f*sigmafun(x)
c 
        return
        end
c 
c 
c 
c 
c 
        function sigmafun(x)
        implicit real *8 (a-h,o-z)
c 
cccc        sigmafun=sin(x*3)
        save
  
        sigmafun=(x+1)*log(x+1)
  
        sigmafun=sigmafun*cos(6*x)
  
  
ccccc        sigmafun=sqrt(x+1)*log(x+1)*sin(x*3)
cccc        sigmafun=log(x+1)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine matvec(a,x,n,y)
        implicit real *8 (a-h,o-z)
        save
        dimension x(n),a(n,n),y(n)
c 
        do 1400 i=1,n
        d=0
        do 1200 j=1,n
        d=d+x(j)*a(i,j)
 1200 continue
        y(i)=d
 1400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine interact(j,i,xaux,xsup,f)
        implicit real *8 (a-h,o-z)
        save
  
        d=(xaux-xsup)**2
        f=log(d)
  
  
cccc        f=f*xaux
  
cccc        f=d
  
  
        return
        end
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c       This si the end of the debugging code and the beginning of the
c       quadrature code proper
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine unisegm(ier,ir,ifinit,ndigits,n,naux,interact,
     1      totmat,w)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),totmat(n,n)
c 
        external interact
  
c 
c        This subroutine constructs the matrix totmat of interactions
c        of a single segment with itself. The potential is assumed
c        to be of the form
c 
c                 G(x,y)=\phi(x,y)*log(||x-y||) +\psi(x,y),            (1)
c 
c        with \phi, \psi smooth functions. From the point of view of
c        this subroutine, the segment is parametrized by the interval
c        [-1,1], and it is the user's responsibility to scale all of
c        the elements of the matrix totmat (see below) accordingly,
c        and to design the subroutine interact (see below) in such a
c        way that with this in mind.
c 
c        The subroutine constructs an n-point discretization of the
c        interval [-1,1], as well as an auxiliary naux-point discre-
c        tization, and uses the user-supplied subroutine interact
c        to evaluate the interactions between the points of the
c        discretization and the auxiliary points. The locations on
c        the interval [-1,1] of both sets of points are highly
c        non-equispaced, and are read by the subroutine from the
c        file number ir. Certain other types of data are read from
c        the file ir, and it is fervently hoped that the user has
c        stored the required data on that file by preceding call(s)
c        to the subroutine oneseg. The locations of both the nodes
c        of the discretization and the support nodes can (and should)
c        be obtained by the user via a preceding call to the subroutine
c        uniptret (see); uniptret also provides the weights corresponding
c        to the discretization nodes.
c 
c        ATTENHUT! IMPORTANT NOTE! In the absence of the appropriate
c        file number ir, this subroutine is completely useless; it has
c        no use as a stand-alone device. Plese note that normally, the
c        file it is created long before they are used by this subroutine.
c        the contents of ir do not depend on the Green's function or on
c        much else; they are fairly "universal".
c 
c                      Input parameters:
c 
c  ir - the FORTRAN file from which the data are to be read.
c        Hopefully, these data have been stored there by previous
c        calls to the the subroutine oneseg
c  ifinit - the parameter telling the subroutine whether the data have
c        to be read from the file ir, or have been read from that file
c        during a preceding call.
c    ifinit=1 tells the subroutine to read the data
c    ifinit=0 tells the subroutine that the data are already in the
c        array w (see below). In this case, the first
c 
c                      2*n+naux+2*n*naux+20                               (2)
c 
c        elements of the array w have to be unchanged since the last
c        call to this subroutine
c  ndigits - the number of digits used in the design of the discretization
c        to be read from disk
c  n - the number of nodes in the discretization of the segment to be
c        constructed
c  naux - the number of auxiliary nodes to be used in the construction
c        of the discretization
c  intercat - the subroutine (user-supplied) evaluation the interaction
c        function (1). The calling sequence of interact is
c 
c               interact(j,i,xaux,xsup,f).                                (3)
c 
c****************************************************************************
c 
c            Input parameters of the subroutine interact:
c 
c  j - the sequence number of the auxiliary node at which the
c        source is located
c  i - the sequence number of the discretization node at which the
c        potential is to be evaluated
c  xaux - the location on the interval [-1,1] of the auxiliary node
c        at which the source is located
c  xsup -  the location on the interval [-1,1] of the discretization
c          node at which the potential is to be evaluated
c 
c            Output parameters of the subroutine interact:
c 
c  f - the potential created at the ponint xsup by the unit charge at
c        the point xaux
c 
c****************************************************************************
c 
c  w - the array in which the data are stored read during a preceding call
c        to this subroutine. Plese note that this is only an input parameter
c        if ifinit (see above) has been set to 0; otherwise, this is an
c        output parameter (or a work array,if this is the first and only
c        call to this subroutine). Also, please note that the length of w
c        should not be less than 2*n+naux+2*n*naux+20
c 
c                              Output parameters:
c 
c  ier - error return code.
c     ier=0 means successful execution;
c     ier=4 means that the file number IR does not contain data
c        corresponding to the three values (ndigit,n,naux) specified by
c        the user. The user should either select a different combuination,
c        or rerun the subroutine oneseg (currently, a part of the file
c        test35.f) in order to add these data to the file IR: Do you think
c        the subroutine oneseg is a conversation piece, or what?
c  totmat - the matrix of interactions of nodes xs with themselves
c 
c        . . . allocate memory for various arrays
c 
        ixs=1
        lxs=n+2
c 
        iwhts=ixs+lxs
        lwhts=n+2
c 
        ixsaux=iwhts+lwhts
        lxsaux=naux+2
c 
        iwhtsaux=ixsaux+lxsaux
        lwhtsaux=naux+2
c 
        icorrs=iwhtsaux+lwhtsaux
        lcorrs=n*naux+2
c 
        iainterp=icorrs+lcorrs
        lainterp=n*naux+2
c 
        iautosup=iainterp+lainterp
        lautosup=n*naux+2
c 
        ltot=iautosup+lautosup
c 
c       if the user so requested - retrieve from disk the
c       various arrays and matrices
c 
        if(ifinit .ne. 0)
     1      call unicoret(ier,ir,ndigits,n,naux,
     2           w(ixs),w(ixsaux),w(iwhts),w(iwhtsaux),
     3      w(icorrs),w(iainterp) )
c 
        if(ier .ne. 0) return
c 
c        construct the matrix of interactions
c 
        call onematr(n,naux,w(icorrs),w(iainterp),
     1      w(ixs),w(ixsaux),totmat,w(iautosup),interact)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine uniptret(ier,ir,ndigits,n,
     1      naux,xs,xsaux,whts,whtsaux)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),xsaux(1),whts(1),whtsaux(1)
c 
        character *1 card(70),template(1)
        character *100 image
        character * 26 templ20
c 
        equivalence(templ20,template(1))
c 
        data templ20/'      In this run, ndigits'/
c 
c       This subroutine reads from the file ir the data believed
c       necessary to the user in dealing with the subroutine oesegm
c       (see). The data are the nodes of the "universal" nodes on
c       the interval [-1,1], the weights whts of the quadrature
c       associated with these "universal" nodes, and the auxiliary
c       nodes xsaux.
c 
c 
c                      Input parameters:
c 
c  ir - the FORTRAN file from which the data are to be read.
c        Hopefully, these data have been stored there by previous
c        calls to the the subroutine oneseg
c  ndigits - the number of digits used in the design of the discretization
c        to be read from disk
c  n - the number of nodes in the discretization of the segment to be
c        constructed
c  naux - the number of auxiliary nodes to be used in the construction
c        of the discretization
c 
c                              Output parameters:
c 
c  ier - error return code.
c     ier=0 means successful execution;
c     ier=4 means that the file number IR does not contain data
c        corresponding to the three values (ndigit,n,naux) specified by
c        the user. The user should either select a different combuination,
c        or rerun the subroutine oneseg (currently, a part of the file
c        test35.f) in order to add these data to the file IR: Do you think
c        the subroutine oneseg is a conversation piece, or what?
c  xs - the nodes of the "universal" n-point discretization of the
c        interval [-1,1]
c  xsaux - the nodes of the "universal" naux-point discretization of the
c        interval [-1,1]
c  whts - the weights of the quadrature formula with the nodes xs
c  whtsaux - the weights of the quadrature formula with the nodes xsaux
c 
c 
c       . . . find the header of the record corresponding to the
c             user-specified ndigits and n
c 
        ier=0
        rewind (ir)
c 
        do 2000 i=1,10000000
c 
        ncard=i
c 
c       read the next card
c 
        read(ir,1400,end=2200) card
 1400 format(70a1)
c 
c       check if this is the header
c 
        ifheader=1
        do 1600 j=1,25
c 
        if(card(j) .ne. template(j)) ifheader=0
 1600 continue
c 
        if(ifheader .eq. 0) goto 2000
c 
c       this is the header. Check if this is the header
c       of the record we want
c 
        write(image,1400) card
c 
cccc 1800 format(29x,i2,6x,i4,12x,i4)
 1800 format(29x,i2,4x,i4,13x,i4)
c 
        read(image,1800) ndigits7,n7,naux7
c 
cccc        call prinf('ndigits7=*',ndigits7,1)
cccc        call prinf('n7=*',n7,1)
cccc        call prinf('naux7=*',naux7,1)
        ifwant=1
        if(ndigits7 .ne. ndigits) ifwant=0
        if(n7 .ne. n) ifwant=0
        if(naux7 .ne. naux) ifwant=0
c 
        if(ifwant .eq. 1) goto 3000
c 
c       this was the header of a wrong record. skip most
c       of this record without reading
c 
        do 1900 jj=1,n7*naux7*2/3
  
        read(ir,1800)
 1900 continue
c 
 2000 continue
        ier=8
c 
 2200 continue
c 
        ier=4
        return
  
 3000 continue
c 
c       the desired header has been found; read data
c 
 2400 format(1a1)
 2600 format(1x,e21.15,1x,e21.15,1x,e21.15)
c 
        read(ir,2400)
        read(ir,2400)
c 
        read(ir,2600) (xs(i),i=1,n)
c 
c 
        read(ir,2400)
        read(ir,2400)
c 
        read(ir,2600) (xsaux(i),i=1,naux)
c 
c 
        read(ir,2400)
        read(ir,2400)
c 
        read(ir,2600) (whts(i),i=1,n)
c 
c 
        read(ir,2400)
        read(ir,2400)
c 
        read(ir,2600) (whtsaux(i),i=1,naux)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine onematr(n,naux,corrs,ainterp,
     1      xs,xsaux,totmat,auxtosup,interact)
        implicit real *8 (a-h,o-z)
c 
        save
        real *8 ainterp(n,naux),corrs(naux,n),auxtosup(n,naux),
     1          totmat(n,n),xs(1),xsaux(1)
c 
        external interact
c 
c       construct the matrix of values of potentials created
c       at the support nodes by charges at the auxiliary nodes
c 
        amatmax=0
  
        do 2400 i=1,n
        do 2200 j=1,naux
c 
        call interact(j,i,xsaux(j),xs(i),cd)
c 
        auxtosup(i,j)=corrs(j,i)*cd
c 
        if (abs(auxtosup(i,j)) .gt. amatmax)
     1      amatmax=abs(auxtosup(i,j))
 2200 continue
 2400 continue
c 
c        construct the matrix connecting the values of the charge
c        at the support nodes with the values of the potential at the
c        same nodes
c 
        call unimatmu(auxtosup,ainterp,n,naux,n,totmat)
  
        return
        end
  
c 
c 
c 
c 
c 
        subroutine unimatmu(a,b,n,m,k,c)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),b(m,k),c(n,k)
c 
        do 1600 i=1,n
        do 1400 j=1,k
        d=0
        do 1200 ii=1,m
c 
        d=d+a(i,ii)*b(ii,j)
 1200 continue
c 
        c(i,j)=d
 1400 continue
 1600 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine unicoret(ier,ir,ndigits,n,
     1      naux,xs,xsaux,whts,whtsaux,corrs,ainterp)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),xsaux(1),whts(1),corrs(1),
     1      ainterp(1),whtsaux(1)
c 
        character *1 card(70),template(1)
        character *100 image
        character * 26 templ20
c 
        equivalence(templ20,template(1))
c 
        data templ20/'      In this run, ndigits'/
c 
c       This subroutine reads from the file ir the all of the data
c       needed by the subroutine unisegm (see), and possibly by the
c       user dealing with the subroutine unisegm. The data are the
c       nodes and weights of the "universal" quadratures (n-point
c       and naux-point) on the interval [-1,1], the matrix corrs of
c       "universal" quadrature weights with nodes xsaux for the
c       evaluation of "universal" integrals at the nodes xs, and the
c       "universal" interpolation matrix connecting the values of
c       fairly general functions at the nodes xs to the values of
c       such functions at the nodes xsaux.
c 
c                      Input parameters:
c 
c  ir - the FORTRAN file from which the data are to be read.
c        Hopefully, these data have been stored there by previous
c        calls to the the subroutine oneseg
c  ndigits - the number of digits used in the design of the discretization
c        to be read from disk
c  n - the number of nodes in the discretization of the segment to be
c        constructed
c  naux - the number of auxiliary nodes to be used in the construction
c        of the discretization
c 
c                              Output parameters:
c 
c  ier - error return code.
c     ier=0 means successful execution;
c     ier=4 means that the file number IR does not contain data
c        corresponding to the three values (ndigit,n,naux) specified by
c        the user. The user should either select a different combuination,
c        or rerun the subroutine oneseg (currently, a part of the file
c        test35.f) in order to add these data to the file IR: Do you think
c        the subroutine oneseg is a conversation piece, or what?
c  xs - the nodes of the "universal" n-point discretization of the
c        interval [-1,1]
c  xsaux - the nodes of the "universal" naux-point discretization of the
c        interval [-1,1]
c  whts - the weights of the quadrature formula with the nodes xs
c  whtsaux - the weights of the quadrature formula with the nodes xsaux
c  corrs - the array dimensionsd (naux,n) of "universal" quadrature
c        weights with nodes xsaux for the evaluation of "universal"
c        integrals at the nodes xs
c  ainterp - the array dimensioned (n,naux) of the "universal"
c        interpolation coefficients connecting the values of
c        fairly general functions at the nodes xs to the values of
c        such functions at the nodes xsaux.
c 
c       . . . find the header of the record corresponding to the
c             user-specified ndigits and n
c 
        ier=0
        rewind (ir)
c 
        do 2000 i=1,10000000
c 
        ncard=i
c 
c       read the next card
c 
        read(ir,1400,end=2200) card
 1400 format(70a1)
c 
c       check if this is the header
c 
        ifheader=1
        do 1600 j=1,25
c 
        if(card(j) .ne. template(j)) ifheader=0
 1600 continue
c 
        if(ifheader .eq. 0) goto 2000
c 
c       this is the header. Check if this is the header
c       of the record we want
c 
        write(image,1400) card
c 
cccc 1800 format(29x,i2,6x,i4,12x,i4)
cccc 1800 format(29x,i2,5x,i5,12x,i4)
 1800 format(29x,i2,4x,i4,13x,i4)
c 
        read(image,1800) ndigits7,n7,naux7
c 
cccc        call prinf('ndigits7=*',ndigits7,1)
cccc        call prinf('n7=*',n7,1)
cccc        call prinf('naux7=*',naux7,1)
        ifwant=1
        if(ndigits7 .ne. ndigits) ifwant=0
        if(n7 .ne. n) ifwant=0
        if(naux7 .ne. naux) ifwant=0
c 
        if(ifwant .eq. 1) goto 3000
c 
c       this was the header of a wrong record. skip most
c       of this record without reading
c 
        do 1900 jj=1,n7*naux7*2/3
        read(ir,1800)
 1900 continue
  
  
  
c 
 2000 continue
        ier=8
  
 2200 continue
c 
        ier=4
        return
  
 3000 continue
c 
c       the desired header has been found; read data
c 
 2400 format(1a1)
 2600 format(1x,e21.15,1x,e21.15,1x,e21.15)
c 
        read(ir,2400)
        read(ir,2400)
c 
        read(ir,2600) (xs(i),i=1,n)
c 
c 
        read(ir,2400)
        read(ir,2400)
c 
        read(ir,2600) (xsaux(i),i=1,naux)
c 
c 
        read(ir,2400)
        read(ir,2400)
c 
        read(ir,2600) (whts(i),i=1,n)
c 
c 
        read(ir,2400)
        read(ir,2400)
c 
        read(ir,2600) (whtsaux(i),i=1,naux)
c 
c 
        read(ir,2400)
        read(ir,2400)
c 
        read(ir,2600) (ainterp(i),i=1,n*naux)
c 
c 
        read(ir,2400)
        read(ir,2400)
c 
        read(ir,2600) (corrs(i),i=1,n*naux)
  
        return
        end
  
  
  
  
  
  
