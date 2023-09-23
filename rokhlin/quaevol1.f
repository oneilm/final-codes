        implicit real *8 (a-h,o-z)
        dimension w(9 000 000)
c 
        external funuser,funwht,funuser2
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
C 
        PRINT *, 'ENTER c'
        READ *,c
        CALL PRIN2('c=*',c,1 )
c 
c        initialize the function computer
c 
        epsbell=1.0d-12
        call funwhtini(epsbell)
  
        epssmall=1.0d-14
        call funuserini(c,nvects,epssmall)
c 
        nfuns=nvects-2
c 
        ax=-1
cccc        ax=0
        bx=1
c 
        k=24
        lw=9 000 000
c 
        call prinf('originally,  lw=*',lw,1)
c 
        call prinf('before quaevol2, lw=*',lw,1)
c 
  
        epssvd=1.0d-10
        epsquadr=1.0d-10
        iw=71
        nsings=nfuns
        npols=10
  
  
        call quaevol2(ier,ax,bx,k,funwht,funuser,funuser2,
     1      nsings,npols,epssvd,epsquadr,w,lw,lused,iw,ncols,
     2      nfintegr)
  
  
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine funwht(x,weight)
        implicit real *8 (a-h,o-z)
c 
c 
        save
        dimension wbell(100 000),rkhibell(10000)
c 
        call prol0eva(x,wbell,psi0,derpsi0)
c 
        weight=psi0
  
  
        weight=weight*(1+1.0d-10+cos(20*x))
cccc        weight=weight*(1+1.0d-10+sin(x))
  
        return
c 
c 
c 
c 
        entry funwhtini(eps)
c 
c        construct the bell corresponding to this eps
c 
        call provcget(ier,eps,cbell)
  
        call prin2('in funuserini, cbell=*',cbell,1)
c 
        lenwbell=100 000
        call prol0ini(ier,cbell,wbell,rlam20,rkhibell,
     1      lenwbell,keep2,ltot2)
c 
        call prinf('in funuserini after prol0ini, ier=*',ier,1)
  
        xend=1
        call prol0eva(xend,wbell,psi0,derpsi0)
c 
        call prin2('in funuserini, psi0(xend)=*',psi0,1)
c 
        return
        end
c 
c 
c 
c 
c 
  
        subroutine funuser(x,i,f)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension w(4000 000),rkhi(10 000),rlams2(10 000)
  
        complex *16 rlams(10 000)
  
c 
        dimension wbell(100 000),rkhibell(10000)
  
cccc        call prin2('in funuser, x=*',x,1)
  
        call prolevv(i-1,x,w,val)
c 
  
        f=val * sqrt(rlams2(i))
  
  
cccc           call prin2('in funuser, f=*',f,1)
  
c 
        a=30
        return
c 
c 
c 
c 
        entry funuserini(c7,nvects7,eps)
c 
c        construct the PSWFs corresponding to this c
c 
        lenw=4000 000
        c=c7
        call prolcrea(ier,c,w,lenw,nvects,nhigh,
     1    rkhi,rlams,rlams2,keep,lused)
  
  
        call prin2('in funuserini after prolcrea, rlams=*',
     1      rlams,nvects*2)
  
        nvects7=nvects
c 
  
c 
c        construct the bell corresponding to this eps
c 
        call provcget(ier,eps,cbell)
  
        call prin2('in funuserini, cbell=*',cbell,1)
  
  
  
        lenwbell=100 000
        call prol0ini(ier,cbell,wbell,rlam20,rkhibell,
     1      lenwbell,keep2,ltot2)
  
c 
        call prinf('in funuserini after prol0ini, ier=*',ier,1)
  
        xend=1
        call prol0eva(xend,wbell,psi0,derpsi0)
  
        call prin2('in funuserini, psi0(xend)=*',psi0,1)
  
  
  
        return
        end
  
c 
c 
c 
c 
c 
        subroutine funuser2(x,i,f)
        implicit real *8 (a-h,o-z)
c 
        save
        call legepol(x,i-1,f,der)
c 
        return
        end
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code and the beginning of the
c        quadrature code proper
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine quaevol2(ier,a,b,k,funwht,funuser1,funuser2,
     1      nsings,npols,epssvd,epsquadr,www,lw,lused,iw,ncols,
     2      nfintegr)
        implicit real *8 (a-h,o-z)
        save
        dimension www(1)
c 
        external quaevosu,funwht,funuser,funuser2
c 
ccccc        external quaevosu,funwht,funuser
c 
c       This subroutine constructs quadratures for a class of functions.
c       The subroutine starts with constructing an SVD of a set of functions,
c       and continues by constructing (or rather, attempting to construct)
c       a Generalized Gaussian quadrature for the obtained singular
c       vectors. In fact, the subroutine builds a Chebychev-type quadrature
c       (i.e. an n-point quadrature for n functions), which is trivial, and
c       attempts to reduce the number of nodes, one node at a time. If the
c       subroutine failed to obtain a Generalized Gaussian Quadrature, it
c       still has (usually) achieved something; in practice, the subroutine
c       tends to get within one or two nodes from a Gaussian quadrature.
c 
c       This subroutine assumes that the user-specified functions to be
c       integrated are of the form
c 
c       f(x)=omega(x) * phi_i(x) * psi_j(x),
c 
c       with omega a (non-negative) weight function, and phi_1, phi_2,...
c       ...phi_{nsings}, and  psi_1, psi_2,...,psi{npols} are two families
c       of functions. It is assumed that the rank of the family {phi_i}
c       is considerably lower than nsings (their number). Otherwise, the
c       user is advised to use the subroutine quaeval1 (see)
c 
c   IMPORTANT NOTE: This subroutine is a someahat specialized version of the
c       subroutine quaevol1. The difference is that this subroutine uses
c       the subroutine allsvcmb to compress the input functions, while
c       quaevol1 uses allsvbld. In other words, in most cases when the user
c       wants to use this subroutine, quaevol1 should be used instead!!!
c 
c                            Input parameters:
c 
c  a,b - the ends of the interval on which the functions are to be
c       integrated
c  k - the number of Gaussian nodes to be used on each subinterval
c  funwht - the weight function, specified by a subroutine (user-supplied)
c       The calling sequence of funwht has to be
c 
c                       funwht(x,weight),
c 
c       with x the argument (input) and weight the value of the weight
c       function (output). Please note that weight must be non-negative
c       for all x \in [a,b]
c  funuser - the user-supplied subroutine, defining the functions to
c       be integrated. The calling sequance of funuser has to be
c 
c                       funuser(x,i,f),
c 
c       with x the argument (input), i the sequence number of the
c       function to be evaluated (input), and f the value of the
c       function (output). Please note that i should be between 1
c       and nfuns (see below)
c  funuser2 - the the other user-supplied subroutine, defining the
c       functions to be integrated. The calling sequance of funuser2
c       has to be
c 
c                       funuser2(x,i,f),
c 
c       with x the argument (input), i the sequence number of the
c       function to be evaluated (input), and f the value of the function
c       (output). Please note that i should be between 1 and nfuns (see below)
c  nfuns - the number of functions to be integrated
c  epssvd - the accuracy the SVD to be performed in the construction of the
c       "basis" functions for which the quadratures are to be designed
c  epsquadr - the accuracy to which the subroutine will attempt to design
c       the quadratures.
C 
C        PLEASE NOTE THAT THE PRECISE RELATIONSHIP BETWEEN EPSSVD, EPSQUADR
C        IS TOO COMPLICATED TO BE DESCRIBED HERE IN DETAIL, BUT THE RATIO
C        EPSQUADR/EPSSVD SHOULD BE AT LEAST 100; THE GREATER THIS RATIO,
C        THE MORE STABLY THE ALGORITHM WILL PERFORM.
c 
c  lw - the length of the user-supplied work array w, in real *8 elements.
c       Should be large
c  iw - the FORTRAN unit number on which the obtained quadratures will be
c       stored
c 
c                             Output parameters:
c 
c  ier - error return code.
c     ier=0 means that quadratures were found for all numbers of
c             nodes up to nfintegr/2 (see below for the definition
c             of nfintegr)
c     ier=2 means that the subroutine failed to get the smallest
c             possible number of nodes (which is equal to nfintegr/2).
c             Otherwise, the execution was successful
c     ier=1024 means that during one of calls to adapgaus, the latter
c             failed to obtain the required accuracy. Normally, this
c             means that the requested accuracy eps (user-specified) is
c             WAY too small. This is a fatal error.
c     ier=64 means that the amount of memory was supplied in the work
c             array w is insufficient (also, see parameter lw above).
c             This is a fatal error.
c 
c  lused - the amount of memory in array w (in real *8 words) actually
c       utilized by this subroutine
c  ncols - the number of singular values obtained by allsvbld inside this
c        subroutine while obtaining the SVD of the user-supplied set of
c        functions, to precision epssvd
c  nfintegr - the number of obtained singular values that are greater
c        than epsintegr (ralative to the first singular value); this is
c        the number of singular vectors that are going to be integrated
c        by the quadratures stored by this subroutine on the FORTRAN
c        file iw
c 
c        . . . allocate memory for the construction of basis functions
c 
        nfuns=nsings*npols
  
        irlams=1
        lrlams=nsings+npols+100 000
c 
        irints=irlams+lrlams
        lrints=nsings+npols+100 000
c 
        iwww=irints+lrints
        lwww=lw-iwww-10
c 
c        construct the basis functions
c 
        ier=0
c 
        call prinf('before allsvcmb, lwww=*',lwww,1)
c 
        igraph=0
c 
        call allsvcmb(jer,a,b,k,quaevosu,funwht,funuser,funuser2,
     1      nsings,npols,epssvd,www(irlams),lrlams,www(iwww),lwww,
     2      ncols,nn,www(irints),lused0,keep,igraph)
c 
        call prinf('in quaevol2 after allsvbld, lused0=*',lused0,1)
c 
        lused=lused0+iwww
        call prinf('and lused=*',lused,1)
  
  
        call prinf('in quaevol2 after allsvbld, keep=*',keep,1)
        call prinf('in quaevol2 after allsvbld, jer=*',jer,1)
  
        if(jer .ne. 0) then
            ier=64
            return
        endif
c 
        call prinf('in quaevol2 after allsvbld, nn=*',nn,1)
        call prinf('in quaevol2 after allsvbld, k=*',k,1)
        call prinf('in quaevol2 after allsvbld, ncols=*',ncols,1)
c 
c       perform garbage collection and allocate memory for the
c       construction of quadratures proper
c 
        n=k*nn
c 
        irints2=irlams+ncols+1
        call arrmove(www(irints),www(irints2),ncols+1)
c 
        irints=irints2
        lrints=ncols+1
c 
        iwww2=irints+lrints
        call arrmove(www(iwww),www(iwww2),keep)
c 
        iwww=iwww2
        lwww=keep
c 
        iroots=iwww+lwww
        lroots=ncols+2
c 
        iweights=iroots+lroots
        lweights=ncols+2
c 
        iamatr=iweights+lweights
        lamatr=ncols*n+10
c 
        irnorms=iamatr+lamatr
        lrnorms=n+2
c 
        iipivots=irnorms+lrnorms
        lipivots=n+2
c 
        ltot1=iipivots+lipivots
c 
        call prinf('ltot1=*',ltot1,1)
c 
c 
        ixs2=iweights+lweights
        lxs2=ncols+2
c 
        iws2=ixs2+lxs2
        lws2=ncols+2
c 
        iips=iws2+lws2
        lips=ncols+2
c 
        iwwww=iips+lips
        lenwwww=lw-iwwww
  
        call prinf('and lenwwww=*',lenwwww,1)
c 
        ix=www(iwww+1)+iwww-1
c 
c       construct the quadratures
c 
        call quaevol0(ier,a,b,k,funwht,funuser,
     1      nfuns,epsquadr,www,lw,lused2,iw,www(irlams),
     2      www(irints),www(iwww),ncols,www(ix),nn,
     3      www(iamatr),nfintegr,
     4      www(iroots),www(iweights),www(irnorms),www(iipivots),
     5      www(ixs2),www(iws2),www(iips),www(iwwww),lenwwww)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine quaevol1(ier,a,b,k,funwht,funuser,
     1      nfuns,epssvd,epsquadr,www,lw,lused,iw,ncols,nfintegr)
        implicit real *8 (a-h,o-z)
        save
        dimension www(1)
c 
        external quaevosu,funwht,funuser
c 
c       This subroutine constructs quadratures for a class of functions.
c       The subroutine starts with constructing an SVD of a set of functions,
c       and continues by constructing (or rather, attempting to construct)
c       a Generalized Gaussian quadrature for the obtained singular
c       vectors. In fact, the subroutine builds a Chebychev-type quadrature
c       (i.e. an n-point quadrature for n functions), which is trivial, and
c       attempts to reduce the number of nodes, one node at a time. If the
c       subroutine failed to obtain a Generalized Gaussian Quadrature, it
c       still has (usually) achieved something; in practice, the subroutine
c       tends to get within one or two nodes from a Gaussian quadrature.
c 
c                            Input parameters:
c 
c  a,b - the ends of the interval on which the functions are to be
c       integrated
c  k - the number of Gaussian nodes to be used on each subinterval
c  funwht - the weight function, specified by a subroutine (user-supplied)
c       The calling sequence of funwht has to be
c 
c                       funwht(x,weight),
c 
c       with x the argument (input) and weight the value of the weight
c       function (output). Please note that weight must be non-negative
c       for all x \in [a,b]
c  funuser - the sssser-supplied subroutine, defining the functions to
c       be integrated. The calling sequance of funuser has to be
c 
c                       funuser(x,i,f),
c 
c       with x the argument (input), i the sequence number of the
c       function to be evaluated (input), and f the value of the function
c       (output). Please note that i should be between 1 and nfuns (see below)
c  nfuns - the number of functions to be integrated
c  epssvd - the accuracy the SVD to be performed in the construction of the
c       "basis" functions for which the quadratures are to be designed
c  epsquadr - the accuracy to which the subroutine will attempt to design
c       the quadratures.
C 
C        PLEASE NOTE THAT THE PRECISE RELATIONSHIP BETWEEN EPSSVD, EPSQUADR
C        IS TOO COMPLICATED TO BE DESCRIBED HERE IN DETAIL, BUT THE RATIO
C        EPSQUADR/EPSSVD SHOULD BE AT LEAST 100; THE GREATER THIS RATIO,
C        THE MORE STABLY THE ALGORITHM WILL PERFORM.
c 
c  lw - the length of the user-supplied work array w, in real *8 elements.
c       Should be large
c  iw - the FORTRAN unit number on which the obtained quadratures will be
c       stored
c 
c                             Output parameters:
c 
c  ier - error return code.
c     ier=0 means that quadratures were found for all numbers of
c             nodes up to nfintegr/2 (see below for the definition
c             of nfintegr)
c     ier=2 means that the subroutine failed to get the smallest
c             possible number of nodes (which is equal to nfintegr/2).
c             Otherwise, the execution was successful
c     ier=1024 means that during one of calls to adapgaus, the latter
c             failed to obtain the required accuracy. Normally, this
c             means that the requested accuracy eps (user-specified) is
c             WAY too small. This is a fatal error.
c     ier=64 means that the amount of memory was supplied in the work
c             array w is insufficient (also, see parameter lw above).
c             This is a fatal error.
c 
c  lused - the amount of memory in array w (in real *8 words) actually
c       utilized by this subroutine
c  ncols - the number of singular values obtained by allsvbld inside this
c        subroutine while obtaining the SVD of the user-supplied set of
c        functions, to precision epssvd
c  nfintegr - the number of obtained singular values that are greater
c        than epsintegr (ralative to the first singular value); this is
c        the number of singular vectors that are going to be integrated
c        by the quadratures stored by this subroutine on the FORTRAN
c        file iw
c 
c        . . . allocate memory for the construction of basis functions
c 
        irlams=1
        lrlams=nfuns+120
c 
        irints=irlams+lrlams
        lrints=nfuns+120
c 
        iwww=irints+lrints
        lwww=lw-iwww-10
c 
c        construct the basis functions
c 
        ier=0
c 
        call prinf('before allsvbld, lw=*',lw,1)
c 
        igraph=0
c 
        call allsvbld(jer,a,b,k,quaevosu,funwht,funuser,
     1      nfuns,epssvd,www(irlams),lrlams,www(iwww),lwww,
     2      ncols,nn,www(irints),lused0,keep,igraph)
c 
        call prinf('in quaevol2 after allsvbld, lused0=*',lused0,1)
c 
        lused=lused0+iwww
        call prinf('and lused=*',lused,1)
  
  
        call prinf('in quaevol2 after allsvbld, keep=*',keep,1)
        call prinf('in quaevol2 after allsvbld, jer=*',jer,1)
  
        if(jer .ne. 0) then
            ier=64
            return
        endif
c 
        call prinf('in quaevol2 after allsvbld, nn=*',nn,1)
        call prinf('in quaevol2 after allsvbld, k=*',k,1)
        call prinf('in quaevol2 after allsvbld, ncols=*',ncols,1)
c 
c       perform garbage collection and allocate memory for the
c       construction of quadratures proper
c 
        n=k*nn
c 
        irints2=irlams+ncols+1
        call arrmove(www(irints),www(irints2),ncols+1)
c 
        irints=irints2
        lrints=ncols+1
c 
        iwww2=irints+lrints
        call arrmove(www(iwww),www(iwww2),keep)
c 
        iwww=iwww2
        lwww=keep
c 
        iroots=iwww+lwww
        lroots=ncols+2
c 
        iweights=iroots+lroots
        lweights=ncols+2
c 
        iamatr=iweights+lweights
        lamatr=ncols*n+10
c 
        irnorms=iamatr+lamatr
        lrnorms=n+2
c 
        iipivots=irnorms+lrnorms
        lipivots=n+2
c 
        ltot1=iipivots+lipivots
c 
        call prinf('ltot1=*',ltot1,1)
  
c 
c 
        ixs2=iweights+lweights
        lxs2=ncols+2
c 
        iws2=ixs2+lxs2
        lws2=ncols+2
c 
        iips=iws2+lws2
        lips=ncols+2
c 
        iwwww=iips+lips
        lenwwww=lw-iwwww
  
        call prinf('and lenwwww=*',lenwwww,1)
c 
        ix=www(iwww+1)+iwww-1
c 
c       construct the quadratures
c 
        call quaevol0(ier,a,b,k,funwht,funuser,
     1      nfuns,epsquadr,www,lw,lused2,iw,www(irlams),
     2      www(irints),www(iwww),ncols,www(ix),nn,
     3      www(iamatr),nfintegr,
     4      www(iroots),www(iweights),www(irnorms),www(iipivots),
     5      www(ixs2),www(iws2),www(iips),www(iwwww),lenwwww)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine quaevol0(ier,a,b,k,funwht,funuser,
     1      nfuns,eps,www,lw,lused2,iw,rlams,rints,coefs,
     2      ncols,x,nn,amatr,nfuns7,
     3      roots,weights,rnorms,ipivots,xs2,ws2,ips,wwww,lwwww)
c 
        implicit real *8 (a-h,o-z)
        save
        dimension coefs(1),x(1),amatr(1),www(1),rlams(1),
     1      rints(1),roots(1),weights(1),rnorms(1),ipivots(1),
     2      ips(1),xs2(1),ws2(1),wwww(1)
c 
cccc        external quaevosu,quaevofu,funwht,funuser
        external quaevofu,funwht,funuser
c 
c       evaluate the correct integrals of the functions
c 
        mmm=20
        do 1600 i=1,ncols
c 
        eps333=eps/100/rlams(i)*rlams(1)
  
        call prinf('i=*',i,1)
        call prin2('and eps333=*',eps333,1)
c 
        iiii= quaevofi(i)
  
        call adapgaus(jer,a,b,quaevofu,funwht,coefs,mmm,eps333,
     1      rint,maxrec,numint)
c 
        if(jer .ne. 0) then
              call prinf(
     1   'bombing from quaevol0 with error code from adapgaus; jer=*',
     2            jer,1)
            ier=1024
            return
        endif
c 
        rints(i)=rint
 1600 continue
c 
        call prin2('and after adapgaus, rints=*',
     1      rints,ncols)
c 
c       determine the number of functions for which the quadrature
c       is to be constructed
c 
        do 2200 i=1,ncols-1
c 
        nfuns7=i
        if(rlams(i) .lt. eps) goto 2400
 2200 continue
 2400 continue
c 
        i=nfuns7/2
        nfuns7=i*2
cccc        ix=www(iwww+1)+iwww-1
c 
        n=nn*k
c 
c       find the initial (m-point) quadrature using the pivoted
c       Gram-Schmidt
c 
cccc        call prin2('before quaeiniq, x=*',x,n)
  
        call prinf('before quaeiniq, nfuns7=*',nfuns7,1)
c 
        call quaeiniq(amatr,n,nfuns7,coefs,x,
     1      rints,roots,weights,rnorms,ipivots)
c 
        call prin2('after quaeiniq, roots=*',roots,nfuns7)
        call prin2('after quaeiniq, weights=*',weights,nfuns7)
c 
        call quaebubb(roots,weights,nfuns7)
c 
c       write the header on the disk
c 
        call quaestin(iw,ncols,a,b,rlams)
c 
        do 2600 j=1,nfuns7
c 
        call funwht(roots(j),ddd)
c 
        ws2(j)=weights(j)*sqrt(ddd)
 2600 continue
c 
        call quaestor(iw,roots,ws2,nfuns7)
c 
c       find all quadratures for this collection of nfuns7 functions
c 
        call prin2('before quaepnts, roots=*',roots,nfuns7)
        call prin2('before quaepnts, weights=*',weights,nfuns7)
c 
        call quaepnts(ier,nfuns7,rints,roots,weights,nfuns7,
     1      coefs,eps*1000,ips,xs2,ws2,wwww,lwwww,iw,a,b)
c 
        return
        end
  
  
c 
c 
c 
c 
        subroutine quaepnts(ier,nfuns,rints,xs,ws,npts,
     1      coefs,epsout,ips,xs2,ws2,w,lw,iw,a,b)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),ws(1),ips(1),
     1      xs2(1),ws2(1),rints(1),w(1)
c 
c       Starting with an npts-node quadrature formula integrating
c       exactly nfuns functions, this subroutine successively
c       attempts to omit one node after another, and to deform the
c       remaining nodes and weights into quadrature formulae
c       integrating exactly the same nfuns functions. The process
c       is terminated when it fails.
c 
c                      Input parameters:
c 
c  nfuns - the number of functions that the input quadrature formula
c       integrates exactly (and, hopefully, that the output formula
c       will also integrate exactly)
c  rints - the integrals of the nfuns functions to be integrated
c  xs,ws - the nodes and weights of the n-point quadrature formula
c  npts - the number of nodes in the input quadrature formula
c  coefs,ab,nn,k - the parameters produced by the soubroutine allsvbld
c       (or allsvcmb) and used by this subroutine to evaluate the
c       functions to be integrated
c  epsout - the value of discrepancy at wich the subroutine will declare
c       victory. Should be somewhat loose (the subroutine will make two
c       extra steps of Newton after this accuracy has been reached).
c  iw - the FORTRAN unit number on which the obtained quadratures will
c       be stored
c 
c                      Output parameters:
c 
c                      Work array:
c 
c  ips,xs2,ws2,w - TO BE DOCUMENTED
c 
c 
c       construct the array of pointers to the nodes to be thrown out
c 
        call prin2('in quaepnts, a=*',a,1)
        call prin2('in quaepnts, b=*',b,1)
  
        ier=0
        n=npts
c 
        i00=n/2
        do 1200 i=1,n
c 
        ii=i+i00
        call quaewrap(n,ii,j)
        ips(i)=j
 1200 continue
c 
        call prinf('ips as created before the loop*',ips,n)
c 
        iips=n/2
        imin=ips(iips)
c 
        if( (imin .gt. n) .or. (imin .le. 0) ) then
            call prinf('bombing from quaepnts with imin=*',imin,1)
c 
            ier=2
            return
        endif
c 
c       conduct the main loop
c 
        ifodd=1
        do 2800 i=1,100000
c 
        ifodd=-ifodd
c 
c       if this is an odd loop - check if one of weights is negative
c 
cccc        goto 1300
c 
        if(ifodd .eq. -1) goto 1300
c 
cccc        do 1250 j=1,n-1
        do 1250 j=1,n
c 
        call quaewrap(n,j,jjj)
        if(ws(jjj) .lt. 0) then
            imin=j
            iips=iips-1
            goto 1300
        endif
c 
 1250 continue
 1300 continue
c 
        call quaepnt1(jer,xs,ws,n,nfuns,rints,xs2,ws2,imin,
     1      coefs,epsout,w,lw,a,b)
c 
          if(jer .eq. 0) goto 1350
c 
c       the search failed. try throwing out the next node and starting
c       again
c 
          iips=iips+1
          imin=ips(iips)
c 
        if( iips .gt. n ) then
            call prinf('bombing from quaepnts with iips=*',iips,1)
c 
        ier=2
        return
        endif
c 
        goto 2800
c 
 1350 continue
c 
c       the search has been successful. Store the results on disk,
c       and go to the next (smaller) n
c 
        n=n-1
        iips=1
c 
        call mvecmove(xs2,xs,n)
        call mvecmove(ws2,ws,n)
c 
        do 1370 j=1,n
c 
        call funwht(xs2(j),ddd)
c 
        ws2(j)=ws2(j)*sqrt(ddd)
 1370 continue
c 
        call quaestor(iw,xs2,ws2,n)
c 
c       construct the new sequence of elements that we will try
c       throwing out
c 
        i00=imin-1
c 
        do 1400 j=1,n
c 
        ii=j+i00
        call quaewrap(n,ii,jj)
        ips(j)=jj
 1400 continue
c 
        j1=ips(1)
        j2=ips(2)
        j3=ips(3)
c 
        ips(1)=j2
        ips(2)=j3
        ips(3)=j1
c 
        call prinf('ips as created*',ips,n)
        call prinf('while n= *',n,1)
        call prinf('and imin= *',imin,1)
c 
        jj=nfuns/2
        if(jj .ge. n) goto 3200
c 
 2800 continue
c 
 3200 continue
c 
        call prinf('at the end, n=*',n,1)
c 
        call prin2('for this n, xs=*',xs,n)
        call prin2('for this n, ws=*',ws,n)
c 
        ier=0
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine quaewrap(n,i,j)
        implicit real *8 (a-h,o-z)
c 
        save
        j=i
        if( (i .ge. 1) .and. (i .le. n) ) return
c 
        if(i .gt. n) j=i-n
        if(i .lt. 1) j=n+i
        return
        end
c 
c 
c 
c 
c 
        subroutine quaepnt1(ier,xs,ws,n,nfuns,rints,xs2,ws2,imin,
     1      coefs,epsout,w,lw,a,b)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),ws(1),w(1),rints(1),xs2(1),ws2(1)
c 
c       Starting with an n-node quadrature formula integrating
c       exactly nfuns functions, this subroutine attempts to
c       omit the node number imin (input) and starting with the
c       resulting n-1 points, to construct an n-1-node quadrature
c       formula integrating exactly the same nfuns functions.
c 
c                      Input parameters:
c 
c  xs,ws - the nodes and weights of the n-point quadrature formula
c  n - the number of nodes in the input formula
c  nfuns - the number of functions that the input quadrature formula
c       integrates exactly (and, hopefully, that the output formula
c       will also integrate exactly)
c  rints - the integrals of the nfuns functions to be integrated
c  imin - the sequence number of the node to be omitted
c  coefs,ab,nn,k - the parameters produced by the soubroutine allsvbld
c       (or allsvcmb) and used by this subroutine to evaluate the
c       functions to be integrated
c  epsout - the value of discrepancy at wich the subroutine will declare
c       victory. Should be somewhat loose (the subroutine will make two
c       extra steps of Newton after this accuracy has been reached).
c 
c                      Output parameters:
c 
c  ier - error return code
c  xs2, ws2 - the nodes and weights of the n-1-node quadrature formula
c 
c                      Work array:
c 
c  w - must be sufficiently large
c 
c 
c       . . . copy the nodes and weights, omitting the one
c             selected at this step
c 
        ier=0
        j=0
        do 1400 i=1,n
c 
        if(i .eq. imin) goto 1400
c 
        j=j+1
        xs2(j)=xs(i)
        ws2(j)=ws(i)
 1400 continue
c 
c       find the n-1-point quadrature
c 
        npts=n-1
        kold=-7
c 
        ifout=0
c 
        isol=1
        lsol=npts*2+2
c 
        irhs=isol+lsol
        lrhs=npts*2+nfuns
c 
        ixs22=irhs+lrhs
        lxs22=npts*2+nfuns
c 
        iws22=ixs22+lxs22
        lws22=npts*2+nfuns
c 
        isums=iws22+lws22
        lsums=npts*2+nfuns
c 
        iww=isums+lsums
        lenww=lw-iww-2
  
        call prinf('before quaenewt, lenww=*',lenww,1)
  
        lneed=nfuns*npts*10+nfuns+npts*2 +100
  
        call prinf('before quaenewt, lneed=*',lneed,1)
c 
        if(lneed .gt. lenww) then
            ier=2048
            call prinf('in quaepnt1, memory needed is*',lneed,1)
            call prinf('in quaepnt1, memory available is*',lenww,1)
            stop
        endif
c 
        do 2400 i=1,120
  
        call prinf('i=*',i,1)
c 
        call quaenewt(jer,xs2,ws2,npts,nfuns,w(isums),rints,dnew,
     1      dold,coefs,kcontr,w(iww),
     2      w(isol),w(irhs),w(ixs22),w(iws22),a,b)
  
  
        call prinf('after quaenewt, jer=*',jer,1)
        call prinf('after quaenewt, kcontr=*',kcontr,1)
c 
        if( (kcontr .eq. 4) .and. (kold .eq. 4) ) jer=3
        kold=kcontr
c 
        if(jer .eq. 0) goto 2200
c 
        call prinf('in quaepnt1 after quaenewt, jer=*',jer,1)
        ier=4
        return
 2200 continue
  
  
        call prin2('in quaepnt1, epsout=*',epsout,1)
  
cccc          stop
c 
        if(ifout .eq. 3) goto 2600
        if(ifout .ge. 1) ifout=ifout+1
        if( (dnew .lt. epsout) .and. (ifout .eq. 0) ) ifout=1
c 
        if( (i .eq. 30) .and. (dold/dnew .lt. 1.01) ) goto 2500
        if( (i .eq. 60) .and. (dold/dnew .lt. 1.01) ) goto 2500
        if( (i .eq. 90) .and. (dold/dnew .lt. 1.01) ) goto 2500
 2400 continue
c 
 2500 continue
  
        ier=4
        return
c 
 2600 continue
c 
        call quaebubb(xs2,ws2,npts)
c 
        call prin2('and xs2=*',xs2,npts)
        call prin2('and ws2=*',ws2,npts)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine quaenewt(ier,xs,ws,npts,nfuns,sums,
     1      rints,dnew,dold,coefs,kcontr,w,
     2      sol,rhs,xs2,ws2,a,b)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),ws(1),sums(1)
        dimension w(1),sol(1),rints(1),rhs(1),xs2(1),ws2(1)
c 
c       This subroutine conducts a single Newton iteration for the
c       determination of the npts-point quadrature for nfuns
c       functions.
c 
c                          Input parameters:
c 
c  xs - starting nodes
c  ws - starting weights
c  npts - the number of nodes in the quadrature being sought
c  nfuns - the number of functions we want to integrate
c  rints - the integrals of the functions for which we are designing
c       the quadratures
c  coefs, ab - the parameters produced by a prior call to the
c       subroutine allsvbld (or allsvcmb); will be used by this
c       subroutine to evaluate the functions for which the quadrature
c       is being designed
c  kk - the parameter k as given to the subroutine allsvbld (see);
c       will be used by this subroutine to evaluate the functions
c       for which the quadrature
c 
c                          Output parameters:
c 
c  ier - error return code. Ier=0 means successful conclusion.
c       ..............
c  xs - updated nodes
c  ws - updated weights
c  dnew - the discrepancy at the new nodes
c  kcontr - the value of the control parameter. EXPLANATION: the
c       subroutine attempts to make a step of the Newton process. If
c       such a step reduces the discrepancy, the subroutine sets the
c       coefficient kcontr to 1, and exits. If the step fails to reduce
c       the discrepancy, the subroutine divides it by 4, sets kconts to
c       2, and exits. If after scaling by 1/4 the step fails to reduce
c       the discrepancy, the subroutine scales the step by 1/16, sets
c       kcontr to 3, and exits,...
c 
c                          Work arrays:
c 
c  sums, sol, rhs, xs2, ws2
c  w - must be at least nfuns*npts*10 + nfuns+npts*2+100
c 
c       . . . allocate memory
c 
  
        call prin2('in quaenewt, a=*',a,1)
        call prin2('in quaenewt, b=*',b,1)
  
        np2=npts*2
c 
        ia=1
        la=nfuns*np2+10
c 
        iu=ia+la
        lu=nfuns*np2+10
c 
        iw=iu+lu
        lw=nfuns*np2+10
c 
        itt=iw+lw
        ltt=nfuns*np2+10
c 
        iv=itt+ltt
        lv=nfuns*np2+10
c 
        irnorms=iv+lv
        lrnorms=nfuns+np2
  
  
        ltot2=irnorms+lrnorms
c 
c       construct the linearized problem to be solved
c 
        ier=0
        itype=1
        call quaematr(xs,ws,npts,nfuns,w(ia),sums,itype,
     1      coefs)
c 
        do 1200 i=1,nfuns
c 
        rhs(i)=sums(i)-rints(i)
 1200 continue
c 
        dold=0
        rnor22=0
        do 1400 i=1,nfuns
c 
        dold=dold+(sums(i)-rints(i))**2
        rnor22=rnor22+rints(i)**2
 1400 continue
c 
c       construct the Gram-schmidt decomposition of the matrix a
c 
        eps=1.0d-12
        call leastsq(w(ia),w(iu),w(iw),w(itt),
     1      nfuns,np2,ncols,w(irnorms),eps,w(iv) )
c 
c       construct the least squares solution of the Newton equation
c 
        call leastsq2(w(iu),w(iw),w(itt),nfuns,np2,ncols,rhs,sol,w(iv))
c 
c       perform step-length control
c 
        scale=1
        icount=0
        do 4000 k=1,50
c 
        kcontr=k
  
        call prinf('in quaenewt, k=*',k,1)
c 
c       . . . recompute the right-hand sides
c 
        do 2400 i=1,npts
c 
        xs2(i)=xs(i)-sol(i)           *scale
  
cccc        if(xs2(i) .gt. 1) xs2(i)=1 -1.0d-8
cccc        if(xs2(i) .lt. -1) xs2(i)=-1 +1.0d-8
c 
        if(xs2(i) .gt. b) xs2(i)=b *(1-1.0d-14)
        if(xs2(i) .lt. a) xs2(i)=a*(1+1.0d-14)
c 
        ws2(i)=ws(i)-sol(npts+i)      *scale
 2400 continue
c 
        itype=2
        call quaematr(xs2,ws2,npts,nfuns,w(ia),sums,itype,
     1      coefs )
c 
        dnew=0
        do 2600 i=1,nfuns
c 
        dnew=dnew+(sums(i)-rints(i))**2
 2600 continue
c 
        call prin2('and dold=*',dold,1)
        call prin2('and dnew=*',dnew,1)
        call prin2('and rnor22=*',rnor22,1)
  
        if( (dold .gt. dnew) .and. (k .eq. 1) ) goto 4200
c 
        if(dold .gt. dnew) icount=icount+1
        if(icount .ge. 1) goto 4200
        if(dnew .lt. rnor22*1.0d-22) icount=icount+1
        scale=scale/2
        scale=scale/2
c 
        if(k .gt. 5) ier=4
        if(k .gt. 5) return
c 
 4000 continue
c 
 4200 continue
c 
        call prinf('kcontr=*',kcontr,1)
        call prinf('and npts=*',npts,1)
        call prin2('and dnew=*',dnew,1)
        call prin2('and dold=*',dold,1)
c 
        if(kcontr .gt. 4) ier=4
c 
        call mvecmove(xs2,xs,npts)
        call mvecmove(ws2,ws,npts)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine quaematr(xs,ws,npts,nfuns,a,sums,itype,
     1      coefs)
        implicit real *8 (a-h,o-z)
        save
        dimension a(nfuns,npts*2),xs(1),ws(1),sums(1)
c 
c        construct the matrix of the Jacobian of the mapping
c        connection the points xs and the weights ws with the
c        integrals of the functions
c 
        if(itype. eq. 2) goto 2100
  
        do 2000 i=1,nfuns
c 
        d=0
        do 1400 j=1,npts
c 
        call nesteva(ier,coefs,xs(j),i,pol,der)
        d=d+ws(j)*pol
c 
        a(i,j)=der*ws(j)
        a(i,j+npts)=pol
 1400 continue
        sums(i)=d
 2000 continue
c 
        return
c 
 2100 continue
c 
c       calculate the approximations to the integrals
c 
        do 2400 i=1,nfuns
c 
        d=0
        do 2200 j=1,npts
c 
        call nestev2(ier,coefs,xs(j),i,pol)
c 
        d=d+ws(j)*pol
 2200 continue
c 
        sums(i)=d
 2400 continue
c 
       return
       end
c 
c 
c 
c 
c 
        subroutine quaeiniq(amatr,n,ncols,coefs,xs,
     1      rints,roots,weights,rnorms,ipivots)
        implicit real *8 (a-h,o-z)
        save
        dimension amatr(ncols,n),coefs(1),xs(1),
     1      rnorms(1),roots(1),weights(1),
     2      ipivots(1),rints(1)
c 
c        This subroutine constructs a CHEBYCHEV quadrature for
c        a collection of functions. The quadrature it constructs
c        (usually) requires twice as many nodes as the Gaussian
c        quadrature, but is much more easily and reliably
c        constructed. It is also a preliminary step in the
c        construction of Gaussian quadratures. Please note
c        that this subroutine uses as its input the autput of the
c        subroutine allsvbld (or allsvcmb); this subroutine has
c        no known uses as a stand-alone device. Also, please note
c        that the resulting weights are not guaranteed to be
c        positive, though they tend to be so.
c 
c                      Input parameters:
c 
c  n - the number of nodes in the discretization of the functions
c  ncols - the number of functions for which the quadrature is to
c        be designed
c  coefs - Legendre expansions of the functions to be integrated,
c        as produced by allsvbld or allsvcmb (see)
c  ab, nn,k,xs,rints - produced by allsvbld or allsvcmb (see)
c 
c                      Output parameters:
c 
c  roots - the ncols quadrature nodes
c  weights - the weights corresponding to the nodes roots
c 
c                      Work arrays:
c 
c  rnorms - must be at least n real *8 elements long
c  ipivots - must be at least n integer *4 elements long
c 
c       . . . construct the matrix of values of the functions for which
c             the quadratures are to be designed
c 
        do 1400 i=1,ncols
        do 1200 j=1,n
c 
        call nestev2(ier,coefs,xs(j),i,val)
c 
        amatr(i,j)=val
 1200 continue
 1400 continue
c 
c       conduct the pivoted Gram-Schmidt on amatr
c 
        eps=1.0d-30
        call quaeppiv(amatr,ncols,n,rnorms,eps,ncols,ipivots)
c 
c       sort the array ipivots
c 
        call allsbubb(ipivots,ncols)
c 
c       construct the nodes of the rudimentary quadrature
c 
        do 1600 i=1,ncols
c 
        roots(i)=xs(ipivots(i))
 1600 continue
c 
c       construct the corresponding weights
c 
        do 2400 i=1,ncols
        do 2200 j=1,ncols
c 
        call nestev2(ier,coefs,roots(j),i,val)
        amatr(i,j)=val
 2200 continue
c 
        weights(i)=rints(i)
 2400 continue
c 
        call qrsolv(amatr,ncols,weights,rcond)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine quaebubb(xs,ws,n)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),ws(1)
c 
        do 1400 i=1,n
        do 1200 j=1,n-1
        if(xs(j) .lt. xs(j+1)) goto 1200
c 
        d=xs(j)
        xs(j)=xs(j+1)
        xs(j+1)=d
c 
        d=ws(j)
        ws(j)=ws(j+1)
        ws(j+1)=d
 1200 continue
 1400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine quaeppiv(b,n,m,rnorms,eps,ncols,ipivots)
        implicit real *8 (a-h,o-z)
        save
        dimension b(n,m),rnorms(1),ipivots(1)
c 
c       This subroutine applies a pivoted double gram-schmidt
c       procedure to the matrix b. THIS SUBROUTINE IS DIFFERENT
C       FROM MOST OF ITS RELATIVES IN THAT IS TERMINATES THE
C       GRAM-SCHMIDT PROCESS AFTER NCOLS PIVOTS, WITH NCOLS THE
C       USER-SPECIFIED INTEGER PARAMETER. The number of the resulting
c       orthogonal vectors vectors is (obviously) equal to ncols.
c 
c                    input paramneters:
c 
c  b - the matrix to be gram-schmidt'ed. it is destroyed by this
c       subroutine
c  n,m - dimensionalities of the matrix b
c  eps - the machine accuracy
c  ncols - the number of pivots to use
c 
c                     output parameters:
c 
c  b - the matrix of gram-schmidt vectors of the matrix a. note
c        that on exit from this subroutine only the first
c        ncols columns of  b  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit b is dimensioned
c        b(n,ncols)
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c  ipivots - the sequence numbers of the pivots used, in the order in
c        which they were used
c 
        do 1100 i=1,m
        ipivots(i)=i
 1100 continue
c 
c       . . . prepare the array of values of norms
c             of columns
c 
        done=1
        dtot=0
        do 1400 i=1,m
        d=0
        do 1200 j=1,n
        d=d+b(j,i)**2
 1200 continue
        rnorms(i)=d
        dtot=dtot+d
 1400 continue
c 
c       . . . conduct gram-schmidt iterations
c 
        thresh=dtot*eps**2
        do 4000 i=1,ncols
c 
c       find the pivot
c 
        ipivot=i
        rn=rnorms(i)
c 
        do 2200 j=i+1,m
        if(rnorms(j) .le. rn) goto 2200
        rn=rnorms(j)
        ipivot=j
 2200 continue
 2400 continue
c 
c       put the column number ipivot in the i-th place
c 
        do 2600 j=1,n
        d=b(j,i)
        b(j,i)=b(j,ipivot)
        b(j,ipivot)=d
 2600 continue
c 
        iijj=ipivots(i)
        ipivots(i)=ipivots(ipivot)
        ipivots(ipivot)=iijj
c 
        d=rnorms(ipivot)
        rnorms(ipivot)=rnorms(i)
        rnorms(i)=d
c 
c       orthogonalize the i-th column to all preceding ones
c 
        if(i .eq. 1) goto 2790
        do 2780 j=1,i-1
c 
        call mvecscap(b(1,i),b(1,j),n,d)
c 
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*d
 2770 continue
 2780 continue
 2790 continue
c 
c       normalize the i-th column
c 
        call mvecscap(b(1,i),b(1,i),n,d)
c 
c       if the norm of the remaining part of the matrix
c       is sufficiently small - terminate the process
c 
        if(d .lt. thresh) return
c 
        d=done/dsqrt(d)
        do 2800 j=1,n
        b(j,i)=b(j,i)*d
 2800 continue
c 
c        orthogonalize everything else to it
c 
        do 3200 j=i+1,m
c 
        if(rnorms(j) .lt. thresh) goto 3200
c 
        call mvecscap(b(1,i),b(1,j),n,d)
c 
        rrn=0
        do 3000 l=1,n
        b(l,j)=b(l,j)-b(l,i)*d
        rrn=rrn+b(l,j)**2
 3000 continue
        rnorms(j)=dsqrt(rrn)
 3200 continue
 3400 continue
 4000 continue
c 
 4200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine quaestin(iw,ncols,a,b,rlams)
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
        entry quaestor(iw,xs,ws,ns)
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
	i0=mod(i,10)
        if(i0 .eq. 0) i0=11
        write(iw,2600) i0,xs(i*2-1),xs(i*2)
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
	i0=mod(i,10)
        if(i0 .eq. 0) i0=11
        write(iw,2600) i0,ws(i*2-1),ws(i*2)
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
        subroutine quaeretr(ier,ir,n,x,ws)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),ws(1)
        character *1 card(100),Nodes(5),aaa1(100),weights(7)
        character *100 aaa
c 
        equivalence(aaa1(1),aaa)
c 
        data Nodes/'N','o','d','e','s'/,
     1      weights/'W','e','i','g','h','t','s'/
c 
c       form the line to be found in trhe file
c 
 1020 format('c  Data for ',i3,' nodes')
c 
        write(aaa,1020) n
c 
c      keep reading the file till we get to the right line
c 
 1200 format(80a1)
c 
        ier=8
        rewind(ir)
        numrec=0
        do 1350 i=1,10000
c 
        read(ir,1200,end=3000) (card(j),j=1,30)
c 
c       determine if this is the correct line
c 
        ifcorr=1
c 
        do 1300 j=1,20
c 
        if(card(j) .ne. aaa1(j)) ifcorr=0
 1300 continue
c 
        if(ifcorr .eq. 1) ier=0
        if(ifcorr .eq. 1) goto 1360
c 
 1350 continue
c 
 1360 continue
c 
c       find the nodes
c 
        do 1500 i=1,1000
c 
        read(ir,1200,end=3000) (card(j),j=1,30)
c 
        ifline=1
        do 1400 j=1,5
c 
        if(card(j+9) .ne. Nodes(j)) ifline=0
 1400 continue
c 
        if(ifline .eq. 1) goto 1600
 1500 continue
c 
        ier=16
        return
c 
 1600 continue
c 
        read(ir,1200,end=3000) (card(j),j=1,30)
c 
c       read the nodes
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
c       find the weights
c 
        do 3500 i=1,1000
c 
        read(ir,1200,end=3000) (card(j),j=1,30)
c 
cccc        call prina('card as read*',card,20)
  
        ifline=1
        do 3400 j=1,5
c 
        if(card(j+9) .ne. weights(j)) ifline=0
 3400 continue
c 
        if(ifline .eq. 1) goto 3600
 3500 continue
c 
 3600 continue
c 
        read(ir,1200,end=3000) (card(j),j=1,30)
c 
c       read the nodes
c 
        do 4600 i=1,n/2
c 
 4400 format(8x,e22.16,1x,e22.16)
c 
        read(ir,4400) ws(2*i-1),ws(2*i)
 4600 continue
c 
        n2=n/2
        j=n-n2*2
        if(j .ne. 0) read(ir,2400) ws(n)
c 
 3000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine quaevosu(x,i,funwht,funuser,f)
        implicit real *8 (a-h,o-z)
c 
c       evaluate the function
c 
        save
        call funuser(x,i,f)
c 
c       evaluate the weight and scale the function by its square root
c 
        call funwht(x,weight)
c 
        f=f*sqrt(weight)
c 
        return
        end
c 
c 
c 
c 
c 
        function quaevofu(x,funwht,coefs)
        implicit real *8 (a-h,o-z)
        save
        dimension coefs(1)
c 
c       evaluate the integrand
c 
        call nestev2(ier,coefs,x,i,f)
c 
c       evaluate the weight and scale the function by its square root
c 
        call funwht(x,weight)
        quaevofu=f*sqrt(weight)
c 
        return
c 
c 
c 
c 
        entry quaevofi(i7)
        i=i7
        quaevofi=0
        return
        end
