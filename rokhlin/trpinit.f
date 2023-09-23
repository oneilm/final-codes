        implicit real *8 (a-h,o-z)
        dimension x0y0(20),vert1(3),vert2(3),
     1      vert3(3),rints(1000),rints2(1000),
     2      wplot(10 000),diff(100),www(50 000),work(50000)
c 
        external fun6
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER nord'
        READ *,nord
        CALL PRINf('nord=*',nord,1 )
c 
        x0=1.1
        y0=-0.7
        z0=0.01
cccc        z0=2.01
  
cccc        x0=0.3
        y0=-2
        y0=-1.5
        y0=-1.2
c 
        x0=0.7
        y0=-0.6
  
  
        x0=5
        y0=0.1
        z0=1.0d-7
  
c 
        x0y0(1)=x0
        x0y0(2)=y0
        x0y0(3)=z0
c 
c        design the computer of interals of orthogonal polynomials
c 
        call trpinit(nord,www,keep)
c 
        call prinf('after trpinit, keep=*',keep,1)
c 
c       to evaluate the INTEGRALS of the form
c 
c      \int P_{i,j} / sqrt( (x-x0)^2 + (y-y0)^2) dx dy
c 
c       over the standard triangle
c 
        call trpinteg(x0,y0,z0,www,rints,work)
c 
        n=(nord+1)*(nord+2)/2
c 
        call prin2('after trpinteg, rints=*',rints,n)
c 
c       evaluate the integral via the two-dimensional quadrature
c 
c 
c       . . . create the standard triangle
c 
        vert1(1)=0
        vert1(2)=sqrt(3.0d0)
c 
        vert2(1)=-1
        vert2(2)=0
c 
        vert3(1)=1
        vert3(2)=0
c 
        d=3
        d=1/sqrt(d)
        vert1(2)=vert1(2)-d
        vert2(2)=vert2(2)-d
        vert3(2)=vert3(2)-d
c 
c        . . . plot the whole picture
c 
        call trplopen(wplot)
c 
        call trplline(vert1(1),vert1(2),vert2(1),vert2(2),wplot)
        call trplline(vert1(1),vert1(2),vert3(1),vert3(2),wplot)
        call trplline(vert2(1),vert2(2),vert3(1),vert3(2),wplot)
c 
        call trplpnt(x0,y0,wplot)
c 
        iw=21
        call trplwrt(iw,wplot,'triangle as created*')
c 
c       . . . calculate the integrals via the two-dimensional quadrature
c 
        m=19
        eps=1.0d-8
c 
        x0y0(3)=0
  
  
        call triaadam(ier,vert1,vert2,vert3,fun6,n,
     1      x0y0,www,m,eps,rints2,maxrec,numfunev,work)
c 
        call prinf('after triaadam, ier=*',ier,1)
        call prin2('after triaadam, rints2=*',rints2,n)
        call prin2('while rints=*',rints,n)
        call prinf('and numfunev=*',numfunev,1)
c 
        do 2600 i=1,n
        diff(i)=rints2(i)-rints(i)
 2600 continue
c 
        call prin2('and diff=*',diff,n)
c 
        stop
        end
c 
c 
c 
c 
c 
        subroutine fun6(x,y,z0,www,f)
        implicit real *8 (a-h,o-z)
        save
        dimension f(1),z0(3),www(1)
c 
c       evaluate all orthogonal polynomials at this point
c 
        call trpeval(x,y,www,f)
c 
        nord=www(3)
        x0=z0(1)
        y0=z0(2)
        z00=z0(3)
c 
        d=sqrt( (x-x0)**2+(y-y0)**2+z00**2 )
c 
        nfuns=(nord+1)*(nord+2)/2
c 
        do 1400 i=1,nfuns
        f(i)=f(i)/d
 1400 continue
c 
        return
        end
c 
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code and the beginning of the
c        integration code proper
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        subroutine trpinteg(x0,y0,z0,w,rints,work)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension w(1),rints(1),work(1)
c 
c       For a user-specified point (x0,y0,z0) in R^3, this subroutine
c       evaluates the integrals
c 
c    \int P_{i,j} (x,y) / sqrt( (x-x0)^2+(y-y0)^2+(z-z0)^2) dx dy     (1)
c 
c       over the standard triangle vith the vertices
c 
c       (0, 2/sqrt(3)), (-1,-1/sqrt(3)), (1,-1/sqrt(3)),          (2)
c 
c       and i+j .le. nord.                                        (3)
c 
c       In (1), P_{i,j} denotes th (i,j) -th orthogonal polynomial
c       on the standard triangle. Obviously, the selection of such
c       polynomials is not unique, and our choice is somewhat
c       random. For our purposes, it is sufficient that the total
c       order of orthogonal polynomials in (1) is non-decreasing,
c       and that the output of this subroutine is consistent with
c       that of trpeval (see). In other words, the m-th element of
c       the array rints produced by this subroutine contains the
c       integral of the polynomial whose value is returned by the
c       subroutine trpeval in the m-th element of the array pols.
c 
c       IMPORTANT: This subroutine uses the array w that (hopefully)
c       has been created by a prior call to the subroutine trpinit
c       (see). The parameter nord in (3) is among a fairly large
c       amount of data stored by trpinit in the array w. Of course,
c       contents of the array w do not depend on the point (x0,y0,z0),
c       depending only on nord; again, see trpinit for details.
c 
c       Please note that nord must be .le. 12.
c 
c             input parameters:
c 
c  x0,y0,z0 - the coordinates of the target point in (1)
c  w - an array that has been (hopefully) created by a prior call to
c       the subroutine trpinit
c 
c             output parameters:
c 
c  rints - integrals (1) with the limits on i, j defined by (3)
c 
c             work arrays:
c 
c  work - ideally, must be at least 101 *(nord+1)*(nord+2) + 20
c         real *8 eelements long. However, it is only used as a work
c         array for the subroutine adapgaum; see the latter for details
c 
c 
c       . . . construct the memory map for the evaluation of the
c             orthogonal polynomials at the point (x,y)
c 
        iicontr=w(1)
        iccontr=w(2)
        nord=w(3)
        nsteps=w(5)
c 
c 
c       evaluate the integrals of functions of the form
c       x**i *y**j /sqrt((x-x0)**2+(y-y0)**2+z0**2)
c 
        n=(nord+1)*(nord+2)/2
c 
        call trptrint(x0,y0,z0,nord,rints,work)
c 
c       linearly combine the obtained integrals to convert
c       them into integrals of functions of the form
c       p_{n,m} (x,y) /sqrt((x-x0)**2+(y-y0)**2+z0**2)
c 
        call trpfunev(w(iicontr),w(iccontr),nsteps,rints)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine trpeval(x,y,w,pols)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension w(1),pols(1)
c 
c       This subroutine evaluates at the user-specified point
c       (x,y) \ in R^2 a collection of orthogonal polynomials
c       of orders 0,1,...,nord, with nord .leq. 12. The values
c       of these polynomials are returned in the array pols;
c       they are ordered in non-decreasing, order, so that
c 
c    pols(1) contains the value of the polynomial of order 0
c    pols(2), pols(3)  contain the values of the polynomials
c       of order 1
c    pols(4), pols(5), pols(6)  contain the values of the
c       polynomials of order 2, etc.
c 
c       Obviously, such a collection of polynomials is not unsique.
c       For our purposes, it is sufficient that the total order of
c       orthogonal polynomials whose values are returned is
c       non-decreasing, and that the output of this subroutine is
c       consistent with that of trpinteg (see). In other words,
c       the m-th element of the array rints produced by trpinteg
c       contains the integral of the polynomial whose value is
c       returned by this subroutine trpinteg in the m-th element of
c       the array pols.
c 
c       IMPORTANT: This subroutine uses the array w that (hopefully)
c       Has been created by a prior call to the subroutine trpinit
c       (see).
c 
c             input parameters:
c 
c  x,y - the coordinates of the point where the polynomials are to
c       be evaluated
c  w - an array that has been (hopefully) created by a prior call to
c       the subroutine trpinit
c 
c             output parameters:
c 
c  pols - the values of the polynomials
c 
c             work arrays: none
c 
c       . . . construct the memory map for the evaluation of the
c             orthogonal polynomials at the point (x,y)
c 
c 
c       construct the memory map for the evaluation of the
c       orthogonal polynomials at the point (x,y)
c 
        iicontr=w(1)
        iccontr=w(2)
        nord=w(3)
        nsteps=w(5)
c 
        call trppolev(x,y,nord,w(iicontr),w(iccontr),nsteps,pols)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine trpinit(nord,w,keep)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension w(1),zswand(2,200),wswand(200),kks(30)
c 
        data kks/1,3,6,6,7,12,15,16,19,25,28,36,40,46,54,58,
     1       66,73,82,85,93,100,106,118,126,138,145,154,166,175/
c 
c       This is an initialization subroutine for the subroutines
c       trpeval, trpinteg that evaluate and integrate with a weight
c 
c        1/sqrt( (x-x0)**2+ (y-y0)**2+ (z-z0)**2)                 (1)
c 
c       a set of polynomials of (x,y) orthogonal on the standard
c       triangle  with the vertices
c 
c       (0, 2/sqrt(3)), (-1,-1/sqrt(3)), (1,-1/sqrt(3)),          (2)
c 
c       and the cumulative order from 0 to nord; in (1), (x0,y0,z0)
c       is a user-specified point in R^3. The actual evaluation
c       is performed by the subroutine trpeval (see), and the actual
c       integrals are evaluated by trpinteg (see). Please note that
c       nord must be .le. 12. Also, please note that the principal
c       output of this subroutine is the array w (see description
c       below), to be used by the subroutines trpeval, trpinteg
c       (see); this subroutine has no use as a stand-alone device.
c 
c             input parameters:
c 
c  nord - the maximum order of the polynomials to be evaluated and
c       integrated; must be between 1 and 12.
c 
c             output parameters:
c 
c  w - an array to be used by the subroutines trpeval, trpinteg;
c       its first keep elements should not be changed between the
c       call to this subroutine and subsequent calls to trpeval,
c       trpinteg. Has to be at least 20 000 real *8 elements long.
c  keep - the length of the chunk of array w that should not be
c       changed between the call to this subroutine and subsequent
c       calls to trpeval, trpinteg. It is always much less than
c       the 35 000 elements used by this subroutine as work space.
c 
c        . . . allocate memory for the arrays icontr, ccontr
c 
        ninire=2
        iicontr=10
        licontr=(9000*3)/ninire+100
c 
        iccontr=iicontr+licontr
        lccontr=9000
c 
        iamatr=iccontr+lccontr
c 
        nfuns=(nord+1)*(nord+2)/2
        nwand=26
        kk=kks(nwand)
c 
        lamatr=kk*nfuns
c 
cccc        call prinf('lamatr=*',lamatr,1)
c 
c        construct the control tables to be used in the
c        evaluation of orthogonal polynomials and of corresponding
c        integrals
c 
        ifpivot=0
cccc        call trportpo(nord,nwand,amatr,kk,ifpivot,w(iicontr),
        call trportpo(nord,nwand,w(iamatr),kk,ifpivot,w(iicontr),
     1          w(iccontr),nsteps,zswand,wswand)
c 
cccc        call prinf('after trportpo, w(iicontr)=*',w(iicontr),3*nsteps)
cccc        call prin2('after trportpo, w(iccontr)=*',w(iccontr),nsteps)
cccc        call prinf('after trportpo, nsteps=*',nsteps,1)
cccc        call prinf('after trportpo, kk=*',kk,1)
cccc        call prinf('after trportpo, nord=*',nord,1)
c 
        w(1)=iicontr+0.01
        w(2)=iccontr+0.01
        w(3)=nord+0.01
        w(4)=kk+0.01
        w(5)=nsteps+0.01
        keep=300+(nsteps+3)/ninire+nsteps
        w(6)=keep+0.01
        return
        end
c 
c 
c 
c 
c 
        subroutine trpfunev(icontr,ccontr,nsteps,funs)
        implicit real *8 (a-h,o-z)
        save
        dimension icontr(3,1),ccontr(1),funs(1)
c 
c        evaluate the orthogonalized functions at the point x by using
c        the control table icontr and the corresponding coefficinets
c        in ccontr
c 
        do 3000 k=1,nsteps
c 
c       if this operation is a swap - perform it
c 
        if(icontr(1,k) .ne. 1) goto 2200
c 
        i=icontr(2,k)
        j=icontr(3,k)
        d=funs(i)
        funs(i)=funs(j)
        funs(j)=d
c 
        goto 3000
c 
 2200 continue
c 
c       if this operation is addition of one element of funs
c       to another - act accordingly
c 
        if(icontr(1,k) .ne. 2) goto 2400
c 
        i=icontr(2,k)
        j=icontr(3,k)
        d=ccontr(k)
        funs(j)=funs(j)-d*funs(i)
c 
        goto 3000
c 
 2400 continue
c 
c       if this operation is scaling - act accordingly
c 
        i=icontr(2,k)
        d=ccontr(k)
        funs(i)=d*funs(i)
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
        subroutine trppolev(x,y,nord,icontr,ccontr,nsteps,pols)
        implicit real *8 (a-h,o-z)
        save
        dimension icontr(3,1),ccontr(1),pols(1)
c 
c        evaluate the products of powers of x, y
c 
        ii=0
        do 1400 i1=1,nord+1
        do 1200 j1=1,nord+1
c 
        i=i1-1
        j=j1-1
c 
        if(i+j .gt. nord) goto 1200
c 
        ii=ii+1
        pols(ii)=x**i*y**j
 1200 continue
 1400 continue
c 
         nfuns=(nord+1)*(nord+2)/2
c 
c        evaluate the orthogonal polynomials at the point x by using
c        the control table icontr and the corresponding coefficinets
c        in ccontr
c 
        call trpfunev(icontr,ccontr,nsteps,pols)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine trportpo(nord,nwand,amatr,kk,ifpivot,
     1      icontr,ccontr,nsteps,zs,ws)
        implicit real *8 (a-h,o-z)
        save
        dimension zs(2,200),ws(200),icontr(3,1),ccontr(1),
     1      amatr(kk,1),rnorms(1000)
c 
c        retrieve from disk the wandzura's weights
c 
        eps=1.0d-14
c 
        n=(nord+1)*(nord+2)/2
        call trpwand(ier,nwand,zs,ws,kk1)
c 
        d=3
        d=sqrt(d)
        d=sqrt(d)
        do 1200 i=1,kk
        ws(i)=ws(i)/d
 1200 continue
c 
c        construct the table of powers of x,y at the nodes to be
c        linearly combined into orthogonal polynomials
c 
        do 2000 k=1,kk
c 
        ii=0
        do 1600 i1=1,nord+1
        do 1400 j1=1,nord+1
c 
        i=i1-1
        j=j1-1
c 
        if(i+j .gt. nord) goto 1400
c 
        ii=ii+1
        amatr(k,ii)=zs(1,k)**i * zs(2,k)**j *sqrt(ws(k))
 1400 continue
 1600 continue
c 
 2000 continue
c 
c       conduct the gram-schmidt process
c 
        call trpgrm(amatr,kk,n,rnorms,rcond,ifpivot,
     1      icontr,ccontr,nsteps)
  
        call prin2('after trpgrm, rnorms=*',rnorms,n)
  
        dmin=1.0d39
        dmax=0
        do 2200 i=1,n
        if(abs(rnorms(i)) .lt. dmin) dmin=abs(rnorms(i))
        if(abs(rnorms(i)) .gt. dmax) dmax=abs(rnorms(i))
 2200 continue
  
        call prin2('and condition number is*',dmax/dmin,1)
c 
c        scale the matrix amatr back to compensate for the weights
c 
        do 2600 i=1,n
        do 2400 j=1,kk
c 
        amatr(j,i)=amatr(j,i)/sqrt(ws(j))
 2400 continue
 2600 continue
c 
c       scan the arrays icontr, ccontr, eliminating those
c       operations that do not do anything
c 
cccc        call prinf('before compression, nsteps=*',nsteps,1)
c 
        ii=0
        do 4000 i=1,nsteps
c 
c       if this is a swap of an element with itself - skip it
c 
        if(icontr(1,i) .ne. 1) goto 2800
        if(icontr(2,i) .eq. icontr(3,i)) goto 4000
c 
 2800 continue
c 
c       if this is a subtraction of a zero from something - skip it
c 
        if(icontr(1,i) .ne. 2) goto 3000
        if(abs(ccontr(i)) .lt. eps) goto 4000
 3000 continue
c 
        ii=ii+1
        icontr(1,ii)=icontr(1,i)
        icontr(2,ii)=icontr(2,i)
        icontr(3,ii)=icontr(3,i)
c 
        ccontr(ii)=ccontr(i)
c 
 4000 continue
c 
        nsteps=ii
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine trpgrm(b,m,n,rnorms,rcond,ifpivot,
     1      icontr,ccontr,nsteps)
        implicit real *8 (a-h,o-z)
        save
        real *8 b(m,n),cd
        dimension rnorms(1),icontr(3,1),ccontr(1)
c 
c       . . . prepare the array of values of norms
c             of columns
c 
        done=1
        dtot=0
        do 1400 i=1,n
        d=0
        call trpscap(b(1,i),b(1,i),m,d)
c 
        rnorms(i)=d
        dtot=dtot+d
 1400 continue
c 
c       . . . conduct gram-schmidt iterations
c 
        ii=0
c 
        do 4000 i=1,n
c 
c       find the pivot
c 
        ipivot=i
        rn=rnorms(i)
c 
        do 2200 j=i+1,n
        if(rnorms(j) .le. rn) goto 2200
        rn=rnorms(j)
        ipivot=j
 2200 continue
 2400 continue
c 
        if(ifpivot .eq. 0) ipivot=i
c 
c       put the column number ipivot in the i-th place
c 
        do 2600 j=1,m
        cd=b(j,i)
        b(j,i)=b(j,ipivot)
        b(j,ipivot)=cd
 2600 continue
c 
        d=rnorms(ipivot)
        rnorms(ipivot)=rnorms(i)
        rnorms(i)=d
c 
        ii=ii+1
        icontr(1,ii)=1
        icontr(2,ii)=ipivot
        icontr(3,ii)=i
        ccontr(ii)=0
c 
 2700 continue
c 
c       orthogonalize the i-th column to all preceding ones
c 
        if(i .eq. 1) goto 2790
c 
        do 2780 j=1,i-1
c 
        call trpscap(b(1,i),b(1,j),m,cd)
c 
        do 2770 l=1,m
        b(l,i)=b(l,i)-b(l,j)*cd
 2770 continue
  
cccc        call prin2('orthogonalizing, cd=*',cd,1)
c 
        ii=ii+1
        icontr(1,ii)=2
        icontr(2,ii)=j
        icontr(3,ii)=i
        ccontr(ii)=cd
c 
 2780 continue
 2790 continue
c 
c       normalize the i-th column
c 
        call trpscap(b(1,i),b(1,i),m,cd)
c 
        d=cd
c 
        d=done/dsqrt(d)
c 
        do 2800 j=1,m
        b(j,i)=b(j,i)*d
 2800 continue
c 
        ii=ii+1
        icontr(1,ii)=3
        icontr(2,ii)=i
        icontr(3,ii)=0
c 
        ccontr(ii)=d
c 
c        orthogonalize everything else to it
c 
        do 3200 j=i+1,n
c 
        call trpscap(b(1,i),b(1,j),m,cd)
c 
        rrn=0
        do 3000 l=1,m
        b(l,j)=b(l,j)-b(l,i)*cd
        rrn=rrn+b(l,j)*b(l,j)
 3000 continue
  
cccc        call prin2('orthogonalizing, cd=*',cd,1)
  
c 
        ii=ii+1
        icontr(1,ii)=2
        icontr(2,ii)=i
        icontr(3,ii)=j
        ccontr(ii)=cd
  
        rnorms(j)=dsqrt(rrn)
 3200 continue
 4000 continue
  
 4200 continue
c 
        rcond=rnorms(1)/rnorms(n)
        nsteps=ii
        return
        end
c 
c 
c 
c 
c 
        subroutine trpscap(x,y,n,prod)
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
        subroutine trptrint(x0,y0,z0,nord,rints,work)
        implicit real *8 (a-h,o-z)
        save
        dimension rints(1),zs(2,200),ws(200),work(1),x0y0(4)
c 
        external trpvert
c 
        data ifcalled/0/
c 
c       for a user-specified point (x0,y0,z0) in R^3, this subroutine
c       evaluates the integrals
c 
c      \int x^i * y^j / sqrt( (x-x0)^2 + (y-y0)^2) dx dy          (1)
c 
c       over the standard triangle vith the vertices
c 
c       (0, 2/sqrt(3)), (-1,-1/sqrt(3)), (1,-1/sqrt(3)),          (2)
c 
c       and i+j .le. nord;                                        (3)
c 
c       please note that nord must be .le. 12.
c 
c             input parameters:
c 
c  x0,y0 - the coordinates of the target point in (1)
c  nord - the upper limit in (3)
c 
c             output parameters:
c 
c  rints - integrals (1) with the limits on i, j defined by (3)
c 
c             work arrays:
c 
c  work - ideally, must be at least 101 *(nord+1)*(nord+2) + 20
c         real *8 eelements long. However, it is only used as a work
c         array for the subroutine adapgaum; see the latter for details
c 
c 
c       . . . if this is the first call to this subroutine
c             - precompute everything that needs to be precomputed
c 
        if(ifcalled .eq. 1) goto 1200
c 
        ym=-1/sqrt(3.0d0)
        yp=2/sqrt(3.0d0)
c 
        nwand=30
        call trpwand(ier,nwand,zs,ws,kk)
c 
        d=3
        d=sqrt(d)
        d=sqrt(d)
        do 1100 i=1,kk
        ws(i)=ws(i)/d
 1100 continue
c 
        ifcalled=1
 1200 continue
c 
c       determine whether the integrals should be evaluated
c       via adaptive quadrature in y combined with analytic
c       formula in x, or by a high-order wandzura schems
c 
        d=x0**2+y0**2 + z0**2
cccc        if(d .gt. 4) goto 2000
cccc        if(d .gt. 2.3) goto 2000
        if(d .gt. 3.3) goto 2000
c 
c       the point is close to the triangle; use adaptive quadrature
c 
        x0y0(1)=x0
        x0y0(2)=y0
        x0y0(3)=z0
c 
        eps=1.0d-10
        m=18
        n=(nord+1)*(nord+2)/2
c 
        call adapgaum(ier,ym,yp,trpvert,n,x0y0,nord,m,eps,
     1      rints,maxrec,numint,work)
c 
        return
c 
 2000 continue
c 
c       the point is far from the triangle; use wandzura's formula
c 
        n=(nord+1)*(nord+2)/2
        do 2200 i=1,n
        rints(i)=0
 2200 continue
c 
        do 3000 k=1,kk
c 
        d=sqrt( (zs(1,k)-x0)**2+(zs(2,k)-y0)**2+z0**2)
        d=1/d
c 
        ii=0
        do 2600 i1=1,nord+1
        do 2400 j1=1,nord+1
c 
        i=i1-1
        j=j1-1
        if(i+j .gt. nord) goto 2400
c 
        ii=ii+1
        dd=d* zs(1,k)**i * zs(2,k)**j
        rints(ii)=rints(ii)+dd*ws(k)
 2400 continue
 2600 continue
 3000 continue
  
  
  
  
        return
        end
c 
c 
c 
c 
c 
        subroutine trpwand(ier,n,z,w,kk)
        implicit real *8 (a-h,o-z)
        save
        dimension z(2,1),w(1)
c 
c        this subroutine constructs (or rather, retrieves)
c        Wandzura's quadrature formulae for smooth functions
c        on the triangle with vertices
c 
c 
c          (0,2/sqrt(3)), (-1,-1/sqrt(3)), (-1,-1/sqrt(3).
c 
c                        input parameters:
c 
c  n - the order of the quadrature (the maximum order of
c        the polynomials of two variables that are integrated
c        exactly. must be between 0 and 30
c 
c                        output parameters:
c 
c  ier - error return code.
c          ier=0 means successful execution
c          ier=4 means that n is outside the interval [0,30]
c  z - the nodes in the plane (all inside the treiangle (1) above
c  w - weights corresponding to nodes z(i) (all positive)
c 
        ier=0
        if( (n .gt. 0) .and. (n .lt. 31) ) goto 1200
        ier=4
        return
 1200 continue
c 
c        construct the wights if  0 .le. n .le. 5
c 
        if(n .gt. 4) goto 1400
c 
        call wandwht0(n+1,z,w,kk)
        goto 2000
 1400 continue
c 
c        construct the wights if  5 .le. n .le. 30
c 
        call wandwhts(ier,n+1,z,w,kk)
c 
        three=3
        coef=3*dsqrt(three)/4
        do 1600 i=1,kk
        w(i)=w(i)*coef
 1600 continue
c 
 2000 continue
  
        call trpwandr(z,w,kk)
        return
        end
c 
c 
c 
c 
c 
        subroutine trpwandr(zs,ws,n)
        implicit real *8 (a-h,o-z)
        save
        dimension zs(2,1),ws(1)
c 
c       shrink wandzura's nodes
c 
        do 1200 i=1,n
        zs(1,i)=zs(1,i)*sqrt(3.0d0)/1.5d0
        zs(2,i)=zs(2,i)*sqrt(3.0d0)/1.5d0
        ws(i)=3*ws(i)/2.25d0*sqrt(sqrt(3.0d0))
 1200 continue
c 
c       rotate the wandzura nodes
c 
        done=1
        pi=atan(done)*4
        theta=-pi/2
c 
        u11=cos(theta)
        u22=cos(theta)
        u12=sin(theta)
        u21=-sin(theta)
c 
        do 1400 i=1,n
        d1=u11*zs(1,i)+u12*zs(2,i)
        d2=u21*zs(1,i)+u22*zs(2,i)
c 
        zs(1,i)=d1
        zs(2,i)=d2
 1400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine trpvert(y,n,x0y0,nord,vals)
        implicit real *8 (a-h,o-z)
        save
        dimension x0y0(1),rintsm(15),rintsp(15),vals(1)
c 
c        determine the limits of horizontal integration integration
c 
        xm=y/sqrt(3.0d0)-2/3.0d0
        xp=-xm
c 
c       evaluate all horizontal integrals
c 
        x0=x0y0(1)
        y0=x0y0(2)
        z0=x0y0(3)
c 
        a=x0
        b=sqrt( (y-y0)**2+z0**2)
c 
        call trpintho(a,b,xm,rintsm,nord)
        call trpintho(a,b,xp,rintsp,nord)
c 
        ii=0
        do 1400 i1=1,nord+1
        do 1200 j1=1,nord+1
c 
        i=i1-1
        j=j1-1
c 
        if(i+j .gt. nord) goto 1200
c 
        ii=ii+1
        vals(ii)=(rintsp(i1)-rintsm(i1)) * y**j
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
        subroutine trpintho(a,b,x,rints,nord)
        implicit real *8 (a-h,o-z)
        save
        dimension rints(1)
c 
c        this subroutine evaluates the indefinite integrals
c        of the form
c 
c        \int  x^k/Sqrt(b**2 + (-a + x)**2) dx               (1)
c 
c        for k=0, 1, . . . , nord                            (2)
c 
c        Input parameters:
c 
c  a, b - parameters in (1) above
c  x - the value at which (1) will be evaluated
c  rints - the integrals (1)
c  nord - the upper limit in (2)
c 
c 
c 
c 
        b2=b*b
        b3=b2*b
        b4=b3*b
        b5=b4*b
        b6=b5*b
        b7=b6*b
        b8=b7*b
        b9=b8*b
        b10=b9*b
        b11=b10*b
        b12=b11*b
c 
        a2=a*a
        a3=a2*a
        a4=a3*a
        a5=a4*a
        a6=a5*a
        a7=a6*a
        a8=a7*a
        a9=a8*a
        a10=a9*a
        a11=a10*a
        a12=a11*a
c 
        x2=x*x
        x3=x2*x
        x4=x3*x
        x5=x4*x
        x6=x5*x
        x7=x6*x
        x8=x7*x
        x9=x8*x
        x10=x9*x
        x11=x10*x
        x12=x11*x
c 
        sq22=Sqrt(b2 + (a - x)**2)
        dl22=log(sq22+x-a)
c 
        rints(1)=Log(-a + sq22 + x)
c 
        if(nord. eq. 0) return
c 
        rints(2)=sq22 + a*dl22
c 
        if(nord. le. 1) return
c 
        rints(3)=(sq22*(3*a + x) +
     1      (2*a2 - b2)*dl22)/2
c 
        if(nord. le. 2) return
c 
        rints(4)=(sq22*(11*a2 - 4*b2 +
     1        5*a*x + 2*x2))/6 + (a3 - (3*a*b2)/2)
     2      *dl22
c 
        if(nord. le. 3) return
c 
        rints(5)= (-(sq22*(-50*a3 + 55*a*b2
     1      - 26*a2*x + 9*b2*x - 14*a*x2 - 6*x3)) +
     2      3*(8*a4 - 24*a2*b2 + 3*b4)*dl22)/24
c 
        if(nord. le. 4) return
c 
        rints(6)=(sq22*(274*a4 - 607*a2*b2
     1      + 64*b4 + 154*a3*x - 161*a*b2*x + 94*a2*x2 -
     2      32*b2*x2 + 54*a*x3 + 24*x4) + 15*a*(8*a4 -
     3      40*a2*b2 + 15*b4)*dl22)/120
c 
        if(nord. le. 5) return
c 
        rints(7)=(sq22*(588*a5 + 348*a4*x
     1      + 75*b4*x - 50*b2*x3 + 40*x5 + a3*(-2184*b2
     2      + 228*x2) - 4*a2*(183*b2*x - 37*x3) + a*(693*b4
     3      - 234*b2*x2 + 88*x4)) + 15*(16*a6 - 120*a4*b2
     4      + 90*a2*b4 - 5*b6)*dl22)/240
c 
        if(nord. le. 6) return
c 
        rints(8)=-(sq22*(-4356*a6 - 2676*a5*x
     1      + 12*a4*(2033*b2 - 153*x2) + 4*a3*(2358*b2*x
     2      - 319*x3) + a2*(-15525*b4 + 3786*b2*x2 - 856*x4)
     3      + a*(-2907*b4*x + 1298*b2*x3 - 520*x5) + 48*(16*b6
     4      - 8*b4*x2 + 6*b2*x4 - 5*x6)))/1680 +
     5      (a*(16*a6 - 168*a4*b2 + 210*a2*b4 - 35*b6)*
     6      dl22)/16
c 
        if(nord. le. 7) return
c 
        rints(9)=-(sq22*(-36528*a7 -
     1      23088*a6*x + 24*a5*(11989*b2 - 682*x2) + 8*a4*
     2      (15333*b2*x - 1486*x3) - 2*a3*(152967*b4 -
     3      28248*b2*x2 + 4264*x4) - 2*a2*(37899*b4*x -
     4      12136*b2*x3 + 2920*x5) + a*(45477*b6 -
     5      17226*b4*x2 + 8632*b2*x4 - 3600*x6) +
     6      35*(105*b6*x - 70*b4*x3 + 56*b2*x5 -
     7      48*x7)))/13440 + ((128*a8 - 1792*a6*b2 +
     8      3360*a4*b4 - 1120*a2*b6 + 35*b8)*dl22)/128
c 
        if(nord. le. 8) return
c 
        rints(10)=(sq22*(114064*a8 +
     1      73744*a7*x + a6*(-1202984*b2 + 53584*x2) +
     2      a5*(-550968*b2*x + 40144*x3) + 2*a4*(961437*b4
     3      - 139272*b2*x2 + 15032*x4) + a3*(568722*b4*x
     4      - 137072*b2*x3 + 22000*x5) + a2*(-572519*b6 +
     5      170190*b4*x2 - 61032*b2*x4 + 15280*x6) +
     6      a*(-82841*b6*x + 41574*b4*x3 - 22200*b2*x5 +
     7      9520*x7) + 128*(128*b8 - 64*b6*x2 + 48*b4*x4
     8      - 40*b2*x6 + 35*x8)) + 315*a*(128*a8 -
     9      2304*a6*b2 + 6048*a4*b4 - 3360*a2*b6 +
     a      315*b8)*dl22)/40320
c 
        if(nord. le. 9) return
c 
        rints(11)= ( sq22*(236192*a9
     1      + 155552*a8*x - 32*a7*(100463*b2 - 3601*x2)
     2      - 32*a6*(48624*b2*x - 2761*x3) +
     3      4*a5*(1802163*b4 - 210444*b2*x2 +
     4      17048*x4) + 4*a4*(603555*b4*x - 113500*b2*x3
     5      + 13016*x5) - 4*a3*(895510*b6 - 214695*b4*x2
     6      + 57840*b2*x4 - 9656*x6) - 4*a2*(176065*b6*x
     7      - 68955*b4*x3 + 26328*b2*x5 - 6776*x7) +
     8      a*(307835*b8 - 124150*b6*x2 + 69960*b4*x4 -
     9      38896*b2*x6 + 17024*x8) + 63*(315*b8*x -
     a      210*b6*x3 + 168*b4*x5 - 144*b2*x7 + 128*x9))
     b      + 315*(256*a10 - 5760*a8*b2 + 20160*a6*b4 -
     c      16800*a4*b6 + 3150*a2*b8 - 63*b10)*
     d      dl22 )/80640
c 
        if(nord. le. 10) return
c 
        rints(12)=-(sq22 *(-2678752*a10 -
     1      1791712*a9*x + 32*a8*(1429148*b2 - 42131*x2)
     2      + 32*a7*(722839*b2*x - 32891*x3) - 4*a6*
     3      (34245973*b4 - 3294524*b2*x2 + 207688*x4) -
     4      4*a5*(12585885*b4*x - 1898860*b2*x3 +
     5      163336*x5) + 4*a4*(25557485*b6 - 5050365*b4*x2
     6      + 1059800*b2*x4 - 126376*x6) + 4*a3*(6126770*b6*x
     7      - 1925545*b4*x3 + 551888*b2*x5 - 94696*x7) +
     8      a2*(-17587235*b8 + 5847110*b6*x2 -
     9      2566920*b4*x4 + 1020016*b2*x6 - 267904*x8)
     a      + a*(-2073565*b8*x + 1109310*b6*x3 -
     b      666264*b4*x5 + 380912*b2*x7 - 169344*x9) +
     c      1280*(256*b10 - 128*b8*x2 + 96*b6*x4 -
     d      80*b4*x6 + 70*b2*x8 - 63*x10)))/887040 +
     e      (a*(256*a10 - 7040*a8*b2 + 31680*a6*b4 -
     f      36960*a4*b6 + 11550*a2*b8 - 693*b10)*
     g      dl22 )/256
c 
        if(nord. le. 11) return
c 
        rints(13)= -(sq22 *(-11010688*a11
     1      - 7462528*a10*x + 64*a9*(3601247*b2 - 88882*x2)
     2      + 64*a8*(1888329*b2*x - 70402*x3) - 128*a7*
     3      (6946290*b4 - 560376*b2*x2 + 28271*x4) -
     4      64*a6*(5489859*b4*x - 679436*b2*x3 + 45454*x5)
     5      + 32*a5*(29063632*b6 - 4818501*b4*x2 +
     6      809982*b2*x4 - 72428*x6) + 32*a4*(8026750*b6*x
     7      - 2069895*b4*x3 + 461778*b2*x5 - 56588*x7) -
     8      2*a3*(133447535*b8 - 37228960*b6*x2 +
     9      13094760*b4*x4 - 3904576*b2*x6 + 683648*x8) -
     a      6*a2*(7250345*b8*x - 3182640*b6*x3 +
     b      1488120*b4*x5 - 607936*b2*x7 + 162176*x9) +
     c      3*a*(4976075*b10 - 2087830*b8*x2 + 1254600*b6*x4
     d      - 784624*b4*x6 + 457856*b2*x8 - 206080*x10) +
     e      231*(3465*b10*x - 2310*b8*x3 + 1848*b6*x5 -
     f      1584*b4*x7 + 1408*b2*x9 - 1280*x11)))/3.54816d6
     g      + ((1024*a12 - 33792*a10*b2 + 190080*a8*b4 -
     h      295680*a6*b6 + 138600*a4*b8 - 16632*a2*b10
     i      + 231*b12)*dl22 )/1024
c 
        return
        end
