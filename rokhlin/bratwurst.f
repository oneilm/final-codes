        implicit real *8 (a-h,o-z)
        real *8 a(5 000 000),c(1000 000),
     1      u(1000000),v(1000000),uout(1000000),
     2      u2(1000000),errs(1000000),u4(1000000),
     3      prods(1000 000)
c 
        CALL PRINI(6,13)
C 
c       . . . construct the testing parameters
c 
         PRINT *, 'ENTER n'
         READ *,n
         CALL PRInf('n=*',n,1)
c 
c        construct the test matrix and plot it
c 
        call mattest(a,n,c,u,v,uout,u2,errs,u4,prods)
  
  
         stop
         end
c 
c 
c 
c 
c 
        subroutine mattest(a,n,c,u,v,uout,u2,errs,u4,prods)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(n,n),x(100 000),y(100 000),
     1      diff(10000),aa(10000),bb(10000),xs(10000),
     2      pols(10000),fs(10000),c(n*2,n*2),
     3      u(n,n),v(n,n),uout(n*2,n*2),u2(n,n),errs(n,n),
     4      rhs(10000),u4(n,n),rlams(10000),
     5      w(1000 000),prods(n,1)
  
c 
c        construct the matrix
c 
        call BIDIAG_RAND(N,Y)
        call BIDIAG_RAND(N,Y)
        call BIDIAG_RAND(N,Y)
  
  
        call BIDIAG_RAND(N,Y)
c 
        done=1
        do 1400 i=1,n
        do 1200 j=1,n
c 
        a(j,i)=0
 1200 continue
c 
        a(i,i)=1
        if(i .lt. n) a(i+1,i)=i
        if(i .lt. n) a(i+1,i)=1.0d-6
c 
        a(i,i)=a(i,i)*(1+y(i))
        if(i .lt. n) a(i+1,i)=a(i+1,i)*(1+y(i))
c 
 1400 continue
c 
c       display the matrix
c 
  
        iw=21
        call matplot(a,n,iw,'a as created*')
c 
c        construct the matrix in the bidiagonal form
c 
        do 1600 i=1,n
c 
        aa(i)=a(i,i)
        if(i .eq. n) goto 1600
c 
        bb(i)=a(i+1,i)
 1600 continue
c 
c       construct the singular values of a via the Sturm
c       procedure
c 
  
        eps=1.0d-13
        lenw=1000 000
  
  
        call prin2('before bidiag_svd, eps=*',eps,1)
  
        call bidiag_svd(ier,aa,bb,n,eps,u4,rlams,nlams,
     1      w,lenw,lused)
  
        call prinf('after bidiag_svd, ier=*',ier,1)
        call prinf('and lused=*',lused,1)
        call prin2('and rlams=*',rlams,nlams)
  
  
cccc        stop
cccc        call prin2('and rlams=*',rlams,n)
c 
c       construct the SVD of A
c 
        eps=1.0d-12
        ifuv=1
        call d2mudv(a,n,u,v,eps,ifuv,numrows)
  
        call prin2('after d2mudv, singular values are*',a,n)
        call prinf('after d2mudv, numrows=*',numrows,1)
c 
c       check the accuracy of u4
c 
  
        do 4200 i=1,n
c 
        do 4000 j=1,n
c 
        d1=u(j,i)-u4(j,i)
        errs(j,i)=u(j,i)+u4(j,i)
        if(abs(errs(j,i)) .gt. abs(d1)) errs(j,i)=d1
c 
 4000 continue
c 
 4200 continue
c 
        call prin2('and errs in u4 are*',errs,n*n)
c 
c       check the orthonormality of the columns of u4
c 
        do 4600 i=1,nlams
        do 4400 j=1,nlams
c 
        call bidiagscap(u4(1,i),u4(1,j),n,prods(j,i))
c 
 4400 continue
 4600 continue
c 
        call prin2('and prods=*',prods,nlams**2)
  
        eps=1.0d-12
        ifuv=0
        call d2mudv(prods,n,u,v,eps,ifuv,numrows)
  
        call prin2('and eigenvalues of prods are*',prods,nlams)
  
        do 4800 i=1,nlams
c 
        prods(i,1)=prods(i,1)-1
 4800 continue
c 
        call prin2('and eigenvalues of prods -1 are*',prods,nlams)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine matplot(a,n,iw,title)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(n,n)
        dimension x(100000),y(100000),ww(1000 000)
        character *4 title(1)
c 
c       this subroutine plots the geometry of a mockup of
c       a square matrix
c 
c 
c                input parameters:
c 
c  a - the mockup of the matrix to be plotted. it is an INTEGER *4
c       n \times n matrix; non-zero elements will be plotted.
c  n - the dimensionality of the matrix to be plotted
c  iw - the FORTRAN unit number on which the plot will be written
c  title - the title to be put on the plot
c 
c                output paramaters: none
c 
c                work arrays:
c 
c  x,y - must be at least n**2 real *8 elements each
c 
c       . . . construct the array of coordinates
c 
  
        call trplopen(ww)
  
  
        m=1
        k=n
        eps=1.0d-10
        num=0
        do 1400 i=1,n
        do 1200 j=1,n
        if(abs(a(j,i)) .lt. eps) goto 1200
        num=num+1
        x(num)=i
        y(num)=-j
c 
        xnum=i
        ynum=-j
        call trplpnt(xnum,ynum,ww)
c 
 1200 continue
 1400 continue
cccc          call prin2('before rsplot, x=*',x,num)
cccc          call prin2('before rsplot, y=*',y,num)
  
c 
c       construct the frame
c 
        call plotframe(iw,n,m,k,ww)
c 
c       plot the junk
c 
  
  
        call trplwrt(iw,ww,title)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine plotframe(iw,n,m,k,w)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1)
c 
c        create the vertical boundaries of blocks
c 
cccc         call prinf('in plotframe, n=*',n,1)
cccc         call prinf('in plotframe, m=*',m,1)
cccc         call prinf('in plotframe, k=*',k,1)
  
        do 2000 i=1,k+1
        y1= - ((i-1)*m+0.5)
        y2=y1
        x1=0.5
        x2=n+0.5
  
        call trplline(x1,y1,x2,y2,w)
  
 1200 format(4x,e11.5,4x,e11.5)
cccc        write(iw,1200) x1,y1
cccc        write(iw,1200) x2,y2
 1400 format('        join')
cccc        write(iw,1400)
 2000 continue
  
c 
c        create the horizontal boundaries of blocks
c 
        do 3000 i=1,k+1
        x1=  ((i-1)*m+0.5)
        x2=x1
        y1=-0.5
        y2=-(n+0.5)
  
        call trplline(x1,y1,x2,y2,w)
  
  
cccc        write(iw,1200) x1,y1
cccc        write(iw,1200) x2,y2
cccc        write(iw,1400)
 3000 continue
        return
        end
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code and the beginning
c        of the code finding the SVD of a bidiagonal matrix
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c        This file contains two user-callable subroutines: bidiag_svd
c        and charpol_roots. Following is a brief description of them
c        things.
c 
c 
c   bidiag_svd - constructs the singular values and the LEFT singular
c        vectors of the user-supplied LOWER bidiagonal matrix A. The
c        singular values are determined  via a version of the Sturm
c        method; subsequently a version of the inverse power method
c        is used to find the corresponding singular vectors.
c 
c   charpol_roots - constructs the eigenvalues of a user-supplied
c        symmetric tridiagonal matrix located between the limits
c        xs,x2 supplied by the user. If some of the eigenvalues are
c        so close to each other that separation is difficult, the
c        subroutine returns them in a lump.
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        subroutine bidiag_svd(ier,aa,bb,n,eps,u,rlams,nlams,
     1      w,lenw,lused)
        implicit real *8 (a-h,o-z)
        save
        dimension aa(1),bb(1),w(1),u(n,n),rlams(1)
c 
c        This subroutine constructs the singular values and the
c        LEFT singular vectors of the user-supplied LOWER bidiagonal
c        matrix A. The singular values are determined  via a version
c        of the Sturm method; subsequently a version of the inverse
c        power method is used to find the corresponding singular
c        vectors. Obviously, this approach becomes difficult when
c        the matrix A has multiple singular values; the subroutine
c        makes no attempt to deal competently with this situation.
c        It checks for multiple singular values; if some are found,
c        there are two possibilities: either they are all less than
c        10*eps (with eps supplied by the user), or some of them are
c        greater than 10*eps. In the first case, the problem is
c        ignored; in the second case, the error return code ier is
c        set to 4, and the subroutine bombs.
c 
c                      Input parameters:
c 
c  aa - the diagonal of the bidiagonal matrix A
c  bb - the subdiagonal of the bidiagonal matrix A
c  n - the dimensionality of the matrix A
c  eps - the machine zero. More specifically, any two singular values
c        s_1, s_2 of A such that
c 
c               |s1-s_2| < s_1 * eps,                               (1)
c 
c        AND
c 
c               |s1-s_2| < s_2 * eps,                               (2)
c 
c        are considered to be too close (please note "AND" above, as
c        opposed to "OR"). Furthermore, any singular value s such that
c        s < 10*eps is considered ignorable.
c  lenw - the amount of space (in real *8 words) in the work array w
c        supplied by the user
c 
c                       Output parameters:
c 
c  ier - error return code.
c    ier = 0 means successful execution
c    ier = 4 means that the matrix A has multiple singular values
c        (according to the definition above) that are greater than
c        eps*10. This error is not quite fatal, in the sense that
c        the subroutine returns the singular values it has found
c    ier > 4 means that the amount of space provided by the user
c        in the work array w is insufficient. This is a fatal error
c  u - the (n,nlams) matrix of left singular vectors of A
c  rlams - the singular values of A
c  nlams - the number of simple (non-multiple) singular values found;
c        also the number of columns returned in u
c  lused - the size of the part of the work array w that was actually
c        used by the subroutine.
c 
c 
c        . . . allocate memory
c 
        ier=0
c 
        irhs=1
        lrhs=n*2+10
c 
        ia=irhs+lrhs
        la=n*2+10
c 
        ib=ia+la
        lb=n*2+10
c 
        ia2=ib+lb
        la2=n*2+10
c 
        ib2=ia2+la2
        lb2=n*2+10
c 
        ic2=ib2+lb2
        lc2=n*2+10
c 
        iuu=ic2+lc2
        luu=n*2+10
c 
        ivv=iuu+luu
        lvv=n*2+10
c 
        iww=ivv+lvv
        lww=n*2+10
c 
        iroots=iww+lww
        lroots=n*2+10
c 
        inums_bad=iroots+lroots
        lnums_bad=n*2+10
c 
        iends_bad=inums_bad+lnums_bad
        lends_bad=n*4+10
c 
        iw2=iends_bad+lends_bad
        lenw2=lenw-iw2-10
c 
        if(lenw2 .lt. 3*n+100) then
            ier=32
            return
        endif
c 
        call bidiag_svd0(jer,aa,bb,n,w(irhs),u,eps,
     1      w(iw2),lenw2,w(ia),w(ib),w(ia2),w(ib2),w(ic2),
     2      w(iuu),w(ivv),w(iww),w(iroots),w(inums_bad),
     3      w(iends_bad),rlams,lused2,nlams)
c 
        lused=lused2+iw2
c 
        ier=jer
        if(ier .gt. 4) ier=32
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine bidiag_svd0(ier,aa,bb,n,rhs,u4,eps,
     1      w,lenw,a,b,a2,b2,c2,uu,vv,ww,roots,nums_bad,
     2      ends_bad,rlams,lused,nlams)
        implicit real *8 (a-h,o-z)
        save
        dimension aa(1),bb(1),a(1),b(1),roots(1),
     1      ends_bad(2,1),nums_bad(1),w(1),a2(1),b2(1),
     2      c2(1),uu(1),vv(1),ww(1),rhs(1),u4(n,n),rlams(1)
c 
c        construct the tridiagonal matrix
c 
        ier=0
c 
        do 1200 i=1,n*2
        a(i)=0
 1200 continue
c 
        ii=0
        d=0
        do 2200 i=1,n
c 
        ii=ii+1
        b(ii)=aa(i)
        d=d+b(ii)**2
        ii=ii+1
c 
        b(ii)=bb(i)
        d=d+b(ii)**2
 2200 continue
c 
c       use Sturm procedure to find all eigenvalues of the
c       tridiagonal matrix
c 
        x1=-sqrt(d)
        x2=-x1
        x2=0
c 
        call charpol_roots(jer,a,b,n*2,x1,x2,eps,
     1      nroots,roots,ends_bad,nums_bad,nbad,w,lenw,
     2      lused)
c 
        if(jer .ne. 0) then
            ier=16
            return
        endif
c 
cccc        call prin2('after charpol_roots, roots=*',roots,nroots)
cccc        call prinf('after charpol_roots, nbad=*',nbad,1)
cccc        call prin2('after charpol_roots, eps=*',eps,1)
c 
        nlams=nroots
        do 2220 i=1,nlams
c 
        rlams(i)=-roots(i)
 2220 continue
c 
c       check whether there are any two relatively large
c       singular values that are too close to each other.
c       If so - bomb
c 
        if(nbad .eq. 0) goto 2300
  
        call prin2('ends_bad=*',ends_bad,nbad*2)
  
c 
        if (ends_bad(1,1) .lt. -10*eps) then
            ier=4
            return
        endif
c 
 2300 continue
c 
c        find the eigenvectors corresponding to the last nroots
c        eigenvalues
c 
        do 2800 i=1,nroots
c 
        d1=1.0d20
        if(i .ne. 1) d1=abs(roots(i-1)-roots(i))
        d2=1.0d20
        if(i .ne. nroots) d2=abs(roots(i+1)-roots(i))
c 
        d=d2
        if(d1 .lt. d2) d=d1
        delta=d/100
        if(d1 .gt. d2) delta=-d/100
c 
        call prin2('delta=*',delta,1)
        do 2400 j=1,n*2
c 
        a2(j)=roots(i)+delta
        b2(j)=b(j)
        c2(j+1)=b(j)
c 
 2400 continue
c 
        call BIDIAG_RAND(N*2,rhs)
c 
        call bidiag_fact(a2,b2,c2,n*2,uu,vv,ww)
c 
        do 2600 j=1,10
  
        call bidiag_solv(uu,vv,ww,n*2,rhs)
c 
        call bidiagscap(rhs,rhs,n*2,d)
c 
        do 2500 l=1,n*2
        rhs(l)=rhs(l)/sqrt(d/2)
 2500 continue
c 
 2600 continue
c 
c        retrieve from the rhs the elements of the appropriate
c        singular vector
c 
        d=0
        do 2700 j=1,n
c 
        u4(j,i)=rhs(1+(j-1)*2)
        d=d+u4(j,i)**2
 2700 continue
c 
        d=1/sqrt(d)
        do 2750 j=1,n
c 
        u4(j,i)=u4(j,i)*d
 2750 continue
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
        subroutine bidiag_fact(a,b,c,n,u,v,w)
        implicit real *8 (a-h,o-z)
        save
        dimension a(1),b(1),c(1),u(1),v(1),w(1),rhs(1)
c 
c        eliminate down
c 
        do 1200 i=1,n-1
        d=c(i+1)/a(i)
        a(i+1)=a(i+1)-b(i)*d
        u(i)=d
 1200 continue
c 
c        eliminate up
c 
        do 1400 i=n,2,-1
        d=b(i-1)/a(i)
        v(i)=d
 1400 continue
c 
c       scale the diagonal
c 
        done=1
        do 1600 i=1,n
        w(i)=done/a(i)
 1600 continue
c 
        return
c 
c 
c 
c 
        entry bidiag_solv(u,v,w,n,rhs)
c 
c        eliminate down
c 
        do 2400 i=1,n-1
        rhs(i+1)=rhs(i+1)-u(i)*rhs(i)
 2400 continue
c 
c        eliminate up
c 
        do 2600 i=n,2,-1
        rhs(i-1)=rhs(i-1)-rhs(i)*v(i)
 2600 continue
c 
c       scale
c 
        do 2800 i=1,n
        rhs(i)=rhs(i)*w(i)
 2800 continue
        return
        end
c 
c 
c 
c 
c 
        SUBROUTINE BIDIAG_RAND(N,Y)
        IMPLICIT REAL *8 (A-H,O-Z)
        save
        DIMENSION Y(1)
C 
C       GENERATE PSEUDO-RANDOM NUMBERS
C 
        done=1
        pi=atan(done)*4
        add=sqrt(7*done)
        x=1.1010101010101
c 
        do 1200 i=1,n
c 
        phi=x*100000*pi+add
cccc        phi=x*1000*pi+add
        j=phi
        phi=(phi-j)
        x=phi
c 
        y(i)=x
 1200 continue
c 
        RETURN
        END
c 
c 
c 
c 
c 
        subroutine bidiagscap(x,y,n,prod)
        implicit real*8 (a-h,o-z)
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
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code and the beginning
c        of the code finding the spectrum of a tridiagonal symmetric
c        matrix
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine charpol_roots(ier,a,b,n,x1,x2,eps,
     1      nroots,roots,ends_bad,nums_bad,nbad,w,lenw,
     2      lused)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(n),b(n-1),roots(1),ends_bad(2,1),w(1)
        integer *4 nums_bad(1)
c 
c        This subroutine constructs the eigenvalues of a user
c        -supplied symmetric tridiagonal matrix located between
c        the user-supplied limits x1,x2. If some of the eigenvalues
c        are so close to each other that separation is difficult,
c        the subroutine returns them in a lump.
c 
c                 Input parameters:
c 
c  a - diagonal of the matrix whose eigenvalues are to be
c        constructed
c  b - subdiagonal (also superdiagonal) of the matrix whose
c        eigenvalues are to be constructed. EXPLANATION:
c        The matrix looks as follows:
c 
c 
c            a(1) b(1) 0    0    0    0  . . .
c            b(1) a(2) b(2) 0    0    0  . . .
c            0    b(2) a(3) b(3) 0    0  . . .
c            0    0    b(3) a(4) b(4) 0  . . .
c            .    .    .    .    .    .  . . .
c            .    .    .    .    .    .  . . .
c            .    .    .    .    .    .  . . .
c  n - the dimensionality of the matrix. Please note that the
c        diagonal a contains n elements, while the subdiagonal
c        b contains n-1 elements
c  x1,x2 - the limits between which the subroutine will be
c        looking for the roots (and God willing, will find
c        all of themn things)
c  eps - the clustering parameter. EXPLANATION: If the subroutine
c        encounters tow or more roots on an interval of length
c        eps or less, it assumes that it can not separate these
c        eigenvalues, and returns them as a lump (see parameters
c        ends_bad, nums_bad,nbad below)
c  lenw - the length of the user-supplied work array w (in
c        real *8 words)
c 
c                 Output parameters:
c  ier - error return code:
c     ier=0 means successful conclusion
c     ier=4 means that the amount of memory supplied by the user
c        in the work array w (and communicated to the subroutine
c        via parameter lenw) is insufficient
c     ier=16 means that the amount of memory supplied by the user
c        in the work array w (and communicated to the subroutine
c        via parameter lenw) is grossly insufficient
c  nroots - the number of roots returned in the array roots
c  roots - the simple roots of the characteristic polynomial found
c        by the subroutine on the interval [x1,x2]
c  ends_bad - the array of ends of "bad" intervals, i.e. intervals
c        of length less than eps containing more than one root
c  nums_bad - integer array containing numbers of roots on each
c        of the "bad" intervals
c  nbads - the number of "bad" intervals
c 
c 
c        . . . allocate memory
c 
        ipols=1
        lpols=n+10
c 
        imap=ipols+lpols
        lenmap=lenw-imap
c 
c        if the amount of space in the array map is
c        clearly insufficient - bomb
c 
        nnneed=120
        if(nnneed .gt. lenmap) then
           ier=16
           return
        endif
c 
c        find the roots
c 
        call charpol_roots0(jer,a,b,n,x1,x2,ii,w(ipols),
     1      w(imap),lenmap,roots,ends_bad,nums_bad,nbad,eps,
     2      nn)
c 
        lused=imap+nn*6+20
c 
        if(jer .ne. 0) then
            ier=4
            return
        endif
c 
        nroots=ii
        call charpol_bubble2(ends_bad,nums_bad,nbad)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine charpol_roots0(ier,a,b,n,x1,x2,ii,pols,
     1      map,lenmap,roots,ends_bad,nums_bad,nbad,eps,nn)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(1),b(1),pols(1),map(6,1),roots(1),
     1      ends_bad(2,1)
        integer *4 nums_bad(1)
c 
c       calculate the sturm index at the ends x1,x2 of the
c       user-supplied interval
c 
        ier=0
        call charpol_sturm(a,b,n,x1,pols(2),ns1)
        call charpol_sturm(a,b,n,x2,pols(2),ns2)
c 
c        initialize the recursive subdivisoon
c 
        map(1,1)=0+0.1
        map(2,1)=x1
        map(3,1)=x2
        map(4,1)=ns1+0.1
        map(5,1)=ns2+0.1
        map(6,1)=0.1
c 
c        construct the subdivision threshold
c 
        nn=1
c 
c        repeatedly scan the map of the interval [x1,x2],
c        subdividing whatever needs to be subdivided
c 
        ifdone=1
        ifbad=0
        do 2000 ijk=1,10000
c 
        ifdone=1
        nn2=nn
c 
        do 1800 i=1,nn2
c 
c       if this chunk has been subdivided - skip it
c 
        ifkids=map(6,i)
        if(ifkids .ne. 0) goto 1800
c 
c       if this chunk has only one root on it - skip it
c 
        nrts=map(5,i)-map(4,i)+0.1
        if(nrts .eq. 1) goto 1800
c 
c        this chunk contains more than one root. if it is too
c        short to be split - mark it and continue
c 
        rl=abs(map(3,i)-map(2,i))
        if(rl .lt. eps*abs(map(3,i)) ) then
            ifdone=0
            map(6,i)=nrts+0.1
            goto 1800
        endif
c 
c        this chunk contains more than one root, and is not too short.
c        split it
c 
        ifdone=0
c 
        x3=(map(3,i)+map(2,i))/2
        call charpol_sturm(a,b,n,x3,pols(2),ns3)
c 
        ns1=map(4,i)
        ns2=map(5,i)
c 
        map(6,i)=1.1
c 
c       process the left subchunk
c 
        nrtsl=ns3-ns1
        if(nrtsl .eq. 0) goto 1400
c 
        nn=nn+1
c 
c        if the amount of space in the array map is
c        insufficient - bomb
c 
        nnneed=nn*6+18
        if(nnneed .gt. lenmap) then
           ier=16
           return
        endif
c 
        map(1,nn)=map(1,i)+1
        map(2,nn)=map(2,i)
        map(3,nn)=x3
        map(4,nn)=ns1+0.1
        map(5,nn)=ns3+0.1
        map(6,nn)=0.1
c 
 1400 continue
c 
c       process the right subchunk
c 
        nrtsr=ns2-ns3
        if(nrtsr .eq. 0) goto 1600
c 
        nn=nn+1
c 
        map(1,nn)=map(1,i)+1
        map(2,nn)=x3
        map(3,nn)=map(3,i)
        map(4,nn)=ns3+0.1
        map(5,nn)=ns2+0.1
        map(6,nn)=0.1
c 
 1600 continue
c 
 1800 continue
c 
        if(ifdone .eq. 1) goto 2100
 2000 continue
c 
 2100 continue
c 
c        scan the map of the interval [x1,x2]. If a subinterval
c        contains exactly one root, find the said root. If a
c        subinterval contains more than one root, and is too
c        short to be subdivided, inform the user accordingly.
c 
        ii=0
        nbad=0
        do 2600 i=1,nn
c 
        ifkids=map(6,i)
        if(ifkids .eq. 1) goto 2600
        if(ifkids .ne. 0) goto 2400
c 
        ii=ii+1
c 
        d1=map(2,i)
        d2=map(3,i)
c 
        call charpol_bisect(a,b,n,d1,d2,x,pols)
c 
        roots(ii)=x
c 
        goto 2600
c 
 2400 continue
c 
c        this subinterval is too short to be subdivided, and
c        contains more than one node. Inform the user accordingly
c 
         nbad=nbad+1
         ends_bad(1,nbad)=map(2,i)
         ends_bad(2,nbad)=map(3,i)
         nums_bad(nbad)=map(6,i)
 2600 continue
c 
c        . . . sort them things
c 
        call charpol_bubble(roots,ii)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine charpols_eval(a,b,n,x,pols)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(1),b(1),pols(1)
c 
c        start the recursion
c 
        pkm2=1
        pkm1=a(1)-x
c 
        i0=0
        pols(i0)=pkm2
        pols(1)=pkm1
c 
        do 1200 i=2,n
c 
        pk=(a(i)-x)*pkm1-b(i-1)**2*pkm2
c 
        pols(i)=pk
c 
        pkm2=pkm1
        pkm1=pk
c 
        d=pkm1**2+pkm2**2
c 
        if(d .gt. 1.0d12) then
            pkm2=pkm2*1.0d-24
            pkm1=pkm1*1.0d-24
        endif
c 
        if(d .lt. 1.0d-12) then
            pkm2=pkm2*1.0d+24
            pkm1=pkm1*1.0d+24
        endif
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
        subroutine charpol_bisect(a,b,n,x1,x2,x,pols)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(1),b(1),pols(1)
c 
c        conduct bisection
c 
        call charpols_eval(a,b,n,x1,pols(2))
        f1=pols(n+1)
c 
        call charpols_eval(a,b,n,x2,pols(2))
        f2=pols(n+1)
c 
cccc        do 2000 k=1,120
        do 2000 k=1,52
c 
c       if x1 or x2 is so close to being a root that never mind,
c       set the output to it and get out
c 
        x3=(x1+x2)/2
        call charpols_eval(a,b,n,x3,pols(2))
        f3=pols(n+1)
c 
        if(f3*f2 .gt. 0) then
            x2=x3
            f2=f3
            goto 2000
        endif
c 
        x1=x3
        f1=f3
c 
 2000 continue
c 
        x=(x1+x2)/2
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine charpol_sturm(a,b,n,x,pols,nsturm)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(1),b(1),pols(1)
c 
c        calculate the sturm sequence
c 
        call charpols_eval(a,b,n,x,pols)
c 
c       calculate the sturm index
c 
        nsturm=0
c 
        do 1400 i=1,n
c 
        dd=pols(i)*pols(i-1)
        if(dd .le. 0) nsturm=nsturm+1
c 
        if(pols(i) .eq. 0) nsturm=nsturm-1
 1400 continue
  
        return
        end
c 
c 
c 
c 
c 
        subroutine charpol_bubble2(ends_bad,nums_bad,nbad)
        implicit real *8 (a-h,o-z)
        save
        dimension ends_bad(2,1),nums_bad(1)
c 
        do 1400 i=1,nbad
        do 1200 j=1,nbad-1
c 
        if(ends_bad(1,j) .le. ends_bad(1,j+1) ) goto 1200
c 
        d=ends_bad(1,j)
        ends_bad(1,j)=ends_bad(1,j+1)
        ends_bad(1,j+1)=d
c 
        d=ends_bad(2,j)
        ends_bad(2,j)=ends_bad(2,j+1)
        ends_bad(2,j+1)=d
c 
        id=nums_bad(j)
        nums_bad(j)=nums_bad(j+1)
        nums_bad(j+1)=id
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
        subroutine charpol_bubble(a,n)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(1)
c 
        do 1400 i=1,n
        do 1200 j=1,n-1
c 
        if(a(j) .le. a(j+1) ) goto 1200
c 
        d=a(j)
        a(j)=a(j+1)
        a(j+1)=d
  
 1200 continue
 1400 continue
c 
        return
        end
