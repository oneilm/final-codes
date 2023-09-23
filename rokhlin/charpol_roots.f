        implicit real *8 (a-h,o-z)
        real *8 a(100 000),b(100 000),roots(100000)
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
         call mattest(a,n,b,roots)
  
  
         stop
         end
c 
c 
c 
c 
c 
        subroutine mattest(a,n,b,roots)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(n,n),x(100 000),y(100 000),b(n,n),
     1      aa(10000),bb(10000),
     2      roots(1),diffs(10000),rands(10000),
     3      arr(2,10000),ends_bad(2,10000),w(1000000)
  
  
        integer *4 nums_bad(10000)
c 
c        construct the matrix
c 
        call RIRAND(N*2,rands)
        call RIRAND(N*2,rands)
  
        call prin2('rands as constructed*',rands,n)
  
        done=1
        do 1400 i=1,n
        do 1200 j=1,n
c 
        a(j,i)=0
 1200 continue
c 
        a(i,i)=1000
        a(i,i)=1
        if(i .lt. n) a(i+1,i)=i
        if(i .lt. n) a(i+1,i)=done/i
c 
        do 1300 j=1,n
c 
        a(j,i)=a(j,i)*(rands(i)+1)
cccc        a(j,i)=a(j,i)*sqrt(sqrt(done*i))
 1300 continue
  
  
 1400 continue
c 
        call d2mcopy(a,b,n*n)
  
cccc        call prin2('a as created*',a,n*n)
c 
c       display the matrix
c 
  
        iw=21
        call matplot(a,n,iw,'a as created*',x,y)
c 
c        multiply a by its adjoint, and display the result
c 
        call matmulta(a,a,b,n)
c 
c        normalize b
c 
  
        call d2mscap(b,b,n**2,d)
        d=1/sqrt(d)
  
        d=d*100
c 
        do 1800 i=1,n
        do 1600 j=1,n
c 
        b(j,i)=b(j,i)*d
 1600 continue
 1800 continue
c 
cccc        call prin2('b as created*',b,n*n)
c 
        iw=22
        call matplot(b,n,iw,'a after squaring*',x,y)
c 
c       construct the vectors aa, bb
c 
        do 2200 i=1,n-1
c 
        aa(i)=b(i,i)
        bb(i)=b(i+1,i)
 2200 continue
c 
        aa(n)=b(n,n)
c 
        call prin2('aa=*',aa,n)
        call prin2('bb=*',bb,n)
c 
c        diagonalize b
c 
        eps=1.0d-14
        ifvect=0
        call jaceig(b,n,eps,ifvect,uout,nelems)
  
        call prin2('eigenvalues of b are*',b,n)
  
  
        x1=10
        x1=-0.00001
        x2=100
        lenw=1000000
        eps=1.0d-3
        call charpol_roots(ier,aa,bb,n,x1,x2,eps,ii,roots,
     1      ends_bad,nums_bad,nbad,w,lenw,lused)
  
         call prinf('nums_bad=*',nums_bad,nbad)
         call prin2('ends_bad=*',ends_bad,nbad*2)
         call prinf('nbad=*',nbad,1)
  
  
        call prin2('after charpol_roots0, roots=*',roots,ii)
  
  
  
        call charpol_bubble(b,n)
  
        call prin2('and again, eigenvalues*',b,n)
c 
        do 2800 i=1,n
c 
        diffs(i)=b(i,1)-roots(i)
 2800 continue
c 
        call prin2('and diffs=*',diffs,n)
  
        call prinf('and lused=*',lused,1)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine matplot(a,n,iw,title,x,y)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(n,n)
        dimension x(1),y(1),ww(1000 000)
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
c 
c 
        subroutine d2mcopy(a,b,n)
        implicit real *8 (a-h,o-z)
        save
        dimension a(1),b(1)
c 
        do 1200 i=1,n
        b(i)=a(i)
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine d2mscap(x,y,n,prod)
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
        subroutine matmulta(a,x,c,n)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),x(n,n),c(n,n)
c 
        do 1600 i=1,n
        do 1400 j=1,n
c 
        cd=0
        do 1200 ll=1,n
c 
cccc        cd=cd+a(i,ll)*x(ll,j)
        cd=cd+a(i,ll)*x(j,ll)
 1200 continue
c 
        c(i,j)=cd
 1400 continue
 1600 continue
        return
        end
c 
c 
c 
c 
c 
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
c 
c 
c 
c 
        entry rirandin(ifcall7)
        ifcall=0
        return
        END
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
c        -supplied symmetric positive definite tridiagonal
c        matrix located between the user-supplied limits x1,x2
c        If some of the eigenvalues are so close to each other
c        that separation is difficult, the subroutine returns
c        them in a lump.
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
