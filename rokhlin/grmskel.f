        implicit real *8 (a-h,o-z)
        dimension a(10 000 000)
c
        call prini(6,13)
C
C       SET ALL PARAMETERS
C
         PRINT *, 'ENTER n'
         READ *,n
         CALL PRINF('n=*',n,1 )
c        
ccc        m=n*2
        m=n/2
        call testmatr(n,m,a)

       
c
        stop
        end
c
c
c
c
c
        subroutine testmatr(n,m,a)
        implicit real *8 (a-h,o-z)
        save
c
        real *8 ts(10000),a(n,m),pols(10000),w(20 000 000),
     1      rnorms(10 000),b(10 000 000),rnorms3(10000)
        integer *4 icols(10 000),irows(10 000),
     1      ipivots(10 000),ipivots3(10 000)
c
c        construct the test nodes
c
        done=1
        h=2*done/(n-1)
        do 1200 i=1,n
c
        ts(i)=(i-1)*h-1
 1200 continue
c
        call prin2('ts=*',ts,n)
c
c       construct the matrix to be skeletonized
c
        do 1600 i=1,n
c
cccc        call legepols(ts(i),m+2,pols)
c

        do 1400 j=1,m
        d=j
c
cc        a(i,j)=pols(j)
        a(i,j)=ts(i)**d

        if(a(i,j) .lt. 1.0d-20) a(i,j)=0


cccc        a(i,j)=0
 1400 continue
c
 1600 continue


cc        call bidiag_rand(n*n,a)
cc        call bidiag_rand(n*n,a)
cc        call bidiag_rand(n*n,a)

        do 1800 i=1,n
       
cccc        a(i,i)=1
 1800 continue


        call rarrcopy(a,b,n*m)

c
cccc        call prin2('a as created*',a,n*m)

        eps=1.0d-10
c
c        skeletonize the matrix via jacskel
c


        call rleastra(a,n,m,w)
        call rarrcopy(w,a,n*m)


        t1=clotatim()


        call jacskel(a,m,n,eps,ipivots,ncols,rnorms)
cccc        call jacskel(a,n,m,eps,ipivots,ncols,rnorms)


        t2=clotatim()

        call prin2('time from jacskel is*',t2-t1,1)




cccc        call prin2('after jacskel, a=*',a,n*n)
        call prin2('after jacskel, rnorms=*',rnorms,ncols)
        call prinf('after jacskel, ipivots=*',ipivots,ncols)
        call prinf('after jacskel, ncols=*',ncols,1)

c
c        skeletonize the matrix via rskel_basic
c
        call rarrcopy(b,a,n*m)


        t3=clotatim()

cccc        call rskel_piv3(a,n,m,eps,rnorms3,ipivots3,ncols3)
ccccc        call rskel_piv3(a,n,m,eps,ipivots3,rnorms3,ncols3)
        call grmskel(a,n,m,eps,ipivots3,ncols3,rnorms3)

        t4=clotatim()

        call prin2('time from rskel_basic is*',t4-t3,1)

        call prinf('after rskel_basic, ipivots3=*',ipivots3,ncols3)

        call prin2('after rskel_basic, rnorms3=*',rnorms3,ncols3)
c
        return
        end
c 
c 
c 
c 
c 
        subroutine rarrcopy(x,y,n)
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c 
c         This is the end of the debugging code, and the beginning
c         of the skeletonizing code proper.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c 
c 
        subroutine grmskel(b,n,m,eps,ipivots,ncols,rnorms)
        implicit real *8 (a-h,o-z)
        save
        dimension rnorms(1),ipivots(1)
        real *8 b(n,m),cd
c
c        This subroutine skeletonizes the columns of the user-supplied
c        matrix a via the standard double Gram-Schmidt procedure. At 
c        the present time (10.03.04) this is the fastest general-purpose 
c        skeletonizer we have; it is also mercifully short.
c
c                       Input parameters:
c
c  a - the matrix whose columns are to be skeletonized
c  (n,m) - dimensionalities of the matrix a
c  eps - the precision to which a is to be skeletonized
c  
c                       Output parameters:
c
c  ipivots - the row skeleton of a (please note that ipivots HAS 
c        to be at least n integer *4 elements long)
c  ncols - the number of elements returned in ipivots 
c  rnorms - the norms associated with the pivots in array ipivots
c        (in the ideal worls, these would be singular values, but
c        this is NOT the ideal world)
c 
c 
c        . . . initialize the array of pivots
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
c 
        d=0
        do 1200 j=1,n
        d=d+b(j,i)*(b(j,i))
 1200 continue
        rnorms(i)=sqrt(d)
        dtot=dtot+d
 1400 continue
c 
        thresh=dtot*eps**2
  
        thresh=sqrt(thresh)
c 
c       . . . conduct gram-schmidt iterations
c 
        do 4000 i=1,n
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
        cd=b(j,i)
        b(j,i)=b(j,ipivot)
        b(j,ipivot)=cd
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
        call grmprod(b(1,i),b(1,j),n,cd)
c 
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*cd
 2770 continue
 2780 continue
 2790 continue
c 
c       normalize the i-th column
c 
        call grmprod(b(1,i),b(1,i),n,cd)
c 
        d=cd
        if(d .lt. thresh**2 ) return
c 
        ncols=i
c 
        d=done/dsqrt(d)
        do 2800 j=1,n
        b(j,i)=b(j,i)*d
 2800 continue
c 
        if(i .eq. m) goto 3400
c 
c        orthogonalize everything else to it
c 
        do 3200 j=i+1,m
c 
        if(rnorms(j) .lt. thresh) goto 3200
c 
        call grmprod(b(1,i),b(1,j),n,cd)
c 
        cd=(cd)
c 
        rrn=0
        do 3000 l=1,n
        b(l,j)=b(l,j)-b(l,i)*cd
        rrn=rrn+b(l,j)*(b(l,j))
 3000 continue
        rnorms(j)=dsqrt(rrn)
 3200 continue
 3400 continue
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
        subroutine grmprod(x,y,n,prod)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1)
c 
        prod=0
        do 1200 i=1,n
        prod=prod+x(i)*(y(i))
 1200 continue
        return
        end
