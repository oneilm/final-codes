        implicit real *8 (a-h,o-z)
        dimension a(4000 000)
c
        call prini(6,13)
C
C       SET ALL PARAMETERS
C
         PRINT *, 'ENTER n'
         READ *,n
         CALL PRINF('n=*',n,1 )
c        
        m=n*2
        m=n
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
        real *8 ts(10000),a(n,m),pols(10000),w(1000 000),
     1      rnorms(1000),b(1000 000)
        integer *4 icols(1000),irows(1000),ipivots(1000)
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

        eps=1.0d-14
 

        ifrelskel=1
        call jacskel(a,n,m,eps,ifrelskel,ipivots,ncols,rnorms)


cccc        call prin2('after jacskel, a=*',a,n*n)
        call prin2('after jacskel, rnorms=*',rnorms,ncols)
        call prinf('after jacskel, ipivots=*',ipivots,ncols)
        call prinf('after jacskel, ncols=*',ncols,1)
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
        subroutine jacskel(a,n,m,eps,ifrelskel,ipivots,ncols,rnorms)
        implicit real *8 (a-h,o-z)
        save
        dimension aa(2),uu(2,2),ipivots(1),rnorms(1),a(n,m)
c
c        This subroutine skeletonizes the ROWS of the user-supplied
c        matrix a. Please note that this is a Jacobi-based code, as
c        opposed to a Gram-Schmidt-based one. Due to the author's 
c        lack of experience in this environment, omissions are 
c        possible! Use at your own risk! - 03.06.04
c
c                       Input parameters:
c
c  a - the matrix whose ROWS are to be skeletonized
c  (n,m) - dimensionalities of the matrix a
c  eps - the precision to which a is to be skeletonized
c  ifrelskel - set to 0 when a is to be skeletonized
c              to absolute precision eps;
c              set to 1 when a is to be skeletonized
c              to relative precision eps,
c              relative to the Frobenius norm of a
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
c        . . . initialize the array of pivots
c 
        do 1100 i=1,n
        ipivots(i)=i
 1100 continue
c
c
c        construct the threshold precision
c        as absolute when ifrelskel=0 and
c        as relative when ifrelskel=1
c
        if(ifrelskel .eq. 0) thresh=eps**2
c
        if(ifrelskel .eq. 1) then
c
            dtot=0
            do 1400 i=1,m
            do 1200 j=1,n
            dtot=dtot+a(j,i)**2
 1200 continue
 1400 continue
c
            thresh=dtot*eps**2
c
        endif
c
c 
c       . . . conduct pivoted Jacobi elimination
c 
        do 4000 i=1,n
c
c       recompute the array rnorms
c
        n0=i-2
        if(i .le. 2) n0=1
        call jacskel_rnorms(a,n,n0,m,i,rnorms)
c
        d=0
        do 1500 j=1,n
c
        d=d+rnorms(j)
 1500 continue
c
        if(d .lt. thresh) goto 4100
c
        ncols=i
c
c       find the row with the biggest norm and put it in the 
c       i-th position
c
        d=0
        do 1600 j=i,n
c     
        if(d .ge. rnorms(j)) goto 1600
        d=rnorms(j)
        ipiv=j
 1600 continue
c
        do 1800 j=1,m
c
        d=a(ipiv,j)
        a(ipiv,j)=a(i,j)
        a(i,j)=d
 1800 continue
c
        d=rnorms(ipiv)
        rnorms(ipiv)=rnorms(i)
        rnorms(i)=d
c
        j=ipivots(ipiv)
        ipivots(ipiv)=ipivots(i)
        ipivots(i)=j
c
c        put the largest element of the i-th row in the i-th place
c
        jpiv=i
        d=abs(a(i,i))
        do 1850 j=i+1,m
c
        dd=abs(a(i,j))
        if(dd .gt. d) then
            jpiv=j
            d=dd
        endif
c
 1850 continue
c
cccc        jpiv=i
cccc        call prinf('jpiv=*',jpiv,1)
c
        if(jpiv .eq. i) goto 1950
c
        do 1900 j=1,n
c
        d=a(j,i)
        a(j,i)=a(j,jpiv)
        a(j,jpiv)=d
 1900 continue
 1950 continue
c
c       eliminate the i-th row
c
        n0=1
        size22=thresh/1000
        do 2200 j=m,i+1,-1
c
        aa(1)=a(i,i)
        aa(2)=a(i,j)
c
        i0=i-1
        if(i .eq. 1) i0=1
        call jacskel_rotfnd(aa,uu,size22)
        call jacskel_rotate(uu,a(i0,i),a(i0,j),n,n0)

cccc        call prin2('a=*',a,n*m)

c
 2200 continue
c
 4000 continue
c
 4100 continue
c
c       return to the user the diagonal terms
c
        do 4200 i=1,ncols
c
        rnorms(i)=a(i,i)
 4200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine jacskel_rnorms(a,n,n0,m,m0,rnorms)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),rnorms(n)
c 
cccc        do 1400 i=1,n
        do 1400 i=n0,n
c 
        d=0
        do 1200 j=m0,m

        d=d+a(i,j)**2
 1200 continue
c 
        rnorms(i)=d
 1400 continue
c
        return
        end
c 
c 
c 
c 
c 
        subroutine jacskel_rotate(a,x,y,n,n0)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(2,2),x(1),y(1)
c 
        do 1200 i=n0,n
c 
        d1=a(1,1)*x(i)+a(1,2)*y(i)
        d2=a(2,1)*x(i)+a(2,2)*y(i)
c 
        x(i)=d1
        y(i)=d2
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine jacskel_rotfnd(a,u,size22)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(2),u(2,2)
c 
        u22=-a(1)
        u12=a(2)
c 
        d=u22**2+u12**2
c 
        if(d .lt. size22*1.0d-66) then
c 
            u(2,2)=1
            u(1,2)=0
            u(1,1)=1
            u(2,1)=0
            return
        endif
c 
        d=sqrt(d)
c 
        u(2,2)=u22/d
        u(1,2)=u12/d
        u(1,1)=-u(2,2)
        u(2,1)=u(1,2)
        return
        end
