        integer *4 a(10000),b(10000),c(10000),
     1    w(21000)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
         PRINT *, 'ENTER n'
        READ *,n
        CALL PRINF('n=*',n,1 )
c 
c       test the merging routine
c 
         if(2 .ne.3) goto 2100
        na=10
        nb=14
        do 1200 i=1,nb
        a(i)=i*3
        b(i)=i*2
 1200 continue
        call prinf('a=*',a,na)
        call prinf('b=*',b,nb)
c 
c       . . . merge
c 
        call sortmerg(a,na,b,nb,c)
        call prinf('after merge, c=*',c,na+nb)
c 
c        now, test the sorting routine on arrays of length
c        2**k
c 
        k=7
        n=2**k
        do 2000 i=1,n
        a(i)=-i*(-1)**i
 2000 continue
        call prinf('a before sortpwr2,*',a,n)
c 
        call sortpwr2(a,k,w)
         call prinf('a after sortwowr2*',a,n)
  
        n=1+4+8+32+128
        n=128
 2100 continue
c 
        do 2200 i=1,n
        a(i)=-i*(-1)**i
        b(i)=a(i)
        c(i)=a(i)
 2200 continue
  
cccc        call sortanyn(c,n,w)
        call sortanyn(a,n,w)
        call prinf('after sortanyn, a=*',a,10)
c 
        call bubble(b,n)
        call prinf('after bubble, b=*',b,10)
c 
        ifbad=0
        do 2400 i=1,n
        c(i)=a(i)-b(i)
         if(c(i) .ne. 0) ifbad=1
 2400 continue
ccccc        call prinf('and difference is*',c,n)
         call prinf('and ifbad=*',ifbad,1)
  
        stop
        end
c 
c 
c 
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c            this is the end of the debugging code and the
c            beginning of the actual sorting routine
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c 
c 
        subroutine sortanyn(a,n,w)
        save
        integer *4 ls(30),a(1),w(1),starts(30),lengths(30),
     1      lenlog(30)
        data nold/0/
c 
c        this subroutine sorts the integer array a of length n
c 
c               input parameters:
c 
c  a - the array to be sorted
c  n - the length of the array a
c 
c               output parameters:
c 
c  a - the input array after sorting
c 
c               work arrays:
c 
c  w - must be at least n+1 integer *4 elements long
c 
c        . . . if n is less than 16 - use the bubble algorithm
c 
        if(n .ge. 16) goto 1100
        call bubble(a,n)
        return
 1100 continue
c 
c       construct the representation of n in base 2
c 
        if(n .eq. nold) goto 1900
        j1=n
        do 1200 i=1,100
        nls=i
        j2=j1/2
        l=j1-j2*2
        ls(i)=l
        j1=j2
        if(j1 .eq. 0) goto 1400
 1200 continue
 1400 continue
        nold=n
cccc        call prinf('binary representation of n is*',ls,nls)
c 
c       construct the map for consecutive sorts and merges
c 
        nsort=0
        starts(1)=1
        do 1800 i=1,nls
        if(ls(i) .eq. 0) goto 1800
        ki=i
        nsort=nsort+1
        lengths(nsort)=2**(i-1)
        starts(nsort+1)=starts(nsort)+lengths(nsort)
        lenlog(nsort)=i-1
 1800 continue
 1900 continue
cccc        call prinf('lengths=*',lengths,nsort)
cccc        call prinf('lenlog=*',lenlog,nsort)
cccc        call prinf('starts=*',starts,nsort)
c 
c       if n is a power of two - do a simple shell sort
c 
        if(nsort .ne. 1) goto 1950
        call sortpwr2(a,ki-1,w)
        return
 1950 continue
c 
c       sort each of the pieces of length 2**k
c 
        do 2000 i=1,nsort
c 
cccc        call prinf('i=*',i,1)
        ia=starts(i)
  
        la=lenlog(i)
        lla=lengths(i)
cccc        call prinf('ia=*',ia,1)
cccc        call prinf('la=*',la,1)
c 
cccc        if(la .ne. 1) call sortpwr2(a(ia),la,w)
cccc        call prinf('during preliminary sort, a=*',a,n)
        if(lla .eq. 1) goto 2000
        if(lla .gt. 10)  call sortpwr2(a(ia),la,w)
        if(lla .lt. 10)  call bubble(a(ia),lla)
 2000 continue
cccc        call prinf('after preliminary sort, a=*',a,n)
c 
c        perform the final merge
c 
        ia=1
        na=lengths(1)
        do 2400 i=2,nsort
        ib=starts(i)
        nb=lengths(i)
        call sortmerg(a,na,a(ib),nb,w)
c 
        do 2200 j=1,na+nb
        a(j)=w(j)
 2200 continue
c 
        na=na+nb
 2400 continue
cccc        call prinf('after final sort, a=*',a,n)
        return
        end
c 
c 
c 
c 
c 
        subroutine sortpwr2(a,k,w)
        save
        integer *4 a(1),w(1)
c 
c       sort an integer array of  length 2**k
c 
c           work array:
c 
c w - must be at least 2**(k-1)+1 integer *4 elements long
c 
        n=2**k
c 
c       perform bubble sort on the pairs of points
c 
        do 1200 i=1,n-1,2
        if(a(i) .le. a(i+1) ) goto 1200
        j=a(i)
        a(i)=a(i+1)
        a(i+1)=j
 1200 continue
c 
c       perform merging recursively
c 
        do 2000 i=1,k-1
        length=2**i
        nint=n/length
        do 1600 j=1,nint/2
        j1=(j-1)*length*2+1
        j2=j1+length
        call sortmerg(a(j1),length,a(j2),length,w)
c 
        do 1400 l=1,length*2
        a(j1+l-1)=w(l)
 1400 continue
 1600 continue
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine sortmerg(a,na,b,nb,c)
        save
        integer *4 a(1),b(1),c(1)
c 
c       merge the sorted arrays ia, ib obtaining the sorted array ic
c 
        ia=1
        ib=1
        ic=1
c 
        do 1600 i=1,na+nb
        ic=i
c 
        if(a(ia) .gt. b(ib) ) goto 1200
        c(i)=a(ia)
        ia=ia+1
        goto 1400
c 
 1200 continue
        c(i)=b(ib)
        ib=ib+1
 1400 continue
c 
        if(ib .gt. nb) goto 1800
        if(ia .gt. na) goto 2200
 1600 continue
c 
 1800 continue
c 
c       we have run out of elemnets in b. finish off the array a
c 
        do 2000 i=1,na-ia+1
        c(ic+i)=a(ia+i-1)
 2000 continue
        return
c 
 2200 continue
c 
c       we have run out of elements in a. finish off the array b
c 
        do 2400 i=1,nb-ib+1
        c(ic+i)=b(ib+i-1)
 2400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine bubble(a,n)
        save
        integer *4 a(n)
c 
        do 1400 i=1,n
        do 1200 j=1,n-i
        if(a(j+1) .ge. a(j) ) goto 1200
        k=a(j)
        a(j)=a(j+1)
        a(j+1)=k
 1200 continue
 1400 continue
        return
        end
