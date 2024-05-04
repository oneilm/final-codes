c 
c 
c (c) 2007 Vladimir Rokhlin from his final codes
c 
c 
        subroutine rsortanyn_vec(a, ind, lda, n, w)
        implicit real *8 (a-h,o-z)
        save
        double precision  a(*), w(*)

cccc        data nold/0/
c 
c        This subroutine sorts the COLUMNS of the real matrix a of
c        dimension lda X n according to the ROW ind; it is an exact real
c        version of the subroutine sortanyn.
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
c        . . . if n is less than 16 - use the rbubble algorithm
c 
        integer ls(30),starts(30),lengths(30)
        integer lenlog(30)


        if(n .ge. 16) goto 1100
        call rbubble(a,n)
        return
 1100 continue
c 
c       construct the representation of n in base 2
c 
cccc        if(n .eq. nold) goto 1900

        j1=n
        do 1200 i=1,100
           nls=i
           j2=j1/2
           l=j1-j2*2
           ls(i)=l
           j1=j2
           if(j1 .eq. 0) goto 1400
 1200   continue
 1400   continue
        
cccc        nold=n
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
        call rsortpwr2(a,ki-1,w)
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
cccc        if(la .ne. 1) call rsortpwr2(a(ia),la,w)
cccc        call prinf('during preliminary sort, a=*',a,n)
        if(lla .eq. 1) goto 2000
        if(lla .gt. 10)  call rsortpwr2(a(ia),la,w)
        if(lla .lt. 10)  call rbubble(a(ia),lla)
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
        call rsortmerg(a,na,a(ib),nb,w)
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
c
      subroutine rsortanyn(a,n,w)
        implicit real *8 (a-h,o-z)
        save
        integer *4 ls(30),starts(30),lengths(30),
     1      lenlog(30)
c 
        real *8 a(1),w(1)
  
  
        data nold/0/
c 
c        this subroutine sorts the real array a of length n;
c        it is an exact real version of the subroutine sortanyn
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
c        . . . if n is less than 16 - use the rbubble algorithm
c 
        if(n .ge. 16) goto 1100
        call rbubble(a,n)
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
        call rsortpwr2(a,ki-1,w)
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
cccc        if(la .ne. 1) call rsortpwr2(a(ia),la,w)
cccc        call prinf('during preliminary sort, a=*',a,n)
        if(lla .eq. 1) goto 2000
        if(lla .gt. 10)  call rsortpwr2(a(ia),la,w)
        if(lla .lt. 10)  call rbubble(a(ia),lla)
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
        call rsortmerg(a,na,a(ib),nb,w)
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
        subroutine rsortpwr2(a,k,w)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(1),w(1)
c 
c       sort a real array of  length 2**k
c 
c           work array:
c 
c w - must be at least 2**(k-1)+1 integer *4 elements long
c 
        n=2**k
c 
c       perform rbubble sort on the pairs of points
c 
        do 1200 i=1,n-1,2
        if(a(i) .le. a(i+1) ) goto 1200
        rj=a(i)
        a(i)=a(i+1)
        a(i+1)=rj
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
        call rsortmerg(a(j1),length,a(j2),length,w)
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
        subroutine rsortmerg(a,na,b,nb,c)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(1),b(1),c(1)
c 
c       merge the sorted arrays a, b obtaining the sorted array c
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
        subroutine rbubble(a,n)
        implicit real *8 (a-h,o-z)
        real *8 a(n)
c 
        do 1400 i=1,n
           do 1200 j=1,n-i
              if(a(j+1) .ge. a(j) ) goto 1200
              rk=a(j)

              a(j)=a(j+1)
              a(j+1)=rk
 1200      continue
 1400   continue
        return
        end
