        implicit real *8 (a-h,o-z)
        real *8 a(100 000),b(100 000),c(100 000),
     1    w(410 000)
c 
        external elcomp2
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER n'
        READ *,n
        CALL PRINF('n=*',n,1 )
c 
        PRINT *, 'ENTER m'
        READ *,m
        CALL PRINF('m=*',m,1 )
  
  
c 
c       generate the random vectors to be sorted
c 
        call RIRAND(N*m,a)
  
cccc        call prin2('a as created*',a,n*2)
        call prin2('a as created*',a,10)
  
c 
c       test the merging routine
c 
        call mergetest(a,m,n,b,c,elcomp2,w)
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine mergetest(a,m,n,b,c,elcomp,w)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(m,n),b(m,n),c(m,n),w(1)
        external elcomp
c 
c 
  
c 
c       this subroutine tests the merging routine
c 
c        . . . copy array a into arrays b,c
c 
  
        call sortrcpy(a,b,n*m)
        call sortrcpy(a,c,n*m)
c 
c       merge-sort the array a
c 
        call sortrear(a,m,n,elcomp,w)
cccc        call prin2('after sortrear, a=*',a,n*m)
        call prin2('after sortrear, a=*',a,10)
c 
c 
        call sortrbub(b,m,n,elcomp,w)
        call prin2('after bubble, b=*',b,10)
c 
        ifbad=0
        do 2400 i=1,n
        do 2200 j=1,m
c 
        c(j,i)=a(j,i)-b(j,i)
c 
        if(c(j,i) .ne. 0) ifbad=1
 2200 continue
 2400 continue
ccccc        call prinf('and difference is*',c,n)
         call prinf('and ifbad=*',ifbad,1)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine elcomp2(a,b,m,ifagtb)
        implicit real *8 (a-h,o-z)
        save
        integer a(2),b(2)
c 
        ifagtb=0
        if(a(1) .eq. b(1)) goto 1200
c 
        if(a(1) .gt. b(1)) ifagtb=1
        return
c 
 1200 continue
c 
        if(a(2) .gt. b(2) ) ifagtb=1
  
  
cccc        call prinf('exiting elcomp, ifabtb=*',ifagtb,1)
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
        subroutine sortrear(a,m,n,elcomp,w)
        implicit real *8 (a-h,o-z)
c 
        save
        integer *4 ls(30),starts(30),lengths(30),
     1      lenlog(30)
c 
        real *8 a(m,n),w(m,n)
        external elcomp
c 
        data nold/0/
c 
c        this subroutine sorts the columns of the user-provided
c        real two-dimensional array a. The sorting is performed with
c        respect to the user-provided comparizon subroutine elcomp.
c 
c               input parameters:
c 
c  a - the array to be sorted
c  m - the number of rows in a
c  n - the number of columns in a
c  elcomp - the subroutine performing the comparizon between two
c       columns of a. The calling sequence of elcomp is
c 
c               elcomp(a,b,m,ifagtb).                               (1)
c 
c              The input parameters in (1) above are:
c 
c    a,b - the two arrays to be compared
c    m   - the length of the arrays a,b
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c              The output parameters in (1) above are:
c 
c     ifagtb - an integer parameter;
c 
c           if a < b, ifagtb must be set to 0
c           if a > b, ifagtb must be set to 1
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c               output parameters:
c 
c  a - the input array after sorting
c 
c               work arrays:
c 
c  w - must be at least n*m+1  real *8 elements long
c 
c        . . . if n is less than 16 - use the bubble algorithm
c 
        if(n .ge. 16) goto 1100
        call sortrbub(a,m,n,elcomp,w)
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
        call sortrep2(a,m,ki-1,elcomp,w)
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
c 
        if(lla .eq. 1) goto 2000
        if(lla .gt. 10)  call sortrep2(a(1,ia),m,la,elcomp,w)
        if(lla .le. 10)  call sortrbub(a(1,ia),m,lla,elcomp,w)
 2000 continue
c 
c        perform the final merge
c 
        ia=1
        na=lengths(1)
        do 2400 i=2,nsort
        ib=starts(i)
        nb=lengths(i)
        call sortmep2(a,m,na,a(1,ib),nb,elcomp,w)
c 
        do 2200 j=1,na+nb
c 
        call sortrcpy(w(1,j),a(1,j),m)
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
        subroutine sortrep2(a,m,k,elcomp,w)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(m,1),w(m,1)
        external elcomp
c 
c       sort a real array of  length 2**k
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
c 
        call elcomp(a(1,i),a(1,i+1),m,ifagtb)
c 
        if(ifagtb .le. 0) goto 1200
c 
        call sortrcpy(a(1,i),w,m)
        call sortrcpy(a(1,i+1),a(1,i),m)
        call sortrcpy(w,a(1,i+1),m)
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
        call sortmep2(a(1,j1),m,length,a(1,j2),length,elcomp,w)
c 
        do 1400 l=1,length*2
c 
        call sortrcpy(w(1,l),a(1,j1+l-1),m)
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
        subroutine sortmep2(a,m,na,b,nb,elcomp,c)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(m,1),b(m,1),c(m,1)
  
c 
c       merge the sorted arrays ia, ib obtaining the sorted array ic
c 
        ia=1
        ib=1
        ic=1
c 
cccc        call prinf('in sortmep2 before loop, a=*',a,na*2)
cccc        call prinf('in sortmep2 before loop, b=*',b,nb*2)
  
  
        do 1600 i=1,na+nb
        ic=i
c 
cccc        call prinf('in sortmep2 before elcomp, i=*',i,1)
  
  
        call elcomp(a(1,ia),b(1,ib),m,ifagtb)
  
cccc        call prinf('in sortmep2 after ifagtb=*',ifagtb,1)
  
        if(ifagtb .gt. 0) goto 1200
  
cccc        if(a(ia) .gt. b(ib) ) goto 1200
  
  
        call sortrcpy(a(1,ia),c(1,i),m)
  
  
cccc        c(1,i)=a(1,ia)
cccc        c(2,i)=a(2,ia)
  
        ia=ia+1
        goto 1400
c 
 1200 continue
  
        call sortrcpy(b(1,ib),c(1,i),m)
  
cccc        c(1,i)=b(1,ib)
cccc        c(2,i)=b(2,ib)
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
cccc        c(1,ic+i)=a(1,ia+i-1)
cccc        c(2,ic+i)=a(2,ia+i-1)
  
        call sortrcpy(a(1,ia+i-1),c(1,ic+i),m)
  
 2000 continue
        return
c 
 2200 continue
c 
c       we have run out of elements in a. finish off the array b
c 
        do 2400 i=1,nb-ib+1
  
cccc        c(1,ic+i)=b(1,ib+i-1)
cccc        c(2,ic+i)=b(2,ib+i-1)
  
        call sortrcpy(b(1,ib+i-1),c(1,ic+i),m)
  
 2400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine sortrbub(a,m,n,elcomp,w)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(m,1),w(1)
        external elcomp
  
cccc    call prin2('in bubble, a=*',a,n*m)
  
c 
        do 1400 i=1,n
        do 1200 j=1,n-i
c 
        call elcomp(a(1,j),a(1,j+1),m,ifagtb)
c 
        if(ifagtb .eq. 0) goto 1200
  
        call sortrcpy(a(1,j),w,m)
        call sortrcpy(a(1,j+1),a(1,j),m)
        call sortrcpy(w,a(1,j+1),m)
 1200 continue
 1400 continue
  
cccc        call prin2('exiting sortrbub, a=*',a,n*m)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine sortrcpy(a,b,n)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension a(1),b(1)
c 
        do 1200 i=1,n
        b(i)=a(i)
 1200 continue
        return
        end
