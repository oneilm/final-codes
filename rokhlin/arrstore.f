c 
c 
c 
c 
c 
        subroutine arrstore(iw,arr,n,mess)
        implicit real *8 (a-h,o-z)
        save
        character *1 mess(1),ast
        dimension arr(2,1)
        data ast/'*'/
c 
  
c        This subroutine stores on the FORTRAN unit iw the
c        user-supplied real *8 array arr of length n in a
c        format easily converted into a set of FORTRAN data
c        statements. If the array arr is longish, it will be
c        broken up into pieces of length 180 elements (so
c        that each data statement is no longer than 90 cards)
c 
c                      Input parameters:
c 
c  iw - the FORTRAN unit number on which the data are to be stored
c  arr - the real *8 array to be stored
c  n - the number of elements in array arr
c  mess - a message to be written on the file iw before the stored
c        data. Must be ended with an asterisk.
c 
c                      Output parameters: none
c 
c 
        i2=0
        l=0
        do 1100 i=1,100
        if(mess(i) .eq. ast) goto 1150
        l=i
 1100 continue
 1150 continue
  
  
        write(iw,1250)
 1200 format(80a1)
        write(iw,1200) (mess(i),i=1,l)
 1250 format(20x)
        write(iw,1250)
 1270 format(20x,'total number of elements is',i5)
        write(iw,1270) n
        write(iw,1250)
c 
        n2=n/2
        ifodd=n-n2*2
c 
        nblocks=n2/90
  
  
        call prinf('nblocks=*',nblocks,1)
  
        nleft=n2-nblocks*90
c 
c       store on disk the array arr in blocks of 180 elements;
c       if there is a tail left that is smaller than 180 elements,
c       it will be handled later
c 
        if(nblocks .eq. 0) goto 2000
c 
        do 1800 i=1,nblocks
c 
 1300 format('c   ',/,'c        following is block number ',i3,
     1    ' of length 180 elements',/,'c    ')
c 
        write(iw,1300) i
c 
 1400 format(5x,i1,2x,e22.16,',',e22.16,',')
c 
        i1=(i-1)*90+1
        i2=i1+90-1
        icount=0
        do 1600 j=i1,i2
c 
        icount=icount+1
        if(icount .eq. 11) icount=1
c 
        write(iw,1400) icount,arr(1,j),arr(2,j)
 1600 continue
c 
 1800 continue
c 
 2000 continue
c 
c       if there is a tail - store it
c 
        if(nleft .eq. 0) goto 3000
c 
 2200 format('c   ',/,'c        following is the last block ',
     1    '(block number ',i3,')',/,'c        of length ',i3,
     2    ' elements',/,'c    ')
c 
        i=n-nblocks*180
        j=nblocks+1
        write(iw,2200) j,i
c 
        i1=nblocks*90+1
        i2=i2+nleft
  
        do 2600 j=i1,i2
c 
        icount=icount+1
        if(icount .eq. 11) icount=1
c 
        write(iw,1400) icount,arr(1,j),arr(2,j)
 2600 continue
 3000 continue
c 
c       if the number of elements in array arr is actually odd,
c       print the last element
c 
        icount=icount+1
        if(ifodd .eq. 1) write(iw,1400) icount,arr(1,n2+1)
  
 32000 format(20x)
        write(iw,3200)
  
        return
        end
  
