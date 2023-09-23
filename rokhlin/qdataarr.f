        implicit real *16 (a-h,o-z)
        real *16 arr(10000)
c 
c 
        call prini(6,13)
c 
C 
C       SET ALL PARAMETERS
C 
C 
        PRINT *, 'ENTER n'
        READ *,n
        CALL PRINF('n=*',n,1 )
c 
c       construct the test array
c 
        do 1200 i=1,n
        arr(i)=i
 1200 continue
c 
c       test the writing routine
c 
        iw=31
        call qdataarr(n,iw,'arr',arr)
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine qdataarr(n,iw,name,arr)
        implicit real *16 (a-h,o-z)
        save
        dimension arr(1)
        character *3 name
        character *1 ch1(20)
        character *20 ch20
        equivalence (ch1(1),ch20)
        data ch20/'123456789a123456789a'/
c 
c       this subroutine prints out a real *16 user-supplied
c       array in a form of array declarations and data statements
c       that can be incorporated into a FORTRAN code. The DATA
c       statements will contain 32-digit number, and can be
c       used to construct extended precision codes.
c 
c                     input parameters:
c 
c  n - the number of real *16 elements to be incorporated
c       into the DATA statements. also the number of elements
c       in array arr.
c  iw - the fortran unit umber on which the output will be
c       printed
c  name - three characters that will be the prefix in the names
c       of all FORTRAN arrays that will be DATAed by the
c       subroutine.
c  arr - a real *16 array that will be stored in the data
c       statements
c 
c                     output parameters: none
c 
c       . . . write on disk the declarations of arrays
c             to be DATAed
c 
        ndata=n/16
        ntail=n-ndata*16
        nn=ndata
        if(ntail .ne. 0) nn=ndata+1
c 
        call arrequ(nn,name,iw,ntail)
c 
c       . . . write on disk the declarations of biffer arrays
c             and the equivalence statements
c 
        call arrbequ(nn-1,name,iw)
c 
c        . . . write on disk the data statements
c 
c 
c       . . . the main ones
c 
        ndata=n/16
        ntail=n-ndata*16
        if(ndata .eq. 0) goto 3000
        do 2000 i=1,ndata
c 
 1200 format(9x, 'data ',a3,i4,'/')
c 
        jj=1000+i
        write(iw,1200) name,jj
c 
 1400 format(5x,a1,6x,e38.32,',')
 1600 format(5x,a1,6x,e38.32,'/')
c 
        do 1800 j=1,16
        if(j .ne. 16) write(iw,1400) ch1(j), arr((i-1)*16+j)
        if(j .eq. 16) write(iw,1600) ch1(j), arr((i-1)*16+j)
 1800 continue
 2000 continue
 3000 continue
c 
c       write the data statements for the tail
c 
        if(ntail .eq. 0) goto 4000
c 
        jj=1000+ndata+1
        write(iw,1200) name,jj
        do 3200 j=1,ntail
        if(j .ne. ntail) write(iw,1400) ch1(j), arr(ndata*16+j)
        if(j .eq. ntail) write(iw,1600) ch1(j), arr(ndata*16+j)
 3200 continue
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
        subroutine arrequ(n,name,iw,nnt)
        implicit real *16 (a-h,o-z)
        save
        character *3 name
        character *1 ch1(20)
        character *20 ch20
        equivalence (ch1(1),ch20)
        data ch20/'123456789a123456789a'/
c 
c       write the array declarations
c 
        call prinf('in arrequ, n=*',n,1)
        call prinf('in arrequ, nnt=*',nnt,1)
 2200 format(9x, 'real *16')
c 
        ncards=n/4
        ntail=n-ncards*4
c 
        if(ncards .eq. 0) goto 3000
        j=0
        do 2600 i=1,ncards
        j=j+1
        if(j .eq. 11) j=1
c 
        if(j .eq. 1) write(iw,2200)
        jj=1000+i
        i16=16
        ilast=i16
 2400 format(5x,a1,2x,3(a3,i4,'(',i2,'),'),a3,i4,'(',i2,'),')
 2500 format(5x,a1,2x,3(a3,i4,'(',i2,'),'),a3,i4,'(',i2,')')
cccc 2500 format(5x,a1,2x,3(a3,i4,'(16),'),a3,i4,'(16)')
        i1=(i-1)*4+1+1000
        i2=i1+1
        i3=i2+1
        i4=i3+1
c 
        if( (i .ne. ncards) .and. (j .ne. 10) )
     1  write(iw,2400) ch1(j), name,i1,i16, name,i2,i16,
     2     name,i3,i16, name,i4,i16
c 
        ilast=16
        if( (i .eq. ncards) .and. (nnt .ne. 0) .and.
     1      (ntail .eq. 0) ) ilast=nnt
          call prinf('and ilast=*',ilast,1)
        if( (i .eq. ncards) .or. (j .eq. 10) )
     1  write(iw,2500) ch1(j), name,i1,i16, name,i2,i16,
     2     name,i3,i16, name,i4,ilast
 2600 continue
 3000 continue
c 
c      now, declare the tail
c 
        if(ntail .eq. 0) goto 4000
         write(iw,2200)
c 
 3200 format(5x,a1,2x,2(a3,i4,'(16),'),a3,i4,'(',i2,')')
 3400 format(5x,a1,2x,1(a3,i4,'(16),'),a3,i4,'(',i2,')')
 3600 format(5x,a1,2x,a3,i4,'(',i2,')')
c 
        ilast=16
        if(nnt .ne. 0) ilast=nnt
c 
        i1=ncards*4+1+1000
        i2=i1+1
        i3=i2+1
c 
        if(ntail .eq. 3)
     1  write(iw,3200) ch1(1),name,i1,name,i2,name,i3,ilast
c 
        if(ntail .eq. 2)
     1  write(iw,3400) ch1(1),name,i1,name,i2,ilast
c 
        if(ntail .eq. 1)
     1  write(iw,3600) ch1(1),name,i1,ilast
 4000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine arrbequ(n,name,iw)
        implicit real *16 (a-h,o-z)
        save
        character *3 name
        character *1 ch1(20)
        character *20 ch20
        equivalence (ch1(1),ch20)
        data ch20/'123456789a123456789a'/
c 
c       write the array declarations for the biffer arrays
c 
 2200 format(9x, 'real *16')
c 
        ncards=n/4
        ntail=n-ncards*4
c 
        if(ncards .eq. 0) goto 3000
        j=0
        do 2600 i=1,ncards
        j=j+1
        if(j .eq. 11) j=1
c 
        if(j .eq. 1) write(iw,2200)
        jj=1000+i
 2400 format(5x,a1,2x,3(a3,i4,'(2),'),a3,i4,'(2),')
 2500 format(5x,a1,2x,3(a3,i4,'(2),'),a3,i4,'(2)')
        i1=(i-1)*4+1+2000
        i2=i1+1
        i3=i2+1
        i4=i3+1
c 
        if( (i .ne. ncards) .and. (j .ne. 10) )
     1  write(iw,2400) ch1(j),name,i1,name,i2,name,i3,
     2      name,i4
c 
        if( (i .eq. ncards) .or. (j .eq. 10) )
     1  write(iw,2500) ch1(j),name,i1,name,i2,name,i3,
     2      name,i4
 2600 continue
 3000 continue
c 
c      now, declare the tail
c 
        if(ntail .eq. 0) goto 4000
         write(iw,2200)
c 
 3200 format(5x,a1,2x,2(a3,i4,'(2),'),a3,i4,'(2)')
 3400 format(5x,a1,2x,1(a3,i4,'(2),'),a3,i4,'(2)')
 3600 format(5x,a1,2x,a3,i4,'(2)')
c 
        i1=ncards*4+1+2000
        i2=i1+1
        i3=i2+1
c 
        if(ntail .eq. 3)
     1  write(iw,3200) ch1(1),name,i1,name,i2,name,i3
c 
        if(ntail .eq. 2)
     1  write(iw,3400) ch1(1),name,i1,name,i2
c 
        if(ntail .eq. 1)
     1  write(iw,3600) ch1(1),name,i1
 4000 continue
c 
c       write the equivalence statements, stitching the
c       small arrays into a long one
c 
 4200 format(9x,'equivalence')
cccc
c 
 4400 format(5x,a1,2x,'(',a3,i4,'(16),',a3,i4,'(1)),',
     1    '(',a3,i4,'(1),',a3,i4,'(2)),'  )
c 
 4500 format(5x,a1,2x,'(',a3,i4,'(16),',a3,i4,'(1)),',
     1    '(',a3,i4,'(1),',a3,i4,'(2))'  )
c 
        j=0
        do 5000 i=1,n
c 
        j=j+1
        if(j .eq. 11) j=1
c 
        if(j .eq. 1) write(iw,4200)
c 
        i1=i+1000
        i2=i1+1
        ii1=i1+1000
        ii2=i2+1000
c 
        if( (i .ne. n) .and. (j .ne. 10) )
     1         write(iw,4400) ch1(j),
     2      name,i1,name,ii1,name,i2,name,ii1
c 
        if( (i .eq. n) .or. (j .eq. 10) )
     1      write(iw,4500) ch1(j),
     2      name,i1,name,ii1,name,i2,name,ii1
c 
 5000 continue
 5200 continue
  
  
  
  
        return
        end
  
  
