        IMPLICIT REAL *8 (A-H,O-Z)
        dimension x(10 000),y(10 000),x2(10 000),y2(10 000),
     1      x3(10 000),y3(10 000)
c 
        CALL PRINI(6,13)
  
         PRINT *, 'ENTEr n'
         READ *,n
         call prinf('n=*',n,1)
  
  
cccc        if(2 .ne. 3) goto 2000
  
c 
c       construct the test arrays to be plotted
c 
        done=1
        h=done/(n-1)
c 
        do 1200 i=1,n
        x(i)=(i-1)*h
        y(i)=x(i)
 1200 continue
c 
        iw=21
c 
        call lotagraph(iw,x,y,n,'from lotagraph*')
  
        done=1
        n2=n*2
        h=done/(n2-1)
c 
        do 1400 i=1,n2
        x2(i)=(i-1)*h
        y2(i)=x2(i)**2
 1400 continue
c 
        iw=22
c 
        call lotagraph2(iw,x,y,n,x2,y2,n2,'from lotagraph2*')
  
  
  
        done=1
        n3=n*3
        h=done/(n3-1)
c 
        do 1600 i=1,n3
        x3(i)=(i-1)*h
        y3(i)=x3(i)**3 *3
 1600 continue
c 
        iw=23
c 
        call lotagraph3(iw,x,y,n,x2,y2,n2,x3,y3,n3,
     1      'from lotagraph3*')
  
  
  
 2000 continue
  
c 
c       now, test the "plot" routines
c 
        done=1
        pi=datan(done)*4
c 
        h=2*pi/n
        do 2200 i=1,n
        x(i)=dcos((i-1)*h)
        y(i)=dsin((i-1)*h)
 2200 continue
  
        iw=31
        call lotaplot(iw,x,y,n,'from lotaplot*')
  
  
c 
        done=1
        pi=datan(done)*4
        n2=n*2
  
c 
        h=2*pi/n2
        do 2400 i=1,n2
        x2(i)=dcos((i-1)*h)
        y2(i)=2*dsin((i-1)*h)
 2400 continue
  
        iw=32
        call lotaplot2(iw,x,y,n,x2,y2,n2,'from lotaplot2*')
  
  
        done=1
        pi=datan(done)*4
        n3=n*3
  
c 
        h=2*pi/n3
        do 2600 i=1,n3
        x3(i)=dcos((i-1)*h)
        y3(i)=3*dsin((i-1)*h)
 2600 continue
  
        iw=33
        call lotaplot3(iw,x,y,n,x2,y2,n2,x3,y3,n3,
     1      'from lotaplot3*')
  
  
  
  
  
  
  
        STOP
        END
c 
c 
c 
c 
c 
        subroutine lotaplot(iw,x,y,n,title)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 title(1)
c 
        inumgr=1
c 
        call lotaplo0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
     1      title,inumgr)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine lotaplot2(iw,x,y,n,x2,y2,n2,title)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 title(1)
c 
        inumgr=2
c 
        call lotaplo0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
     1      title,inumgr)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine lotaplot3(iw,x,y,n,x2,y2,n2,x3,y3,
     1      n3,title)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 title(1)
c 
        inumgr=3
c 
        call lotaplo0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
     1      title,inumgr)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine lotagraph(iw,x,y,n,title)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 title(1)
c 
        inumgr=1
c 
        call lotagrap0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
     1      title,inumgr)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine lotagraph2(iw,x,y,n,x2,y2,n2,title)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 title(1)
c 
        inumgr=2
c 
        call lotagrap0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
     1      title,inumgr)
c 
        return
        end
c 
c 
c 
c 
c 
c 
        subroutine lotagraph3(iw,x,y,n,x2,y2,n2,x3,y3,
     1      n3,title)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 title(1)
c 
        inumgr=3
c 
        call lotagrap0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
     1      title,inumgr)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine lotagrap0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
     1      title,inumgr)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 title(1),gn(2),file1(8),anum1(8),line(82),
     1      blank,quo,backslash
        character *8 dummy,anum8,file8
c 
        equivalence (file1,file8),(anum1,anum8)
        data gn/'g','n'/,blank/' '/,quo/'"'/,backslash/'\\'/
c 
c        convert the user-specified Fortran unit number to
c        character format
c 
        if( (iw .ge. 0) .and. (iw .le. 9) )  write(dummy,1100) iw
 1100 format(i1)
c 
        if( (iw .ge. 10) .and. (iw .le. 99) ) write(dummy,1200) iw
 1200 format(i2)
c 
        if( (iw .ge. 100) .and. (iw .le. 999) ) write(dummy,1300) iw
 1300 format(i3)
c 
        if( (iw .ge. 1000) .and. (iw .le. 9999) ) write(dummy,1400) iw
 1400 format(i4)
c 
 2000 format(1a8)
        read(dummy,2000) anum8
c 
c        construct the file name on which the Gnuplot instructions
c        are to be written
c 
        file1(1)=gn(1)
        file1(2)=gn(2)
        do 2200 i=1,6
        file1(i+2)=anum1(i)
 2200 continue
c 
c        open the fortran file with the unit 87 and name file8
c 
        iun=87
        open(unit=iun,file=file8)
c 
 2250 format('   # set terminal postscript',/,
     1    '   # set output "plot.ps"')
c 
        write(iun,2250)
c 
c        generate the title for the plot
c 
        line(1)=blank
        line(2)=blank
        line(3)=quo
c 
        call MESslen(title,nchar,line(4))
c 
        line(nchar+4)=quo
 2300 format('  set title ',80a1)
        write(iun,2300) (line(i),i=1,nchar+4)
c 
 2350 format('   show title')
        write(iun,2350)
c 
c        write the instructions
c 
 2400 format(i6)
        iun2=iw+100000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c 
        do 2600 i=1,6
        file1(i+2)=anum1(i)
 2600 continue
c 
 2800 format('plot ','"',1a8,'"     ','notitle  with lines')
 2830 format('plot ','"',1a8,'"     ','notitle  with lines, ',1a1)
c 
        if(inumgr .eq. 1) write(iun,2800) file8
        if(inumgr .ne. 1) write(iun,2830) file8,backslash
c 
c        write the first data file to be plotted
c 
        iun22=88
        open(unit=iun22,file=file8)
c 
        write(iun22,3000) (x(i),y(i),i=1,n)
c 
        close(iun22)
c 
c       if the user so requested - write the instructions for the
c       plotting the second data file
c 
        if(inumgr .eq. 1) close(iun)
        if(inumgr .eq. 1) return
  
        iun2=iw+200000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c 
        do 2840 i=1,6
        file1(i+2)=anum1(i)
 2840 continue
c 
        if(inumgr .eq. 2) write(iun,2850) file8
        if(inumgr .eq. 3) write(iun,2855) file8,backslash
 2850 format('     ','"',1a8,'"     ','notitle  with lines')
 2855 format('     ','"',1a8,'"     ','notitle  with lines, ',1a1)
c 
c        write the second data file to be plotted
c 
        iun22=88
        open(unit=iun22,file=file8)
c 
        write(iun22,3000) (x2(i),y2(i),i=1,n2)
c 
        close(iun22)
c 
c       if the user so requested - write the instructions for the
c       plotting the third data file
c 
        if(inumgr .eq. 2) close(iun)
        if(inumgr .eq. 2) return
  
        iun2=iw+300000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c 
        do 2860 i=1,6
        file1(i+2)=anum1(i)
 2860 continue
c 
        write(iun,2870) file8
 2870 format('     ','"',1a8,'"     ','notitle  with lines')
c 
c        write the third data file to be plotted
c 
        iun22=88
        open(unit=iun22,file=file8)
c 
 3000 format(2x,e11.5,2x,e11.5)
c 
        write(iun22,3000) (x3(i),y3(i),i=1,n3)
c 
        close(iun22)
        close(iun)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine lotaplo0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
     1      title,inumgr)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 title(1),gn(2),file1(8),anum1(8),line(82),
     1      blank,quo,backslash
        character *8 dummy,anum8,file8
c 
        equivalence (file1,file8),(anum1,anum8)
        data gn/'g','n'/,blank/' '/,quo/'"'/,backslash/'\\'/
c 
c        convert the user-specified Fortran unit number to
c        character format
c 
        if( (iw .ge. 0) .and. (iw .le. 9) )  write(dummy,1100) iw
 1100 format(i1)
c 
        if( (iw .ge. 10) .and. (iw .le. 99) ) write(dummy,1200) iw
 1200 format(i2)
c 
        if( (iw .ge. 100) .and. (iw .le. 999) ) write(dummy,1300) iw
 1300 format(i3)
c 
        if( (iw .ge. 1000) .and. (iw .le. 9999) ) write(dummy,1400) iw
 1400 format(i4)
c 
 2000 format(1a8)
        read(dummy,2000) anum8
c 
c        construct the file name on which the Gnuplot instructions
c        are to be written
c 
        file1(1)=gn(1)
        file1(2)=gn(2)
        do 2200 i=1,6
        file1(i+2)=anum1(i)
 2200 continue
c 
c        open the fortran file with the unit 87 and name file8
c 
        iun=87
        open(unit=iun,file=file8)
c 
 2250 format('   # set terminal postscript',/,
     1    '   # set output "plot.ps"')
c 
        write(iun,2250)
c 
c        generate the title for the plot
c 
        line(1)=blank
        line(2)=blank
        line(3)=quo
c 
        call MESslen(title,nchar,line(4))
c 
        line(nchar+4)=quo
 2270 format('  set title ',80a1)
c 
        write(iun,2270) (line(i),i=1,nchar+4)
c 
 2280 format('   show title')
        write(iun,2280)
c 
c        find the limits for both x and y
c 
        xmin=1.0d30
        ymin=1.0d30
        xmax=-1.0d20
        ymax=-1.0d20
c 
        do 2300 i=1,n
        if(x(i) .lt. xmin) xmin=x(i)
        if(y(i) .lt. ymin) ymin=y(i)
        if(x(i) .gt. xmax) xmax=x(i)
        if(y(i) .gt. ymax) ymax=y(i)
 2300 continue
c 
        if(inumgr .eq. 1) goto 2340
c 
        do 2310 i=1,n2
        if(x2(i) .lt. xmin) xmin=x2(i)
        if(y2(i) .lt. ymin) ymin=y2(i)
        if(x2(i) .gt. xmax) xmax=x2(i)
        if(y2(i) .gt. ymax) ymax=y2(i)
 2310 continue
        if(inumgr .eq. 2) goto 2340
c 
        do 2320 i=1,n3
        if(x3(i) .lt. xmin) xmin=x3(i)
        if(y3(i) .lt. ymin) ymin=y3(i)
        if(x3(i) .gt. xmax) xmax=x3(i)
        if(y3(i) .gt. ymax) ymax=y3(i)
 2320 continue
c 
 2340 continue
c 
        xcenter=(xmin+xmax)/2
        ycenter=(ymax+ymin)/2
c 
        xsize=(xmax-xmin)
        ysize=(ymax-ymin)
        size=xsize
        if(ysize .gt. size) size=ysize
        size=size*1.1
c 
        xmin=xcenter-size/2
        xmax=xcenter+size/2
        ymin=ycenter-size/2
        ymax=ycenter+size/2
c 
c        set the size of the stupid thing
c 
 2350 format(2x,' set size 0.75,1.0')
        write(iun,2350)
c 
c        write the instructions
c 
 2400 format(i6)
        iun2=iw+100000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c 
        do 2600 i=1,6
        file1(i+2)=anum1(i)
 2600 continue
c 
 2800 format('plot ', '[',  e8.2, ':', e8.2,'] ',
     1    '[',  e8.2, ':', e8.2,'] ',
     2    '"',1a8,'" ','notitle with lines')
c 
 2830 format('plot ', '[',  e8.2, ':', e8.2,'] ',
     1    '[',  e8.2, ':', e8.2,'] ',
     2    '"',1a8,'" ','notitle with lines, ',1a1)
c 
        if(inumgr .eq. 1) write(iun,2800) xmin,xmax,ymin,ymax,file8
        if(inumgr .ne. 1) write(iun,2830) xmin,xmax,ymin,ymax,
     1      file8,backslash
c 
c        write the first data file to be plotted
c 
        iun22=88
        open(unit=iun22,file=file8)
c 
        write(iun22,3000) (x(i),y(i),i=1,n)
c 
        close(iun22)
c 
c       if the user so requested - write the instructions for the
c       plotting the second data file
c 
        if(inumgr .eq. 1) close(iun)
        if(inumgr .eq. 1) return
c 
        iun2=iw+200000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c 
        do 2840 i=1,6
        file1(i+2)=anum1(i)
 2840 continue
c 
 2850 format('     ', '"',1a8,'" ','notitle with lines')
c 
 2855 format('     ', '"',1a8,'" ','notitle with lines, ',1a1)
c 
         if(inumgr .eq. 2) write(iun,2850) file8
         if(inumgr .ne. 2) write(iun,2855) file8,backslash
c 
c        write the second data file to be plotted
c 
        iun22=88
        open(unit=iun22,file=file8)
c 
        write(iun22,3000) (x2(i),y2(i),i=1,n2)
c 
        close(iun22)
c 
c       if the user so requested - write the instructions for the
c       plotting the third data file
c 
        if(inumgr .eq. 2) close(iun)
        if(inumgr .eq. 2) return
  
        iun2=iw+300000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c 
        do 2860 i=1,6
        file1(i+2)=anum1(i)
 2860 continue
c 
        write(iun,2870) file8
 2870 format('     ','"',1a8,'"','notitle  with lines')
c 
c        write the third data file to be plotted
c 
        iun22=88
        open(unit=iun22,file=file8)
c 
 3000 format(2x,e11.5,2x,e11.5)
c 
        write(iun22,3000) (x3(i),y3(i),i=1,n3)
c 
        close(iun22)
        close(iun)
c 
        return
        end
c 
c 
c 
c 
c 
        SUBROUTINE MESslen(MES,nchar,line)
        save
        CHARACTER *1 MES(1),AST,line(1)
        DATA AST/'*'/
C 
C         DETERMINE THE LENGTH OF THE MESSAGE
C 
        I=0
        DO 1400 I=1,10000
        IF(MES(I).EQ.AST) GOTO 1600
        I1=I
 1400 CONTINUE
 1600 CONTINUE
c 
        nchar=i1
        do 1800 i=1,nchar
        line(i)=mes(i)
 1800 continue
         RETURN
         END
  
  
  
