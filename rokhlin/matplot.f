c        
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code and the beginning of the
c        plotting routines proper
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This file contains 6 user callable subroutines: matplot,
c        trplopen, trplpnt, trpldot, trplline, trplwrt, matplot. 
c        Folllowing is a brief description of these subroutines.
c
c   trplopen - initializes the array w to be used as storage by 
c        the subroutines trpldot, trplpnt, trplline, trplwrt (see), 
c        used for plotting data in GNUPLOT-readable form.
c
c   trplpnt - plots one point, with the coordinates (x,y)
c
c   trpldot - plots one dot, with the coordinates (x,y)
c
c   trplline - plots one line, whose ends have coordinates (x1,y1), 
c        (x2,y2). 
c
c   trplwrt - interprets the array w that has been created by prior 
c        calls to the subroutines trplopen, trpldot, trplpoint, 
c        trplline (see), and produces a collection of files readable 
c        by GNUPLOT.
c
c   matplot - plots the geometry of a mockup of a rectangular matrix
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
        subroutine matplot(a,n,m,eps,iw,title,ww)
        implicit real *8 (a-h,o-z)
        real *8 a(n,m)
        dimension ww(1)
        character *1 title(1)
c
c       this subroutine plots the geometry of a mockup of 
c       a rectangular matrix
c
c
c                input parameters:
c
c  a - the mockup of the matrix to be plotted. it is a real *8 
c       n \times m matrix; elements greater than eps will be plotted.
c  n,m - the dimensionalities of the matrix to be plotted
c  iw - the FORTRAN unit number on which the plot will be written
c  title - the title to be put on the plot
c
c                output paramaters: none
c
c                work arrays: 
c
c  ww - must be at least n**2 *5/2+20 real *8 elements 
c
c       . . . construct the array of coordinates 
c

        call trplopen(ww)
c
        num=0
        do 1400 i=1,m
        do 1200 j=1,n
        if(abs(a(j,i)) .lt. eps) goto 1200
        num=num+1
c
        xnum=i
        ynum=-j
        call trplpnt(xnum,ynum,ww)
c
 1200 continue
 1400 continue
c
c       construct the frame
c
        call trpl_plotframe(iw,n,m,ww)
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
        subroutine trpl_plotframe(iw,n,m,w)
        implicit real *8 (a-h,o-z)
        dimension w(1)
c
c        create the horizontal boundaries of blocks
c
        do 2000 i=1,n+1
        y1= - ((i-1)+0.5)
        y2=y1
        x1=0.5
        x2=m+0.5
c
        call trplline(x1,y1,x2,y2,w)
 2000 continue
c
c        create the vertical boundaries of blocks
c
        do 3000 i=1,m+1
        x1=  ((i-1)+0.5)
        x2=x1
        y1=-0.5
        y2=-(n+0.5)
c
        call trplline(x1,y1,x2,y2,w)
 3000 continue
        return
        end
c
c
c
c 
c 
        subroutine trplopen(w)
        implicit real *8 (a-h,o-z)
        save
        real *4 w(3,1)
c 
c       this subroutine initializes the array w to be used as storage
c       by the subroutines trpldot, trplpnt, trplline,
c       trplwrt (see) , used for plotting data in GNUPLOT-readable form.
c       The array w must be sufficiently long, since each plotting
c       operation (dot, point, line) takes 5 real *4 locations to
c       memorize
c 
        w(1,1)=1
        return
        end
c 
c 
c 
c 
c 
        subroutine trpldot(x,y,w)
        implicit real *8 (a-h,o-z)
        save
        real *4 w(5,1)
c 
c       This subroutine plots one dot, with the coordinates (x,y).
c       Actually, the information is stored in array w, to be
c       plotted later by a call to the subroutine trplwrt (see)
c 
        w(1,1)=w(1,1)+1
c 
        i=w(1,1)+0.1
        w(1,i)=1
        w(2,i)=x
        w(3,i)=y
        return
        end
c 
c 
c 
c 
c 
        subroutine trplpnt(x,y,w)
        implicit real *8 (a-h,o-z)
        save
        real *4 w(5,1)
c 
c       This subroutine plots one point, with the coordinates (x,y).
c       Actually, the information is stored in array w, to be
c       plotted later by a call to the subroutine trplwrt (see)
c 
        w(1,1)=w(1,1)+1
c 
        i=w(1,1)+0.1
        w(1,i)=2
        w(2,i)=x
        w(3,i)=y
        return
        end
c 
c 
c 
c 
c 
        subroutine trplline(x1,y1,x2,y2,w)
        implicit real *8 (a-h,o-z)
        save
        real *4 w(5,1)
c 
c       This subroutine plots one line, whose ends have
c       coordinates (x1,y1), (x2,y2). Actually, the information
c       is stored in array w, to be plotted later by a call to
c       the subroutine trplwrt (see)
c 
        w(1,1)=w(1,1)+1
c 
        i=w(1,1)+0.1
        w(1,i)=3
        w(2,i)=x1
        w(3,i)=y1
        w(4,i)=x2
        w(5,i)=y2
        return
        end
c 
c 
c 
c 
c 
        subroutine trplwrt(iw,w,title)
        implicit real *4 (a-h,o-z)
        save
        real *4 w(5,1)
        character *1 title(1),gn(2),file1(8),anum1(8),line(82),
     1      blank,quo,backslash
        character *8 dummy,anum8,file8,file81,file82
c 
        equivalence (file1,file8),(anum1,anum8)
        data gn/'g','n'/,blank/' '/,quo/'"'/,backslash/'\\'/
c 
c        This subroutine interprets the array w that has been
c        created by prior calls to the subroutines trplopen,
c        trpldot, trplpoint, trplline (see), and produces a
c        collection of files readable by GNUPLOT. The subroutine
c        creates from 2 to 4 files, depending on which of the
c        subroutines trpldot, trplpoint, trplline have been
c        previously called. The names of the files this subroutine
c        will created is determined (to some extent) by the value
c        of the parameter iw. EXPLANATION: Suppose that iw=21.
c        Then the main GNUPLOT-readable file will be named gn21;
c        the auxiliary file gn100021 will contain data describing
c        dots to be plotted; the auxiliary file gn200021 will
c        contain data describing points to be plotted;
c        dots to be plotted; and the auxiliary file gn300021 will
c        contain data describing lines to be plotted.
  
c 
c                Input parameters:
c 
c  iw - integer parameter determining the names of the files to be
c        produced
c  w - the array produced by prior calls to the subroutines trplopen,
c        trpldot, trplpoint, trplline (see).
c  title - the title of the plot. MUST end with the asterisk, which
c        is not plotted
c 
c        . . . convert the user-specified Fortran unit number to
c              character format
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
        call trpl_messlen(title,nchar,line(4))
c 
        line(nchar+4)=quo
 2270 format('  set title ',80a1)
c 
        write(iun,2270) (line(i),i=1,nchar+4)
c 
 2280 format('   show title')
        write(iun,2280)
c 
c        reorder the array w, putting the instructions to write
c        dots first, instructions to write points second, and
c        instructions to write lines last
c 
        call trplreor(w,idots,ndots,ipoints,npoints,
     1      ilines,nlines)
c 
c        find the limits for both x and y
c 
        xmin=1.0d30
        ymin=1.0d30
        xmax=-1.0d20
        ymax=-1.0d20
c 
        n=ndots+npoints+nlines
c
        do 2300 i=2,n+1
c 
        if(w(2,i) .gt. xmax) xmax=w(2,i)
        if(w(2,i) .lt. xmin) xmin=w(2,i)
        if(w(3,i) .gt. ymax) ymax=w(3,i)
        if(w(3,i) .lt. ymin) ymin=w(3,i)
c 
        j=w(1,i)+0.1
        if(j .ne. 3) goto 2300
c 
        if(w(4,i) .gt. xmax) xmax=w(4,i)
        if(w(4,i) .lt. xmin) xmin=w(4,i)
        if(w(5,i) .gt. ymax) ymax=w(5,i)
        if(w(5,i) .lt. ymin) ymin=w(5,i)
 2300 continue
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
 2350 format(2x,' set size 1.0,1.0')
        write(iun,2350)
c 
c        write the instructions and the data file in the case
c        when only the dots are to be plotted
c 
c        . . . instructions
c 
        if(ndots .eq. 0) goto 3000
        if(npoints .ne. 0) goto 3000
        if(nlines .ne. 0) goto 3000
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
 2800 format('plot ', '[',  e11.5, ':', e11.5,'] ',
     1    '[',  e11.5, ':', e11.5,'] ',
     2    '"',1a8,'" ',1a1,/,6x,'notitle with dots')
c 
        write(iun,2800) xmin,xmax,ymin,ymax,file8
c 
c        . . . the data file
c 
        iun22=88
        open(unit=iun22,file=file8)
c 
 2900 format(2x,e11.5,2x,e11.5)
        write(iun22,2900) (w(2,i),w(3,i),i=idots,idots+ndots-1)
c 
        close(iun22)
        return
  
 3000 continue
c 
c        write the instructions and the data file in the case
c        when only the points are to be plotted
c 
c        . . . instructions
c 
        if(npoints .eq. 0) goto 4000
        if(ndots .ne. 0) goto 4000
        if(nlines .ne. 0) goto 4000
c  
 3400 format(i6)
        iun2=iw+200000
        write (dummy,3400) iun2
        read(dummy,2000) anum8
c 
        do 3600 i=1,6
        file1(i+2)=anum1(i)
 3600 continue
c 
 3800 format('plot ', '[',  e11.5, ':', e11.5,'] ',
     1    '[',  e11.5, ':', e11.5,'] ',
     2    '"',1a8,'" ',1a1,/,6x,'notitle with points')
c 
        write(iun,3800) xmin,xmax,ymin,ymax,file8
c 
c        . . . the data file
c 
        iun22=88
        open(unit=iun22,file=file8)
c 
 3900 format(2x,e11.5,2x,e11.5)
        write(iun22,3900) (w(2,i),w(3,i),i=ipoints,ipoints+npoints-1)
c 
        close(iun22)
        return
 4000 continue
c 
c        write the instructions and the data file in the case
c        when only the lines are to be plotted
c 
c        . . . instructions
c 
        if(nlines .eq. 0) goto 5000
        if(ndots .ne. 0) goto 5000
        if(npoints .ne. 0) goto 5000
c  
 4400 format(i6)
        iun2=iw+300000
        write (dummy,4400) iun2
        read(dummy,2000) anum8
c 
        do 4600 i=1,6
        file1(i+2)=anum1(i)
 4600 continue
c 
 4800 format('plot ', '[',  e11.5, ':', e11.5,'] ',
     1    '[',  e11.5, ':', e11.5,'] ',
     2    '"',1a8,'" ',1a1,/,6x,'notitle with lines')
c 
        write(iun,4800) xmin,xmax,ymin,ymax,file8,backslash
c 
c        . . . the data file
c 
        iun22=88
        open(unit=iun22,file=file8)
c 
 4900 format(2x,e11.5,2x,e11.5)
c 
        do 4950 i=1,nlines
c  
        j=i+ilines-1
        write(iun22,4900) w(2,j),w(3,j)
        write(iun22,4900) w(4,j),w(5,j)
c 
 4940 format(20x)
        write(iun22,4940)
c 
 4950 continue
c 
        close(iun22)
        return
 5000 continue
c 
c        write the instructions and the data files in the case
c        when points and dots are to be plotted
c 
c        . . . instructions
c 
        if(ndots .eq. 0) goto 6000
        if(npoints .eq. 0) goto 6000
        if(nlines .ne. 0) goto 6000
c 
 5400 format(i6)
        iun2=iw+100000
        write (dummy,5400) iun2
        read(dummy,2000) anum8
c 
        do 5600 i=1,6
        file1(i+2)=anum1(i)
 5600 continue
c 
        file81=file8
c 
 5800 format('plot ', '[',  e11.5, ':', e11.5,'] ',
     1    '[',  e11.5, ':', e11.5,'] ',
     2    '"',1a8,'" ',1a1,/,6x,'notitle with dots, ',1a1)
c 
        write(iun,5800) xmin,xmax,ymin,ymax,file8,backslash,
     1      backslash
c 
        iun2=iw+200000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c 
        do 5820 i=1,6
        file1(i+2)=anum1(i)
 5820 continue
c 
 5850 format('     ', '"',1a8,'" ','notitle with points')
        write(iun,5850) file8
c 
c        . . . the data files
c 
        iun22=88
        open(unit=iun22,file=file8)
c 
 5900 format(2x,e11.5,2x,e11.5)
c 
        iun22=88
        open(unit=iun22,file=file8)
c 
        write(iun22,5900) (w(2,i),w(3,i),i=ipoints,ipoints+npoints-1)
c 
        close(iun22)
c 
        iun22=88
        open(unit=iun22,file=file81)
c 
        write(iun22,5900) (w(2,i),w(3,i),i=idots,idots+ndots-1)
        close(iun22)
 6000 continue
c 
c        write the instructions and the data files in the case
c        when dots and lines are to be plotted
c 
c        . . . instructions
c 
        if(ndots .eq. 0) goto 7000
        if(nlines .eq. 0) goto 7000
        if(npoints .ne. 0) goto 7000
c 
 6400 format(i6)
        iun2=iw+100000
        write (dummy,5400) iun2
        read(dummy,2000) anum8
c 
        do 6600 i=1,6
        file1(i+2)=anum1(i)
 6600 continue
c 
        file81=file8
c 
 6800 format('plot ', '[',  e11.5, ':', e11.5,'] ',
     1    '[',  e11.5, ':', e11.5,'] ',
     2    '"',1a8,'" ',1a1,/,'notitle with dots, ',1a1)
c 
        write(iun,6800) xmin,xmax,ymin,ymax,file8,
     1          backslash,backslash
c 
        iun2=iw+300000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c 
        do 6820 i=1,6
        file1(i+2)=anum1(i)
 6820 continue
c 
 6850 format('     ', '"',1a8,'" ','notitle with lines')
        write(iun,6850) file8
c 
c        . . . the data files
c 
        iun22=88
        open(unit=iun22,file=file81)
c 
        write(iun22,6900) (w(2,i),w(3,i),i=idots,idots+ndots-1)
c 
        close(iun22)
c 
 6900 format(2x,e11.5,2x,e11.5)
c 
        iun22=88
        open(unit=iun22,file=file8)
c 
        do 6950 i=1,nlines
c  
        j=i+ilines-1
        write(iun22,4900) w(2,j),w(3,j)
        write(iun22,4900) w(4,j),w(5,j)
c 
 6940 format(20x)
        write(iun22,6940)
c 
 6950 continue
c  
        close(iun22)
 7000 continue
c 
c        write the instructions and the data files in the case
c        when points and lines are to be plotted
c 
c        . . . instructions
c 
        if(npoints .eq. 0) goto 8000
        if(nlines .eq. 0) goto 8000
        if(ndots .ne. 0) goto 8000
c 
 7400 format(i6)
        iun2=iw+200000
        write (dummy,7400) iun2
        read(dummy,2000) anum8
c 
        do 7600 i=1,6
        file1(i+2)=anum1(i)
 7600 continue
c 
        file81=file8
c 
 7800 format('plot ', '[',  e11.5, ':', e11.5,'] ',
     1    '[',  e11.5, ':', e11.5,'] ',
     2    '"',1a8,'" ',1a1,/,6x,'notitle with points, ',1a1)
c 
        write(iun,7800) xmin,xmax,ymin,ymax,file8,backslash,
     1      backslash
c 
        iun2=iw+300000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c 
        do 7820 i=1,6
        file1(i+2)=anum1(i)
 7820 continue
c 
 7850 format('     ', '"',1a8,'" ','notitle with lines')
        write(iun,7850) file8
c 
c        . . . the data files
c 
        iun22=88
        open(unit=iun22,file=file81)
c 
        write(iun22,7900) (w(2,i),w(3,i),i=ipoints,ipoints+npoints-1)
c 
        close(iun22)
 7900 format(2x,e11.5,2x,e11.5)
c 
        iun22=88
        open(unit=iun22,file=file8)
c 
        do 7950 i=1,nlines
  
        j=i+ilines-1
        write(iun22,7900) w(2,j),w(3,j)
        write(iun22,7900) w(4,j),w(5,j)
c 
 7940 format(20x)
        write(iun22,7940)
c 
 7950 continue
c  
        close(iun22)
 8000 continue
c 
c        write the instructions and the data files in the case
c        when dots, points, and lines are to be plotted
c 
c        . . . instructions
c 
        if(npoints .eq. 0) goto 9000
        if(nlines .eq. 0) goto 9000
        if(ndots .eq. 0) goto 9000
c 
 8400 format(i6)
        iun2=iw+100000
        write (dummy,8400) iun2
        read(dummy,2000) anum8
c 
        do 8600 i=1,6
        file1(i+2)=anum1(i)
 8600 continue
c 
        file81=file8
c 
 8800 format('plot ', '[',  e11.5, ':', e11.5,'] ',
     1    '[',  e11.5, ':', e11.5,'] ',
     2    '"',1a8,'" ',1a1,/,6x,'notitle with dots, ',1a1)
c 
        write(iun,8800) xmin,xmax,ymin,ymax,file8,backslash,
     1      backslash
c 
        iun2=iw+200000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c 
        do 8820 i=1,6
        file1(i+2)=anum1(i)
 8820 continue
c 
        file82=file8
c 
 8850 format('     ', '"',1a8,'" ','notitle with points, ',1a1)
        write(iun,8850) file8,backslash
c 
        iun2=iw+300000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c 
        do 8860 i=1,6
        file1(i+2)=anum1(i)
 8860 continue
c 
 8870 format('     ', '"',1a8,'" ','notitle with lines')
        write(iun,8870) file8
c 
c        . . . the data files
c 
        iun22=88
        open(unit=iun22,file=file81)
c 
        write(iun22,8900) (w(2,i),w(3,i),i=idots,idots+ndots-1)
c 
        close(iun22)
c 
        iun22=88
        open(unit=iun22,file=file82)
c 
        write(iun22,8900) (w(2,i),w(3,i),i=ipoints,ipoints+npoints-1)
c 
        close(iun22)
c 
 8900 format(2x,e11.5,2x,e11.5)
c 
        iun22=88
        open(unit=iun22,file=file8)
c 
        do 8950 i=1,nlines
        j=i+ilines-1
        write(iun22,8900) w(2,j),w(3,j)
        write(iun22,8900) w(4,j),w(5,j)
c 
 8940 format(20x)
        write(iun22,8940)
c 
 8950 continue
        close(iun22)
 9000 continue
        return
        end
c 
c 
c 
c 
c 
        SUBROUTINE trpl_messlen(MES,nchar,line)
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
c 
c 
c 
c 
c 
        subroutine trplreor(w,idots,ndots,ipoints,npoints,
     1      ilines,nlines)
c 
        save
        real *4 w(5,1)
c 
c        reorder the array w, putting the instructions to write
c        dots first, instructions to write points second, and
c        instructions to write lines last
c 
c       . . . dots
c 
        n=w(1,1)+0.1
c 
        i1=2
        n1=1
        ndots=0
        do 1400 i=i1,n
c 
        j=w(1,i)+0.1
        if(j .ne. 1) goto 1400
        n1=n1+1
        ndots=ndots+1
c 
        call trplswtc(w(1,i),w(1,n1),5)
 1400 continue
c 
c       . . . points
c 
        i2=n1+1
        n2=i2-1
        npoints=0
        do 1800 i=i2,n
c 
        j=w(1,i)+0.1
        if(j .ne. 2) goto 1800
        n2=n2+1
        npoints=npoints+1
c 
        call trplswtc(w(1,i),w(1,n2),5)
 1800 continue
c 
        i3=n2+1
        nlines=n-ndots-npoints-1
c 
        idots=2
        ipoints=idots+ndots
        ilines=ipoints+npoints
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine trplswtc(x,y,n)
        implicit real *4 (a-h,o-z)
        save
        dimension x(1),y(1)
c 
        do 1200 i=1,n
        d=y(i)
        y(i)=x(i)
        x(i)=d
 1200 continue
        return
        end
