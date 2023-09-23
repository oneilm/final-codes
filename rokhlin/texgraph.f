        implicit real *8 (a-h,o-z)
        dimension f(1000),t(1000)
c 
        CALL PRINI(6,13)
C 
c       . . . construct the testing parameters
c 
c 
         PRINT *, 'ENTER n'
         READ *,n
         CALL PRINf('n=*',n,1)
c 
c       construct the test curve to be plotted
c 
        done=1
        h=done/n
        h=h*2
        do 1200 i=1,n
        t(i)=i
        f(i)=t(i)**2/n *2
 1200 continue
c 
        t(4)=1.0d21
c 
c 
        ifpad=1
        iw=21
        call texgraph(iw,t,f,n,ifpad,
     1      'graphing in the new packaging*')
  
  
        iw=23
        call texplot(iw,t,f,n,ifpad,
     1      'plotting in the new packaging*')
  
  
  
  
  
c 
        stop
        end
c 
c 
c 
c 
c 
        subroutine texgraph(iw,x,y,n,ifpad,mes)
        implicit real *8 (a-h,o-z)
        save
        character *1 mes(1)
        dimension x(1),y(1)
c 
c       this subroutine creates a LATEX-readable file
c       plotting a user-specified curve. It does not force
c       the x and y scales to be equal, unlike its entry
c       texplot (see below).
c 
c                       input parameters:
c 
c  iw - the Fortran file on which the output will be written
c  x - the values of the x-coordinate
c  y - the values of the y-coordinate
c  n - number of elements in each of the arrays t, f
c  ifpad - tells the subroutine how whether to pad the spaces
c           between the user-specified points. ifpadit=0 will
c           cause the subroutine to put dots at the points
c           (t(i),f(i)). ifpadit=1 will  cause the subroutine to
c           connect consecutive points with straight lines.
c  mes - the title of the plot. must end with * (asterisk), which
c           is not plotted.
c 
c                      output parameters: none
c 
c       . . . graph the curve specified by the nodes (x(i), y(i))
c 
        scalex=1
        scaley=1
        ifequal=0
c 
        call texgrph(iw,x,y,n,scalex,scaley,
     1      ifequal,mes,ifpad)
c 
        return
c 
c 
c 
c 
        entry texplot(iw,x,y,n,ifpad,mes)
c 
c       this subroutine creates a LATEX-readable file
c       plotting a user-specified curve. It forces
c       the x and y scales to be equal, unlike its entry
c       texgraph (see above).
c 
c                       input parameters:
c 
c  iw - the Fortran file on which the output will be written
c  x - the values of the x-coordinate
c  y - the values of the y-coordinate
c  n - number of elements in each of the arrays t, f
c  ifpad - tells the subroutine how whether to pad the spaces
c           between the user-specified points. ifpadit=0 will
c           cause the subroutine to put dots at the points
c           (t(i),f(i)). ifpadit=1 will  cause the subroutine to
c           connect consecutive points with straight lines.
c  mes - the title of the plot. must end with * (asterisk), which
c           is not plotted.
c 
c                      output parameters: none
c 
c       . . . plot the curve specified by the nodes (x(i), y(1))
c 
        scalex=1
        scaley=1
        ifequal=1
c 
        call texgrph(iw,x,y,n,scalex,scaley,
     1      ifequal,mes,ifpad)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine texgrph(iw,t,f,n,scalex,scaley,
     1      ifequal,mes,ifpadit)
        implicit real *8 (a-h,o-z)
        save
        dimension t(1),f(1),rmarks(20),numbers(20)
        character *1 mes(1)
c 
c       this subroutine creates a LATEX-readable file
c       plotting a user-specified curve.
c 
c                       input parameters:
c 
c  iw - the Fortran file on which the output will be written
c  t - the values of the x-coordinate
c  f - the values of the y-coordinate
c  n - number of elements in each of the arrays t, f
c  scalex - the coefficient by which the x-axis is compressed.
c           it has only two legal values: 1 and 0.5.
c  scaley - the coefficient by which the y-axis is compressed.
c           it has only two legal values: 1 and 0.5.
c  ifequal - the parameter telling the subroutine whether the
c           scales in t, f will be equal. if ifequal=1, the
c           equality will be enforced. Otherwise, the subroutine
c           will feel free to define at will the maximum and
c            munimum values to be plotted.
c  mes - the title of the plot. must end with * (asterisk), which
c           is not plotted.
c  ifpadit - tells the subroutine how whether to pad the spaces
c           between the user-specified points. ifpadit=0 will
c           cause the subroutine to put dots at the points
c           (t(i),f(i)). ifpadit=1 will  cause the subroutine to
c           connect consecutive points with straight lines.
c 
c 
c                      output parameters: none
c 
c           . . . construct the header of the file
c 
 1200 format(
     1  '  \\documentstyle [fleqn,11pt]{article}    ',/,
     2  '  \\textwidth 450 pt                       ',/,
     3  '  \\textheight 575 pt                      ',/,
     4  '  \\topmargin 0 pt                         ',/,
     5  '  \\oddsidemargin 0 pt                     ',/,
     6  '  \\evensidemargin 0 pt                    ',/,
     7  '  \\mathindent 72pt                        ',/,
     8  '  \\pagestyle{empty}                       ',/,
     9  '  \\begin{document}                        ',/,
     a  '                                          ',/,
     b  '  \\setlength{\\unitlength}{5.0mm}         ')
c 
        write(iw,1200)
c 
 1300 format('  \\begin{picture}(' F20.10,',',F20.10,') (-2,0)   '/)
c 
        dx=24*scalex
        dy=24*scaley
        write(iw,1300) dx, dy
c 
c       find the minimum and maximum of x, y
c 
        tmin=1.0d20
        tmax=-1.0d20
        fmin=1.0d20
        fmax=-1.0d20
c 
        do 1400 i=1,n
c 
        if(dabs(t(i)) .gt. 1.0d20) goto 1400
c 
        if(tmin. gt. t(i)) tmin=t(i)
        if(tmax. lt. t(i)) tmax=t(i)
c 
        if(fmin. gt. f(i)) fmin=f(i)
        if(fmax. lt. f(i)) fmax=f(i)
c 
 1400 continue
c 
        call prinf('ifequal=*',ifequal,1)
        if(ifequal .eq. 0) goto 1600
        d1=tmin
        if(fmin .lt. d1) d1=fmin
c 
        d2=tmax
        if(fmax .gt. d2) d2=fmax
c 
        tmin=d1
        fmin=d1
        tmax=d2
        fmax=d2
c 
 1600 continue
c 
c        determine the linear mapping that will convert the
c        user-supplied values of t, f into coordinates
c        on the screen
c 
        dx=tmax-tmin
        dy=fmax-fmin
c 
        tleft=tmin-dx/10
        tright=tmax+dx/10
c 
        fbottom=fmin-dy/10
        ftop=fmax+dy/10
c 
c 
        xleft=0
        xright=24*scalex
        ybottom=0
        ytop=24*scaley
c 
        alpha=(xright-xleft)/(tright-tleft)
        beta=xleft-alpha*tleft
c 
        gamma=(ytop-ybottom)/(ftop-fbottom)
        delta=ybottom-gamma*fbottom
  
 2000 format(' %',/,
     1    ' % plotting the user-specified function',/,' %')
        write(iw,2000)
c 
c       plot the user-specified function
c 
        do 3000 i=1,n
c 
        if(dabs(t(i)) .gt. 1.0d20) goto 3000
c 
 2200 format('   \\put (',f11.5, ',',f11.5,
     1    ') { \\line (1,0)   {0.025 } }')
c 
        x=alpha*t(i)+beta
        y=gamma*f(i)+delta
c 
c       if the need be - pad this interval
c 
        ifpad=0
        if(i .eq. 1) goto 2800
        if(dabs(t(i-1)) .gt. 1.0d20) goto 2800
        if(ifpadit .ne. 0)
     1      call padit1(iw,xold,yold,x,y,ifpad)
 2800 continue
c 
        if(ifpad .eq. 0) write(iw,2200) x,y
c 
        xold=x
        yold=y
 3000 continue
c 
c        create the coordinate system
c 
 3100 format(' %',/,
     1    ' % plotting the system of coordinates',/,' %')
        write(iw,3100)
c 
cccc 3200 format('   \\put (1,1) { \\vector (1,0)   {25 } }',/,
 3200 format('   \\put (1,1) { \\vector (1,0)   {',f20.10,' } }',/,
cccc     1       '    \\put (1,1) { \\vector (0,1)   {25 } } ')
     1     '  \\put (1,1) { \\vector (0,1)  {', f20.10,' } } ')
c 
        dx=25*scalex
        dy=25*scaley
        write(iw,3200) dx, dy
c 
c       annotate the x axis
c 
 3250 format(' %',/,
     1    ' % annotating the x axis',/,' %')
        write(iw,3250)
c 
        call fndmarks(ier,tmin,tmax,nmarks,rmarks,numbers)
         call prinf('after first fndmarks, numbers=*',
     1       numbers,nmarks)
c 
        do 3500 i=1,nmarks
c 
 3300 format('   \\put (',f11.5, ',',f11.5,
     1    ') { \\line (0,1)   {0.2 } }')
c 
 3350 format('   \\put (',f11.5, ',',f11.5,
     1    ') { \\line (0,1)   {0.4 } }')
c 
        x=alpha*rmarks(i)+beta
        y=0.9
        if(numbers(i) .eq. 0) write(iw,3300) x,y
        if(numbers(i) .ne. 0) write(iw,3350) x,y
c 
 3400 format('   \\put (',f11.5, ',',f11.5,
     1    ') {', e8.2,'}')
c 
 3450 format('   \\put (',f11.5, ',',f11.5,
     1    ') {\\tiny ', e8.2,'}')
c 
        x=x-1
        y=-0.5
c 
        if( (numbers(i) .ne. 0) .and. (scalex .gt. 0.75) )
     1      write(iw,3400) x, y, rmarks(i)
c 
        if( (numbers(i) .ne. 0) .and. (scalex .lt. 0.75) )
     1      write(iw,3450) x, y, rmarks(i)
c 
 3500 continue
c 
c      annotate the y axis
c 
 3600 format(' %',/,
     1    ' % annotating the y axis',/,' %')
        write(iw,3600)
c 
        call fndmarks(ier,fmin,fmax,nmarks,rmarks,numbers)
         call prinf('after second fndmarks, numbers=*',
     1       numbers,nmarks)
c 
        do 4500 i=1,nmarks
c 
 4300 format('   \\put (',f11.5, ',',f11.5,
     1    ') { \\line (1,0)   {0.2 } }')
c 
 4350 format('   \\put (',f11.5, ',',f11.5,
     1    ') { \\line (1,0)   {0.4 } }')
c 
        y=gamma*rmarks(i)+delta
        x=0.9
cccc        x=1
cccc        x=0.95
        if(numbers(i) .eq. 0) write(iw,4300) x,y
        if(numbers(i) .ne. 0) write(iw,4350) x,y
c 
 4400 format('   \\put (',f11.5, ',',f11.5,
     1    ') {', e8.2,'}')
c 
 4450 format('   \\put (',f11.5, ',',f11.5,
     1    ') {\\tiny ', e8.2,'}')
c 
        x=x-4
        if(scalex .lt. 0.75) x=x+1
c 
        if( (numbers(i) .ne. 0) .and. (scalex .gt. 0.75) )
     1      write(iw,4400) x, y, rmarks(i)
c 
        if( (numbers(i) .ne. 0) .and. (scalex .lt. 0.75) )
     1      write(iw,4450) x, y, rmarks(i)
c 
 4500 continue
c 
c       create the title
c 
 6250 format(' %',/,
     1    ' % creating the title',/,' %')
        write(iw,6250)
c 
 6280 format('  \\put(1,-4) { ')
        if(scalex .gt. 0.75) write(iw,6280)
 6290 format('  \\put(1,-2.5) { ')
        if(scalex .lt. 0.75) write(iw,6290)
c 
 6300 format('  \\begin{minipage}[t]{',f10.5,' mm}')
        d=120
        if(scalex .lt. 0.75) d=d/2
        write(iw,6300) d
c 
 6400 format('  {  \\begin{center}')
c 
 6500 format('  {\\small \\begin{center}')
       if(scalex .lt. 0.75) write(iw,6500)
c 
       if(scalex .gt. 0.75) write(iw,6400)
c 
        call titletex(MES,Iw)
c 
 6600 format('   \\end{center}   }')
        write(iw,6600)
c 
 6800 format('   \\end{minipage}   }')
        write(iw,6800)
c 
c       construct the tail of the file
c 
 7200 format('   \\end{picture}        ')
        write(iw,7200)
c 
 7800 format('   \\end{document}             ')
        write(iw,7800)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine padit1(iw,xold,yold,x,y,ifpad)
        implicit real *8 (a-h,o-z)
c 
c       find the distance between this dot and the
c       preceding one, in mm
c 
        save
        thresh=0.1
c 
        d=(x-xold)**2+(y-yold)**2
        d=dsqrt(d)
        d=d*5
cccc        call prin2('inside padit1, d=*',d,1)
c 
        ifpad=0
        if(d .gt. thresh) ifpad=1
        if(ifpad .eq. 0) return
cccc        call prinf('inside padit1, ifpad=*',ifpad,1)
c 
c        check how much to pad this interval
c 
        n=d/thresh
        n=n+1
c 
cccc        call prinf('inside padit1, n=*',n,1)
        hx=(x-xold)/n
        hy=(y-yold)/n
c 
c       . . . pad
c 
 1200 format('   \\multiput (',f8.5, ',',f8.5,
     1    ')(',f8.5, ',',f8.5, '){',i5,
     2    '}{\\line(1,0){0.025}}')
c 
cccc 1200 format('   \\multiput (',f20.10, ',',f20.10,
cccc     1    ')(',f20.10, ',',f20.10, '){',i5,
cccc     2    '}{\\line(1,0){0.025}}')
c 
cccc        write(iw,1200) xold+hx,yold+hy,hx,hy,n
        write(iw,1200) xold,yold,hx,hy,n+1
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine fndmarks(ier,a,b,nmarks,rmarks,numbers)
        implicit real *8 (a-h,o-z)
        save
        dimension rmarks(1),numbers(1)
c 
c       find the rough scale on which marks will be set
c 
        ier=0
        c=b-a
        ten=10
        do 1200 i=-20,20
        d1=ten**i
        d2=d1*ten
        if(( c .ge. d1) .and. (c. le. d2) ) goto 1400
 1200 continue
        ier=4
        return
 1400 continue
         call prin2('d1 as constructed*',d1,1)
c 
c       determine the number of marks and their spacing
c 
        n=c/d1
         call prinf('n as constructed*',n,1)
c 
        if(n .eq. 1) step=d1/10
        if(n .eq. 2) step=d1/5
        if(n .eq. 3) step=d1/2
        if(n .eq. 4) step=d1/2
        if(n .eq. 5) step=d1/2
        if(n .eq. 6) step=d1
        if(n .eq. 7) step=d1
        if(n .eq. 8) step=d1
        if(n .eq. 9) step=d1
        if(n .eq. 10) step=d1
c 
        if(c/step .gt. 13) step=step*2
c 
         call prin2('step as constructed*',step,1)
c 
c       construct the locations of the marks
c 
        j=0
        t=a+step
        k=t/step
        t=k*step
        do 1600 i=1,20
        j=j+1
        t=t+step
        rmarks(j)=t
c 
        if(t .gt. b) goto 1800
 1600 continue
 1800 continue
c 
        nmarks=j
        if(rmarks(j) .gt. b) nmarks=nmarks-1
c 
        call prinf('nmarks=*',nmarks,1)
c 
c       construct the array controllong the printing of numbers
c 
        do 2200 i=1,nmarks
        numbers(i)=0
 2200 continue
c 
        if(nmarks .ne. 5) goto 3105
        numbers(1)=1
        numbers(3)=1
        numbers(5)=1
 3105 continue
c 
        if(nmarks .ne. 6) goto 3106
        numbers(2)=1
        numbers(4)=1
        numbers(6)=1
 3106 continue
c 
        if(nmarks .ne. 7) goto 3107
        numbers(2)=1
        numbers(4)=1
        numbers(6)=1
 3107 continue
c 
        if(nmarks .ne. 8) goto 3108
        numbers(2)=1
        numbers(4)=1
        numbers(6)=1
        numbers(8)=1
 3108 continue
c 
        if(nmarks .ne. 9) goto 3109
        numbers(2)=1
        numbers(5)=1
        numbers(8)=1
 3109 continue
c 
        if(nmarks .ne. 10) goto 3110
        numbers(1)=1
        numbers(4)=1
        numbers(7)=1
        numbers(10)=1
 3110 continue
c 
        if(nmarks .ne. 11) goto 3111
        numbers(2)=1
        numbers(5)=1
        numbers(8)=1
        numbers(11)=1
 3111 continue
c 
        if(nmarks .ne. 12) goto 3112
        numbers(2)=1
        numbers(5)=1
        numbers(8)=1
        numbers(11)=1
 3112 continue
c 
        if(nmarks .ne. 13) goto 3113
        numbers(3)=1
        numbers(6)=1
        numbers(9)=1
        numbers(12)=1
 3113   continue
c 
         call prinf('numbers as constructed*',numbers,nmarks)
        return
        end
c 
c 
c 
c 
c 
        SUBROUTINE titletex(MES,IP)
        save
        CHARACTER *1 MES(1),AST,card(80),blank,title(5),quo
        DATA AST/'*'/,quo/' '/,blank/' '/
     1      title/' ',' ',' ',' ',' '/
C 
C         DETERMINE THE LENGTH OF THE MESSAGE
C 
        I=0
        DO 1400 I=1,10000
        IF(MES(I).EQ.AST) GOTO 1600
        I1=I
 1400 CONTINUE
 1600 CONTINUE
 1800 FORMAT(1X, 80A1)
         do 2000 i=1,80
         card(i)=blank
 2000 continue
         do 2200 i=1,5
         card(i+1)=title(i)
 2200 continue
         card(8)=quo
c 
         do 2400 i=1,i1
         card(i+8)=mes(i)
 2400 continue
         card(8+i1+1)=quo
c 
         write(ip,1800) card
         return
         end
  
