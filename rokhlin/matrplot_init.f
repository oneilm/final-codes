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

        xmin=-1
        xmax=1
        ymin=-1
        ymax=1
 
        iw=21  
        call matrplot_init(iw,xmin,xmax,ymin,ymax)

        x1=-1
        x2=1
        y=-1

        call matrplot_hor(iw,x1,x2,y)
  


        x1=-1
        x2=1
        y=1

        call matrplot_hor(iw,x1,x2,y)
  


        x=-1
        y1=-1
        y2=1

        call matrplot_vert(iw,x,y1,y2)
  
        x=1
        y1=-1
        y2=1

        call matrplot_vert(iw,x,y1,y2)
  
        x=0
        y=0
        call matrplot_text(iw,x,y,
     1      '$see IT again$ \ \ and again and again*')
cccc        call matrplot_text(iw,x,y,'$ IT$*')

        x=0.3
        y=0.3
        width=0.1 *2
        hight=0.1 *2
        
        call matrplot_dark(iw,x,y,width,hight)


        call prin2('done!*',x,0)
        call matrplot_close(iw)

        x1=10
        x2=20
        y=2

        x1=0
        x2=10
        y=0  
  
c 
        stop
        end
c 
c 
c 
c 
c 
        subroutine matrplot_init(iw,xmin,xmax,ymin,ymax)
        implicit real *8 (a-h,o-z)
        save
        character *1 text(1),ast,beg_line1(9),card(100),
     1      comma,brac(3),brac1
        character *10 beg_line
        character *80 card2
c
        equivalence(beg_line1(1),beg_line)
c
        data ast/'*'/,comma/','/,brac/')',' ','{'/,brac1/'}'/
        data beg_line/'   \\put ( '/

c 
c       this entry initializes the construction of 
c       LATEX-readable file plotting a collection of vertical 
c       and horizontal intervals, pieces of text, and 
c       solid rectangles
c 
c                       input parameters:
c 
c  iw - the Fortran file on which the output will be written
c  xmin - the expected minimum value of x 
c  xmax - the expected maximum value of x 
c  ymin - the expected minimum value of y 
c  ymax - the expected maximum value of y 
c
c                       output parameters: none
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
cccc     b  '  \\setlength{\\unitlength}{2pt}         ')
c 
        write(iw,1200)
c 
 1300 format('  \\begin{picture}(' F20.10,',',F20.10,') (0,0)   '/)
c 
        dx=24
        dy=24

        write(iw,1300) dx, dy
c 
c        determine the linear mapping that will convert the
c        user-supplied values of x, y into coordinates
c        on the screen
c 
        alpha=24/(xmax-xmin)
        beta=-alpha*xmin
c
        gamma=24/(ymax-ymin)
        delta=-ymin*gamma
c 
        return
c
c
c
c
        entry matrplot_hor(iw,x1,x2,y)
c 
c       this entry plots a horizontal line
c 
c                       input parameters:
c 
c  iw - the Fortran file on which the output will be written
c  x1 - the x-coordinate of the first end of the interval to be plotted
c  x2 - the x-coordinate of the second end of the interval to be plotted
c  y - the y-coordinate of all points on the interval to be plotted
c
 2200 format('   \\put (',f11.5, ',',f11.5,
     1    ') { \\line (1,0) {',f11.5,' } }')
c
        xx1=alpha*x1+beta
        xx2=alpha*x2+beta
        yy=gamma*y+delta
c
        d=xx2-xx1
        write(iw,2200) xx1,yy,d
c
        return
c
c
c
c
        entry matrplot_vert(iw,x,y1,y2)
c 
c       this entry plots a vertical line
c 
c                       input parameters:
c 
c  iw - the Fortran file on which the output will be written
c  x - the x-coordinate of all points on the interval to be plotted
c  y1 - the y-coordinate of the first end of the interval to be plotted
c  y2 - the y-coordinate of the second end of the interval to be plotted
c
c
c        determine the parameters of the interval to be plotted
c
 2400 format('   \\put (',f11.5, ',',f11.5,
     1    ') { \\line (0,1) {',f11.5,' } }')
c
        yy1=gamma*y1+delta
        yy2=gamma*y2+delta
        xx=alpha*x+beta
c
        d=yy2-yy1
        write(iw,2400) xx,yy1,d
c
        return
c
c
c
c
        entry matrplot_text(iw,x,y,text)
c 
c       this entry prints a user-supplied text
c 
c                       input parameters:
c 
c  iw - the Fortran file on which the output will be written
c  (x,y) - the x-coordinates of the left lower corner of the text
c  text - the text to be printed; must be shorter than 100 
c       characters and end with the asterisk
c
c        . . . determine the length of the text
c
        lentext=0
        do 2600 i=1,100
c
        if(text(i) .eq. ast) then
            lentext=i-1
            goto 2700
        endif
c
 2600 continue
 2700 continue
c
c        build the line to be written
c
        yy=gamma*y+delta
        xx=alpha*x+beta
c
        do 2750 i=1,10
c
        card(i)=beg_line1(i)
 2750 continue
c
 2770 format(f11.5)
 2790 format(11a1)
        write(card2,2770) xx
        read(card2,2790) (card(i+10),i=1,11)
c
        card(22)=comma
c
        write(card2,2770) yy
        read(card2,2790) (card(i+22),i=1,11)
c
        card(34)=brac(1)
        card(35)=brac(2)
        card(36)=brac(3)
c
        do 2800 i=1,lentext
c
        card(36+i)=text(i)
 2800 continue
c
        card(36+lentext+1)=brac1
 2820 format(100a1)
c
        write(iw,2820) (card(i),i=1,36+lentext+1)
c
        return
c
c
c
c
        entry matrplot_dark(iw,x,y,width,hight)
c 
c       this entry prints a black rectangle 
c 
c                       input parameters:
c 
c  iw - the Fortran file on which the output will be written
c  (x,y) - the x-coordinates of the left lower corner of the 
c       rectangle
c  width - the width of the rectangle
c  hight - the hight of the rectangle
c
        yy=gamma*y+delta
        xx=alpha*x+beta

        hhight=gamma*hight
        wwidth=alpha*width

       call prin2('width=*',width,1)
       call prin2('wwidth=*',wwidth,1)
       call prin2('alpha=*',alpha,1)


 3000 format('   \\put (',f11.5, ',',f11.5,
     1    ')',/,'   { \\rule {', f11.5,' \\unitlength}',
     2    '{',f11.5,' \\unitlength } }')
c
        write(iw,3000) xx,yy,wwidth,hhight
c
        return
c
c
c
c
        entry matrplot_close(iw)
c 
c       this entry closes the file iw on which a LATEX picture has 
c       been written
c
 3200 format(' \\end{picture}')
 3400 format(' \\end{document}')

        write(iw,3200) 
        write(iw,3400) 

        return
        end











