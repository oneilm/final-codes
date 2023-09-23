        implicit real *8 (a-h,o-z)
c 
        call prini(6,13)
  
  
c 
c       the following file will run the executable code
c       "savemfort" on one or more fortran files.
c 
  
C        while
C        test -n "$1"
C        do
  
C          cp $1 t9654236
C          savemfort
C          cp t9654237 $1
C          rm t9654236
C          rm t9654237
  
C          shift
C        done
  
  
  
  
  
        call savem('t9654236','t9654237')
  
  
  
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine savem(name1,name2)
        implicit real *8 (a-h,o-z)
c 
        save
        character *1 card1(80),com,blank,com2
        character *8 name1,name2
        character *12 card_save
c 
c        read the FORTRAN file "name1"; find subroutines and
c        functions that have not been declared "save",
c        declare them things "save", and store the obtained
c        (much improved) code in the file "name2"
c 
        data com/'c'/,blank/' '/,com2/'C'/
c 
        data card_save/'        save'/
C 
c       open the input and output files
c 
        open (unit=11,file=name1)
        open (unit=12,file=name2)
c 
c       read cards one after another
c 
        ifsub=0
        ifsave=0
        do 2000 i=1,10000000
c 
        do 1100 j=1,80
        card1(j)=blank
 1100 continue
c 
 1200 format(80a1)
 1400 format(a12)
c 
c        . . . read the card
c 
  
        read(11,1200,end=2200) card1
c 
c       if this is a "save" card - skip it
c 
        call findword(card1,'save','SAVE',4,ifsame)
c 
        if(ifsame .eq. 1) goto 2000
c 
c       if this is a "comment"  card - copy it without further ado
c 
        if(card1(1) .eq. com ) goto 1600
        if(card1(1) .eq. com2) goto 1600
c 
c       if this is a "subroutine"  card - prepare to insert
c       the "save" card
c 
        call findword(card1,'subroutine','SUBROUTINE',10,ifsame)
        call findword(card1,'function','FUNCTION',8,ifsame2)
c 
        if(ifsame2 .eq. 1) ifsame=1
c 
        if(ifsame .eq. 1) then
            ifsub=1
            goto 1600
        endif
c 
cccc        if( (ifsub .eq. 1) .and. (card1(6) .ne. blank))
cccc     1      goto 1600
c 
  
c 
        if(ifsub .eq. 1) then
            if(card1(6) .ne. blank) goto 1600
c 
c       . . . if this is an "implicit" card - delay
c             dumping the "save" statement
c 
            call findword(card1,'implicit','IMPLICIT',8,ifsame3)
            if(ifsame3 .eq. 1) goto 1600
        endif
c 
        if( (ifsub .eq. 1) .and. (card1(6) .eq. blank)) then
        ifsave=1
        ifsub=0
            ifsub=0
        endif
c 
 1600 continue
c 
c        If there are trailing blanks on this card
c        - eliminate them things
c 
        nj=2
        do 1800 j=3,80
c 
        if(card1(j) .ne. blank) nj=j
 1800 continue
c 
        if(ifsave .eq. 1) then
            write(12,1400) card_save
            ifsave=0
        endif
c 
        write(12,1200) (card1(j),j=1,nj)
 2000 continue
 2200 continue
  
        close(11)
        close(12)
  
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine findword(card,word,word_upper,nword,iffound)
        implicit real *8 (a-h,o-z)
c 
        save
        character *1 blank,card(1),word(1),word_upper(1)
c 
        data blank/' '/
c 
        i=0
        j=1
c 
  
        iffound=1
        ifnonblank=0
        do 1400 ii=1,40
c 
        i=i+1
        if(card(i) .eq. blank) goto 1400
        ifnonblank=1
c 
c        if this location contains the proper letter
c        advance both i and j and continue
c 
         if( (card(i) .ne. word(j)) .and.
     1       (card(i) .ne. word_upper(j)) ) then
             iffound=0
             return
         endif
c 
         if( (card(i) .eq. word(j)) .or.
     1       (card(i) .eq. word_upper(j)) ) then
             j=j+1
c 
             if(j .gt. nword) then
                 iffound=iffound*ifnonblank
                 return
             endif
c 
             goto 1400
         endif
c 
 1400 continue
c 
         iffound=iffound*ifnonblank
c 
        return
        end
