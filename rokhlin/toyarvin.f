  
c 
c 
c 
c 
c 
        subroutine foolcosi(x,extcos,extsin)
        implicit real *16 (a-h,o-z)
        save
        dimension factin(0:60)
c 
c        compute, slowly but correctly, the sine and cosine
c        of x with extended precision
c 
        done=1
        two=2
        half=done/2
        n=20
        small=0.01
c 
        nhalf=0
        xx=x
ccccc         if(2 .ne. 3) goto 1500
        do 1200 i=1,100
        if(qabs(xx) .lt. small) goto 1400
        xx=xx*half
        nhalf=nhalf+1
 1200 continue
 1400 continue
 1500 continue
cccc         call prinf('nhalf=*',nhalf,1)
cccc         call prinq('xx=*',xx,1)
c 
c       calculate sine and cosine of xx
c 
        extcos=done
        extsin=0
c 
        d1=xx
        d2=xx**2
        dd=d2
        do 1600 i=1,n
        extcos=extcos+d2*factin(2*i)
        extsin=extsin+d1*factin(2*i-1)
        d1=d1*dd
        d2=d2*dd
 1600 continue
        extsin=-extsin
cccc          if(2 .ne. 3) return
c 
c       now, square the thing nhalf times
c 
        if(nhalf .eq. 0) return
        do 1800 i=1,nhalf
cccc          call prinf('i=*',i,1)
cccc          call prinq('extcos=*',extcos,1)
cccc          call prinq('extsin=*',extsin,1)
        d2=extcos**2-extsin**2
        d1=two*extcos*extsin
        extsin=d1
        extcos=d2
 1800 continue
        return
c 
c 
c 
c 
        entry foolcoi(ier)
c 
c       this is initialization entry point
c 
        coe=1
        factin(0)=1
        do 2200 i=1,55
        coe=-coe
        d=i
        factin(i)=factin(i-1)/d*coe
 2200 continue
         call prinq('factin as created*',factin(0),55)
        return
        end
