c
c    (c) 2007 Vladimir Rokhlin from his final codes
c
        implicit real *8 (a-h,o-z)
        dimension zs(2,100 000),xs(100 000),w(100 000)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
c 
        PRINT *, 'ENTER num'
        READ *,num
        call prinf('num=*',num,1 )
  
c 
c       construct the geometry of the test
c 
  
        call newrand(num,xs)
  
        call prin2('xs as constructed*',xs,num)
  
  
        iw=21
        call zquaplot(iw,xs,num/2,1,'random points*')
  
  
        do 1200 i=1,num/2
        zs(1,i)=xs(i)
        zs(2,i)=xs(i+num/2)
 1200 continue
  
  
        iw=22
        call zquaplot(iw,zs,num/2,1,'other random points*')
  
        iw=23
        call zquaplot2(iw,xs,num/2,1,zs,num/2,1,'both*')
  
  
  
c 
c       sort the obtained elements
c 
        call rsortanyn(xs,num,w)
  
        call prin2('and xs after sorting*',xs,num)
  
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine newrand(num,xs)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1)
c 
        integer  *4 alpha,beta,gamma
c 
        done=1
        m=37333
c 
        alpha=2**7+9
        beta=1
        gamma=2**22
c 
        do 2200 i=1,num
c 
        m1=m*alpha+beta
        j=m1/gamma
        m=m1-j*gamma
c 
        xs(i)=m*done/gamma
 2200 continue
        return
        end
