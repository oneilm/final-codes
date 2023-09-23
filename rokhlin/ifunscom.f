        implicit real *8 (a-h,o-z)
        dimension funs(10 000)
c 
        CALL PRINI(6,13)
C 
c       . . . construct the testing parameters
c 
c 
cccc         PRINT *, 'ENTER m'
cccc         READ *,m
cccc         CALL PRINf('m=*',m,1)
c 
cccc         PRINT *, 'ENTER eps'
cccc         READ *,eps
cccc         CALL PRIN2('eps=*',eps,1)
c 
         PRINT *, 'ENTER x'
         READ *,x
         CALL PRIN2('x=*',x,1)
c 
c        test the i-function computer
c 
         eps=1.0d-17
         call ifunscom(x,eps,m,funs(2))
         call prin2('after ifunscom, funs=*',funs,10)
c 
        stop
        end
c 
c 
c 
c 
c 
        subroutine ifunscom(x,eps,n,funs)
        implicit real *8 (a-h,o-z)
        save
        dimension funs(1)
c 
c       this subroutine evaluates all scaled I-type (modified)
c       Bessel functions of the argument x that are greater
c       than eps. More specifically, on exit, the i-the element
c       of the array funs is
c 
c              e^{-x} \cdot I_{i-1} (x).                            (1)
c 
c                       input parameters:
c 
c  x - the argument of which the modified Bessel functions are
c       to be evaluated.
c  eps - the threshold after which the functions (1) will be viewed
c        as zero.
c 
c                        output parameters:
c 
c  n - the number of scaled Bessel functions (1) returned. In other
c         words, e^{-x} \cdot I_{i-1} (x) < eps for all i>n.
c  funs - the array of functions (1). Note that it must have
c         element number 0 !!)
c 
c       . . . recurse up till the thing becomes fairly large
c 
        funs(1)=1
        funs(2)=0
        rlarge=1/eps*1.0d17
c 
        do 1200 i=1,1000
        ni=i
        funs(i+1)=-2*i/x*funs(i)+funs(i-1)
        if(funs(i+1) .gt. rlarge) goto 1400
 1200 continue
 1400 continue
         call prin2('funs after firt recursion*',funs,ni)
c 
c       now, starting with the position ni, recurse down
c 
        funs(ni+1)=0
        funs(ni)=1
c 
        do 1600 i=ni,1,-1
        funs(i-1)=2*i/x*funs(i)+funs(i+1)
 1600 continue
c 
c       sum up the squares of these recursed values, to get
c       the normalizing coefficient
c 
        i0=0
        sum=funs(i0)/2
        do 1800 i=1,ni
        sum=sum+funs(i)
 1800 continue
        sum=1/(sum*2)
c 
c       . . . scale
c 
        do 2000 i=i0,ni
        funs(i)=funs(i)*sum
 2000 continue
        n=ni
        return
        end
