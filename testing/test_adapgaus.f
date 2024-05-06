        implicit real *8 (a-h,o-z)
        external fun
c
        call prini(6,13)
C
C       SET ALL PARAMETERS
C
c
        PRINT *, 'ENTER m '
        READ *,m
        CALL PRINF('m=*',m,1 )
c
c        integrate
c
        jtest=4
        a=0
        b=1
        eps=1.0d-15
c
        call adapgaus(ier,a,b,fun,jtest,par2,m,eps,
     1      rint,maxrec,numint)
        call prin2('after adapgaus, rint=*',rint,1)
        call prinf('after adapgaus, ier=*',ier,1)
        call prinf('after adapgaus, maxrec=*',maxrec,1)
        call prinf('after adapgaus, numint=*',numint,1)
 2200 format(' and more accurately, rint= ',e22.16)
        write(6,2200) rint
        write(13,2200) rint


        if(2 .ne. 2) stop
c
c        evaluate the integral via trapezoidal rule
c
        trap=0
        n=10000
        h=(b-a)/(n+1)
        do 2000 i=1,n
        t=a+i*h
        trap=trap+fun(t,jtest,dummy)
 2000 continue
        trap=trap*h
        call prin2('and trap=*',trap,1)



c
        stop
        end
c
c
c
c
c
c
        function fun(x,jtest,dummy)
        implicit real *8 (a-h,o-z)
c
        save
        fun=x**jtest
         half=1
         half=half/2
         done=1
         fun=dlog(x)+dlog(1-x)+dlog(dabs(x-half))
cccc         fun=dlog(dabs(x))+dlog(dabs(done-x))
cccc          fun=half/x
ccc         call prin2('inside fun, fun=*',fun,1)
        return
        end
