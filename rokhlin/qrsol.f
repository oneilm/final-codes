        implicit real *8 (a-h,o-z)
        dimension a(1000 000),b(1000 000),a2(1000 000)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
         PRINT *, 'ENTER n'
         READ *,n
         CALL PRINF('n=*',n,1 )
c 
        call testit(a,n,b,a2)
  
c 
        stop
        end
c 
c 
c 
c 
c 
        subroutine testit(a,n,b,a2)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(n,n),b(n,n),
     1      x(1000),y(1000),x2(1000),a2(1),diffs(1000),
     2      work(1000 000),rhs(1000)
c 
c 
c        create the test matrix
c 
        do 1400 i=1,n
        do 1200 j=1,n
c 
        a(j,i)=i**2-j
 1200 continue
  
        a(i,i)=a(i,i)+i**3
 1400 continue
  
        call arrmove(a,a2,n*n*2)
  
  
  
  
  
cccc        call prin2('test matrix as created *',a,n*n)
  
  
  
c 
c       construct the test vector
c 
        done=1
c 
        do 2200 i=1,n
c 
        x(i)=done/(i**2+1)
 2200 continue
  
        call qrmatvec(a,n,x,y)
  
        call prin2('x as created*',x,n)
        call prin2('y as created*',y,n)
  
  
        call arrmove(y,rhs,n)
  
  
c 
c       construct its QR decomposition
c 
  
  
        t1=clotatim()
  
  
  
        call qrsol(a,n,rhs,rcond)
  
        call prin2('after qrsol, rhs=*',rhs,n)
  
        t2=clotatim()
  
        call prin2('cpu time for qrsol=*',t2-t1,1)
        call prin2('after qrsol, rcond=*',rcond,1)
  
  
        t3=clotatim()
  
        call grmsol(a2,n,y,x2,rcond2)
  
  
cccc        call orthom(a2,n,work,cond)
  
        t4=clotatim()
  
  
        call prin2('and CPU time for grmsol=*',t4-t3,1)
  
  
        call prin2('after grmsol, rcond2=*',rcond2,1)
  
  
        stop
        end
  
c 
c 
c 
c 
c 
        subroutine qrmatvec(a,n,x,y)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(n,n),x(1),y(1),cd
c 
c        apply the matrix a to the vector x obtaining y
c 
        do 1400 i=1,n
        cd=0
        do 1200 j=1,n
        cd=cd+a(i,j)*x(j)
 1200 continue
        y(i)=cd
 1400 continue
        return
        end
  
  
c 
c 
c 
c 
c 
        subroutine arrmove(x,y,n)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension x(1),y(1)
c 
        do 1200 i=1,n
c 
        y(i)=x(i)
 1200 continue
c 
        return
        end
c 
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c        This is the end of the testing code and the begtinning of the
c        linear solver proper
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccvc
c 
c 
c 
        subroutine qrsol(a,n,rhs,rcond)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(n,n),u(2,2),aa(2),rhs(1)
c 
c        This subroutine solves the linear system with the real matrix
c        a and the right-hand side rhs via the standard QR-decomposition.
c        In the process, both rthe matrix and rhe right-had side are
c        destroyed.
c 
c                     Input parameters:
c 
c  a - the input matrix
c  n - the dimensionality of the system to be solved
c  rhs - the right-hand side
c 
c                     output parameters:
c 
c  rhs - the solution
c  rcond - a fairly crude estimate of the condition number
c 
c       . . . eliminate the upper right triangle
c 
        do 2000 i=1,n-1
c 
        do 1600 j=n,i+1,-1
c 
        aa(1)=a(j-1,i)
        aa(2)=a(j,i)
        call qrrotfn0(aa,u)
c 
        do 1400 k=1,n
c 
        d1=u(1,1)*a(j-1,k)+u(1,2)*a(j,k)
        d2=u(2,1)*a(j-1,k)+u(2,2)*a(j,k)
c 
        a(j-1,k)=d1
        a(j,k)=d2
 1400 continue
c 
c 
        d1=u(1,1)*rhs(j-1)+u(1,2)*rhs(j)
        d2=u(2,1)*rhs(j-1)+u(2,2)*rhs(j)
c 
        rhs(j-1)=d1
        rhs(j)=d2
c 
 1600 continue
 2000 continue
c 
c       estimate the condition number
c 
        dmax=-1
        dmin=abs(a(1,1))
        do 2200 i=1,n
c 
        if(dmax .lt. abs(a(i,i)) ) dmax=abs(a(i,i))
        if(dmin .gt. abs(a(i,i)) ) dmin=abs(a(i,i))
 2200 continue
c 
        if(dmin .eq. 0) then
            rcond=-1
            return
        endif
c 
        rcond=dmax/dmin
c 
c        now, apply the inverse of the triangular matrix to the
c        transformed rhs
c 
        do 2600 i=n,2,-1
c 
        do 2400 j=1,i-1
c 
        rhs(j)=rhs(j)-a(j,i)/a(i,i)*rhs(i)
 2400 continue
 2600 continue
c 
        do 2800 i=1,n
        rhs(i)=rhs(i)/a(i,i)
 2800 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine qrrotfn0(a,u)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(2),u(2,2)
c 
        u22=-a(1)
        u12=a(2)
c 
        d=sqrt(u22**2+u12**2)
c 
        if(d .eq. 0) then
c 
            u(2,2)=1
            u(1,2)=0
            u(1,1)=1
            u(2,1)=0
            return
        endif
c 
        u(2,2)=u22/d
        u(1,2)=u12/d
        u(1,1)=-u(2,2)
        u(2,1)=u(1,2)
        return
        end
  
