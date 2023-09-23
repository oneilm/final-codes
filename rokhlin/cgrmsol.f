        implicit real *8 (a-h,o-z)
        complex *16 a(1100),c(1000)
c 
        call prini(6,13)
c 
c       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER n '
        READ *,n
        CALL PRINF('n=*',n,1 )
C 
        call invtes(a,n,c)
c 
        stop
        end
c 
c 
c 
c 
c 
         subroutine invtes(a,n,c)
         implicit real *8 (a-h,o-z)
        save
         complex *16 a(n,n),ima,
     1        c(n,n),x(100),y(100),z(100),xx(100)
        data ima/(0.0d0,1.0d0)/
c 
c       construct the test matrix
c 
        do 1400 i=1,n
        do 1200 j=1,n
        a(i,j)=i**2-j**2*ima
c 
cccc        a(i,j)=0
 1200 continue
        a(i,i)=i*100-ima*i**2
 1400 continue
c 
        do 1800 i=1,n
        do 1600 j=1,n
        c(j,i)=a(j,i)
 1600 continue
 1800 continue
c 
         call prin2('a as created*',a,n*n*2)
c 
c      . . . create the test vector
c 
        do 2200 i=1,n
        x(i)=i+i**2*ima
c 
        x(i)=i**3+ima*i
        xx(i)=x(i)
 2200 continue
         call prin2('x as created *',x,n*2)
c 
c        apply the inverse of a to x
c 
        call cgrmsol(a,n,x,y,rcond)
  
        call prin2('after cgrmsol, rcond=*',rcond,1)
  
        call prin2('after cgrmsol, y=*',y,n*2)
  
        call cmatvec(c,n,y,z)
  
          call prin2('a y =*',z,n*2)
c 
        do 2400 i=1,n
        z(i)=z(i)-xx(i)
 2400 continue
c 
          call prin2('and the difference*',z,n*2)
  
  
        return
        end
c 
c 
c 
c 
c 
        subroutine cmatvec(a,n,x,y)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,n),x(1),y(1),cd
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c        This is the end of the debugging code and the beginning of the
c        linear solver subroutine proper
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine cgrmsol(b,n,x,y,rcond)
        implicit real *8 (a-h,o-z)
        save
        complex *16 b(n,n),cd,x(1),y(1)
c 
c        This subroutine uses the pivoted Gram-schmidt procedure
c        to solve a linear system with complex coefficients.
c        This is a primitive routine in that it solves the system
c        with a single right-hand side, and will perform the
c        Gram-Schmidt again and again for each new right-hand
c        side. It should not be used in this regime!
c 
c                 Input parameters:
c 
c  b - the matrix of the system; destroyed by the subroutine
c  n - dimensionality of the system
c  x - the right-hand side; destroyed by the subroutine
c 
c 
c                 Output parameters:
c 
c  y - the solution of the system
c  rcond - estimate of the condition number of the matrix b
c 
c        . . . transpose the matrix b
c 
        do 1400 i=1,n
        do 1200 j=1,i
        cd=b(i,j)
        b(i,j)=b(j,i)
        b(j,i)=cd
 1200 continue
 1400 continue
c 
c        apply the Gram-Schmidt process to the matrix b
c 
        call cgrmsol0(b,n,x,y,rcond)
c 
c       apply the adjoint of the Gram-Schmidt'ed matrix to the vector
c 
        do 2400 i=1,n
        cd=0
        do 2200 j=1,n
        cd=cd+dconjg(b(i,j))*x(j)
 2200 continue
        y(i)=cd
 2400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine cgrmsol0(b,n,x,rnorms,rcond)
        implicit real *8 (a-h,o-z)
        save
        complex *16 b(n,n),cd,x(1)
        dimension rnorms(1)
c 
c       . . . prepare the array of values of norms
c             of columns
c 
        done=1
        dtot=0
        do 1400 i=1,n
        d=0
        do 1200 j=1,n
        d=d+b(j,i)*dconjg(b(j,i))
 1200 continue
        rnorms(i)=d
        dtot=dtot+d
 1400 continue
c 
c       . . . conduct gram-schmidt iterations
c 
        do 4000 i=1,n
c 
c       find the pivot
c 
        ipivot=i
        rn=rnorms(i)
c 
        do 2200 j=i+1,n
        if(rnorms(j) .le. rn) goto 2200
        rn=rnorms(j)
        ipivot=j
 2200 continue
 2400 continue
c 
c       put the column number ipivot in the i-th place
c 
        do 2600 j=1,n
        cd=b(j,i)
        b(j,i)=b(j,ipivot)
        b(j,ipivot)=cd
 2600 continue
c 
        cd=x(i)
        x(i)=x(ipivot)
        x(ipivot)=cd
c 
        d=rnorms(ipivot)
        rnorms(ipivot)=rnorms(i)
        rnorms(i)=d
 2700 continue
c 
c       orthogonalize the i-th column to all preceeding ones
c 
        if(i .eq. 1) goto 2790
c 
        do 2780 j=1,i-1
c 
        call cgrmscap(b(1,i),b(1,j),n,cd)
c 
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*cd
 2770 continue
c 
        x(i)=x(i)-x(j)*cd
c 
 2780 continue
 2790 continue
c 
c       normalize the i-th column
c 
        call cgrmscap(b(1,i),b(1,i),n,cd)
c 
c       if the norm of the remaining part of the matrix
c       is sufficiently small - terminate the process
c 
        d=cd
c 
        d=done/dsqrt(d)
        do 2800 j=1,n
        b(j,i)=b(j,i)*d
 2800 continue
c 
        x(i)=x(i)*d
c 
c        orthogonalize everything else to it
c 
        do 3200 j=i+1,n
c 
        call cgrmscap(b(1,i),b(1,j),n,cd)
        cd=dconjg(cd)
c 
         rrn=0
        do 3000 l=1,n
        b(l,j)=b(l,j)-b(l,i)*cd
        rrn=rrn+b(l,j)*dconjg(b(l,j))
 3000 continue
        x(j)=x(j)-x(i)*cd
        rnorms(j)=dsqrt(rrn)
 3200 continue
 4000 continue
  
 4200 continue
c 
        rcond=rnorms(1)/rnorms(n)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine cgrmscap(x,y,n,prod)
        implicit complex *16 (a-h,o-z)
        save
        dimension x(1),y(1)
c 
        prod=0
        do 1200 i=1,n
        prod=prod+x(i)*dconjg(y(i))
 1200 continue
        return
        end
