        implicit real *8 (a-h,o-z)
        dimension funs(1000),vals(1000 000),prods(100 000),
     1      w(30 000)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER x'
        READ *,x
        CALL PRIN2('x=*',x,1 )
C 
        PRINT *, 'ENTER nmax'
        READ *,nmax
        CALL PRINf('nmax=*',nmax,1 )
C 
cccc        nmax=6
        ifinit=1
        call herfuns(x,nmax,funs,ifinit,w)
  
  
        call prin2('funs as created*',funs,nmax)
c 
c       test the orthonormality of the obtained Hermite functions
c 
  
cccc        stop
  
  
        ntest=1000
        call herstest(nmax,ntest,vals,prods,w)
c 
        stop
        end
  
c 
c 
c 
c 
c 
        subroutine herstest(nmax,ntest,vals,prods,w)
        implicit real *8 (a-h,o-z)
        save
        dimension funs(1000),xs(10 000),whts(10 000),
     1      vals(ntest,nmax),prods(nmax,nmax),w(1)
c 
c        find the interval on which the inner products are to
c        be evaluated
c 
        a=1
        eps=1.0d-15
        do 1200 i=1,100
c 
cccc        call prinf('i=*',i,1)
  
        ifinit=0
cccc        call herfs0(a,nmax,funs,ifinit)
        call herfuns(a,nmax,funs,ifinit,w)
  
  
cccc        call prin2('and funs (nmax)=*',funs(nmax),1)
  
        if(abs(funs(nmax)) .lt. eps) goto 1400
        a=a*1.1
 1200 continue
  
 1400 continue
c 
c       construct the equispaced nodes on the interval [-a,a]
c 
        h=2*a/(ntest-1)
        do 1600 i=1,ntest
c 
        xs(i)=-a+(i-1)*h
        whts(i)=h
 1600 continue
c 
c       construct the matrix of values of Hermite functions
c 
        call prini(0,0)
  
        do 2000 i=1,ntest
c 
        ifinit=0
cccc        call herfs0(xs(i),nmax,funs,ifinit)
        call herfuns(xs(i),nmax,funs,ifinit,w)
c 
        do 1800 j=1,nmax
c 
        vals(i,j)=funs(j)
 1800 continue
 2000 continue
c 
        call prini(6,13)
  
  
cccc        call prin2('and vals are*',vals,ntest*nmax)
c 
c        calculate the inner products
c 
        do 2400 i=1,nmax
        do 2200 j=1,nmax
c 
        call scapro(vals(1,i),vals(1,j),whts,ntest,prod)
        prods(i,j)=prod
 2200 continue
 2400 continue
c 
        call prin2('and prods are*',prods,nmax*nmax)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine scapro(f,g,whts,n,prod)
        implicit real *8 (a-h,o-z)
        save
        dimension whts(1),f(1),g(1)
c 
        prod=0
        do 1200 i=1,n
        prod=prod+f(i)*g(i)*whts(i)
 1200 continue
        return
        end
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code and the beginning
c        of the Hermite function code proper
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine herfuns(x,nmax,funs,ifinit,w)
        implicit real *8 (a-h,o-z)
        save
        dimension funs(1),w(1)
c 
c        This subroutine evaluates the first nmax+1 Hermite
c        functions of the user-supplied argument x. The n-th
c        Hermite function is defined by the formula
c 
c        fun_n=c_n \cdot H_n(x) \cdot e^(-x^2/2),                      (1)
c 
c        where H_n is the n-th Hermite polynomial, and the
c        coefficients c_n are chosen in such a way that
c        the functions fun_0, fun_1, ... are orthonormal on R^1.
c 
c                        Input parameters:
c 
c  x - the point on R^1 where the Hermite functions are to
c        be evaluated
c  nmax - the maximum order n for which the function will be evaluated
c        (well, maybe it evaluates nmax+2 of them or so)
c  ifinit - the integer parameter telling the subroutine whether it should
c        perform certain initialization:
c     ifinit=1 tells the subroutine to initialize the array w;
c     ifinit=0 will cause the subroutine to skip the initialization
c        ifinit=0 will save some 50% of cpu time
c  w - array containing the initialization data. This is only an input
c        parameter if ifinit has been set to 0; in this case, the first
c        3*nmax+30 elements of w must be unchanged after the preceeding
c        call to this subroutine with ifinit-1
c 
c                        Output parameters:
c 
c  funs - the values of the first nmax+2 Hermite functions defined
c        by (1) above at the suer-specified point x
c  w - array in which the initialization data are return; must be at
c        least 3*nmax+30 real *8 locations. This is only an output
c        parameter if ifinit has been set to 1; otherwise, it is an
c        input parameter
c 
c 
c       . . . if the user so requested, allocate memory for the
c             arrays coefs1, coefs2(10000), scales in the user-provided
c             array w
c 
        if(ifinit .eq. 0) goto 1200
c 
        icoefs1=1
        lcoefs1=nmax+10
c 
        icoefs2=icoefs1+lcoefs1
        lcoefs2=nmax+10
c 
        iscales=icoefs2+lcoefs2
c 
 1200 continue
c 
c        evaluate the hermite functions
c 
        call herfs0(x,nmax,funs,ifinit,
     1      w(icoefs1),w(icoefs2),w(iscales))
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine herfs0(x,nmax,funs,ifinit,coefs1,coefs2,
     1      scales)
        implicit real *8 (a-h,o-z)
        save
        dimension funs(1),coefs1(1),coefs2(1),scales(1)
c 
c 
c       if the user so requested - initialize the array of
c       coefficients of the recursion
c 
        if(ifinit .eq. 0) goto 1060
c 
        nthresh=20
        two=2
        done=1
        pi=atan(done)*4
        sqsqpi=1/sqrt(sqrt(pi))
        big=exp(nthresh*done)
c 
        do 1050 n=1,nmax+2
c 
        n1=n+1
c 
        coefs1(n1)=-sqrt(n/(n+done))
        coefs2(n1)=sqrt(two/n1)
 1050 continue
 1060 continue
c 
c        initialize the array of scales
c 
        scales(1)=0
        scales(2)=0
c 
c       evaluate the Hermite functions of orders 0, 1
c 
        funs(1)=1
        funs(2)=2*x/sqrt(two)
c 
c       conduct the recursion
c 
        do 1200 n=1,nmax
c 
        n1=n+1
        d=coefs1(n1)*funs(n1-1)+ coefs2(n1)*x*funs(n1)
c 
        funs(n1+1)=d
c 
c       adjust the array of scales
c 
        scales(n1+1)=scales(n1)
c 
        if(abs(funs(n1+1)) .lt. big) goto 1200
c 
        funs(n1+1)=funs(n1+1)/big
        funs(n1)=funs(n1)/big
c 
        scales(n1+1)=scales(n1+1)+1
        scales(n1)=scales(n1)+1
 1200 continue
c 
        do 1400 i=1,nmax+2
c 
        d=exp(-x**2/2+scales(i)*nthresh)
        funs(i)=funs(i)*d*sqsqpi
 1400 continue
c 
        return
        end
