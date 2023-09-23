        implicit real *8 (a-h,o-z)
        complex *16 rints2(1000),rhs(1000),
     1      diffs(1000),z,rint,ima,w(1000 000)
c 
        external fun2
c 
        data ima/(0.0d0,1.0d0)/
c 
        call prini(6,13)
c 
C 
        PRINT *, 'ENTER nmax'
        READ *,nmax
        CALL PRINf('nmax=*',nmax,1 )
c 
c       test the recursive evaluator of integrals
c 
        z=-3+ima*200
cccc        z=3
        z=0.4
  
  
cccc        z=1
c 
        a=0.3
        b=0.4
ccccc        b=0.300000001d0
  
  
        x0=0.7
cccc        a=0
  
        do 1400 i=1,nmax
  
        aa=a
        bb=b
        eps=1.0d-14
        m=20
c 
        call cadapgau(ier,aa,bb,fun2,i-1,z,m,eps,
     1      rint,maxrec,numint)
c 
        rints2(i)=rint/exp(z*x0)
c 
 1400 continue
  
        call prin2('and rints2=*',rints2,nmax*2)
  
  
cccc        stop
  
c 
c        construct  the dense matrix
c 
        lenw=1000 000
        call expolege_prod(ier,a,b,x0,z,nmax,rhs,w,lenw,lused)
  
c 
        call prin2('after expolege_prod0, rhs=*',rhs,nmax*2)
  
  
  
cccc        stop
c 
c       compare them things
c 
        do 1600 i=1,nmax
c 
        diffs(i)=rhs(i)-rints2(i)
 1600 continue
c 
        call prin2('and diffs=*',diffs,nmax*2)
c 
  
  
        stop
        end
c 
c 
c 
c 
        function fun2(x,n,z)
        implicit real *8 (a-h,o-z)
        save
        complex *16 z,d,fun2
c 
        call legepol(x,n,poln,dern)
c 
        d=exp(z*x)*poln
        fun2=d
c 
        return
        end
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code, and the beginning
c        of the code for the evaluation of Legendre-exponential
c        integrals
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine expolege_prod(ier,a,b,x0,z,nn,rhs,w,lenw,lused)
        implicit real *8 (a-h,o-z)
        save
        complex *16 z
        dimension w(1),rhs(1)
c 
c 
c       This subroutines evaluates integrals of the form
c 
c       \int_a^b P_n(x) \cdot e^{z \cdot (x-x0)} dx                     (1)
c 
c       for n=0,1,...,nn, where P_n denotes the n-th Legendre
c       polynomial, and a,z,nn are user-supplied parameters. The
c       integrals are evaluated to full double precision or a little
c       better.
C 
C                   Input parameters:
c 
c  a - the lower integration limit in (1) above
c  b - the upper integration limit in (1) above
c  x0 - the paremeter in the exponent in (1) above
c  z - the complex coefficient in the exponential in (1)
c  nn - the number of integrals (1) to be returned (actually, nn+1
c       integrals are returned)
c  lenw - the amount of space provided in the work array w (in
c       real *8 words)
c 
c                   Output parameters:
c 
c  ier - error return code.
c      ier=0 means successful execution
c      ier=8 means that the amount of space provided by the user
c            in the work array w (as specified by the parameter
c            lenw) is insufficient. Obviously, this is a fatal
c            error
c  rhs - the (complex) array containing the integrals (1). The name
c       is incongruous; if the user does not like it, he should
c       feel free to write his own code.
c  lused - the amount of space (in real *8 elements) in array w
c       actually used by the code
c 
c 
c        . . . depending on the values of the user-supplied
c              parameters z, nn, determine the size nnn of
c              the linear system to be solved
c 
        nadd=60
c 
c       PLEASE NOTE THAT THE PARAMETER NADD DETERMINES THE
C       ACCURACY OF THE EVALUATION OF THE FUNCTIONS. NADD=60
C       INSURES FULL DOUBLE PRECISION. NADD=100 OR SO IS
c       SUFFICIENT FOR THE EXTENDED PRECISION
c 
        nnn=abs(z)+nadd
        nnn2=nn+nadd
        if(nnn2 .gt. nnn) nnn=nnn2
c 
c        allocate memory
c 
        irhs7=1
        lrhs7=nnn*2+30
c 
        iaa=irhs7+lrhs7
        laa=nnn*2+30
c 
        ibb=iaa+laa
        lbb=nnn*2+30
c 
        icc=ibb+lbb
        lcc=nnn*2+30
c 
        iuu=icc+lcc
        luu=nnn*2+30
c 
        ivv=iuu+luu
        lvv=nnn*2+30
c 
        iww=ivv+lvv
        lww=nnn*2+30
c 
        ipols=iww+lww
        lpols=nnn+20
c 
        ipols0=ipols+lpols
        lpols0=nnn+20
c 
        lused=ipols0+lpols0
c 
        call prinf('and lused=*',lused,1)
c 
        ier=0
        if(lused .gt. lenw) then
            ier=8
            return
        endif
c 
        call expolege_prod0(a,b,x0,z,nnn,w(irhs7),w(iaa),
     1      w(ibb),w(icc),w(iuu),w(ivv),w(iww),
     2      w(ipols),w(ipols0) )
  
        do 1400 i=1,nn*2+2
c 
        rhs(i)=w(i)
 1400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine expolege_prod0(a,b,x0,z,nn,rhs,
     1     aa,bb,cc,uu,vv,ww,pols,pols0)
        implicit real *8 (a-h,o-z)
        save
        complex *16 rhs(1),rints(4),z,rexa,rexb,cd,
     1      aa(1),bb(1),cc(1),uu(1),vv(1),ww(1)
c 
        dimension pols(1),pols0(1)
c 
c        construct the Legendre polynomials of a and of 0
c 
        call legepols(b,nn+5,pols)
c 
        call legepols(a,nn+5,pols0)
c 
        do 2400 n=2,nn
c 
        rexb=exp((b-x0)*z)
        rexa=exp((a-x0)*z)
c 
        tnp1=pols(n+2)
        tnm1=pols(n)
c 
        rhs(n)=( tnp1/(2*n+1)-tnm1/(2*n+1)) * rexb
c 
        tnp1=pols0(n+2)
        tnm1=pols0(n)
c 
        cd=( tnp1/(2*n+1)-tnm1/(2*n+1)) * rexa
        rhs(n)=rhs(n)-cd
c 
        bb(n)=z/(2*n+1)
        cc(n)=-z/(2*n+1)
c 
 2400 continue
c 
        bb(1)=0
        cc(1)=0
c 
        rints(1)=(exp(z*(b-x0))-exp(z*(a-x0)))/z
c 
        rints(2)=b/z*exp((b-x0)*z)-a/z*exp((a-x0)*z)-
     1      (exp((b-x0)*z)-exp((a-x0)*z))/z**2
c 
        if(abs((b-a)*z) .lt. 0.01) then
c 
            az=(a-b)*z
c 
             cd=
     1       1+az/2+az**2/6+az**3/24+az**4/120+az**5/720+
     2       az**6/5040+az**7/40320+az**8/362880+az**9/3 628 800
c 
            rints(1)=exp(z*(b-x0))* cd *(b-a)
c 
        endif
c 
        if(abs(z) .lt. 0.01) then
c 
        rints(2)=
     a      0.5d0*(a - x0)**2 - 0.5D0*(b - x0)**2 +
     1     b*(b - x0) + a*(-a
     $     + x0) + (-0.33333333333333333333D0*a**3 + 0
     $     .33333333333333333333D0*b**3 + 0.5D0*a**2*x0 -
     4      0.5*b**2*x0 )*z + (0
     $     .041666666666666666667D0*(a - x0)**4 - 0
     $     .041666666666666666667D0*(b - x0)**4 - 0
     $     .16666666666666666667D0*a*(a - x0)**3 + 0
     $     .16666666666666666667D0*b*(b - x0)**3)*z**2 + (0
     $     .0083333333333333333333D0*(a - x0)**5 - 0
     $     .0083333333333333333333D0*(b - x0)**5 - 0
     $     .041666666666666666667D0*a*(a - x0)**4 + 0
     $     .041666666666666666667D0*b*(b - x0)**4)*z**3 + (0
     $     .0013888888888888888889D0*(a - x0)**6 - 0
     $     .0013888888888888888889D0*(b - x0)**6 - 0
     $     .0083333333333333333333D0*a*(a - x0)**5 + 0
     $     .0083333333333333333333D0*b*(b - x0)**5)*z**4 + (0
     $     .0001984126984126984127D0*(a - x0)**7 - 0
     $     .0001984126984126984127D0*(b - x0)**7 - 0
     $     .0013888888888888888889D0*a*(a - x0)**6 + 0
     $     .0013888888888888888889D0*b*(b - x0)**6)*z**5 + (0
     $     .000024801587301587301587D0*(a - x0)**8 - 0
     $     .000024801587301587301587D0*(b - x0)**8 - 0
     $     .0001984126984126984127D0*a*(a - x0)**7 + 0
     $     .0001984126984126984127D0*b*(b - x0)**7)*z**6 + (2
     $     .7557319223985890653D-6*(a - x0)**9 - 2
     $     .7557319223985890653D-6*(b - x0)**9 - 0
     $     .000024801587301587301587D0*a*(a - x0)**8 + 0
     $     .000024801587301587301587D0*b*(b - x0)**8)*z**7 + (2
     $     .7557319223985890653D-7*(a - x0)**10 - 2
     $     .7557319223985890653D-7*(b - x0)**10 - 2
     $     .7557319223985890653D-6*a*(a - x0)**9 + 2
     $     .7557319223985890653D-6*b*(b - x0)**9)*z**8 + (2
     $     .5052108385441718775D-8*(a - x0)**11 - 2
     $     .5052108385441718775D-8*(b - x0)**11 - 2
     $     .7557319223985890653D-7*a*(a - x0)**10 + 2
     $     .7557319223985890653D-7*b*(b - x0)**10)*z**9 + (2
     $     .0876756987868098979D-9*(a - x0)**12 - 2
     $     .0876756987868098979D-9*(b - x0)**12 - 2
     $     .5052108385441718775D-8*a*(a - x0)**11 + 2
     $     .5052108385441718775D-8*b*(b - x0)**11)*z**10 + (1
     $     .6059043836821614599D-10*(a - x0)**13 - 1
     $     .6059043836821614599D-10*(b - x0)**13 - 2
     $     .0876756987868098979D-9*a*(a - x0)**12 + 2
     $     .0876756987868098979D-9*b*(b - x0)**12)*z**11 + (1
     $     .1470745597729724714D-11*(a - x0)**14 - 1
     $     .1470745597729724714D-11*(b - x0)**14 - 1
     $     .6059043836821614599D-10*a*(a - x0)**13 + 1
     $     .6059043836821614599D-10*b*(b - x0)**13)*z**12
c 
        endif
c 
        rints(3)=0
c 
        rhs(1)=rints(2)+cc(1)*rints(3)
c 
        do 3200 i=1,nn
c 
        aa(i)=1
 3200 continue
c 
c       . . . solve the sparse system
c 
        call expolege_prod_fact(aa,bb,cc,nn,uu,vv,ww)
        call expolege_prod_solve(uu,vv,ww,nn,rhs)
c 
c       shift them things
c 
        do 3600 i=nn,1,-1
c 
        rhs(i+1)=rhs(i)
 3600 continue
c 
        rhs(1)=rints(1)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine expolege_prod_fact(a,b,c,n,u,v,w)
        implicit complex *16 (a-h,o-z)
        save
        dimension a(1),b(1),c(1),u(1),v(1),w(1),rhs(1)
c 
c        eliminate down
c 
        do 1200 i=1,n-1
        d=c(i+1)/a(i)
        a(i+1)=a(i+1)-b(i)*d
        u(i)=d
 1200 continue
c 
c        eliminate up
c 
        do 1400 i=n,2,-1
        d=b(i-1)/a(i)
        v(i)=d
 1400 continue
c 
c       scale the diagonal
c 
        done=1
        do 1600 i=1,n
        w(i)=done/a(i)
 1600 continue
c 
        return
c 
c 
c 
c 
        entry expolege_prod_solve(u,v,w,n,rhs)
c 
c        eliminate down
c 
        do 2400 i=1,n-1
        rhs(i+1)=rhs(i+1)-u(i)*rhs(i)
 2400 continue
c 
c        eliminate up
c 
        do 2600 i=n,2,-1
        rhs(i-1)=rhs(i-1)-rhs(i)*v(i)
 2600 continue
c 
c       scale
c 
        do 2800 i=1,n
        rhs(i)=rhs(i)*w(i)
 2800 continue
        return
        end
