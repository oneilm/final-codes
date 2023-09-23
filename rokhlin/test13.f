        implicit real *8 (a-h,o-z)
        complex *16 a,b,c,z1,z2,z3,f1,f2,f3,
     1      root1,root2,ftest1,ftest2,root,zs(200),
     2      roots(200)
        complex *16 ima
        dimension errs(200)
        integer *4 niters(200)
        external funeva,funeva2,funeva3,funeva4
        data ima/(0.0d0,1.0d0)/
c
        call prini(6,13)
c
c       now, find all roots at once
c
        n=5

        z1=1
        z2=0
        z3=-1

        call fndroots(funeva4,n,roots,z1,z2,z3,
     1      eps,niters,errs)
         call prin2('after fndroots, roots=*',roots,n*2)
         call prinf('and niters=*',niters,n)
         call prin2('and errs=*',errs,n)
        stop
        end
c
c        
c
c
c
        subroutine funeva3(z,f)
        implicit real *8 (a-h,o-z)
        complex *16 z,f,z1,z2,z3,z4,ima
        data ima/(0.0d0,1.0d0)/


        f=12+22*z+23*z**2+2*z**3+z**4


        return
        end
c
c        
c
c
c

        subroutine funeva4(z,f)
        implicit real *8 (a-h,o-z)
        complex *16 z,f,z1,z2,z3,z4,ima
        data ima/(0.0d0,1.0d0)/

c
        f=60+158*z+115*z**2+80*z**3+5*z**4+2*z**5
c
        return
        end
c
c
c
c
c
        subroutine funeva(z,f)
        implicit real *8 (a-h,o-z)
        complex *16 z,f,z1,z2,z3,z4,ima
        data ima/(0.0d0,1.0d0)/

c
        f=4808+3470*z+669*z**2-10*z**3-z**4
        return
        end
c
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        this is the end of the debugging code and the beginning of the 
c        actual root-finder.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine fndroots(funcom,n,roots,z01,z02,z03,
     1      eps,niters,errs)
        implicit real *8 (a-h,o-z)
        complex *16 roots(1),zs(200),z01,z02,z03
        integer *4 niters(1)
        dimension errs(1)
c
c        this subroutine finds all roots of a polynomial via 
c        muller's process with numerical deflation
c
c                input parameters:
c
c  funcom - the name of the subroutine evaluating the polynomial
c        its calling sequence should be 
c
c             funcom(z,f)
c
c  n - the order of the polynomial (must be less than 100)
c  z01,z02,z03 - the starting values for the muller process
c        (the same values are used for all roots)
c  eps - the accuracy to which the roots will be computed 
c        (at least, an attempt to get them to this accuracy 
c        will be made)
c
c               output parameters:
c
c  roots - the roots of the polynomial (n of them)
c  niters - the numbers of iterations taken by the muller 
c        scheme for each of the roots (niter(i) is the 
c        number of iterations the algorithm took to obtain
c        the i-th root)
c  errs - the accuracy with which each root has been 
c        computed (errs(i) is the value of the polynomial 
c        at the putative i-th root)
c
c        . . . one root after the other, find them things 
c              via muller's process
c
        do 2000 i=1,n
          call prinf('in fndroots, i=*',i,1)
c
c       . . . find this root
c
         nzs=i-1
         maxit=200
        call muller(ier,z01,z02,z03,funcom,
     1      eps,roots(i),niters(i),zs,nzs,errs(i),maxit)
c
        zs(i)=roots(i)
          call prinf('in fndroots, niters(i)=*',niters(i),1)
 2000 continue
        return
        end
c
c
c
c
c
        subroutine funcomp(z,funcom,zs,m,f)
        implicit real *8 (a-h,o-z)
        complex *16 z,zs(1),f
c
c        calculate the undeflated function
c
        call funcom(z,f)
c
c        . . . deflate
c
         if(m .eq. 0) return
c
        do 1200 i=1,m
        f=f/(z-zs(i))
 1200 continue
        return
        end
c
c
c
c
c
        subroutine muller(ier,z01,z02,z03,funcom,
     1      eps,root,niter,zs,nzs,err,maxit)
        implicit real *8 (a-h,o-z)
        complex *16 z01,z02,z03,z1,z2,z3,f1,f2,f3,
     1      root1,root2,z4,f4,root,zs(1),ima
        data ima/(0.0d0,1.0d0)/
c
c        this subroutine finds a root of a polynomial via 
c        muller's process. the algorithm assumes that a 
c        number of roots has already been found, and avoids
c        these roots by the numerical deflation process.
c
c                input parameters:
c
c  z01,z02,z03 - the starting values for the muller process
c  funcom - the name of the subroutine evaluating the polynomial
c        its calling sequence should be 
c
c             funcom(z,f)
c
c  eps - the accuracy to which the root will be computed 
c        (at least, an attempt to get them to this accuracy 
c        will be made)
c  n - the order of the polynomial (must be less than 100)
c  zs - the array of roots already found 
c  nzs - the number of elements in array zs
c  maxit - the maximum number of muller iterations permitted
c
c               output parameters:
c
c  ier - the error return code.
c        ier=0 means a successful conclusion
c        ier=4 means that the user-specified pricision 
c          has not been achieved after maxit iterations
c  root - the root of the polynomial
c  niter - the number of iterations taken by the muller 
c        scheme for the root 
c  errs - the accuracy with which each root has been 
c        computed (errs(i) is the value of the polynomial 
c        at the putative i-th root)
c
c        . . . one root after the other, find them things 
c              via muller's process
c
c        initialize the muller iterations
c
        z1=z01        
        z2=z02        
        z3=z03        
         done=1
        call funcomp(z1,funcom,zs,nzs,f1)
        call funcomp(z2,funcom,zs,nzs,f2)
        call funcomp(z3,funcom,zs,nzs,f3)
        call reord(z1,z2,z3,f1,f2,f3)
c
c       conduct muller iterations
c
        ier=0
        dd2=1.0d37
        do 1200 i=1,maxit
        niter=i
c
c       . . . find the roots on this step
c
        call fndroot2(z1,z2,z3,f1,f2,f3,
     1      root1,root2)
c
c       . . . find the root nearest z3
c
        d1=cdabs(root1-z3)
        d2=cdabs(root2-z3)
c
        z4=root1
       if(d2 .lt. d1) z4=root2
c
c       shift the whole stack
c
        call funcomp(z4,funcom,zs,nzs,f4)
c
        z1=z2
        z2=z3
        z3=z4
        f1=f2
        f2=f3
        f3=f4
c
        call reord(z1,z2,z3,f1,f2,f3)
c
c        check if the desired accuracy has been achieved
c
        dd=cdabs(f3)
        if(dd .le. eps) goto 2000
 1200 continue
        ier=4
 2000 continue
        root=z3
        err=dd
        return
        end
c
c
c
c
c
        subroutine reord(z1,z2,z3,f1,f2,f3)
        implicit real *8 (a-h,o-z)
        complex *16 z1,z2,z3,f1,f2,f3,cd1,cd2
c
        if(cdabs(f1) .ge. cdabs(f2) ) goto 1200
        cd1=f1
        f1=f2
        f2=cd1
        cd2=z1
        z1=z2
        z2=cd2
 1200 continue
        if(cdabs(f2) .ge. cdabs(f3) ) goto 1400
        cd1=f2
        f2=f3
        f3=cd1
        cd2=z2
        z2=z3
        z3=cd2
 1400 continue
c
        if(cdabs(f1) .ge. cdabs(f2) ) goto 1600
        cd1=f1
        f1=f2
        f2=cd1
        cd2=z1
        z1=z2
        z2=cd2
 1600 continue
        return
        end
c
c
c
c
c
        subroutine fndroot2(z71,z72,z73,f1,f2,f3,
     1      root1,root2)
        implicit real *8 (a-h,o-z)
        complex *16 z1,z2,z3,f1,f2,f3,root1,root2,a,b,c,
     1      discr,z71,z72,z73,zcent
c
c        this subroutine finds the two roots root1, root2
c        of the second order polynomial taking the 
c        user-specified values f1, f2, f3 at the 
c        user-specified nodes z71, z72, z73
c
c        . . . shift and scale the points
c
        zcent=z71+z72+z73
        zcent=zcent/3
        z1=z71-zcent
        z2=z72-zcent
        z3=z73-zcent
c
        dd=cdabs(z1)+cdabs(z2)+cdabs(z3)
        z1=z1/dd
        z2=z2/dd
        z3=z3/dd
c
c       find the coefficients of the quadratic polynomial
c       defined by its values at the polint z1,z2,z3
c
         eps=1.0d-10
        a=( (f2-f1)/(z2-z1) - (f3-f1)/(z3-z1) )/(z2-z3)
        b=(f2-f1)/(z2-z1)-a*(z2+z1)
        c=f1-a*z1**2-b*z1
cccc         call prin2('in fndroot2, a=*',a,2)
cccc         call prin2('in fndroot2, b=*',b,2)
cccc         call prin2('in fndroot2, c=*',c,2)
c
c         if this is in reality a linear equation - act accordingly
c
        if (cdabs(a) .lt. eps) goto 1400
c
c       solve the quadratic equations with the coefficinets
c       a, b, c
c
        discr=b**2-4*a*c
        discr=cdsqrt(discr)
c
        root1=(-b+discr)/(2*a)
        root2=(-b-discr)/(2*a)
        goto 2000
 1400 continue
        root1=-c/b
        root2=1.0d20
 2000 continue
c
c       shift and scale back
c
        root1=root1*dd
        root2=root2*dd
        root1=root1+zcent
        root2=root2+zcent
        return
        end


