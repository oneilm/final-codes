        implicit real *8 (a-h,o-z)
        dimension a(100 000),alphas(1000 000),
     1      rhss(1000 000)
  
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER n'
        READ *,n
        CALL PRInf('n=*',n,1 )
  
  
        m=n/2
        m=n*2
        nalphas=30
        call eigtest(a,n,m,alphas,nalphas,rhss)
  
  
  
  
  
        stop
        end

c 
c 
c 
c 
        subroutine eigtest(a,n,m,alphas,nalphas,rhss)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,m),alphas(m,nalphas),ww(10000),
     1      ima,rhss(n,nalphas),tt(10000)
c 
        complex *16 w(2000 000),x(10000)
c
        data ima/(0.0d0,1.0d0)/
c 
c       construct the matrix a
c 
  
        done=1
cccc        ima=1
        eps=2.3
        eps=0.01
        do 1400 i=1,n
        do 1200 j=1,m
c 
        a(i,j)=eps *(i**2 +j**2) +ima*j*eps
 1200 continue
c 
        a(i,i)=-i+ima/i
 1400 continue
c
c       construct the solution
c
        call qrsolve_rand(m*2+2,x)

        call prin2('x as created*',x,m*2)
c
c        construct the vectors alphas
c
        call qrsolve_rand(m*2*nalphas+2,alphas)

        call prin2('alphas=*',alphas,n*nalphas*2)
c
c        construct the right-hand sides
c
        do 1800 i=1,nalphas
c
        do 1600 j=1,m
c
        ww(j)=x(j)*alphas(j,i)
 1600 continue
c
        call cleast_trans_matvec2(a,ww,rhss(1,i),n,m)

        call prinf('i=*',i,1)
        call prin2('rhss(i)=*',rhss(1,i),n*2)
c
 1800 continue

cccc        stop
c
        call prin2('a as created*',a,n*m*2)  
c
c       solve the system in the least squares sense
c
   
        eps=1.0d-9
        lw=2000 000
        ifnew=1
        call cleast_trans(ier,a,n,m,alphas,rhss,nalphas,
cccc     1      tt,eps,w,lw)
     1      ifnew,tt,eps,w,lw)

        call prin2('after cleast_trans, tt*',tt,m*2)
        call prin2('and x=*',x,m*2)


        do 2200 i=1,m
c
        ww(i)=tt(i)-x(i)
 2200 continue

        call prin2('and difference is*',ww,m*2)

        return
        end
  
c 
c 
c 
c 
c 
        subroutine qrsolve_rand(n,y)
        implicit real *8 (a-h,o-z)
        save
        dimension y(1)
        data ifcalled/0/
c 
c       generate pseudo-random numbers
c 
        if(ifcalled .ne. 0) goto 1100
c 
        done=1
        pi=atan(done)*4
        add=sqrt(7*done)
        x=1.1010101010101
        ifcalled=1
 1100 continue
c 
        do 1200 i=1,n
c 
        phi=x*100000*pi+add
        j=phi
        phi=(phi-j)
        x=phi
c 
        y(i)=x
 1200 continue
c 
        return
c 
c 
c 
c 
        entry qrsolve_rand_reset(x7)
c 
        done=1
        pi=atan(done)*4
        add=sqrt(7*done)
        x=x7
        ifcalled=1
        return
        end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        This is the end of the debugging code, and the beginning of
c        the actual least squares code
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine cleast_trans(ier,a,n,m,alphas,rhss,nalphas,
cccc     1      tt,eps,w,lw)
     1      ifnew,tt,eps,w,lw)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,m),alphas(m,nalphas),rhss(n,nalphas),
     1      tt(1)
c
        real *8 w(1)
c
c        This subroutine solves in the least squares sense a
c        collection of systems of linear equations
c
c        A * alpha_j (T) \sim f_j,                              (1)
c
c        with j=1,2,...nalphas. 
c
c        In (1) above, alpha_j is the j-th diagonal 
c        m \times m - matrix, rhs_j is the j-th right hand side, 
c        and T is the vector to be found (via a least squares 
c        procedure). 
c        
c        Thus, this subroutine solves in the least squares sense
c        a set of linesr systems (1), in which the systems differ
c        by the right-hand sides rhs_j and the disgonal matrices
c        alpha_j, while the matrix A is the same.
c
c 
c                         Input parameters:
c
c  a - the matrix in (1)
c  (n,m) - the dimensionalities of a
c  alphas - the array of nalphas vectors of length m; viewed as 
c        diagonal m \times m matrices (nalphas of them things)
c  rhss - the array of nalphas vectors of length n; viewed as
c        the right-hand sides in (1)
c  nalphas - the number of right-hand sides in (1)
c  ifnew - integer parameter telling the subroutine whether it
c        has to construct the SVD of the input matrix a. ifnew=1
c        means that the SVD has to be constructed. ifnew=0 means
c        that the decomposition has been constructed during a 
c        preceding call to this subroutine. Needless to say, this
c        option can only be used when the matrix a is the same
c        for these two calls; furthermore, the first 
c
c                   2*(n*m+m**2)+30
c
c        elements of the array w have to be unchanged between the 
c        two calls. PLEASE NOTE THAT THIS IS A CPU-TIME SAVING 
C        OPTION. IF YOU DO NOT WANT TO MESS WITH IT - SET IFNEW=1
C        AND DON'T WORRY (UNTIL YOYU NEED TO SAVE CPU TIME)
c  eps - the accuracy to which the least squares will be minimized
c  lw - the amount of space supplied to the subroutine in the work
c        array w (in real *8 elements)
c
c                          Output parameters:
c
c  ier - error return code;
c      ier=0 means successful execution
c      ier=64 means that the amount of memory supplied by the user
c        in the work array w is insufficient. This is a fatal error
c      ier=128 means that the amount of memory supplied by the user
c        in the work array w is grossly insufficient. This is a 
c        fatal error
c  tt - the complex vector solving equations (1) in the least squares
c       sense (all nalphas of them things)
c
c        allocate memory
c
        ier=0
        eps0=1.0d-14
c
        ivsvstar=1
        lvsvstar=m*m*2+10
c
        ivustar=ivsvstar+lvsvstar
        lvustar=m*n*2+10
c
        ibmatr=ivustar+lvustar
        lbmatr=m*m*2+10        
c
        ibmatr0=ibmatr+lbmatr
        lbmatr0=m*m*2+10        
c
        iu=ibmatr0+lbmatr0
        lu=m*n*2+10
c
        iv=iu+lu
        lv=n*m*2+10
c
        irhs=iv+lv
        lrhs=m*2+4
c
        irhs0=irhs+lrhs
        lrhs0=m*2+4
c
        iww=irhs0+lrhs0
        lww=m*2+4
c
        is=iww+lww
        ls=m
        if(m .lt. n) ls=n
        ls=ls*2+4
c
        iw1=is+ls
        lw1=lw-iw1
c
        if(lw1 .lt. 10 000) then
            ier=128
            return
        endif
c
        call cleast_trans0(jer,a,n,m,alphas,rhss,nalphas,
     1      w(iu),w(iv),w(iw1),lw1,w(ivsvstar),w(ibmatr),w(ibmatr0),
     2      tt,eps,w(ivustar),
     3      w(is),w(irhs0),w(irhs),w(iww),eps0,ifnew)
c
        if(jer .ne. 0) ier=64
c
        return
        end
c
c
c
c
c
        subroutine cleast_trans0(ier,a,n,m,alphas,rhss,nalphas,
     1      u,v,w,lw,vsvstar,bmatr,bmatr0,tt,eps,vustar,
     2      s,rhs0,rhs,ww,eps0,ifnew)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,m),alphas(m,nalphas),rhss(n,nalphas)
c
        complex *16 u(n,1),v(m,1),bmatr(m,m),w(m,1),vsvstar(m,m),
     1      rhs(1),ww(1),tt(1),bmatr0(m,m),rhs0(1),s(1),vustar(1)
c
c        construct the decomposition of the left-hand side and the rhs
c
        ier=0
        do 1200 i=1,m
        do 1100 j=1,m
c
        bmatr(j,i)=0
 1100 continue
c
        rhs(i)=0
 1200 continue
c
        ifsvd=1
        if(ifnew .ne. 1) ifsvd=0
        do 1800 i=1,nalphas
c        
        call cleast_trans_dec(jer,a,n,m,alphas(1,i),rhss(1,i),
     1      u,v,w,vsvstar,bmatr0,rhs0,ifsvd,
     2      s,vustar,eps0,lw)
c
        if(jer .ne. 0) then
            ier=32
            return
        endif
c        
        ifsvd=0
        do 1600 ii=1,m
        do 1400 j=1,m
c
        bmatr(j,ii)=bmatr(j,ii)+bmatr0(j,ii)
 1400 continue
c
        rhs(ii)=rhs(ii)+rhs0(ii)
 1600 continue
c
 1800 continue
c
c        solve the obtained system of linear equations in
c        the least squares sense
c
        eps2=eps0*10
        ifvect=1
        call cjaceig(bmatr,m,eps2,ifvect,bmatr0,nelems,w)
c
        call prin2('after cjaceig, eigenvalues are*',
     1      bmatr,m*2)
c        
        call cleast_trans_matavec(bmatr0,rhs,ww,m,m)
c
        do 2200 i=1,m
c
        ww(i)=ww(i)/bmatr(i,1)
c
        if(abs(bmatr(1,i)) .lt. eps) ww(i)=0
 2200 continue
c
        call cleast_trans_matvec(bmatr0,ww,tt,m,m)
c
        return
        end
c
c
c
c
c
        subroutine cleast_trans_dec(ier,a,n,m,alpha,f,
     1      u,v,w,vsvstar,bmatr,rhs,ifsvd,
     2      s,vustar,eps,lw)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,m),alpha(1),f(1)
c
        complex *16 u(n,1),v(m,1),s(1),bmatr(m,m),
     1      w(m,1),vsvstar(m,m),vustar(1),rhs(1)
c
c        SVD the matrix a
c
        ier=0
        if(ifsvd .eq. 0) goto 1800
c
        call csvdpiv(jer,a,n,m,u,v,s,ncols,eps,
     1      w,lw,ltot)
c
        if(jer .ne.0) then
            ier=16
            return
        endif
c        
        call prin2('after csvdpiv, s=*',s,ncols*2)
        call prinf('while ncols=*',ncols,1)
c
c        calculate V * S * V^*
c
        do 1600 i=1,ncols
        do 1400 j=1,m
c
        w(j,i)=v(j,i)*s(i)
 1400 continue
 1600 continue
c
        call cleast_trans_matmua(w,v,m,m,ncols,vsvstar)
c
c        construct the matrix vustar
c
        call cleast_trans_matmua(v,u,m,n,ncols,vustar)
c
 1800 continue
c
c        construct the whole matrix in the left side
c
        do 2400 i=1,m
        do 2200 j=1,m
c
        bmatr(j,i)=vsvstar(j,i)*conjg(alpha(j))*alpha(i)
 2200 continue
 2400 continue
c
c        construct the RHS
c
        call cleast_trans_matvec(vustar,f,rhs,m,n)
c
        do 2600 i=1,m
c
        rhs(i)=rhs(i)*conjg(alpha(i))
 2600 continue
c
        return
        end
c 
c 
c 
c 
c 
        subroutine cleast_trans_matmua(a,b,n,m,ncols,c)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,ncols),b(m,ncols),c(n,m),cd
c 
        do 1600 i=1,n
        do 1400 j=1,m
c 
        cd=0
        do 1200 k=1,ncols
        cd=cd+a(i,k)*dconjg(b(j,k))
 1200 continue
        c(i,j)=cd
 1400 continue
 1600 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine cleast_trans_matavec(a,x,y,n,m)
        implicit  complex *16 (a-h,o-z)
        save
        dimension a(n,m),x(n),y(m)
c 
        do 1400 i=1,n
        d=0
        do 1200 j=1,m
        d=d+conjg(a(j,i))*x(j)
 1200 continue
        y(i)=d
 1400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine cleast_trans_matvec2(a,x,y,n,m)
        implicit  complex *16 (a-h,o-z)
        save
        dimension a(n,m),x(m),y(n)
c 
        do 1400 i=1,n
        d=0
        do 1200 j=1,m
        d=d+a(i,j)*x(j)
 1200 continue
        y(i)=d
 1400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine cleast_trans_matvec(a,x,y,n,m)
        implicit  complex *16 (a-h,o-z)
        save
        dimension a(n,m),x(n),y(m)
c 
        do 1400 i=1,n
        d=0
        do 1200 j=1,m
        d=d+a(i,j)*x(j)
 1200 continue
        y(i)=d
 1400 continue
        return
        end






