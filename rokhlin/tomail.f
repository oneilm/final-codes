       implicit real *8 (a-h,o-z)
        dimension x(10 000),a(10 000),b(10 000),
     1     c(10 000),ua(10 000),ub(10 000),va(10 000),
     2     vb(10 000)
c 
        call prini(6,13)
c 
C 
C       SET ALL PARAMETERS
C 
C 
        PRINT *, 'ENTER k'
        READ *,k
        CALL PRINF('k=*',k,1 )
C 
        PRINT *, 'ENTER l'
        READ *,l
        CALL PRINF('l=*',l,1 )
C 
        PRINT *, 'ENTER m'
        READ *,m
        CALL PRINF('m=*',m,1 )
C 
        PRINT *, 'ENTER n'
        READ *,n
        CALL PRINF('n=*',n,1 )
C 
c        test the matrix least squares subroutine
c 
        delta=0
        eps=1.0d-12
  
        call mattest(a,b,c,k,l,m,n,eps,delta,x,
     1      ua,va,ub,vb)
c 
c      test the right matrix least squares subroutine
c 
       call mattest2(a,c,k,m,n,eps,delta,x,
     1      ua,va)
c 
c        test the left matrix least squares subroutine
c 
        call mattest3(a,c,k,m,n,eps,delta,x,ua,va)
  
  
c 
        stop
        end
c 
c 
c 
c 
c 
        subroutine mattest3(a,c,k,m,n,eps,delta,x,
     1      u,v)
        implicit real *8 (a-h,o-z)
        save
        dimension a(k,m),x(m,n),c(k,n),
     1      u(1),v(1),w(100 000)
c 
c        construct the test matrix a
c 
        done=1
        do 1400 i=1,m
        do 1200 j=1,k
        a(j,i)=done/(i**2+j**2+1)
 1200 continue
cccc        if(i .lt. k) a(i,i)=10+a(i,i)
 1400 continue
c 
        call prin2('in mattest3, a as created*',a,m*k)
c 
c       create the matrix c
c 
        done=1
        do 2400 i=1,n
        do 2200 j=1,k
cccc        c(j,i)=done/(i**2+j**2+1)
        c(j,i)=done/(i**3+j**4+1)
 2200 continue
ccc        if(i .lt. n) c(i,i)=10+c(i,i)
 2400 continue
c 
        call prin2('in mattest3, c as created*',c,k*n)
c 
c        construct the matrix x such that
c 
c                      A * A = C,                                 (1)
c 
        lw=100 000
c 
        call leamatll(ier,a,c,k,m,n,eps,delta,x,w,lw,ltot)
c 
        call prinf('after leamatr2,  ier=*',ier,1)
        call prinf('after leamatr2,  ltot=*',ltot,1)
        call prin2('after leamatr2,  X=*',x,k*m)
  
  
cccc        if(2 .ne. 3) stop
c 
c       multiply the matrices X, A  back together and compare
c       them with C
c 
cccc        call leamat(ier,x,k,m,a,m,n,u)
        call leamat(ier,a,k,m,x,m,n,u)
c 
        call prin2('and C as recomputed from A, X*',u,k*n)
        call prin2('again, original C=*',c,k*n)
  
        call diffprin(u,c,k*n,'and the difference*')
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine mattest2(a,c,k,m,n,eps,delta,x,
     1      u,v)
        implicit real *8 (a-h,o-z)
        save
        dimension a(m,n),x(k,m),c(k,n),
     1      u(1),v(1),w(100 000)
c 
c        construct the test matrix a
c 
        done=1
        do 1400 i=1,n
        do 1200 j=1,m
        a(j,i)=done/(i**2+j**2+1)
 1200 continue
cccc        if(i .lt. k) a(i,i)=10+a(i,i)
 1400 continue
c 
        call prin2('a as created*',a,m*n)
c 
c       create the matrix c
c 
        done=1
        do 2400 i=1,n
        do 2200 j=1,k
cccc        c(j,i)=done/(i**2+j**2+1)
        c(j,i)=done/(i**3+j**4+1)
 2200 continue
ccc        if(i .lt. n) c(i,i)=10+c(i,i)
 2400 continue
c 
        call prin2('and c as created*',c,k*n)
c 
c        construct the matrix x such that
c 
c                     X B = C,                                 (1)
c 
        lw=100 000
c 
        call leamatrr(ier,a,c,k,m,n,eps,delta,x,w,lw,ltot)
c 
        call prinf('after leamatr2,  ier=*',ier,1)
        call prinf('after leamatr2,  ltot=*',ltot,1)
        call prin2('after leamatr2,  X=*',x,k*m)
c 
c       multiply the matrices X, A  back together and compare
c       them with C
c 
        call leamat(ier,x,k,m,a,m,n,u)
c 
        call prin2('and C as recomputed from A, X*',u,k*n)
        call prin2('again, original C=*',c,k*n)
  
        call diffprin(u,c,k*n,'and the difference*')
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine mattest(a,b,c,k,l,m,n,eps,delta,x,
     1      ua,va,ub,vb)
        implicit real *8 (a-h,o-z)
        save
        dimension a(k,l),x(l,m),b(m,n),c(k,n),
     1      ua(k,1),va(l,1),ub(m,1),vb(n,1),
     2      w(100 000)
c 
c        construct the test matrices a,b
c 
        done=1
        do 1400 i=1,l
        do 1200 j=1,k
        a(j,i)=done/(i**2+j**2+1)
 1200 continue
cccc        if(i .lt. k) a(i,i)=10+a(i,i)
 1400 continue
c 
        call prin2('a as created*',a,l*k)
c 
        done=1
        do 1800 i=1,n
        do 1600 j=1,m
        b(j,i)=done/(i**2+j**2+1)
 1600 continue
cccc        if(i .lt. m) b(i,i)=10+b(i,i)
 1800 continue
c 
        call prin2('and b as created*',b,m*n)
c 
c       create the matrix c
c 
        done=1
        do 2400 i=1,n
        do 2200 j=1,k
        c(j,i)=done/(i**2+j**2+1)
 2200 continue
ccc        if(i .lt. n) c(i,i)=10+c(i,i)
 2400 continue
c 
        call prin2('and c as created*',c,k*n)
c 
c        construct the matrix x such that
c 
c                     A X B = C,                                 (1)
c 
        lw=100 000
  
        call leamatr(ier,a,b,c,k,l,m,n,eps,delta,x,
     1      w,lw,ltot)
  
  
        call prinf('after leamatr,  ier=*',ier,1)
        call prinf('after leamatr,  ltot=*',ltot,1)
        call prin2('after leamatr,  X=*',x,l*m)
c 
c       multiply the matrices A, X, B back together and compare
c       them with C
c 
        call leamat(ier,a,k,l,x,l,m,ua)
        call leamat(ier,ua,k,m,b,m,n,ub)
c 
        call prin2('and C as recomputed from A, B, X*',ub,k*n)
cccc        call prin2('again, original C=*',c,k*n)
c 
        call diffprin(ub,c,k*n,'and the difference*')
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine diffprin(a,b,n,mes)
        implicit real *8 (a-h,o-z)
        save
        dimension a(1),b(1),c(10 000)
        character *1 mes(1)
c 
        difmax=0
        do 1200 i=1,n
        c(i)=b(i)-a(i)
        dd=dabs(b(i)-a(i))
        if(difmax .lt. dd) difmax=dd
 1200 continue
        call prin2(mes,c,n)
        call prin2('and maximum difference is*',difmax,1)
        return
        end
c 
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of thedebugging code and the  beginning of
c        the actual matrix least squares subroutine
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        subroutine leamatll(ier,a,c,k,m,n,eps,delta,x,w,lw,ltot)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension x(m,n),c(k,n),a(k,m),w(1)
c 
c        this subroutine solves in the least squares sense the
c        equation
c 
c                       A(k,m) * X(m,n)  = C(k,n),               (1)
c 
c       where A, C are user-specified matrices, and X is the
c       matrix to be found. Note that the dimensionalities of
c       the matrices A,C, X are as general as could be. Note
c       that this subroutine does is only an interface; it
c       transposes the matrices A, C, and calls the subroutine
c       leamatrr (see). Finally, note that this subroutine
c       DOES NOT change the input matrices A, C.
c 
c 
c                          input parameters:
c 
c  a,c - matrices in (1)
c  k,m,n - dimensionalities in (1)
c  eps - the parameter telling the subroutine at which point a
c       a singular value of the matrix is to be declared to be zero
c       and ignored (see subroutine leasasr for details)
c  delta - the "Tikhonov constant" (see subroutine leasasr for details)
c  lw - the length of the user-supplied work array
c 
c                          output parameters:
c 
c  ier - error return code
c      ier=0 means successful execution
c      ier .ne. 0 means that the amount of space lw allocated supplied
c            in the work array w is insufficient
c  x - the unknown matrix in (1)
c  ltot - the length of the work array w actually used by the subroutine
c 
c                          work array:
c 
c  w - must be sufficiently long
c 
c       . . . transpose the arrays A, C to convert the problem into
c             one solvable by the subroutine leamatrr
c 
        call leatrans(a,k,m,w)
        call leatrans(c,k,n,w)
c 
c        solve the problem
c 
c               (X(m,n))^*  *  (A(k,m))^*  = (C(k,n))^*
c 
c        in the least squares sense
c 
        call leamatrr(ier,a,c,n,m,k,eps,delta,x,w,lw,ltot)
c 
c       transpose the matrices a, c, x
c 
        call leatrans(a,m,k,w)
        call leatrans(c,n,k,w)
        call leatrans(x,n,m,w)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine leatrans(a,m,n,w)
        implicit real *8 (a-h,o-z)
        save
        dimension a(m,n),w(n,m)
c 
        do 1400 i=1,m
        do 1200 j=1,n
c 
        w(j,i)=a(i,j)
 1200 continue
 1400 continue
c 
        call leamem(w,a,n*m)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine leamatrr(ier,a,c,k,m,n,eps,delta,x,w,lw,ltot)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension x(k,m),c(k,n),a(m,n),w(1)
c 
c        this subroutine solves in the least squares sense the
c        equation
c 
c                      X(k,m) * A(m,n) = C(k,n),               (1)
c 
c       where A, C are user-specified matrices, and X is the
c       matrix to be found. Note that the dimensionalities of
c       the matrices A,C, X are as general as could be. Note
c       that this is only a memory management routine; all the
c       work is done by the subroutine leamatr2 (see).
c 
c 
c                          input parameters:
c 
c  a,c - matrices in (1)
c  k,m,n - dimensionalities in (1)
c  eps - the parameter telling the subroutine at which point a
c       a singular value of the matrix is to be declared to be zero
c       and ignored (see subroutine leasasr for details)
c  delta - the "Tikhonov constant" (see subroutine leasasr for details)
c  lw - the length of the user-supplied work array
c 
c                          output parameters:
c 
c  ier - error return code
c      ier=0 means successful execution
c      ier .ne. 0 means that the amount of space lw allocated supplied
c            in the work array w is insufficient
c  x - the unknown matrix in (1)
c  ltot - the length of the work array w actually used by the subroutine
c 
c                          work array:
c 
c  w - must be sufficiently long
c 
c        . . . allocate memory for the subroutine
c              leamatr2
c 
c 
        iu=1
        lu=m*n
c 
        iv=iu+lu
        lv=m*n
c 
        is=iv+lv
        ls=m+n
c 
        iw=is+ls
c 
        lenw=lw-iw
c 
c        solve the equation (1)
c 
        call leamatr2(ier,a,c,k,m,n,eps,delta,x,
     1      w(iu),w(iv),w(is),w(iw),lenw,ltot1)
c 
        ltot=iw+ltot1+10
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine leamatr2(ier,a,c,k,m,n,eps,delta,x,
     1      u,v,s,w,lw,ltot)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension x(k,m),c(k,n),a(m,n),
     1      u(1),v(1),s(1),w(1)
c 
c        this subroutine solves in the least squares sense the
c        equation
c 
c                      X(k,m) * A(m,n) = C(k,n),               (1)
c 
c       where A, C are user-specified matrices, and X is the
c       matrix to be found. Note that the dimensionalities of
c       the matrices A,C, X are as general as could be.
c 
  
c 
c 
c        . . . construct the SVD of A
c 
        eps2=1.0d-13
        call prin2('before svdpivot, a=*',a,m*n)
        call svdpivot(ier,a,m,n,u,v,s,ncols,eps2,
     1      w,lw,ltot)
c 
        if(ier .ne. 0) return
c 
        call prin2('in leamatr2 after svdpivot, s=*',s,ncols)
cccc        call prin2('in leamatr2 after svdpivot, u=*',u,ncols*m)
cccc        call prin2('in leamatr2 after svdpivot, v=*',v,ncols*n)
c 
c        multiply c from the right by v
c 
        call leamat(ier,c,k,n,v,n,ncols,w)
c 
c        multiply CV from the right by S^{-1}
c 
cccc         call prin2('before leasasr, w=*',w,k*ncols)
        delta=0
        call leasasr(w,k,ncols,s,eps,delta)
c 
cccc         call prin2('after leasasr, w=*',w,k*ncols)
c 
c       multiply C V S^{-1} from the right by U^*
c 
        call leamar(ier,w,k,ncols,u,m,ncols,x)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine leasasr(a,m,n,s,eps,delta)
        implicit real *8 (a-h,o-z)
        save
        dimension a(m,n),s(n)
c 
c       multiply the matrix a by the generalized inverse of the
c       diagonal matrix with the vector s on the diagonal
c 
        do 2000 i=1,m
        do 1600 j=1,n
c 
        if(s(j) .gt. eps) goto 1400
        a(i,j)=0
        goto 1600
c 
 1400 continue
c 
        a(i,j)=a(i,j)/(delta+s(j))
c 
 1600 continue
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine leamatr(ier,a,b,c,k,l,m,n,eps,delta,x,
     1      w,lw,ltot)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension a(k,l),x(l,m),b(m,n),c(k,n),w(1)
c 
c        this subroutine solves in the least squares sense the
c        equation
c 
c                     A(k,l) * X(l,m) * B(m,n) = C(k,n),               (1)
c 
c       where A, B, C are user-specified matrices, and X is the
c       matrix to be found. Note that the dimensionalities of
c       the matrices A,B,C, X are as general as could be
c 
c 
c                          input parameters:
c 
c  a,b,c - matrices in (1)
c  k,l,m,n - dimensionalities in (1)
c  eps - the parameter telling the subroutine at which point a
c       a singular value of the matrix is to be declared to be zero
c       and ignored (see subroutine leasas for details)
c  delta - the "Tikhonov constant" (see subroutine leasas for details)
c  lw - the length of the user-supplied work array
c 
c                          output parameters:
c 
c  ier - error return code
c      ier=0 means successful execution
c      ier .ne. 0 means that the amount of space lw allocated supplied
c            in the work array w is insufficient
c  x - the unknown matrix in (1)
c  ltot - the length of the work array w actually used by the subroutine
c 
c                          work array:
c 
c  w - must be sufficiently long
c 
c 
c        . . . construct the SVD of A
c 
        ier=0
        isa=1
        lsa=min(k,l)+2
c 
        iua=isa+lsa
        lua=k*l+10
c 
        iva=iua+lua
        lva=k*l+10
c 
        iwork=iva+lva
c 
        eps2=1.0d-13
        lw1=lw-iwork-2
c 
        call prinf('in leamatr, lw1=*',lw1,1)
c 
        call svdpivot(ier,a,k,l,w(iua),w(iva),w(isa),ncolsa,eps2,
     1      w(iwork),lw1,ltot1)
c 
        ltot=ltot1+lwork
c 
        if(ier .ne. 0) return
c 
        call prinf('after first svdpivot, ncolsa=*',ncolsa,1)
        call prin2('after first svdpivot, sa=*',w(isa),ncolsa)
c 
c        perform garbage collection
c 
        lua2=k*ncolsa+2
        lva2=l*ncolsa+2
        iva2=iua+lua2
        lua=lua2
c 
        call leamem(w(iva),w(iva2),lva2)
c 
        iva=iva2
        lva=lva2
c 
c        . . . construct the SVD of B
c 
        isb=iva+lva
        lsb=min(m,n)+2
c 
        iub=isb+lsb
        lub=m*n+10
c 
        ivb=iub+lub
        lvb=m*n+10
c 
        iwork=ivb+lvb
        lw2=lw-iwork
cccc        call prinf('in leamatr, lw2=*',lw2,1)
c 
        eps2=1.0d-13
c 
        call svdpivot(ier,b,m,n,w(iub),w(ivb),w(isb),ncolsb,eps2,
     1      w(iwork),lw2,ltot2)
c 
        ii=iwork+ltot2
        if(ii .gt. ltot) ltot=ii
c 
cccc        call prinf('after second svdpivot, ltot=*',ltot,1)
c 
        if(ier .ne. 0) return
c 
        call prinf('after second svdpivot, ncolsb=*',ncolsb,1)
        call prin2('after second svdpivot, sb=*',w(isb),ncolsb)
c 
c        perform garbage collection
c 
        lub2=m*ncolsb+2
        lvb2=n*ncolsb+2
        ivb2=iub+lub2
c 
        call leamem(w(ivb),w(ivb2),lvb2)
c 
        ivb=ivb2
        lvb=lvb2
c 
        iwork=ivb+lvb
        lwork=ncolsa*n+10
c 
        ii=iwork+lwork
        call prinf('in leamatr, final ii=*',ii,1)
        if(ii .gt. ltot) ltot=ii
        if(ltot .lt. lw) goto 2200
        ier=2
        return
 2200 continue
c 
c         perform the remainder of the solution
c 
        call leamat0(c,k,l,m,n,eps,delta,x,
     1      w(isa),w(iua),w(iva),w(isb),w(iub),w(ivb),
     2      w(iwork),lw,ltot,ncolsa,ncolsb)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine leamat0(c,k,l,m,n,eps,delta,x,
     1      sa,ua,va,sb,ub,vb,work,lw,ltot,ncolsa,ncolsb)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension x(l,m),c(k,n),work(1),
     1      ua(k,1),va(l,1),ub(m,1),vb(n,1),sa(1),sb(1)
c 
c        this subroutine solves in the least squares sense the
c        equation
c 
c                     A(k,l) * X(l,m) * B(m,n) = C(k,n),               (1)
c 
c       where A, B, C are user-specified matrices, and X is the
c       matrix to be found. Note that the dimensionalities of
c       the matrices A,B,C, X are as general as could be.
c       Actually, this subroutine uses as its input the SVDs of
c       the matrices A,B, provided by the calling subroutine leamatr
c       (see), so that the equation (1) has the form
c 
c                    UA * SA * VA^* * X * UB * SB * VB^* = C       (2)
c 
c       . . . multiply the matrix c from the left by UA^* ;
c             note that the resulting matrix UAC is dimensioned
c             uac(ncolsa,n)
c 
        call leamatl(ier,ua,k,ncolsa,c,k,n,work)
c 
c       . . . multiply the matrix UAC by VB from the right;
c             note that the resulting matrix uacvb is
c             dimensioned uacvb(ncolsa,ncolsb)
c 
        call leamat(ier,work,ncolsa,n,vb,n,ncolsb,x)
c 
c        now, the equation (2) has assumed the form
c 
c                   SA * VA^* * X * UB * SB  = UACVB              (3)
c 
c       . . . multiply (3) by SA^{-1} from the left and
c             by SB^{-1} from the right in the appropriate least
c             squares sense
c 
         call leasas(sa,ncolsa,sb,ncolsb,x,eps,delta)
c 
c        now, the equation (3) has assumed the form
c 
c                  VA^* * X * UB  = UACVB                     (4)
c 
c       . . . multiply (4) by VA from the left and
c             by UB^* from the right, obtaining the
c             solution X
c 
        call leamat(ier,va,l,ncolsa,x,ncolsa,ncolsb,work)
c 
        call leamar(ier,work,l,ncolsb,ub,m,ncolsb,x)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine leamat(ier,a,k,l,b,m,n,c)
        implicit real *8 (a-h,o-z)
        save
        dimension a(k,l),b(m,n),c(k,n)
c 
c        this subroutine multiplies the matrix a by the matrix b,
c        getting the matrix c
c 
c        . . . check that the dimensionalities of the input
c              matrices agree
c 
        ier=4
        if(l .eq. m) ier=0
        if(ier .ne. 0) return
c 
c       multiply the matrix a by the matrix b getting c
c 
        do 2000 i=1,k
        do 1800 j=1,n
        d=0
        do 1600 kk=1,l
        d=d+a(i,kk)*b(kk,j)
 1600 continue
        c(i,j)=d
 1800 continue
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine leamatl(ier,a,k,l,b,m,n,c)
        implicit real *8 (a-h,o-z)
        save
        dimension a(k,l),b(m,n),c(l,n)
c 
c        this subroutine multiplies the matrix a^* by the matrix b,
c        getting the matrix c
c 
c        check that the dimensionalities of the input matrices agree
c 
        ier=4
        if(k .eq. m) ier=0
        if(ier .ne. 0) return
c 
c       multiply the matrix a by the matrix b getting c
c 
        do 2000 i=1,l
        do 1800 j=1,n
        d=0
        do 1600 kk=1,k
        d=d+a(kk,i)*b(kk,j)
 1600 continue
        c(i,j)=d
 1800 continue
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine leamar(ier,a,k,l,b,m,n,c)
        implicit real *8 (a-h,o-z)
        save
        dimension a(k,l),b(m,n),c(k,m)
c 
c        this subroutine multiplies the matrix a by the matrix b^*,
c        getting the matrix c
c 
c        check that the dimensionalities of the input matrices agree
c 
        ier=4
        if(l .eq. n) ier=0
        if(ier .ne. 0) return
c 
c       multiply the matrix a by the matrix b getting c
c 
        do 2000 i=1,k
        do 1800 j=1,m
        d=0
        do 1600 kk=1,l
        d=d+a(i,kk)*b(j,kk)
 1600 continue
        c(i,j)=d
 1800 continue
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine leasas(sa,ncolsa,sb,ncolsb,uacvb,eps,delta)
        implicit real *8 (a-h,o-z)
        save
        dimension sa(1),sb(1),uacvb(ncolsa,ncolsb)
c 
c       this subroutine multiplies the user-supplied matrix UACVB
c       by SA^{-1} from the left and by SB^{-1} from the right in
c       the appropriate least squares sense
c 
        do 1800 i=1,ncolsb
        do 1600 j=1,ncolsa
c 
        d=sa(j)*sb(i)
        if(d .gt. eps) goto 1400
        uacvb(j,i)=0
        goto 1600
c 
 1400 continue
c 
        uacvb(j,i)=uacvb(j,i)/(delta+d)
 1600 continue
 1800 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine leamem(x,y,n)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1)
c 
        do 1200 i=1,n
        y(i)=x(i)
 1200 continue
        return
        end
c 
c 
c 
c 
c 
  
c 
c 
c 
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c          this is the end of the debugging output and the beginning
c          of the actual SVD code
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c 
c 
        subroutine svdpivot(ier,a,n,m,u,v,s,ncols,eps,
     1      w,lw,ltot)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),u(n,1),v(m,1),s(1),w(1)
c 
c        this subroutine constructs the singular value decomposition
c        of the usert-specified matrix a, under the assumption
c        that the rank of a is considerably smaller than its
c        dimensionalities n,m. the output of this subroutine
c        consists of two matrices u,v, and a vector s such that
c 
c                 a = u d v^*                                   (1)
c 
c        with t the diagonal matrix with the vector s on the
c        diagonal, and u, v two matrices whose columns are
c        orthonormal. the reason for the existence of this
c        subroutine is the fact that on exit, the matrices u,v
c        have dimensionalities (n,ncols), (m,ncols), respectively,
c        with ncols the numerical rank of a.
c 
c                   input parameters:
c 
c  a - the matrix to be svd-decomposed; note that it is NOT damaged
c       by this subroutine in any way
c  n,m - dimensionalities of a
c  eps - the accuracy with which the calculation will be performed;
c        should not be less than the machine precision; recommended
c        value: 1.0d-16
c  lw - the number of real *8 words provided by the user in the work
c        array w. must be at least m+n+10 (otherwise, the subroutine
c        can not even estimate the actual memory requirements).
c 
c                   output parameters:
c 
c  ier - error return code
c      ier=0 means successful execution
c      ier=4 means that the amount of space lw allocated supplied
c            in the work array w is insufficient
c      ier=1024 means that the amount of space lw allocated supplied
c            in the work array w is grossly insufficient
c  u - the matrix in (1). note that the matrix u has to be dimensioned
c        u(n,m), but on exit from this subroutine only the first
c        ncols columns of  b  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit b is dimensioned
c        b(n,ncols)
c  v - the second matrix in (1). it is dimensioned v(m,ncols) on
c        exit. on entry the parameter ncols is unknown,
c        so v should be dimensioned v(m,n) by the user. however,
c        on exit only the first ncols columns of v are meaningful
c  s - the array of singular values of a
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices b, v on exit
c  ltot - the total amount of space (in terms of real *8 words)
c        needed by the subroutine in the array w
c 
c                   work arrays:
c 
c  w - must be sufficiently long. 3*n*m+ 2*(min(m,n))**2 +2*n+2*m+100
c        are always sufficient. however, must of the time it is much
c        less than that, and depends on the rank of a. the actual
c        amount of space available in w is communicated to the
c        subroutine via the parameter lw (see above), and if this
c        amount is insufficient, the output error code ier is set
c        accordingly, and the memory requirements of the subroutine
c        are returned in the parameter ltot (see above).
c 
c        . . . using gram-schmidt process with pivoting, decompose
c              the matrix a in the form
c 
c                    a=U  V^*,
c 
c              with  u an orthogonal matrix of minimum rank,
c              and v whatever it wants to be
c 
        if( (lw .gt. m) .and. (lw. gt. n)) goto 1200
        ier=1024
        return
 1200 continue
c 
        ifpivot=1
        call svdpiv0(a,n,m,u,v,ncols,w,eps,ifpivot)
cccc         call prin2('after first svdpiv0, rnorms=*',w,ncols)
c 
c       allocate memory for the subsequent operations
c 
        iw=1
        liw=m*ncols+10
c 
        iu2=iw+liw
        lu2=ncols**2+10
c 
        iv2=iu2+lu2
        lv2=ncols**2+10
c 
        it=iv2+lv2
        lt=lv2
c 
        irnorms=it+lt
        lrnorms=n
        if(m .gt. n) lrnorms=m
        lrnorms=lrnorms+10
c 
        iww=irnorms+lrnorms
        lww=n
        if(m .gt. n) lww=m
        lww=lww+10
c 
        ltot=iww+lww
        if(ltot .le. lw) goto 1400
        ier=4
        return
 1400 continue
c 
c       conclude the construction of the singular value decomposition
c 
        iffirst=0
        call svdpiv1(a,u,v,w(iw),w(it),n,m,ncols,
     1    w(irnorms),eps,w(iv2),w(iu2),w(iww),iffirst)
cccc         call prin2('after svdpiv1, w(it)=*',w(it),ncols)
c 
c       copy the right singluar matrix into the array v
c 
        jj=0
        do 1800 i=1,ncols
        do 1600 j=1,m
        jj=jj+1
        v(j,i)=w(jj)
 1600 continue
 1800 continue
c 
c       now, copy the singular values into the array s
c 
        do 2000 i=1,ncols
        s(i)=w(it+i-1)
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine svdpiv1(a,u,v,w,t,n,m,ncols,rnorms,
     1    eps2,v2,u2,ww,iffirst)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),u(n,1),v(m,1),w(m,1),t(1),rnorms(1),
     1      v2(1),u2(1),ww(1)
c 
c        using gram-schmidt process with pivoting, decompose
c        the matrix a in the form
c 
c          a=U  V^*,
c 
c        with  u an orthogonal matrix of minimum rank,
c        and v whatever it wants to be
c 
        if(iffirst .eq. 0) goto 1200
        ifpivot=1
        call svdpiv0(a,n,m,u,v,ncols,rnorms,eps2,ifpivot)
cccc         call prin2('after first svdpiv0, rnorms=*',rnorms,ncols)
 1200 continue
c 
c        using gram-schmidt process without pivoting, decompose
c        the matrix v in the form
c 
c          v=w t^*,
c 
c        with  w an orthogonal matrix of minimum rank,
c        and t a triangular matrix of dimensionality ncols * ncols
c 
        ifpivot=0
        call svdpiv0(v,m,ncols,w,t,ncols2,rnorms,eps2,ifpivot)
cccc         call prin2('after second svdpiv0, rnorms=*',rnorms,ncols)
c 
c        use jacobi rotations to sv-decompose the matrix t
c 
        ifuv=1
        call d2mudv(t,ncols,u2,v2,eps2,ifuv,numrows)
        call prin2('after d2mudv, t=*',t,ncols)
        call prinf('after d2mudv, numrows=*',numrows,1)
c 
c       multiply the matrix u from the right by u2, and the
c       matrix w from the right by v2^*
c 
        call svdfinm(u,u2,w,v2,n,m,ncols,ww)
        return
        end
c 
c 
c 
c 
c 
        subroutine svdfinm(u,u2,w,v2,n,m,ncols,ww)
        implicit real *8 (a-h,o-z)
        save
        dimension u(n,ncols),w(m,ncols),u2(ncols,ncols),
     1      v2(ncols,ncols),ww(1)
c 
c       multiply u from the right by u2
c 
        do 1600 i=1,n
        do 1400 j=1,ncols
c 
        d=0
        do 1200 k=1,ncols
        d=d+u(i,k)*u2(k,j)
 1200 continue
        ww(j)=d
 1400 continue
        do 1500 j=1,ncols
        u(i,j)=ww(j)
 1500 continue
 1600 continue
c 
c         multiply w from the right by the adjoint of v2
c 
        do 2600 i=1,m
        do 2400 j=1,ncols
c 
        d=0
        do 2200 k=1,ncols
        d=d+w(i,k)*v2(j,k)
 2200 continue
        ww(j)=d
 2400 continue
        do 2500 j=1,ncols
        w(i,j)=ww(j)
 2500 continue
 2600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine svdpiv0(a,n,m,b,v,ncols,rnorms,eps,ifpivot)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),b(n,m),v(m,1),rnorms(1)
c 
c       this matrix applies the compressing gram-schmidt
c       process to the matrix a, obtaining its decomposition
c       in the form
c 
c             a=b v^T,                                        (1)
c 
c       with the matrices b, v having dimensionalities
c       b(n,ncols), v(m,ncols), respectively. the reason for the
c       existence of this subroutine is the hope that the
c       dimensionality ncols, determined by this subroutine,
c       is comsiderably lower than either n or m
c 
c                     input parameters:
c 
c  a - the matrix to be decomposed (note that a is NOT damaged
c       by this subroutine in any way)
c  n,m - the dimensionalities of a
c  eps - the accuracy to which the decomposition (1) will be computed.
c 
c                     output parameters:
c 
c  b - the matrix in (1). note that the matrix b has to be dimensioned
c        b(n,m), but on exit from this subroutine only the first
c        ncols columns of  b  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit b is dimensioned
c        b(n,ncols)
c  v - the second matrix in (1). it is dimensioned v(m,ncols) on
c        exit. on entry the parameter ncols is unknown,
c        so v should be dimensioned v(m,n) by the user. however,
c        on exit only the first ncols columns of v are meaningful
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices b, v on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c 
c        . . . copy the user-supplied matrix a into b
c 
        do 1400 i=1,m
        do 1200 j=1,n
        b(j,i)=a(j,i)
 1200 continue
 1400 continue
c 
c        apply the gram-schmidt proces (with pivoting) to b
c 
         call svdgrm(b,n,m,rnorms,eps,ncols,ifpivot)
c 
c        project the original matrix on the obtained orthogonal
c        columns
c 
        do 2400 i=1,ncols
        do 2200 j=1,m
c 
        call svdscapr(a(1,j),b(1,i),n,prod)
        v(j,i)=prod
 2200 continue
 2400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine svdgrm(b,n,m,rnorms,eps,ncols,ifpivot)
        implicit real *8 (a-h,o-z)
        save
        dimension b(n,m),rnorms(1)
c 
c       this subroutine applies a pivoted double gram-schmidt
c       procedure to the matrix b, resulting in a collection of
c       orthogonal vectors. the number of these vectors is
c       equal to the numerical rank of b, given the precision
c       of computations eps, which is supplied by the user.
c 
c                    input paramneters:
c 
c  b - the matrix to be gram-schmidt'ed. it is destroyed by this
c       subroutine
c  n,m - dimensionalities of the matrix b
c  eps - the machine accuracy
c 
c                     output parameters:
c 
c  b - the matrix of gram-schmidt vectors of the matrix a. note
c        that on exit from this subroutine only the first
c        ncols columns of  b  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit b is dimensioned
c        b(n,ncols)
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrix b on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c 
c       . . . prepare the array of values of norms
c             of columns
c 
        done=1
        dtot=0
        do 1400 i=1,m
        d=0
        do 1200 j=1,n
        d=d+b(j,i)**2
 1200 continue
        rnorms(i)=d
        dtot=dtot+d
 1400 continue
c 
c       . . . conduct gram-schmidt iterations
c 
        thresh=dtot*eps**2
        do 4000 i=1,m
c 
c       find the pivot
c 
         if(ifpivot .eq. 0) goto 2700
        ipivot=i
        rn=rnorms(i)
c 
        do 2200 j=i+1,m
        if(rnorms(j) .le. rn) goto 2200
        rn=rnorms(j)
        ipivot=j
 2200 continue
 2400 continue
c 
c       put the column number ipivot in the i-th place
c 
        do 2600 j=1,n
        d=b(j,i)
        b(j,i)=b(j,ipivot)
        b(j,ipivot)=d
 2600 continue
c 
        d=rnorms(ipivot)
        rnorms(ipivot)=rnorms(i)
        rnorms(i)=d
 2700 continue
c 
c       orthogonalize the i-th column to all preceeding ones
c 
        if(i .eq. 1) goto 2790
        do 2780 j=1,i-1
c 
        call svdscapr(b(1,i),b(1,j),n,d)
c 
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*d
 2770 continue
 2780 continue
 2790 continue
c 
c       normalize the i-th column
c 
        call svdscapr(b(1,i),b(1,i),n,d)
c 
c       if the norm of the remaining part of the matrix
c       is sufficiently small - terminate the process
c 
        if(d .lt. thresh) return
        ncols=i
c 
        d=done/dsqrt(d)
        do 2800 j=1,n
        b(j,i)=b(j,i)*d
 2800 continue
c 
c        orthogonalize everything else to it
c 
        do 3200 j=i+1,m
c 
        if(rnorms(j) .lt. thresh) goto 3200
c 
        call svdscapr(b(1,i),b(1,j),n,d)
c 
         rrn=0
        do 3000 l=1,n
        b(l,j)=b(l,j)-b(l,i)*d
        rrn=rrn+b(l,j)**2
 3000 continue
        rnorms(j)=dsqrt(rrn)
 3200 continue
 3400 continue
 4000 continue
  
 4200 continue
  
  
        return
        end
c 
c 
c 
c 
c 
        subroutine svdscapr(x,y,n,prod)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1)
c 
        prod=0
        do 1200 i=1,n
        prod=prod+x(i)*y(i)
 1200 continue
        return
        end
c 
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c            this is the end of the debuging code and the beginning
c            of the actual SVD routines
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c 
        subroutine d2mudv(a,n,u,v,eps,ifuv,numrows)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),u0(2,2),v0(2,2),u(n,n),v(n,n),
     1      b(2,2),rlam(2)
c 
c        this subroutine uses the classical jacobi iteration
c        to construct the singular value decomposition of a
c        general real matrix, of the form
c 
c          a = u s v,                                        (1)
c 
c        with u,v orthogonal matrices, and s diagonal, with
c        the diagonal elements monotonically decreasing.
c 
c                     input parameters:
c 
c  a - the n * n matrix to be sv-decomposed, is destroyed
c        by the subroutine
c  n - the dimensionality of a
c  eps - the accuracy to which the computations will be performed
c  ifuv - tells the subroutine whether the singular vectors are
c        to be computed, or only singular values
c 
c                     output parameters:
c 
c  a - the first column of a contains the (ordered) singular values
c  u,v - orthogonal matrices in (1)
c  numrows - the number of jacobi iterations (i.e. single rotations)
c        that the algorithm has taken
c 
c       . . . one sweep after another, conduct jacobi
c             iterations
c 
        maxsweep=10
        if(ifuv .eq. 0) goto 1150
        zero=0
        do 1100 i=1,n
        do 1050 j=1,n
        u(j,i)=zero
        v(j,i)=zero
 1050 continue
        u(i,i)=1
        v(i,i)=1
 1100 continue
 1150 continue
        numrows=0
        thresh=dsqrt(eps)
        thresh=dsqrt(thresh)
        thresh=dsqrt(thresh)
        do 3400 iijjkk=1,4
c 
        do 3000 ijk=1,maxsweep
        ifany=0
c 
        do 2000 i=1,n
        do 1800 j=1,n
         if(j .eq. i) goto 1800
         if( (dabs(a(i,j)) .lt. thresh) .and.
     1       (dabs(a(j,i)) .lt. thresh) ) goto 1800
        ifany=1
        numrows=numrows+1
c 
c       construct the svd of the (i,j)-th 2 * 2 submatrix
c 
        b(1,1)=a(i,i)
        b(1,2)=a(i,j)
        b(2,1)=a(j,i)
        b(2,2)=a(j,j)
        call oneuv(b,rlam,u0,v0)
c 
c       transform the big matrix
c 
c       . . . the rows
c 
        do 1200 l=1,n
        di=u0(1,1)*a(i,l)-u0(1,2)*a(j,l)
        dj=-u0(2,1)*a(i,l)+u0(2,2)*a(j,l)
        a(i,l)=di
        a(j,l)=dj
 1200 continue
c 
c       . . . and the columns
c 
        do 1400 l=1,n
        di=v0(1,1)*a(l,i)+v0(1,2)*a(l,j)
        dj=v0(2,1)*a(l,i)+v0(2,2)*a(l,j)
        a(l,i)=di
        a(l,j)=dj
 1400 continue
c 
c 
c       if the user so requested, adjust the matrices u,v
c 
        if(ifuv .eq. 0) goto 1800
c 
        do 1600 l=1,n
c 
c       the matrix u
c 
        di=u0(1,1)*u(l,i)-u0(1,2)*u(l,j)
        dj=-u0(2,1)*u(l,i)+u0(2,2)*u(l,j)
        u(l,i)=di
        u(l,j)=dj
c 
c       . . . and the columns
c 
        di=v0(1,1)*v(i,l)+v0(1,2)*v(j,l)
        dj=v0(2,1)*v(i,l)+v0(2,2)*v(j,l)
        v(i,l)=di
        v(j,l)=dj
 1600 continue
c 
 1800 continue
 2000 continue
c 
      if(ifany .eq. 0) goto 3200
 3000 continue
 3200 continue
        thresh=thresh**2
 3400 continue
c 
c        now, reorder the singular values of a to put them
c        in increasing order
c 
        do 4200 i=1,n
        a(i,1)=a(i,i)
 4200 continue
c 
       call d2msvsrt(a(1,1),n,u,v,a(1,2),ifuv)
        return
        end
c 
c 
c 
c 
c 
        subroutine d2msvsrt(s,n,u,v,iw,ifuv)
        implicit real *8 (a-h,o-z)
        save
        dimension u(n,n),v(n,n),s(1),iw(1)
c 
c        if some of the alleged singular values
c        are negative - change their sign, and that
c        of the corresponding left eigenvectors
c 
        do 1100 i=1,n
        if(s(i) .ge. 0) goto 1100
        s(i)=-s(i)
        if(ifuv .eq. 0) goto 1100
        do 1050 j=1,n
        v(i,j)=-v(i,j)
 1050 continue
 1100 continue
  
c        using the bubble sort, reorder the elements
c        of array s
c 
        do 1200 i=1,n
        iw(i)=i
 1200 continue
c 
        do 2000 i=1,n-1
c 
        ifmoved=0
        do 1600 j=1,n-i
        if(s(j). ge. s(j+1)) goto 1600
        jj=iw(j)
        iw(j)=iw(j+1)
        iw(j+1)=jj
c 
        d=s(j)
        s(j)=s(j+1)
        s(j+1)=d
 1600 continue
 2000 continue
        if(ifuv .eq. 0) return
c 
c       now, reorder the rows of u and columns of v
c 
        do 2400 i=1,n-1
        ii=iw(i)
c 
c       exchange the rows number i and ii in the matrix u
c 
        do 2200 j=1,n
        d=v(i,j)
        v(i,j)=v(ii,j)
        v(ii,j)=d
c 
c       exchange columns number i and ii in the matrix v
c 
        d=u(j,i)
        u(j,i)=u(j,ii)
        u(j,ii)=d
 2200 continue
c 
        do 2300 j=i+1,n
        if(iw(j) .ne. i) goto 2300
        iw(j)=ii
        goto 2350
 2300 continue
 2350 continue
c 
 2400 continue
        return
        end
  
c 
c 
c 
c 
c 
        subroutine oneuv(a,rlam,u,v)
        implicit real *8 (a-h,o-z)
        save
        dimension a(2,2),u(2,2),v(2,2),rlam(2),
     1      b(2,2),w(2,2),vstar(2,2),z(2,2),rlams(2,2)
        data eps/1.0d-10/
c 
c       this subroutine produces the singular value decomposition
c       of a 2 * 2 real matrix a, so that on exit,
c 
c       a = u d v,
c 
c       with u and v orthogonal, and d a diagonal matrix with the
c       vector rlam=(rlam(1),rlam(2)) on the diagonal.
c 
c       . . . simmetrize a
c 
        den=a(2,2)+a(1,1)
        rnum=a(2,1)-a(1,2)
c 
c       if denominator is not too small, use exact formulae
c 
        if(dabs(den) .lt. eps*dabs(rnum)) goto 1200
        tgphi=rnum/den
        phi=datan(tgphi)
        alpha=dcos(phi)
        beta=dsin(phi)
        goto 1400
 1200 continue
c 
c       denominator is too small. use taylor series
c 
        alpha=den/rnum
        beta=1
 1400 continue
c 
c      calculate the simmetrizing matrix and the simmetric one
c 
        w(1,1)=alpha
        w(1,2)=beta
        w(2,1)=-beta
        w(2,2)=alpha
        call prod2(w,a,b)
c 
c       now, diagonalize the simmetrized matrix
c 
        den=b(2,2)-b(1,1)
        rnum=b(1,2)*2
c 
c       if denominator is not too small, use exact formulae
c 
        if(dabs(den) .lt. eps*dabs(rnum)) goto 1600
        tg2phi=-rnum/den
        phi=datan(tg2phi)/2
        alpha=dcos(phi)
        beta=dsin(phi)
        goto 1800
 1600 continue
c 
c       denominator is too small. use taylor expansion
c 
        alpha=dsqrt((1-den/rnum)/2)
        beta=dsqrt((1+den/rnum)/2)
 1800 continue
c 
c       construct the diagonalizing matrix
c       and the resulting diagonal
c 
        v(1,1)=alpha
        v(1,2)=beta
        v(2,1)=-beta
        v(2,2)=alpha
c 
        vstar(1,1)=v(1,1)
        vstar(1,2)=v(2,1)
        vstar(2,1)=v(1,2)
        vstar(2,2)=v(2,2)
c 
        call prod2(v,w,u)
c 
c       finally, compute the resulting diagonal elements
c 
        call prod2(u,a,z)
        call prod2(z,vstar,rlams)
          rlam(1)=rlams(1,1)
          rlam(2)=rlams(2,2)
c 
        u(1,2)=-u(1,2)
        u(2,1)=-u(2,1)
        return
        end
c 
c 
c 
c 
c 
        subroutine prod2(a,b,c)
        implicit real *8 (a-h,o-z)
        save
        dimension a(2,2),b(2,2),c(2,2)
c 
        c(1,1)=a(1,1)*b(1,1)+a(1,2)*b(2,1)
        c(1,2)=a(1,1)*b(1,2)+a(1,2)*b(2,2)
        c(2,1)=a(2,1)*b(1,1)+a(2,2)*b(2,1)
        c(2,2)=a(2,1)*b(1,2)+a(2,2)*b(2,2)
        return
        end
c 
c 
c 
c 
c 
C 
C 
C 
C 
        SUBROUTINE PRINI(IP1,IQ1)
        save
        CHARACTER *1 MES(1), AA(1)
        REAL *4 A(1)
        REAL *8 A2(1)
        REAL *16 A4(1)
        INTEGER *4 IA(1)
        INTEGER *2 IA2(1)
        IP=IP1
        IQ=IQ1
        RETURN
  
C 
C 
C 
C 
C 
        ENTRY PRIN(MES,A,N)
        CALL  MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1200)(A(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1200)(A(J),J=1,N)
 1200 FORMAT(6(2X,E11.5))
         RETURN
C 
C 
C 
C 
        ENTRY PRIN2(MES,A2,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1400)(A2(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1400)(A2(J),J=1,N)
 1400 FORMAT(6(2X,E11.5))
        RETURN
C 
C 
C 
C 
        ENTRY PRINQ(MES,A4,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1500)(A4(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1500)(A4(J),J=1,N)
 1500 FORMAT(6(2X,e11.5))
        RETURN
C 
C 
C 
C 
        ENTRY PRINF(MES,IA,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1600)(IA(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1600)(IA(J),J=1,N)
 1600 FORMAT(10(1X,I7))
        RETURN
C 
C 
C 
C 
        ENTRY PRINF2(MES,IA2,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1600)(IA2(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1600)(IA2(J),J=1,N)
        RETURN
C 
C 
C 
C 
        ENTRY PRINA(MES,AA,N)
        CALL MESSPR(MES,IP,IQ)
 2000 FORMAT(1X,80A1)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,2000)(AA(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,2000)(AA(J),J=1,N)
        RETURN
        END
c 
c 
c 
c 
c 
        SUBROUTINE MESSPR(MES,IP,IQ)
        save
        CHARACTER *1 MES(1),AST
        DATA AST/'*'/
C 
C         DETERMINE THE LENGTH OF THE MESSAGE
C 
        I=0
        DO 1400 I=1,10000
        IF(MES(I).EQ.AST) GOTO 1600
        I1=I
 1400 CONTINUE
 1600 CONTINUE
         IF ( (I1.NE.0) .AND. (IP.NE.0) )
     1     WRITE(IP,1800) (MES(I),I=1,I1)
         IF ( (I1.NE.0) .AND. (IQ.NE.0) )
     1     WRITE(IQ,1800) (MES(I),I=1,I1)
 1800 FORMAT(1X,80A1)
         RETURN
         END
C 
C 
C 
C 
C 
        SUBROUTINE ZTIME(I)
        save
        J=1
        J=7-I+J
CCCC    I=MRUN(J)
        RETURN
        END
c 
c 
c 
c 
c 
        subroutine msgmerge(a,b,c)
        save
        character *1 a(1),b(1),c(1),ast
        data ast/'*'/
c 
        do 1200 i=1,1000
c 
        if(a(i) .eq. ast) goto 1400
        c(i)=a(i)
        iadd=i
 1200 continue
c 
 1400 continue
c 
        do 1800 i=1,1000
c 
        c(iadd+i)=b(i)
        if(b(i) .eq. ast) return
 1800 continue
        return
        end
c 
