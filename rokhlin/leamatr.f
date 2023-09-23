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
        call prin2('before svdpivot, a=*',a,m*n)
        call svdpivot(ier,a,m,n,u,v,s,ncols,eps,
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
        lw1=lw-iwork-2
c 
        call prinf('in leamatr, lw1=*',lw1,1)
c 
        call svdpivot(ier,a,k,l,w(iua),w(iva),w(isa),ncolsa,eps,
     1      w(iwork),lw1,ltot1)
c 
        ltot=ltot1+iwork
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
        call svdpivot(ier,b,m,n,w(iub),w(ivb),w(isb),ncolsb,eps,
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
