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
c 
  
  
        stop
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
        complex *16 a(k,l),x(l,m),b(m,n),c(k,n),
     1      ua(k,1),va(l,1),ub(m,1),vb(n,1),
     2      w(100 000)
c 
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c 
c        construct the test matrices a,b
c 
        call prin2('ima=*',ima,2)
        done=1
        do 1400 i=1,l
        do 1200 j=1,k
        a(j,i)=done/(i**2+j**2+1)
  
        a(j,i)=a(j,i)+done/(ima+j**2-i**2)
  
 1200 continue
cccc        if(i .lt. k) a(i,i)=10+a(i,i)
 1400 continue
c 
        call prin2('a as created*',a,l*k*2)
c 
        done=1
        do 1800 i=1,n
        do 1600 j=1,m
        b(j,i)=done/(i**2+j**2+1)
c 
       b(j,i)=b(j,i)*(2-ima)
 1600 continue
cccc        if(i .lt. m) b(i,i)=10+b(i,i)
 1800 continue
c 
        call prin2('and b as created*',b,m*n*2)
c 
c       create the matrix c
c 
c 
        done=1
        do 2400 i=1,n
        do 2200 j=1,k
        c(j,i)=done/(i**2+j**2+1)
  
         c(j,i)=c(j,i)+c(j,i)**2 * ima
 2200 continue
ccc        if(i .lt. n) c(i,i)=10+c(i,i)
 2400 continue
c 
        call prin2('and c as created*',c,k*n*2)
c 
c        construct the matrix x such that
c 
c                     A X B = C,                                 (1)
c 
        lw=100 000
  
        call cleamatr(ier,a,b,c,k,l,m,n,eps,delta,x,
     1      w,lw,ltot)
  
  
        call prinf('after cleamat0,  ier=*',ier,1)
        call prinf('after cleamat0,  ltot=*',ltot,1)
        call prin2('after cleamat0,  X=*',x,l*m*2)
  
c 
c       multiply the matrices A, X, B back together and compare
c       them with C
c 
        call cleamat(ier,a,k,l,x,l,m,ua)
        call cleamat(ier,ua,k,m,b,m,n,ub)
c 
        call prin2('and C as recomputed from A, B, X*',ub,k*n*2)
        call prin2('again, original C=*',c,k*n*2)
  
        call diffprin(ub,c,k*n,'and the difference*')
  
  
  
  
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
        complex *16 a(1),b(1),c(10 000)
        character *1 mes(1)
c 
        difmax=0
        do 1200 i=1,n
        c(i)=b(i)-a(i)
        dd=cdabs(b(i)-a(i))
        if(difmax .lt. dd) difmax=dd
 1200 continue
        call prin2(mes,c,n*2)
        call prin2('and maximum difference is*',difmax,1)
        return
        end
c 
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code and the  beginning of
c        the actual matrix least squares subroutine
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        subroutine cleamatr(ier,a,b,c,k,l,m,n,eps,delta,x,
     1      w,lw,ltot)
        implicit real *8 (a-h,o-z)
c 
        save
        complex *16  a(k,l),x(l,m),b(m,n),c(k,n),w(1)
c 
c        this subroutine solves in the least squares sense the
c        equation
c 
c                     A(k,l) * X(l,m) * B(m,n) = C(k,n),               (1)
c 
c       where A, B, C are user-specified complex matrices, and X
c       is the matrix to be found. Note that the dimensionalities
c       of the matrices A,B,C, X are as general as could be
c 
c 
c                          input parameters:
c 
c  a,b,c - matrices in (1)
c  k,l,m,n - dimensionalities in (1)
c  eps - the parameter telling the subroutine at which point a
c       a singular value of the matrix is to be declared to be zero
c       and ignored (see subroutine cleasas for details)
c  delta - the "Tikhonov constant" (see subroutine cleasas for details)
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
        call prinf('in cleamatr, lw1=*',lw1,1)
c 
  
cccc        subroutine csvdpiv(ier,a,n,m,u,v,s,ncols,eps,
cccc     1      w,lw,ltot)
  
  
        call csvdpiv(ier,a,k,l,w(iua),w(iva),w(isa),ncolsa,eps2,
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
        call cleamem(w(iva),w(iva2),lva2)
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
cccc        call prinf('in cleamatr, lw2=*',lw2,1)
c 
        eps2=1.0d-13
c 
        call csvdpiv(ier,b,m,n,w(iub),w(ivb),w(isb),ncolsb,eps2,
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
        call cleamem(w(ivb),w(ivb2),lvb2)
c 
        ivb=ivb2
        lvb=lvb2
c 
        iwork=ivb+lvb
        lwork=ncolsa*n+10
c 
        ii=iwork+lwork
        call prinf('in cleamatr, final ii=*',ii,1)
        if(ii .gt. ltot) ltot=ii
        if(ltot .lt. lw) goto 2200
        ier=2
        return
 2200 continue
c 
c         perform the remainder of the solution
c 
        call cleamat0(c,k,l,m,n,eps,delta,x,
     1      w(isa),w(iua),w(iva),w(isb),w(iub),w(ivb),
     2      w(iwork),lw,ltot,ncolsa,ncolsb)
  
  
        return
        end
  
c 
c 
c 
c 
        subroutine cleamat0(c,k,l,m,n,eps,delta,x,
     1      sa,ua,va,sb,ub,vb,work,lw,ltot,ncolsa,ncolsb)
        implicit real *8 (a-h,o-z)
c 
        save
        complex *16 x(l,m),c(k,n),work(1),
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
c       the matrices A,B, provided by the calling subroutine cleamatr
c       (see), so that the equation (1) has the form
c 
c                    UA * SA * VA^* * X * UB * SB * VB^* = C       (2)
c 
c       . . . multiply the matrix c from the left by UA^* ;
c             note that the resulting matrix UAC is dimensioned
c             uac(ncolsa,n)
c 
        call cleamatl(ier,ua,k,ncolsa,c,k,n,work)
c 
c       . . . multiply the matrix UAC by VB from the right;
c             note that the resulting matrix uacvb is
c             dimensioned uacvb(ncolsa,ncolsb)
c 
        call cleamat(ier,work,ncolsa,n,vb,n,ncolsb,x)
c 
c        now, the equation (2) has assumed the form
c 
c                   SA * VA^* * X * UB * SB  = UACVB              (3)
c 
c       . . . multiply (3) by SA^{-1} from the left and
c             by SB^{-1} from the right in the appropriate least
c             squares sense
c 
         call cleasas(sa,ncolsa,sb,ncolsb,x,eps,delta)
c 
c        now, the equation (3) has assumed the form
c 
c                  VA^* * X * UB  = UACVB                     (4)
c 
c       . . . multiply (4) by VA from the left and
c             by UB^* from the right, obtaining the
c             solution X
c 
        call cleamat(ier,va,l,ncolsa,x,ncolsa,ncolsb,work)
c 
        call cleamar(ier,work,l,ncolsb,ub,m,ncolsb,x)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cleamat(ier,a,k,l,b,m,n,c)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(k,l),b(m,n),c(k,n),d
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
        subroutine cleamatl(ier,a,k,l,b,m,n,c)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(k,l),b(m,n),c(l,n),d
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
c       multiply the adjoint of the matrix a by the matrix b getting c
c 
        do 2000 i=1,l
        do 1800 j=1,n
        d=0
        do 1600 kk=1,k
cccc        d=d+a(kk,i)*b(kk,j)
        d=d+dconjg(a(kk,i))*b(kk,j)
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
        subroutine cleamar(ier,a,k,l,b,m,n,c)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(k,l),b(m,n),c(k,m),d
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
c       multiply the matrix a by the the adjoint of the matrix b getting c
c 
        do 2000 i=1,k
        do 1800 j=1,m
        d=0
        do 1600 kk=1,l
cccc        d=d+a(i,kk)*b(j,kk)
        d=d+a(i,kk)*dconjg(b(j,kk))
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
        subroutine cleasas(sa,ncolsa,sb,ncolsb,uacvb,eps,delta)
        implicit real *8 (a-h,o-z)
        save
        complex *16 sa(1),sb(1),uacvb(ncolsa,ncolsb)
c 
c       this subroutine multiplies the user-supplied matrix UACVB
c       by SA^{-1} from the left and by SB^{-1} from the right in
c       the appropriate least squares sense
c 
        do 1800 i=1,ncolsb
        do 1600 j=1,ncolsa
c 
        d=abs(sa(j)*sb(i))
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
        subroutine cleamem(x,y,n)
        implicit real *8 (a-h,o-z)
        save
        complex *16 x(1),y(1)
c 
        do 1200 i=1,n
        y(i)=x(i)
 1200 continue
        return
        end
