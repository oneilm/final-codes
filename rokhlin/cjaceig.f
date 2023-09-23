        implicit real *8 (a-h,o-z)
        dimension a(100 000),w3(100000)
  
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER n'
        READ *,n
        CALL PRInf('n=*',n,1 )
  
  
  
        call eigtest(a,n,w3)
  
  
  
  
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine eigtest(a,n,w3)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,n),b(100000),u(100 000),
     1      b2(100 000),ima,w4(100000),w5(100000),
     2      w1(10000),w(100000),w3(n,n)
c 
        data ima/(0.0d0,1.0d0)/
c 
c       construct the matrix a to be diagonalized
c 
  
        done=1
cccc        ima=1
        eps=2.3
        eps=0.01
        do 1400 i=1,n
        do 1200 j=1,n
c 
        a(j,i)=eps *(i**2 +j**2) +ima*j*eps
cccc        a(j,i)=eps *(i**2 -j**2) +11*j -7*i
 1200 continue
c 
        a(i,i)=-i+ima/i
cccc        a(i,i)=0
 1400 continue
  
        call squarmat(a,n,b)
        call squarmat(a,n,b2)

        call prin2('b as constructed*',b,n*n*2)  
        call prin2('while a is*',a,n*n*2)  

        iw=21
        call matplot(b,n,iw,'as as created*')
  
c 
c       eigendecompose the test matrix
c 

        ifuout=1
        eps=1.0d-12
        call cjaceig(b,n,eps,ifuout,u,nelems,w)
c
        call prin2('after cjaceig, b=*',b,n*2)

cccc        stop

c
c        test them things
c
        do 2400 i=1,n
        do 2200 j=1,n
c
        w3(j,i)=0
 2200 continue
c
cccc        w3(i,i)=w2(i)
        w3(i,i)=b(i)
 2400 continue
c
        call matmul(u,w3,n,n,n,w4)
        call matmua(w4,u,n,n,n,w5)
c
        call prin2('w5=*',w5,n*n*2)
        call prin2('while b2=*',b2,n*n*2)
c
        do 2600 i=1,n*n
c
        w1(i)=w5(i)-b2(i)
 2600 continue
c
        call prin2('and errors are*',w1,n*n*2)
c
        return
        end
c 
c 
c 
c 
c 
        subroutine matmul(a,b,n,m,ncols,c)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,ncols),b(m,ncols),c(n,m),cd
c 
        do 1600 i=1,n
        do 1400 j=1,m
c 
        cd=0
        do 1200 k=1,ncols
        cd=cd+a(i,k)*b(k,j)
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
        subroutine matmua(a,b,n,m,ncols,c)
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
        subroutine squarmat(a,n,b)
        implicit complex *16 (a-h,o-z)
        save
        dimension a(n,n),b(n,n)
c 
        do 1600 i=1,n
        do 1400 j=1,n
        d=0
        do 1200 k=1,n
        d=d+conjg(a(k,i))*a(k,j)
 1200 continue
c 
        b(i,j)=conjg(d)
 1400 continue
 1600 continue
        return
        end

c 
c 
c 
c 
c 
        subroutine matplot(a,n,iw,title)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,n)
        dimension x(100000),y(100000),ww(1000 000)
        character *4 title(1)
c 
c       this subroutine plots the geometry of a mockup of
c       a square matrix
c 
c 
c                input parameters:
c 
c  a - the mockup of the matrix to be plotted. it is an INTEGER *4
c       n \times n matrix; non-zero elements will be plotted.
c  n - the dimensionality of the matrix to be plotted
c  iw - the FORTRAN unit number on which the plot will be written
c  title - the title to be put on the plot
c 
c                output paramaters: none
c 
c                work arrays:
c 
c  x,y - must be at least n**2 real *8 elements each
c 
c       . . . construct the array of coordinates
c 
  
        call trplopen(ww)
  
  
        m=1
        k=n
        eps=1.0d-10
        num=0
        do 1400 i=1,n
        do 1200 j=1,n
        if(abs(a(j,i)) .lt. eps) goto 1200
        num=num+1
        x(num)=i
        y(num)=-j
c 
        xnum=i
        ynum=-j
        call trplpnt(xnum,ynum,ww)
c 
 1200 continue
 1400 continue
cccc          call prin2('before rsplot, x=*',x,num)
cccc          call prin2('before rsplot, y=*',y,num)
  
c 
c       construct the frame
c 
        call plotframe(iw,n,m,k,ww)
c 
c       plot the junk
c 
  
  
        call trplwrt(iw,ww,title)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine plotframe(iw,n,m,k,w)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1)
c 
c        create the vertical boundaries of blocks
c 
cccc         call prinf('in plotframe, n=*',n,1)
cccc         call prinf('in plotframe, m=*',m,1)
cccc         call prinf('in plotframe, k=*',k,1)
  
        do 2000 i=1,k+1
        y1= - ((i-1)*m+0.5)
        y2=y1
        x1=0.5
        x2=n+0.5
  
        call trplline(x1,y1,x2,y2,w)
  
 1200 format(4x,e11.5,4x,e11.5)
cccc        write(iw,1200) x1,y1
cccc        write(iw,1200) x2,y2
 1400 format('        join')
cccc        write(iw,1400)
 2000 continue
  
c 
c        create the horizontal boundaries of blocks
c 
        do 3000 i=1,k+1
        x1=  ((i-1)*m+0.5)
        x2=x1
        y1=-0.5
        y2=-(n+0.5)
  
        call trplline(x1,y1,x2,y2,w)
  
  
cccc        write(iw,1200) x1,y1
cccc        write(iw,1200) x2,y2
cccc        write(iw,1400)
 3000 continue
        return
        end


c
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code and the beginning of the
c        code for the evaluation of eigenvalues and eigenvectors
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c
        subroutine cjaceig(a,n,eps,ifvect,uout,nelems,w)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,n),uout(n,n)
        real *8 w(n,1)
c 
c        This subroutine uses Jacobi iterations to find the
c        spectrum of a Hermotian matrix. It is a primitive
c        and robust code. It is also mercifully short.
c        IMPORTANT! The dimensionality of the matrix must be
c        .ge. 3. Otherwise, a terrible bomb will happen!
c 
c                 Input parameters:
c 
c  a - the matrix whose spectrum is to be found; destroyed by the
c        subroutine
c  n - the dimensionality of a; must be .ge. 3
c  eps - the accuracy (relative) to which the calculations are to
c        be conducted
c  ifvect - the integer parameter telling the subroutine whether it
c        should compute the eigenvectors (or only the eigenvalues)
c        ifvect=1 will cause the eigenvectors to be computed
c        ifvect=0 will suppress the computation of the eigenvectors
c 
c                 Output parameters:
c 
c  a - the first n elements of a contain the (ordered) spectrum of a;
c        please note that the spectrum (while, obviously, real) is
c        stored in the form of complex numbers
c  uout - the eigenvectors of a corresponding to the eigenvalues
c        returned in the first column of a. These are only returned if
c        on input, ifvect was set to 1
c  nelems - the number of jacobi rotations performed by the subroutine;
c        provided solely for the amusement of a sophisticated reader
c 
c
c        . . . tridiagonalize the matrix a
c
        call cjaceig_tridiag(a,n,ifvect,uout)
c
c       reformat the matrix w into a real matrix w2
c
        call cjaceig_convit(a,a,n)
c
c        construct the eigendecomposition of the matrix w2
c
        call jaceig(a,n,eps,ifvect,w,nelems)
c
c        if the user requested that only the eigenvalues
c        be constructed, copy them things to the beginning
c        of array a and get out
c
        if(ifvect .eq. 0) then
            call cjaceig_convit2(a,a,n)
            return
        endif
c
        call cjaceig_mem(a,w(1,n+1),n)
c
c        multiply the complex unitary matrix obtained from the
c        cjaceig_tridiagonalization by the real orthogonal one obtained 
c        from the diagonalization of the real cjaceig_tridiagonal
c        matrix
c
        call cjaceig_matmul(w,uout,n,a)
c
        call cjaceig_mem(a,uout,n*n)

        do 1400 i=1,n
c
        a(i,1)=w(i,n+1)
 1400 continue
c 
        return
        end
c
c
c
c
c
        subroutine cjaceig_convit2(a,b,n)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(n,n)
        complex *16 b(n,n)
c
        do 1200 i=1,n
c
        b(i,4)=a(i,1)
 1200 continue
c
        do 1400 i=1,n
c
        b(i,1)=b(i,4)
 1400 continue
        return
        end
c
c
c
c
c
        subroutine cjaceig_convit(a,b,n)
        implicit real *8 (a-h,o-z)
        save
        real *8 b(n)
        complex *16 a(n)
c
        do 1200 i=1,n*n
c
        b(i)=a(i)
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine cjaceig_matmul(a,b,n,c)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(n,n)
        complex *16 b(n,n),c(n,n),cd
c 
        do 1600 i=1,n
        do 1400 j=1,n
c 
        cd=0
        do 1200 k=1,n
        cd=cd+b(i,k)*a(k,j)
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
        subroutine cjaceig_mem(x,y,n)
        implicit real *8 (a-h,o-z)
        save
        complex *16 x(1),y(1)
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
        subroutine cjaceig_tridiag(a,n,ifuout,uout)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,n),uout(n,n),aa(3),uu(2,2),vv(2,2),
     1      zz
c
c        initialize the matrix uout
c
        if(ifuout .eq. 0) goto 1600
c
        do 1400 i=1,n
        do 1200 j=1,n
c
        uout(j,i)=0
 1200 continue
c
        uout(i,i)=1
 1400 continue
c
 1600 continue
c
        do 2700 i=1,n
        do 2650 j=1,n
c
        d=d+a(j,i)*conjg(a(j,i))
 2650 continue
 2700 continue
c
        d=sqrt(d)
        eps=d*1.0d-37
c
c        use jacobi rotations to tridiagonalize the 
c        matrix A
c
        iw=20
        do 2000 k=1,n-2
c
        do 1800 i=n,k+2,-1
c
        aa(1)=a(k,i-1)
        aa(2)=a(k,i)
c
        call csvd_rotfnd(aa,uu)
c
        call csvd_cols_rotate(uu,a,n,i-1,i)
        if(ifuout .ne. 0) call csvd_cols_rotate(uu,uout,n,i-1,i)
c
        vv(1,1)=conjg(uu(1,1))
        vv(2,2)=conjg(uu(2,2))
        vv(1,2)=conjg(uu(1,2))
        vv(2,1)=conjg(uu(2,1))
c
        call csvd_rows_rotate(vv,a,n,i-1,i)
 1800 continue
 2000 continue
c
c        convert the complex tridiagonal matrix into
c        a real one
c
        do 2600 i=1,n-1
c
        if(abs(a(i+1,i)) .lt. eps) goto 2400
c
        zz=a(i+1,i)/abs(a(i+1,i))
c
        call csvd_colmult(a,n,i+1,zz)
        call csvd_rowmult(a,n,i+1,1/zz)

        if(ifuout .ne. 0) call csvd_colmult(uout,n,i+1,zz)
c     
 2400 continue
 2600 continue
c 
        return
        end
