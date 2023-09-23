        implicit real *8 (a-h,o-z)
        dimension z(40 000),f(20 000),coefs(20 000),
     1      f2(20 000),f3(20 000)
c 
        call prini(6,13)
c 
C 
C       SET ALL PARAMETERS
C 
C 
        PRINT *, 'ENTER n'
        READ *,n
        CALL PRINF('n=*',n,1 )
C 
c       construct the matrix
c 
        a=1
        b=3
        c=1
c 
        call cexpatest(a,b,c,n,z,f,coefs,f2,f3)
c 
        stop
        end
c 
c 
c 
c 
c 
        subroutine expatest(a,b,c,n,z,f,coefs,f2,f3)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension z(2,n,n),f(n,n),coefs(n,n),
     2      f2(n,n),www(100 000),f3(n,n),wwww(100 000)
c 
c        construct the Chebychev discretization of the
c        user-supplied square
c 
        ifinit=1
        call ch2dinit(a,b,c,n,ifinit,z,www)
c 
c 
        call prin2('z as constructed*',z,n*n*2)
  
c 
c       evaluate the test function at all these points
c 
        done=1
        do 1600 i=1,n
        do 1400 j=1,n
        f(j,i)=done/(z(1,j,i)-z(2,j,i))
c 
        f(j,i)=dlog((z(1,j,i)-z(2,j,i))**2)
        f(j,i)=f(j,i)*z(1,j,i)
c 
 1400 continue
 1600 continue
c 
        call prin2('test function as constructed*',f,n**2)
c 
c        decompose the thing into the two-dimensional Chebychev
c        series
c 
        call ch2dexp(f,n,coefs,www)
c 
        call prin2('after ch2dexp, coefs=*',coefs,n**2)
c 
c        evaluate the expansion at the nodes of the Chebichev expansion
c 
c        truncate the expansion after n/2 terms
c 
        do 2400 i=1,n
c 
        if(i .le. n/2) goto 2400
c 
        do 2200 j=1,n
c 
        if(j .le. n/2) goto 2200
c 
        coefs(j,i)=0
 2200 continue
 2400 continue
c 
        call ch2deva(coefs,n,f2,www)
c 
        call prin2('and f2 evaluated from expansion=*',f2,n**2)
c 
c        evaluate the maximum error
c 
        errmax=0
        do 2800 i=1,n
        do 2600 j=1,n
c 
        d=dabs(f(j,i)-f2(j,i))
        if(d .gt. errmax) errmax=d
 2600 continue
 2800 continue
c 
        call prin2('and maximum error is*',errmax,1)
c 
c       now, test the entry evaluating the expansion at a single
c       point
c 
        do 3400 i=1,n
        do 3200 j=1,n
c 
        call ch2deva1(coefs,n,z(1,i,j),z(2,i,j),
     1      f3(i,j),derx1,dery1,wwww)
c 
 3200 continue
 3400 continue
c 
        call prin2('and f recomputed one value at a time*',f3,n*n)
c 
c        evaluate the maximum error
c 
        errmax=0
        do 3800 i=1,n
        do 3600 j=1,n
c 
        d=dabs(f(j,i)-f3(j,i))
        if(d .gt. errmax) errmax=d
 3600 continue
 3800 continue
c 
        call prin2('and maximum error is*',errmax,1)
c 
c        now, test the evaluation of the edrivatives
c 
        i=n/2
        j=n-1
c 
        call ch2deva1(coefs,n,z(1,i,j),z(2,i,j),
     1      f1,derx1,dery1,wwww)
c 
        call prin2('derx1 as evaluated is*',derx1,1)
        call prin2('dery1 as evaluated is*',dery1,1)
c 
        h=0.001
c 
        call ch2deva1(coefs,n,z(1,i,j)+h,z(2,i,j),
     1      f1p,derx1p,dery1p,wwww)
c 
        call ch2deva1(coefs,n,z(1,i,j)-h,z(2,i,j),
     1      f1m,derx1m,dery1m,wwww)
c 
        call prin2('and derx numerically*',(f1p-f1m)/2/h,1)
  
c 
c 
        call ch2deva1(coefs,n,z(1,i,j),z(2,i,j)+h,
     1      f1p,derx1p,dery1p,wwww)
c 
        call ch2deva1(coefs,n,z(1,i,j),z(2,i,j)-h,
     1      f1m,derx1m,dery1m,wwww)
c 
        call prin2('and dery numerically*',(f1p-f1m)/2/h,1)
  
  
  
  
        return
        end
c 
c 
c 
c 
c 
        subroutine cexpatest(a,b,c,n,z,f,coefs,f2,f3)
        implicit real *8 (a-h,o-z)
c 
        save
        dimension z(2,n,n),
     2      www(100 000),wwww(100 000)
c 
        complex *16  f(n,n),coefs(n,n),ima,derx1,dery1,
     2      f2(n,n),f3(n,n),errs(1000)
c 
        data ima/(0.0d0,1.0d0)/
  
c 
c        construct the Chebychev discretization of the
c        user-supplied square
c 
        ifinit=1
        call ch2dinit(a,b,c,n,ifinit,z,www)
c 
c 
        call prin2('z as constructed*',z,n*n*2)
  
c 
c       evaluate the test function at all these points
c 
        done=1
        do 1600 i=1,n
        do 1400 j=1,n
        f(j,i)=z(1,j,i)**2+z(2,j,i)**3
        f(j,i)=z(1,j,i)*ima
        f(j,i)=z(1,j,i)
cccc        f(j,i)=z(1,i,j)
cccc        f(j,i)=z(2,j,i)
c 
  
        f(j,i)=z(1,j,i)**2+z(2,j,i)**3 * ima
  
 1400 continue
 1600 continue
c 
        call prin2('test function as constructed*',f,n**2*2)
c 
c        decompose the thing into the two-dimensional Chebychev
c        series
c 
        call cch2dexp(f,n,coefs,www)
c 
        call prin2('after cch2dexp, coefs=*',coefs,n**2*2)
c 
c        evaluate the expansion at the nodes of the Chebichev expansion
c 
c        truncate the expansion after n/2 terms
c 
        do 2400 i=1,n
c 
        if(i .le. n/2) goto 2400
c 
        do 2200 j=1,n
c 
        if(j .le. n/2) goto 2200
c 
        coefs(j,i)=0
 2200 continue
 2400 continue
c 
        call cch2deva(coefs,n,f2,www)
c 
        call prin2('and f2 evaluated from expansion=*',
     1      f2,n**2*2)
  
cccc        stop
  
c 
c        evaluate the maximum error
c 
        errmax=0
        ii=0
        do 2800 i=1,n
        do 2600 j=1,n
c 
        d=abs(f(j,i)-f2(j,i))
        if(d .gt. errmax) errmax=d
  
        ii=ii+1
        errs(ii)=f(j,i)-f2(j,i)
  
  
 2600 continue
 2800 continue
  
  
        call prin2('and errs=*',errs,n*n*2)
        call prin2('and www=*',www,20)
  
ccc        stop
  
  
c 
c       now, test the entry evaluating the expansion at a single
c       point
c 
        do 3400 i=1,n
        do 3200 j=1,n
c 
        call cch2deva1(coefs,n,z(1,i,j),z(2,i,j),
     1      f3(i,j),derx1,dery1,www)
c 
 3200 continue
 3400 continue
c 
        call prin2('and f recomputed one value at a time*',f3,n*n)
c 
c        evaluate the maximum error
c 
        errmax=0
        do 3800 i=1,n
        do 3600 j=1,n
c 
        d=abs(f(j,i)-f3(j,i))
        if(d .gt. errmax) errmax=d
 3600 continue
 3800 continue
c 
        call prin2('and maximum error is*',errmax,1)
  
        return
        end
c 
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c               this is the end of the debugging code and the beginning
c               of the actual chebychev expansion routine
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c        This file contains 7 user-callable subroutines (or rather,
c        entries). Following is a brief description of the said
c        subrouines.
c 
c 
c   ch2dexp - evaluates the coefficients of the two-dimensional
c        Chebychev series of the user supplied function f
c 
c   ch2deva - evaluates the user-supplied Chebychev expansion
c        at the (tensor product) Chebychev nodes on the square
c 
c   ch2deva1 - evaluates the user-supplied two-dimensional
c        Chebychev expansion at the point (x1,y1), as well as both
c        partial derivatives
c 
c   cch2dexp - evaluates the coefficients of the COMPLEX
c        two-dimensional Chebychev series of the user supplied
c        COMPLEX function f. This is as complex version of the
c        subroutine ch2dexp.
c 
c   cch2deva - evaluates the user-supplied COMPLEX Chebychev
c        expansion at the (tensor product) Chebychev nodes on the
c        square. This is a complex version of the subroutine ch2deva
c 
c   cch2deva1 - evaluates the user-supplied two-dimensional
c        COMPLEX Chebychev expansion at the point (x1,y1), as well
c        as both partial derivatives. This is a complex version
c        of the subroutine ch2deva1
c 
c   ch2dinit - initialization entry point. It prepares data (matrices
c        of transformations, etc.) to be used by the other 6 entries.
c 
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        subroutine cch2dexp(f,n,coefs,w)
        implicit real *8 (a-h,o-z)
        save
        dimension w(1)
        complex *16 f(n,n),coefs(n,n),f1,derx1,dery1
c 
c        this entry evaluates the coefficients of the COMPLEX
c        two-dimensional Chebychev series of the user supplied
c        COMPLEX function f. This is as complex version of the
c        subroutine ch2dexp.
c 
c                       input parameters:
c 
c  f - the values of the function whose two-dimensional Chebychev
c        transform is to be calculated, evaluated at the appropriate
c        (tensor product Chebychev) grid in two dimensions; the points
c        where the values f have been evaluated are returned by the
c        entry ch2dinit of this subroutine (see above). Note that the
c        array f is to be constructed and supplied by the user
c 
c  n - the dimensionality of array f (f is dimensioned n \times n)
c  w - the array that must have been constructed by the entry
c        ch2dinit of this subroutine (see above).
c 
c                       output parameters:
c 
c  coefs - the coefficients of the two-dimensional Chebychev expansion
c        of the user-supplied function f. Note that the array is
c        dimensioned n \times n; in other words, the expansion is of
c        order (n-1)
c 
c 
        iu=w(1)
        iv=w(2)
        ix=w(3)
        iy=w(4)
c 
        call cch2dexp0(n,w(iu),f,coefs,w(ix),w(iy))
c 
        return
c 
c 
c 
c 
        entry cch2deva(coefs,n,f,w)
c 
c        this entry evaluates the user-supplied COMPLEX Chebychev
c        expansion at the (tensor product) Chebychev nodes on the
c        square. This is a complex version of the subroutine ch2deva
c 
c                       input parameters:
c 
c  coefs - the coefficients of the two-dimensional Chebychev expansion
c        to be evaluated. Note that the array is dimensioned n \times n;
c        in other words, the expansion is of order (n-1)
c  n - the dimensionality of array f (f is dimensioned n \times n)
c  w - the array that must have been constructed by the entry
c        ch2dinit of this subroutine (see above).
c 
c                       output parameters:
c 
c  f - the values of the function whose two-dimensional Chebychev
c        expansion is given by the array coefs, evaluated at the
c        appropriate (tensor product Chebychev) grid in two dimensions;
c        the points where the values f have been evaluated are returned
c        by the entry ch2dinit of this subroutine (see above). Note that
c        the array f is to be constructed and supplied by the user
c 
        iu=w(1)
        iv=w(2)
        ix=w(3)
        iy=w(4)
c 
        call cch2deva0(n,w(iv),coefs,f,w(ix),w(iy) )
c 
        return
c 
c 
c 
c 
        entry cch2deva1(coefs,n,x1,y1,f1,derx1,dery1,w)
c 
c 
c        this entry evaluates the user-supplied two-dimensional
c        COMPLEX Chebychev expansion at the point (x1,y1), as well
c        as both partial derivatives. This is a complex version
c        of the subroutine ch2deva1
c 
c                       input parameters:
c 
c  coefs - the coefficients of the two-dimensional Chebychev expansion
c        to be evaluated. Note that the array is
c        dimensioned n \times n; in other words, the expansion is of
c        order (n-1)
c  n - the dimensionality of array coefs (coefs is dimensioned n \times n)
c  (x1,y1) - the coordinates of the point in R^2 where the expansion
c        is to be evaluated
c 
c                       output parameters:
c 
c  f1 - the value of the expansion at the point (x1,y1)
c  derx1 - the x-derivative of the expansion at the point (x1,y1)
c  dery1 - the y-derivative of the expansion at the point (x1,y1)
c 
        ix=w(3)
        iy=w(4)
        idersx=w(5)
        idersy=w(6)
c 
        uua=w(7)
        vva=w(8)
        uub=w(9)
        vvb=w(10)
c 
        xx=x1*uua+vva
        yy=y1*uub+vvb
c 
        call cch2dev10(coefs,xx,yy,w(ix+2),w(iy+2),
     1      w(idersx+2),w(idersy+2),f1,derx1,dery1,n)
c 
        derx1=derx1*uua
        dery1=dery1*uub
        return
        end
c 
c 
c 
c 
c 
        subroutine cch2dev10(coefs,x,y,tsx,tsy,
     1      dersx,dersy,fun,derx,dery,n)
        implicit real *8 (a-h,o-z)
        save
        dimension tsx(1),dersx(1),tsy(1),dersy(1)
        complex *16 coefs(n,n),fun,derx,dery
c 
c        construct the values of chebychev polynomials and
c        their derivatives, of both x and y
c 
        i=0
        tsx(i)=1
        tsy(i)=1
c 
        tsx(1)=x
        tsy(1)=y
c 
        dersx(i)=0
        dersy(i)=0
c 
        dersx(1)=1
        dersy(1)=1
c 
        do 1200 i=1,n
c 
        tsx(i+1)=2*x*tsx(i)-tsx(i-1)
c 
        dersx(i+1)=2*tsx(i)+2*x*dersx(i)-dersx(i-1)
c 
c 
        tsy(i+1)=2*y*tsy(i)-tsy(i-1)
c 
        dersy(i+1)=2*tsy(i)+2*y*dersy(i)-dersy(i-1)
c 
 1200 continue
  
cccc        call prin2('in ch2dev10, tsx =*',tsx(0),n+1)
cccc        call prin2('in ch2dev10, tsy =*',tsy(0),n+1)
  
c 
c        evaluate the expansion and its partial derivatives
c 
        fun=0
        derx=0
        dery=0
c 
        do 1600 i=1,n
        do 1400 j=1,n
c 
        fun=fun+coefs(i,j)*tsx(i-1)*tsy(j-1)
c 
        derx=derx+coefs(i,j)*dersx(i-1)*tsy(j-1)
c 
        dery=dery+coefs(i,j)*tsx(i-1)*dersy(j-1)
 1400 continue
 1600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cch2deva0(n,v,coefs,f,w,w2)
        implicit real *8 (a-h,o-z)
        save
        dimension v(n,n)
        complex *16 f(n,n),coefs(n,n),w(1),w2(1)
c 
c      evaluste a user-specified two-dimensional chebychev series
c      on a chebychev grid
c 
c        . . . horizontally
c 
        do 1200 i=1,n
c 
        call crmatvec(v,coefs(1,i),f(1,i),n)
 1200 continue
c 
c        . . . vertically
c 
        do 1800 i=1,n
c 
        do 1400 j=1,n
        w(j)=f(i,j)
 1400 continue
        call crmatvec(v,w,w2,n)
c 
        do 1600 j=1,n
        f(i,j)=w2(j)
 1600 continue
c 
 1800 continue
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine cch2dexp0(n,u,f,coefs,w,w2)
        implicit real *8 (a-h,o-z)
        save
        dimension u(n,n)
        complex *16 f(n,n),coefs(n,n),w(1),w2(1)
c 
c      decompose the user-supplied function f into a double
c      chebychev series
c 
c        . . . horizontally
c 
        do 1200 i=1,n
c 
        call crmatvec(u,f(1,i),coefs(1,i),n)
 1200 continue
c 
c        . . . vertically
c 
        do 1800 i=1,n
c 
        do 1400 j=1,n
        w(j)=coefs(i,j)
 1400 continue
c 
        call crmatvec(u,w,w2,n)
c 
        do 1600 j=1,n
        coefs(i,j)=w2(j)
 1600 continue
c 
 1800 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine crmatvec(a,x,y,n)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n)
        complex *16 x(1),y(1),d
c 
        do 2000 i=1,n
        d=0
        do 1400 j=1,n
        d=d+a(i,j)*x(j)
 1400 continue
        y(i)=d
 2000 continue
        return
        end
c 
c 
c 
c 
  
        subroutine ch2dinit(a,b,c,n,ifinit,z,w)
        implicit real *8 (a-h,o-z)
        save
        dimension z(2,n,n),w(1),f(n,n),coefs(n,n)
c 
c        this subroutine has four entries. the present entry is
c        for initialization only. The entry ch2dexp (see below)
c        calculates the coefficients of the two-dimensional
c        Chebychev transform of a function (user-supplied) tabulated
c        at (tensor product) Chebychev nodes on a square in R^2.
c        The entry ch2deva (see below) evaluates at the Chebychev nodes
c        on a square the function given by its two-dimensional
c        Chebychev series (obviously, ch2dexp and ch2deva are inverses
c        of each other). Finally, the entry ch2deva1 (see below)
c        evaluates a user-supplied Chebychev expansion at a single
c        point in R^2.
c 
c   IMPORTANT NOTE: This is a "slow" subroutine. It uses explicit
c        matrix forms of both the forward and inverse one-dimensional
c        Chebychev transforms; THIS SUBROUTINE SHOULD NOT BE USED
c        WITH LARGE N!!!!!
c 
c                       input paramaters:
c 
c  (a,b) - the coordinates of the lower left corner of the square in R^2
c  c - the size (length of side) of the square.
c  n - the number of nodes in each direction (note that the total number
c        of nodes on the square is n**2)
c  ifinit - the parameter telling the subroutine whether it should
c        initialize the Chebychev transform routine, which can be pricey;
c     ifinit=1 means that the subroutine will be initialized
c     ifinit=0 means that the subroutine has already been initialized;
c        also note that the first 2*n**2+6*n+100 elements of array w
c        must not have been changed between this call and the preceding
c        call to this subroutine (when ifinit was equal to 1)
c 
c                       output parameters:
c 
c  z - the n**2 (tensor product) Chebychev nodes on the user-specified
c        square
c  w - the array to be used by the entries ch2deva, ch2dexp, ch2deva1
c        (see below). The subroutine uses 2*n**2+6*n+100 elements of
c        the array w. Note that these elements should not be changed
c        between the call to this subroutine and the subsequent calls
c        to ch2deva, ch2dexp, ch2deva1.
c 
c       . . . allocate memory for the subroutine cheb2d0 that
c             will discretiza the square into an n \times n
c             chebychev grid, and initialize the Chebychev transform
c             subroutine
c 
        it=21
        lt=n+1
c 
        ix=it+lt
        lx=n*2+10
c 
        iy=ix+lx
        ly=n*2+10
c 
        iwhts=iy+ly
        lwhts=n+1
c 
        iu=iwhts+lwhts
        lu=n*n+10
c 
        iv=iu+lu
        lv=n*n+10
c 
        idersx=iv+lv
        ldersx=n*2+10
c 
        idersy=idersx+ldersx
        ldersy=n*2+10
c 
        w(1)=iu+0.1
        w(2)=iv+0.1
        w(3)=ix+0.1
        w(4)=iy+0.1
        w(5)=idersx+0.1
        w(6)=idersy+0.1
c 
c        construct the discretization of the square and
c        initialize the Chebychev transform
c 
        call cheb2d0(a,b,c,n,w(it),w(iu),w(iv),
     1      w(iwhts),w(ix),w(iy),z,ifinit)
c 
        a1=a
        a2=a1+c
c 
        b1=b
        b2=b1+c
c 
        uua=2/(a2-a1)
        vva=1-uua*a2
c 
        uub=2/(b2-b1)
        vvb=1-uub*b2
c 
        w(7)=uua
        w(8)=vva
        w(9)=uub
        w(10)=vvb
c 
        return
c 
c 
c 
c 
        entry ch2dexp(f,n,coefs,w)
c 
c        this entry evaluates the coefficients of the two-dimensional
c        Chebychev series of the user supplied function f
c 
c                       input parameters:
c 
c  f - the values of the function whose two-dimensional Chebychev
c        transform is to be calculated, evaluated at the appropriate
c        (tensor product Chebychev) grid in two dimensions; the points
c        where the values f have been evaluated are returned by the
c        entry ch2dinit of this subroutine (see above). Note that the
c        array f is to be constructed and supplied by the user
c 
c  n - the dimensionality of array f (f is dimensioned n \times n)
c  w - the array that must have been constructed by the entry
c        ch2dinit of this subroutine (see above).
c 
c                       output parameters:
c 
c  coefs - the coefficients of the two-dimensional Chebychev expansion
c        of the user-supplied function f. Note that the array is
c        dimensioned n \times n; in other words, the expansion is of
c        order (n-1)
c 
        call ch2dexp0(n,w(iu),f,coefs,w(ix),w(iy))
c 
        return
c 
c 
c 
c 
        entry ch2deva(coefs,n,f,w)
c 
c        this entry evaluates the user-supplied Chebychev expansion
c        at the (tensor product) Chebychev nodes on the square
c 
c                       input parameters:
c 
c  coefs - the coefficients of the two-dimensional Chebychev expansion
c        to be evaluated. Note that the array is dimensioned n \times n;
c        in other words, the expansion is of order (n-1)
c  n - the dimensionality of array f (f is dimensioned n \times n)
c  w - the array that must have been constructed by the entry
c        ch2dinit of this subroutine (see above).
c 
c                       output parameters:
c 
c  f - the values of the function whose two-dimensional Chebychev
c        expansion is given by the array coefs, evaluated at the
c        appropriate (tensor product Chebychev) grid in two dimensions;
c        the points where the values f have been evaluated are returned
c        by the entry ch2dinit of this subroutine (see above). Note that
c        the array f is to be constructed and supplied by the user
c 
        call ch2deva0(n,w(iv),coefs,f,w(ix),w(iy) )
c 
        return
c 
c 
c 
c 
        entry ch2deva1(coefs,n,x1,y1,f1,derx1,dery1,w)
c 
c        this entry evaluates the user-supplied two-dimensional
c        Chebychev expansion at the point (x1,y1), as well as both
c        partial derivatives
c 
c                       input parameters:
c 
c  coefs - the coefficients of the two-dimensional Chebychev expansion
c        to be evaluated. Note that the array is
c        dimensioned n \times n; in other words, the expansion is of
c        order (n-1)
c  n - the dimensionality of array coefs (coefs is dimensioned n \times n)
c  (x1,y1) - the coordinates of the point in R^2 where the expansion
c        is to be evaluated
c 
c                       output parameters:
c 
c  f1 - the value of the expansion at the point (x1,y1)
c  derx1 - the x-derivative of the expansion at the point (x1,y1)
c  dery1 - the y-derivative of the expansion at the point (x1,y1)
c 
        xx=x1*uua+vva
        yy=y1*uub+vvb
c 
        call ch2dev10(coefs,xx,yy,w(ix+2),w(iy+2),
     1      w(idersx+2),w(idersy+2),f1,derx1,dery1,n)
c 
        derx1=derx1*uua
        dery1=dery1*uub
        return
        end
c 
c 
c 
c 
c 
        subroutine ch2dev10(coefs,x,y,tsx,tsy,
     1      dersx,dersy,fun,derx,dery,n)
        implicit real *8 (a-h,o-z)
        save
        dimension coefs(n,n),tsx(1),dersx(1),tsy(1),dersy(1)
c 
c        construct the values of chebychev polynomials and
c        their derivatives, of both x and y
c 
        i=0
        tsx(i)=1
        tsy(i)=1
c 
        tsx(1)=x
        tsy(1)=y
c 
        dersx(i)=0
        dersy(i)=0
c 
        dersx(1)=1
        dersy(1)=1
c 
        do 1200 i=1,n
c 
        tsx(i+1)=2*x*tsx(i)-tsx(i-1)
c 
        dersx(i+1)=2*tsx(i)+2*x*dersx(i)-dersx(i-1)
c 
c 
        tsy(i+1)=2*y*tsy(i)-tsy(i-1)
c 
        dersy(i+1)=2*tsy(i)+2*y*dersy(i)-dersy(i-1)
c 
 1200 continue
  
cccc        call prin2('in ch2dev10, tsx =*',tsx(0),n+1)
cccc        call prin2('in ch2dev10, tsy =*',tsy(0),n+1)
  
c 
c        evaluate the expansion and its partial derivatives
c 
        fun=0
        derx=0
        dery=0
c 
        do 1600 i=1,n
        do 1400 j=1,n
c 
        fun=fun+coefs(i,j)*tsx(i-1)*tsy(j-1)
c 
        derx=derx+coefs(i,j)*dersx(i-1)*tsy(j-1)
c 
        dery=dery+coefs(i,j)*tsx(i-1)*dersy(j-1)
 1400 continue
 1600 continue
c 
        return
        end
c 
c 
c 
c 
c 
  
        subroutine cheb2d0(a,b,c,n,t,u,v,whts,x,y,z,ifinit)
        implicit real *8 (a-h,o-z)
        save
        dimension t(1),u(n,n),v(n,n),whts(1),x(1),y(1),
     1      z(2,n,n)
c 
c       construct the Chebychev nodes and the matrices
c       for the forward and inverse Chebychev transforms
c 
        if(ifinit .eq. 0) goto 1100
        itype=2
        call chebexps(itype,n,t,u,v,whts)
c 
 1100 continue
c 
c       construct the discretization of the square
c 
        a1=a
        a2=a1+c
c 
        b1=b
        b2=b1+c
c 
        ua=(a2-a1)/2
        va=(a2+a1)/2
c 
        ub=(b2-b1)/2
        vb=(b2+b1)/2
c 
        do 1200 i=1,n
c 
        x(i)=ua*t(i)+va
        y(i)=ub*t(i)+vb
 1200 continue
c 
c 
c       evaluate the test function at all these points
c 
        done=1
        do 1600 i=1,n
        do 1400 j=1,n
        z(1,j,i)=x(j)
        z(2,j,i)=y(i)
 1400 continue
 1600 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine ch2dexp0(n,u,f,coefs,w,w2)
        implicit real *8 (a-h,o-z)
        save
        dimension u(n,n),f(n,n),coefs(n,n),w(1),w2(1)
c 
c      decompose the user-supplied function f into a double
c      chebychev series
c 
c        . . . horizontally
c 
        do 1200 i=1,n
c 
        call matvec(u,f(1,i),coefs(1,i),n)
 1200 continue
c 
c        . . . vertically
c 
        do 1800 i=1,n
c 
        do 1400 j=1,n
        w(j)=coefs(i,j)
 1400 continue
        call matvec(u,w,w2,n)
c 
        do 1600 j=1,n
        coefs(i,j)=w2(j)
 1600 continue
c 
 1800 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine ch2deva0(n,v,coefs,f,w,w2)
        implicit real *8 (a-h,o-z)
        save
        dimension v(n,n),f(n,n),coefs(n,n),w(1),w2(1)
c 
c      evaluste a user-specified two-dimensional chebychev series
c      on a chebychev grid
c 
c        . . . horizontally
c 
        do 1200 i=1,n
c 
        call matvec(v,coefs(1,i),f(1,i),n)
 1200 continue
c 
c        . . . vertically
c 
        do 1800 i=1,n
c 
        do 1400 j=1,n
        w(j)=f(i,j)
 1400 continue
        call matvec(v,w,w2,n)
c 
        do 1600 j=1,n
        f(i,j)=w2(j)
 1600 continue
c 
 1800 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine matvec(a,x,y,n)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),x(1),y(1)
c 
        do 2000 i=1,n
        d=0
        do 1400 j=1,n
        d=d+a(i,j)*x(j)
 1400 continue
        y(i)=d
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine chebexps(itype,n,x,u,v,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),whts(1),u(n,n),v(n,n)
c 
c         this subroutine constructs the chebychev nodes
c         on the interval [-1,1], and the weights for the
c         corresponding order n quadrature. it also constructs
c         the matrix converting the coefficients
c         of a chebychev expansion into its values at the n
c         chebychev nodes. no attempt has been made to
c         make this code efficient, but its speed is normally
c         sufficient, and it is mercifully short.
c 
c                 input parameters:
c 
c  itype - the type of the calculation to be performed
c          itype=0 means that only the chebychev nodes are
c                  to be constructed.
c          itype=1 means that only the nodes and the weights
c                  are to be constructed
c          itype=2 means that the nodes, the weights, and
c                  the matrices u, v are to be constructed
c  n - the number of chebychev nodes and weights to be generated
c 
c                 output parameters:
c 
c  x - the order n chebychev nodes - computed independently
c          of the value of itype.
c  u - the n*n matrix converting the  values at of a polynomial of order
c         n-1 at n chebychev nodes into the coefficients of its
c         chebychev expansion - computed only in itype=2
c  v - the n*n matrix converting the coefficients
c         of an n-term chebychev expansion into its values at
c         n chebychev nodes (note that v is the inverse of u)
c          - computed only in itype=2
c  whts - the corresponding quadrature weights - computed only
c         if itype .ge. 1
c 
c       . . . construct the chebychev nodes on the interval [-1,1]
c 
        ZERO=0
        DONE=1
        pi=datan(done)*4
        h=pi/(2*n)
ccc        do 1200 i=n,1,-1
        do 1200 i=1,n
        t=(2*i-1)*h
        x(n-i+1)=dcos(t)
1200  CONTINUE
cccc        call prin2('chebychev nodes as constructed*',x,n)
c 
        if(itype .eq. 0) return
c 
c        construct the weights of the quadrature
c        formula based on the chebychev nodes,
c        and also the matrix of the chebychev transform
c 
c        . . . construct the first two rows of the matrix
c 
         if(itype .le. 1) goto 1350
        do 1300 i=1,n
        u(1,i)=1
        u(2,i)=x(i)
 1300 continue
 1350 continue
c 
c       construct all quadrature weights and the rest of the rows
c 
        do 2000 i=1,n
c 
c       construct the weight for the i-th node
c 
        Tjm2=1
        Tjm1=x(i)
        whts(i)=2
c 
        ic=-1
        do 1400 j=2,n-1
c 
c       calculate the T_j(x(i))
c 
        Tj=2*x(i)*Tjm1-Tjm2
c 
        if(itype .eq. 2) u(j+1,i)=tj
c 
        tjm2=tjm1
        tjm1=tj
c 
c       calculate the contribution of this power to the
c       weight
c 
  
        ic=-ic
        if(ic .lt. 0) goto 1400
        rint=-2*(done/(j+1)-done/(j-1))
        whts(i)=whts(i)-rint*tj
ccc        whts(i)=whts(i)+rint*tj
 1400 continue
        whts(i)=whts(i)/n
 2000 continue
           if(itype .ne. 2) return
c 
c        now, normalize the matrix of the chebychev transform
c 
        do 3000 i=1,n
c 
        d=0
        do 2200 j=1,n
        d=d+u(i,j)**2
 2200 continue
        d=done/dsqrt(d)
        do 2400 j=1,n
        u(i,j)=u(i,j)*d
 2400 continue
 3000 continue
c 
c        now, rescale the matrix
c 
        ddd=2
        ddd=dsqrt(ddd)
        dd=n
        dd=done/dsqrt(dd/2)
        do 3400 i=1,n
        do 3200 j=1,n
        u(j,i)=u(j,i)*dd
 3200 continue
        u(1,i)=u(1,i)/ddd
 3400 continue
c 
c        finally, construct the matrix v, converting the values at the
c        chebychev nodes into the coefficients of the chebychev
c        expansion
c 
        dd=n
        dd=dd/2
        do 4000 i=1,n
        do 3800 j=1,n
        v(j,i)=u(i,j)*dd
 3800 continue
 4000 continue
c 
        do 4200 i=1,n
        v(i,1)=v(i,1)*2
 4200 continue
        return
        end
c 
c 
c 
c 
c 
      SUBROUTINE CHFUNDER(X,VAL,der,TEXP,N)
C 
C     This subroutine computes the value and the derivative
c     of a chebychev expansion with coefficients TEXP
C     at point X in interval [-1,1]
C 
c                input parameters:
c 
C     X = evaluation point
C     TEXP = expansion coefficients
C     N  = order of expansion
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c 
c                output parameters:
c 
C     VAL = computed value
C     der = computed value of the derivative
C 
      IMPLICIT REAL *8 (A-H,O-Z)
        save
      REAL *8 TEXP(1)
C 
        done=1
        tjm2=1
        tjm1=x
        derjm2=0
        derjm1=1
c 
        val=texp(1)*tjm2+texp(2)*tjm1
        der=texp(2)
c 
      DO 600 J = 2,N
c 
        tj=2*x*tjm1-tjm2
        val=val+texp(j+1)*tj
c 
        derj=2*tjm1+2*x*derjm1-derjm2
        der=der+texp(j+1)*derj
c 
        tjm2=tjm1
        tjm1=tj
        derjm2=derjm1
        derjm1=derj
 600   CONTINUE
c 
      RETURN
      END
  
  
