      implicit real *8 (a-h,o-z)
      dimension vert1(2),vert2(2),vert3(2)
      dimension ts(2,1000000), whts(1000000)
      dimension pols1(10 000000),  pols2(10 000000), pols3(10 000000)
      dimension dersx2(10 000000),  dersy2(10 000000)
      dimension rnorms(1000), z(2), z1(2), diffs(1000), errs(1000)
      dimension w(10000)
  
      call prini(6,13)
  
C 
C       SET ALL PARAMETERS
C 
      PRINT *, 'ENTER mmax'
      READ *,mmax
      CALL PRINf('mmax=*',mmax,1)
  
c      PRINT *, 'ENTER z'
c      READ *,z(1),z(2)
c      CALL PRIN2('z=*',z,2)
  
      z(1)=-1.0d0/3
      z(2)=-1.0d0/3
  
      z(1)=0
      z(2)=0
  
      CALL PRIN2('z=*',z,2)
  
      n=10 000
      n=1
      do i=1,n
      call ortho2eva(mmax,z,pols2,w)
      enddo
      call prin2('pols=*', pols2, (mmax+1)*(mmax+2)/2)
  
       call ortho2eva3(mmax,z,pols2,dersx2,dersy2,w)
c      call prin2('pols=*', pols2, (mmax+1)*(mmax+2)/2)
c      call prin2('dersx=*', dersx2, (mmax+1)*(mmax+2)/2)
c      call prin2('dersy=*', dersy2, (mmax+1)*(mmax+2)/2)
  
      h=1d-14
      h=1d-7
  
c 
c     ... check the derivative in x
c 
      z1(1)=z(1)-h
      z1(2)=z(2)
      call ortho2eva3(mmax,z1,pols1,dersx2,dersy2,w)
ccc      call prin2('pols=*', pols1, (mmax+1)*(mmax+2)/2)
      z1(1)=z(1)
      z1(2)=z(2)
      call ortho2eva3(mmax,z1,pols2,dersx2,dersy2,w)
ccc      call prin2('pols=*', pols2, (mmax+1)*(mmax+2)/2)
      z1(1)=z(1)+h
      z1(2)=z(2)
      call ortho2eva3(mmax,z1,pols3,dersx2,dersy2,w)
ccc      call prin2('pols=*', pols3, (mmax+1)*(mmax+2)/2)
      z1(1)=z(1)
      z1(2)=z(2)
      call ortho2eva3(mmax,z1,pols2,dersx2,dersy2,w)
cccc      call prin2('dersx=*', dersx2, (mmax+1)*(mmax+2)/2)
c 
      do i=1,(mmax+1)*(mmax+2)/2
         diffs(i) = (pols3(i)-pols1(i))/(2*h)
         errs(i)=diffs(i)-dersx2(i)
      enddo
c 
cccc      call prin2('aprox dersx=*', diffs, (mmax+1)*(mmax+2)/2)
      call prin2('and error in dersx=*', errs, (mmax+1)*(mmax+2)/2)
  
c 
c     ... check the derivative in y
c 
      z1(1)=z(1)
      z1(2)=z(2)-h
      call ortho2eva3(mmax,z1,pols1,dersx2,dersy2,w)
ccc      call prin2('pols=*', pols1, (mmax+1)*(mmax+2)/2)
      z1(1)=z(1)
      z1(2)=z(2)
      call ortho2eva3(mmax,z1,pols2,dersx2,dersy2,w)
ccc      call prin2('pols=*', pols2, (mmax+1)*(mmax+2)/2)
      z1(1)=z(1)
      z1(2)=z(2)+h
      call ortho2eva3(mmax,z1,pols3,dersx2,dersy2,w)
ccc      call prin2('pols=*', pols3, (mmax+1)*(mmax+2)/2)
      z1(1)=z(1)
      z1(2)=z(2)
      call ortho2eva3(mmax,z1,pols2,dersx2,dersy2,w)
cccc      call prin2('dersy=*', dersy2, (mmax+1)*(mmax+2)/2)
c 
      do i=1,(mmax+1)*(mmax+2)/2
         diffs(i) = (pols3(i)-pols1(i))/(2*h)
         errs(i)=diffs(i)-dersy2(i)
      enddo
c 
cccc      call prin2('approx dersy=*', diffs, (mmax+1)*(mmax+2)/2)
      call prin2('and error in dersy=*', errs, (mmax+1)*(mmax+2)/2)
  
c      stop
 5000 continue
  
c 
c     right triangle
c 
c      vert1(1)=-1
c      vert1(2)=-1
c      vert2(1)=-1
c      vert2(2)= 1
c      vert3(1)= 1
c      vert3(2)=-1
c 
c     standard triangle
c 
      vert1(1)=-1
      vert1(2)=-1/dsqrt(3.0d0)
      vert2(1)= 0
      vert2(2)= 2/dsqrt(3.0d0)
      vert3(1)= 1
      vert3(2)=-1/dsqrt(3.0d0)
  
  
ccc      np=10+1
ccc      np=30+1
      np=40+1
ccc     np=50+1
ccc     np=100+1
  
  
      ifinit=1
ccc      call triagauc(np,vert1,vert2,vert3,ts,whts,ifinit,w)
      call triagauc(np,vert2,vert3,vert1,ts,whts,ifinit,w)
ccc      call triagauc(np,vert3,vert1,vert2,ts,whts,ifinit,w)
  
      kk=np**2
  
c      do i=1,kk
c         write(99,*) ts(1,i), ts(2,i), whts(i)
c      enddo
  
c 
c     ... check the orthogonality
c 
      nvals=(mmax+1)*(mmax+2)/2
      call test1(kk,ts,whts,pols1,pols2,mmax,nvals)
  
      stop
      end
  
  
  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  
      subroutine test1(kk,zs,ws,pols1,pols2,mmax,nvals)
      implicit real *8 (a-h,o-z)
        save
      dimension arr(1), w(10000)
      dimension zs(2,kk), ws(kk), pols1(nvals,kk), pols2(kk,nvals)
c 
c       ... evaluate the orthogonal polynomials up to porder mmax
c 
      do i=1,kk
         call ortho2eva(mmax,zs(1,i),pols1(1,i),w)
      enddo
      do i=1,kk
         do j=1,nvals
            pols2(i,j) = pols1(j,i)
         enddo
      enddo
      n = nvals
c 
c       ... test the orthogonality
c 
      s1 = 0
      s2 = 0
      do i=1,nvals
      do j=1,nvals
      call ortscap(pols2(1,i),pols2(1,j),ws,kk,d)
c      write(*,*) i,j,d
c      if(i.eq.j) write(*,*) i,j,d
      if(i.eq.j) s1 = s1 + (d-1)**2
      if(i.ne.j) s2 = s2 + (d)**2
      enddo
      enddo
      s1 = dsqrt(s1/n)
      s2 = dsqrt(s2/(n*(n-1)))
      call prin2('and error in norms=*', s1, 1)
      call prin2('and error in products=*', s2, 1)
      return
      end
  
  
      subroutine ortscap(x,y,w,n,d)
      implicit real *8 (a-h,o-z)
        save
      dimension x(1),y(1),w(1)
c 
      d=0
      do 2000 i=1,n
         d=d+x(i)*y(i)*w(i)
 2000 continue
      return
      end
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c       This is the end of the debugging routines and the beginning
c       of the code for the evaluation of orthogonal polynomials
c       on the standard triangle
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c       This file contains a set of subroutines for the handling
c       of the evaluation of othogonal polynomials on the standard
c       triangle.  It contains 2 user-callable subroutines.
c 
c   ortho2eva - evaluates at the user-supplied point z
c       a collection of polynomials (of x,y) orthogonal on the
c       standard triangle
c 
c   ortho2eva3 - evaluates at the user-supplied point z a collection
c       of polynomials (of x,y) orthogonal on the standard triangle,
c       together with their derivatives with respect to x and y.
  
c 
c 
      subroutine ortho2eva(mmax,z,pols,w)
c 
c     This subroutine evaluates at the user-supplied point z
c     a collection of polynomials (of x,y) orthogonal on the
c     standard triangle.
c 
c     The "standard" triangle is the triangle with the
c     vertices
c 
c     (0,2/sqrt(3)), (-1,sqrt(1/sqrt(3)), (1,sqrt(1/sqrt(3)).       (1)
c 
c     The polynomials evaluated by this subroutine are all
c     orthogonal polynomials up to order mmax, arranged in the
c     increasing order.
c 
c     This subroutine is based on the Koornwinder's representation
c     of the orthogonal polynomials on the right triangle
c 
c     (-1,-1), (-1,1), (1,-1)                                       (2)
c 
c     given by
c 
c     K_mn(x,y) = P_m ((2*x+1+y)/(1-y)) ((1-y)/2)^m P_n^{2m+1,0} (y)
c 
c     where P_m are the Legendge polynomials or order m
c     and P_n^{2m+1,0} are the Jacobi polynomials of order n with
c     the parameters alpha=2*m+1 and beta=0.
c 
c                   Input parameters:
c 
c  mmax - the maximum order to which the polynomials are to be evaluated;
c  z - the location in R^2 where the polynomials are to be evaluated;
c       normally, expected to be inside (including boundary) the
c       standard triangle (1) above.
c  w - the work array. must be at least mmax+1+(mmax+1)**2 long
c 
c                   Output parameters:
c 
c  pols - the orthogonal polynomials evaluated at the point z
c       ( (mmax+1)*(mmax+2)/2 of them things)
c 
c       . . . allocate memory for the work arrays
c 
      implicit real *8 (a-h,o-z)
        save
      dimension z(2), pols(1), w(1)
      iw1=1
      iw2=iw1+mmax+1
      iw3=iw2+(mmax+1)**2
      call ortho2eva0(mmax,z,pols,w(iw1),w(iw2))
      return
      end
c 
c 
c 
c 
c 
      subroutine ortho2eva0(mmax,z,pols,f1,f2)
c 
c     evaluate the orthonormal polynomials on the triangle
c 
      implicit real *8 (a-h,o-z)
        save
      dimension z(2), pols(1), f1(1), f2(1)
      data init/1/
c 
      if( init.eq.1 ) then
         zero=0
         sqrt2=dsqrt(2.0d0)
         sqrt3=dsqrt(3.0d0)
         init=0
         r11=-1.0d0/3.0d0
         r12=-1.0d0/sqrt3
         r21=1.0d0/9.0d0*(-sqrt3)*sqrt3
         r22=1.0d0/9.0d0*(6.0d0)*sqrt3
      endif
c 
      a=z(1)
      b=z(2)
c 
c     ... map the standard triangle to the right
c     triangle with the vertices (-1,-1), (1,-1), (-1,1)
c 
      x=r11+r12*b+a
      y=r21+r22*b
c 
c     ... evaluate the Koornwinder's polynomials
c     the via three terms recursion
c 
      par1=(2*x+1+y)/2
      par2=(1-y)/2
      call klegeypols(par1,par2,mmax,f1)
      do 1200 m=0,mmax
         par1=2*m+1
ccc         call kjacopols(y,par1,zero,mmax,f2(1+m*(mmax+1)))
         call kjacopols(y,par1,zero,(mmax-m),f2(1+m*(mmax+1)))
 1200 continue
	 
      kk=0
      do 2200 m=0,mmax
         do 2000 n=0,m
         kk=kk+1
c 
c     ... evaluate the polynomial (m-n, n)
c 
         pols(kk)=f1(m-n+1)*f2(n+1+(m-n)*(mmax+1))
c 
c     ... and normalize it
c 
         scale=dsqrt((1.0d0+(m-n)+n)*(1.0d0+(m-n)+(m-n))/sqrt3)
cccc         scale=1
         pols(kk)=pols(kk)*scale
 2000	 continue
 2200 continue
      return
      end
c 
c 
c 
c 
c 
      subroutine ortho2eva3(mmax,z,pols,dersx,dersy,w)
c 
c     This subroutine evaluates at the user-supplied point z
c     a collection of polynomials (of x,y) orthogonal on the
c     standard triangle, together with their derivatives with
c     respect to x and y.
c 
c     The "standard" triangle is the triangle with the
c     vertices
c 
c     (0,2/sqrt(3)), (-1,sqrt(1/sqrt(3)), (1,sqrt(1/sqrt(3)).       (1)
c 
c     The polynomials evaluated by this subroutine are all
c     orthogonal polynomials up to order mmax, arranged in the
c     increasing order.
c 
c     This subroutine is based on the Koornwinder's representation
c     of the orthogonal polynomials on the right triangle
c 
c     (-1,-1), (-1,1), (1,-1)                                       (2)
c 
c     given by
c 
c     K_mn(x,y) = P_m ((2*x+1+y)/(1-y)) ((1-y)/2)^m P_n^{2m+1,0} (y)
c 
c     where P_m are the Legendge polynomials or order m
c     and P_n^{2m+1,0} are the Jacobi polynomials of order n with
c     the parameters alpha=2*m+1 and beta=0.
c 
c                   Input parameters:
c 
c  mmax - the maximum order to which the polynomials are to be evaluated;
c  z - the location in R^2 where the polynomials are to be evaluated;
c       normally, expected to be inside (including boundary) the
c       standard triangle (1) above.
c  w - the work array. must be at least mmax+1+(mmax+1)**2 long
c 
c                   Output parameters:
c 
c  pols - the orthogonal polynomials evaluated at the point z
c       ( (mmax+1)*(mmax+2)/2 of them things)
c  dersx - the derivatives with respect to x of the polynomials
c       returned in array pols
c  dersy - the derivatives with respect to y of the polynomials
c       returned in array pols
c 
c     NOTE: this subroutine will fail to evaluate the derivatives
c           at the vertex (0,2/sqrt(3)).
c 
c       . . . allocate memory for the work arrays
c 
      implicit real *8 (a-h,o-z)
        save
      dimension z(2), pols(1), w(1)
      iw1=1
      iw2=iw1+mmax+1
      iw3=iw2+(mmax+1)**2
      iw4=iw3+mmax+1
      iw5=iw4+mmax+1
      iw6=iw5+(mmax+1)**2
      iw7=iw6+(mmax+1)**2
      call ortho2eva30(mmax,z,pols,dersx,dersy,
     $     w(iw1),w(iw2),w(iw3),w(iw4),w(iw5),w(iw6))
      return
      end
c 
c 
c 
c 
c 
      subroutine ortho2eva30(mmax,z,pols,dersx,dersy,
     $     f1,f2,f3,f4,f5,f6)
c 
c     evaluate the orthonormal polynomials on the triangle,
c     together with their derivatives with respect to x and y.
c 
      implicit real *8 (a-h,o-z)
        save
      dimension z(2), pols(1), dersx(1), dersy(1)
      dimension f1(1), f2(1)
      dimension f3(1), f4(1), f5(1), f6(1)
      data init/1/
c 
      if( init.eq.1 ) then
         zero=0
         sqrt2=dsqrt(2.0d0)
         sqrt3=dsqrt(3.0d0)
         init=0
         r11=-1.0d0/3.0d0
         r12=-1.0d0/sqrt3
         r21=1.0d0/9.0d0*(-sqrt3)*sqrt3
         r22=1.0d0/9.0d0*(6.0d0)*sqrt3
      endif
c 
      a=z(1)
      b=z(2)
c 
c     ... map the standard triangle to the right
c     triangle with the vertices (-1,-1), (1,-1), (-1,1)
c 
      x=r11+r12*b+a
      y=r21+r22*b
c 
c     ... evaluate the Koornwinder's polynomials
c     the via three terms recursion
c 
      par1=(2*x+1+y)/2
      par2=(1-y)/2
      call klegeypols(par1,par2,mmax,f1)
c 
      do 1200 m=0,mmax
         par1=2*m+1
ccc         call kjacopols(y,par1,zero,mmax,f2(1+m*(mmax+1)))
         call kjacopols(y,par1,zero,(mmax-m),f2(1+m*(mmax+1)))
 1200 continue
c 
      u=(2*x+1+y)/(1-y)
      v=y
      call klegepols2(u,mmax,f3,f4)
      do 1400 m=0,mmax
         par1=2*m+1
         call kjacopols2(v,par1,zero,mmax-m,
     $        f5(1+m*(mmax+1)),f6(1+m*(mmax+1)))
 1400 continue
c 
      kk=0
      do 2200 m=0,mmax
c 
         do 2000 n=0,m
         kk=kk+1
c 
c     ... evaluate the polynomial (m-n, n)
c 
         pols(kk)=f1(m-n+1)*f2(n+1+(m-n)*(mmax+1))
c 
c     ... and normalize it
c 
         scale=dsqrt((1.0d0+(m-n)+n)*(1.0d0+(m-n)+(m-n))/sqrt3)
         pols(kk)=pols(kk)*scale
c 
c     .. evaluate the derivatives
c 
         s11=2/(1-y)
         s12=2*(1+x)/(1-y)**2
c 
         pol11=f3(1+(m-n))
         der11=f4(1+(m-n))
  
c 
         pol12=f5(n+1+(m-n)*(mmax+1))
         der12=f6(n+1+(m-n)*(mmax+1))
c 
         pol13=((1-y)/2)**(m-n)
         der13=-((m-n)/(1-y))*((1-y)/2)**(m-n)
c 
         dersx(kk)=
     $        (der11*s11)*pol12*pol13
         dersy(kk)=
     $        (der11*s12)*pol12*pol13+pol11*(der12*pol13+der13*pol12)
c 
         t1=dersx(kk)
         t2=dersy(kk)
         dersx(kk)=t1
         dersy(kk)=t1*r12+t2*r22
         dersx(kk)=dersx(kk)*scale
         dersy(kk)=dersy(kk)*scale
c 
 2000	 continue
 2200 continue
      return
      end
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        this is the end of the debugging code and the beginning
c        of the jacobi expansion routines
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        subroutine kjacopols(x,a,b,n,pols)
        implicit real *8 (a-h,o-z)
        save
        dimension pols(1)
c 
c       evaluates a bunch of Jacobi polynomials
c       at the user-provided point x
c 
c       ... if n=0 or n=1 - exit
c 
        if(n .ge. 2) goto 1200
        pols(1)=1
        if(n .eq. 0) return
c 
        pols(2)=(a/2-b/2)+(1+a/2+b/2)*x
        return
 1200 continue
c 
        pols(1)=1
        pols(2)=(a/2-b/2)+(1+a/2+b/2)*x
c 
c       n is greater than 2. conduct recursion
c 
        pkm1=1
        pk=(a/2-b/2)+(1+a/2+b/2)*x
        pk=1
        pkp1=(a/2-b/2)+(1+a/2+b/2)*x
c 
        do 2000 k=2,n
c 
        pkm1=pk
        pk=pkp1
        alpha=(2*k+a+b-1)*(a**2-b**2+(2*k+a+b-2)*(2*k+a+b)*x)
        beta=2*(k+a-1)*(k+b-1)*(2*k+a+b)
        pkp1=(alpha*pk-beta*pkm1)/(2*k*(k+a+b)*(2*k+a+b-2))
        pols(k+1)=pkp1
 2000 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine kjacopol(x,a,b,n,pol,der)
        implicit real *8 (a-h,o-z)
c 
c       evaluates a single Jacobi polynomial (together
c       with its derivative) at the user-provided point x
c 
c       ...
c 
        save
        pkm1=1
        dkm1=0
        pk=(a/2-b/2)+(1+a/2+b/2)*x
        dk=(1+a/2+b/2)
c 
        pk=1
        dk=0
        pkp1=(a/2-b/2)+(1+a/2+b/2)*x
        dkp1=(1+a/2+b/2)
c 
c 
c        if n=0 or n=1 - exit
c 
        if(n .ge. 2) goto 1200
        pol=1
        der=0
        if(n .eq. 0) return
c 
        pol=(a/2-b/2)+(1+a/2+b/2)*x
        der=(1+a/2+b/2)
        return
 1200 continue
c 
        pol=1
        der=0
        pol=(a/2-b/2)+(1+a/2+b/2)*x
        der=(1+a/2+b/2)
c 
c       n is greater than 2. conduct recursion
c 
        do 2000 k=2,n
c 
        pkm1=pk
        pk=pkp1
        dkm1=dk
        dk=dkp1
        alpha1=(2*k+a+b-1)*(a**2-b**2)
        alpha2=(2*k+a+b-1)*((2*k+a+b-2)*(2*k+a+b))
        beta=2*(k+a-1)*(k+b-1)*(2*k+a+b)
        gamma=(2*k*(k+a+b)*(2*k+a+b-2))
        pkp1=((alpha1+alpha2*x)*pk-beta*pkm1)/gamma
        dkp1=((alpha1+alpha2*x)*dk-beta*dkm1+alpha2*pk)/gamma
        pol=pkp1
        der=dkp1
 2000 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine kjacopols2(x,a,b,n,pols,ders)
        implicit real *8 (a-h,o-z)
        save
        dimension pols(1), ders(1)
c 
c       evaluates a bunch of Jacobi polynomials (together
c       with their derivatives) at the user-provided point x
c 
c       ... if n=0 or n=1 - exit
c 
        if(n .ge. 2) goto 1200
        pols(1)=1
        ders(1)=0
        if(n .eq. 0) return
c 
        pols(2)=(a/2-b/2)+(1+a/2+b/2)*x
        ders(2)=(1+a/2+b/2)
        return
 1200 continue
c 
        pols(1)=1
        ders(1)=0
        pols(2)=(a/2-b/2)+(1+a/2+b/2)*x
        ders(2)=(1+a/2+b/2)
c 
c       n is greater than 2. conduct recursion
c 
        pkm1=1
        dkm1=0
        pk=(a/2-b/2)+(1+a/2+b/2)*x
        dk=(1+a/2+b/2)
c 
        pk=1
        dk=0
        pkp1=(a/2-b/2)+(1+a/2+b/2)*x
        dkp1=(1+a/2+b/2)
c 
        do 2000 k=2,n
c 
        pkm1=pk
        pk=pkp1
        dkm1=dk
        dk=dkp1
        alpha1=(2*k+a+b-1)*(a**2-b**2)
        alpha2=(2*k+a+b-1)*((2*k+a+b-2)*(2*k+a+b))
        beta=2*(k+a-1)*(k+b-1)*(2*k+a+b)
        gamma=(2*k*(k+a+b)*(2*k+a+b-2))
        pkp1=((alpha1+alpha2*x)*pk-beta*pkm1)/gamma
        dkp1=((alpha1+alpha2*x)*dk-beta*dkm1+alpha2*pk)/gamma
        pols(k+1)=pkp1
        ders(k+1)=dkp1
 2000 continue
c 
        return
        end
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        this is the end of the debugging code and the beginning
c        of the evaluation of scaled Legendre polynomials
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine klegeypols(x,y,n,pols)
        implicit real *8 (a-h,o-z)
        save
        dimension pols(1)
c 
c       evaluate a sequence of scaled Legendre polynomials
c       P_n(x/y) y^n, with the parameter y \in [0..1]
c 
c       ...
c 
        pkm1=1
        pk=x
c 
        pk=1
        pkp1=x
c 
c 
c        if n=0 or n=1 - exit
c 
        if(n .ge. 2) goto 1200
        pols(1)=1
        if(n .eq. 0) return
c 
        pols(2)=x
        return
 1200 continue
c 
        pols(1)=1
        pols(2)=x
c 
c       n is greater than 2. conduct recursion
c 
        do 2000 k=1,n-1
c 
        pkm1=pk
        pk=pkp1
        pkp1=( (2*k+1)*x*pk-k*pkm1*y**2 )/(k+1)
        pols(k+2)=pkp1
 2000 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine klegepols2(x,n,pols,ders)
        implicit real *8 (a-h,o-z)
        save
        dimension pols(1),ders(1)
c 
c       evaluates a bunch of Legendre polynomials (together
c       with their derivatives) at the user-provided point x
c 
c       ... if n=0 or n=1 - exit
c 
        if(n .ge. 2) goto 1200
        pols(1)=1
        ders(1)=0
        if(n .eq. 0) return
c 
        pols(2)=x
        ders(2)=1
        return
 1200 continue
c 
c       n is greater than 1. conduct recursion
c 
        pols(1)=1
        ders(1)=0
        pols(2)=x
        ders(2)=1
c 
        pkm1=1
        pk=x
        dkm1=0
        dk=1
c 
        pk=1
        pkp1=x
        dk=0
        dkp1=1
c 
        do 2000 k=1,n-1
c 
        pkm1=pk
        pk=pkp1
        dkm1=dk
        dk=dkp1
c 
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
        dkp1=( (2*k+1)*(x*dk+pk)-k*dkm1 )/(k+1)
        pols(k+2)=pkp1
        ders(k+2)=dkp1
 2000 continue
c 
        return
        end
  
  
  
  
  
  
  
  
  
  
  
  
  
