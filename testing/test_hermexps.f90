program test_hermexps
  !
  ! DEPENDENCIES: prini.f
  !               hermexps.f
  !               pplot.f
  !               gammanew_eval.f
  !

  implicit real *8 (a-h,o-z)
  dimension :: funs(1000),ders(1000),diffs(1000), &
       funsp(1000),funsm(1000),dersp(1000),dersm(1000)

  dimension :: w(10000),rints(1000),rintsp(1000), &
       rintsm(1000),funs2(1000)

  double precision :: ts(10000), vals(10000), vals2(10000)

  dimension :: xs(1000),uu(100000),vv(100000),whts(1000), &
           ww(100000),coefs(1000), dints(10000)



  
  call prini(6,13)
  !
  ! SET ALL PARAMETERS
  !
  PRINT *, 'ENTER n'
  READ *,n
  CALL PRINf('n=*',n,1 )

  nmax=n

  h = 0.00001d0
  x=2

  call prin2('test point, x = *', x, 1)
  
  ifinit=1
  call herminte(x, funs, ders, rints,n,w,ifinit)

  call prin2('funs=*',funs,nmax)
  call prin2('ders=*',ders,nmax)
  call prin2('rints=*',rints,2)

  call herminte(x+h,funsp,dersp,rintsp,n,w,ifinit)
  call herminte(x-h,funsm,dersm,rintsm,n,w,ifinit)

  funs2(1) = (rintsp(1)-rintsm(1))/h/2

  !cccc        call prin2('and rintsp=*',rintsp,1)
  !cccc        call prin2('and rintsm=*',rintsm,1)

  funs2(2)=(rintsp(2)-rintsm(2))/h/2

  call prin2('and funs2=*',funs2,2)
  call prin2('and funs2(1)/funs(1)=*',funs2(1)/funs(1),1)
  call prin2('and funs2(2)/funs(2)=*',funs2(2)/funs(2),1)

  do  i=1,n
     funs2(i)=(rintsp(i)-rintsm(i))/h/2
     diffs(i)=funs2(i)-funs(i)
  end do


  call prin2('and funs2=*',funs2,n)
  call prin2('while funs=*',funs,n)
  call prin2('and diffs=*',diffs,n)

  print *
  print *
  print *
  
  !c
  !c       check that the matrices u,v are inverses of each other
  !c
  itype=2
  call hermexps(itype,n,xs,uu,vv,whts,w)


  call prin2('uu=*',uu,n*n)
  !cccc        call matmat(uu,vv,n,ww)

  call matmat(vv,uu,n,ww)


  !call prin2('uu vv =*',ww,n*n)
  ddd = 0
  do i=1,n
     ww(i+n*(i-1)) = ww(i+n*(i-1)) - 1
  end do

  do i = 1,n**2
     ddd = ddd + ww(i)**2
  end do
  ddd = sqrt(ddd)
  call prin2('frobenius norm of I - UV = *', ddd, 1)
  
  
  !c
  !c       construct a test function and decompose it
  !c
  do  i=1,n
     funs(i)=exp(-xs(i)**2/2) * xs(i)
     !cccc        funs(i)=1
  end do


  call prin2('funs as created*',funs,n)

  call matvec(uu,funs,coefs,n)
  call prin2('and coefs=*',coefs,n)


  !
  ! test (debug) the quadratures
  !
  call prin2('quadrature weights = *', whts, n)
  call prin2('quadrature nodes = *', xs, n)
  print *
  do i =1,n
     print *, 'node ', i, '  ',xs(i)
  end do
  
  print *

  dint = 0
  p = 4.0d0
  do i = 1,n
     dint = dint + whts(i) * xs(i)**p * exp(-xs(i)**2)
     !dint = dint + whts(i) * cos(xs(i)) * exp(-xs(i)**2)
  enddo

  xxx = (1+p)/2
  call gammanew_eval(xxx, gam)
  dexact = 0.5d0 * (1.0d0 + (-1.0d0)**p) *gam 

  !pi = 4*atan(1.0d0)
  !dexact = sqrt(pi)/exp(0.25d0)
  


  call prin2('integrating x^p, p = *', p, 1)
  call prin2('from quad, integral = *', dint, 1)
  call prin2('exact integral = *', dexact, 1)
  call prin2('error = *', dint - dexact, 1)
  call prin2('ratio, dint/dexact = *', dint/dexact, 1)
  print *
  print *

  dint = 0
  do i = 1,n
     dint = dint + whts(i) * cos(xs(i)) * exp(-xs(i)**2)
  enddo

  pi = 4*atan(1.0d0)
  dexact = sqrt(pi)/exp(0.25d0)
  
  call prin2('integrating cos(x)*', p, 0)
  call prin2('from quad, integral = *', dint, 1)
  call prin2('exact integral = *', dexact, 1)
  call prin2('error = *', dint - dexact, 1)
  call prin2('ratio, dint/dexact = *', dint/dexact, 1)
  print *
  print *

  
  !
  ! plot a specific hermite function
  !
  m = n
  a = xs(1)
  b = xs(n)
  nplot = 400

  ifinit = 1
  do i = 1,nplot
     ts(i) = a + (b-a)/(nplot-1)*(i-1)
     call hermfuns(ts(i) , m+10, funs, ders, ifinit, w)
     vals(i) = funs(m+1)
  end do

  iw = 21
  itype = 1
  call pyplot(iw, ts, vals, nplot, itype, 'hermite polynomial*')



  !
  ! now try out the scaled routines
  !
  print *, '----------- testing the scaled function routines ----------'
  sc = 2.5d0
  print *
  print *
  call prin2('scale factor, sc = *', sc, 1)

  ! plot the scaled function on the same axes
  ifinit = 1
  do i = 1,nplot
     call hermfuns_scaled(sc, ts(i) , m+10, funs, ders, ifinit, w)
     vals2(i) = funs(m+1)
  end do

  !call prin2('vals2 = *', vals2, nplot)
  
  iw = 22
  itype = 1
  call pyplot2(iw, ts, vals, nplot, itype, &
       ts, vals2, nplot, itype, &
       'hermite polynomial and its scaled version*')

  
  !
  ! create a big grid, apply trapezoidal rule and make sure they are
  ! unit functions
  !
  a = -30
  b = 30
  nquad = 500
  h = (b-a)/(nquad-1)

  do j = 1,n
     dints(j) = 0
     do i = 1,nquad
        t = a + h*(i-1)
        call hermfuns_scaled(sc, t, 2*n+20, funs, ders, ifinit, w)
        dints(j) = dints(j) + h*funs(j)**2
     end do
  end do

  call prin2('norm of scaled h_0 through h_{n-1} = *', dints, n)

  do i = 1,n
     dints(i) = dints(i) - 1
  end do
  call prin2('and errors in norms = *', dints, n)
  

  print *
  print *
  print *
  
  !
  ! now construct the quadrature, and integrate to see if we get the
  ! same answer using the quadrature weights
  !
  itype = 2
  call hermexps_scaled(sc, itype, n, xs, uu, vv, whts, w)

  dint = 0
  do i = 1,n
     dint = dint + whts(i)*xs(i)**p * exp(-xs(i)**2 / sc**2)
  end do

  xxx = (1+p)/2
  call gammanew_eval(xxx, gam)
  dexact = 0.5d0 * (1.0d0 + (-1.0d0)**p) * (sc**(1+p))*gam 

  
  call prin2('integrating x^p, p = *', p, 1)
  call prin2('from quad, integral = *', dint, 1)
  call prin2('exact integral = *', dexact, 1)
  call prin2('error = *', dint - dexact, 1)
  call prin2('ratio, dint/dexact = *', dint/dexact, 1)
  print *
  print *


  !
  ! check the values to coefficients matrix
  !
  do i = 1,n
     vals(i) = xs(i) * exp(-xs(i)**2 / 2/sc**2) &
          / (5*sqrt(5.0d0)*pi**(0.25d0)/4)
  end do

  call matvec(uu, vals, coefs, n)
  call prin2('coefs in scaled expansion = *', coefs, n)
  
  

  ! check that uu * vv = I
  call matmat(vv, uu, n, ww)
  !call prin2('uu vv =*',ww,n*n)

  ddd = 0
  do i=1,n
     ww(i+n*(i-1)) = ww(i+n*(i-1)) - 1
  end do

  do i = 1,n**2
     ddd = ddd + ww(i)**2
  end do
  ddd = sqrt(ddd)
  call prin2('frobenius norm of scaled I - UV = *', ddd, 1)
   

  stop
end program test_hermexps






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   these are the local testing subroutines
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine matvec(a,x,y,n)
  implicit real *8 (a-h,o-z)
  dimension a(n,n)
  real *8 x(n),y(n),d

  do  i=1,n
     d=0
     do  j=1,n
        d=d+a(i,j)*x(j)
     end do

     y(i)=d
  end do

  return
end subroutine matvec





subroutine matmat(a,b,n,c)
  implicit real *8 (a-h,o-z)
  save
  dimension a(n,n),b(n,n),c(n,n)

  do  i=1,n
     do  j=1,n
        d=0
        do  k=1,n
           d=d+a(i,k)*b(k,j)
        end do

        c(i,j)=d
     end do
  end do

  return
end subroutine matmat
