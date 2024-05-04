program test_gaussq
  IMPLICIT REAL *8 (A-H,O-Z)
  integer :: kind, kpts
  double precision ::  t(2000), w(2000), b(2000), endpts(2), f(2000)
  double precision :: done, rints(10000), exacts(10000), diffs(10000)
  double precision :: pars(10)

  external :: fjacobi

  call prini(6,13)

  !C 
  !C       CREATE THE INPUT PARAMETERS
  !C 

  !
  ! set the kind of quadrature to generate
  ! 1: legendre
  ! 2: chebyshev of first kind
  ! 3: chebyshev of second kind
  ! 4: hermite
  ! 5: jacobi with alpha,beta > -1
  ! 6: generalized laguerre, with alpha > -1
  !

  kind = 5

  if (kind .eq. 1) then
     print *, '- - - generating gauss-legendre quadrature - - -'
     print *
  else if (kind .eq. 2) then
     print *, '- - - generating gauss-cheby 1st quadrature - - -'
     print *     
  else if (kind .eq. 5) then
     print *, '- - - generating gauss-jacobi quadrature - - -'
     print *     
  else
     print *, ' kind not set to 1 -> 6'
     stop
  end if
  
  
  
  PRINT *, 'ENTER n'
  READ *,n
  CALL PRINF('n=*',n,1 )
  
  !PRINT *, 'ENTER k'
  !READ *,k
  !CALL PRINF('k=*',k,1 )
 
  alpha = 1
  beta = 1
  kpts = 0

  
  ! construct the gaussian quadrature
 
  call prinf('before gaussq, n=*',n,1)
  call gaussq(kind, n, alpha, beta, kpts, endpts, b, t, w)
  call prin2('the nodes, t =*', t, n)
  call prin2('the weights, w =*', w, n)



  done = 1
  
  if (kind .eq. 1) then
     
     ! test gauss-legendre
     do k = 0,2*n-1
        rints(k+1) = 0
        do i = 1,n
           rints(k+1) = rints(k+1) + w(i) * (t(i)**k)
        end do
        exacts(k+1) = done/(k+1) * ( done**(k+1) - (-done)**(k+1) )
        diffs(k+1) = rints(k+1) - exacts(k+1)
     end do

  else if (kind .eq. 2) then
     ! test gauss-chebyshev of first kind
     print *, 'testing code for chebyshev not yet implemented'
     stop
     
  else if (kind .eq. 3) then
     ! test gauss-chebyshev of second kind
     print *, 'testing code for chebyshev 2nd kind not yet implemented'
     stop

  else if (kind .eq. 4) then
     ! test gauss-hermite
     print *, 'testing code for hermite not yet implemented'
     stop


  else if (kind .eq. 5) then
     ! test gauss-jacobi, comparing against adaptive integration

     pars(1) = alpha
     pars(2) = beta

     do k = 0,2*n-1
        rints(k+1) = 0
        do i = 1,n
           rints(k+1) = rints(k+1) + w(i) * (t(i)**k)
        end do

        ! calculate the exact integral
        a = -1
        b = 1
        m = 16
        eps = 1.0d-13

        p = k
        call  adapgaus(ier, a, b, fjacobi, p, pars, m, eps, &
             exacts(k+1), maxrec, numint)

        diffs(k+1) = rints(k+1) - exacts(k+1)
     end do
     


     
  else if (kind .eq. 6) then
     ! test gauss-laguerre
     print *, 'testing code for laguerre not yet implemented'
     stop

     
  end if


  call prin2('exact integrals = *', exacts, 2*n)
  call prin2('rints integrals = *', rints, 2*n)
  call prin2('and diffs = *', diffs, 2*n)

  errmax = -1
  do i = 1,2*n
     if (abs(diffs(i)) .gt. errmax) errmax = abs(diffs(i))
  end do
  
  call prin2('and maximum abs quad error is errmax = *', errmax, 1)
  
end program test_gaussq






function fjacobi(x, p, pars)
  implicit double precision (a-h,o-z)
  double precision :: x, pars(*)

  alpha = pars(1)
  beta = pars(2)
  fjacobi = x**p * (1-x)**alpha * (1+x)**beta

  return
end function fjacobi
