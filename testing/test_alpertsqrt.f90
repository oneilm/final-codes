program test_alpertsqrt

  use iso_fortran_env

  implicit real( real64 ) (a-h,o-z)
  integer( int32 ) :: kind, kpts
  real( real64 ) ::  t(2000), w(2000), b(2000), endpts(2), f(2000)
  real( real64 ) :: done, rints(10000), exacts(10000), diffs(10000)
  real( real64 ) :: pars(10)

  real( real64 ) :: xnodes(10000), whts(10000)
  
  real( real128 ) :: aq, bq, qint, epsq, pq, qpars(10), alphaq, betaq
  real( real128 ) :: hq, tq, fq, q1, q2
  
  real( real128 ), external :: fsmooth, fleft, fright

  call prini(6,13)



  print *, '- - - - - testing smooth alpert quadrature - - - - -'

  aq = -1
  bq = 1
  a = aq
  b = bq
  norder = 16
  n = 2*norder+1+5
  call prinf('norder = *', norder, 1)
  call prin2('a = *', a, 1)
  call prin2('b = *', b, 1)
  call prinf('number of nodes, n = *', n, 1)
  call alpert_smooth(a, b, n, norder, xnodes, whts)

  call prin2('after alpert_smooth, xnodes = *', xnodes, n)
  call prin2('after alpert_smooth, whts = *', whts, n)

  dint = 0
  do i = 1,n
     tq = xnodes(i)
     fq = fsmooth(tq, q1, q2)
     dint = dint + whts(i)*fq
  end do

  call prin2('from alpert, integral = *', dint, 1)

  ! compute the integral using adaptive integration and compare
  ier = 0
  qint = 0
  m = 16
  epsq = 1.0q-15
  call  adapgaus_quad(ier, aq, bq, fsmooth, q1, q2, m, epsq, &
       qint, maxrec, numint)

  call prinq('after adapgaus_quad, integral = *', qint, 1)

  err = (qint-dint)/qint
  call prin2('and relative error = *', err, 1)





  print *
  print *
  print *
  print *, '- - - - - testing left sqrt alpert quadrature - - - - -'

  aq = .23q0
  bq = 4.3q0
  a = aq
  b = bq
  norder = 16
  n = 2*norder+1+5
  call prinf('norder = *', norder, 1)
  call prin2('a = *', a, 1)
  call prin2('b = *', b, 1)
  call prinf('number of nodes, n = *', n, 1)
  call alpert_sqrta(a, b, n, norder, xnodes, whts)

  call prin2('after alpert_smooth, xnodes = *', xnodes, n)
  call prin2('after alpert_smooth, whts = *', whts, n)

  dint = 0
  do i = 1,n
     tq = xnodes(i)
     fq = fleft(tq, aq, q2)
     dint = dint + whts(i)*fq
  end do

  call prin2('from alpert, integral = *', dint, 1)

  ! compute the integral using adaptive integration and compare
  ier = 0
  qint = 0
  m = 16
  epsq = 1.0q-15
  call  adapgaus_quad(ier, aq, bq, fleft, aq, q2, m, epsq, &
       qint, maxrec, numint)

  call prinq('after adapgaus_quad, integral = *', qint, 1)

  err = (qint-dint)/qint
  call prin2('and relative error = *', err, 1)


  
  
  print *
  print *
  print *
  print *, '- - - - - testing right sqrt alpert quadrature - - - - -'

  aq = 0
  bq = 2.3q0
  a = aq
  b = bq
  norder = 16
  n = 2*norder+1+5
  call prinf('norder = *', norder, 1)
  call prin2('a = *', a, 1)
  call prin2('b = *', b, 1)
  call prinf('number of nodes, n = *', n, 1)
  call alpert_sqrtb(a, b, n, norder, xnodes, whts)

  call prin2('after alpert_smooth, xnodes = *', xnodes, n)
  call prin2('after alpert_smooth, whts = *', whts, n)

  dint = 0
  do i = 1,n
     tq = xnodes(i)
     fq = fright(tq, bq, q2)
     dint = dint + whts(i)*fq
  end do

  call prin2('from alpert, integral = *', dint, 1)

  ! compute the integral using adaptive integration and compare
  ier = 0
  qint = 0
  m = 16
  epsq = 1.0q-15
  call  adapgaus_quad(ier, aq, bq, fright, bq, q2, m, epsq, &
       qint, maxrec, numint)

  call prinq('after adapgaus_quad, integral = *', qint, 1)

  err = (qint-dint)/qint
  call prin2('and relative error = *', err, 1)


  
  
  
end program test_alpertsqrt





function fsmooth(t, p1, p2)
  use iso_fortran_env
  implicit real( real128 ) (a-h,o-z)
  real( real128 ) :: t, p1, p2, fsmooth

  fsmooth = cos(10*t)
  
  return
end function fsmooth





function fleft(t, a, p2)
  use iso_fortran_env
  implicit real( real128 ) (a-h,o-z)
  real( real128 ) :: t, a, p2, fleft

  fleft = cos(10*t)/sqrt(t-a)
  
  return
end function fleft





function fright(t, b, p2)
  use iso_fortran_env
  implicit real( real128 ) (a-h,o-z)
  real( real128 ) :: t, b, p2, fright

  fright = cos(10*t)/sqrt(b-t)
  
  return
end function fright






