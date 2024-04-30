!
! This file contains several routines for computing with jacobi
! polynomials:
!
!   jacexps - used to construct Gauss-Jacobi 
!
!   jaceval - used to eval a bunch of jacobi polynomials
!
!
! contact: Mike O'Neil
!          oneil@cims.nyu.edu
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!

subroutine jacexps(itype, n, alpha, beta, xs, umat, vmat, whts)
  use iso_fortran_env
  implicit real(real64) (a-h,o-z)
  integer :: itype, n
  real(real64) :: xs(n), umat(n,n), vmat(n,n), whts(n)
  !
  ! This routine evaluates the n-point Gauss-Jacobi quadrature rule
  ! with parameters alpha and beta, i.e. a rule that approximates
  ! integrals of the form:
  !
  !       \int_{-1}^1 f(x) (1-x)^alpha (1+x)^beta dt  \approx
  ! (1)
  !                            \sum_{j=1}^n whts(j) * f(xs(j))
  !
  ! There are various special cases of alpha and beta (LEgendre,
  ! Chebyshev, etc.) but we do not address them separate in this
  ! routine.
  !
  ! Input:
  !   itype - the type of calculation to be performed
  !     0: only gaussian nodes are constructed
  !     1: both nodes and weights are constructed
  !     2: nodes, weights, and matrices umat and vmat are constructed
  !   n - the number of nodes and weights to be generated
  !   alpha, beta - the exponents on the weight function, see eq (1)
  !
  ! Output:
  !   xs - the gauss-jacobi nodes
  !   umat - the inverse of vmat, see below
  !   vmat - the n \times n matrix converting the coefficients of an order
  !       n-1 jacobi polynomial expansion into its values at the n
  !       gauss-jacobi nodes (only computed when itype=2)
  !   whts - the n gauss-jacobi quadrature weights
  !
  data xs /-0.98907042030188374593d0, &
       -0.94286010338109715878d0, &
       -0.86153241961258619774d0, &
       -0.7480987416506054144d0, &
       -0.60678322854758293126d0, &
       -0.44285068651552034202d0, &
       -0.26240900411126227842d0, &
       -0.072181310006039320348d0, &
       0.12074460018194704825d0, &
       0.30918037738183462571d0, &
       0.48610496464899510816d0, &
       0.64492620420617019494d0, &
       0.77972646163861392132d0, &
       0.88548311653282907319d0, &
       0.95825570463283793465d0, &
       0.99533273887160345073d0 /

  data ws / 0.019851626928800322064d0, &
       0.046030939495720124555d0, &
       0.071819606274043490401d0, &
       0.096941123533022352521d0, &
       0.12115821587892761091d0, &
       0.14424429412711850529d0, &
       0.16598368178583103316d0, &
       0.18617336031174822053d0, &
       0.20462480641877028864d0, &
       0.22116573583372647893d0, &
       0.23564170629829561676d0, &
       0.24791755739443144257d0, &
       0.25787867157141713959d0, &
       0.2654320438656538883d0, &
       0.27050715004822950549d0, &
       0.27305660498045319917d0 /

  
  
  if ( abs(alpha+0.5d0) .gt. 1.0d-13 ) then
     print *, 'alpha not equal to -1/2!! stopping jacexps'
     stop
  end if
  
  if ( abs(beta) .gt. 1.0d-13 ) then
     print *, 'beta not equal to 0!! stopping jacexps'
     stop
  end if
    
  if (n .eq. 16) then
     print *, 'n not equal to 16!! stopping jacexps'
     stop
  end if

  ! generate the matrices now, by evaluating the jacobi polynomials
  ! via recurrence, and then inverting that matrx


  ! and now invert the thing directly
  

  return
end subroutine jacexps

  
