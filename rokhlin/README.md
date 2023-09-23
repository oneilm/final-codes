#Vladimir Rokhlin code library circa 2007


1. The muller polynomial root-finder is in the file test13.f. An exact copy of it can be found in the file muller.f. The file test13.f can be run by means of the command file test13.

2. A Gram-Schmidt inverter for complex *16 matrices is in corthom.f, together with the debugging code. It can be run by means of the command file corthom.

3. A Gram-Schmidt inverter for real *8 matrices is in orthom.f, together with the debugging code. It can be run by means of the command file orthom.

4. A subroutine constructing the Chebychev quadratures is in the 
file chebwhts.f, together with the debugging code. 

5. A subroutine constructing the matrices converting 
values at the Chebychev nodes of a Chebychev expansion 
into its coefficients (and back) is in the file chebexps.f, 
together with the debugging code. This subroutine
also produces the weights of the Chebychev quadrature;
the other subroutine in this file calculates the value
and the derivative of a chebychev expansion at a point.

6. Page's adaptive gaussian quadrature routine is in the file
gausspk.f, together with the debugging (testing) code.
It can be run by calling the routine gausspk

7. A pair of subroutines alptrap, alptrap2 are in the file test102.f,
together with their debugging code. These subroutines implement 
Alpert's end-point corrected trapezoidal quadrature rules for 
smooth functions on the interval. They implement rules of orders
2,4,6,8,10,12.

8. A subroutine funstep is in the file funstep.f. It constructs a smooth 
(c^{\infty}) elementary function f such that 

		c^{\infty} on the real line. 
		f(x)=1 for all x .le. 0
		f(x)=0 for all x .ge. \pi/2

9. My adaptive Chebychev quadrature code is in 
the file adapcheb.f, together with the debugging code.

10. My adaptive Gaussian quadrature code  for real functions is 
in the file adapgaus.f, together with the debugging code.

11. My adaptive Gaussian quadrature code  for complex functions is 
in the file cadapgau.f, together with the debugging code.

12. A complex conjugate gradient routine is in the file ccongr.f,
together with the debugging code.

13. A real conjugate gradient algorithm is in the file congra.f,
together with the debugging routine.

14. A complex conjugate residual algorithm is in the file cmconr.f,
together with the debugging routine.

15. A subroutine implementing the end-point corrected trapezoidal rule 
for the evaluation of the time integration for the potential 
theory for the wave equation in two dimensions is in the file 
wavgrcor.f, together with the debugging code.

16. A subroutine for the rapid evaluation of the Bessel functions 
j0, j1, y0, y1 of the real argument is in the file jy01rea.f, 
together with the debugging routines.

17. A slow but accurate function for the evaluation of the 
gamma function is in the file gammas4.f, together with the 
debugging routines.

18. The code constructing the 6 classical gaussian quadratures is in
the file gaussq.f (obtained from Golub years and years ago). 
This version uses a portable gamma function routine, unlike the 
original.
> 18a - The code constructing the 6 classical gaussian
> 		quadratures in extended precision is in the file gauq16.f
> 		(obtained from 
> 		Golub and modified years and years ago). This version uses a
> 		portable gamma function routine, unlike the original.

19. The code constructing the gaussian weights and nodes on the 
interval [-1,1] is in the file gauswhts.f, together with 
the debugging code. The code is less general than that in 
gaussq.f, but it is very short, transparent, and portable.

20. The subroutine resampling (in an equispaced manner) an analytically
specified curve is in the file anaresa.f, together with the
debugging code.

21. The subroutine resampling (in an equispaced manner) a user-specified
curve (the user specifies it as a collection of points in the plane)
is in the file rsresa.f, together with the debugging code.

22. The subroutine sortanyn sorting an array of integers is in the 
file sortanyn.f, together with the debugging code.

23. The subroutine d2mallb constructing the quad-tree structure
for the FMM in two dimensions is in the file d2mallb.f, together
with the debugging code.

24. The subroutine d2mudv uses the Jacobi iterations to construct the 
singular value decomposition for a dense REAL *8  matrix. It is 
in the file d2mudv.f, together with the debugging code. 

25. The subroutine d2mcudv uses the Jacobi iterations to construct the 
singular value decomposition for a dense COMPLEX *16 matrix. It is 
in the file d2mudv.f, together with the debugging code. 

26. The subroutine d2mcomp for the recursive compression of a matrix can
be found in the file d2mcomp.f, together with the debugging codes.
It avoids the construction of the whole matrix by doing a recursive
QR decomposition of a matrix of large dimensionality and low rank.

27. The subroutine d2mstrcr for the construction of the complete 
adaptive multipole logical structure can be found in the file
d2mstrcr.f, together with the debugging code. 

28. The subroutine (written by Grishka) for the accurate and reasonably
rapid evaluation of the Hankel functions H_0, H_1 of the complex 
argument can be found in the file hank71.f. This file contains no
debugging code.

29. The subroutine for the chebychev approximation of functions tabulated
at non-chebychev nodes can be found in the file chebnon.f, together
with the debugging code.

30. The subroutine svdpivot calculating the singular value decomposition
of a real *8 matrix can be found in the file svdpivot.f, together with
the debugging code. It uses a combination of the pivoted Gram-Schmidt
and the Jacobi rotations, and is supposed to be quite efficient when
the matrix is far from being full-rank. This subroutine is quite a
glutton for memory, and should only be used for problems that are not
too large. Otherwise, see D2MCOMP.f

30a.
The subroutine qsvdpiv calculating the singular value decomposition
of a real *16 matrix can be found in the file qsvdpiv.f, together with
the debugging code. It uses a combination of the pivoted Gram-Schmidt
and the Jacobi rotations, and is supposed to be quite efficient when
the matrix is far from being full-rank. This is a real *16 version of
the subroutine svdpivot (see above). 

31.
A REAL *16 subroutine constructing the matrices converting 
values at the Chebychev nodes of a Chebychev expansion 
into its coefficients (and back) is in the file qchfunde.f,
together with the debugging code. This subroutine
also produces the weights of the Chebychev quadrature;
the other subroutine in this file calculates the value
and the derivative of a chebychev expansion at a point.
Note that this is the REAL *16 version of the subroutine
chebexps (see above). Note also that there are differences 
between this routine and chebexps, since SUN has a bug in the 
calculation of the sine and cosine functions in REAL *16.

32.
A REAL *16 subroutine qerrfun evaluating the error function erf is 
in the file qerrfun.f, together with the debugging code. The other 
user-accessible subroutine in that file is qerrfunc, evaluating the
second error function erfc(x)=1-erf(x).

33.
A subroutine constructing the matrices converting values at 
the Gaussian nodes of a Gaussian expansion into its coefficients 
(and back) is in the file legeexps.f, together with the debugging 
code. This subroutine also produces the nodes and weights of the 
Gaussian quadrature; the other subroutine in this file calculates
the value and the derivative of a gaussian expansion at a point.
Another pair of subroutines produces the Legendre expansions of the
indefinite integral and the derivative of the user-specified Legendre
expansion.

34.
The subroutines that pack FORTRAN-produced data into FORTRAN
DATA atatements are in the files qdataarr.f, qdataard.f. The 
subroutine qdataarr prodices data with the extended precision,
and qdataard produced data with the double precision. Both use
extended precision data for input.

35.
The extended precision subroutine qpsquare calculating the 
orthogonal polynomials on the standard square is in the file 
qpsquare.f, together with the debugging code. It calculates the 
first 100 orthogonal polynomials.

36.
The double precision subroutine dpsquare calculating the 
orthogonal polynomials on the standard square is in the file 
dpsquare.f, together with the debugging code. It calculates the 
first 100 orthogonal polynomials.

37.
The extended precision subroutine qpsquain calculating the 
orthogonal polynomials of 1/z on the standard square is in the file 
qpsquain.f, together with the debugging code. It calculates the 
first 100 orthogonal polynomials.

38.
The double precision subroutine dpsquain calculating the 
orthogonal polynomials of 1/z on the standard square is in the file 
dpsquain.f, together with the debugging code. It calculates the 
first 100 orthogonal polynomials.

39.
The subroutine csvdpiv calculating the singular value decomposition
of a complex matrix can be found in the file csvdpiv.f, together with
the debugging code. it uses a combination of the pivoted Gram-Schmidt
and the Jacobi rotations, and is supposed to be quite efficient when
the matrix is far from being full-rank. This subroutine is quite a
glutton for memory, and should only be used for problems that are not
too large.

40.
The subroutine hank101 calculating the Hankel functions H_0, H_1
of a complex argument is located in the file hank101.f, together
with the debugging code. hank101 evaluates H_0(z), H_1(z) for 
arbitrary z in the complex plane, unless the exponent overflow 
occurs. It gives relative precision of 0.5E-14 everywhere, and is
reasonably fast. However, the code has not been optimized, and 
a simple removal of double calls to functions, etc. should speed
it up by a factor of 2 or so.

41. The subroutine d2msvd in the file d2msvd.f constructs the 
compressing singular value decomposition of a real *8 matrix. This
subroutine should be used for large matrices whose rank is relatively
small. d2msvd avoids the construction of the whole big matrix, 
performing the whole process recursively. The file d2msvd.f also
contains debugging subroutines.

42. The subroutine legexrts in the file legexrts. it attempts to find
all roots of a user-specified function located on the interval [a,b]
the function is specified by its values tabulated at the gaussian
nodes on that same interval [a,b]. evaluating this function at these
nodes is the user's responsibility; however, the gaussian nodes on the
interval [-1,1] can be obtained by calling the subroutine legeexps
(see). The file legexrts.f also contains debugging code.

43. The subroutine clagreva in the file clagreva lagrange interpolates 
a function f tabulated at the nodes t(i) in the complex plane to obtain
the value of f at the point x. The file clagreva  also contains 
debugging code.

44. The subroutine lagrcoef in the file lagrcoef.f constructs 
interpolation coiefficients connecting the values of a function
f tabulated at the nodes t(i) on the real line  with the value of f at
the real point x. The file lagrcoef.f does not contain ant debugging 
code, with all that follows fro this fact.

45. The subroutine gaubeam in the file gaubeam.f evaluates hankel 
functions H_0(z), H_1(z) scaled by an exponential. The purpose of
this subroutine is to prevent overflow in the calculations with
Gaussian beams. In reality, it uses a slightly changed (and renamed,
thank God) version of the subroutine hank101 (see).

46. The subroutine ifunscom in the file ifunscom.f evaluates all 
scaled modified Bessel functions e^{-x} \cdot I_m (x) of an argument x.
It uses the recursuion (up and down, and scaling), and is reasonably 
efficient, but not too efficeint.

46. The subroutine jbfun in the file jb.f evaluates all 
Bessel functions J_m (z) of a complex argument z.
It uses the recursuion (up and down, and scaling), and is reasonably 
efficient, but not too efficeint.

47. The subroutine fourtmat constructing explicitly the matrix
of DFT can be found in the file fourtmat.f. The same file contains
the subroutine fourtint constructing the interpolating matrix
from an equispaced grid to Gaussian one, and a few other (possibly)
useful subroutines.

48. The subroutine prolcrea evaluating the prolate spheroidal wave
functions and various associated eigenvalues is in the file prolcrea.f,
together with the debugging code. 

49. The subroutine prolall builds the machinery necessary for the 
extension of a band-limited function outside the interval on which
it has been specified by the user. It is in the file prolall.f,
together with the debugging code. Please note that this subroutine 
needs the subroutines found in the file proleva.f (see).

50. The subroutines leastsq1,leastsq2 are in the file leastsq.f,
together with the debugging code. They use a version of the Gram-Schmidt
process to solve a linear least squares problem (with real coefficients).

51. The file triaadap.f contains the subroutine triaadap that evaluates 
the integral of a function on a triangle via an adaptive nested tensor
product Gaussian rule. The file also contains the debugging code.

52. The files hjtran30.f, hjtran31.f contain a relatively final 
version of the code designing the beam version of the translation
operator in three dimensions. The code uses the numerically obtained 
beam, in the more or less final version. The files also contains the 
necessary debugging code. Note that the file hjtran31.f contains various 
service subroutines used by the file hjtran30.f, that contains all the 
ideologically sensitive material.

53. The file ch2dinit.f contains a code for the evaluation of 
two-dimensional Chebychev coefficients of a function tabulated at
(tensor product) Chebychev nodes on a square in two dimensions.
It also contains a code for the evaluation of a two-dimensional 
Chebychev expansion at the Chebychev (tensor product) nodes on a 
square in R^2, and for the evaluation of such an expansion at an
individual point in R^2.
    
54. The file wandquad.f contains the subroutine wandquad providing
wandzura's quadratures for functions on a triangle in the plane. 
The quadratures are nearly Gaussian. The file also contains the
necessary debugging code.

55. The subroutines cleastsq, cleasts2 are in the file cleastsq.f,
together with the debugging code. They use a version of the Gram-Schmidt
process to solve a linear least squares problem (with complex 
coefficients).

56.
The subroutine d3mstrcr for the construction of the complete 
adaptive multipole logical structure in three dimensaions can be 
found in the file d3mstrcr.f, together with the debugging code. 

57.
The subroutine ccauadap solving (adaptively) a system of non-stiff 
complex ODEs is in the file ccauadap.f, together with the debugging code.
The subroutine uses a spectral defferred correction scheme (of 
arbitrary user-specifeid order) driven by the explicit Euler method,
trapezoidally corrected to second order.  - 12.26.99

58.

The file stauadap.f contains two subroutines stauadap, stauexp, 
solving the Cauchy problem for a stiff system of complex ODEs 
via an arbitrary order defferred correction scheme (driven by the
implicit Euler iteration). The subroutine stauadap does not use 
the extrapolation to insure the proper behavior at infinity. Thus, 
it is limited to relatively low order schemes (for stiff problems; 
for non-stiff ones, the subroutine ccauadap should be used). The 
subroutine stauexp is the ultimate defferred correction subroutine 
for stiff problems. It uses two defferred correction schemes, and 
linearly combines ("extrapolates") them to obtain strong stiff 
stability at infinity. The file cauadap.f also contains debugging 
codes for both subroutines stauadap, stauexp.

59.

The file bellret.f contains the subroutine bellret evaluatin
at the point x \in [0,1] a smooth step function, i.e. function 
bell(x) such that 

            | bell(x) | < eps,                                       (1)
                                                                     
for all x \leq 0, and

            | bell(1) -1 | < eps,                                    (2)
         
for all x \geq 1; furthermore, 

          bell(x)= \sum_{i=1}^{ncoss}
          coefs(i) * cos((i-1)*(x*belllen+shift))                    (3)

The claim to fame of the function bell :[0,1] \to R^1 is that 
it has (approximately) the lowest frequency content of all
functions satisfying the conditions (1),(2).

The coefficients coefs (as well as parameters ncoss, 
shift, bellen) must have been generated by a prior call to 
the subroutine bellret, also in the file bellret.f


60.
The file pnmrts.f contains subroutines pnmini, pnmeva 
for the evaluation of NORMALIZED Associated Legendre Functions 
P_n^m at arbitrary points on the interval [-1,1], and the 
subroutine pnmrts that finds the roots of functions P_n^m.

61. 
The file leamatr.f contains the subroutines leamatr, leamatll,
leamatrr solving the least squares approximation problems for 
matrices from both sides, from the left, and from the right, 
respectively.

62.
The file qneval.f contains the subroutine qneval evaluating 
(for the user-specified real x and integer n) the Legendre 
functions Q_0, Q_1, ..., Q_n.

63. 
The file cgrmsol.f contains the subroutine cgrmsol, solving a 
system of complex linear equations (with a signle right-hand
side). It is a very short, primitive, and robust code.

64.
The file qrdpivot.f contains the subroutine qrdpivot constructing
the left singular vectors (and also, singular values) of the 
user-supplied matrix. The claim to fame of this subroutine is that
it requires less storage and (usually) less CPU time than 
svdpivot.

65.  
The file allqrbld.f contains the subroutine allqrbld compressing
a bunch of user-supplied functions on the interval [a,b].  The
functions are supplied in the form of a subroutine evaluating them.
The subroutine creates an appropriate discretization of the interval
[a,b] (a nested Legendre discretization), and evaluates the uaer-
supplied functions at the nodes of the discretization. Then it 
applies the SVD to the result (via the subroutine qrdpivot), and
declares the left singular vectors to be the compressed version of 
the user-supplied functions.

66. 
The file grmsol.f contains the subroutine grmsol, solving a 
system of real linear equations (with a signle right-hand
side). It is a very short, primitive, and robust code.

67.
The file trplopen.f contains the subroutines trplopen, trpldot,
trplpoint, trplline, trplwrt, used for plotting data in GNUPLOT
-readable format. The subroutine in this file are much more
flexible than those in the file lotapklot.

68.
The file p2inret.f contains subroutine for the construction of
interpolation nodes (up to order 10) on the standard equilateral
triangle, and also for the evaluation of orthogonal polynomials
(also up to order 10) on that triangle. Pleae note that the 
subroutines in this file use data they retrieve from the fairly 
large subroutine p2data3 stored separately in the file p2data3.f

69. 
The file p2eval contains the subroutine p2eval, evaluating at 
points in R^2 various polynomials orthogonal on the standard 
equilateral triangle. This is a self-contained file, in the sense 
that it contains the subroutine p2data3, storing the precomputed 
coefficients to be used in the evaluation of orthogonal polynomials.
Please note that this file does NOT contain a testing code.

70. 
The file triaadam.f contains the subroutine triaadam that evaluates 
the integral of a vector-valued function on a triangle via an adaptive 
nested tensor product Gaussian rule. The file also contains the debugging
code.

71.
The file adapgaum.f contains the subroutine adapgaum that evaluates 
the integral of a vector-valued function on an interval via an adaptive 
nested Gaussian rule. The file also contains the debugging
code.

72.  
The file alltrint.f contains the subroutine alltrint that for a
user-specified point (x,y,z) in R^3, evaluates the integrals

      \int x^i * y^j / sqrt( (x-x0)^2 + (y-y0)^2) dx dy          (1)

       over the standard triangle vith the vertices

       (0, 2/sqrt(3)), (-1,-1/sqrt(3)), (1,-1/sqrt(3)),          (2)

       and i+j .le. 12

73. The file mveccomp.f contains the subroutine mevccomp that 
compresses via recursive Gram-Schmidt the user-provided collection
of vectors; the vectors can (and should be, for the sake of efficiency)
blocked into matrices of reasonable size. The subroutine does not
create the whole mess at once; the user provides the subroutine
that builds the vectors to be compressed, one block at a time. The
principal intended use of this subroutine is in the least squares 
environment.

74.
The file mvecsvd.f contains a fairly civilized version of the 
subroutine mvecsvd that constructs the left singular vectors 
and the singular values of a user-supplied matrix. The matrix is
supplied via a subroutine, and never constructed explicitly.
When the matrix is large but the rank is small, the subroutine
eats very little memory. This file supersedes the file mveccomp.f

75.
The file allsvbld.f contains the subroutime allsvbld constructing the
left singular functions of the set of of user-specified functions on
the interval [a,b], and returns to the user nested Legendre expansions
for the first ncols singular functions; all elements with higher
number than ncols have weights associated with them that are less than
eps (user-specified). The subroutine also returns to the user the
singular values associated with the said singular functions. Once this
subroutine has been called, any one of the obtained singular functions
(together with its derivative) can be evaluated by a call to the
subroutine nesteva (see).

76.  
The file allsvcmb.f contains the subroutines allsvcmb, allsvbld
constructing the left singular functions of the set of of
user-specified functions on the interval [a,b], and return to the user
nested Legendre expansions for the first ncols singular functions; all
elements with higher number than ncols have weights associated with
them that are less than eps (user-specified). The subroutines also
return to the user the singular values associated with the said
singular functions. Once one of these subroutines has been called, any
one of the obtained singular functions (together with its derivative)
can be evaluated by a call to the subroutine nesteva (see).

Please note that the subroutine allsvcmb is a highly specialized
version of the subroutine allsvbld (see). In most situations where one
would think of using this subroutine, one should actually use
allsvbld!!! The subroutine allsvcmb should be used when the input
functions \phi_i to be SVD'ed are naturally represented as products

            \phi_i (x) = \psi_i (x) + \khi_i (x),                    (1)

and the number of functions \psi_i is much greater than the
(numerical) rank of their set. The subroutine starts with constructing
the SVD of the collection of functions \psi_i; then it multiplies the
obtained singular vectors by the functions \khi_i, and constructs the
SVDs of the product (scaled by the corresponding singular values, of
course). When the number of functions \psi_i is much greater than the
(numerical) rank of their set, this subroutine tends to be
considerably faster than allsvbld.

Also, note that the file allsvcmb.f makes the fime allsvbld obsolete.

77.
The file homot1.f contains the subroutines homot0, homot1,
designing Generalized Gaussian quadratures
for a user-specified set of basis functions. The user 
supplies the functions to be integrated via the subroutine
fun3 (see below). The exact conditions under which this
subroutines are guaranteed to work are complicated and given
by Markov's theorems (totally positive kernels, Chebychev
systems, and other nonsense). On the other hand, these
subroutines often work when there appear to be no known
theorems to guarantee their successful operation. In other 
terms, THIS IS NOT A BLACK BOX; HUMAN INTERVENTION IS 
OFTEN NEEDED.

78.
The file qhlogev2.f contains the subroutine qhlogev2 that for the
user-supplied real x and integer n, constructs the coefficients
hilcoefs, coefsqua, coefslog of the linear forms converting the values
of a function f at n Gaussian nodes on the interval [a,b] into the
values at the point x of the Hilbert, quadrupole, and logarithm
transforms of f, respectively. This subroutine uses as its input the
array w, constructed by the initialization subroutine qhlogini.

79.
The file perilagr.f contains the subroutines perilagr and eqlagrev 
interpolating (via standard Lagrange interpolation) a user-specified 
function from the user-provided equispaced grid to the user-supplied 
point x on the line. perilagr is the version of the scheme designed
for the periodic case; in the non-periodic case, eqlagrev should be
used.


80.
The file jaceig.f contains the subroutine jaceig that uses Jacobi 
iterations to find the spectrum of a real symmetric matrix. It is a 
primitive and robust code. It is also mercifully short. 

81.
The file cjaceig.f contains the subroutine cjaceig that uses Jacobi 
iterations to find the spectrum of a Hermitian matrix. It is a 
primitive and robust code. It is also mercifully short. 

82.
The file proquadr.f contains the subroutine proquadr constructing
 quadratures for band-limited functions on the interval [-1,1]. 
Specifically, it constructs an npts-point quadrature integrating 
exactly the first npts*2 prolate functions corresponding to the 
user-specified band-limit c; npts is chosen to be the smallest 
integer such that the prolate eigenvalue number npts*2+1 is less 
than eps.


83.  
The file slepapin.f contains the user-callable subroutines slepapin,
slepapev, slepextr, sgainleg for the construction and evaluation of
band-limited approximations to functions on the interval [-1,1],
having (approximately) the user-specified supergain. The other obvious
use of these subroutines is to extend functions defined on the
interval [-1,1] from that interval to a bigger one. Please note that
the approximation constructed is always in the standard L^2 sense (no
user- specified weight functions, etc), since it uses standard prolate
functions, as opposed to taylor-designed basis function (that will
come later). Also, please note that the code in slepapin is fairly
civilized, but not very carefully tested.

84.
The file lesqinit.f contains the subroutines lesqinit, lesqeval
constructing supergain-controlled least squares approximations of
functions (vectors, whatever) by other functions (vectors, etc.).  It
is a fairly civilized code, but at this time the testing performed has
been fairly perfunctory. 11.9.99

The subroutine lesqevac has been added to the mix, permitting complex
functions to be approximated (though the approximating functions
have to be real); several service subroutines have also be added in
order to implement lesqevac. All of these have undergone only 
very preliminary testing.  7.12.00

85.
The file sllinit.f contains the user-callable subroutines sllinit,
slleval approximating a user-supplied function by a band-limited 
function, having the user-prescribed amount of supergain. 12.11.99

86.
The file prolbelin contains the subroutines prolbelin, prolbell
constructing the prolate bell (possessing all the optimal properties).

87.  
The file prolinte contains the subroutine prolinte constructing
prolate function-based interpolation formulae with user-supplied
nodes.

88.  
The file prolfili.f contains a set of subroutines for the
interpolation and filtering of band-limited functions on both the
interval and the square.


89.
The file hank102.f contains the subroutine evaluating the Hankel 
functions H_0^1, H_1^1 for an arbitrary user-specified complex 
number z. This is a slightly accelerated and prettified version 
of the old hank101, whose principal claim to fame is that it is valid 
on the whole complex plane, and is reasonably accurate (14-digit 
relative accuracy) and reasonably fast. -6.12.00

90.
The file hank103.f contains the subroutines hank103, hanks103.
The subroutine hank103 is a modification of the subroutine hank102, 
providing to the user the additional option of evaluating the 
functions h0, h1 scaled by the (complex) coefficient e^{-i \cdot z}. 
The subroutine hanks103 combines the subroutine hank103 with the
standard recursion to evaluate the first n+1 Hankel functions of
the argument z, with both n and z user-provided parameters. -6.13.00

91.
The file projfilt.f contains the subroutine projfilt applying a
projecting filter to a user-supplied (complex-valued) function
fin. Fin is assumed to be periodic in the standard sense, i.e. fin is
tabulated at n nodes, and it is assumed that fin(n+1)=fin(1)

92.
The file cjseval.f contains the subroutines cjseval, rjseval for
the evaluation of Bessel functions of complex and real arguments, 
respectively.

93.
The file protaper contains various subroutines constructing taper
functions, "universal" quadratures, etc. based on the prolate
function of order zero. This file is self-sufficient, in the sense
that it calculates its own prolate functions, and does not use 
any subroutines from the file prolcrea.f 

94.
The file unidiscr.f contains the subroutine unidiscr creating 
"universal" discretizations of the interval [-1,1], together with 
accompanying machinery. Please note that the principal output of
this subroutine is a FORTRAN-written file, to be read by the 
subroutines unicoret, uniptret and otherwise used by the subroutine
usisegm, to be found in the file unisegm.f. The other user-callable
subroutine in this file is uninodes, producing the nodes of the
universal discretization on the interval [-1,1], together with the
corresponding quadrature weights.

95.
The file unisegm.f contains the subroutines unisegm, unicoret, 
uniptret, reading from disk the data stored there previously by the 
subroutine unidiscr stored in the file unidiscr.f, and using these 
data to create discretizations of integral operators.

96.
The file hermexps contains the subroutine hermexps constructing
the Hermite nodes on R^1, and the weights for the corresponding 
order n quadrature. It also constructs the matrix v converting 
the coefficients of an Hermite expansion into its values at the 
n Hermite nodes and its inverse u, converting the values of a 
function at n Hermite nodes into the coefficients of the 
corresponding Hermite series. No attempt has been made to make 
this code efficient, but its speed is normally sufficient, and it 
is mercifully short.


97.
The file prq2side contains the subroutines prq2side, prq1side,
returning to the user the two and one-sided universal quadratures, 
respectively. These subroutines do not construct the quadratures, 
but have several of these quadratures stored  in them; thus, the
selection of such quadratures is limited. For subroutines constructing
such quadratures for reasonably general combinations of accuracies 
and numbers of ndes, see the file probexps.f

98.
The file probexps.f contains the subroutine probexps that is a
"universal " analogue of the subroutine legeexps for the Legendre 
polynomials, together with various companion subroutines. 

99.
The file adapgaub contains the subroutine adapgaub that is a very
similar to the subroutine adapgaus. However, adapgaub can (at the
user's discretion) replace Gaussian quadratures with "universal"
ones on the subintervals near the ends of the interval of integration.
For integrands with singularities at (or near) the ends of the interval
of integration, this results in a dramatic reduction in the required
number of function evaluations. THIS SUBROUTINE IS EXPECTED TO SUPERSEDE
THE SUBROUTINE ADAPGAUS COMPLETELY.

100.
The file qu10by20 contains the subroutines qu10by20 returning to the 
user the nodes and weights of 10 "universal" quadratures on the 
interval [-1,1], corresponding to the 10 Gaussian nodes. The i-th 
quadrature integrates exactly all functions f of the form

         f(x)=phi(t)+psi(t) * log(|x_i-t|),                     (1)

with phi, psi two polynomials of degree 19, and x_i the i-th Gaussian 
node. The file also contains several companion subroutines for qu10by20.


101.
The file fourlog contains the subroutines fourlog, sinint, evaluating
the integral 

        int_{-rlim}^{rlim} cos(a*x) * log(|x|) dx               (1)

and the integral 

        sinint(x)=\int_0^x sin(x)/x dx,                          (2)

respectively.


102.
The file o1r20qua.f contains a collection of subroutines for the 
handling of integral integral equations of potential theory in two
dimensions on short sides of polygons. In this file, the side of
the polygon is discretized into 20 nodes, and assumes that on this
side everything is a polynomial of order 9. Under such conditions,
the resulting accuracy is about 7 digits.

103.
The file o3rrdall.f contains the machinery for the construction of 
"near" quadratures for the handling of integral equations 
of potential theory in two dimensions on short sides of polygons. 
In this file, the side of the polygon is discretized into 20 nodes, 
and it is assumed that on this side everything is a polynomial of 
order 9. Under such conditions, the resulting accuracy is better
than 6 digits. Plese note that this file contains a lot of tabulated
data and no test driver.

104.
The file wholmatr.f contains the subroutine wholmatr constructing 
the matrix of interactions corresponding to the user-specified 
boundary and interaction law. The subroutine uses (inter alia) 
the code found in the file o3rrdall.f, and is supposed to contain 
many useful subroutines in an essentially final form.

105.
The file wholmatc.f contains the complex version of the code in the
file wholmatr.f. The subroutine wholmatc constructs
the matrix of interactions corresponding to the user-specified 
boundary and interaction law. The subroutine uses (inter alia) 
the code found in the file o3rrdall.f, and is supposed to contain 
many useful subroutines in an essentially final form.

106.
The file qrsolve.f contains the subroutines qrsolv,qrdecom, qrsolve, 
qrsolve_c,qrsolve_adj, qrcond qrsolv, qrdecom, qrsolve using a 
version of QR-decomposition to solve user-supplied systems of linear 
algebraic equations with real coefficients and right-hand that are
either real or complex.

107.
The file cqrsolve.f contains the subroutines cqrsolv, cqrdecom, cqrsolve 
that use a version of QR-decomposition to solve user-supplied systems
of linear algebraic equations with complex coefficients and right-hand 
sides.

108.
The file homot3.f contains the subroutines homot2, homot3, designing
Generalized Gaussian quadratures for a user-specified set of basis
functions. The user supplies the functions to be integrated via the
subroutines. The exact conditions under which this subroutines are
guaranteed to work are complicated and given by Markov's theorems
(totally positive kernels, Chebychev systems, and other nonsense). On
the other hand, these subroutines often work when there appear to be
no known theorems to guarantee their successful operation. In other
terms, THIS IS NOT A BLACK BOX; HUMAN INTERVENTION IS OFTEN NEEDED.
%%%%%
ALSO, PLEASE NOTE THAT THIS FILE IS INTENDED TO SUPERSEDE THE FILE
HOMOT1.F COMPLETELY!!!!
 
109.  
The file orthopol.f contains the subroutine orthopol constructing the
machinery for the evaluation of orthogonal polynomials on the standard
triangle in the plane, and the subroutine orthoeva that actually
evaluates such polynomials. The other user-callable subroutines in the
file are ortretr, ortconv, ortread. Plese note that the subroutine
orthopol is extremely crude: it does not even have memory management,
being able to evaluate orthogonal polynomials of orders up to 40. On
the other hand, the subroutine orthoeva is reasonably civilized.

110.
The file reprinit.f contains the subroutines reprinit, reprcon1, 
reprcon2, repreval that, for a user-specified collection of functions 
choose a (much smaller) subset of functions that span (to the precision 
eps) all of the original functions. Then, a set of nodes is chosen 
such that the collocation matrix of the chosen functions at the chosen
nodes admits a decent generalized (in the appropriate sense) inverse.
Subsequently, any linear combination of the original functions can be 
efficiently expressed in the form as a linear combnation of the functions
in the chosen subset. The only data required by the procedure are the 
values of the function to be expanded at the chosen nodes (of whom there
are not too many).

111. 
The file cmvecsvd.f contains the subroutine cmvecsvd constructing
the left singular vectors and the singular values of the user-specified 
rectangular matrix. The matrix is provided via the user-supplied 
subroutine. The principal output of the subroutine are the compressed 
version of the matrix and the vector of its singular values. The principal 
anticipated uses for this subroutine are in the Least squares and 
quadrature environments. This is a complex version of the subroutine
mvecsvd (see).

112. 
The file cnestexs.f contains the subroutine cnestxs finding a subdivision 
of the interval [a,b], such that on every subinterval, the user-provided 
functions (nfuns of them things) given by the subroutine fun is 
interpolated to precision eps by a k-point Legendre interpolation formula; 
then it evaluates on each subinterval the coefficients of the Legendre 
series approximation of f. The resulting approximation is gauaranteed 
to accuracy eps, and the output of this subroutine can be used by the 
subroutines cnesteva, cnestev2, cnestem (see) to evaluate the interpolated 
functions at arbitrary points on the interval [a,b]. The file contains 
several other user-callable subroutine, whose functions are closely
related to the function of cnestexs.

113.
The file calsvbld.f contains a collection of subroutines for the 
construction of singular value decompositions (or rather, of left 
singular vectors) of families of complex functions of a real variable, 
and for the performance of certain associated operations. The assumption 
is that the number of functions to be compressed is too great for them 
to be discretized and evaluated simultaneously; the user specifies them 
things via subroutines that evaluate them. TO A LARGE DEGREE, THIS IS A 
COMPLEX VERSION OF THE FILE allsvcmb (see)

114.
The file calsvcmb.f contains the subroutines calsvcmb, calsvbld
constructing the left singular functions of the set of user-specified 
complex functions on the interval [a,b], and return to the user
nested Legendre expansions for the first ncols singular functions; all
elements with higher number than ncols have weights associated with
them that are less than eps (user-specified). The subroutines also
return to the user the singular values associated with the said
singular functions. Once one of these subroutines has been called, any
one of the obtained singular functions (together with its derivative)
can be evaluated by a call to a subroutine of the cnesteva family (see).

Please note that the subroutine calsvcmb is a highly specialized
version of the subroutine calsvbld (see). In most situations where one
would think of using this subroutine, one should actually use
calsvbld!!! The subroutine calsvcmb should be used when the input
functions \phi_i to be SVD'ed are naturally represented as products

            \phi_i (x) = \psi_i (x) + \khi_i (x),                    (1)

and the number of functions \psi_i is much greater than the
(numerical) rank of their set. The subroutine starts with constructing
the SVD of the collection of functions \psi_i; then it multiplies the
obtained singular vectors by the functions \khi_i, and constructs the
SVDs of the product (scaled by the corresponding singular values, of
course). When the number of functions \psi_i is much greater than the
(numerical) rank of their set, calsvcmb tends to be considerably faster 
than calsvbld.

Also, note that the file calsvcmb.f makes the fime calsvbld obsolete.

115.
The file expoquad.f contains the subroutine expoquad designing special
-purpose quadrature for the integration of a certain type of functions 
frequently encountered in the numerical complex analysis. The functions
are defined on intervals located in C^2 and parallel to the real exis, 
and have the form 

 f(z)=exp(rlam*x) * phi(x),                                   (1)
   
with the parameter rlam on the (user-supplied) interval  [rlmin,rlmax], 
and x on the (user-supplied) interval [a,b], and phi a reasonably 
general integrable (complex-valued) function on the interval [a,b].

116.
The file quaevol1.f containe the subroutine quaevol1, constructing 
quadratures for a user-specified class of functions. The subroutine 
starts with constructing an SVD of a set of functions, and continues 
by constructing (or rather, attempting to construct) a Generalized 
Gaussian quadrature for the obtained singular vectors. In fact, the 
subroutine builds a Chebychev-type quadrature (i.e. an n-point 
quadrature for n functions), which is trivial, and attempts to reduce 
the number of nodes, one node at a time. If the subroutine failed to 
obtain a Generalized Gaussian Quadrature, it still has (usually) achieved 
something; in practice, the subroutine tends to get within one or two 
nodes from a Gaussian quadrature.

117. 
The file fileflush.f contains the subroutine fileflush, which flushes
a FORTRAN formatted output file. Specifically, it closes the file, 
opens it as a read/write file, and positions it at the end of the 
last record.

118.
The file cadapgaum.f contains the subroutine cadapgaum using the adaptive 
Gaussian quadrature to evaluate the integral of the user-supplied 
vector-valued function fun on the user-specified interval [a,b].
PLEASE NOTE THAT THIS IS A COMPLEX-VALUED VERSION OF THE SUBROUTINE 
ADAPGAUM (SEE).

119.
The file calsquads.f contains the subroutine calsquads producing 
quadratures for the class of complex functions (that has been compressed 
by "calsvcmb" or "calsvbld" previously) multiplied by a user-specified 
weight function "fun".
 
120.
The file ortho2eva.f contains a set of subroutines for the evaluation 
of othogonal polynomials on the standard triangle. It contains 2 
user-callable subroutines: ortho2eva (evaluates the polynomials)
and ortho2eva3 (evaluates the polynomials together with their 
derivatives). This code has been written by Zydrunas.

121.
The file quaevol2.f contains the subroutines quaevol1, quaevol2,
constructing quadratures for a user-specified class of functions. 
The subroutines starts with constructing an SVD of a set of functions, 
and continue by constructing (or rather, attempting to construct) a 
Generalized Gaussian quadrature for the obtained singular vectors. In 
fact, the subroutines build a Chebychev-type quadrature (i.e. an n-point 
quadrature for n functions), which is trivial, and attempt to reduce 
the number of nodes, one node at a time. If one of these subroutines 
failed to obtain a Generalized Gaussian Quadrature, it still has (usually) 
achieved something; in practice, the subroutines tend to get within one or 
two nodes from a Gaussian quadrature. The difference between the two
subroutines is that quaevol1 uses allsvbld to construct the SVD, while
quaevol2 uses allsvcmb.

122. 
The file arbgauss contains the subroutine arbgauss constructing
gaussian quadratures on the interval [a,b] with the weight given by the 
user-specified function. The weight function is fairly arbitrary. PLEASE
NOTE THAT THIS IS AN OLD CODE. IT IS ONLY ENTERED IN THE DOCUMENTATION
FILE NOW. 7.8.01.

123.
The file tables14.f contans a fairly final version of the subroutine 
tables14 returning to the user quadratures for the evaluation of
        integrals of the form

              int_{-1}^1 p(x) * s(x) dx,                                  (1)

        where p is a polynomial of order up to 19, and s a linear 
        combination of functions

              s_0(x) = 1,                                                 (2)
              s_1(x) = log(x**2+y**2),                                    (3)
              s_2(x) = d/dx (log(x**2+y**2))= x/(x**2+y**2),              (4)
              s_3(x) = d/dy (log(x**2+y**2))= y/(x**2+y**2),              (5)

        and 

              y \in [1.0d-6, \infty).                                     (6)

124.
The file tables14b.f contans a fairly final version of the subroutine 
tables14 returning to the user quadratures for the evaluation of
        integrals of the form

              int_{-1}^1 p(x) * s(x) dx,                                  (1)

        where p is a polynomial of order up to 19, and s a linear 
        combination of functions

              s_0(x) = 1,                                                 (2)
              s_1(x) = log(x**2+y**2),                                    (3)
              s_2(x) = d/dx (log(x**2+y**2))= x/(x**2+y**2),              (4)
              s_3(x) = d/dy (log(x**2+y**2))= y/(x**2+y**2),              (5)
              s_4(x) = d^2/dy^2 (log(x**2+y**2))= 
                     =  (x**2-y**2)/(x**2+y**2)**2                        (6)

        and 

              y \in [1.0d-6, \infty).                                     (7)

PLEASE NOTE THAT THESE QUADRATURES ARE DIFFERET FROM THOSE IN THE FILE
TABLES14.F IN THAT THEY HANDLE THE QUADRUPOLE POTENTIALS.

125.
The file univinterp.f contains machinery for the interpolation and 
integration on the interval [-1,1] of functions of the form (1) below. 
Specifically, this subroutine constructs a 30-node interpolation 
formulae for functions 

        P_n(x), 

        P_n (x) \cdot (1-x) \cdot log(1-x),
                                                                    (1)
        P_n (x) \cdot  \cdot(1+x) log(1+x),
                 
with n=0,1,2,...,9.  The interpolation is exact in exact arithmetic; 
in practice, 14-15 digits are obtained with double precision 
computations. 

126.
The file arrstore.f contains the subroutine arrstore that stores on 
the FORTRAN unit iw the user-supplied real *8 array arr of length n 
in a format easily converted into a set of FORTRAN data statements. 
If the array arr is longish, it will be broken up into pieces of 
length 180 elements (so that each data statement is no longer than 
90 cards). Please note that this file does NOT contain any debugging
testing code.

127.
The file univintall.f contains machinery for the interpolation and 
integration on the interval [-1,1] of functions of the form (1) below. 
Specifically, this subroutine constructs a 30-node interpolation 
formulae for functions 

        P_n(x), 

        P_n (x) \cdot (1-x) \cdot log(1-x),
                                                                    (1)
        P_n (x) \cdot  \cdot(1+x) log(1+x),
                 
with n=0,1,2,...,9.  The interpolation is exact in exact arithmetic; 
in practice, 14-15 digits are obtained with double precision 
computations. PLEASE NOTE THAT THIS FILE CONTAINS ALL OF THE ANALYTICAL
MACHINERY NEEDED TO CONSTRUCT SELF-INTERACTIONS OF A CHUNK IN THE
SISCRETIZATION OF INTEGRAL EQUATIONS OF SCATTERING THEORY IN TWO 
DIMENSIONS. This file is intended to supercede the file univinterp.f
completely.



128.
The file prolfltinte.f contains the subroutine prolfltinte producing
filtering and interpolation matrices for (non-periodic) band-limited 
functions on the interval [-1,1]. It also generates the requisite
discretization nodes and the corresponding weights. This is a 
companion file for the file prolfili (see).

129.
The file hank106.f contains the subroutine hank106 evaluating the 
Hankel functions H_0^1, H^1_1 of a complex argument, the argument 
living on a ray. The subroutine you are looking at is the 
initialization subroutine; the evaluation subroutine is hank106 
(see). THIS IS A HIGHLY SPECIALIZED SUBSTITUTE FOR HANK103 (SEE);
IT IS MUCH CRANKIER THAN HENK103, BUT IS ABOUT 4 TIMES FASTER THAN
THE LATTER.

130.
The file curvefilt.f contains the subroutine curvefilt filtering
and resampling the user-supplied curve. The curve is NOT assumed 
to be closed, and is supplied by the user in the form of a collection 
of points in the plane.

131.
This file cnestex_cheb.f contains the numerical machinery for the 
construction and evaluation of nested Chebychev expansions of complex 
functions of one real variable. There are 3 user-callable subroutines 
in the file: cnestex_cheb, cneste_cheb, cneste2_cheb. 

132.
The file expoapp20.f contains 5 user-callable subroutines and entries:
expoapp20, expo_ttww_ret20, expoapp63, expo_ttww_ret63,
expoextr. Following is a brief description of these entries.

expoapp20 - constructs an expansion of the user-supplied function 
    into a linear combination of 38 exponentials, of the form

       f(x)=\sum_{j=1}^{38} coesexp_j * exp(ima*c*tt_j),             (1)

    with c=20.

expo_ttww_ret20 - returns to the user the 38 prolate nodes tt 
    and their corresponding weights ww, discretizing to 16-digit
    (or so) accuracy band-limited functions with c=20 on the
    interval [-1,1]. 

expoapp63 - constructs an expansion of the user-supplied 
    function into a linear combination of 74
    exponentials, of the form

       f(x)=\sum_{j=1}^{74} coesexp_j * exp(ima*c*tt_j),             (2)

    with c=63.

expo_ttww_ret63 - returns to the user the 74 prolate nodes tt 
    and their corresponding weights ww, discretizing to 16-digit
    (or so) accuracy band-limited functions with c=63 on the
    interval [-1,1]. 

expoextr - given a collection of user-provided frequencies ttt 
    and (complex) coefficients vals(j),this subroutine evaluates 
    the expansion of the form 

       fout = sum_{j=1}^nnn 
              vals(i) \cdot e^{i \cdot c \cdot ttt(j) \cdot x}       (3)

    at the point x \in R^2; the subroutine also returns the 
    derivative of fout with respect to x

133.
The file expoexpand.f contains the subroutine expoexpand constructing
the coefficients coesexp of nn exponential expansion of a 
user-specified (complex-valued) function. The function is supplied to 
the subroutine via collection fs of its values tabulated at 40 nodes 
tt on the interval [-1,1]; the the user can obtain the latter by calling 
the entry exporetr of this subroutine. This is a simplified, improved,
and much more robust version of the code in the file expoapp20.

134.
The file exponexpand.f contains five user-callable subroutines:
exponexpand, exponeval, expoexpand, exporetr, expoeval (actually,
exporetr is an entry in the subroutine expoexpand). Together, these
entries handle approximation of functions on R^1 by exponentials.
This is a fairly civilized code, inter alia meant to supercede that
in the file expoexpand.f.

135.
The file inteeval135 contains the subroutine inteeval135 evaluating
the integral

        \int_{0}^{\pi}   {cos(m*theta) \over                           (1)
            ( sqrt(1+q*(sin (theta/2) )**2 ) ^k },

with k=1, 3, 5, and m \in [0,35]. PLEASE NOTE THAT THIS CODE SHOULD BE
TESTED WITH THE FILE inteeval135f90 and inteeval135, rather than with 
the file inteeval135 alone!.
 
136. The file chebtorat.f contains the subroutine chebtorat that
For a Chebychev polynomial T of order n, and a (complex ) number 
delta, constructs a set of 7*n points zsout in the complex plane, 
and the coefficients deltas(1), deltas(1),...,deltas(7), such that 
on the interval [-1,1] on the real axis, 

          delta * T(x) = 1-1 / Q(x) + O(|delta|**8).                   (1)

The file also contains the subroutine tappreval evaluating the 
expansion (1) at arbitrary points in C.

137.
The file blind_bld.f contains the subroutine blind_bld constructing
a Helmholtz potential in the plane having a user-prescribed far-field
signiture in certain (user-prescribed) directions. The potential to be 
constructed is generated by a collection of charges located at 
user-prescribed positions in the plane, and the subroutine uses a
least-squares scheme to find the appropriate (complex) intensities. 
PLEASE NOTE THAT THIS SUBROUTINE IS NOT DESIGNED TO BE VERY FAST OR TO
USE MEMORY VERY EFFICIENTLY

138.
The file ncleastsq.f contains the user-callable subroutines: 
ncleastsq, ncleastsq2, cleamatr, ncleamatll, ncleamatrl. Following 
is a brief description of these subroutines. 

   ncleastsq - constructs a QR-type (not exactly QR) decomposition 
         of the input matrix A, to be used by the subroutine 
         ncleastsq2 for the solution of least squares problems of 
         the form 

                    A X \sim Y,                                        (1)

        and by the subroutine ncleamatrl for the least squares
        solution of the matrix equation

                    A(n,m) * X(m,k) = C(n,k).                          (2)

   ncleastsq2 - uses the QR-type (not exactly QR) decomposition
         of the matrix A produced by a prior call to nsleastsq (see)
         to solve in the least squares sense the linear system 

                    A X = RHS.                                         (3)

   cleamatr - solves in the least squares sense the matrix
        equation 

                     A(k,l) * X(l,m) * B(m,n) = C(k,n),                (3a)

       where A, B, C are user-specified complex matrices, and X 
       is the matrix to be found. Note that the dimensionalities 
       of the matrices A,B,C, X are as general as could be 


   nclearmatll - solves in the least squares sense the matrix 
         equation

                    A(n,m) * X(m,k) = C(n,k),                          (4)

         without using any additional data (i.e. it performs all 
         factorizations itself).

   nclearmatrr - solves in the least squares sense the matrix 
         equation

                    X(n,m) * A(m,k) = C(n,k),                          (5)

         without using any additional data (i.e. it performs all 
         factorizations itself).

   ncleamatrl - solves in the least squares sense the 
        matrix equation

                    A(n,m) * X(m,k) = C(n,k),                          (6)

        using as input the matrix c and the array w obtained 
        via a preceding call to the subroutine ncleastsq

   PLEASE NOTE THAT THIS FILE IS EXPECTED TO SUPERCEDE THE FILES
   CLEASTSQ, CLEAMATR



139.
The subroutine rsortanyn sorting an array of integers is in the 
file rsortanyn.f, together with the debugging code. It is an exact 
real *8 version of the real code to be found in the file sortanyn.f.

140. 

        Thie file cskeleton.f contains four user-callable 
        subroutines: cskeleton, cskeleton3, cskel_basic, 
        cskel_utv. Following is a brief description of 
        these four subroutines

   cskeleton - for the user-specified complex n \times m - matrix a, 
       this subroutine finds a subset of ncols rows and a subset 
       of ncols columns, such that each of the three submatrices
       dimensined rowsout(ncols,m), colsout(n,ncols), 
       aout(ncols,ncols) provides a good approximation to the 
       original matrix a. Then, it constructs the matrices
       expand, eval, (approximately) solving the equations 

        aout(ncols,ncols) * expand(ncols,m) = rowsout(ncols,m),     (1)

        colsout(n,ncols) * expand(ncols,m) = a(n,m),                (2)

        and

        eval(n,ncols) * aout(ncols,ncols)  = colsout(n,ncols),      (3)

        eval(n,k) * colsout(k,m) = a,                               (4)

        respectively
        
   cskeleton3 - has the purpose in life identical to that of
        cskeleton; has an additional input parameter iwhich,
        permitting the user to suppress some of the calculations
        that are always performed by cskeleton. Such suppression
        can be used to save both CPU time and storage.

   cskel_basic - this is a simplified and accelerated version
        of the subroutines cskeleton3, cskeleton. the integer 
        vectors irows, icols are its only output: it does not 
        create any of the matrices rowsout, colsout, eval, expand, 
        whatever.
      
   cskel_utv - decomposes the user-supplied matrix a(n,m) in 
        the form

        a(n,m)=u(n,ncols) * t(ncols,ncols) * ( v(m,ncols) )^*,       (1)

       with the columns of the matrices u, v orthonormal, and 
       t a lower triangular matrix

141.
The file merges_1chunk.f contains three user-callable subroutines:
bandgeo, skels_1chunk, merges_1chunk; the subroutine bandgeo
performs preliminary geometrical processing of a band to be 
simulated; the subroutine skels_1chunk skeletonizes one chunk on
the finest level of subdivision (both incoming and outgoing 
skeletons are produced); the subroutine merges_1chunk merges
skeletons of two borother-chunks, producing the skeleton for
the daddy. 6.20.02

142.
The file scattmatr_merge2.f contains the subroutines scatmatr0, 
scattmatr_merge2; the subroutine scatmatr0 constructs the scattering 
matrix (in its four incarnations) for one chunk on level 0. It uses 
data (hopefully) created by a prior call to the subroutines bandgeo, 
skels_1chunk.

The subroutine scattmatr_merge2 merges the user-supplied scattering 
matrices s1, s2 of two brother-chunks, obtaining the scattering 
matrix sout of the daddy chunk. It also produces the splitting 
matrices split1, split2, "exchange matrices" exchg11, exchg12, 
exchg21, exchg22, and the "merging matrices" cmerge1, cmerge2. 

143.
The file emseval.f contains three user-callable subroutines:
e1eval, emseval, alphaseval. Following is a brief 
description of all three.

  e1eval - for the user-provided argument x \in (0,\infty), 
        this subroutine evaluates (to relative precision
        13 - 14 digits when run in double precision, to
        20-21 digits when run in extended precision) the
        exponential integral 

 
        E_1(x) = \int_x^{\infty} e^{-t}/t dt               (1)

  emseval - for the user-provided argument x \in (0,\infty), 
        this subroutine evaluates (to relative precision
        13 - 14 digits when run in double precision, to
        20-21 digits when run in extended precision) the
        exponential integral s
 
        E_m(x) = \int_1^{\infty} e^{-t*x}/t^m dt,          (2)

        with m=1,2,...,n.

  alphaseval - for the user-provided argument x \in (0,\infty), 
        this subroutine evaluates (to the machine precision)
        the exponential integrals

 
        alpha_m(x) = \int_1^{\infty} e^{-t*x}*t^m dt,      (3)

        with m=0,1,2,...,n.

144.

The file d2nest_expand contains two user-callable subroutines: 
        d2nest_expand, d2nest_eval. Following is a brief 
        description of these two subroutines.

  d2nest_expand - constructs a nested Chebychev expansion
        of the user-supplied function in the user-supplied
        rectangle. The expansion is expected to be subesequently
        evaluated at arbitrary points within the rectangle via 
        the subroutine d2nest_eval (see). 
       
   d2nest_eval - uses the data generated by a preceding 
        call to the subroutine d2nest_expand to (approximately)
        evaluate the user-supplied function at the user-supplied
        point in R^2. 


145.
The file cd2nest_expand.f contains two user-callable subroutines: 
        cd2nest_expand, cd2nest_eval. Following is a brief 
        description of these two subroutines. PLEASE NOTE 
        THAT THESE ARE COMPLEX VERSIONS OF THE SUBROUTINES 
        D2NEST_EXPAND, D2NEST_EVAL (SEE)

  cd2nest_expand - constructs a nested Chebychev expansion
        of the user-supplied function in the user-supplied
        rectangle. The expansion is expected to be subesequently
        evaluated at arbitrary points within the rectangle via 
        the subroutine d2nest_eval (see). 
       
   cd2nest_eval - uses the data generated by a preceding 
        call to the subroutine d2nest_expand to (approximately)
        evaluate the user-supplied function at the user-supplied
        point in R^2. 

146.
The file r2dstrcr.f contains three user-callable subroutines: 
r2dstrcr, r2dretr, and r2d_points_in. Following is a brief 
description of these subroutines.

   r2dstrcr - constructs the "bastardized" binary
        tree, obtained from the standard quad-tree by the 
        introduction of "bastard" boxes, so that every box
        is subdivided into two sons (one of which might not
        exist). Thus, there are two types of boxes: square 
        and rectangular; the rectangular boxes are twice as 
        high as they are wide. Please note that square boxes
        have even-numbered levels of subdivision: 0 (main 
        box), 2, 4,... . Rectangular boxes have odd-numbered 
        levels of subdivision. 

   r2dretr - retrieves from array w (hopefully)
        produced via a preceding call to the subroutine 
        r2dstrcr (see) the data corresponding to the 
        user-specified box in the structure. 

   r2d_points_in - for the user-specified conventionally oriented
        rectangle (defined by the coordinates of its sides),
        this subroutine finds all particles from the 
        user-specified collection that are inside the soid
        rectangle; the subroutine uses the "bastardized" 
        tree structure constructed during the preceding call 
        to the subroutine r2dstrcr (see). 


147.
The file square_basis.f contains the subroutine square_basis.
For the user-supplied subroutine poteval, this subroutine
constructs the singular values and the coresponding singular
vectors of the operator of interactions between two "concentric" 
squares in the plane. The subroutine also returns the nodes 
on the two squares, the corresponding quadrature weights, and 
several other parameters.


148.
The file square_skeleton.f contains the subroutine square_skeleton
skeletonizing the interactions between two "concentric" squares 
in R^2. Its expected use is for the Helmholtz interactions, and 
this is the regime where it has been tested.


149.
The file rects_skeleton.f contains the subroutine rects_skeleton
skeletonizing the interactions between two "concentric" rectangles 
in R^2. Its expected use is for the Helmholtz interactions, and 
this is the regime where it has been tested.

150.
The file r2d_outer_skel.f contains two user-callable subroutines:
r2dstrcr_neighb and r2d_outer_skel. Following is a brief description 
of these subroutines.

  r2dstrcr_neighb - an augmentation of the subroutine r2dstrcr,
       constructing (in addition to the machinery built by r2dstrcr)
       lists of neighbors for all boxes in the structure.
  
  r2d_outer_skel - for the box number ibox, this subroutine uses the 
       various data constructed by a prior call to the
       subroutine r2dstrcr_neighb (see) to obtain the outer
       skeleton for the box number ibox. The outer skeleton
       consists of a bunch of boxes on various levels (higher
       than the level of ibox), plus a collection of "free"
       points. 

151.
The file r2d_outer_skels contains four user-callable subroutines:
        r2dstrcr_neighb, r2d_outer_skel, r2d_outer_skels,
        r2d_skels_retr.  Following is a brief description 
        of these subroutines.

   r2dstrcr_neighb - an augmentation of the subroutine r2dstrcr,
        constructing (in addition to the machinery built by r2dstrcr)
        lists of neighbors for all boxes in the structure.
   
   r2d_outer_skel - for the box number ibox, this subroutine 
        uses the various data constructed by a prior call to the
        subroutine r2dstrcr_neighb (see) to obtain the 
        UNCOMPRESSED outer skeleton for the box number ibox. 
        The outer skeleton consists of a bunch of boxes on various 
        levels (higher than the level of ibox), plus a collection of 
        "free" points. 

   r2d_outer_skels - constructs outer and inner skeletons for all
        "standard" boxes on all levels. The resulting skeletons 
        are expected to be used in the construction of actual 
        skeletons for actual boxes

   r2d_skels_retr - retrieves from the storage area a single 
        skeleton for a single box; hopefully, the said skeleton 
        had been stored in the said area via a preceding call
        to r2d_outer_skels

PLEASE NOTE THAT THIS FILE COMPLETELY SUPERSEDES THE FILE 
r2d_outer_skel.f

152.
The file r2d_allskels.f contains 4 user-callable subroutines:
        r2d_allskels, r2d_evalexpands, r2d_store, r2d_retr. 
        Following is a brief description of these subroutines.

   r2d_allskels - first, this subroutine constructs the "bastardized" 
        binary tree, obtained from the standard quad-tree by 
        the introduction of "bastard" boxes, so that every 
        box is subdivided into two sons (one of which might 
        not exist). 

        Second, this subroutine constructs all skeletons on all 
        levels. The principal output of this subroutine is the 
        array w, containing the structure created by the subroutine:
        definitions of boxes, their corners and centers, various 
        types of control information, and inner and outer skeletons

   r2d_evalexpands - constructs the matrices eval, expand 
        for all boxes on all levels. It also uses the skeletons
        of all boxes constructed via a preceding call to the
        subroutine r2d_allskels (see) to construct the arrays 
        irows, icols.

   r2d_store - a call to this subroutine stores in the storage 
       area store one of (many) possible arrays pertaining to the
       chunk number ibox. 

   r2d_retr - a call to this subroutine retrieves from the storage 
       area store one of (many) possible arrays pertaining to 
       the chunk number ibox. 
        
154.
The file r2d_allskels.f contains the subroutine r2d_allskels,
        This subroutine constructs all of the matrices involved 
        in the sparse inversion of the two-dimensional "scattered
        points" scattering problem. There are two user-callable 
        subroutines in this file: r2d_allskels, r2d_finretr.
        Following is a brief description of these subroutines.

   r2d_allskels - first, this subroutine constructs the 
        "bastardized" binary tree, obtained from the standard 
        quad-tree by the introduction of "bastard" boxes, so that 
        every box is subdivided into two sons (one of which might 
        not exist). Thus, there are two types of boxes: square 
        and rectangular; the rectangular boxes are twice as 
        high as they are wide. Please note that square boxes
        have even-numbered levels of subdivision: 0 (main 
        box), 2, 4,... . Rectangular boxes have odd-numbered 
        levels of subdivision. The obtained structure is accessed
        via the subroutine r2dretr (see) returning to the user
        various types of data associated with the user-specified
        box in the structure.

        Second, this subroutine constructs all skeletons on all 
        levels. The said skeletons are also stored in the array
        w, from which array they are retrieved via the subroutine
        r2d_finskels_retr (see)

        Third, this subroutine constructs the sparse inverse of
        the interaction matrix on the whole structure, in the
        form of scattering, merging, splitting, and exchange
        matrices. All of these matrices are stored in the array
        w, from where they can be retrieved by calling the 
        subroutine r2d_finretr (see).

        PLEASE NOTE THAT WHILE THE SUBROUTINES R2D_FINSKELS_RETR,
        R2D_FINRETR ARE (AS OF 10.25.02) A PART OF THIS FILE, THE 
        SUBROUTINE R2DRETR IS A PART OF THE FILE r2dstrcr.f. 
        ATTEMPTS TO FIND A SUBROUTINE IN A FILE WHERE IT IS NOT 
        ARE UNLIKELY TO BE SUCCESSFUL.

   r2d_finretr -  retrieves from the main storage area w one 
       of (many) arrays pertaining to the chunk number ibox; 
       it is fervently hoped that the information to be 
       retrieved had been stores in the array w during a 
       preceding call to the subroutine r2d_allmatrs (see). 
       The chunk number ibox can be a chunk on any level; the 
       type of the array to be retrieved is determined by the 
       parameter kind.


155. 
The file lensemake.f contains two user-callable subroutines:
       lensemake and lenseeval. Following is a brief description
       of these subroutines.

  lensemake - constructs a solution of the Helmholtz
       equation resembling a potential focused by a lense. 
       More specifically, it constructs the coefficients of an 
       H-expansion such that the potential represented by the 
       said expansion (together with its gradient) is small on 
       the circle of radius r, except for a small spot around
       theta=pi/2. The obtaioned expansion and its derivative 
       with respect to r can be evaluated (inter alia) by calling
       the subroutine lenseeval (see)

  lenseeval - evaluates at the user-supplied point
        (x,y) in R^2 the potential represented by an H-expansion 
        with the (user-supplied) coefficients coefs; the
        derivative of the said potential with respect to the 
        radius sqrt(x**2+y**2) is also retrned.

156. 
The file hole_skel.f contains the subroutine hole_skel, skeletonizing
        one hole. It returns both the inner and the outer skeletons; 
        the skeletons are returned both as nodes in R^2 and as maps 
        of their locations in array z of user-supplied nodes. Ths 
        subroutine also returns the matrices expand, eval. The matrix 
        expand "expands"  a charge distribution defined on all 
        elements of which the hole consists into an equivalent charge 
        distributions on the skeleton. The matrix eval converts 
        values of an incoming potentials defined on the skeleton 
        into the values of the same incoming potential on elements 
        of the hole. Finally, the subroutine returns the artificial 
        matrix of interactions artif. EXPLANATION: if the hole is 
        replaced with its skeleton, and the matrix of interactions 
        between the elements of the hole is replaced with the matrix 
        artif, the scattering matrix of the hole (with respect to the 
        rest of the structure) is not changed, to precision eps.

157.
The file cmvecpsvd.f contains the subroutine cmvecpsvd constructing 
        the left singular vectors, the singular values and the integer 
        array of the most "significant" columns of the user-specified 
        rectangular matrix

        The matrix is provided via the user-supplied subroutine
        matfun (see below). The principal output of the subroutine 
        are the compressed version of the matrix (returned in the 
        array a, and the vector of its singular values. The principal 
        anticipated uses for this subroutine are in the Least squares 
        and quadrature environments.

        PLEASE NOTE THAT THIS SUBROUTINE IS ALMOST IDENTICAL TO THE 
        SUBROUTINE CMVECSVD; THE ONLY DIFFERENCE IS THE PRESENCE OF 
        THE INTEGER OUTPUT ARRAY IPIVOTS, CONTAINING THE SEQUENCE 
        NUMBERS OF THE "SIGNIFICANT" COLUMNS


158.
The file mvecpsvd.f contains the subroutine mvecpsvd constructing 
        the left singular vectors, the singular values and the integer 
        array of the most "significant" columns of the user-specified 
        rectangular matrix

        The matrix is provided via the user-supplied subroutine
        matfun (see below). The principal outputs of the subroutine 
        are the compressed version of the matrix (returned in the 
        array a, and the vector of its singular values. The principal 
        anticipated uses for this subroutine are in the Least squares 
        and quadrature environments.

        PLEASE NOTE THAT THIS SUBROUTINE IS ALMOST IDENTICAL TO THE 
        SUBROUTINE MVECSVD; THE ONLY DIFFERENCE IS THE PRESENCE OF 
        THE INTEGER OUTPUT ARRAY IPIVOTS, CONTAINING THE SEQUENCE 
        NUMBERS OF THE "SIGNIFICANT" COLUMNS
       
159.
The file rskeleton.f contains four user-callable subroutines: 
        rskeleton, rskeleton3, rskel_basic, rskel_utv. Following 
        is a brief description of these four subroutines

   rskeleton - for the user-specified real n \times m - matrix a, 
       this subroutine finds a subset of ncols rows and a subset 
       of ncols columns, such that each of the three submatrices
       dimensined rowsout(ncols,m), colsout(n,ncols), 
       aout(ncols,ncols) provides a good approximation to the 
       original matrix a. Then, it constructs the matrices
       expand, eval, (approximately) solving the equations 

        aout(ncols,ncols) * expand(ncols,m) = rowsout(ncols,m),     (1)

        colsout(n,ncols) * expand(ncols,m) = a(n,m),                (2)

        and

        eval(n,ncols) * aout(ncols,ncols)  = colsout(n,ncols),      (3)

        eval(n,k) * colsout(k,m) = a,                               (4)

        respectively
        
   rskeleton3 - has the purpose in life identical to that of
        rskeleton; has an additional input parameter iwhich,
        permitting the user to suppress some of the calculations
        that are always performed by rskeleton. Such suppression
        can be used to save both CPU time and storage.

   rskel_basic - this is a simplified and accelerated version
        of the subroutines rskeleton3, rskeleton. the integer 
        vectors irows, icols are its only output: it does not 
        create any of the matrices rowsout, colsout, eval, expand, 
        whatever.
      
   rskel_utv - decomposes the user-supplied matrix a(n,m) in 
        the form

        a(n,m)=u(n,ncols) * t(ncols,ncols) * ( v(m,ncols) )^*,       (1)

       with the columns of the matrices u, v orthonormal, and 
       t a lower triangular matrix

160.
The file nrleastsq.f contains 7 user-callable subroutines: 
         nrleastsq, nrleastsq2, nrleamatlr, rleamatr, nrleamatrr, 
         nrleamatll, nrleamatrl. Following is a brief description 
         of these subroutines. 

    NOTE: PLEASE NOTE THAT THE SUBROUTINE nrleamatrl IS ALMOST NOT 
         A USER-CALLABLE ONE. WHILE IT IS CONCIEVABLE THAT SOME 
         USERS WILL FIND IT USEFUL, MOST ARE EXPECTED TO BE HAPPY
         WITH THE SUBROUTINES nrleamatrr, nrleamatll, TO THE EXTENT 
         THAT THEY ARE HAPPY WITH ANYTHING IN THIS FILE.
       
   nrleastsq - constructs a QR-type (not exactly QR) decomposition 
         of the real input matrix A, to be used by the subroutine 
         nrleastsq2 for the solution of least squares problems of 
         the form 

               A X \sim Y,                                             (1)

        and by the subroutine nrleamatrl for the least squares
        solution of the matrix equation

         A(n,m) * X(m,k) = C(n,k).                                     (2)

   nrleastsq2 - uses the QR-type (not exactly QR) decomposition
         of the matrix A produced by a prior call to nrleastsq (see)
         to solve in the least squares sense the linear system 

                    A X = RHS.                                         (3)

   rleamatr - solves in the least squares sense the matrix
        equation 

                     A(k,l) * X(l,m) * B(m,n) = C(k,n),                (3a)

       where A, B, C are user-specified real matrices, and X 
       is the matrix to be found. Note that the dimensionalities 
       of the matrices A,B,C, X are as general as could be 

   nrleamatlr - solves in the least squares sense the matrix
        equation 

                     A(k,l) * X(l,m) * B(m,n) = C(k,n),                (3b)

       where A, B, C are user-specified real matrices, and X 
       is the matrix to be found. Note that the dimensionalities 
       of the matrices A,B,C, X are as general as could be 

       PLEASE NOTE THAT THIS SUBROUTINE CONSTRUCTS GENERALIZED     
       INVERSES OF THE MATRICES A, B INDEPENDENTLY FROM EACH OTHER.
       THUS, A CONSPIRACY IS POSSIBLE THAT WILL CAUSE SOME OF THE
       ELEMENTS OF THE MATRIX X TO BE OF THE ORDER ||C||/EPS^2,
       AS OPPOSED TO THE PROPER ||C||/EPS. CLEARLY, THIS WILL NOT
       HAPPEN IF THE MATRICES A,B (WHILE RECTANGULAR) HAVE NO SMALL
       SINGULAR VALUES. THE PURPOSE OF THIS SHORT-CUT IS TO AVOID 
       THE CONSTRUCTION OF SINGULAR VALUE DECOMPOSITIONS OF A, B;
       THE LATTER PROCEDURE IS KNOWN TO BE POTENTIALLY EXPENSIVE 
       AND EVEN UNRELIABLE. THE PROPER PROCEDURE USING THE SVDs
       IS PERFORMED BY THE SUBROUTINE RLEAMATR (SEE)

   nrlearmatll - solves in the least squares sense the matrix 
         equation

                    A(n,m) * X(m,k) = C(n,k),                          (4)

         without using any additional data (i.e. it performs all 
         factorizations itself).

   nrlearmatrr - solves in the least squares sense the matrix 
         equation

                    X(n,m) * A(m,k) = C(n,k),                          (5)

         without using any additional data (i.e. it performs all 
         factorizations itself).

   nrleamatrl - solves in the least squares sense the 
        matrix equation

                    A(n,m) * X(m,k) = C(n,k),                          (6)

        using as input the matrix c and the array w obtained 
        via a preceding call to the subroutine nrleastsq
        (see).

161.
The file quad_retr6.f contains the subroutine quad_retr6, returning
        to the user the 6-digit quadrature for the low-frequency 
        Helmholtz FMM with the Helmholtz coefficient rk in the 
        rectangle

        Re(rk) \in [5,64], Im(rk) \in [0,1]                         (1)
       
        in the complex plane. Actually, the quadratures are
        applicable in the region

        Re(rk) \in [5,64], Im(rk) \in [0,\infty],                   (1a)
       
        but are not very optimal in it. More specifically, the 
        quarature returned by this subroutine evaluates integral 
        of the form

         p * \int_0^{\infty} t*exp(-p*x*cd)/cd*J_0(p*y*t) dt        (2)

        with 

                   cd=sqrt(t^2-(rk/p)^2),                           (3)

        and p=Re(rk).

        PLEASE NOTE THAT BOTH THE NODES AND THE WEIGHTS OF THE
        RETURNED QUADRATURE ARE COMPLEX!! 

162.
The file quad_retr3.f contains the subroutine quad_retr3, returning
        to the user the 3-digit quadrature for the low-frequency 
        Helmholtz FMM with the Helmholtz coefficient rk in the 
        rectangle

        Re(rk) \in [5,64], Im(rk) \in [0,1]                         (1)
       
        in the complex plane. Actually, the quadratures are
        applicable in the region

        Re(rk) \in [5,64], Im(rk) \in [0,\infty],                   (1a)
       
        but are not very optimal in it. More specifically, the 
        quarature returned by this subroutine evaluates integral 
        of the form

         p * \int_0^{\infty} t*exp(-p*x*cd)/cd*J_0(p*y*t) dt        (2)

        with 

                   cd=sqrt(t^2-(rk/p)^2),                           (3)

        and p=Re(rk).

        PLEASE NOTE THAT BOTH THE NODES AND THE WEIGHTS OF THE
        RETURNED QUADRATURE ARE COMPLEX!! 


163.
The file quad_retr9.f contains the subroutine quad_retr9, returning
        to the user the 9-digit quadrature for the low-frequency 
        Helmholtz FMM with the Helmholtz coefficient rk in the 
        rectangle

        Re(rk) \in [5,64], Im(rk) \in [0,1]                         (1)
       
        in the complex plane. Actually, the quadratures are
        applicable in the region

        Re(rk) \in [5,64], Im(rk) \in [0,\infty],                   (1a)
       
        but are not very optimal in it. More specifically, the 
        quarature returned by this subroutine evaluates integral 
        of the form

         p * \int_0^{\infty} t*exp(-p*x*cd)/cd*J_0(p*y*t) dt        (2)

        with 

                   cd=sqrt(t^2-(rk/p)^2),                           (3)

        and p=Re(rk).

        PLEASE NOTE THAT BOTH THE NODES AND THE WEIGHTS OF THE
        RETURNED QUADRATURE ARE COMPLEX!! 


163a.
The file quad_retr12.f contains the subroutine quad_retr12, returning
        to the user the 12-digit quadrature for the low-frequency 
        Helmholtz FMM with the Helmholtz coefficient rk in the 
        rectangle

        Re(rk) \in [5,64], Im(rk) \in [0,1]                         (1)
       
        in the complex plane. Actually, the quadratures are
        applicable in the region

        Re(rk) \in [5,64], Im(rk) \in [0,\infty],                   (1a)
       
        but are not very optimal in it. More specifically, the 
        quarature returned by this subroutine evaluates integral 
        of the form

         p * \int_0^{\infty} t*exp(-p*x*cd)/cd*J_0(p*y*t) dt        (2)

        with 

                   cd=sqrt(t^2-(rk/p)^2),                           (3)

        and p=Re(rk).

        PLEASE NOTE THAT BOTH THE NODES AND THE WEIGHTS OF THE
        RETURNED QUADRATURE ARE COMPLEX!! 

164.
The file helmholtz_quadrs_get.f contains the subroutine
helmholtz_quadrs_get whose purpose in life is the construction 
of quadratures to be used in the FMM for the Helmholtz equation
at low frequencies.

165.
The file expocheb_prod.f contains the subroutines expocheb_prod 
evaluating integrals of the form

       \int_a^b T_n(x) \cdot e^{z \cdot (x-x0)} dx                     (1)


       for n=0,1,...,nn, where T_n denotes the n-th Chebychev 
       polynomial, and a,z,nn,x0 are user-supplied parameters. The
       integrals are evaluated to full double precision or a little
       better.

166.
The file expolege_prod.f contains the subroutines expolege_prod 
evaluating integrals of the form

       \int_a^b P_n(x) \cdot e^{z \cdot (x-x0)} dx                     (1)

       for n=0,1,...,nn, where P_n denotes the n-th Legendre
       polynomial, and a,z,nn,x0 are user-supplied parameters. The
       integrals are evaluated to full double precision or a little
       better.

167.
The file jfuns_asymp.f contains a code documenting the asymptotic 
         expansion for Bessel functions (J_{nu}) for large nu, 
         valid to 1/(nu)^5. The said expansion is obtained by 
         elementary analysis from the formula on page 227 of 
         Watson's "Bessel Functions".

168.
The file charpol_roots.f contains the subroutine charpol_roots 
constructing the eigenvalues of a user-supplied symmetric positive 
definite tridiagonal matrix located between the user-supplied 
limits x1,x2. If some of the eigenvalues are so close to each 
other that separation is difficult, the subroutine returns
them in a lump.

169.
The file bidiag_svd.f contains the subroutines bidiag_svd,
charpol_roots, constructing the singular values and the LEFT 
singular vectors of the user-supplied LOWER bidiagonal matrix A. 
The singular values are determined via a version of the Sturm 
method; subsequently a version of the inverse power method 
is used to find the corresponding singular vectors.

charpol_roots - constructs the eigenvalues of a user-supplied 
symmetric tridiagonal matrix located between the limits xs,x2 
supplied by the user. If some of the eigenvalues are so close
to each other that separation is difficult, the subroutine 
returns them in a lump.

Please note that, inter alia, this file sepercedes the file 
charpol_roots.f

170.
The file fvdpivot,f contains the subroutine fvdpivot, 
constructing a distorted singular value decomposition (to be 
referred to as the FVD) of the uer-specified matrix A. THIS
IS AN EARLY-STAGE CODE. I DO NOT EXPECT IT TO DEVELOPE ANY 
FURTHER. 4.20.03.

171. 
The file rskel_twomatr.f contains two user-callable subroutines: 
        rskel_twomatr and rskel_twomatr_basic. Following is 
        a brief description of these subroutines.


   rskel_twomatr - starts with constructing a two-sided 
        skeletonization of the matrix 

                c(n,k)=a(n,m)*b(m,k).                              (1)

        In other words, it finds the rank (to precision eps)
        ncols of the matrix a*b and a subset of ncols coordinates
        in R^m such that the matrix

        a      * b         = c                                     (2)
         |skel    |skel       |skel 

        has the smallest singular value that is reasonably close
        to the smallest singular value of c.

        Then, this subroutine constructs three mappings
        fillin(ncols,ncols), eval(m,ncols), expand(ncols,m)
        such that

        a(n,ncols) * fillin (ncols,ncols) * b(ncols,k)  = a*b,    (3)
         |skel                               |skel       

        a(n,m) * eval(m,ncols) * b(ncols,k)  = a*b,               (4)
                                  |skel       

        a(n,m) * expand(ncols,m) * b(m,k)    = a*b,               (5)
         |skel       

        all of the above equalities to be understood to precision 
        eps.
        
   rskel_twomatr_basic - constructs a two-sided skeletonization
        of the matrix 

                c(n,k)=a(n,m)*b(m,k).                              (1)

        In other words, it finds the rank (to precision eps)
        ncols of the matrix a*b and a subset of ncols coordinates
        in R^m such that the matrix

        a      * b         = c                                     (2)
         |skel    |skel       |skel 

        has the smallest singular value that is reasonably close
        to the smallest singular value of c.

        PLEASE NOTE THAT BOTH MATRICES A,B ARE UTTERLY DESTROYED 
        BY THIS SUBROUTINE!!!!

172.
The file cubeadap.f contains the subroutine cubeadap using the 
adaptive tensor product gaussian quadrature to integrate the
user-supplied function on a cube in R^3 (also user-supplied).

173.
The file ccubeadap.f contains the subroutine cubeadap using the 
adaptive tensor product gaussian quadrature to integrate the
user-supplied COMPLEX function on a cube in R^3 (also user-supplied).

174.
The file constcoef_crea.f contains 6 user-callable subroutines (or,
rather, entries): constcoef_crea, constcoef_eval, 
constcoef_convert_init, constcoef_convert_eval,
constcoef_convert_evalc, constcoef_discr0. Following 
is a brief description of these subroutines.


   constcoef_crea - constructs the various types of data 
        to be used in the (numerical) evaluation of integrals 
        of the form
       
        -1/(2*rlam)* \int_a^b e^{-rlam*|x-t|} * sigma(t) dt,          (1)

        with rlam a complex number, [a,b] a user specified 
        interval, and sigma a user supplied function to be 
        tabulated at the nodes xs provided by this subroutine.

   constcoef_eval - uses various types of data constructed
        via a preceding call to the entry constcoef_crea 
        to evaluate the expression (1) with x=x(1), x=x(2), 
        ... x=x(n*k) constructed by the entry constcoef_crea 

   constcoef_convert_init - constructs various types of data
        to be used (by the entries constcoef_convert_eval, 
        cconstcoef_convert_eval) to interpolate functions 
        from one composite Legendre discretization of the 
        interval [a,b] to another. Please note that all 
        subintervals within each of the discretizations are 
        of the same length.

   constcoef_convert_eval - interpolates a real function from 
        one composite Gaussian discretization of the interval 
        [a,b] to another. It uses data generated via a 
        preceding call to the entry constcoef_convert_init.

   constcoef_convert_evalc - interpolates a complex function 
        from one composite Gaussian discretization of the 
        interval [a,b] to another. It uses data generated via 
        a preceding call to the entry constcoef_convert_init.

   constcoef_discr0 - constructs a composite Gaussian discre-
        tization of the interval [-1,1]


175.
The file abskel.f contains two user-callable subroutines: abskel 
and abskel_basic. Following is a brief description of these two 
subroutines.

   abskel - skeletonizes the columns of the input matrix a IN 
        TERMS OF THE COLUMNS OF THE INPUT MATRIX B; then, it uses 
        the subroutine abskel_basic to approximate all columns of 
        b with the skeleton columns of a.

   abskel_basic - skeletonizes the columns of the input matrix a 
        IN TERMS OF THE COLUMNS OF THE INPUT MATRIX B.

176.
The file cskel_twomatr.f contains two user-callable subroutines: 
cskel_twomatr and cskel_twomatr_basic. Following is a brief 
description of these subroutines.


   cskel_twomatr - starts with constructing a two-sided 
        skeletonization of the complex matrix 

                c(n,k)=a(n,m)*b(m,k).                              (1)

        In other words, it finds the rank (to precision eps)
        ncols of the matrix a*b and a subset of ncols coordinates
        in R^m such that the matrix

        a      * b         = c                                     (2)
         |skel    |skel       |skel 

        has the smallest singular value that is reasonably close
        to the smallest singular value of c.

        Then, this subroutine constructs three mappings
        fillin(ncols,ncols), eval(m,ncols), expand(ncols,m)
        such that

        a(n,ncols) * fillin (ncols,ncols) * b(ncols,k)  = a*b,    (3)
         |skel                               |skel       

        a(n,m) * eval(m,ncols) * b(ncols,k)  = a*b,               (4)
                                  |skel       

        a(n,m) * expand(ncols,m) * b(m,k)    = a*b,               (5)
         |skel       

        all of the above equalities to be understood to precision 
        eps.
        
   cskel_twomatr_basic - constructs a two-sided skeletonization
        of the complex matrix 

                c(n,k)=a(n,m)*b(m,k).                              (1)

        In other words, it finds the rank (to precision eps)
        ncols of the matrix a*b and a subset of ncols coordinates
        in R^m such that the matrix

        a      * b         = c                                     (2)
         |skel    |skel       |skel 

        has the smallest singular value that is reasonably close
        to the smallest singular value of c.

        PLEASE NOTE THAT BOTH MATRICES A,B ARE UTTERLY DESTROYED 
        BY THIS SUBROUTINE!!!!

177.
The file savem.f contains the subroutine savem that reads the 
FORTRAN file "name1"; finds subroutines and functions in it 
that have not been declared "save", declares them things "save", 
and stores the obtained (much improved) code in the file "name2".


178.
The file cabskel.f contains two user-callable subroutines: cabskel 
and cabskel_basic. Following is a brief description of these two 
subroutines.

   cabskel - skeletonizes the columns of the input complex matrix 
        a IN TERMS OF THE COLUMNS OF THE complex INPUT MATRIX B; 
        then, it uses the subroutine abskel_basic to approximate 
        all columns of b with the skeleton columns of a.

   abskel_basic - skeletonizes the columns of the complex input 
        matrix a IN TERMS OF THE COLUMNS OF THE complex INPUT 
        MATRIX B.

179.
The file cleast_trans.f contains the subroutine cleast_trans 
solving in the least squares sense a collection of systems of 
linear equations

        A * alpha_j (T) \sim f_j,                              (1)

with j=1,2,...nalphas. 

In (1) above, alpha_j is the j-th diagonal m \times m - matrix, 
rhs_j is the j-th right hand side, and T is the vector to be 
found (via a least squares procedure).
     
Thus, this subroutine solves in the least squares sense a set of 
linesr systems (1), in which the systems differ by the right-hand 
sides rhs_j and the diagonal matrices alpha_j, while the matrix 
A is the same.


180.
The file prol_fir_bld.f contains the subroutine prol_fir_bld,
constructing a prolate function-based FIR filter approximating 
the user-supplied one; the nodes at which the signal is tabulated 
are equispaced, and the desired frequency response of the filter 
is supplied by the user via a subroutine evaluating the said filter.


181.
The file onedfmm_init.f is file contains five user-callable
subroutines: onedfmm_init, onedfmm_eval, onedfmm_eval2, 
onedfmm_evalc, onedfmm_eval2c. An adventurous user might 
also consider callable the entries onedfmm_eval_oneside, 
onedfmm_eval_oneside2, onedfmm_eval_onesidec, onedfmm_eval_oneside2c 
of the subroutine onedfmm_init. Following is a brief description of 
the first five of these subroutines.

   onedfmm_init - prepares the various data for the rapid application 
        of the one-sided FMM on the line, i.e. of the sums of the form

           \sum_{j=1}^n_i  sigma(j)/(ys(i)-xs(j)),               (1)

        for all j .leq. n_i. The actual evaluation of the sums
        (1) is performed by the entries onedfmm_eval_oneside,
        onedfmm_eval_oneside2, onedfmm_eval_onesidec,
        onedfmm_eval_oneside2c of this subroutine (see below).

   onedfmm_eval - applies the FMM operator on the line to the 
        user-provided real vector. In other words, it evaluates 
        the sum (1) for all pairs (xs(j),ys(i)) such that 
        
            xs(j) \neq ys(i).                                    (2)

   onedfmm_eval2 - applies the FMM operator on the line to TWO 
        user-provided real vectors. In other words,
        it evaluates the sums

           \sum_{j=1}^n sigma(j)/(ys(i)-xs(j)),                  (3)
 

           \sum_{j=1}^n sigma2(j)/(ys(i)-xs(j)),                 (4)
 
        subject to the condition (2)
   
   onedfmm_evalc - the same as onedfmm_eval, but with complex sigma
   onedfmm_eval2c - the same as onedfmm_evalc, but with complex 
        sigma, sigma2


182.
The file pnmproj_init.f contains three user-callable subroutines:
pnmproj_init, pnmproj_eval, pnmproj_evalc. Following is a brief 
description of these three subroutines.
        

   pnmproj_init - initializes the application of the spherical
         filter to user-specified functions; the latter
        procedure is performed by the subroutine pnmproj_eval
        (see). Specifically, given the function 

          sigma(x)=\sum{k=m}^N \alpha_k * P_k^m (x)            (1)

       tabulated at the nodes xs(1),xs(2),...,xs(nxs),
       the subroutine pnmproj_eval produces the function

          sigma(x)=\sum{k=m}^{nlege} \alpha_k * P_k^m (x),     (2)

       tabulated at the points ys(1),ys(2),...,ys(nys).

   pnmproj_eval - applies the sphrical filter to the 
        user-specified function sigma. Specifically, given 
        the function 

          sigma(x)=\sum{k=m}^N \alpha_k * P_k^m (x)            (1)

       tabulated at the nodes xs(1),xs(2),...,xs(nxs),
       the subroutine pnmproj_eval produces the function

          sigma(x)=\sum{k=m}^{nlege} \alpha_k * P_k^m (x),     (2)

       tabulated at the points ys(1),ys(2),...,ys(nys).

       PLEASE NOTE THAT THIS SUBROUTINE USES THE ARRAY W THAT
       (hopefully) has been created by a preceding call to
       the subroutine pnmproj_init (see). This subroutine
       has no uses known to man as a stand-alone device.       

   pnmproj_evalc - a complex *16 version of the subroutine 
       pnmproj_eval


183.
The file pnmlagr_init.f contains three user-callable subroutines: 
pnmlagr_init, pnmlagr_eval_fast, and pnmlagr_eval_fastc. 
Following is a brief description of these subroutines.

   pnmlagr_init - prepares data for the use of the subroutines 
        pnmlagr_eval, pnmlagr_evalc (see), performing fast 
        Lagrange interpolation. 

   pnmlagr_eval_fast - performs "fast" Lagrange interpolation
        of a real-valued function fs from the points xs to 
        the points ys. Please note that this subroutine uses 
        data prepared by a preceding call to the subroutine 
        pnmlagr_init. 

   pnmlagr_eval_fastc - performs "fast" Lagrange interpolation
        of a real-valued function fs from the points xs to 
        the points ys. Please note that this subroutine uses 
        data prepared by a preceding call to the subroutine 
        pnmlagr_init. 

184.
The file sphe_interp_init.f contains four user-callable subroutines: 
sphecre, spheint, sphe_interp_init, sphe_interp_eval. Following is
       a brief description of these subroutines.

   sphecre - discretizes the two-dimensional sphere, generating the
       standard grid, equispaced in phi and Gaussian in z. In 
       addition, coefficients of the associated quadrature formula
       are created. 

   spheint - integrates the user-specified function on the sphere. 
       It is assumed that the function is tabulated on a standard
       grid, constructed by the subroutine sphecre (see). In addition, 
       coefficients of the quadrature formula are used, also created 
       by the subroutine sphecre.
   sphe_interp_init - initializes the "fast" interpolation scheme on 
       the sphere from a courser "standard" grid to a finer one. 
       The actual interpolation is performed by the subroutine 
       sphe_interp_eval 

   sphe_interp_eval - applies the "fast" interpolation scheme on the 
       sphere from a courser standard grid to a finer one. It uses 
       the array w constructed via a preceding call to the subroutine 
       sphe_interp_init (see)

185.
The file sphe_filter_init.f contains two user-callable subroutines:
sphe_filter_init, sphe_filter_eval. Following is a brief description 
of these subroutines.

 
   sphe_filter_init - initializes the "fast" filtering scheme on 
       the sphere from a finer standard grid to a coarser one. 
       The actual filtering is performed by the subroutine 
       sphe_filter_eval 

   sphe_filter_eval - applies the "fast" filtering scheme on the 
       sphere from a finer standard grid to a coarser one. It uses 
       the array w constructed via a preceding call to the subroutine 
       sphe_filter_init (see)

186.  
The file prol_fir_dsgn.f contains the subroutines prol_fir_bld,
prol_fir_dsgn, constructing a prolate function-based FIR filter
approximating the user-supplied one; the nodes at which the signal is
tabulated are equispaced, and the desired frequency response of the
filter is supplied by the user via a subroutine evaluating the said
filter. PLEASE NOTE THAT THIS FILE SUPERSEDES THE FILE PROL_FIR_BLD.F.

187.
The file symm_grn12ders.f contains the subroutine symm_grn12ders
This subroutine evaluates the 2-jet of the Green's function for the 
cylindrically symmetrized Helmholtz equation; the accuracy produced
by this subroutine is about 10 digits (detailed tests have not been
performed as of 2.8.04). The Green's function is given by the formula

        F(rk,r0,r,z0,z)=
        \int_{0}^{2*\pi} e^{i * rk * rr}/rr *cos(m*t) dt,             (1)

        with 

        rr=sqrt((z-z0)^2+(r*cos(t)-r0)^2+(r*sin(t))^2 )               (2)


        Please note that the subroutine is designed to perform for 
        complex values of rk, but only when the argument of rk is 
        between 0 and pi/6. Also, note that 

        b=sqrt(r0*r)*abs(rk) 

        should not be greater than 50, at least if z and z0 can
        be close to each other

188.
The file appr_build.f contains the subroutines appr_build, appr_eval, 
appr_eval_c for the approximation of user-supplied functions by 
polynomials in their pristine polynomial form. Following is a brief 
description of these subroutines.

   appr_build - constructs the polynomial approximation of order n to 
        the user-supplied real- or complex-valued function; the 
        function is assumed to be supplied at the n Chebychev nodes ts 
        on the interval [-1,1]. The approximation is constructed in the 
        most straightforward way imaginable: by solving the linear 
        system

        \sigma_{j=1}^n coefs(j)*ts(i)^{j-1}=fs(i),                 (1)

        with i=1,2,...,n, and fs the values of the user-supplied
        function at the Gaussian nodes ts(1),ts(2),...,ts(n).

   appr_eval - evaluates the polynomial of degree n-1 with real 
        coefficients coefs at the real point t, returning the value 
        fout

   appr_eval_c - evaluates the polynomial of degree n-1 with complex
        coefficients coefs at the real point t, returning the value 
        fout

189.
The file d1mstrcr.f contains three user-callable subroutines:
d1mstrcr, d1m_onebox, d1m_onelist. Following is a brief 
description of these three subroutines.

   d1mstrcr - constructs the logical structure for the fully 
        adaptive FMM in one dimension and stores it in the
        array w. After that, the user can obtain the information 
        about various boxes and lists in it by calling the 
        subroutines d1m_onebox, d1m_onelist (see).

   d1m_onebox - the standard data for the user-specified box number 
        ibox (see the subroutine d1m_onebox for a complete description)

   d1m_onelist - returns to the user a single user-specified list of 
        the type itype (the possible types are 1,11,2,3,4) for the 
        user-specified box number ibox.

190.
The file jacskel.f contains the subroutine jacskel skeletonizing the 
ROWS of the user-supplied matrix a. Please note that this is a 
Jacobi-based code, as opposed to a Gram-Schmidt-based one. Due to the 
author's lack of experience in this environment, omissions are 
possible. Use at your own risk! - 03.06.04

191.
The file symm_grn12ders_new.f contains the subroutine symm_grn12ders.
This subroutine evaluates the 2-jet of the Green's function for the 
cylindrically symmetrized Helmholtz equation; the accuracy produced
by this subroutine is about 12 digits (detailed tests have not been
performed as of 2.8.04). The Green's function is given by the formula

        F(rk,r0,r,z0,z)=
        \int_{0}^{2*\pi} e^{i * rk * rr}/rr *cos(m*t) dt,             (1)

        with 

        rr=sqrt((z-z0)^2+(r*cos(t)-r0)^2+(r*sin(t))^2 )               (2)


        Please note that the subroutine is designed to perform for 
        complex values of rk, but only when the argument of rk is 
        between 0 and pi/6. Also, note that 

        b=sqrt(r0*r)*abs(rk) 

        should not be greater than 50, at least if z and z0 can
        be close to each other

PLEASE NOTE THAT THIS IS A HIGHER ACCURACY VERSION OF THE 
SUBROUTINE GRN12DERS FOUND IN THE FILE grn12ders.f. PLEASE
NOTE THAT THIS VERSION IS ALSO SLIGHTLY SLOWER.


192.
Thi file rnestex_cheb.f contains the numerical machinery for the 
construction and evaluation of nested Chebychev expansions of real
functions of one real variable. There are 3 user-callable subroutines 
in the file: rnestex_cheb, rneste_cheb, rneste2_cheb. 

193.
The file grmskel skeletonizes the columns of the user-supplied
matrix via the standard double Gram-Schmidt procedure. At the 
present time (10.03.04) this is the fastest general-purpose 
skeletonizer we have; it is also mercifully short.
 
