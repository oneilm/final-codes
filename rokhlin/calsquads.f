        implicit real *8 (a-h,o-z)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER k '
ccc        READ *,k
	k=20
        CALL PRINF('k=*',k,1 )
C 
        PRINT *, 'ENTER nsings '
ccc        READ *,nsings
	nsings=11
        CALL PRINF('nsings=*',nsings,1 )
  
  
cccc        call testex(k,nsings)
  
        call testexs(k,nsings)
  
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine testexs(k,nsings)
        implicit real *8 (a-h,o-z)
        save
        dimension par1(100),ab(20 000),x(10000),whts(10000),
     1      rlams(10 000),rints(10 000),rnorms(10 000),work(1000000),
     2      xsout(1000),uu(100 000),vv(100 000)
c 
        complex *16 w(4 000 000),fs2(1000 000),prods2(1000 000),
     1      ftest(10000),der,coefs(10000),errs(10000),ftest2(10000),
     2      whtsout(1000),forminte(1000),wsout2(10000)
  
        complex *16 fi,zint,ztru,zint2,fi2
  
        external fun1,fun2,fun3
c 
c      construct all parameters
c 
        call fun1ini(nsings)
  
        a=1.0d-6
        a=1.0d-2
        a=0
ccc        a=1.0d-10
        b=0.5
cccc        b=1-a
  
        eps=1.0d-13
        lab=10 000
c 
  
  
        nfuns=nsings
  
        lw=4 000 000
        ifab=1
        lrlams=10000
        igraph=0
  
        call calsvbld(ier,a,b,k,fun1,par1,par2,
     1      nfuns,eps,ab,lab,ifab,nn,x,whts,n,rlams,lrlams,
     2      w,lw,ncols,rints,lused,igraph)
  
  
        call prinf('after calsvbld, nn=*',nn,1)
        call prin2('after calsvbld, rlams=*',rlams,ncols*2)
        call prin2('after calsvbld, rints=*',rints,ncols*2)
c 
c       check the obtained singular functions
c 
ccc        call allsvchk(ier,w,ncols,eps,x,whts,n,rlams,
ccc     1      fs2,prods2,ab,nn,k)
c 
c 
c       find a collection of points on the interval [a,b]
c       such that the matrix of values of singular vectors at these
c       points is well-conditioned
c 
        ifuv=1
        ncols2=ncols-5
        ncols2=ncols
  
        call calspnts(ier,ncols2,x,whts,ab,k,nn,w,
     1      xsout,whtsout,rints,uu,vv,ifuv,rcond,w(4*n*ncols+10) )
c 
        call prin2('after calspnts, rcond=*',rcond,1)
        call prin2('after calspnts, xsout=*',xsout,ncols2)
        call prin2('after calspnts, vv=*',vv,2*ncols2**2)
  
  
  
cccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccc
  
c 
  
	do 1600 iii=1,6
c 
	lwork=1000000
	eps=1.0d-8
	mm=10
	call calsquads(ier,a,b,fun2,par2,mm,eps,
     1	   ncols,x,whts,ab,k,nn,w,xsout,whtsout,wsout2,
     2	   uu,work,lwork)
ccc     2	   uu,rnorms,rints,work)
c 
        call prin2('after calsquads, xsout=*',xsout,ncols2)
        call prin2('after calsquads, wsout=*',whtsout,2*ncols2)
        call prin2('after calsquads, wsout2=*',wsout2,2*ncols2)
c 
c	test the quadrature:
c 
ccc	iii=6
c 
	zint=0
	zint2=0
	do 1300 i=1,ncols
	   t=xsout(i)
           call fun1(t,iii,par1,par2,fi)
	   zint=zint+wsout2(i)*fi
c 
	   call fun2(t,par1,fi2)
	   zint2=zint2+whtsout(i)*fi*fi2
 1300 continue
  
	call prin2('after quads, zint=*',zint,2)
	call prin2('after quads, zint2=*',zint2,2)
c 
c	do it use cadapgau
c 
	m=20
	 
	call cadapgau(ier,a,b,fun3,iii,par2,m,eps,
     1      ztru,maxrec,numint)
  
  
	call prin2('after cadapgau, ztru=*',ztru,2)
	call prin2('the difference is=*',abs(ztru-zint),1)
	call prin2('the 2nd difference is=*',abs(ztru-zint2),1)
  
 1600 continue
  
        return
        end
c 
c 
c 
c 
c 
c 
        function fun3(t,iii,par2)
        implicit real *8 (a-h,o-z)
        save
        dimension par1(1),par2(1)
c 
        complex *16 f,ima,fun3,f2
c 
	external fun1,fun2
c 
        data ima/(0.0d0,1.0d0)/
c 
ccc	iii=6
        call fun1(t,iii,par1,par2,f)
        call fun2(t,par1,f2)
c 
ccc	fun3=f
	fun3=f*f2
  
        return
	end
c 
c 
c 
c 
c 
        subroutine fun2(t,par1,f)
        implicit real *8 (a-h,o-z)
        save
        dimension par1(1),par2(1)
c 
        complex *16 f,ima
c 
        data ima/(0.0d0,1.0d0)/
c 
        f=exp(ima*t)
  
        return
	end
c 
c 
c 
c 
c 
c 
        subroutine apprtest(xsout,ncols,v,itest,ntest,a,b,
     1      ab,nn,k,coefs)
        implicit real *8 (a-h,o-z)
        save
        dimension xsout(1),xstest(10000),v(1),ab(2,1),w(1)
c 
        complex *16 fsout(1000),fstest(1000),forminte(1000),cd,
     1      errs(1000)
c 
c       evaluate the test function at the support nodes
c 
        do 1100 i=1,ncols
c 
        call fun1(xsout(i),itest,par1,par2,fsout(i))
 1100 continue
c 
        call prin2('in apprtest, xsout=*',xsout,ncols)
        call prin2('in apprtest, fsout=*',fsout,ncols*2)
c 
c       create the test nodes and the test function at them
c 
        h=(b-a)/ntest
        do 1200 i=1,ntest
c 
        xstest(i)=a+(i-1)*h+h/2
        call fun1(xstest(i),itest,par1,par2,fstest(i))
 1200 continue
c 
        call prin2('in apprtest, xstest=*',xstest,ntest)
        call prin2('in apprtest, fstest=*',fstest,ntest*2)
c 
c       for each of the test points, interpolate the function
c       from the support nodes to the test point
c 
        do 1600 i=1,ntest
c 
        call cnesinte(ab,nn,k,coefs,xsout,v,ncols,xstest(i),
     1      forminte)
c 
        call prinf('i=*',i,1)
        call prin2('and forminte=*',forminte,ncols*2)
c 
        cd=0
        do 1400 j=1,ncols
c 
        cd=cd+forminte(j)*fsout(j)
 1400 continue
c 
        call prin2('and cd=*',cd,2)
        call prin2('while fstest(i)=*',fstest(i),2)
        call prin2('and fstest(i)-cd=*',fstest(i)-cd,2)
c 
        errs(i)=fstest(i)-cd
 1600 continue
  
        call prin2('and errs=*',errs,ntest*2)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine fun1(x,i,par1,par2,f)
        implicit real *8 (a-h,o-z)
        save
        dimension par1(1),par2(1),sings(10000)
c 
        complex *16 f,ima
c 
        data ima/(0.0d0,1.0d0)/
c 
        f=x**sings(i)+ima*(1-x)**sings(i)
ccccc        f=x**sings(i)
cccc        f=x**(sings(i)+1)+ima*(1-x)**sings(i)
  
        return
c 
c 
c 
        entry fun1ini(nsings7)
c 
        nsings=nsings7
        done=1
        h=10.3*done/nsings
        do 2200 j=1,nsings
        sings(j)=j*h
 2200 continue
c 
        call prin2('in fun1ini, sings=*',sings,nsings)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine allsvchk(ier,coefs,ncols,eps,x,whts,n,rlams,
     1      fs,prods,ab,nn,k)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),whts(1),rlams(1),ab(2)
        complex *16 fs(n,ncols),prods(ncols,ncols),cd,f,cds(1000)
c 
c       evaluate all singular functions at all nodes in array x
c 
        do 1400 i=1,ncols
c 
        do 1200 j=1,n
c 
        call cnestev2(ier,coefs,ab,nn,k,x(j),i,fs(j,i))
 1200 continue
 1400 continue
  
cccc         call prin2('in allsvchk, fs=*',fs,n*ncols*2)
cccc         call prin2('in allsvchk, whts=*',whts,n)
  
cccc           stop
  
  
c 
c       evaluate all inner products of all functions
c 
        do 2400 i=1,ncols
        do 2200 j=1,ncols
c 
        cd=0
        do 1800 ll=1,n
c 
        cd=cd+fs(ll,i)*conjg(fs(ll,j))*whts(ll)
 1800 continue
  
cccc         call prinf('in allsvchk, i=*',i,1)
cccc         call prinf('in allsvchk, j=*',j,1)
cccc         call prin2('in allsvchk, cd=*',cd,2)
c 
         prods(j,i)=cd
 2200 continue
 2400 continue
c 
        call prin2('prods in allsvchk are*',prods,ncols**2*2)
c 
c       project one of the original functions on the singular
c       vectors
c 
        ipar1=3
c 
        do 3400 i=1,ncols
c 
        cd=0
        do 3200 j=1,n
c 
        call fun1(x(j),ipar1,par1,par2,f)
c 
        cd=cd+f*conjg(fs(j,i))*whts(j)
 3200 continue
  
ccc        call prinf('i=*',i,1)
ccc        call prin2('and cd=*',cd,2)
c 
        cds(i)=cd
 3400 continue
  
        call prin2('and cds=*',cds,ncols*2)
  
  
        return
        end
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c   In between these two sets of 4 "ccc..." lines is the subroutine
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c 
  
  
	subroutine calsquads(ier,a,b,fun,par2,mm,eps,ncols,
     1	   xs,whts,ab,k,nn,w,xsout,wsout,wsout2,uu,work,lwork)
	implicit real *8 (a-h,o-z)
c 
c 
c Purpose:
c 
c	This subroutine produces quadratures for the class of
c	functions (that has been compressed by "calsvcmb" previously)
c	multiplied by a user-specified weight function "fun".
c	The weight function can be accessed via the calling sequence:
c 
c	   call fun(t,par1,f)
c 
c	So this code basically generates:
c 
c   	     \int_a^b fun(t) * \phi_j (t) dt =
c 
c	\sum_{k=1}^{ncols} wsout(k) * fun(xsout(k)) * \phi_j (xsout(k))
c	
c Inputs:
c 
c   a,b --- the ends of integration interval [a,b]
c   mm --- the order of Gaussian quads used on subintervals
c   eps --- the accuracy of the quadrature
c   ncols --- the number of singular functions from alsvcmb
c   xs,whts,ab,k,nn,w --- the info produced by "alsvcmb" during
c		compressing process.
c 
c Outputs:
c 
c   ier --- error return codes
c		=0 OK
c		=8 not enough memory in workspace
c   xsout --- the nodes for the qudratures
c   wsout --- the weights of the quadratures
c   wsout2 --- another set of weights, if this set of weights is used
c		instead of the "wsout" above, one need not compute
c		the weight function at the nodes to calculate the
c		the integral, i.e. this set of weights already include
c		the weight function "fun"
c 
c Workspace:
c 
c   uu --- used to store the inverse of the matrix that produces
c		the quad weights, must be (2*ncols**2) long,
c		AND VERY IMPORTANTLY, is should not be overwritten
c		between calls to this subroutine, if one wants ti
c		produce more than ONE set of quadratures
c 
c   work --- the work space that needs to be
c		(nn*k+ncols*2+ncols*nn*k*2+20) long.
c 
c 
c 
c 
        dimension xs(1),whts(1),ab(1),w(1),xsout(1),work(1)
        complex *16 uu(1),wsout(1),wsout2(1),fw
c 
	ier=0
c 
	irnorms=1
	lrnorms=nn*k+10
c 
	irints=irnorms+lrnorms
	lrints=2*ncols+10
c 
	iwork=irints+lrints
c 
	if (lwork.lt.iwork+ncols*nn*k*2) then
	   ier=8
	   return
	endif
c 
	call calsquads0(jer,a,b,fun,par2,mm,eps,
     1	   ncols,xs,whts,ab,k,nn,w,xsout,wsout,wsout2,
     2	   uu,work(irnorms),work(irints),work(iwork))
c 
        return
        end
c 
c 
c 
c 
c 
c 
c 
c 
	subroutine calsquads0(ier,a,b,fun,par2,mm,eps,
     1	   ncols,xs,whts,ab,k,nn,w,xsout,wsout,wsout2,
     2	   uu,rnorms,rints,work)
	implicit real *8 (a-h,o-z)
c 
        dimension xs(1),whts(1),ab(1),w(1),xsout(1)
        complex *16 uu(1),wsout(1),wsout2(1),fw
c 
        dimension rnorms(1),work(1)
        complex *16 rints(1)
c 
	external calsrhs,fun
c 
	data ifinit/0/
c 
c ... Initialization:
c       find a collection of points on the interval [a,b]
c       such that the matrix of values of singular vectors at these
c       points is well-conditioned as well as the inverse of this
c	matrix
c 
	if (ifinit.eq.0) then
c 
           call calsmat(ncols,xs,whts,ab,k,nn,w,
     1         xsout,uu,rcond,work,rnorms)
  
ccc        call prin2('after calsmat, xsout=*',xsout,ncols)
ccc        call prin2('after calsmat, uu=*',uu,2*ncols**2)
ccc        call prin2('after calsmat, rcond=*',rcond,1)
  
           call calsrini(w,ab,nn,k,ncols)
	 
	   ifinit=1
c 
	endif
  
c 
c ... CALCULATE the integrals with the user-specified weight funs:
c 
        nfuns=ncols
  
        call cadapgaum(ier,a,b,calsrhs,nfuns,fun,par2,mm,eps,
     1      rints,maxrec,numint,work)
  
ccc        call prin2('after adapgaum, rints=*',rints,ncols*2)
  
c 
c ... Get the weights:
c 
        call calqrma5(uu,rints,wsout2,ncols)
  
ccc        call prin2('after calqrma5, wsout2=*',wsout2,ncols*2)
c 
	do 1400 i=1,ncols
	   t=xsout(i)
	   call fun(t,par2,fw)
c 
	   if (abs(fw).lt.1.0d-20) then
	      wsout(i)=0
	   else
	      wsout(i)=wsout2(i)/fw
	   endif
 1400 continue
ccc        call prin2('after calqrma5, wsout=*',wsout,ncols*2)
  
c 
	return
	end
c 
c 
c 
c 
c 
c 
c 
c 
c 
c 
c 
c 
	subroutine calsrhs(t,n,fun,par2,vals)
	implicit real *8 (a-h,o-z)
	save
c 
c	This function is used for generating the integrals of
c	the left singular functions (from calsvcmb) multiplied by
c	a user-specified function
c 
	dimension w7(1),ab7(1)
	dimension w(1000000),ab(200000)
c 
	complex *16 vals(1),fw
c 
	call cnestem(ier,w,ab,nn,k,ncols,t,vals)
	call fun(t,par2,fw)
c 
        do 1200 i=1,ncols
	   vals(i)=vals(i)*fw
 1200 continue
c 
	return
c 
c 
c 
c 
c 
c 
c 
	entry calsrini(w7,ab7,nn7,k7,ncols7)
c 
c	initialize the stuff that comes out of "calsvcmb"
c 
	ncols=ncols7
	k=k7
	nn=nn7
c 
	do 2020 i=1,ncols7*nn7*k7*4
	   w(i)=w7(i)
 2020 continue
c 
	do 2040 i=1,nn7*2
	   ab(i)=ab7(i)
 2040 continue
c 
        return
        end
c 
c 
c 
c 
c 
c 
c 
c 
c 
        subroutine calsmat(ncols,xs,whts,ab,k,nn,w,
     1      xsout,uu,rcond,uuadj,rnorms)
        implicit real *8 (a-h,o-z)
c 
c Purpose:
c 
c	This subroutine generates the nodes (good for interpolation
c	as well as quadrature for the functions that has been
c	compressed via alsvbld or alsvcmb) and the inverse
c	of the matrix that will produce quadrature weights.
c 
c 
c Inputs:
c 
c	standard inputs (see calsquad above)
c 
c Outputs:
c 
c    xsout --- the interpolation or quadrature nodes
c    uu --- the inverse matrix that generates qudrature weights
c 
c 
        save
        dimension ipivots(10 000),whts(1),xs(1),w(1),xsout(1),
     1      ab(1)
c 
        complex *16 uuadj(ncols,nn*k),rnorms(1),uu(ncols,ncols)
c 
	n=nn*k
c 
c       evaluate all of the singular vectors at the nodes of
c       the discretization xs of the interval of definition
c 
        do 1450 j=1,n
c 
        call cnestem(ier,w,ab,nn,k,ncols,xs(j),uuadj(1,j) )
c 
        d=sqrt(whts(j))
        do 1430 i=1,ncols
c 
        uuadj(i,j)=uuadj(i,j)*d
 1430 continue
 1450 continue
  
        eps2=1.0d-30
c 
        call calsppiv(uuadj,ncols,n,rnorms,eps2,ncols3,ipivots)
ccc	call prinf(' in calsquads, ncols3=*',ncols3,1)
c 
c       sort the array ipivots
c 
        call calsbubb(ipivots,ncols3)
c 
c       pick out of the array xs the points at which the matrix
c       of values of the singular vectors is well-conditioned
c 
        do 1500 i=1,ncols
c 
        ii=ipivots(i)
        xsout(i)=xs(ii)
 1500 continue
c 
c       construct the matrix of values of singular vectors at the
c       selected nodes
c 
c 
        do 1600 j=1,ncols
c 
        call cnestem(ier,w,ab,nn,k,ncols,xsout(j),uu(1,j) )
 1600 continue
c 
c       invert the matrix uu
c 
	call matrans2(uu,ncols)
c 
        call corthom(uu,ncols,uuadj,rcond)
c 
        return
        end
c 
c 
c 
c 
c 
c 
c 
c 
        subroutine matrans2(uu,n)
        implicit real *8 (a-h,o-z)
c 
        save
        complex *16 uu(n,n),ztmp
c 
        do 1040 i=1,n-1
        do 1020 j=i+1,n
           ztmp=uu(i,j)
           uu(i,j)=uu(j,i)
           uu(j,i)=ztmp
 1020 continue
 1040 continue
c 
        return
        end
