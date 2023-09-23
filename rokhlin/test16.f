        implicit real *8 (a-h,o-z)
        dimension par1(100),rlams(10 000),
     1      ftest(10 000),w(4000 000),vmatr(100 000)
c 
        dimension coefs(10 000),errs(10 000),s(1000),usout(100)
c 
        dimension xs(1000),ws(1000),u(100 000),v(100 000),
     1      coefs2(10000)
        external fun1,fun2
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER npols '
        READ *,npols
        CALL PRINF('npols=*',npols,1 )
c 
  
c 
c       retrieve the nodes and weights from the user-supplied
c       subroutine xsandws
c 
        call xsandws2(xs,ws,npts1)
  
  
  
        call prinf('npts1 as read from subroutine*',npts1,1)
        call prin2('xs as read from subroutine*',xs,npts1)
        call prin2('ws as read from subroutine*',ws,npts1)
  
  
cccc          stop
  
        eps=1.0d-14
        eps2=1.0d-13
c 
  
cccc        eps=1.0d-14
        eps2=1.0d-7
c 
cccc        eps=1.0d-32
cccc        eps2=1.0d-26
c 
c       construct the SVD of the user-specified functions,
c       and use it to construct the matrix vmatr
c 
        a=-1
        b=1
        k=24
        lw=4000 000
        igraph=0
  
        nfuns=npols*3
  
        nfuns=npols*2
        call fun1ini(npols)
cccc        call fun2ini(npols)
  
        iwrite=31
        call bigsub(ier,a,b,k,fun1,par1,par2,
     1       nfuns,eps,eps2,w,lw,lused,ncols,rlams,nn,
     2      vmatr,xs,ws,npts1,u,v,errs,s,iwrite)
  
  
        n=nn*k
  
        call prinf('after bigsub, ier=*',ier,1)
        call prin2('and after bigsub, errs=*',errs,ncols)
        call prin2('and after bigsub, s=*',s,ncols)
  
        n=npts1
c 
c       . . . test u
c 
cccc          stop
c 
c       test the matrix v
c 
  
c 
        i=3
        do 3200 ii=1,npts1
c 
        call nesteva(ier,w,xs(ii),i,ftest(ii),der)
c 
 3200 continue
  
        call prin2('ftest as obtained while testing v*',ftest,npts1)
        call prin2('while xs=*',xs,npts1)
        call prin2('while v=*',xs,npts1)
  
        call uvmatvec2(v,ftest,coefs,ncols,npts1)
  
        call prin2('and coefs after uvmatvec2*',coefs,ncols)
  
  
  
  
        x=xs(10)
  
        call ufuneval(fun1,nfuns,ncols,vmatr,x,usout)
  
        call prin2('after ufuneval, usout=*',usout,ncols)
  
  
  
ccccc            stop
  
  
        call ufunsinterp(fun1,vmatr,nfuns,ncols,v,npts1,
     1      x,coefs2)
  
        call prin2('after ufunsinterp, coefs2=*',coefs2,npts1)
  
        d=0
        do 3400 i=1,20
c 
        d=d+ftest(i)*coefs2(i)
 3400 continue
  
        call prin2('and d=*',d,1)
        call prin2('and d-ftest(10)=*',d-ftest(10),1)
  
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine xsandws2(xs,ws,npts)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),ws(1)
        dimension xs0(20),ws0(20)
  
c 
        data xs0/
     1  -.9999219386549200588D+00,-.9991275778562456269D+00,
     2  -.9958847972714327735D+00,-.9870507654600293429D+00,
     3  -.9682235882727754893D+00,-.9343306635306926492D+00,
     4  -.8804434022267430932D+00,-.8025841146730825185D+00,
     5  -.6983746898933523478D+00,-.5674661041648135342D+00,
     6  -.4117399983159040265D+00,-.2352952734665198287D+00,
     7  -.4424120722671363372D-01,0.1536756076188653947D+00,
     8  0.3495769666991975487D+00,0.5340515181579167453D+00,
     9  0.6978009400925435295D+00,0.8322860528975736109D+00,
     *  0.9303262529811095438D+00,0.9866229602021953277D+00/
c 
        data ws0/
     1  0.2562641557039653290D-03,0.1615681421487374685D-02,
     2  0.5396278774612809361D-02,0.1301582524338190343D-01,
     3  0.2549814314465490130D-01,0.4312293995802031334D-01,
     4  0.6532546890661728074D-01,0.9079662496606951401D-01,
     5  0.1176837195144482431D+00,0.1438195699911100382D+00,
     6  0.1669459788486692028D+00,0.1849189604743250196D+00,
     7  0.1958895493744871426D+00,0.1984542969496517426D+00,
     8  0.1917689777886093850D+00,0.1756195226107159696D+00,
     9  0.1504460385325518634D+00,0.1173187882395179351D+00,
     *  0.7787042972344134866D-01,0.3423694138192404703D-01/
c 
c       flip the nodes and weights
c 
        npts=20
        do 1200 i=1,npts
c 
        xs(i)=xs0(i)
        ws(i)=ws0(i)
 1200 continue
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
        dimension par1(1),par2(1)
c 
c 
        j=i
        if(i .gt. n) j=i-n
        if(i .gt. n*2) j=i-n*2
  
        call legepol(x,j-1,d,der)
  
cccc        d=x**(i-1)
  
        f=d
c 
cccc        if(i .gt. n) f=d*log(1+x)
cccc        if(i .gt. n*2) f=d*log(1-x)
  
        if(i .gt. n) f=d*log(1+x)*(1+x)
        if(i .gt. n*2) f=d*log(1-x)*(1-x)
  
cccc        call prin2('f=*',f,1)
c 
        return
c 
c 
c 
c 
        entry fun1ini(n7)
c 
        n=n7
        call prinf('in fun1ini, n=*',n,1)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine fun2(x,i,par1,par2,f)
        implicit real *8 (a-h,o-z)
        save
        dimension par1(1),par2(1)
c 
c 
        j=i
        if(i .gt. n) j=i-n
        if(i .gt. n*2) j=i-n*2
  
        call legepol(x,j-1,d,der)
  
cccc        d=x**(i-1)
  
        f=d
c 
cccc        if(i .gt. n) f=d*log(1+x)
cccc        if(i .gt. n*2) f=d*log(1-x)
  
        if(i .gt. n) f=d*log(1+x)*(1+x)
        if(i .gt. n*2) f=d*log(1-x)*(1-x)
  
cccc        call prin2('f=*',f,1)
c 
        return
c 
c 
c 
c 
        entry fun2ini(n7)
c 
        n=n7
        call prinf('in fun1ini, n=*',n,1)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine uvmatmul(a,n,m,b,k,c)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),b(m,k),c(n,k)
c 
        do 1600 i=1,n
        do 1400 j=1,k
        d=0
        do 1200 l=1,m
c 
        d=d+a(i,l)*b(l,j)
 1200 continue
c 
        c(i,j)=d
 1400 continue
 1600 continue
        return
        end
  
c 
c 
c 
c 
c 
        subroutine uvmatvec2(a,x,y,n,m)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),x(m),y(n)
c 
        do 1400 i=1,n
        d=0
        do 1200 j=1,m
c 
        d=d+a(i,j)*x(j)
 1200 continue
c 
        y(i)=d
 1400 continue
        return
        end
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code and the beginning of
c        the interpolation code proper
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine bigsub(ier,a,b,k,fun1,par1,par2,
     1      nfuns,eps,eps2,w,lw,lused,ncols,rlams,nn,
     2      vmatr,xs,ws,npts1,u,v,errs,s,iwrite)
        implicit real *8 (a-h,o-z)
        save
        dimension par1(1),par2(1),rlams(1),w(1),vmatr(1),
     1      xs(1),ws(1),u(1),v(1),s(1),errs(1)
c 
        external fun1
c 
c        This subroutine builds (primarily) the matrices vmatr,
c        u,v to be used for the design of interpolation formulae
c        in the "special function" environment. It uses the SVD
c        to construct the basis functions from the user-supplied
c        function evaluator fun1, and constructs the array vmatr
c        that can be used by the subroutine ufuneval (see) to
c        evaluate all of the constructed basis functions at
c        arbitrary points on the interval [a,b].
c 
c        Then, it constructs the matrix u of values of the basis
c        functions at the USER-SUPPLIED nodes xs, and the inverse
c        v of u. The pair (vmatr,v) can now be used to construct
c        interpolation formulae on the interval [a,b] for functions
c        defined by the user-supplied subroutine fun1. Finally. it
c        stores the matrices vmatr,u,v and the arrays xs, ws on the
c        FORTRAN unit iwrite (user-supplied) in a format easily
c        convertible into FORTRAN DATA statements.
c 
c        One the matrixs vmatr has been constructed, the basis
c        functions can (and should) be evaluated via the subroutine
c        ufuneval (see), with the aid of the user-supplied subroutine
c        fun1.
c 
c 
c                        Input parameters:
c 
c  a,b - the ends of the interval on which the functions live
c  k - the number of nodes in the Legendre interpolation to be used
c        by the subroutine allsvbld
c  fun1 - the subroutine evaluating the functions to be compressed;
c 
c      The calling sequence of fun1 must be:
c 
c      fun(x,i,par1,par2,f),
c 
c       where x is the input argument where the function is to be
c       evaluated, i the sequence number of function, f (output)
c       is the value of the i-th function at the point x, and par1,
c       par2 are input parameters, carrying whatever information fun1
c       might need.
c 
c  par1,par2 - whatever the user-supplied subroutine fun needs
c  nfuns - the number of functions to be compressed
c  eps - the accuracy to which the legendre expansion must represent
c      the user-supplied function for the subroutine to decide that
c      the representation is accurate
c  lw - the length of the user-provided array w. It has to be
c      large. If it is insufficient, ier os set to 64, and the
c      execution of the subroutine terminated
c  xs - the nodes at which the basis functions are to be evaluated
c      in order to construct the matrix v
c  ws - the quadrature weights corresponding to the nodes xs
c  npts1 - the number of elements in each of the arrays xs, ws
c  iwrite - the FORTRAN unit number on which the arrays xs, ws,
c      vmatr, u,v are to be stored.
c 
c                       Output parameters:
c 
c  ier -error return code;
c          ier=0 means normal return
c          ier=4 means that the amount of space in the array ab
c                 (given by the parameter lab) is insufficient. very
c                 often, the real cause of the problem is eps that
c                 is too small.
c          ier=64 means that the amount of space in the user-provided
c                 array w (given by the parameter lw) is insufficient.
c  w - the first n*ncols elements of array w contain the coefficients of
c        the nested Legendre expansions of the singular functions;
c 
c        IMPORTANT NOTE: the first n*ncols elements of array w should
c                        not be changed between the call to this
c                        subroutine and subsequent calls to the
c                        subroutine nesteva or its relatives
c  lused - the length of the piece of the array w that was actually
c          used by the subroutine
c  ncols - the number of QR-functions in the QR-decomposition of the
c      nfuns functions given by the subroutine fun (see above) whose
c      corresponding weights are greater than eps (see above). It is
c      also the number of meaningful columns of arrays fs, coefs on
c      exit from this subroutine.
c  rlams - the singular values of the matrix of functions given by the
c        function fun. Note that though only ncols elements of array rlams
c        are meaningful, nfuns elements have to be allocated by the
c        user.
c  nn - the number of subintervals in the nested discretization of
c      the interval [a,b].
c  vmatr - the nfuns \times ncols matrix to be used by the subroutine
c        ufuneval for the evaluation of the singular functions
c  u - the npts1 \times ncols - matrix of values at the nodes xs supplied
c        by the userof the singular functions (obtained via SVD)
c  v - the generalized inverse of u
c  errs - errs(i) is the maximum of errors of the evaluation of the
c        i-th singular functions, evluated at 100 equispaced points
c        on the interval [a,b]
c  s - the singular values of the matrix uu given by the formula
c 
c                    uu(j,i)=u(j,i)*sqrt(ws(j))                       (1)
c 
c        the condition number of uu is a good estimate of the quality
c        of interpolation
c 
c                            Work array:
c 
c  w - must be large
c 
c       . . . allocate memory
c 
        ier=0
c 
        irints=1
        lrints=nfuns+10
c 
        iw=irints+lrints
        lenw=lw-iw-5
c 
        lrlams=nfuns+100
c 
        call allsvbld(jer,a,b,k,fun1,par1,par2,
     1      nfuns,eps,rlams,lrlams,w(iw),lenw,ncols,nn,
     2      w(irints),lused,keep,igraph)
c 
        n=nn*k
        call prinf('after allsvbld, ier=*',ier,1)
        call prinf('after allsvbld, lused=*',lused,1)
        call prinf('after allsvbld, ncols=*',ncols,1)
        call prin2('after allsvbld, rlams=*',rlams,ncols)
c 
        if(jer .ne. 0) then
            ier=16
            return
        endif
c 
c       calculate the appropriate ncols
c 
        do 1200 i=1,ncols
c 
        ncols2=i
        if(rlams(i)/rlams(1) .lt. eps2) goto 1300
 1200 continue
 1300 continue
c 
        ncols=ncols2
c 
c        compress the array w
c 
        call mvecmove(w(iw),w,keep+10)
        iw=1
c 
        ifmatr=keep+10
        lfmatr=nfuns*n+10
c 
        iumatr=ifmatr+lfmatr
        lumatr=n*ncols+10
c 
        ixlarge=iumatr+lumatr
        lxlarge=n+10
c 
        iwlarge=ixlarge+lxlarge
        lwlarge=n+10
c 
        iab=iwlarge+lwlarge
        lab=2*nn+10
c 
        iwork=iab+lab
        lenwork=lw-iwork-10
c 
        call bigsub0(ier,a,b,n,nfuns,w(iw),ncols,rlams,
     1      w(ifmatr),w(iumatr),vmatr,xs,ws,npts1,u,v,
     2      w(iwork),lenwork,errs,w(ixlarge),w(iwlarge),
     3      w(iab),fun1,eps)
c 
c        calculate the singular values of the obtained matrix u
c 
        iuu=ifmatr
        luu=npts1*ncols+10
c 
        iuuu=iuu+luu
        luuu=npts1*ncols+10
c 
        ivvv=iuuu+luuu
        lvvv=npts1*ncols+10
c 
        iwww=ivvv+lvvv
c 
        call uvmattst(ws,npts1,ncols,u,w(iuu),w(iuuu),
     1        w(ivvv),w(iwww),s)
c 
c        store the arrays xs, ws
c 
 2200 format('        following is the ',i3,'-point array xs',
     1    '  as supplied by the user')
c 
         write(iwrite,2200) npts1
c 
        call arrstore(iwrite,xs,npts1,'      array xs*')
  
c 
 2400 format('        following is the ',i3,'-point array ws',
     1    '  as supplied by the user')
c 
         write(iwrite,2200) npts1
c 
        call arrstore(iwrite,ws,npts1,'      array ws*')
c 
c        store the obtained matrix vmatr
c 
 2600 format('        following is the matrix vmatr, with ',
     1    'nfuns=',i3,', and ncols=',i3,/,20x)
c 
        write (iwrite,2600) nfuns,ncols
        call arrstore(iwrite,vmatr,ncols*nfuns,
     1      '        matrix vmatr*')
c 
c 
c       store matrices u,v
c 
        iwrite=31
c 
 3600 format('        following is the matrix u, with ',
     1    'npts=',i3)
c 
        write (iwrite,3600) npts1
        call arrstore(iwrite,u,npts1**2,'        matrix u*')
c 
 3800 format('        following is the matrix v, with ',
     1    'npts=',i3)
c 
        write (iwrite,3800) npts1
        call arrstore(iwrite,v,npts1**2,'        matrix v*')
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine uvmattst(ws,npts,ncols,u,uu,uuu,vvv,w,s)
        implicit real *8 (a-h,o-z)
        save
        dimension u(npts,ncols),uu(npts,ncols),uuu(1),vvv(1),
     1      s(1),w(1),ws(1)
c 
c       scale u by square roots of weights
c 
        do 1400 i=1,n   cols
        do 1200 j=1,npts
c 
        uu(j,i)=u(j,i)*sqrt(ws(j))
 1200 continue
 1400 continue
c 
c       construct the SVD of UU
c 
        eps=1.0d-10
        lw=100 000
        call svdpivot(ier,uu,npts,ncols,uuu,vvv,s,ncols3,eps,
     1      w,lw,ltot)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine bigsub0(ier,a,b,n,nfuns,w,ncols,rlams,
     1      fmatr,umatr,vmatr,xs,ws,npts1,u,v,work,lenwork,
     2      errs,xlarge,wlarge,ab,fun1,eps)
        implicit real *8 (a-h,o-z)
        save
        dimension rlams(1),xtest(110),w(1),fmatr(1),
     1      umatr(1),vmatr(1),xs(1),ws(1),u(1),v(1),
     2      usout(1000),errs(1),work(1),xlarge(1),wlarge(1),
     3      ab(1)
c 
        external fun1
c 
c        construct the matrix v in the SVD
c 
        call vmatcrea(w,n,nfuns,ncols,rlams,
     1      fmatr,umatr,vmatr,xlarge,wlarge,ab,fun1)
c 
c       test the accuracy of the evaluation of singular vectors
c 
        ntest=100
        h=(b-a)/ntest
        do 1400 i=1,ntest
        xtest(i)=a+h/2+(i-1)*h
 1400 continue
c 
        do 1500 i=1,ncols
c 
        errs(i)=0
 1500 continue
c 
        do 2000 i=1,ntest
c 
        call ufuneval(fun1,nfuns,ncols,vmatr,xtest(i),usout)
c 
        do 1600 j=1,ncols
c 
        call nesteva(ier,w,xtest(i),j,dj,der)
c 
        d=abs(usout(j)-dj)
        if(d .gt. errs(j)) errs(j)=d
c 
 1600 continue
c 
 2000 continue
c 
c       construct the interpolation matrices u,v
c 
        irhs=1
        lrhs=ncols**2+10
c 
        iwork=irhs+lrhs
        lenwork2=lenwork-iwork
c 
        call uvmatrbld(w,xs,ws,n,u,v,ncols,npts1,work(iwork),
     1      lenwork2,work(irhs),eps)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine uvmatrbld(w,xs,ws,n,u,v,ncols,npts1,work,
     1      lenwork,rhs,eps)
        implicit real *8 (a-h,o-z)
        save
        dimension u(npts1,ncols),v(ncols,npts1),xs(1),ws(1),
     1      work(1),rhs(ncols,ncols)
c 
        npts=npts1
        n=npts1
c 
        do 1200 i=1,ncols
        do 1100 j=1,ncols
c 
        rhs(j,i)=0
 1100 continue
c 
        rhs(i,i)=1
 1200 continue
c 
c       construct the matrix ov values of basis functions at
c       the support nodes
c 
        do 1600 i=1,ncols
        do 1400 j=1,npts
c 
        call nesteva(ier,w,xs(j),i,f,der)
c 
        u(j,i)=f
 1400 continue
c 
 1600 continue
c 
c       construct the "generalized" inverse v of the matrix u
c 
        delta=0
        lw=100 000
c 
        call leamatrr(ier,u,rhs,ncols,npts1,ncols,eps,delta,
     1      v,work,lenwork,ltot)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine vmatcrea(w,n,nfuns,ncols,rlams,
     1      fmatr,umatr,vmatr,x,whts,ab,fun1)
        implicit real *8 (a-h,o-z)
        save
        dimension par1(1),par2(1),fmatr(n,nfuns),umatr(n,ncols),
     1      vmatr(nfuns,ncols),w(1),rlams(1),x(1),whts(1),ab(1)
c 
c        retrieve from array w the nodes and weights of the
c        big quadrature
c 
        call alldecod(w,k,nn,n,ab,x,whts)
c 
c        construct the values of the original functions at the
c        nodes of the discretization
c 
        do 1400 i=1,nfuns
        do 1200 j=1,n
c 
        call fun1(x(j),i,par1,par2,f)
        fmatr(j,i)=f
 1200 continue
 1400 continue
c 
c       construct the matrix u
c 
        do 1800 i=1,ncols
        do 1600 j=1,n
c 
        call nesteva(ier,w,x(j),i,f,der)
c 
        umatr(j,i)=f
 1600 continue
 1800 continue
c 
c       evaluate the elements of the matrix v
c 
        do 2400 i=1,ncols
        do 2200 j=1,nfuns
c 
        d=0
        do 2100 jj=1,n
c 
        d=d+umatr(jj,i)*fmatr(jj,j)*whts(jj)
        vmatr(j,i)=d/rlams(i)**2
 2100 continue
c 
 2200 continue
c 
 2400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine ufuneval(fun1,nfuns,ncols,vmatr,x,usout)
        implicit real *8 (a-h,o-z)
        save
        dimension vmatr(nfuns,ncols),usout(1)
c 
c        This subroutine uses the matrix vmatr (hopefully, constructed
c        via a prior call to the subroutine bigsub (see) to evaluate
c        the ncols basis functions at the point x.
c 
c                      Input parameters:
c 
c  fun1 - the subroutine evaluating the functions to be compressed;
c        must be EXACTLY the same as in the preceding call to bigsub
c  nfuns - must be EXACTLY the same as in the preceding call to bigsub
c  ncols - the number of QR-functions in the QR-decomposition of the
c      nfuns functions given by the subroutine fun1, as produced by
c      the prior call to bigsub
c  vmatr - the nfuns \times ncols matrix to be used by the subroutine
c        ufuneval for the evaluation of the singular functions
c  x - the point where the ncols singular functions are to be evaluated
c 
c                      Output parameters:
c 
c  usout - the values of the ncols singular functions at the point x
c 
c       . . . evaluate all of the original functions at the user-specified
c             point x
c 
        do 1100 i=1,ncols
c 
        usout(i)=0
 1100 continue
c 
 1200 continue
c 
c       apply the matrix vmatr to the obtrained values
c 
        d=0
        do 1600 j=1,nfuns
c 
        call fun1(x,j,par1,par2,d)
c 
        do 1400 i=1,ncols
c 
        usout(i)=usout(i)+vmatr(j,i)*d
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
        subroutine ufunsinterp(fun1,vmatr,nfuns,ncols,v,npts,
     1      x,coefs)
        implicit real *8 (a-h,o-z)
        save
        dimension usout(22),vmatr(nfuns,ncols),coefs(1),
     1      v(ncols,npts)
c 
c       This subroutine returns to the user the coefficients
c       coefs of the interpolation formula expressing the
c       values of a function at the user-specified point x
c       via that function's values at the npts support nodes.
c 
c        The function to be interpolated is assumed to be a linear
c        combination of functions generated by a prior call to
c        the suibroutine bigsub (see). The input parameters vmatr,
c        ncols,v of this subroutine are supposed to have been
c        produced by a prior call to bigsub; it is critically
c        important that the subroutine fun1 and the parameters
c        nfuns, npts passed to this subroutine be exactly
c        identical to the subroutine fun1 and the parameters nfuns,
c        npts passed to the subroutine bigsub during the call that
c        produced the parameters vmatr, ncols, v.
c 
c        PLESE NOTE THAT THIS SUBROUTINE USES SEVERAL TYPES OF DATA
C        CREATED BY A PRIOR CALL TO THE SUBROUTINE BIGSUB (SEE);
C        THIS SUBROUTINE HAS NO KNOWN USES AS A STAND-ALONE DEVICE.
c 
c                   Input parameters:
c 
c  fun1 - the subroutine evaluating the functions to be compressed;
c        must be EXACTLY the same as in the preceding call to bigsub
c  vmatr - the nfuns \times ncols matrix produced by a prior call
  
c        to the subroutine bigsub
c  nfuns - must be EXACTLY the same as in the preceding call to bigsub
c  ncols - the number of QR-functions in the QR-decomposition of the
c      nfuns functions given by the subroutine fun1, as produced by
c      the prior call to bigsub
c  x - the point where the ncols singular functions are to be evaluated
c 
c                   Output parameters:
c 
c  coefs - the coefficients of the interpolation formula
c        (20 of them things)
c 
c        . . . if this is the first call to this subroutine, retrieve
c              the matrices vmatr, v from the subroutine onesideretr
c 
c 
c       . . . Evaluate all basis functions at the point x
c 
         call ufuneval(fun1,nfuns,ncols,vmatr,x,usout)
c 
c        apply the adjoint of the matrix v to the vector
c        of values of basis functions at the user-specified
c        point x, obtaining the interpolation formula
c 
        do 1400 i=1,npts
c 
        d=0
        do 1200 j=1,ncols
c 
        d=d+v(j,i)*usout(j)
 1200 continue
c 
        coefs(i)=d
 1400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine arrstore(iw,arr,n,mess)
        implicit real *8 (a-h,o-z)
        save
        character *1 mess(1),ast
        dimension arr(2,1)
        data ast/'*'/
c 
        i2=0
        l=0
        do 1100 i=1,100
        if(mess(i) .eq. ast) goto 1150
        l=i
 1100 continue
 1150 continue
c 
        write(iw,1250)
 1200 format(80a1)
        write(iw,1200) (mess(i),i=1,l)
 1250 format(20x)
        write(iw,1250)
 1270 format(20x,'total number of elements is',i5)
        write(iw,1270) n
        write(iw,1250)
c 
        n2=n/2
        ifodd=n-n2*2
c 
        nblocks=n2/90
c 
        call prinf('nblocks=*',nblocks,1)
c 
        nleft=n2-nblocks*90
c 
c       store on disk the array arr in blocks of 180 elements;
c       if there is a tail left that is smaller than 180 elements,
c       it will be handled later
c 
        if(nblocks .eq. 0) goto 2000
c 
        do 1800 i=1,nblocks
c 
 1300 format('c   ',/,'c        following is block number ',i3,
     1    ' of length 180 elements',/,'c    ')
c 
        write(iw,1300) i
c 
 1400 format(5x,i1,2x,e22.16,',',e22.16,',')
c 
        i1=(i-1)*90+1
        i2=i1+90-1
        icount=0
        do 1600 j=i1,i2
c 
        icount=icount+1
        if(icount .eq. 11) icount=1
c 
        write(iw,1400) icount,arr(1,j),arr(2,j)
 1600 continue
c 
 1800 continue
c 
 2000 continue
c 
c       if there is a tail - store it
c 
        if(nleft .eq. 0) goto 3000
c 
 2200 format('c   ',/,'c        following is the last block ',
     1    '(block number ',i3,')',/,'c        of length ',i3,
     2    ' elements',/,'c    ')
c 
        i=n-nblocks*180
        j=nblocks+1
        write(iw,2200) j,i
c 
        i1=nblocks*90+1
        i2=i2+nleft
  
        do 2600 j=i1,i2
c 
        icount=icount+1
        if(icount .eq. 11) icount=1
c 
        write(iw,1400) icount,arr(1,j),arr(2,j)
 2600 continue
 3000 continue
c 
c       if the number of elements in array arr is actually odd,
c       print the last element
c 
        icount=icount+1
        if(ifodd .eq. 1) write(iw,1400) icount,arr(1,n2+1)
  
 3200 format(20x)
        write(iw,3200)
c 
        return
        end
  
