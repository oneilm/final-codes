        implicit real *8 (a-h,o-z)
        real *8 amatr(100 000),u(1000 000),v(1000 000),
     1      w(1000 000)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER n'
        READ *,n
        CALL PRINf('n=*',n,1 )
C 
        PRINT *, 'ENTER m'
        READ *,m
        CALL PRINf('m=*',m,1 )
C 
        PRINT *, 'ENTER coef'
        READ *,coef
        CALL PRIN2('coef=*',coef,1 )
C 
        call lea2test(n,m,amatr,coef,u,v,w)
  
  
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine lea2test(n,m,amatr,coef,u,v,w)
        implicit real *8 (a-h,o-z)
        save
        dimension amatr(n,m),ts(1000),whts(1000)
c 
        dimension u(n,1),v(m,1),rea(2),fima2(10 000),
     1      frea1(10 000),fima1(10 000),frea2(10 000)
  
        complex *16 rhs(1000),coefs(1000),rhsout(1000),
     1      prods(1000),cd,com,
     2      apprrhs(1000),rat(1000),diff(1000),ima
c 
        equivalence (rea(1),com)
c 
        data ima/(0.0d0,1.0d0)/
C 
c       discretize the interval [-1,1] and evaluate the
c       Legendre polynomials on it
c 
        itype=1
        call legeexps(itype,n,ts,u,v,whts)
cccc        call prin2('Legendre nodes as constructed*',ts,n)
cccc        call prin2('Legendre weights as constructed*',whts,n)
  
  
c 
c        construct the test matrix
c 
        do 1400 i=1,n
c 
        do 1200 j=1,m
c 
cccc        amatr(i,j)=sin(coef*j*ts(i))
c 
        amatr(i,j)=ts(i)**(j-1)
 1200 continue
c 
        rhs(i)=ts(i)**6
        rhs(i)=sin(coef*ts(i))+1.0d-12 +ima*cos(ts(i))
c 
 1400 continue
c 
c       evaluate the SVD of the test matrix to be used for the
c       construction of the least squares approximation
c 
        lw=1000 000
c 
        call lesqinit(ier,amatr,n,m,whts,w,lw,lused,keep)
  
        call prinf('after lesqinit, ier=*',ier,1)
        call prinf('after lesqinit, lused=*',lused,1)
        call prinf('after lesqinit, keep=*',keep,1)
c 
c       approximate the right-hand sied with a linear combinations
c       of singular vectors
c 
        sgain=382.05
        sgain=382000
  
  
c 
c       construct the least squares approximation to rhs with
c       linear combinations of columns of amatr
c 
        ifsecant=1
cccc        ifsecant=0
c 
        call lesqevac(ier,amatr,w,rhs,sgain,ifsecant,
     1      prods,sgainout,apprrhs,err,errnosu)
c 
        call prinf('after lesqeval, ier=*',ier,1)
c 
        call prin2('after lesqeval, err=*',err,1)
        call prin2('after lesqeval, errnosu=*',errnosu,1)
c 
        call prin2('after lesqeval, prods=*',prods,m*2)
        call prin2('after lesqeval, apprrhs=*',apprrhs,n*2)
        call prin2('while rhs=*',rhs,n*2)
  
c 
        do 2200 i=1,n
        rat(i)=apprrhs(i)/rhs(i)
        diff(i)=apprrhs(i)-rhs(i)
 2200 continue
c 
        call prin2('and rat=*',rat,n)
c 
        call prin2('and diff=*',diff,n)
  
  
cccc         stop
c 
c       plot the exact right-hand side and the obtained
c       approximation together
c 
        do 2400 i=1,n
c 
        com=rhs(i)
        frea1(i)=rea(1)
        fima1(i)=rea(2)
c 
        com=apprrhs(i)
        frea2(i)=rea(1)
        fima2(i)=rea(2)
 2400 continue
c 
        iw=21
        call lotagraph2(iw,ts,frea1,n,ts,frea2,n,
     1      'real parts*')
c 
        iw=22
        call lotagraph2(iw,ts,fima1,n,ts,fima2,n,
     1      'imaginary parts*')
c 
        call prin2('and sgainout=*',sgainout,1)
        call prin2('while sgain=*',sgain,1)
  
        call lesqarrm(prods,coefs,m*2)
  
        call prin2('and coefs=*',coefs,m*2)
c 
c       recompute the approximation to the rhs from the
c       obtained final coefficients
c 
        err=0
c 
        do 3400 i=1,n
        cd=0
        do 3200 j=1,m
c 
        cd=cd+coefs(j)*amatr(i,j)
 3200 continue
  
        rhsout(i)=cd
c 
        err=(rhsout(i)-rhs(i))**2*whts(i)+err
 3400 continue
  
  
        err=sqrt(err)
        call prin2('and error is*',err,1)
  
  
        return
        end
c 
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c       This is the end of the debugging code and the beginning of the
c       least squares codse proper
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        subroutine lesqevac(ier,a,w,rhs,sgain,ifsecant,
     1      coefsout,sgainout,apprrhs,err,errnosu)
c 
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),errssec(20)
        complex *16 rhs(1),coefsout(1),apprrhs(1)
c 
c       This subroutine uses a version of the slepian scheme to
c       construct a least squares approximation of the user-supplied
c       vector rhs with a linear combination of (also user-supplied)
c       vectors (provided as columns of the matrix a). The subroutine
c       uses as input the array w produced by the subroutine lesqinit
c       (see) and has no uses as a stand-alone device. The scheme
c       constructs an approximation with the user-prescribed
c       supergain factor.
c 
c   PLEASE NOTE THAT THIS SUBROUTINE USES real BASIS FUNCTIONS
c   (PREPROCESSED BY A PRIOR CALL TO lesqinit) AND complex RIGHT-HAND
c   SIDE. IN THIS SENSE, THIS SUBROUTINE IS A SEMI-COMPLEXIFIED VERSION
c       OF THE SUBROUTINE lesqeval (SEE).
c 
c                         Input parameters:
c 
c  a - the array containing the vectors with which the approximation
c       is to be constructed; not damaged by this subroutine n any
c       way
c  rhs - the vector to be approximated
c  sgain - the supergain of the output vector
c  ifsecant - tells the subroutine which of the two algorithms to use.
c    ifsecant=0 will cause the subroutine to use the simple, efficient,
c       and reliable algorithm that does not guarantee that the
c       approximation will have the user-prescribed supergain with a
c       particularly high precision.
c    ifsecant=1 will cause the subroutine to construct an approximation
c       with the user-prescribed supergain sgain, with something like
c       15-digit precision. In this version, the subroutine is somewhat
c       slower (especially, for small-scale problems). More importantly,
c       there is no theorem that it will always work, though it has
c       never failed in the author's experience.
c 
c                         Output parameters:
c 
c  ier - error return code.
c         ier=0 means execution successful in all respects
c         ier=4 means that the supergain condition has not been activated;
c               the user-specified supergain was too high (the problem was
c               solved optimally with lower supergain); this is a very mild
c               condition
c        ier=16 means that the Newton performed 200 iterations (or
c               attempted iterations, dedicated to step-length control)
c               without achieving the requested accuracy; perhaps the
c               user-specifed eps is too small
c        ier=32 means that the process wanted to perform steps that are
c               shorter than 1.0d-5; this problem has never been
c               encountered so far; serious programming trouble likely
c        ier=64 means serious programming trouble (exceeded array?)
c        ier=256 means serious programming trouble (exceeded array?)
c 
c  coefsout - the projections of the approximation on the columns of a
c       (m of them things)
c  sgainout - supergain factor of the obtained approximation
c  apprrhs - the obtained approximation (n values)
c  err - the discrepancy between the obtained approximation (n the L^2
c        sense, and properly scaled by the weights in the array whts)
c  errnosu - the discrepancy between the approximation obtained
c        without supergain control (n the L^2 sense, and properly scaled
c        by the weights in the array whts). Needless to say, it should be
c        smaller than or equal to err.
c 
c       . . . unpack the first 10 elements of array w to obtain
c             various forms of discrete information to be used
c             by the subroutine
c 
        n=w(1)
        m=w(2)
        ncols=w(3)
c 
        iu=w(6)
        iwhts=w(7)
        icm=w(8)
        is=w(9)
        iv=w(10)
        iwwww=w(11)
        irprod=w(12)
c 
c        construct the least squares approximation with the supergain
c        specified by the user
c 
        call lesq2evc(ier,w(iu),n,ncols,rhs,w(iwhts),w(icm),sgain,
     1       coefsout,sgainout,ifsecant,errnosu,
     2      errssec,niter,w(irprod),w(iwwww) )
c 
c       construct the linear combination of original functions
c       that is equal (approximately) to the rhs
c 
        do 2400 i=1,ncols
c 
        coefsout(i)=coefsout(i)/w(is+i-1)
 2400 continue
c 
        call lesqmatc(w(iv),m,ncols,coefsout,apprrhs)
        call lesqarrm(apprrhs,coefsout,m*2)
c 
c        recompute the rhs from the final coefficients and
c        compare it to the original rhs
c 
        call lesqmatc(a,n,m,coefsout,apprrhs)
c 
c        evaluate the difference between the recomputed rhs and the
c        original one
c 
        err2=0
        do 2600 i=1,n
        err2=abs(apprrhs(i)-rhs(i))**2*w(iwhts+i-1) +err2
 2600 continue
c 
        err2=sqrt(err2)
        err=err2
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine lesq2evc(ier,amatr,n,m,rhs,whts,cm,sgain,
     1      prods,sgainout,ifsecant,errnosu,
     2      errssec,niter,rprods,wwww)
        implicit real *8 (a-h,o-z)
        save
        dimension amatr(n,m),whts(1),cm(1),
     1      errssec(1),rprods(1),wwww(1)
c 
        complex *16 rhs(1),prods(1),cd
c 
c       This subroutine uses a version of the slepian scheme to
c       construct a least squares approximation of the user-supplied
c       vector rhs with a linear combination of (also user-supplied)
c       vectors (provided as columns of the matrix a). The subroutine
c       uses as input the array w produced by the subroutine lesqinit
c       (see) and has no uses as a stand-alone device. The scheme
c       constructs an approximation with the user-prescribed
c       supergain factor.
c 
c   PLEASE NOTE THAT THIS SUBROUTINE USES real BASIS FUNCTIONS
c   (PREPROCESSED BY A PRIOR CALL TO lesqinit) AND complex RIGHT-HAND
c   SIDE. IN THIS SENSE, THIS SUBROUTINE IS A SEMI-COMPLEXIFIED VERSION
c       OF THE SUBROUTINE lesqeval (SEE).
c 
c                         Input parameters:
c 
c  a - the array containing the vectors with which the approximation
c       is to be constructed; not damaged by this subroutine n any
c       way
c  rhs - the vector to be approximated
c  sgain - the supergain of the output vector
c  ifsecant - tells the subroutine which of the two algorithms to use.
c    ifsecant=0 will cause the subroutine to use the simple, efficient,
c       and reliable algorithm that does not guarantee that the
c       approximation will have the user-prescribed supergain with a
c       particularly high precision.
c    ifsecant=1 will cause the subroutine to construct an approximation
c       with the user-prescribed supergain sgain, with something like
c       15-digit precision. In this version, the subroutine is somewhat
c       slower (especially, for small-scale problems). More importantly,
c       there is no theorem that it will always work, though it has
c       never failed in the author's experience.
c 
c                         Output parameters:
c 
c  ier - error return code.
c         ier=0 means execution successful in all respects
c         ier=4 means that the supergain condition has not been activated;
c               the user-specified supergain was too high (the problem was
c               solved optimally with lower supergain); this is a very mild
c               condition
c        ier=16 means that the Newton performed 200 iterations (or
c               attempted iterations, dedicated to step-length control)
c               without achieving the requested accuracy; perhaps the
c               user-specifed eps is too small
c        ier=32 means that the process wanted to perform steps that are
c               shorter than 1.0d-5; this problem has never been
c               encountered so far; serious programming trouble likely
c        ier=64 means serious programming trouble (exceeded array?)
c        ier=256 means serious programming trouble (exceeded array?)
c 
c  coefsout - the projections of the approximation on the columns of a
c       (m of them things)
c  sgainout - supergain factor of the obtained approximation
c  apprrhs - the obtained approximation (n values)
c  err - the discrepancy between the obtained approximation (n the L^2
c        sense, and properly scaled by the weights in the array whts)
c  errnosu - the discrepancy between the approximation obtained
c        without supergain control (n the L^2 sense, and properly scaled
c        by the weights in the array whts). Needless to say, it should be
c        smaller than or equal to err.
c 
c        . . . calculate the inner products of all of the basis
c              elements with the right-hand side
c 
        ier=0
        eps=1.0d-11
c 
        errnosu=0
        do 1200 i=1,m
c 
        call lesqscac(rhs,amatr(1,i),n,prods(i),whts)
c 
        rprods(i)=abs(prods(i))
 1200 continue
c 
c       evaluate the error of the approximation without supergain
c       control
c 
        do 1500 i=1,n
c 
        cd=0
        do 1400 j=1,m
c 
        cd=cd+prods(j)*amatr(i,j)
 1400 continue
        errnosu=errnosu+abs((rhs(i)-cd))**2*whts(i)
 1500 continue
c 
        errnosu=sqrt(errnosu)
c 
c       evaluate the norm of the rhs
c 
        call lesqsccc(rhs,rhs,n,cd,whts)
        d=cd
c 
        totpowr1=d*sgain
c 
c 
c       construct the approximation to the right-hand side
c       subject to the approximate supergain condition
c 
        call projfal(jer,cm,rprods,m,totpowr1,
     1      rnum,denom,rr,drrds,rlout,wwww,sout)
c 
        call prinf('after first projfal, jer=*',jer,1)
c 
c       if disaster has happened - bomb
c 
         if(jer .eq. 16) ier=256
         if(jer .eq. 16) return
c 
c       if the supergain condition has been enforced automatically,
c       or the user requested that it not be enforced exactly -
c       clean up and get out
c 
        if(jer .eq. 4) ier=4
        if( (jer .ne. 4) .and. (ifsecant .eq. 1) ) goto 1550
c 
        sgainout=rr
        do 1540 i=1,m
        prods(i)=prods(i)/(1+rlout*cm(i))
 1540 continue
        return
 1550 continue
c 
c       the user has requested that the supergain condition be
c       enforced exactly, and the enforcement has not happened
c       automatically. Enforce the thing explicitly.
c 
        call lesqnewt(ier,m,cm,sgain,
     1      rprods,sgainout,niter,totpowr1,rr,drrds,eps,rlout,wwww)
c 
        if(ier .ne. 0) return
c 
        do 2200 i=1,m
cccc        rprods(i)=rprods(i)/(1+rlout*cm(i))
        prods(i)=prods(i)/(1+rlout*cm(i))
 2200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine lesqeval(ier,a,w,rhs,sgain,ifsecant,
     1      coefsout,sgainout,apprrhs,err,errnosu)
c 
        implicit real *8 (a-h,o-z)
        save
        dimension w(1),rhs(1),coefsout(1),apprrhs(1),errssec(20)
c 
c       This subroutine uses a version of the slepian scheme to
c       construct a least squares approximation of the user-supplied
c       vector rhs with a linear combination of (also user-supplied)
c       vectors (provided as columns of the matrix a). The subroutine
c       uses as input the array w produced by the subroutine lesqinit
c       (see) and has no uses as a stand-alone device. The scheme
c       constructs an approximation with the user-prescribed
c       supergain factor.
c 
c                         Input parameters:
c 
c  a - the array containing the vectors with which the approximation
c       is to be constructed; not damaged by this subroutine n any
c       way
c  rhs - the vector to be approximated
c  sgain - the supergain of the output vector
c  ifsecant - tells the subroutine which of the two algorithms to use.
c    ifsecant=0 will cause the subroutine to use the simple, efficient,
c       and reliable algorithm that does not guarantee that the
c       approximation will have the user-prescribed supergain with a
c       particularly high precision.
c    ifsecant=1 will cause the subroutine to construct an approximation
c       with the user-prescribed supergain sgain, with something like
c       15-digit precision. In this version, the subroutine is somewhat
c       slower (especially, for small-scale problems). More importantly,
c       there is no theorem that it will always work, though it has
c       never failed in the author's experience.
c 
c                         Output parameters:
c 
c  ier - error return code.
c         ier=0 means execution successful in all respects
c         ier=4 means that the supergain condition has not been activated;
c               the user-specified supergain was too high (the problem was
c               solved optimally with lower supergain); this is a very mild
c               condition
c        ier=16 means that the Newton performed 200 iterations (or
c               attempted iterations, dedicated to step-length control)
c               without achieving the requested accuracy; perhaps the
c               user-specifed eps is too small
c        ier=32 means that the process wanted to perform steps that are
c               shorter than 1.0d-5; this problem has never been
c               encountered so far; serious programming trouble likely
c        ier=64 means serious programming trouble (exceeded array?)
c        ier=256 means serious programming trouble (exceeded array?)
c 
c  coefsout - the projections of the approximation on the columns of a
c       (m of them things)
c  sgainout - supergain factor of the obtained approximation
c  apprrhs - the obtained approximation (n values)
c  err - the discrepancy between the obtained approximation (n the L^2
c        sense, and properly scaled by the weights in the array whts)
c  errnosu - the discrepancy between the approximation obtained
c        without supergain control (n the L^2 sense, and properly scaled
c        by the weights in the array whts). Needless to say, it should be
c        smaller than or equal to err.
c 
c       . . . unpack the first 10 elements of array w to obtain
c             various forms of discrete information to be used
c             by the subroutine
c 
c       . . . unpack the first 10 elements of array w to obtain
c             various forms of discrete information to be used
c             by the subroutine
c 
        n=w(1)
        m=w(2)
        ncols=w(3)
c 
        iu=w(6)
        iwhts=w(7)
        icm=w(8)
        is=w(9)
        iv=w(10)
        iwwww=w(11)
c 
c        construct the least squares approximation with the supergain
c        specified by the user
c 
        call lesq2ev(ier,w(iu),n,ncols,rhs,w(iwhts),w(icm),sgain,
     1      coefsout,sgainout,apprrhs,ifsecant,errnosu,
     2      errssec,niter,w(iwwww) )
c 
c       construct the linear combination of original functions
c       that is equal (approximately) to the rhs
c 
        do 2400 i=1,ncols
c 
        coefsout(i)=coefsout(i)/w(is+i-1)
 2400 continue
c 
        call lesqmatv(w(iv),m,ncols,coefsout,apprrhs)
        call lesqarrm(apprrhs,coefsout,m)
c 
c        recompute the rhs from the final coefficients and
c        compare it to the original rhs
c 
        call lesqmatv(a,n,m,coefsout,apprrhs)
c 
c        evaluate the difference between the recomputed rhs and the
c        original one
c 
        err2=0
        do 2600 i=1,n
        err2=(apprrhs(i)-rhs(i))**2*w(iwhts+i-1) +err2
 2600 continue
c 
        err2=sqrt(err2)
  
        call prin2('and err2=*',err2,1)
  
        err=err2
  
  
        return
        end
c 
c 
c 
c 
c 
        subroutine lesq2ev(ier,amatr,n,m,rhs,whts,cm,sgain,
     1      prods,sgainout,apprrhs,ifsecant,errnosu,
     2      errssec,niter,wwww)
        implicit real *8 (a-h,o-z)
        save
        dimension amatr(n,m),rhs(1),whts(1),cm(1),prods(1),
     1      apprrhs(1),errssec(1),wwww(1)
c 
c       This subroutine uses a version of the slepian scheme to
c       construct a least squares approximation of the user-supplied
c       vector rhs with a linear combination of (also user-supplied)
c       vectors (provided as columns of the matrix amatr). The
c       scheme constructs an approximation with the user-prescribed
c       supergain factor, and does so subject to the weights in the
c       (again, user-supplied) array whts.
c 
c                         Input parameters:
c 
c  amatr - the array containing the vectors with which the approximation
c       is to be constructed
c  n - the length of the vector rhs to be approximated, and also of
c       the approximating vectors (columns of matrix amatr)
c  m - the number of approximating vectors (columns of matrix amatr)
c  rhs - the vector to be approximated
c  cm - the supergains associated with vectors stored in matrix amatr
c       (m of them things)
c  sgain - the supergain of the output vector
c  ifsecant - tells the subroutine which of the two algorithms to use.
c    ifsecant=0 will cause the subroutine to use the simple, efficient,
c       and reliable algorithm that does not guarantee that the
c       approximation will have the user-prescribed supergain with a
c       particularly high precision.
c    ifsecant=1 will cause the subroutine to construct an approximation'
c       with the user-prescribed supergain sgain, with something like
c       15-digit precision. In this version, the subroutine is somewhat
c       slower (especially, for small-scale problems). More importantly,
c       there is no theorem that it will always work, though it has
c       never failed in the author's experience.
c 
c                         Output parameters:
c 
c  ier - error return code.
c         ier=0 means execution successful in all respects
c         ier=4 means that the supergain condition has not been activated;
c               the user-specified supergain was too high (the problem was
c               solved optimally with lower supergain); this is a very mild
c               condition
c        ier=16 means that the Newton performed 200 iterations (or
c               attempted iterations, dedicated to step-length control)
c               without achieving the requested accuracy; perhaps the
c               user-specifed eps is too small
c        ier=32 means that the process wanted to perform steps that are
c               shorter than 1.0d-5; this problem has never been
c               encountered so far; serious programming trouble likely
c        ier=64 means serious programming trouble (exceeded array?)
c        ier=256 means serious programming trouble (exceeded array?)
c 
c  prods - the projections of the approximation on the columns of amatr
c       (m of them things)
c  sgainout - supergain factor of the obtained approximation
c  apprrhs - the obtained approximation (n coordinates)
c  err - the discrepancy between the obtained approximation (n the L^2
c        sense, and properly scaled by the weights in the array whts)
c  errnosu - the discrepancy between the approximation obtained
c        without supergain control (n the L^2 sense, and properly scaled
c        by the weights in the array whts). Needless to say, it should be
c        smaller than or equal to err.
c  errsec - the array (no more that 13 elements) containing the history
c        (errors obtained on various steps) of the secant scheme used
c        by the subroutine. It is only provided as a testing and debugging
c        tool, to be used when the parameter ifsecant had been set to 1
c  niter - the number of iterations taken by the secant; also the number
c        of elements returned in the array errsec. It is only provided as
c        a testing and debugging tool, to be used when the parameter
c        ifsecant had been set to 1
c 
c        . . . calculate the inner products of all of the basis
c              elements with the right-hand side
c 
        ier=0
        eps=1.0d-11
c 
        errnosu=0
        do 1200 i=1,m
c 
        call lesqsca(amatr(1,i),rhs,n,prods(i),whts)
 1200 continue
c 
c       evaluate the error of the approximation without supergain
c       control
c 
        do 1500 i=1,n
c 
        d=0
        do 1400 j=1,m
c 
        d=d+prods(j)*amatr(i,j)
 1400 continue
        errnosu=errnosu+(rhs(i)-d)**2*whts(i)
 1500 continue
c 
        errnosu=sqrt(errnosu)
c 
c       evaluate the norm of the rhs
c 
        call lesqsca(rhs,rhs,n,d,whts)
c 
        totpowr1=d*sgain
c 
c       construct the approximation to the right-hand side
c       subject to the approximate supergain condition
c 
        call projfal(jer,cm,prods,m,totpowr1,
     1      rnum,denom,rr,drrds,rlout,wwww,sout)
c 
        call prinf('after first projfal, jer=*',jer,1)
c 
c       if disaster has happened - bomb
c 
         if(jer .eq. 16) ier=256
         if(jer .eq. 16) return
c 
c       if the supergain condition has been enforced automatically,
c       or the user requested that it not be enforced exactly -
c       clean up and get out
c 
        if(jer .eq. 4) ier=4
        if( (jer .ne. 4) .and. (ifsecant .eq. 1) ) goto 1550
c 
        sgainout=rr
        do 1540 i=1,m
        prods(i)=prods(i)/(1+rlout*cm(i))
 1540 continue
        return
 1550 continue
c 
c       the user has requested that the supergain condition be
c       enforced exactly, and the enforcement has not happened
c       automatically. Enforce the thing explicitly.
c 
        call lesqnewt(ier,m,cm,sgain,
     1      prods,sgainout,niter,totpowr1,rr,drrds,eps,rlout,wwww)
c 
        if(ier .ne. 0) return
c 
        do 2200 i=1,m
        prods(i)=prods(i)/(1+rlout*cm(i))
 2200 continue
  
        return
        end
c 
c 
c 
c 
c 
        subroutine lesqnewt(ier,m,cm,sgain,
     1      prods,sgainout,niter,s0,rr0,drrds0,eps,rlout,wwww)
        implicit real *8 (a-h,o-z)
        save
        dimension cm(1),prods(1),wwww(1)
c 
c        This subroutine uses Newton iterations to determine the
c        value of the parameter s (total gain) that will result in the
c        user-requested  value of the supregain. Specifically, this
c        subroutine solves the equation
c 
c       { \sum_{m=1}^n  cm(m) * alph2(m) \over
c             \sum_{m=1}^n  alph2(m) }          = sgain                (1)
c 
c       (see subroutine projfal for the definitions of the parameters
c       alph2, etc.).
c 
c                        Input parameters:
c 
c  m - the numbers of elements in arrays cm,prods
c  cm - the array of supergain factors in (1) above
c  sgain - the right-hand side in (1) above
c  prods - the array of inner products in (1) above
c  s0 - the initial guess for the parameter s (gain)
c  rr0 - supergain corresponding to the initial guess s0
c  drrds0 - derivative of rr with respect to s at s=s0
c  eps - the precision to which the Newton will be conducted
c 
c                        Output parameters:
c 
c  ier - error return code
c        ier=0 means normal conclusion
c        ier=16 means that the Newton performed 200 iterations (or
c               attempted iterations, dedicated to step-length control)
c               without achieving the requested accuracy; perhaps the
c               user-specifed eps is too small
c        ier=32 means that the process wanted to perform steps that are
c               shorter than 1.0d-5; this problem has never been
c               encountered so far; serious programming trouble likely
c        ier=64 means serious programming trouble (exceeded array?)
c  sgainout - should be very close to sgain
c  niter - the number of iterations (including attempted iterations,
c        dedicated to step-length control) performed by the process
c  rlout - the value of the Lagrange multiplier at which the supergain
c        sgain (or sgainout, which is the same) is achieved
c 
c        . . . use newton to determine the value of s that will lead to
c              the user-requested supergain
c 
        ier=0
        ddold=1.0d20
c 
        s=s0
        rr=rr0
        drrds=drrds0
        stepcoe=1
c 
        do 3000 ijk=1,200
c 
        niter=ijk
c 
        call prinf('ijk=*',ijk,1)
        call prin2('stepcoe=*',stepcoe,1)
        call prin2('rr-sgain=*',rr-sgain,1)
  
        sold=s
        rrold=rr
        drrdsold=drrds
c 
        s=s-(rr-sgain)/drrds * stepcoe
        ddold=abs(rr-sgain)/sgain
c 
c       if the obtained s is negative - reduce steplength and try again
c 
        if(s .gt. 0) goto 2200
c 
        stepcoe=stepcoe/2
        s=sold
        rr=rrold
        drrds=drrdsold
        goto 3000
 2200 continue
c 
        call projfal(jer,cm,prods,m,s,
     1      rnum,denom,rr,drrds,rlout,wwww,sout)
c 
        dd=abs(rr-sgain)/sgain
c 
c       if the step failed to reduce the discrepancy - reduce the
c       steplength and try again
c 
        call prin2('dd=*',dd,1)
        call prin2('ddold=*',ddold,1)
        call prin2('and ddold-dd=*',ddold-dd,1)
  
        if( (ddold*1.00001d0 .ge. dd) .and. (jer .eq. 0) ) goto 2400
c 
        if(jer .eq. 16) ier=64
        if(jer .eq. 16) return
c 
        stepcoe=stepcoe/2
c 
c       if the step-length is less than 1.0e-5, bomb
c 
        sgainout=rr
        if(stepcoe .lt. 1.0d-5) ier=32
        if(stepcoe .lt. 1.0d-5) return
c 
        s=sold
        rr=rrold
        drrds=drrdsold
        goto 3000
 2400 continue
c 
        stepcoe=stepcoe*2
        if(stepcoe .gt. 1) stepcoe=1
c 
c       if the precision eps has been achieved - get out
c 
        call prin2('ddold=*',ddold,1)
        call prin2('dd=*',dd,1)
c 
        if(ddold .gt. eps) goto 3000
        sgainout=rr
        return
c 
 3000 continue
c 
        sgainout=rr
        ier=16
        return
        end
c 
c 
c 
c 
c 
        subroutine projfal(ier,cm,prods,n,s,
     1      rnum,denom,rr,drrds,rl,w,sout)
        implicit real *8 (a-h,o-z)
        save
        dimension cm(1),prods(1),w(1)
  
        ialph2=1
        lalph2=n+2
c 
        ialph2d=ialph2+lalph2
        lalph2d=n+2
c 
        irddral=ialph2d+lalph2d
c 
        call projfal0(ier,cm,prods,n,s,w(ialph2),w(ialph2d),
     1      rr,drrds,rl,w(irddral),rnum,denom,sout)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine projfal0(ier,cm,prods,n,s,alph2,alph2der,
     1      rr,drrds,rl,drrdalph,rnum,den,sout)
        implicit real *8 (a-h,o-z)
        save
        dimension cm(1),prods(1),alph2(1),alph2der(1),
     1      drrdalph(1)
c 
c       For user-specified parameters cm,prods,s,n, this subroutine
c       evaluates a number of expressions, and returns them to
c       the user. These expressions are discussed below, more or less in
c       the order of their evaluation by the subroutine.
c 
c       First, it finds (via a Newton procedure) the (unique) value
c       of rl, such that
c 
c           \sum_{m=0}^n { cm(m) * prods(m)^2  \over
c                        (1+cm(m)*rl)^2 }            =s,                (1)
c 
c       and the values alph2(m) defined by the formula
c 
c           alph2(m)={ prods(m)^2 \over (1+cm(m)*rl)^2 }.                (2)
c 
c       Then, it evaluates the array of alph2der of derivatives of the
c       coefficients alp2 with respect to s.
c 
c       Finally, it uses the obtained data to evaluate the ratio
c 
c       rr= { \sum_{m=1}^n  cm(m) * alph2(m) \over
c             \sum_{m=1}^n  alph2(m) },                                  (3)
c 
c       and the derivative drrds of rr with respect to s.
c 
c       This subroutine also returns the numerator rnum and the
c       denominator den in the formula (3).
c 
c                        Input parameters:
c 
c  cm - the array of coefficients in (2), (2),(3)
c  prods - the array of coefficients in (2), (2),(3)
c  n - the number of elements in each of arrays cm, prods
c  s - the right-hand side in (1)
c 
c                        Output parameters:
c 
c  ier - error return code
c  alph2 - given by the expression (2) above
c  alph2der - derivatives of alph2 with respect to s
c  rr - defined by the expression (3) above
c  drrds - derivative of rr with respect to s
c  drrdalph - array of partial derivatives of rr with respect to
c       coefficients in array alph2
c  rl - solution of equation (1) above
c  rnum - numerator in (3)
c  den - denominator in (3)
c 
c       . . . construct the Lagrange multiplier corresponding to
c             the user-specified parameters
c 
        ier=0
c 
        call projfbis(ier,cm,prods,n,s,rl,sout)
c 
        call prinf('after projfbis, ier=*',ier,1)
        if(ier .eq. 16) return
c 
c       evaluate the derivative of S with respect to rl
c 
        call projfder(cm,prods,n,rl,dsdrl)
c 
c       evaluate the |alphas**2| and their derivatives with respect to rl
c 
        do 1400 i=1,n
c 
        alph2(i)=prods(i)**2/(1+rl*cm(i))**2
c 
        alph2der(i)=-2*cm(i)*prods(i)**2/(1+rl*cm(i))**3
 1400 continue
c 
c       evaluate the derivatives of |alphas**2| with respect to s
c 
        do 1600 i=1,n
c 
        alph2der(i)=alph2der(i)/dsdrl
 1600 continue
c 
c       . . . now alph2der is the array of  derivatives of
c             |alphas**2| with respect to s
c 
c       . . . evaluate rr
c 
        rnum=0
        den=0
        do 1800 i=1,n
c 
        rnum=rnum+cm(i)*alph2(i)
        den=den+alph2(i)
 1800 continue
c 
        rr=rnum/den
c 
c       evaluate the derivatives of rr with respect to alph2
c 
        do 2000 i=1,n
c 
        drrdalph(i)=(cm(i)*den-rnum)/den**2
 2000 continue
c 
c       finally, evluate the derivative of rr with respect to s
c 
        drrds=0
        do 2200 i=1,n
c 
        drrds=drrds+drrdalph(i)*alph2der(i)
 2200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine projfbis(ier,cm,prods,n,s,rlout,sout)
        implicit real *8 (a-h,o-z)
        save
        dimension cm(1),prods(1)
c 
c       For the user-specified positive number s and arrays cm, prods,
c       this subroutine finds the coefficient rlout such that
c 
c           \sum_{m=0}^n { cm(m) * prods(m)^2  \over
c                        (1+cm(m)*rlout)^2 }            =s.             (1)
c 
c 
c                        Input parameters:
c 
c  cm - the array in (1)
c  prods - the array in (1)
c  n - the length of the arrays cm, prods
c  s - the right-hand side in (1)
c 
c                        output parameters:
c 
c  ier - error return code;
c        ier=0 means normal conclusion
c        ier=4 means that the parameter s is too big: with rlout=0,
c              the left-hand side in (1) is still smaller than the
c              right-hand side. This is a very mild condition, and
c              the subroutine will in this case set the output
c              parameter rlout to zero, and the output parameter
c              sout to the corresponding value of the left-hand side
c        ier=16 means that with rlout=1.0d20, the left-hand side of
c              (1) above is still larger than the user-specified s. It
c              is a very serious condition. Either the scaling used is
c              utterly unreasonable (sufficiently poor choice of units
c              will break any computer), or there is serious programming
c              trouble. This should be viewed as a fatal condition.
c  rlout - the solution of equation (1) above
c  sout -  the value of s corresponding toe the parameter rlout
c        being returned; when ier=0, sout is equal to s with very high
c        precision. When ier=4, sout is equal to the value of the left-hand
c        side of (1) corresponding to rlout=0.
c 
c 
c       . . . initialize the bisection process
c 
        ier=0
        rl1=0
        rl2=1.0d20
        rlout=0
c 
        call projfsum(cm,prods,n,rl1,f1)
        call projfsum(cm,prods,n,rl2,f2)
c 
        call prin2('f1=*',f1,1)
        call prin2('f2=*',f2,1)
        call prin2('s=*',s,1)
  
        if( (f1 .gt. s) .and. (f2 .lt. s) ) goto 1200
c 
        if(f1 .gt. s) goto 1100
c 
            ier=4
            rlout=0
            sout=f1
            return
 1100 continue
c 
c 
        if(f2 .lt. s) goto 1150
c 
        ier=16
        rlout=1.0d20
        sout=f2
        return
 1150 continue
c 
 1200 continue
c 
c        reduce rl2 by decimal scaling as far as it will co
c 
        do 1300 i=1,100
c 
        rl22=rl2/10
c 
        call projfsum(cm,prods,n,rl22,f22)
c 
        if(f22 .ge. s) goto 1350
c 
        rl2=rl22
        f2=f22
 1300 continue
 1350 continue
c 
        do 2000 i=1,10
c 
        rl3=(rl1+rl2)/2
c 
        call projfsum(cm,prods,n,rl3,f3)
c 
        if(f3 .lt. s) goto 1400
c 
        f1=f3
        rl1=rl3
        goto 2000
 1400 continue
c 
        f2=f3
        rl2=rl3
c 
 2000 continue
c 
        rlout=(rl1+rl2)/2
c 
c       use Newton to tighten the obtained value of rlout
c 
        do 1600 i=1,5
c 
        call projfsum(cm,prods,n,rlout,sum)
        call projfder(cm,prods,n,rlout,der)
c 
        ff=sum-s
c 
        rlout=rlout-ff/der
 1600 continue
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine projfder(cm,prods,n,rlagr,der)
        implicit real *8 (a-h,o-z)
        save
        dimension cm(1)
        real *8 prods(1)
c 
c       evaluate the sum that is supposed to be equal to the
c       total power
c 
        der=0
        do 1200 i=1,n
c 
        der=der + cm(i)**2*prods(i)**2/(1+rlagr*cm(i))**3
 1200 continue
c 
        der=-der*2
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine projfsum(cm,prods,n,rlagr,sum)
        implicit real *8 (a-h,o-z)
        save
        dimension cm(1)
        real *8 prods(1)
c 
c       evaluate the sum that is supposed to be equal to the
c       total power
c 
        sum=0
        do 1200 i=1,n
c 
        sum=sum + cm(i)*prods(i)**2/(1+rlagr*cm(i))**2
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine lesqinit(ier,a,n,m,whts,w,lw,lused,keep)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),whts(n),w(1)
c 
c        this subroutine prepares data to be used by the subroutine
c        lesqeval (see) for the construction of controlled-supergain
c        solutions of least squares problems. This subroutine has no
c        uses as a stand-alone device.  Please note that this is merely
c        a memory-management routine; all work is performed by the
c        subroutine lesq2in (see)
c 
c                         Input parameters:
c 
c  a - the matrix whose columns are the vectors with which the least
c        squares approximations are to be constructed. This matrix is not
c        changed by the subroutine PRACTICALLY SPEAKING. More precisely,
c        the subroutine scales rows of this matrix by the corresponding
c        elements of the array whts, and later rescales them back.
c  n - the lengths of approximating vectors
c  m - the number of appeoximating vectors
c  whts - the array of POSITIVE numbers that are to be used as weights
c        for the least squares approximation. Please note that elements
c        of whts can be positive and very small, without causing stability
c        problems; they can not be negative or zero.
c  lw - the length of the user-supplied work array w; if it is too short,
c       the error return code ier is set to 4 (the array is somewhat
c       to short) or t0 1024 (the array is grossly too short), and the
c       execution of the subroutine is terminated
c 
c                         Output parameters:
c 
c  ier - error return code.
c           ier=0   means successful conclusion
c           ier > 0 means that the subroutine svdpivot used by this
c                 subroutine has run out of space in the work array w.
c                 This is a fatal error
c  w - the array in which on exit the data have been stored to be used
c       by the subroutine lesqeval (see) The subroutine is informed of
c       its length via the parameter lw (see above). It should be
c       sufficiently long. If this array is too short, the parameter ier
c       (see above) is set to something other than zero, and the execution
c       of the subroutine is terminated.
c  lused - the size (in real *8 words) of the chunk of the array w actually
c       used by the subroutine. This is a slightly conservative number,
c       likely to be a little larger than the amount actually used
c  keep - the chunk of the array w that has to be unchanged between the
c       call to this subroutine and subsequent calls to the subroutine
c       lesqueval; likely to be much less than lused
c 
c       . . . allocate memory for the subroutine preparing data for
c             the least square approximation
c 
        iu=15
        lu=n*m+10
c 
        is=iu+lu
        ls=n+m
c 
        iv=is+ls
        lv=n*m+10
c 
        iw=iv+lv
        lenw=lw-iw
c 
        lused= iw+ 3*n*m+ 2*(min(m,n))**2 +2*n+2*m+100
c 
c       construct the data for the least square approximation
c 
        eps=1.0d-15
        call lesq2in(ier,a,n,m,whts,eps,w(iu),w(iv),w(is),
     1      ncols,w(iw),lenw)
c 
        if(ier .ne. 0) return
c 
c        perform garbage collection and store the pertinent
c        information in the first 10 elements of array w
c 
        iu2=iu
        lu2=n*ncols+4
c 
        is2=iu2+lu2
        ls2=ncols+2
c 
        call lesqarrm(w(is),w(is2),ncols)
c 
        iv2=is2+ls2
        lv2=m*ncols+4
c 
        call lesqarrm(w(iv),w(iv2),lv2)
c 
        iwhts=iv2+lv2
        lwhts=n+2
c 
        call lesqarrm(whts,w(iwhts),n)
c 
        icm=iwhts+lwhts
        lcm=ncols+2
c 
c       construct the array of penalty coefficients to be used by
c       the subroutine lesqieval in the design of a supergain-controlled
c       least squares approximation
c 
        do 1600 i=1,ncols
c 
        w(i+icm-1)=1/w(is2+i-1)**2
 1600 continue
c 
        w(1)=n+0.1
        w(2)=m+0.1
        w(3)=ncols+0.1
c 
        w(6)=iu2+0.1
        w(7)=iwhts+0.1
        w(8)=icm+0.1
        w(9)=is2+0.1
        w(10)=iv2+0.1
c 
        iwwww=icm+lcm
        lwwww=3*m+30
c 
        w(11)=iwwww+0.1
c 
        irprod=iwwww+lwwww
        lrprod=m+10
c 
        keep=irprod+lrprod
        w(12)=irprod
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine lesqarrm(x,y,n)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),y(1)
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
        subroutine lesqsca(cx,y,n,prod,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension y(1),whts(1)
        real *8 cx(1),prod
c 
        prod=0
        do 1800 i=1,n
        prod=prod+cx(i)*y(i)*whts(i)
 1800 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine lesqmatv(v,m,ncols,coefsout,coefs)
        implicit real *8 (a-h,o-z)
        save
        dimension v(m,ncols),coefsout(ncols),coefs(m)
c 
        do 2800 i=1,m
c 
        d=0
        do 2600 j=1,ncols
c 
        d=d+v(i,j)*coefsout(j)
 2600 continue
c 
        coefs(i)=d
 2800 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine lesq2in(ier,a,n,m,whts,eps,u,v,s,ncols,w,lw)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,m),whts(1),u(n,1),v(m,1),s(1),w(1)
c 
c        this subroutine prepares data to be used by the subroutine
c        lesq2ev (see) for the construction of controlled-supergain
c        solutions of least squares problems. This subroutine has no
c        uses as a stand-alone device
c 
c                         Input parameters:
c 
c  a - the matrix whose columns are the vectors with which the least
c        squares approximations are to be constructed. This matrix is not
c        changed by the subroutine PRACTICALLY SPEAKING. More precisely,
c        the subroutine scales rows of this matrix by the corresponding
c        elements of the array whts, and later rescales them back.
c  n - the lengths of approximation vectors
c  m - the number of appeoximating vectors
c  whts - the array of POSITIVE numbers that are to be used as weights
c        for the least squares approximation. Please note that elements
c        of whts can be positive and very small, without causing stability
c        problems; they can not be negative or zero.
c  lw - the length of the user-supplied work array w; if it is too short,
c       the error return code ier is set to 4 (the array is somewhat
c       to short) or t0 1024 (the array is grossly too short), and the
c       execution of the subroutine is terminated
c  eps - the precision to which the underlying SVD of the matrix a
c       is calculated. A good value is 1.0d-15
c 
c                         Output parameters:
c 
c  ier - error return code.
c           ier=0   means successful conclusion
c           ier > 0 means that the subroutine svdpivot used by this
c                 subroutine has run out of space in the work array w.
c                 This is a fatal error
c  u,s,v - the three matrices in the SVD of the user-provided matrix a
c  ncols - the rank of the matrix a as determined by the SVD subroutine;
c       also the numbers of columns returned in the matrices u,s, and
c       the number of elements returned in the array s
c 
c                         Work arrays:
c 
c  w - the work array to be used by the subroutine svdpivot used by this
c       subroutine. The subroutine is informed of its length via the
c       parameter lw (see above). The value
c       lw=3*n*m+ 2*(min(m,n))**2 +2*n+2*m+100 is always sufficient. If this
c       array is too short, the parameter ier (see above) is set to
c       something other than zero, and the execution of the subroutine
c       is terminated.
c 
c 
c       . . . premultiply the user-supplied functions by the square roots
c             of weights to prepare it for the SVD
c 
        ier=0
        do 1200 i=1,n
        w(i)=sqrt(whts(i))
 1200 continue
c 
        do 1600 j=1,n
        do 1400 i=1,m
        a(j,i)=a(j,i)*w(j)
 1400 continue
 1600 continue
c 
c       construct the SVD of  the user-supplied vectors
c 
        call svdpivot(ier,a,n,m,u,v,s,ncols,eps,
     1      w,lw,ltot)
c 
        call prinf('after svdpivot, ncols=*',ncols,1)
        call prinf('after svdpivot, n=*',n,1)
        call prinf('after svdpivot, m=*',m,1)
c 
        if(ier .ne. 0) goto 2500
c 
c       compensate the obtained singular vectors and the original
c       matrix a for the scaling we introduced previously
c 
        do 2100 i=1,n
        w(i)=sqrt(whts(i))
 2100 continue
c 
        do 2400 i=1,ncols
        do 2200 j=1,n
c 
        u(j,i)=u(j,i)/w(j)
 2200 continue
 2400 continue
c 
 2500 continue
c 
        do 2800 i=1,m
        do 2600 j=1,n
c 
        a(j,i)=a(j,i)/w(j)
 2600 continue
 2800 continue
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine lesqscac(cx,y,n,prod,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension y(1),whts(1)
        complex *16 cx(1),prod
c 
        prod=0
        do 1800 i=1,n
        prod=prod+cx(i)*y(i)*whts(i)
 1800 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine lesqsccc(cx,y,n,prod,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension whts(1)
        complex *16 cx(1),prod,y(1)
c 
        prod=0
        do 1800 i=1,n
        prod=prod+cx(i)*conjg(y(i))*whts(i)
 1800 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine lesqmatc(v,m,ncols,coefsout,coefs)
        implicit real *8 (a-h,o-z)
        save
        dimension v(m,ncols)
        complex *16 coefsout(ncols),coefs(m),d
c 
        do 2800 i=1,m
c 
        d=0
        do 2600 j=1,ncols
c 
        d=d+v(i,j)*coefsout(j)
 2600 continue
c 
        coefs(i)=d
 2800 continue
        return
        end
  
  
  
  
  
  
