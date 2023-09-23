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

C
C
C
C
        SUBROUTINE PRINI(IP1,IQ1)
        CHARACTER *1 MES(1), AA(1)
         save
        REAL *4 A(1)
        REAL *8 A2(1)
        REAL *16 A4(1)
        INTEGER *4 IA(1)
        INTEGER *2 IA2(1)
        IP=IP1
        IQ=IQ1
        RETURN

C
C
C
C
C
        ENTRY PRIN(MES,A,N)
        CALL  MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1200)(A(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1200)(A(J),J=1,N)
 1200 FORMAT(6(2X,E11.5))
         RETURN
C
C
C
C
        ENTRY PRIN2(MES,A2,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1400)(A2(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1400)(A2(J),J=1,N)
 1400 FORMAT(6(2X,E11.5))
        RETURN
C
C
C
C
        ENTRY PRINQ(MES,A4,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1500)(A4(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1500)(A4(J),J=1,N)
 1500 FORMAT(6(2X,e11.5))
        RETURN
C
C
C
C
        ENTRY PRINF(MES,IA,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1600)(IA(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1600)(IA(J),J=1,N)
 1600 FORMAT(10(1X,I7))
        RETURN
C
C
C
C
        ENTRY PRINF2(MES,IA2,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1600)(IA2(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1600)(IA2(J),J=1,N)
        RETURN
C
C
C
C
        ENTRY PRINA(MES,AA,N)
        CALL MESSPR(MES,IP,IQ)
 2000 FORMAT(1X,80A1)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,2000)(AA(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,2000)(AA(J),J=1,N)
        RETURN
        END
c
c
c
c
c
        SUBROUTINE MESSPR(MES,IP,IQ)
        CHARACTER *1 MES(1),AST
        DATA AST/'*'/
C
C         DETERMINE THE LENGTH OF THE MESSAGE
C
        I=0
        DO 1400 I=1,10000
        IF(MES(I).EQ.AST) GOTO 1600
        I1=I
 1400 CONTINUE
 1600 CONTINUE
         IF ( (I1.NE.0) .AND. (IP.NE.0) )
     1     WRITE(IP,1800) (MES(I),I=1,I1)
         IF ( (I1.NE.0) .AND. (IQ.NE.0) )
     1     WRITE(IQ,1800) (MES(I),I=1,I1)
 1800 FORMAT(1X,80A1)
         RETURN
         END
C
C
C
C
C
        SUBROUTINE ZTIME(I)
        J=1
        J=7-I+J
CCCC    I=MRUN(J)
        RETURN
        END
c
c
c
c
c
        subroutine msgmerge(a,b,c)
        character *1 a(1),b(1),c(1),ast
        data ast/'*'/
c
        do 1200 i=1,1000
c
        if(a(i) .eq. ast) goto 1400
        c(i)=a(i)       
        iadd=i
 1200 continue
c
 1400 continue
c
        do 1800 i=1,1000
c
        c(iadd+i)=b(i)
        if(b(i) .eq. ast) return
 1800 continue
        return
        end
c

c
c
c
c
c
        subroutine lotaplot(iw,x,y,n,title)
        implicit real *8 (a-h,o-z)
        dimension x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 title(1)
c
        inumgr=1
c
        call lotaplo0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
     1      title,inumgr)
c
        return
        end
c
c
c
c
c
        subroutine lotaplot2(iw,x,y,n,x2,y2,n2,title)
        implicit real *8 (a-h,o-z)
        dimension x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 title(1)
c
        inumgr=2
c
        call lotaplo0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
     1      title,inumgr)
c
        return
        end
c
c
c
c
c
        subroutine lotaplot3(iw,x,y,n,x2,y2,n2,x3,y3,
     1      n3,title)
        implicit real *8 (a-h,o-z)
        dimension x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 title(1)
c
        inumgr=3
c
        call lotaplo0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
     1      title,inumgr)
c
        return
        end
c
c
c
c
c
        subroutine lotagraph(iw,x,y,n,title)
        implicit real *8 (a-h,o-z)
        dimension x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 title(1)
c
        inumgr=1
c
        call lotagrap0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
     1      title,inumgr)
c
        return
        end
c
c
c
c
c
        subroutine lotagraph2(iw,x,y,n,x2,y2,n2,title)
        implicit real *8 (a-h,o-z)
        dimension x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 title(1)
c
        inumgr=2
c
        call lotagrap0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
     1      title,inumgr)
c
        return
        end
c
c
c
c
c
c
        subroutine lotagraph3(iw,x,y,n,x2,y2,n2,x3,y3,
     1      n3,title)
        implicit real *8 (a-h,o-z)
        dimension x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 title(1)
c
        inumgr=3
c
        call lotagrap0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
     1      title,inumgr)
c
        return
        end
c
c
c
c
c
        subroutine lotagrap0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
     1      title,inumgr)
        implicit real *8 (a-h,o-z)
        dimension x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 title(1),gn(2),file1(8),anum1(8),line(82),
     1      blank,quo,backslash
        character *8 dummy,anum8,file8
c
        equivalence (file1,file8),(anum1,anum8)
        data gn/'g','n'/,blank/' '/,quo/'"'/,backslash/'\\'/
c
c        convert the user-specified Fortran unit number to 
c        character format
c
        if( (iw .ge. 0) .and. (iw .le. 9) )  write(dummy,1100) iw
 1100 format(i1)
c
        if( (iw .ge. 10) .and. (iw .le. 99) ) write(dummy,1200) iw
 1200 format(i2)
c
        if( (iw .ge. 100) .and. (iw .le. 999) ) write(dummy,1300) iw
 1300 format(i3)
c
        if( (iw .ge. 1000) .and. (iw .le. 9999) ) write(dummy,1400) iw
 1400 format(i4)
c
 2000 format(1a8)
        read(dummy,2000) anum8
c
c        construct the file name on which the Gnuplot instructions
c        are to be written
c
        file1(1)=gn(1)
        file1(2)=gn(2)
        do 2200 i=1,6
        file1(i+2)=anum1(i)
 2200 continue
c
c        open the fortran file with the unit 87 and name file8
c 
        iun=87
        open(unit=iun,file=file8)
c
 2250 format('   # set terminal postscript',/,
     1    '   # set output "plot.ps"')
c
        write(iun,2250)
c
c        generate the title for the plot
c
        line(1)=blank
        line(2)=blank
        line(3)=quo
c
        call MESslen(title,nchar,line(4))
c
        line(nchar+4)=quo
 2300 format('  set title ',80a1)
        write(iun,2300) (line(i),i=1,nchar+4)
c
 2350 format('   show title')
        write(iun,2350)
c
c        write the instructions 
c
 2400 format(i6)
        iun2=iw+100000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c
        do 2600 i=1,6
        file1(i+2)=anum1(i)
 2600 continue
c
 2800 format('plot ','"',1a8,'"     ','notitle  with lines')
 2830 format('plot ','"',1a8,'"     ','notitle  with lines, ',1a1)
c
        if(inumgr .eq. 1) write(iun,2800) file8
        if(inumgr .ne. 1) write(iun,2830) file8,backslash
c
c        write the first data file to be plotted
c
        iun22=88
        open(unit=iun22,file=file8)
c
        write(iun22,3000) (x(i),y(i),i=1,n)        
c
        close(iun22)
c
c       if the user so requested - write the instructions for the
c       plotting the second data file
c
        if(inumgr .eq. 1) close(iun)
        if(inumgr .eq. 1) return

        iun2=iw+200000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c
        do 2840 i=1,6
        file1(i+2)=anum1(i)
 2840 continue
c
        if(inumgr .eq. 2) write(iun,2850) file8
        if(inumgr .eq. 3) write(iun,2855) file8,backslash
cccc 2850 format('     ','"',1a8,'"     ','notitle  with points')
cccc 2855 format('     ','"',1a8,'"     ','notitle  with points, ',1a1)
c
 2850 format('     ','"',1a8,'"     ','notitle  with lines')
 2855 format('     ','"',1a8,'"     ','notitle  with lines, ',1a1)
c
c        write the second data file to be plotted
c
        iun22=88
        open(unit=iun22,file=file8)
c
        write(iun22,3000) (x2(i),y2(i),i=1,n2)
c
        close(iun22)
c
c       if the user so requested - write the instructions for the
c       plotting the third data file
c
        if(inumgr .eq. 2) close(iun)
        if(inumgr .eq. 2) return

        iun2=iw+300000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c
        do 2860 i=1,6
        file1(i+2)=anum1(i)
 2860 continue
c
        write(iun,2870) file8
 2870 format('     ','"',1a8,'"     ','notitle  with lines')
cccc 2870 format('     ','"',1a8,'"     ','notitle  with points')
c
c        write the third data file to be plotted
c
        iun22=88
        open(unit=iun22,file=file8)
c
 3000 format(2x,e11.5,2x,e11.5)
c
        write(iun22,3000) (x3(i),y3(i),i=1,n3)
c
        close(iun22)
        close(iun)
c
        return
        end
c
c
c
c
c
        subroutine lotaplo0(iw,x,y,x2,y2,x3,y3,n,n2,n3,
     1      title,inumgr)
        implicit real *8 (a-h,o-z)
        dimension x(1),y(1),x2(1),y2(1),x3(1),y3(1)
        character *1 title(1),gn(2),file1(8),anum1(8),line(82),
     1      blank,quo,backslash
        character *8 dummy,anum8,file8
c
        equivalence (file1,file8),(anum1,anum8)
        data gn/'g','n'/,blank/' '/,quo/'"'/,backslash/'\\'/
c
c        convert the user-specified Fortran unit number to 
c        character format
c
        if( (iw .ge. 0) .and. (iw .le. 9) )  write(dummy,1100) iw
 1100 format(i1)
c
        if( (iw .ge. 10) .and. (iw .le. 99) ) write(dummy,1200) iw
 1200 format(i2)
c
        if( (iw .ge. 100) .and. (iw .le. 999) ) write(dummy,1300) iw
 1300 format(i3)
c
        if( (iw .ge. 1000) .and. (iw .le. 9999) ) write(dummy,1400) iw
 1400 format(i4)
c
 2000 format(1a8)
        read(dummy,2000) anum8
c
c        construct the file name on which the Gnuplot instructions
c        are to be written
c
        file1(1)=gn(1)
        file1(2)=gn(2)
        do 2200 i=1,6
        file1(i+2)=anum1(i)
 2200 continue
c
c        open the fortran file with the unit 87 and name file8
c 
        iun=87
        open(unit=iun,file=file8)
c
 2250 format('   # set terminal postscript',/,
     1    '   # set output "plot.ps"')
c
        write(iun,2250)
c
c        generate the title for the plot
c
        line(1)=blank
        line(2)=blank
        line(3)=quo
c
        call MESslen(title,nchar,line(4))
c
        line(nchar+4)=quo
 2270 format('  set title ',80a1)
c
        write(iun,2270) (line(i),i=1,nchar+4)
c
 2280 format('   show title')
        write(iun,2280)
c
c        find the limits for both x and y
c
        xmin=1.0d30
        ymin=1.0d30
        xmax=-1.0d20
        ymax=-1.0d20
c
        do 2300 i=1,n
        if(x(i) .lt. xmin) xmin=x(i)
        if(y(i) .lt. ymin) ymin=y(i)
        if(x(i) .gt. xmax) xmax=x(i)
        if(y(i) .gt. ymax) ymax=y(i)
 2300 continue
c
        if(inumgr .eq. 1) goto 2340
c
        do 2310 i=1,n2
        if(x2(i) .lt. xmin) xmin=x2(i)
        if(y2(i) .lt. ymin) ymin=y2(i)
        if(x2(i) .gt. xmax) xmax=x2(i)
        if(y2(i) .gt. ymax) ymax=y2(i)
 2310 continue
        if(inumgr .eq. 2) goto 2340
c
        do 2320 i=1,n3
        if(x3(i) .lt. xmin) xmin=x3(i)
        if(y3(i) .lt. ymin) ymin=y3(i)
        if(x3(i) .gt. xmax) xmax=x3(i)
        if(y3(i) .gt. ymax) ymax=y3(i)
 2320 continue
c
 2340 continue
c
        xcenter=(xmin+xmax)/2 
        ycenter=(ymax+ymin)/2
c
        xsize=(xmax-xmin)
        ysize=(ymax-ymin)
        size=xsize
        if(ysize .gt. size) size=ysize
        size=size*1.1
c
        xmin=xcenter-size/2
        xmax=xcenter+size/2
        ymin=ycenter-size/2
        ymax=ycenter+size/2
c
c        set the size of the stupid thing
c
 2350 format(2x,' set size 0.75,1.0')
        write(iun,2350) 
c
c        write the instructions 
c
 2400 format(i6)
        iun2=iw+100000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c
        do 2600 i=1,6
        file1(i+2)=anum1(i)
 2600 continue
c
 2800 format('plot ', '[',  e8.2, ':', e8.2,'] ',
     1    '[',  e8.2, ':', e8.2,'] ',
     2    '"',1a8,'" ','notitle with lines')
cccc     2    '"',1a8,'" ','notitle with points')
c
 2830 format('plot ', '[',  e8.2, ':', e8.2,'] ',
     1    '[',  e8.2, ':', e8.2,'] ',
     2    '"',1a8,'" ','notitle with lines, ',1a1)
ccccc     2    '"',1a8,'" ','notitle with points, ',1a1)
c
        if(inumgr .eq. 1) write(iun,2800) xmin,xmax,ymin,ymax,file8
        if(inumgr .ne. 1) write(iun,2830) xmin,xmax,ymin,ymax,
     1      file8,backslash
c
c        write the first data file to be plotted
c
        iun22=88
        open(unit=iun22,file=file8)
c
        write(iun22,3000) (x(i),y(i),i=1,n)        
c
        close(iun22)
c
c       if the user so requested - write the instructions for the
c       plotting the second data file
c
        if(inumgr .eq. 1) close(iun)
        if(inumgr .eq. 1) return
c
        iun2=iw+200000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c
        do 2840 i=1,6
        file1(i+2)=anum1(i)
 2840 continue
c
 2850 format('     ', '"',1a8,'" ','notitle with lines')
cccc 2850 format('     ', '"',1a8,'" ','notitle with points')
c
 2855 format('     ', '"',1a8,'" ','notitle with lines, ',1a1)
cccc 2855 format('     ', '"',1a8,'" ','notitle with points, ',1a1)
c
         if(inumgr .eq. 2) write(iun,2850) file8
         if(inumgr .ne. 2) write(iun,2855) file8,backslash
c
c        write the second data file to be plotted
c
        iun22=88
        open(unit=iun22,file=file8)
c
        write(iun22,3000) (x2(i),y2(i),i=1,n2)
c
        close(iun22)
c
c       if the user so requested - write the instructions for the
c       plotting the third data file
c
        if(inumgr .eq. 2) close(iun)
        if(inumgr .eq. 2) return

        iun2=iw+300000
        write (dummy,2400) iun2
        read(dummy,2000) anum8
c
        do 2860 i=1,6
        file1(i+2)=anum1(i)
 2860 continue
c
        write(iun,2870) file8
 2870 format('     ','"',1a8,'"','notitle  with lines')
c
c        write the third data file to be plotted
c
        iun22=88
        open(unit=iun22,file=file8)
c
 3000 format(2x,e11.5,2x,e11.5)
c
        write(iun22,3000) (x3(i),y3(i),i=1,n3)
c
        close(iun22)
        close(iun)
c
        return
        end
c
c
c
c
c
        SUBROUTINE MESslen(MES,nchar,line)
        CHARACTER *1 MES(1),AST,line(1)
        DATA AST/'*'/
C
C         DETERMINE THE LENGTH OF THE MESSAGE
C
        I=0
        DO 1400 I=1,10000
        IF(MES(I).EQ.AST) GOTO 1600
        I1=I
 1400 CONTINUE
 1600 CONTINUE
c
        nchar=i1
        do 1800 i=1,nchar
        line(i)=mes(i)
 1800 continue
         RETURN
         END

c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        this is the end of the debugging code and the beginning 
c        of the legendre expansion routines
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c        This file contains a set of subroutines for the handling 
c        of Legendre expansions. It contains 19 subroutines that are 
c        user-callable. Following is a brief description of these 
c        subroutines.
c
c   legeexps - constructs Legendre nodes, and  corresponding Gaussian
c        weights. Also constructs the matrix v converting the 
c         coefficients of a legendre expansion into its values at 
c         the n Gaussian nodes, and its inverse u, converting the
c         values of a function at n Gaussian nodes into the
c         coefficients of the corresponding Legendre series.
c
c   legepol - evaluates a single Legendre polynomial (together
c         with its derivative) at the user-provided point
c
c   legepols - evaluates a bunch of Legendre polynomials
c         at the user-provided point
c   legepls2 - an accelerated version of legepols, evaluating a 
c         Legendre polynomials at the user-provided point; maximum
c         order of the polynomial to be evaluated is 290
c
c   legeinmt - for the user-specified n, constructs the matrices of 
c        spectral indefinite integration differentiation on the n 
c        Gaussian nodes on the interval [-1,1]. 
c
c   legeinte - computes the indefinite integral of the legendre 
c        expansion polin getting the expansion polout
c
c   legediff -  differentiates the legendre expansion polin getting 
c        the expansion polout
c
c   legefder - computes the value and the derivative of a Legendre 
c        expansion at point X in interval [-1,1]; this subroutine 
c        is not designed to be very efficient, but it does not
c        use any exdternally supplied arrays
c
c   legefde2 - the same as legefder, except it is desigmed to be 
c        fairly efficient; it uses externally supplied arrays
c        that are precomputed
c   
c   legeexev - computes the value of a Legendre expansion with 
c        at point X in interval [-1,1]; same as legefder, but does
c        not compute the derivative of the expansion
c
c   legeexe2 - the same as legeexev, except it is desigmed to be 
c        fairly efficient; it uses externally supplied arrays
c        that are precomputed
c
c   lematrin - constructs the matrix interpolating functions from 
c        the n-point Gaussian grid on the interval [-1,1] to an 
c        arbitrary m-point grid (the nodes of the latter are
c        user-provided)
c   
c   levecin - constructs the coefficients of the standard 
c        interpolation formula connecting the values of a 
c        function at n Gaussian nodes on the interval [a,b] with
c        its value at the point x \in R^1
c
c   legeodev - evaluates at the point x a Legendre expansion
c        having only odd-numbered elements; this is a fairly 
c        efficient code, using external arrays that are 
c        precomputed
c
c   legeevev - evaluates at the point x a Legendre expansion
c        having only even-numbered elements; this is a fairly 
c        efficient code, using external arrays that are 
c        precomputed
c
c   legepeven - evaluates even-numbered Legendre polynomials 
c        of the argument x; this is a fairly efficient code, 
c        using external arrays that are precomputed
c
c   legepodd - evaluates odd-numbered Legendre polynomials 
c        of the argument x; this is a fairly efficient code, 
c        using external arrays that are precomputed
c
C   legefdeq - computes the value and the derivative of a
c        Legendre Q-expansion with coefficients coefs
C     at point X in interval (-1,1); please note that this is
c     the evil twin of the subroutine legefder, evaluating the
c     proper (P-function) Legendre expansion; this subroutine 
c        is not designed to be very efficient, but it does not
c        use any exdternally supplied arrays
c
c   legeq - calculates the values and derivatives of a bunch 
c        of Legendre Q-functions at the user-specified point 
c        x on the interval (-1,1)
c
c   legeqs - calculates the value and the derivative of a single
c        Legendre Q-function at the user-specified point 
c        x on the interval (-1,1)
c
c   legecfde - computes the value and the derivative of a Legendre 
c        expansion with complex coefficients at point X in interval 
c        [-1,1]; this subroutine is not designed to be very efficient, 
c        but it does not use any exdternally supplied arrays. This is
c        a complex version of the subroutine legefder.
c
c   legecfd2 - the same as legecfde, except it is designed to be 
c        fairly efficient; it uses externally supplied arrays
c        that are precomputed. This is a complex version of the 
c        subroutine legefde2.
c
c   legecva2 - the same as legecfd2, except it is does not evaluate
c        the derivative of the function
c
c
        subroutine legeexps(itype,n,x,u,v,whts)
        implicit real *8 (a-h,o-z)
        dimension x(1),whts(1),u(n,n),v(n,n)
c
c         this subroutine constructs the gaussiaqn nodes 
c         on the interval [-1,1], and the weights for the 
c         corresponding order n quadrature. it also constructs
c         the matrix v converting the coefficients
c         of a legendre expansion into its values at the n
c         gaussian nodes, and its inverse u, converting the
c         values of a function at n gaussian nodes into the
c         coefficients of the corresponding legendre series.
c         no attempt has been made to make this code efficient, 
c         but its speed is normally sufficient, and it is 
c         mercifully short.
c
c                 input parameters:
c
c  itype - the type of the calculation to be performed
c          itype=0 means that only the gaussian nodes are 
c                  to be constructed. 
c          itype=1 means that only the nodes and the weights 
c                  are to be constructed
c          itype=2 means that the nodes, the weights, and
c                  the matrices u, v are to be constructed
c  n - the number of gaussian nodes and weights to be generated
c  
c                 output parameters:
c
c  x - the order n gaussian nodes - computed independently
c          of the value of itype.
c  u - the n*n matrix converting the  values at of a polynomial of order
c         n-1 at n legendre nodes into the coefficients of its 
c         legendre expansion - computed only in itype=2
c  v - the n*n matrix converting the coefficients
c         of an n-term legendre expansion into its values at
c         n legendre nodes (note that v is the inverse of u)
c          - computed only in itype=2
c  whts - the corresponding quadrature weights - computed only 
c         if itype .ge. 1
c
c       . . . construct the nodes and the weights of the n-point gaussian 
c             quadrature
c
        ifwhts=0
        if(itype. gt. 0) ifwhts=1
        call legewhts(n,x,whts,ifwhts)
c
c       construct the matrix of values of the legendre polynomials
c       at these nodes        
c
        if(itype .ne. 2) return
        do 1400 i=1,n
c
        call legepols(x(i),n-1,u(1,i) )
 1400 continue
c
        do 1800 i=1,n
        do 1600 j=1,n
        v(i,j)=u(j,i)
 1600 continue
 1800 continue
c
c       now, v converts coefficients of a legendre expansion
c       into its values at the gaussian nodes. construct its 
c       inverse u, converting the values of a function at 
c       gaussian nodes into the coefficients of a legendre 
c       expansion of that function
c
        do 2800 i=1,n
        d=1
        d=d*(2*i-1)/2
        do 2600 j=1,n
        u(i,j)=v(j,i)*whts(j)*d
 2600 continue
 2800 continue
        return
        end
c
c
c
c
c
        subroutine legewhts(n,ts,whts,ifwhts)
        implicit real *8 (a-h,o-z)
        dimension ts(1),whts(1)
c
c        this subroutine constructs the nodes and the
c        weights of the n-point gaussian quadrature on 
c        the interval [-1,1]
c
c                input parameters:
c
c  n - the number of nodes in the quadrature
c
c                output parameters:
c
c  ts - the nodes of the n-point gaussian quadrature
c  w - the weights of the n-point gaussian quadrature
c
c       . . . construct the array of initial approximations
c             to the roots of the n-th legendre polynomial
c
        eps=1.0d-14
        ZERO=0
        DONE=1
        pi=datan(done)*4
        h=pi/(2*n) 
        do 1200 i=1,n
        t=(2*i-1)*h
        ts(n-i+1)=dcos(t)
1200  CONTINUE
c
c         use newton to find all roots of the legendre polynomial
c
        ts(n/2+1)=0
        do 2000 i=1,n/2
c
        xk=ts(i)
        deltold=1
        do 1400 k=1,10
        call legepol(xk,n,pol,der)
        delta=-pol/der
cccc         call prin2('delta=*',delta,1)
        xk=xk+delta
        if(deltold .le. eps) goto 1600
        deltold=dabs(delta)
 1400 continue
 1600 continue
        ts(i)=xk
        ts(n-i+1)=-xk
 2000 continue
c
c       now, use the explicit integral formulae 
c       to obtain the weights
c
        if(ifwhts .eq. 0) return
        a=-1
        b=1
        do 2200 i=1,n/2+1
        call prodend(a,ts,n,i,fm)
        call prodend(b,ts,n,i,fp)
        whts(i)=fp-fm
        whts(n-i+1)=whts(i)
 2200 continue
        return
        end
c
c
c
c
c
        subroutine legepol(x,n,pol,der)
        implicit real *8 (a-h,o-z)
c
        pkm1=1
        pk=x
c
        pk=1
        pkp1=x
c
c        if n=0 or n=1 - exit
c
        if(n .ge. 2) goto 1200
        pol=1
        der=0
        if(n .eq. 0) return
c
        pol=x
        der=1
        return
 1200 continue
c
c       n is greater than 1. conduct recursion
c
        do 2000 k=1,n-1
        pkm1=pk
        pk=pkp1
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
 2000 continue
c
c        calculate the derivative
c
        pol=pkp1
        der=n*(x*pkp1-pk)/(x**2-1)
        return
        end
c
c
c
c
c
        subroutine prodend(x,xs,n,i,f)
        implicit real *8 (a-h,o-z)
        dimension xs(1)
c
c      evaluate the product
c
        f=1
        dlarge=1.0d20
        dsmall=f/dlarge
c
        large=0
        do 2000 j=1,n
        dd=dabs(f)
        if( dd .gt. dsmall) goto 1200
         f=f*10000
         large=large-1
 1200 continue
c
        if( dd .lt. dlarge) goto 1400
        f=f/10000
        large=large+1
 1400 continue
        if(j .eq. i) goto 2000
        f=f*(x-xs(j))/(xs(i)-xs(j))
 2000 continue
        d10000=10000
        f=f*d10000**large
        f=f**2*(x-xs(i))
        return
        end
c
c
c
c
c
        subroutine legepols(x,n,pols)
        implicit real *8 (a-h,o-z)
        dimension pols(1)
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
        pkm1=pk
        pk=pkp1
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
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
        subroutine legepls2(x,n,pols)
        implicit real *8 (a-h,o-z)
        dimension pols(1),pjcoefs1(2000),pjcoefs2(300)
        save
        data ifcalled/0/
c
c        if need be - initialize the arrays pjcoefs1, pjcoefs2
c
        if(ifcalled .eq. 1) goto 1100
c
        done=1
        ninit=290
        do 1050 j=2,ninit
c
        pjcoefs1(j)=(2*j-done)/j
        pjcoefs2(j)=-(j-done)/j
c
 1050 continue
c 
        ifcalled=1
 1100 continue

        pkm1=1
        pk=x
c
        pk=1
        pkp1=x
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
        pkm1=pk
        pk=pkp1
        pkp1=x*pk*pjcoefs1(k+1)+pkm1*pjcoefs2(k+1)
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
        subroutine legeinmt(n,ainte,adiff,x,whts,endinter,
     1      itype,w)
        implicit real *8 (a-h,o-z)
        dimension ainte(1),w(1),x(1),whts(1),adiff(1),endinter(1)
c
c
c        for the user-specified n, this subroutine constructs
c        the matrices of spectral indefinite integration and/or
c        spectral differentiation on the n Gaussian nodes 
c        on the interval [-1,1]. Actually, this is omnly a 
c        memory management routine. All the actual work is done
c        by the subroutine legeinm0 (see)
c
c                           input parameters:
c
c  n - the number of Gaussian nodes on the interval [-1,1]
c  itype - the type of the calculation to be performed
c          EXPLANATION: 
c       itype=1 means that only the matrix ainte will 
c               be constructed
c       itype=2 means that only the matrix adiff will 
c               be constructed
c       itype=3 means that both matrices ainte and adiff
c               will be constructed
c
c                           output paramaters:
c
c  ainte - the matrix of spectral indefinite integration on 
c          the Gaussian nodes
c  adiff - the matrix of spectral differentiation on 
c          the Gaussian nodes
c  x - the n Gaussian nodes on the intervl [-1,1]
c  whts - the n Gaussian weights on the interval [-1,1]
c  endinter - the interpolation coefficients converting the 
c          values of a function at n Gaussian nodes into its
c          value at 1 (the right end of the interval)
c
c                           work arrays:
c
c  w - must be 3* n**2 + 2*n +50 *8 locations long
c
c        . . . allocate memory for the construction of the integrating
c              matrix
c
        ipolin=1
        lpolin=n+5
c
        ipolout=ipolin+lpolin
        lpolout=n+5
c
        iu=ipolout+lpolout
        lu=n**2+1
c
        iv=iu+lu
        lv=n**2+1
c
        iw=iv+lv
        lw=n**2+1
c
        ltot=iw+lw
c
c        construct the integrating matrix
c
        call legeinm0(n,ainte,adiff,w(ipolin),w(ipolout),
     1      x,whts,w(iu),w(iv),w(iw),itype,endinter)
c
        return
        end
c
c
c
c
c
        subroutine legeinm0(n,ainte,adiff,polin,polout,
     1      x,whts,u,v,w,itype,endinter)
        implicit real *8 (a-h,o-z)
        dimension ainte(n,n),u(n,n),v(n,n),w(n,n),
     1      endinter(1),x(n),whts(n),polin(n),polout(n),
     2      adiff(n,n)
c
c        for the user-specified n, this subroutine constructs
c        the matrices of spectral indefinite integration and/or
c        spectral differentiation on the n Gaussian nodes 
c        on the interval [-1,1]
c
c                           input parameters:
c
c  n - the number of Gaussian nodes on the interval [-1,1]
c  itype - the type of the calculation to be performed
c          EXPLANATION: 
c       itype=1 means that only the matrix ainte will 
c               be constructed
c       itype=2 means that only the matrix adiff will 
c               be constructed
c       itype=3 means that both matrices ainte and adiff
c               will be constructed
c
c                           output paramaters:
c
c  ainte - the matrix of spectral indefinite integration on 
c          the Gaussian nodes
c  adiff - the matrix of spectral differentiation on 
c          the Gaussian nodes
c  x - the n Gaussian nodes on the intervl [-1,1]
c  whts - the n Gaussian weights on the interval [-1,1]
c
c                           work arrays:
c
c  polin, polout - must be n+3 real *8 locations each
c
c  u, v, w - must be n**2+1 real *8 locations each
c
c        . . . construct the matrices of the forward and inverse 
c              Legendre transforms
c
        itype2=2
        call legeexps(itype2,n,x,u,v,whts)
c
cccc         call prin2('after legeexps, u=*',u,n*n)
c
c        if the user so requested,
c        construct the matrix converting the coefficients of
c        the Legendre series of a function into the coefficients
c        of the indefinite integral of that function
c
        if(itype. eq. 2) goto 2000
c
        do 1600 i=1,n
c
        do 1200 j=1,n+2
        polin(j)=0
 1200 continue
c
        polin(i)=1
c
        call legeinte(polin,n,polout)
c
        do 1400 j=1,n
        ainte(j,i)=polout(j)
 1400 continue
c
 1600 continue
c
cccc         call prin2('ainte initially is*',ainte,n*n)
c
c        multiply the three, obtaining the integrating matrix
c
        call matmul(ainte,u,w,n)
        call matmul(v,w,ainte,n)
c
 2000 continue
c
c        if the user so requested,
c        construct the matrix converting the coefficients of
c        the Legendre series of a function into the coefficients
c        of the derivative of that function
c
        if(itype. eq. 1) goto 3000
c
        do 2600 i=1,n
c
        do 2200 j=1,n
        polin(j)=0
 2200 continue
c
        polin(i)=1
c
        call legediff(polin,n,polout)
c
        do 2400 j=1,n
        adiff(j,i)=polout(j)
cccc        ainte(i,j)=polout(j)
 2400 continue
c
 2600 continue
c
cccc         call prin2('adiff initially is*',adiff,n*n)
c
c        multiply the three, obtaining the integrating matrix
c
        call matmul(adiff,u,w,n)
        call matmul(v,w,adiff,n)
c
 3000 continue
c
c        construct the vector of interpolation coefficients
c        converting the values of a polynomial at the Gaussian
c        nodes into its value at the right end of the interval
c
        do 3400 i=1,n
c
        d=0
        do 3200 j=1,n
        d=d+u(j,i)
 3200 continue
        endinter(i)=d
 3400 continue
c
        return
        end
c
c
c
c
c
        subroutine legeinte(polin,n,polout)
        implicit real *8 (a-h,o-z)
        dimension polin(1),polout(1)
c
c       this subroutine computes the indefinite integral of the 
c       legendre expansion polin getting the expansion polout
c
c
c                       input parameters:
c
c  polin - the legendre expansion to be integrated
c  n - the order of the expansion polin 
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c         also nothe that the order of the integrated expansion is
c         n+1 (who could think!)
c
c                       output parameters:
c
c  polout - the legendre expansion of the integral of the function 
c         represented by the expansion polin
c
        do 1200 i=1,n+2
        polout(i)=0
 1200 continue
c
        do 2000 k=2,n+1
        j=k-1
c
cccc        polout(k+1)=polin(k)/(2*j+1)+polout(k+1)
        polout(k+1)=polin(k)/(2*j+1)
        polout(k-1)=-polin(k)/(2*j+1)+polout(k-1)
c
 2000 continue
c
        polout(2)=polin(1)+polout(2)
c
        dd=0
        sss=-1
        do 2200 k=2,n+1
c
        dd=dd+polout(k)*sss
        sss=-sss
 2200 continue
c
        call prin2('dd=*',dd,1)
        polout(1)=-dd
c
        return
        end
c
c
c
c
c
        subroutine legediff(polin,n,polout)
        implicit real *8 (a-h,o-z)
        dimension polin(1),polout(1)
c
c       this subroutine differentiates the legendre 
c       expansion polin getting the expansion polout
c
c
c                       input parameters:
c
c  polin - the legendre expansion to be differentiated
c  n - the order of the expansion polin 
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c         also nothe that the order of the integrated expansion is
c         n+1 (who could think!)
c
c                       output parameters:
c
c  polout - the legendre expansion of the derivative of the function 
c         represented by the expansion polin
c
        do 1200 k=1,n+1
        polout(i)=0
 1200 continue
c
        pk=polin(n+1)
        pkm1=polin(n)
        pkm2=0
        do 2000 k=n+1,2,-1
c
        j=k-1
c         
        polout(k-1)=pk*(2*j-1)
        if(k .ge. 3) pkm2=polin(k-2)+pk
c
        pk=pkm1
        pkm1=pkm2
c
 2000 continue
         return
         end
c
c
c
c
c
      SUBROUTINE legeFDER(X,VAL,der,PEXP,N)
      IMPLICIT REAL *8 (A-H,O-Z)
      REAL *8 PEXP(1)
C
C     This subroutine computes the value and the derivative
c     of a gaussian expansion with coefficients PEXP
C     at point X in interval [-1,1].
c
c                input parameters:
c
C     X = evaluation point
C     PEXP = expansion coefficients
C     N  = order of expansion 
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c
c                output parameters:
c
C     VAL = computed value
C     der = computed value of the derivative
C
C


        done=1
        pjm2=1
        pjm1=x
        derjm2=0
        derjm1=1
c
        val=pexp(1)*pjm2+pexp(2)*pjm1      
        der=pexp(2)
c
        DO 600 J = 2,N
c
        pj= ( (2*j-1)*x*pjm1-(j-1)*pjm2 ) / j
        val=val+pexp(j+1)*pj
c
        derj=(2*j-1)*(pjm1+x*derjm1)-(j-1)*derjm2
c
        derj=derj/j
        der=der+pexp(j+1)*derj
c 
        pjm2=pjm1
        pjm1=pj
        derjm2=derjm1
        derjm1=derj
 600   CONTINUE
c
      RETURN
      END


c
c
c
c
c
      SUBROUTINE legeFDE2(X,VAL,der,PEXP,N,
     1    pjcoefs1,pjcoefs2,ninit)
      IMPLICIT REAL *8 (A-H,O-Z)
      REAL *8 PEXP(1),pjcoefs1(1),pjcoefs2(1)
c
C     This subroutine computes the value and the derivative
c     of a gaussian expansion with coefficients PEXP
C     at point X in interval [-1,1].
c
c                input parameters:
c
C  X - evaluation point
C  PEXP - expansion coefficients
C  N  - order of expansion 
c  pjcoefs1, pjcoefs2 - two arrays precomputed on a previous call 
c      on a previous call to this subroutine. Please note that this
c      is only an input parameter if the parameter ninit (see below) 
c      has been set to 0; otherwise, these are output parameters
c  ninit - tells the subroutine whether and to what maximum order the 
c       arrays coepnm1,coepnp1,coexpnp1 should be initialized.
c     EXPLANATION: The subroutine will initialize the first ninit
c       elements of each of the arrays pjcoefs1, pjcoefs2. On the first 
c       call to this subroutine, ninit should be set to the maximum 
c       order n for which this subroutine might have to be called; 
c       on subsequent calls, ninit should be set to 0. PLEASE NOTE 
c       THAT THAT THESE ARRAYS USED BY THIS SUBROUTINE
c       ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES USED BY THE 
c       SUBROUTINE LEGEEXE2. If these arrays have been initialized
c       by one of these two subroutines, they do not need to be 
c       initialized by the other one.
c
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c
c                output parameters:
c
C  VAL - computed value
C  der - computed value of the derivative
C
C
        if(ninit .eq. 0) goto 1400
c
        done=1
        do 1200 j=2,ninit
c
        pjcoefs1(j)=(2*j-done)/j
        pjcoefs2(j)=-(j-done)/j
c
 1200 continue
c 
        ifcalled=1
 1400 continue
c
        pjm2=1
        pjm1=x
        derjm2=0
        derjm1=1
c
        val=pexp(1)*pjm2+pexp(2)*pjm1      
        der=pexp(2)
c
        DO 1600 J = 2,N
c
cccc        pj= ( (2*j-1)*x*pjm1-(j-1)*pjm2 ) / j
c
        pj= pjcoefs1(j)*x*pjm1+pjcoefs2(j)*pjm2


        val=val+pexp(j+1)*pj
c
cccc        derj=(2*j-1)*(pjm1+x*derjm1)-(j-1)*derjm2
        derj=pjcoefs1(j)*(pjm1+x*derjm1)+pjcoefs2(j)*derjm2

ccc         call prin2('derj=*',derj,1)


cccc        derj=derj/j
        der=der+pexp(j+1)*derj
c 
        pjm2=pjm1
        pjm1=pj
        derjm2=derjm1
        derjm1=derj
 1600   CONTINUE
c
      RETURN
      END
c
c
c
c
c
      SUBROUTINE legeexev(X,VAL,PEXP,N)
      IMPLICIT REAL *8 (A-H,O-Z)
      REAL *8 PEXP(1)
C
C     This subroutine computes the value o a Legendre
c     expansion with coefficients PEXP at point X in interval [-1,1]
C
c                input parameters:
c
C     X = evaluation point
C     PEXP = expansion coefficients
C     N  = order of expansion 
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c
c                output parameters:
c
C     VAL = computed value
C
        done=1
        pjm2=1
        pjm1=x
c
        val=pexp(1)*pjm2+pexp(2)*pjm1      
        der=pexp(2)
c
        DO 600 J = 2,N
c
        pj= ( (2*j-1)*x*pjm1-(j-1)*pjm2 ) / j
        val=val+pexp(j+1)*pj
c
        pjm2=pjm1
        pjm1=pj
 600   CONTINUE
c
        RETURN
        END
c
c
c
c
c
      SUBROUTINE legeexe2(X,VAL,PEXP,N,
     1      pjcoefs1,pjcoefs2,ninit)
      IMPLICIT REAL *8 (A-H,O-Z)
      REAL *8 PEXP(1),pjcoefs1(1),pjcoefs2(1)
c
C     This subroutine computes the value o a Legendre
c     expansion with coefficients PEXP at point X in interval [-1,1]
C
c                input parameters:
c
C     X = evaluation point
C     PEXP = expansion coefficients
C     N  = order of expansion 
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c
c                output parameters:
c
C     VAL = computed value
C
        done=1
        if(ninit .eq. 0) goto 1400
c
        done=1
        do 1200 j=2,ninit
c
        pjcoefs1(j)=(2*j-done)/j
        pjcoefs2(j)=-(j-done)/j
c
 1200 continue
c 
        ifcalled=1
 1400 continue
c
        pjm2=1
        pjm1=x
c
        val=pexp(1)*pjm2+pexp(2)*pjm1      
        der=pexp(2)
c
        DO 600 J = 2,N
c
        pj= pjcoefs1(j)*x*pjm1+pjcoefs2(j)*pjm2

cccc        pj= ( (2*j-1)*x*pjm1-(j-1)*pjm2 ) / j
        val=val+pexp(j+1)*pj
c
        pjm2=pjm1
        pjm1=pj
 600   CONTINUE
c
        RETURN
        END
c
c
c
c
c
        subroutine lematrin(n,m,xs,amatrint,ts,w)
        implicit real *8 (a-h,o-z)
        dimension amatrint(m,n),xs(1),w(1),ts(1)
c
c
c        This subroutine constructs the matrix interpolating
c        functions from the n-point Gaussian grid on the interval [-1,1]
c        to an arbitrary m-point grid (the nodes of the latter are
c        user-provided)
c
c                 Input parameters:
c
c  n - the number of interpolation nodes
c  m - the number of nodes to which the functions will be interpolated
c  xs - the points at which the function is to be interpolated
c
c                  Output parameters:
c
c  amatrint - the m \times n matrix conerting the values of a function
c        at the n Legendre nodes into its values at m user-specified
c        (arbitrary) nodes
c  ts - the n Gaussian nodes on the interval [-1,1]
c
c                  Work arrays:
c
c  w - must be at least 2*n**2+n + 100 real *8 locations long 
c

        icoefs=1
        lcoefs=n+2
c
        iu=icoefs+lcoefs
        lu=n**2+10
c
        iv=iu+lu
c
        ifinit=1
        do 2000 i=1,m
c
        call levecin(n,xs(i),ts,w(iu),w(iv),w(icoefs),ifinit)
c
        do 1400 j=1,n
        amatrint(i,j)=w(j)
 1400 continue
c
        ifinit=0
 2000 continue
c
        return
        end

c
c
c
c
c
        subroutine levecin(n,x,ts,u,v,coefs,ifinit)
        implicit real *8 (a-h,o-z)
        dimension u(n,n),v(n,n),ts(1),coefs(1)
c
c        This subroutine constructs the coefficients of the 
c        standard interpolation formula connecting the values of a 
c        function at n Gaussian nodes on the interval [a,b] with
c        its value at the point x \in R^1
c
c                 Input parameters:
c
c  n - the number of interpolation nodes
c  x - the points at which the function is to be interpolated
c  ts - the n Gaussian nodes on the interval [-1,1]; please note that
c        it is an input parameter only if the parameter ifinit (see 
c        below) has been set to 1; otherwise, it is an output parameter
c  u - the n*n matrix converting the  values at of a polynomial of order
c         n-1 at n legendre nodes into the coefficients of its 
c        legendre expansion; please note that
c        it is an input parameter only if the parameter ifinit (see 
c        below) has been set to 1; otherwise, it is an output parameter
c  ifinit - an integer parameter telling the subroutine whether it should
c        initialize the Legendre expander; 
c     ifinit=1 will cause the subroutine to perform the initialization
c     ifinit=0 will cause the subroutine to  skip the initialization
c
c                  Output parameters:
c
c  coefs - the interpolation coefficients
c
c                 Work arrays: 
c
c  v - must be at least n*n real *8 locations long
c
c       . . . construct the n Gausian nodes on the interval [-1,1]; 
c             also the corresponding Gaussian expansion-evaluation 
c             matrices
c
        itype=2
        if(ifinit .ne.0) call legeexps(itype,n,ts,u,v,coefs)
c
c       evaluate the n Legendre polynomials at the point where the
c       functions will have to be interpolated
c
        call legepols(x,n+1,v)
c
c       apply the interpolation matrix to the ector of values 
c       of polynomials from the right 
c
        call lematvec(u,v,coefs,n)
        return
        end
c
c
c
c
c
        subroutine lematvec(a,x,y,n)
        implicit real *8 (a-h,o-z)
        dimension a(n,n),x(n),y(n)
c
        do 1400 i=1,n
        d=0
        do 1200 j=1,n
        d=d+a(j,i)*x(j)
 1200 continue
        y(i)=d
 1400 continue
        return
        end
c
c
c
c
c
        subroutine matmul(a,b,c,n)
        implicit real *8 (a-h,o-z)
        dimension a(n,n),b(n,n),c(n,n)
c
        do 2000 i=1,n
        do 1800 j=1,n
        d=0
        do 1600 k=1,n
        d=d+a(i,k)*b(k,j)
 1600 continue
        c(i,j)=d
 1800 continue
 2000 continue
        return
c
c
c
c
        entry matmua(a,b,c,n)
ccc          call prin2('in matmua, a=*',a,n**2)
ccc          call prin2('in matmua, b=*',b,n**2)
        do 3000 i=1,n
        do 2800 j=1,n
        d=0
        do 2600 k=1,n
        d=d+a(i,k)*b(j,k)
 2600 continue
        c(i,j)=d
 2800 continue
 3000 continue
ccc          call prin2('exiting, c=*',c,n**2)
        return
        end

c
c
c
c
c
        subroutine legeodev(x,nn,coefs,val,ninit,
     1      coepnm1,coepnp1,coexpnp1)
        implicit real *8 (a-h,o-z)
        dimension coepnm1(1),coepnp1(1),
     1            coexpnp1(1),coefs(1)
c
c
c       This subroutine evaluates at the point x a Legendre expansion
c       having only odd-numbered elements
c
c                  Input parameters:
c
c  x - point on the interval [-1,1] at which the Legendre expansion 
c       is to be evaluated
c  nn - order of the expansion to be evaluated
c  coefs - odd-numbered coefficients of the Legendre expansion
c       to be evaluated at the point x (nn/2+2 of them things)
c  ninit - tells the subroutine whether and to what maximum order the 
c       arrays coepnm1,coepnp1,coexpnp1 should be initialized.
c     EXPLANATION: The subroutine will initialize the first ninit/2+2
c                  (or so) elements of each of the arrays  coepnm1,
c       coepnp1, coexpnp1. On the first call to this subroutine, ninit 
c       should be set to the maximum order nn for which this subroutine 
c       might have to be called; on subsequent calls, ninit should be 
c       set to 0. PLEASE NOTE THAT THAT THESE ARRAYS USED BY THIS SUBROUTINE
c       ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES USED BY THE 
c       SUBROUTINE LEGEPODD. IF these arrays have been initialized
c       by one of these two subroutines, they do not need to be 
c       initialized by the other one.
c  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 real *8 elements long 
c       each. Please note that these are input arrays only if ninit 
c       (see above) has been set to 0; otherwise, these are output arrays.
c 
c                  Output parameters:
c
c  val - the value at the point x of the Legendre expansion with 
c       coefficients coefs (see above) 
c    EXPLANATION: On exit from the subroutine, pols(1) = P_0 (x),
c       pols(2) = P_2(x), pols(3) =P_4(x),  etc.
c  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 real *8 elements long 
c       each. Please note that these are output parameters only if ninit 
c       (see above) has not been set to 0; otherwise, these are input 
c       parameters
c 
c       
        if(ninit .eq. 0) goto 1400
        done=1
        n=0
        i=0
c
        do 1200 nnn=2,ninit,2
c
        n=n+2
        i=i+1
c
        coepnm1(i)=-(5*n+7*(n*done)**2+2*(n*done)**3)
        coepnp1(i)=-(9+24*n+18*(n*done)**2+4*(n*done)**3)
        coexpnp1(i)=15+46*n+36*(n*done)**2+8*(n*done)**3
c
        d=(2+n*done)*(3+n*done)*(1+2*n*done)
        coepnm1(i)=coepnm1(i)/d
        coepnp1(i)=coepnp1(i)/d
        coexpnp1(i)=coexpnp1(i)/d
c
 1200 continue
c
 1400 continue
c
        x22=x**2
c
        pi=x
        pip1=x*(2.5d0*x22-1.5d0)
c
        val=coefs(1)*pi+coefs(2)*pip1

        do 2000 i=1,nn/2-2
c
        pip2 = coepnm1(i)*pi +
     1      (coepnp1(i)+coexpnp1(i)*x22)*pip1
c
        val=val+coefs(i+2)*pip2
c
        pi=pip1
        pip1=pip2


 2000 continue
c
        return
        end
c
c
c
c
c
        subroutine legeevev(x,nn,coefs,val,ninit,
     1      coepnm1,coepnp1,coexpnp1)
        implicit real *8 (a-h,o-z)
        dimension coepnm1(1),coepnp1(1),
     1            coexpnp1(1),coefs(1)
c
c
c       This subroutine evaluates at the point x a Legendre expansion
c       having only even-numbered elements
c
c                  Input parameters:
c
c  x - point on the interval [-1,1] at which the Legendre expansion 
c       is to be evaluated
c  nn - order of the expansion to be evaluated
c  coefs - even-numbered coefficients of the Legendre expansion
c       to be evaluated at the point x (nn/2+2 of them things)
c  ninit - tells the subroutine whether and to what maximum order the 
c       arrays coepnm1,coepnp1,coexpnp1 should be initialized.
c     EXPLANATION: The subroutine will initialize the first ninit/2+2
c                  (or so) elements of each of the arrays  coepnm1,
c       coepnp1, coexpnp1. On the first call to this subroutine, ninit 
c       should be set to the maximum order nn for which this subroutine 
c       might have to be called; on subsequent calls, ninit should be 
c       set to 0. PLEASE NOTE THAT THAT THESE ARRAYS USED BY THIS SUBROUTINE
c       ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES USED BY THE 
c       SUBROUTINE LEGEPEVEN. IF these aqrrays have been initialized
c       by one of these two subroutines, they do not need to be 
c       initialized by the other one.
c  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 real *8 elements long 
c       each. Please note that these are input arrays only if ninit 
c       (see above) has been set to 0; otherwise, these are output arrays.
c 
c                  Output parameters:
c
c  val - the value at the point x of the Legendre expansion with 
c       coefficients coefs (see above) 
c    EXPLANATION: On exit from the subroutine, pols(1) = P_0 (x),
c       pols(2) = P_2(x), pols(3) =P_4(x),  etc.
c  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 real *8 elements long 
c       each. Please note that these are output parameters only if ninit 
c       (see above) has not been set to 0; otherwise, these are input 
c       parameters
c       
        if(ninit .eq. 0) goto 1400
c
        done=1
        n=-1
        i=0
        do 1200 nnn=1,ninit,2
c
        n=n+2
        i=i+1
c
        coepnm1(i)=-(5*n+7*(n*done)**2+2*(n*done)**3)
        coepnp1(i)=-(9+24*n+18*(n*done)**2+4*(n*done)**3)
        coexpnp1(i)=15+46*n+36*(n*done)**2+8*(n*done)**3
c
        d=(2+n*done)*(3+n*done)*(1+2*n*done)
        coepnm1(i)=coepnm1(i)/d
        coepnp1(i)=coepnp1(i)/d
        coexpnp1(i)=coexpnp1(i)/d
c
 1200 continue
c
 1400 continue
c
        x22=x**2
c
        pi=1
        pip1=1.5d0*x22-0.5d0
c
        val=coefs(1)+coefs(2)*pip1
c
c       n is greater than 2. conduct recursion
c
        do 2000 i=1,nn/2-2
c
        pip2 = coepnm1(i)*pi +
     1      (coepnp1(i)+coexpnp1(i)*x22) *pip1
        val=val+coefs(i+2)*pip2
c
        pi=pip1
        pip1=pip2
c
 2000 continue
c
        return
        end
c
c
c
c
c
        subroutine legepeven(x,nn,pols,ninit,
     1      coepnm1,coepnp1,coexpnp1)
        implicit real *8 (a-h,o-z)
        dimension pols(1),coepnm1(1),coepnp1(1),
     1            coexpnp1(1)
c
c       This subroutine evaluates even-numbered Legendre polynomials 
c       of the argument x, up to order nn+1
c
c                  Input parameters:
c
c  x - the argument for which the Legendre polynomials are 
c       to be evaluated
c  nn - the maximum order for which the Legendre polynomials are
c       to be evaluated
c  ninit - tells the subroutine whether and to what maximum order the 
c       arrays coepnm1,coepnp1,coexpnp1 should be initialized.
c     EXPLANATION: The subroutine ill initialize the first ninit/2+2
c                  (or so) elements of each of the arrays  coepnm1,
c       coepnp1, coexpnp1. On the first call to this subroutine, ninit 
c       should be set to the maximum order nn for which this subroutine 
c       might have to be called; on subsequent calls, ninit should be 
c       set to 0.
c  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 real *8 elements long 
c       each. Please note that these are input arrays only if ninit 
c       (see above) has been set to 0; otherwise, these are output arrays.
c       PLEASE NOTE THAT THAT THESE ARRAYS USED BY THIS SUBROUTINE
c       ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES USED BY THE 
c       SUBROUTINE LEGEEVEV. IF these aqrrays have been initialized
c       by one of these two subroutines, they do not need to be 
c       initialized by the other one.
c 
c                  Output parameters:
c
c  pols - even-numbered Legendre polynomials of the input parameter x
c         (nn/2+2 of them things)
c    EXPLANATION: On exit from the subroutine, pols(1) = P_0 (x),
c       pols(2) = P_2(x), pols(3) =P_4(x),  etc.
c  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 real *8 elements long 
c       each. Please note that these are output parameters only if ninit 
c       (see above) has not been set to 0; otherwise, these are input 
c       parameters. PLEASE NOTE THAT THAT THESE ARRAYS USED BY THIS SUBROUTINE
c       ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES USED BY THE 
c       SUBROUTINE LEGEEVEV. If these arrays have been initialized
c       by one of these two subroutines, they do not need to be 
c       initialized by the other one.
c 
c       
        if(ninit .eq. 0) goto 1400
c
        done=1
        n=-1
        i=0
        do 1200 nnn=1,ninit,2
c
        n=n+2
        i=i+1
c
        coepnm1(i)=-(5*n+7*(n*done)**2+2*(n*done)**3)
        coepnp1(i)=-(9+24*n+18*(n*done)**2+4*(n*done)**3)
        coexpnp1(i)=15+46*n+36*(n*done)**2+8*(n*done)**3
c
        d=(2+n*done)*(3+n*done)*(1+2*n*done)
        coepnm1(i)=coepnm1(i)/d
        coepnp1(i)=coepnp1(i)/d
        coexpnp1(i)=coexpnp1(i)/d
c
 1200 continue
c
 1400 continue
c
        x22=x**2
c
        pols(1)=1
        pols(2)=1.5d0*x22-0.5d0
c
c       n is greater than 2. conduct recursion
c
        do 2000 i=1,nn/2
c
        pols(i+2) = coepnm1(i)*pols(i) +
     1      (coepnp1(i)+coexpnp1(i)*x22) *pols(i+1)
c
 2000 continue
c
        return
        end
c
c
c
c
c
        subroutine legepodd(x,nn,pols,ninit,
     1      coepnm1,coepnp1,coexpnp1)
        implicit real *8 (a-h,o-z)
        dimension pols(1),coepnm1(1),coepnp1(1),
     1            coexpnp1(1)
c
c       This subroutine evaluates odd-numbered Legendre polynomials 
c       of the argument x, up to order nn+1
c
c                  Input parameters:
c
c  x - the argument for which the Legendre polynomials are 
c       to be evaluated
c  nn - the maximum order for which the Legendre polynomials are
c       to be evaluated
c  ninit - tells the subroutine whether and to what maximum order the 
c       arrays coepnm1,coepnp1,coexpnp1 should be initialized.
c     EXPLANATION: The subroutine will initialize the first ninit/2+2
c                  (or so) elements of each of the arrays  coepnm1,
c       coepnp1, coexpnp1. On the first call to this subroutine, ninit 
c       should be set to the maximum order nn for which this subroutine 
c       might have to be called; on subsequent calls, ninit should be 
c       set to 0.
c  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 real *8 elements long 
c       each. Please note that these are input arrays only if ninit 
c       (see above) has been set to 0; otherwise, these are output arrays.
c 
c                  Output parameters:
c
c  pols - the odd-numbered Legendre polynomials of the input parameter x
c         (nn/2+2 of them things)
c    EXPLANATION: On exit from the subroutine, pols(1) = P_1(x),
c       pols(2) = P_3(x), pols(3) = P_5 (x), etc.
c  coepnm1,coepnp1,coexpnp1 - should be nn/2+4 real *8 elements long 
c       each. Please note that these are output parameters only if ninit 
c       (see above) has not been set to 0; otherwise, these are input 
c       parameters. PLEASE NOTE THAT THAT THESE ARRAYS USED BY THIS 
c       SUBROUTINE ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES 
c       SUSED BY THE UBROUTINE LEGEODEV. IF these arrays have been 
c       initialized by one of these two subroutines, they do not need 
c       to be initialized by the other one.
c       
        if(ninit .eq. 0) goto 1400
        done=1
        n=0
        i=0
c
        do 1200 nnn=2,ninit,2
c
        n=n+2
        i=i+1
c
        coepnm1(i)=-(5*n+7*(n*done)**2+2*(n*done)**3)
        coepnp1(i)=-(9+24*n+18*(n*done)**2+4*(n*done)**3)
        coexpnp1(i)=15+46*n+36*(n*done)**2+8*(n*done)**3
c
        d=(2+n*done)*(3+n*done)*(1+2*n*done)
        coepnm1(i)=coepnm1(i)/d
        coepnp1(i)=coepnp1(i)/d
        coexpnp1(i)=coexpnp1(i)/d
c
 1200 continue
c
 1400 continue
c
        x22=x**2
c
        pols(1)=x
        pols(2)=x*(2.5d0*x22-1.5d0)
c
        do 2000 i=1,nn/2
c
        pols(i+2) = coepnm1(i)*pols(i) +
     1      (coepnp1(i)+coexpnp1(i)*x22)*pols(i+1)
c
 2000 continue
c
        return
        end
c
c
c
c
c
        subroutine legefdeq(x,val,der,coefs,n)
        implicit real *8 (a-h,o-z)
        dimension coefs(1)
C
C     This subroutine computes the value and the derivative
c     of a Legendre Q-expansion with coefficients coefs
C     at point X in interval (-1,1); please note that this is
c     the evil twin of the subroutine legefder, evaluating the
c     proper (P-function) Legendre expansion
c
c                input parameters:
c
C  X = evaluation point
C  coefs = expansion coefficients
C  N  = order of expansion 
c
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c
c                output parameters:
c
C     VAL = computed value
C     der = computed value of the derivative
C
c
        val=0
        der=0
c
        d= log( (1+x) /(1-x) ) /2
        pkm1=d
        pk=d*x-1
c
        pk=d
        pkp1=d*x-1

        derk=(1/(1+x)+1/(1-x)) /2
        derkp1=d + derk *x 
c
        val=coefs(1)*pk+coefs(2)*pkp1
        der=coefs(1)*derk+coefs(2)*derkp1
c
c        if n=0 or n=1 - exit
c
        if(n .ge. 2) goto 1200
c
        if(n .eq. 0) return
c
        return
 1200 continue
c
c       n is greater than 2. conduct recursion
c
        do 2000 k=1,n-1
        pkm1=pk
        pk=pkp1
c
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
c
        derkm1=derk
        derk=derkp1
c
        derkp1= ( (2*k+1)*pk+(2*k+1)*x*derk - k*derkm1 )/(k+1)
c
        val=val+coefs(k+2)*pkp1
        der=der+coefs(k+2)*derkp1
c
 2000 continue
c
        return
        end
c
c
c
c
c
        subroutine legeqs(x,n,pols,ders)
        implicit real *8 (a-h,o-z)
        dimension pols(1),ders(1)
c
c       This subroutine calculates the values and derivatives of 
c       a bunch of Legendre Q-functions at the user-specified point 
c       x on the interval (-1,1)
c
c                     Input parameters:
c
c  x - the point on the interval [-1,1] where the Q-functions and 
c       their derivatives are to be evaluated
c  n - the highest order for which the functions are to be evaluated
c  
c                     Output parameters:
c
c  pols - the values of the Q-functions (the evil twins of the 
c       Legeendre polynomials) at the point x (n+1 of them things)
c  ders - the derivatives of the Q-functions (the evil twins of the 
c       Legeendre polynomials) at the point x (n+1 of them things)
c  
c
        d= log( (1+x) /(1-x) ) /2
        pkm1=d
        pk=d*x-1
c
        pk=d
        pkp1=d*x-1

        derk=(1/(1+x)+1/(1-x)) /2
        derkp1=d + derk *x 
c
c        if n=0 or n=1 - exit
c
        if(n .ge. 2) goto 1200
        pols(1)=pk
        ders(1)=derk
        if(n .eq. 0) return
c
        pols(2)=pkp1
        ders(2)=derkp1
        return
 1200 continue
c
        pols(1)=pk
        pols(2)=pkp1
c
c       n is greater than 2. conduct recursion
c
        ders(1)=derk
        ders(2)=derkp1
c
        do 2000 k=1,n-1
        pkm1=pk
        pk=pkp1
c
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
        pols(k+2)=pkp1
c
        derkm1=derk
        derk=derkp1
c
        derkp1= ( (2*k+1)*pk+(2*k+1)*x*derk - k*derkm1 )/(k+1)
        ders(k+2)=derkp1
 2000 continue
c
        return
        end
c
c
c
c
c


        subroutine legeq(x,n,pol,der)
        implicit real *8 (a-h,o-z)
c
c       This subroutine calculates the value and derivative of 
c       a Legendre Q-function at the user-specified point 
c       x on the interval (-1,1)
c
c
c                     Input parameters:
c
c  x - the point on the interval [-1,1] where the Q-functions and 
c       their derivatives are to be evaluated
c  n - the order for which the function is to be evaluated
c  
c                     Output parameters:
c
c  pol - the value of the n-th Q-function (the evil twin of the 
c       Legeendre polynomial) at the point x 
c  ders - the derivatives of the Q-function at the point x 
c  
c
        d= log( (1+x) /(1-x) ) /2
        pk=d
        pkp1=d*x-1
c
c        if n=0 or n=1 - exit
c
        if(n .ge. 2) goto 1200
        pol=d

        der=(1/(1+x)+1/(1-x)) /2

        if(n .eq. 0) return
c
        pol=pkp1
        der=d + der *x 
        return
 1200 continue
c
c       n is greater than 1. conduct recursion
c
        do 2000 k=1,n-1
        pkm1=pk
        pk=pkp1
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
 2000 continue
c
c        calculate the derivative
c
        pol=pkp1
        der=n*(x*pkp1-pk)/(x**2-1)
        return
        end

c
c
c
c
c
      SUBROUTINE legecFDE(X,VAL,der,PEXP,N)
      IMPLICIT REAL *8 (A-H,O-Z)
      complex *16 PEXP(1),val,der
C
C     This subroutine computes the value and the derivative
c     of a gaussian expansion with complex coefficients PEXP
C     at point X in interval [-1,1].
c
c                input parameters:
c
C     X = evaluation point
C     PEXP = expansion coefficients
C     N  = order of expansion 
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c
c                output parameters:
c
C     VAL = computed value
C     der = computed value of the derivative
C
C
        done=1
        pjm2=1
        pjm1=x
        derjm2=0
        derjm1=1
c
        val=pexp(1)*pjm2+pexp(2)*pjm1      
        der=pexp(2)
c
        DO 600 J = 2,N
c
        pj= ( (2*j-1)*x*pjm1-(j-1)*pjm2 ) / j
        val=val+pexp(j+1)*pj
c
        derj=(2*j-1)*(pjm1+x*derjm1)-(j-1)*derjm2
c
        derj=derj/j
        der=der+pexp(j+1)*derj
c 
        pjm2=pjm1
        pjm1=pj
        derjm2=derjm1
        derjm1=derj
 600   CONTINUE
c
      RETURN
      END
c
c
c
c
c
      SUBROUTINE legecFD2(X,VAL,der,PEXP,N,
     1    pjcoefs1,pjcoefs2,ninit)
      IMPLICIT REAL *8 (A-H,O-Z)
      REAL *8 pjcoefs1(1),pjcoefs2(1)
      complex *16 PEXP(1),val,der
c
C     This subroutine computes the value and the derivative
c     of a Legendre expansion with complex coefficients PEXP
C     at point X in interval [-1,1].
c
c                input parameters:
c
C  X - evaluation point
C  PEXP - expansion coefficients
C  N  - order of expansion 
c  pjcoefs1, pjcoefs2 - two arrays precomputed on a previous call 
c      on a previous call to this subroutine. Please note that this
c      is only an input parameter if the parameter ninit (see below) 
c      has been set to 0; otherwise, these are output parameters
c  ninit - tells the subroutine whether and to what maximum order the 
c       arrays coepnm1,coepnp1,coexpnp1 should be initialized.
c     EXPLANATION: The subroutine will initialize the first ninit
c       elements of each of the arrays pjcoefs1, pjcoefs2. On the first 
c       call to this subroutine, ninit should be set to the maximum 
c       order n for which this subroutine might have to be called; 
c       on subsequent calls, ninit should be set to 0. PLEASE NOTE 
c       THAT THAT THESE ARRAYS USED BY THIS SUBROUTINE
c       ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES USED BY THE 
c       SUBROUTINE LEGEEXE2. If these arrays have been initialized
c       by one of these two subroutines, they do not need to be 
c       initialized by the other one.
c
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c
c                output parameters:
c
C  VAL - computed value
C  der - computed value of the derivative
C
        if(ninit .eq. 0) goto 1400
c
        done=1
        do 1200 j=2,ninit
c
        pjcoefs1(j)=(2*j-done)/j
        pjcoefs2(j)=-(j-done)/j
c
 1200 continue
c 
        ifcalled=1
 1400 continue
c
        pjm2=1
        pjm1=x
        derjm2=0
        derjm1=1
c
        val=pexp(1)*pjm2+pexp(2)*pjm1      
        der=pexp(2)
c
        DO 1600 J = 2,N
c
cccc        pj= ( (2*j-1)*x*pjm1-(j-1)*pjm2 ) / j
c
        pj= pjcoefs1(j)*x*pjm1+pjcoefs2(j)*pjm2


        val=val+pexp(j+1)*pj
c
cccc        derj=(2*j-1)*(pjm1+x*derjm1)-(j-1)*derjm2
        derj=pjcoefs1(j)*(pjm1+x*derjm1)+pjcoefs2(j)*derjm2

ccc         call prin2('derj=*',derj,1)


cccc        derj=derj/j
        der=der+pexp(j+1)*derj
c 
        pjm2=pjm1
        pjm1=pj
        derjm2=derjm1
        derjm1=derj
 1600   CONTINUE
c
      RETURN
      END
c
c
c
c
c
      SUBROUTINE legecva2(X,VAL,PEXP,N,
     1    pjcoefs1,pjcoefs2,ninit)
      IMPLICIT REAL *8 (A-H,O-Z)
      REAL *8 pjcoefs1(1),pjcoefs2(1)
      complex *16 PEXP(1),val
c
C     This subroutine computes the value of a Legendre expansion 
c     with complex coefficients PEXP at point X in interval [-1,1].
c
c                input parameters:
c
C  X - evaluation point
C  PEXP - expansion coefficients
C  N  - order of expansion 
c  pjcoefs1, pjcoefs2 - two arrays precomputed on a previous call 
c      on a previous call to this subroutine. Please note that this
c      is only an input parameter if the parameter ninit (see below) 
c      has been set to 0; otherwise, these are output parameters
c  ninit - tells the subroutine whether and to what maximum order the 
c       arrays coepnm1,coepnp1,coexpnp1 should be initialized.
c     EXPLANATION: The subroutine will initialize the first ninit
c       elements of each of the arrays pjcoefs1, pjcoefs2. On the first 
c       call to this subroutine, ninit should be set to the maximum 
c       order n for which this subroutine might have to be called; 
c       on subsequent calls, ninit should be set to 0. PLEASE NOTE 
c       THAT THAT THESE ARRAYS USED BY THIS SUBROUTINE
c       ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES USED BY THE 
c       SUBROUTINE LEGEEXE2. If these arrays have been initialized
c       by one of these two subroutines, they do not need to be 
c       initialized by the other one.
c
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c
c                output parameters:
c
C  VAL - computed value
C
        if(ninit .eq. 0) goto 1400
c
        done=1
        do 1200 j=2,ninit
c
        pjcoefs1(j)=(2*j-done)/j
        pjcoefs2(j)=-(j-done)/j
c
 1200 continue
c 
        ifcalled=1
 1400 continue
c
        pjm2=1
        pjm1=x
c
        val=pexp(1)*pjm2+pexp(2)*pjm1      
c
        DO 1600 J = 2,N
c
        pj= pjcoefs1(j)*x*pjm1+pjcoefs2(j)*pjm2
        val=val+pexp(j+1)*pj
c
        pjm2=pjm1
        pjm1=pj
 1600   CONTINUE
c
      RETURN
      END

c
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c          this is the end of the debugging output and the beginning
c          of the actual SVD code
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
        subroutine svdpivot(ier,a,n,m,u,v,s,ncols,eps,
     1      w,lw,ltot)
        implicit real *8 (a-h,o-z)
        dimension a(n,m),u(n,1),v(m,1),s(1),w(1)
c
c        this subroutine constructs the singular value decomposition 
c        of the usert-specified matrix a, under the assumption
c        that the rank of a is considerably smaller than its 
c        dimensionalities n,m. the output of this subroutine
c        consists of two matrices u,v, and a vector s such that
c
c                 a = u d v^*                                   (1)
c
c        with t the diagonal matrix with the vector s on the 
c        diagonal, and u, v two matrices whose columns are 
c        orthonormal. the reason for the existence of this
c        subroutine is the fact that on exit, the matrices u,v
c        have dimensionalities (n,ncols), (m,ncols), respectively,
c        with ncols the numerical rank of a. 
c
c                   input parameters:
c
c  a - the matrix to be svd-decomposed; note that it is NOT damaged
c       by this subroutine in any way
c  n,m - dimensionalities of a
c  eps - the accuracy with which the calculation will be performed;
c        should not be less than the machine precision; recommended
c        value: 1.0d-16
c  lw - the number of real *8 words provided by the user in the work 
c        array w. must be at least m+n+10 (otherwise, the subroutine
c        can not even estimate the actual memory requirements).
c
c                   output parameters:
c
c  ier - error return code
c      ier=0 means successful execution
c      ier=4 means that the amount of space lw allocated supplied
c            in the work array w is insufficient
c      ier=1024 means that the amount of space lw allocated supplied
c            in the work array w is grossly insufficient
c  u - the matrix in (1). note that the matrix u has to be dimensioned
c        u(n,m), but on exit from this subroutine only the first 
c        ncols columns of  b  are meaningful, the rest making no 
c        sense whatsoever, so that effectively on exit b is dimensioned
c        b(n,ncols)
c  v - the second matrix in (1). it is dimensioned v(m,ncols) on
c        exit. on entry the parameter ncols is unknown, 
c        so v should be dimensioned v(m,n) by the user. however,
c        on exit only the first ncols columns of v are meaningful
c  s - the array of singular values of a
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices b, v on exit
c  ltot - the total amount of space (in terms of real *8 words)
c        needed by the subroutine in the array w
c
c                   work arrays:
c
c  w - must be sufficiently long. 3*n*m+ 2*(min(m,n))**2 +2*n+2*m+100
c        are always sufficient. however, must of the time it is much 
c        less than that, and depends on the rank of a. the actual 
c        amount of space available in w is communicated to the 
c        subroutine via the parameter lw (see above), and if this 
c        amount is insufficient, the output error code ier is set
c        accordingly, and the memory requirements of the subroutine
c        are returned in the parameter ltot (see above).
c
c        . . . using gram-schmidt process with pivoting, decompose 
c              the matrix a in the form
c
c                    a=U  V^*,
c
c              with  u an orthogonal matrix of minimum rank, 
c              and v whatever it wants to be
c
        if( (lw .gt. m) .and. (lw. gt. n)) goto 1200
        ier=1024
        return
 1200 continue
c
        ifpivot=1
        call svdpiv0(a,n,m,u,v,ncols,w,eps,ifpivot)
cccc         call prin2('after first svdpiv0, rnorms=*',w,ncols)
c
c       allocate memory for the subsequent operations
c
        iw=1
        liw=m*ncols+10
c
        iu2=iw+liw
        lu2=ncols**2+10
c
        iv2=iu2+lu2
        lv2=ncols**2+10
c
        it=iv2+lv2
        lt=lv2
c
        irnorms=it+lt
        lrnorms=n
        if(m .gt. n) lrnorms=m
        lrnorms=lrnorms+10
c
        iww=irnorms+lrnorms
        lww=n
        if(m .gt. n) lww=m
        lww=lww+10
c
        ltot=iww+lww
        if(ltot .le. lw) goto 1400
        ier=4
        return
 1400 continue
c
c       conclude the construction of the singular value decomposition
c
        iffirst=0
        call svdpiv1(a,u,v,w(iw),w(it),n,m,ncols,
     1    w(irnorms),eps,w(iv2),w(iu2),w(iww),iffirst)
cccc         call prin2('after svdpiv1, w(it)=*',w(it),ncols)
c
c       copy the right singluar matrix into the array v
c
        jj=0
        do 1800 i=1,ncols
        do 1600 j=1,m        
        jj=jj+1
        v(j,i)=w(jj)
 1600 continue
 1800 continue
c
c       now, copy the singular values into the array s
c
        do 2000 i=1,ncols
        s(i)=w(it+i-1)
 2000 continue
        return
        end
c
c
c
c
c
        subroutine svdpiv1(a,u,v,w,t,n,m,ncols,rnorms,
     1    eps2,v2,u2,ww,iffirst)
        implicit real *8 (a-h,o-z)
        dimension a(n,m),u(n,1),v(m,1),w(m,1),t(1),rnorms(1),
     1      v2(1),u2(1),ww(1)
c
c        using gram-schmidt process with pivoting, decompose 
c        the matrix a in the form
c
c          a=U  V^*,
c
c        with  u an orthogonal matrix of minimum rank, 
c        and v whatever it wants to be
c
        if(iffirst .eq. 0) goto 1200
        ifpivot=1
        call svdpiv0(a,n,m,u,v,ncols,rnorms,eps2,ifpivot)
cccc         call prin2('after first svdpiv0, rnorms=*',rnorms,ncols)
 1200 continue
c
c        using gram-schmidt process without pivoting, decompose 
c        the matrix v in the form
c
c          v=w t^*,
c
c        with  w an orthogonal matrix of minimum rank, 
c        and t a triangular matrix of dimensionality ncols * ncols
c
        ifpivot=0
        call svdpiv0(v,m,ncols,w,t,ncols2,rnorms,eps2,ifpivot)
cccc         call prin2('after second svdpiv0, rnorms=*',rnorms,ncols)
c
c        use jacobi rotations to sv-decompose the matrix t
c
        ifuv=1
        call d2mudv(t,ncols,u2,v2,eps2,ifuv,numrows)
cccc        call prin2('after d2mudv, t=*',t,ncols)
        call prinf('after d2mudv, numrows=*',numrows,1)
c
c       multiply the matrix u from the right by u2, and the 
c       matrix w from the right by v2^*
c
        call svdfinm(u,u2,w,v2,n,m,ncols,ww)
        return
        end
c
c 
c
c
c
        subroutine svdfinm(u,u2,w,v2,n,m,ncols,ww)
        implicit real *8 (a-h,o-z)
        dimension u(n,ncols),w(m,ncols),u2(ncols,ncols),
     1      v2(ncols,ncols),ww(1)
c
c       multiply u from the right by u2
c
        do 1600 i=1,n
        do 1400 j=1,ncols
c
        d=0
        do 1200 k=1,ncols
        d=d+u(i,k)*u2(k,j)
 1200 continue
        ww(j)=d
 1400 continue
        do 1500 j=1,ncols
        u(i,j)=ww(j)
 1500 continue
 1600 continue
c
c         multiply w from the right by the adjoint of v2
c
        do 2600 i=1,m
        do 2400 j=1,ncols
c
        d=0
        do 2200 k=1,ncols
        d=d+w(i,k)*v2(j,k)
 2200 continue
        ww(j)=d
 2400 continue
        do 2500 j=1,ncols
        w(i,j)=ww(j)
 2500 continue
 2600 continue
c    
        return
        end
c
c
c
c
c
        subroutine svdpiv0(a,n,m,b,v,ncols,rnorms,eps,ifpivot)
        implicit real *8 (a-h,o-z)
        dimension a(n,m),b(n,m),v(m,1),rnorms(1)
c
c       this matrix applies the compressing gram-schmidt 
c       process to the matrix a, obtaining its decomposition
c       in the form
c
c             a=b v^T,                                        (1)
c
c       with the matrices b, v having dimensionalities
c       b(n,ncols), v(m,ncols), respectively. the reason for the
c       existence of this subroutine is the hope that the 
c       dimensionality ncols, determined by this subroutine,
c       is comsiderably lower than either n or m
c      
c                     input parameters: 
c
c  a - the matrix to be decomposed (note that a is NOT damaged
c       by this subroutine in any way)
c  n,m - the dimensionalities of a
c  eps - the accuracy to which the decomposition (1) will be computed.
c
c                     output parameters:
c
c  b - the matrix in (1). note that the matrix b has to be dimensioned
c        b(n,m), but on exit from this subroutine only the first 
c        ncols columns of  b  are meaningful, the rest making no 
c        sense whatsoever, so that effectively on exit b is dimensioned
c        b(n,ncols)
c  v - the second matrix in (1). it is dimensioned v(m,ncols) on
c        exit. on entry the parameter ncols is unknown, 
c        so v should be dimensioned v(m,n) by the user. however,
c        on exit only the first ncols columns of v are meaningful
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrices b, v on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c
c        . . . copy the user-supplied matrix a into b
c
        do 1400 i=1,m
        do 1200 j=1,n
        b(j,i)=a(j,i)
 1200 continue
 1400 continue
c
c        apply the gram-schmidt proces (with pivoting) to b
c
         call svdgrm(b,n,m,rnorms,eps,ncols,ifpivot)
c
c        project the original matrix on the obtained orthogonal
c        columns
c
        do 2400 i=1,ncols
        do 2200 j=1,m
c
        call svdscapr(a(1,j),b(1,i),n,prod)
        v(j,i)=prod
 2200 continue
 2400 continue
        return
        end
c
c
c
c
c
        subroutine svdgrm(b,n,m,rnorms,eps,ncols,ifpivot)
        implicit real *8 (a-h,o-z)
        dimension b(n,m),rnorms(1)
c
c       this subroutine applies a pivoted double gram-schmidt 
c       procedure to the matrix b, resulting in a collection of
c       orthogonal vectors. the number of these vectors is 
c       equal to the numerical rank of b, given the precision
c       of computations eps, which is supplied by the user.
c
c                    input paramneters:
c
c  b - the matrix to be gram-schmidt'ed. it is destroyed by this
c       subroutine
c  n,m - dimensionalities of the matrix b
c  eps - the machine accuracy 
c  
c                     output parameters:
c
c  b - the matrix of gram-schmidt vectors of the matrix a. note 
c        that on exit from this subroutine only the first 
c        ncols columns of  b  are meaningful, the rest making no 
c        sense whatsoever, so that effectively on exit b is dimensioned
c        b(n,ncols)
c  ncols - the rank of the matrix a to the precision eps. also the
c        second dimension of the matrix b on exit
c  rnorms - the normalizing factors in the gram-schmidt process.
c        only the first ncols of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        real *8 elements long.
c
c       . . . prepare the array of values of norms 
c             of columns
c
        done=1
        dtot=0
        do 1400 i=1,m
        d=0
        do 1200 j=1,n
        d=d+b(j,i)**2
 1200 continue
        rnorms(i)=d
        dtot=dtot+d
 1400 continue
c
c       . . . conduct gram-schmidt iterations
c
        thresh=dtot*eps**2
        do 4000 i=1,m
c
c       find the pivot  
c
         if(ifpivot .eq. 0) goto 2700
        ipivot=i
        rn=rnorms(i)
c
        do 2200 j=i+1,m
        if(rnorms(j) .le. rn) goto 2200
        rn=rnorms(j)
        ipivot=j
 2200 continue
 2400 continue
c
c       put the column number ipivot in the i-th place
c
        do 2600 j=1,n
        d=b(j,i)
        b(j,i)=b(j,ipivot)
        b(j,ipivot)=d
 2600 continue
c
        d=rnorms(ipivot)
        rnorms(ipivot)=rnorms(i)
        rnorms(i)=d
 2700 continue
c
c       orthogonalize the i-th column to all preceeding ones
c
        if(i .eq. 1) goto 2790
        do 2780 j=1,i-1
c
        call svdscapr(b(1,i),b(1,j),n,d)
c
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*d
 2770 continue
 2780 continue
 2790 continue
c
c       normalize the i-th column
c
        call svdscapr(b(1,i),b(1,i),n,d)
c
c       if the norm of the remaining part of the matrix 
c       is sufficiently small - terminate the process
c
        if( (ifpivot .eq. 1) .and.
     1      (d .lt. thresh) ) return
        ncols=i
c
        d=done/dsqrt(d)
        do 2800 j=1,n
        b(j,i)=b(j,i)*d
 2800 continue
c
c        orthogonalize everything else to it 
c
        do 3200 j=i+1,m
c
        if(rnorms(j) .lt. thresh) goto 3200
c
        call svdscapr(b(1,i),b(1,j),n,d)
c
         rrn=0
        do 3000 l=1,n
        b(l,j)=b(l,j)-b(l,i)*d
        rrn=rrn+b(l,j)**2
 3000 continue
        rnorms(j)=dsqrt(rrn)
 3200 continue 
 3400 continue
 4000 continue

 4200 continue


        return
        end
c
c
c
c
c
        subroutine svdscapr(x,y,n,prod)
        implicit real *8 (a-h,o-z)
        dimension x(1),y(1)
c
        prod=0
        do 1200 i=1,n
        prod=prod+x(i)*y(i)
 1200 continue
        return
        end
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c            this is the end of the debuging code and the beginning
c            of the actual SVD routines
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine d2mudv(a,n,u,v,eps,ifuv,numrows)
        implicit real *8 (a-h,o-z)
        dimension a(n,n),u0(2,2),v0(2,2),u(n,n),v(n,n),
     1      b(2,2),rlam(2)
c
c        this subroutine uses the classical jacobi iteration
c        to construct the singular value decomposition of a 
c        general real matrix, of the form
c
c          a = u s v,                                        (1)
c         
c        with u,v orthogonal matrices, and s diagonal, with
c        the diagonal elements monotonically decreasing.
c
c                     input parameters:        
c
c  a - the n * n matrix to be sv-decomposed, is destroyed 
c        by the subroutine
c  n - the dimensionality of a
c  eps - the accuracy to which the computations will be performed
c  ifuv - tells the subroutine whether the singular vectors are 
c        to be computed, or only singular values
c  
c                     output parameters:
c
c  a - the first column of a contains the (ordered) singular values 
c  u,v - orthogonal matrices in (1)
c  numrows - the number of jacobi iterations (i.e. single rotations)
c        that the algorithm has taken
c
c       . . . one sweep after another, conduct jacobi 
c             iterations
c
        maxsweep=10
        if(ifuv .eq. 0) goto 1150
        zero=0
        do 1100 i=1,n
        do 1050 j=1,n
        u(j,i)=zero
        v(j,i)=zero
 1050 continue
        u(i,i)=1
        v(i,i)=1
 1100 continue
 1150 continue
        numrows=0
        thresh=dsqrt(eps)
        thresh=dsqrt(thresh)
        thresh=dsqrt(thresh)
        do 3400 iijjkk=1,4
c
        do 3000 ijk=1,maxsweep
        ifany=0
c
        do 2000 i=1,n
        do 1800 j=1,n
         if(j .eq. i) goto 1800
         if( (dabs(a(i,j)) .lt. thresh) .and. 
     1       (dabs(a(j,i)) .lt. thresh) ) goto 1800
        ifany=1
        numrows=numrows+1
c
c       construct the svd of the (i,j)-th 2 * 2 submatrix
c
        b(1,1)=a(i,i)
        b(1,2)=a(i,j)
        b(2,1)=a(j,i)
        b(2,2)=a(j,j)
        call oneuv(b,rlam,u0,v0)
c
c       transform the big matrix
c
c       . . . the rows
c
        do 1200 l=1,n
        di=u0(1,1)*a(i,l)-u0(1,2)*a(j,l)
        dj=-u0(2,1)*a(i,l)+u0(2,2)*a(j,l)
        a(i,l)=di
        a(j,l)=dj
 1200 continue
c
c       . . . and the columns
c
        do 1400 l=1,n
        di=v0(1,1)*a(l,i)+v0(1,2)*a(l,j)
        dj=v0(2,1)*a(l,i)+v0(2,2)*a(l,j)
        a(l,i)=di
        a(l,j)=dj
 1400 continue
c
c
c       if the user so requested, adjust the matrices u,v
c
        if(ifuv .eq. 0) goto 1800
c
        do 1600 l=1,n
c
c       the matrix u
c
        di=u0(1,1)*u(l,i)-u0(1,2)*u(l,j)
        dj=-u0(2,1)*u(l,i)+u0(2,2)*u(l,j)
        u(l,i)=di
        u(l,j)=dj
c
c       . . . and the columns
c
        di=v0(1,1)*v(i,l)+v0(1,2)*v(j,l)
        dj=v0(2,1)*v(i,l)+v0(2,2)*v(j,l)
        v(i,l)=di
        v(j,l)=dj
 1600 continue
c
 1800 continue
 2000 continue
c
      if(ifany .eq. 0) goto 3200
 3000 continue
 3200 continue
        thresh=thresh**2
 3400 continue
c
c        now, reorder the singular values of a to put them
c        in increasing order
c
        do 4200 i=1,n
        a(i,1)=a(i,i)
 4200 continue
c
       call d2msvsrt(a(1,1),n,u,v,a(1,2),ifuv)
        return
        end
c
c
c
c
c
        subroutine d2msvsrt(s,n,u,v,iw,ifuv)
        implicit real *8 (a-h,o-z)
        dimension u(n,n),v(n,n),s(1),iw(1)
c
c        if some of the alleged singular values 
c        are negative - change their sign, and that
c        of the corresponding left eigenvectors
c
        do 1100 i=1,n
        if(s(i) .ge. 0) goto 1100
        s(i)=-s(i)
        if(ifuv .eq. 0) goto 1100
        do 1050 j=1,n
        v(i,j)=-v(i,j)
 1050 continue
 1100 continue

c        using the bubble sort, reorder the elements 
c        of array s
c
        do 1200 i=1,n
        iw(i)=i
 1200 continue
c
        do 2000 i=1,n-1
c
        ifmoved=0
        do 1600 j=1,n-i
        if(s(j). ge. s(j+1)) goto 1600
        jj=iw(j)
        iw(j)=iw(j+1)
        iw(j+1)=jj
c
        d=s(j)
        s(j)=s(j+1)
        s(j+1)=d
 1600 continue
 2000 continue
        if(ifuv .eq. 0) return
c
c       now, reorder the rows of u and columns of v
c
        do 2400 i=1,n-1
        ii=iw(i)
c
c       exchange the rows number i and ii in the matrix u
c
        do 2200 j=1,n
        d=v(i,j)
        v(i,j)=v(ii,j)
        v(ii,j)=d
c
c       exchange columns number i and ii in the matrix v
c
        d=u(j,i)
        u(j,i)=u(j,ii)
        u(j,ii)=d
 2200 continue
c
        do 2300 j=i+1,n
        if(iw(j) .ne. i) goto 2300
        iw(j)=ii
        goto 2350
 2300 continue
 2350 continue
c
 2400 continue
        return
        end

c
c
c
c
c
        subroutine oneuv(a,rlam,u,v)
        implicit real *8 (a-h,o-z)
        dimension a(2,2),u(2,2),v(2,2),rlam(2),
     1      b(2,2),w(2,2),vstar(2,2),z(2,2),rlams(2,2)
        data eps/1.0d-10/,eps2/1.0d-1/,done/1.0d0/
c
c       this subroutine produces the singular value decomposition
c       of a 2 * 2 real matrix a, so that on exit,
c
c       a = u d v,
c
c       with u and v orthogonal, and d a diagonal matrix with the 
c       vector rlam=(rlam(1),rlam(2)) on the diagonal.
c
c       . . . simmetrize a
c
        den=a(2,2)+a(1,1)
        rnum=a(2,1)-a(1,2)
c
c       if denominator is not too small, use exact formulae
c
        if(dabs(den) .lt. eps*dabs(rnum)) goto 1200
        tgphi=rnum/den
c
        if(dabs(tgphi) .lt. eps2) goto 1100
c
        dd=done/(done+tgphi**2)
        alpha=dsqrt(dd)
        beta=dsqrt(done-dd)
        if(tgphi .lt. 0) beta=-beta
        goto 1400
 1100 continue
c
        dd=tgphi**2/(1+tgphi**2)
c
        beta=tgphi-tgphi**3/2
        alpha=dsqrt(done-beta**2)
c
cccc        beta=dsqrt(dd)
cccc        alpha=dsqrt(done-dd)
cccc        if(tgphi .lt. 0) beta=-beta
c        
cccc        phi=datan(tgphi)
cccc        alpha=dcos(phi)
cccc        beta=dsin(phi)
        goto 1400
 1200 continue
c
c       denominator is too small. use taylor series
c
        alpha=den/rnum
        beta=1
 1400 continue
c
c      calculate the simmetrizing matrix and the simmetric one
c
        w(1,1)=alpha
        w(1,2)=beta
        w(2,1)=-beta
        w(2,2)=alpha
        call prod2(w,a,b)
c
c       now, diagonalize the simmetrized matrix
c
        den=b(2,2)-b(1,1)
        rnum=b(1,2)*2
c
c       if denominator is not too small, use exact formulae
c
        if(dabs(den) .lt. eps*dabs(rnum)) goto 1600
        tg2phi=-rnum/den
        if(dabs(tg2phi) .lt. 1.0d-3) goto 1500
c
        dd=dsqrt(4/tg2phi**2+4)
        if(tg2phi .lt. 0) tgphi=(-2/tg2phi-dd)/2
        if(tg2phi .gt. 0) tgphi=(-2/tg2phi+dd)/2
c
        phi=datan(tgphi)
        alpha=dcos(phi)
        beta=dsin(phi)
        goto 1800
 1500 continue
c
        beta=tg2phi/2
        alpha=dsqrt(done-beta**2)

cccc        phi=datan(tg2phi)/2
cccc        alpha=dcos(phi)
cccc        beta=dsin(phi)
        goto 1800
 1600 continue
c
c       denominator is too small. use taylor expansion
c
        alpha=dsqrt((1-den/rnum)/2)
        beta=dsqrt((1+den/rnum)/2)
 1800 continue
c
c       construct the diagonalizing matrix
c       and the resulting diagonal
c
        v(1,1)=alpha
        v(1,2)=beta
        v(2,1)=-beta
        v(2,2)=alpha
c
        vstar(1,1)=v(1,1)
        vstar(1,2)=v(2,1)      
        vstar(2,1)=v(1,2)
        vstar(2,2)=v(2,2)       
c       
        call prod2(v,w,u)
c
c       finally, compute the resulting diagonal elements
c
        call prod2(u,a,z)
        call prod2(z,vstar,rlams)
        rlam(1)=rlams(1,1)
        rlam(2)=rlams(2,2)
c
        u(1,2)=-u(1,2)
        u(2,1)=-u(2,1)
        return
        end
c
c
c
c
c
        subroutine prod2(a,b,c)
        implicit real *8 (a-h,o-z)
        dimension a(2,2),b(2,2),c(2,2)
c
        c(1,1)=a(1,1)*b(1,1)+a(1,2)*b(2,1)
        c(1,2)=a(1,1)*b(1,2)+a(1,2)*b(2,2)
        c(2,1)=a(2,1)*b(1,1)+a(2,2)*b(2,1)
        c(2,2)=a(2,1)*b(1,2)+a(2,2)*b(2,2)
        return
        end
