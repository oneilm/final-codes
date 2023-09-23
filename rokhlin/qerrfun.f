        implicit real *16(a-h,o-z)
        dimension errs(300),
     1     coefs(1000),f2(1000),x0(1000)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
         PRINT *, 'ENTER a'
         READ *,a
         CALL PRINQ('a=*',a,1 )
c 
c       check the fast error function routine against the slow one
c 
        n=240
        call erf0(a,n,erf00)
        call prinq('after erf0, erf00=*',erf00,1)
c 
        call qerrfun(a,erfun)
  
        call prinq('and after qerrfun, erfun=*',erfun,1)
        call prinq('and the ratio-1=*',erfun/erf00-1,1)
c 
c       perform massive testing
c 
        nn=200
        h=10
        h=h/nn
        errmax=0
        do 1200 i=1,nn
        t=(i-1)*h+h/2
        call erf0(t,n,erf00)
c 
        call qerrfun(t,erfun)
c 
        errs(i)=erfun/erf00-1
        call prinq('errs(i)=*',errs(i),1)
        if(errmax .lt. qabs(errs(i)) ) errmax=qabs(errs(i))
 1200 continue
        call prinq('and errs=*',errs,200)
        call prinq('and errmax=*',errmax,1)
  
  
  
        stop
        end
  
c 
c 
c 
c 
c 
       subroutine erf0(x,n,erf00)
       implicit real *16 (a-h,o-z)
        save
        dimension t(1000),w(1000),tt(1000),ww(1000),f(1000)
        data n7/-9654236/
c 
c        construct the chebychev nodes and weights
c        on the interval [-1,1]
c 
        if(n7 .eq.n) goto 1100
        n7=n
        itype=1
        call qchebexp(itype,n,t,u,v,w)
 1100 continue
c 
c       construct the chebychev nodes on the interval 0,x
c 
        do 1200 i=1,n
        tt(i)=(t(i)+1)/2*x
        ww(i)=w(i)/2*x
 1200 continue
c 
c        evaluate exp(-x**2) at all these nodes
c 
        do 1400 i=1,n
        f(i)=qexp(-tt(i)**2)
 1400 continue
cccc        call prinq('and f=*',f,n)
c 
c        . . . integrate
c 
        erf00=0
        do 1600 i=1,n
        erf00=erf00+ww(i)*f(i)
 1600 continue
        done=1
        pi=qatan(done)*4
        erf00=erf00*2/qsqrt(pi)
  
  
        return
        end
c 
c 
c 
c 
c 
       subroutine erfc0(x,n,erfc00)
       implicit real *16 (a-h,o-z)
        save
        dimension t(1000),w(1000),tt(1000),ww(1000),f(1000)
        data n7/-9654236/
c 
c        construct the chebychev nodes and weights
c        on the interval [-1,1]
c 
        if(n7 .eq.n) goto 1100
        n7=n
        itype=1
        call qchebexp(itype,n,t,u,v,w)
 1100 continue
c 
c       construct the chebychev nodes on the interval (x,x+20)
c 
        do 1200 i=1,n
        tt(i)=(t(i)+1)/2*20 +x
        ww(i)=w(i)/2*20
 1200 continue
c 
c        evaluate exp(-x**2) at all these nodes
c 
        do 1400 i=1,n
        f(i)=qexp(-tt(i)**2)
 1400 continue
cccc        call prinq('and f=*',f,n)
c 
c        . . . integrate
c 
        erfc00=0
        do 1600 i=1,n
        erfc00=erfc00+ww(i)*f(i)
 1600 continue
c 
        done=1
        pi=qatan(done)*4
        erfc00=erfc00*2/qsqrt(pi)
  
  
  
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
        subroutine qchebexp(itype,n,x,u,v,whts)
        implicit real *16(a-h,o-z)
        save
        dimension x(1),whts(1),u(n,n),v(n,n)
c 
c         this subroutine constructs the chebychev nodes
c         on the interval [-1,1], and the weights for the
c         corresponding order n quadrature. it also constructs
c         the matrix converting the coefficients
c         of a chebychev expansion into its values at the n
c         chebychev nodes. no attempt has been made to
c         make this code efficient, but its speed is normally
c         sufficient, and it is mercifully short.
c 
c                 input parameters:
c 
c  itype - the type of the calculation to be performed
c          itype=0 means that only the chebychev nodes are
c                  to be constructed.
c          itype=1 means that only the nodes and the weights
c                  are to be constructed
c          itype=2 means that the nodes, the weights, and
c                  the matrices u, v are to be constructed
c  n - the number of chebychev nodes and weights to be generated
c 
c                 output parameters:
c 
c  x - the order n chebychev nodes - computed independently
c          of the value of itype.
c  u - the n*n matrix converting the  values at of a polynomial of order
c         n-1 at n chebychev nodes into the coefficients of its
c         chebychev expansion - computed only in itype=2
c  v - the n*n matrix converting the coefficients
c         of an n-term chebychev expansion into its values at
c         n chebychev nodes (note that v is the inverse of u)
c          - computed only in itype=2
c  whts - the corresponding quadrature weights - computed only
c         if itype .ge. 1
c 
c       . . . construct the chebychev nodes on the interval [-1,1]
c 
        ZERO=0
        DONE=1
        pi=qatan(done)*4
        h=pi/(2*n)
        call foolcoi(ier)
ccc        do 1200 i=n,1,-1
        do 1200 i=1,n
        t=(2*i-1)*h
cccc        x(n-i+1)=qcos(t)
        call foolcosi(t,x(n-i+1),extsin)
1200  CONTINUE
ccc        call prinq('chebychev nodes as constructed*',x,n)
c 
        if(itype .eq. 0) return
c 
c        construct the weights of the quadrature
c        formula based on the chebychev nodes,
c        and also the matrix of the chebychev transform
c 
c        . . . construct the first two rows of the matrix
c 
         if(itype .le. 1) goto 1350
        do 1300 i=1,n
        u(1,i)=1
        u(2,i)=x(i)
 1300 continue
 1350 continue
c 
c       construct all quadrature weights and the rest of the rows
c 
        do 2000 i=1,n
c 
c       construct the weight for the i-th node
c 
        Tjm2=1
        Tjm1=x(i)
        whts(i)=2
c 
        ic=-1
        do 1400 j=2,n-1
c 
c       calculate the T_j(x(i))
c 
        Tj=2*x(i)*Tjm1-Tjm2
c 
        if(itype .eq. 2) u(j+1,i)=tj
c 
        tjm2=tjm1
        tjm1=tj
c 
c       calculate the contribution of this power to the
c       weight
c 
        ic=-ic
        if(ic .lt. 0) goto 1400
        rint=-2*(done/(j+1)-done/(j-1))
        whts(i)=whts(i)-rint*tj
ccc        whts(i)=whts(i)+rint*tj
 1400 continue
        whts(i)=whts(i)/n
 2000 continue
           if(itype .ne. 2) return
c 
c        now, normalize the matrix of the chebychev transform
c 
        do 3000 i=1,n
c 
        d=0
        do 2200 j=1,n
        d=d+u(i,j)**2
 2200 continue
        d=done/qsqrt(d)
        do 2400 j=1,n
        u(i,j)=u(i,j)*d
 2400 continue
 3000 continue
c 
c        now, rescale the matrix
c 
        ddd=2
        ddd=qsqrt(ddd)
        dd=n
        dd=done/qsqrt(dd/2)
        do 3400 i=1,n
        do 3200 j=1,n
        u(j,i)=u(j,i)*dd
 3200 continue
        u(1,i)=u(1,i)/ddd
 3400 continue
c 
c        finally, construct the matrix v, converting the values at the
c        chebychev nodes into the coefficients of the chebychev
c        expansion
c 
        dd=n
        dd=dd/2
        do 4000 i=1,n
        do 3800 j=1,n
        v(j,i)=u(i,j)*dd
 3800 continue
 4000 continue
c 
        do 4200 i=1,n
        v(i,1)=v(i,1)*2
 4200 continue
        return
        end
c 
c 
c 
c 
c 
      SUBROUTINE QCHFUNDE(X,VAL,der,TEXP,N)
C 
C     This subroutine computes the value and the derivative
c     of a chebychev expansion with coefficients TEXP
C     at point X in interval [-1,1]
C 
c                input parameters:
c 
C     X = evaluation point
C     TEXP = expansion coefficients
C     N  = order of expansion (0 ,..., N-1)
c 
c                output parameters:
c 
C     VAL = computed value
C     der = computed value of the derivative
C 
      IMPLICIT REAL *16(A-H,O-Z)
        save
      REAL *16TEXP(1)
C 
        done=1
        tjm2=1
        tjm1=x
        derjm2=0
        derjm1=1
c 
        val=texp(1)*tjm2+texp(2)*tjm1
        der=texp(2)
c 
      DO 600 J = 2,N
c 
        tj=2*x*tjm1-tjm2
        val=val+texp(j+1)*tj
c 
        derj=2*tjm1+2*x*derjm1-derjm2
        der=der+texp(j+1)*derj
c 
        tjm2=tjm1
        tjm1=tj
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
        subroutine foolcosi(x,extcos,extsin)
        implicit real *16 (a-h,o-z)
        save
        dimension factin(0:30)
c 
c        compute, slowly but correctly, the sine and cosine
c        of x with extended precision
c 
        done=1
        two=2
        half=done/2
        n=20
        small=0.01
c 
        nhalf=0
        xx=x
ccccc         if(2 .ne. 3) goto 1500
        do 1200 i=1,100
        if(qabs(xx) .lt. small) goto 1400
        xx=xx*half
        nhalf=nhalf+1
 1200 continue
 1400 continue
 1500 continue
cccc         call prinf('nhalf=*',nhalf,1)
cccc         call prinq('xx=*',xx,1)
c 
c       calculate sine and cosine of xx
c 
        extcos=done
        extsin=0
c 
        d1=xx
        d2=xx**2
        dd=d2
        do 1600 i=1,n
        extcos=extcos+d2*factin(2*i)
        extsin=extsin+d1*factin(2*i-1)
        d1=d1*dd
        d2=d2*dd
 1600 continue
        extsin=-extsin
cccc          if(2 .ne. 3) return
c 
c       now, square the thing nhalf times
c 
        if(nhalf .eq. 0) return
        do 1800 i=1,nhalf
cccc          call prinf('i=*',i,1)
cccc          call prinq('extcos=*',extcos,1)
cccc          call prinq('extsin=*',extsin,1)
        d2=extcos**2-extsin**2
        d1=two*extcos*extsin
        extsin=d1
        extcos=d2
 1800 continue
        return
c 
c 
c 
c 
        entry foolcoi(ier)
c 
c       this is initialization entry point
c 
        coe=1
        factin(0)=1
        do 2200 i=1,22
        coe=-coe
        d=i
        factin(i)=factin(i-1)/d*coe
 2200 continue
cccc         call prinq('factin as created*',factin(0),30)
        return
        end
c 
c 
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c       this is the end of the debugging code and the beginning
c       of the error function routines
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c 
c 
        subroutine qerrfun(x7,erfun)
        implicit real *16 (a-h,o-z)
        save
        dimension coefs(18),coefs2(13),bufcoe(2)
        equivalence (coefs(18),bufcoe(1)),(coefs2(1),bufcoe(2))
c 
c        this subroutine evaluates the error function defined
c        by the formula
c 
c        erf(x)=\int_0^x exp(-t**2) dt
c 
c       with extended precision accuracy.
c 
c 
c                      input parameters:
c 
c  x7 - the argument
c                      output parameters:
c 
c  erfun - error function of x
c 
        data coefs/
     1     0.46963291789115009827412922472047D+00,
     2     0.42582445804381043569204735290848D+00,
     3     -.49552626796204340403576830798593D-01,
     4     -.44929348876838274955800124201535D-02,
     5     0.12919410465849695349422476112098D-02,
     6     0.18363892921493962704169788672122D-04,
     7     -.22111147040995262915385554699197D-04,
     8     0.52333748523425713467369319779734D-06,
     9     0.27818478883353788538253099109576D-06,
     a     -.14115809274881311456031659199125D-07,
     1     -.27257129633056169998453886261514D-08,
     2     0.20634390487207062940639851585784D-09,
     3     0.21427399199678536792416443571655D-10,
     4     -.22299025553935820457965847853047D-11,
     5     -.13625007465069828057004944758513D-12,
     6     0.19514401092229309192942348364070D-13,
     7     0.68562716923170458929126778282939D-15,
     8     -.14450649286969997803041431643477D-15/
         data coefs2/
     9     -.24593530646054176842074301775172D-17,
     a     0.92959956122049262809396442021631D-18,
     1     0.28743670971555272132905892954155D-20,
     2     -.52871778042148820515979592731087D-20,
     3     0.44807385096647751700661175657819D-22,
     4     0.26922569397181593320864331505555D-22,
     5     -.50472040668074925808180499798587D-24,
     6     -.12386920454438148903803863927236D-24,
     7     0.35241362346072537624979164129445D-26,
     8     0.51841341153297714329403112468134D-27,
     9     -.19798231818357131118309229653530D-28,
     a     -.19862928018160851717691808512573D-29,
     1     0.64298976752690731987886063592519D-31/
c 
        x=qabs(x7)
c 
c        if x .gt. 0.1 - evaluate erf via erfc
c 
        if(x .lt. 1.0q0) goto 1400
        call qerrfunc(x,erfunc)
        erfun=1-erfunc
        if(x7 .lt. 0) erfun=-erfun
        return
 1400 continue
c 
c        x .lt. 0.1. use chebychev series to evaluate erfun
c 
        t=x*2-1
        n=30
        call qerrexev(coefs,n,t,erfun)
        if(x7 .lt. 0) erfun=-erfun
        return
        end
c 
c 
c 
c 
c 
        subroutine qerrfunc(x7,erfunc)
        implicit real * 16 (a-h,o-z)
        save
        dimension c02n1(10),c02n2(10),c02n3(10),c02n4(8),
     1     c24n1(10),c24n2(10),c24n3(12),
     2     c46n1(10),c46n2(10),c46n3(8),
     3     c68n1(10),c68n2(10),c68n3(6),
     4      c810n1(10),c810n2(15),
     5     c02(30),c24(30),c46(30),c68(30),c810(30)
        equivalence
     1   (c02(1),c02n1(1)),(c02(11),c02n2(1)),(c02(21),c02n3(1)),
     2   (c24(1),c24n1(1)),(c24(11),c24n2(1)),(c24(21),c24n3(1)),
     3   (c46(1),c46n1(1)),(c46(11),c46n2(1)),(c46(21),c46n3(1)),
     4   (c68(1),c68n1(1)),(c68(11),c68n2(1)),(c68(21),c68n3(1)),
     5   (c810(1),c810n1(1)),(c810(11),c810n2(1)),
     6   (c02(31),c02n4(1))
c 
c        this subroutine evaluates the second error function defined
c        by the formula
c 
c        erfc(x)=\int_x^{\infty}  exp(-t**2) dt = 1 - erf(x)
c 
c       with extended precision accuracy
c 
c 
c                      input parameters:
c 
c  x7 - the argument
c                      output parameters:
c 
c  erfunc - second error function of x
  
c 
c   coefficients for the interval [0, 2]
c 
        data c02n1/
     1     0.52136087049960619781849508312991Q+00,
     2     -.34477074547252352449342439419592Q+00,
     3     0.99759779268200963479718559480923Q-01,
     4     -.26063665256478758984350198490020Q-01,
     5     0.62684634943190106999675369238913Q-02,
     6     -.14061223858870650278211581416723Q-02,
     7     0.29698557246548996692018339283689Q-03,
     8     -.59485553901067240043910010542263Q-04,
     9     0.11362960761135521370839946482502Q-04,
     a     -.20794078634147761077310996754027Q-05/
c 
        data c02n2/
  
     1     0.36590821201669711438642556301882Q-06,
     2     -.62107261363621591843172860865957Q-07,
     3     0.10195316699270051335580482106724Q-07,
     4     -.16223227802469610796062165768710Q-08,
     5     0.25073406746662080536115811286194Q-09,
     6     -.37703813593742712974168897019594Q-10,
     7     0.55248768982320522072097020012181Q-11,
     8     -.78999129790218407822939938977739Q-12,
     9     0.11036213151407694552179330257083Q-12,
     a     -.15079931632503800581334842519161Q-13/
c 
        data c02n3/
     1     0.20174311859216235287859448030821Q-14,
     2     -.26449521072902580164793348211251Q-15,
     3     0.34011233662469910562850493803972Q-16,
     4     -.42928773664251293839650309006466Q-17,
     5     0.53223804776089466438685653147602Q-18,
     6     -.64860644065140435856158116762429Q-19,
     7     0.77739252371807837876370560162840Q-20,
     8     -.91691812750698378690597135480207Q-21,
     9     0.10648355585514911368237797170254Q-21,
     a     -.12181814947819617722300168602698Q-22/
c 
        data c02n4/
     1     0.13734755945084779867462110216846Q-23,
     2     -.15268575573384197852242343081360Q-24,
     3     0.16742427117522907125140060745447Q-25,
     4     -.18114692950611147766479753260373Q-26,
     5     0.19349079790342965313845605002607Q-27,
     6     -.20546472152712031825029589284034Q-28,
     7     0.23436671585903654146398399582971Q-29,
     8     -.51064772345483477151570712662071Q-30/
c 
c   coefficients for the interval [2, 4]
        data c24n1/
     1     0.18742956902622409306925144479872Q+00,
     2     -.57946528836690967785413378444440Q-01,
     3     0.85952473873965326996308911399058Q-02,
     4     -.12284613640242593736724370562102Q-02,
     5     0.16974366686180486953367670784488Q-03,
     6     -.22738329337044562795469514814350Q-04,
     7     0.29598443788048780996722287708018Q-05,
     8     -.37513575581379711371503352440765Q-06,
     9     0.46372828761738247083074639709199Q-07,
     a     -.55994553925640879241202334507183Q-08/
c 
        data c24n2/
     1     0.66131608966799159755426938686925Q-09,
     2     -.76482715221998167746578867813974Q-10,
     3     0.86709043258765937469529030691067Q-11,
     4     -.96454562743821085144258450212591Q-12,
     5     0.10536827959360625660700528989713Q-12,
     6     -.11312630906760667269865441308807Q-13,
     7     0.11945180669174576517515458382592Q-14,
     8     -.12413107525804078010209734851108Q-15,
     9     0.12702463231854429281216324752301Q-16,
     a     -.12807258736639517232877482476199Q-17/
c 
        data c24n3/
     1     0.12729426443963546096673617369866Q-18,
     2     -.12478254325373568523046115484300Q-19,
     3     0.12069382415080905790174026421391Q-20,
     4     -.11523481357557426836576718890519Q-21,
     5     0.10864746405132768883292074304440Q-22,
     6     -.10119339333958407890345558606455Q-23,
     7     0.93138978055373967457176232654886Q-25,
     8     -.84741890625381197262746465981092Q-26,
     9     0.76236749091460841177609634990828Q-27,
     a     -.67781247299713562765608554360546Q-28,
     1     0.58847289741147701455688797529953Q-29,
     2     -.43217491793857309410573595035914Q-30/
c 
c   coefficients for the interval [4, 6]
c 
        data c46n1/
     1     0.11277819743166245319088851431478Q+00,
     2     -.21913497155272348854403231709269Q-01,
     3     0.20915397845687494704395083143931Q-02,
     4     -.19628624819122982453937138730095Q-03,
     5     0.18126654238718201380392151756897Q-04,
     6     -.16483628246572165758417591008001Q-05,
     7     0.14769699327777976107481836575099Q-06,
     8     -.13047490253442587780462366837498Q-07,
     9     0.11369753471043668813639541094894Q-08,
     a     -.97781802266404980741610121580943Q-10/
c 
        data c46n2/
     1     0.83032123480629757872312378628952Q-11,
     2     -.69646508425907134583559420164119Q-12,
     3     0.57728321648271661335653788247210Q-13,
     4     -.47301485583698741521430457595131Q-14,
     5     0.38327149744119831383914884472204Q-15,
     6     -.30720233149987745126274517435573Q-16,
     7     0.24364677192448132179513790068467Q-17,
     8     -.19126665070736228731222529219220Q-18,
     9     0.14865432561241565654130238348592Q-19,
     a     -.11441591030872881659849514692371Q-20/
c 
        data c46n3/
     1     0.87230998571718318806405344232366Q-22,
     2     -.65891651936999413334501938612018Q-23,
     3     0.49324158413158222073023864676717Q-24,
     4     -.36597300135688308474178347999623Q-25,
     5     0.26920767745242577505300186270589Q-26,
     6     -.19637445933108154996248853884989Q-27,
     7     0.14216229691161037924544550733445Q-28,
     8     -.10136635733468428578316773906482Q-29/
c 
c   coefficients for the interval [6, 8]
c 
        data c68n1/
     1     0.80586726359943425276197738287254Q-01,
     2     -.11340871814321887442002116443655Q-01,
     3     0.79038930590881203560170174393107Q-03,
     4     -.54574581010721115202842178559706Q-04,
     5     0.37342298952738347988298869910226Q-05,
     6     -.25326406802593954485737513767179Q-06,
     7     0.17029548888897999275400134926829Q-07,
     8     -.11354810726045605056112156438796Q-08,
     9     0.75091260369026600390440269335786Q-10,
     a     -.49262100700803286822366796153529Q-11/
c 
        data c68n2/
     1     0.32064751085119214716490647595778Q-12,
     2     -.20711328658873924858613392100474Q-13,
     3     0.13277721023908436484100441115886Q-14,
     4     -.84497101884086011030267767173465Q-16,
     5     0.53385959541556925857276665429843Q-17,
     6     -.33491977405720461950783951234504Q-18,
     7     0.20866134081155284916762621651560Q-19,
     8     -.12911813267013401127435041289982Q-20,
     9     0.79365332146565188043514089609902Q-22,
     a     -.48464641692687441074256035233402Q-23/
c 
        data c68n3/
     1     0.29404919865774879091294502596371Q-24,
     2     -.17728209305496042549085176793617Q-25,
     3     0.10622256869982957270388264515729Q-26,
     4     -.63281696042167086852717865443297Q-28,
     5     0.37634510526886755284206278477043Q-29,
     6     -.23276123507863826776390035921707Q-30/
c 
c   coefficients for the interval [8,10]
c 
        data c810n1/
     1     0.62684290243453396610825314277179Q-01,
     2     -.69014792436725362656233676014409Q-02,
     3     0.37767451908518504573804627058453Q-03,
     4     -.20547526563664360267957542649865Q-04,
     5     0.11115026063587423010624134024677Q-05,
     6     -.59787671676642149912916275307666Q-07,
     7     0.31981785361812047436886673174224Q-08,
     8     -.17014609183903509133469574199210Q-09,
     9     0.90033981117773496809677415990439Q-11,
     a     -.47390564596554209182714435555536Q-12/
c 
        data c810n2/
     1     0.24814917515262421373897053355326Q-13,
     2     -.12927149377630089157277771898625Q-14,
     3     0.67002972077760909937376758563788Q-16,
     4     -.34555574320072334655326037049359Q-17,
     5     0.17733943812920036876032543058533Q-18,
     6     -.90570174500582326678144704124431Q-20,
     7     0.46034811426248772172003339213424Q-21,
     8     -.23288236840799981234004091731405Q-22,
     9     0.11726351904475800826720593142939Q-23,
     a     -.58774876038377951007211157322197Q-25,
     1     0.29325940462141413026313215565472Q-26,
     2     -.14569212683159290298634440553278Q-27,
     3     0.72274922213980504030355407245095Q-29,
     4     -.37605098717833700992926894313892Q-30,
     5     0.31974386737846878046756008801675Q-31/
c 
        x=qabs(x7)
c 
c       if x .gt. 10 - set the function to zero and exit
c 
        if(x .lt. 9.9999q0) goto 1100
        erfunc=0
        return
 1100 continue
c 
c        if x is on the interval [0,2] - act accordingly
c 
        if(x .gt. 2) goto 1200
        n=37
        t=x-1
        call qerrexev(c02,n,t,erfunc)
        erfunc=erfunc*qexp(-x**2)
        if(x7 .lt. 0) erfunc=2-erfunc
        return
 1200 continue
c 
        if(x .gt. 4) goto 1400
        n=31
        t=x-3
        call qerrexev(c24,n,t,erfunc)
        erfunc=erfunc*qexp(-x**2)
        if(x7 .lt. 0) erfunc=2-erfunc
        return
 1400 continue
c 
        if(x .gt. 6) goto 1600
        n=27
        t=x-5
        call qerrexev(c46,n,t,erfunc)
        erfunc=erfunc*qexp(-x**2)
        if(x7 .lt. 0) erfunc=2-erfunc
        return
 1600 continue
c 
        if(x .gt. 8) goto 1800
        n=25
        t=x-7
        call qerrexev(c68,n,t,erfunc)
        erfunc=erfunc*qexp(-x**2)
        if(x7 .lt. 0) erfunc=2-erfunc
        return
 1800 continue
c 
        if(x .gt. 10) goto 2000
        n=24
        t=x-9
        call qerrexev(c810,n,t,erfunc)
        erfunc=erfunc*qexp(-x**2)
        if(x7 .lt. 0) erfunc=2-erfunc
        return
 2000 continue
  
  
  
        return
        end
c 
c 
c 
c 
c 
        subroutine qerrexev(c,n,x,fun)
        implicit real *16 (a-h,o-z)
        save
        dimension c(1)
c 
c        evaluate the chebychev expansion with coefficients
c        c at the point x
c 
cccc        call prinq('in qerrexev, x=*',x,1)
cccc        call prinq('in qerrexev, c=*',c,10)
        x2=x+x
        tkm1=1
        tk=x
        fun=c(1)*tkm1
        do 1400 k=1,n
        tkp1=x2*tk-tkm1
c 
        fun=fun+c(k+1)*tk
c 
        tkm1=tk
        tk=tkp1
 1400 continue
cccc        call prinq('and fun=*',fun,1)
        return
        end
