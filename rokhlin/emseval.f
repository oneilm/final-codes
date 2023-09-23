        implicit real *8 (a-h,o-z)
        dimension ems(1000)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER n '
        READ *,n
        CALL PRINf('n=*',n,1 )
c 
        eps=1.0d-20
        x=2
c 
c       evaluate the exponential integrals via the subroutine
c 
        call emseval(x,n,ems)
  
        call prin2('after emseval, ems=*',ems,n)
c 
c       evaluate the same integrals via the adaptive integration
c 
  
        call emfun_adap(x,n,en,eps)
  
        call prin2('after emfun_adap, en=*',en,1)
        call prin2('and en-ems(n)=*',en-ems(n),1)
  
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine emfun_adap(x,n,en,eps)
        implicit real *8 (a-h,o-z)
        save
        external fevalm
c 
        a=1
        b=-log(eps/100)
        mm=24
  
        call adapgaus(ier,a,b,fevalm,n,x,mm,eps,
     1      en,maxrec,numint)
  
        return
        end
c 
c 
c 
c 
c 
        function fevalm(t,n,x)
        implicit real *8 (a-h,o-z)
c 
        save
        fevalm=exp(-x*t)/t**n
  
        return
        end
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code, and the beginning of
c        the Exponential Integral code proper
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c        This file contains three user-callable subroutines:
c        e1eval, emseval, alphaseval. Following is a brief
c        description of all three.
c 
c   e1eval - for the user-provided argument x \in (0,\infty),
c        this subroutine evaluates (to relative precision
c        13 - 14 digits when run in double precision, to
c        20-21 digits when run in extended precision) the
c        exponential integral
c 
c 
c        E_1(x) = \int_x^{\infty} e^{-t}/t dt               (1)
c 
c   emseval - for the user-provided argument x \in (0,\infty),
c        this subroutine evaluates (to relative precision
c        13 - 14 digits when run in double precision, to
c        20-21 digits when run in extended precision) the
c        exponential integral s
c 
c        E_m(x) = \int_1^{\infty} e^{-t*x}/t^m dt,          (2)
c 
c        with m=1,2,...,n.
c 
c   alphaseval - for the user-provided argument x \in (0,\infty),
c        this subroutine evaluates (to the machine precision)
c        the exponential integrals
c 
c 
c        alpha_m(x) = \int_1^{\infty} e^{-t*x}*t^m dt,      (3)
c 
c        with m=0,1,2,...,n.
c 
c 
c   emsexp_eval - evaluates the expansion of the form
c 
c        f=\sum_{k=1}^n coefs(k) * \int_a^{\infty} exp(-p*x) dx   (1)
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        subroutine alphaseval(x,n,alphas)
        implicit real *8 (a-h,o-z)
        save
        dimension alphas(1)
c 
c        For the user-provided argument x \in (0,\infty),
c        this subroutine evaluates (to the machine precision)
c        the exponential integrals
c 
c 
c        alpha_m(x) = \int_1^{\infty} e^{-t*x}*t^m dt,               (1)
c 
c        with m=0,1,2,...,n.
c 
c                      Input parameters:
c 
c  x - the argument in (1) above
c  n - the highest order m in (1) above
c 
c                      Output parameterds:
c 
c  alphas - the integrals alpha_m(x) in (1) above
c 
c 
c        . . . evaluate alpha_0 (x)
c 
        d=exp(-x)
        alphas(1)=d/x
c 
c       use recursion to evaluate the higher order
c       exponential integrals
c 
        do 1200 i=1,n
c 
        i1=i+1
        alphas(i1)=(d+i*alphas(i))/x
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine emsexp_eval(a,p,coefs,n,f,ems)
        implicit real *8 (a-h,o-z)
        save
        dimension ems(1000)
        complex *16 cd,coefs(1),f
c 
c        This subroutine evaluates the expansion of the
c        form
c 
c    f=\sum_{k=1}^n coefs(k) * \int_a^{\infty} exp(-p*x) /x^k dx   (1)
c 
c                    Input parameters:
c 
c  a - the left integration integration limit in (1)
c  p - the parameter in (1)
c  coefs - the coefficients in (1)
c  n - the number of terms in (1)
c 
c                    Output parameters:
c 
c  f - the sum in (1)
c 
c                    Work arrays:
c 
c  ems - must be at least n complex *16 locations long
c 
c 
c       . . . evaluate the exponential integrals
c 
        b=p*a
        call emseval(b,n,ems)
c 
        call prin2('after emseval, ems=*',ems,n)
c 
c       evaluate the sum (1)
c 
        f=0
        do 1200 i=1,n
c 
        d=ems(i)/a**(i-1)
c 
        f=f+coefs(i)*d
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine emseval(x,n,ems)
        implicit real *8 (a-h,o-z)
        save
        dimension ems(1)
c 
c        For the user-provided argument x \in (0,\infty),
c        this subroutine evaluates (to relative precision
c        13 - 14 digits when run in double precision, to
c        20-21 digits when run in extended precision) the
c        exponential integral s
c 
c 
c        E_m(x) = \int_1^{\infty} e^{-t*x}/t^m dt,               (1)
c 
c        with m=1,2,...,n.
c 
c                      Input parameters:
c 
c  x - the argument in (1) above
c  ems - the integrals E_m(x) in (1) above
c 
c 
c        . . . evaluate E_1 (x)
c 
        call e1eval(x,ems(1))
c 
c       use recursion to evaluate the higher order
c       exponential integrals
c 
        d=exp(-x)
c 
        do 1200 i=1,n-1
c 
        ems(i+1)=(d-x*ems(i))/i
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine e1eval(x,e1)
        implicit real *8 (a-h,o-z)
        save
        dimension coefs1(20),coefs2(20),coefs3(20)
c 
c        For the user-provided argument x \in (0,\infty),
c        this subroutine evaluates (to relative precision
c        13 - 14 digits when run in double precision, to
c        20-21 digits when run in extended precision) the
c        exponential integral
c 
c 
c        E_1(x) = \int_x^{\infty} e^{-t}/t dt               (1)
c 
c 
c                      Input parameters:
c 
c  x - the argument in (1) above
c  e1 - the integral E_1(x) in (1) above
c 
c 
c      . . . Coefficients for aa=0, bb=1/15
c 
        data coefs1/
     1  0.96954526666896356766D+00,-.29562457063550202274D-01,
     2  0.85495670053301315140D-03,-.35331370699233228427D-04,
     3  0.18613050094651754209D-05,-.11755179722996717467D-06,
     4  0.85668866793235894309D-08,-.70205510624291400882D-09,
     5  0.63504494788988946650D-10,-.62527170191982832598D-11,
     6  0.66293491405666360904D-12,-.75037307569501107324D-13,
     7  0.90044000800329928498D-14,-.11389284240928611247D-14,
     8  0.15111213380271024670D-15,-.20944862985693284537D-16,
     9  0.30220206234056770979D-17,-.45250304036640517271D-18,
     *  0.70077246146712664722D-19,-.10909616144053653581D-19/
c 
c      Coefficients for aa=1/15, 1/5.3d0
c 
        data coefs2/
     1  0.89805987077586343707D+00,-.41032532576611054362D-01,
     2  0.16209747524170531119D-02,-.85036347080121688414D-04,
     3  0.53546335050405476057D-05,-.38419939114265847633D-06,
     4  0.30449601248472359744D-07,-.26118923597442445666D-08,
     5  0.23905606264920304191D-09,-.23104541665515644027D-10,
     6  0.23395093084178672016D-11,-.24666860289090228845D-12,
     7  0.26948639829947749980D-13,-.30385350820983721391D-14,
     8  0.35242579224695154966D-15,-.41932532005257319720D-16,
     9  0.51062370910638810255D-17,-.63510785818162754924D-18,
     *  0.80522071301475187570D-19,-.10218151895178754782D-19/
c 
c      Coefficients for aa=1/5.3d0, 1/2.5d0
c 
        data coefs3/
     1  0.80624135836386946839D+00,-.49725826975926992323D-01,
     2  0.24393578230521649195D-02,-.15017585534319530860D-03,
     3  0.10646254877701552133D-04,-.83279579878924993739D-06,
     4  0.70122405079781471004D-07,-.62552614242547852133D-08,
     5  0.58474339053445056855D-09,-.56832816430267113131D-10,
     6  0.57093663632819514041D-11,-.59014224936643652001D-12,
     7  0.62538024065268930625D-13,-.67747085492694378598D-14,
     8  0.74844957538995627386D-15,-.84158806229556451611D-16,
     9  0.96156010914349796671D-17,-.11147407649020231352D-17,
     *  0.13093973951286023137D-18,-.15349144923597118126D-19/
c 
        n=19
c 
c       evaluate the expansion in the case when 15 < x < \infty
c 
        if(x .lt. 15) goto 2000
c 
        done=1
        aa=0
        bb=done/15
c 
        alpha=2/(bb-aa)
        beta=1-alpha*bb
        t=alpha/x+beta
c 
        call CHFUNDER(t,e1,der,coefs1,N)
c 
        e1=e1*exp(-x)/x
c 
        return
 2000 continue
c 
c       evaluate the expansion in the case when 5.3 < x < 15
c 
        if( (x .gt. 15) .or. (x. lt. 5.3d0) ) goto 3000
c 
        done=1
        aa=done/15
        bb=done/5.3d0
c 
        alpha=2/(bb-aa)
        beta=1-alpha*bb
c 
        t=alpha/x+beta
c 
        call CHFUNDER(t,e1,der,coefs2,N)
c 
        e1=e1*exp(-x)/x
c 
        return
c 
 3000 continue
c 
c       evaluate the expansion in the case when 2.5 < x < 5.3
c 
        if( (x .gt. 5.3d0) .or. (x. lt. 2.5d0) ) goto 4000
c 
        done=1
        aa=done/5.3d0
        bb=done/2.5d0
c 
        alpha=2/(bb-aa)
        beta=1-alpha*bb
c 
        t=alpha/x+beta
c 
        call CHFUNDER(t,e1,der,coefs3,N)
c 
        e1=e1*exp(-x)/x
c 
        return
c 
 4000 continue
c 
 4200 continue
c 
c       x is less than 2.5. Use the taylor series
c 
        n2=26
        n2=29
        call e1fun_tayl(x,n2,e1)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine e1fun_tayl(x,n,e1)
        implicit real *8 (a-h,o-z)
        save
        dimension coefs(100)
c 
c 
        data ifcalled/0/
c 
        data gamma/0.577215664901532860606512090082D+00/
c 
c       if this is the first call tothis subroutine,
c       construct the table of factorials
c 
        if(ifcalled .eq. 1) goto 1400
c 
        coefs(1)=1
        do 1200 i=2,45
        coefs(i)=coefs(i-1)*i
c 
        coefs(i-1)=1/(coefs(i-1)*(i-1))
 1200 continue
c 
        ifcalled=1
c 
 1400 continue
c 
c        evaluate the exponential integral function E_1
c        using the Taylor series
c 
        xk=-x
        sum=0
        do 1600 k=1,n
c 
        sum=sum+xk*coefs(k)
c 
        xk=-xk*x
 1600 continue
c 
        e1=-(sum+log(x)+gamma)
c 
        return
        end
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
