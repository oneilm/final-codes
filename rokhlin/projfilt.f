        implicit real *8 (a-h,o-z)
        real *8 cm(100 000),ts(100 000),finrea(100 000),
     1          finima(100 000),foutrea(100 000),
     2      foutima(100 000),rea(2)
  
        complex *16 ff(100 000),ff2(100 000),fout2(100 000)
c 
        dimension w(100 000)
c 
        complex *16 ima,fin(100 000),fout(100 000),com
c 
        equivalence (rea,com)
c 
        data ima/(0.0d0,1.0d0)/
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER n'
        READ *,n
        CALL PRINf('n=*',n,1 )
  
        eps2=1.0d-12
        sgain=1000000
        sgain=1.1
        sgain=100
        aa=30
cccc        aa=4
        m1=0
        m2=2
c 
c       construct the function to be filtered
c 
        done=1
        pi=atan(done)*4
        h=2*pi/n
        do 1200 i=1,n
c 
        ts(i)=(i-1)*h
cccc        fin(i)=exp(ima*ts(i)*m1)+exp(ima*ts(i)*m2)
        fin(i)=exp(ima*aa*cos(ts(i))) +10
 1200 continue
c 
c       construct the coefficients of the "filter"
c 
        eps=1.0d-20
        do 1400 i=1,n
c 
        cm(i)=exp(-(i-1)**2* 0.1)
c 
        if(cm(i) .gt. eps) cm(i)=1/cm(i)
        if(cm(i) .le. eps) cm(i)=1/eps
 1400 continue
  
        call prin2('cm as created*',cm,n)
c 
c        filter the test function
c 
        ifinit=1
        ifnewton=1
  
        call projfilt(ier,fin,n,cm,sgain,ifnewton,eps2,fout,
     1      ff,ff2,w,ifinit)
c 
c       plot both input and output signals
c 
        do 2200 i=1,n
c 
        com=fin(i)
c 
        finrea(i)=rea(1)
        finima(i)=rea(2)
c 
        com=fout(i)
c 
        foutrea(i)=rea(1)
        foutima(i)=rea(2)
 2200 continue
c 
        iw=21
        call lotagraph2(iw,ts,finrea,n,ts,foutrea,n,
     1       'real parts*')
c 
c 
        iw=22
        call lotagraph2(iw,ts,finima,n,ts,foutima,n,
     1       'imaginary parts*')
c 
  
        sgain2=sgain-0.1
        sgain2=sgain
  
        sgain2=sgain-1.0d
  
  
        ifinit=0
  
        call projfilt(ier,fout,n,cm,sgain2,ifnewton,eps2,fout2,
     1      ff,ff2,w,ifinit)
  
  
        diffmax=0
        do 2400 i=1,n
c 
        ddd=abs(fout2(i)-fout(i))
        if(ddd .gt. diffmax) diffmax=ddd
 2400 continue
c 
        call prin2('and diffmax=*',diffmax,1)
  
  
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine projfilt(ier,fin,n,cm,sgain,ifnewton,eps,fout,
     1      ff,ff2,w,ifinit)
        implicit real *8 (a-h,o-z)
        save
        complex *16 fin(1),ff(1),ff2(1),fout(1)
        dimension cm(1),w(1)
c 
c       This subroutine applies a projecting filter to a user-supplied
c       (complex) function fin. Fin is assumed to be periodic in
c       the standard sense, i.e. fin is tabulated at n nodes, and
c       it is assumed that fin(n+1)=fin(1). The user supplies the
c       filter to the subroutine via the real array cm.
c       Depending on  the value of the input parameter
c       ifnewton (see below) the actual supergain of the constructed
c       function will either be (more or less) exactly equal to the
c       user-specified value, or be only a reasonable approximation to
c       the requested supergain (this way, the scheme is noticeably
c       simpler, and who needs the supergain exactly, anyway!).
c 
c                    Input parameters:
c 
c  fin - the function to be filtered
c  n - the number of elements in array fin. Since fin is going to be
c       Fourier transformed, n should be a product of small primes
c       (unless the user likes order n**2 algorithms!)
c  cm - the array of penalty coefficients corresponding to the first n/2+1
c       modes.
c 
c   Explanation: The filter the subroutine applies is a "zero-phase" one,
c       in the sense that all frequencies are scaled by real coefficients,
c       and the function exp(i \cdot m \cdot t) is scaled by the same
c       coefficient aas the function exp(-i \cdot m \cdot t). Thus, only
c       the first n/2+1 elements of array cm are used
c  sgain - the amount of supergain the user would like the appreoximation
c       to possess
c  ifnewton - tells the subroutine whether to enforce the supergain
c       condition exactly, or to settle for an approximate result;
c     ifnewton=0 will cause the subroutine to settle for the approximate
c       result
c     ifnewton=1 will cause the subroutine to enforce the supergain
c       condition exactly
c  eps - the accuracy to which the subroutine will conduct the Newton
c       process; this parameter is only used if the parameter ifnewton
c     (see above) had been set to 1.
c  ifinit - the parameter telling the soubroutine whether it should
c       initialize the FFT routine it uses.
c     ifinit=1 will cause the subroutine to initialize the FFT routine
c     ifinit=0 will cause the subroutine to skip the initialization;
c       in this case, the first 4*n+30 elements of the array w should
c       not have been changed since the last call to this subroutine
c       when the parameter ifinit had been set to 1. Also, please note
c       that ifinit should be set to 1 each time this subroutine is
c       used with a new n
c  w - this is only an input parameter if ifinit has been set to 0;
c       otherwise, it is an output parameter (the first n*4+30 elements
c       of it; the remainder is treated as work space)
c 
c                     Output arrays:
c 
c  fout - the filtered function
c  ff - Fourier transform of fin
c  ff2 - Fourier transform of out
c 
c                      Work arrays:
c 
c  w - must be at least 9*n+100 real *8 elements long. Please note that
c       if ifinit (see above) has been set to 0, the first 4*n+30
c       elements of w are treated as an input array, and will be used
c       as such by the FFT routine called by this subroutine. In this
c       case, these elements should not have been changed since the
c       preceding call by this subroutine with the same n, and ifinit
c       set to 1.
c 
c 
c       . . . allocate memory for the filtering procedure
c 
        iwsave=1
        lwsave=4*n+30
c 
        icm2=iwsave+lwsave
        lcm2=n+2
c 
        iprods=icm2+lcm2
        lprods=n+2
c 
        iw=iprods+lprods
        lw=3*n+10
c 
c       filter the user-specified signal
c 
        call prin2('in prolfilt before prolflt, sgain=*',sgain,1)
  
        call projflt(ier,fin,n,cm,sgain,ifnewton,eps,fout,
     1        ff,ff2,w(iwsave),w(icm2),w(iprods),ifinit,w(iw))
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine projflt(ier,fin,n,cm,sgain,ifnewton,eps,fout,
     1        ff,ff2,wsave,cm2,prods,ifinit,w)
        implicit real *8 (a-h,o-z)
        save
        complex *16 fin(1),wsave(1),ff(1),
     1      ff2(1),fout(1),w(1)
c 
        dimension cm(1),cm2(1),prods(1)
c 
        ier=0
  
        call prin2('in projflt, eps=*',eps,1)
c 
c       if the user so requested - initialize the FFT code
c 
        if(ifinit .eq. 1) call DCFFTI (N,WSAVE)
c 
c        Fourier transform the user-provided signal
c 
        call projfarm(fin,ff,n*2)
c 
        call DCFFTF(N,ff,WSAVE)
c 
        done=1
        d=done/n
        do 1200 i=1,n
        ff(i)=ff(i)*d
 1200 continue
c 
c       construct the array of absolute values of Fourier coefficients
c 
        cm2(1)=cm(1)
        do 1300 i=1,n/2
        cm2(i+1)=cm(i+1)
        cm2(n-i+1)=cm(i+1)
 1300 continue
c 
        totpowr=0
        do 1400 i=1,n
c 
        prods(i)=abs(ff(i))
        totpowr=totpowr+prods(i)**2
 1400 continue
c 
c       determine the Lagrange multiplier
c 
        totpowr=totpowr*sgain
c 
        call projfal(jer,cm2,prods,n,totpowr,
     1      rnum,denom,rr,drrds,rlout,w)
c 
        call prinf('after first projfal, jer=*',jer,1)
  
        if(jer .eq. 4) ier=4
        if(jer .eq. 16) ier=512
        if(jer .eq. 16) return
c 
        sgainout=rr
c 
c       if the user has requested that the supergain condition
c       be enforced approximately - bypass the Newton process
c 
        if(ifnewton .eq. 0) goto 3200
        if(ier .ne. 0) goto 3200
c 
c       if the user so requested, use Newton iterations
c       to obtain exact supergain condition
c 
        rloldold=rlout
c 
        totpowr3=totpowr
        sgainout1=sgainout
c 
        ddold=1.0d20
        dd=1.0d20
  
        do 2500 ijk=1,10
c 
        niter=ijk
  
        call prinf('ijk=*',ijk,1)
  
        ddold=abs(rr-sgain)/sgain
c 
        totpowr3=totpowr3-(rr-sgain)/drrds
        rlold=rlout
c 
        call prin2('rr-sgain=*',rr-sgain,1)
        call prin2('ddold=*',ddold,1)
  
        call projfal(ier,cm2,prods,n,totpowr3,
     1      rnum,denom,rr,drrds,rlout,w)
c 
        dd=abs(rr-sgain)/sgain
c 
c       if the user-prescribed accuracy has been achieved -
c       terminate the Newton process
c 
        if(ddold .lt. eps) goto 3200
c 
c       if the error has failed to decrease on this step
c       bomb, unless the desired accuracy had been achieved
c       on the preceding step
c 
        if(dd .le. ddold) goto 2500
c 
        if(dd .lt. eps) rlout=rlold
        if(dd .lt. eps) goto 3200
c 
        ier=64
        rlout=rloldold
        goto 3200
c 
 2500 continue
c 
        ier=32
        rlout=rloldold
c 
 3200 continue
c 
        do 3400 i=1,n
c 
        ff2(i)=ff(i)/(1+rlout*cm2(i))
 3400 continue
c 
c       Fourier transform the thing back
c 
        call projfarm(ff2,fout,n*2)
c 
        call DCFFTb(N,fout,WSAVE)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine projfal(ier,cm,prods,n,s,
     1      rnum,denom,rr,drrds,rl,w)
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
     1      rr,drrds,rl,w(irddral),rnum,denom)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine projfal0(ier,cm,prods,n,s,alph2,alph2der,
     1      rr,drrds,rl,drrdalph,rnum,den)
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
c       Then, it evaloates the array of alph2der of derivatives of the
c       coefficients alp2 with respect to s.
c 
c       Finally, it uses the obtained data to evaluate the ratio
c 
c       rr= { \sum_{m=1}^n  cm(m) * prods(m)^2 \over
c             \sum_{m=1}^n  prods(m)^2 },                               (3)
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
        call projfbis(ier,cm,prods,n,s,rl)
c 
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
        subroutine projfarm(x,y,n)
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
        subroutine projfbis(ier,cm,prods,n,totpowr,rlout)
        implicit real *8 (a-h,o-z)
        save
        dimension cm(1),prods(1)
c 
c       initialize the bisection process
c 
        ier=0
        rl1=0
        rl2=1.0d10
        rlout=0
c 
        call projfsum(cm,prods,n,rl1,f1)
        call projfsum(cm,prods,n,rl2,f2)
c 
        if( (f1 .gt. totpowr) .and. (f2 .lt. totpowr) ) goto 1200
c 
        if(f1 .lt. totpowr) ier=4
        if(f2 .gt. totpowr) ier=16
        return
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
        if(f22 .ge. totpowr) goto 1350
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
        if(f3 .lt. totpowr) goto 1400
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
        ff=sum-totpowr
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
