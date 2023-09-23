        implicit real *8 (a-h,o-z)
        complex *16 ima,zs(100 000)
c 
        data ima/(0.0d0,1.0d0)/
  
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
         PRINT *, 'ENTER m'
         READ *,m
         CALL PRINf('m=*',m,1)
  
         a=10
         eps=1.0d-12
        call testit(a,m,zs,eps)
  
  
        stop
        end
  
  
  
  
  
  
        subroutine testit(a,m,zs,eps)
        implicit real *8 (a-h,o-z)
        save
        complex *16 ima,z,fjs(10 000),rjs(10 000),
     1              rats(10 000),diffs(10 000),zs(m,m)
c 
        data ima/(0.0d0,1.0d0)/
c 
c        construct the grid on which the subroutine is to be
c        tested
c 
         h=2*a/(m-1)
         do 1400 i=1,m
         do 1200 j=1,m
c 
         zs(j,i)=-a+(j-1)*h +ima*(-a+(i-1)*h)
c 
  
 1200 continue
 1400 continue
c 
cccc         call prin2('zs as created*',zs,m*m*2)
c 
c        evaluate the Bessel functions at all points on the grid
c        using both subroutines, and calculate the maximum error
c 
  
         eps2=eps**2
c 
         errsmax=0
  
         icount=0
         do 2400 i=1,m
         do 2200 j=1,m
c 
cccc        call prinf('i=*',i,1)
cccc        call prinf('j=*',j,1)
  
  
  
        call cjseval(zs(j,i),eps,n1,rjs,nmax)
c 
  
cccc        call prin2('after cjseval, rjs=*',rjs,n1*2)
  
  
        call JBFUN(Zs(j,i),fJs(2),n2,EPS2)
c 
  
        errmax=0
c 
        do 1600 k=1,n1
        d=(abs(fjs(k)-rjs(k))/rjs(k))
c 
        if(d .gt. errmax) errmax=d
 1600 continue
c 
cccc        call prin2('errmax=*',errmax,1)
  
        if(errsmax .lt. errmax) errsmax=errmax
  
  
        goto 2200
        icount=icount+1
        if(icount .eq. 20) icount=1
  
        if(icount .ne. 1) goto 2200
  
c 
 2200 continue
  
        call prinf('i=*',i,1)
        call prin2('errsmax=*',errsmax,1)
 2400 continue
  
        call prin2('errsmax=*',errsmax,1)
        call prinf('nmax=*',nmax,1)
  
  
        return
        end
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code and the beginning of the
c        Bessel function code proper.
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine cjseval(z,eps,n,rjs,nmax)
        implicit real *8 (a-h,o-z)
c 
        save
        complex *16 z,rjs(1)
c 
c       This subroutine evaluates all J-Bessel functions of the
c       user-supplied complex argument z that are greater than eps.
c 
c                   Input parameters:
c 
c  z - the argument of the Bessel functions to be evaluated
c  eps - the accuracy
c 
c                   Output parameters:
c 
c  n - the number of Bessel functions that are greater than eps
c  rjs - the array of Bessel functions
c  rjs - the total number of elements of array rjs actually used
c       by this subroutine
c 
        call cjseva0(z,eps,n,rjs,rjs,nmax)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine cjseva0(z,eps,n,rjs,rrjs,nmax)
        implicit real *8 (a-h,o-z)
c 
        save
        complex *16 rjs(1),z,cd,cd2,scale,ima,
     1      sum,cdd,com,coeima,sum1,sum2,sum3,sum4
        dimension rea(2),rrjs(1)
        equivalence (rea(1),com)
c 
        data ima/(0.0d0,1.0d0)/
c 
c       find n corresponding to these rl and x
c 
        rl=1/eps
c 
        rjs(1)=0
        rjs(2)=1
c 
        com=z
c 
        cd=2/z
        cd2=cd
c 
cccc        rrl=rl**2*exp(abs(rea(2)))
        rrl=rl*exp(abs(rea(2))) * (abs(rea(1))+abs(rea(2)))
cccc        rrl2=rrl**2*sqrt(sqrt(sqrt(rl) ))
cccc        rrl2=rrl**2* ( sqrt(rl) )
        rrl2=rrl**2 *rl
c 
        do 1200 i1=2,10 000 000
c 
        i=i1-1
        rjs(i1+1)=cd*rjs(i1)-rjs(i)
c 
        cd=cd+cd2
        n=i
c 
        i3=i1+i1+1
        dddd=rrjs(i3)**2+rrjs(i3+1)**2
        if(dddd .gt. rrl2) goto 1400
 1200 continue
 1400 continue
c 
        n=n+2
        nmax=n
        cd=cd+2*cd2
  
cccc        call prin2('rjs as constructed*',rjs,n*2)
  
c 
c       starting from n-th position, conduct recursion down
c 
        rjs(n)=0
        rjs(n-1)=1
c 
        cd=cd-3*cd2
c 
        i=n-1
        sum1=0
        sum2=0
        sum3=0
        sum4=0
        do 1700 ijk=1,n
c 
        i1=i
        i=i1-1
        rjs(i)=cd*rjs(i1)-rjs(i1+1)
c 
        cd=cd-cd2
        sum1=sum1+rjs(i1)
c 
        if(i .eq. 1) goto 1800
c 
        i1=i
        i=i1-1
        rjs(i)=cd*rjs(i1)-rjs(i1+1)
c 
        cd=cd-cd2
        sum2=sum2+rjs(i1)
c 
        if(i .eq. 1) goto 1800
c 
        i1=i
        i=i1-1
        rjs(i)=cd*rjs(i1)-rjs(i1+1)
c 
        cd=cd-cd2
        sum3=sum3+rjs(i1)
        if(i .eq. 1) goto 1800
c 
        i1=i
        i=i1-1
        rjs(i)=cd*rjs(i1)-rjs(i1+1)
c 
        cd=cd-cd2
        sum4=sum4+rjs(i1)
c 
        if(i .eq. 1) goto 1800
 1700 continue
 1800 continue
c 
        n4=(n-2)/4
        i=n-2-n4*4
        if(i .eq. 0) cdd=2
        if(i .eq. 1) cdd=2*ima
        if(i .eq. 2) cdd=-2
        if(i .eq. 3) cdd=-2*ima
c 
        coeima=-ima
c 
        if(rea(2) .le. 0) goto 1850
c 
        coeima=-coeima
        if(i .eq. 1) cdd=-2*ima
        if(i .eq. 3) cdd=2*ima
 1850 continue
c 
        sum=( ( ( (sum4*coeima)+sum3)*coeima+sum2)
     1      *coeima+sum1)*cdd+rjs(1)
c 
        if(rea(2) .le. 0)
     1      scale=exp(z*ima)/sum
c 
        if(rea(2) .gt. 0)
     1      scale=exp(-z*ima)/sum
c 
c       . . . scale
c 
        nnn=0
        ddd=scale*conjg(scale)*
     1      (rjs(1)*conjg(rjs(1))+rjs(2)*conjg(rjs(2)) )
c 
        eps22=ddd*eps**2
c 
        do 2000 i=1,n
c 
        rjs(i)=rjs(i)*scale
c 
        i3=i+i
        d=rrjs(i3-1)**2+rrjs(i3)**2
        if(d .gt. eps22) nnn=i
 2000 continue
c 
        n=nnn
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine rjseval(x7,eps,rjs,nstop,nmax)
        implicit real *8 (a-h,o-z)
        save
        dimension coefs1(3,50),rjs(1),coefs2(3,50)
c 
c       This subroutine evaluates Bessel (J) functions for the
c       user-specified real argument x. This subroutine has been
c       tested for |x| \in [1.0d-20, 100 000]
c 
c                       Input parameters:
c 
c  x7 - the point where the Bessel functions are to be evaluated
c  eps - the size of the highest order Bessel function to be evaluated
c       accurately
c 
c                       Output parameters:
c 
c  rjs - the Bessel functions
c  nstop - the number of the element in the array rjseval after which all
c       elements are less than eps; please note that it tends to be on
c       the conservative side, in the sense that several elements BEFORE
c       the element number nstop might also be smaller than eps
c  nmax - the total number of elements in array rjs that have been used
c       by this subroutine
c 
        data coefs2/
     1     0.485161590E+01,-.1533E+02,0.664E+02,
     2     0.557479200E+01,-.1416E+02,0.619E+02,
     3     0.625643784E+01,-.1349E+02,0.603E+02,
     4     0.689232299E+01,-.1256E+02,0.562E+02,
     5     0.750445331E+01,-.1199E+02,0.533E+02,
     6     0.809769445E+01,-.1183E+02,0.545E+02,
     7     0.866076542E+01,-.1124E+02,0.506E+02,
     8     0.920747303E+01,-.1080E+02,0.484E+02,
     9     0.973283423E+01,-.1029E+02,0.463E+02,
     *     0.102674828E+02,-.1072E+02,0.520E+02,
     1     0.107638292E+02,-.1005E+02,0.461E+02,
     2     0.112552191E+02,-.9771E+01,0.452E+02,
     3     0.117400439E+02,-.9684E+01,0.462E+02,
     4     0.122131046E+02,-.9548E+01,0.462E+02,
     5     0.126749006E+02,-.9358E+01,0.470E+02,
     6     0.131292750E+02,-.9155E+01,0.453E+02,
     7     0.135752300E+02,-.8913E+01,0.442E+02,
     8     0.140178440E+02,-.8854E+01,0.439E+02,
     9     0.144478335E+02,-.8637E+01,0.437E+02,
     *     0.148793032E+02,-.8786E+01,0.461E+02,
     1     0.152948194E+02,-.8414E+01,0.439E+02,
     2     0.157161866E+02,-.8594E+01,0.451E+02,
     3     0.161193485E+02,-.8204E+01,0.434E+02,
     4     0.165234217E+02,-.8085E+01,0.422E+02,
     5     0.169256656E+02,-.8157E+01,0.445E+02,
     6     0.173176690E+02,-.7952E+01,0.432E+02,
     7     0.177059893E+02,-.7824E+01,0.429E+02,
     8     0.180918885E+02,-.7694E+01,0.424E+02,
     9     0.184765030E+02,-.7763E+01,0.439E+02,
     *     0.188575190E+02,-.7824E+01,0.446E+02,
     1     0.192287687E+02,-.7652E+01,0.450E+02,
     2     0.195996729E+02,-.7663E+01,0.462E+02,
     3     0.199721860E+02,-.7815E+01,0.473E+02,
     4     0.203136715E+02,-.6779E+01,0.395E+02,
     5     0.206806958E+02,-.6981E+01,0.417E+02,
     6     0.210426011E+02,-.7173E+01,0.448E+02,
     7     0.213951929E+02,-.7009E+01,0.434E+02,
     8     0.217508944E+02,-.7219E+01,0.470E+02,
     9     0.221018888E+02,-.7262E+01,0.474E+02,
     *     0.224470460E+02,-.7166E+01,0.470E+02,
     1     0.227865308E+02,-.7011E+01,0.470E+02,
     2     0.231219476E+02,-.6682E+01,0.444E+02,
     3     0.234611544E+02,-.6669E+01,0.449E+02,
     4     0.237988923E+02,-.6761E+01,0.468E+02,
     5     0.241322868E+02,-.6723E+01,0.465E+02,
     6     0.244623627E+02,-.6619E+01,0.461E+02,
     7     0.247883271E+02,-.6500E+01,0.465E+02,
     8     0.251209727E+02,-.6755E+01,0.498E+02,
     9     0.254401624E+02,-.6480E+01,0.476E+02,
     *     0.257623322E+02,-.6418E+01,0.470E+02/
c 
        data coefs1/
     1     0.427413535E+01,-.3349E+01,0.273E+01,
     2     0.502544757E+01,-.2925E+01,0.240E+01,
     3     0.572223994E+01,-.2609E+01,0.233E+01,
     4     0.638983543E+01,-.2453E+01,0.236E+01,
     5     0.702755808E+01,-.2364E+01,0.260E+01,
     6     0.763437008E+01,-.2279E+01,0.285E+01,
     7     0.819679197E+01,-.2002E+01,0.258E+01,
     8     0.875971135E+01,-.1932E+01,0.281E+01,
     9     0.931099024E+01,-.1937E+01,0.317E+01,
     *     0.984156676E+01,-.1910E+01,0.353E+01,
     1     0.103166856E+02,-.1473E+01,0.300E+01,
     2     0.108289737E+02,-.1495E+01,0.332E+01,
     3     0.113081940E+02,-.1331E+01,0.331E+01,
     4     0.117971436E+02,-.1372E+01,0.374E+01,
     5     0.122488845E+02,-.1140E+01,0.359E+01,
     6     0.127193238E+02,-.1189E+01,0.402E+01,
     7     0.131742042E+02,-.1158E+01,0.429E+01,
     8     0.135966029E+02,-.8787E+00,0.398E+01,
     9     0.140322900E+02,-.8197E+00,0.427E+01,
     *     0.144431279E+02,-.5903E+00,0.422E+01,
     1     0.148814159E+02,-.6986E+00,0.487E+01,
     2     0.152936276E+02,-.6033E+00,0.505E+01,
     3     0.156980072E+02,-.4458E+00,0.493E+01,
     4     0.160946391E+02,-.3004E+00,0.498E+01,
     5     0.164917605E+02,-.2255E+00,0.525E+01,
     6     0.168752327E+02,-.5469E-01,0.528E+01,
     7     0.172828680E+02,-.1421E+00,0.582E+01,
     8     0.176252242E+02,0.3692E+00,0.500E+01,
     9     0.180235004E+02,0.2490E+00,0.568E+01,
     *     0.184123783E+02,0.2194E+00,0.604E+01,
     1     0.187785237E+02,0.3685E+00,0.608E+01,
     2     0.191348832E+02,0.5594E+00,0.605E+01,
     3     0.194883340E+02,0.7464E+00,0.602E+01,
     4     0.198528220E+02,0.7489E+00,0.657E+01,
     5     0.201975720E+02,0.9743E+00,0.650E+01,
     6     0.205420819E+02,0.1148E+01,0.650E+01,
     7     0.209215340E+02,0.9597E+00,0.730E+01,
     8     0.212504013E+02,0.1271E+01,0.686E+01,
     9     0.215955692E+02,0.1332E+01,0.720E+01,
     *     0.219204337E+02,0.1635E+01,0.681E+01,
     1     0.222742886E+02,0.1541E+01,0.747E+01,
     2     0.225804435E+02,0.1935E+01,0.707E+01,
     3     0.229174710E+02,0.1968E+01,0.745E+01,
     4     0.232625937E+02,0.1930E+01,0.792E+01,
     5     0.235946826E+02,0.2022E+01,0.801E+01,
     6     0.239184511E+02,0.2188E+01,0.784E+01,
     7     0.242375327E+02,0.2305E+01,0.810E+01,
     8     0.245690213E+02,0.2303E+01,0.847E+01,
     9     0.248411680E+02,0.2837E+01,0.792E+01,
     *     0.251856702E+02,0.2681E+01,0.866E+01/
c 
        x=abs(x7)
        if(x .lt.1) x=1
c 
        dd=log10(eps)
        jj=1.01-dd
        j=jj+15
c 
c       find n corresponding to this x
c 
        d=x**0.3333333333d0
c 
        if(x .gt. 500) goto 1100
c 
        nmax=x+coefs1(1,j)*d+coefs1(2,j)+coefs1(3,j)/d
        nstop=x+coefs1(1,jj)*d+coefs1(2,jj)+coefs1(3,jj)/d
c 
        goto 1190
 1100 continue
c 
        nmax=x+coefs2(1,j)*d+coefs2(2,j)+coefs2(3,j)/d
        nstop=x+coefs2(1,jj)*d+coefs2(2,jj)+coefs2(3,jj)/d
c 
 1190 continue
c 
        if(abs(x7) .gt. 1.0d-3) goto 1193
        nmax=12
        nstop=10
 1193   continue
c 
        if(abs(x7) .gt. 1.0d-5) goto 1195
        nmax=6
        nstop=5
 1195 continue
c 
c       conduct recursion down
c 
        x=abs(x7)
        rjs(nmax+2)=0
        rjs(nmax+1)=1
        d=0
c 
        rr=2*nmax/x
        rrr=2/x
c 
        do 1200 i=nmax,1,-1
c 
        rjs(i)=rr*rjs(i+1)-rjs(i+2)
        rr=rr-rrr
c 
        d=d+rjs(i)**2
 1200 continue
c 
        d=2*d-rjs(1)**2
        dd=1/sqrt(d)
c 
c       normalize
c 
        sss=1
        if(x7 .lt. 0) sss=-1
c 
        do 1400 i=1,nmax
c 
        rjs(i)=rjs(i)*dd
c 
        dd=dd*sss
 1400 continue
c 
        return
        end
  
