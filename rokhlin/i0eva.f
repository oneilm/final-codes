        implicit real *8 (a-h,o-z)
        real *16 i0,dt
        real *8 errs(500)
c 
        CALL PRINI(6,13)
C 
c       . . . construct the testing parameters
c 
c 
         PRINT *, 'ENTER t'
         READ *,t
         CALL PRIN2('t=*',t,1)
  
        do 2000 ijk=1,200
        t=ijk
  
c 
c        calculate i0(t) via trapezoidal quadrature
c 
        dt=t
        call i0comp01(ier,qabs(dt),i0)
        call prinq('i0=*',i0,1)
        call prinq('i0 e^(-t) =*',i0*dexp(-t),1)
c 
c       . . . and via expansion
c 
        call i0eva(t,fi0)
        call prin2('and fi0 via expansion*',fi0,1)
        call prinq('and fi0-i0*',fi0-i0,1)
        call prinq('and (fi0-i0)/i0=*',(fi0-i0)/i0,1)
         errs(ijk)=(fi0-i0)/i0
 2000 continue
         call prin2('and errs=*',errs,200)
        stop
        end
c 
c 
c 
c 
c 
        subroutine i0comp01(ier,x,i0)
        implicit real *16 (a-h,o-z)
        save
        real *16  i0
c 
c       this subroutine evaluates the function
c       i0(x) for values of x between 0 and 250.
c 
c                    input parameters:
c 
c  x - the values of the argument at which the functioon i0
c      is to be evaluated
c 
c                    output parameters:
c 
c  ier - error return code; ier=0 means successful conclusion.
c                           ier=8 means that x > 250
c                           ier=16 means that x < 0
c  i0 - the value of the bessel function (real *16).
c 
c       . . . determine the number of nodes in the trapezoidal
c             approximation
c 
        ier=0
        if(x .lt. 0) ier=16
        if(x .gt. 250) ier=8
        if(ier .ne. 0) return
c 
        n=32
        if(x .gt. 3) n=64
        if(x .gt. 20) n=128
        if(x .gt. 100) n=256
c 
c        use trapezoidal approximation to compute i0(x)
c 
        done=1
        pi=qatan(done)*4
c 
        h=2*pi/n
        dd=0
        hcos=qcos(h)
        hsin=qsin(h)
cccc         hcos=qsqrt(done-hsin**2)
        costhet=1
        sinthet=0
c 
        do 1200 i=1,n
        theta=(i-1)*h
cccc        d=qcos(theta) * x
        d=costhet*x
        ddd=qexp(d)
cccc        dd=dd+qexp(d)+qexp(-d)
        dd=dd+ddd+done/ddd
c 
        costhet2=costhet*hcos-sinthet*hsin
        sinthet2=costhet*hsin+sinthet*hcos
        costhet=costhet2
        sinthet=sinthet2
c 
 1200 continue
        i0=dd*h/4/pi
        return
        end
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c         this is the end of the debugging code and the beginning
c         of the actual routine for the evaluation of the function
c         i0 - modified bessel function
c 
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine i0eva(x7,fi0)
        implicit real *8 (a-h,o-z)
        save
        dimension p0to3(21)
        data p0to3/
     1  0.3674336090541583E+00, -.1483942216332326E+00,
     2  0.7538109249292410E-01, -.3429202889217209E-01,
     3  0.1335283653302592E-01, -.4461174022644684E-02,
     4  0.1294920713528336E-02, -.3309968933824789E-03,
     5  0.7542634832396139E-04, -.1548671523094681E-04,
     6  0.2891239479464111E-05, -.4946276244165349E-06,
     7  0.7806369735032947E-07, -.1143166990263007E-07,
     8  0.1561212097141843E-08, -.1997725692654135E-09,
     9  0.2403196202942698E-10, -.2710668355445937E-11,
     *  0.2916600802142503E-12, -.3325583782158792E-13,
     *  0.3199318882156914E-14/
c 
        dimension p3to6(19)
        data p3to6/
     1  0.1941982776283822D+00, -.2323945533089018D-01,
     2  0.4244030631168826D-02, -.8759373293717004D-03,
     3  0.1909924203989986D-03, -.4222271728163107D-04,
     4  0.9160432785148500D-05, -.1904704016240210D-05,
     5  0.3739731123576933D-06, -.6879564187569590D-07,
     6  0.1182061524819559D-07, -.1896547862750232D-08,
     7  0.2845293860689127D-09, -.4001288511461210D-10,
     8  0.5284459250803188D-11, -.6526404831544231D-12,
     9  0.7680294564558936D-13, -.9627486914965938D-14,
     *  0.1006473132584426D-14/
c 
        dimension p6to9(13)
        data p6to9/
     1  0.3989284896164599D+00, 0.5092106731858055D-01,
     2  -.6127221009097870D-02, 0.6131729304562125D+00,
     3  -.4628528030678180D+01, -.1205084525934594D+02,
     4  0.7340780560072762D+03, -.8729685960965613D+04,
     5  0.5831335775767627D+05, -.2430536532475642D+06,
     6  0.6287065304211588D+06, -.9292779115166327D+06,
     7  0.6032651250101362D+06/
c 
        dimension p9to12(12)
        data p9to12/
     1  0.3989407763246097D+00, 0.5000863927874731D-01,
     2  0.2226512800214346D-01, 0.1652919141248679D+00,
     3  -.1957907568960066D+01, 0.1894748754972728D+02,
     4  -.1103884918016903D+03, 0.3682201782413656D+03,
     5  -.5202693680397354D+03, 0.7939108402881759D+01,
     6  0.6287065304211588D+06, -.9292779115166327D+06/
c 
        dimension p12to24(12)
        data p12to24/
     1  0.3989422838292772D+00, 0.4986710003332280D-01,
     2  0.2811267232744196D-01, 0.2585760370970408D-01,
     3  0.1658198437867097D+00, -.2952322171370709D+01,
     4  0.5466563227185643D+02, -.6921662882875580D+03,
     5  0.6144958292895270D+04, -.3605643906301946D+05,
     6  0.1258751734607877D+06, -.1960251421440677D+06/
c 
        dimension p24to48(8)
        data p24to48/
     1  0.3989422803956930D+00, 0.4986778658145504D-01,
     2  0.2805045279642165D-01, 0.2923081851609006D-01,
     3  0.4428963589957261D-01, 0.1017678215019001D+00,
     4  0.6423782560416375D-01, 0.1927258316093262D+01/
c 
        dimension p48to96(8)
        data p48to96/
     1  0.3989422804014209D+00, 0.4986778505649491D-01,
     2  0.2805062761844382D-01, 0.2921959904479739D-01,
     3  0.4472650398122760D-01, 0.9140365337986182D-01,
     4  0.2035846167608939D+00, 0.1104329030954804D+01/
        dimension p96to192(6)
        data p96to192/
     1  0.3989422804013822D+00, 0.4986778509112728D-01,
     2  0.2805061543528473D-01, 0.2922179806303123D-01,
     3  0.4451046595352375D-01, 0.1022708740517523D+00/
c 
c       this subroutine evaluates the modified bessel
c       function I_0 of a real argument x7
c 
c               input parametedrs:
c 
c  x7 - the argumant of the bessel function to be evaluated
c 
c               output parameters:
c 
c  fi0 - the modified bessel function i0 of the argument x7
c 
c           . . . evaluate the function i0
c 
        x=dabs(x7)
        done=1
c 
c        . . . for  0< x <3
c 
        if(x .gt. 3) goto 1200
c 
        n=21
        shft=3
        shft=shft/2
        call polev(p0to3(2),x-shft,n-1,fi0)
        fi0=fi0*dexp(x)
        return
 1200 continue
c 
c        . . . for  3< x <6
c 
        if(x .gt. 6) goto 1400
c 
        n=19
        shft=3+6
        shft=shft/2
c 
        call polev(p3to6(2),x-shft,n-1,fi0)
        fi0=fi0*dexp(x)
        return
 1400 continue
c 
c        . . . for  6< x <9
c 
        if(x .gt. 9) goto 1600
c 
        n=13
        xinv=done/x
        call polev(p6to9(2),xinv,n-1,fi0)
        fi0=fi0*dexp(x)*dsqrt(xinv)
        return
 1600 continue
c 
c        . . . for  9< x <12
c 
        if(x .gt. 12) goto 1800
c 
        n=10
        xinv=done/x
        call polev(p9to12(2),xinv,n-1,fi0)
        fi0=fi0*dexp(x)*dsqrt(xinv)
        return
 1800 continue
c 
c        . . . for  12< x <24
c 
        if(x .gt. 24) goto 2000
c 
        n=12
        xinv=done/x
        call polev(p12to24(2),xinv,n-1,fi0)
        fi0=fi0*dexp(x)*dsqrt(xinv)
        return
 2000 continue
c 
c        . . . for  24< x <48
c 
        if(x .gt. 48) goto 2200
c 
        n=8
        xinv=done/x
        call polev(p24to48(2),xinv,n-1,fi0)
        fi0=fi0*dexp(x)*dsqrt(xinv)
        return
 2200 continue
c 
c        . . . for  48< x <96
c 
        if(x .gt. 96) goto 2400
c 
        n=8
        xinv=done/x
        call polev(p48to96(2),xinv,n-1,fi0)
        fi0=fi0*dexp(x)*dsqrt(xinv)
        return
 2400 continue
c 
c        . . . for  48< x <96
c 
        if(x .gt. 192) goto 2600
c 
        n=6
        xinv=done/x
        call polev(p96to192(2),xinv,n-1,fi0)
        fi0=fi0*dexp(x)*dsqrt(xinv)
        return
 2600 continue
c 
c        for x>192
c 
        d=done+done/(8*x)+9*done/(2*(8*x)**2)
     1      +9*25*done/(6*(8*x)**3)
     2     + 9*25*49*done/(24*(8*x)**4 )
     3      +9*25*49*81/(120*(8*x)**5)
     4      +9*25*49*81*121/(720*(8*x)**6)
  
c 
        pi=datan(done)*4
cccc         d=done
        fi0=dexp(x)/dsqrt(2*pi*x)*d
        return
        end
c 
c 
c 
c 
c 
        subroutine polev(p,x,n,ff)
        implicit real *8 (a-h,o-z)
        save
        dimension p(1)
        ff=0
        d=1
        do 1200 i=0,n
        ff=ff+d*p(i)
        d=d*x
 1200 continue
        return
        end
  
