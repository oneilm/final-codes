       implicit real *8 (a-h,o-z)
        dimension p(100),a(10000),err(10000),
     1     b(10000),erry(10000),err1(1000),erry1(1000)
c 
        CALL PRINI(6,13)
C 
c       . . . construct the testing parameters
c 
c 
         PRINT *, 'ENTER t'
         READ *,t
         CALL PRIN2('t=*',t,1)
c 
        call jy0rea(t,fj,fy)
        call prin2('after jy0rea, fj=*',fj,1)
        call prin2('after jy0rea, fy=*',fy,1)
c 
c       test the subroutine against the direct one
c 
        nnn=1000
        eps = 1.0d-18
        h=0.25
ccc        h=h/nnn
        dj=0
        dy=0
        dj1=0
        dy1=0
        ddj=0
        ddy=0
        dy1=0
        ddy1=0
        do 1400 i=1,nnn
        xx=(i-1)*h+0.001
c 
cccc       call jy0rea(xx,fj,fy)
        call jy01rea(xx,fj,fy,fj1,fy1)
cccc        call j01rea(xx,fj,fj1)
cccc        call j0rea(xx,fj)
        call jbq(xx,eps,a(2),n)
        err(i)=a(1)-fj
c 
        err1(i)=a(2)-fj1
        if( (i/20)*20 .eq. i)
     1      call prin2('xx=*',xx,1)
ccc         call prin2('fy1=*',fy1,1)
c 
        call ybq(xx,y0,a(2),b(2),y1)
ccc         call prin2('y1=*',y1,1)
ccc         call prin2('fj1=*',fj1,1)
ccc         call prin2('a(2)=*',a(2),1)
        erry(i)=fy-y0
c 
          dy1=dy1+(fy1-y1)**2
          ddy1=ddy1+y1**2
        erry1(i)=fy1-y1
ccc           call prin2('xx=*',xx,1)
ccc           call prin2('fj1=*',fj1,1)
ccc           call prin2('a(2)=*',a(2),1)
ccc           call prin2('fy1=*',fy1,1)
ccc           call prin2('a(2)=*',a(2),1)
  
        if(dj .lt. dabs(fj-a(1)) ) xj=xx
        if(dy .lt. dabs(fy-y0) ) xy=xx
  
        if(dj .lt. dabs(fj-a(1)) ) dj=dabs(fj-a(1))
        if(dy .lt. dabs(fy-y0) ) dy=dabs(fy-y0)
c 
        if(dj1 .lt. dabs(err1(i) ) ) dj1=dabs(err1(i) )
 1400 continue
        call  prin2('err=*',err,nnn)
        call  prin2('erry=*',erry,nnn)
        call  prin2('err1=*',err1,nnn)
        call  prin2('erry1=*',erry1,nnn)
  
        call prin2('dj=*',dj,1)
        call prin2('dy=*',dy,1)
        call prin2('dj1=*',dj1,1)
  
        call prin2('xj=*',xj,1)
        call prin2('xy=*',xy,1)
          call prin2('finally, l2 error in y1 is*',
     1        dsqrt(dy1/ddy1),1)
          if(2 .ne. 3) stop
c 
        nnn=10000
        t1=second(jt)
        do 2400 i=1,nnn
         if(i .ge. nnn/4) h=5.2
        xx=0.001+(i-1)*h
cccc         xx=4
c 
        call jy0rea(xx,fj,fy)
ccc        call j0rea(xx,fj)
 2400 continue
        t2=second(jt)
  
        call  prin2('t2-t1=*',t2-t1,1)
        stop
        end
C 
C 
C 
C 
c 
        subroutine jbq(x,eps,fjs,n)
        implicit real *8 (a-h,o-z)
        save
        dimension fjs(1)
c 
c       determine how far we have to go
c 
        i0=0
        done=1
        fjs(i0)=1
        fjs(1)=0
        dd=1.0d50
c 
        do 1200 i=1,100000
         nn=i
        fjs(i+1)=(2*i/x)*fjs(i)-fjs(i-1)
        if( dabs(fjs(i)) .ge. dd) goto 1400
 1200 continue
 1400 continue
ccccc          call prinf('nn=*',nn,1)
c 
c        conduct recursion down
c 
        fjs(nn)=0
        fjs(nn-1)=1
        do 2200 i=nn-1,1,-1
        fjs(i-1)=(2*i/x)*fjs(i)-fjs(i+1)
 2200 continue
c 
c        determine the scaling coefficients
c 
        ddd=0
        do 2400 i=1,nn
        ddd=ddd+fjs(i)**2
 2400 continue
        ddd=2*ddd+fjs(i0)**2
        ddd=done/dsqrt(ddd)
c 
c       . . . scale
c 
        do 2600 i=0,nn
        fjs(i)=fjs(i)*ddd
        if(dabs(fjs(i)) .ge. eps) n=i
 2600 continue
        n=n+2
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine ybq(x,y0,fjs,ders,y1)
        implicit real *8 (a-h,o-z)
        save
        dimension fjs(1),ders(1000)
        data gamma/0.5772156649015328606d+00/
c 
c       calculate the j-functions
c 
        eps=1.0d-24
        call jbq(x,eps,fjs,n)
c 
c       calculate y0
c 
        done=1
        pi=datan(done)*4
        i0=0
c 
ccc        call prinq('in ybq, fjs(0)=*',fjs(i0),n)
ccc        call prinq('in ybq, gamma=*',gamma,1)
        y0=dlog(x/2)+gamma
        y0=2/pi*y0*fjs(i0)
c 
c        evaluate y0
c 
        dd=0
        d=-1
        do 1200 i=1,n/2
        dd=dd+d*fjs(2*i)/i
        d=-d
 1200 continue
ccc         call prinq('y0 before addition =*',y0,1)
          dd=4*dd/pi
ccc         call prinq('dd=*',dd,1)
c 
         four=4
        y0=y0-dd
c 
ccc         call prinq('inside ybd, y0=*',y0,1)
c 
c        . . . calculate the derivatives of j-functions
c 
        do 1400 i=1,n
        ders(i)=(fjs(i-1)-fjs(i+1))/2
 1400 continue
c 
c       evaluate y1
c 
        y1=dlog(x/2)+gamma
        y1=-y1*fjs(1)*2/pi
        y1=y1+fjs(i0)/(pi*x)*2
c 
        dd=0
        d=-1
        do 1600 i=1,n/2
        dd=dd+d*ders(2*i)/i
        d=-d
 1600 continue
c 
        y1=y1-4/pi*dd
        y1=-y1
ccc        call prinq('in ybq, y1=*',y1,1)
        return
        end
c 
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        this is the end of the debugging routines and the
c        beginning of the actual code for the evaluation of
c        bessel functions j0, j1, y0, y1 of a real argument
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c 
c 
        subroutine jy0rea(x,fj,fy)
        implicit real *8 (a-h,o-z)
        save
        dimension p1(10),p2(10),p3(10)
        dimension q1(10),q2(10),q3(10),p01(14),pa(14),
     1      pcoe(4),qcoe(3),ppcoe(4),qqcoe(3)
  
        dimension pp1(10),pp2(10),pp3(10)
        dimension qq1(10),qq2(10),qq3(10),pp01(14),ppa(14)
c 
c        this subroutine evaluates the bessel functions
c        j0, y0, j1, j1  of the real argument x. it is a double
c        precision code, and it produces better than
c        14-digit accuracy.
c 
c        the subroutine has 3 entries:
c 
c   1. The main entry jy0rea evaluates the functions j0, y0.
c 
c   2. the second entry j0rea evaluates the function j0
c 
c   3. the third entry jy01rea evaluates all four functions
c        j0, y0, j1, y1.
c 
c   4. tha last entry j01rea evaluates the functions j0, j1.
c 
c    note:
c        as a matter of fact, the difference in speed between
c        these entries is minor, except when the argument x
c        is smaller than 5. in most cases, it is cosher to
c        call jy0rea and use whatever functions are needed.
c 
c    see individual entries for the description  of parameters
c 
c 
c            input parameters (for the main entry jy0rea):
c 
c  x - the point at which the bessel functions are to be
c        computed
c 
c            output parameters:
c 
c  fj - j_0(x)
c  fy - y_0(x)
c 
        data n/9/,done/1.0d0/,
     1        pi/0.31415926535897932d+01/,
     2      pi4/ 0.78539816339744831d+00/,
     3      twopi/0.63661977236758134d+00/,half/0.5d0/
  
        data gamma/0.5772156649015328606d+00/
c 
c        following are the data for the evaluation of j0, y1
c 
c       polynomials for the interval[  5,  8]
c       coefficients u, v are
c 
        data u1/-.26666666666666666D+02/,
     1       v1/0.43333333333333333D+01/
        data p1/
     1   0.55188630535417531D+00,0.29842611342375205D-02,
     2   -.22253610298042293D-04,-.24059562114410511D-05,
     3   -.62519562186941359D-07,0.24934055666238108D-08,
     4   0.41206115755368409D-09,0.25528739624711840D-10,
     5   0.25523813472930658D-12,-.14105831621855898D-12/
        data q1/
     1   -.57447651351824614D+00,0.20863543714478222D-02,
     2   0.67640390399059588D-04,0.21557363995613308D-06,
     3   -.93127374605599322D-07,-.59030870188416595D-08,
     4   -.10310872406889678D-09,0.18457676936265269D-10,
     5   0.25407175781017002D-11,0.17219044905026530D-12/
c 
c       polynomials for the interval[  8, 16]
c       coefficients u, v are
c 
       data u2/-.32000000000000000D+02/,
     1      v2/ 0.30000000000000000D+01/
        data p2/
     1   0.55726719109133432D+00,0.23975519140921681D-02,
     2   -.25356091519412556D-04,-.15741000190497896D-05,
     3   -.66082369527225043D-08,0.34568315096067513D-08,
     4   0.21202642743659065D-09,-.21382820333999236D-11,
     5   -.15199234054449315D-11,-.13265185344533817D-12/
        data q2/
     1   -.57042403185595782D+00,0.19449230725145111D-02,
     2   0.46234660955849426D-04,-.38639510641770226D-06,
     3   -.70414485151628035D-07,-.18596728814008043D-08,
     4   0.15098408425498177D-09,0.19693743665260420D-10,
     5   0.69185850796570623D-12,-.81099121304145928D-13/
c 
c       polynomials for the interval[ 16,200]
c       coefficients u, v are
c 
       data u3/-.33032258064516129D+02/,
     1      v3/0.10645161290322580D+01/
        data p3/
     1   0.56187709478339232D+00,0.22082674022637160D-02,
     2   -.32377039027116634D-04,-.13310422741973540D-05,
     3   0.32696371624251529D-07,0.41987947660104233D-08,
     4   -.38698138851351260D-10,-.28629742386609819D-10,
     5   -.89701868091845451D-12,0.26904152943996875D-12/
        data q3/
     1   -.56641981030146299D+00,0.20539683046510319D-02,
     2   0.39622470761182690D-04,-.89007635266569935D-06,
     3   -.65556699881311379D-07,0.13238721366291133D-08,
     4   0.32440195816636586D-09,0.27845179856889703D-11,
     5   -.26876647471195694D-11,-.17902792592772566D-12/
c 
c       polynomial for the interval[ -5,  5]
c       approximating j0
c       coefficients u, v are
c 
        data u01/0.20000D+00/,v01/0.00D+01/
        data p01/
     1   0.99999999999999999D+00,-.62499999999999999D+01,
     2   0.97656249999999971D+01,-.67816840277777046D+01,
     3   0.26490953233497169D+01,-.66227383082984269D+00,
     4   0.11497809559088356D+00,-.14665573279047978D-01,
     5   0.14321845520995797D-02,-.11050748338748342D-03,
     6   0.69060109259921222D-05,-.35612800636151098D-06,
     7   0.15137695102037179D-07,-.45766751791858648D-09/
c 
c       polynomial for the interval[ -5,  5]
c       approximating the auxiliary series
c       for y0.
c 
c       coefficients u, v are
c 
        data ua/0.200D+00/,va/0.00D+01/
        data pa/
     1   0.23080374166345875D-18,0.39788735772973833D+01,
     2   -.93254849467907364D+01,0.79151492603932262D+01,
     3   -.35134717952720334D+01,0.96269127188894342D+00,
     4   -.17933345613328207D+00,0.24207935611640996D-01,
     5   -.24780251676542572D-02,0.19902129007039959D-03,
     6   -.12877070467925354D-04,0.68457290371870522D-06,
     7   -.29874878428850324D-07,0.92184735587807664D-09/
c 
c        . . . and the polynomials for the far-field
c              expansion
c 
         data  pcoe/
     1  0.100000d+01, -0.703125000d-01,
     2  0.112152099609375000d+00,-0.572501420974731445d+00/
         data qcoe /
     1  -0.125000000000000000d+00,0.732421875000000000d-01,
     2    -0.227108001708984375d+00/
c 
c        following are the data for the evaluation of j1, y1
c 
c       . . . polynomials for the interval[  5,  8]
c             coefficients u, v are
c 
         data uu1/ -.26666666666666666D+02/,
     1     vv1/0.43333333333333333D+01/
c 
        data pp1/
     1   -.53173716842357832D+00,-.70178172669864168D-02,
     2   -.11141466843623608D-03,-.63344921233317975D-07,
     3   0.12752038123508872D-06,0.71606938492205609D-08,
     4   0.93727788222353480D-10,-.23557767481351388D-10,
     5   -.29624451496572487D-11,-.18816170183253387D-12/
        data qq1/
     1   -.60003166328076073D+00,0.85394947411512136D-02,
     2   -.46097186926304210D-04,-.34453219531793409D-05,
     3   -.73622918800211383D-07,0.35853041907877011D-08,
     4   0.50514959096077148D-09,0.28823625074304308D-10,
     5   0.12955236601959938D-12,-.17206089234769226D-12/
c 
c       polynomials for the interval[  8, 16]
c       coefficients u, v are
c 
        data uu2/-.3200d+02/,vv2/0.300000D+01/
        data pp2/
     1   -.54497644374938988D+00,-.61861877569424901D-02,
     2   -.75519539239281760D-04,0.64276660968355753D-06,
     3   0.91625624280644285D-07,0.21146401912860349D-08,
     4   -.19056156225658613D-09,-.22870252895459067D-10,
     5   -.72685987817955162D-12,0.97310924711447901D-13/
        data qq2/
     1   -.58455282719872705D+00,0.69451579319909565D-02,
     2   -.45980909043254114D-04,-.21942261041953223D-05,
     3   -.46261694080917413D-08,0.43598630512244850D-08,
     4   0.24573711261656784D-09,-.34057690674508506D-11,
     5   -.17729655060214402D-11,-.14665688721652535D-12/
c 
c       polynomials for the interval[ 16,200]
c       coefficients u, v are
c 
         data uu3/-.33032258064516129D+02/,
     1     vv3/  0.10645161290322580D+01/
        data pp3/
     1   -.55744188279451250D+00,-.62708736488248012D-02,
     2   -.65213613041335236D-04,0.12770306722443576D-05,
     3   0.83712703670230413D-07,-.17016400001870425D-08,
     4   -.38327669173427131D-09,-.27018853980163044D-11,
     5   0.30698129249049068D-11,0.19617321436013121D-12/
        data qq3/
     1   -.57107443930370983D+00,0.65282322527159425D-02,
     2   -.55056346423860902D-04,-.18451794882073581D-05,
     3   0.43428771082098162D-07,0.51132086696947351D-08,
     4   -.51845748226931597D-10,-.33144493956423250D-10,
     5   -.97340298914856351D-12,0.30457057625214024D-12/
c 
c       polynomial for the interval [-5,  5]
c       approximating j1
c 
        data pp01/
     1   0.24999999999999999D+01,-.78124999999999999D+01,
     2   0.81380208333333328D+01,-.42385525173610987D+01,
     3   0.13245476616751823D+01,-.27594742951437636D+00,
     4   0.41063605575835601D-01,-.45829916703240591D-02,
     5   0.39782908072775114D-03,-.27626920623275312D-04,
     6   0.15695917451666021D-05,-.74218359709148273D-07,
     7   0.29194624036124853D-08,-.82969140436250033D-10/
c 
c       polynomial for the interval [-5,5]
c       approximating the auxiliary series
c       for y1
c 
        data ppa/
     1   0.15915494309189533D+01,-.74603879574325938D+01,
     2   0.94981791124720513D+01,-.56215548724384520D+01,
     3   0.19253825438098315D+01,-.43040029491828996D+00,
     4   0.67782220525766800D-01,-.79296828112744362D-02,
     5   0.71648105630007065D-03,-.51514211657500308D-04,
     6   0.30175384873850692D-05,-.14660942256048710D-06,
     7   0.59057445601699370D-08,-.17100368998812547D-09/
c 
c        . . . and the polynomials for the far-field
c              expansion
c 
         data  ppcoe/
     1  0.10d+01,0.1171875000d+00,
     2  -0.1441955566406250d+00,0.676592588424682617d+00/
         data qqcoe /
     1   0.37500d+00,  -0.102539062500d+00,
     2   0.277576446533203125d+00/
c 
c        if x is less than 5 - use local expansion
c 
        if(x .lt. 5.0d0) goto 2000
c 
c       if x is greater than 200 - use far-field expansion
c 
        if(x .gt. 200) goto 2200
c 
c        x is between 5 and 200 - use the main regime
c 
        if(x .gt. 8.0d0) goto 1400
c 
c       . . . x is between 5 and 8 - use first expansion
c 
        ttt=done/x
        t=u1*ttt+v1
        call polev2(p1(2),q1(2),t,n,f1,f2)
c 
        xcos=dcos(x)
        xsin=dsin(x)
        tsq=dsqrt(ttt)
        fj=xcos*f1-xsin*f2
        fj=fj*tsq
        fy=xcos*f2+xsin*f1
        fy=fy*tsq
         return
 1400 continue
c 
        if(x .gt. 16.0d0) goto 1600
c 
c       . . . x is between 8 and 16 - use second expansion
c 
        ttt=done/x
        t=u2*ttt+v2
        call polev2(p2(2),q2(2),t,n,f1,f2)
c 
        xcos=dcos(x)
        xsin=dsin(x)
        tsq=dsqrt(ttt)
        fj=xcos*f1-xsin*f2
        fj=fj*tsq
        fy=xcos*f2+xsin*f1
        fy=fy*tsq
         return
 1600 continue
c 
c       . . . x is between 16 and 200 - use third expansion
c 
        ttt=done/x
        t=u3*ttt+v3
        call polev2(p3(2),q3(2),t,n,f1,f2)
c 
        xcos=dcos(x)
        xsin=dsin(x)
        tsq=dsqrt(ttt)
        fj=xcos*f1-xsin*f2
        fj=fj*tsq
        fy=xcos*f2+xsin*f1
        fy=fy*tsq
         return
 2000 continue
c 
c        x is on the interval [0,5]. use the local
c        expansion
c 
        t=u01*x+v01
        t=t*t
        n0=13
        call polev2(p01(2),pa(2),t,n0,fj,ser)
ccc        call polev(p01(2),t,n0,fj)
ccc        call polev(pa(2),t,n0,ser)
        fy=dlog(x*half)+gamma
        fy=fy*fj*twopi +ser
        return
c 
 2200 continue
c 
c       x is greater than 200. use far-field expansion
c 
        t=done/x
        tt=t**2
        pd=((pcoe(4)*tt+pcoe(3))*tt+pcoe(2))*tt+pcoe(1)
        qd=(qcoe(3)*tt+qcoe(2))*tt+qcoe(1)
        qd=qd*t
        ttt=x-pi4
c 
        xcos=dcos(ttt)
        xsin=dsin(ttt)
        tsq=dsqrt(twopi*t)
        fj=xcos*pd-xsin*qd
        fy=xcos*qd+xsin*pd
        fj=fj*tsq
        fy=fy*tsq
c 
        return
c 
c 
c 
c 
        entry j0rea(x,fj)
c 
c        this entry evaluates the bessel function
c        j0 of the real argument x. it is a double
c        precision code, and it produces better than
c        14-digit accuracy.
c 
c            input parameters:
c 
c  x - the point at which the bessel functions are to be
c        computed
c 
c            output parameters:
c 
c  fj - j_0(x)
c 
c        if x is less than 5 - use local expansion
c 
        if(x .lt. 5.0d0) goto 4000
c 
c       if x is greater than 200 - use far-field expansion
c 
        if(x .gt. 200) goto 4200
c 
c        x is between 5 and 100 - use the main regime
c 
        if(x .gt. 8.0d0) goto 3400
c 
c       . . . x is between 5 and 8 - use first expansion
c 
        ttt=done/x
        t=u1*ttt+v1
        call polev2(p1(2),q1(2),t,n,f1,f2)
c 
        xcos=dcos(x)
        xsin=dsin(x)
        tsq=dsqrt(ttt)
        fj=xcos*f1-xsin*f2
        fj=fj*tsq
ccc        fy=xcos*f2+xsin*f1
ccc        fy=fy*tsq
         return
 3400 continue
c 
        if(x .gt. 16.0d0) goto 3600
c 
c       . . . x is between 8 and 16 - use second expansion
c 
        ttt=done/x
        t=u2*ttt+v2
        call polev2(p2(2),q2(2),t,n,f1,f2)
c 
        xcos=dcos(x)
        xsin=dsin(x)
        tsq=dsqrt(ttt)
        fj=xcos*f1-xsin*f2
        fj=fj*tsq
         return
 3600 continue
c 
c       . . . x is between 16 and 200 - use third expansion
c 
        ttt=done/x
        t=u3*ttt+v3
        call polev2(p3(2),q3(2),t,n,f1,f2)
c 
        xcos=dcos(x)
        xsin=dsin(x)
        tsq=dsqrt(ttt)
        fj=xcos*f1-xsin*f2
        fj=fj*tsq
         return
 4000 continue
c 
c        x is on the interval [0,5]. use the local
c        expansion
c 
        t=u01*x+v01
        t=t*t
        n0=13
        call polev(p01(2),t,n0,fj)
        return
 4200 continue
c 
c       x is greater than 200. use far-field expansion
c 
        t=done/x
        tt=t**2
        pd=((pcoe(4)*tt+pcoe(3))*tt+pcoe(2))*tt+pcoe(1)
        qd=(qcoe(3)*tt+qcoe(2))*tt+qcoe(1)
        qd=qd*t
        ttt=x-pi4
c 
        xcos=dcos(ttt)
        xsin=dsin(ttt)
        tsq=dsqrt(twopi*t)
        fj=xcos*pd-xsin*qd
        fj=fj*tsq
        return
c 
c 
c 
c 
        entry jy01rea(x,fj,fy,fj1,fy1)
c 
c        this entry evaluates the bessel functions
c        j0, y0, j1, y1  of the real argument x. it is a double
c        precision code, and it produces better than
c        14-digit accuracy.
c 
c            input parameters:
c 
c  x - the point at which the bessel functions are to be
c        computed
c 
c            output parameters:
c 
c  fj - j_0(x)
c  fy - y_0(x)
c  fj1 - j_1(x)
c  fy1 - y_1(x)
c 
c 
c        if x is less than 5 - use local expansion
c 
        if(x .lt. 5.0d0) goto 5000
c 
c       if x is greater than 200 - use far-field expansion
c 
        if(x .gt. 200) goto 5200
c 
c        x is between 5 and 200 - use the main regime
c 
        if(x .gt. 8.0d0) goto 4400
c 
c       . . . x is between 5 and 8 - use first expansion
c 
        ttt=done/x
        t=u1*ttt+v1
        call polev4(p1(2),q1(2),pp1(2),qq1(2),
     1      t,n,f1,f2,ff1,ff2)
c 
        xcos=dcos(x)
        xsin=dsin(x)
        tsq=dsqrt(ttt)
        fj=xcos*f1-xsin*f2
        fj=fj*tsq
        fy=xcos*f2+xsin*f1
        fy=fy*tsq
c 
        fj1=xcos*ff1-xsin*ff2
        fj1=fj1*tsq
        fy1=xcos*ff2+xsin*ff1
        fy1=fy1*tsq
c 
         return
 4400 continue
c 
        if(x .gt. 16.0d0) goto 4600
c 
c       . . . x is between 8 and 16 - use second expansion
c 
        ttt=done/x
        t=u2*ttt+v2
c 
        call polev4(p2(2),q2(2),pp2(2),qq2(2),
     1      t,n,f1,f2,ff1,ff2)
c 
        xcos=dcos(x)
        xsin=dsin(x)
        tsq=dsqrt(ttt)
        fj=xcos*f1-xsin*f2
        fj=fj*tsq
        fy=xcos*f2+xsin*f1
        fy=fy*tsq
c 
        fj1=xcos*ff1-xsin*ff2
        fj1=fj1*tsq
        fy1=xcos*ff2+xsin*ff1
        fy1=fy1*tsq
c 
         return
 4600 continue
c 
c       . . . x is between 16 and 200 - use third expansion
c 
        ttt=done/x
        t=u3*ttt+v3
c 
        call polev4(p3(2),q3(2),pp3(2),qq3(2),
     1      t,n,f1,f2,ff1,ff2)
c 
        xcos=dcos(x)
        xsin=dsin(x)
        tsq=dsqrt(ttt)
        fj=xcos*f1-xsin*f2
        fj=fj*tsq
        fy=xcos*f2+xsin*f1
        fy=fy*tsq
c 
        fj1=xcos*ff1-xsin*ff2
        fj1=fj1*tsq
        fy1=xcos*ff2+xsin*ff1
        fy1=fy1*tsq
c 
         return
 5000 continue
c 
c        x is on the interval [0,5]. use the local
c        expansion
c 
        ttt=u01*x+v01
        t=ttt*ttt
        n0=13
        call polev4(p01(2),pa(2),pp01(2),ppa(2),
     1      t,n0,fj,ser,fj1,ser2)
        fy=dlog(x*half)+gamma
         ddd=fy
        fy=fy*fj*twopi +ser
        fj1=fj1*ttt
  
c 
         ser2=ser2*ttt
       fy1=(ddd*fj1-fj/x)*twopi-ser2
        return
 5200 continue
c 
c       x is greater than 200. use far-field expansion
c 
        t=done/x
        tt=t**2
        pd=((pcoe(4)*tt+pcoe(3))*tt+pcoe(2))*tt+pcoe(1)
        qd=(qcoe(3)*tt+qcoe(2))*tt+qcoe(1)
        qd=qd*t
        ttt=x-pi4
c 
        xcos=dcos(ttt)
        xsin=dsin(ttt)
        tsq=dsqrt(twopi*t)
        fj=xcos*pd-xsin*qd
        fy=xcos*qd+xsin*pd
        fj=fj*tsq
        fy=fy*tsq
c 
        ppd=((ppcoe(4)*tt+ppcoe(3))*tt+ppcoe(2))*tt+ppcoe(1)
        qqd=(qqcoe(3)*tt+qqcoe(2))*tt+qqcoe(1)
        qqd=qqd*t
c 
        fy1=xcos*ppd-xsin*qqd
        fj1=xcos*qqd+xsin*ppd
        fj1=fj1*tsq
        fy1=-fy1*tsq
        return
c 
c 
c 
c 
        entry j01rea(x,fj,fj1)
c 
c        this entry evaluates the bessel functions
c        j0, j1 of the real argument x. it is a double
c        precision code, and it produces better than
c        14-digit accuracy.
c 
c            input parameters:
c 
c  x - the point at which the bessel functions are to be
c        computed
c 
c            output parameters:
c 
c  fj - j_0(x)
c  fj1 - j_1(x)
c 
c 
c        if x is less than 5 - use local expansion
c 
        if(x .lt. 5.0d0) goto 6000
c 
c       if x is greater than 200 - use far-field expansion
c 
        if(x .gt. 200) goto 6200
c 
c        x is between 5 and 200 - use the main regime
c 
        if(x .gt. 8.0d0) goto 5400
c 
c       . . . x is between 5 and 8 - use first expansion
c 
        ttt=done/x
        t=u1*ttt+v1
        call polev4(p1(2),q1(2),pp1(2),qq1(2),
     1      t,n,f1,f2,ff1,ff2)
c 
        xcos=dcos(x)
        xsin=dsin(x)
        tsq=dsqrt(ttt)
        fj=xcos*f1-xsin*f2
        fj=fj*tsq
c 
        fj1=xcos*ff1-xsin*ff2
        fj1=fj1*tsq
         return
 5400 continue
c 
        if(x .gt. 16.0d0) goto 5600
c 
c       . . . x is between 8 and 16 - use second expansion
c 
        ttt=done/x
        t=u2*ttt+v2
        call polev4(p2(2),q2(2),pp2(2),qq2(2),
     1      t,n,f1,f2,ff1,ff2)
c 
        xcos=dcos(x)
        xsin=dsin(x)
        tsq=dsqrt(ttt)
        fj=xcos*f1-xsin*f2
        fj=fj*tsq
c 
        fj1=xcos*ff1-xsin*ff2
        fj1=fj1*tsq
         return
 5600 continue
c 
c       . . . x is between 16 and 200 - use third expansion
c 
        ttt=done/x
        t=u3*ttt+v3
        call polev4(p3(2),q3(2),pp3(2),qq3(2),
     1      t,n,f1,f2,ff1,ff2)
c 
        xcos=dcos(x)
        xsin=dsin(x)
        tsq=dsqrt(ttt)
        fj=xcos*f1-xsin*f2
        fj=fj*tsq
c 
        fj1=xcos*ff1-xsin*ff2
        fj1=fj1*tsq
         return
 6000 continue
c 
c        x is on the interval [0,5]. use the local
c        expansion
c 
        ttt=u01*x+v01
        t=ttt*ttt
        n0=13
        call polev2(p01(2),pp01(2),t,n0,fj,fj1)
        fj1=fj1*ttt
        return
 6200 continue
c 
c       x is greater than 200. use far-field expansion
c 
        t=done/x
        tt=t**2
        pd=((pcoe(4)*tt+pcoe(3))*tt+pcoe(2))*tt+pcoe(1)
        qd=(qcoe(3)*tt+qcoe(2))*tt+qcoe(1)
        qd=qd*t
        ttt=x-pi4
c 
        xcos=dcos(ttt)
        xsin=dsin(ttt)
        tsq=dsqrt(twopi*t)
        fj=xcos*pd-xsin*qd
        fj=fj*tsq
c 
        ppd=((ppcoe(4)*tt+ppcoe(3))*tt+ppcoe(2))*tt+ppcoe(1)
        qqd=(qqcoe(3)*tt+qqcoe(2))*tt+qqcoe(1)
        qqd=qqd*t
c 
        fj1=xcos*qqd+xsin*ppd
        fj1=fj1*tsq
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
c 
c 
c 
c 
c 
        subroutine polev2(p,q,x,n,f1,f2)
        implicit real *8 (a-h,o-z)
        save
        dimension p(1),q(1)
        f1=0
        f2=0
        d=1
        do 1200 i=0,n
        f1=f1+d*p(i)
        f2=f2+d*q(i)
        d=d*x
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine polev4(p,q,pp,qq,x,n,f1,f2,ff1,ff2)
        implicit real *8 (a-h,o-z)
        save
        dimension p(1),q(1),pp(1),qq(1)
c 
        f1=0
        f2=0
        ff1=0
        ff2=0
        d=1
        do 1200 i=0,n
        f1=f1+d*p(i)
        f2=f2+d*q(i)
c 
        ff1=ff1+d*pp(i)
        ff2=ff2+d*qq(i)
        d=d*x
 1200 continue
        return
        end
  
  
  
  
