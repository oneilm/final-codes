         implicit real *8 (a-h,o-z)
         dimension points(10000),weights(10000)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER k'
        READ *,k
        CALL PRINf('k=*',k,1 )
c 
c       construct the test arrays to be summed by parts
c 
        nn=10
        n=10
        call qua10test(k,nn,n,points,weights)
  
  
  
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine qua10test(k,nn,n,points,weights)
        implicit real *8 (a-h,o-z)
        save
        dimension xgauss(100),points(nn,n),weights(nn,n),
     1      rints(100),rints2(100),diffs(100)
c 
c       construct the Gaussian nodes
c 
  
  
        itype=0
        call legeexps(itype,n,xgauss,u,v,whts)
  
        call prin2('xgauss=*',xgauss,n)
  
        if(nn .eq. 20) call qu10by20(points,weights)
        if(nn .eq. 10) call qu10by10(points,weights)
  
c 
        call prin2('points(1,k) as retrieved *',points(1,k),nn)
        call prin2('weights(1,k) as retrieved *',weights(1,k),nn)
  
c 
c       integrate the first several polynomials via each of the
c       columns
c 
  
  
        two=2
  
  
        do 2200 i=1,n
c 
        rints(i)=0
 2200 continue
  
  
  
  
c 
c       integrate the first several polynomials multiplied by logs
c       via each of the columns
c 
  
  
        two=2
  
        iord=1
  
        do 3200 i=1,n
c 
        rints(i)=0
 3200 continue
  
c 
        do 3600 i=1,n
c 
        do 3400 j=1,nn
c 
        rints(i)=rints(i)+weights(j,i)*
     1           fun2(points(j,i),xgauss(i),iord)
  
cccc     1      log(abs(xgauss(i)-points(j,i)) )
 3400 continue
  
c 
c        evaluate the same integral via adapgaus
c 
        a=-1
        b=1
        m=20
        eps=1.0d-14
c 
        call adapgaus(ier,a,b,fun2,xgauss(i),iord,m,eps,
     1      rint2,maxrec,numint)
  
        rints2(i)=rint2
  
        diffs(i)=rints2(i)-rints(i)
 3600 continue
  
        call prin2('and integrals are*',rints,n)
        call prin2('and rints2=*',rints2,n)
        call prin2('and diffs=*',diffs,n)
c 
  
  
  
        return
        end
c 
c 
c 
c 
c 
        function fun2(point,xgauss,iord)
        implicit real *8 (a-h,o-z)
        save
  
        fun2=abs(point-xgauss)**8
cccc        fun2=log(abs(point-xgauss))*point**7
  
cccc        fun2=1
        return
        end
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        This is the end of the debugging code and the beginning of
c        the quadrature code proper
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c        This file contains four user-callable subroutines:
c 
c 
c  qu10by20 - returns to the user the nodes and weights
c       of 10 "universal" quadratures on the interval [-1,1],
c       corresponding to the 10 Gaussian nodes. The i-th quadrature
c       integrates exactly all functions f of the form
c 
c         f(x)=phi(t)+psi(t) * log(|x_i-t|),                     (1)
c 
c       with phi, psi two polynomials of degree 19, and x_i
c       the i-th Gaussian node
c 
c  qu10by10 - returns to the user the nodes and weights
c       of 10 "universal" quadratures on the interval [-1,1],
c       corresponding to the 10 Gaussian nodes. The i-th quadrature
c       integrates exactly all functions f of the form (1) (see above),
c       with phi, psi two polynomials of degree 9, and x_i
c       the i-th Gaussian node
c 
c  qu10near - returns to the user the nodes and weights
c       of a 25-point quadrature evaluating integrals on the interval
c       [-1,1]  of functions f of the form
c 
c        phi(x)+psi(x)*log(|x-x_0),                                (2)
c 
c       with phi, psi polynomials of order 19, and x_0 on the
c       interval [1.002, \infty]
c 
c  qu10nea2 - returns to the user the nodes and weights
c       of a 34-point quadrature evaluating integrals on the interval
c       [-1,1]  of functions f of the form
c 
c        phi(x)+psi(x)*log(|x-x_0),                                (3)
c 
c       with phi, psi polynomials of order 19, and x_0 on the
c       interval [1.0000001, \infty]
c 
c 
c 
c 
c 
        subroutine qu10by20(points,weights)
        implicit real *8 (a-h,o-z)
        save
        dimension points(20,10),weights(20,10)
c 
        dimension pts1(20),wts1(20),pts2(20),wts2(20),
     1      pts3(20),wts3(20),pts4(20),wts4(20),
     2      pts5(20),wts5(20)
c 
c       This subroutine returns to the user the nodes and weights
c       of 10 "universal" quadratures on the interval [-1,1],
c       corresponding to the 10 Gaussian nodes. The i-th quadrature
c       integrates exactly all functions f of the form
c 
c         f(x)=phi(t)+psi(t) * log(|x_i-t|),                     (1)
c 
c       with phi, psi two polynomials of degree 19, and x_i
c       the i-th Gaussian node
c 
c                    Input parameters:    NONE
c 
c                    Ootput parameters:
c 
c  points - the nodes of all 10 quadratures: the i-th column of the
c       matrix points contains the 20 nodes of the quadrature
c       corresponding to the -th node
c  weights - the nodes of all 10 quadratures: the i-th column of the
c       matrix points contains the 20 weights of the quadrature
c       corresponding to the -th node
c 
        data pts1/
     1  -.9981634982072590E+00,-.9915538426155355E+00,
     2  -.9832830597261916E+00,-.9717182247177029E+00,
     3  -.9767808048940554E+00,-.9510745407650473E+00,
     4  -.9076092073466108E+00,-.8383218262873255E+00,
     5  -.7409531298650585E+00,-.6149009548967510E+00,
     6  -.4616959693629515E+00,-.2851700678259385E+00,
     7  -.9137517586662577E-01,0.1117191859035582E+00,
     8  0.3147183396123194E+00,0.5074418050864351E+00,
     9  0.6796547909428153E+00,0.8218291223445154E+00,
     *  0.9258690628914792E+00,0.9857550031566521E+00/
c 
        data wts1/
     1  0.4549530276129239E-02,0.8061972356831065E-02,
     2  0.7846491509221885E-02,0.1020921859922621E-01,
     3  0.4376668743055731E-02,0.3155637007552215E-01,
     4  0.5589853336302420E-01,0.8306764492704936E-01,
     5  0.1117777015165651E+00,0.1400742356811925E+00,
     6  0.1656956764325668E+00,0.1863423599940235E+00,
     7  0.1999108602395829E+00,0.2047012957108593E+00,
     8  0.1995879219279501E+00,0.1841403709262305E+00,
     9  0.1586853343160054E+00,0.1243033740307224E+00,
     *  0.8276342604950238E-01,0.3645101332473923E-01/
c 
        data pts2/
     1  -.9988058420291902E+00,-.9872843436634010E+00,
     2  -.9636524670903349E+00,-.9334253562249570E+00,
     3  -.9032141694291794E+00,-.8791242404522657E+00,
     4  -.8659919778214799E+00,-.8441518232383018E+00,
     5  -.7830122293627457E+00,-.6822223115182355E+00,
     6  -.5443651540245961E+00,-.3741596361213601E+00,
     7  -.1787362969799313E+00,0.3250764940420090E-01,
     8  0.2483970127657333E+00,0.4567153332148303E+00,
     9  0.6450928664129666E+00,0.8019509591260683E+00,
     *  0.9174123181749484E+00,0.9841102587843858E+00/
c 
        data wts2/
     1  0.4613165480267725E-02,0.1828807253690649E-01,
     2  0.2800711607426281E-01,0.3132125713737401E-01,
     3  0.2807689733067758E-01,0.1953894632857636E-01,
     4  0.8989706381844279E-02,0.4084585343113525E-01,
     5  0.8124629229974090E-01,0.1199153310922528E+00,
     6  0.1550170278157623E+00,0.1842029675489165E+00,
     7  0.2050691298166857E+00,0.2155437225112376E+00,
     8  0.2141839011938948E+00,0.2003751573220286E+00,
     9  0.1744347717161011E+00,0.1376217543725583E+00,
     *  0.9206104077797359E-01,0.4064788883180330E-01/
c 
        data pts3/
     1  -.9930122613344544E+00,-.9643941805858567E+00,
     2  -.9175869557574051E+00,-.8596474179143115E+00,
     3  -.7990442705578241E+00,-.7443700669795086E+00,
     4  -.7031684479180525E+00,-.6811221147338119E+00,
     5  -.6579449967122333E+00,-.4893032825170349E+00,
     6  -.5949471708817615E+00,-.3441659269876558E+00,
     7  -.1665388360359957E+00,0.3344207235048554E-01,
     8  0.2434356234005800E+00,0.4498696841532687E+00,
     9  0.6389777503450650E+00,0.7978632869142587E+00,
     *  0.9155180699594585E+00,0.9837258757112682E+00/
c 
        data wts3/
     1  0.1779185047295351E-01,0.3870503130688927E-01,
     2  0.5371120504006760E-01,0.6073467935312451E-01,
     3  0.5901993368079349E-01,0.4905519952687902E-01,
     4  0.3249237024565949E-01,0.1335394647359688E-01,
     5  0.4151626285032320E-01,0.1262522598897083E+00,
     6  0.8451456029785398E-01,0.1628408262182312E+00,
     7  0.1907085688234115E+00,0.2071802235592876E+00,
     8  0.2105274840046563E+00,0.2000282919611018E+00,
     9  0.1760212452196401E+00,0.1399000910241113E+00,
     *  0.9402669113438403E-01,0.4161927891732695E-01/
c 
        data pts4/
     1  -.9903478871133064E+00,-.9504025146897737E+00,
     2  -.8834986023815021E+00,-.7974523551287402E+00,
     3  -.7022255002503293E+00,-.6087194789244770E+00,
     4  -.5275278952351443E+00,-.4677586540799004E+00,
     5  -.4360689210457631E+00,-.4121945474876214E+00,
     6  -.3494226766912416E+00,-.2425993523587580E+00,
     7  -.9646839923921976E-01,0.7921243716755010E-01,
     8  0.2715178194483623E+00,0.4658440358656128E+00,
     9  0.6472213975763011E+00,0.8015601619414561E+00,
     *  0.9168056007307856E+00,0.9839468743284697E+00/
c 
        data wts4/
     1  0.2462513260640958E-01,0.5449201732063149E-01,
     2  0.7799498604905830E-01,0.9241688894090970E-01,
     3  0.9619882322938871E-01,0.8902783806613932E-01,
     4  0.7181973054765565E-01,0.4663017060125383E-01,
     5  0.1794303974049506E-01,0.4061799823409467E-01,
     6  0.8507517518442898E-01,0.1277525783356953E+00,
     7  0.1628510773009285E+00,0.1863323765408475E+00,
     8  0.1958227701928090E+00,0.1903138548150774E+00,
     9  0.1700731513382045E+00,0.1365784674773715E+00,
     *  0.9239595239694548E-01,0.4103797108165556E-01/
c 
        data pts5/
     1  -.9883561797860962E+00,-.9398305159297058E+00,
     2  -.8572399919019391E+00,-.7482086250804680E+00,
     3  -.6228514167093103E+00,-.4928317114329242E+00,
     4  -.3702771193724618E+00,-.2666412108172462E+00,
     5  -.1916083010783278E+00,-.1521937160593462E+00,
     6  -.1233125650067164E+00,-.5257959675044444E-01,
     7  0.2012559739993003E+00,0.5877314311857769E-01,
     8  0.3627988191760868E+00,0.5297121321076324E+00,
     9  0.6878399330187783E+00,0.8237603202215137E+00,
     *  0.9259297297557394E+00,0.9856881498392896E+00/
c 
        data wts5/
     1  0.2974603958509253E-01,0.6657945456889161E-01,
     2  0.9731775484182560E-01,0.1190433988432928E+00,
     3  0.1297088242013777E+00,0.1282900896966494E+00,
     4  0.1148917968875341E+00,0.9074932908233869E-01,
     5  0.5818196361216744E-01,0.2224697059733436E-01,
     6  0.4788826761346366E-01,0.9237500180593535E-01,
     7  0.1541960911507042E+00,0.1287410543031414E+00,
     8  0.1665885274544506E+00,0.1648585116745725E+00,
     9  0.1491408089644011E+00,0.1207592726093190E+00,
     *  0.8212177982524418E-01,0.3657506268226379E-01/
c 
c       build the first columns of matrices points, weights
c 
        do 1600 i=1,20
c 
        points(i,1)=pts1(i)
        points(i,2)=pts2(i)
        points(i,3)=pts3(i)
        points(i,4)=pts4(i)
        points(i,5)=pts5(i)
c 
        weights(i,1)=wts1(i)
        weights(i,2)=wts2(i)
        weights(i,3)=wts3(i)
        weights(i,4)=wts4(i)
        weights(i,5)=wts5(i)
 1600 continue
c 
  
        call prin2('points to start with*',points,20)
  
  
        n=20
        do 1800 i=1,20
c 
        points(i,10)=-points(n-i+1,1)
        weights(i,10)=weights(n-i+1,1)
c 
        points(i,9)=-points(n-i+1,2)
        weights(i,9)=weights(n-i+1,2)
c 
        points(i,8)=-points(n-i+1,3)
        weights(i,8)=weights(n-i+1,3)
c 
        points(i,7)=-points(n-i+1,4)
        weights(i,7)=weights(n-i+1,4)
c 
        points(i,6)=-points(n-i+1,5)
        weights(i,6)=weights(n-i+1,5)
  
 1800 continue
c 
  
  
        do 2000 i=1,10
c 
        call bubble(points(1,i),weights(1,i),n)
  
 2000 continue
  
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine qu10by10(points,weights)
        implicit real *8 (a-h,o-z)
        save
        dimension points(10,10),weights(10,10)
c 
        dimension pts1(10),wts1(10),pts2(10),wts2(10),
     1      pts3(10),wts3(10),pts4(10),wts4(10),
     2      pts5(10),wts5(10)
c 
c       This subroutine returns to the user the nodes and weights
c       of 10 "universal" quadratures on the interval [-1,1],
c       corresponding to the 10 Gaussian nodes. The i-th quadrature
c       integrates exactly all functions f of the form
c 
c         f(x)=phi(t)+psi(t) * log(|x_i-t|),                     (1)
c 
c       with phi, psi two polynomials of degree 9, and x_i
c       the i-th Gaussian node
c 
c                    Input parameters:    NONE
c 
c                    Ootput parameters:
c 
c  points - the nodes of all 10 quadratures: the i-th column of the
c       matrix points contains the 10 nodes of the quadrature
c       corresponding to the -th node
c  weights - the nodes of all 10 quadratures: the i-th column of the
c       matrix points contains the 10 weights of the quadrature
c       corresponding to the -th node
c 
        data pts1/
     1  -.9775265978924164E+00,-.9927130008947402E+00,
     2  -.9628113134850536E+00,-.8746495860203248E+00,
     3  -.6821339364616946E+00,-.3834417428981817E+00,
     4  -.7449652732270578E-02,0.3879589161915740E+00,
     5  0.7291092524068535E+00,0.9461148381878638E+00/
  
c 
        data wts1/
     1  0.8299949884569751E-02,0.1532128141680198E-01,
     2  0.4162022974256374E-01,0.1379116907862554E+00,
     3  0.2476744444947763E+00,0.3448552394449961E+00,
     4  0.3972500342673991E+00,0.3809468640076146E+00,
     5  0.2893572427647569E+00,0.1367630231902661E+00/
c 
        data pts2/
     1  -.9845969809489035E+00,-.9349433881396992E+00,
     2  -.8870413885114848E+00,-.8492563405908384E+00,
     3  -.7041082990212924E+00,-.4259170601972792E+00,
     4  -.5079433144193205E-01,0.3563999933656665E+00,
     5  0.7136441715103631E+00,0.9428798320552685E+00/
c 
        data wts2/
     1  0.3722795870469910E-01,0.5475110055262285E-01,
     2  0.3357646969727819E-01,0.7404497504606667E-01,
     3  0.2152092060261831E+00,0.3352215726635366E+00,
     4  0.4040031892573215E+00,0.3963426329646647E+00,
     5  0.3047432791192180E+00,0.1448796159684093E+00/
c 
        data pts3/
     1  -.9738920297674256E+00,-.8789060022103482E+00,
     2  -.7629886381456456E+00,-.6859819895167676E+00,
     3  -.6062524208127935E+00,-.3913446525037240E+00,
     4  -.5247439421702681E-01,0.3437179523871567E+00,
     5  0.7046916223373251E+00,0.9407467189116827E+00/
c 
        data wts3/
     1  0.6497814625107976E-01,0.1158394942782285E+00,
     2  0.1057742370610585E+00,0.4873777565096340E-01,
     3  0.1418296907147840E+00,0.2845019134462200E+00,
     4  0.3817248076397339E+00,0.3947432139533795E+00,
     5  0.3117853367819897E+00,0.1500853842225628E+00/
c 
        data pts4/
     1  -.9641626350443678E+00,-.8287675882140083E+00,
     2  -.6475108693713116E+00,-.4971588163794493E+00,
     3  -.4148596373849821E+00,-.2300352756653790E+00,
     4  0.8725358899128882E-01,0.4499048353734389E+00,
     5  0.7626246751946268E+00,0.9540485100712768E+00/
c 
        data wts4/
     1  0.8986015426853379E-01,0.1706946146862437E+00,
     2  0.1784284930016252E+00,0.1095051902759451E+00,
     3  0.1011946610928487E+00,0.2632884950700838E+00,
     4  0.3562706727172412E+00,0.3525757971974298E+00,
     5  0.2606065288583410E+00,0.1175753928317077E+00/
c 
        data pts5/
     1  -.9558295208928536E+00,-.7842118273308901E+00,
     2  -.5384583459526444E+00,-.3038725818772138E+00,
     3  -.1617355128130771E+00,-.6242317102246658E-01,
     4  0.4602799258277928E+00,0.1596361057568556E+00,
     5  0.7524999924246040E+00,0.9498336474156162E+00/
c 
        data wts5/
     1  0.1113635024916333E+00,0.2216627940962435E+00,
     2  0.2550957370683026E+00,0.2001613963896193E+00,
     3  0.8283574479929280E-01,0.1579781227372286E+00,
     4  0.3117957006034794E+00,0.2745833512528460E+00,
     5  0.2577570782145362E+00,0.1267665723468184E+00/
c 
c       build the first columns of matrices points, weights
c 
        do 1600 i=1,10
c 
        points(i,1)=pts1(i)
        points(i,2)=pts2(i)
        points(i,3)=pts3(i)
        points(i,4)=pts4(i)
        points(i,5)=pts5(i)
c 
        weights(i,1)=wts1(i)
        weights(i,2)=wts2(i)
        weights(i,3)=wts3(i)
        weights(i,4)=wts4(i)
        weights(i,5)=wts5(i)
 1600 continue
c 
        n=10
        do 1800 i=1,10
c 
        points(i,10)=-points(n-i+1,1)
        weights(i,10)=weights(n-i+1,1)
c 
        points(i,9)=-points(n-i+1,2)
        weights(i,9)=weights(n-i+1,2)
c 
        points(i,8)=-points(n-i+1,3)
        weights(i,8)=weights(n-i+1,3)
c 
        points(i,7)=-points(n-i+1,4)
        weights(i,7)=weights(n-i+1,4)
c 
        points(i,6)=-points(n-i+1,5)
        weights(i,6)=weights(n-i+1,5)
  
 1800 continue
c 
        do 2000 i=1,n
c 
        call bubble(points(1,i),weights(1,i),n)
  
 2000 continue
  
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine qu10near(points,weights,n)
        implicit real *8 (a-h,o-z)
        save
        dimension points(25),weights(25),pts(25),whts(25)
c 
c       this subroutine returns to the user the nodes and weights
c       of a 25-point quadrature evaluating integrals on the interval
c       [-1,1]  of functions f of the form
c 
c        phi(x)+psi(x)*log(|x-x_0),                                (2)
c 
c       with phi, psi polynomials of order 19, and x_0 on the
c       interval [1.002, \infty]
c 
c       PLEASE NOTE THAT THIS SET OF QUADRTAURE NODES BECOMES EFFECTIVE
C       RELATIVELY CLOSE TO THE RIGHT END OF THE INTERVAL [-1,1]!
C       FOR A SET OF NODES (34 OF THEM) THAT BECOMES EFFECTIVE AT THE
C       POINT 1.000001, USE THE SUBROUTINE QU10NEA2
c 
        data pts/
     1  -.9868079513084912E+00,-.9313713985637742E+00,
     2  -.8351474330482964E+00,-.7038230006310025E+00,
     3  -.5449889659962071E+00,-.3674797037434158E+00,
     4  -.1806333747770838E+00,0.6469995759472734E-02,
     5  0.1857161993639606E+00,0.3505340282335115E+00,
     6  0.4962759155016478E+00,0.6203584049552849E+00,
     7  0.7221545437898833E+00,0.8026745130945755E+00,
     8  0.8641104537033558E+00,0.9093403505295012E+00,
     9  0.9414764029375800E+00,0.9635131868568672E+00,
     *  0.9780963820516231E+00,0.9874061851916429E+00,
     *  0.9931336721028822E+00,0.9965208518059550E+00,
     *  0.9984340218257851E+00,0.9994449905472546E+00,
     *  0.9999030522998785E+00/
c 
        data whts/
     1  0.3375398278085781E-01,0.7658900778650240E-01,
     2  0.1149011473014309E+00,0.1464716008520721E+00,
     3  0.1697115816496150E+00,0.1837376980421925E+00,
     4  0.1884296423985280E+00,0.1844149486887468E+00,
     5  0.1729758036227731E+00,0.1558828955555725E+00,
     6  0.1351721843243647E+00,0.1128938150467626E+00,
     7  0.9087465083521377E-01,0.7053912550864444E-01,
     8  0.5282008001565255E-01,0.3816432298995675E-01,
     9  0.2661071367522212E-01,0.1790568895894973E-01,
     *  0.1162438953748007E-01,0.7276320151441085E-02,
     *  0.4384370278100738E-02,0.2532978865932456E-02,
     *  0.1387439296179184E-02,0.6918255471516529E-03,
     *  0.2537862906569214E-03/
c 
        n=25
        do 1200 i=1,n
c 
        points(i)=pts(i)
        weights(i)=whts(i)
 1200 continue
c 
        call bubble(points,weights,n)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine qu10nea2(points,weights,n)
        implicit real *8 (a-h,o-z)
        save
        dimension points(34),weights(34),pts(34),whts(34)
c 
c       this subroutine returns to the user the nodes and weights
c       of a 34-point quadrature evaluating integrals on the interval
c       [-1,1]  of functions f of the form
c 
c        phi(x)+psi(x)*log(|x-x_0),                                (2)
c 
c       with phi, psi polynomials of order 19, and x_0 on the
c       interval [1.0000001, \infty]
c 
c       PLEASE NOTE THAT A SMALLER SET OF 25 QUADRTAURE NODES IS
C       EFFECTIVE STARTING AT THE POINT FARTHER FROM THE RIGHT END
C       OF THE INTERVAL [-1,1] (SPECIFICALLY, 1.002), AND SHOULD BE
C       USED IN THAT ENVIRONMENT. SEE THE SUBROUTINE QU10NEAR
C 
        data pts/
     1  -.9849805433760962E+00,-.9219410357745987E+00,
     2  -.8128292799190389E+00,-.6646128633084092E+00,
     3  -.4865840834313151E+00,-.2895430097661786E+00,
     4  -.8487299410726155E-01,0.1164259192266738E+00,
     5  0.3046585524271848E+00,0.4722214196879411E+00,
     6  0.6141664068675482E+00,0.8158420311881339E+00,
     7  0.7284566224512793E+00,0.8793242240454028E+00,
     8  0.9524615978034928E+00,0.9711755000852873E+00,
     9  0.9828491374171867E+00,0.9899682249534673E+00,
     *  0.9232800780039728E+00,0.9942255748805810E+00,
     *  0.9967265386062437E+00,0.9981714067947011E+00,
     *  0.9989930857452949E+00,0.9994534043047958E+00,
     *  0.9997075770580029E+00,0.9998459394519500E+00,
     *  0.9999201878070153E+00,0.9999594475293126E+00,
     *  0.9999798861098794E+00,0.9999903486994319E+00,
     *  0.9999956033896319E+00,0.9999981795311328E+00,
     *  0.9999993916581695E+00,0.9999998971836818E+00/
c 
        data whts/
     1  0.3842092808267745E-01,0.8700589004861498E-01,
     2  0.1300442467176496E+00,0.1648262904407643E+00,
     3  0.1894156638807060E+00,0.2027529374248529E+00,
     4  0.2047381789830294E+00,0.1962304599285781E+00,
     5  0.1789579603273322E+00,0.1553386415928526E+00,
     6  0.1282138251321751E+00,0.7478556672577059E-01,
     7  0.1005006587821183E+00,0.5292578109350610E-01,
     8  0.2329582331146040E-01,0.1470378563551970E-01,
     9  0.9054904583994560E-02,0.5461804952976727E-02,
     *  0.3579215807321622E-01,0.3234043346027102E-02,
     *  0.1882319504974399E-02,0.1077987513486747E-02,
     *  0.6079562026039394E-03,0.3378716045832491E-03,
     *  0.1851050184055118E-03,0.9997671706614300E-04,
     *  0.5321976912692718E-04,0.2790399536319090E-04,
     *  0.1439552694561194E-04,0.7295049502162178E-05,
     *  0.3619616300469218E-05,0.1742911820829311E-05,
     *  0.7863150996170421E-06,0.2711908988157842E-06/
c 
        n=34
c 
        do 1200 i=1,n
c 
        points(i)=pts(i)
        weights(i)=whts(i)
 1200 continue
c 
        call bubble(points,weights,n)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine bubble(pts,wts,n)
        implicit real *8 (a-h,o-z)
        save
        dimension pts(n),wts(n)
c 
        do 1400 j=1,n
        do 1200 i=1,n-1
c 
        if(pts(i+1) .gt. pts(i)) goto 1200
c 
        d=pts(i+1)
        pts(i+1)=pts(i)
        pts(i)=d
c 
        d=wts(i+1)
        wts(i+1)=wts(i)
        wts(i)=d
 1200 continue
 1400 continue
c 
        return
        end
  
  
  
  
  
  
  
  
  
  
  
  
