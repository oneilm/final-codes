        implicit real *8 (a-h,o-z)
        real *8 coefs(10 000),ttest(10 000),ftest(10 000)
        complex *16 fc(20 000),wsave(100 000),com
        dimension rea(2),tfour(20 000),frea(20 000),fima(20 000),
     1      fabs(20 000)
        equivalence(rea(1),com)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
C 
        PRINT *, 'ENTER ntest'
        READ *,ntest
        CALL PRINf('ntest=*',ntest,1 )
C 
        PRINT *, 'ENTER ndigits'
        READ *,ndigits
        CALL PRINf('ndigits=*',ndigits,1 )
c 
c        retrieve the bell corresponding to this epsilon
c 
        call bellret(ier,ndigits,shift,coefs,ncoss,
     1    belllen,waves)
c 
        call prinf('after bellrt, ier=*',ier,1)
        if(ier .ne. 0) stop
  
        call prin2('after bellrt, shift=*',shift,1)
        call prin2('after bellrt, belllen=*',belllen,1)
        call prin2('after bellrt, waves=*',waves,1)
        call prinf('after bellrt, ncoss=*',ncoss,1)
        call prin2('after bellrt, coefs=*',coefs,ncoss)
c 
c       construct the test values of the argument and
c       evaluate the bell at them
c 
        t0=0
        t1=1
c 
cccc        h=(t1-t0)/(ntest-1)
        h=(t1-t0)/ntest
c 
        do 1200 i=1,ntest+1
        ttest(i)=(i-1)*h
        call belleva(ttest(i),coefs,ncoss,shift,belllen,ftest(i))
 1200 continue
c 
        call prin2('ttest as created*',ttest,ntest)
        call prin2('ftest as created*',ftest,ntest)
c 
c        flip the obtained bell, double it, and Fourier-transform
c 
        do 2200 i=1,ntest
        fc(i)=ftest(i)
  
        if(i .ne. 1) fc(2*ntest-i+2)=ftest(i)
 2200 continue
        call prin2('fc as constructed*',fc,ntest*4)
c 
        call DCFFTI(Ntest*2,WSAVE)
c 
        call DCFFTF(Ntest*2,fC,WSAVE)
c 
        call prin2('and fc after FFT*',fc,ntest*4)
c 
        do 2300 i=2,ntest*2
        fc(i)=fc(i)/fc(1)
 2300 continue
        fc(1)=1
c 
        call prin2('and fc after scaling*',fc,ntest*4)
c 
c       plot the bell
c 
        iw=31
        call lotagraph(iw,ttest,ftest,ntest,
     1      'bell as constructed*')
c 
c       plot the real part of the Fourier transform of the bell
c 
        com=fc(1)
        offset=log10(dabs(rea(1)))
  
        call prin2('offset*',offset,1)
        delta=1.0d-15
        do 2400 i=1,ntest*2
        tfour(i)=i
        com=fc(i)
        frea(i)=rea(1)
        fima(i)=rea(2)
c 
cccc        frea(i)=log10( abs(rea(1))+delta)-offset
cccc        fima(i)=log10( abs(rea(2))+delta)-offset
c 
        fabs(i)=log10(abs(fc(i))+delta)-offset
 2400 continue
c 
        iw=32
        call lotagraph(iw,tfour,frea,ntest*2,
     1      'real part of the bell s Fourier transform*')
c 
        iw=33
        call lotagraph(iw,tfour,fima,ntest*2,
     1      'imaginary part of the bell s Fourier transform*')
  
c 
        iw=34
        call lotagraph(iw,tfour,fabs,ntest*2,
     1      'absolute value of the bell s Fourier transform*')
  
        stop
        end
c 
c 
c 
c 
c 
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        this is the end of the debugging code and the
c        start of the actual bell-generating routines.
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c 
c 
c 
        subroutine belleva(x,coefs,ncoss,shift,belllen,bell)
        implicit real *8 (a-h,o-z)
        save
        dimension coefs(1)
c 
c        this subroutine evaluates at the point x \in [0,1]
c        a smooth step function, i.e. function bell(x)
c        such that
c 
c            | bell(x) | < eps,                                       (1)
c 
c         for all x \leq 0, and
c 
c            | bell(1) -1 | < eps,                                    (2)
c 
c         for all x \geq 1.
c 
c        and bell(x)= \sum_{i=1}^{ncoss}
c          coefs(i) * cos((i-1)*(x*belllen+shift))                    (3)
c 
c        The claim to fame of the function bell :[0,1] \to R^1 is that
c        it has (approximately) the lowest frequency content of all
c        functions satisfying the conditions (1),(2).
c 
c        The coefficients coefs (as well as parameters ncoss,
c        shift, bellen) must have been generated by a prior call to
c        the subroutine bellret (see).
c 
c                          Input parameters:
c 
c  shift - the parameter in (3)
c  coefs - the coefficients in (3)
c  ncoss - the number of terms in the expansion (3)
c  belllen - the scaling coefficient in (3)
c 
c                          Output parameters:
c 
c  bell - the value of the step-function at the point x
c 
c       . . . evaluate the bell at the point x
c 
        bell=0
        do 1200 i=1,ncoss
        bell=bell+coefs(i)*cos((i-1)*(x*belllen+shift))
 1200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine bellret(ier,ndigits,shift,coefs,ncoss,
     1    belllen,waves)
        implicit real *8 (a-h,o-z)
c 
c        this subroutine generates the coefficients of the cosine
c        expansion of a smooth step function, i.e. function bell(x)
c        such that
c 
c            | bell(x) | < eps,                                       (1)
c 
c         for all x \leq 0, and
c 
c            | bell(1) -1 | < eps,                                    (2)
c 
c         for all x \geq 1.
c 
c        and bell(x)= \sum_{i=1}^{ncoss}
c          coefs(i) * cos((i-1)*(x*belllen+shift))                    (3)
c 
c        The claim to fame of the function bell :[0,1] \to R^1 is that
c        it has (approximately) the lowest frequency content of all
c        functions satisfying the conditions (1),(2).
c 
c                          Input parameters:
c 
c  ndigits - specifies the accuracy eps in the formulae (1),(2).
c        specifically, the subroutine sets eps to 10**(-ndigits);
c        The permitted values of ndigits are 2,3,...,14.
c 
c                          Output parameters:
c 
c  ier - error return code.
c              ier=0 means successful conclusion
c              ier=8 means that ndigts .lt. 2 or ndigits .gt. 14
c  shift - the parameter in (3)
c  coefs - the coefficients in (3)
c  ncoss - the number of terms in the expansion (3)
c  belllen - the scaling coefficient in (3)
c  waves - the length of the transition region in terms of the wavelengths
c              of the highest frequency present in (2)
c 
        save
        dimension coefs(1)
         real *8
     1  a021001(10),a031001(12),a041001(13),a051001(15)
c 
         real *8 a061001(16),a061002(2),a062001(2)
         equivalence
     1  (a061001(16),a062001(1)),(a061002(1),a062001(2))
c 
         real *8
     1  a071001(16),a071002(4),a072001(2)
         equivalence
     1  (a071001(16),a072001(1)),(a071002(1),a072001(2))
c 
         real *8
     1  a081001(16),a081002( 3)
         real *8
     1  a082001(2)
         equivalence
     1  (a081001(16),a082001(1)),(a081002(1),a082001(2))
c 
         real *8
     1  a101001(16),a101002(10)
         real *8
     1  a102001(2)
         equivalence
     1  (a101001(16),a102001(1)),(a101002(1),a102001(2))
c 
         real *8
     1  a111001(16),a111002(10)
         real *8
     1  a112001(2)
         equivalence
     1  (a111001(16),a112001(1)),(a111002(1),a112001(2))
c 
         real *8
     1  a121001(16),a121002(14)
         real *8
     1  a122001(2)
         equivalence
     1  (a121001(16),a122001(1)),(a121002(1),a122001(2))
c 
         real *8
     1  a131001(16),a131002(16)
         real *8
     1  a132001(2)
         equivalence
     1  (a131001(16),a132001(1)),(a131002(1),a132001(2))
c 
         real *8
     1  a141001(16),a141002(16),a141003( 2)
         real *8
     1  a142001(2),a142002(2)
         equivalence
     1  (a141001(16),a142001(1)),(a141002(1),a142001(2)),
     2  (a141002(16),a142002(1)),(a141003(1),a142002(2))
c 
         real *8
     1  a091001(16),a091002( 5)
         real *8
     1  a092001(2)
         equivalence
     1  (a091001(16),a092001(1)),(a091002(1),a092001(2))
         data a091001/
     1       .7643557767651380D+00, .4229226564853893D+00,
     2      -.2983863978861431D+00, .1476878104856316D+00,
     3      -.2319375153124024D-01,-.4520552787499382D-01,
     4       .5845614045295102D-01,-.3874588346924604D-01,
     5       .1192728883518737D-01, .5796569036137246D-02,
     6      -.1100645818875433D-01, .8276175544584536D-02,
     7      -.3558340363024086D-02, .3053282521448297D-03,
     8       .8885873844173872D-03,-.8493407732223650D-03/
         data a091002/
     1       .4683173384363421D-03,-.1794749082947602D-03,
     2       .4812874820548043D-04,-.8313252373568385D-05,
     3       .7089148119090272D-06/
  
          data shift09/ .1572532561484233D+01/,
     1   bellen09/ .1488862585026712D+01/,ncoss09/21/,
     2   waves09/ .4739196799831570D+01/
c 
         data a101001/
     1       .7363054352965710D+00, .4636439530097259D+00,
     2      -.3026615934088443D+00, .1166983409979529D+00,
     3       .2228652266941917D-01,-.7964160884981655D-01,
     4       .6722956881004725D-01,-.2428080622623282D-01,
     5      -.1184007282785399D-01, .2442161428883887D-01,
     6      -.1765871037474883D-01, .4788498684812595D-02,
     7       .3908667878890873D-02,-.5776256212098999D-02,
     8       .3514949149217457D-02,-.7374469296848851D-03/
         data a101002/
     1      -.7099058622005291D-03, .8518481552128887D-03,
     2      -.4539880582188875D-03, .9921238179910405D-04,
     3       .4746887809157340D-04,-.5760720131400754D-04,
     4       .2977394304601655D-04,-.9562962113975990D-05,
     5       .1883940043528357D-05,-.1792003483638115D-06/
  
          data shift10/ .1570879004637246D+01/,
     1   bellen10/ .1400231938028232D+01/,ncoss10/26/,
     2   waves10/ .5571345860308438D+01/
c 
         data a111001/
     1       .7619342333283250D+00, .4276646446738760D+00,
     2      -.3022032758445941D+00, .1487758339246545D+00,
     3      -.2005265033410434D-01,-.5201955206023457D-01,
     4       .6599829690772643D-01,-.4338160500831945D-01,
     5       .1173090576943680D-01, .9916184175032931D-02,
     6      -.1623661483894504D-01, .1192910958427234D-01,
     7      -.4561133944605161D-02,-.7223949902631700D-03,
     8       .2595281927218540D-02,-.2173796196379104D-02/
         data a111002/
     1       .1059081222024067D-02,-.2166851421185096D-03,
     2      -.1390329593824966D-03, .1779675181922006D-03,
     3      -.1099685046446558D-03, .4700722401713838D-04,
     4      -.1465519149369965D-04, .3253192710535821D-05,
     5      -.4675421593587349D-06, .3310923435176988D-07/
  
          data shift11/ .1571540427376041D+01/,
     1   bellen11/ .1485555471332739D+01/,ncoss11/26/,
     2   waves11/ .5910837412495394D+01/
c 
         data a121001/
     1       .7484097397932418D+00, .4474820158295624D+00,
     2      -.3047044024410046D+00, .1340747817278568D+00,
     3       .2381752802604024D-02,-.6984621345076295D-01,
     4       .7139430809177576D-01,-.3657254016875451D-01,
     5      -.8627475994993295D-03, .2068855891547915D-01,
     6      -.2079677365006402D-01, .1042526701520445D-01,
     7      -.5224178798050015D-04,-.4937288472300977D-02,
     8       .4771385409647879D-02,-.2339804235225448D-02/
         data a121002/
     1       .1827934736335808D-03, .7456700230884726D-03,
     2      -.7174033472459562D-03, .3556887547909816D-03,
     3      -.6818763474620012D-04,-.4706116471511307D-04,
     4       .5314577051114207D-04,-.2794919083646614D-04,
     5       .8941517229589557D-05,-.1333193652603514D-05,
     6      -.2598325839848416D-06, .2023562001511036D-06,
     7      -.5056472958636079D-07, .5253908145664064D-08/
c 
          data shift12/ .1570217581898451D+01/,
     1   bellen12/ .1460090695889145D+01/,ncoss12/30/,
     2   waves12/ .6739038896784041D+01/
c 
         data a131001/
     1      0.7458497704884032D+00,0.4515234301673326D+00,
     2      -.3059431412774121D+00,0.1319660320259424D+00,
     3      0.6866854542687525D-02,-.7457994918062562D-01,
     4      0.7413539582569672D-01,-.3617929838275591D-01,
     5      -.3821266588878148D-02,0.2436632160297016D-01,
     6      -.2325759554327966D-01,0.1075244244008057D-01,
     7      0.1336015103253427D-02,-.6840941383771983D-02,
     8      0.6095051045400239D-02,-.2667134885015925D-02/
         data a131002/
     1      -.2321472693336557D-03,0.1375425725801576D-02,
     2      -.1162069227453696D-02,0.5011225049860961D-03,
     3      -.8969637981642418D-05,-.1657241958695529D-03,
     4      0.1402363253761760D-03,-.6345251723383736D-04,
     5      0.1113465115411384D-04,0.7356840070227565D-05,
     6      -.7827360146441762D-05,0.4035838307736876D-05,
     7      -.1391289349189919D-05,0.3297145490672966D-06,
     8      -.4978809338067784D-07,0.3685260576616166D-08/
c 
          data shift13/0.1570217581898451D+01/,
     1   bellen13/0.1472988439295640D+01/,ncoss13/32/,
     2   waves13/0.7267435128164640D+01/
c 
         data a141001/
     1      0.7528122226141722D+00,0.4420516588453662D+00,
     2      -.3067184152515317D+00,0.1417948035513841D+00,
     3      -.5110145722998855D-02,-.6790997906103520D-01,
     4      0.7581842785066114D-01,-.4359137689445185D-01,
     5      0.3775795775693809D-02,0.2092247776264682D-01,
     6      -.2470353392800002D-01,0.1476130646488070D-01,
     7      -.2181833191942940D-02,-.5524836482684357D-02,
     8      0.6810602086990080D-02,-.4180785142151859D-02/
         data a141002/
     1      0.9332869504977770D-03,0.9942340146217241D-03,
     2      -.1367184659497391D-02,0.8795031575408905D-03,
     3      -.2744206237802587D-03,-.8030775511255443D-04,
     4      0.1681419356494381D-03,-.1192557296387221D-03,
     5      0.4908225538636047D-04,-.6470847902049055D-05,
     6      -.7620555202695809D-05,0.7441622205658471D-05,
     7      -.3977351705950825D-05,0.1502876901265932D-05,
     8      -.4174790655665884D-06,0.8293335894405784D-07/
         data a141003/
     1      -.1070532797923474D-07,0.6840759766720436D-09/
c 
          data shift14/0.1570879004637245D+01/,
     1   bellen14/0.1464059232321913D+01/,ncoss14/34/,
     2   waves14/0.7689404705510816D+01/
c 
         data a021001/
     1       .6053956340755589D+00, .5938743987911688D+00,
     2      -.1852431984032939D+00,-.1023007464216133D+00,
     3       .1232895809926354D+00,-.7438078127054907D-02,
     4      -.5745285188905543D-01, .2952483752197335D-01,
     5       .1413738277994295D-01,-.1924976746198741D-01/
c 
          data shift02/ .1558388772870766D+01/,
     1   bellen02/ .6882056576584277D+00/,ncoss02/10/,
     2   waves02/ .9857819905213272D+00/
c 
         data a031001/
     1       .6554053094484179D+00, .5524916085576551D+00,
     2      -.2459153191712111D+00,-.1891767273854600D-01,
     3       .1107884971703763D+00,-.6222404666944975D-01,
     4      -.1079013412585298D-01, .3437851026800580D-01,
     5      -.1640351395589483D-01,-.3492991610874020D-02,
     6       .7305625792748465D-02,-.2706623631567897D-02/
c 
          data shift03/ .1571623497056505D+01/,
     1   bellen03/ .9462827792803379D+00/,ncoss03/12/,
     2   waves03/ .1656661400737230D+01/
c 
         data a041001/
     1       .6825158178941713D+00, .5247521077509075D+00,
     2      -.2689116049989937D+00, .2678510521867087D-01,
     3       .8742236404291184D-01,-.7530599597646512D-01,
     4       .1579990079395468D-01, .2125303576110858D-01,
     5      -.2104426890780379D-01, .5960304490656361D-02,
     6       .3198560915643166D-02,-.3466779827443726D-02,
     7       .1087300517522000D-02/
c 
          data shift04/ .1568314816010071D+01/,
     1   bellen04/ .1128260236834249D+01/,ncoss04/13/,
     2   waves04/ .2154818325434440D+01/
c 
         data a051001/
     1       .7104120476196825D+00, .4932308249889199D+00,
     2      -.2863691172460001D+00, .7216553698864252D-01,
     3       .5564270965593014D-01,-.7791142323906103D-01,
     4       .3958965941197767D-01, .1934804285314537D-02,
     5      -.1814359554866653D-01, .1325243774933803D-01,
     6      -.3375856491773806D-02,-.1788495667446429D-02,
     7       .2087831594162704D-02,-.8881524596016997D-03,
     8       .1612443859740161D-03/
c 
          data shift05/ .1571623497056505D+01/,
     1   bellen05/ .1263916159738074D+01/,ncoss05/15/,
     2   waves05/ .2816219062664560D+01/
c 
         data a061001/
     1       .7038076239269709D+00, .5027971862498118D+00,
     2      -.2869713292325525D+00, .6358554463674791D-01,
     3       .6800038877510228D-01,-.8631562021170471D-01,
     4       .3940493329060428D-01, .9036956056858798D-02,
     5      -.2613027754114167D-01, .1695208641073167D-01,
     6      -.2059898589282143D-02,-.5237160213810244D-02,
     7       .4588674809866285D-02,-.1465753196906162D-02,
     8      -.3733634333176731D-03, .5975885276909976D-03/
         data a061002/
     1      -.2670381936622739D-03, .4982882163019569D-04/
c 
          data shift06/ .1568314816010071D+01/,
     1   bellen06/ .1247372754505900D+01/,ncoss06/18/,
     2   waves06/ .3374934175882043D+01/
c 
         data a071001/
     1       .7216950017224211D+00, .4812317807555351D+00,
     2      -.2954469371576006D+00, .9236026565118963D-01,
     3       .4337221267665335D-01,-.8260796098797379D-01,
     4       .5372489915912662D-01,-.7860785581244152D-02,
     5      -.1901095836772260D-01, .2043439606724708D-01,
     6      -.9129333967203329D-02,-.9314098086262938D-03,
     7       .4338195455608017D-02,-.3077593866672819D-02,
     8       .9083300428580979D-03, .2398409719481314D-03/
         data a071002/
     1      -.3945107170602797D-03, .2042448937307646D-03,
     2      -.5694742291527392D-04, .7269875281701878D-05/
c 
          data shift07/ .1571623497056505D+01/,
     1   bellen07/ .1323472418573899D+01/,ncoss07/20/,
     2   waves07/ .4002106371774618D+01/
c 
         data a081001/
     1       .7519456846037513D+00, .4401493226946239D+00,
     2      -.2982323692426962D+00, .1322434378832087D+00,
     3      -.3746407667494753D-02,-.5705579041673717D-01,
     4       .5834937289748275D-01,-.3051108406716311D-01,
     5       .2788481709794087D-02, .1078786811133925D-01,
     6      -.1119252405956546D-01, .5933471834731667D-02,
     7      -.1201008571163611D-02,-.9011159780546036D-03,
     8       .1070372457792187D-02,-.5987745588242490D-03/
         data a081002/
     1       .2132491235665612D-03,-.4733875070317052D-04,
     2       .5155338990740036D-05/
c 
          data shift08/ .1571209716006643D+01/,
     1   bellen08/ .1457114293564570D+01/,ncoss08/19/,
     2   waves08/ .4174324964471815D+01/
c 
        ier=0
        if(( ndigits .gt. 14) .or. (ndigits .lt. 2) ) ier=8
        if(ier .eq. 8) return
c 
c        return to the user the array of cosine coefficients of the
c        bell corresponding to the precision eps
c 
        if(ndigits .ne. 2) goto 1250
c 
         shift=shift02
         ncoss=ncoss02
         belllen=bellen02
         waves=waves02
c 
        do 1200 i=1,ncoss
        coefs(i)=a021001(i)
 1200 continue
 1250 continue
c 
        if(ndigits .ne. 3) goto 1350
c 
         shift=shift03
         ncoss=ncoss03
         belllen=bellen03
         waves=waves03
c 
        do 1300 i=1,ncoss
        coefs(i)=a031001(i)
 1300 continue
 1350 continue
c 
        if(ndigits .ne. 4) goto 1450
c 
         shift=shift04
         ncoss=ncoss04
         belllen=bellen04
         waves=waves04
c 
        do 1400 i=1,ncoss
        coefs(i)=a041001(i)
 1400 continue
 1450 continue
c 
        if(ndigits .ne. 5) goto 1550
c 
         shift=shift05
         ncoss=ncoss05
         belllen=bellen05
         waves=waves05
c 
        do 1500 i=1,ncoss
        coefs(i)=a051001(i)
 1500 continue
 1550 continue
c 
        if(ndigits .ne. 6) goto 1650
c 
         shift=shift06
         ncoss=ncoss06
         belllen=bellen06
         waves=waves06
c 
        do 1600 i=1,ncoss
        coefs(i)=a061001(i)
 1600 continue
 1650 continue
c 
        if(ndigits .ne. 7) goto 1750
c 
         shift=shift07
         ncoss=ncoss07
         belllen=bellen07
         waves=waves07
c 
        do 1700 i=1,ncoss
        coefs(i)=a071001(i)
 1700 continue
 1750 continue
c 
        if(ndigits .ne. 8) goto 1850
c 
         shift=shift08
         ncoss=ncoss08
         belllen=bellen08
         waves=waves08
c 
        do 1800 i=1,ncoss
        coefs(i)=a081001(i)
 1800 continue
 1850 continue
c 
        if(ndigits .ne. 9) goto 1950
c 
         shift=shift09
         ncoss=ncoss09
         belllen=bellen09
         waves=waves09
c 
        do 1900 i=1,ncoss
        coefs(i)=a091001(i)
 1900 continue
 1950 continue
c 
        if(ndigits .ne. 10) goto 2050
c 
         shift=shift10
         ncoss=ncoss10
         belllen=bellen10
         waves=waves10
c 
        do 2000 i=1,ncoss
        coefs(i)=a101001(i)
 2000 continue
 2050 continue
c 
        if(ndigits .ne. 11) goto 2150
c 
         shift=shift11
         ncoss=ncoss11
         belllen=bellen11
         waves=waves11
c 
        do 2100 i=1,ncoss
        coefs(i)=a111001(i)
 2100 continue
 2150 continue
c 
        if(ndigits .ne. 12) goto 2250
c 
         shift=shift12
         ncoss=ncoss12
         belllen=bellen12
         waves=waves12
c 
        do 2200 i=1,ncoss
        coefs(i)=a121001(i)
 2200 continue
 2250 continue
c 
        if(ndigits .ne. 13) goto 2350
c 
         shift=shift13
         ncoss=ncoss13
         belllen=bellen13
         waves=waves13
c 
        do 2300 i=1,ncoss
        coefs(i)=a131001(i)
 2300 continue
 2350 continue
c 
        if(ndigits .ne. 14) goto 2450
c 
         shift=shift14
         ncoss=ncoss14
         belllen=bellen14
         waves=waves14
c 
        do 2400 i=1,ncoss
        coefs(i)=a141001(i)
 2400 continue
 2450 continue
        return
        end
  
  
  
  
  
  
  
