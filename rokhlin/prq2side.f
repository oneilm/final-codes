        implicit real *8 (a-h,o-z)
c 
        external fun1
c 
        dimension xs(10 000),whts(10 000),fs(10 000)
c 
        call prini(6,13)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER n'
        READ *,n
        CALL PRINf('n=*',n,1 )
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER ndigits'
        READ *,ndigits
        CALL PRINf('ndigits=*',ndigits,1 )
c 
c       construct the sl0 nodes and weights
c 
  
  
        call prq1side(ier,ndigits,n,xs,whts)
cccc        call prq2side(ier,ndigits,n,xs,whts)
  
        call prin2('after prq1side, xs=*',xs,n)
        call prin2('after prq1side, whts=*',whts,n)
  
cccc        stop
c 
c       construct the values of the test function at the
c       sl0 nodes, and evaluate the integral
c 
        rint=0
        do 1200 i=1,n
c 
        fs(i)=fun1(xs(i),par1,par2)
c 
        rint=rint+whts(i)*fs(i)
c 
 1200 continue
c 
        call prin2('fs as calculated*',fs,n)
        call prin2('rint as calculated*',rint,1)
  
c 
c        evaluate the same integral via adapgaus
  
c 
        a=-1
        b=0
        m=20
        eps=1.0d-15
  
  
c 
        call adapgaus(ier,a,b,fun1,x,n,m,eps,
     1      rint2,maxrec,numint)
  
        call prinf('after adapgaus, ier=*',ier,1)
        call prin2('after adapgaus, eps=*',eps,1)
  
        call prin2('after adapgaus, rint2=*',rint2,1)
        call prinf('after adapgaus, numint=*',numint,1)
  
        call prin2('and rint2-rint=*',rint2-rint,1)
  
  
  
  
  
  
        stop
        end
  
  
c 
c 
c 
c 
c 
        function fun1(x,par1,par2)
        implicit real *8 (a-h,o-z)
c 
        save
        bb=1+1.0d-12
cccc        bb=1+1.0d-20
cccc        bb=1+1
  
  
  
cccc        fun1=(1-x**2)*log(1-x**2)
        fun1=sqrt(bb-x**2)*log(bb-x**2)
        fun1=log(bb-x**2)
        fun1=1/sqrt(sqrt((bb-x**2)))
c 
  
  
  
cccc        fun1=(bb-x**2)*log(bb-x**2) *cos(120*x)
        fun1=(bb-x**2)*log(bb-x**2)
  
        fun1=fun1+2
  
        nnn=40
        nnn=0
        call legepol(x,nnn,pol,der)
  
        fun1=fun1*pol
  
cccc        fun1=pol
  
        fun1=3
  
  
        fun1=log(bb-x**2)
cccc        fun1=1/sqrt(bb-x**2)
cccc        fun1=1/sqrt(bb+x)
  
cccc        fun1=(x+1)*log(x+1)
  
        fun1=log(x+1)
  
        return
        end
c 
c 
c 
c 
c 
        subroutine prq2side(ier,ndigits,n,xs,ws)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),ws(1),xs16by26(13),ws16by26(13),
     1      xs16by40(20),ws16by40(20),xs16by60(30),ws16by60(30)
        dimension xs08by12(6),ws08by12(6),
     1      xs08by20(10),ws08by20(10),xs08by30(15),ws08by30(15)
c 
c        This subroutine returns to the user the nodes and weights of
c        a two-sided "universal" quadrature on the interval [-1,1],
c        for one of six combinations of the number of nodes and accuracy.
c        Following are the combinations:
c 
c           number of digits       number of nodes
c 
c                 8                      12
c                 8                      20
c                 8                      30
c                 16                     26
c                 16                     40
c                 16                     60
c 
c        Any other combination will cause the error return code to be set
c        to 16, and the execution of the subroutine to be terminated
c 
c 
c        . . . Nodes and weights for ndigits=  8, n= 12
c 
         data xs08by12/
     1   -.999999729099894785360E+00,-.999882537881617381592E+00,
     2   -.995839834709499859166E+00,-.953802580359643591765E+00,
     3   -.761118739369036024704E+00,-.303294195149091290567E+00/
c 
         data ws08by12/
     1   0.240739668555841537543E-05,0.518020851953743067901E-03,
     2   0.120235397593903208400E-01,0.916426359802621402394E-01,
     3   0.318306868570953636374E+00,0.577506528025773038583E+00/
c 
c        Nodes and weights for ndigits=  8, n= 20
c 
         data xs08by20/
     1   -.999999964872134119756E+00,-.999994844904069693730E+00,
     2   -.999882537881617381592E+00,-.998809723505868785420E+00,
     3   -.992742383035922991671E+00,-.969404890519814048263E+00,
     4   -.903573853138300038524E+00,-.761118739369036024704E+00,
     5   -.517935428865798223578E+00,-.184805969116997921534E+00/
c 
         data ws08by20/
     1   0.251587893553888154008E-06,0.192545748864176907772E-04,
     2   0.310812511172245840741E-03,0.241712033604581953853E-02,
     3   0.116619191758766702119E-01,0.392626982307954266603E-01,
     4   0.984240297536079824420E-01,0.190984121142572181824E+00,
     5   0.293759567762540749298E+00,0.363160224363420309086E+00/
c 
c        Nodes and weights for ndigits=  8, n= 30
c 
         data xs08by30/
     1   -.999999990737835239040E+00,-.999999371713544705206E+00,
     2   -.999990590707727395938E+00,-.999925257817662816394E+00,
     3   -.999596427226425295256E+00,-.998341565948008431353E+00,
     4   -.994475995636242065392E+00,-.984485197591920577018E+00,
     5   -.962258540745565921136E+00,-.918931228170263342446E+00,
     6   -.843977614969984244175E+00,-.727837611918054062723E+00,
     7   -.565566369710069592439E+00,-.360149162683674750126E+00,
     8   -.123803865822690590760E+00/
c 
         data ws08by30/
     1   0.557803472916037625096E-07,0.200537005099555983871E-05,
     2   0.219023189877352023219E-04,0.138566634343566069036E-03,
     3   0.620070585904679240280E-03,0.215632125870179281401E-02,
     4   0.614637042177997180867E-02,0.148493063720519318423E-01,
     5   0.310988230551658960603E-01,0.573550212532691449529E-01,
     6   0.942098389139264551187E-01,0.138953963482541395347E+00,
     7   0.185115589954657878136E+00,0.223652112665888611536E+00,
     8   0.245680052072321023614E+00/
c 
c        Nodes and weights for ndigits= 16, n= 26
c 
         data xs16by26/
     1   -.999999999999989539936E+00,-.999999999964481145620E+00,
     2   -.999999993486306540062E+00,-.999999662822186650009E+00,
     3   -.999992039204416856712E+00,-.999893112888412695021E+00,
     4   -.999076549322726720683E+00,-.994451794944780465762E+00,
     5   -.975558228798968795274E+00,-.917943890999820829340E+00,
     6   -.783646314990568770145E+00,-.540409566144105107587E+00,
     7   -.194424901960809770104E+00/
c 
         data ws16by26/
     1   0.120440511024722947568E-12,0.219025735350353874427E-09,
     2   0.290536878392669352995E-07,0.118171285448640879129E-05,
     3   0.227148190492242018357E-04,0.252361906647554310628E-03,
     4   0.181393652351249849832E-02,0.904681836995841248192E-02,
     5   0.327859798178494176642E-01,0.890952486471638038240E-01,
     6   0.185510473511829905650E+00,0.300268937354182075819E+00,
     7   0.381202318064118602749E+00/
c 
c        Nodes and weights for ndigits= 16, n= 40
c 
         data xs16by40/
     1   -.999999999999998901423E+00,-.999999999999079161306E+00,
     2   -.999999999924640003850E+00,-.999999997687125287226E+00,
     3   -.999999960809510804131E+00,-.999999562390947190985E+00,
     4   -.999996445409194116976E+00,-.999977654305572177075E+00,
     5   -.999886629393802540057E+00,-.999521682889993170071E+00,
     6   -.998283951201297561751E+00,-.994673932704008498804E+00,
     7   -.985503365281910622828E+00,-.965010732231033324977E+00,
     8   -.924412638868528279965E+00,-.852667981892730692388E+00,
     9   -.739026815345513683388E+00,-.577090481232282542352E+00,
     *   -.368929683053911208564E+00,-.127097463955535675797E+00/
c 
         data ws16by40/
     1   0.104087843940600518980E-13,0.474794646362121774626E-11,
     2   0.287758174510788713952E-09,0.714034229361019358610E-08,
     3   0.101914738613885592920E-06,0.980757408229216018358E-06,
     4   0.696073542250644715838E-05,0.385554760783768391682E-04,
     5   0.173181620510733715269E-03,0.648278084303752303068E-03,
     6   0.206358312337005357885E-02,0.567179597641819644551E-02,
     7   0.136200762911666841514E-01,0.288401683968840231011E-01,
     8   0.542391562248948211732E-01,0.911131295034249543507E-01,
     9   0.137309986441574331112E+00,0.186255582848380180136E+00,
     *   0.227944881578662733173E+00,0.252073573593902988865E+00/
c 
c        Nodes and weights for ndigits= 16, n= 60
c 
         data xs16by60/
     1   -.999999999999999809996E+00,-.999999999999948033819E+00,
     2   -.999999999997880803965E+00,-.999999999960649989691E+00,
     3   -.999999999542896839841E+00,-.999999996161247666847E+00,
     4   -.999999974740322689606E+00,-.999999863042333403129E+00,
     5   -.999999366961603982906E+00,-.999997444389978396352E+00,
     6   -.999990824816626311351E+00,-.999970296518791103107E+00,
     7   -.999912337477079221838E+00,-.999762068634743313586E+00,
     8   -.999401796210486579376E+00,-.998598442719560571390E+00,
     9   -.996924362860624118108E+00,-.993651063468102921668E+00,
     *   -.987625317012845137620E+00,-.977151520715363595547E+00,
     1   -.959919669944219721134E+00,-.933028916022993298555E+00,
     2   -.893154833555614310337E+00,-.836888206549400774624E+00,
     3   -.761233548679631130204E+00,-.664203841411112896133E+00,
     4   -.545399049279050177383E+00,-.406428596279089634054E+00,
     5   -.251047640319134245168E+00,-.849284655246994799414E-01/
c 
         data ws16by60/
     1   0.149837722802727978318E-14,0.224332283739078914796E-12,
     2   0.685743900199185124759E-11,0.104450918003149856417E-09,
     3   0.103906169965088241734E-08,0.766128597327812554009E-08,
     4   0.449756318073493549704E-07,0.219948329853291930877E-06,
     5   0.924076395403929520865E-06,0.340984021017146020314E-05,
     6   0.112345985556718291797E-04,0.334721221614215356087E-04,
     7   0.910849628216576560077E-04,0.228206858504618142287E-03,
     8   0.529867731948734860251E-03,0.114632148727471144612E-02,
     9   0.232113287746209519639E-02,0.441561425952663555711E-02,
     *   0.791715465965982735898E-02,0.134157090347964865888E-01,
     1   0.215341805971397509482E-01,0.328070671999727130023E-01,
     2   0.475177610922899805285E-01,0.655249294376756109806E-01,
     3   0.861255463856479021955E-01,0.108007909676503470366E+00,
     4   0.129335909473458633292E+00,0.147975319719175865770E+00,
     5   0.161831158207605126340E+00,0.169225811965369987400E+00/
c 
        ier=16
  
        call prosicp2(ier,12,n,8,ndigits,
     1    xs08by12,ws08by12,xs,ws)
c 
        call prosicp2(ier,20,n,8,ndigits,
     1    xs08by20,ws08by20,xs,ws)
c 
        call prosicp2(ier,30,n,8,ndigits,
     1    xs08by30,ws08by30,xs,ws)
c 
c 
        call prosicp2(ier,26,n,16,ndigits,
     1    xs16by26,ws16by26,xs,ws)
c 
        call prosicp2(ier,40,n,16,ndigits,
     1    xs16by40,ws16by40,xs,ws)
c 
        call prosicp2(ier,60,n,16,ndigits,
     1    xs16by60,ws16by60,xs,ws)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine prq1side(ier,ndigits,n,xs,ws)
        implicit real *8 (a-h,o-z)
        save
        dimension xs(1),ws(1),xs16by20(20),ws16by20(20),
     1      xs16by40(40),ws16by40(40),xs16by30(30),ws16by30(30)
        dimension xs08by10(10),ws08by10(10),
     1      xs08by20(20),ws08by20(20),xs08by30(30),ws08by30(30)
c 
c        This subroutine returns to the user the nodes and weights of
c        a one-sided "universal" quadrature on the interval [-1,0],
c        with -1 the "universal" ens. There are six permitted
c        combinations of the number of nodes and accuracy, as follows
c 
c           number of digits       number of nodes
c 
c                 8                      10
c                 8                      20
c                 8                      30
c                 16                     20
c                 16                     30
c                 16                     40
c 
c        Any other combination will cause the error return code to be set
c        to 16, and the execution of the subroutine to be terminated
c 
c        . . . Nodes and weights for ndigits=  8, n= 10
c 
         data xs08by10/
     1   -.999999813253686052721E+00,-.999936027322580423400E+00,
     2   -.998099696363856976756E+00,-.981899818930478575180E+00,
     3   -.916564335254533525877E+00,-.768054460201517174278E+00,
     4   -.549342198123048273888E+00,-.318591178447270259526E+00,
     5   -.133724867538290458686E+00,-.256166076288578759556E-01/
c 
         data ws08by10/
     1   0.159861691049885575343E-05,0.269943559504672507390E-03,
     2   0.517837529014611016204E-02,0.333963956339820469200E-01,
     3   0.103994654289530340740E+00,0.190796531645971936310E+00,
     4   0.235963236024975059393E+00,0.215426177418941163764E+00,
     5   0.149329785852711042756E+00,0.656433014153022155621E-01/
c 
c        Nodes and weights for ndigits=  8, n= 20
c 
         data xs08by20/
     1   -.999999985178525707345E+00,-.999998688630542192801E+00,
     2   -.999977628956435711015E+00,-.999812437778915853091E+00,
     3   -.998989949745294563844E+00,-.996043672369446504118E+00,
     4   -.987888599298756567813E+00,-.969613522433065286273E+00,
     5   -.935323916239201115728E+00,-.880036416724001907430E+00,
     6   -.801794013937641060009E+00,-.702863760766393918393E+00,
     7   -.589403744636303880666E+00,-.469881191704544743029E+00,
     8   -.353092560034175018693E+00,-.246596289406667127731E+00,
     9   -.155941193127595732810E+00,-.846352570102547745221E-01,
     *   -.345679065101535395476E-01,-.657084283752559494041E-02/
c 
         data ws08by20/
     1   0.951617252416327163463E-07,0.442217165148626740574E-05,
     2   0.540327394955019034583E-04,0.352115606061683203289E-03,
     3   0.152471440118519379927E-02,0.487682169786703975742E-02,
     4   0.122657715165627243288E-01,0.253006076364886357707E-01,
     5   0.441475303556881750615E-01,0.667752235788721423896E-01,
     6   0.893253481994801231693E-01,0.107480134178066431866E+00,
     7   0.117993601783447180578E+00,0.119558407062832116736E+00,
     8   0.112746487552020164416E+00,0.993215045080989642825E-01,
     9   0.814118407725382677797E-01,0.608950893237661455033E-01,
     *   0.391086697884789506201E-01,0.168575819407681147076E-01/
c 
c        Nodes and weights for ndigits=  8, n= 30
c 
         data xs08by30/
     1   -.999999995465714322841E+00,-.999999799404069663014E+00,
     2   -.999997691200755087420E+00,-.999984817122633934606E+00,
     3   -.999928896509343355481E+00,-.999737899463030867122E+00,
     4   -.999195707406048233936E+00,-.997871157375605554298E+00,
     5   -.995019293574516435596E+00,-.989512419427454494965E+00,
     6   -.979846976877212489805E+00,-.964260475419441155227E+00,
     7   -.940960279327969583651E+00,-.908424803946445188503E+00,
     8   -.865706071986740882458E+00,-.812655630301328223785E+00,
     9   -.750016771193762611301E+00,-.679365268890433219114E+00,
     *   -.602922045473970307953E+00,-.523289672933326941603E+00,
     1   -.443173461865516009096E+00,-.365138828702329405016E+00,
     2   -.291437015152121536735E+00,-.223909564736322221222E+00,
     3   -.163964553967948687332E+00,-.112607253172683783678E+00,
     4   -.705044479397954686606E-01,-.380631333685766556405E-01,
     5   -.155081673614736694971E-01,-.294514553981481018844E-02/
c 
         data ws08by30/
     1   0.245589913889518821601E-07,0.576339670502366130913E-06,
     2   0.486799535386316330919E-05,0.256773972243492231646E-04,
     3   0.100430228168185289330E-03,0.316023863201754114245E-03,
     4   0.838481219915332025140E-03,0.193351062028678006662E-02,
     5   0.395819039793741467855E-02,0.730809134241219819173E-02,
     6   0.123204533418547805400E-01,0.191565220436689250480E-01,
     7   0.277027080061892675756E-01,0.375301478978371540515E-01,
     8   0.479336381884889228578E-01,0.580423139771682986314E-01,
     9   0.669695208936204152916E-01,0.739587862308428661769E-01,
     *   0.784891227515316539515E-01,0.803206163622822573622E-01,
     1   0.794813088872250640313E-01,0.762111781435689194108E-01,
     2   0.708851729932990067080E-01,0.639356108660177847019E-01,
     3   0.557879845166148402876E-01,0.468168343057179276637E-01,
     4   0.373222649060204306710E-01,0.275239906628699990484E-01,
     5   0.175687850669437217945E-01,0.755716621838480397313E-02/
c 
c        Nodes and weights for ndigits= 16, n= 20
c 
         data xs16by20/
     1   -.999999999999989042809E+00,-.999999999962436500168E+00,
     2   -.999999993337472276826E+00,-.999999679015401981769E+00,
     3   -.999993183058839692777E+00,-.999919908145364500307E+00,
     4   -.999405877215916315461E+00,-.996961269861555610803E+00,
     5   -.988569320307808907800E+00,-.966729658527527102279E+00,
     6   -.921864155521278294040E+00,-.846551401502783250423E+00,
     7   -.740090232829234028561E+00,-.609990119027047067383E+00,
     8   -.469452265745686840230E+00,-.332902986440549074133E+00,
     9   -.212320674716339636460E+00,-.115691433685623394370E+00,
     *   -.473170586401741268554E-01,-.899663996888450576492E-02/
c 
         data ws16by20/
     1   0.126559429536872747920E-12,0.231142783636367864838E-09,
     2   0.293633767777849588644E-07,0.109615129184568890365E-05,
     3   0.186264278174928835072E-04,0.177520224008959110747E-03,
     4   0.107282654346218332325E-02,0.446418581040610910818E-02,
     5   0.135911641032624387279E-01,0.317620986486005580962E-01,
     6   0.592917309076419032905E-01,0.914778812534188849421E-01,
     7   0.120144754639921162717E+00,0.137773349084684830770E+00,
     8   0.140833590428404167957E+00,0.130243957461256181716E+00,
     9   0.109597441127431570220E+00,0.829660047578474531778E-01,
     *   0.535032036197896105466E-01,0.230805392161085236360E-01/
c 
c        Nodes and weights for ndigits= 16, n= 30
c 
         data xs16by30/
     1   -.999999999999998684762E+00,-.999999999998773463379E+00,
     2   -.999999999894530225181E+00,-.999999996715751022109E+00,
     3   -.999999945112416743364E+00,-.999999410160672868331E+00,
     4   -.999995487207796362452E+00,-.999973771078563061238E+00,
     5   -.999878865719757809933E+00,-.999540451671060290689E+00,
     6   -.998530465085129086129E+00,-.995955330833475821047E+00,
     7   -.990249429747843568940E+00,-.979100023137368264462E+00,
     8   -.959639131868156385208E+00,-.928946210707455787404E+00,
     9   -.884750435073266184793E+00,-.826097412601983706184E+00,
     *   -.753736007516020338766E+00,-.670094029977421494360E+00,
     1   -.578876648863245531500E+00,-.484448280611535121145E+00,
     2   -.391196738191278166195E+00,-.303034828883388701600E+00,
     3   -.223112109227824592313E+00,-.153732256127736284737E+00,
     4   -.964239815400234459146E-01,-.520983454001167738804E-01,
     5   -.212324500285282040718E-01,-.403245430203905538153E-02/
c 
         data ws16by30/
     1   0.126734956240622335462E-13,0.641159388075627990598E-11,
     2   0.405936561916933120182E-09,0.101378822029028627915E-07,
     3   0.141275385339877406803E-06,0.129269244960258003161E-05,
     4   0.852284542303239885723E-05,0.429804051786219138198E-04,
     5   0.172843766354573556838E-03,0.571711526736654301304E-03,
     6   0.159341367520349177672E-02,0.381617279617319989713E-02,
     7   0.798414880590996760489E-02,0.148013167345523108794E-01,
     8   0.246195403543120910178E-01,0.371565213351993551908E-01,
     9   0.514008769701991180613E-01,0.657798527205917474465E-01,
     *   0.785312064698290295427E-01,0.881244102227865098761E-01,
     1   0.935719669713184832590E-01,0.945449913465447066345E-01,
     2   0.913011791898088219540E-01,0.844964739787426626738E-01,
     3   0.749663699548167884144E-01,0.635413188360784566208E-01,
     4   0.509265529661300802142E-01,0.376480323373051428303E-01,
     5   0.240510203987872869555E-01,0.103471308739398951547E-01/
c 
c        Nodes and weights for ndigits= 16, n= 40
c 
         data xs16by40/
     1   -.999999999999999638213E+00,-.999999999999850990523E+00,
     2   -.999999999992186550978E+00,-.999999999827124757114E+00,
     3   -.999999997727871019351E+00,-.999999979263791480363E+00,
     4   -.999999856621124858991E+00,-.999999206545693172313E+00,
     5   -.999996351787170751847E+00,-.999985681269457789468E+00,
     6   -.999951034877759847348E+00,-.999851743339178315595E+00,
     7   -.999597329306835988124E+00,-.999008237700004006109E+00,
     8   -.997764371468272913548E+00,-.995350461097497918423E+00,
     9   -.991015085374625768967E+00,-.983764483687498345634E+00,
     *   -.972409335006481187259E+00,-.955671238715704403634E+00,
     1   -.932338571443200562382E+00,-.901444800414161659578E+00,
     2   -.862432764581415088179E+00,-.815269764207316844432E+00,
     3   -.760489988127791768482E+00,-.699158430372988125010E+00,
     4   -.632767944051531988853E+00,-.563093354764438499422E+00,
     5   -.492031154779559888297E+00,-.421450653742502302812E+00,
     6   -.353074986777930032128E+00,-.288401137975065100288E+00,
     7   -.228659712288983866694E+00,-.174809114334128649934E+00,
     8   -.127555508163444178377E+00,-.873890822901681860001E-01,
     9   -.546279769530087384795E-01,-.294629258638504824256E-01,
     *   -.119974454227582840418E-01,-.227781718570115190455E-02/
c 
         data ws16by40/
     1   0.306323159953996293057E-14,0.688802148757552364188E-12,
     2   0.269033149926221370397E-10,0.484265815622118312601E-09,
     3   0.539651338414544324975E-08,0.427424888191308276147E-07,
     4   0.260216687163476989726E-06,0.127987751960248555570E-05,
     5   0.526289764203711167827E-05,0.185491092719621790396E-04,
     6   0.571078103728377118989E-04,0.155897400539893271980E-03,
     7   0.381977637458643632856E-03,0.848603156477204212050E-03,
     8   0.172427114008955818331E-02,0.322863637784786086469E-02,
     9   0.560850295466324093409E-02,0.909272000271185977379E-02,
     *   0.138333088403933531919E-01,0.198478844161633891520E-01,
     1   0.269813650737590886660E-01,0.349010136783211678237E-01,
     2   0.431292288541331445114E-01,0.511071090134078287838E-01,
     3   0.582731681433076104215E-01,0.641386361403872542852E-01,
     4   0.683437969038799700087E-01,0.706867743339613677743E-01,
     5   0.711240694337088694613E-01,0.697484097820044878554E-01,
     6   0.667527702143166846811E-01,0.623897510768036949789E-01,
     7   0.569337220873735894009E-01,0.506504051417116620242E-01,
     8   0.437758578245055005959E-01,0.365047265616776777117E-01,
     9   0.289863743573970395161E-01,0.213270101039698862548E-01,
     *   0.135963366403290116612E-01,0.584516414634265575243E-02/
c 
        ier=16
c 
        call prosicp1(ier,10,n,8,ndigits,
     1    xs08by10,ws08by10,xs,ws)
c 
        call prosicp1(ier,20,n,8,ndigits,
     1    xs08by20,ws08by20,xs,ws)
c 
        call prosicp1(ier,30,n,8,ndigits,
     1    xs08by30,ws08by30,xs,ws)
c 
        call prosicp1(ier,20,n,16,ndigits,
     1    xs16by20,ws16by20,xs,ws)
c 
        call prosicp1(ier,30,n,16,ndigits,
     1    xs16by30,ws16by30,xs,ws)
c 
        call prosicp1(ier,40,n,16,ndigits,
     1    xs16by40,ws16by40,xs,ws)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine prosicp1(ier,n0,n,ndigits0,ndigits,
     1    xs0,ws0,xs,ws)
c 
        implicit real *8 (a-h,o-z)
        save
        dimension xs0(1),xs(1),ws0(1),ws(1)
c 
        if(n0 .ne. n) return
        if(ndigits0 .ne. ndigits) return
c 
        ier=0
        do 1200 i=1,n
c 
        ws(i)=ws0(i)
        xs(i)=xs0(i)
 1200 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine prosicp2(ier,n0,n,ndigits0,ndigits,
     1    xs0,ws0,xs,ws)
c 
        implicit real *8 (a-h,o-z)
        save
        dimension xs0(1),xs(1),ws0(1),ws(1)
c 
        if(n0 .ne. n) return
        if(ndigits0 .ne. ndigits) return
c 
        ier=0
        do 1200 i=1,n/2
c 
        ws(i)=ws0(i)
        xs(i)=xs0(i)
c 
        ws(n-i+1)=ws0(i)
        xs(n-i+1)=-xs0(i)
c 
 1200 continue
c 
        return
        end
  
  
  
