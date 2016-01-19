*     -*-fortran-*-
*
*     Tabutalet values to check APFEL
*
      double precision Q20,Q2(4)
      double precision xlha(9)
      double precision Ref(1049)
      data Q20  / 2.000E+00/
      data Q2   / 1.000E+01, 1.000E+02, 1.000E+03, 1.000E+04/
      data xlha / 1.000E-05, 1.000E-04, 1.000E-03, 1.000E-02, 1.000E-01,
     1            3.000E-01, 5.000E-01, 7.000E-01, 9.000E-01/
      data Ref /
     1    2.4423490E-01,
     1    1.1939223E-03, 7.4858485E-04, 6.8684503E+00, 1.8321320E+00,
     1    2.6560823E+01,
     1    6.0824859E-03, 3.6740014E-03, 4.3111983E+00, 9.8276287E-01,
     1    1.7089551E+01,
     1    3.1542453E-02, 1.8628394E-02, 2.6356330E+00, 4.8396049E-01,
     1    9.8707682E+00,
     1    1.6006304E-01, 9.2860680E-02, 1.5065607E+00, 2.0301457E-01,
     1    4.7524677E+00,
     1    5.9241427E-01, 3.0960038E-01, 4.6979232E-01, 3.8074420E-02,
     1    1.2535894E+00,
     1    5.5326046E-01, 2.2278383E-01, 6.6416151E-02, 4.0412765E-03,
     1    2.1485760E-01,
     1    2.6098665E-01, 7.4586190E-02, 6.2465416E-03, 3.5570529E-04,
     1    2.9829282E-02,
     1    6.2066608E-02, 1.0588029E-02, 2.0331899E-04, 1.2618583E-05,
     1    1.8974303E-03,
     1    2.0655154E-03, 1.1692735E-04, 1.7884365E-07, 1.8676692E-08,
     1    9.0467836E-06,
     1    2.0592026E+00, 6.9947860E-01, 2.0938702E-02, 2.7796199E+00,
     1    3.8760232E-01, 1.0065547E-01, 6.7038941E-04, 4.8892818E-01,
     1    1.2569250E-06, 3.4806578E-08, 0.0000000E+00, 1.2917316E-06,
     1    1.2613240E+00, 3.4949260E-01, 8.8538245E-03, 1.6196704E+00,
     1    2.6645817E-01, 5.5516869E-02, 3.0464059E-04, 3.2227968E-01,
     1    5.5984641E-06, 6.5956467E-08, 0.0000000E+00, 5.6644205E-06,
     1    7.6943795E-01, 1.5130053E-01, 2.7249652E-03, 9.2346345E-01,
     1    1.6779289E-01, 2.7004284E-02, 1.0439254E-04, 1.9490157E-01,
     1    2.4917866E-05, 5.0592241E-08, 0.0000000E+00, 2.4968458E-05,
     1    5.0422391E-01, 5.0221697E-02, 3.2382347E-04, 5.5476943E-01,
     1    9.4535761E-02, 1.0084555E-02, 1.4146659E-05, 1.0463446E-01,
     1    1.0741616E-04,-8.5261622E-08, 0.0000000E+00, 1.0733090E-04,
     1    4.0012177E-01, 4.7207136E-03, 7.5314148E-11, 4.0484248E-01,
     1    3.9493949E-02, 8.9703145E-04, 7.4202211E-13, 4.0390980E-02,
     1    3.5252575E-04,-8.8034173E-08, 0.0000000E+00, 3.5243771E-04,
     1    2.8104960E-01, 3.7331706E-04, 0.0000000E+00, 2.8142292E-01,
     1    1.3919444E-02, 2.5879799E-05, 0.0000000E+00, 1.3945324E-02,
     1    3.4567099E-04,-1.5289564E-08, 0.0000000E+00, 3.4565570E-04,
     1    1.4229595E-01, 5.4688994E-05, 0.0000000E+00, 1.4235064E-01,
     1    3.8774229E-03, 1.4371700E-06, 0.0000000E+00, 3.8788600E-03,
     1    1.8504903E-04,-2.0793803E-09, 0.0000000E+00, 1.8504695E-04,
     1    4.2159990E-02, 3.4219816E-06, 0.0000000E+00, 4.2163412E-02,
     1    5.3043971E-04, 3.3349712E-08, 0.0000000E+00, 5.3047306E-04,
     1    5.4847207E-05,-1.4171546E-10, 0.0000000E+00, 5.4847065E-05,
     1    2.3832936E-03, 1.0751953E-08, 0.0000000E+00, 2.3833043E-03,
     1    6.6110260E-06, 2.1917172E-11, 0.0000000E+00, 6.6110479E-06,
     1    3.0380922E-06,-6.7433231E-13, 0.0000000E+00, 3.0380916E-06,
     1    6.3800495E+00, 4.4951112E+00, 1.4422570E-03, 1.0876603E+01,
     1    1.1708981E+00, 1.3519116E+00, 1.2397063E-03, 2.5240494E+00,
     1    3.5987694E-02,-7.2542390E-01, 3.0659886E-04,-6.8912961E-01,
     1    3.9478452E+00, 2.5020682E+00, 7.6199065E-04, 6.4506754E+00,
     1    8.0692221E-01, 8.3625116E-01, 6.6766155E-04, 1.6438410E+00,
     1    4.0747619E-02,-4.4856639E-01, 1.8935775E-04,-4.0762942E-01,
     1    2.4530538E+00, 1.2843720E+00, 3.0311301E-04, 3.7377290E+00,
     1    5.1116811E-01, 4.6848232E-01, 2.7761429E-04, 9.7992805E-01,
     1    8.7169553E-02,-2.7621777E-01, 1.0335295E-04,-1.8894486E-01,
     1    1.6886636E+00, 5.6295094E-01, 4.0626558E-05, 2.2516552E+00,
     1    2.9530734E-01, 2.1456022E-01, 4.7258074E-05, 5.0991482E-01,
     1    3.0479657E-01,-1.6345448E-01, 4.4823934E-05, 1.4138692E-01,
     1    1.5161514E+00, 1.1164973E-01, 1.1818769E-05, 1.6278129E+00,
     1    1.3985638E-01, 4.0309014E-02, 7.1568211E-06, 1.8017255E-01,
     1    9.7717957E-01,-5.2170689E-02, 8.1484305E-06, 9.2501703E-01,
     1    1.1375627E+00, 1.0088756E-02, 1.3688489E-07, 1.1476516E+00,
     1    5.5054523E-02, 3.4398979E-03, 8.6701629E-08, 5.8494508E-02,
     1    1.0037954E+00,-5.6530750E-03, 4.5880182E-08, 9.9814242E-01,
     1    5.9396456E-01, 4.6308796E-04, 0.0000000E+00, 5.9442765E-01,
     1    1.6044182E-02, 1.8524195E-04, 0.0000000E+00, 1.6229424E-02,
     1    5.6708421E-01,-2.1030671E-04, 0.0000000E+00, 5.6687390E-01,
     1    1.8069608E-01, 3.1362838E-06, 0.0000000E+00, 1.8069921E-01,
     1    2.2593923E-03, 2.4962933E-06, 0.0000000E+00, 2.2618886E-03,
     1    1.7786561E-01, 3.1852933E-06, 0.0000000E+00, 1.7786879E-01,
     1    1.0489216E-02, 7.0095229E-09, 0.0000000E+00, 1.0489223E-02,
     1    2.8885235E-05, 1.4811944E-09, 0.0000000E+00, 2.8886716E-05,
     1    1.0459270E-02, 5.7564148E-09, 0.0000000E+00, 1.0459275E-02,
     1    1.7629258E-01,
     1    1.9542327E-03, 1.2169593E-03, 1.4804139E+01, 5.8506311E+00,
     1    7.8563463E+01,
     1    9.1796241E-03, 5.4938721E-03, 7.9694728E+00, 2.8521250E+00,
     1    4.0303920E+01,
     1    4.3139382E-02, 2.5161340E-02, 4.0483431E+00, 1.2236123E+00,
     1    1.8076741E+01,
     1    1.9233533E-01, 1.0966349E-01, 1.8562089E+00, 4.1242790E-01,
     1    6.4104169E+00,
     1    5.8096893E-01, 2.9533327E-01, 4.3848677E-01, 5.4976145E-02,
     1    1.1044630E+00,
     1    4.5616133E-01, 1.7700381E-01, 5.0625554E-02, 4.4480562E-03,
     1    1.4037824E-01,
     1    1.8726949E-01, 5.1271850E-02, 4.0973716E-03, 3.1932692E-04,
     1    1.6312206E-02,
     1    3.7929220E-02, 6.1704483E-03, 1.1337548E-04, 9.2082694E-06,
     1    9.1911057E-04,
     1    9.4813132E-04, 5.1002889E-05, 7.7645242E-08, 1.0328532E-08,
     1    3.9209800E-06,
     1    4.5864870E+00, 2.3926065E+00, 2.2579028E-01, 7.2048838E+00,
     1    8.6455735E-01, 4.7802491E-01, 3.4593101E-02, 1.3771754E+00,
     1    1.7461134E-05, 8.7350539E-07, 5.6260205E-08, 1.8390900E-05,
     1    2.4374374E+00, 1.1433220E+00, 9.8566289E-02, 3.6793257E+00,
     1    4.5658419E-01, 2.3855658E-01, 1.5946480E-02, 7.1108726E-01,
     1    7.4976902E-05, 1.5346565E-06, 9.8854832E-08, 7.6610414E-05,
     1    1.2385597E+00, 4.7409717E-01, 3.5921692E-02, 1.7485785E+00,
     1    2.1341437E-01, 1.0281964E-01, 6.1867919E-03, 3.2242080E-01,
     1    3.1936217E-04, 9.2577611E-07, 5.9796242E-08, 3.2034774E-04,
     1    6.4231756E-01, 1.5229228E-01, 8.8649943E-03, 8.0347484E-01,
     1    8.3360619E-02, 3.3281381E-02, 1.6530887E-03, 1.1829509E-01,
     1    1.2891439E-03,-2.1877770E-06,-1.4107339E-07, 1.2868151E-03,
     1    3.9953919E-01, 1.7920563E-02, 4.3975956E-04, 4.1789951E-01,
     1    2.3015973E-02, 3.4638656E-03, 7.8718217E-05, 2.6558556E-02,
     1    3.6361932E-03,-1.6907060E-06,-1.0916009E-07, 3.6343934E-03,
     1    2.3548935E-01, 1.2814665E-03, 1.5970227E-05, 2.3678679E-01,
     1    6.7725826E-03, 1.5736030E-04, 1.0849835E-06, 6.9310279E-03,
     1    2.9298533E-03,-2.3936038E-07,-1.5447361E-08, 2.9295985E-03,
     1    9.9686720E-02, 9.0118465E-05, 1.6984178E-06, 9.9778537E-02,
     1    1.6325595E-03, 5.3798387E-06, 4.3297988E-08, 1.6379826E-03,
     1    1.2935091E-03,-2.7219136E-08,-1.7567988E-09, 1.2934801E-03,
     1    2.3352933E-02, 2.8422167E-06, 9.0633920E-08, 2.3355865E-02,
     1    1.8683340E-04, 4.8872849E-08, 8.9289861E-10, 1.8688316E-04,
     1    3.0120435E-04,-1.4820982E-09,-9.5736073E-11, 3.0120278E-04,
     1    8.6128657E-04, 4.7425789E-09, 3.7162006E-10, 8.6129169E-04,
     1    1.6955810E-06, 8.1136599E-12, 7.3967101E-13, 1.6955898E-06,
     1    1.0862367E-05,-5.0338246E-12,-3.0339479E-13, 1.0862362E-05,
     1    1.3947267E+01, 1.1911910E+01, 1.0151370E-02, 2.5869329E+01,
     1    2.5922548E+00, 2.5150605E+00, 2.9418572E-03, 5.1102571E+00,
     1    3.7605244E-02,-7.7772409E-01, 5.1474537E-03,-7.3497139E-01,
     1    7.4681354E+00, 5.8903275E+00, 4.7619298E-03, 1.3363225E+01,
     1    1.3714480E+00, 1.3023910E+00, 1.4671310E-03, 2.6753061E+00,
     1    4.6144634E-02,-5.3487388E-01, 2.5569829E-03,-4.8617226E-01,
     1    3.8587081E+00, 2.6065621E+00, 1.9256657E-03, 6.4671959E+00,
     1    6.4418710E-01, 5.8659472E-01, 6.2185390E-04, 1.2314037E+00,
     1    1.0692654E-01,-3.4788329E-01, 1.1145611E-03,-2.3984219E-01,
     1    2.1103292E+00, 9.4254040E-01, 5.8092374E-04, 3.0534505E+00,
     1    2.5792852E-01, 2.0257508E-01, 1.8955078E-04, 4.6069316E-01,
     1    3.6686373E-01,-2.0342332E-01, 3.8537362E-04, 1.6382578E-01,
     1    1.5115259E+00, 1.5318849E-01, 7.2060118E-05, 1.6647864E+00,
     1    8.2009692E-02, 2.5299917E-02, 1.8757600E-05, 1.0732837E-01,
     1    1.0243472E+00,-6.2866865E-02, 6.4724272E-05, 9.6154509E-01,
     1    9.5611656E-01, 1.5274278E-02, 1.4351689E-05, 9.7140519E-01,
     1    2.6889108E-02, 1.6179666E-03, 2.5687425E-06, 2.8509644E-02,
     1    8.6832007E-01,-8.5832556E-03, 1.2656386E-05, 8.5974947E-01,
     1    4.1774783E-01, 1.2254238E-03, 3.4963743E-06, 4.1897675E-01,
     1    6.7342938E-03, 8.8459137E-05, 6.0262636E-07, 6.8233556E-03,
     1    4.0438196E-01,-7.5975396E-04, 2.8782640E-06, 4.0362509E-01,
     1    1.0042336E-01, 3.4049402E-05, 2.5523370E-07, 1.0045767E-01,
     1    7.8943839E-04, 1.7726372E-06, 3.9073389E-08, 7.9125010E-04,
     1    9.9366963E-02,-2.0591291E-05, 2.1267813E-07, 9.9346584E-02,
     1    3.7954084E-03, 1.4528691E-08, 5.2409711E-09, 3.7954281E-03,
     1    7.3136501E-06, 8.0153882E-10, 3.2255033E-11, 7.3144839E-06,
     1    3.7877733E-03,-2.7918604E-10, 5.2037488E-09, 3.7877783E-03,
     1    1.3946000E-01,
     1    2.6048820E-03, 1.6070835E-03, 2.4245063E+01, 1.0606264E+01,
     1    1.4501628E+02,
     1    1.1756958E-02, 6.9806384E-03, 1.1784443E+01, 4.7875990E+00,
     1    6.4678854E+01,
     1    5.2299806E-02, 3.0258211E-02, 5.3029905E+00, 1.8743871E+00,
     1    2.4801767E+01,
     1    2.1501936E-01, 1.2124582E-01, 2.1001344E+00, 5.5915055E-01,
     1    7.2865648E+00,
     1    5.6565146E-01, 2.8215104E-01, 4.1024135E-01, 6.1789287E-02,
     1    9.6493989E-01,
     1    3.9250455E-01, 1.4841896E-01, 4.1107944E-02, 4.2901853E-03,
     1    1.0188520E-01,
     1    1.4572636E-01, 3.8720096E-02, 2.9931511E-03, 2.7504108E-04,
     1    1.0640176E-02,
     1    2.6263640E-02, 4.1334847E-03, 7.3943243E-05, 7.0912032E-06,
     1    5.5237309E-04,
     1    5.3320550E-04, 2.7682475E-05, 4.2194620E-08, 6.4920290E-09,
     1    2.0204476E-06,
     1    7.7458765E+00, 4.5478960E+00, 7.8099016E-01, 1.3074763E+01,
     1    1.2750472E+00, 8.0948541E-01, 1.6252984E-01, 2.2470625E+00,
     1    1.9802779E-04, 1.0079834E-05, 1.4508922E-06, 2.0955851E-04,
     1    3.7319977E+00, 2.0373747E+00, 3.3186609E-01, 6.1012385E+00,
     1    5.7728792E-01, 3.5882745E-01, 7.0278899E-02, 1.0063943E+00,
     1    8.3284169E-04, 1.6690957E-05, 2.4031494E-06, 8.5193579E-04,
     1    1.6816112E+00, 7.7952812E-01, 1.1896913E-01, 2.5801085E+00,
     1    2.2661215E-01, 1.3502153E-01, 2.5619717E-02, 3.8725340E-01,
     1    3.4417126E-03, 8.0601390E-06, 1.1642238E-06, 3.4509369E-03,
     1    7.4992995E-01, 2.2758256E-01, 3.0588307E-02, 1.0081008E+00,
     1    7.1678983E-02, 3.6965592E-02, 6.6750952E-03, 1.1531967E-01,
     1    1.3174098E-02,-2.5111281E-05,-3.6195626E-06, 1.3145367E-02,
     1    3.9475134E-01, 2.4125499E-02, 2.4139506E-03, 4.2129079E-01,
     1    1.5766850E-02, 3.1802490E-03, 5.0107889E-04, 1.9448178E-02,
     1    3.2975581E-02,-1.6079391E-05,-2.3203938E-06, 3.2957181E-02,
     1    2.0520079E-01, 1.6433573E-03, 1.1663270E-04, 2.0696078E-01,
     1    4.1317709E-03, 1.4707353E-04, 1.7925473E-05, 4.2967699E-03,
     1    2.3179603E-02,-1.9907052E-06,-2.8720101E-07, 2.3177325E-02,
     1    7.6986299E-02, 1.0708545E-04, 5.6997679E-06, 7.7099085E-02,
     1    8.9439816E-04, 6.2026706E-06, 5.2129601E-07, 9.0112212E-04,
     1    9.0218749E-03,-2.0131628E-07,-2.9057223E-08, 9.0216445E-03,
     1    1.5443963E-02, 2.8629588E-06, 1.3377623E-07, 1.5446960E-02,
     1    8.9680153E-05, 9.1379951E-08, 3.8080033E-09, 8.9775341E-05,
     1    1.7968448E-03,-9.4432716E-09,-1.3650368E-09, 1.7968339E-03,
     1    4.2605976E-04, 2.7468992E-09, 2.5527743E-10, 4.2606276E-04,
     1    6.4499630E-07, 8.5880521E-12, 4.8638899E-13, 6.4500538E-07,
     1    4.8574108E-05,-2.2217553E-11,-3.2702832E-12, 4.8574083E-05,
     1    2.3094560E+01, 2.0968429E+01, 2.8711771E-02, 4.4091701E+01,
     1    3.7724412E+00, 3.6907432E+00, 5.7969252E-03, 7.4689814E+00,
     1    3.8689371E-02,-8.0674591E-01, 5.5620250E-03,-7.6249452E-01,
     1    1.1193874E+01, 9.5828669E+00, 1.2598162E-02, 2.0789339E+01,
     1    1.7103471E+00, 1.6520657E+00, 2.5691928E-03, 3.3649820E+00,
     1    5.0567545E-02,-5.5299266E-01, 2.7124875E-03,-4.9971263E-01,
     1    5.1202903E+00, 3.8581451E+00, 4.7508572E-03, 8.9831863E+00,
     1    6.7440132E-01, 6.3032398E-01, 9.6387937E-04, 1.3056892E+00,
     1    1.2310829E-01,-3.5513343E-01, 1.1694890E-03,-2.3085565E-01,
     1    2.4112240E+00, 1.2314780E+00, 1.3363304E-03, 3.6440383E+00,
     1    2.1895130E-01, 1.7709803E-01, 2.6161967E-04, 3.9631094E-01,
     1    4.1168225E-01,-2.0246746E-01, 4.0440085E-04, 2.0961919E-01,
     1    1.4787137E+00, 1.6718313E-01, 1.4368394E-04, 1.6460405E+00,
     1    5.6427208E-02, 1.6322265E-02, 2.2574275E-05, 7.2772047E-02,
     1    1.0280923E+00,-5.7705799E-02, 6.7557915E-05, 9.7045405E-01,
     1    8.2807776E-01, 1.4629239E-02, 1.8110174E-05, 8.4272511E-01,
     1    1.6542912E-02, 8.2076886E-04, 1.3557751E-06, 1.7365036E-02,
     1    7.6193336E-01,-7.1634642E-03, 1.3965844E-05, 7.5478387E-01,
     1    3.2113011E-01, 1.0872466E-03, 4.3733359E-06, 3.2222173E-01,
     1    3.7198072E-03, 3.7659385E-05, 1.5856914E-07, 3.7576251E-03,
     1    3.1275450E-01,-6.0321667E-04, 4.0644131E-06, 3.1215534E-01,
     1    6.6189200E-02, 2.9933563E-05, 7.5547755E-07, 6.6219889E-02,
     1    3.8189478E-04, 6.4738893E-07, 1.9207858E-08, 3.8256138E-04,
     1    6.5649137E-02,-1.7106858E-05, 7.3210669E-07, 6.5632762E-02,
     1    1.8748561E-03, 2.1930447E-08, 1.4128944E-08, 1.8748922E-03,
     1    2.7966043E-06, 2.7994021E-10, 2.5605549E-10, 2.7971403E-06,
     1    1.8719141E-03,-1.0303972E-08, 1.3862818E-08, 1.8719177E-03,
     1    1.1560477E-01,
     1    3.1909185E-03, 1.9532932E-03, 3.4733122E+01, 1.5877446E+01,
     1    2.2012545E+02,
     1    1.4024062E-02, 8.2753090E-03, 1.5617475E+01, 6.7257637E+00,
     1    8.8805560E+01,
     1    6.0020239E-02, 3.4519592E-02, 6.4173019E+00, 2.4498614E+00,
     1    3.0404496E+01,
     1    2.3243375E-01, 1.2999499E-01, 2.2778615E+00, 6.6750544E-01,
     1    7.7914884E+00,
     1    5.4992465E-01, 2.7035099E-01, 3.8526613E-01, 6.4426108E-02,
     1    8.5269107E-01,
     1    3.4622762E-01, 1.2833111E-01, 3.4601234E-02, 4.0075020E-03,
     1    7.8903427E-02,
     1    1.1868570E-01, 3.0812615E-02, 2.3199306E-03, 2.3698518E-04,
     1    7.6405937E-03,
     1    1.9491222E-02, 2.9901608E-03, 5.2394816E-05, 5.6053533E-06,
     1    3.7085772E-04,
     1    3.3526317E-04, 1.6935391E-05, 2.5451660E-08, 4.1835755E-09,
     1    1.1724106E-06,
     1    1.3888328E+01, 7.5072828E+00, 2.1262250E+00, 2.3521836E+01,
     1    1.9673290E+00, 1.1336161E+00, 3.7603035E-01, 3.4769755E+00,
     1    1.2170098E-03, 5.8419648E-05, 1.0326032E-05, 1.2857555E-03,
     1    6.1962000E+00, 3.1627604E+00, 8.5845287E-01, 1.0217413E+01,
     1    7.9910962E-01, 4.5500074E-01, 1.4921925E-01, 1.4033296E+00,
     1    5.0370433E-03, 9.2099295E-05, 1.6285345E-05, 5.1454279E-03,
     1    2.5437513E+00, 1.1445575E+00, 2.8999281E-01, 3.9783016E+00,
     1    2.7703462E-01, 1.5288877E-01, 4.9365561E-02, 4.7928896E-01,
     1    2.0306780E-02, 3.5566820E-05, 6.3129732E-06, 2.0348660E-02,
     1    1.0037757E+00, 3.0865879E-01, 7.0310754E-02, 1.3827452E+00,
     1    7.5196811E-02, 3.6411854E-02, 1.1474693E-02, 1.2308336E-01,
     1    7.4362435E-02,-1.4373018E-04,-2.5444363E-05, 7.4193261E-02,
     1    4.4969141E-01, 2.9305829E-02, 5.4275991E-03, 4.8442483E-01,
     1    1.3699847E-02, 2.5840765E-03, 7.7980229E-04, 1.7063726E-02,
     1    1.6883333E-01,-7.9932598E-05,-1.4164747E-05, 1.6873923E-01,
     1    2.0735496E-01, 1.8440344E-03, 2.8197762E-04, 2.0948098E-01,
     1    3.1862818E-03, 1.0735031E-04, 3.0793520E-05, 3.3244257E-03,
     1    1.0657116E-01,-8.9351430E-06,-1.5827221E-06, 1.0656065E-01,
     1    7.0123292E-02, 1.1259421E-04, 1.4991106E-05, 7.0250878E-02,
     1    6.2406264E-04, 4.4210565E-06, 1.1890735E-06, 6.2967277E-04,
     1    3.7638049E-02,-8.2719163E-07,-1.4649797E-07, 3.7637075E-02,
     1    1.2394596E-02, 2.7899685E-06, 3.3317361E-07, 1.2397719E-02,
     1    5.5724442E-05, 6.9752456E-08, 1.6225341E-08, 5.5810420E-05,
     1    6.6521639E-03,-3.4753529E-08,-6.1472905E-09, 6.6521230E-03,
     1    2.7180980E-04, 2.2872644E-09, 2.7401280E-10, 2.7181236E-04,
     1    3.2928942E-07, 2.5026453E-11, 1.1837319E-12, 3.2931563E-07,
     1    1.4415290E-04,-6.6848576E-11,-1.1873423E-11, 1.4415282E-04,
     1    3.3337681E+01, 3.1196145E+01, 4.7283311E-02, 6.4581110E+01,
     1    4.6936420E+00, 4.6288712E+00, 7.7669888E-03, 9.3302803E+00,
     1    3.9604605E-02,-8.0880658E-01, 5.0727437E-03,-7.6412923E-01,
     1    1.4965327E+01, 1.3358527E+01, 1.9562642E-02, 2.8343416E+01,
     1    1.9092541E+00, 1.8636133E+00, 3.1056228E-03, 3.7759730E+00,
     1    5.4501296E-02,-5.4993441E-01, 2.4878059E-03,-4.9294531E-01,
     1    6.2481200E+00, 4.9901309E+00, 6.9036422E-03, 1.1245154E+01,
     1    6.6529673E-01, 6.3039578E-01, 1.0376593E-03, 1.2967302E+00,
     1    1.3707935E-01,-3.5023479E-01, 1.0792868E-03,-2.1207615E-01,
     1    2.6351863E+00, 1.4493662E+00, 1.7917500E-03, 4.0863443E+00,
     1    1.8655876E-01, 1.5277470E-01, 2.4500173E-04, 3.3957846E-01,
     1    4.4673266E-01,-1.9651124E-01, 3.7349055E-04, 2.5059492E-01,
     1    1.4380477E+00, 1.7104203E-01, 1.7127855E-04, 1.6092610E+00,
     1    4.1817795E-02, 1.1399123E-02, 1.7476809E-05, 5.3234395E-02,
     1    1.0179112E+00,-5.2166912E-02, 6.0701124E-05, 9.6580500E-01,
     1    7.3200992E-01, 1.3291595E-02, 1.8159480E-05, 7.4531967E-01,
     1    1.1250000E-02, 4.9352399E-04, 8.5073228E-07, 1.1744375E-02,
     1    6.7908693E-01,-5.8457434E-03, 1.2226787E-05, 6.7325341E-01,
     1    2.5898975E-01, 9.0044413E-04, 3.7718237E-06, 2.5989397E-01,
     1    2.3218354E-03, 2.0537930E-05, 6.3388733E-08, 2.3424368E-03,
     1    2.5316309E-01,-4.4898240E-04, 3.4260734E-06, 2.5271754E-01,
     1    4.7492682E-02, 2.2575363E-05, 6.0484650E-07, 4.7515862E-02,
     1    2.1461785E-04, 3.1874276E-07, 4.2478598E-09, 2.1494084E-04,
     1    4.7173218E-02,-1.1521084E-05, 5.9395224E-07, 4.7162291E-02,
     1    1.0795287E-03, 1.5153679E-08, 1.2988184E-08, 1.0795568E-03,
     1    1.3040619E-06, 1.0973624E-10, 3.9989157E-11, 1.3042116E-06,
     1    1.0781466E-03,-6.4030869E-09, 1.2942186E-08, 1.0781532E-03,
     1    0.0000000E+00/