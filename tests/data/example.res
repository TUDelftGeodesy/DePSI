===============================================
MASTER RESULTFILE:                      
Created by:                             G.Mulder TU Delft
DVersion:                               Version (2015) (For TOPSAR)
FFTW library:                           used
VECLIB library:                         not used
LAPACK library:                         not used
Compiled at:                            XXXXXXXX
By GUN gcc:                             XXXXXXXX
===============================================



Start_process_control
readfiles:		1
precise_orbits:		1
crop:		1
sim_amplitude:		0
master_timing:		0
oversample:		0
resample:		1
filt_azi:		0
filt_range:		0
NOT_USED:		0
End_process_control



******************************************************************* 
*_Start_readfiles:
******************************************************************* 
Volume_file:                                              dummy
Volume_ID:                                                162245
Volume_identifier:                                        dummy
Volume_set_identifier:                                    dummy
Number of records in ref. file:                           dummy
SAR_PROCESSOR:                                            Sentinel-1B
SWATH:                                                    IW2
PASS:                                                     Descending
IMAGE_MODE:                                               IW
polarisation:                                             VV
Product type specifier:                                   S1B
Logical volume generating facility:                       dummy
Location and date/time of product creation:               dummy
Number_of_lines_Swath:                                    13680
number_of_pixels_Swath:                                   27028
rangePixelSpacing:                                        2.329562e+00
azimuthPixelSpacing:                                      1.389183e+01
total_Burst:                                              9
Burst_number_index:                                       4
RADAR_FREQUENCY (HZ):                                     5.405000454334350e+09
Scene identification:                                     Orbit: 20888
Scene location:                                           lat: 5.285133749958482e+01 lon:6.883094329390478e+00
Sensor platform mission identifer:                        S1B
Scene_center_heading:                                     -1.645706245495673e+02
Scene_centre_latitude:                                    52.34636275554789
Scene_centre_longitude:                                   5.952591327951736
Radar_wavelength (m):                                     0.055465760
Azimuth_steering_rate (deg/s):                            9.798633249999998e-01
Pulse_Repetition_Frequency_raw_data(TOPSAR):              1.451627112193990e+03
First_pixel_azimuth_time (UTC):                           2020-Mar-28 05:49:19.933511
Pulse_Repetition_Frequency (computed, Hz):                4.864863102995529e+02
Azimuth_time_interval (s):                                2.055556299999998e-03
Total_azimuth_band_width (Hz):                            3.130000000000000e+02
Weighting_azimuth:                                        Hamming
Range_time_to_first_pixel (2way) (ms):                    5.641912410057507
Range_sampling_rate (computed, MHz):                      64.345238126
Total_range_band_width (MHz):                             48.300000000
Weighting_range:                                          Hamming
DC_reference_azimuth_time:                                2020-Mar-28 05:49:35.293439
DC_reference_range_time:                                  5.357127211302541e-03
Xtrack_f_DC_constant (Hz, early edge):                    -4.190346e+00
Xtrack_f_DC_linear (Hz/s, early edge):                    5.519143e+04
Xtrack_f_DC_quadratic (Hz/s/s, early edge):               -4.465147e+07
FM_reference_azimuth_time:                                2020-Mar-28 05:49:35.285819
FM_reference_range_time:                                  5.641912410057507e-03
FM_polynomial_constant_coeff (Hz, early edge):            -2.190052980863128e+03
FM_polynomial_linear_coeff (Hz/s, early edge):            4.027293936151257e+05
FM_polynomial_quadratic_coeff (Hz/s/s, early edge):       -6.554460524600121e+07
Datafile:                                                 s1b-iw2-slc-vv-20200328t054925-20200328t054950-020888-0279c5-005.tiff
Dataformat:                                               tiff
Number_of_lines_original:                                 14065
Number_of_pixels_original:                                49843
Scene_ul_corner_latitude:                                 52.35671298100559
Scene_ur_corner_latitude:                                 52.51325612534782
Scene_lr_corner_latitude:                                 52.326518849251606
Scene_ll_corner_latitude:                                 52.17013287457889
Scene_ul_corner_longitude:                                6.729190990864142
Scene_ur_corner_longitude:                                5.275280947023443
Scene_lr_corner_longitude:                                5.223952756994337
Scene_ll_corner_longitude:                                6.6714418521468675
deramp:                                                   0
reramp:                                                   0
ESD_correct:                                              0
First_line (w.r.t. output_image):                         1
Last_line (w.r.t. output_image):                          14065
First_pixel (w.r.t. output_image):                        1
Last_pixel (w.r.t. output_image):                         49843
Number_of_pixels_output_image:                            49843
Number_of_lines_output_image:                             14065
******************************************************************* 
* End_readfiles:_NORMAL
******************************************************************* 



******************************************************************* 
*_Start_precise_orbits:
******************************************************************* 
 t(s)    X(m)                 Y(m)                 Z(m)                
NUMBER_OF_DATAPOINTS:       200
 20865   3626008.9058349766   1060535.1214947468   5973108.5284061     
 20866   3632520.9210890783   1060068.1662619063   5969242.514157453   
 20867   3639028.7999242027   1059599.0743025467   5965369.765174575   
 20868   3645532.5347308316   1059127.8467443837   5961490.285810932   
 20869   3652032.117899445    1058654.4847151318   5957604.0804199865  
 20870   3658527.541820523    1058178.9893425072   5953711.153355206   
 20871   3665018.7988845482   1057701.3617542249   5949811.508970056   
 20872   3671505.8814819995   1057221.6030780002   5945905.151618      
 20873   3677988.7820113543   1056739.7144473724   5941992.085665138   
 20874   3684467.4929030617   1056255.697019175    5938072.315528107   
 20875   3690942.0065955706   1055769.5519560657   5934145.845636174   
 20876   3697412.3155273246   1055281.2804207022   5930212.680418607   
 20877   3703878.412136771    1054790.8835757417   5926272.824304678   
 20878   3710340.288862356    1054298.362583842    5922326.281723655   
 20879   3716797.9381425246   1053803.7186076604   5918373.057104806   
 20880   3723251.352415724    1053306.9528098544   5914413.154877401   
 20881   3729700.5241204007   1052808.066353082    5910446.579470709   
 20882   3736145.445695       1052307.0604         5906473.335314      
 20883   3742586.1095861145   1051803.936119055    5902493.426849097   
 20884   3749022.50827292     1051298.6947018448   5898506.858568038   
 20885   3755454.634242741    1050791.337345756    5894513.634975417   
 20886   3761882.4799828986   1050281.8652481749   5890513.76057582    
 20887   3768306.037980716    1049770.2796064883   5886507.239873844   
 20888   3774725.3007235164   1049256.5816180827   5882494.077374081   
 20889   3781140.260698622    1048740.7724803446   5878474.277581124   
 20890   3787550.9103933554   1048222.8533906601   5874447.844999562   
 20891   3793957.242295041    1047702.8255464166   5870414.784133989   
 20892   3800359.2488909997   1047180.6901449999   5866375.099488999   
 20893   3806756.9226768473   1046656.4483895535   5862328.795581659   
 20894   3813150.256181361    1046130.101506247    5858275.876978951   
 20895   3819539.241941615    1045601.6507270071   5854216.348260333   
 20896   3825923.872494677    1045071.0972837596   5850150.214005262   
 20897   3832304.140377618    1044538.4424084312   5846077.478793196   
 20898   3838680.0381275103   1044003.6873329486   5841998.147203597   
 20899   3845051.5582814245   1043466.8332892377   5837912.22381592    
 20900   3851418.693376429    1042927.8815092247   5833819.713209621   
 20901   3857781.435949598    1042386.8332248369   5829720.619964164   
 20902   3864139.7785380003   1041843.689668       5825614.948659002   
 20903   3870493.713687152    1041298.452076362    5821502.703885988   
 20904   3876843.233976346    1040751.1217104553   5817383.890286535   
 20905   3883188.331993322    1040201.6998365341   5813258.512514454   
 20906   3889529.000325815    1039650.1877208519   5809126.5752235465  
 20907   3895865.231561565    1039096.5866296625   5804988.083067623   
 20908   3902197.0182883106   1038540.8978292206   5800843.040700491   
 20909   3908524.353093789    1037983.1225857794   5796691.452775955   
 20910   3914847.2285657367   1037423.2621655926   5792533.323947822   
 20911   3921165.637291895    1036861.3178349149   5788368.658869902   
 20912   3927479.57186        1036297.29086        5784197.462196      
 20913   3933789.0248663886   1035731.1825127867   5780019.738592234   
 20914   3940093.9889417915   1035162.9940879558   5775835.492773978   
 20915   3946394.4567255382   1034592.726885874    5771644.729468911   
 20916   3952690.420856956    1034020.3822069062   5767447.4534047125  
 20917   3958981.8739753747   1033445.961351419    5763243.669309068   
 20918   3965268.808720124    1032869.4656197784   5759033.381909658   
 20919   3971551.2177305324   1032290.8963123502   5754816.595934166   
 20920   3977829.093645928    1031710.2547295      5750593.316110268   
 20921   3984102.4291056413   1031127.5421715949   5746363.547165654   
 20922   3990371.2167490004   1030542.7599390001   5742127.293828001   
 20923   3996635.4492240753   1029955.9093377313   5737884.560837217   
 20924   4002895.1192138973   1029366.9916964008   5733635.352982106   
 20925   4009150.219410239    1028776.0083492714   5729379.675063699   
 20926   4015400.7425048705   1028182.9606306046   5725117.531883024   
 20927   4021646.681189565    1027587.8498746626   5720848.9282411095  
 20928   4027888.0281560947   1026990.6774157077   5716573.868938988   
 20929   4034124.7760962313   1026391.4445880019   5712292.358777686   
 20930   4040356.9177017454   1025790.152725807    5708004.402558234   
 20931   4046584.445664411    1025186.8031633858   5703710.005081661   
 20932   4052807.3526759995   1024581.3972349998   5699409.171148999   
 20933   4059025.631437168    1023973.9362805228   5695101.905573415   
 20934   4065239.2746841162   1023364.4216622727   5690788.213216632   
 20935   4071448.275161931    1022752.8547481797   5686468.098952519   
 20936   4077652.6256156936   1022139.2369061726   5682141.567654935   
 20937   4083852.3187904926   1021523.5695041814   5677808.624197748   
 20938   4090047.3474314124   1020905.8539101359   5673469.273454821   
 20939   4096237.704283538    1020286.0914919653   5669123.520300019   
 20940   4102423.3820919534   1019664.2836175992   5664771.369607204   
 20941   4108604.3736017463   1019040.4316549677   5660412.826250244   
 20942   4114780.671558       1018414.5369719999   5656047.895103      
 20943   4120952.2687148326   1017786.6009421963   5651676.581051387   
 20944   4127119.157862484    1017156.6249613394   5647298.889029507   
 20945   4133281.3318002294   1016524.6104307835   5642914.823983513   
 20946   4139438.783327341    1015890.5587518814   5638524.390859552   
 20947   4145591.505243091    1015254.4713259866   5634127.594603778   
 20948   4151739.4903467554   1014616.349554453    5629724.440162343   
 20949   4157882.731437605    1013976.1948386341   5625314.932481397   
 20950   4164021.221314912    1013334.0085798831   5620899.076507089   
 20951   4170154.952777954    1012689.792179554    5616476.877185574   
 20952   4176283.918626       1012043.5470389999   5612048.339463001   
 20953   4182408.111667505    1011395.2745651099   5607613.46829749    
 20954   4188527.5247476287   1010744.9761869125   5603172.2686950285  
 20955   4194642.150720717    1010092.653338973    5598724.745673577   
 20956   4200751.982441109    1009438.3074558546   5594270.904251089   
 20957   4206857.012763147    1008781.9399721222   5589810.749445522   
 20958   4212957.234541174    1008123.5523223401   5585344.286274833   
 20959   4219052.640629531    1007463.1459410726   5580871.51975698    
 20960   4225143.223882559    1006800.7222628836   5576392.454909917   
 20961   4231228.9771546      1006136.282722338    5571907.096751606   
 20962   4237309.893299999    1005469.8287539999   5567415.4503        
 20963   4243385.965182407    1004801.3617979272   5562917.520584935   
 20964   4249457.185702727    1004130.8833161537   5558413.312683762   
 20965   4255523.547771178    1003458.3947762066   5553902.831685713   
 20966   4261585.04429797     1002783.8976456127   5549386.082680012   
 20967   4267641.668193322    1002107.3933918996   5544863.07075589    
 20968   4273693.412367447    1001428.8834825947   5540333.801002579   
 20969   4279740.269730563    1000748.369385225    5535798.278509305   
 20970   4285782.2331928825   1000065.8525673177   5531256.508365296   
 20971   4291819.295664623    999381.3344964004    5526708.495659785   
 20972   4297851.450056       998694.8166400002    5522154.2454820005  
 20973   4303878.689286696    998006.3004710989    5517593.762932958   
 20974   4309901.006314264    997315.7874844964    5513027.053160823   
 20975   4315918.3941057315   996623.2791804478    5508454.121325554   
 20976   4321930.845628115    995928.7770592076    5503874.972587101   
 20977   4327938.3538484415   995232.2826210299    5499289.61210542    
 20978   4333940.9117337335   994533.7973661701    5494698.045040466   
 20979   4339938.512251013    993833.3227948827    5490100.276552193   
 20980   4345931.148367301    993130.8604074217    5485496.311800553   
 20981   4351918.813049623    992426.4117040427    5480886.155945506   
 20982   4357901.499265       991719.978185        5476269.814147001   
 20983   4363879.199990065    991011.5613559601    5471647.291576694   
 20984   4369851.908239886    990301.162744237     5467018.593453034   
 20985   4375819.617039141    989588.7838825574    5462383.725006172   
 20986   4381782.319412508    988874.4263036469    5457742.691466256   
 20987   4387740.008384664    988158.0915402316    5453095.498063436   
 20988   4393692.676980286    987439.7811250382    5448442.15002786    
 20989   4399640.318224053    986719.4965907921    5443782.652589681   
 20990   4405582.925140643    985997.2394702195    5439117.010979042   
 20991   4411520.490754733    985273.0112960468    5434445.2304261     
 20992   4417453.008091001    984546.813601        5429767.316161      
 20993   4423380.47018388     983818.6479231772    5425083.273425499   
 20994   4429302.870106832    983088.5158221633    5420393.10750778    
 20995   4435220.200943074    982356.4188629155    5415696.8237076355  
 20996   4441132.455775818    981622.35861039      5410994.427324852   
 20997   4447039.627688284    980886.3366295442    5406285.923659219   
 20998   4452941.709763688    980148.3544853351    5401571.318010531   
 20999   4458838.695085245    979408.413742719     5396850.615678575   
 21000   4464730.576736171    978666.5159666532    5392123.82196314    
 21001   4470617.347799684    977922.6627220946    5387390.942164019   
 21002   4476499.001359       977176.855574        5382651.981581001   
 21003   4482375.53050723     976429.0960926536    5377906.945525376   
 21004   4488246.928377065    975679.3858696461    5373155.839354444   
 21005   4494113.188111095    974927.7265018965    5368398.668437008   
 21006   4499974.302851904    974174.1195863227    5363635.438141864   
 21007   4505830.2657420775   973418.5667198428    5358866.153837814   
 21008   4511681.0699242065   972661.0694993759    5354090.820893661   
 21009   4517526.708540877    971901.62952184      5349309.444678201   
 21010   4523367.174734673    971140.2483841529    5344522.030560236   
 21011   4529202.461648187    970376.9276832335    5339728.583908571   
 21012   4535032.562424       969611.6690159999    5334929.110092001   
 21013   4540857.470214752    968844.4739846535    5330123.614490735   
 21014   4546677.178213274    968075.344212525     5325312.102530613   
 21015   4552491.67962245     967304.2813282284    5320494.579648876   
 21016   4558300.967645158    966531.286960377     5315671.051282766   
 21017   4564105.035484281    965756.3627375851    5310841.52286953    
 21018   4569903.876342704    964979.5102884667    5306005.999846413   
 21019   4575697.483423306    964200.7312416353    5301164.487650657   
 21020   4581485.849928966    963420.0272257042    5296316.991719506   
 21021   4587268.9690625705   962637.3998692879    5291463.517490205   
 21022   4593046.834027       961852.8508009999    5286604.0704        
 21023   4598819.438035309    961066.381654694     5281738.655897432   
 21024   4604586.774341254    960277.9940851827    5276867.279476252   
 21025   4610348.836208762    959487.6897525189    5271989.946641516   
 21026   4616105.616901758    958695.4703167543    5267106.6628982695  
 21027   4621857.109684173    957901.3374379418    5262217.433751567   
 21028   4627603.307819935    957105.2927761342    5257322.26470646    
 21029   4633344.20457297     956307.3379913839    5252421.161267998   
 21030   4639079.793207206    955507.4747437427    5247514.128941232   
 21031   4644810.066986575    954705.704693264     5242601.173231217   
 21032   4650535.019174999    953902.0294999998    5237682.299643      
 21033   4656254.643046729    953096.4508291963    5232757.513692834   
 21034   4661968.931917276    952288.9703668719    5227826.8209417565  
 21035   4667677.879112476    951479.5898042393    5222890.226962009   
 21036   4673381.477958157    950668.3108325104    5217947.737325828   
 21037   4679079.721780152    949855.1351428975    5212999.357605451   
 21038   4684772.6039042985   949040.0644266128    5208045.09337312    
 21039   4690460.117656423    948223.1003748685    5203084.95020107    
 21040   4696142.25636236     947404.2446788766    5198118.93366154    
 21041   4701819.013347941    946583.4990298498    5193147.049326771   
 21042   4707490.381939001    945760.865119        5188169.302769      
 21043   4713156.355471818    944936.3446426855    5183185.699571556   
 21044   4718816.9273244655   944109.939317849     5178196.245362129   
 21045   4724472.090885463    943281.6508665793    5173200.9457795005  
 21046   4730121.839543333    942451.4810109645    5168199.806462449   
 21047   4735766.1666865945   941619.4314730937    5163192.8330497565  
 21048   4741405.065703769    940785.5039750553    5158180.031180205   
 21049   4747038.529983376    939949.7002389382    5153161.406492572   
 21050   4752666.552913935    939112.0219868306    5148136.96462564    
 21051   4758289.127883971    938272.4709408218    5143106.711218189   
 21052   4763906.248282       937431.048823        5138070.651909001   
 21053   4769517.9075071355   936587.7573605565    5133028.7923478475  
 21054   4775124.099000837    935742.5983010919    5127981.138228473   
 21055   4780724.816215156    934895.5733973098    5122927.695255617   
 21056   4786320.052602145    934046.6844019126    5117868.469134014   
 21057   4791909.801613852    933195.9330676039    5112803.465568401   
 21058   4797494.0567023335   932343.321147087     5107732.690263517   
 21059   4803072.811319638    931488.8503930651    5102656.148924097   
 21060   4808646.058917816    930632.5225582408    5097573.84725488    
 21061   4814213.79294892     929774.3393953182    5092485.790960601   
 21062   4819776.006865001    928914.3026569999    5087391.985746      
 21063   4825332.694128812    928052.4141010405    5082292.437326704   
 21064   4830883.848245909    927188.6755053991    5077187.151461916   
******************************************************************* 
* End_precise_orbits:_NORMAL
******************************************************************* 



******************************************************************* 
*_Start_crop:
******************************************************************* 
Data_output_file:                          2020-03-28.raw
Data_output_format:                        complex_short
First_line (w.r.t. original_image):        1
Last_line (w.r.t. original_image):         14065
First_pixel (w.r.t. original_image):       1
Last_pixel (w.r.t. original_image):        49843
******************************************************************* 
* End_crop:_NORMAL
******************************************************************* 



******************************************************************* 
*_Start_resample:
******************************************************************* 
Normalization_Lines:                        1 1469
Normalization_Pixels:                       1 23567
Shifted azimuth spectrum:                   0
Data_output_file:                           slave_rsmp_reramped.raw
Data_output_format:                         complex_real4
Interpolation kernel:                       12 point raised cosine kernel
First_line (w.r.t. original_master):        1
Last_line (w.r.t. original_master):         14065
First_pixel (w.r.t. original_master):       1
Last_pixel (w.r.t. original_master):        49843
******************************************************************* 
* End_resample:_NORMAL
******************************************************************* 
