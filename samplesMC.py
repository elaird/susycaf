import samples
from configuration import srm
mc = samples.SampleHolder()

#PY 6 QCD
tag1 = "/bm409/ICF/automated/2010_10_11_23_08_00/"
tag2 = "/bbetchar/ICF/automated/2010_10_18_17_55_22/"
tag3 = "/bbetchar/ICF/automated/2010_10_17_20_12_27/"
tag4 = "/bm409/ICF/automated/2010_10_11_23_03_46/"
tagB = "_TuneZ2_7TeV_pythia6.Fall10-START38_V12-v1.GEN-SIM-RECO/"
tagC = "_TuneZ2_7TeV_pythia6.Fall10-START38_V12-v2.GEN-SIM-RECO/"

mc.add('qcd_py6_pt_0to5'      , '%s/%s/QCD_Pt_0to5%s")'      %(srm,tag1,tagB), xs = 4.844e+10)
mc.add('qcd_py6_pt_5to15'     , '%s/%s/QCD_Pt_5to15%s")'     %(srm,tag1,tagB), xs = 3.675e+10)
mc.add('qcd_py6_pt_15to30'    , '%s/%s/QCD_Pt_15to30%s")'    %(srm,tag1,tagB), xs = 8.159e+08)
mc.add('qcd_py6_pt_30to50'    , '%s/%s/QCD_Pt_30to50%s")'    %(srm,tag1,tagB), xs = 5.312e+07)
mc.add('qcd_py6_pt_50to80'    , '%s/%s/QCD_Pt_50to80%s")'    %(srm,tag1,tagB), xs = 6.359e+06)
mc.add('qcd_py6_pt_80to120'   , '%s/%s/QCD_Pt_80to120%s")'   %(srm,tag1,tagB), xs = 7.843e+05)
mc.add('qcd_py6_pt_120to170'  , '%s/%s/QCD_Pt_120to170%s")'  %(srm,tag1,tagB), xs = 1.151e+05)
mc.add('qcd_py6_pt_170to300'  , '%s/%s/QCD_Pt_170to300%s")'  %(srm,tag1,tagB), xs = 2.426e+04)
mc.add('qcd_py6_pt_300to470'  , '%s/%s/QCD_Pt_300to470%s", alwaysUseLastAttempt = True)'%(srm,tag2,tagB), xs = 1.168e+03)
mc.add('qcd_py6_pt_470to600'  , '%s/%s/QCD_Pt_470to600%s")'  %(srm,tag2,tagB), xs = 7.022e+01)
mc.add('qcd_py6_pt_600to800'  , '%s/%s/QCD_Pt_600to800%s")'  %(srm,tag3,tagC), xs = 1.555e+01)
mc.add('qcd_py6_pt_800to1000' , '%s/%s/QCD_Pt_800to1000%s")' %(srm,tag3,tagB), xs = 1.844e+00)
mc.add('qcd_py6_pt_1000to1400', '%s/%s")'                    %(srm,tag3)     , xs = 3.321e-01)
mc.add('qcd_py6_pt_1400to1800', '%s/%s/QCD_Pt_1400to1800%s")'%(srm,tag4,tagB), xs = 1.087e-02)
mc.add('qcd_py6_pt_1800'      , '%s/%s/QCD_Pt_1800%s")'      %(srm,tag4,tagB), xs = 3.575e-04)
            
#PY 8 QCD
tag1 = "bbetchar/ICF/automated/2010_10_06_02_22_48"
tag2 = "Tune1_7TeV_pythia8.Fall10-START38_V12-v1.GEN-SIM-RECO"

mc.add('qcd_py8_pt0to15'     , '%s/%s/QCD_Pt_0to15_%s/")'     %(srm,tag1,tag2), xs = 2.119e+12 )
mc.add('qcd_py8_pt15to30'    , '%s/%s/QCD_Pt_15to30_%s/")'    %(srm,tag1,tag2), xs = 7.883e+08 )
mc.add('qcd_py8_pt30to50'    , '%s/%s/QCD_Pt_30to50_%s/")'    %(srm,tag1,tag2), xs = 5.033e+07 )
mc.add('qcd_py8_pt50to80'    , '%s/%s/QCD_Pt_50to80_%s/")'    %(srm,tag1,tag2), xs = 6.027e+06 )
mc.add('qcd_py8_pt80to120'   , '%s/%s/QCD_Pt_80to120_%s/")'   %(srm,tag1,tag2), xs = 7.510e+05 )
mc.add('qcd_py8_pt120to170'  , '%s/%s/QCD_Pt_120to170_%s/")'  %(srm,tag1,tag2), xs = 1.119e+05 )
mc.add('qcd_py8_pt170to300'  , '%s/%s/QCD_Pt_170to300_%s/")'  %(srm,tag1,tag2), xs = 2.409e+04 )
mc.add('qcd_py8_pt300to470'  , '%s/%s/QCD_Pt_300to470_%s/")'  %(srm,tag1,tag2), xs = 1.201e+03 )
mc.add('qcd_py8_pt470to600'  , '%s/%s/QCD_Pt_470to600_%s/")'  %(srm,tag1,tag2), xs = 7.416e+01 )
mc.add('qcd_py8_pt600to800'  , '%s/%s/QCD_Pt_600to800_%s/")'  %(srm,tag1,tag2), xs = 1.672e+01 )
mc.add('qcd_py8_pt800to1000' , '%s/%s/QCD_Pt_800to1000_%s/")' %(srm,tag1,tag2), xs = 2.002e+00 )
mc.add('qcd_py8_pt1000to1400', '%s/%s/QCD_Pt_1000to1400_%s/")'%(srm,tag1,tag2), xs = 3.625e-01 )
mc.add('qcd_py8_pt1400to1800', '%s/%s/QCD_Pt_1400to1800_%s/")'%(srm,tag1,tag2), xs = 1.180e-02 )
mc.add('qcd_py8_pt1800'      , '%s/%s/QCD_Pt_1800_%s/")'      %(srm,tag1,tag2), xs = 3.802e-04 )

##PY6 mu-enriched QCD
#tag1 = "/bbetchar//ICF/automated/2010_10_29_20_17_55/"
#mc.add('qcd_py6_mu_pt20to30'  , '%s/%s/QCD_Mu_Pt20to30.Summer10-START36_V9_S09-v1.GEN-SIM-RECO/'  %(srm,tag1), xs = None)
#mc.add('qcd_py6_mu_pt30to50'  , '%s/%s/QCD_Mu_Pt30to50.Summer10-START36_V9_S09-v1.GEN-SIM-RECO/'  %(srm,tag1), xs = None)
#mc.add('qcd_py6_mu_pt50to80'  , '%s/%s/QCD_Mu_Pt50to80.Summer10-START36_V9_S09-v1.GEN-SIM-RECO/'  %(srm,tag1), xs = None)
#mc.add('qcd_py6_mu_pt80to120' , '%s/%s/QCD_Mu_Pt80to120.Summer10-START36_V9_S09-v1.GEN-SIM-RECO/' %(srm,tag1), xs = None)
#mc.add('qcd_py6_mu_pt120to170', '%s/%s/QCD_Mu_Pt120to170.Summer10-START36_V9_S09-v1.GEN-SIM-RECO/'%(srm,tag1), xs = None)
#mc.add('qcd_py6_mu_pt170toInf', '%s/%s/QCD_Mu_Pt170toInf.Summer10-START36_V9_S09-v1.GEN-SIM-RECO/'%(srm,tag1), xs = None)
       
###############
##V12 samples##
###############

#PY 6 QCD
mc.add("v12_qcd_py6_pt15",  '%s/gouskos//ICF/automated/2010_06_24_18_00_11/")'%srm, xs = 8.762e+08, ptHatMin =   15 )
mc.add("v12_qcd_py6_pt30",  '%s/gouskos//ICF/automated/2010_06_24_18_09_51/")'%srm, xs = 6.041e+07, ptHatMin =   30 )
mc.add("v12_qcd_py6_pt80",  '%s/gouskos//ICF/automated/2010_07_06_00_55_17/")'%srm, xs = 9.238e+05, ptHatMin =   80 )
mc.add("v12_qcd_py6_pt170", '%s/gouskos//ICF/automated/2010_07_06_01_33_23/")'%srm, xs = 2.547e+04, ptHatMin =  170 )
mc.add("v12_qcd_py6_pt300", '%s/gouskos//ICF/automated/2010_07_09_19_13_09/")'%srm, xs = 1.256e+03, ptHatMin =  300 )
mc.add("v12_qcd_py6_pt470", '%s/gouskos//ICF/automated/2010_07_10_04_22_06/")'%srm, xs = 8.798e+01, ptHatMin =  470 )
mc.add("v12_qcd_py6_pt800", '%s/gouskos//ICF/automated/2010_07_10_04_37_56/")'%srm, xs = 2.186e+00, ptHatMin =  800 )
mc.add("v12_qcd_py6_pt1400",'%s/gouskos//ICF/automated/2010_07_10_04_47_48/")'%srm, xs = 1.122e-02, ptHatMin = 1400 )
mc.adjustOverlappingSamples( ["v12_qcd_py6_pt%d"%i for i in [15,30,80,170,300,470,800,1400] ] )

#PY 8 QCD
py8Dir = "/bm409//ICF/automated/2010_07_26_10_23_33/"
py8Gt1  = "Summer10-START36_V10_S09-v1.GEN-SIM-RECO"
py8Gt2  = "Summer10-START36_V10_S09-v2.GEN-SIM-RECO"
mc.add("v12_qcd_py8_pt0to15",       '%s/%s/QCD_Pt-0to15_7TeV-pythia8.%s/")'     %(srm,py8Dir,py8Gt1), xs = 2.117e+12)
mc.add("v12_qcd_py8_pt15to20",      '%s/%s/QCD_Pt-15to20_7TeV-pythia8.%s/")'    %(srm,py8Dir,py8Gt1), xs = 5.638e+08)
mc.add("v12_qcd_py8_pt20to30",      '%s/%s/QCD_Pt-20to30_7TeV-pythia8.%s/")'    %(srm,py8Dir,py8Gt1), xs = 2.264e+08)
mc.add("v12_qcd_py8_pt30to50",      '%s/%s/QCD_Pt-30to50_7TeV-pythia8.%s/")'    %(srm,py8Dir,py8Gt2), xs = 5.018e+07)
mc.add("v12_qcd_py8_pt50to80",      '%s/%s/QCD_Pt-50to80_7TeV-pythia8.%s/")'    %(srm,py8Dir,py8Gt1), xs = 6.035e+06)
mc.add("v12_qcd_py8_pt80to120",     '%s/%s/QCD_Pt-80to120_7TeV-pythia8.%s/")'   %(srm,py8Dir,py8Gt1), xs = 7.519e+05)
mc.add("v12_qcd_py8_pt120to170",    '%s/%s/QCD_Pt-120to170_7TeV-pythia8.%s/")'  %(srm,py8Dir,py8Gt1), xs = 1.120e+05)
mc.add("v12_qcd_py8_pt170to230",    '%s/%s/QCD_Pt-170to230_7TeV-pythia8.%s/")'  %(srm,py8Dir,py8Gt2), xs = 1.994e+04)
mc.add("v12_qcd_py8_pt230to300",    '%s/%s/QCD_Pt-230to300_7TeV-pythia8.%s/")'  %(srm,py8Dir,py8Gt2), xs = 4.123e+03)
mc.add("v12_qcd_py8_pt300to380",    '%s/%s/QCD_Pt-300to380_7TeV-pythia8.%s/")'  %(srm,py8Dir,py8Gt1), xs = 9.593e+02)
mc.add("v12_qcd_py8_pt380to470",    '%s/%s/QCD_Pt-380to470_7TeV-pythia8.%s/")'  %(srm,py8Dir,py8Gt1), xs = 2.434e+02)
mc.add("v12_qcd_py8_pt470to600",    '%s/%s/QCD_Pt-470to600_7TeV-pythia8.%s/")'  %(srm,py8Dir,py8Gt1), xs = 7.410e+01)
mc.add("v12_qcd_py8_pt600to800",    '%s/%s/QCD_Pt-600to800_7TeV-pythia8.%s/")'  %(srm,py8Dir,py8Gt1), xs = 1.657e+01)
mc.add("v12_qcd_py8_pt800to1000",   '%s/%s/QCD_Pt-800to1000_7TeV-pythia8.%s/")' %(srm,py8Dir,py8Gt2), xs = 1.997e+00)
mc.add("v12_qcd_py8_pt1000to1400",  '%s/%s/QCD_Pt-1000to1400_7TeV-pythia8.%s/")'%(srm,py8Dir,py8Gt2), xs = 3.621e-01)
mc.add("v12_qcd_py8_pt1400to1800",  '%s/%s/QCD_Pt-1400to1800_7TeV-pythia8.%s/")'%(srm,py8Dir,py8Gt2), xs = 1.179e-02)
mc.add("v12_qcd_py8_pt1800to2200",  '%s/%s/QCD_Pt-1800to2200_7TeV-pythia8.%s/")'%(srm,py8Dir,py8Gt2), xs = 3.743e-04)
mc.add("v12_qcd_py8_pt2200to2600",  '%s/%s/QCD_Pt-2200to2600_7TeV-pythia8.%s/")'%(srm,py8Dir,py8Gt2), xs = 7.590e-06)
mc.add("v12_qcd_py8_pt2600to3000",  '%s/%s/QCD_Pt-2600to3000_7TeV-pythia8.%s/")'%(srm,py8Dir,py8Gt2), xs = 5.458e-08)
mc.add("v12_qcd_py8_pt3000to3500",  '%s/%s/QCD_Pt-3000to3500_7TeV-pythia8.%s/")'%(srm,py8Dir,py8Gt2), xs = 3.283e-11)

mgKFactor = 3048.0/2400.0 #Z+jets NNLO/LO

#MG QCD
mgDir = "/as1604/ICF/automated/2010_08_15_23_34_33/"
mc.add("v12_qcd_mg_ht_50_100",  '%s/%s/QCD_Pt-50To100_7TeV-madgraph.Spring10-START3X_V26-v1.GEN-SIM-RECO/")'%(srm,mgDir),
       xs = {"LO":30e+06, "fakeNLO":30e+06*mgKFactor}["fakeNLO"])
mc.add("v12_qcd_mg_ht_100_250", '%s/%s/QCD_Pt100to250-madgraph.Spring10-START3X_V26_S09-v2.GEN-SIM-RECO/ ")'%(srm,mgDir),
       xs = {"LO":7e+06, "fakeNLO":7e+06*mgKFactor}["fakeNLO"])
mc.add("v12_qcd_mg_ht_250_500", '%s/%s/QCD_Pt250to500-madgraph.Spring10-START3X_V26_S09-v1.GEN-SIM-RECO/ ")'%(srm,mgDir),
       xs = {"LO":171e+03, "fakeNLO":171e+03*mgKFactor}["fakeNLO"])
mc.add("v12_qcd_mg_ht_500_1000",'%s/%s/QCD_Pt500to1000-madgraph.Spring10-START3X_V26_S09-v1.GEN-SIM-RECO/")'%(srm,mgDir),
       xs = {"LO":5200, "fakeNLO":5200*mgKFactor}["fakeNLO"])
mc.add("v12_qcd_mg_ht_1000_inf",'%s/%s/QCD_Pt1000toInf-madgraph.Spring10-START3X_V26_S09-v1.GEN-SIM-RECO/")'%(srm,mgDir),
       xs = {"LO":83, "fakeNLO":83*mgKFactor}["fakeNLO"])

#ALPGEN QCD
#efficiences from https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/851.html
agDir = "/elaird/ICF/automated/2010_08_26_19_33_35/"
mc.add("v12_qcd_ag_2jets_pt40to120"  ,'%s/%s/QCD2Jets_Pt40to120-alpgen.Summer10-START36_V9_S09-v1.GEN-SIM-RECO/")'  %(srm,agDir), xs = 0.571 * 1.67E+07 )
mc.add("v12_qcd_ag_2jets_pt120to280" ,'%s/%s/QCD2Jets_Pt120to280-alpgen.Summer10-START36_V9_S09-v1.GEN-SIM-RECO/")' %(srm,agDir), xs = 0.279 * 1.23E+05 )
mc.add("v12_qcd_ag_2jets_pt280to500" ,'%s/%s/QCD2Jets_Pt280to500-alpgen.Summer10-START36_V9_S09-v1.GEN-SIM-RECO/")' %(srm,agDir), xs = 0.204 * 1.57E+03 )
mc.add("v12_qcd_ag_2jets_pt500to5000",'%s/%s/QCD2Jets_Pt500to5000-alpgen.Summer10-START36_V9_S09-v1.GEN-SIM-RECO/")'%(srm,agDir), xs = 0.192 * 5.23E+01 )

mc.add("v12_qcd_ag_3jets_pt40to120"  ,'%s/%s/QCD3Jets_Pt40to120-alpgen.Summer10-START36_V9_S09-v1.GEN-SIM-RECO/")'  %(srm,agDir), xs = 0.299 * 1.13E+07 )
mc.add("v12_qcd_ag_3jets_pt120to280" ,'%s/%s/QCD3Jets_Pt120to280-alpgen.Summer10-START36_V9_S09-v1.GEN-SIM-RECO/")' %(srm,agDir), xs = 0.207 * 2.22E+05 )
mc.add("v12_qcd_ag_3jets_pt280to500" ,'%s/%s/QCD3Jets_Pt280to500-alpgen.Summer10-START36_V9_S09-v1.GEN-SIM-RECO/")' %(srm,agDir), xs = 0.144 * 3.52E+03 )
mc.add("v12_qcd_ag_3jets_pt500to5000",'%s/%s/QCD3Jets_Pt500to5000-alpgen.Summer10-START36_V9_S09-v1.GEN-SIM-RECO/")'%(srm,agDir), xs = 0.134 * 1.23E+02 )

mc.add("v12_qcd_ag_4jets_pt40to120"  ,'%s/%s/QCD4Jets_Pt40to120-alpgen.Summer10-START36_V9_S09-v1.GEN-SIM-RECO/")'  %(srm,agDir), xs = 0.184 * 2.34E+06 )
mc.add("v12_qcd_ag_4jets_pt120to280" ,'%s/%s/QCD4Jets_Pt120to280-alpgen.Summer10-START36_V9_S09-v1.GEN-SIM-RECO/")' %(srm,agDir), xs = 0.147 * 1.47E+05 )
mc.add("v12_qcd_ag_4jets_pt280to500" ,'%s/%s/QCD4Jets_Pt280to500-alpgen.Summer10-START36_V9_S09-v1.GEN-SIM-RECO/")' %(srm,agDir), xs = 0.106 * 3.30E+03 )
mc.add("v12_qcd_ag_4jets_pt500to5000",'%s/%s/QCD4Jets_Pt500to5000-alpgen.Summer10-START36_V9_S09-v1.GEN-SIM-RECO/")'%(srm,agDir), xs = 0.103 * 1.27E+02 )

mc.add("v12_qcd_ag_5jets_pt40to120"  ,'%s/%s/QCD5Jets_Pt40to120-alpgen.Summer10-START36_V9_S09-v1.GEN-SIM-RECO/")'  %(srm,agDir), xs = 0.115 * 4.71E+05 )
mc.add("v12_qcd_ag_5jets_pt120to280" ,'%s/%s/QCD5Jets_Pt120to280-alpgen.Summer10-START36_V9_S09-v1.GEN-SIM-RECO/")' %(srm,agDir), xs = 0.101 * 6.30E+04 )
mc.add("v12_qcd_ag_5jets_pt280to500" ,'%s/%s/QCD5Jets_Pt280to500-alpgen.Summer10-START36_V9_S09-v1.GEN-SIM-RECO/")' %(srm,agDir), xs = 0.075 * 2.01E+03 )
mc.add("v12_qcd_ag_5jets_pt500to5000",'%s/%s/QCD5Jets_Pt500to5000-alpgen.Summer10-START36_V9_S09-v1.GEN-SIM-RECO/")'%(srm,agDir), xs = 0.076 * 8.49E+01 )

mc.add("v12_qcd_ag_6jets_pt40to120"  ,'%s/%s/QCD6Jets_Pt40to120-alpgen.Summer10-START36_V9_S09-v2.GEN-SIM-RECO/")'  %(srm,agDir), xs = 0.095 * 8.33E+04 )
mc.add("v12_qcd_ag_6jets_pt120to280" ,'%s/%s/QCD6Jets_Pt120to280-alpgen.Summer10-START36_V9_S09-v2.GEN-SIM-RECO/")' %(srm,agDir), xs = 0.102 * 2.02E+04 )
mc.add("v12_qcd_ag_6jets_pt280to500" ,'%s/%s/QCD6Jets_Pt280to500-alpgen.Summer10-START36_V9_S09-v2.GEN-SIM-RECO/")' %(srm,agDir), xs = 0.099 * 9.05E+02 )
mc.add("v12_qcd_ag_6jets_pt500to5000",'%s/%s/QCD6Jets_Pt500to5000-alpgen.Summer10-START36_V9_S09-v2.GEN-SIM-RECO/")'%(srm,agDir), xs = 0.099 * 4.20E+01 )

#MG GAMMA + JETS
mc.add("v12_g_jets_mg_pt40_100", '%s/arlogb//ICF/automated/2010_07_26_15_14_40//PhotonJets_Pt40to100-madgraph.Spring10-START3X_V26_S09-v1.GEN-SIM-RECO/")'%srm,
       xs = {"LO":23620, "fakeNLO":23620*mgKFactor}["fakeNLO"])
mc.add("v12_g_jets_mg_pt100_200",'%s/arlogb/ICF/automated/2010_07_26_15_14_40/PhotonJets_Pt100to200-madgraph.Spring10-START3X_V26_S09-v1.GEN-SIM-RECO/")'%srm,
       xs = {"LO":3476, "fakeNLO":3476*mgKFactor}["fakeNLO"])
mc.add("v12_g_jets_mg_pt200",    '%s/arlogb/ICF/automated/2010_07_26_15_14_40/PhotonJets_Pt200toInf-madgraph.Spring10-START3X_V26_S09-v1.GEN-SIM-RECO/")'%srm,
       xs = {"LO":485, "fakeNLO":485*mgKFactor}["fakeNLO"])

#PY GAMMA + JETS INCLUSIVE
rnString = 'utils.fileListFromDisk(location = "/vols/cms02/rnandi/PhotonJet_Pt%d.root", isDirectory = False)'
mc.add("v12_g_jets_py6_pt15"  , rnString%15  , xs = 1.922e+05, ptHatMin =   15) 
mc.add("v12_g_jets_py6_pt30"  , rnString%30  , xs = 2.007e+04, ptHatMin =   30) 
mc.add("v12_g_jets_py6_pt80"  , rnString%80  , xs = 5.565e+02, ptHatMin =   80) 
mc.add("v12_g_jets_py6_pt170" , rnString%170 , xs = 2.437e+01, ptHatMin =  170)
mc.add("v12_g_jets_py6_pt300" , rnString%300 , xs = 1.636e+00, ptHatMin =  300)
mc.add("v12_g_jets_py6_pt470" , rnString%470 , xs = 1.360e-01, ptHatMin =  470)
mc.add("v12_g_jets_py6_pt800" , rnString%800 , xs = 3.477e-03, ptHatMin =  800)
mc.add("v12_g_jets_py6_pt1400", rnString%1400, xs = 1.286e-05, ptHatMin = 1400)
mc.add("v12_g_jets_py6_pt2200", rnString%2200, xs = 4.035e-09, ptHatMin = 2200)
mc.add("v12_g_jets_py6_pt3000", rnString%3000, xs = 1.779e-14, ptHatMin = 3000)
mc.adjustOverlappingSamples( ["v12_g_jets_py6_pt%d"%i for i in [15,30,80,170,300,470,800,1400,2200,3000] ] )

#MG TT/EWK
mc.add("tt_tauola_mg_v12",'utils.fileListFromDisk(location = "/vols/cms01/mstoye/ttTauola_madgraph_V11tag/SusyCAF_Tree*.root", isDirectory = False, itemsToSkip = ["_4_2"] )',
       xs = {"LO":95.0,"NLO":157.5}["NLO"] )
mc.add("z_inv_mg_v12",'%s/zph04/ICF/automated/2010_07_14_11_52_58/",itemsToSkip=["14_3.root"])'%srm,
       xs = {"LO":4500.0,"fakeNLO":4500.0*mgKFactor}["fakeNLO"] )
mc.add("z_jets_mg_v12",'%s/jad/ICF/automated//2010_07_05_22_43_20/", pruneList=False)'%srm,
       xs = {"LO":2400.0,"NNLO":3048.0}["NNLO"] )
mc.add("w_jets_mg_v12",'%s/jad/ICF/automated//2010_06_18_22_33_23/")'%srm,
       xs = {"LO":24170.0,"NNLO":31314.0}["NNLO"] )
mc.add("zg_nunu", '%s/bbetchar/ICF/automated/2010_10_09_21_40_58/")'%srm, xs = 2.68 )

mc.add("tt_tauola_mg",'%s/bbetchar//ICF/automated/2010_11_17_05_40_17/")'%srm, xs = {"LO":121, "BurtGuessNLO":157.5}["BurtGuessNLO"])

#PYTHIA
mc.add("w_enu",       '%s/bbetchar//ICF/automated/2010_10_11_21_21_39/")'%srm, xs = {"LO":7899,"BurtGuessNNLO":10234}["BurtGuessNNLO"]) 
mc.add("w_munu",      '%s/trommers/ICF/automated/2010_10_06_18_40_51/")'%srm,  xs = {"LO":7899,"BurtGuessNNLO":10234}["BurtGuessNNLO"]) 
mc.add("w_taunu",     '%s/bbetchar//ICF/automated/2010_10_11_21_21_39/")'%srm, xs = {"LO":7899,"BurtGuessNNLO":10234}["BurtGuessNNLO"])

#SUSY LM
mc.add("lm0_v12" ,'%s/bainbrid/ICF/automated/2010_07_16_12_54_00/LM0.Spring10-START3X_V26_S09-v1.GEN-SIM-RECO/")'%srm,  xs = 38.93   )
mc.add("lm1_v12" ,'%s/bainbrid/ICF/automated/2010_07_12_17_52_54/LM1.Spring10-START3X_V26_S09-v1.GEN-SIM-RECO/")'%srm,  xs = 4.888   )
mc.add("lm2_v12" ,'%s/bainbrid/ICF/automated/2010_07_16_12_54_00/LM2.Spring10-START3X_V26_S09-v1.GEN-SIM-RECO/")'%srm,  xs = 0.6027  )
mc.add("lm3_v12" ,'%s/bainbrid/ICF/automated/2010_07_12_17_52_54/LM3.Spring10-START3X_V26_S09-v1.GEN-SIM-RECO/")'%srm,  xs = 3.438   )
mc.add("lm4_v12" ,'%s/bainbrid/ICF/automated/2010_07_16_12_54_00/LM4.Spring10-START3X_V26_S09-v1.GEN-SIM-RECO/")'%srm,  xs = 1.879   )
mc.add("lm5_v12" ,'%s/bainbrid/ICF/automated/2010_07_16_12_54_00/LM5.Spring10-START3X_V26_S09-v1.GEN-SIM-RECO/")'%srm,  xs = 0.4734  )
mc.add("lm6_v12" ,'%s/bainbrid/ICF/automated/2010_07_16_12_54_00/LM6.Spring10-START3X_V26_S09-v1.GEN-SIM-RECO/")'%srm,  xs = 0.3104  )
mc.add("lm7_v12" ,'%s/bainbrid/ICF/automated/2010_07_16_12_54_00/LM7.Spring10-START3X_V26_S09-v1.GEN-SIM-RECO/")'%srm,  xs = 1.209   )
mc.add("lm8_v12" ,'%s/bainbrid/ICF/automated/2010_07_16_12_54_00/LM8.Spring10-START3X_V26_S09-v1.GEN-SIM-RECO/")'%srm,  xs = 0.7300  )
mc.add("lm9_v12" ,'%s/bainbrid/ICF/automated/2010_07_12_17_52_54/LM9.Spring10-START3X_V26_S09-v1.GEN-SIM-RECO/")'%srm,  xs = 7.134   )
mc.add("lm10_v12",'%s/bainbrid/ICF/automated/2010_07_12_17_52_54/LM10.Spring10-START3X_V26_S09-v1.GEN-SIM-RECO/")'%srm, xs = 0.04778 )
mc.add("lm11_v12",'%s/bainbrid/ICF/automated/2010_07_16_12_54_00/LM11.Spring10-START3X_V26_S09-v1.GEN-SIM-RECO/")'%srm, xs = 0.8236  )
mc.add("lm12_v12",'%s/bainbrid/ICF/automated/2010_07_12_17_52_54/LM12.Spring10-START3X_V26_S09-v1.GEN-SIM-RECO/")'%srm, xs = 4.414   )
mc.add("lm13_v12",'%s/bainbrid/ICF/automated/2010_07_16_12_54_00/LM13.Spring10-START3X_V26_S09-v1.GEN-SIM-RECO/")'%srm, xs = 6.899   )

#SMS
mc.add("t1", '%s/bbetchar/ICF/automated/2011_01_24_18_09_42/")'%srm, xs = 1.0) #dummy xs
mc.add("t2", '%s/bbetchar/ICF/automated/2011_01_24_18_14_21/")'%srm, xs = 1.0) #dummy xs
mc.add("t3", '%s/bbetchar/ICF/automated/2011_01_25_21_33_22/")'%srm, xs = 1.0) #dummy xs
mc.add("t4", '%s/henning//ICF/automated/2011_01_27_12_39_13/")'%srm, xs = 1.0) #dummy xs

#MSUGRA SCANS
mc.add("scan_tanbeta3_tanja1",'%s/trommers//ICF/automated/2010_12_03_15_28_38/")'%srm, xs = 100000.0 ) #dummy xs
mc.add("scan_tanbeta3_tanja2",'%s/trommers//ICF/automated/2010_12_17_20_09_06/")'%srm, xs = 100000.0 ) #dummy xs
dead = str(["_43_1_vfJ", "_109_1_K7V", "_648_1_SNq", "_19_1_P2A", '_106_1_5y9', '_107_1_xWu', '_116_1_GYJ', '_157_1_P2J', '_163_1_zG0', '_171_1_mkQ', '_173_1_gew', '_177_1_8Ew', '_191_1_Kdf', '_226_1_21E', '_266_1_ygR', '_295_1_BNB', '_311_1_Thz', '_319_1_UoK', '_328_2_Pn2', '_333_1_6tK', '_369_1_NG4', '_378_1_Nhp', '_413_1_IrU', '_42_1_h0i', '_446_1_5vP', '_449_1_a7D', '_487_1_VjZ', '_509_1_SiQ', '_547_1_hG9', '_559_1_Jjs', '_579_1_GaN', '_599_1_7Fs', '_644_1_1Zq', '_67_1_glp']).replace("'",'"')
mc.add("scan_tanbeta3_burt1", '%s/bbetchar//ICF/automated/2010_11_10_19_34_17/", itemsToSkip = %s)'%(srm, dead), xs = 100000.0 ) #dummy xs
mc.add("scan_tanbeta3_skim200", 'utils.fileListFromDisk(location = "/vols/cms02/elaird1/11_skims/33_m12_200/")', xs = 1.0 ) #dummy xs
mc.add("scan_tanbeta3_skim100", 'utils.fileListFromDisk(location = "/vols/cms02/elaird1/11_skims/34_m12_100/")', xs = 1.0 ) #dummy xs

#MG EWK SKIMS
mc.add("w_jets_mg_v12_skim", 'utils.fileListFromDisk(location = "/vols/cms02/elaird1/14_skims/WJets-madgraph.Spring10-START3X_V26_S09-v1.GEN-SIM-RECO/")',
       xs = 4.462162e-03 * 3.131400e+04)
mc.add("z_inv_mg_v12_skim",  'utils.fileListFromDisk(location = "/vols/cms02/elaird1/14_skims/ZinvisibleJets-madgraph.Spring10-START3X_V26_S09-v1.GEN-SIM-RECO/")',
       xs = 2.067155e-03 * 5.715000e+03)
mc.add("z_jets_mg_v12_skim", 'utils.fileListFromDisk(location = "/vols/cms02/elaird1/14_skims/ZJets-madgraph.Spring10-START3X_V26_S09-v1.GEN-SIM-RECO/", itemsToSkip = %s)'%str(["v12_1%d_skim.root"%i for i in range(1,10) ]), xs = 1.185708e-02 * 3.048000e+03)

#TEST SKIMS
mc.add("v12_g_jets_test_skim",
       'utils.fileListFromDisk(location = "/vols/cms02/elaird1/tmp//photonLook/pythia6_ge2_caloAK5_photonLoose/config00/g_jets_py6_pt*skim.root", isDirectory = False)',
       xs = 1.0)
mc.add("v12_qcd_test_skim",
       'utils.fileListFromDisk(location = "/vols/cms02/elaird1/tmp//photonLook/pythia6_ge2_caloAK5_photonLoose/config00/qcd_py6_pt*skim.root", isDirectory = False)',
       xs = 1.0)
mc.add("markus1",'["dcap://gfe02.grid.hep.ph.ic.ac.uk:22128/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/arlogb/ICF/automated/2010_07_14_17_25_26//PhotonJets_Pt200toInf-madgraph.Spring10-START3X_V26_S09-v1.GEN-SIM-RECO//SusyCAF_Tree_4_2.root"]', xs = 1.0)
