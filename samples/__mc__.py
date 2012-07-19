from supy.samples import SampleHolder
from supy.sites import srm,pnfs
mc = SampleHolder()
pnfs = pnfs()
mgKFactor = 3048.0/2400.0 #Z+jets NNLO/LO
mgKFactor2 = 3048.0/2475.0#Z+jets NNLO/LO

######## QCD ##########
mgQcdLoc = '/bbetchar//ICF/automated/2011_04_07_20_30_16/'
mgQcdDset = "QCD_TuneD6T_HT-%s_7TeV-madgraph.Spring11-PU_S1_START311_V1G1-v1.AODSIM"

mc.add("qcd_mg_ht_100_250",  '%s/%s/%s")'%(srm,mgQcdLoc,mgQcdDset%"100To250"), xs = {"LO":7e+06, "fakeNLO":7e+06*mgKFactor}["fakeNLO"])
mc.add("qcd_mg_ht_250_500",  '%s/%s/%s")'%(srm,mgQcdLoc,mgQcdDset%"250To500"), xs = {"LO":171e+03, "fakeNLO":171e+03*mgKFactor}["fakeNLO"])
mc.add("qcd_mg_ht_500_1000", '%s/%s/%s")'%(srm,mgQcdLoc,mgQcdDset%"500To1000"),xs = {"LO":5200, "fakeNLO":5200*mgKFactor}["fakeNLO"])
mc.add("qcd_mg_ht_1000_inf", '%s/%s/%s")'%(srm,mgQcdLoc,mgQcdDset%"1000"),     xs = {"LO":83, "fakeNLO":83*mgKFactor}["fakeNLO"])

bins,xss = zip(*tuple(
    [(  15, 8.159e+08),
     (  30, 5.312e+07),
     (  50, 6.359e+06),
     (  80, 7.843e+05),
     ( 120, 1.151e+05),
     ( 170, 2.426e+04),
     ( 300, 1.168e+03),
     ( 470, 7.022e+01),
     ( 600, 1.555e+01),
     ( 800, 1.844e+00),
     (1000, 3.321e-01),
     (1400, 1.087e-02),
     (1800, 3.575e-04),
     (None,None)]))

py6Loc = '/bbetchar//ICF/automated/2011_04_07_19_50_45/'
py6Dset = "/QCD_Pt_%s_TuneZ2_7TeV_pythia6.Spring11-PU_S1_START311_V1G1-v1.AODSIM"
for low,high,xs in zip(bins[:-1],bins[1:],xss) :
    mc.add("qcd_py6_pt_%s"%("%d_%d"%(low,high) if high else "%d"%low),
           '%s/%s/%s")'%(srm,py6Loc,py6Dset%("%dto%d"%(low,high) if high else "%d"%low)),
           xs = xs)

py6FJLoc = '/bbetchar//ICF/automated/2011_06_12_05_20_30/'
py6Summer11Dset = "/QCD_Pt-%s_TuneZ2_7TeV_pythia6.Summer11-PU_S3_START42_V11-v2.AODSIM"
for low,high,xs in zip(bins[:-1],bins[1:],xss[:-1]) :
    mc.add("qcd_py6fj_pt_%s"%("%d_%d"%(low,high) if high else str(low)),
           '%s/%s/%s")'%(srm,py6FJLoc,py6Summer11Dset%("%dto%d"%(low,high) if high else "%d"%low)),
           xs = xs)

#### QCD mu enriched ####
bins,xss,forms = zip(*tuple(
    [(15,  5.792e+08 * 0.00254, 0), # xs * filter efficiency
     (20,  2.363e+08 * 0.00518, 0),
     (30,  5.307e+07 * 0.01090, 0),
     (50,  6.351e+06 * 0.02274, 1),
     (80,  7.851e+05 * 0.03700, 1),
     (120, 9.295e+04 * 0.04777, 1),
     (150, 4.758e+04 * 0.05964, 1),
     (None,None,None)]))
py6FJmuLoc = '/bbetchar//ICF/automated/2011_07_15_16_44_18/'

formats = ["QCD_Pt-%s_MuPt5Enriched_TuneZ2_7TeV-pythia6.Summer11-PU_S3_START42_V11-v2.AODSIM/",
           "QCD_Pt-%s_MuPt5Enriched_TuneZ2_7TeV-pythia6.Summer11-PU_S4_START42_V11-v1.AODSIM/"]
for low,high,xs,form in zip(bins[:-1],bins[1:],xss[:-1],forms[:-1]) :
    mc.add("qcd_py6fjmu_pt_%s"%("%d_%d"%(low,high) if high else str(low)),
           '%s/%s/%s")'%(srm,py6FJmuLoc,formats[form]%("%dto%d"%(low,high) if high else str(low))),
           xs = xs)

######### TT / EWK ############

pyFJLoc = '/bbetchar//ICF/automated/2011_07_21_00_58_10/'
summer11pu = '.Summer11-PU_S3_START42_V11-v2.AODSIM/'
mc.add("tt_tauola_fj", '%s/arlogb//ICF/automated/2011_07_11_17_17_07/")'%srm,    xs = {"LO":94,    "BurtGuessNLO":157.5}["BurtGuessNLO"])   
mc.add("w_enu_fj",     '%s/%s/WToENu_TuneZ2_7TeV-pythia6%s")'%(srm,pyFJLoc,summer11pu),       xs = {"LO": 7899, "BurtGuessNNLO": 15639}["BurtGuessNNLO"])
mc.add("w_munu_fj",    '%s/%s/WToMuNu_TuneZ2_7TeV-pythia6%s")'%(srm,pyFJLoc,summer11pu),       xs = {"LO": 7899, "BurtGuessNNLO": 15639}["BurtGuessNNLO"])
mc.add("w_taunu_fj",   '%s/%s/WToTauNu_TuneZ2_7TeV-pythia6-tauola%s")'%(srm,pyFJLoc,summer11pu),xs = {"LO": 7899, "BurtGuessNNLO": 15639}["BurtGuessNNLO"])

# https://twiki.cern.ch/twiki/bin/view/CMS/MadGraphSummer11Production
mc.add("tt_tauola_fj_mg",'%s/bbetchar//ICF/automated/2011_07_20_22_27_52/")'%srm, xs = 319.18)
mc.add("w_jets_fj_mg", '%s/gouskos//ICF/automated/2011_07_18_17_43_04/")'%srm, xs = 55854)

burt_ttbar = '/bbetchar//ICF/automated/2011_04_07_19_30_01/'
burt_ewk = '/bbetchar//ICF/automated/2011_04_07_19_40_51/'
spring11pu = "Spring11-PU_S1_START311_V1G1-v1.AODSIM"

#MG Spring '11
mc.add("tt_tauola_mg",'%s%s/TTJets_TuneZ2_7TeV-madgraph-tauola.%s",itemsToSkip = ["SusyCAF_Tree_49_1_0fC.root","SusyCAF_Tree_24_1_cuB.root"])'%(srm,burt_ttbar,spring11pu), xs = {"LO":121, "BurtGuessNLO":157.5}["BurtGuessNLO"])
mc.add("w_jets_mg",'%s%s/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola.%s")'%(srm,burt_ewk,spring11pu), xs = {"LO":24640, "BurtGuessNNLO": 31924}["BurtGuessNNLO"])
mc.add("zinv_jets_mg",'%s/henning//ICF/automated/2011_04_13_12_16_30/")'%srm, xs = {"LO":4500.0,"fakeNLO":4500.0*mgKFactor}["fakeNLO"])

mc.add("dyll_jets_mg", '%s/gouskos//ICF/automated/2011_05_02_12_39_53/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola.Spring11-PU_S1_START311_V1G1-v1.AODSIM ")'%srm,
       xs = {"LO":2321.0, "NNLO":3048}["NNLO"])

mc.add("zinv_jets_mg_spring11_skim", 'utils.fileListFromDisk(location = "/vols/cms02/elaird1/29_skims/07_zinv/v1/zinv_jets_mg_*_skim.root", isDirectory = False)', xs = 2.353024e-03 * 5.715000e+03)

#PYTHIA
mc.add("tt_tauola_pythia",'%s%s/TT_TuneZ2_7TeV-pythia6-tauola.%s")'%(srm,burt_ttbar,spring11pu), xs = {"LO":94, "BurtGuessNLO":122}["BurtGuessNLO"])
mc.add("w_enu", '%s%s/WToENu_TuneZ2_7TeV-pythia6.%s)'%(srm,burt_ewk,spring11pu),            xs = {"LO": 7899, "BurtGuessNNLO": 10234}["BurtGuessNNLO"])
mc.add("w_munu", '%s%s/WToMuNu_TuneZ2_7TeV-pythia6.%s)'%(srm,burt_ewk,spring11pu),          xs = {"LO": 7899, "BurtGuessNNLO": 10234}["BurtGuessNNLO"])
mc.add("w_taunu", '%s%s/WToTauNu_TuneZ2_7TeV-pythia6-tauola.%s)'%(srm,burt_ewk,spring11pu), xs = {"LO": 7899, "BurtGuessNNLO": 10234}["BurtGuessNNLO"])
mc.add("z_nunu",  '%s/henning//ICF/automated/2011_04_15_10_55_57/")'%srm, xs = 4292)


## MG Summer '11 samples (CMSSW_4_2_5 / V15-03-14 / L1FastJetL2L3)
mc.add("dyll_jets_mg_summer11", '%s/mstoye//ICF/automated/2011_08_06_17_06_28/")'%srm, xs = {"LO":2475.0, "NNLO":3048}["NNLO"])
mc.add("tt_jets_mg_tauola_summer11", '%s/mstoye//ICF/automated/2011_07_30_02_06_09/", alwaysUseLastAttempt = True)'%srm, xs = {"LO":94.76, "guessNLO":157.5}["guessNLO"])

l = 'utils.fileListFromDisk(isDirectory = False, location = "/vols/cms02/elaird1/29_skims/08_mumu/v2'
mc.add("dyll_jets_mg_summer11_mumuSkim",      '%s/dyll_jets_mg_summer11_*_skim.root")'     %l, xs = 1.581657e-01 * 3.048000e+03)
mc.add("tt_jets_mg_tauola_summer11_mumuSkim", '%s/tt_jets_mg_tauola_summer11_*_skim.root")'%l, xs = 1.276058e-02 * 1.575000e+02)

l = 'utils.fileListFromDisk(isDirectory = False, location = "/vols/cms02/elaird1/29_skims/08_mumu/v3'
mc.add("dyll_jets_mg_summer11_mumuHtSkim",      '%s/dyll_jets_mg_summer11_mumuSkim_*_skim.root")'     %l, xs = 6.685125e-03 * 4.820891e+02)
mc.add("tt_jets_mg_tauola_summer11_mumuHtSkim", '%s/tt_jets_mg_tauola_summer11_mumuSkim_*_skim.root")'%l, xs = 2.900146e-01 * 2.009791e+00)

# Z->nunu (HT binned)
zNunuLoc = "/dburton//ICF/automated/2011_10_10_13_33_29/"
tag = "7TeV-madgraph.Summer11-PU_S4_START42_V11-v1.AODSIM"
mc.add("znunu_jets_mg_ht_50_100_summer11",  '%s/%s/ZJetsToNuNu_50_HT_100_%s_")'%(srm, zNunuLoc, tag), xs = {"LO": 309.5, "fakeNLO": 309.5*mgKFactor}["fakeNLO"])
mc.add("znunu_jets_mg_ht_100_200_summer11", '%s/%s/ZJetsToNuNu_100_HT_200_%s")'%(srm, zNunuLoc, tag), xs = {"LO": 125.2, "fakeNLO": 125.2*mgKFactor}["fakeNLO"])
mc.add("znunu_jets_mg_ht_200_inf_summer11", '%s/%s/ZJetsToNuNu_200_HT_inf_%s")'%(srm, zNunuLoc, tag), xs = {"LO":  32.9, "fakeNLO":  32.9*mgKFactor}["fakeNLO"])

dir = "/vols/cms02/elaird1/29_skims/07_zinv/v2"
l = 'utils.fileListFromDisk(isDirectory = False, location = '
mc.add("znunu_jets_mg_ht_50_100_summer11_skim",   '%s"%s/znunu_jets_mg_ht_50_100_*_skim.root")'%(l, dir), xs = 1.288242e-03 * 309.5*mgKFactor)
mc.add("znunu_jets_mg_ht_100_200_summer11_skim", '%s"%s/znunu_jets_mg_ht_100_200_*_skim.root")'%(l, dir), xs = 9.922905e-03 * 125.2*mgKFactor)
mc.add("znunu_jets_mg_ht_200_inf_summer11_skim", '%s"%s/znunu_jets_mg_ht_200_inf_*_skim.root")'%(l, dir), xs = 5.364460e-01 *  32.9*mgKFactor)

# W->lnu (HT binned)
wFactor = 31314./27770.0
mc.add("w_jets_mg_tauola_ht_300_inf_summer11", '%s/gouskos//ICF/automated/2011_08_19_18_17_37/")'%srm, xs = {"LO":48.49, "guessNLO":48.49*wFactor}["guessNLO"])

# G + jets (HT binned)
gJetMg = "/elaird//ICF/automated/2011_10_14_14_40_23/"
tag = "7TeV-madgraph.Summer11-PU_S4_START42_V11-v1.AODSIM"
mc.add("g_jets_mg_ht_40_100_summer11",   '%s/%s/GJets_TuneZ2_40_HT_100_%s")'%(srm, gJetMg, tag), xs = {"LO":25690,   "fakeNLO":25690  *mgKFactor}["fakeNLO"])
mc.add("g_jets_mg_ht_100_200_summer11", '%s/%s/GJets_TuneZ2_100_HT_200_%s")'%(srm, gJetMg, tag), xs = {"LO": 5213,   "fakeNLO": 5213  *mgKFactor}["fakeNLO"])
mc.add("g_jets_mg_ht_200_inf_summer11", '%s/%s/GJets_TuneZ2_200_HT_inf_%s")'%(srm, gJetMg, tag), xs = {"LO":  798.3, "fakeNLO":  798.3*mgKFactor}["fakeNLO"])

# G + jets (HT binned) skims
dir = "/vols/cms02/elaird1/29_skims/04_photons/v6"
l = 'utils.fileListFromDisk(isDirectory = False, location = '
mc.add("g_jets_mg_ht_40_100_summer11_skim",   '%s"%s/g_jets_mg_ht_40_100_summer11_*_skim.root")'%(l, dir), xs = 1.114472e-04 * 25690  *mgKFactor)
mc.add("g_jets_mg_ht_100_200_summer11_skim", '%s"%s/g_jets_mg_ht_100_200_summer11_*_skim.root")'%(l, dir), xs = 3.303680e-02 *  5213  *mgKFactor)
mc.add("g_jets_mg_ht_200_inf_summer11_skim", '%s"%s/g_jets_mg_ht_200_inf_summer11_*_skim.root")'%(l, dir), xs = 8.159796e-02 *   798.3*mgKFactor)

# QCD (HT binned)
qcdDset1 = "QCD_TuneZ2_HT-%s_7TeV-madgraph.Summer11-PU_S4_START42_V11-v1.AODSIM"
qcdDset3 = "QCD_TuneZ2_HT-%s_7TeV-madgraph.Summer11-PU_S4_START42_V11-v3.AODSIM"
qcdDir1 = "/dburton//ICF/automated/2011_10_26_12_39_58/"
qcdDir2 = "/dburton//ICF/automated/2011_11_24_13_30_10/"

mc.add("qcd_mg_ht_100_250_summer11",  '%s/%s/%s")'%(srm, qcdDir1, qcdDset1%"100To250"), xs = {"LO":4194000.0, "fakeNLO":4194000.0*mgKFactor}["fakeNLO"])
mc.add("qcd_mg_ht_250_500_summer11",  '%s/%s")'   %(srm, qcdDir2,                    ), xs = {"LO": 198500.0, "fakeNLO": 198500.0*mgKFactor}["fakeNLO"])
mc.add("qcd_mg_ht_500_1000_summer11", '%s/%s/%s")'%(srm, qcdDir1, qcdDset1%"500To1000"),xs = {"LO":   5856.0, "fakeNLO":   5856.0*mgKFactor}["fakeNLO"])
mc.add("qcd_mg_ht_1000_inf_summer11", '%s/%s/%s")'%(srm, qcdDir1, qcdDset1%"1000"),     xs = {"LO":    122.6, "fakeNLO":    122.6*mgKFactor}["fakeNLO"])

# QCD (HT binned) skims
dir = "/vols/cms02/elaird1/29_skims/04_photons/v8"
l = 'utils.fileListFromDisk(isDirectory = False, location = '
mc.add("qcd_mg_ht_100_250_summer11_skim",   '%s"%s/qcd_mg_ht_100_250_summer11_*_skim.root")'%(l, dir), xs = 3.180463e-06 * 5.326380e+06)
mc.add("qcd_mg_ht_250_500_summer11_skim",   '%s"%s/qcd_mg_ht_250_500_summer11_*_skim.root")'%(l, dir), xs = 2.491150e-04 * 2.520950e+05)
mc.add("qcd_mg_ht_500_1000_summer11_skim", '%s"%s/qcd_mg_ht_500_1000_summer11_*_skim.root")'%(l, dir), xs = 2.405978e-04 * 7.437120e+03)
mc.add("qcd_mg_ht_1000_inf_summer11_skim", '%s"%s/qcd_mg_ht_1000_inf_summer11_*_skim.root")'%(l, dir), xs = 1.675973e-04 * 1.557020e+02)

## MG  Spring '11 (L2L3)
#gJetMg = "/elaird//ICF/automated/2011_04_04_11_45_04/"
gJetMg = "/bm409//ICF/automated/2011_06_08_16_42_46/"
mc.add("g_jets_mg_ht_40_100_spring11",'%s/%s/GJets_TuneD6T_HT-40To100_7TeV-madgraph.Spring11-PU_S1_START311_V1G1-v1.AODSIM")'%(srm, gJetMg),
       xs = {"LO":23620, "fakeNLO":23620*mgKFactor}["fakeNLO"])
mc.add("g_jets_mg_ht_100_200_spring11",'%s/%s/GJets_TuneD6T_HT-100To200_7TeV-madgraph.Spring11-PU_S1_START311_V1G1-v1.AODSIM")'%(srm, gJetMg),
       xs = {"LO":3476, "fakeNLO":3476*mgKFactor}["fakeNLO"])
mc.add("g_jets_mg_ht_200_inf_spring11"   ,'%s/%s/GJets_TuneD6T_HT-200_7TeV-madgraph.Spring11-PU_S1_START311_V1G1-v1.AODSIM")'%(srm, gJetMg),
       xs = {"LO":485, "fakeNLO":485*mgKFactor}["fakeNLO"])

#MG (L1OffsetL2L3)
gJetMg = "/bm409//ICF/automated/2011_06_16_17_15_39/"
mc.add("g_jets_mg_ht_40_100.L1",'%s/%s/_.GJets_TuneD6T_HT-40To100_7TeV-madgraph.Spring11-PU_S1_START311_V1G1-v1.AODSIM")'%(srm, gJetMg),
       xs = {"LO":23620, "fakeNLO":23620*mgKFactor}["fakeNLO"])
mc.add("g_jets_mg_ht_100_200.L1",'%s/%s/GJets_TuneD6T_HT-100To200_7TeV-madgraph.Spring11-PU_S1_START311_V1G1-v1.AODSIM")'%(srm, gJetMg),
       xs = {"LO":3476, "fakeNLO":3476*mgKFactor}["fakeNLO"])
mc.add("g_jets_mg_ht_200_inf.L1"   ,'%s/%s/GJets_TuneD6T_HT-200_7TeV-madgraph.Spring11-PU_S1_START311_V1G1-v1.AODSIM")'%(srm, gJetMg),
       xs = {"LO":485, "fakeNLO":485*mgKFactor}["fakeNLO"])

#PYTHIA
GPyLoc = "/dburton//ICF/automated/2011_04_15_14_36_55/"
GPyDset = "G_Pt_%s_TuneZ2_7TeV_pythia6.Spring11-PU_S1_START311_V1G1-v1.AODSIM"
mc.add("g_jets_py6_pt_0_15",      '%s/%s/%s")'%(srm,GPyLoc,GPyDset%"0to15"),      xs = 8.420e+07)
mc.add("g_jets_py6_pt_15_30",     '%s/%s/%s")'%(srm,GPyLoc,GPyDset%"15to30"),     xs = 1.717e+05)
mc.add("g_jets_py6_pt_30_50",     '%s/%s/%s")'%(srm,GPyLoc,GPyDset%"30to50"),     xs = 1.669e+04)
mc.add("g_jets_py6_pt_50_80",     '%s/%s/%s")'%(srm,GPyLoc,GPyDset%"50to80"),     xs = 2.722e+03)
mc.add("g_jets_py6_pt_80_120",    '%s/%s/%s")'%(srm,GPyLoc,GPyDset%"80to120"),    xs = 4.472e+02)
mc.add("g_jets_py6_pt_120_170",   '%s/%s/%s")'%(srm,GPyLoc,GPyDset%"120to170"),   xs = 8.417e+01)
mc.add("g_jets_py6_pt_170_300",   '%s/%s/%s")'%(srm,GPyLoc,GPyDset%"170to300"),   xs = 2.264e+01)
mc.add("g_jets_py6_pt_300_470",   '%s/%s/%s")'%(srm,GPyLoc,GPyDset%"300to470"),   xs = 1.493e+00)
mc.add("g_jets_py6_pt_470_800",   '%s/%s/%s")'%(srm,GPyLoc,GPyDset%"470to800"),   xs = 1.323e-01)
mc.add("g_jets_py6_pt_800_1400",  '%s/%s/%s")'%(srm,GPyLoc,GPyDset%"800to1400"),  xs = 3.481e-03)
mc.add("g_jets_py6_pt_1400_1800", '%s/%s/%s")'%(srm,GPyLoc,GPyDset%"1400to1800"), xs = 1.270e-05)
mc.add("g_jets_py6_pt_1800_inf",  '%s/%s/%s")'%(srm,GPyLoc,GPyDset%"1800"),       xs = 2.936e-07)

##### SUSY LM  ##########
lmDir = "/henning//ICF/automated/2011_04_12_14_13_06/"
lmTag = "_SUSY_sftsht_7TeV-pythia6.Spring11-PU_S1_START311_V1G1-v1.AODSIM"
mc.add("lm0" ,'%s/%s/LM0%s")' %(srm, lmDir, lmTag), xs = 38.93  )
mc.add("lm1" ,'%s/%s/LM1%s")' %(srm, lmDir, lmTag), xs = 4.888  )
mc.add("lm2" ,'%s/%s/LM2%s")' %(srm, lmDir, lmTag), xs = 0.6027 )
mc.add("lm3" ,'%s/%s/LM3%s")' %(srm, lmDir, lmTag), xs = 3.438  )
mc.add("lm4" ,'%s/%s/LM4%s")' %(srm, lmDir, lmTag), xs = 1.879  )
mc.add("lm5" ,'%s/%s/LM5%s")' %(srm, lmDir, lmTag), xs = 0.4734 )
mc.add("lm6" ,'%s/%s/LM6%s")' %(srm, lmDir, lmTag), xs = 0.3104 )
mc.add("lm7" ,'%s/%s/LM7%s")' %(srm, lmDir, lmTag), xs = 1.209  )
mc.add("lm8" ,'%s/%s/LM8%s")' %(srm, lmDir, lmTag), xs = 0.7300 )
mc.add("lm9" ,'%s/%s/LM9%s")' %(srm, lmDir, lmTag), xs = 7.134  )
mc.add("lm11",'%s/%s/LM11%s")'%(srm, lmDir, lmTag), xs = 0.8236 )
mc.add("lm12",'%s/%s/LM12%s")'%(srm, lmDir, lmTag), xs = 4.414  )
mc.add("lm13",'%s/%s/LM13%s")'%(srm, lmDir, lmTag), xs = 6.899  )

#SMS
mc.add("t1.ted", '%s/elaird//ICF/automated/2011_03_09_13_08_06/")'%srm, xs = 1.0) #dummy xs
mc.add("t1.yos", '%s/yeshaq//ICF/automated/2012_07_02_17_24_20/")'%pnfs, xs = {"LO":107.5}["LO"]) #dummy xs
mc.add("t2tt.yos", '%s/yeshaq//ICF/automated/2012_07_10_16_39_39/")'%pnfs, xs = {"LO":107.5}["LO"]) #dummy xs
mc.add("t2bb.yos", '%s/yeshaq//ICF/automated/2012_07_10_16_50_17/")'%pnfs, xs = {"LO":107.5}["LO"]) #dummy xs
mc.add("t1tttt.chr", '%s/clucas//ICF/automated/2012_06_04_18_26_52/")'%pnfs, xs = 1.0) #dummy xs
mc.add("t2.ted", '%s/elaird/ICF/automated/2011_03_08_23_04_26/")'%srm, lumi = 1.0) #dummy xs

#MSUGRA SCANS
mc.add("scan_tanbeta10_burt1",'%s/bbetchar//ICF/automated/2011_07_10_04_13_24/",)'%srm, lumi = 1.0) #dummy lumi
mc.add("scan_tanbeta3_tanja1",'%s/trommers//ICF/automated/2010_12_03_15_28_38/")'%srm, xs = 100000.0 ) #dummy xs
mc.add("scan_tanbeta3_tanja2",'%s/trommers//ICF/automated/2010_12_17_20_09_06/")'%srm, xs = 100000.0 ) #dummy xs
mc.add("scan_tanbeta10",('%s/trommers//ICF/automated/2011_02_24_16_41_52/")'%srm).replace("user","user_old"), xs = 1.0) #fake xs
dead = str(["_43_1_vfJ", "_109_1_K7V", "_648_1_SNq", "_19_1_P2A", '_106_1_5y9', '_107_1_xWu', '_116_1_GYJ', '_157_1_P2J', '_163_1_zG0', '_171_1_mkQ', '_173_1_gew', '_177_1_8Ew', '_191_1_Kdf', '_226_1_21E', '_266_1_ygR', '_295_1_BNB', '_311_1_Thz', '_319_1_UoK', '_328_2_Pn2', '_333_1_6tK', '_369_1_NG4', '_378_1_Nhp', '_413_1_IrU', '_42_1_h0i', '_446_1_5vP', '_449_1_a7D', '_487_1_VjZ', '_509_1_SiQ', '_547_1_hG9', '_559_1_Jjs', '_579_1_GaN', '_599_1_7Fs', '_644_1_1Zq', '_67_1_glp']).replace("'",'"')
mc.add("scan_tanbeta3_burt1", '%s/bbetchar//ICF/automated/2010_11_10_19_34_17/", itemsToSkip = %s)'%(srm, dead), xs = 100000.0 ) #dummy xs
mc.add("scan_tanbeta3_skim200", 'utils.fileListFromDisk(location = "/vols/cms02/elaird1/11_skims/33_m12_200/")', xs = 1.0 ) #dummy xs
mc.add("scan_tanbeta3_skim100", 'utils.fileListFromDisk(location = "/vols/cms02/elaird1/11_skims/34_m12_100/")', xs = 1.0 ) #dummy xs
