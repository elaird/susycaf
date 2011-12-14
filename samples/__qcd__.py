from supy.samples import SampleHolder
from supy.sites import srm
qcd = SampleHolder()

mgKFactor = 3048.0/2400.0 #Z+jets NNLO/LO
mgKFactor2 = 3048.0/2475.0#Z+jets NNLO/LO

mgQcdLoc = '/bbetchar//ICF/automated/2011_04_07_20_30_16/'
mgQcdDset = "QCD_TuneD6T_HT-%s_7TeV-madgraph.Spring11-PU_S1_START311_V1G1-v1.AODSIM"

qcd.add("qcd_mg_ht_100_250",  '%s/%s/%s")'%(srm,mgQcdLoc,mgQcdDset%"100To250"), xs = {"LO":7e+06, "fakeNLO":7e+06*mgKFactor}["fakeNLO"])
qcd.add("qcd_mg_ht_250_500",  '%s/%s/%s")'%(srm,mgQcdLoc,mgQcdDset%"250To500"), xs = {"LO":171e+03, "fakeNLO":171e+03*mgKFactor}["fakeNLO"])
qcd.add("qcd_mg_ht_500_1000", '%s/%s/%s")'%(srm,mgQcdLoc,mgQcdDset%"500To1000"),xs = {"LO":5200, "fakeNLO":5200*mgKFactor}["fakeNLO"])
qcd.add("qcd_mg_ht_1000_inf", '%s/%s/%s")'%(srm,mgQcdLoc,mgQcdDset%"1000"),     xs = {"LO":83, "fakeNLO":83*mgKFactor}["fakeNLO"])

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
    qcd.add("qcd_py6_pt_%s"%("%d_%d"%(low,high) if high else "%d"%low),
           '%s/%s/%s")'%(srm,py6Loc,py6Dset%("%dto%d"%(low,high) if high else "%d"%low)),
           xs = xs)

py6FJLoc = '/bbetchar//ICF/automated/2011_06_12_05_20_30/'
py6Summer11Dset = "/QCD_Pt-%s_TuneZ2_7TeV_pythia6.Summer11-PU_S3_START42_V11-v2.AODSIM"
for low,high,xs in zip(bins[:-1],bins[1:],xss[:-1]) :
    qcd.add("qcd_py6fj_pt_%s"%("%d_%d"%(low,high) if high else str(low)),
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
    qcd.add("qcd_py6fjmu_pt_%s"%("%d_%d"%(low,high) if high else str(low)),
           '%s/%s/%s")'%(srm,py6FJmuLoc,formats[form]%("%dto%d"%(low,high) if high else str(low))),
           xs = xs)

loc = '/bbetchar/ICF/automated/2011_12_08_09_09_32/'
form = 'QCD_Pt-%s_MuPt5Enriched_TuneZ2_7TeV-pythia6.Fall11-PU_S6_START42_V14B-v1.AODSIM")'
for low,high,xs,_ in zip(bins[:-1],bins[1:],xss[:-1],forms[:-1]) :
    qcd.add("qcd_mu_%s"%("%d_%d"%(low,high) if high else str(low)),
            srm+loc+form%(('%dto%d'%(low,high)) if high else str(low)),
            xs = xs)

# QCD (HT binned)
qcdDset1 = "QCD_TuneZ2_HT-%s_7TeV-madgraph.Summer11-PU_S4_START42_V11-v1.AODSIM"
qcdDset3 = "QCD_TuneZ2_HT-%s_7TeV-madgraph.Summer11-PU_S4_START42_V11-v3.AODSIM"
qcdDir1 = "/dburton//ICF/automated/2011_10_26_12_39_58/"
qcdDir2 = "/dburton//ICF/automated/2011_11_24_13_30_10/"

qcd.add("qcd_mg_ht_100_250_summer11",  '%s/%s/%s")'%(srm, qcdDir1, qcdDset1%"100To250"), xs = {"LO":4194000.0, "fakeNLO":4194000.0*mgKFactor}["fakeNLO"])
qcd.add("qcd_mg_ht_250_500_summer11",  '%s/%s")'   %(srm, qcdDir2,                    ), xs = {"LO": 198500.0, "fakeNLO": 198500.0*mgKFactor}["fakeNLO"])
qcd.add("qcd_mg_ht_500_1000_summer11", '%s/%s/%s")'%(srm, qcdDir1, qcdDset1%"500To1000"),xs = {"LO":   5856.0, "fakeNLO":   5856.0*mgKFactor}["fakeNLO"])
qcd.add("qcd_mg_ht_1000_inf_summer11", '%s/%s/%s")'%(srm, qcdDir1, qcdDset1%"1000"),     xs = {"LO":    122.6, "fakeNLO":    122.6*mgKFactor}["fakeNLO"])

# QCD (HT binned) skims
dir = "/vols/cms02/elaird1/29_skims/04_photons/v8"
l = 'utils.fileListFromDisk(isDirectory = False, location = '
qcd.add("qcd_mg_ht_100_250_summer11_skim",   '%s"%s/qcd_mg_ht_100_250_summer11_*_skim.root")'%(l, dir), xs = 3.180463e-06 * 5.326380e+06)
qcd.add("qcd_mg_ht_250_500_summer11_skim",   '%s"%s/qcd_mg_ht_250_500_summer11_*_skim.root")'%(l, dir), xs = 2.491150e-04 * 2.520950e+05)
qcd.add("qcd_mg_ht_500_1000_summer11_skim", '%s"%s/qcd_mg_ht_500_1000_summer11_*_skim.root")'%(l, dir), xs = 2.405978e-04 * 7.437120e+03)
qcd.add("qcd_mg_ht_1000_inf_summer11_skim", '%s"%s/qcd_mg_ht_1000_inf_summer11_*_skim.root")'%(l, dir), xs = 1.675973e-04 * 1.557020e+02)

