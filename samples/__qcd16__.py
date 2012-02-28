from supy.samples import SampleHolder
from supy.sites import srm
qcd16 = SampleHolder()


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

#py6Loc = '/bbetchar//ICF/automated/2011_04_07_19_50_45/'
#py6Dset = "/QCD_Pt_%s_TuneZ2_7TeV_pythia6.Spring11-PU_S1_START311_V1G1-v1.AODSIM"
#for low,high,xs in zip(bins[:-1],bins[1:],xss) :
#    qcd16.add("qcd_py6_pt_%s"%("%d_%d"%(low,high) if high else "%d"%low),
#           '%s/%s/%s")'%(srm,py6Loc,py6Dset%("%dto%d"%(low,high) if high else "%d"%low)),
#           xs = xs)

#### QCD mu enriched ####
bins,xss,corrupt = zip(*tuple(
    [(15,  5.792e+08 * 0.00254, []), # xs * filter efficiency
     (20,  2.363e+08 * 0.00518, []),
     (30,  5.307e+07 * 0.01090, []),
     (50,  6.351e+06 * 0.02274, []),
     (80,  7.851e+05 * 0.03700, []),
     (120, 9.295e+04 * 0.04777, []),
     (150, 4.758e+04 * 0.05964, [41]),
     (None,None,None)]))
loc = '/bbetchar/ICF/automated/2012_02_27_01_45_05/'
form = 'QCD_Pt-%s_MuPt5Enriched_TuneZ2_7TeV-pythia6.Fall11-PU_S6_START44_V9B-v1.AODSIM", itemsToSkip = %s)'
for low,high,xs,cor in zip(bins[:-1],bins[1:],xss[:-1],corrupt[:-1]) :
    cstr = ','.join('"SusyCAF_Tree_%d"'%d for d in cor).join(['[',']'])
    qcd16.add("qcd_mu_%s"%("%d_%d"%(low,high) if high else str(low)),
              srm+loc+form%(('%dto%d'%(low,high),cstr) if high else (str(low),cstr)),
              xs = xs)
