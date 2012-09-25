from supy.samples import SampleHolder
from supy.sites import pnfs
pnfs = pnfs()
ewk17 = SampleHolder()

four = "WJetsToLNu_HT-400ToInf_8TeV-madgraph.Summer12-PU_S7_START52_V9-v1.AODSIM/"
incl = "WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball.Summer12-PU_S7_START52_V9-v1.AODSIM/"
ewk17.add("wj_lv_mg_ht_0_250_other_reqs",   '%s/dburton/WJet_Skim_v2/")'%pnfs, xs = 36140.2*4135747.0/18273090.0) #wrong XS!
ewk17.add("wj_lv_mg_ht_250_300", '%s/clucas/ICF/automated/2012_06_20_16_29_08/")'%pnfs,           xs = 48.01)
ewk17.add("wj_lv_mg_ht_300_400", '%s/clucas/ICF/automated/2012_06_21_02_04_54/")'%pnfs,           xs = 38.30)
ewk17.add("wj_lv_mg_ht_400_inf", '%s/clucas/ICF/automated/2012_06_05_09_30_46/%s")'%(pnfs, four), xs = 25.22)
ewk17.add("wj_lv_mg_ht_incl",    '%s/clucas/ICF/automated/2012_06_05_09_30_46/%s")'%(pnfs, incl), xs = {"LO":30400.0, "NLO":36257.2}["NLO"])

ewk17.add("zinv_mg_ht_50_100.job214", '%s/clucas//ICF/automated/2012_06_10_22_27_16/")'%pnfs, xs = {"NNLO":452.75, "LO":381.2}["NNLO"])
ewk17.add("zinv_mg_ht_100_200.job234", '%s/clucas//ICF/automated/2012_06_17_22_29_52/")'%pnfs, xs = {"NNLO":190.39, "LO":160.3}["NNLO"])
ewk17.add("zinv_mg_ht_200_400.job233", '%s/clucas//ICF/automated/2012_06_15_00_32_54/")'%pnfs, xs = {"NNLO":49.2776, "LO":41.49}["NNLO"])
ewk17.add("zinv_mg_ht_400_inf.job213", '%s/clucas//ICF/automated/2012_06_10_22_08_24/")'%pnfs, xs = {"NNLO":6.2639, "LO":5.274}["NNLO"])

vv = '%s/clucas//ICF/automated/2012_06_05_09_30_46/'%pnfs
ewk17.add("ww_py.job188", '%s/WW_TuneZ2star_8TeV_pythia6_tauola.Summer12-PU_S7_START52_V9-v1.AODSIM/")'%vv, xs = {"NNLO":57.1097, "LO":33.61}["NNLO"])
ewk17.add("wz_py.job188", '%s/WZ_TuneZ2star_8TeV_pythia6_tauola.Summer12-PU_S7_START52_V9-v1.AODSIM/")'%vv, xs = {"NNLO":32.3161, "LO":12.63}["NNLO"])
ewk17.add("zz_py.job188", '%s/ZZ_TuneZ2star_8TeV_pythia6_tauola.Summer12-PU_S7_START52_V9-v1.AODSIM/")'%vv, xs = {"NNLO":8.25561, "LO":5.196}["NNLO"])

#https://twiki.cern.ch/twiki/bin/view/CMS/HiggsMCProductionSummer12#VBF_H_bb_POWHEG
ewk17.add("zinv_hbb_125_powheg.job342", '%s/yeshaq//ICF/automated/2012_08_31_15_26_30/")'%pnfs, xs = 0.0361)
