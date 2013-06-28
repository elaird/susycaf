from supy.samples import SampleHolder
from supy.sites import pnfs
pnfs = pnfs()
ewk17 = SampleHolder()

### Wjets ###
ewk17.add("wj_lv_mg_ht_5_250_skim",     '%s/yeshaq/ICF/supy-output/wj_skim/")'%pnfs, xs = {"LO":30400.0, "NLO":36257.2}["NLO"])
ewk17.add("wj_lv_mg_ht_0_250_other_reqs",   '%s/dburton/WJet_Skim_v2/")'%pnfs, xs = 36140.2*4135747.0/18273090.0) #wrong XS!

S10_incl = "WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM"
S10_excl = "/WJetsToLNu_HT-%s_8TeV-madgraph_v2.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM"
ewk17.add("wj_lv_mg_ht_10_150",    '%s/yeshaq/ICF/supy-output/wj_skim_v2/")'%pnfs, xs = {"NLO":9090.43}["NLO"])
ewk17.add("wj_lv_mg_ht_150_200.job663", '%s/agapitos//ICF/automated/2013_05_15_18_11_50/")'%pnfs, xs = 290.69475) # see comment 1
ewk17.add("wj_lv_mg_ht_200_250.job672", '%s/agapitos//ICF/automated/2013_05_15_17_18_45/")'%pnfs, xs = 123.3417) # see comment 2
ewk17.add("wj_lv_mg_ht_250_300.job498", '%s/karage//ICF/automated/2012_11_29_20_11_51/%s")'%(pnfs, S10_excl)%"250To300", xs = {"LO":48.01, "NLO":57.26}["NLO"])
ewk17.add("wj_lv_mg_ht_300_400.job498", '%s/karage//ICF/automated/2012_11_29_20_11_51/%s")'%(pnfs, S10_excl)%"300To400", xs = {"LO":38.30, "NLO":45.68}["NLO"])
ewk17.add("wj_lv_mg_ht_400_inf.job498", '%s/karage//ICF/automated/2012_11_29_20_11_51/%s")'%(pnfs, S10_excl)%"400ToInf", xs = {"LO":25.22, "NLO":30.08}["NLO"])
ewk17.add("wj_lv_mg_ht_incl.job363",    '%s/clucas/ICF//automated/2012_09_21_09_36_56/%s")'%(pnfs, S10_incl), xs = {"LO":30400.0, "NLO":36257.2}["NLO"])


for part in [1,2,3,4,5] :
    ewk17.add("wj_lv_mg_ht_incl_v2.job673_part%i"%part, '%s/yeshaq//ICF/automated/2013_05_23_14_13_07/part%i")'%(pnfs,part), xs = {"LO":30400.0, "NLO":36257.2}["NLO"])

#comment 1
# CrossSection = ( NNLO_W_inclusive_from[1] / LO_W_inclusive_from_PREP[2] ) * LO_W_exclusive_from_PREP: [3]
                                                        
#[1] https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat8TeV
#[2] http://cms.cern.ch/iCMS/prep/requestmanagement?dsn=WJetsToLNu_TuneZ2Star_8TeV-madgraph*&campid=Summer12
#[3] http://cms.cern.ch/iCMS/prep/requestmanagement?dsn=WJetsToLNu_HT*200*madgraph*&campid=Summer12
                                                        
#comment 2 
#CrossSection = 123.3417 #NNLO correct with the disconnection correction factor see: https://twiki.cern.ch/twiki/pub/CMS/SusyRA1Material/FixSmapleProblems_24May2013.pdf, correction factor is 1.1074
#[1] https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat8TeV
#[2] http://cms.cern.ch/iCMS/prep/requestmanagement?dsn=WJetsToLNu_TuneZ2Star_8TeV-madgraph*&campid=Summer12
#[3] http://cms.cern.ch/iCMS/prep/requestmanagement?dsn=WJetsToLNu_HT*200*madgraph*&campid=Summer12


### Zinv ###

#https://twiki.cern.ch/twiki/bin/view/CMS/HiggsMCProductionSummer12#VBF_H_bb_POWHEG
ewk17.add("zinv_hbb_125_powheg.job342", '%s/yeshaq//ICF/automated/2012_08_31_15_26_30/")'%pnfs, xs = 0.0361)

Zinv_S10_excl = "ZJetsToNuNu_%s_HT_%s_TuneZ2Star_8TeV_madgraph.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM"
ewk17.add("zinv_mg_ht_50_100.job407", '%s/clucas//ICF/automated/2012_09_23_19_53_51/")'%pnfs, xs = {"NNLO":452.75, "LO":381.2}["NNLO"])
ewk17.add("zinv_mg_ht_100_200.job365", '%s/clucas//ICF/automated/2012_09_21_10_32_59/%s")'%(pnfs, Zinv_S10_excl)%(100,200), xs = {"NNLO":190.39, "LO":160.3}["NNLO"])
ewk17.add("zinv_mg_ht_200_400.job365", '%s/clucas//ICF/automated/2012_09_21_10_32_59/%s")'%(pnfs, Zinv_S10_excl)%(200,400), xs = {"NNLO":49.2776, "LO":41.49}["NNLO"])
ewk17.add("zinv_mg_ht_400_inf.job365", '%s/clucas//ICF/automated/2012_09_21_10_32_59/%s")'%(pnfs, Zinv_S10_excl)%(400,"inf"), xs = {"NNLO":6.2639, "LO":5.274}["NNLO"])

