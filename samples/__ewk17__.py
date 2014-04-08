from supy.samples import SampleHolder
from supy.sites import pnfs, eos
pnfs = pnfs()
eos = eos()
ewk17 = SampleHolder()

### Wjets ###
ewk17.add("wj_lv_mg_ht_5_250_skim",     '%s/yeshaq/ICF/supy-output/wj_skim/")'%pnfs, xs = {"LO":30400.0, "NLO":36257.2}["NLO"])
ewk17.add("wj_lv_mg_ht_0_250_other_reqs",   '%s/dburton/WJet_Skim_v2/")'%pnfs, xs = 36140.2*4135747.0/18273090.0) #wrong XS!

S10_incl = "WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM"
S10_excl = "/WJetsToLNu_HT-%s_8TeV-madgraph_v2.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM"
ewk17.add("wj_lv_mg_ht_10_150", '%s/yeshaq/ICF/supy-output/wj_skim_v2", pruneList = False)'%pnfs, xs = {"NLO":9090.43}["NLO"])
ewk17.add("wj_lv_mg_ht_150_200.job663", '%s/agapitos//ICF/automated/2013_05_15_18_11_50/")'%pnfs, xs = 290.69475) # see comment 1
ewk17.add("wj_lv_mg_ht_200_250.job672", '%s/agapitos//ICF/automated/2013_05_15_17_18_45/")'%pnfs, xs = 123.3417) # see comment 2
ewk17.add("wj_lv_mg_ht_250_300.job498", '%s/karage//ICF/automated/2012_11_29_20_11_51/%s")'%(pnfs, S10_excl)%"250To300", xs = {"LO":48.01, "NLO":57.26}["NLO"])
ewk17.add("wj_lv_mg_ht_300_400.job498", '%s/karage//ICF/automated/2012_11_29_20_11_51/%s")'%(pnfs, S10_excl)%"300To400", xs = {"LO":38.30, "NLO":45.68}["NLO"])
ewk17.add("wj_lv_mg_ht_400_inf.job498", '%s/karage//ICF/automated/2012_11_29_20_11_51/%s")'%(pnfs, S10_excl)%"400ToInf", xs = {"LO":25.22, "NLO":30.08}["NLO"])
ewk17.add("wj_lv_mg_ht_incl.job363",    '%s/clucas/ICF//automated/2012_09_21_09_36_56/%s")'%(pnfs, S10_incl), xs = {"LO":30400.0, "NLO":36257.2}["NLO"])


for part in [1,2,3,4,5] :
    ewk17.add("wj_lv_mg_ht_incl_v2.job673_part%i"%part, '%s/yeshaq//ICF/automated/2013_05_23_14_13_07/part%i")'%(pnfs,part), xs = {"LO": 30400.0, "NNLO":37509}["NNLO"])

for part in [1,2,3,4] :
    ewk17.add("wj_lv_mg_ht_incl_v2.job50_part%i"%part, '%s/zmeng/ICF/automated/2013_09_02_15_42_12/part%i")'%(pnfs,part), xs = {"LO":30400.0, "NNLO":37509}["NNLO"])

##The LO xs are from PREP
skimFactor = 13992013./57645905.
ewk17.add("wj_lv_mg_ht_incl_LO"    , '%s/clucas/Parked13/WJets_inc/")'      % eos, xs = {"LO":30400}["LO"])
ewk17.add("wj_lv_mg_ht_10To150_LO" , '%s/clucas/Parked13/WJets_10to150/")'  % eos, xs = {"LO":30400*skimFactor}["LO"])
ewk17.add("wj_lv_mg_ht_150To200_LO", '%s/clucas/Parked13/WJets_150to200/")' % eos, xs = {"LO":235.6}["LO"])
ewk17.add("wj_lv_mg_ht_200To250_LO", '%s/clucas/Parked13/WJets_200to250/")' % eos, xs = {"LO":90.27}["LO"])
ewk17.add("wj_lv_mg_ht_250To300_LO", '%s/clucas/Parked13/WJets_250to300/")' % eos, xs = {"LO":48.01}["LO"])
ewk17.add("wj_lv_mg_ht_300To400_LO", '%s/clucas/Parked13/WJets_300to400/")' % eos, xs = {"LO":38.30}["LO"])
ewk17.add("wj_lv_mg_ht_400ToInf_LO", '%s/clucas/Parked13/WJets_400toinf/")' % eos, xs = {"LO":25.22}["LO"])

kFactor = 37509./30400.
ewk17.add("wj_lv_mg_ht_incl"    , '%s/clucas/Parked13/WJets_inc/")'      % eos, xs = {"NNLO":30400*kFactor}["NNLO"])
ewk17.add("wj_lv_mg_ht_10To150" , '%s/clucas/Parked13/WJets_10to150/")'  % eos, xs = {"NNLO":30400*kFactor*skimFactor}["NNLO"])
ewk17.add("wj_lv_mg_ht_150To200", '%s/clucas/Parked13/WJets_150to200/")' % eos, xs = {"NNLO":235.6*kFactor}["NNLO"])
ewk17.add("wj_lv_mg_ht_200To250", '%s/clucas/Parked13/WJets_200to250/")' % eos, xs = {"NNLO":90.27*kFactor}["NNLO"])
ewk17.add("wj_lv_mg_ht_250To300", '%s/clucas/Parked13/WJets_250to300/")' % eos, xs = {"NNLO":48.01*kFactor}["NNLO"])
ewk17.add("wj_lv_mg_ht_300To400", '%s/clucas/Parked13/WJets_300to400/")' % eos, xs = {"NNLO":38.30*kFactor}["NNLO"]) 
ewk17.add("wj_lv_mg_ht_400ToInf", '%s/clucas/Parked13/WJets_400toinf/")' % eos, xs = {"NNLO":25.22*kFactor}["NNLO"])

ewk17.add("wj_lv_mg_ht_N1", '%s/clucas/HCP_transfer/W1JetsToLNu/")' % eos, xs = {"IC":6381.2, "PREP": 5500, "Other":5400, "NNLO":5500*kFactor}["NNLO"])
ewk17.add("wj_lv_mg_ht_N2", '%s/clucas/HCP_transfer/W2JetsToLNu/")' % eos, xs = {"IC":2039.8, "PREP": 1800, "Other":1750, "NNLO":1800*kFactor}["NNLO"])
ewk17.add("wj_lv_mg_ht_N3", '%s/clucas/HCP_transfer/W3JetsToLNu/")' % eos, xs = {"IC":612.5, "PREP": 520, "Other":519., "NNLO":520*kFactor}["NNLO"])
ewk17.add("wj_lv_mg_ht_N4", '%s/clucas/HCP_transfer/W4JetsToLNu/")' % eos, xs = {"IC":251.0, "PREP": 210, "Other":214., "NNLO":210*kFactor}["NNLO"])

ewk17.add("wj_lv_mg_ht_N1_LO", '%s/clucas/HCP_transfer/W1JetsToLNu/")' % eos, xs = {"IC":6381.2, "PREP": 5500, "Other":5400, "NNLO":5500*kFactor}["PREP"])
ewk17.add("wj_lv_mg_ht_N2_LO", '%s/clucas/HCP_transfer/W2JetsToLNu/")' % eos, xs = {"IC":2039.8, "PREP": 1800, "Other":1750, "NNLO":1800*kFactor}["PREP"])
ewk17.add("wj_lv_mg_ht_N3_LO", '%s/clucas/HCP_transfer/W3JetsToLNu/")' % eos, xs = {"IC":612.5, "PREP": 520, "Other":519., "NNLO":520*kFactor}["PREP"])
ewk17.add("wj_lv_mg_ht_N4_LO", '%s/clucas/HCP_transfer/W4JetsToLNu/")' % eos, xs = {"IC":251.0, "PREP": 210, "Other":214., "NNLO":210*kFactor}["PREP"])
#
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

ext_suffix = "ZJetsToNuNu_%s_HT_%s_TuneZ2Star_8TeV_madgraph_ext.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM/"
ewk17.add("zinv_mg_ht_50_100_ext.job500", '%s/karage//ICF/automated/2012_11_29_23_55_17/%s")' % (pnfs, ext_suffix % (50, 100)), xs={"NNLO":452.75, "LO":381.2}["NNLO"])
ewk17.add("zinv_mg_ht_100_200_ext.job680", '%s/clucas//ICF/automated/2013_06_17_17_08_07/%s")' % (pnfs,
                                                                                                  ext_suffix.replace("V7A","V7C") % (100, 200)),xs={"NNLO":190.39,
                                                                                                                                                    "LO":160.3}["NNLO"])
ewk17.add("zinv_mg_ht_200_400_ext.job500", '%s/karage//ICF/automated/2012_11_29_23_55_17/%s")' % (pnfs, ext_suffix % (200, 400)), xs={"NNLO":49.2776, "LO":41.49}["NNLO"])
ewk17.add("zinv_mg_ht_400_inf_ext.job500", '%s/karage//ICF/automated/2012_11_29_23_55_17/%s")' % (pnfs, ext_suffix % (400, "inf")), xs={"NNLO":6.2639, "LO":5.274}["NNLO"])

ht = [(50,100),(100,200),(200,400),(400,"inf")]
xss = [{"NNLO":452.75, "LO":381.2}["NNLO"],
       {"NNLO":190.39, "LO":160.3}["NNLO"],
       {"NNLO":49.2776, "LO":41.49}["NNLO"],
       {"NNLO":6.2639, "LO":5.274}["NNLO"]]
for ht,xs in zip(ht,xss):
    ewk17.add("zinv_mg_ht_%s_%s" % ht, '%s/clucas/Parked13/ZJets_%sto%s_Combined/")' % (eos, ht[0], ht[1]), xs=xs)

diBos = "%s_TuneZ2star_8TeV_pythia6_tauola.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM/"
diBosons ={"ZZ":8.25561, "WZ":32.3161, "WW":57.1097}
for key in diBosons:
    ewk17.add("%s_pythia6.job370" % key, '%s/zmeng/ICF/automated/2012_09_21_17_43_41/%s")' % (pnfs, diBos % key), xs = diBosons[key])
    ewk17.add("%s_py6" % key, '%s/clucas/Parked13/%s")' % (eos, key), xs = diBosons[key])
