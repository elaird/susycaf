from supy.samples import SampleHolder
from supy.sites import pnfs, eos
pnfs = pnfs()
eos = eos()
dyll17 = SampleHolder()

dyChris = "clucas//ICF/automated/2012_09_21_09_51_06" 
dyll17.add("dyll_M-10To50_mg", '%s/%s/DYJetsToLL_M-10To50filter_8TeV-madgraph.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM/")'%(pnfs,dyChris), xs = 11050.0)
dyll17.add("dyll_HT-200To400_mg", '%s/%s/DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM/")'%(pnfs,dyChris), xs = 2.826)
dyll17.add("dyll_HT-400ToInf_mg", '%s/%s/DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM/")'%(pnfs,dyChris), xs = 19.73 )
dyll17.add("dyll_M-50_mg", '%s/%s/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM/")'%(pnfs,dyChris), xs = 2950.0)
skimFactor = 8111524./30171503.
dyll17.add("dyll_HT_10To200_M-10To50", '%s/clucas/Parked13/DYJets_10to200_M_10to50/")' % eos, xs = 13124.07)
dyll17.add("dyll_HT_10To200_M-50", '%s/clucas/Parked13/DYJets_10to200/")'% eos, xs = 3503.71*skimFactor)
dyll17.add("dyll_HT_200To400_M-50", '%s/clucas/Parked13/DYJets_200to400_Combined/")'% eos, xs = 24.29)
dyll17.add("dyll_HT_400ToInf_M-50", '%s/clucas/Parked13/DYJets_400toinf_Combined/")'% eos, xs = 3.3564)           

dyNChris = "/clucas/ICF/automated/2013_07_30_14_43_01/"
dyll17.add("dy1ll_M-50", '%s/%s/DY1JetsToLL_M-50_TuneZ2Star_8TeV-madgraph.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM/")'   %(pnfs,dyNChris), xs= 629.19)
dyll17.add("dy2ll_M-50", '%s/%s/_.DY2JetsToLL_M-50_TuneZ2Star_8TeV-madgraph.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM/")' %(pnfs,dyNChris), xs = 223.32)
dyll17.add("dy3ll_M-50", '%s/%s/_.DY3JetsToLL_M-50_TuneZ2Star_8TeV-madgraph.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM/")' %(pnfs,dyNChris), xs = 62.93)
dyll17.add("dy4ll_M-50", '%s/%s/_.DY4JetsToLL_M-50_TuneZ2Star_8TeV-madgraph.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM/")' %(pnfs,dyNChris), xs = 28.43)
