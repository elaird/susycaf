from supy.samples import SampleHolder
from supy.sites import pnfs
pnfs = pnfs()
dyll17 = SampleHolder()

dyChris = "clucas//ICF/automated/2012_09_21_09_51_06" 
dyll17.add("dyll_M-10To50_mg", '%s/%s/DYJetsToLL_M-10To50filter_8TeV-madgraph.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM/")'%(pnfs,dyChris), xs = 11050.0)
dyll17.add("dyll_HT-200To400_mg", '%s/%s/DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM/")'%(pnfs,dyChris), xs = 2.826)
dyll17.add("dyll_HT-400ToInf_mg", '%s/%s/DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM/")'%(pnfs,dyChris), xs = 19.73 )
dyll17.add("dyll_M-50_mg", '%s/%s/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM/")'%(pnfs,dyChris), xs = 2950.0)


