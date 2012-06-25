from supy.samples import SampleHolder
from supy.sites import pnfs
pnfs = pnfs()
ht17 = SampleHolder()

a = ", alwaysUseLastAttempt = True"
#L1FJL2L3Residual
ht17.add("JetHT.Run2012B-PromptReco-v1.AOD.job217",  '%s/zmeng//ICF/automated/2012_06_13_00_19_12/JetHT.Run2012B-PromptReco-v1.AOD/")'%pnfs, lumi = 1.0) #job 217,   226/226 comp.

#ht.add("HT.Run2011B-PromptReco-v1.AOD.job570",  '%s/bm409//ICF/automated/2011_10_17_12_55_58/HT.Run2011B-PromptReco-v1.AOD")'%srm,         lumi = 1.0) #job 570,   82/432  comp.
