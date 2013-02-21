from supy.samples import SampleHolder
from supy.sites import pnfs
pnfs = pnfs()
ht17 = SampleHolder()

a = ", alwaysUseLastAttempt = True"
#L1FJL2L3Residual
ht17.add("HT.Run2012A-PromptReco-v1.AOD.job229",  '%s/zmeng//ICF/automated/2012_06_14_14_26_12/HT.Run2012A-PromptReco-v1.AOD/")'%pnfs, lumi = 1.0) #job 217,   226/226 comp.

ht17.add("HTMHT.Run2012B-PromptReco-v1.AOD.job228", '%s/zmeng//ICF/automated/2012_06_14_11_22_04/HTMHT.Run2012B-PromptReco-v1.AOD/")'%pnfs, lumi = 1.0)
ht17.add("HTMHT.Run2012B-PromptReco-v1.AOD.job238", '%s/zmeng//ICF/automated/2012_06_22_14_25_16/HTMHT.Run2012B-PromptReco-v1.AOD/")'%pnfs, lumi = 1.0)

#53X
ht17.add("HTMHT.Run2012B-13Jul2012-v1.AOD.job358", '%s/yeshaq//ICF/automated/2012_09_19_01_03_03/HTMHT.Run2012B-13Jul2012-v1.AOD/")'%pnfs, lumi = 1.0)
ht17.add("HTMHT.Run2012C-24Aug2012-v1.AOD.job361", '%s/yeshaq//ICF/automated/2012_09_19_01_59_39/HTMHT.Run2012C-24Aug2012-v1.AOD/")'%pnfs, lumi = 1.0)
ht17.add("HTMHT.Run2012C-PromptReco-v2.AOD.job360", '%s/yeshaq//ICF/automated/2012_09_19_01_39_12/HTMHT.Run2012C-PromptReco-v2.AOD/")'%pnfs, lumi = 1.0)

ht17.add("HT.Run2012A-13Jul2012-v1.AOD.job358", '%s/yeshaq//ICF/automated/2012_09_19_01_03_03/HT.Run2012A-13Jul2012-v1.AOD/")'%pnfs, lumi = 1.0)
ht17.add("HT.Run2012A-recover-06Aug2012-v1.AOD.job359", '%s/yeshaq//ICF/automated/2012_09_19_01_25_46/HT.Run2012A-recover-06Aug2012-v1.AOD/")'%pnfs, lumi = 1.0)







#ht.add("HT.Run2011B-PromptReco-v1.AOD.job570",  '%s/bm409//ICF/automated/2011_10_17_12_55_58/HT.Run2011B-PromptReco-v1.AOD")'%srm,         lumi = 1.0) #job 570,   82/432  comp.
