from supy.samples import SampleHolder
from supy.sites import pnfs
pnfs = pnfs()
ht17 = SampleHolder()

a = ", alwaysUseLastAttempt = True"

#53X
ht17.add("HTMHT.Run2012B-13Jul2012-v1.AOD.job358", '%s/yeshaq//ICF/automated/2012_09_19_01_03_03/HTMHT.Run2012B-13Jul2012-v1.AOD/")'%pnfs, lumi = 1.0)
ht17.add("HTMHT.Run2012C-24Aug2012-v1.AOD.job361", '%s/yeshaq//ICF/automated/2012_09_19_01_59_39/HTMHT.Run2012C-24Aug2012-v1.AOD/")'%pnfs, lumi = 1.0)
ht17.add("HTMHT.Run2012C-PromptReco-v2.AOD.job360", '%s/yeshaq//ICF/automated/2012_09_19_01_39_12/HTMHT.Run2012C-PromptReco-v2.AOD/")'%pnfs, lumi = 1.0)

ht17.add("HT.Run2012A-13Jul2012-v1.AOD.job358", '%s/yeshaq//ICF/automated/2012_09_19_01_03_03/HT.Run2012A-13Jul2012-v1.AOD/")'%pnfs, lumi = 1.0)
ht17.add("HT.Run2012A-recover-06Aug2012-v1.AOD.job359", '%s/yeshaq//ICF/automated/2012_09_19_01_25_46/HT.Run2012A-recover-06Aug2012-v1.AOD/")'%pnfs, lumi = 1.0)

for era,lum in zip(["B","C","D"],[4427.0,6893.0,7263.0])  :
    ht17.add("HTMHTParked.Run2012%s-22Jan2013-v1.job649"%era, '%s/yeshaq//ICF/automated/2013_04_19_17_30_49/HTMHTParked.Run2012%s-22Jan2013-v1.AOD/")'%(pnfs,era), lumi=lum) 

ht17.add("HTMHTParked_ICF_sync_test", '%s/yeshaq//ICF/supy-output/ICF_sync_test/")'%pnfs, lumi = 1.0)
