from supy.samples import SampleHolder
from supy.sites import pnfs
pnfs = pnfs()
mumu17 = SampleHolder()


#L1FJL2L3Residual
a = ", alwaysUseLastAttempt = True"
mumu17.add("DoubleMu.Run2012A-PromptReco-v1.AOD.job229", '%s/zmeng//ICF/automated/2012_06_14_14_26_12/DoubleMu.Run2012A-PromptReco-v1.AOD")'%pnfs, lumi = 1)
mumu17.add("DoubleMu.Run2012B-PromptReco-v1.AOD.job228", '%s/zmeng//ICF/automated/2012_06_14_11_22_04/DoubleMu.Run2012B-PromptReco-v1.AOD")'%pnfs, lumi = 1)
mumu17.add("DoubleMu.Run2012B-PromptReco-v1.AOD.job239", '%s/zmeng//ICF/automated/2012_06_22_14_25_16/DoubleMu.Run2012B-PromptReco-v1.AOD")'%pnfs, lumi = 1)




