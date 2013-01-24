from supy.samples import SampleHolder
from supy.sites import pnfs
pnfs = pnfs()
mumu17 = SampleHolder()


#L1FJL2L3Residual
a = ", alwaysUseLastAttempt = True"
mumu17.add("DoubleMu.Run2012A-PromptReco-v1.AOD.job229", '%s/zmeng//ICF/automated/2012_06_14_14_26_12/DoubleMu.Run2012A-PromptReco-v1.AOD")'%pnfs, lumi = 1)
mumu17.add("DoubleMu.Run2012B-PromptReco-v1.AOD.job228", '%s/zmeng//ICF/automated/2012_06_14_11_22_04/DoubleMu.Run2012B-PromptReco-v1.AOD")'%pnfs, lumi = 1)
mumu17.add("DoubleMu.Run2012B-PromptReco-v1.AOD.job239", '%s/zmeng//ICF/automated/2012_06_22_14_25_16/DoubleMu.Run2012B-PromptReco-v1.AOD")'%pnfs, lumi = 1)


mumu17.add("DoubleMu.Run2012A-13Jul2012-v1.AOD.job375", '%s/karage/ICF/automated/2012_09_22_23_23_49/")'%pnfs, lumi = 1)
mumu17.add("DoubleMu.Run2012A-recover-06Aug2012-v1.AOD.job389", '%s/karage/ICF/automated/2012_09_22_23_09_55/")'%pnfs, lumi = 1)

mumu17.add("DoubleMu.Run2012B-13Jul2012-v4.AOD.job408", '%s/karage/ICF/automated/2012_09_24_00_17_01/")'%pnfs, lumi = 1)

mumu17.add("DoubleMu.Run2012C-24Aug2012-v1.AOD.job401", '%s/karage/ICF/automated/2012_09_22_12_21_59/")'%pnfs, lumi = 1)
mumu17.add("DoubleMu.Run2012C-PromptReco-v2.AOD.job395", '%s/karage/ICF/automated/2012_09_23_00_37_41/")'%pnfs, lumi = 1)

mumu17.add("DoubleMu.Run2012D-PromptReco-v1.AOD.job508", '%s/yeshaq/ICF/automated/2012_12_04_16_51_31/DoubleMu.Run2012D-PromptReco-v1.AOD")'%pnfs, lumi = 1)
mumu17.add("DoubleMu.Run2012D-PromptReco-v1.AOD.job529", '%s/yeshaq/ICF/automated/2013_01_18_15_15_54/")'%pnfs, lumi = 1) #catchup

