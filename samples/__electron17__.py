from supy.samples import SampleHolder
from supy.sites import pnfs
pnfs = pnfs()
electron17 = SampleHolder()

#2012
electron17.add("SingleElectron.Run2012B-PromptReco-v1.AOD.job228", '%s/zmeng//ICF/automated/2012_06_14_11_22_04/SingleElectron.Run2012B-PromptReco-v1.AOD/", alwaysUseLastAttempt = True)'%pnfs, lumi = 1 ) #/pb
electron17.add("SingleElectron.Run2012A-PromptReco-v1.AOD.job229", '%s/zmeng//ICF/automated/2012_06_14_14_26_12/SingleElectron.Run2012A-PromptReco-v1.AOD/", alwaysUseLastAttempt = True)'%pnfs, lumi = 1 ) #/pb
electron17.add("SingleElectron.Run2012B-PromptReco-v1.AOD.job238", '%s/zmeng//ICF/automated/2012_06_22_14_25_16/SingleElectron.Run2012B-PromptReco-v1.AOD/", alwaysUseLastAttempt = True)'%pnfs, lumi = 1 ) #/pb


