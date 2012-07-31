from supy.samples import SampleHolder
from supy.sites import pnfs
pnfs = pnfs()
muon17 = SampleHolder()

muon17.add("SingleMu.Run2012A-PromptReco-v1.AOD.job228", '%s/zmeng//ICF/automated/2012_06_14_11_22_04/SingleMu.Run2012A-PromptReco-v1.AOD/", alwaysUseLastAttempt = True)'%pnfs, lumi = 1 ) #/pb
muon17.add("SingleMu.Run2012B-PromptReco-v1.AOD.job229", '%s/zmeng//ICF/automated/2012_06_14_14_26_12/SingleMu.Run2012B-PromptReco-v1.AOD/", alwaysUseLastAttempt = True)'%pnfs, lumi = 1 ) #/pb
muon17.add("SingleMu.Run2012B-PromptReco-v1.AOD.job238", '%s/zmeng//ICF/automated/2012_06_22_14_25_16/SingleMu.Run2012B-PromptReco-v1.AOD/", alwaysUseLastAttempt = True)'%pnfs, lumi = 1 ) #/pb

