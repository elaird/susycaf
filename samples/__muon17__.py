from supy.samples import SampleHolder
from supy.sites import pnfs, eos
pnfs = pnfs()
eos = eos()
muon17 = SampleHolder()


muon17.add("SingleMu.Run2012A-PromptReco-v1.AOD.job228", '%s/zmeng//ICF/automated/2012_06_14_11_22_04/SingleMu.Run2012A-PromptReco-v1.AOD/", alwaysUseLastAttempt = True)'%pnfs, lumi = 1 ) #/pb
muon17.add("SingleMu.Run2012B-PromptReco-v1.AOD.job229", '%s/zmeng//ICF/automated/2012_06_14_14_26_12/SingleMu.Run2012B-PromptReco-v1.AOD/", alwaysUseLastAttempt = True)'%pnfs, lumi = 1 ) #/pb
muon17.add("SingleMu.Run2012B-PromptReco-v1.AOD.job238", '%s/zmeng//ICF/automated/2012_06_22_14_25_16/SingleMu.Run2012B-PromptReco-v1.AOD/", alwaysUseLastAttempt = True)'%pnfs, lumi = 1 ) #/pb

muon17.add("SingleMu.Run2012A-13Jul2012-v1.AOD.job358", '%s/yeshaq/ICF/automated/2012_09_19_01_03_03/SingleMu.Run2012A-13Jul2012-v1.AOD/")'%pnfs, lumi = 1 ) #/pb
muon17.add("SingleMu.Run2012A-13Jul2012-v1.AOD.job467", '%s/yeshaq/ICF/automated/2012_10_15_19_19_13/")'%pnfs, lumi = 1 ) #/pb

muon17.add("SingleMu.Run2012A-recover-06Aug2012-v1.AOD.job359",'%s/yeshaq/ICF/automated/2012_09_19_01_25_46/SingleMu.Run2012A-recover-06Aug2012-v1.AOD/")'%pnfs, lumi = 1 ) #/pb
muon17.add("SingleMu.Run2012A-recover-06Aug2012-v1.AOD.job477", '%s/yeshaq/ICF/automated/2012_10_15_20_12_27/")'%pnfs, lumi = 1 ) #/pb

muon17.add("SingleMu.Run2012B-13Jul2012-v1.AOD.job358", '%s/yeshaq/ICF/automated/2012_09_19_01_03_03/SingleMu.Run2012B-13Jul2012-v1.AOD/")'%pnfs, lumi = 1 ) #/pb
muon17.add("SingleMu.Run2012B-13Jul2012-v1.AOD.job461", '%s/yeshaq/ICF/automated/2012_10_15_18_00_32/")'%pnfs, lumi = 1 ) #/pb

muon17.add("SingleMu.Run2012C-24Aug2012-v1.AOD.job361", '%s/yeshaq/ICF/automated/2012_09_19_01_59_39/SingleMu.Run2012C-24Aug2012-v1.AOD/")'%pnfs, lumi = 1 ) #/pb
muon17.add("SingleMu.Run2012C-24Aug2012-v1.AOD.job470", '%s/yeshaq/ICF/automated/2012_10_15_19_40_29/")'%pnfs, lumi = 1 ) #/pb

muon17.add("SingleMu.Run2012C-PromptReco-v2.AOD.job360", '%s/yeshaq/ICF/automated/2012_09_19_01_39_12/SingleMu.Run2012C-PromptReco-v2.AOD/")'%pnfs, lumi = 1 ) #/pb
muon17.add("SingleMu.Run2012C-PromptReco-v2.AOD.job474", '%s/yeshaq/ICF/automated/2012_10_15_20_00_38/")'%pnfs, lumi = 1 ) #/pb

muon17.add("SingleMu.Run2012D-PromptReco-v1.AOD.job508", '%s/yeshaq/ICF/automated/2012_12_04_16_51_31/SingleMu.Run2012D-PromptReco-v1.AOD/")'%pnfs, lumi = 1 ) #/pb
muon17.add("SingleMu.Run2012D-PromptReco-v1.AOD.job525", '%s/yeshaq/ICF/automated/2013_01_18_14_00_11/")'%pnfs, lumi = 1 ) #/pb

eras = {"A":887.6540,"B":4419.0,"C":7119.0,"D":7295.0} 
for era in eras:
    muon17.add("SingleMu.Run2012%s-22Jan2013" % era, '%s/clucas/Parked13/SingleMu_Run2012%s_22Jan2013/")'% (eos, era), lumi = eras[era]) #/pb

