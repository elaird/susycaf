from supy.samples import SampleHolder
from supy.sites import pnfs
pnfs = pnfs()
photon17 = SampleHolder()

#photon17.add("Photon.Run2012A-PromptReco-v1.AOD.job29", '%s/zmeng//ICF/automated/2012_04_13_19_13_29/")'%pnfs, lumi = 1.0)
#photon17.add("Photon.Run2012A-PromptReco-v1.AOD.job44", '%s/zmeng//ICF/automated/2012_04_16_09_55_55/")'%pnfs, lumi = 1.0)
#photon17.add("Photon.Run2012A-PromptReco-v1.AOD.job57", '%s/zmeng//ICF/automated/2012_04_22_16_58_56/")'%pnfs, lumi = 1.0)
#photon17.add("Photon.Run2012A-PromptReco-v1.AOD.job69", '%s/zmeng//ICF/automated/2012_04_26_16_22_08/")'%pnfs, lumi = 1.0)
#photon17.add("Photon.Run2012A-PromptReco-v1.AOD.job74", '%s/zmeng//ICF/automated/2012_04_27_00_07_16/")'%pnfs, lumi = 1.0)
#photon17.add("Photon.Run2012A-PromptReco-v1.AOD.job81", '%s/clucas//ICF/automated/2012_05_06_13_23_30/")'%pnfs,lumi = 1.0)

a = "%s/yeshaq/ICF/automated/2012_05_25_16_29_39"%pnfs
photon17.add("Photon.2012A.job171",       '%s/Photon.Run2012A-PromptReco-v1.AOD")'%a,       lumi = 1.0)
photon17.add("SinglePhoton.2012B.job171", '%s/SinglePhoton.Run2012B-PromptReco-v1.AOD")'%a, lumi = 1.0)

photon17.add("GJets_HT400.job92",    '%s/clucas//ICF/automated/2012_05_08_11_07_51/")'%pnfs, xs = {"LO":107.5}["LO"])
photon17.add("GJets_HT400.job174",  '%s/yeshaq//ICF/automated/2012_05_30_22_41_51/")'%pnfs, xs = {"LO":107.5}["LO"])
