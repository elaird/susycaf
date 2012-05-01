from supy.samples import SampleHolder
from supy.sites import pnfs
pnfs = pnfs()
print pnfs
photon17 = SampleHolder()

photon17.add("Photon.Run2012A-PromptReco-v1.AOD.job29", '%s/zmeng//ICF/automated/2012_04_13_19_13_29/")'%pnfs, lumi = 1.0)
photon17.add("Photon.Run2012A-PromptReco-v1.AOD.job44", '%s/zmeng//ICF/automated/2012_04_16_09_55_55/")'%pnfs, lumi = 1.0)
photon17.add("Photon.Run2012A-PromptReco-v1.AOD.job57", '%s/zmeng//ICF/automated/2012_04_22_16_58_56/")'%pnfs, lumi = 1.0)
photon17.add("Photon.Run2012A-PromptReco-v1.AOD.job69", '%s/zmeng//ICF/automated/2012_04_26_16_22_08/")'%pnfs, lumi = 1.0)
photon17.add("Photon.Run2012A-PromptReco-v1.AOD.job74", '%s/zmeng//ICF/automated/2012_04_27_00_07_16/")'%pnfs, lumi = 1.0)
