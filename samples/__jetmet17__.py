from supy.samples import SampleHolder
from supy.sites import pnfs
jetmet17 = SampleHolder()
pnfs = pnfs()
#2012
        
jetmet17.add("JetHT.Run2012B-PromptReco-v1.AOD.job228",  '%s/zmeng//ICF/automated/2012_06_14_11_22_04/JetHT.Run2012B-PromptReco-v1.AOD/")'%pnfs, lumi = 1.0) #job 217,   226/226 comp.
jetmet17.add("JetHT.Run2012B-PromptReco-v1.AOD.job238",  '%s/zmeng//ICF/automated/2012_06_22_14_25_16/JetHT.Run2012B-PromptReco-v1.AOD/")'%pnfs, lumi = 1.0) #job 217,   226/226 comp.

