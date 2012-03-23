from supy.samples import SampleHolder
from supy.sites import srm
muon16 = SampleHolder()

muon16.add("SingleMu.2011A", '%s/bbetchar//ICF/automated/2012_02_29_14_08_41/", alwaysUseLastAttempt = True)'%srm, lumi = 2282.8165 ) #/pb
muon16.add("SingleMu.2011B", '%s/bbetchar//ICF/automated/2012_02_29_14_22_43/")'%srm, lumi = 2725.4681 ) #/pb
