from supy.samples import SampleHolder
from supy.sites import srm
muon16 = SampleHolder()

muon16.add("SingleMu.2011A", '%s/bbetchar//ICF/automated/2012_02_29_14_08_41/", alwaysUseLastAttempt = True)'%srm, lumi = 2331 ) #/pb 
muon16.add("SingleMu.2011B", '%s/bbetchar//ICF/automated/2012_02_29_14_22_43/")'%srm, lumi = 2765 ) #/pb

muon16.add("MuHad.2011A", '%s/bbetchar//ICF/automated/2012_06_25_16_46_14/")'%srm, lumi = 2279.7997 ) #/pb 
muon16.add("MuHad.2011B", '%s/bbetchar//ICF/automated/2012_06_25_16_56_29/")'%srm, lumi = 2691.4278 ) #/pb 
