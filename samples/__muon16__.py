from supy.samples import SampleHolder
from supy.sites import srm
muon16 = SampleHolder()

muon16.add("SingleMu.2011A.1", '%s/bbetchar//ICF/automated/2012_02_09_05_54_46/")'%srm, lumi = 9999 )
muon16.add("SingleMu.2011A.2", '%s/bbetchar//ICF/automated/2012_02_09_06_02_10/")'%srm, lumi = 9999 )
muon16.add("SingleMu.2011B",   '%s/bbetchar//ICF/automated/2012_02_09_06_09_26/")'%srm, lumi = 9999 )
