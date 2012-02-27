from supy.samples import SampleHolder
from supy.sites import srm
muon16 = SampleHolder()

corrupt = ','.join('"SusyCAF_Tree_%d"'%d for d in [1125,290,339,428,579,606,742,721]).join(['[',']'])
muon16.add("SingleMu.2011A", '%s/bbetchar//ICF/automated/2012_02_23_02_37_07/", itemsToSkip=%s)'%(srm,corrupt), lumi = 2281.1693 ) # 2287.3459 including corrupt
muon16.add("SingleMu.2011B", '%s/bbetchar//ICF/automated/2012_02_22_00_30_32/")'%srm, lumi = 2738.0194 )
