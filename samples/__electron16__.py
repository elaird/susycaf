from supy.samples import SampleHolder
from supy.sites import srm
electron16 = SampleHolder()

corrupt = ','.join('"SusyCAF_Tree_%d"'%d for d in [493]).join(['[',']'])
electron16.add("SingleEl.2011A", '%s/bbetchar//ICF/automated/2012_02_18_01_31_30/", itemsToSkip = %s)'%(srm,corrupt), lumi = 2287.1116 )
electron16.add("SingleEl.2011B", '%s/bbetchar//ICF/automated/2012_02_18_02_01_43/")'%srm, lumi = 2735.2943 )
