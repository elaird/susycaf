from supy.samples import SampleHolder
from supy.sites import srm
electron16 = SampleHolder()

electron16.add("SingleEl.Run2011A.1", '%s/bbetchar//ICF/automated/2012_02_09_06_16_46/")'%srm, lumi = 9999 )
electron16.add("SingleEl.Run2011A.2", '%s/bbetchar//ICF/automated/2012_02_09_06_22_00/")'%srm, lumi = 9999 )
electron16.add("SingleEl.Run2011B",   '%s/bbetchar//ICF/automated/2012_02_09_06_27_23/")'%srm, lumi = 9999 )
