from supy.samples import SampleHolder
from supy.sites import srm
electron16 = SampleHolder()

electron16.add("SingleEl.Run2011A", '%s/bbetchar//ICF/automated/2012_02_18_01_31_30/")'%srm, lumi = 9999 )
electron16.add("SingleEl.Run2011B", '%s/bbetchar//ICF/automated/2012_02_18_02_01_43/")'%srm, lumi = 9999 )
