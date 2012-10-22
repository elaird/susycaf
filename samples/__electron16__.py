from supy.samples import SampleHolder
from supy.sites import srm
electron16 = SampleHolder()

electron16.add("SingleEl.2011A", '%s/bbetchar//ICF/automated/2012_03_02_08_42_57/")'%srm, lumi = 2311.5019 ) #/pb 
electron16.add("SingleEl.2011B", '%s/bbetchar//ICF/automated/2012_03_02_08_49_19/")'%srm, lumi = 2736.7709 ) #/pb 

electron16.add("EleHad.2011A", '%s/bbetchar//ICF/automated/2012_03_15_18_08_05/")'%srm, lumi = 2333 ) #/pb
electron16.add("EleHad.2011B", '%s/bbetchar//ICF/automated/2012_03_15_18_19_18/")'%srm, lumi = 2767 ) #/pb
