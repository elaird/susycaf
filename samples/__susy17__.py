from supy.samples import SampleHolder
from supy.sites import pnfs
pnfs = pnfs()
susy17 = SampleHolder()

susy17.add("T2tt_8.job351", '%s/yeshaq//ICF/automated/2012_09_04_23_14_16/")'%pnfs, xs = 1.0) #dummy XS
