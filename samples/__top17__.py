from supy.samples import SampleHolder
from supy.sites import pnfs
top17 = SampleHolder()

top17.add("tt_8_mg.job188",  '%s/clucas//ICF/automated/2012_06_05_09_30_46")'%pnfs, xs = 234.) #from XS twiki; in PREP, xs = 491.3
top17.add("ttz_8_mg.job269", '%s/zmeng//ICF/automated/2012_07_20_15_24_50")'%pnfs, xs = 0.172) #xs from PREP
top17.add("ttz_8_mg.job269_1", '["/home/elaird/susycaf/SusyCAF_Tree_11_1_7LX.root"]', xs = 0.172) #xs from PREP
