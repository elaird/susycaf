from supy.samples import SampleHolder
from supy.sites import pnfs
pnfs = pnfs()
susy17 = SampleHolder()

susy17.add("T2tt_8.job351", '%s/yeshaq//ICF/automated/2012_09_04_23_14_16/")'%pnfs, xs = 1.0) #dummy XS
susy17.add("T2tt_8.job351_1", '["/uscms/home/elaird/141_ntuples/job351/SusyCAF_Tree_1000_1_Bim.root"]', xs = 1.0) #dummy XS; [(500.0, 100.0), (775.0, 200.0)]
susy17.add("T2tt_500_100", '["/uscms/home/elaird/141_ntuples/job351/T2tt_8.job351_500_100_partial.root"]', xs = 0.0855) #/pb
susy17.add("T2tt_500_300", '["/uscms/home/elaird/141_ntuples/job351/T2tt_8.job351_500_300.root"]', xs = 0.0855) #/pb

