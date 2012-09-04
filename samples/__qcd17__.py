from supy.samples import SampleHolder
from supy.sites import pnfs
qcd17 = SampleHolder()
pnfs = pnfs()

binsXs = [ (   0, 4.9586188E10  ),
           (   5, 4.2639499E10  ),
           (  15, 9.8828742E8   ),
           (  30, 6.6285328E7   ),
           (  50, 8148778.0     ),
           (  80, 1033680.0     ),
           ( 120, 156293.3      ),
           ( 170, 34138.15      ),
           ( 300, 1759.549      ),
           ( 470, 113.8791      ),
           ( 600, 26.9921       ),
           ( 800, 3.550036      ),
           (1000, 0.737844      ),
           (1400, 0.03352235    ),
           (1800, 0.001829005   ),
           (None, None),
           ]

for i,(low,xs) in enumerate(binsXs[:-1]) :
    hi = binsXs[i+1][0]
    dir = "QCD_Pt-%d%s_TuneZ2star_8TeV_pythia6.Summer12-PU_S7_START52_V9-v1.AODSIM"%(low, "to%d"%hi if hi else "")
    qcd17.add("qcd_py6_pt_%d"%low, '%s/clucas//ICF/automated/2012_07_03_11_15_20/%s")'%(pnfs, dir), xs = xs)
