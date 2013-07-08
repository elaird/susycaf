from supy.samples import SampleHolder
from supy.sites import pnfs
import itertools
pnfs = pnfs()
top17 = SampleHolder()

top17.add("tt_8_mg.job315_1", '["/uscms/home/elaird/141_ntuples/job315/SusyCAF_Tree_100_1_h3D.root"]', xs = 234.)
top17.add("tt_8_mg.job315_zeroNu_50Met", '["/uscms/home/elaird/141_ntuples/job315/tt_8_mg.job315_zeroNu_50Met.root"]', xs = 234.*(44005./6923750))
top17.add("tt_8_mg.job315_1_oneNu_100Met", '["/uscms/home/elaird/141_ntuples/job315/tt_8_mg.job315_oneNu_100Met.root"]', xs = 234.*(1465./20000.))
top17.add("tt_8_mg.job315_1_twoNu_100Met", '["/uscms/home/elaird/141_ntuples/job315/tt_8_mg.job315_twoNu_100Met.root"]', xs = 234.*(510./20000.))

loc = "/uscms/home/elaird/141_ntuples/job315"
top17.add("tt_8_mg.job315_oneNu_100Met", 'utils.fileListFromDisk(location = "%s/oneNu_100Met/tt_8_mg.job315_*_skim.root", isDirectory = False)'%loc,
          xs = 7.037891e-02 * 2.340000e+02)
top17.add("tt_8_mg.job315_twoNu_100Met", 'utils.fileListFromDisk(location = "%s/twoNu_100Met/tt_8_mg.job315_*_skim.root", isDirectory = False)'%loc,
          xs = 2.542509e-02 * 2.340000e+02)
top17.add("ttz_8_mg.job269_1", '["/uscms/home/elaird/141_ntuples/job269/SusyCAF_Tree_11_1_7LX.root"]', xs = 0.172) #xs from PREP

top17.add("ttbar_powheg_v1.job410", '%s/clucas//ICF/automated/2012_09_26_14_02_42/")'%pnfs, xs = 211.)
top17.add("ttbar_powheg_v2.job404", '%s/zmeng//ICF/automated/2012_09_21_17_07_53/")'%pnfs, xs = 211.)

sTop = "-channel_TuneZ2star_8TeV-powheg-tauola.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM/"
combs = {"T_s":{"NLO":3.79, "LO":2.82},"T_t":{"NLO":56.4, "LO":47.0},"T_tW":{"NLO":11.1, "LO":10.7},
      "Tbar_s":{"NLO":1.76, "LO":1.57},"Tbar_t": {"NLO":30.7, "LO":25.0},"Tbar_tW":{"NLO":11.1, "LO":10.7}}

for key in combs: 
    top17.add("%s_powheg.job368"%key,
              '%s/zmeng//ICF/automated/2012_09_21_17_19_15/%s")'%(pnfs, (key + sTop).replace("tW-channel","tW-channel-DR")),combs[key]["NLO"])
