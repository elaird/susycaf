from supy.samples import SampleHolder
from supy.sites import pnfs
pnfs = pnfs()
top17 = SampleHolder()

top17.add("tt_8_mg.job188", '%s/clucas//ICF/automated/2012_06_05_09_30_46/TTJets_TuneZ2star_8TeV-madgraph-tauola.Summer12-PU_S7_START52_V9-v1.AODSIM")'%pnfs, xs = 234.)

top17.add("tt_8_mg.job315", '%s/clucas/ICF/automated/2012_08_17_14_42_02/")'%pnfs, xs = 234.)
top17.add("tt_8_mg.job315_1", '["/uscms/home/elaird/141_ntuples/job315/SusyCAF_Tree_100_1_h3D.root"]', xs = 234.)
top17.add("tt_8_mg.job315_zeroNu_50Met", '["/uscms/home/elaird/141_ntuples/job315/tt_8_mg.job315_zeroNu_50Met.root"]', xs = 234.*(44005./6923750))
top17.add("tt_8_mg.job315_1_oneNu_100Met", '["/uscms/home/elaird/141_ntuples/job315/tt_8_mg.job315_oneNu_100Met.root"]', xs = 234.*(1465./20000.))
top17.add("tt_8_mg.job315_1_twoNu_100Met", '["/uscms/home/elaird/141_ntuples/job315/tt_8_mg.job315_twoNu_100Met.root"]', xs = 234.*(510./20000.))

loc = "/uscms/home/elaird/141_ntuples/job315"
top17.add("tt_8_mg.job315_oneNu_100Met", 'utils.fileListFromDisk(location = "%s/oneNu_100Met/tt_8_mg.job315_*_skim.root", isDirectory = False)'%loc, xs = 7.037891e-02 * 2.340000e+02)
top17.add("tt_8_mg.job315_twoNu_100Met", 'utils.fileListFromDisk(location = "%s/twoNu_100Met/tt_8_mg.job315_*_skim.root", isDirectory = False)'%loc, xs = 2.542509e-02 * 2.340000e+02)

top17.add("ttz_8_mg.job269", '%s/zmeng//ICF/automated/2012_07_20_15_24_50")'%pnfs, xs = 0.172) #xs from PREP
top17.add("ttz_8_mg.job269_1", '["/uscms/home/elaird/141_ntuples/job269/SusyCAF_Tree_11_1_7LX.root"]', xs = 0.172) #xs from PREP


y = '%s/yeshaq//ICF/automated/2012_06_05_02_07_12/'%pnfs
top17.add("t_s_powheg.job200", '%s/clucas//ICF/automated/2012_06_06_18_49_32/")'%pnfs, xs = {"NNLO":5.55, "LO":2.82}["NNLO"])
top17.add("t_t_powheg.job187", '%s/T_t-channel_TuneZ2star_8TeV-powheg-tauola.Summer12-PU_S7_START52_V9-v1.AODSIM")'%y, xs = {"NNLO":56.4, "LO":47.0}["NNLO"])
top17.add("t_tw_powheg.job187", '%s/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola.Summer12-PU_S7_START52_V9-v1.AODSIM")'%y, xs = {"NNLO":11.1, "LO":10.7}["NNLO"])

top17.add("tbar_t_powheg.job194", '%s/clucas//ICF/automated/2012_06_06_11_53_47/")'%pnfs, xs = {"NNLO":30.7, "LO":25.0}["NNLO"])
top17.add("tbar_tw_powheg.job187", '%s/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola.Summer12-PU_S7_START52_V9-v1.AODSIM")'%y, xs = {"NNLO":11.1, "LO":10.7}["NNLO"])
