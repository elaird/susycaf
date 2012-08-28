from supy.samples import SampleHolder
from supy.sites import pnfs
pnfs = pnfs()
ewk17 = SampleHolder()

incl = "WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball.Summer12-PU_S7_START52_V9-v1.AODSIM/"
ewk17.add("wj_lv_mg_ht_0_250",   '%s/dburton/WJet_Skim_v2/")'%pnfs, xs = 36140.2)
ewk17.add("wj_lv_mg_ht_250_300", '%s/clucas/ICF/automated/2012_06_20_16_29_08/")'%pnfs,           xs = 48.01)
ewk17.add("wj_lv_mg_ht_300_400", '%s/clucas/ICF/automated/2012_06_21_02_04_54/")'%pnfs,           xs = 38.30)
ewk17.add("wj_lv_mg_ht_400_inf", '%s/clucas/ICF/automated/2012_06_05_09_30_46/")'%pnfs,           xs = 25.22)
ewk17.add("wj_lv_mg_ht_incl",    '%s/clucas/ICF/automated/2012_06_05_09_30_46/%s")'%(pnfs, incl), xs = 36257.2)
