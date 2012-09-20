from supy.samples import SampleHolder
from supy.sites import pnfs
pnfs = pnfs()
ewk17 = SampleHolder()

#mgKFactor = 3048.0/2400.0 #Z+jets NNLO/LO
#mgKFactor2 = 3048.0/2475.0#Z+jets NNLO/LO

wj_lv_mgLoc = '/clucas/ICF/automated/'
wj_lvDset= "WJetsToLNu_HT-%s_8TeV-madgraph.Summer12-PU_S7_START52_V9-v1.AODSIM/"

250_300 = '2012_06_20_16_29_08/'
300_400 = '2012_06_21_02_04_54/'
400_inf = '2012_06_05_09_30_46/'

ewk17.add("wj_lv_mg_ht_250_300", '%s/%s/%s/%s")'%(pnfs, wj_lv_mgLoc, 250_300, wj_lvDset%"250To300"), xs = {"LO":7e+06, "fakeNLO":7e+06*mgKFactor}["fakeNLO"])
ewk17.add("wj_lv_mg_ht_300_400", '%s/%s/%s/%s")'%(pnfs, wj_lv_mgLoc, 300_400, wj_lvDset%"300To400"), xs = {"LO":171e+03, "fakeNLO":171e+03*mgKFactor}["fakeNLO"])
ewk17.add("wj_lv_mg_ht_400_inf", '%s/%s/%s/%s")'%(pnfs, wj_lv_mgLoc, 400_inf, wj_lvDset%"400ToInf"),xs = {"LO":5200, "fakeNLO":5200*mgKFactor}["fakeNLO"])

##### Incluive W+jet Sample ######

wj_lvInclDset = "WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball.Summer12-PU_S7_START52_V9-v1.AODSIM/"
incl =  '2012_06_05_09_30_46/'

ewk17.add("wj_lv_mg_ht_incl", '%s/%s/%s/%s")'%(pnfs, wj_lv_mgLoc, incl, wj_lvInclDset),xs = {"LO":30400.0, "NLO":36257.2})



