from supy.samples import SampleHolder
from supy.sites import pnfs
pnfs = pnfs()
photon17 = SampleHolder()

a = "%s/yeshaq/ICF/automated/2012_05_25_16_29_39"%pnfs
photon17.add("Photon.Run2012A-PromptReco-v1.AOD.job171", '%s/Photon.Run2012A-PromptReco-v1.AOD")'%a,       lumi = 1.0)
photon17.add("SinglePhoton.Run2012B-PromptReco-v1.AOD.job171", '%s/SinglePhoton.Run2012B-PromptReco-v1.AOD")'%a, lumi = 1.0)

photon17.add("Photon.Run2012A-PromptReco-v1.AOD.job228", '%s/zmeng//ICF/automated/2012_06_14_11_22_04/Photon.Run2012A-PromptReco-v1.AOD")'%pnfs, lumi = 1.0)
photon17.add("SinglePhoton.Run2012B-PromptReco-v1.AOD.job229", '%s/zmeng//ICF/automated/2012_06_14_14_26_12/SinglePhoton.Run2012B-PromptReco-v1.AOD")'%pnfs, lumi = 1.0)
photon17.add("SinglePhoton.Run2012B-PromptReco-v1.AOD.job238", '%s/zmeng//ICF/automated/2012_06_22_14_25_16/SinglePhoton.Run2012B-PromptReco-v1.AOD")'%pnfs, lumi = 1.0)

gJet = "GJets_HT-%s_8TeV-madgraph_v2.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM/"

photon17.add("g_jets_mg_ht_200_400.job501", '%s/agapitos//ICF/automated/2012_11_29_18_15_22/%s")' % (pnfs, gJet % "200To400"),
             xs={"LO": 960.5, "NLO":1140.78}["NLO"])
photon17.add("g_jets_mg_ht_400_inf.job501", '%s/agapitos//ICF/automated/2012_11_29_18_15_22/%s")' % (pnfs, gJet % "400ToInf"),
             xs={"LO": 107.5, "NLO":124.68}["NLO"])
