from supy.samples import SampleHolder
from supy.sites import pnfs, eos
qcd17 = SampleHolder()
pnfs = pnfs()
eos = eos()

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
    qcd17.add("qcd_py6_pt_%d_dCache"%low, '%s/clucas//ICF/automated/2012_07_03_11_15_20/%s")'%(pnfs, dir), xs = xs)

bEn = '%s/yeshaq//ICF/automated/2012_09_03_15_56_26/'%pnfs
qcd17.add("qcd_b_mg_ht_250", '%s/BJets_HT-250To500_8TeV-madgraph.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM")'%bEn, xs = 5828.0)
qcd17.add("qcd_b_mg_ht_500", '%s/BJets_HT-500To1000_8TeV-madgraph.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM")'%bEn, xs = 217.6)

qcd17.add("qcd_b_py6_pt_30", '%s/QCD_Pt-30To50_bEnriched_TuneZ2star_8TeV-pythia6-evtgen.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM")'%bEn, xs = 6.677E7)
qcd17.add("qcd_b_py6_pt_50", '%s/QCD_Pt-50To150_bEnriched_TuneZ2star_8TeV-pythia6-evtgen.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM")'%bEn, xs = 9355000.0)
qcd17.add("qcd_b_py6_pt_150", '%s/QCD_Pt-150_bEnriched_TuneZ2star_8TeV-pythia6-evtgen.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM")'%bEn, xs = 67340.0)

for i,(low,xs) in enumerate(binsXs[:-1]) :
    hi = binsXs[i+1][0]
    dir = "_.QCD_Pt-%d%s_TuneZ2star_8TeV_pythia6.Summer12-PU_S7_START52_V9-v1.AODSIM_/"%(low, "to%d"%hi if hi else "")
    if i < 1 : dir = dir.replace("_.QCD", "QCD").replace("AODSIM_","AODSIM")
    if low == 1800.: dir = dir.replace("AODSIM_","AODSIM") 
    qcd17.add("qcd_py6_pt_v2_%d"%low, '%s/agapitos//ICF/automated/2012_09_22_14_28_53/%s")'%(pnfs, dir), xs = xs)

for i,(low,xs) in enumerate(binsXs[4:-1]) :
    hi = binsXs[i+1][0]
    dir = "QCD_%d%s/"%(low, "to%d"%hi if hi else "")
    if low in (300,170) : dir = dir.replace("/","_merged/") 
    qcd17.add("qcd_py6_pt_%d"% low, '%s/clucas/Parked13/%s")'%(eos, dir), xs = xs)

