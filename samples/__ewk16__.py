from supy.samples import SampleHolder
from supy.sites import srm
ewk16 = SampleHolder()

srm_burt = srm+'/bbetchar/ICF/automated'

#Fall 11
# ewk16.add("wbb_mg", ) #pending
ewk16.add("wj_lv_mg", '%s/2012_02_29_14_33_59/")'%srm_burt, xs = {"LO":27770.0, "guessNLO": None}['LO'] )
ewk16.add("dyj_ll_mg", '%s/2012_02_14_21_04_44/")'%srm_burt, xs = {"LO":2475.0, "guessNLO":None}['LO'] )

wnj = "W%dJets_TuneZ2_7TeV-madgraph-tauola.Fall11-PU_S6_START44_V9B-v%d.AODSIM"
ewk16.add("w2j_mg", '%s/2012_06_15_23_47_40/%s")'%(srm_burt,wnj%(2,1)), xs = 1435.0 )
ewk16.add("w3j_mg", '%s/2012_06_15_23_47_40/%s")'%(srm_burt,wnj%(3,2)), xs =  304.2 )
ewk16.add("w4j_mg", '%s/2012_06_15_23_47_40/%s")'%(srm_burt,wnj%(4,1)), xs =  172.6 )
