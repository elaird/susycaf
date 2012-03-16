from supy.samples import SampleHolder
from supy.sites import srm
ewk16 = SampleHolder()

srm_burt = srm+'/bbetchar/ICF/automated'

#Fall 11
# ewk16.add("wbb_mg", ) #pending
ewk16.add("wj_lv_mg", '%s/2012_02_29_14_33_59/")'%srm_burt, xs = {"LO":27770.0, "guessNLO": 37700}['guessNLO'] )
ewk16.add("dyj_ll_mg", '%s/2012_02_14_21_04_44/")'%srm_burt, xs = {"LO":2475.0, "guessNLO":None}['LO'] )

