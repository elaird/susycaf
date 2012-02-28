from supy.samples import SampleHolder
from supy.sites import srm
ewk16 = SampleHolder()

srm_burt = srm+'/bbetchar/ICF/automated'

#Fall 11
# ewk16.add("wbb_mg", ) #pending
corrupt = ','.join('"SusyCAF_Tree_%d"'%d for d in [1287,1603,1286,1092,1237,1509,1415,1279,1385,1266,1150,1483,288]).join(['[',']'])
ewk16.add("wj_lv_mg", '%s/2012_02_18_01_38_16/", itemsToSkip = %s)'%(srm_burt,corrupt), xs = {"LO":27770.0, "guessNLO":0.9 * 55854}['guessNLO'] )
ewk16.add("dyj_ll_mg", '%s/2012_02_14_21_04_44/")'%srm_burt, xs = {"LO":2475.0, "guessNLO":None}['LO'] )
