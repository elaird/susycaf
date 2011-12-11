from supy.samples import SampleHolder
from supy.sites import srm
ewk = SampleHolder()

srm_burt = srm+'/bbetchar/ICF/automated'

#Fall 11
# ewk.add("wbb_mg", ) #pending
ewk.add("wj_lv_mg", '%s/2011_12_08_08_51_04/")'%srm_burt, xs = {"LO":27770.0, "guessNLO":0.9 * 55854}['guessNLO'] )
ewk.add("dyj_ll_mg", '%s/2011_12_08_09_00_34/")'%srm_burt, xs = {"LO":2475.0, "guessNLO":None}['LO'] )


#Summer 11
ewk.add("w_jets_fj_mg", '%s/gouskos//ICF/automated/2011_07_18_17_43_04/")'%srm, xs = 0.9 * 55854)
