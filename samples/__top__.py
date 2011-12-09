from supy.samples import SampleHolder
from supy.sites import srm
top = SampleHolder()

srm_burt = srm + '/bbetchar/ICF/automated'
singleT = srm_burt + '/2011_12_08_08_39_30/%s_TuneZ2_%s_7TeV-powheg-tauola.Fall11-PU_S6_START42_V14B-v1.AODSIM/")'

# Fall2011 reprocessing of Summer 2011
top.add("ttj_mg", '%s/2011_12_08_08_31_06/")'%srm_burt , xs = {"LO":94.76, "guessNLO":157.5 }['guessNLO'] )
top.add("top_s_ph", singleT%('T','s-channel'), xs = 2.341 )
top.add("top_t_ph", singleT%('T','t-channel'), xs = 35.72 )
top.add("top_tW_ph", singleT%('T','tW-channel-DS'), xs = 7.104 )
top.add("tbar_s_ph", singleT%('Tbar','s-channel'), xs = 1.265 )
top.add("tbar_t_ph", singleT%('Tbar','t-channel'), xs = 18.43 )
top.add("tbar_tW_ph", singleT%('Tbar','tW-channel-DS'), xs = 7.108 )

# Summer2011
top.add("tt_tauola_fj_mg",'%s/2011_07_20_22_27_52/")'%srm_burt, xs = { "LO":94.76, "guessNLO":157.5}["guessNLO"] )
top.add("tt_tauola_fj", '%s/arlogb//ICF/automated/2011_07_11_17_17_07/")'%srm, xs = { "LO":94.76, "guessNLO":157.5}["guessNLO"])
