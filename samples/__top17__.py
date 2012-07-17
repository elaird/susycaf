from supy.samples import SampleHolder
from supy.sites import srm
top16 = SampleHolder()

srm_burt = srm + '/bbetchar/ICF/automated'
singleT = srm_burt + '/2012_02_17_17_53_18/%s_TuneZ2_%s_7TeV-powheg-tauola.Fall11-PU_S6_START44_V9B-v1.AODSIM/")'

# Fall2011 reprocessing of Summer 2011
top16.add("ttj_mg", '%s/2012_02_14_20_57_14/")'%srm_burt, xs = {"LO":94.76, "guessNLO":157.5 }['guessNLO'] )
top16.add("ttj_ph", '%s/2012_02_24_06_31_07/")'%srm_burt , xs = 149.6 )
top16.add("top_s_ph", singleT%('T','s-channel'), xs = 2.341 )
top16.add("top_t_ph", singleT%('T','t-channel'), xs = 35.72 )
top16.add("top_tW_ph", singleT%('T','tW-channel-DS'), xs = 7.104 )
top16.add("tbar_s_ph", singleT%('Tbar','s-channel'), xs = 1.265 )
top16.add("tbar_t_ph", singleT%('Tbar','t-channel'), xs = 18.43 )
top16.add("tbar_tW_ph", singleT%('Tbar','tW-channel-DS'), xs = 7.108 )

