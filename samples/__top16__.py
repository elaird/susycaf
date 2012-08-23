from supy.samples import SampleHolder
from supy.sites import srm
top16 = SampleHolder()

srm_burt = srm + '/bbetchar/ICF/automated'
singleT = srm_burt + '/2012_02_17_17_53_18/%s_TuneZ2_%s_7TeV-powheg-tauola.Fall11-PU_S6_START44_V9B-v1.AODSIM/")'

skip = ['SusyCAF_Tree_255_1_bpq.root',
        'SusyCAF_Tree_19_1_I7H.root',
        'SusyCAF_Tree_2259_1_0sQ.root',
        'SusyCAF_Tree_2224_1_ZqC.root',
        'SusyCAF_Tree_1042_1_6ID.root',
        'SusyCAF_Tree_1531_1_LfW.root',
        'SusyCAF_Tree_2323_1_cym.root',
        'SusyCAF_Tree_858_1_kol.root',
        'SusyCAF_Tree_1028_1_5M6.root',
        'SusyCAF_Tree_2405_1_oOe.root',
        'SusyCAF_Tree_6_1_o3A.root',
        'SusyCAF_Tree_1944_1_0nc.root',
        'SusyCAF_Tree_2287_1_lGR.root',
        'SusyCAF_Tree_2146_1_UTy.root',
        ]

# Fall2011 reprocessing of Summer 2011
top16.add("ttj_mg", '%s/2012_02_14_20_57_14/", itemsToSkip=["%s"])'%(srm_burt,'","'.join(skip)), xs = {"LO":94.76, "guessNLO":157.5 }['guessNLO'] )
top16.add("ttj_ph", '%s/2012_02_24_06_31_07/")'%srm_burt, xs = 149.6 )
top16.add("ttj_mn", '%s/2012_08_01_00_43_36/")'%srm_burt, xs = 147.4 )
top16.add("top_s_ph", singleT%('T','s-channel'), xs = 2.341 )
top16.add("top_t_ph", singleT%('T','t-channel'), xs = 35.72 )
top16.add("top_tW_ph", singleT%('T','tW-channel-DS'), xs = 7.104 )
top16.add("tbar_s_ph", singleT%('Tbar','s-channel'), xs = 1.265 )
top16.add("tbar_t_ph", singleT%('Tbar','t-channel'), xs = 18.43 )
top16.add("tbar_tW_ph", singleT%('Tbar','tW-channel-DS'), xs = 7.108 )

