import supy,steps,calculables,samples
import os,math,copy,ROOT as r, numpy as np

class topAsymm(supy.analysis) :
    ''' Analysis for measurement of production asymmetry in top/anti-top pairs
    
    This analysis contains several secondary calculables, which need
    to be primed in the following order:
    
    1. Reweightings: run all samples in both tags, inverting label "selection";  --update
    2. Prime the b-tagging variable for [BQN]-type jets: top samples only, inverting label "top reco"; --update
    3. Prime the discriminants: [ top.tt, top.w_jet, QCD.SingleMu ] samples only, inverting after discriminants if at all; --update
    4. Run the analysis, all samples, no label inversion
    '''

    def parameters(self) :

        reweights = {
            'label' :['nvr',                'rrho',          'pu',        ],
            'abbr'  :['nvr',                'rrho',          'pu',        ],
            'func'  :['ratio_vertices',     'ratio_rho',     'pileup' ],
            'var'   :['nVertexRatio',       'rhoRatio',      'pileupTrueNumInteractionsBX0Target'],
            }

        objects = {
            'label'       :[               'pf' ,               'pat' ],
            'muonsInJets' :[               True ,               False ],
            'jet'         :[('xcak5JetPF','Pat'),   ('xcak5Jet','Pat')],
            'muon'        :[       ('muon','PF'),       ('muon','Pat')],
            'electron'    :[   ('electron','PF'),   ('electron','Pat')],
            'photon'      :[    ('photon','Pat'),     ('photon','Pat')],
            'met'         :[          'metP4PF' ,    'metP4AK5TypeII' ],
            'sumP4'       :[          'pfSumP4' ,           'xcSumP4' ],
            'sumPt'       :[       'metSumEtPF' ,           'xcSumPt' ],
            }

        leptons = {
            'name'     : [                  'muon',                'electron'],
            'ptMin'    : [                   20.0 ,                     30.0 ],
            'etaMax'   : [                    2.1 ,                      2.5 ],
            'iso'      : [    'CombinedRelativeIso',    'IsoCombinedAdjusted'],
            'triggers' : [       self.mutriggers(),         self.eltriggers()],
            'isoNormal': [             {"max":0.10},             {'max':0.09}],
            'isoInvert': [ {"min":0.15, "max":0.55}, {'min':0.14, 'max':0.55}]
            }

        #btag working points: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP
        csvWP = {"CSVL" : 0.244, "CSVM" : 0.679, "CSVT" : 0.898 }
        bCut = {"normal"   : {"index":0, "min":csvWP['CSVT']},
                "inverted" : {"index":0, "min":csvWP['CSVL'], "max":csvWP['CSVM']}}

        return { "vary" : ['selection','lepton','objects','reweights'],
                 "discriminant2DPlots": True,
                 "nJets" :  {"min":4,"max":None},
                 "unreliable": self.unreliableTriggers(),
                 "bVar" : "CSV", # "Combined Secondary Vertex"
                 "objects": self.vary([ ( objects['label'][index], dict((key,val[index]) for key,val in objects.iteritems())) for index in range(2) if objects['label'][index] in ['pf']]),
                 "lepton" : self.vary([ ( leptons['name'][index], dict((key,val[index]) for key,val in leptons.iteritems())) for index in range(2) if leptons['name'][index] in ['muon','electron'][:2]]),
                 "reweights" : self.vary([ ( reweights['label'][index], dict((key,val[index]) for key,val in reweights.iteritems())) for index in range(3) if reweights['label'][index] in ['pu']]),
                 "selection" : self.vary({"top" : {"bCut":bCut["normal"],  "iso":"isoNormal"},
                                          "QCD" : {"bCut":bCut["normal"],  "iso":"isoInvert"}
                                          #"Wlv" : {"bCut":bCut["inverted"],"iso":"isoNormal"},
                                          }),
                 "topBsamples": { "pythia"   : ("tt_tauola_fj",["tt_tauola_fj.wNonQQbar.tw.%s","tt_tauola_fj.wTopAsymP00.tw.%s"]),
                                  "madgraph" : ("ttj_mg",['ttj_mg.wNonQQbar.tw.%s','ttj_mg.wTopAsymP00.tw.%s']),
                                  }["madgraph"]
                 }

    @staticmethod
    def mutriggers() :                 # L1 prescaling evidence
        ptv = { 12   :(1,2,3,4,5),     # 7,8,11,12
                15   :(2,3,4,5,6,8,9), # 12,13
                20   :(1,2,3,4,5,7,8),
                24   :(1,2,3,4,5,7,8,11,12),
                24.21:(1,),
                30   :(1,2,3,4,5,7,8,11,12),
                30.21:(1,),
                40   :(1,2,3,5,6,9,10),
                40.21:(1,4,5),
                }
        return sum([[("HLT_Mu%d%s_v%d"%(int(pt),"_eta2p1" if type(pt)!=int else "",v),int(pt)+1) for v in vs] for pt,vs in sorted(ptv.iteritems())],[])
    
    @staticmethod
    def unreliableTriggers() :
        '''Evidence of L1 prescaling at these ostensible prescale values'''
        return { "HLT_Mu15_v9":(25,65),
                 "HLT_Mu20_v8":(30,),
                 "HLT_Mu24_v8":(20,25),
                 "HLT_Mu24_v11":(35,),
                 "HLT_Mu24_v12":(35,),
                 "HLT_Mu30_v8":(4,10),
                 "HLT_Mu30_v11":(4,20),
                 "HLT_Mu30_v12":(8,20),
                 "HLT_Mu40_v6":(4,),
                 "HLT_Mu40_v9":(10,),
                 "HLT_Mu40_v10":(10,)
                 }

    @staticmethod
    def eltriggers() :
        return zip(["%s_v%d"%(s,i) for s in
                    ["HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30",
                     "HLT_Ele25_CaloIdVT_TrkIdT_TriCentralJet30",
                     "HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30",
                     "HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30"]
                    for i in range(20) ], 20*4*[25])

    ########################################################################################

    def listOfSampleDictionaries(self) : return [getattr(samples,item) for item in ['muon16', 'electron16', 'top16', 'ewk16', 'qcd16']]

    def data(self,pars) :
        return { "muon" : supy.samples.specify( names = ['SingleMu.2011A',
                                                         'SingleMu.2011B'], weights = 'tw'),
                 "electron" : supy.samples.specify( names = ['EleHad.2011A',
                                                             'EleHad.2011B'], weights = 'tw')
                 }[pars['lepton']['name']]

    @staticmethod
    def single_top() :
        return ['top_s_ph','top_t_ph','top_tW_ph','tbar_s_ph','tbar_t_ph','tbar_tW_ph']

    def listOfSamples(self,pars) :
        rw = pars['reweights']['abbr']
        lname = pars['lepton']["name"]

        def qcd_py6_mu(eL = None) :
            if True or lname == "electron" : return []
            return supy.samples.specify( names = ["qcd_mu_15_20",
                                                  "qcd_mu_20_30",
                                                  "qcd_mu_30_50",
                                                  "qcd_mu_50_80",
                                                  "qcd_mu_80_120",
                                                  "qcd_mu_120_150",
                                                  "qcd_mu_150"],  weights = ['tw',rw], effectiveLumi = eL )     if 'Wlv' not in pars['tag'] else []
        def ewk(eL = None) :
            return supy.samples.specify( names = [#"wj_lv_mg",
                                                  "dyj_ll_mg",
                                                  "w2j_mg",
                                                  "w3j_mg",
                                                  "w4j_mg"], effectiveLumi = eL, color = 28, weights = ['tw',rw] ) if "QCD" not in pars['tag'] else []

        def single_top(eL = None) :
            return supy.samples.specify( names = self.single_top(),
                                         effectiveLumi = eL, color = r.kGray, weights = ['tw',rw]) if "QCD" not in pars['tag'] else []
        
        def ttbar_mg(eL = None) :
            return (supy.samples.specify(names = "ttj_mg", effectiveLumi = eL, color = r.kBlue, weights = ["wNonQQbar",'tw',rw] ) +
                    sum( [supy.samples.specify( names = "ttj_mg", effectiveLumi = eL, 
                                                color = color, weights = [ calculables.top.wTopAsym(asym, R_sm = -0.05), 'tw',rw ] )
                          for asym,color in [(0.0,r.kOrange),
                                             (-0.3,r.kGreen),(0.3,r.kRed),
                                             #(-0.6,r.kYellow),(0.6,r.kYellow),
                                             #(-0.5,r.kYellow),(0.5,r.kYellow),
                                             #(-0.4,r.kYellow),(0.4,r.kYellow),
                                             #(-0.2,r.kYellow),
                                             (0.2,r.kYellow),
                                             #(-0.1,r.kYellow),
                                             (0.1,r.kYellow),
                                             ]], [])
                    )[: 2 if "QCD" in pars['tag'] else 2 if 'Wlv' in pars['tag'] else None]
        
        return  ( self.data(pars) + qcd_py6_mu() + ewk() + ttbar_mg(8e4) + single_top() )


    ########################################################################################
    def listOfCalculables(self, pars) :
        obj = pars["objects"]
        lepton = obj[pars["lepton"]["name"]]
        calcs = sum( [supy.calculables.zeroArgs(module) for module in [calculables, supy.calculables]]+
                     [supy.calculables.fromCollections(getattr(calculables,item), [obj[item]])
                      for item in  ['jet','photon','electron','muon']], [])
        calcs += supy.calculables.fromCollections(calculables.top,[('genTop',""),('fitTop',"")])
        calcs += [
            calculables.jet.IndicesBtagged(obj["jet"],pars["bVar"]),

            calculables.jet.Indices(       obj["jet"],      ptMin = 20, etaMax = 3.1, flagName = "JetIDloose"),
            calculables.jet.Indices(       obj["jet"],      ptMin = 25, etaMax = 2.4, flagName = "JetIDloose", extraName = "triggering"),

            calculables.electron.Indices_TopPAG(  obj["electron"], ptMin = 30, absEtaMax = 2.5, id = "ID70"),
            calculables.muon.Indices(             obj["muon"],     ptMin = 20, absEtaMax = 2.1, ID = "ID_TOPPAG",

                                                  isoMax = 0.25, ISO = "CombinedRelativeIso"), #these two are kind of irrelevant, since we use IndicesAnyIsoIsoOrder

            calculables.electron.IndicesIsoLoose( obj["electron"], ptMin = 15, absEtaMax = 2.5, iso = "IsoCombinedAdjusted", isoMax = 0.15),
            calculables.muon.IndicesIsoLoose( obj["muon"], ptMin = 10, absEtaMax = 2.5, iso = "CombinedRelativeIso", isoMax = 0.20 ),

            calculables.electron.IsoCombinedAdjusted(obj["electron"], barrelCIso = 0.09, endcapCIso = 0.06 ), # VBTF 0.85
            calculables.muon.IndicesTriggering(lepton),
            calculables.muon.IndicesAnyIsoIsoOrder(lepton, pars["lepton"]["iso"]),
            calculables.muon.MinJetDR(lepton, obj["jet"], jetPtMin = 30),
            calculables.xclean.xcJet_SingleLepton( obj["jet"], leptons = lepton, indices = "IndicesAnyIsoIsoOrder" ),

            calculables.trigger.lowestUnPrescaledTrigger(zip(*pars["lepton"]["triggers"])[0]),

            calculables.vertex.ID(),
            calculables.vertex.Indices(),

            calculables.top.TopJets( obj['jet'] ),
            calculables.top.TopLeptons( lepton ),
            calculables.top.mixedSumP4( transverse = obj["met"], longitudinal = obj["sumP4"] ),
            calculables.top.TopComboQQBBLikelihood( pars['bVar'] ),
            calculables.top.OtherJetsLikelihood( pars['bVar'] ),
            calculables.top.TopRatherThanWProbability( priorTop = 0.5 ),
            calculables.top.IndicesGenTopPQHL( obj['jet'] ),
            calculables.top.IndicesGenTopExtra (obj['jet'] ),
            calculables.top.genTopRecoIndex(),
            calculables.top.TopReconstruction(),
            calculables.top.TTbarSignExpectation(nSamples = 16, qDirFunc = "qDirExpectation_EtaSum"),

            calculables.other.Mt( lepton, "mixedSumP4", allowNonIso = True, isSumP4 = True),
            calculables.other.Covariance(('met','PF')),
            calculables.other.TriDiscriminant( LR = "DiscriminantWQCD", LC = "DiscriminantTopW", RC = "DiscriminantTopQCD"),
            calculables.other.KarlsruheDiscriminant( obj['jet'], obj['met'] ),

            calculables.jet.pt( obj['jet'], index = 0, Btagged = True ),
            calculables.jet.absEta( obj['jet'], index = 3, Btagged = False),

            supy.calculables.other.pt( "mixedSumP4" ),
            supy.calculables.other.size( "Indices".join(obj['jet']) ),
            supy.calculables.other.abbreviation( "TrkCountingHighEffBJetTags", "TCHE", fixes = calculables.jet.xcStrip(obj['jet']) ),
            supy.calculables.other.abbreviation( "CombinedSecondaryVertexBJetTags", "CSV", fixes = calculables.jet.xcStrip(obj['jet']) ),
            supy.calculables.other.abbreviation( pars['reweights']['var'], pars['reweights']['abbr'] ),
            supy.calculables.other.abbreviation( {'muon':'muonTriggerWeightPF','electron':"CrossTriggerWeight"}[pars['lepton']['name']], 'tw' ),
            supy.calculables.other.abbreviation( "xcak5JetPFCorrectedP4Pat","xcak5JetPFCorrectedP4Pattriggering"),
            ]
        return calcs
    ########################################################################################

    def listOfSteps(self, pars) :
        obj = pars["objects"]
        lname = pars["lepton"]["name"]
        lepton = obj[lname]
        otherLepton = obj[{'electron':'muon','muon':'electron'}[lname]]
        lPtMin = pars["lepton"]["ptMin"]
        lIndices = 'IndicesAnyIsoIsoOrder'.join(lepton)
        lIso = pars['lepton']['iso'].join(lepton)
        lIsoMinMax = pars["lepton"][pars['selection']['iso']]
        topTag = pars['tag'].replace("Wlv","top").replace("QCD","top")
        bVar = pars["bVar"].join(obj["jet"])[2:]
        rw = pars['reweights']['abbr']
        
        ssteps = supy.steps
        
        return (
            [ssteps.printer.progressPrinter()
             , ssteps.histos.value("genpthat",200,0,1000,xtitle="#hat{p_{T}} (GeV)").onlySim()
             , ssteps.histos.value("genQ",200,0,1000,xtitle="#hat{Q} (GeV)").onlySim()
             ####################################
             , ssteps.filters.label('data cleanup'),
             ssteps.filters.multiplicity("vertexIndices",min=1),
             ssteps.filters.value('physicsDeclared',min=1).onlyData(),
             ssteps.filters.value('trackingFailureFilterFlag',min=1),
             ssteps.filters.value('beamHaloCSCTightHaloId', max=0),
             ssteps.filters.value('hbheNoiseFilterResult',min=1),
             ssteps.filters.value('hcalLaserEventFilterFlag', min=1).onlyData(),
             steps.trigger.l1Filter("L1Tech_BPTX_plus_AND_minus.v0").onlyData(),
             steps.filters.monster()

             ####################################

             , ssteps.filters.label('reweighting')
             , getattr(self,pars['reweights']['func'])(pars)
             , ssteps.histos.value(obj["sumPt"],50,0,1500)
             , ssteps.histos.value("rho",100,0,40)

             , ssteps.filters.label('trigger reweighting'),
             self.triggerWeight(pars, [ss.weightedName for ss in self.data(pars)])
             , steps.trigger.lowestUnPrescaledTriggerHistogrammer().onlyData()
             
             ####################################
             , ssteps.filters.label('selection'),

             ssteps.filters.multiplicity( min=1, var = "IndicesAnyIso".join(lepton) )
             , ssteps.histos.multiplicity("Indices".join(obj["jet"])),
             ssteps.filters.multiplicity("Indices".join(obj["jet"]), **pars["nJets"])
             
             , ssteps.histos.value( lIso, 100, 0, 1, indices = "IndicesIsoLoose".join(lepton), index = 1),
             ssteps.filters.multiplicity( max=1, var = "IndicesIsoLoose".join(lepton) ),
             ssteps.filters.multiplicity( max=0, var = "IndicesIsoLoose".join(otherLepton) )
             
             , ssteps.histos.absEta("P4".join(lepton), 100,0,4, indices = lIndices, index = 0)
             , ssteps.histos.pt("P4".join(lepton), 200,0,200, indices = lIndices, index = 0)
             , ssteps.histos.value("MinJetDR".join(lepton), 100,0,2, indices = lIndices, index = 0),
             ssteps.filters.value("MinJetDR".join(lepton), min = 0.3, indices = lIndices, index = 0)
 
             , ssteps.histos.pt("mixedSumP4",100,0,300),
             ssteps.filters.value('ecalDeadCellTPFilterFlag',min=1),
             steps.jet.failedJetVeto( obj["jet"], ptMin = 20, id = "PFJetIDloose")
             , ssteps.histos.pt("mixedSumP4",100,0,300)

             , ssteps.histos.value( lIso, 55,0,1.1, indices = lIndices, index=0),
             ssteps.filters.value( lIso, indices = lIndices, index = 0, **lIsoMinMax)

             , calculables.jet.ProbabilityGivenBQN(obj["jet"], pars['bVar'], binning=(51,-0.02,1), samples = (pars['topBsamples'][0],[s%rw for s in pars['topBsamples'][1]]), tag = topTag)
             , ssteps.histos.value("TopRatherThanWProbability", 100,0,1)
             , ssteps.histos.value(bVar, 51,-0.02,1, indices = "IndicesBtagged".join(obj["jet"]), index = 0)
             , ssteps.histos.value(bVar, 51,-0.02,1, indices = "IndicesBtagged".join(obj["jet"]), index = 1)
             , ssteps.histos.value(bVar, 51,-0.02,1, indices = "IndicesBtagged".join(obj["jet"]), index = 2),
             ssteps.filters.value(bVar, indices = "IndicesBtagged".join(obj["jet"]), **pars["selection"]["bCut"])
             
             , ssteps.filters.label('top reco'),
             ssteps.filters.multiplicity("TopReconstruction",min=1)
             , ssteps.filters.label("selection complete")

             , steps.top.channelClassification().onlySim()
             , steps.top.combinatorialFrequency().onlySim()
             #, steps.displayer.ttbar(jets=obj["jet"], met=obj['met'], muons = obj['muon'], electrons = obj['electron'])
             ####################################
             , ssteps.filters.label('discriminants')
             , ssteps.histos.value("KarlsruheDiscriminant", 28, -320, 800 )
             , ssteps.histos.value("TriDiscriminant",50,-1,1)
             , ssteps.filters.label('qq:gg'),   self.discriminantQQgg(pars)
             , ssteps.filters.label('qq:gg:4j'),self.discriminantQQgg4Jet(pars)
             , ssteps.filters.label('top:W'),   self.discriminantTopW(pars)
             , ssteps.filters.label('top:QCD'), self.discriminantTopQCD(pars)
             , ssteps.filters.label('W:QCD'),   self.discriminantWQCD(pars)
             , calculables.gen.qDirExpectation_('fitTopSumP4Eta', 8, 'top_muon_pf_%s'%rw, 'ttj_mg.wTopAsymP00.tw.%s'%rw)
             , calculables.gen.qDirExpectation_SumRapidities('top_muon_pf_%s'%rw, 'ttj_mg.wTopAsymP00.tw.%s'%rw)
             , calculables.gen.qDirExpectation_EtaSum('top_muon_pf_%s'%rw, 'ttj_mg.wTopAsymP00.tw.%s'%rw)
             , calculables.gen.qDirExpectation_RapiditySum('top_muon_pf_%s'%rw, 'ttj_mg.wTopAsymP00.tw.%s'%rw)
             ####################################
             , ssteps.filters.label('gen top kinfit ')
             , steps.top.kinFitLook('genTopRecoIndex')
             #, steps.top.kinematics('genTop')
             , steps.top.resolutions('genTopRecoIndex')
             , ssteps.filters.label('reco top kinfit ')
             , steps.top.kinFitLook('fitTopRecoIndex')
             , steps.top.kinematics('fitTop')
             , steps.top.resolutions('fitTopRecoIndex')
             ####################################

             , ssteps.histos.value("M3".join(obj['jet']), 20,0,800)
             , ssteps.histos.multiplicity("Indices".join(obj["jet"]))
             , ssteps.filters.label('object pt')
             , ssteps.histos.pt("mixedSumP4",100,1,201)
             , ssteps.histos.pt("P4".join(lepton), 100,1,201, indices = lIndices, index = 0)
             , ssteps.histos.pt("CorrectedP4".join(obj['jet']), 100,1,201, indices = "Indices".join(obj['jet']), index = 0)
             , ssteps.histos.pt("CorrectedP4".join(obj['jet']), 100,1,201, indices = "Indices".join(obj['jet']), index = 1)
             , ssteps.histos.pt("CorrectedP4".join(obj['jet']), 100,1,201, indices = "Indices".join(obj['jet']), index = 2)
             , ssteps.histos.pt("CorrectedP4".join(obj['jet']), 100,1,201, indices = "Indices".join(obj['jet']), index = 3)
             , ssteps.filters.label('object eta')
             , ssteps.histos.absEta("P4".join(lepton), 100,0,4, indices = lIndices, index = 0)
             , ssteps.histos.absEta("CorrectedP4".join(obj['jet']), 100,0,4, indices = "Indices".join(obj['jet']), index = 0)
             , ssteps.histos.absEta("CorrectedP4".join(obj['jet']), 100,0,4, indices = "Indices".join(obj['jet']), index = 1)
             , ssteps.histos.absEta("CorrectedP4".join(obj['jet']), 100,0,4, indices = "Indices".join(obj['jet']), index = 2)
             , ssteps.histos.absEta("CorrectedP4".join(obj['jet']), 100,0,4, indices = "Indices".join(obj['jet']), index = 3)
             
             , ssteps.filters.label('extended jets')
             , ssteps.histos.value('FourJetPtThreshold'.join(obj['jet']), 50,0,100)
             , ssteps.histos.value('fitTopJetPtMin', 50,0,100)
             , ssteps.histos.value('FourJetAbsEtaThreshold'.join(obj['jet']), 40,0,4)
             , ssteps.histos.value('fitTopJetAbsEtaMax', 40,0,4)
             
             , ssteps.filters.label('signal distributions'), steps.top.Asymmetry(('fitTop',''), bins = 41)#640)
             , ssteps.filters.label('spin distributions'),   steps.top.Spin(('fitTop',''))
             , ssteps.filters.label('sign check'),           steps.top.signCheck()

             #steps.histos.value('fitTopSumP4Eta', 12, -6, 6),
             #steps.filters.absEta('fitTopSumP4', min = 1),
             #steps.histos.value('fitTopSumP4Eta', 12, -6, 6),
             #steps.top.Asymmetry(('fitTop',''), bins = 640),
             #steps.top.Spin(('fitTop','')),
             
             #steps.top.kinFitLook("fitTopRecoIndex"),
             #steps.filters.value("TriDiscriminant",min=-0.64,max=0.8),
             #steps.histos.value("TriDiscriminant",50,-1,1),
             #steps.top.Asymmetry(('fitTop',''), bins = 640),
             #steps.top.Spin(('fitTop','')),
             #steps.filters.value("TriDiscriminant",min=-.56,max=0.72),
             #steps.histos.value("TriDiscriminant",50,-1,1),
             #steps.top.Asymmetry(('fitTop','')),
             #steps.top.Spin(('fitTop','')),
             ])
    ########################################################################################

    @staticmethod
    def lepIso(index,pars) :
        lepton = pars["objects"][pars["lepton"]["name"]]
        return 

    @staticmethod
    def triggerWeight(pars,samples) :
        triggers,thresholds = zip(*pars['lepton']['triggers'])
        jets = pars['objects']['jet']
        return { 'electron' : calculables.trigger.CrossTriggerWeight( samples = samples,
                                                                      triggers = triggers,
                                                                      jets = (jets[0],jets[1]+'triggering')),
                 'muon'     : calculables.trigger.TriggerWeight( samples,
                                                                 unreliable = pars['unreliable'],
                                                                 triggers = triggers,
                                                                 thresholds = thresholds )
                 }[pars['lepton']['name']]
    
    @classmethod
    def pileup(cls,pars) :
        rw = pars['reweights']['abbr']
        return supy.calculables.other.Target("pileupTrueNumInteractionsBX0", thisSample = pars['baseSample'],
                                             target = ("data/pileup_true_Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON.root","pileup"),
                                             groups = [('qcd_mu',[]),('wj_lv_mg',[]),('dyj_ll_mg',[]),
                                                       ("w2j_mg",[]),("w3j_mg",[]),("w4j_mg",[]),
                                                       ('single_top', ['%s.tw.%s'%(s,rw) for s in cls.single_top()]),
                                                       ('ttj_mg',['ttj_mg%s.tw.%s'%(s,rw) for s in ['',
                                                                                                    '.wNonQQbar',
                                                                                                    '.wTopAsymP00']])]).onlySim()
    @classmethod
    def ratio_vertices(cls,pars) :
        rw = pars['reweights']['abbr']
        return supy.calculables.other.Ratio("nVertex", binning = (20,-0.5,19.5), thisSample = pars['baseSample'],
                                            target = ("SingleMu",[]), groups = [('qcd_mu',[]),('wj_lv_mg',[]),('dyj_ll_mg',[]),
                                                                                ("w2j_mg",[]),("w3j_mg",[]),("w4j_mg",[]),
                                                                                ('single_top', ['%s.tw.%s'%(s,rw) for s in cls.single_top()]),
                                                                                ('ttj_mg',['ttj_mg%s.tw.'%(s,rw) for s in ['',
                                                                                                                           '.wNonQQbar',
                                                                                                                           '.wTopAsymP00']])])
    @classmethod
    def ratio_rho(cls,pars) :
        rw = pars['reweights']['abbr']
        return supy.calculables.other.Ratio("rho", binning = (90,0,30), thisSample = pars['baseSample'],
                                            target = ("SingleMu",[]), groups = [('qcd_mu',[]),('wj_lv_mg',[]),('dyj_ll_mg',[]),
                                                                                ("w2j_mg",[]),("w3j_mg",[]),("w4j_mg",[]),
                                                                                ('single_top', ['%s.tw.s'%(s,rw) for s in cls.single_top()]),
                                                                                ('ttj_mg',['ttj_mg%s.tw.%s'%(s,rw) for s in ['',
                                                                                                                             '.wNonQQbar',
                                                                                                                             '.wTopAsymP00']])])
    @staticmethod
    def discriminantQQgg(pars) :
        rw = pars['reweights']['abbr']
        lname = pars['lepton']['name']
        return supy.calculables.other.Discriminant( fixes = ("","QQgg"),
                                                    left = {"pre":"gg", "tag":"top_%s_pf_%s"%(lname,rw), "samples":['ttj_mg.wNonQQbar.tw.%s'%rw]},
                                                    right = {"pre":"qq", "tag":"top_%s_pf_%s"%(lname,rw), "samples":['ttj_mg.wTopAsymP00.tw.%s'%rw]},
                                                    dists = {"M3".join(pars["objects"]['jet']) : (20,0,600), # 0.047
                                                             "fitTopNtracksExtra" : (20,0,160),             
                                                             "tracksCountwithPrimaryHighPurityTracks" :  (20,0,250),               # 0.049
                                                             "fitTopFifthJet": (2,-0.5, 1.5),              # 0.052
                                                             "fitTopAbsSumRapidities" : (20, 0, 4),
                                                             #"fitTopPartonXlo" : (20,0,0.12),              # 0.036
                                                             #"fitTopBMomentsSum2" : (20,0,0.2),           # 0.004
                                                             #"fitTopPartonXhi" : (20,0.04,0.4),           # 0.003
                                                             },
                                                    correlations = pars['discriminant2DPlots'],
                                                    bins = 14)
    @staticmethod
    def discriminantQQgg4Jet(pars) :
        rw = pars['reweights']['abbr']
        obj = pars['objects']
        lname = pars['lepton']['name']
        return supy.calculables.other.Discriminant( fixes = ("","QQgg4Jet"),
                                                    left = {"pre":"gg", "tag":"top_%s_pf_%s"%(lname,rw), "samples":['ttj_mg.wNonQQbar.tw.%s'%rw]},
                                                    right = {"pre":"qq", "tag":"top_%s_pf_%s"%(lname,rw), "samples":['ttj_mg.wTopAsymP00.tw.%s'%rw]},
                                                    dists = {'FourJetPtThreshold'.join(obj['jet']) : (25,0,100),
                                                             'FourJetAbsEtaThreshold'.join(obj['jet']) : (20,0,4),
                                                             },
                                                    correlations = pars['discriminant2DPlots'],
                                                    bins = 14)
    @staticmethod
    def discriminantTopW(pars) :
        rw = pars['reweights']['abbr']
        lname = pars['lepton']['name']
        return supy.calculables.other.Discriminant( fixes = ("","TopW"),
                                                    left = {"pre":"wj", "tag":"top_%s_pf_%s"%(lname,rw), "samples":['w%dj_mg.tw.%s'%(n,rw) for n in [2,3,4]]},
                                                    right = {"pre":"ttj_mg", "tag":"top_%s_pf_%s"%(lname,rw), "samples": ['ttj_mg.%s.tw.%s'%(s,rw) for s in ['wNonQQbar','wTopAsymP00']]},
                                                    correlations = pars['discriminant2DPlots'],
                                                    dists = {"TopRatherThanWProbability" : (20,0.5,1),          # 0.185
                                                             "B0pt".join(pars['objects']["jet"]) : (20,0,200),  # 0.128
                                                             "fitTopHadChi2"     : (20,0,20),                   # 0.049
                                                             #"mixedSumP4.pt"     : (30,0,180),                 # 0.043
                                                             #"Kt".join(pars['objects']["jet"]) : (25,0,150),   # 0.018
                                                             #"fitTopDeltaPhiLNu" : (20,0,math.pi),             # 0.019
                                                             #"fitTopLeptonPt"    : (20,0,180),                 # Find out
                                                             })
    @staticmethod
    def discriminantTopQCD(pars) :
        rw = pars['reweights']['abbr']
        lname = pars['lepton']['name']
        tops = ['ttj_mg.%s.tw.%s'%(s,rw) for s in ['wNonQQbar','wTopAsymP00']]
        datas = {"muon" : ["SingleMu.2011A.tw","SingltMu.2011B.tw"],
                 "electron": ["EleHad.2011A.tw","EleHad.2011B.tw"]}[lname]
        lumi = 5008 # FIXME HACK !!!
        return supy.calculables.other.Discriminant( fixes = ("","TopQCD"),
                                                    left = {"pre":"Multijet", "tag":"QCD_%s_pf_%s"%(lname,rw), "samples":datas+tops, 'sf':[1,1,-lumi,-lumi]},
                                                    right = {"pre":"ttj_mg", "tag":"top_%s_pf_%s"%(lname,rw), "samples": tops},
                                                    correlations = pars['discriminant2DPlots'],
                                                    dists = {"Mt".join(pars['objects'][lname])+"mixedSumP4" : (10,0,100), # 0.329
                                                             "Kt".join(pars['objects']["jet"]) : (12,0,120),              # 0.079
                                                             #"B0pt".join(pars['objects']["jet"]) : (30,0,300),
                                                             })
    @staticmethod
    def discriminantWQCD(pars) :
        rw = pars['reweights']['abbr']
        lname = pars['lepton']['name']
        tops = ['ttj_mg.%s.tw.%s'%(s,rw) for s in ['wNonQQbar','wTopAsymP00']]
        datas = {"muon" : ["SingleMu.2011A.tw","SingltMu.2011B.tw"],
                 "electron": ["EleHad.2011A.tw","EleHad.2011B.tw"]}[lname]
        lumi = 5008 # FIXME HACK !!!
        return supy.calculables.other.Discriminant( fixes = ("","WQCD"),
                                                    left = {"pre":"wj", "tag":"top_%s_pf_%s"%(lname,rw), "samples":['w%dj_mg.tw.%s'%(n,rw) for n in [2,3,4]]},
                                                    right = {"pre":"Multijet", "tag":"QCD_%s_pf_%s"%(lname,rw), "samples":datas+tops, 'sf':[1,1,-lumi,-lumi]},
                                                    correlations = pars['discriminant2DPlots'],
                                                    dists = {"Mt".join(pars['objects'][lname])+"mixedSumP4" : (10,0,100), # 0.395
                                                             "B0pt".join(pars['objects']["jet"]) : (10,0,200),            # 0.122
                                                             })
    
    ########################################################################################
    def concludeAll(self) :
        self.rowcolors = 2*[13] + 2*[45]
        super(topAsymm,self).concludeAll()
        #self.meldWpartitions()
        #self.meldQCDpartitions()
        for rw,lname in set([(pars['reweights']['abbr'],pars['lepton']['name']) for pars in self.readyConfs]) :
            self.meldScale(rw,lname)
            self.plotMeldScale(rw,lname)
            #self.ensembleTest(rw,lname)
        #self.PEcurves()

    def conclude(self,pars) :
        rw = pars['reweights']['abbr']
        org = self.organizer(pars, verbose = True )
        org.mergeSamples(targetSpec = {"name":"SingleMu.2011", "color":r.kBlack, "markerStyle":20}, allWithPrefix="SingleMu")
        org.mergeSamples(targetSpec = {"name":"EleHad.2011", "color":r.kBlack, "markerStyle":20}, allWithPrefix="EleHad")
        org.mergeSamples(targetSpec = {"name":"multijet", "color":r.kBlue}, allWithPrefix="qcd_mu")
        org.mergeSamples(targetSpec = {"name":"t#bar{t}", "color":r.kViolet}, sources=["ttj_mg.wNonQQbar.tw.%s"%rw,"ttj_mg.wTopAsymP00.tw.%s"%rw], keepSources = True)
        org.mergeSamples(targetSpec = {"name":"t#bar{t}.q#bar{q}.N30", "color":r.kRed}, sources = ["ttj_mg.wTopAsymN30.tw.%s"%rw,"ttj_mg.wNonQQbar.tw.%s"%rw][:1])
        org.mergeSamples(targetSpec = {"name":"t#bar{t}.q#bar{q}.P30", "color":r.kGreen}, sources = ["ttj_mg.wTopAsymP30.tw.%s"%rw,"ttj_mg.wNonQQbar.tw.%s"%rw][:1])
        #org.mergeSamples(targetSpec = {"name":"W+jets", "color":28}, allWithPrefix="wj_lv_mg")
        org.mergeSamples(targetSpec = {"name":"W+jets", "color":28}, sources = ["w%dj_mg.tw.%s"%(n,rw) for n in [2,3,4]])
        org.mergeSamples(targetSpec = {"name":"DY+jets", "color":r.kYellow}, allWithPrefix="dyj_ll_mg")
        org.mergeSamples(targetSpec = {"name":"Single top", "color":r.kGray}, sources = ["%s.tw.%s"%(s,rw) for s in self.single_top()])
        org.mergeSamples(targetSpec = {"name":"Standard Model", "color":r.kGreen+2}, sources = ["multijet","t#bar{t}","W+jets","DY+jets","Single top"], keepSources = True)
        #for ss in filter(lambda ss: 'ttj_mg' in ss['name'],org.samples) : org.drop(ss['name'])

        orgpdf = copy.deepcopy(org)
        orgpdf.scale( toPdf = True )
        org.scale( lumiToUseInAbsenceOfData = 1.1e3 )

        names = [ss["name"] for ss in org.samples]
        kwargs = {"detailedCalculables": False,
                  "blackList":["lumiHisto","xsHisto","nJobsHisto"],
                  "samplesForRatios" : next(iter(filter(lambda x: x[0] in names and x[1] in names, [("SingleMu.2011","Standard Model"),
                                                                                                    ("EleHad.2011","Standard Model")])), ("","")),
                  "sampleLabelsForRatios" : ("data","s.m."),
                  "detailedCalculables" : True,
                  "rowColors" : self.rowcolors,
                  "rowCycle" : 100,
                  "omit2D" : True,
                  }
        
        supy.plotter(org, pdfFileName = self.pdfFileName(org.tag+"_log"),  doLog = True, pegMinimum = 0.01, **kwargs ).plotAll()
        supy.plotter(org, pdfFileName = self.pdfFileName(org.tag+"_nolog"), doLog = False, **kwargs ).plotAll()

        kwargs.update({"samplesForRatios":("",""),
                       "omit2D" : False,
                       "dependence2D" : True})
        supy.plotter(orgpdf, pdfFileName = self.pdfFileName(org.tag+"_pdf"), doLog = False, **kwargs ).plotAll()

    #def meldWpartitions(self,pars) :
    #    rw = pars['reweights']['abbr']
    #    samples = {"top_muon_pf_%s"%rw : ["w_"],
    #               "Wlv_muon_pf_%s"%rw : ["w_","SingleMu"],
    #               "QCD_muon_pf_%s"%rw : []}
    #    organizers = [supy.organizer(tag, [s for s in self.sampleSpecs(tag) if any(item in s['name'] for item in samples[tag])])
    #                  for tag in [p['tag'] for p in self.readyConfs]]
    #    if len(organizers)<2 : return
    #    for org in organizers :
    #        org.mergeSamples(targetSpec = {"name":"Data 2011", "color":r.kBlack, "markerStyle":20}, allWithPrefix="SingleMu")
    #        org.mergeSamples(targetSpec = {"name":"w_mg", "color":r.kRed if "Wlv" in org.tag else r.kBlue, "markerStyle": 22}, sources = ["wj_lv_mg.tw.%s"%rw])
    #        org.scale(toPdf=True)
    #
    #    melded = supy.organizer.meld("wpartitions",filter(lambda o: o.samples, organizers))
    #    pl = supy.plotter(melded,
    #                      pdfFileName = self.pdfFileName(melded.tag),
    #                      doLog = False,
    #                      blackList = ["lumiHisto","xsHisto","nJobsHisto"],
    #                      rowColors = self.rowcolors,
    #                      rowCycle = 100,
    #                      omit2D = True,
    #                      ).plotAll()
    #    
    #def meldQCDpartitions(self) :
    #    samples = {"top_muon_pf_%s"%rw : ["qcd_mu"],
    #               "Wlv_muon_pf_%s"%rw : [],
    #               "QCD_muon_pf_%s"%rw : ["qcd_mu","SingleMu"]}
    #    organizers = [supy.organizer(tag, [s for s in self.sampleSpecs(tag) if any(item in s['name'] for item in samples[tag])])
    #                  for tag in [p['tag'] for p in self.readyConfs]]
    #    if len(organizers)<2 : return
    #    for org in organizers :
    #        org.mergeSamples(targetSpec = {"name":"Data 2011", "color":r.kBlack, "markerStyle":20}, allWithPrefix="SingleMu")
    #        org.mergeSamples(targetSpec = {"name":"qcd_mu", "color":r.kRed if "QCD" in org.tag else r.kBlue, "markerStyle": 22}, allWithPrefix="qcd_mu")
    #        org.scale(toPdf=True)
    #
    #    melded = supy.organizer.meld("qcdPartitions",filter(lambda o: o.samples, organizers))
    #    pl = supy.plotter(melded,
    #                      pdfFileName = self.pdfFileName(melded.tag),
    #                      doLog = False,
    #                      blackList = ["lumiHisto","xsHisto","nJobsHisto"],
    #                      rowColors = self.rowcolors,
    #                      rowCycle = 100,
    #                      omit2D = True,
    #                      ).plotAll()
        
    def plotMeldScale(self, rw, lname) :
        if not hasattr(self,"orgMelded") : print "run meldScale() before plotMeldScale()"; return
        melded = copy.deepcopy(self.orgMelded)
        for ss in filter(lambda ss: 'ttj_mg' in ss['name'], melded.samples) : melded.drop(ss['name'])
        for log,label in [(False,""),(True,"_log")] : 
            pl = supy.plotter(melded, pdfFileName = self.pdfFileName(melded.tag + label),
                              doLog = log,
                              blackList = ["lumiHisto","xsHisto","nJobsHisto"],
                              rowColors = self.rowcolors,
                              samplesForRatios = ("top.Data 2011","S.M."),
                              sampleLabelsForRatios = ('data','s.m.'),
                              rowCycle = 100,
                              omit2D = True,
                              pageNumbers = False,
                              ).plotAll()
        
    def meldScale(self,rw,lname) :
        meldSamples = {"top_%s_pf_%s"%(lname,rw) : [{ 'muon':"SingleMu",
                                                      'electron':'EleHad'}[lname],
                                                    "ttj_mg",
                                                    #"wj_lv_mg",
                                                    "dyj_ll_mg"]+self.single_top()+["w%dj_mg"%n for n in [2,3,4]],
                       #"Wlv_%s_pf_%s"%(lname,rw) : ["w_jets"],
                       "QCD_%s_pf_%s"%(lname,rw) : [{ 'muon':"SingleMu",
                                                      'electron':'EleHad'}[lname],
                                                    "ttj_mg"]}
        
        organizers = [supy.organizer(tag, [s for s in self.sampleSpecs(tag) if any(item in s['name'] for item in meldSamples[tag])])
                      for tag in [p['tag'] for p in self.readyConfs if p["tag"] in meldSamples]]
        if len(organizers) < len(meldSamples) : return
        for org in organizers :
            org.mergeSamples(targetSpec = {"name":"t#bar{t}", "color":r.kViolet}, sources=["ttj_mg.wNonQQbar.tw.%s"%rw,"ttj_mg.wTopAsymP00.tw.%s"%rw], keepSources = True)
            #org.mergeSamples(targetSpec = {"name":"W", "color":r.kRed}, allWithPrefix = "wj_lv_mg")
            org.mergeSamples(targetSpec = {"name":"W", "color":r.kRed}, sources = ["w%dj_mg.tw.%s"%(n,rw) for n in [2,3,4]] )
            org.mergeSamples(targetSpec = {"name":"DY", "color":28}, allWithPrefix = "dyj_ll_mg")
            org.mergeSamples(targetSpec = {"name":"Single", "color":r.kGray}, sources = ["%s.tw.%s"%(s,rw) for s in self.single_top()], keepSources = False )
            org.mergeSamples(targetSpec = {"name":"Data 2011", "color":r.kBlack, "markerStyle":20}, allWithPrefix={'muon':"SingleMu",'electron':"EleHad"}[lname])
            org.scale()
            if "QCD_" in org.tag :
                org.mergeSamples(targetSpec = {"name":"multijet","color":r.kBlue},
                                 sources=["Data 2011",'t#bar{t}'],
                                 scaleFactors = [1,-1],
                                 force=True, keepSources = False)

        self.orgMelded = supy.organizer.meld(organizers = organizers)

        templateSamples = ['top.t#bar{t}','top.W','QCD.multijet']
        baseSamples = ['top.Single','top.DY']

        mfCanvas = r.TCanvas()
        def printMF(open = False, close = False) :
            name = "%s/measuredFractions.pdf"%self.globalStem
            mfCanvas.Print(name + ("[" if open else "]" if close else ""), "pdf" if open else "")
            if close : print supy.utils.hyphens, '\n', name, " has been written."
        printMF(open=True)
        
        def measureFractions(dist, rebin = 1) :
            before = next(self.orgMelded.indicesOfStep("label","selection complete"))
            distTup = self.orgMelded.steps[next(iter(filter(lambda i: before<i, self.orgMelded.indicesOfStepsWithKey(dist))))][dist]

            templates = [None] * len(templateSamples)
            bases = []
            for ss,hist in zip(self.orgMelded.samples,distTup) :
                contents = [hist.GetBinContent(i) for i in range(hist.GetNbinsX()+2)]
                if rebin!=1 : contents = contents[:1] + [sum(bins) for bins in zip(*[contents[1:-1][i::rebin] for i in range(rebin)])] + contents[-1:]
                if ss['name'] == "top.Data 2011" :
                    observed = contents
                    nEventsObserved = sum(observed)
                elif ss['name'] in templateSamples :
                    templates[templateSamples.index(ss['name'])] = contents
                elif ss['name'] in baseSamples :
                    bases.append(contents)
                else : pass

            from supy.utils.fractions import componentSolver,drawComponentSolver
            cs = componentSolver(observed, templates, 1e4, base = np.sum(bases, axis=0) )
            stuff = drawComponentSolver( cs, mfCanvas, distName = dist,
                                         templateNames = [t.replace("top.ttj_mg.wTopAsymP00.tw.%s"%rw,"q#bar{q}-->t#bar{t}").replace("top.ttj_mg.wNonQQbar.tw.%s"%rw,"gg-->t#bar{t}").replace("QCD.Data 2011","Multijet").replace("top.W","W+jets").replace('top.',"") for t in  templateSamples])
            printMF()
            return distTup,cs

        def mf2(dist) : return measureFractions(dist,2)

        distTup,cs = map(mf2,["TriDiscriminant","KarlsruheDiscriminant"])[0]

        fractions = dict(zip(templateSamples,cs.fractions))        
        for iSample,ss in enumerate(self.orgMelded.samples) :
            if ss['name'] in fractions : self.orgMelded.scaleOneRaw(iSample, fractions[ss['name']] * sum(cs.observed) / distTup[iSample].Integral(0,distTup[iSample].GetNbinsX()+1))

        self.orgMelded.mergeSamples(targetSpec = {"name":"bg", "color":r.kWhite}, sources = set(baseSamples + templateSamples) - set(['top.t#bar{t}']), keepSources = True, force = True)
        templateSamples = ['top.ttj_mg.wTopAsymP00.tw.%s'%rw,'top.ttj_mg.wNonQQbar.tw.%s'%rw]
        baseSamples = ['bg']
        distTup,cs = map(measureFractions,["tracksCountwithPrimaryHighPurityTracks","fitTopNtracksExtra","fitTopFifthJet","xcak5JetPFM3Pat","fitTopAbsSumRapidities"])[0]
                                           #"xcak5JetPFFourJetPtThresholdPat","xcak5JetPFFourJetAbsEtaThresholdPat","fitTopPartonXlo"])[0]
        printMF(close=True)

        #fractions = dict(zip(templateSamples,cs.fractions))
        #for iSample,ss in enumerate(self.orgMelded.samples) :
        #    if ss['name'] in fractions : self.orgMelded.scaleOneRaw(iSample, fractions[ss['name']] * sum(cs.observed) / distTup[iSample].Integral(0,distTup[iSample].GetNbinsX()+1))

        #self.orgMelded.drop("top.t#bar{t}")
        #self.orgMelded.mergeSamples(targetSpec = {"name":"top.t#bar{t}", "color":r.kViolet}, sources = templateSamples, keepSources = True, force = True)
        

        templateSamples = ['top.t#bar{t}'] # hack !!
        self.orgMelded.mergeSamples(targetSpec = {"name":"S.M.", "color":r.kGreen+2}, sources = templateSamples + baseSamples , keepSources = True, force = True)
        #for ss in filter(lambda ss: 'ttj_mg' in ss['name'],self.orgMelded.samples) : self.orgMelded.drop(ss['name'])
        self.orgMelded.drop('bg')
        
    def PEcurves(self) :
        if not hasattr(self, 'orgMelded') : return
        specs = ([{'var' : "ak5JetPFCSVPat[i[%d]]:xcak5JetPFIndicesBtaggedPat"%bIndex, 'left':True, 'right':False} for bIndex in [0,1,2]] +
                 [{'var' : "TopRatherThanWProbability",                                      'left':True, 'right':False},
                  {'var' : "TriDiscriminant",                                                'left':True, 'right':True}])
        pes = {}
        for spec in specs :
            dists = dict(zip([ss['name'] for ss in self.orgMelded.samples ],
                             self.orgMelded.steps[next(self.orgMelded.indicesOfStepsWithKey(spec['var']))][spec['var']] ) )
            contours = supy.utils.optimizationContours( [dists['top.t#bar{t}']],
                                                        [dists[s] for s in ['QCD.multijet','top.W','top.Single','top.DY']],
                                                        **spec
                                                        )
            supy.utils.tCanvasPrintPdf(contours[0], "%s/PE_%s"%(self.globalStem,spec['var']))
            if spec['left']^spec['right'] : pes[spec['var']] = contours[1]
        c = r.TCanvas()
        leg = r.TLegend(0.5,0.8,1.0,1.0)
        graphs = []
        for i,(var,pe) in enumerate(pes.iteritems()) :
            pur,eff = zip(*pe)
            g = r.TGraph(len(pe), np.array(eff), np.array(pur))
            g.SetTitle(";efficiency;purity")
            g.SetLineColor(i+2)
            leg.AddEntry(g,var,'l')
            graphs.append(g)
            g.Draw('' if i else 'AL')
        leg.Draw()
        c.Update()
        supy.utils.tCanvasPrintPdf(c, "%s/purity_v_efficiency"%self.globalStem)
        return


    def templates(self, iStep, dist, qqFrac, rw) :
        if not hasattr(self,'orgMelded') : print 'run meldScale() before asking for templates()'; return
        topQQs = [s['name'] for s in self.orgMelded.samples if 'wTopAsym' in s['name']]
        asymm = [eval(name.replace("top.ttj_mg.wTopAsym","").replace(".tw.%s"%rw,"").replace("P",".").replace("N","-.")) for name in topQQs if 'QCD' not in name]
        distTup = self.orgMelded.steps[iStep][dist]
        hist = distTup[ self.orgMelded.indexOfSampleWithName("S.M.") ]
        edges = supy.utils.edgesRebinned( hist, targetUncRel = 0.015, offset = 0 )

        def nparray(name, scaleToN = None) :
            hist_orig = distTup[ self.orgMelded.indexOfSampleWithName(name) ]
            hist = hist_orig.Rebin(len(edges)-1, "%s_rebinned"%hist_orig.GetName(), edges)
            bins = np.array([hist.GetBinContent(j) for j in range(hist.GetNbinsX()+2)])
            if scaleToN : bins *= (scaleToN / sum(bins))
            return bins

        nTT = sum(nparray('top.t#bar{t}'))
        observed = nparray('top.Data 2011')
        base = ( nparray('QCD.multijet') +
                 nparray('top.W') +
                 nparray('top.DY') +
                 nparray('top.Single') +
                 nparray('top.ttj_mg.wNonQQbar.tw.%s'%rw, scaleToN = (1-qqFrac) * nTT )
                 )
        templates = [base +  nparray(qqtt, qqFrac*nTT ) for qqtt in topQQs]
        return zip(asymm, templates), observed
    

    def ensembleFileName(self, iStep, dist, qqFrac, suffix = '.pickleData') :
        return "%s/ensembles/%d_%s_%.3f%s"%(self.globalStem,iStep,dist,qqFrac,suffix)

    def ensembleTest(self,rw) :
        qqFracs = sorted([0.10, 0.12, 0.15, 0.20, 0.25, 0.30, 0.40, 0.60, 1.0])
        dists = [#'lHadtDeltaY',
                 'ttbarDeltaAbsY',
                 #'leptonRelativeY',
                 #'ttbarSignedDeltaY',
                 'ttbarSignExpectation'
                ]
        args = sum([[(iStep, dist, qqFrac, rw) for iStep in list(self.orgMelded.indicesOfStepsWithKey(dist))[:None] for qqFrac in qqFracs] for dist in dists],[])
        supy.utils.operateOnListUsingQueue(6, supy.utils.qWorker(self.pickleEnsemble), args)
        ensembles = dict([(arg[:-1],supy.utils.readPickle(self.ensembleFileName(*arg[:-1]))) for arg in args])

        for iStep in sorted(set([iStep for iStep,dist,qqFrac in ensembles])) :
            canvas = r.TCanvas()
            dists = sorted(set([dist for jStep,dist,qqFrac in ensembles if jStep==iStep]))
            legend = r.TLegend(0.7,0.5,0.9,0.9)
            graphs = {}
            for iDist,dist in enumerate(dists) :
                points = sorted([(qqFrac,ensemble.sensitivity) for (jStep, jDist, qqFrac),ensemble in ensembles.iteritems() if jStep==iStep and jDist==dist])
                qqs,sens = zip(*points)
                graphs[dist] = r.TGraph(len(points),np.array(qqs),np.array(sens))
                graphs[dist].SetLineColor(iDist+1)
                graphs[dist].Draw('' if iDist else "AL")
                graphs[dist].SetMinimum(0)
                graphs[dist].SetTitle("Sensitivity @ step %d;fraction of t#bar{t} from q#bar{q};expected uncertainty on R"%iStep)
                legend.AddEntry(graphs[dist],dist,'l')
            legend.Draw()
            supy.utils.tCanvasPrintPdf(canvas, '%s/sensitivity_%d'%(self.globalStem,iStep))
                
    def pickleEnsemble(self, iStep, dist, qqFrac, rw) :
        supy.utils.mkdir(self.globalStem+'/ensembles')
        templates,observed = self.templates(iStep, dist, qqFrac, rw)
        ensemble = supy.utils.templateFit.templateEnsembles(2e3, *zip(*templates) )
        supy.utils.writePickle(self.ensembleFileName(iStep,dist,qqFrac), ensemble)

        name = self.ensembleFileName(iStep,dist,qqFrac,'')
        canvas = r.TCanvas()
        canvas.Print(name+'.ps[')
        stuff = supy.utils.templateFit.drawTemplateEnsembles(ensemble, canvas)
        canvas.Print(name+".ps")
        import random
        for i in range(20) :
            par, templ = random.choice(zip(ensemble.pars,ensemble.templates)[2:-2])
            pseudo = [np.random.poisson(mu) for mu in templ]
            tf = supy.utils.templateFit.templateFitter(pseudo,ensemble.pars,ensemble.templates, 1e3)
            stuff = supy.utils.templateFit.drawTemplateFitter(tf,canvas, trueVal = par)
            canvas.Print(name+".ps")
            for item in sum([i if type(i) is list else [i] for i in stuff[1:]],[]) : supy.utils.delete(item)

        canvas.Print(name+'.ps]')
        os.system('ps2pdf %s.ps %s.pdf'%(name,name))
        os.system('rm %s.ps'%name)
        
    def templateFit(self, iStep, dist, qqFrac = 0.15) :
        print "FIXME"; return
        if not hasattr(self,'orgMelded') : print 'run meldScale() before templateFit().'; return

        outName = self.globalStem + '/templateFit_%s_%d'%(dist,qqFrac*100)
        #TF = templateFit.templateFitter(observed, *zip(*templates) )
        #print supy.utils.roundString(TF.value, TF.error , noSci=True)
        
        #stuff = templateFit.drawTemplateFitter(TF, canvas)
        #canvas.Print(outName+'.ps')
        #for item in sum([i if type(i) is list else [i] for i in stuff[1:]],[]) : supy.utils.delete(item)
        
        #canvas.Print(outName+'.ps]')
        #os.system('ps2pdf %s.ps %s.pdf'%(outName,outName))
        #os.system('rm %s.ps'%outName)
