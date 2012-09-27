import supy,steps,calculables,samples
import os,math,copy,itertools,ROOT as r, numpy as np

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
            'abbr'  :'pu',
            'func'  :'pileup',
            'var'   :'pileupTrueNumInteractionsBX0Target'
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

        return { "vary" : ['selection','lepton','toptype'],
                 "discriminant2DPlots": True,
                 "nJets" :  {"min":4,"max":None},
                 "unreliable": self.unreliableTriggers(),
                 "bVar" : "CSV", # "Combined Secondary Vertex"
                 "objects" : dict((key,val[0]) for key,val in objects.iteritems()),
                 "lepton" : self.vary([ ( leptons['name'][index], dict((key,val[index]) for key,val in leptons.iteritems())) for index in range(2) if leptons['name'][index] in ['muon','electron'][:2]]),
                 "reweights" : reweights,
                 "selection" : self.vary({"top" : {"bCut":bCut["normal"],  "iso":"isoNormal"},
                                          "QCD" : {"bCut":bCut["normal"],  "iso":"isoInvert"}
                                          }),
                 "toptype" : self.vary({"mg":"mg","ph":"ph"}),
                 "topBsamples": ("ttj_%s",['ttj_%s.wGG.tw.%s','ttj_%s.wQG.tw.%s','ttj_%s.wAG.tw.%s','ttj_%s.wQQ.tw.%s'])
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

        def ewk(eL = None) :
            return supy.samples.specify( names = ["dyj_ll_mg",
                                                  "w2j_mg",
                                                  "w3j_mg",
                                                  "w4j_mg"], effectiveLumi = eL, color = 28, weights = ['tw',rw] ) if "QCD" not in pars['tag'] else []

        def single_top(eL = None) :
            return supy.samples.specify( names = self.single_top(),
                                         effectiveLumi = eL, color = r.kGray, weights = ['tw',rw]) if "QCD" not in pars['tag'] else []
        
        def ttbar(eL = None) :
            tt = pars['toptype']
            return (supy.samples.specify(names = "ttj_%s"%tt, effectiveLumi = eL, color = r.kBlue, weights = ["wGG",'tw',rw] ) +
                    supy.samples.specify(names = "ttj_%s"%tt, effectiveLumi = eL, color = r.kCyan, weights = [ "wQG", 'tw', rw ] ) +
                    supy.samples.specify(names = "ttj_%s"%tt, effectiveLumi = eL, color = r.kCyan+1, weights = [ "wAG", 'tw', rw ] ) +
                    supy.samples.specify(names = "ttj_%s"%tt, effectiveLumi = eL, color = r.kOrange, weights = [ "wQQ", 'tw', rw ] )
                    )
        
        return  ( self.data(pars) + ewk() + ttbar(8e4) + single_top() )

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

            calculables.gen.genIndicesHardPartons(),
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
            calculables.other.KarlsruheDiscriminant( obj['jet'], obj['met'] ),

            calculables.jet.pt( obj['jet'], index = 0, Btagged = True ),
            calculables.jet.pt( obj['jet'], index = 3, Btagged = False ),
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
        topTag = pars['tag'].replace("QCD","top")
        bVar = pars["bVar"].join(obj["jet"])[2:]
        rw = pars['reweights']['abbr']
        tt = pars['toptype']
        
        ssteps = supy.steps

        saDisable = 'ttj' not in pars['sample'] or not any(w in pars['sample'].split('.') for w in ['wQQ','wQG','wAG','wGG'])
        saWeights = []
        return (
            [ssteps.printer.progressPrinter()
             , ssteps.histos.value("genpthat",200,0,1000,xtitle="#hat{p_{T}} (GeV)").onlySim()
             , ssteps.histos.value("genQ",200,0,1000,xtitle="#hat{Q} (GeV)").onlySim()
             , getattr(self,pars['reweights']['func'])(pars),
             supy.calculables.other.SymmAnti(pars['sample'],"genTopCosPhiBoost",1, inspect=True, nbins=160, weights = saWeights,
                                             funcEven = r.TF1('phiboost',"[0]*(1+[1]*x**2)/sqrt(1-x**2)",-1,1),
                                             funcOdd = r.TF1('phiboostodd','[0]*x/sqrt(1-x**2)',-1,1)).disable(saDisable),
             supy.calculables.other.SymmAnti(pars['sample'],"genTopCosThetaBoostAlt",1, inspect=True, weights = saWeights,
                                             funcEven = '++'.join('x**%d'%(2*d) for d in range(5)),
                                             funcOdd = '++'.join('x**%d'%(2*d+1) for d in range(5))).disable(saDisable),
             supy.calculables.other.SymmAnti(pars['sample'],"genTopDeltaBetazRel",1, inspect=True, weights = saWeights,
                                             funcEven = '++'.join(['(1-abs(x))']+['x**%d'%d for d in [0,2,4,6,8,10,12,14,16,18]]),
                                             funcOdd = '++'.join(['x**%d'%d for d in [1,3,5,7,9,11,13]])).disable(saDisable)
             , ssteps.filters.label('symm anti')
             , ssteps.filters.value('tw',min=1e-8)
             , ssteps.histos.symmAnti('genTopCosPhiBoost','genTopCosPhiBoost',100,-1,1).disable(saDisable)
             , ssteps.histos.symmAnti('genTopDeltaBetazRel','genTopDeltaBetazRel',100,-1,1).disable(saDisable)
             , ssteps.histos.symmAnti('genTopCosThetaBoostAlt','genTopCosThetaBoostAlt',100,-1,1).disable(saDisable)
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

             , ssteps.filters.label('trigger reweighting')
             , self.triggerWeight(pars, [ss.weightedName for ss in self.data(pars)])
             , steps.trigger.lowestUnPrescaledTriggerHistogrammer().onlyData()

             ####################################
             , ssteps.filters.label('selection'),

             ssteps.filters.multiplicity( min=1, var = "IndicesAnyIso".join(lepton) )
             , ssteps.histos.multiplicity("Indices".join(obj["jet"])),
             ssteps.filters.multiplicity("Indices".join(obj["jet"]), **pars["nJets"]),
             ssteps.filters.pt("CorrectedP4".join(obj['jet']), min = 40, indices = "Indices".join(obj['jet']), index=0),
             ssteps.filters.pt("CorrectedP4".join(obj['jet']), min = 30, indices = "Indices".join(obj['jet']), index=1)
             
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

             , calculables.jet.ProbabilityGivenBQN(obj["jet"], pars['bVar'], binning=(51,-0.02,1), samples = (pars['topBsamples'][0]%tt,[s%(tt,rw) for s in pars['topBsamples'][1]]), tag = topTag)
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
             , steps.top.leptonSigned('TridiscriminantWTopQCD', (60,-1,1))
             , steps.top.leptonSigned('KarlsruheDiscriminant', (28,-320,800) )
             , ssteps.filters.label('discriminants')
             , ssteps.histos.value("KarlsruheDiscriminant", 28, -320, 800 )
             , self.tridiscriminant(pars)
             , self.tridiscriminant2(pars)
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
             
             , ssteps.filters.label('signal distributions')
             , ssteps.histos.symmAnti('genTopCosPhiBoost','genTopCosPhiBoost',100,-1,1).disable(saDisable)
             , ssteps.histos.symmAnti('genTopDeltaBetazRel','genTopDeltaBetazRel',100,-1,1).disable(saDisable)
             , ssteps.histos.symmAnti('genTopCosThetaBoostAlt','genTopCosThetaBoostAlt',100,-1,1).disable(saDisable)

             , ssteps.histos.symmAnti('genTopCosPhiBoost','fitTopCosPhiBoost',100,-1,1)
             , ssteps.histos.symmAnti('genTopCosThetaBoostAlt','fitTopCosThetaBoostAlt',100,-1,1)
             , ssteps.histos.symmAnti('genTopDeltaBetazRel','fitTopDeltaBetazRel',100,-1,1)

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
        tt = pars['toptype']
        return supy.calculables.other.Target("pileupTrueNumInteractionsBX0", thisSample = pars['baseSample'],
                                             target = ("data/pileup_true_Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON.root","pileup"),
                                             groups = [('qcd_mu',[]),('wj_lv_mg',[]),('dyj_ll_mg',[]),
                                                       ("w2j_mg",[]),("w3j_mg",[]),("w4j_mg",[]),
                                                       ('single_top', ['%s.tw.%s'%(s,rw) for s in cls.single_top()]),
                                                       ('ttj_%s'%tt,['ttj_%s.%s.tw.%s'%(tt,s,rw) for s in ['',
                                                                                                           'wGG','wQG','wAG','wQQ']])]).onlySim()
    @staticmethod
    def tridiscriminant(pars) :
        rw = pars['reweights']['abbr']
        lname = pars['lepton']['name']
        tt = pars['toptype']
        lumi = 5008 # FIXME HACK !!!
        tops = ['ttj_%s.%s.tw.%s'%(tt,s,rw) for s in ['wGG','wQG','wAG','wQQ']]
        datas = {"muon" : ["SingleMu.2011A.tw","SingleMu.2011B.tw"],
                 "electron": ["EleHad.2011A.tw","EleHad.2011B.tw"]}[lname]
        return supy.calculables.other.Tridiscriminant( fixes = ("","WTopQCD"),
                                                       zero = {"pre":"ttj_%s"%tt, "tag":"top_%s_%s"%(lname,tt), "samples": tops},
                                                       pi23 = {"pre":"Multijet", "tag":"QCD_%s_%s"%(lname,tt), "samples":datas+tops, 'sf':[1,1,-lumi,-lumi,-lumi]},
                                                       pi43 = {"pre":"wj", "tag":"top_%s_%s"%(lname,tt), "samples":['w%dj_mg.tw.%s'%(n,rw) for n in [2,3,4]]},
                                                       correlations = pars['discriminant2DPlots'],
                                                       dists = {"TopRatherThanWProbability" : (20,0.5,1),
                                                                "B0pt".join(pars['objects']["jet"]) : (20,0,100),
                                                                "fitTopHadChi2"     : (20,0,20),
                                                                "Mt".join(pars['objects'][lname])+"mixedSumP4" : (20,0,100),
                                                                })
    @staticmethod
    def tridiscriminant2(pars) :
        rw = pars['reweights']['abbr']
        lname = pars['lepton']['name']
        tt = pars['toptype']
        jets = pars["objects"]['jet']
        return supy.calculables.other.Tridiscriminant( fixes = ("","GGqqGq"),
                                                       zero = {"pre":"gg", "tag":"top_%s_%s"%(lname,tt), "samples":['ttj_%s.wGG.tw.%s'%(tt,rw)]},
                                                       pi23 = {"pre":"qq", "tag":"top_%s_%s"%(lname,tt), "samples":['ttj_%s.wQQ.tw.%s'%(tt,rw)]},
                                                       pi43 = {"pre":"qg", "tag":"top_%s_%s"%(lname,tt), "samples":['ttj_%s.%s.tw.%s'%(tt,s,rw) for s in ['wQG','wAG']]},
                                                       correlations = pars['discriminant2DPlots'],
                                                       dists = { "fitTopPtPlusSumPt" : (20,0,600),
                                                                 "fitTopPtOverSumPt" : (20,0,1),
                                                                 "fitTopSumP4AbsEta" : (20,0,6),
                                                                 #"fitTopAbsSumRapidities" : (20, 0, 4),
                                                                 #"M3".join(jets) : (20,0,600),
                                                                 #"fitTopPtSum" : (30, 0, 150), # extra jet
                                                                 #"fitTopPartonLo" : (20,-0.2,0.2),
                                                                 #"fitTopPartonHi" : (20,0,0.4),
                                                                 #"fitTopMassSum" : (30, 300, 900), # pdf sum
                                                                 #"fitTopRapiditySum" : (20, 0, 2), # pdf difference
                                                                 #"fitTopNtracksExtra" : (20,0,160),
                                                                 #"tracksCountwithPrimaryHighPurityTracks" :  (20,0,250),               # 0.049
                                                                 #"fitTopPartonXlo" : (20,0,0.12),              # 0.036
                                                                 #"fitTopBMomentsSum2" : (20,0,0.2),           # 0.004
                                                                 #"fitTopPartonXhi" : (20,0.04,0.4),           # 0.003
                                                                 })
    ########################################################################################
    def concludeAll(self) :
        self.orgMelded = {}
        self.sensitivityPoints = {}
        self.rowcolors = 2*[13] + 2*[45]
        super(topAsymm,self).concludeAll()
        for tt,rw,lname in set([(pars['toptype'],pars['reweights']['abbr'],pars['lepton']['name']) for pars in self.readyConfs]) :
            self.meldScale(rw,lname,tt)
            self.plotMeldScale(rw,lname,tt)
            self.PEcurves(rw,lname, tt)
        for tt,rw in set([(pars['toptype'],pars['reweights']['abbr']) for pars in self.readyConfs]) :
            self.measureAmplitudes(rw,tt)
        #self.sensitivity_graphs()
        #self.grant_proposal_plots()

    def conclude(self,pars) :
        rw = pars['reweights']['abbr']
        tt = pars['toptype']
        org = self.organizer(pars, verbose = True )

        if pars['lepton']['name']=='muon' :
            org.mergeSamples(targetSpec = {"name":"SingleMu.2011", "color":r.kBlack, "markerStyle":20}, allWithPrefix="SingleMu")
        else:
            org.mergeSamples(targetSpec = {"name":"EleHad.2011", "color":r.kBlack, "markerStyle":20}, allWithPrefix="EleHad")
            
        org.mergeSamples(targetSpec = {"name":"t#bar{t}", "color":r.kViolet}, sources=["ttj_%s.%s.tw.%s"%(tt,s,rw) for s in ['wQQ','wQG','wAG','wGG']], keepSources = True)
        org.mergeSamples(targetSpec = {"name":"W+jets", "color":28}, sources = ["w%dj_mg.tw.%s"%(n,rw) for n in [2,3,4]])
        org.mergeSamples(targetSpec = {"name":"DY+jets", "color":r.kYellow}, allWithPrefix="dyj_ll_mg")
        org.mergeSamples(targetSpec = {"name":"Single top", "color":r.kGray}, sources = ["%s.tw.%s"%(s,rw) for s in self.single_top()])
        org.mergeSamples(targetSpec = {"name":"Standard Model", "color":r.kGreen+2}, sources = ["multijet","t#bar{t}","W+jets","DY+jets","Single top"], keepSources = True)

        org.scale( lumiToUseInAbsenceOfData = 5008 )

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

    def plotMeldScale(self,rw,lname,tt) :
        if (lname,rw, tt) not in self.orgMelded : print "run meldScale() before plotMeldScale()"; return
        melded = copy.deepcopy(self.orgMelded[(lname,rw,tt)])
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

    def meldScale(self,rw,lname,tt) :
        meldSamples = {"top_%s_%s"%(lname,tt) : [{ 'muon':"SingleMu",
                                                   'electron':'EleHad'}[lname],
                                                 "ttj_%s"%tt,
                                                 "dyj_ll_mg"]+self.single_top()+["w%dj_mg"%n for n in [2,3,4]],
                       "QCD_%s_%s"%(lname,tt) : [{ 'muon':"SingleMu",
                                                   'electron':'EleHad'}[lname],
                                                 "ttj_%s"%tt]}

        organizers = [supy.organizer(tag, [s for s in self.sampleSpecs(tag) if any(item in s['name'] for item in meldSamples[tag])])
                      for tag in [p['tag'] for p in self.readyConfs if p["tag"] in meldSamples]]

        if len(organizers) < len(meldSamples) : return
        for org in organizers :
            org.mergeSamples(targetSpec = {"name":"t#bar{t}", "color":r.kViolet}, sources=["ttj_%s.%s.tw.%s"%(tt,s,rw) for s in ['wQQ','wQG','wAG','wGG']], keepSources = 'top' in org.tag)
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

        self.orgMelded[(lname,rw,tt)] = supy.organizer.meld(organizers = organizers)
        org = self.orgMelded[(lname,rw,tt)]
        templateSamples = ['top.t#bar{t}','top.W','QCD.multijet']
        baseSamples = ['top.Single','top.DY']

        mfCanvas = r.TCanvas()
        mfFileName = "%s/%s_measuredFractions"%(self.globalStem, lname )
        supy.utils.tCanvasPrintPdf( mfCanvas, mfFileName, option = '[', verbose = False)
        with open(mfFileName+'.txt','w') as file : print >> file, ""
        
        def measureFractions(dist, rebin = 1) :
            before = next(org.indicesOfStep("label","selection complete"))
            distTup = org.steps[next(iter(filter(lambda i: before<i, org.indicesOfStepsWithKey(dist))))][dist]

            templates = [None] * len(templateSamples)
            bases = []
            for ss,hist in zip(org.samples,distTup) :
                contents = supy.utils.binValues(hist)
                if rebin!=1 : contents = contents[:1] + [sum(bins) for bins in zip(*[contents[1:-1][i::rebin] for i in range(rebin)])] + contents[-1:]
                if ss['name'] == "top.Data 2011" :
                    observed = contents
                elif ss['name'] in templateSamples :
                    templates[templateSamples.index(ss['name'])] = contents
                elif ss['name'] in baseSamples :
                    bases.append(contents)
                else : pass

            from supy.utils.fractions import componentSolver,drawComponentSolver
            cs = componentSolver(observed, templates, 1e4, base = np.sum(bases, axis=0) )
            stuff = drawComponentSolver( cs, mfCanvas, distName = dist,
                                         templateNames = [t.replace("top.ttj_%s.wQQ.tw.%s"%(tt,rw),"q#bar{q}-->t#bar{t}").replace("top.ttj_%s.wQG.tw.%s"%(tt,rw),"qg-->t#bar{t}").replace("top.ttj_%s.wAG.tw.%s"%(tt,rw),"#bar{q}g-->t#bar{t}").replace("top.ttj_%s.wGG.tw.%s"%(tt,rw),"gg-->t#bar{t}").replace("QCD.Data 2011","Multijet").replace("top.W","W+jets").replace('top.',"") for t in  templateSamples])
            supy.utils.tCanvasPrintPdf( mfCanvas, mfFileName, verbose = False)
            with open(mfFileName+'.txt','a') as file : print >> file, "\n",dist+"\n", cs
            return distTup,cs

        def mf2(dist) : return measureFractions(dist,1)

        distTup,cs = map(mf2,["KarlsruheDiscriminant","TridiscriminantWTopQCD"])[-1]

        iTT = next(i for i,ss in enumerate(org.samples) if ss['name']=='top.t#bar{t}')
        nTT = distTup[iTT].Integral(0,distTup[iTT].GetNbinsX()+1)
        fractions = dict(zip(templateSamples,cs.fractions))
        for iSample,ss in enumerate(org.samples) :
            if ss['name'] in fractions :
                f = fractions[ss['name']]
                n = distTup[iSample].Integral(0,distTup[iSample].GetNbinsX()+1)
            elif ss['name'] in ['top.ttj_%s.%s.tw.%s'%(tt,s,rw) for s in ['wQQ','wQG','wAG','wGG']] :
                f = fractions['top.t#bar{t}']
                n = nTT
            else : continue
            org.scaleOneRaw(iSample, f * sum(cs.observed) / n )

        org.mergeSamples(targetSpec = {"name":"bg", "color":r.kBlack,"fillColor":r.kGray, "markerStyle":1, "goptions":"hist"}, sources = set(baseSamples + templateSamples) - set(['top.t#bar{t}']), keepSources = True, force = True)
        templateSamples = ['top.ttj_%s.%s.tw.%s'%(tt,s,rw) for s in ['wQQ','wQG','wAG','wGG']]
        baseSamples = ['bg']
        distTup,cs = map(measureFractions,["fitTopPtOverSumPt","fitTopPtPlusSumPt","fitTopSumP4AbsEta","TridiscriminantGGqqGq"])[-1]
        org.mergeSamples(targetSpec = {"name":'qgqqbar'}, sources = templateSamples[:-1], keepSources = True, force = True)
        templateSamples = ['qgqqbar', templateSamples[-1]]
        distTup,cs = map(measureFractions,["fitTopPtOverSumPt","fitTopPtPlusSumPt","fitTopSumP4AbsEta","TridiscriminantGGqqGq"])[-1]
        supy.utils.tCanvasPrintPdf( mfCanvas, mfFileName, option = ']')
        org.drop('qgqqbar')

        templateSamples = ['top.t#bar{t}'] # hack !!
        org.mergeSamples(targetSpec = {"name":"S.M.", "color":r.kGreen+2}, sources = templateSamples + baseSamples , keepSources = True, force = True)
        org.drop('bg')

    def PEcurves(self, rw, lname, tt) :
        if (lname,rw,tt) not in self.orgMelded : return
        org = self.orgMelded[(lname,rw,tt)]
        c = r.TCanvas()
        fileName = "%s/pur_eff_%s_%s"%(self.globalStem,lname,rw)
        supy.utils.tCanvasPrintPdf(c, fileName, option = '[', verbose = False)
        specs = ([{'var' : "TopRatherThanWProbability",                                "canvas":c, 'left':True, 'right':False}] +
                 [{'var' : "ak5JetPFCSVPat[i[%d]]:xcak5JetPFIndicesBtaggedPat"%bIndex, "canvas":c, 'left':True, 'right':False} for bIndex in [0,1,2]]
                 )
        pes = {}
        for spec in specs :
            dists = dict(zip([ss['name'] for ss in org.samples ],
                             org.steps[next(org.indicesOfStepsWithKey(spec['var']))][spec['var']] ) )
            contours = supy.utils.optimizationContours( [dists['top.t#bar{t}']],
                                                        [dists[s] for s in ['QCD.multijet','top.W','top.Single','top.DY']],
                                                        **spec
                                                        )
            supy.utils.tCanvasPrintPdf(c, fileName, verbose = False)
            if spec['left']^spec['right'] : pes[spec['var']] = contours[1]
            c.Clear()
        leg = r.TLegend(0.5,0.8,1.0,1.0)
        graphs = []
        for i,(var,pe) in enumerate(pes.items()) :
            pur,eff = zip(*pe)
            g = r.TGraph(len(pe), np.array(eff), np.array(pur))
            g.SetTitle(";efficiency;purity")
            g.SetLineColor(i+2)
            leg.AddEntry(g,var,'l')
            graphs.append(g)
            g.Draw('' if i else 'AL')
        leg.Draw()
        c.Update()
        supy.utils.tCanvasPrintPdf(c, fileName, option = ')' )
        return


    def measureAmplitudes(self,rw,tt) :
        lnames = ['muon','electron']
        if not all((lname,rw, tt) in self.orgMelded for lname in lnames):
            print "run meldScale() before measureAmplitudes()"
            return
        orgMuEl = [self.orgMelded[(lname,rw,tt)] for lname in lnames]
        omitSamples = ["S.M.","top.Data 2011","top.t#bar{t}"]
        for org in orgMuEl :
            print ", ".join(ss["name"] for ss in org.samples if ss["name"] not in omitSamples)
        canvas = r.TCanvas()
        maFileName = "%s/measuredAmplitudes"%(self.globalStem)
        supy.utils.tCanvasPrintPdf( canvas, maFileName, option = '[', verbose = False)
        with open(maFileName+'.txt','w') as file : print >> file, ""

        def BV(h, rebin=4) :
            contents = supy.utils.binValues(h)
            if rebin!=1 : contents = contents[:1] + [sum(bins) for bins in zip(*[contents[1:-1][i::rebin] for i in range(rebin)])] + contents[-1:]
            return contents

        SA = supy.utils.symmAnti
        def measure(fitvar,genvar, antisamples, sumTemplates = False) :
            steps = [org.steps[ next(org.indicesOfStep('symmAnti','%s in (anti)symm parts of %s'%(fitvar,genvar))) ] for org in orgMuEl]
            templatess = [[ BV( SA(step[fitvar+'_anti'][org.indexOfSampleWithName(sample)])[1]) for sample in antisamples ] for step,org in zip(steps,orgMuEl)]
            basess = [ [ BV( SA(step[fitvar+'_symm'][org.indexOfSampleWithName(sample)])[0]) for sample in antisamples ] +
                       [ BV( step[fitvar][i]) for i,ss in enumerate(org.samples) if ss['name'] not in antisamples+omitSamples ]
                       for step,org in zip(steps,orgMuEl)]
            observeds = [ BV( step[fitvar][org.indexOfSampleWithName("top.Data 2011")] ) for step,org in zip(steps,orgMuEl) ]

            def stack(Lists) : return [sum(lists,[]) for lists in zip(*Lists)]
            templates = stack(templatess)
            bases = stack(basess)
            observed = sum(observeds,[])

            from supy.utils.fractions import componentSolver,drawComponentSolver
            if sumTemplates : templates = [np.sum(templates, axis=0)]
            cs = componentSolver(observed, templates, 1e4, base = np.sum(bases, axis=0) , normalize = False )
            stuff = drawComponentSolver(cs, canvas, distName = fitvar, showDifference = True,
                                        templateNames = ["antisymmetric %s-->t#bar{t}"%('(qg|q#bar{q})' if sumTemplates else "qg" if 'QG' in t else "q#bar{q}" if 'QQ' in t else '') for t in antisamples])
            supy.utils.tCanvasPrintPdf( canvas, maFileName, verbose = False)
            with open(maFileName+'.txt','a') as file : print >> file, "\n",fitvar+"\n", cs
            return steps,cs

        samples = ['top.ttj_%s.%s.tw.%s'%(tt,w,rw) for w in ['wQQ','wQG','wAG']]
        measure('fitTopCosPhiBoost','genTopCosPhiBoost', samples[1:], sumTemplates=True) #qg only
        measure('fitTopDeltaBetazRel','genTopDeltaBetazRel', samples, sumTemplates=True)
        #measure('fitTopDeltaBetazRel','genTopDeltaBetazRel', samples)
        measure('fitTopCosThetaBoostAlt','genTopCosThetaBoostAlt', samples, sumTemplates=True)
        #measure('fitTopCosThetaBoostAlt','genTopCosThetaBoostAlt', samples)
        supy.utils.tCanvasPrintPdf( canvas, maFileName, option = ']')
