import supy,steps,calculables,samples
import os,math,copy,ROOT as r, numpy as np

class topAsymm(supy.analysis) :

    def parameters(self) :

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
            'name'     : [               'muon',      'electron'],
            'ptMin'    : [                20.0 ,           30.0 ],
            'etaMax'   : [                 2.1 ,            9.0 ],
            'iso'      : ['CombinedRelativeIso',   'IsoCombined'],
            'triggers' : [    self.mutriggers(),           None ]
            }

        bCut = {"normal"   : {"index":1, "min":2.0},
                "inverted" : {"index":1, "max":2.0}}
        lIso = {"normal":  {"N":1, "indices":"Indices"},
                "inverted":{"N":0, "indices":"IndicesNonIso"}}

        return { "vary" : ['selection','lepton','objects'],
                 "discriminant2DPlots": True,
                 "nJets" :  {"min":4,"max":None},
                 "unreliable": self.unreliableTriggers(),
                 "bVar" : "NTrkHiEff", # "TrkCountingHighEffBJetTags"
                 "objects": self.vary([ ( objects['label'][index], dict((key,val[index]) for key,val in objects.iteritems())) for index in range(2) if objects['label'][index] in ['pf']]),
                 "lepton" : self.vary([ ( leptons['name'][index], dict((key,val[index]) for key,val in leptons.iteritems())) for index in range(2) if leptons['name'][index] in ['muon']]),
                 "selection" : self.vary({"top" : {"bCut":bCut["normal"],  "lIso":lIso["normal"]},
                                          "QCD" : {"bCut":bCut["normal"],  "lIso":lIso["inverted"]}
                                          #"Wlv" : {"bCut":bCut["inverted"],"lIso":lIso["normal"]},
                                          }),
                 "topBsamples": { "pythia"   : ("tt_tauola_fj",["tt_tauola_fj.wNonQQbar.tw.nvr","tt_tauola_fj.wTopAsymP00.tw.nvr"]),
                                  "madgraph" : ("FIXME",[]),
                                  }["pythia"]
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

    ########################################################################################

    def listOfSampleDictionaries(self) : return [getattr(samples,item) for item in ['muon', 'top', 'ewk', 'qcd']]

    def data(self,pars) :
        return supy.samples.specify( names = ['SingleMu.2011B-PR1.1b',
                                              'SingleMu.2011B-PR1.1a',
                                              'SingleMu.2011A-Oct.1',
                                              'SingleMu.2011A-Aug.1',
                                              'SingleMu.2011A-PR4.1',
                                              'SingleMu.2011A-May.1'], weights = 'tw')
    def listOfSamples(self,pars) :

        def qcd_py6_mu(eL = None) :
            return supy.samples.specify( names = ["qcd_py6fjmu_pt_15_20",
                                                  "qcd_py6fjmu_pt_20_30",
                                                  "qcd_py6fjmu_pt_30_50",
                                                  "qcd_py6fjmu_pt_50_80",
                                                  "qcd_py6fjmu_pt_80_120",
                                                  "qcd_py6fjmu_pt_120_150",
                                                  "qcd_py6fjmu_pt_150"],  weights = ['tw','nvr'], effectiveLumi = eL )     if 'Wlv' not in pars['tag'] else []
        def qcd_mg(eL = None) :
            return supy.samples.specify( names = ["qcd_mg_ht_50_100",
                                                  "qcd_mg_ht_100_250",
                                                  "qcd_mg_ht_250_500",
                                                  "qcd_mg_ht_500_1000",
                                                  "qcd_mg_ht_1000_inf"],  weights = ['tw','nvr'], effectiveLumi = eL )     if 'Wlv' not in pars['tag'] else []
        def ewk(eL = None) :
            return supy.samples.specify( names = ["w_jets_fj_mg","dyj_ll_mg"], effectiveLumi = eL, color = 28, weights = ['tw','nvr'] ) if "QCD" not in pars['tag'] else []

        def ttbar_mg(eL = None) :
            return (supy.samples.specify( names = "tt_tauola_mg", effectiveLumi = eL, color = r.kBlue, weights = ['wNonQQbar','tw','nvr']) +
                    sum([supy.samples.specify( names = "tt_tauola_mg", effectiveLumi = eL, 
                                               color = color, weights = [calculables.top.wTopAsym( asym, R_sm = -0.05), 'nvr' ])
                         for asym,color in [(0.0,r.kOrange),
                                            (-0.3,r.kGreen),
                                            (0.3,r.kRed)
                                            ]], [])
                    )[: 0 if "QCD" in pars['tag'] else 2 if 'Wlv' in pars['tag'] else None]
        
        def ttbar_py(eL = None) :
            return (supy.samples.specify(names = "tt_tauola_fj", effectiveLumi = eL, color = r.kBlue, weights = ["wNonQQbar",'tw','nvr']) +
                    sum( [supy.samples.specify( names = "tt_tauola_fj", effectiveLumi = eL, 
                                                color = color, weights = [ calculables.top.wTopAsym(asym), 'tw','nvr' ] )
                          for asym,color in [(0.0,r.kOrange),
                                             (-0.3,r.kGreen),(0.3,r.kRed),
                                             #(-0.6,r.kYellow),(0.6,r.kYellow),
                                             #(-0.5,r.kYellow),(0.5,r.kYellow),
                                             #(-0.4,r.kYellow),(0.4,r.kYellow),
                                             #(-0.2,r.kYellow),(0.2,r.kYellow),
                                             #(-0.1,r.kYellow),(0.1,r.kYellow),
                                             ]], [])
                    )[: 0 if "QCD" in pars['tag'] else 2 if 'Wlv' in pars['tag'] else None]
        
        return  ( self.data(pars) + qcd_py6_mu(100) + ewk() + ttbar_py() )


    ########################################################################################
    def listOfCalculables(self, pars) :
        obj = pars["objects"]
        lepton = obj[pars["lepton"]["name"]]
        calcs  = supy.calculables.zeroArgs(supy.calculables)
        calcs += supy.calculables.zeroArgs(calculables)
        for item in ['jet','photon','electron','muon'] :
            calcs += supy.calculables.fromCollections( getattr(calculables, item), [obj[item]])
        calcs += [
            calculables.jet.IndicesBtagged(obj["jet"],pars["bVar"]),
            calculables.jet.Indices(       obj["jet"],      ptMin = 20, etaMax = 3.5, flagName = "JetIDloose"),
            calculables.photon.Indices(    obj["photon"],   ptMin = 25, flagName = "photonIDLooseFromTwikiPat"),
            calculables.electron.Indices(  obj["electron"], ptMin = 10, simpleEleID = "80", useCombinedIso = True),
            calculables.muon.Indices(      obj["muon"],     ptMin = 10, combinedRelIsoMax = 0.15),
            calculables.muon.IndicesTriggering(obj["muon"]),
            calculables.muon.IndicesAnyIsoIsoOrder(obj[pars["lepton"]["name"]], pars["lepton"]["iso"]),

            calculables.xclean.IndicesUnmatched(collection = obj["photon"], xcjets = obj["jet"], DR = 0.5),
            calculables.xclean.IndicesUnmatched(collection = obj["electron"], xcjets = obj["jet"], DR = 0.5),
            calculables.xclean.xcJet(obj["jet"], applyResidualCorrectionsToData = False,
                                     gamma    = obj["photon"],      gammaDR = 0.5,
                                     electron = obj["electron"], electronDR = 0.5,
                                     muon     = obj["muon"],         muonDR = 0.5, correctForMuons = not obj["muonsInJets"]),

            calculables.vertex.ID(),
            calculables.vertex.Indices(),

            calculables.top.mixedSumP4(transverse = obj["met"], longitudinal = obj["sumP4"]),
            calculables.top.SemileptonicTopIndex(lepton),
            calculables.top.fitTopLeptonCharge(lepton),
            calculables.top.TopReconstruction(lepton,obj["jet"],"mixedSumP4"),

            calculables.top.TopComboQQBBLikelihood(pars['objects']['jet'], pars['bVar']),
            calculables.top.OtherJetsLikelihood(pars['objects']['jet'], pars['bVar']),
            calculables.top.TopRatherThanWProbability(priorTop=0.5),
            calculables.top.RadiativeCoherence(('fitTop',''),pars['objects']['jet']),

            calculables.other.Mt(lepton,"mixedSumP4", allowNonIso=True, isSumP4=True),
            calculables.other.Covariance(('met','PF')),
            calculables.other.PtSorted(obj['muon']),
            calculables.other.lowestUnPrescaledTrigger(zip(*pars["lepton"]["triggers"])[0]),

            calculables.other.TriDiscriminant(LR = "DiscriminantWQCD", LC = "DiscriminantTopW", RC = "DiscriminantTopQCD"),
            calculables.other.KarlsruheDiscriminant(pars['objects']['jet'], pars['objects']['met']),

            calculables.jet.pt( obj['jet'], index = 0, Btagged = True ),
            calculables.jet.absEta( obj['jet'], index = 3, Btagged = False),
            supy.calculables.other.size("%sIndices%s"%pars['objects']['jet']),
            supy.calculables.other.pt("mixedSumP4"),
            
            supy.calculables.other.abbreviation( "TrkCountingHighEffBJetTags", "NTrkHiEff", fixes = calculables.jet.xcStrip(obj['jet']) ),
            supy.calculables.other.abbreviation( "nVertexRatio", "nvr" ),
            supy.calculables.other.abbreviation('muonTriggerWeightPF','tw'),
            ]
        calcs += supy.calculables.fromCollections(calculables.top,[('genTop',""),('fitTop',"")])
        return calcs
    ########################################################################################

    def listOfSteps(self, pars) :
        obj = pars["objects"]
        lepton = obj[pars["lepton"]["name"]]
        lPtMin = pars["lepton"]["ptMin"]
        lEtaMax = pars["lepton"]["etaMax"]
        lIso = ('%s'+pars['lepton']['iso']+'%s') % lepton
        lIsoIndices = ("%s"+pars["selection"]["lIso"]["indices"]+"%s")%lepton
        topTag = pars['tag'].replace("Wlv","top").replace("QCD","top")
        bVar = ("%s"+pars["bVar"]+"%s")%calculables.jet.xcStrip(obj["jet"])
        
        ssteps = supy.steps
        
        return (
            [ssteps.printer.progressPrinter()
             , ssteps.histos.value("genpthat",200,0,1000,xtitle="#hat{p_{T}} (GeV)").onlySim()
             
             ####################################
             , ssteps.filters.label('data cleanup'),
             ssteps.filters.multiplicity("vertexIndices",min=1),
             ssteps.filters.value('physicsDeclared',min=1).onlyData(),
             ssteps.filters.value('hbheNoiseFilterResult',min=1).onlyData(),
             steps.trigger.l1Filter("L1Tech_BPTX_plus_AND_minus.v0").onlyData(),
             steps.filters.monster()

             ####################################

             , ssteps.filters.label('nvertex reweighting')
             , self.ratio(pars)
             , ssteps.histos.value(obj["sumPt"],50,0,1500)
             , ssteps.histos.value("rho",100,0,40)

             , ssteps.filters.label('trigger reweighting')
             , self.triggerWeight(pars, [ss.weightedName for ss in self.data(pars)])
             , steps.trigger.lowestUnPrescaledTriggerHistogrammer()
             
             ####################################
             , ssteps.filters.label('cross-cleaning'),
             steps.jet.forwardFailedJetVeto( obj["jet"], ptAbove = 50, etaAbove = 3.5),
             ssteps.filters.multiplicity( max=0, var = "%sIndicesOther%s"%obj['jet']),
             ssteps.filters.multiplicity( max=0, var = "%sIndicesOther%s"%obj['muon']),
             ssteps.filters.multiplicity( max=0, var = "%sIndicesUnmatched%s"%obj['photon']),
             ssteps.filters.multiplicity( max=0, var = "%sIndicesUnmatched%s"%obj['electron']),
             steps.jet.uniquelyMatchedNonisoMuons(obj["jet"])
             
             ####################################
             , ssteps.filters.label('selection'),
             ssteps.filters.OR( [ssteps.filters.multiplicity( max=1, var = '%sIndicesAnyIso%s'%lepton),
                                 ssteps.filters.value( lIso, min=0.25, indices = "%sIndicesAnyIsoIsoOrder%s"%lepton, index = 1)]),
             ssteps.filters.multiplicity( max=0, var = "%sIndices%s"%obj['photon']),
             ssteps.filters.multiplicity( max=0, var = "%sIndices%s"%obj['electron'])

             , ssteps.histos.value("%sTriggeringPt%s"%lepton, 200,0,200),
             ssteps.filters.value("%sTriggeringPt%s"%lepton, min = lPtMin)

             , ssteps.histos.multiplicity("%sIndices%s"%obj["jet"]),
             ssteps.filters.multiplicity("%sIndices%s"%obj["jet"], **pars["nJets"])

             , ssteps.histos.pt("mixedSumP4",100,0,300),
             ssteps.filters.pt("mixedSumP4",min=20)

             , ssteps.histos.value( lIso, 55,0,1.1, indices = "%sIndicesAnyIsoIsoOrder%s"%lepton, index=0),
             ssteps.filters.value( lIso, max = 1.0, indices = "%sIndicesAnyIsoIsoOrder%s"%lepton, index = 0),
             ssteps.filters.multiplicity("%sIndices%s"%lepton, min=pars["selection"]["lIso"]["N"],max=pars["selection"]["lIso"]["N"]),
             ssteps.filters.pt("%sP4%s"%lepton, min = lPtMin, indices = lIsoIndices, index = 0),
             ssteps.filters.absEta("%sP4%s"%lepton, max = lEtaMax, indices = lIsoIndices, index = 0)

             , ssteps.histos.value(bVar, 60,0,15, indices = "%sIndicesBtagged%s"%obj["jet"], index = 0)
             , ssteps.histos.value(bVar, 60,0,15, indices = "%sIndicesBtagged%s"%obj["jet"], index = 1)
             , ssteps.histos.value(bVar, 60,0,15, indices = "%sIndicesBtagged%s"%obj["jet"], index = 2)
             , calculables.jet.ProbabilityGivenBQN(obj["jet"], pars['bVar'], binning=(64,-1,15), samples = pars['topBsamples'], tag = topTag)
             , ssteps.histos.value("TopRatherThanWProbability", 100,0,1),
             ssteps.filters.value(bVar, indices = "%sIndicesBtagged%s"%obj["jet"], index = 1, min = 0.0),
             ssteps.filters.value(bVar, indices = "%sIndicesBtagged%s"%obj["jet"], **pars["selection"]["bCut"])
             
             , ssteps.filters.label('top reco'),
             ssteps.filters.multiplicity("TopReconstruction",min=1)
             
             ####################################
             , ssteps.filters.label("selection complete")

             , ssteps.histos.multiplicity("%sIndices%s"%obj["jet"])
             , ssteps.histos.value("%sM3%s"%obj['jet'], 20,0,800)
             , ssteps.histos.value("fitTopRadiativeCoherence", 100,-1,1)
             
             ####################################
             , ssteps.filters.label('discriminants')
             , ssteps.histos.value("KarlsruheDiscriminant", 28, -320, 800 )
             , ssteps.histos.value("TriDiscriminant",50,-1,1)
             #, ssteps.filters.label('qq:gg'),   self.discriminantQQgg(pars)
             , ssteps.filters.label('top:W'),   self.discriminantTopW(pars)
             , ssteps.filters.label('top:QCD'), self.discriminantTopQCD(pars)
             , ssteps.filters.label('W:QCD'),   self.discriminantWQCD(pars)
             , calculables.gen.qDirProbPlus('fitTopSumP4Eta', 10, 'top_muon_pf', 'tt_tauola_fj.wTopAsymP00.tw.nvr', path = self.globalStem)
             
             , ssteps.filters.label('signal distributions')
             , steps.top.Asymmetry(('fitTop',''), bins = 640)
             , steps.top.Spin(('fitTop',''))
             
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
        return calculables.trigger.TriggerWeight( samples,
                                                  unreliable = pars['unreliable'],
                                                  **dict(zip(['triggers','thresholds'], zip(*pars['lepton']['triggers']))))
    
    @staticmethod
    def ratio(pars) : 
        return supy.calculables.other.Ratio("nVertex", binning = (20,-0.5,19.5), thisSample = pars['baseSample'],
                                            target = ("SingleMu",[]), groups = [('qcd_mg',[]),('qcd_py6',[]),('w_jets_fj_mg',[]),('dyj_ll_mg',[]),
                                                                                ('tt_tauola_fj',['tt_tauola_fj%s.tw.nvr'%s for s in ['',
                                                                                                                                     '.wNonQQbar',
                                                                                                                                     '.wTopAsymP00']])])
    @staticmethod
    def discriminantQQgg(pars) :
        return supy.calculables.other.Discriminant( fixes = ("","TopQqQg"),
                                                    left = {"pre":"qq", "tag":"top_muon_pf", "samples":['tt_tauola_fj.wNonQQbar.tw.nvr']},
                                                    right = {"pre":"qg", "tag":"top_muon_pf", "samples":['tt_tauola_fj.wTopAsymP00.tw.nvr']},
                                                    dists = {"fitTopPtOverSumPt" : (20,0,1),
                                                             "fitTopCosThetaDaggerTT" : (40,-1,1),
                                                             "fitTopSumP4AbsEta" : (21,0,7),
                                                             "%sIndices%s.size"%pars['objects']["jet"] : (10,-0.5,9.5)
                                                             },
                                                    correlations = pars['discriminant2DPlots'],
                                                    bins = 14)
    @staticmethod
    def discriminantTopW(pars) :
        return supy.calculables.other.Discriminant( fixes = ("","TopW"),
                                                    left = {"pre":"w_jets_fj_mg", "tag":"top_muon_pf", "samples":[]},
                                                    right = {"pre":"tt_tauola_fj", "tag":"top_muon_pf", "samples": ['tt_tauola_fj.%s.tw.nvr'%s for s in ['wNonQQbar','wTopAsymP00']]},
                                                    correlations = pars['discriminant2DPlots'],
                                                    dists = {"%sKt%s"%pars['objects']["jet"] : (25,0,150),
                                                             "%sB0pt%s"%pars['objects']["jet"] : (30,0,300),
                                                             "%s3absEta%s"%pars['objects']["jet"] : (20,0,4),
                                                             "fitTopHadChi2"     : (20,0,100),
                                                             "mixedSumP4.pt"     : (30,0,180),
                                                             #"fitTopLeptonPt"    : (30,0,180),  # not so powerful?
                                                             "fitTopDeltaPhiLNu" : (20,0,math.pi),
                                                             "TopRatherThanWProbability" : (20,0,1),
                                                             })
    @staticmethod
    def discriminantTopQCD(pars) :
        return supy.calculables.other.Discriminant( fixes = ("","TopQCD"),
                                                    left = {"pre":"SingleMu", "tag":"QCD_muon_pf", "samples":[]},
                                                    right = {"pre":"tt_tauola_fj", "tag":"top_muon_pf", "samples": ['tt_tauola_fj.%s.tw.nvr'%s for s in ['wNonQQbar','wTopAsymP00']]},
                                                    correlations = pars['discriminant2DPlots'],
                                                    dists = {"%sKt%s"%pars['objects']["jet"] : (25,0,150),
                                                             "%sB0pt%s"%pars['objects']["jet"] : (30,0,300),
                                                             "%s3absEta%s"%pars['objects']["jet"] : (20,0,4),
                                                             "%sMt%s"%pars['objects']['muon']+"mixedSumP4" : (30,0,180),
                                                             "%sDeltaPhiB01%s"%pars['objects']["jet"] : (20,0,math.pi),
                                                             #"mixedSumP4.pt"     : (30,0,180),
                                                             #"fitTopLeptonPt"    : (30,0,180),
                                                             #"fitTopDeltaPhiLNu" : (20,0,math.pi),
                                                             })
    @staticmethod
    def discriminantWQCD(pars) :
        return supy.calculables.other.Discriminant( fixes = ("","WQCD"),
                                                    left = {"pre":"w_jets_fj_mg", "tag":"top_muon_pf", "samples":[]},
                                                    right = {"pre":"SingleMu", "tag":"QCD_muon_pf", "samples":[]},
                                                    correlations = pars['discriminant2DPlots'],
                                                    dists = {"%sB0pt%s"%pars['objects']["jet"] : (30,0,300),
                                                             "%sMt%s"%pars['objects']['muon']+"mixedSumP4" : (30,0,180),
                                                             "%sDeltaPhiB01%s"%pars['objects']["jet"] : (20,0,math.pi),
                                                             "fitTopCosHelicityThetaL": (20,-1,1),
                                                             })
    
    ########################################################################################
    def concludeAll(self) :
        self.rowcolors = 2*[13] + 2*[45]
        super(topAsymm,self).concludeAll()
        #self.meldWpartitions()
        #self.meldQCDpartitions()
        self.meldScale()
        #self.dilutions()
        #self.measureQQbarComponent()
        self.plotMeldScale()
        #self.ensembleTest()
        #self.PEcurves()

    def conclude(self,pars) :
        org = self.organizer(pars)
        org.mergeSamples(targetSpec = {"name":"Data 2011", "color":r.kBlack, "markerStyle":20}, allWithPrefix="SingleMu")
        org.mergeSamples(targetSpec = {"name":"qcd_py6", "color":r.kBlue}, allWithPrefix="qcd_py6")
        org.mergeSamples(targetSpec = {"name":"t#bar{t}", "color":r.kViolet}, sources=["tt_tauola_fj.wNonQQbar.tw.nvr","tt_tauola_fj.wTopAsymP00.tw.nvr"], keepSources = True)
        org.mergeSamples(targetSpec = {"name":"t#bar{t}.q#bar{q}.N30", "color":r.kRed}, sources = ["tt_tauola_fj.wTopAsymN30.tw.nvr","tt_tauola_fj.wNonQQbar.tw.nvr"][:1])
        org.mergeSamples(targetSpec = {"name":"t#bar{t}.q#bar{q}.P30", "color":r.kGreen}, sources = ["tt_tauola_fj.wTopAsymP30.tw.nvr","tt_tauola_fj.wNonQQbar.tw.nvr"][:1])
        org.mergeSamples(targetSpec = {"name":"standard_model", "color":r.kGreen+2}, sources = ["qcd_py6","t#bar{t}","w_jets_fj_mg.tw.nvr"], keepSources = True)
        for ss in filter(lambda ss: 'tt_tauola' in ss['name'],org.samples) : org.drop(ss['name'])

        orgpdf = copy.deepcopy(org)
        orgpdf.scale( toPdf = True )
        org.scale( lumiToUseInAbsenceOfData = 1.1e3 )

        names = [ss["name"] for ss in org.samples]
        kwargs = {"detailedCalculables": False,
                  "blackList":["lumiHisto","xsHisto","nJobsHisto"],
                  "samplesForRatios" : next(iter(filter(lambda x: x[0] in names and x[1] in names, [("Data 2011","standard_model")])), ("","")),
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

    def meldWpartitions(self) :
        samples = {"top_muon_pf" : ["w_"],
                   "Wlv_muon_pf" : ["w_","SingleMu"],
                   "QCD_muon_pf" : []}
        organizers = [supy.organizer(tag, [s for s in self.sampleSpecs(tag) if any(item in s['name'] for item in samples[tag])])
                      for tag in [p['tag'] for p in self.readyConfs]]
        if len(organizers)<2 : return
        for org in organizers :
            org.mergeSamples(targetSpec = {"name":"Data 2011", "color":r.kBlack, "markerStyle":20}, allWithPrefix="SingleMu")
            org.mergeSamples(targetSpec = {"name":"w_mg", "color":r.kRed if "Wlv" in org.tag else r.kBlue, "markerStyle": 22}, sources = ["w_jets_fj_mg.tw.nvr"])
            org.scale(toPdf=True)

        melded = supy.organizer.meld("wpartitions",filter(lambda o: o.samples, organizers))
        pl = supy.plotter(melded,
                          pdfFileName = self.pdfFileName(melded.tag),
                          doLog = False,
                          blackList = ["lumiHisto","xsHisto","nJobsHisto"],
                          rowColors = self.rowcolors,
                          rowCycle = 100,
                          omit2D = True,
                          ).plotAll()
        
    def meldQCDpartitions(self) :
        samples = {"top_muon_pf" : ["qcd_py6fjmu"],
                   "Wlv_muon_pf" : [],
                   "QCD_muon_pf" : ["qcd_py6fjmu","SingleMu"]}
        organizers = [supy.organizer(tag, [s for s in self.sampleSpecs(tag) if any(item in s['name'] for item in samples[tag])])
                      for tag in [p['tag'] for p in self.readyConfs]]
        if len(organizers)<2 : return
        for org in organizers :
            org.mergeSamples(targetSpec = {"name":"Data 2011", "color":r.kBlack, "markerStyle":20}, allWithPrefix="SingleMu")
            org.mergeSamples(targetSpec = {"name":"qcd_py6mu", "color":r.kRed if "QCD" in org.tag else r.kBlue, "markerStyle": 22}, allWithPrefix="qcd_py6fjmu")
            org.scale(toPdf=True)

        melded = supy.organizer.meld("qcdPartitions",filter(lambda o: o.samples, organizers))
        pl = supy.plotter(melded,
                          pdfFileName = self.pdfFileName(melded.tag),
                          doLog = False,
                          blackList = ["lumiHisto","xsHisto","nJobsHisto"],
                          rowColors = self.rowcolors,
                          rowCycle = 100,
                          omit2D = True,
                          ).plotAll()
        
    def plotMeldScale(self) :
        if not hasattr(self,"orgMelded") : print "run meldScale() before plotMeldScale()"; return
        melded = copy.deepcopy(self.orgMelded)
        for ss in filter(lambda ss: 'tt_tauola_fj' in ss['name'], melded.samples) : melded.drop(ss['name'])
        pl = supy.plotter(melded, pdfFileName = self.pdfFileName(melded.tag),
                          doLog = False,
                          blackList = ["lumiHisto","xsHisto","nJobsHisto"],
                          rowColors = self.rowcolors,
                          samplesForRatios = ("top.Data 2011","S.M."),
                          sampleLabelsForRatios = ('data','s.m.'),
                          rowCycle = 100,
                          omit2D = True,
                          ).plotAll()
        
    def meldScale(self) :
        meldSamples = {"top_muon_pf" : ["SingleMu","tt_tauola_fj","w_jets","dyj_ll_mg"],
                       #"Wlv_muon_pf" : ["w_jets"],
                       "QCD_muon_pf" : ["SingleMu"]}
        
        organizers = [supy.organizer(tag, [s for s in self.sampleSpecs(tag) if any(item in s['name'] for item in meldSamples[tag])])
                      for tag in [p['tag'] for p in self.readyConfs if p["tag"] in meldSamples]]
        if len(organizers) < len(meldSamples) : return
        for org in organizers :
            org.mergeSamples(targetSpec = {"name":"t#bar{t}", "color":r.kViolet}, sources=["tt_tauola_fj.wNonQQbar.tw.nvr","tt_tauola_fj.wTopAsymP00.tw.nvr"], keepSources = True)
            org.mergeSamples(targetSpec = {"name":"w_jets", "color":r.kRed}, allWithPrefix = "w_jets")
            org.mergeSamples(targetSpec = {"name":"dy", "color":r.kGray}, allWithPrefix = "dyj_ll_mg")
            org.mergeSamples(targetSpec = {"name":"Data 2011",
                                           "color":r.kBlue if "QCD_" in org.tag else r.kBlack,
                                           "markerStyle":(20 if "top" in org.tag else 1)}, allWithPrefix="SingleMu")

        self.orgMelded = supy.organizer.meld(organizers = organizers)

        templateSamples = ['top.t#bar{t}','top.w_jets','QCD.Data 2011']

        def measureFractions(dist) :
            before = next(self.orgMelded.indicesOfStep("label","selection complete"))
            distTup = self.orgMelded.steps[next(iter(filter(lambda i: before<i, self.orgMelded.indicesOfStepsWithKey(dist))))][dist]
            #distTup = self.orgMelded.steps[next(self.orgMelded.indicesOfStepsWithKey(dist))][dist]

            templates = [None] * len(templateSamples)
            for ss,hist in zip(self.orgMelded.samples,distTup) :
                contents = [hist.GetBinContent(i) for i in range(hist.GetNbinsX()+2)]
                if ss['name'] == "top.Data 2011" :
                    observed = contents
                    nEventsObserved = sum(observed)
                elif ss['name'] in templateSamples :
                    templates[templateSamples.index(ss['name'])] = contents
                else : pass
        
            from supy.utils.fractions import componentSolver,drawComponentSolver
            cs = componentSolver(observed, templates, 1e4)
            csCanvas = drawComponentSolver(cs)
            name = "measuredFractions_%s"%dist
            supy.utils.tCanvasPrintPdf(csCanvas[0], "%s/%s"%(self.globalStem,name))
            with open(self.globalStem+"/%s.txt"%name,"w") as file :
                print >> file, cs
                print >> file, cs.components
            return distTup,cs

        distTup,cs = map(measureFractions,["KarlsruheDiscriminant","TriDiscriminant"])[-1]

        fractions = dict(zip(templateSamples,cs.fractions))        
        for iSample,ss in enumerate(self.orgMelded.samples) :
            if ss['name'] in fractions : self.orgMelded.scaleOneRaw(iSample, fractions[ss['name']] * sum(cs.observed) / distTup[iSample].Integral(0,distTup[iSample].GetNbinsX()+1))
        self.orgMelded.mergeSamples(targetSpec = {"name":"S.M.", "color":r.kGreen+2}, sources = ['top.w_jets','top.t#bar{t}','QCD.Data 2011'], keepSources = True, force = True)

    def dilutions(self) :
        import itertools
        fileName = '%s/dilutions.txt'%self.globalStem
        with open(fileName, "w") as file :
            names = [ss['name'] for ss in self.orgMelded.samples]
            iSamples = [names.index(n) for n in ['top.t#bar{t}','top.w_jets','QCD.Data 2011','top.tt_tauola_fj.wNonQQbar.tw.nvr','top.tt_tauola_fj.wTopAsymP00.tw.nvr']]
            for i,iS in enumerate(iSamples) : print >> file, i,names[iS]
            print >> file
            print >> file, ''.rjust(40), ''.join(("[%d,%d]"%pair).rjust(10) for pair in itertools.combinations(range(len(iSamples)), 2))
            for step in self.orgMelded.steps :
                if not any("Discriminant" in item for item in step.nameTitle) : continue
                print >> file
                print >> file
                print >> file, supy.utils.hyphens
                print >> file, step.name, step.title
                print >> file, supy.utils.hyphens
                print >> file
                for hname,hists in step.iteritems() :
                    if issubclass(type(next(iter(filter(None,hists)))),r.TH2) : continue
                    aHists = [[hists[i].GetBinContent(j) for j in range(0,hists[i].GetNbinsX()+2)] for i in iSamples]
                    print >> file, hname.rjust(40), ''.join([(("%.3f"%round(supy.utils.dilution(A,B),3))).rjust(10) for A,B in itertools.combinations(aHists,2)])
        print "Output file: %s"%fileName

    def PEcurves(self) :
        if not hasattr(self, 'orgMelded') : return
        specs = ([{'var' : "ak5JetPFNTrkHiEffPat[i[%d]]:xcak5JetPFIndicesBtaggedPat"%bIndex, 'left':True, 'right':False} for bIndex in [0,1,2]] +
                 [{'var' : "TopRatherThanWProbability",                                      'left':True, 'right':False},
                  {'var' : "TriDiscriminant",                                                'left':True, 'right':True}])
        pes = {}
        for spec in specs :
            dists = dict(zip([ss['name'] for ss in self.orgMelded.samples ],
                             self.orgMelded.steps[next(self.orgMelded.indicesOfStepsWithKey(spec['var']))][spec['var']] ) )
            contours = supy.utils.optimizationContours( [dists['top.t#bar{t}']],
                                                        [dists[s] for s in ['QCD.Data 2011','top.w_jets']],
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

    def measureQQbarComponent(self) :
        dist = "DiscriminantTopQqQg"
        dists = dict(zip([ss['name'] for ss in self.orgMelded.samples ],
                         self.orgMelded.steps[next(self.orgMelded.indicesOfStepsWithKey(dist))][dist] ) )
        def contents(name) : return np.array([dists[name].GetBinContent(i) for i in range(dists[name].GetNbinsX()+2)])

        from supy.utils.fractions import componentSolver, drawComponentSolver
        cs = componentSolver(observed = contents('top.Data 2011'),
                             components = [ contents('top.tt_tauola_fj.wTopAsymP00.tw.nvr'), contents('top.tt_tauola_fj.wNonQQbar.tw.nvr')],
                             ensembleSize = 1e4,
                             base = contents('top.w_jets') + contents('QCD.Data 2011')
                             )
        csCanvas = drawComponentSolver(cs)
        supy.utils.tCanvasPrintPdf(csCanvas[0], "%s/measuredQQFractions"%self.globalStem)
        with open(self.globalStem+"/measuredQQFractions.txt","w") as file :  print >> file, cs


    def templates(self, iStep, dist, qqFrac) :
        if not hasattr(self,'orgMelded') : print 'run meldScale() before asking for templates()'; return
        topQQs = [s['name'] for s in self.orgMelded.samples if 'wTopAsym' in s['name']]
        asymm = [eval(name.replace("top.tt_tauola_fj.wTopAsym","").replace(".tw.nvr","").replace("P",".").replace("N","-.")) for name in topQQs]
        distTup = self.orgMelded.steps[iStep][dist]
        edges = supy.utils.edgesRebinned( distTup[ self.orgMelded.indexOfSampleWithName("S.M.") ], targetUncRel = 0.015, offset = 2 )

        def nparray(name, scaleToN = None) :
            hist_orig = distTup[ self.orgMelded.indexOfSampleWithName(name) ]
            hist = hist_orig.Rebin(len(edges)-1, "%s_rebinned"%hist_orig.GetName(), edges)
            bins = np.array([hist.GetBinContent(j) for j in range(hist.GetNbinsX()+2)])
            if scaleToN : bins *= (scaleToN / sum(bins))
            return bins

        nTT = sum(nparray('top.t#bar{t}'))
        observed = nparray('top.Data 2011')
        base = ( nparray('QCD.Data 2011') +
                 nparray('top.w_jets') +
                 nparray('top.tt_tauola_fj.wNonQQbar.tw.nvr', scaleToN = (1-qqFrac) * nTT )
                 )
        templates = [base +  nparray(qqtt, qqFrac*nTT ) for qqtt in topQQs]
        return zip(asymm, templates), observed
    

    def ensembleFileName(self, iStep, dist, qqFrac, suffix = '.pickleData') :
        return "%s/ensembles/%d_%s_%.3f%s"%(self.globalStem,iStep,dist,qqFrac,suffix)

    def ensembleTest(self) :
        qqFracs = sorted([0.10, 0.12, 0.15, 0.20, 0.25, 0.30, 0.40, 0.60, 1.0])
        dists = ['lHadtDeltaY',
                 'ttbarDeltaAbsY',
                 'leptonRelativeY',
                 'ttbarSignedDeltaY'
                ]
        args = sum([[(iStep, dist, qqFrac) for iStep in list(self.orgMelded.indicesOfStepsWithKey(dist))[:None] for qqFrac in qqFracs] for dist in dists],[])
        supy.utils.operateOnListUsingQueue(6, supy.utils.qWorker(self.pickleEnsemble), args)
        ensembles = dict([(arg,supy.utils.readPickle(self.ensembleFileName(*arg))) for arg in args])

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
                
    def pickleEnsemble(self, iStep, dist, qqFrac ) :
        supy.utils.mkdir(self.globalStem+'/ensembles')
        templates,observed = self.templates(iStep, dist, qqFrac)
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
