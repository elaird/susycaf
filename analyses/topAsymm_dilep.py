import supy,steps,calculables,samples
import os,math,copy,itertools,ROOT as r, numpy as np

class topAsymm_dilep(supy.analysis) :

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
            'name'     : [                  'muon',                'electron'],
            'ptMin'    : [                   15.0 ,                     15.0 ],
            'etaMax'   : [                    2.1 ,                      2.5 ],
            'iso'      : [    'CombinedRelativeIso',    'IsoCombinedAdjusted'],
            'isoNormal': [             {"max":0.10},             {'max':0.09}],
            'isoInvert': [ {"min":0.15, "max":0.55}, {'min':0.14, 'max':0.55}]
            }

        #btag working points: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP
        csvWP = {"CSVL" : 0.244, "CSVM" : 0.679, "CSVT" : 0.898 }
        bCut = {"normal"   : {"index":0, "min":csvWP['CSVT']},
                "inverted" : {"index":0, "min":csvWP['CSVL'], "max":csvWP['CSVM']}}

        return { "vary" : ['toptype'],
                 "nJets" :  {"min":2,"max":None},
                 "bVar" : "CSV", # "Combined Secondary Vertex"
                 "objects" : dict((key,val[0]) for key,val in objects.iteritems()),
                 "muon"     : dict((key,val[0]) for key,val in leptons.iteritems()),
                 "electron" : dict((key,val[1]) for key,val in leptons.iteritems()),
                 "bCut"    : bCut['normal'],
                 "toptype" : self.vary({"mg":"mg","ph":"ph"}),
                 }

    ########################################################################################

    def listOfSampleDictionaries(self) : return [getattr(samples,item) for item in ['muon16', 'electron16', 'top16', 'ewk16', 'qcd16']]

    def data(self,pars) :
        return { "muon" : supy.samples.specify( names = ['SingleMu.2011A',
                                                         'SingleMu.2011B']),
                 "electron" : supy.samples.specify( names = ['EleHad.2011A',
                                                             'EleHad.2011B'])
                 }['muon']

    @staticmethod
    def single_top() :
        return ['top_s_ph','top_t_ph','top_tW_ph','tbar_s_ph','tbar_t_ph','tbar_tW_ph']

    def listOfSamples(self,pars) :
        def ewk(eL = None) :
            return supy.samples.specify( names = ["dyj_ll_mg",
                                                  "w2j_mg",
                                                  "w3j_mg",
                                                  "w4j_mg"], effectiveLumi = eL, color = 28 ) if "QCD" not in pars['tag'] else []

        def single_top(eL = None) :
            return supy.samples.specify( names = self.single_top(),
                                         effectiveLumi = eL, color = r.kGray) if "QCD" not in pars['tag'] else []
        
        def ttbar(eL = None) :
            tt = pars['toptype']
            return (supy.samples.specify(names = "ttj_%s"%tt, effectiveLumi = eL, color = r.kBlue, weights = ["wGG"] ) +
                    supy.samples.specify(names = "ttj_%s"%tt, effectiveLumi = eL, color = r.kCyan, weights = [ "wQG"] ) +
                    supy.samples.specify(names = "ttj_%s"%tt, effectiveLumi = eL, color = r.kOrange, weights = [ "wQQ"] )
                    )
        
        return  ( self.data(pars)  + ewk() + ttbar(None) + single_top() )

    ########################################################################################
    def listOfCalculables(self, pars) :
        obj = pars["objects"]
        mu = obj['muon']
        el = obj['electron']
        calcs = sum( [supy.calculables.zeroArgs(module) for module in [calculables, supy.calculables]]+
                     [supy.calculables.fromCollections(getattr(calculables,item), [obj[item]])
                      for item in  ['jet','photon','electron','muon']], [])
        calcs += supy.calculables.fromCollections(calculables.top,[('genTop',""),('fitTop',"")])
        calcs += [
            calculables.jet.IndicesBtagged(obj["jet"],pars["bVar"]),

            calculables.jet.Indices(       obj["jet"],      ptMin = 20, etaMax = 3.1, flagName = "JetIDloose"),
            calculables.jet.Indices(       obj["jet"],      ptMin = 25, etaMax = 2.4, flagName = "JetIDloose", extraName = "triggering"),

            #calculables.electron.Indices_TopPAG(  obj["electron"], ptMin = 30, absEtaMax = 2.5, id = "ID70"),
            calculables.electron.Indices_TopPAG(  obj["electron"], ptMin = 15, absEtaMax = 2.5, id = "ID95"),
            calculables.muon.Indices(             obj["muon"],     ptMin = 15, absEtaMax = 2.1, ID = "ID_TOPPAG",
                                                  isoMax = 0.25, ISO = "CombinedRelativeIso"), #these two are kind of irrelevant, since we use IndicesAnyIsoIsoOrder

            calculables.muon.IndicesAnyIsoIsoOrder( el, pars["electron"]["iso"]),
            calculables.muon.IndicesAnyIsoIsoOrder( mu, pars["muon"]["iso"]),
            calculables.electron.IndicesIsoLoose( obj["electron"], ptMin = 15, absEtaMax = 2.5, iso = "IsoCombinedAdjusted", isoMax = 0.15),
            calculables.muon.IndicesIsoLoose( obj["muon"], ptMin = 10, absEtaMax = 2.5, iso = "CombinedRelativeIso", isoMax = 0.20 ),

            calculables.electron.IsoCombinedAdjusted(obj["electron"], barrelCIso = 0.09, endcapCIso = 0.06 ), # VBTF 0.85
            calculables.xclean.xcJet_DoubleLepton( obj["jet"], leptons1 = mu, leptons2 = el, indices = "IndicesAnyIsoIsoOrder" ),

            calculables.vertex.ID(),
            calculables.vertex.Indices(),

            calculables.gen.genIndicesHardPartons(),
            calculables.top.TopJets( obj['jet'] ),
            calculables.top.mixedSumP4( transverse = obj["met"], longitudinal = obj["sumP4"] ),
            calculables.top.TopComboQQBBLikelihood( pars['bVar'] ),
            calculables.top.OtherJetsLikelihood( pars['bVar'] ),
            calculables.top.TopRatherThanWProbability( priorTop = 0.5 ),
            calculables.top.IndicesGenTopPQHL( obj['jet'] ),
            calculables.top.IndicesGenTopExtra (obj['jet'] ),
            calculables.top.DiLeptonsPM( mu, el),
            calculables.top.TopsReconstructionDiLepton(obj['met']),
            calculables.other.Covariance(('met','PF')),

            calculables.jet.pt( obj['jet'], index = 0, Btagged = True ),
            calculables.jet.pt( obj['jet'], index = 3, Btagged = False ),
            calculables.jet.absEta( obj['jet'], index = 3, Btagged = False),

            supy.calculables.other.size( "Indices".join(obj['jet']) ),
            supy.calculables.other.abbreviation( "TrkCountingHighEffBJetTags", "TCHE", fixes = calculables.jet.xcStrip(obj['jet']) ),
            supy.calculables.other.abbreviation( "CombinedSecondaryVertexBJetTags", "CSV", fixes = calculables.jet.xcStrip(obj['jet']) ),
            supy.calculables.other.abbreviation( "xcak5JetPFCorrectedP4Pat","xcak5JetPFCorrectedP4Pattriggering"),
            ]
        return calcs
    ########################################################################################

    def listOfSteps(self, pars) :
        obj = pars["objects"]
        mu = obj['muon']
        el = obj['electron']
        bVar = pars["bVar"].join(obj["jet"])[2:]
        tt = pars['toptype']
        
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
             steps.filters.monster(),
             ssteps.filters.value('ecalDeadCellTPFilterFlag',min=1),
             steps.jet.failedJetVeto( obj["jet"], ptMin = 20, id = "PFJetIDloose")

             ####################################

             , ssteps.filters.label('selection'),

             ssteps.filters.multiplicity("Indices".join(obj["jet"]), **pars["nJets"]),
             ssteps.filters.multiplicity( min=1, var = "IndicesAnyIso".join(el) ),
             ssteps.filters.multiplicity( min=1, var = "IndicesAnyIso".join(mu) ),
             ssteps.filters.value( pars['electron']['iso'].join(el), indices = "IndicesAnyIsoIsoOrder".join(el), index = 0, **pars['electron']['isoNormal']),
             ssteps.filters.value( pars['muon']['iso'].join(mu), indices = "IndicesAnyIsoIsoOrder".join(mu), index = 0, **pars['muon']['isoNormal']),
             ssteps.filters.value( 'DiLeptonsPM', min = 0)
             
             , ssteps.histos.value( pars['muon']['iso'].join(mu), 55,0,1.1, indices = "IndicesAnyIsoIsoOrder".join(mu), index=0)
             , ssteps.histos.value( pars['electron']['iso'].join(el), 55,0,1.1, indices = "IndicesAnyIsoIsoOrder".join(el), index=0)
             , ssteps.histos.absEta("P4".join(mu), 100,0,4, indices = "IndicesAnyIsoIsoOrder".join(mu), index = 0)
             , ssteps.histos.absEta("P4".join(el), 100,0,4, indices = "IndicesAnyIsoIsoOrder".join(el), index = 0)
             , ssteps.histos.pt("P4".join(mu), 200,0,200, indices = "IndicesAnyIsoIsoOrder".join(mu), index = 0)
             , ssteps.histos.pt("P4".join(el), 200,0,200, indices = "IndicesAnyIsoIsoOrder".join(el), index = 0)
 
             , ssteps.histos.value(bVar, 51,-0.02,1, indices = "IndicesBtagged".join(obj["jet"]), index = 0)
             , ssteps.histos.value(bVar, 51,-0.02,1, indices = "IndicesBtagged".join(obj["jet"]), index = 1),
             ssteps.filters.value(bVar, indices = "IndicesBtagged".join(obj["jet"]), **pars["bCut"]),
             ssteps.filters.multiplicity('TopsReconstructionDiLepton',min=1)
             
             , ssteps.filters.label("selection complete")
             , ssteps.histos.multiplicity('TopsReconstructionDiLepton')
             , ssteps.histos.value('LeptonsSumAbsRapidity',50,0,5)
             , ssteps.histos.value('TopsSumAbsY', 50, 0, 5)
             , ssteps.histos.value('TopsPtOverSumPt', 50, 0, 1)
             , ssteps.histos.value('DiLeptonsPtOverSumPt', 50, 0, 1)
             , ssteps.histos.value('TopsDeltaAbsY', 50, -3,3)
             , ssteps.histos.value('genTopsDDeltaY', 50, -2, 2)
             , ssteps.histos.value('genTopsDDeltaAbsY', 50, -2, 2)

             , steps.top.channelClassification().onlySim()
             #, steps.displayer.ttbar(jets=obj["jet"], met=obj['met'], muons = obj['muon'], electrons = obj['electron'])
             ####################################
             #, steps.top.resolutions('genTopRecoIndex')
             #, steps.top.kinematics('fitTop')
             #, steps.top.resolutions('fitTopRecoIndex')
             ####################################

             , ssteps.histos.multiplicity("Indices".join(obj["jet"]))
             , ssteps.filters.label('object pt')
             , ssteps.histos.pt("P4".join(el), 100,1,201, indices = "IndicesAnyIsoIsoOrder".join(el), index = 0)
             , ssteps.histos.pt("P4".join(mu), 100,1,201, indices = "IndicesAnyIsoIsoOrder".join(mu), index = 0)
             , ssteps.histos.pt("CorrectedP4".join(obj['jet']), 100,1,201, indices = "Indices".join(obj['jet']), index = 0)
             , ssteps.histos.pt("CorrectedP4".join(obj['jet']), 100,1,201, indices = "Indices".join(obj['jet']), index = 1)
             , ssteps.filters.label('object eta')
             , ssteps.histos.absEta("P4".join(el), 100,0,4, indices = "IndicesAnyIsoIsoOrder".join(el), index = 0)
             , ssteps.histos.absEta("P4".join(mu), 100,0,4, indices = "IndicesAnyIsoIsoOrder".join(mu), index = 0)
             , ssteps.histos.absEta("CorrectedP4".join(obj['jet']), 100,0,4, indices = "Indices".join(obj['jet']), index = 0)
             , ssteps.histos.absEta("CorrectedP4".join(obj['jet']), 100,0,4, indices = "Indices".join(obj['jet']), index = 1)
             
             , ssteps.filters.label('signal distributions')
             #, ssteps.histos.symmAnti('genTopCosPhiBoost','genTopCosPhiBoost',100,-1,1).disable(saDisable)
             #, ssteps.histos.symmAnti('genTopDeltaBetazRel','genTopDeltaBetazRel',100,-1,1).disable(saDisable)
             #, ssteps.histos.symmAnti('genTopCosThetaBoostAlt','genTopCosThetaBoostAlt',100,-1,1).disable(saDisable)

             #, ssteps.histos.symmAnti('genTopCosPhiBoost','fitTopCosPhiBoost',100,-1,1)
             #, ssteps.histos.symmAnti('genTopCosThetaBoostAlt','fitTopCosThetaBoostAlt',100,-1,1)
             #, ssteps.histos.symmAnti('genTopDeltaBetazRel','fitTopDeltaBetazRel',100,-1,1)

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
                                                                                                           'wGG','wQG','wQQ']])]).onlySim()
    @staticmethod
    def tridiscriminant(pars) :
        rw = pars['reweights']['abbr']
        lname = pars['lepton']['name']
        tt = pars['toptype']
        lumi = 5008 # FIXME HACK !!!
        tops = ['ttj_%s.%s.tw.%s'%(tt,s,rw) for s in ['wGG','wQG','wQQ']]
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
                                                       pi43 = {"pre":"qg", "tag":"top_%s_%s"%(lname,tt), "samples":['ttj_%s.wQG.tw.%s'%(tt,rw)]},
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
#    def concludeAll(self) :
#        self.orgMelded = {}
#        self.sensitivityPoints = {}
#        self.rowcolors = 2*[13] + 2*[45]
#        super(topAsymm,self).concludeAll()
#        for tt,rw,lname in set([(pars['toptype'],pars['reweights']['abbr'],pars['lepton']['name']) for pars in self.readyConfs]) :
#            self.meldScale(rw,lname,tt)
#            self.plotMeldScale(rw,lname,tt)
#            self.PEcurves(rw,lname, tt)
#        for tt,rw in set([(pars['toptype'],pars['reweights']['abbr']) for pars in self.readyConfs]) :
#            self.measureAmplitudes(rw,tt)
#        #self.sensitivity_graphs()
#        #self.grant_proposal_plots()

    def conclude(self,pars) :
        self.rowcolors = 2*[13] + 2*[45]
        tt = pars['toptype']
        org = self.organizer(pars, verbose = True )

#        if pars['lepton']['name']=='muon' :
#            org.mergeSamples(targetSpec = {"name":"SingleMu.2011", "color":r.kBlack, "markerStyle":20}, allWithPrefix="SingleMu")
#        else:
#            org.mergeSamples(targetSpec = {"name":"EleHad.2011", "color":r.kBlack, "markerStyle":20}, allWithPrefix="EleHad")
#            
        org.mergeSamples(targetSpec = {"name":"t#bar{t}", "color":r.kViolet}, sources=["ttj_%s.%s"%(tt,s) for s in ['wQQ','wQG','wGG']], keepSources = True)
        #org.mergeSamples(targetSpec = {"name":"W+jets", "color":28}, sources = ["w%dj_mg.tw.%s"%(n,rw) for n in [2,3,4]])
        #org.mergeSamples(targetSpec = {"name":"DY+jets", "color":r.kYellow}, allWithPrefix="dyj_ll_mg")
        #org.mergeSamples(targetSpec = {"name":"Single top", "color":r.kGray}, sources = ["%s.tw.%s"%(s,rw) for s in self.single_top()])
        #org.mergeSamples(targetSpec = {"name":"Standard Model", "color":r.kGreen+2}, sources = ["multijet","t#bar{t}","W+jets","DY+jets","Single top"], keepSources = True)

        #org.scale( lumiToUseInAbsenceOfData = 5008 )
        org.scale( toPdf = True )

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

    def plotQQGGmeasured(self, org, rw, tt) :
        melded = copy.deepcopy(org)
        assert "FAIL"
        melded.mergeSamples( targetSpec = {"name":"qqggBgStack", "fillColor":r.kRed,"color":r.kRed,"markerStyle":1, "goptions" : "hist"}, sources = ['bg','top.ttj_mg.wQQbar.tw.%s'%rw,'top.ttj_mg.wNonQQbar.tw.%s'%rw], keepSources = True, force = True )
        melded.mergeSamples( targetSpec = {"name":"ggBgStack", "fillColor":r.kGreen,"color":r.kGreen,"markerStyle":1, "goptions" : "hist"}, sources = ['bg','top.ttj_mg.wNonQQbar.tw.%s'%rw], keepSources = True, force = True )
        new = {"ggBgStack" : "gg->t#bar{t}", "qqggBgStack":"q#bar{q}->t#bar{t}","bg":"background","top.Data 2011":"Data 2011"}
        [melded.drop(ss['name']) for ss in melded.samples if ss['name'] not in new]
        spec = {"stepName":"DiscriminantQQgg",
                "stepDesc":None,
                "legendCoords": (0.6, 0.7, 0.82, 0.92),
                "stamp" : True,
                "order" : ["qqggBgStack","ggBgStack","bg","top.Data 2011"]}
        specs = [dict(spec,**{"plotName":"fitTopAbsSumRapidities","newTitle":";|y_{t}+y_{#bar{t}}|;events"})]
        kwargs = {"showStatBox":False, "anMode":True}
        supy.plotter(melded, pdfFileName = self.pdfFileName(org.tag+"_ind"), doLog = False, **kwargs).individualPlots(specs, newSampleNames = new)
        
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
            org.mergeSamples(targetSpec = {"name":"t#bar{t}", "color":r.kViolet}, sources=["ttj_%s.%s.tw.%s"%(tt,s,rw) for s in ['wQQ','wQG','wGG']], keepSources = 'top' in org.tag)
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
            cs = componentSolver(observed, templates, 1e4, base = np.sum(bases, axis=0))
            stuff = drawComponentSolver( cs, mfCanvas, distName = dist,
                                         templateNames = [t.replace("top.ttj_%s.wQQ.tw.%s"%(tt,rw),"q#bar{q}-->t#bar{t}").replace("top.ttj_%s.wQG.tw.%s"%(tt,rw),"qg-->t#bar{t}").replace("top.ttj_%s.wGG.tw.%s"%(tt,rw),"gg-->t#bar{t}").replace("QCD.Data 2011","Multijet").replace("top.W","W+jets").replace('top.',"") for t in  templateSamples])
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
            elif ss['name'] in ['top.ttj_%s.%s.tw.%s'%(tt,s,rw) for s in ['wQQ','wQG','wGG']] :
                f = fractions['top.t#bar{t}']
                n = nTT
            else : continue
            org.scaleOneRaw(iSample, f * sum(cs.observed) / n )

        org.mergeSamples(targetSpec = {"name":"bg", "color":r.kBlack,"fillColor":r.kGray, "markerStyle":1, "goptions":"hist"}, sources = set(baseSamples + templateSamples) - set(['top.t#bar{t}']), keepSources = True, force = True)
        templateSamples = ['top.ttj_%s.%s.tw.%s'%(tt,s,rw) for s in ['wQQ','wQG','wGG']]
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

        samples = ['top.ttj_%s.%s.tw.%s'%(tt,w,rw) for w in ['wQQ','wQG']]
        measure('fitTopCosPhiBoost','genTopCosPhiBoost', samples[1:]) #qg only
        measure('fitTopDeltaBetazRel','genTopDeltaBetazRel', samples, sumTemplates=True)
        measure('fitTopDeltaBetazRel','genTopDeltaBetazRel', samples)
        measure('fitTopCosThetaBoostAlt','genTopCosThetaBoostAlt', samples, sumTemplates=True)
        measure('fitTopCosThetaBoostAlt','genTopCosThetaBoostAlt', samples)
        supy.utils.tCanvasPrintPdf( canvas, maFileName, option = ']')










########################################deprecated
    def sensitivity_graphs(self) :
        qqFracs = {0.07:r.kBlack, 0.14:r.kRed, 0.21:r.kGreen, 0.28:r.kBlue}
        supy.plotter.setupStyle()
        supy.plotter.setupTdrStyle()
        pointkeys = self.sensitivityPoints.keys()
        if not pointkeys : return
        for (step,),keys in itertools.groupby( sorted(self.sensitivityPoints[pointkeys[0]]), lambda key : key[:1] ) :
            canvas = supy.plotter.tcanvas(anMode=True)
            canvas.UseCurrentStyle()
            legend = r.TLegend(0.75,0.7,0.95,0.9)
            legend.SetHeader("#frac{#sigma(q#bar{q}->t#bar{t})}{#sigma(pp=>t#bar{t})}   (%)")
            legend.SetFillStyle(0)
            legend.SetBorderSize(0)
            graphs = []
            for i,key in enumerate(keys) :
                comb_sens = 1./np.sqrt(sum([1./np.square(points[key]) for pointskey,points in self.sensitivityPoints.items()]))
                graph = r.TGraph(len(self.sampleSizeFactor), np.array(self.sampleSizeFactor), comb_sens)
                graph.SetLineColor(qqFracs[key[1]])
                graph.SetMinimum(0)
                graph.SetTitle(";N_{t#bar{t}} / (^{}N_{t#bar{t}} @5fb^{-1});expected precision")
                graph.GetYaxis().SetTitleOffset(1.4)
                graph.Draw('' if i else 'AL')
                legend.AddEntry(graph,"%.2f"%key[1],'l')
                graphs.append(graph)
            legend.Draw()
            fileName = '%s/sensitivity_%d'%(self.globalStem,step) + ".eps"
            canvas.Print(fileName)
            del canvas
            del graphs
            print "Wrote: %s"%fileName
                
    def grant_proposal_plots(self) :
        pars = next((pars for pars in self.readyConfs if "top_" in pars["tag"]),None)
        if not pars : return
        rw = pars['reweights']['abbr']
        tt = pars['toptype']
        names = ["N30","P00","P10","P20","P30"]
        raise "FIXME"
        new = dict([("ttj_mg.wTopAsym%s.tw.%s"%(name,rw), name.replace("P00","  0").replace('P',' +').replace('N','  -')) for name in names])
        org = self.organizer( pars, verbose = True )
        [org.drop(ss['name']) for ss in org.samples if not any(name in ss['name'] for name in names)]
        org.scale( toPdf = True )
        spec = {"stepName":"Asymmetry",
                "stepDesc":"with 41 bins.",
                "legendCoords": (0.18, 0.7, 0.4, 0.92),
                "legendTitle" : "Asymmetry (%)",
                "stamp" : False}
        specs = [dict(spec,**{"plotName":"ttbarDeltaAbsY",
                              "newTitle":";|y_{t}| - |y_{#bar{t}}|;probability density"}),
                 dict(spec,**{"plotName":"ttbarSignExpectation",
                              "newTitle":";#LT sgn(boost) #upoint sgn(#Delta^{}y) #GT;probability density"})
                 ]
        for ss in org.samples : ss['goptions'] = 'hist'
        kwargs = {"showStatBox":False, "anMode":True}
        pl = supy.plotter(org, pdfFileName = self.pdfFileName(org.tag+"_ind"), doLog = False, **kwargs).individualPlots(specs, newSampleNames = new)
        pl = supy.plotter(org, pdfFileName = self.pdfFileName(org.tag+"_ind_log"), doLog = True, pegMinimum=0.001, **kwargs).individualPlots(specs, newSampleNames = new)
