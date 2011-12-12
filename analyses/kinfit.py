import supy,steps,calculables,samples
import os,math,copy,ROOT as r, numpy as np

class kinfit(supy.analysis) :

    def parameters(self) :

        objects = {
            'label'       :[               'pf' ,               'pat' ],
            'muonsInJets' :[               True ,               False ],
            'jet'         :[('ak5JetPF','Pat'),   ('xcak5Jet','Pat')],
            'muon'        :[       ('muon','PF'),       ('muon','Pat')],
            'electron'    :[   ('electron','PF'),   ('electron','Pat')],
            'photon'      :[    ('photon','Pat'),     ('photon','Pat')],
            'met'         :[          'metP4PF' ,    'metP4AK5TypeII' ],
            'sumP4'       :[          'pfSumP4' ,           'xcSumP4' ],
            'sumPt'       :[       'metSumEtPF' ,           'xcSumPt' ],
            }

        return { "bVar" : "NTrkHiEff", # "TrkCountingHighEffBJetTags"
                 "objects": dict((key,val[0]) for key,val in objects.items()),
                 "lepton" : {'name':'muon'},
                 "topBsamples": ("tt_tauola_fj",[]),
                 #"topBsamples": ("ttj_mg",[]),
                 #"lumi" : 2e3,
                 "lumi" : 5e2,
                 "nFilesMax" : self.vary({'test':1, 'full':200})
                 }

    ########################################################################################

    def listOfSampleDictionaries(self) : return [samples.top]

    def listOfSamples(self,pars) :  return supy.samples.specify( names = pars['topBsamples'][0], effectiveLumi = pars['lumi'] , nFilesMax = pars['nFilesMax'])

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

            calculables.xclean.IndicesUnmatched(collection = obj["photon"], xcjets = obj["jet"], DR = 0.5),
            calculables.xclean.IndicesUnmatched(collection = obj["electron"], xcjets = obj["jet"], DR = 0.5),
            calculables.xclean.xcJet(obj["jet"], applyResidualCorrectionsToData = False,
                                     gamma    = obj["photon"],      gammaDR = 0.5,
                                     electron = obj["electron"], electronDR = 0.5,
                                     muon     = obj["muon"],         muonDR = 0.5, correctForMuons = not obj["muonsInJets"]),

            calculables.top.mixedSumP4(transverse = obj["met"], longitudinal = obj["sumP4"]),
            calculables.top.SemileptonicTopIndex(lepton),
            calculables.top.fitTopLeptonCharge(lepton),
            calculables.top.TopReconstruction(lepton,obj["jet"],"mixedSumP4"),
            calculables.top.IndicesGenTopPQHL(obj["jet"]),
            calculables.top.genTopRecoIndex(),
            calculables.top.genTopSemiLeptonicWithinAcceptance(jetPtMin = 20, jetAbsEtaMax = 3.5, lepPtMin = 20, lepAbsEtaMax = 2.1),
            #calculables.top.genTopSemiLeptonicAccepted(obj['jet']),

            calculables.top.TopComboQQBBLikelihood(pars['objects']['jet'], pars['bVar']),
            calculables.top.OtherJetsLikelihood(pars['objects']['jet'], pars['bVar']),

            calculables.other.Covariance(('met','PF')),

            supy.calculables.other.abbreviation( "TrkCountingHighEffBJetTags", "NTrkHiEff", fixes = calculables.jet.xcStrip(obj['jet']) ),
            ]
        calcs += supy.calculables.fromCollections(calculables.top,[('genTop',""),('fitTop',"")])
        return calcs
    ########################################################################################

    def listOfSteps(self, pars) :

        
        obj = pars["objects"]

        el = obj['electron']
        mu = obj['muon']
        jet = obj['jet']
        
        bVar = ("%s"+pars["bVar"]+"%s")%calculables.jet.xcStrip(obj["jet"])
        
        ssteps = supy.steps
        
        return (
            [ssteps.printer.progressPrinter(),
             ssteps.filters.value("genTopSemiMu", min=1),
             ssteps.filters.value("genTopSemiLeptonicWithinAcceptance", min=1)
             
             , ssteps.filters.label('selection'),

             ssteps.filters.pt(obj['met'],min=20),
             ssteps.filters.multiplicity("%sIndices%s"%el, max=0       ),
             ssteps.filters.multiplicity("%sIndices%s"%mu, max=1, min=1),
             ssteps.filters.absEta("%sP4%s"%mu, max = 2.1,  indices = "%sIndices%s"%mu, index = 0),
             ssteps.filters.pt(    "%sP4%s"%mu, min = 20.0, indices = "%sIndices%s"%mu, index = 0),
             ssteps.filters.multiplicity("%sIndices%s"%obj["jet"], min = 4 )

             , calculables.jet.ProbabilityGivenBQN(obj["jet"], pars['bVar'], binning=(64,-1,15), samples = pars['topBsamples']),
             ssteps.filters.value(bVar, indices = "%sIndicesBtagged%s"%obj["jet"], index=1, min=2.0)

             , ssteps.filters.label('top reco'),
             ssteps.filters.multiplicity("TopReconstruction",min=1),

             ssteps.filters.minimum("%sIndicesGenTopPQHL%s"%obj['jet'], 0),
             ssteps.filters.unique("%sIndicesGenTopPQHL%s"%obj['jet'])

             #, steps.top.jetPrinter(obj['jet'])
             
             , ssteps.filters.label("fitTop kinFitLook")          , steps.top.kinFitLook("fitTopRecoIndex")
             , ssteps.filters.label("combinatorial filtering")    , steps.top.combinatorialFiltering(obj['jet'])
             , ssteps.filters.label("combinatorial frequency")    , steps.top.combinatorialFrequency(obj['jet'])
             , ssteps.filters.label("combinatorics look")         , steps.top.combinatoricsLook('genTopRecoIndex',obj['jet'])
             , ssteps.filters.label("combinatorial bg")           , steps.top.combinatorialBG(obj['jet'])
             ])

    ########################################################################################

    def conclude(self,pars) :
        org = self.organizer(pars, verbose = True )
        org.scale( toPdf = True )
        #org.scale( lumiToUseInAbsenceOfData = pars['lumi'] )

        kwargs = {"detailedCalculables": False,
                  "blackList":["lumiHisto","xsHisto","nJobsHisto"],
                  "rowColors" : 2*[13] + 2*[45],
                  "rowCycle" : 100,
                  #"omit2D" : True,
                  }
        
        supy.plotter(org, pdfFileName = self.pdfFileName(org.tag+"_log"),  doLog = True, pegMinimum = 0.01, **kwargs ).plotAll()
        supy.plotter(org, pdfFileName = self.pdfFileName(org.tag+"_nolog"), doLog = False, **kwargs ).plotAll()
        
