import os
import array

import ROOT as r
import calculables
import samples
import steps
import supy

class muonLook(supy.analysis):
    def parameters(self) :
	    objects = self.vary()
	    fields =                              [        "jet",             "met",            "muon",  "electron",
							   "photon",          "rechit",     "muonsInJets",     "jetId",
							   "jetComp",      "jetIdComp", "muonsInJetsComp",   "metComp",
							   "rechitComp",]

	    objects["caloJet"] = dict(zip(fields, [("xcak5Jet","Pat"),  "metP4TypeIPF", ("muon","Pat"), ("electron","Pat"),
						   ("photon","Pat"),           "Calo" ,           True,       "JetIDloose",
						   ("xcak5JetPF", "Pat"), "JetIDtight",           True,          "metP4PF",
						   "PF",]))

	    objects["pfJet"] = dict(zip(fields, [("xcak5JetPF","Pat"), "metP4TypeIPF", ("muon","PF"),      ("electron","PF"),
						 ("photon","Pat"),          "PF" ,          True,       "JetIDloose",
						 ("xcak5Jet", "Pat"),    "JetIDtight",           True,     "metP4PF",
						 "Calo",]))

	    return {"objects": objects,
		    "nJetsMinMax": self.vary({"ge2j": (2, None),
					      #"le3j": (2, 3),
					      #"ge4j": (4, None),
					      }),
		    "thresholds": self.vary(dict([("375",  (375.0, None, 100.0, 50.0)),
                                                  ("200",  (200.0, 325.0 , 73.3, 36.6)),
						  ])),
		    "etRatherThanPt": True,
		    "lowPtThreshold": 30.0,
		    "lowPtName": "lowPt",
		    "highPtThreshold": 50.0,
		    "highPtName": "highPt",
                    #required to be sorted
		    "triggerList": tuple(["HLT_IsoMu24_eta2p1_v%i" % v for v in range(11,16)]),
		    }

    def calcListJet(self, obj, etRatherThanPt, ptMin,
                    lowPtThreshold, lowPtName,
                    highPtThreshold, highPtName, htThreshold):

        def calcList(jet=None, jetId=None, met=None,
                     photon=None, muon=None, muonsInJets=None, electron=None,
                     **_):
            return [calculables.xclean.xcJet(jet,
                                             gamma=photon,
                                             gammaDR=0.5,
                                             muon=muon,
                                             muonDR=0.5,
                                             correctForMuons = not muonsInJets,
                                             electron=electron,
                                             electronDR=0.5),
                    calculables.jet.Indices(jet,
                                            ptMin=ptMin,
                                            etaMax=3.0,
                                            flagName=jetId,
                                            ),
                    calculables.jet.Indices(jet,
                                            ptMin=lowPtThreshold,
                                            etaMax=3.0,
                                            flagName=jetId,
                                            extraName=lowPtName,
                                            ),
                    calculables.jet.Indices(jet,
                                            ptMin=highPtThreshold,
                                            etaMax=3.0,
                                            flagName=jetId,
                                            extraName=highPtName,
                                            ),
                    calculables.jet.IndicesBtagged2(jet,
                                                    tag="CombinedSecondaryVertexBJetTags",
                                                    threshold=0.679,
                                                    ),

                    calculables.jet.SumP4(jet),
                    calculables.jet.SumP4(jet, extraName=lowPtName),
                    calculables.jet.SumP4(jet, extraName=highPtName),
                    calculables.jet.SumP4(jet, extraName="Btagged2"),
                    calculables.jet.DeltaPhiStar(jet, extraName=lowPtName),
                    calculables.jet.DeltaPhiStar(jet),
                    calculables.jet.MaxEmEnergyFraction(jet),
                    calculables.jet.DeltaPseudoJet(jet, etRatherThanPt),
                    calculables.jet.AlphaT(jet, etRatherThanPt),
                    calculables.jet.AlphaTMet(jet, etRatherThanPt, met),
                    calculables.jet.ra1AlphaTCategory(jet, etRatherThanPt),
                    calculables.jet.MhtOverMet((jet[0], highPtName + jet[1]),
                                               met),
                    calculables.jet.MhtOverMet((jet[0], highPtName + jet[1]),
                                               met = "%sPlus%s%s"%(obj["met"], obj["muon"][0], obj["muon"][1])),
                    calculables.jet.MhtOverMet((jet[0], highPtName + jet[1]),
                                               met = "%sPlus%sIndices%s"%(obj["met"], obj["muon"][0], obj["muon"][1])),
                    calculables.jet.DeadEcalDR(jet,
                                               extraName=lowPtName,
                                               minNXtals=10),
                    calculables.jet.IndicesWithOddMuon(collection=jet, muons=muon),
                    calculables.jet.BTagProbability(collection=jet, selection="muon"),
                    supy.calculables.other.fixedValue("%sFixedHtBin%s" % jet,
                                                      htThreshold),
                    ] + supy.calculables.fromCollections(calculables.jet,
                                                         [jet])
        
        outList = calcList(**obj)
        compList = ["jet", "met", "muonsInJets", "jetId"]
        if all([(item+"Comp" in obj) for item in compList]):
            obj2 = obj.copy()
            for item in compList:
                obj2[item] = obj[item+"Comp"]
            outList += calcList(**obj2)
        return outList

    def calcListOther(self, obj, triggers):
        return [calculables.xclean.IndicesUnmatched(collection=obj["photon"],
                                                    xcjets=obj["jet"], DR=0.5,
                                                    ),
                calculables.xclean.IndicesUnmatched(collection=obj["electron"],
                                                    xcjets=obj["jet"], DR=0.5,
                                                    ),
                calculables.muon.Indices(obj["muon"], isoMax=0.12, ptMin=10,
                                         ISO="PfIsolationR04DeltaBCorrected",
                                         ID="IdPog2012Tight",
                                         #requireIsGlobal = False
                                         ),
                calculables.electron.Indices(obj["electron"], ptMin=10,
                                             flag2012="Veto",
                                             ),
                calculables.photon.Indices(obj["photon"], ptMin=25,
                                           flagName="photonSimpleCutBased2012TightPat",
                                           ),
                calculables.photon.CombinedIsoDR03RhoCorrected(obj["photon"]),
                calculables.other.RecHitSumPt(obj["rechit"]),
                calculables.other.RecHitSumP4(obj["rechit"]),
                calculables.other.metPlusIndices(met = obj["met"], collection = obj["muon"]),
                calculables.other.SumP4(obj["muon"]),
                calculables.vertex.ID(),
                calculables.vertex.Indices(),
                calculables.trigger.lowestUnPrescaledTrigger(triggers),
                calculables.other.Mt(obj["muon"], obj["met"]),
                calculables.other.singleIsolatedTrack(ptMin=10., dzMax=.05, relIso=.1,
                                                      muons=obj["muon"], electrons=obj["electron"] ),
                ]

    def listOfCalculables(self, params):
        obj = params["objects"]
        outList = supy.calculables.zeroArgs(supy.calculables)
        outList += supy.calculables.zeroArgs(calculables)
        outList += supy.calculables.fromCollections(calculables.muon,
                                                    [obj["muon"]],
                                                    )
        outList += supy.calculables.fromCollections(calculables.electron,
                                                    [obj["electron"]],
                                                    )
        outList += supy.calculables.fromCollections(calculables.photon,
                                                    [obj["photon"]],
                                                    )
        outList += self.calcListOther(obj, params["triggerList"])
        outList += self.calcListJet(obj,
                                    params["etRatherThanPt"],
                                    params["thresholds"][3],
                                    params["lowPtThreshold"],
                                    params["lowPtName"],
                                    params["highPtThreshold"],
                                    params["highPtName"],
                                    params["thresholds"][0],
                                    )
        outList += [calculables.gen.genIndices(pdgs=[-5, 5],
                                               label="Status3b",
                                               status=[3],
                                               ),
                    calculables.gen.genIndices(pdgs=[-5, 5],
                                               label="Status3bZDaughters",
                                               status=[3],
                                               motherPdgs=[23],
                                               ),
                    calculables.gen.genIndices(pdgs=[-6, 6],
                                               label="Status3t",
                                               status=[3],
                                               ),
                    calculables.gen.topPtWeight(),
                    ]
        return outList

    def stepsEventCount(self, params, label=""):
        _jet = params["objects"]["jet"]
        _et = "Et" if params["etRatherThanPt"] else "Pt"
        out = []
        if params["signalScan"]:
            if "scanBefore" in label:
                out = [supy.steps.filters.label("scanBefore"),
                       steps.gen.scanHistogrammer(htVar="", befOrAf="Before")
                       ]
            elif "scanAfter" in label:
                htVar = "%sSum%s%s" % (_jet[0], _et, _jet[1]),
                out = [supy.steps.filters.label("scanAfter"),
                       steps.gen.scanHistogrammer(htVar=htVar, befOrAf="After")
                       ]
        return out

    def stepsTrigger(self, params):
        return [
            steps.trigger.lowestUnPrescaledTriggerFilter().onlyData(),
            steps.trigger.hltPrescaleHistogrammer(params["triggerList"]).onlyData(),
            steps.trigger.l1Filter("L1Tech_BPTX_plus_AND_minus.v0").onlyData(),
            steps.trigger.physicsDeclaredFilter().onlyData(),
            ]

    def stepsGenValidation(self):
        return [
            #supy.steps.histos.histogrammer("genpthat",200,0,2000,title=";#hat{p_{T}} (GeV);events / bin").onlySim(),
            supy.steps.histos.histogrammer("genPartonHT", 200, 0, 1000, title=";parton H_{T} (GeV);events / bin").onlySim(),
            ]

    def stepsEvent(self, params):
        return [steps.filters.monster(),
                steps.filters.hbheNoise().onlyData(),
                supy.steps.filters.value("beamHaloCSCTightHaloId", max=0),
                supy.steps.filters.value("trackingFailureFilterFlag", min=1),
                supy.steps.filters.value("hcalLaserEventFilterFlag", min=1).onlyData(),
                supy.steps.filters.value("ecalDeadCellTPFilterFlag", min=1),
                supy.steps.filters.value("eeBadScFilterFlag", min=1),
                supy.steps.filters.value("inconsistentMuonPFCandidateFilterFlag", min=1),
                supy.steps.filters.value("greedyMuonPFCandidateFilterFlag", min=1),
                supy.steps.filters.multiplicity("singleIsolatedTrack", max=0),
                #supy.steps.histos.histogrammer("logErrorTooManySeeds",    2, 0.0, 1.0, title = ";logErrorTooManySeeds;events / bin"),
                #supy.steps.histos.histogrammer("logErrorTooManyClusters", 2, 0.0, 1.0, title = ";logErrorTooManyClusters;events / bin"),
                supy.steps.histos.multiplicity("vertexIndices", max=30),
                supy.steps.filters.multiplicity("vertexIndices", min=1),
                ]

    def stepsSlowLaserFilters(self):
        return [supy.steps.filters.value("hcalLaserEvent2012", min=1).onlyData(),
                supy.steps.filters.value("ecalLaserCalibEvent2012", min=1).onlyData(),
                ]

    def stepsHtLeadingJets(self, params):
        jet = params["objects"]["jet"]
        et = "Et" if params["etRatherThanPt"] else "Pt"

        ##when using full scaling
        #steps.jet.htBinFilter(jet, min = params["htBin"], max = params["htBin"]),
        #steps.jet.jetSelector(jet, params["thresholds"][2], 0),
        #steps.jet.jetSelector(jet, params["thresholds"][2], 1),

        #otherwise
        out = [steps.jet.jetEtaSelector(jet, 2.5, 0),
               steps.jet.jetPtSelector(jet, params["thresholds"][2], 0),
               steps.jet.jetPtSelector(jet, params["thresholds"][2], 1),
               supy.steps.filters.value("%sSum%s%s" % (jet[0], et, jet[1]), min=params["thresholds"][0]),
               ]

        if params["thresholds"][1] is not None:
            out.append(supy.steps.filters.value("%sSum%s%s" % (jet[0], et, jet[1]), max=params["thresholds"][1]))
        return out

    def stepsBtagJets(self, params):
        _jet = params["objects"]["jet"]

        out = [steps.jet.bTagEfficiencyHistogrammer(_jet).onlySim(),
               supy.steps.histos.multiplicity("%sIndicesBtagged2%s" % _jet),
#               supy.steps.filters.multiplicity("%sIndicesBtagged2%s" % _jet, min=params["nBTagJets"][0], max=params["nBTagJets"][1]),
               ]
        return out

    def stepsXclean(self, params):
        _muon = params["objects"]["muon"]
        _met = params["objects"]["met"]
        return [
            #steps.other.iterHistogrammer("ecalDeadTowerTrigPrimP4", 256, 0.0, 128.0, title=";E_{T} of ECAL TP in each dead region (GeV);TPs / bin", funcString="lambda x:x.Et()"),
            supy.steps.filters.multiplicity("%sIndices%s" % params["objects"]["jet"], min=params["nJetsMinMax"][0], max=params["nJetsMinMax"][1]),
            supy.steps.filters.multiplicity("%sIndicesOther%s" % params["objects"]["jet"], max=0),
            supy.steps.filters.multiplicity("%sIndicesWithOddMuon%s" % params["objects"]["jet"], max=0),
            steps.jet.uniquelyMatchedNonisoMuons(params["objects"]["jet"]),
            supy.steps.filters.multiplicity("%sIndices%s" % params["objects"]["electron"], max=0),
            supy.steps.filters.multiplicity("%sIndices%s" % params["objects"]["photon"], max=0),
            supy.steps.filters.multiplicity("%sIndicesUnmatched%s" % params["objects"]["electron"], max = 0),
            supy.steps.filters.multiplicity("%sIndicesUnmatched%s" % params["objects"]["photon"], max = 0),
            steps.other.deadEcalFilter(jets=params["objects"]["jet"], extraName=params["lowPtName"], dR=0.3, dPhiStarCut=0.5),
            supy.steps.filters.multiplicity("%sIndices%s" % _muon, min= 1, max=1),
            supy.steps.filters.pt("%sP4%s"%_muon, min = 30.0, indices = "%sIndices%s"%_muon, index = 0),
            supy.steps.filters.eta("%sP4%s"%_muon, min = -2.1, max = 2.1, indices = "%sIndices%s"%_muon, index = 0),
            supy.steps.filters.value("%sMt%s%s"% (_muon[0], _muon[1], params["objects"]["met"]), min = 30.),
            supy.steps.filters.multiplicity("%sDiMuonNonIsoInZMass%s"%_muon, min=0, max=0),
            #supy.steps.filters.multiplicity("%sDiMuonOtherInZMass%s"%_muon, min=0, max=0),
            supy.steps.histos.histogrammer("%sMt%s%s"%(_muon[0], _muon[1], _met),
                                           50, 0, 200, title = ";M_{T} (GeV) of %s%s,%s;events / bin"%(_muon[0], _muon[1], _met)),
            ]

    def stepsPlotsOne(self, params):
        _jet = params["objects"]["jet"]
        _met = params["objects"]["met"]
        _et = "Et" if params["etRatherThanPt"] else "Pt"
        return [
            supy.steps.histos.histogrammer("vertexSumPt", 100, 0.0, 1.0e3, title=";SumPt of 2nd vertex (GeV);events / bin", funcString="lambda x:([0.0,0.0]+sorted(x))[-2]"),
            #steps.trigger.lowestUnPrescaledTriggerHistogrammer(),
            #supy.steps.filters.label("singleJetPlots1"),
            steps.jet.singleJetHistogrammer(_jet),
            supy.steps.filters.label("jetSumPlots1"),
            supy.steps.histos.histogrammer("%sDeltaPhiStar%s%s" % (_jet[0], params["lowPtName"], _jet[1]), 20, 0.0, r.TMath.Pi(), title=";#Delta#phi*;events / bin", funcString='lambda x:x[0][0]'),
            supy.steps.histos.histogrammer("%sDeltaPhiStar%s" % (_jet[0], _jet[1]), 20, 0.0, r.TMath.Pi(), title=";#Delta#phi*;events / bin", funcString='lambda x:x[0][0]'),
            supy.steps.filters.label("kinematicPlots1"),
            steps.jet.alphaHistogrammer(cs=_jet, deltaPhiStarExtraName=params["lowPtName"], etRatherThanPt=params["etRatherThanPt"]),
            ]

    def stepsQcdRejection(self, params):
        _jet = params["objects"]["jet"]
        _met = params["objects"]["met"]
        _et = "Et" if params["etRatherThanPt"] else "Pt"
	outList = []
        outList += [
		supy.steps.filters.value("%sMht%sOver%s" % (_jet[0], params["highPtName"] + _jet[1], _met+"Plus%sIndices%s" % params["objects"]["muon"]), max=1.25),
		steps.jet.cleanJetHtMhtHistogrammer(_jet, params["etRatherThanPt"]),
		steps.jet.alphaHistogrammer(cs=_jet, deltaPhiStarExtraName=params["lowPtName"], etRatherThanPt=params["etRatherThanPt"]),
		]
	if "PF" not in params["objects"]["jet"][0]:
		outList += [supy.steps.filters.value("%sMaxEmEnergyFraction%s"%(params["objects"]["jet"]), min = .1).onlyData(),
			    supy.steps.histos.histogrammer("%sMaxEmEnergyFraction%s" % (params["objects"]["jet"]), 20, 0.0, 1.0,
							   title=";MaxEmEnergyFraction;events / bin").onlyData(),]

	return outList

    def stepsMuonPlots(self, params):
        _muon = params["objects"]["muon"]
        ssteps = supy.steps
        return [
            ssteps.histos.pt("%sP4%s"%_muon, 60, 0, 600, indices = "%sIndices%s"%_muon, index = 0),
            ssteps.histos.eta("%sP4%s"%_muon, 50, -2.5, 2.5, indices = "%sIndices%s"%_muon, index = 0),
            ssteps.histos.phi("%sP4%s"%_muon, 60, 0, r.TMath.Pi(), indices = "%sIndices%s"%_muon, index = 0),]

    def stepsPlotsTwo(self, params):
        _jet = params["objects"]["jet"]
        _met = params["objects"]["met"]
        _et = "Et" if params["etRatherThanPt"] else "Pt"
        ssteps = supy.steps
        return [
            ssteps.histos.multiplicity("vertexIndices" , max = 60),
            steps.jet.cleanJetHtMhtHistogrammer(_jet, params["etRatherThanPt"]),
            ssteps.histos.generic("%sDeltaPhiStar%s%s" % (_jet[0], params["lowPtName"], _jet[1]), 20, 0.0, r.TMath.Pi(),
                                  title=";#Delta#phi*;events / bin", funcString='lambda x:x[0][0]'),
            ssteps.histos.generic("%sDeltaPhiStar%s" % (_jet[0], _jet[1]), 20, 0.0, r.TMath.Pi(), title=";#Delta#phi*;events / bin", funcString='lambda x:x[0][0]'),
            ssteps.histos.generic("%sIndices%s" % _jet, 15, -0.5, 14.5,
                                  title=";number of %s%s passing ID#semicolon p_{T}#semicolon #eta cuts;events / bin" % _jet,
                                  funcString="lambda x:len(x)"),
            ssteps.histos.generic("%sIndices%s" % _jet, 15, -0.5, 14.5,
                                  title=";number of %s%s passing ID#semicolon p_{T}#semicolon #eta cuts;events / bin" % _jet,
                                  funcString="lambda x:len(x)", suffix="ra1Category".join(_jet)),
            ssteps.histos.value("%sMht%sOver%s" % (_jet[0], params["highPtName"]+_jet[1], _met), 100, 0.0, 3.0,
                                xtitle="MHT %s%s / %s" % (_jet[0], params["highPtName"]+_jet[1], _met)),
            ssteps.histos.generic("%sIndicesBtagged2%s" % _jet, 7, -0.5, 6.5,
                                  title=";number of b-tagged %s%s passing ID#semicolon p_{T}#semicolon #eta cuts;events / bin" % _jet, funcString="lambda x:len(x)"),
            ssteps.histos.generic("%sIndicesBtagged2%s" % _jet, 7, -0.5, 6.5,
                                  title=";number of b-tagged %s%s passing ID#semicolon p_{T}#semicolon #eta cuts;events / bin" % _jet,
                                  funcString="lambda x:len(x)", suffix="ra1nJetCategory".join(_jet)),
            ssteps.histos.generic(("%sSumEt%s"%_jet, "%sIndicesBtagged2%s"%_jet), (8, 4), (375.0, -0.5), (1175.0, 3.5),
                                  title = ";H_{T} (GeV);n_{b};events / bin", funcString = "lambda x:(x[0],len(x[1]))"),
            ssteps.histos.generic(("%sSumEt%s"%_jet, "%sIndicesBtagged2%s"%_jet), (8, 4), (375.0, -0.5), (1175.0, 3.5),
                                  title = ";H_{T} (GeV);n_{b};events / bin", funcString = "lambda x:(x[0],len(x[1]))", suffix="ra1nJetCategory".join(_jet)),
            ssteps.histos.generic("%sSumEt%s"%_jet, 9, 150, 375.0, title = ";H_{T} (GeV) from %s%s E_{T}s;events / bin"%_jet, suffix="ra1nBJetCategory".join(_jet)),
            ssteps.histos.generic("%sSumEt%s"%_jet, 9, 150, 375.0, title = ";H_{T} (GeV) from %s%s E_{T}s;events / bin"%_jet,),
            ssteps.histos.generic("%sSumEt%s"%_jet, 8, 375.0, 1175.0, title = ";H_{T} (GeV) from %s%s E_{T}s;events / bin"%_jet),
            ssteps.histos.generic("%sSumEt%s"%_jet, 8, 375.0, 1175.0, title = ";H_{T} (GeV) from %s%s E_{T}s;events / bin"%_jet, suffix=("ra1AlphaTCategory"+_et).join(_jet)),
            ssteps.histos.generic("%sSumEt%s"%_jet, 8, 375.0, 1175.0, title = ";H_{T} (GeV) from %s%s E_{T}s;events / bin"%_jet, suffix="ra1nBJetCategory".join(_jet)),
            ssteps.histos.generic("%sSumEt%s"%_jet, 8, 375.0, 1175.0, title = ";H_{T} (GeV) from %s%s E_{T}s;events / bin"%_jet, suffix="ra1Category".join(_jet))]

    def stepsHtBins(self, params):
        _jet = params["objects"]["jet"]
        _et = "Et" if params["etRatherThanPt"] else "Pt"
        out = []
        out += [supy.steps.filters.value("%sSumEt%s" % _jet, min=bin) for bin in [475, 575, 675, 775, 875, 975, 1075]]
        return out

    def stepsOptional(self, params):
        return [
            #supy.steps.other.skimmer(),
            steps.other.duplicateEventCheck(),
            #steps.other.pickEventSpecMaker(),
            #steps.other.cutBitHistogrammer(self.togglePfJet(_jet), self.togglePfMet(_met)),
            #steps.Print.jetPrinter(_jet),
            #steps.Print.particleP4Printer(_muon),
            #steps.Print.particleP4Printer(_photon),
            #steps.Gen.genParticlePrinter(minPt = 10.0, minStatus = 3),
            ]

    def stepsDisplayer(self, params):
        jet = params["objects"]["jet"]
        return [
            steps.displayer.displayer(jets=jet,
                                      muons=params["objects"]["muon"],
                                      met = "%sPlus%sIndices%s" % (params["objects"]["met"],
                                                                   params["objects"]["muon"][0],
                                                                   params["objects"]["muon"][1]),
                                      electrons=params["objects"]["electron"],
                                      photons=params["objects"]["photon"],
                                      recHits=params["objects"]["rechit"], recHitPtThreshold=1.0,  # GeV
                                      scale=400.0,  # GeV
                                      etRatherThanPt=params["etRatherThanPt"],
                                      deltaPhiStarExtraName=params["lowPtName"],
                                      deltaPhiStarCut=0.5,
                                      deltaPhiStarDR=0.3,
                                      j2Factor=params["thresholds"][2]/params["thresholds"][0],
                                      mhtOverMetName="%sMht%sOver%sPlus%sIndices%s" % (jet[0], params["highPtName"]+jet[1],
                                                                            params["objects"]["met"],
                                                                            params["objects"]["muon"][0],
                                                                            params["objects"]["muon"][1]),
                                      metOtherAlgo=params["objects"]["metComp"],
                                      jetsOtherAlgo=params["objects"]["jetComp"],
                                      prettyMode=False
                                      #doGenJets=True,
                                      ),
            ]

    def listOfSteps(self, params):

        return ([supy.steps.printer.progressPrinter()] +
                #self.stepsEventCount(params, label="scanBefore") +
                self.stepsGenValidation() +
                self.stepsEvent(params) +
                self.stepsTrigger(params) +
                self.stepsHtLeadingJets(params) +
                self.stepsXclean(params) +
                #self.stepsOptional(params) +
                self.stepsQcdRejection(params) +
                self.stepsPlotsOne(params) +
                self.stepsPlotsTwo(params) +
                self.stepsBtagJets(params) +
                self.stepsMuonPlots(params) +
                self.stepsPlotsTwo(params) +
                #self.stepsDisplayer(params) +
                self.stepsHtBins(params) +
                #self.stepsEventCount(params, label="scanAfter") +
                [])

    def listOfSampleDictionaries(self):
        return [samples.muon17, samples.qcd17, samples.top17, samples.ewk17, samples.dyll17]

    def listOfSamples(self, params):
        from supy.samples import specify

        pu = calculables.gen.puWeight(var="pileupTrueNumInteractions", puEra="PU_Parked")
        #w = []
        btag = calculables.jet.BTagWeight(params["objects"]["jet"])

        def data_53X() :
            jw2012 = calculables.other.jsonWeight("cert/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt")
            out = []
            for era in ["A","B","C","D"]:
                out += specify(names="SingleMu.Run2012%s-22Jan2013" % era, weights=jw2012,)# nEventsMax=60000)

            return out

        def dyll(w=[]):
		out = []
		out += specify(names="dyll_HT_10To200_M-50", color=r.kBlue, weights=w  ,)
		out += specify(names="dyll_HT_200To400_M-50", color=r.kBlue, weights=w ,)
		out += specify(names="dyll_HT_400ToInf_M-50",  color=r.kBlue, weights=w,)
		return out

        def w_binned_LO_XS(w=[]):
            out = []
            xs = calculables.gen.xsWeight(file="wj_lv_mg_ht")
            out += specify(names="wj_lv_mg_ht_10To150_LO", color=r.kBlue, weights=w,)
            out += specify(names="wj_lv_mg_ht_150To200_LO", color=r.kOrange+1, weights=w+[xs],)# nEventsMax=50)
            out += specify(names="wj_lv_mg_ht_200To250_LO", color=r.kOrange+3, weights=w+[xs],)# nEventsMax=20)
            out += specify(names="wj_lv_mg_ht_250To300_LO", color=r.kOrange+5, weights=w+[xs], )# nEventsMax=50)
            out += specify(names="wj_lv_mg_ht_300To400_LO", color=r.kOrange+7, weights=w+[xs], )# nEventsMax=50)
            out += specify(names="wj_lv_mg_ht_400ToInf_LO", color=r.kOrange+9, weights=w+[xs], )# nEventsMax=60000)
            return out

        def w_inclusive(w=[]):
            out = []
            out += specify(names="wj_lv_mg_ht_incl_LO", color=r.kBlue, weights=w)
            return out

        def z_binned(w=[]):
            out = []
            out += specify("zinv_mg_ht_50_100", color=r.kRed, weights=w , )
            out += specify("zinv_mg_ht_100_200", color=r.kRed, weights=w, )
            out += specify("zinv_mg_ht_200_400", color=r.kRed, weights=w, )
            out += specify("zinv_mg_ht_400_inf", color=r.kRed, weights=w, )
            return out

        def vv(w=[]):
            out = []
            for diBos in ["ZZ", "WZ", "WW"]:
                out += specify("%s_py6" % diBos, color=r.kGreen, weights=w, )
            return out

        def top(w=[]):
            out = []
            out += specify(names="ttbar_CT10_powheg", color=r.kViolet, weights=w+["topPtWeight"], )
            for sTop in ["T_s", "T_t", "T_tW", "Tbar_s", "Tbar_t", "Tbar_tW"]:
                out += specify("%s_powheg" % sTop, color=r.kCyan, weights=w, )
            return out

        return (
            data_53X() +
            w_binned_LO_XS(w=[pu,btag]) +
            vv(w=[pu,btag]) +
            top(w=[pu,btag])+
	    dyll(w=[pu,btag])+
            []
            )

    def mergeSamples(self, org, withTrigEff=None, withSideBandWeight=None):
        def md(x, y):
            x.update(y)
            return x

        wjetSideBandCorr = 1.0
        ttSideBandCorr = 1.0

        if withSideBandWeight:
            wjetSideBandCorr = .868
            ttSideBandCorr = 1.110
            if "pf" in org.tag:
                wjetSideBandCorr = 1.005
                ttSideBandCorr = 1.122

        muonTrigEff = 1.0

        if withTrigEff :
            muonTrigEff = .9

        weightString = ".puWeight.bTagWeight"
        xsweightString = ".puWeight.bTagWeight.xsWeight"

        weightString2 = ".puWeight"
        xsweightString2 = ".puWeight.xsWeight"


        org.mergeSamples(targetSpec={"name": "2012 Data", "color": r.kBlack, "markerStyle": 20}, allWithPrefix="SingleMu")

        mcOps = {"markerStyle": 1, "lineWidth": 3, "goptions": "hist"}
        ##SingleTop
        org.mergeSamples(targetSpec=md({"name": "SingleTop", "color": r.kBlue+1}, mcOps),
                         sources=[x+weightString for x in ["T_s_powheg", "T_t_powheg",
                                                           "T_tW_powheg", "Tbar_t_powheg", "Tbar_tW_powheg",
                                                           "Tbar_s_powheg"]], scaleFactors=[muonTrigEff]*6)

        ##tt
        ttScaleFactor= 234./245.8
        org.mergeSamples(targetSpec=md({"name": "tt", "color": r.kRed+1}, mcOps),
                         sources=[x+weightString+".topPtWeight" for x in ["ttbar_CT10_powheg"]], scaleFactors=[muonTrigEff*ttSideBandCorr*ttScaleFactor])

        ##wjet
        kFactor = 37509./30400.

        ra1=[x*muonTrigEff for x in [257.73/235.6, 114.22/90.27, 56.39/48.01, 46.60/38.30, 28.97/25.22]]

        org.mergeSamples(targetSpec=md({"name": "W->lv + jets", "color": r.kOrange-3}, mcOps),
                         sources = [m + weightString for m in ["wj_lv_mg_ht_10To150_LO"]] +
                         [x+xsweightString for x in ["wj_lv_mg_ht_%s_LO" % y for y in ["150To200","200To250","250To300","300To400","400ToInf"]]]
                         , scaleFactors=[muonTrigEff*wjetSideBandCorr*kFactor] +[muonTrigEff*wjetSideBandCorr]*5, )


        org.mergeSamples(targetSpec=md({"name": "W->lv + jets", "color": r.kOrange-3}, mcOps),
                                                  sources = [m + weightString for m in ["wj_lv_mg_ht_incl_LO"]]
                                                  , scaleFactors=[muonTrigEff*wjetSideBandCorr])

        ##DY
        org.mergeSamples(targetSpec=md({"name": "Drell-Yan", "color": r.kMagenta-3}, mcOps), allWithPrefix="dyll", scaleFactors=[muonTrigEff]*3)
        ##VV
        org.mergeSamples(targetSpec=md({"name": "VV", "color": r.kOrange+3}, mcOps), sources=[x+weightString for x in ["WW_py6", "ZZ_py6", "WZ_py6"]],
                         scaleFactors=[muonTrigEff]*3)


        ##EWK
        ewkSources = ["tt", "SingleTop", "Z->vv + jets", "W->lv + jets", "VV", "Drell-Yan"]
        ewkSourcesNoTrigEff = [x + " no Trig Eff" for x in ["tt", "SingleTop", "Z->vv + jets", "W->lv + jets", "VV", "Drell-Yan"]]

        org.mergeSamples(targetSpec=md({"name": "Standard Model ", "color": r.kAzure+6}, mcOps), sources=ewkSources, keepSources=True,)
        #org.mergeSamples(targetSpec=md({"name": "Standard Model no btagWeights", "color": r.kAzure+3}, mcOps), sources=ewkSources2, keepSources=False,)


    def conclude(self, conf):
        org = self.organizer(conf)
        self.mergeSamples(org, withTrigEff=True, withSideBandWeight=False)
        org.scale()
        self.makeStandardPlots(org)
        #self.makeRootFiles(org)
        #self.makeAlphaTRootFiles(org)
        #self.makeNJetRootFiles(org)
        #for sample in org.samples:
        #    if sample["name"] in ["Standard Model "]:
        #        self.makeBTagEfficiencyPlots(org, sample["name"])
        #self.makeEfficiencyPlots(org)
        self.wJetHTSidebandCorrectionFactor(org)
        
    def makeStandardPlots(self, org):
        #plot
        pl = supy.plotter(org,
                          pdfFileName=self.pdfFileName(org.tag),
                          samplesForRatios=("2012 Data", "Standard Model "),
                          sampleLabelsForRatios=("data","s.m."),
                          #samplesForRatios=("2012 Data","tt"),
                          #sampleLabelsForRatios=("data","tt"),
                          printRatios=True,
                          showStatBox=True,
                          rowColors=[r.kBlack, r.kViolet+4],
                          #whiteList=["lowestUnPrescaledTrigger"],
                          doLog=True,
                          #pegMinimum=200.,
                          linYAfter=("variableGreaterFilter", "xcak5JetAlphaTEtPat>=0.550 "),
                          blackList=["lumiHisto", "xsHisto", "nJobsHisto"],
                          )
        pl.plotAll()

    def makeEfficiencyPlots(self, org):
        for dct in org.samples:
            self.makeEffPlots(org, sampleName=dct["name"])

    def makeIndividualPlots(self, org):
        ht = "H_{T}^{#color[0]{T}}"
        if "275" in org.tag:
            htLabel = "275 < %s < 325" % ht
        if "325" in org.tag:
            htLabel = "325 < %s < 375" % ht
        if "375" in org.tag:
            htLabel = "375 < %s" % ht

        #plot all
        pl = supy.plotter(org,
                          pdfFileName=self.pdfFileName(org.tag),
                          showStatBox=False,
                          doLog=False,
                          #pegMinimum=0.1,
                          anMode=False,
                          )
        pl.individualPlots(plotSpecs=[{"plotName": "xcak5JetAlphaTRoughPat",
                                       "stepName": "alphaHistogrammer",
                                       "stepDesc": "xcak5JetPat",
                                       "newTitle": ";#alpha_{T};events / bin",
                                       "legendCoords": (0.55, 0.60, 0.85, 0.90),
                                       "stampCoords": (0.75, 0.55)
                                       },
                                      {"plotName": "jetMultiplicity",
                                       "stepName": "singleJetHistogrammer",
                                       "stepDesc": "xcak5JetPat through index 2",
                                       "newTitle": ";N_{jets};events / bin",
                                       "legendCoords": (0.7, 0.7, 0.92, 0.92),
                                       "stampCoords": (0.5, 0.28),
                                       },
                                      {"plotName": "xcak5JetHtPat",
                                       "stepName": "cleanJetHtMhtHistogrammer",
                                       "stepDesc": "xcak5JetPat",
                                       "newTitle": ";H_{T} (GeV);events / bin",
                                       "legendCoords": (0.6, 0.60, 0.92, 0.92),
                                       "stampCoords": (0.45, 0.88)
                                       },
                                      ##after alphaT
                                      {"plotName": "xcak5JetDeltaPhiStarPat",
                                       "stepName": "histogrammer",
                                       "stepDesc": "(lambda x:x[0][0])(xcak5JetDeltaPhiStarPat)",
                                       "index": -1,
                                       "newTitle": ";#Delta#phi*;events / bin",
                                       "legendCoords": (0.6, 0.6, 0.92, 0.92),
                                       "stampCoords": (0.33, 0.88),
                                       },
                                      {"plotName": "xcak5JetHtPlusMhtRoughPat",
                                       "stepName": "cleanJetHtMhtHistogrammer",
                                       "stepDesc": "xcak5JetPat",
                                       "index": -1,
                                       "newTitle": ";M_{eff} (GeV);events / bin",
                                       "legendCoords": (0.7, 0.7, 0.92, 0.92),
                                       "stampCoords": (0.75, 0.4),
                                       },
                                      {"plotName": "xcak5JetIndicesPat",
                                       "stepName": "histogrammer",
                                       "stepDesc": "(lambda x:len(x))(xcak5JetIndicesPat)",
                                       "index": -1,
                                       "newTitle": ";N_{jets};events / bin",
                                       "legendCoords": (0.6, 0.6, 0.92, 0.92),
                                       "stampCoords": (0.6, 0.38),
                                       },
                                      {"plotName": "xcak5JetMbbListPat_2b",
                                       "stepName": "mbbHistogrammer",
                                       "stepDesc": "mbbHistogrammer",
                                       "index": -1,
                                       "stamp": False,
                                       "newTitle": "%s;m_{bb} (GeV) [2 b-jets];events / bin / 5.0 fb^{-1}" % htLabel,
                                       "legendCoords": (0.55, 0.6, 0.8, 0.75),
                                       "reBinFactor": 5,
                                       },
                                      ][-1:],
                           newSampleNames={"tt": "Madgraph t#bar{t}",
                                           "2012 Data": "Data"},
                           #newSampleNames={"qcd_mg_nVtx": "Madgraph QCD",
                           #                  "g_jets_mg_nVtx": "Madgraph #gamma + jets",
                           #                  "2011 Data": "Data",
                           #                  "standard_model_nVtx": "Standard Model",
                           #                  },
                           preliminary=True,
                           tdrStyle=False,
                           )

    def makeEffPlots(self, org, sampleName):

        def sampleIndex(org, name):
            for iSample, sample in enumerate(org.samples):
                if sample["name"] == name:
                    return iSample
            assert False, "could not find sample %s" % name

        def numerAndDenom(org, var):
            d = {}
            for selection in org.steps:
                if "scanBefore" in selection.title:
                    label = "before"
                elif "scanAfter" in selection.title:
                    label = "after"
                if selection.name != "scanHistogrammer":
                    continue
                d[label] = selection[var][sampleIndex(org, sampleName)].Clone(label)
            return d

        keep = []
        file = r.TFile("%s_%s.root" % (sampleName, org.tag), "RECREATE")
        canvas = r.TCanvas()
        canvas.SetRightMargin(0.2)
        canvas.SetTickx()
        canvas.SetTicky()
        psFileName = "%s_%s.ps" % (sampleName, org.tag)
        canvas.Print(psFileName+"[", "Lanscape")

        assert len(self.parameters()["objects"]) == 1
        for key, value in self.parameters()["objects"].iteritems():
            jet = value["jet"]

        for variable in ["nEvents"]:
            histos = numerAndDenom(org, variable)
            if "before" not in histos or "after" not in histos:
                continue
            result = histos["after"].Clone(variable)
            result.Divide(histos["before"])
            result.SetMarkerStyle(20)
            result.SetStats(False)
            if result.ClassName()[2] == "1":
                result.GetYaxis().SetRangeUser(0.0, 0.35)
                result.GetYaxis().SetTitle("efficiency")
                result.Draw()
            else:
                result.GetZaxis().SetRangeUser(0.0, 0.35)
                result.GetZaxis().SetTitle("efficiency")
                result.Draw("colz")
                canvas.Print(psFileName, "Lanscape")
                result.Write()
        canvas.Print(psFileName+"]", "Lanscape")
        temp = psFileName.replace(".chr", "_chr").replace(".ps", ".pdf")
        cmd = "ps2pdfwr " + psFileName + " " + temp
        os.system(cmd)
        os.remove(psFileName)
        file.Close()

    def makeBTagEfficiencyPlots(self, org, sample) :
        r.gStyle.SetNumberContours(40)
        def sampleIndex(org, name) :
            for iSample,sample in enumerate(org.samples) :
                if sample["name"]==name : return iSample
            assert False, "could not find sample %s"%name

        def numerAndDenom(org, var, sample) :
            d = {}
            for selection in org.steps :
                if selection.name != "bTagEfficiencyHistogrammer" : continue
                if var in selection :
                    d["denom"] = selection[var][sampleIndex(org, sample)].Clone("denom")
                    d["numer"] = selection[var+"_btagged"][sampleIndex(org,sample)].Clone("numer")
            return d

        keep = []
        canvas = r.TCanvas()
        canvas.SetRightMargin(0.2)
        canvas.SetTickx()
        canvas.SetTicky()
        psFileName = "data/bTagEff/muon_%s_%s_bTagEff.ps"% (org.tag, sample.replace(" ","_"))
        canvas.Print(psFileName+"[","Lanscape")
        f = r.TFile(psFileName.replace(".ps", ".root"), "RECREATE")
        for variable in [("genB","Bottom Jets",1.1),("genC","Charm Jets",0.8),
                         ("genL","Up, Down, Strange, Gluon Jets", 0.4)]:
            histos = numerAndDenom(org, variable[0], sample)
            if "numer" not in histos or "denom" not in histos : continue
            result = histos["numer"].Clone(variable[0])
            result.Reset()
            result.Divide(histos["numer"], histos["denom"], 1.0, 1.0, "b")
            result.SetMarkerStyle(20)
            result.SetStats(False)
            xtitle = histos["denom"].GetXaxis().GetTitle()
            ytitle = histos["denom"].GetYaxis().GetTitle()
            if result.ClassName()[2]=="1" :
                result.GetYaxis().SetRangeUser(0.0,variable[2])
                result.GetYaxis().SetTitle("b tagging efficiency")
                result.GetXaxis().SetTitle(xtitle)
                result.GetYaxis().SetTitle(ytitle)
                result.SetTitle(variable[1])
                result.Draw()
            else :
                result.GetZaxis().SetRangeUser(0.0,variable[2])
                result.GetXaxis().SetTitle(xtitle)
                result.GetYaxis().SetTitle(ytitle)
                result.GetZaxis().SetTitle("b tagging efficiency")
                result.SetTitle(variable[1])
                result.Draw("colz")
            canvas.Print(psFileName,"Lanscape")
            result.Write()
        canvas.Print(psFileName+"]","Lanscape")
        os.system("ps2pdf "+psFileName)
        os.remove(psFileName)
        f.Close()

    def makeRootFiles(self, org) :

        def sampleIndex(org, name) :
            for iSample,sample in enumerate(org.samples) :
                if sample["name"]==name : return iSample
            assert False, "could not find sample %s"%name

        def histo(name = "", suffix = "", samples = ["2012 Data", "Standard Model "]):
            lst = []
            for selection in org.steps :
                if selection.name != "generic" : continue
                if selection.title!="(lambda x:(x[0],len(x[1])))(%s)"%(name+suffix) : continue
                dct = {}
                for s in samples :
                    dct[s] = {}
                    for nJet in ["le3j","ge4j"]:
                        dct[s][nJet] = selection[name+"_"+nJet][sampleIndex(org, s)]
                lst.append(dct)
            return lst[-1]

        ewkSources = ["Standard Model ", "tt", "SingleTop", "W->lv + jets", "VV", "Drell-Yan"]
        allSources = ["2012 Data"] + ewkSources
        histNames = ["Data", "MCYield", "TTbar", "SingleTop", "WJets", "DiBoson", "DY"]
        if "calo" in org.tag:
            dct = histo(samples=allSources, name = "xcak5JetIndicesBtagged2Pat_vs_xcak5JetSumEtPat",
                        suffix="xcak5Jetra1nJetCategoryPat" )
        else:
            dct = histo(samples=allSources, name="xcak5JetPFIndicesBtagged2Pat_vs_xcak5JetPFSumEtPat",
                        suffix="xcak5JetPFra1nJetCategoryPat" )

        for nJet in ["le3j","ge4j"]:
                if "375" not in org.tag: continue
	        for iBTag in range(4) :
	            f = r.TFile("yields/%s_%s_%db.root"%(org.tag, nJet, iBTag), "RECREATE")
	            f.mkdir("muon")
	            f.cd("muon")
	            for s in ["lumiData", "lumiMc"] :
	                lumi = r.TH1D(s, s, 1, -0.5, 0.5)
	                lumi.SetBinContent(1, org.lumi*1.0e-3)#/fb
	                lumi.Write()
	            for name,key in zip(histNames, allSources):
	                hIn = dct[key][nJet]

	                xMin   = hIn.GetXaxis().GetXmin()
	                xMax   = hIn.GetXaxis().GetXmax()
	                nBinsX = hIn.GetXaxis().GetNbins()

	                assert abs(xMin-375.)<1.0e-6,xMin
	                assert abs(xMax-1175.)<1.0e-6,xMax
	                assert nBinsX==8,nBinsX

	                yMin   = hIn.GetYaxis().GetXmin()
	                yMax   = hIn.GetYaxis().GetXmax()
	                nBinsY = hIn.GetYaxis().GetNbins()

	                assert abs(yMin+0.5)<1.0e-6,yMin
	                assert abs(yMax-3.5)<1.0e-6,yMax
	                assert nBinsY==4,nBinsY

	                h1 = hIn.ProjectionX("%s_projX"%name, 1+iBTag, 1+iBTag)

	                xBinsLo = array.array('d',[275., 325.]+[375.+100*i for i in range(9)])
	                yBinsLo = array.array('d',[55.0, 60.0])
	                hOut = r.TH2D(name, name, len(xBinsLo)-1, xBinsLo, len(yBinsLo)-1, yBinsLo)

	                for iBinX in range(1,1+h1.GetNbinsX()) :
	                    hOut.SetBinContent(2+iBinX, 1, h1.GetBinContent(iBinX))
	                    hOut.SetBinError(2+iBinX, 1, h1.GetBinError(iBinX))

	                hOut.Write()
	            f.Close()

    def makeAlphaTRootFiles(self, org) :

        def sampleIndex(org, name) :
            for iSample,sample in enumerate(org.samples) :
                if sample["name"]==name : return iSample
            assert False, "could not find sample %s"%name

        def histo(name = "", suffix = "", samples = ["2012 Data", "Standard Model "]):
            lst = []
            for selection in org.steps :
                if selection.name != "generic" : continue
                if selection.title!="(%s)"%(name+suffix) : continue
                dct = {}
                for s in samples :
                    dct[s] = {}
                    for nJet in ["le3j","ge4j"]:
                        for aT in ["aTge55","aTl55"]:
                            dct[s][nJet+"_"+aT] = selection[name+"_"+nJet+"_"+aT][sampleIndex(org, s)]
                lst.append(dct)
            return lst[-1]

        ewkSources = ["Standard Model ", "tt", "SingleTop", "W->lv + jets", "VV", "Drell-Yan"]
        allSources = ["2012 Data"] + ewkSources
        histNames = ["Data", "MCYield", "TTbar", "SingleTop", "WJets", "DiBoson", "DY"]
        if "calo" in org.tag:
            dct = histo(samples=allSources, name = "xcak5JetSumEtPat",
                        suffix="xcak5Jetra1AlphaTCategoryEtPat" )
        else:
            dct = histo(samples=allSources, name="xcak5JetPFSumEtPat",
                        suffix="xcak5JetPFra1AlphaTCategoryEtPat" )

        for nJet in ["le3j","ge4j"]:
                if "375" not in org.tag: continue
                f = r.TFile("yields/%s_aT_%s.root"%(org.tag, nJet), "RECREATE")
                f.mkdir("muon")
                f.cd("muon")
                for s in ["lumiData", "lumiMc"] :
                    lumi = r.TH1D(s, s, 1, -0.5, 0.5)
                    lumi.SetBinContent(1, org.lumi*1.0e-3)#/fb
                    lumi.Write()
                for name,key in zip(histNames, allSources):
                    hl55 = dct[key][nJet+"_aTl55"]
                    hge55 = dct[key][nJet+"_aTge55"]

                    xMin   = hl55.GetXaxis().GetXmin()
                    xMax   = hl55.GetXaxis().GetXmax()
                    nBinsX = hl55.GetXaxis().GetNbins()

                    assert abs(xMin-375.)<1.0e-6,xMin
                    assert abs(xMax-1175.)<1.0e-6,xMax
                    assert nBinsX==8,nBinsX

                    xBinsLo = array.array('d',[275., 325.]+[375.+100*i for i in range(9)])
                    yBinsLo = array.array('d',[0,55.0,100.0])
                    hOut = r.TH2D(name, name, len(xBinsLo)-1, xBinsLo, len(yBinsLo)-1, yBinsLo)

                    for iBinX in range(1,1+hl55.GetNbinsX()) :
                        hOut.SetBinContent(2+iBinX, 1, hl55.GetBinContent(iBinX))
                        hOut.SetBinError(2+iBinX, 1, hl55.GetBinError(iBinX))
                        hOut.SetBinContent(2+iBinX, 2, hge55.GetBinContent(iBinX))
                        hOut.SetBinError(2+iBinX, 2, hge55.GetBinError(iBinX))

                    hOut.Write()
                f.Close()

    def makeNJetRootFiles(self, org) :

        def sampleIndex(org, name) :
            for iSample,sample in enumerate(org.samples) :
                if sample["name"]==name : return iSample
            assert False, "could not find sample %s"%name

        def histo(name = "", suffix = "", samples = ["2012 Data", "Standard Model "]):
            lst = []
            for selection in org.steps :
                if selection.name != "generic" : continue
                if selection.title!="(%s)"%(name+suffix) : continue
                dct = {}
                for s in samples :
                    dct[s] = {}
                    for nJet in ["le3j","ge4j"]:
                        for bjet in ["eq0b","eq1b","eq2b"]:
                            dct[s][nJet+"_"+bjet] = selection[name+"_"+nJet+"_"+bjet][sampleIndex(org, s)]
                lst.append(dct)
            return lst[-1]

        ewkSources = ["Standard Model ", "tt", "SingleTop", "W->lv + jets", "VV", "Drell-Yan"]
        allSources = ["2012 Data"] + ewkSources
        histNames = ["Data", "MCYield", "TTbar", "SingleTop", "WJets", "DiBoson", "DY"]
        if "calo" in org.tag:
            dct = histo(samples=allSources, name = "xcak5JetSumEtPat",
                        suffix="xcak5Jetra1CategoryPat" )
        else:
            dct = histo(samples=allSources, name="xcak5JetPFSumEtPat",
                        suffix="xcak5JetPFra1CategoryPat" )

        for nJet in ["le3j","ge4j"]:
                if "375" not in org.tag: continue
                f = r.TFile("yields/%s_nJet_%s.root"%(org.tag, nJet), "RECREATE")
                f.mkdir("muon")
                f.cd("muon")
                for s in ["lumiData", "lumiMc"] :
                    lumi = r.TH1D(s, s, 1, -0.5, 0.5)
                    lumi.SetBinContent(1, org.lumi*1.0e-3)#/fb
                    lumi.Write()
                for name,key in zip(histNames, allSources):
                    h0 = dct[key][nJet+"_eq0b"]
                    h1 = dct[key][nJet+"_eq1b"]
                    h2 = dct[key][nJet+"_eq2b"]
                    hTmp = h0.Clone()
                    hTmp.Add(hTmp,h1)
                    hTmp.Add(hTmp,h2)
                    h = hTmp.Clone()
                    xMin   = h.GetXaxis().GetXmin()
                    xMax   = h.GetXaxis().GetXmax()
                    nBinsX = h.GetXaxis().GetNbins()

                    assert abs(xMin-375.)<1.0e-6,xMin
                    assert abs(xMax-1175.)<1.0e-6,xMax
                    assert nBinsX==8,nBinsX

                    xBinsLo = array.array('d',[275., 325.]+[375.+100*i for i in range(9)])
                    yBinsLo = array.array('d',[0,55.0])
                    hOut = r.TH2D(name, name, len(xBinsLo)-1, xBinsLo, len(yBinsLo)-1, yBinsLo)

                    for iBinX in range(1,1+h.GetNbinsX()) :
                        hOut.SetBinContent(2+iBinX, 1, h.GetBinContent(iBinX))
                        hOut.SetBinError(2+iBinX, 1, h.GetBinError(iBinX))
                    hOut.Write()
                f.Close()

    def wJetHTSidebandCorrectionFactor(self, org) :

        def sampleIndex(org, name) :
            for iSample,sample in enumerate(org.samples) :
                if sample["name"]==name : return iSample
            assert False, "could not find sample %s"%name

        def histo(org = None, btag="",samples = []):
            htName = "xcak5JetPFSumEtPat"
            btagName = "xcak5JetPFra1nBJetCategoryPat"
            if "calo" in org.tag :
                htName = "xcak5JetSumEtPat"
                btagName = "xcak5Jetra1nBJetCategoryPat"
            lst = []
            for selection in org.steps :
                if selection.name != "generic" : continue
                if selection.title!="(%s)"%(htName+btagName) : continue
                dct = {}
                for s in samples :
                    dct[s]= selection[htName+"_"+btag][sampleIndex(org, s)]
                lst.append(dct)
                break
            return lst[-1]

        if "200" not in org.tag : return
        f = open("data/tmp/htSideBandFactors_%s.txt" % org.tag,'w')
        for b in [("eq0b","W->lv + jets"), ("eq2b","tt")]:
           dct = histo(org, btag=b[0], samples =["2012 Data", "Standard Model ", b[1]])
           mc = dct[b[1]].Clone()
           sm = dct["Standard Model "].Clone("sm")
           data = dct["2012 Data"].Clone("data")
           cFactor = dct["2012 Data"].Clone("factor")
           purity = mc.Clone("purity")
           purity.Reset()
           purity.Divide(mc,sm,1,1,"b")
           cFactor.Reset()
           cFactor.Divide(data,sm,1,1,"b")
           f.write("%s \n" %org.tag)
           f.write("%s: \n\nHT bin, Purity, Data, MC, MC error, Correction, Correction error \n" % mc.GetName())
           for ibin in range(data.GetXaxis().GetNbins()):
               ht = data.GetXaxis().GetBinLowEdge(ibin)
               if (ht == 200.0 or ht == 300.0):
                   vals = "%s, %1.3f, %1.3f, %1.3f, %1.3f, %1.3f, %1.3f \n" % ("%i-%i" %(ht,ht+25),
                                                                               purity.GetBinContent(ibin),
                                                                               data.GetBinContent(ibin),
                                                                               mc.GetBinContent(ibin),
                                                                               mc.GetBinError(ibin),
                                                                               cFactor.GetBinContent(ibin),
                                                                               cFactor.GetBinError(ibin))
                   f.write(vals)


        f.close()
