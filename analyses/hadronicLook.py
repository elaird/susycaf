import os

import ROOT as r
import calculables
import samples
import steps
import supy


def triggerTuple(l=[], keys=[]):
    out = []
    for item in l:
        stem = "HLT"
        for key in keys:
            stem += "_%s%s" % (key, str(item[key]).replace(".", "p"))
        for version in item["v"]:
            out.append("%s_v%d" % (stem, version))
    return tuple(out)

triggers_2012 = [{"HT": 200, "AlphaT": 0.57, "v":    [5, 6, 8]},
                 {"HT": 250, "AlphaT": 0.55, "v": range(1,  9)},
                 {"HT": 300, "AlphaT": 0.53, "v": range(1,  4)},
                 {"HT": 350, "AlphaT": 0.52, "v": range(1,  9)},
                 {"HT": 400, "AlphaT": 0.51, "v": range(13,20)},
                 ]


class hadronicLook(supy.analysis):
    def parameters(self):
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
                                          }),
                "thresholds": self.vary(dict([("200", (200.0, 325.0,   73.3, 36.7)),
                                             ("375",  (375.0, None,  100.0,  50.0)),
                                         ])),
                "etRatherThanPt": True,
                "lowPtThreshold": 30.0,
                "lowPtName": "lowPt",
                "highPtThreshold": 50.0,
                "highPtName": "highPt",
                "signalScan": False,
                #required to be sorted
                "triggerList": triggerTuple(l=triggers_2012,
                                            keys=("HT", "AlphaT"),
                                            ),
                }

    def calcListJet(self, obj, etRatherThanPt, ptMin,
                    lowPtThreshold, lowPtName,
                    highPtThreshold, highPtName, htThreshold):

        def calcList(jet=None, jetId=None, met=None,
                     photon=None, muon=None, electron=None,
                     **_):
            return [calculables.xclean.xcJet(jet,
                                             gamma=photon,
                                             gammaDR=0.5,
                                             muon=muon,
                                             muonDR=0.5,
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
                    calculables.jet.MhtOverMet((jet[0], highPtName + jet[1]),
                                               met),
                    calculables.jet.BTagProbability(collection=jet, selection="had"),
                    calculables.jet.DeadEcalDR(jet,
                                               extraName=lowPtName,
                                               minNXtals=10),
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
                calculables.muon.Indices(obj["muon"], ptMin=10, isoMax=0.12,
                                         ISO="PfIsolationR04DeltaBCorrected",
                                         ID="IdPog2012Tight",
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
                calculables.vertex.ID(),
                calculables.vertex.Indices(),
                calculables.trigger.lowestUnPrescaledTrigger(triggers),
                calculables.other.hcalLaserEvent2012(),
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
            #steps.trigger.hltPrescaleHistogrammer(params["triggerList"]).onlyData(),
            steps.trigger.l1Filter("L1Tech_BPTX_plus_AND_minus.v0").onlyData(),
            steps.trigger.physicsDeclaredFilter().onlyData(),
            ]+([] if params["thresholds"][1] != 275.0 else [steps.trigger.lowestUnPrescaledTriggerFilter().onlyData()])  # apply trigger in lowest HT bin

    def stepsGenValidation(self):
        return [
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
        _met = params["objects"]["met"]
        return [
            #steps.other.iterHistogrammer("ecalDeadTowerTrigPrimP4", 256, 0.0, 128.0,
                                         #title=";E_{T} of ECAL TP in each dead region (GeV);TPs / bin", funcString="lambda x:x.Et()"),
            supy.steps.filters.multiplicity("%sIndices%s" % params["objects"]["jet"], min=params["nJetsMinMax"][0], max=params["nJetsMinMax"][1]),
            supy.steps.filters.multiplicity("%sIndicesOther%s" % params["objects"]["jet"], max=0),
            supy.steps.filters.multiplicity("%sIndices%s" % params["objects"]["electron"], max=0),
            supy.steps.filters.multiplicity("%sIndices%s" % params["objects"]["photon"], max=0),
            supy.steps.filters.multiplicity("%sIndices%s" % params["objects"]["muon"], max=0),
            supy.steps.filters.multiplicity("%sIndicesUnmatched%s" % params["objects"]["electron"], max = 0),
            supy.steps.filters.multiplicity("%sIndicesUnmatched%s" % params["objects"]["photon"], max = 0),
            steps.other.deadEcalFilter(jets=params["objects"]["jet"], extraName=params["lowPtName"], dR=0.3, dPhiStarCut=0.5),
            ]

    def stepsPlotsOne(self, params):
        _jet = params["objects"]["jet"]
        _met = params["objects"]["met"]
        _et = "Et" if params["etRatherThanPt"] else "Pt"
        return [
            supy.steps.histos.histogrammer("vertexIndices", 20, -0.5, 19.5,
                                           title=";N vertices;events / bin", funcString="lambda x:len(x)"),
            supy.steps.histos.histogrammer("vertexSumPt", 100, 0.0, 1.0e3,
                                           title=";SumPt of 2nd vertex (GeV);events / bin", funcString="lambda x:([0.0,0.0]+sorted(x))[-2]"),
            steps.jet.singleJetHistogrammer(_jet),
            supy.steps.filters.label("jetSumPlots1"),
            steps.jet.cleanJetHtMhtHistogrammer(_jet, params["etRatherThanPt"]),
            supy.steps.histos.histogrammer("%sDeltaPhiStar%s%s" % (_jet[0], params["lowPtName"], _jet[1]), 20, 0.0, r.TMath.Pi(),
                                           title=";#Delta#phi*;events / bin", funcString='lambda x:x[0][0]'),
            supy.steps.histos.histogrammer("%sDeltaPhiStar%s" % (_jet[0], _jet[1]), 20, 0.0, r.TMath.Pi(),
                                           title=";#Delta#phi*;events / bin", funcString='lambda x:x[0][0]'),
            supy.steps.histos.histogrammer(_met, 100, 0.0, 500.0, title=";"+_met+" (GeV);events / bin", funcString="lambda x: x.pt()"),
            steps.jet.alphaHistogrammer(cs=_jet, deltaPhiStarExtraName=params["lowPtName"], etRatherThanPt=params["etRatherThanPt"]),
            ]

    def stepsQcdRejection(self, params):
        _jet = params["objects"]["jet"]
        _met = params["objects"]["met"]
        _et = "Et" if params["etRatherThanPt"] else "Pt"

        outList = [supy.steps.histos.histogrammer("%sMht%sOver%s" % (_jet[0], params["highPtName"]+_jet[1], _met), 100, 0.0, 3.0,
                                                  title=";MHT %s%s / %s;events / bin" % (_jet[0], params["highPtName"]+_jet[1], _met)),
                   supy.steps.filters.value("%sMht%sOver%s" % (_jet[0], params["highPtName"]+_jet[1], _met), max=1.25),
                   steps.jet.cleanJetHtMhtHistogrammer(_jet, params["etRatherThanPt"]),
                   steps.jet.alphaHistogrammer(cs=_jet, deltaPhiStarExtraName=params["lowPtName"], etRatherThanPt=params["etRatherThanPt"]),
                   supy.steps.histos.histogrammer("%sRecHitSumPt" % params["objects"]["rechit"], 30, 0, 300,
                                                  title=";Sum of HBHE (sev.#geq10), EB,EE (sev.#geq2) RecHit p_{T} (GeV);events / bin"),
                   supy.steps.filters.value("%sRecHitSumPt" % params["objects"]["rechit"], max=30.0),
                   supy.steps.filters.value("%sAlphaT%s%s" % (_jet[0], _et, _jet[1]), min=0.55),]
        if "PF" not in params["objects"]["jet"][0]:
            outList += [supy.steps.filters.value("%sMaxEmEnergyFraction%s"%(_jet[0],_jet[1]), min = .1),
                        supy.steps.histos.histogrammer("%sMaxEmEnergyFraction%s" % (_jet[0], _jet[1]), 20, 0.0, 1.0,
                                                       title=";MaxEmEnergyFraction;events / bin"),]
        return outList
        

    def stepsPlotsTwo(self, params):
        _jet = params["objects"]["jet"]
        _met = params["objects"]["met"]
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
            ssteps.histos.generic("%sSumEt%s"%_jet, 8, 375.0, 1175.0, title = ";H_{T} (GeV) from %s%s E_{T}s;events / bin"%_jet, suffix="ra1nBJetCategory".join(_jet)),
            ssteps.histos.generic("%sSumEt%s"%_jet, 8, 375.0, 1175.0, title = ";H_{T} (GeV) from %s%s E_{T}s;events / bin"%_jet, suffix="ra1Category".join(_jet))]

    def stepsHtBins(self, params):
        _jet = params["objects"]["jet"]
        return [supy.steps.filters.value("%sSumEt%s" % _jet, min=bin) for bin in [475, 575, 675, 775, 875, 975, 1075]]

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
                                      met=params["objects"]["met"],
                                      electrons=params["objects"]["electron"],
                                      photons=params["objects"]["photon"],
                                      recHits=params["objects"]["rechit"], recHitPtThreshold=1.0,  # GeV
                                      scale=400.0,  # GeV
                                      etRatherThanPt=params["etRatherThanPt"],
                                      deltaPhiStarExtraName=params["lowPtName"],
                                      deltaPhiStarCut=0.5,
                                      deltaPhiStarDR=0.3,
                                      j2Factor=params["thresholds"][2]/params["thresholds"][0],
                                      mhtOverMetName="%sMht%sOver%s" % (jet[0], params["highPtName"]+jet[1], params["objects"]["met"]),
                                      metOtherAlgo=params["objects"]["metComp"],
                                      jetsOtherAlgo=params["objects"]["jetComp"],
                                      #doGenJets=True,
                                      ),
            ]

    def stepsMbb(self, params):
        _jet = params["objects"]["jet"]
        return [
            #supy.steps.filters.multiplicity("genIndicesStatus3b", min=4, max=4),
            #steps.printer.eventPrinter(),
            #steps.printer.jetPrinter(_jet),
            #steps.gen.particlePrinter(),
            supy.steps.filters.multiplicity("%sIndicesBtagged2%s" % _jet, min=2, max=2),
            supy.steps.histos.multiplicity("%sIndicesBtagged2%s" % _jet),
            supy.steps.histos.eta("%sCorrectedP4%s" % _jet, 24, -3.0, 3.0, indices="%sIndicesBtagged2%s" % _jet, index=0, xtitle="b jet 0"),
            supy.steps.histos.eta("%sCorrectedP4%s" % _jet, 24, -3.0, 3.0, indices="%sIndicesBtagged2%s" % _jet, index=1, xtitle="b jet 1"),
            supy.steps.histos.mass("%sSumP4Btagged2%s" % _jet, 24, 0.0, 1200.0, xtitle="sum P4 {b jets}"),
            supy.steps.histos.pt("%sSumP4Btagged2%s" % _jet, 24, 0.0, 1200.0, xtitle="sum P4 {b jets}"),
            supy.steps.filters.mass("%sSumP4Btagged2%s" % _jet, min=450.0),
            supy.steps.histos.pt("%sSumP4Btagged2%s" % _jet, 24, 0.0, 1200.0, xtitle="sum P4 {b jets}"),
            #steps.jet.mbbHistogrammer(_jet, drMatch=0.2, bZDaughters="genIndicesStatus3bZDaughters"),
            ]

    def listOfSteps(self, params):
        _jet = params["objects"]["jet"]
        _met = params["objects"]["met"]
        _etRatherThanPt = params["etRatherThanPt"]
        _et = "Et" if _etRatherThanPt else "Pt"

        return ([supy.steps.printer.progressPrinter()] +
                #self.stepsEventCount(params, label="scanBefore") +
                self.stepsGenValidation() +
                self.stepsEvent(params) +
                self.stepsTrigger(params) +
                self.stepsHtLeadingJets(params) +
                self.stepsXclean(params) +
                #self.stepsOptional(params) +
                self.stepsQcdRejection(params) +
                self.stepsBtagJets(params) +
                self.stepsSlowLaserFilters() +
                self.stepsPlotsOne(params) +
                self.stepsPlotsTwo(params) +
                #self.stepsMbb(params) +
                #self.stepsDisplayer(params) +
                self.stepsHtBins(params) +
                #self.stepsEventCount(params, label="scanAfter") +
                [])

    def listOfSampleDictionaries(self):
        return [samples.ht17, samples.top17, samples.ewk17, samples.qcd17, samples.susy17, samples.photon17, samples.dyll17]

    def listOfSamples(self, params):
        from supy.samples import specify

        pu = calculables.gen.puWeight(var="pileupTrueNumInteractions", puEra="PU_Parked")
        #w = []
        btag = calculables.jet.BTagWeight(params["objects"]["jet"])
        w = [pu, btag]

        def data_53X():
            jw2012 = calculables.other.jsonWeight("cert/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt")
            out = []
            #out += specify(names="HT.Run2012A-22Jan2013", weights=jw2012,)
            for era in ["B","C","D"][1:2]:
                out += specify(names="HTMHTParked.Run2012%s-22Jan2013"%era, weights=jw2012, )
            return out

        def dyll():
		out = []
		out += specify(names="dyll_HT_10To200_M-50", color=r.kBlue, weights=w, )
		out += specify(names="dyll_HT_200To400_M-50", color=r.kBlue, weights=w, )
		out += specify(names="dyll_HT_400ToInf_M-50",  color=r.kBlue, weights=w, )
		return out

        def w_binned_LO_XS():
            out = []
            xs = calculables.gen.xsWeight(file="wj_lv_mg_ht")
            out += specify(names="wj_lv_mg_ht_10To150_LO", color=r.kBlue, weights=w, )# )
            out += specify(names="wj_lv_mg_ht_150To200_LO", color=r.kOrange+1, weights=w+[xs],)# nEventsMax=50)
            out += specify(names="wj_lv_mg_ht_200To250_LO", color=r.kOrange+3, weights=w+[xs], )# nEventsMax=50)
            out += specify(names="wj_lv_mg_ht_250To300_LO", color=r.kOrange+5, weights=w+[xs], )# nEventsMax=50)
            out += specify(names="wj_lv_mg_ht_300To400_LO", color=r.kOrange+7, weights=w+[xs], )# nEventsMax=50)
            out += specify(names="wj_lv_mg_ht_400ToInf_LO", color=r.kOrange+9, weights=w+[xs], )#, nEventsMax=60000)
            return out

        def z_binned():
            out = []
            out += specify("zinv_mg_ht_50_100", color=r.kBlue, weights=w,)
            out += specify("zinv_mg_ht_100_200", color=r.kGreen, weights=w,)
            out += specify("zinv_mg_ht_200_400", color=r.kOrange, weights=w,)
            out += specify("zinv_mg_ht_400_inf", color=r.kViolet, weights=w,)
            return out

        def vv():
            out = []
            for diBos in ["ZZ", "WZ", "WW"]:
                out += specify("%s_py6" % diBos, weights=w,)
            return out

        def top():
            out = []
            out += specify(names="ttbar_CT10_powheg", weights=w+["topPtWeight"],)
            for sTop in ["T_s", "T_t", "T_tW", "Tbar_s", "Tbar_t", "Tbar_tW"]:
                out += specify("%s_powheg" % sTop, weights=w,)
            return out

        def susy(eL):
            return specify(names="lm6", effectiveLumi=eL, color=r.kRed)

        def sms():
            out = []
            #out += specify(names="t2bb.job418")#, nFilesMax=1, nEventsMax=500)
            out += specify(names="t1tttt.job442")  # , nFilesMax=1, nEventsMax=500)
            #out += specify(names="t1bbbb.job443", nFilesMax=1, nEventsMax=500)
            #out += specify(names="t1.job444", nFilesMax=1, nEventsMax=500)
            #out += specify(names="t2tt.job445")#, nFilesMax=1, nEventsMax=500)
            #out += specify(names="t2.job446", nFilesMax=1, nEventsMax=500)
            return out

        return (
            data_53X() +
            #dyll() +
            #w_binned_LO_XS() +
            #z_binned() +
            #top() +
            #vv() +
            #sms() +
            []
            )

    def mergeSamples(self, org):
        def md(x, y):
            x.update(y)
            return x

        org.mergeSamples(targetSpec={"name": "2012 Data", "color": r.kBlack, "markerStyle": 20}, allWithPrefix="HT")

        mcOps = {"markerStyle": 1, "lineWidth": 3, "goptions": "hist"}

#        wjetSideBandCorr = .808
#        ttSideBandCorr = 1.038
        wjetSideBandCorr = 1.0
        ttSideBandCorr = 1.0
        if "pf" in org.tag:
            wjetSideBandCorr = 1.0
            ttSideBandCorr = 1.0
#            wjetSideBandCorr = .907
#            ttSideBandCorr = 1.041
        weightString = ".puWeight"
        xsweightString = ".puWeight.xsWeight"


        ##SingleTop
        org.mergeSamples(targetSpec=md({"name": "SingleTop", "color": r.kBlue+1}, mcOps),
                         sources=[x+weightString for x in ["T_s_powheg", "T_t_powheg",
                                                           "T_tW_powheg", "Tbar_t_powheg", "Tbar_tW_powheg",
                                                           "Tbar_s_powheg"]],)
        ##tt
        org.mergeSamples(targetSpec=md({"name": "tt", "color": r.kRed+1}, mcOps),
                         sources=[x+weightString+".topPtWeight" for x in ["ttbar_CT10_powheg"]],)
        ##wjet
        kFactor = 37509./30400.
        org.mergeSamples(targetSpec=md({"name": "W->lv + jets", "color": r.kOrange-3}, mcOps),
                         sources = [m + weightString for m in ["wj_lv_mg_ht_10To150_LO"]] +
                         [x+xsweightString for x in ["wj_lv_mg_ht_%s_LO" % y for y in ["150To200","200To250","250To300","300To400","400ToInf"]]]
                         , scaleFactors=[kFactor]+[1]*6, )
        ##zinv
        org.mergeSamples(targetSpec=md({"name": "Z->vv + jets", "color": r.kMagenta-3}, mcOps), allWithPrefix="zinv")
        ##DY
        org.mergeSamples(targetSpec=md({"name": "Drell-Yan", "color": r.kMagenta-3}, mcOps), allWithPrefix="dyll")
        ##VV
        org.mergeSamples(targetSpec=md({"name": "VV", "color": r.kOrange+3}, mcOps),
                         sources=[x+weightString for x in ["WW_py6", "ZZ_py6", "WZ_py6"]])
        ##EWK
        ewkSources = ["tt", "SingleTop", "Z->vv + jets", "W->lv + jets", "VV", "Drell-Yan"]

        org.mergeSamples(targetSpec=md({"name": "Standard Model ", "color": r.kAzure+6}, mcOps), sources=ewkSources, keepSources=True,)


    def conclude(self, conf):
        org = self.organizer(conf)
        self.mergeSamples(org)
        ##for skimming only
        #utils.printSkimResults(org)
        org.scale(20000) if not self.parameters()["signalScan"] else org.scale(100.0)
        self.makeStandardPlots(org)
        #self.makeIndividualPlots(org)
        #self.makeEfficiencyPlots(org)
        for sample in org.samples:
                if sample["name"] in ["Standard Model "]:
                    self.makeBTagEfficiencyPlots(org, sample["name"])

    def makeStandardPlots(self, org):
        #plot
        pl = supy.plotter(org,
                          pdfFileName=self.pdfFileName(org.tag),
                          samplesForRatios=("2012 Data", "Standard Model "),
                          sampleLabelsForRatios=("data", "s.m."),
                          #samplesForRatios=("2012 Data","tt"),
                          #sampleLabelsForRatios=("data","tt"),
                          printRatios=True,
                          showStatBox=True,
                          rowColors=[r.kBlack, r.kViolet+4],
                          #whiteList=["lowestUnPrescaledTrigger"],
                          doLog=True,
                          #pegMinimum=0.1,
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
        psFileName = "had_%s_%s_bTagEff.ps"% (org.tag, sample.replace(" ","_"))
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
