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
                 {"HT": 250, "AlphaT": 0.55, "v": range(1, 10)},
                 {"HT": 300, "AlphaT": 0.53, "v": range(1,  4)},
                 {"HT": 300, "AlphaT": 0.54, "v": range(1, 10)},
                 {"HT": 350, "AlphaT": 0.52, "v": range(1,  4)},
                 {"HT": 350, "AlphaT": 0.53, "v": range(1, 15)},
                 {"HT": 400, "AlphaT": 0.51, "v": range(1, 15)},
                 {"HT": 400, "AlphaT": 0.52, "v": range(1, 10)},
                 {"HT": 450, "AlphaT": 0.51, "v": range(1, 10)},
                 ]


class hadronicLook(supy.analysis):
    def parameters(self):
        return {"objects": self.vary({"calo":
                                      {"muon": ("muon", "Pat"),
                                       "electron": ("electron", "Pat"),
                                       "photon": ("photon", "Pat"),

                                       "jet": ("xcak5Jet", "Pat"),
                                       "jetId": "JetIDloose",
                                       "muonsInJets": False,
                                       "met": "metP4TypeIPF",
                                       "rechit": "Calo",

                                       "jetComp": ("xcak5JetPF", "Pat"),
                                       "jetIdComp": "JetIDtight",
                                       "muonsInJetsComp": True,
                                       "metComp": "metP4PF",
                                       "rechitComp": "PF",
                                       },
                                      }),
                "nJetsMinMax": self.vary({"ge2j": (2, None),
                                          #"le3j": (2, 3),
                                          #"ge4j": (4, None),
                                          }),
                "nBTagJets": self.vary({"eq0b": (0, 0),
                                        #"eq1b": (1, 1),
                                        #"eq2b": (2, 2),
                                        #"eq3b": (3, 3),
                                        #"ge4b": (4, None),
                                        }),
                "thresholds": self.vary({"200s": (200.0, 275.0,   73.3, 36.7),
                                         #"275s": (275.0, 325.0,  73.3, 36.7),
                                         #"325s": (325.0, 375.0,  86.7, 43.3),
                                         "375":  (375.0, None,  100.0, 50.0),
                                         #"675":  (675.0, None,  100.0, 50.0),
                                         #"875":  (875.0, None,  100.0, 50.0),
                                         }),
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
                calculables.muon.Indices(obj["muon"], ptMin=10, isoMax=0.20,
                                         ISO="PfIsolationR04DeltaBCorrected",
                                         ID="IdPog2012Tight",
                                         ),
                calculables.electron.Indices(obj["electron"], ptMin=10,
                                             flag2012="Veto",
                                             ),
                calculables.photon.Indices(obj["photon"], ptMin=25,
                                           flagName="photonIDRA3Pat",
                                           ),
                calculables.photon.CombinedIsoDR03RhoCorrected(obj["photon"]),
                calculables.other.RecHitSumPt(obj["rechit"]),
                calculables.other.RecHitSumP4(obj["rechit"]),
                calculables.vertex.ID(),
                calculables.vertex.Indices(),
                calculables.trigger.lowestUnPrescaledTrigger(triggers),
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
            #steps.trigger.lowestUnPrescaledTriggerFilter().onlyData(),
            #steps.trigger.hltPrescaleHistogrammer(params["triggerList"]).onlyData(),
            steps.trigger.l1Filter("L1Tech_BPTX_plus_AND_minus.v0").onlyData(),
            steps.trigger.physicsDeclaredFilter().onlyData(),
            ]+([] if params["thresholds"][1] != 275.0 else [steps.trigger.lowestUnPrescaledTriggerFilter().onlyData()])  # apply trigger in lowest HT bin

    def stepsGenValidation(self):
        return [
            #supy.steps.histos.histogrammer("genpthat",200,0,2000,title=";#hat{p_{T}} (GeV);events / bin").onlySim(),
            supy.steps.histos.histogrammer("genPartonHT", 200, 0, 1000, title=";parton H_{T} (GeV);events / bin").onlySim(),
            ]

    def stepsEvent(self, params):
        return [steps.filters.monster(),
                steps.filters.hbheNoise().onlyData(),
                supy.steps.filters.value("beamHaloCSCTightHaloId", max=0).onlyData(),
                supy.steps.filters.value("trackingFailureFilterFlag", min=1).onlyData(),
                supy.steps.filters.value("hcalLaserEventFilterFlag", min=1).onlyData(),
                supy.steps.filters.value("ecalDeadCellTPFilterFlag", min=1).onlyData(),
                #supy.steps.histos.histogrammer("logErrorTooManySeeds",    2, 0.0, 1.0, title = ";logErrorTooManySeeds;events / bin"),
                #supy.steps.histos.histogrammer("logErrorTooManyClusters", 2, 0.0, 1.0, title = ";logErrorTooManyClusters;events / bin"),
                supy.steps.histos.multiplicity("vertexIndices", max=30),
                supy.steps.filters.multiplicity("vertexIndices", min=1),
                ]

    def stepsHtLeadingJets(self, params):
        jet = params["objects"]["jet"]
        et = "Et" if params["etRatherThanPt"] else "Pt"

        ##when using full scaling
        #steps.jet.htBinFilter(jet, min = params["htBin"], max = params["htBin"]),
        #steps.jet.jetSelector(jet, params["thresholds"][2], 0),
        #steps.jet.jetSelector(jet, params["thresholds"][2], 1),

        #otherwise
        out = [steps.jet.jetPtSelector(jet, params["thresholds"][2], 0),
               steps.jet.jetPtSelector(jet, params["thresholds"][2], 1),
               steps.jet.jetEtaSelector(jet, 2.5, 0),
               supy.steps.filters.value("%sSum%s%s" % (jet[0], et, jet[1]), min=params["thresholds"][0]),
               ]

        if params["thresholds"][1] is not None:
            out.append(supy.steps.filters.value("%sSum%s%s" % (jet[0], et, jet[1]), max=params["thresholds"][1]))
        return out

    def stepsBtagJets(self, params):
        _jet = params["objects"]["jet"]
        out = [supy.steps.filters.multiplicity("%sIndicesBtagged2%s" % _jet, min=params["nBTagJets"][0], max=params["nBTagJets"][1]),
               supy.steps.histos.multiplicity("%sIndicesBtagged2%s" % _jet),
               ]
        return out

    def stepsXclean(self, params):
        return [
            #steps.other.iterHistogrammer("ecalDeadTowerTrigPrimP4", 256, 0.0, 128.0, title=";E_{T} of ECAL TP in each dead region (GeV);TPs / bin", funcString="lambda x:x.Et()"),
            supy.steps.filters.multiplicity("%sIndices%s" % params["objects"]["muon"],     max=0),
            supy.steps.filters.multiplicity("%sIndices%s" % params["objects"]["electron"], max=0),
            supy.steps.filters.multiplicity("%sIndices%s" % params["objects"]["photon"],   max=0),
            supy.steps.filters.multiplicity("%sIndicesOther%s" % params["objects"]["jet"], max=0),
            supy.steps.filters.multiplicity("%sIndices%s" % params["objects"]["jet"], min=params["nJetsMinMax"][0], max=params["nJetsMinMax"][1]),
            #steps.jet.uniquelyMatchedNonisoMuons(_jet),
            ]

    def stepsPlotsOne(self, params):
        _jet = params["objects"]["jet"]
        _met = params["objects"]["met"]
        _et = "Et" if params["etRatherThanPt"] else "Pt"
        return [

            supy.steps.histos.multiplicity("%sIndicesBtagged2%s" % _jet),
            supy.steps.histos.histogrammer("%sSum%s%s" % (_jet[0], _et, _jet[1]), 50, 0, 2500,
                                           title=";H_{T} (GeV) from %s%s %ss;events / bin" % (_jet[0], _jet[1], _et)),
            supy.steps.histos.histogrammer("%sSum%s%s" % (_jet[0], _et, _jet[1]), 60, 675, 1275,
                                           title=";H_{T} (GeV) from %s%s %ss;events / bin" % (_jet[0], _jet[1], _et)),
            supy.steps.histos.histogrammer("%sSum%s%s" % (_jet[0], _et, _jet[1]), 100, 0, 1000,
                                           title=";H_{T} (GeV) from %s%s %ss;events / bin" % (_jet[0], _jet[1], _et)),

            supy.steps.histos.histogrammer("%sSumP4%s" % _jet, 50, 0, 500, title=";MHT from %s%s (GeV);events / bin" % _jet, funcString="lambda x:x.pt()"),
            supy.steps.histos.histogrammer("vertexIndices", 20, -0.5, 19.5, title=";N vertices;events / bin", funcString="lambda x:len(x)"),
            supy.steps.histos.histogrammer("vertexSumPt", 100, 0.0, 1.0e3, title=";SumPt of 2nd vertex (GeV);events / bin", funcString="lambda x:([0.0,0.0]+sorted(x))[-2]"),

            #steps.trigger.lowestUnPrescaledTriggerHistogrammer(),
            #supy.steps.filters.label("singleJetPlots1"),
            steps.jet.singleJetHistogrammer(_jet),
            supy.steps.filters.label("jetSumPlots1"),
            steps.jet.cleanJetHtMhtHistogrammer(_jet, params["etRatherThanPt"]),
            supy.steps.histos.histogrammer("%sDeltaPhiStar%s%s" % (_jet[0], params["lowPtName"], _jet[1]), 20, 0.0, r.TMath.Pi(), title=";#Delta#phi*;events / bin", funcString='lambda x:x[0][0]'),
            supy.steps.histos.histogrammer("%sDeltaPhiStar%s" % (_jet[0], _jet[1]), 20, 0.0, r.TMath.Pi(), title=";#Delta#phi*;events / bin", funcString='lambda x:x[0][0]'),
            supy.steps.histos.histogrammer(_met, 100, 0.0, 500.0, title=";"+_met+" (GeV);events / bin", funcString="lambda x: x.pt()"),
            supy.steps.filters.label("kinematicPlots1"),

            steps.jet.alphaHistogrammer(cs=_jet, deltaPhiStarExtraName=params["lowPtName"], etRatherThanPt=params["etRatherThanPt"]),
            ]

    def stepsQcdRejection(self, params):
        _jet = params["objects"]["jet"]
        _met = params["objects"]["met"]
        _et = "Et" if params["etRatherThanPt"] else "Pt"

        return [
            supy.steps.histos.histogrammer("%sMht%sOver%s" % (_jet[0], params["highPtName"]+_jet[1], _met), 100, 0.0, 3.0,
                                           title=";MHT %s%s / %s;events / bin" % (_jet[0], params["highPtName"]+_jet[1], _met)),
            supy.steps.filters.value("%sMht%sOver%s" % (_jet[0], params["highPtName"]+_jet[1], _met), max=1.25),
            steps.other.deadEcalFilter(jets=_jet, extraName=params["lowPtName"], dR=0.3, dPhiStarCut=0.5),

            steps.jet.cleanJetHtMhtHistogrammer(_jet, params["etRatherThanPt"]),
            steps.jet.alphaHistogrammer(cs=_jet, deltaPhiStarExtraName=params["lowPtName"], etRatherThanPt=params["etRatherThanPt"]),
            #steps.jet.alphaMetHistogrammer(cs = _jet, deltaPhiStarExtraName = params["lowPtName"], etRatherThanPt = _etRatherThanPt, metName = _met),

            supy.steps.histos.histogrammer("%sRecHitSumPt" % params["objects"]["rechit"], 30, 0, 300, title=";Sum of HBHE (sev.#geq10), EB,EE (sev.#geq2) RecHit p_{T} (GeV);events / bin"),
            supy.steps.filters.value("%sRecHitSumPt" % params["objects"]["rechit"], max=30.0),

            supy.steps.filters.value("%sAlphaT%s%s" % (_jet[0], _et, _jet[1]), min=0.55),

            supy.steps.histos.histogrammer("%sMaxEmEnergyFraction%s" % (_jet[0], _jet[1]), 20, 0.0, 1.0, title=";MaxEmEnergyFraction;events / bin"),
            #supy.steps.filters.value("%sMaxEmEnergyFraction%s"%(_jet[0],_jet[1]), max = .1),
            ]

    def stepsPlotsTwo(self, params):
        _jet = params["objects"]["jet"]
        _met = params["objects"]["met"]
        return [
            supy.steps.histos.histogrammer("vertexIndices", 20, -0.5, 19.5, title=";N vertices;events / bin", funcString="lambda x:len(x)"),
            supy.steps.histos.histogrammer("%sIndices%s" % _jet, 20, -0.5, 19.5, title=";number of %s%s passing ID#semicolon p_{T}#semicolon #eta cuts;events / bin" % _jet, funcString="lambda x:len(x)"),
            steps.jet.cleanJetHtMhtHistogrammer(_jet, params["etRatherThanPt"]),
            supy.steps.histos.histogrammer("%sDeltaPhiStar%s%s" % (_jet[0], params["lowPtName"], _jet[1]), 20, 0.0, r.TMath.Pi(), title=";#Delta#phi*;events / bin", funcString='lambda x:x[0][0]'),
            supy.steps.histos.histogrammer("%sDeltaPhiStar%s" % (_jet[0], _jet[1]), 20, 0.0, r.TMath.Pi(), title=";#Delta#phi*;events / bin", funcString='lambda x:x[0][0]'),
            supy.steps.histos.histogrammer("%sMht%sOver%s" % (_jet[0], params["highPtName"]+_jet[1], _met), 100, 0.0, 3.0,
                                           title=";MHT %s%s / %s;events / bin" % (_jet[0], params["highPtName"]+_jet[1], _met)),
            ]

    def stepsHtBins(self, params):
        _jet = params["objects"]["jet"]
        return [supy.steps.filters.value("%sSumEt%s" % _jet, min=bin) for bin in [475, 575, 675, 775, 875]]

    def stepsOptional(self, params):
        return [
            supy.steps.other.skimmer(),
            #steps.other.duplicateEventCheck(),
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
                self.stepsOptional(params) +
                self.stepsBtagJets(params) +
                self.stepsQcdRejection(params) +
                self.stepsPlotsOne(params) +
                self.stepsPlotsTwo(params) +
                #self.stepsMbb(params) +
                #self.stepsDisplayer(params) +
                #self.stepsHtBins(params) +
                #self.stepsEventCount(params, label="scanAfter") +
                [])

    def listOfSampleDictionaries(self):
        sh = supy.samples.SampleHolder()
        sh.add("275_ge2b", '["/uscms/home/elaird/08_mbb/02_skim/2012_5fb_275_ge2b.root"]', lumi=5.0e3)
        sh.add("375_ge2b", '["/uscms/home/yeshaq/nobackup/supy-output/hadronicLook/375_calo_ge2/HadronicRegion.root"]', lumi=5.0e3)
        return [samples.ht17, samples.top17, samples.ewk17, samples.qcd17, sh, samples.susy17, samples.photon17]

    def listOfSamples(self, params):
        from supy.samples import specify

        def data_53X():
            jw2012 = calculables.other.jsonWeight("cert/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt")
            out = []
            out += specify(names="HTMHTParked.Run2012B-22Jan2013-v1.job649", weights=jw2012)
            out += specify(names="HTMHTParked.Run2012C-22Jan2013-v1.job649", weights=jw2012)
            out += specify(names="HTMHTParked.Run2012D-22Jan2013-v1.job649", weights=jw2012)
            #out += specify(names="HTMHTParked_ICF_sync_test", weights=jw2012)  # , nFilesMax=1, nEventsMax=20000)
            return out

        def qcd_py6(eL):
            low = map(lambda x: x[0], samples.__qcd17__.binsXs)[:-1]
            out = []
            for pt in low[low.index(80):]:
                out += specify("qcd_py6_pt_v2_%d" % pt, effectiveLumi=eL)
            return out

        def qcd_b_py6(eL):
            out = []
            out += specify("qcd_b_py6_pt_50", effectiveLumi=eL)
            out += specify("qcd_b_py6_pt_150", effectiveLumi=eL)
            return out

        def g_jets_mg(eL):
            out = []
            out += specify("g_jets_mg_ht_200_400.job501", effectiveLumi=eL, color=r.kGreen)
            out += specify("g_jets_mg_ht_400_inf.job501", effectiveLumi=eL, color=r.kGreen)
            return out

        def w_binned():
            out = []
            out += specify(names="wj_lv_mg_ht_10_150", color=r.kBlue)
            out += specify(names="wj_lv_mg_ht_150_200.job663", color=r.kGreen)
            out += specify(names="wj_lv_mg_ht_200_250.job672", color=r.kCyan)
            out += specify(names="wj_lv_mg_ht_250_300.job498", color=r.kOrange)
            out += specify(names="wj_lv_mg_ht_300_400.job498", color=r.kViolet)
            out += specify(names="wj_lv_mg_ht_400_inf.job498", color=r.kAzure)
            return out

        def z_binned():
            out = []
            out += specify("zinv_mg_ht_50_100.job407", color=r.kBlue)
            out += specify("zinv_mg_ht_100_200.job365", color=r.kGreen)
            out += specify("zinv_mg_ht_200_400.job365", color=r.kOrange)
            out += specify("zinv_mg_ht_400_inf.job365", color=r.kViolet)
            out += specify("zinv_mg_ht_50_100_ext.job500", color=r.kBlue)
            out += specify("zinv_mg_ht_100_200_ext.job680", color=r.kGreen)
            out += specify("zinv_mg_ht_200_400_ext.job500", color=r.kOrange)
            out += specify("zinv_mg_ht_400_inf_ext.job500", color=r.kViolet)
            return out

        def w_inclusive():
            out = []
            for part in range(1, 6):
                out += specify(names="wj_lv_mg_ht_incl_v2.job673_part%i" % part)
            return out

        def vv():
            out = []
            #out += specify("zinv_hbb_125_powheg.job342")
            for diBos in ["ZZ", "WZ", "WW"]:
                out += specify("%s_pythia6.job370" % diBos)
            return out

        def top():
            out = []
            out += specify(names="ttbar_powheg_v1.job410")
            for sTop in ["T_s", "T_t", "T_tW", "Tbar_s", "Tbar_t", "Tbar_tW"]:
                out += specify("%s_powheg.job368" % sTop)
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
            #w_binned() +
            #z_binned() +
            #top() +
            #vv() +
            #qcd_py6(30.0e3) +
            #qcd_b_py6(30.0e3) +
            #g_jets_mg(30e3) +
            #w_inclusive() +
            #sms() +
            []
            )

    def mergeSamples(self, org):
        def md(x, y):
            x.update(y)
            return x

        org.mergeSamples(targetSpec={"name": "2012 Data", "color": r.kBlack, "markerStyle": 20}, allWithPrefix="HT")

        mcOps = {"markerStyle": 1, "lineWidth": 3, "goptions": "hist"}

        qcdSources = []

        org.mergeSamples(targetSpec=md({"name": "QCD Multijet", "color": r.kGreen+3}, mcOps), allWithPrefix="qcd_py6")
        qcdSources = ["QCD Multijet"]

        #org.mergeSamples(targetSpec=md({"name":"QCD Multijet (b-en.)", "color":r.kGreen+3}, mcOps), allWithPrefix="qcd_b_py6")
        #qcdSources=["QCD Multijet (b-en.)"]

        #org.mergeSamples(targetSpec=md({"name":"ttz", "color": r.kYellow}, mcOps), allWithPrefix="ttz")
        #org.mergeSamples(targetSpec=md({"name":"tt", "color": r.kRed+1}, mcOps), allWithPrefix="tt_")
        #org.mergeSamples(targetSpec=md({"name":"t", "color": r.kGreen}, mcOps),
        #                 sources=["t_s_powheg.job200", "t_t_powheg.job187", "t_tw_powheg.job187", "tbar_t_powheg.job194", "tbar_tw_powheg.job187"])

        org.mergeSamples(targetSpec=md({"name": "tt/t/ttz", "color": r.kRed+1}, mcOps), sources=[
            "tt_8_mg.job315", "ttz_8_mg.job269",
            "t_s_powheg.job200", "t_t_powheg.job187", "t_tw_powheg.job187", "tbar_t_powheg.job194", "tbar_tw_powheg.job187"])
        org.mergeSamples(targetSpec=md({"name": "Z + jets", "color": r.kBlue}, mcOps), allWithPrefix="zinv_mg_ht")
        org.mergeSamples(targetSpec=md({"name": "W + jets", "color": r.kOrange-3}, mcOps), allWithPrefix="wj_lv_mg_ht_")
        org.mergeSamples(targetSpec=md({"name": "VV", "color": r.kOrange+3}, mcOps), sources=["ww_py.job188", "wz_py.job188", "zz_py.job188"])
        org.mergeSamples(targetSpec=md({"name": "ZH", "color": r.kMagenta}, mcOps), sources=["zinv_hbb_125_powheg.job342"])
        org.mergeSamples(targetSpec=md({"name": "LM6", "color": r.kMagenta}, mcOps), allWithPrefix="lm6")
        ewkSources = ["tt/t/ttz", "Z + jets", "W + jets", "VV"]

        org.mergeSamples(targetSpec=md({"name": "Standard Model ", "color": r.kAzure+6}, mcOps), sources=ewkSources + qcdSources, keepSources=True)

    def conclude(self, conf):
        org = self.organizer(conf)
        self.mergeSamples(org)
        ##for skimming only
        #utils.printSkimResults(org)
        org.scale() if not self.parameters()["signalScan"] else org.scale(100.0)
        self.makeStandardPlots(org)
        #self.makeIndividualPlots(org)
        #self.makeEfficiencyPlots(org)

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
