import os
import array
import copy
import itertools

import ROOT as r
import supy
import steps
import calculables
import samples

class photonLook(supy.analysis) :
    def parameters(self) :
        objects = self.vary()
        fields =                              [       "jet",             "met",        "muon",  "electron",
                                                   "photon",          "rechit", "muonsInJets",   "jetComp",
                                                "jetIdComp", "muonsInJetsComp",     "metComp", "rechitComp",]

        objects["caloJet"] = dict(zip(fields, [("xcak5Jet","Pat"), "metP4TypeIPF", ("muon","Pat"),    ("electron","Pat"),
                                                 ("photon","Pat"),        "Calo" ,           True, ("xcak5JetPF", "Pat"),   
                                                     "JetIDtight",           True,      "metP4PF",                  "PF",]))

        return { "objects": objects,
                 "nJetsMinMax": self.vary({"ge2j": (2, None),
                                           #"le3j": (2, 3),
                                           #"ge4j": (4, None),
                                           }),
                 "nBTagJets": self.vary({"eq0b": (0, 0),
                                         #"eq1b": (1, 1),
                                         }),
                 "photonId" :         ["photonIDTightFromTwikiPat", "photonIDRA3Pat", "photonSimpleCutBased2012TightPat"][2],
                 "zMode" :            self.vary(dict([ ("Z",True), ("g",False) ]                                  [1:]  )),
                 "vertexMode" :       self.vary(dict([ ("vertexMode",True), ("",False) ]                         [1:2] )),
                 "subdet" :           self.vary(dict([ ("barrel", (0.0, 1.4442)), ("endcap", (1.566, 2.5)) ]      [:1 ] )),
                 "jetId" :  ["JetIDloose","JetIDtight"]            [0],
                 "etRatherThanPt" : [True,False]                   [0],
                 "lowPtThreshold": 30.0,
                 "lowPtName":"lowPt",
                 "highPtThreshold" : 50.0,
                 "highPtName" : "highPt",
                 "thresholds": self.vary({"375": (375.0, None,  100.0,  50.0),
                                         }),
                 #required to be sorted
                 "triggerList": tuple(#["HLT_Photon135_v%d"%i for i in range(4,10)]+
                                      #["HLT_Photon150_v%d"%i for i in range(1,5)]#+
                                      ["HLT_Photon150_v%d"%i for i in range(1,20)]+
                                      ["HLT_Photon160_v%d"%i for i in range(1,20)]
                                      ),
                 "tagTriggers": tuple(["HLT_Photon50_CaloIdVL_IsoL_v%d"%i for i in [14,15,16]]+
                                      ["HLT_Photon50_CaloIdVL_v%d"%i for i in [7,8,9]]+
                                      ["HLT_Photon75_CaloIdVL_IsoL_v%d"%i for i in [15,16,17]]+
                                      ["HLT_Photon75_CaloIdVL_v%d"%i for i in [10,11,12]]+
                                      ["HLT_Photon90_CaloIdVL_IsoL_v%d"%i for i in [12,13,14]]+
                                      ["HLT_Photon90_CaloIdVL_v%d"%i for i in [7,8,9]]
                                      ),
                 "possibleTriggers": tuple(["HLT_Photon60_CaloIdL_FJHT300_v%d"%i for i in [1,2,3]]),
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
                    calculables.jet.AlphaTWithPhoton1PtRatherThanMht(jet, photons = photon, etRatherThanPt = etRatherThanPt),
                    calculables.jet.MhtOverMet((jet[0], highPtName + jet[1]),
                                               met),
                  
                    calculables.jet.MhtOverMet((jet[0], highPtName + jet[1]),
                                               met = "%sPlus%sIndices%s"%(met, photon[0], photon[1])),
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


    def calcListOther(self, obj, triggers, photonId):
        return [calculables.xclean.IndicesUnmatched(collection=obj["photon"],
                                                    xcjets=obj["jet"], DR=0.5,
                                                    ),
                calculables.xclean.IndicesUnmatched(collection=obj["electron"],
                                                    xcjets=obj["jet"], DR=0.5,
                                                    ),
                calculables.muon.Indices(obj["muon"], ptMin = 10, isoMax = 0.12,
                                         ISO = "PfIsolationR04DeltaBCorrected",
                                         ID = "IdPog2012Tight"),
                calculables.electron.Indices( obj["electron"], ptMin = 10, etaMax = 2.5,
                                              flag2012 = "Veto"),
                calculables.photon.Indices(obj["photon"], ptMin = 25,
                                           flagName = photonId),
                calculables.photon.CombinedIsoDR03RhoCorrected(obj["photon"]),
                calculables.other.RecHitSumPt(obj["rechit"]),
                calculables.other.RecHitSumP4(obj["rechit"]),
                calculables.other.metPlusIndices(met = obj["met"], collection=obj["photon"]),
                calculables.other.SumP4(obj["photon"]),
                calculables.other.minDeltaRToJet(obj["photon"], obj["jet"]),
                calculables.other.Mt(obj["muon"], obj["met"]),

                calculables.vertex.ID(),
                calculables.vertex.Indices(),
                calculables.trigger.lowestUnPrescaledTrigger(triggers),
                ]

    def calcListGen(self):
        
        return [calculables.gen.genIndices( pdgs = [22], label = "Status3Photon", status = [3]),
                calculables.gen.genMinDeltaRPhotonOther( label = "Status3Photon"),
                calculables.gen.genIndices( pdgs = [22], label = "Status1Photon", status = [1]),
                calculables.gen.genIndices( pdgs = [13], label = "Status1MuPlus", status = [1]),
                calculables.gen.genIndices( pdgs = [-13], label = "Status1MuMinus", status = [1]),
                calculables.gen.genIndices( pdgs = [23], label = "Status3Z", status = [3]),
                calculables.gen.genIsolations(label = "Status1Photon", coneSize = 0.4),
                calculables.gen.genPhotonCategory(label = "Status1Photon"),
                ]


    def listOfCalculables(self, params) :
        if params["vertexMode"] :
            assert params["photonId"] in ["photonIDTightFromTwikiPat"],"In vertexMode but requested %s"%params["photonId"]

        obj = params["objects"]
        outList = supy.calculables.zeroArgs(supy.calculables) 
        outList += supy.calculables.zeroArgs(calculables) 
        outList += supy.calculables.fromCollections(calculables.muon,
                                                    [obj["muon"]]
                                                    ) 
        outList += supy.calculables.fromCollections(calculables.electron,
                                                    [obj["electron"]]
                                                    ) 
        outList += supy.calculables.fromCollections(calculables.photon,
                                                    [obj["photon"]]
                                                    )
        outList += self.calcListOther(obj, params["triggerList"], params["photonId"])
        outList += self.calcListJet(obj,
                                    params["etRatherThanPt"],
                                    params["thresholds"][3],
                                    params["lowPtThreshold"],
                                    params["lowPtName"],
                                    params["highPtThreshold"],
                                    params["highPtName"],
                                    params["thresholds"][0],
                                    )
        outList += self.calcListGen()

        return outList

    def stepsTrigger(self, params):
        return [steps.trigger.l1Filter("L1Tech_BPTX_plus_AND_minus.v0").onlyData(),
                steps.trigger.physicsDeclaredFilter().onlyData(),
                supy.steps.filters.value("lowestUnPrescaledTrigger").onlyData(),
                steps.trigger.lowestUnPrescaledTriggerFilter().onlyData(),
                #steps.trigger.hltPrescaleHistogrammer(params["triggerList"]).onlyData(),
                #steps.trigger.lowestUnPrescaledTriggerHistogrammer().onlyData(),
                #steps.trigger.triggerScan( pattern = r"HLT_Photon\d*_v\d", prescaleRequirement = "prescale==1", tag = "Photon"),
                #steps.trigger.triggerScan( pattern = r"HLT_Photon\d*_v\d", prescaleRequirement = "True", tag = "PhotonAll"),
                ]

    def stepsGenValidation(self, params):
        return [supy.steps.histos.value("genpthat", 200, 0, 1000, xtitle = "#hat{p_{T}} (GeV)").onlySim(),
                supy.steps.histos.histogrammer("genPartonHT", 200, 0, 1000, title=";parton H_{T} (GeV);events / bin").onlySim(),
                steps.gen.genMotherHistogrammer("genIndicesPhoton", specialPtThreshold = 100.0).onlySim(),
                steps.gen.particlePrinter(minPt = 10.0, minStatus = 3).onlySim(),
                steps.gen.particlePrinter(minPt = -1.0, minStatus = 3).onlySim(),
                steps.gen.particlePrinter(minPt=-10.0,minStatus=1).onlySim(),
                ]

    def stepsEvent(self, params):
        return [steps.filters.monster(),
                supy.steps.filters.multiplicity("%sIndices%s" % params["objects"]["jet"], min=params["nJetsMinMax"][0], max=params["nJetsMinMax"][1]),
                steps.filters.hbheNoise().onlyData(),
                supy.steps.filters.value("beamHaloCSCTightHaloId", max=0).onlyData(),
                supy.steps.filters.value("trackingFailureFilterFlag", min=1).onlyData(),
                supy.steps.filters.value("hcalLaserEventFilterFlag", min=1).onlyData(),
                supy.steps.filters.value("ecalDeadCellTPFilterFlag", min=1).onlyData(),
                supy.steps.filters.value("%sMaxEmEnergyFraction%s"%(params["objects"]["jet"]), min = .1),
                supy.steps.histos.histogrammer("%sMaxEmEnergyFraction%s" % (params["objects"]["jet"]), 20, 0.0, 1.0, title=";MaxEmEnergyFraction;events / bin"),
                #supy.steps.histos.histogrammer("logErrorTooManySeeds",    2, 0.0, 1.0, title = ";logErrorTooManySeeds;events / bin"),
                #supy.steps.histos.histogrammer("logErrorTooManyClusters", 2, 0.0, 1.0, title = ";logErrorTooManyClusters;events / bin"),
                supy.steps.histos.multiplicity("vertexIndices", max=60),
                supy.steps.filters.multiplicity("vertexIndices", min=1),
                supy.steps.histos.histogrammer("pileupTrueNumInteractions", 41, -0.5, 40.5, title=";true number of interacions;events / bin",
                                               funcString="lambda x:x[0]").onlySim(),
                ]

    def stepVertexMode(self, params):
        _jet  = params["objects"]["jet"]
        _photon = params["objects"]["photon"]
        
        return [steps.photon.photonPtSelector(_photon, 150.0, 0),
                steps.photon.photonEtaSelector(_photon, 1.45, 0),
                supy.steps.filters.label("vertexDistribution1"),
                supy.steps.histos.histogrammer("vertexIndices", 20, -0.5, 19.5, title=";N vertices;events / bin", funcString="lambda x:len(x)"),
                supy.steps.filters.multiplicity("%sIndices%s"%_jet, min = 1),
                steps.photon.singlePhotonHistogrammer(_photon, _jet),

                ]

    def stepsZMode(self, params):
        return [supy.steps.filters.multiplicity("genIndicesStatus1MuPlus", min = 1).onlySim(),
                supy.steps.filters.multiplicity("genIndicesStatus1MuMinus", min = 1).onlySim(),
                supy.steps.filters.multiplicity("genIndicesStatus3Z", min = 1, max = 1).onlySim(),
                supy.steps.filters.mass("genP4", index = 0, indices = "genIndicesStatus3Z", min = 80, max = 105).onlySim(),
                ]

    def stepsPhoton(self, params):
        _photon = params["objects"]["photon"]
        _jet = params["objects"]["jet"]
        if params["zMode"] :
            return [supy.steps.filters.multiplicity("%sIndices%s"%_photon, max = 0)]
        
        outList = [#supy.steps.filters.label("photonEfficiencyPlots1"),
                   #steps.gen.photonEfficiencyPlots(label = "Status1Photon", ptCut = params["thresholds"]["genPhotonPtMin"],
                   #                                etaCut = 1.4, isoCut = 5.0, deltaRCut = 1.1, jets = _jet, photons = _photon).onlySim(),
                   supy.steps.filters.pt("%sP4%s"%_photon, min = 165.0, indices = "%sIndices%s"%_photon, index = 0),
                   supy.steps.filters.absEta("%sP4%s"%_photon, min = params["subdet"][0], max = params["subdet"][1], indices = "%sIndices%s"%_photon, index = 0),
                   steps.filters.DeltaRGreaterSelector(jets = _jet, particles = _photon, minDeltaR = 1.0, particleIndex = 0),
                   supy.steps.filters.multiplicity("%sIndices%s"%_photon, min = 1, max = 1),
                   #supy.steps.filters.label("purityPlots2"),
                   #steps.gen.photonPurityPlots("Status1Photon", _jet, _photon).onlySim(),
                   #steps.Other.passFilter("photonEfficiencyPlots2"),
                   steps.photon.singlePhotonHistogrammer(_photon, _jet, DR = "03")]
        #outList += [steps.trigger.hltTurnOnHistogrammer("photonLeadingPtPat", (100, 70, 200), t, params["tagTriggers"]) for t in params["triggerList"]]
        #outList += [steps.trigger.hltTurnOnHistogrammer("photonLeadingPtPat", (100, 70, 200), "HLT_Photon135_v4", params["tagTriggers"]),
        #            steps.trigger.hltTurnOnHistogrammer("photonLeadingPtPat", (50, 100, 200), probe = "HLT_Photon135_v5", tags = params["tagTriggers"]),
        #            steps.trigger.hltTurnOnHistogrammer("photonLeadingPtPat", (100, 70, 200), "HLT_Photon135_v6", params["tagTriggers"]),
        #            ]

        return outList

    def stepsHtLeadingJets(self, params):
        _jet = params["objects"]["jet"]
        _et = "Et" if params["etRatherThanPt"] else "Pt"
        htUpper = [supy.steps.filters.value("%sSum%s%s"%(_jet[0], _et, _jet[1]), min = params["thresholds"][1])] if params["thresholds"][1]!=None else []
        
        outList = [steps.jet.jetPtSelector(_jet, params["thresholds"][2], 0),
                   steps.jet.jetPtSelector(_jet, params["thresholds"][2], 1),
                   steps.jet.jetEtaSelector(_jet,2.5,0),
                   supy.steps.histos.value("%sSum%s%s"%(_jet[0], _et, _jet[1]), 50, 0, 2500, xtitle = "H_{T} (GeV) from %s%s %ss"%(_jet[0], _jet[1], _et)),
                   supy.steps.filters.value(var = "%sSum%s%s"%(_jet[0], _et, _jet[1]), min = params["thresholds"][0]),
                   ]
        outList += htUpper
        outList += [supy.steps.filters.label("jetSumPlots"),
                    steps.jet.cleanJetHtMhtHistogrammer(_jet, _et),
                    supy.steps.histos.histogrammer("%sIndices%s"%_jet,10,-0.5,9.5,
                                                   title=";number of %s%s passing ID#semicolon p_{T}#semicolon #eta cuts;events / bin"%_jet,
                                                   funcString="lambda x:len(x)"),
                    steps.jet.cleanJetHtMhtHistogrammer(_jet, _et),
                    supy.steps.histos.value("rho", 100, 0.0, 50.0),]
        
        return outList

    def stepsXclean(self, params):
        _muon = params["objects"]["muon"]
        _met = params["objects"]["met"]
        return [
            #steps.other.iterHistogrammer("ecalDeadTowerTrigPrimP4", 256, 0.0, 128.0, title=";E_{T} of ECAL TP in each dead region (GeV);TPs / bin", funcString="lambda x:x.Et()"),

            supy.steps.filters.multiplicity("%sIndicesOther%s" % params["objects"]["jet"], max=0),
            supy.steps.filters.multiplicity("%sIndicesWithOddMuon%s" % params["objects"]["jet"], max=0),
            steps.jet.uniquelyMatchedNonisoMuons(params["objects"]["jet"]),
            supy.steps.filters.multiplicity("%sIndices%s" % params["objects"]["electron"], max=0),
            supy.steps.filters.multiplicity("%sIndices%s" % params["objects"]["muon"], max=0),
            supy.steps.filters.multiplicity("%sIndicesUnmatched%s" % params["objects"]["electron"], max = 0),
            supy.steps.filters.multiplicity("%sIndicesUnmatched%s" % params["objects"]["photon"], max = 0),
            steps.other.deadEcalFilter(jets=params["objects"]["jet"], extraName=params["lowPtName"], dR=0.3, dPhiStarCut=0.5),
            ]

    def stepsQcdRejection(self, params):
        _jet = params["objects"]["jet"]
        _met = params["objects"]["met"]
        _et = "Et" if params["etRatherThanPt"] else "Pt"
        _photon = params["objects"]["photon"]

        return [
            supy.steps.filters.value("%sMht%sOver%s" %(_jet[0], params["highPtName"]+_jet[1],
                                                       _met if params["zMode"] else _met+ "Plus%sIndices%s" % (_photon[0], _photon[1])), max=1.25),
            supy.steps.filters.value("%sAlphaTEt%s"%_jet, min = 0.55),
            steps.jet.cleanJetHtMhtHistogrammer(_jet, params["etRatherThanPt"]),
            steps.jet.alphaHistogrammer(cs=_jet, deltaPhiStarExtraName=params["lowPtName"], etRatherThanPt=params["etRatherThanPt"]),
            supy.steps.filters.label("kinematicPlots"),
            steps.jet.alphaHistogrammer(_jet, deltaPhiStarExtraName = "", etRatherThanPt = _et),
            supy.steps.histos.histogrammer("%sMht%sOver%s" %(_jet[0], params["highPtName"]+_jet[1], _met if params["zMode"] else _met+"Plus%sIndices%s"%_photon),
                                           100, 0.0, 3.0,
                                           title = ";MHT %s%s / %s;events / bin"%(_jet[0], _jet[1], _met if params["zMode"] else _met+"Plus%sIndices%s"%_photon)),
            ]

    def stepsBtagJets(self, params):
        _jet = params["objects"]["jet"]

        return  [supy.steps.filters.multiplicity("%sIndicesBtagged2%s" % _jet, min=params["nBTagJets"][0], max=params["nBTagJets"][1]),
                 supy.steps.histos.multiplicity("%sIndicesBtagged2%s" % _jet),
                 supy.steps.histos.generic(("%sSumEt%s"%_jet, "%sIndicesBtagged2%s"%_jet), (6, 4), (375.0, -0.5), (975.0, 3.5),
                                           title = ";H_{T} (GeV);n_{b};events / bin", funcString = "lambda x:(x[0],len(x[1]))"),
                 ]

    def stepsDisplayer(self, params):
        return [steps.displayer.displayer(jets=params["objects"]["jet"],
                                  muons=params["objects"]["muon"],
                                  met = "%sPlus%sIndices%s" % (params["objects"]["met"],
                                                               params["objects"]["photon"][0],
                                                               params["objects"]["photon"][1]),
                                  electrons=params["objects"]["electron"],
                                  photons=params["objects"]["photon"],
                                      recHits=params["objects"]["rechit"], recHitPtThreshold=1.0,  # GeV
                                  scale=400.0,  # GeV
                                  etRatherThanPt=params["etRatherThanPt"],
                                  deltaPhiStarExtraName=params["lowPtName"],
                                  deltaPhiStarCut=0.5,
                                  deltaPhiStarDR=0.3,
                                  j2Factor=params["thresholds"][2]/params["thresholds"][0],
                                  mhtOverMetName="%sMht%sOver%sPlus%sIndices%s" % (params["objects"]["jet"][0],
                                                                                   params["highPtName"]+params["objects"]["jet"][1],
                                                                                   params["objects"]["met"],
                                                                                   params["objects"]["photon"][0],
                                                                                   params["objects"]["photon"][1]),
                                  metOtherAlgo=params["objects"]["metComp"],
                                  jetsOtherAlgo=params["objects"]["jetComp"],
                                  prettyMode=False
                                  #doGenJets=True,
                                  ),]

    def stepsHtBins(self, params):
        _jet = params["objects"]["jet"]
        _et = "Et" if params["etRatherThanPt"] else "Pt"

        outList = [supy.steps.histos.histogrammer("%sSumEt%s"%_jet, 8, 375.0, 1175.0, title = ";H_{T} (GeV) from %s%s E_{T}s;events / bin"%_jet)]
        outList += [supy.steps.filters.value("%sSumEt%s" % _jet, min=bin) for bin in [475, 575, 675, 775, 875, 975, 1075]]

        return outList

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


        
    def listOfSteps(self,params) :
        outList = ([supy.steps.printer.progressPrinter()])
        outList += (self.stepsTrigger(params) +
                   self.stepsGenValidation(params) +
                   self.stepsEvent(params))
        

        if params["vertexMode"] :
            outList =+ (stepsVertexMode(params))
            return outList
        
        outList += (self.stepsHtLeadingJets(params))

        if params["zMode"]:
            outList += (self.stepsZMode(params))

        outList += (self.stepsPhoton(params) +
                    self.stepsXclean(params) +
                    self.stepsQcdRejection(params) +
                    self.stepsBtagJets(params) +
                    self.stepsDisplayer(params) +
                    self.stepsOptional(params) +
                    self.stepsHtBins(params) +
                    []
                    )

        return outList

    def listOfSampleDictionaries(self) :
        return [samples.photon17]

    def listOfSamples(self,params) :
        from supy.samples import specify
        jw2012 = calculables.other.jsonWeight("cert/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt")


        data  = []
        #data += specify(names="Photon.Run2012A-22Jan2013", weights=jw2012)#, nFilesMax=1)
        for era in ["B","C","D"][-1:]:
            data += specify(names="SinglePhotonParked.Run2012%s-22Jan2013"%era, weights=jw2012, nFilesMax=10)

        phw = calculables.photon.photonWeight(var = "vertexIndices", weightSet = "ZM")
        pu = calculables.gen.puWeight(var="pileupTrueNumInteractions", puEra="PU_Parked")
        mc = []
#        photon17.add("GJets_HT200to400", '%s/clucas/Parked13/GJets_200to400/")' % eos, xs={"LO": 960.5, "NLO":1140.78}["NLO"])
#        photon17.add("GJets_HT400toinf", '%s/clucas/Parked13/GJets_400toinf/")' % eos, xs={"LO": 107.5, "NLO":124.68}["NLO"])

        mc += specify("GJets_HT200to400", color = r.kRed, weights = pu)#, effectiveLumi=(19708.574*960.5/1140.78))#, nFilesMax=1)
        mc += specify("GJets_HT400toinf", color = r.kBlue, weights = pu)#, nFilesMax=1)#, nEventsMax=5000)
        outList = []

        if not params["zMode"] :
            outList += data
            #outList += mc
        else :
            #outList += specify("dyll_HT_10To200_M-10To50", weights=pu, nFilesMax=1, nEventsMax=5000)
            outList += specify("dyll_HT_10To200_M-50", weights=pu, nFilesMax=1)
            outList += specify("dyll_HT_200To400_M-50", weights=pu, nFilesMax=1)
            outList += specify("dyll_HT_400ToInf_M-50", weights=pu, nFilesMax=1)

            
        return outList

    def mergeSamples(self, org) :
        org.mergeSamples(targetSpec = {"name":"Data", "color":r.kBlack, "markerStyle":20},
                         sources = ["Photon.Run2012A-22Jan2013.jsonWeight"] +
                         ["SinglePhotonParked.Run2012%s-22Jan2013.jsonWeight" % x for x in ["B","C","D"]]
                         )

        #org.mergeSamples(targetSpec = {"name":"SM", "color":r.kBlue, "markerStyle":1, "lineWidth":3, "goptions":"hist"},
        #                 allWithPrefix = "GJets")

        org.mergeSamples(targetSpec = {"name":"SM", "color":r.kRed, "markerStyle":1, "lineWidth":3, "goptions":"hist"},
                         allWithPrefix = "dyll")


#    def concludeAll(self) :
#        #super(photonLook,self).concludeAll()
#
#        for item in ["275","325","375"][-1:] :
#            organizers = [self.organizer(conf) for conf in self.readyConfs if (item in conf["tag"])]
#            for org in organizers :
#                self.mergeSamples(org)
#                if "Z" in org.tag :
#                    lumi = 4529.2
#                    org.scale(lumi)
#                    print "WARNING: HARD-CODED LUMI FOR Z MODE! (%g)"%lumi
#                else :
#                    org.scale()
#            melded = supy.organizer.meld(organizers = organizers)
#            self.makeStandardPlots(melded)
#            #self.makeIndividualPlots(melded)
                                 
    def conclude(self, conf) :
        org = self.organizer(conf)
        
        ##for skimming only
        #utils.printSkimResults(org)

        self.mergeSamples(org)
        if "Z" in org.tag :
            lumi = 4529.2
            org.scale(lumi)
            print "WARNING: HARD-CODED LUMI FOR Z MODE! (%g)"%lumi
        else :
            org.scale()
            
        self.makeStandardPlots(org)
        #self.makeRootFiles(org)
        #self.makeIndividualPlots(org)
        #self.makePurityPlots(org)
        #self.makeEfficiencyPlots(org)
        #self.makeNVertexWeights(org)
        #self.makeMultiModePlots(34.7255)

    def makeStandardPlots(self, org) :
        names = [ss["name"] for ss in org.samples]
        samplesForRatios = filter(lambda x: x[0] in names and x[1] in names,
                                  [("2011 Data","standard_model_nVtx"),
                                   (".2011 Data",".standard_model_nVtx"),
                                   ("g.2011 Data","g.standard_model_nVtx"),
                                   ("g.2011 Data","g.standard_model"),
                                   ("2011 Data","standard_model"),
                                   ("2011 Data","standard_model_py6"),
                                   ("2010 Data","sm_2010"),
                                   ("2012 Data","SM"),
                                   ])

        #plot all
        pl = supy.plotter(org,
                          pdfFileName = self.pdfFileName(org.tag),
                          sampleLabelsForRatios = ("data","s.m."),
                          samplesForRatios = next(iter(samplesForRatios), ("","")),
                          blackList = ["lumiHisto","xsHisto","nJobsHisto",
                                       "deltaRGenReco",
                                       "photonMothergenPt", "photonMotherrecoPt", "photonMothermht",
                                       "quarkMothergenPt",  "quarkMotherrecoPt",  "quarkMothermht",
                                       "otherMothergenPt",  "otherMotherrecoPt",  "otherMothermht",
                                       "nGenPhotonsStatus1Photon","photonEtaStatus1Photon","photonPtStatus1Photon",
                                       "photonPhiVsEtaStatus1Photon", "photonIsoStatus1Photon",
                                       "nJetsStatus1Photon", "jetHtStatus1Photon",
                                       "nJetsPlusnPhotonsStatus1Photon", "jetHtPlusPhotonHtStatus1Photon",
                                       ],
                          doLog = False,
                          printRatios = True,
                          #latexYieldTable = True,
                          rowColors = [r.kBlack, r.kViolet+4],
                          #whiteList = ["xcak5JetIndicesPat",
                          #             #"photonPat1Pt",
                          #             #"photonPat1mhtVsPhotonPt",
                          #             "xcak5JetAlphaTFewBinsPat",
                          #             "xcak5JetAlphaTRoughPat",
                          #             "xcak5JetAlphaTWithPhoton1PtRatherThanMhtPat",
                          #             ],
                          ).plotAll()
            
    def makeIndividualPlots(self, org) :
        #plot all
        pl = supy.plotter(org,
                             pdfFileName = self.pdfFileName(org.tag),
                             showStatBox = False,
                             doLog = False,
                             anMode = True,
                             )
        plots1 = [{"plotName":"xcak5JetAlphaTFewBinsPat",
                   "stepName" :"alphaHistogrammer",
                   "stepDesc" :"xcak5JetPat",
                   "newTitle":";#alpha_{T};events / bin"},
                  
                  {"plotName":"xcak5JetIndicesPat",
                   "stepName" :"histogrammer",
                   "stepDesc" :"(lambda x:len(x))(xcak5JetIndicesPat)",
                   "newTitle":";N_{jets};events / bin"},
                  
                  {"plotName":"photonPat1Pt",
                   "stepName":"singlePhotonHistogrammer",
                   "stepDesc":"photonPat through index 0",
                   "newTitle":";photon p_{T} (GeV);events / bin",
                   "legendCoords": (0.2, 0.60, 0.5, 0.90),
                   "stamp": False,
                   "index":-1,
                   "reBinFactor":5},
                  ]

        plots2 = [{"plotName":"photonPat1mhtVsPhotonPt",
                   "stepName":"singlePhotonHistogrammer",
                   "stepDesc":"photonPat through index 0",
                   "newTitle":";photon p_{T} (GeV);MHT (GeV);events / bin",
                   "stamp": False,
                   "sampleName":"2011 Data",
                   #"index":-1,
                   },
                  ]

        plots3 = [{"plotName":"xcak5JetAlphaTRoughPat",
                   "stepName" :"variablePtGreaterFilter",
                   "stepDesc" :"xcak5JetSumP4Pat.pt()>=140.0 GeV",
                   "newTitle":";#alpha_{T};events / bin / 35 pb^{-1}"},

                  {"plotName":"photonPat1MinDRToJet",
                   "stepName" :"passFilter",
                   "stepDesc" :"singlePhotonPlots2",
                   "newTitle":";#DeltaR(photon, nearest jet);events / bin / 35 pb^{-1}",
                   "reBinFactor":3},
                  
                  {"plotName":"photonPat1SeedTime",
                   "stepName" :"passFilter",
                   "stepDesc" :"singlePhotonPlots2",
                   "newTitle":";time of photon seed crystal hit (ns);events / bin / 35 pb^{-1}",
                   "sampleWhiteList": ["2010 Data"]},
                  
                  {"plotName":"photonPat1sigmaIetaIetaBarrel",
                   "stepName" :"passFilter",
                   "stepDesc" :"singlePhotonPlots2",
                   "newTitle":";#sigma_{i#eta i#eta};events / bin / 35 pb^{-1}"},
                  
                  {"plotName":"photonPat1combinedIsolation",
                   "stepName" :"passFilter",
                   "stepDesc" :"singlePhotonPlots2",
                   "onlyDumpToFile":True},
                  ]

        newSampleNames = {"qcd_mg_nVtx": "Madgraph QCD",
                          "g_jets_mg_nVtx": "Madgraph #gamma + jets",
                          "2011 Data": "Data",
                          "standard_model_nVtx": "Standard Model",
                          }

        pl.individualPlots(plotSpecs = plots2,
                           newSampleNames = newSampleNames,
                           preliminary = True,
                           tdrStyle = False,
                           )

        pl.individualPlots(plotSpecs = plots1,
                           newSampleNames = newSampleNames,
                           preliminary = True,
                           tdrStyle = True,
                           )

    def makeRootFiles(self, org) :
        def sampleIndex(org, name) :
            for iSample,sample in enumerate(org.samples) :
                if sample["name"]==name : return iSample
            assert False, "could not find sample %s"%name


        def histo(name = "", samples = ["Data", "SM"]) :
            lst = []
            for selection in org.steps :
                if selection.name != "generic" : continue
                if selection.title!="(lambda x:(x[0],len(x[1])))(%s)"%name : continue
                dct = {}
                for s in samples :
                    dct[s] = selection[name][sampleIndex(org, s)]
                lst.append(dct)
            return lst[-1]

        dct = histo(name = "xcak5JetIndicesBtagged2Pat_vs_xcak5JetSumEtPat")

        for iBTag in range(4) :
            f = r.TFile("%s_%db.root"%(org.tag, iBTag), "RECREATE")
            f.mkdir("phot")
            f.cd("phot")
            for s in ["lumiData", "lumiMc"] :
                lumi = r.TH1D(s, s, 1, -0.5, 0.5)
                lumi.SetBinContent(1, org.lumi*1.0e-3)#/fb
                lumi.Write()

            for name,key in [("obs","Data"), ("Phot", "SM")] :
                hIn = dct[key]

                xMin   = hIn.GetXaxis().GetXmin()
                xMax   = hIn.GetXaxis().GetXmax()
                nBinsX = hIn.GetXaxis().GetNbins()

                assert abs(xMin-375.)<1.0e-6,xMin
                assert abs(xMax-975.)<1.0e-6,xMax
                assert nBinsX==6,nBinsX

                yMin   = hIn.GetYaxis().GetXmin()
                yMax   = hIn.GetYaxis().GetXmax()
                nBinsY = hIn.GetYaxis().GetNbins()

                assert abs(yMin+0.5)<1.0e-6,yMin
                assert abs(yMax-3.5)<1.0e-6,yMax
                assert nBinsY==4,nBinsY

                h1 = hIn.ProjectionX("%s_projX"%name, 1+iBTag, 1+iBTag)

                xBinsLo = array.array('d',[275., 325.]+[375.+100*i for i in range(7)])
                yBinsLo = array.array('d',[55.0, 60.0])
                hOut = r.TH2D(name, name, len(xBinsLo)-1, xBinsLo, len(yBinsLo)-1, yBinsLo)

                for iBinX in range(1,1+h1.GetNbinsX()) :
                    hOut.SetBinContent(2+iBinX, 1, h1.GetBinContent(iBinX))
                    hOut.SetBinError(2+iBinX, 1, h1.GetBinError(iBinX))

                hOut.Write()
            f.Close()

    def makeNVertexWeights(self, org, chopToOne = False) :
        def sampleIndex(org, name) :
            for iSample,sample in enumerate(org.samples) :
                if sample["name"]==name : return iSample
            assert False, "could not find sample %s"%name

        def numerAndDenom(org, var) :
            d = {}
            for selection in org.steps :
                if selection.name != "histogrammer" : continue
                if selection.title!="(lambda x:len(x))(vertexIndices)" : continue
                if var in selection :
                    sample = "2011 Data";      d["numer"] = selection[var][sampleIndex(org, sample)].Clone("%s_%s_clone"%(var, sample.replace(" ","_")))
                    sample = "standard_model"; d["denom"] = selection[var][sampleIndex(org, sample)].Clone("%s_%s_clone"%(var, sample.replace(" ","_")))
            return d

        def chop(h) :
            for iBin in range(1, 1 + h.GetNbinsX()) :
                if h.GetBinCenter(iBin)!=1.0 :
                    h.SetBinContent(iBin, 0.0)
            return
            
        keep = []
        canvas = r.TCanvas()
        canvas.SetRightMargin(0.2)
        canvas.SetTickx()
        canvas.SetTicky()
        psFileName = "%s.ps"%org.tag
        canvas.Print(psFileName+"[","Lanscape")
        for variable in ["vertexIndices"] :
            histos = numerAndDenom(org, variable)
            if "numer" not in histos or "denom" not in histos : continue

            #get relative bin heights
            result = histos["denom"].Clone("%s_oldDist"%variable)
            result.Reset()
            if chopToOne : chop(histos["numer"])
            result.Divide(histos["numer"], histos["denom"], 1.0/histos["numer"].Integral(), 1.0/histos["denom"].Integral())
            result.SetBinContent(1, 1.0); result.SetBinError(1, 0.0) #hack for zero vertex bin

            #leave MC yield unchanged
            newDist = histos["denom"].Clone("%s_newDist"%variable)
            newDist.Reset()
            newDist.Multiply(histos["denom"], result)
            result.Scale(histos["denom"].Integral()/newDist.Integral())

            #print results
            contents = []
            for iBin in range(1, 1+result.GetNbinsX()) :
                contents.append("%g:%g"%(result.GetBinCenter(iBin), result.GetBinContent(iBin)))
            print "self.weight = {%s}"%", ".join(contents)
            result.GetYaxis().SetTitle("weight to apply to MC")
            result.SetMarkerStyle(20)
            result.SetStats(False)
            result.Draw()
            canvas.Print(psFileName,"Lanscape")

        canvas.Print(psFileName+"]","Lanscape")                
        os.system("ps2pdf "+psFileName)
        os.remove(psFileName)

    def makeEfficiencyPlots(self, org) :
        def sampleIndex(org, name) :
            for iSample,sample in enumerate(org.samples) :
                if sample["name"]==name : return iSample
            assert False, "could not find sample %s"%name

        def numerAndDenom(org, var) :
            d = {}
            for selection in org.selections :
                if selection.name != "passFilter" : continue
                if   "photonEfficiencyPlots1" in selection.title : label = "denom"
                elif "photonEfficiencyPlots2" in selection.title : label = "numer"
                else : continue
                key = var+"Status1Photon"
                if key in selection :
                    d[label] = selection[key][sampleIndex(org,"g_jets_mg_v12")].Clone(label)
                
            return d

        keep = []
        canvas = r.TCanvas()
        canvas.SetRightMargin(0.2)
        canvas.SetTickx()
        canvas.SetTicky()
        psFileName = "%s.ps"%org.tag
        canvas.Print(psFileName+"[","Lanscape")
        for variable in ["photonPt","photonEta","photonIso",
                         "nJets","jetHt","nJetsPlusnPhotons","jetHtPlusPhotonHt","getMinDeltaRPhotonOtherStatus3Photon"] :
            histos = numerAndDenom(org, variable)
            if "numer" not in histos or "denom" not in histos : continue
            result = histos["numer"].Clone(variable)
            result.Reset()
            result.Divide(histos["numer"], histos["denom"], 1.0, 1.0, "b")
            result.SetMarkerStyle(20)
            result.SetStats(False)
            if result.ClassName()[2]=="1" :
                result.GetYaxis().SetRangeUser(0.0,1.0)
                result.GetYaxis().SetTitle("efficiency")
                result.Draw()
            else :
                result.GetZaxis().SetRangeUser(0.0,1.0)
                result.GetZaxis().SetTitle("efficiency")
                result.Draw("colz")
            canvas.Print(psFileName,"Lanscape")

        canvas.Print(psFileName+"]","Lanscape")                
        os.system("ps2pdf "+psFileName)
        os.remove(psFileName)

    def makePurityPlots(self, org) :
        def sampleIndex(org, name) :
            for iSample,sample in enumerate(org.samples) :
                if sample["name"]==name : return iSample
            assert False, "could not find sample %s"%name

        for selection in org.selections :
            if selection.name != "passFilter" : continue
            if "purityPlots3" not in selection.title : continue

            variables = ["genPt"]
            #variables = ["genPt","recoPt","mht"]
            for variable in variables :
                sum = None
                categories = ["photonMother","quarkMother","otherMother"]
                colors = {"photonMother":r.kBlack, "quarkMother":r.kBlue, "otherMother":r.kMagenta}
                labels = {"photonMother":"mother is photon", "quarkMother":"mother is quark", "otherMother":"mother is other"}
                histos = {}
                reBinFactor = 10
                for category in categories :
                    histoOrig = selection[category+variable][sampleIndex(org,"standard_model")]
                    histos[category] = histoOrig.Clone(histoOrig.GetName()+"clone")
                    histos[category].Rebin(reBinFactor)
                    if sum==None :
                        sum = histos[category].Clone()
                    else :
                        sum.Add(histos[category])

                canvas = r.TCanvas()
                canvas.SetRightMargin(0.2)
                null = histos[categories[0]].Clone()
                null.Reset()
                null.SetStats(False)
                null.GetYaxis().SetTitle("fraction of events")
                null.Draw()
                keep = []
                e = 0.02
                legend = r.TLegend(0.80+e, 0.60-e, 0.95+e, 0.90-e)
                for category in categories :
                    result = histos[category].Clone("result"+category)
                    result.Reset()
                    result.Divide(histos[category], sum, 1.0, 1.0, "b")
                    result.SetLineColor(colors[category])
                    result.SetMarkerColor(colors[category])
                    result.SetMarkerStyle(20)
                    result.Draw("same")
                    legend.AddEntry(result,labels[category],"p")
                    keep.append(result)

                legend.Draw()
                eps = "%s_%s_%s.eps"%(org.tag,variable,selection.title)
                pdf = eps.replace(".eps",".pdf")
                canvas.Print(eps)
                os.system("epstopdf "+eps)
                os.remove(eps)

    def makeMultiModePlots(self, lumi) :
        def org(tag) :
            org = organizer.organizer( self.sampleSpecs(tag) )
            self.mergeSamples(org, tag)
            if "Z" in tag : org.scale(1.0)
            else :              org.scale()
            return org

        def fetchPlots(someOrg, plotName) :
            outList = []
            for selection in someOrg.selections :
                if selection.name != plotName[1] : continue
                if plotName[2] not in selection.title : continue
                if plotName[0] not in selection : continue

                histos = selection[plotName[0]]
                for sample,histo in zip(someOrg.samples, histos) :
                    if "color" in sample :
                        histo.SetLineColor(sample["color"])
                        histo.SetMarkerColor(sample["color"])
                    if "markerStyle" in sample :
                        histo.SetMarkerStyle(sample["markerStyle"])
                    outList.append( (sample["name"], histo) )
            return copy.deepcopy(outList)

        def histoDict(plotName, samples) :
            l = []
            for tag in self.sideBySideAnalysisTags() :
                l += fetchPlots(org(tag), plotName)

            d = {}
            names = [item[0] for item in l]
            for item in l :
                name = item[0]
                assert names.count(name)==1, "found %d instances of sample %s"%(names.count(name),name)
                if name in samples :
                    d[name] = item[1]
            return d

        plotName = ("xcak5JetAlphaTEtPat", "variablePtGreaterFilter", "xcak5JetSumP4Pat.pt()>=140.0 GeV")
        samples = ["z_inv_mg_v12_skim", "g_jets_mg_v12", "2010 Data"]
        sampleNames = ["Z -> #nu #bar{#nu}", " #gamma + jets", "Data"]
        styles  = [28,                   29,             20]
        d = histoDict(plotName, samples)

        canvas = r.TCanvas("canvas","canvas",500,500)
        canvas.SetTickx()
        canvas.SetTicky()

        legend = r.TLegend(0.65, 0.65, 0.85, 0.85)
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)

        first = True
        for iSample,sample in enumerate(samples) :
            histo = d[sample]
            histo.SetMarkerStyle(styles[iSample])
            histo.SetStats(False)
            histo.Scale( 1.0/histo.Integral(0,histo.GetNbinsX()+1) )
            legend.AddEntry(histo,sampleNames[iSample],"pl")
            if first :
                histo.Draw()
                histo.GetYaxis().SetRangeUser(0.0,1.0)
                histo.GetYaxis().SetTitle("fraction of events / bin")
                histo.GetYaxis().SetTitleOffset(1.25)
            else :
                histo.Draw("same")
            first = False

        legend.Draw()

        utils.cmsStamp(lumi)
        
        eps = "alphaTCompare.eps"
        canvas.Print(eps)
        os.system("epstopdf "+eps)
        os.remove(eps)
