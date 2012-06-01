#!/usr/bin/env python

import os,copy,ROOT as r
import supy,steps,calculables,samples

class photonLook(supy.analysis) :
    def parameters(self) :
        objects = self.vary()
        fields =                              [ "jet",             "met",            "muon",        "electron",        "photon",       "rechit", "muonsInJets"]
        objects["caloJet"] = dict(zip(fields, [("xcak5Jet","Pat"), "metP4AK5TypeII",("muon","Pat"),("electron","Pat"),("photon","Pat"), "Calo" ,    False,    ]))

        return { "objects": objects,
                 "nJetsMinMax" :      self.vary(dict([ ("ge2",(2,None)),  ("2",(2,2)),  ("ge3",(3,None)) ]       [0:1] )),
                 "photonId" :         ["photonIDTightFromTwikiPat", "photonIDRA3Pat"][1],
                 "zMode" :            self.vary(dict([ ("Z",True), ("g",False) ]                                  [1:]  )),
                 "vertexMode" :       self.vary(dict([ ("vertexMode",True), ("",False) ]                         [1:2] )),
                 "subdet" :           self.vary(dict([ ("barrel", (0.0, 1.444)), ("endcap", (1.566, 2.5)) ]      [:1 ] )),
                 "jetId" :  ["JetIDloose","JetIDtight"]            [0],
                 "etRatherThanPt" : [True,False]                   [0],
                 "lowPtThreshold": 30.0,
                 "lowPtName":"lowPt",
                 "highPtThreshold" : 50.0,
                 "highPtName" : "highPt",
                 "thresholds": (475.0, None,  100.0, 50.0),

                 #required to be sorted
                 "triggerList": tuple(["HLT_Photon135_v%d"%i for i in range(4,10)]+
                                      ["HLT_Photon150_v%d"%i for i in range(1,10)]+
                                      ["HLT_Photon160_v%d"%i for i in range(1,10)]
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

    def listOfCalculables(self, params) :
        if params["vertexMode"] :
            assert params["photonId"] in ["photonIDTightFromTwikiPat"],"In vertexMode but requested %s"%params["photonId"]

        obj = params["objects"]
        _etRatherThanPt = params["etRatherThanPt"]

        return supy.calculables.zeroArgs(supy.calculables) +\
               supy.calculables.zeroArgs(calculables) +\
               supy.calculables.fromCollections(calculables.jet,[obj["jet"]]) +\
               supy.calculables.fromCollections(calculables.muon,[obj["muon"]]) +\
               supy.calculables.fromCollections(calculables.electron,[obj["electron"]]) +\
               supy.calculables.fromCollections(calculables.photon,[obj["photon"]]) +\
               [ calculables.xclean.xcJet( obj["jet"],
                                           gamma = obj["photon"],
                                           gammaDR = 0.5,
                                           muon = obj["muon"],
                                           muonDR = 0.5,
                                           correctForMuons = not obj["muonsInJets"],
                                           electron = obj["electron"], electronDR = 0.5
                                           ),
                 calculables.jet.Indices( obj["jet"], ptMin = params["thresholds"][3], etaMax = 3.0, flagName = params["jetId"]),
                 calculables.jet.Indices( obj["jet"], ptMin = params["lowPtThreshold"], etaMax = 3.0, flagName = params["jetId"], extraName = params["lowPtName"]),
                 calculables.jet.Indices( obj["jet"], ptMin = params["highPtThreshold"], etaMax = 3.0, flagName = params["jetId"], extraName = params["highPtName"]),
                 calculables.jet.IndicesBtagged2(obj["jet"], tag = "CombinedSecondaryVertexBJetTags", threshold = 0.679),

                 
                 calculables.muon.Indices( obj["muon"], ptMin = 10, combinedRelIsoMax = 0.15),
                 calculables.electron.Indices( obj["electron"], ptMin = 10, simpleEleID = "95", useCombinedIso = True),
                 calculables.photon.Indices(obj["photon"], ptMin = 25, flagName = params["photonId"]),

                 calculables.gen.genIndices( pdgs = [22], label = "Status3Photon", status = [3]),
                 calculables.gen.genMinDeltaRPhotonOther( label = "Status3Photon"),
                 
                 calculables.gen.genIndices( pdgs = [22], label = "Status1Photon", status = [1]),
                 calculables.gen.genIsolations(label = "Status1Photon", coneSize = 0.4),
                 calculables.gen.genPhotonCategory(label = "Status1Photon"),

                 calculables.xclean.IndicesUnmatched(collection = obj["photon"], xcjets = obj["jet"], DR = 0.5),
                 calculables.xclean.IndicesUnmatched(collection = obj["electron"], xcjets = obj["jet"], DR = 0.5)
                 ] \
                 + [ calculables.jet.SumP4(obj["jet"]),
                     calculables.jet.SumP4(obj["jet"], extraName = params["lowPtName"]),
                     calculables.jet.SumP4(obj["jet"], extraName = params["highPtName"]),
                     calculables.jet.DeltaPhiStar(obj["jet"], extraName = ""),
                     calculables.jet.DeltaPhiStar(obj["jet"], extraName = params["lowPtName"]),
                     calculables.jet.DeltaPseudoJet(obj["jet"], _etRatherThanPt),
                     calculables.jet.AlphaTWithPhoton1PtRatherThanMht(obj["jet"], photons = obj["photon"], etRatherThanPt = _etRatherThanPt),
                     calculables.jet.AlphaT(obj["jet"], _etRatherThanPt),
                     calculables.jet.AlphaTMet(obj["jet"], _etRatherThanPt, obj["met"]),
                     calculables.jet.MhtOverMet((obj["jet"][0], obj["jet"][1]+params["highPtName"]), met = obj["met"]),
                     calculables.jet.MhtOverMet((obj["jet"][0], obj["jet"][1]+params["highPtName"]), met = "%sPlus%s%s"%(obj["met"], obj["photon"][0], obj["photon"][1])),
                     calculables.other.metPlusParticles(met = obj["met"], particles = obj["photon"]),
                     calculables.other.minDeltaRToJet(obj["photon"], obj["jet"]),
                     calculables.other.SumP4(obj["photon"]),
                     calculables.vertex.ID(),
                     calculables.vertex.Indices(),
                     calculables.jet.deadEcalDR(obj["jet"], extraName = params["lowPtName"], minNXtals = 10),
                     calculables.trigger.lowestUnPrescaledTrigger(params["triggerList"]),
                     ]

    def listOfSteps(self,params) :
        _jet  = params["objects"]["jet"]
        _electron = params["objects"]["electron"]
        _muon = params["objects"]["muon"]
        _photon = params["objects"]["photon"]
        _met  = params["objects"]["met"]
        _etRatherThanPt = params["etRatherThanPt"]
        _et = "Et" if _etRatherThanPt else "Pt"
        
        htUpper = [supy.steps.filters.value("%sSum%s%s"%(_jet[0], _et, _jet[1]), min = params["thresholds"][1])] if params["thresholds"][1]!=None else []
        
        #event and trigger
        outList = [
            supy.steps.printer.progressPrinter(),
            supy.steps.histos.value("genpthat", 200, 0, 1000, xtitle = "#hat{p_{T}} (GeV)").onlySim(),
            steps.trigger.l1Filter("L1Tech_BPTX_plus_AND_minus.v0").onlyData(),
            steps.trigger.physicsDeclaredFilter().onlyData(),
            steps.filters.monster(),
            steps.filters.hbheNoise().onlyData(),
            steps.trigger.hltPrescaleHistogrammer(params["triggerList"]).onlyData(),
            supy.steps.filters.value("lowestUnPrescaledTrigger").onlyData(),
            steps.trigger.lowestUnPrescaledTriggerHistogrammer().onlyData(),
            steps.trigger.lowestUnPrescaledTriggerFilter().onlyData(),
            #steps.trigger.triggerScan( pattern = r"HLT_Photon\d*_v\d", prescaleRequirement = "prescale==1", tag = "Photon"),
            #steps.trigger.triggerScan( pattern = r"HLT_Photon\d*_v\d", prescaleRequirement = "True", tag = "PhotonAll"),
            ]

        if params["vertexMode"] :
            outList += [
                steps.photon.photonPtSelector(_photon, 150.0, 0),
                steps.photon.photonEtaSelector(_photon, 1.45, 0),
                supy.steps.filters.label("vertexDistribution1"),
                supy.steps.histos.histogrammer("vertexIndices", 20, -0.5, 19.5, title=";N vertices;events / bin", funcString="lambda x:len(x)"),
                supy.steps.filters.multiplicity("%sIndices%s"%_jet, min = 1),
                steps.photon.singlePhotonHistogrammer(_photon, _jet),
                ]

        #require vertex
        outList += [supy.steps.filters.multiplicity("vertexIndices", min = 1),
                    supy.steps.histos.multiplicity("vertexIndices", max = 41),
                    ]

        if params["vertexMode"] :
            return outList

        #HT bin and leading jets
        outList += [
            ##when using full scaling
            #steps.jet.htBinFilter(_jet, min = params["htBin"], max = params["htBin"]),
            #steps.jet.jetSelector(_jet, params["referenceThresholds"][0], 0),
            #steps.jet.jetSelector(_jet, params["referenceThresholds"][0], 1),
            #steps.jet.jetEtaSelector(_jet, 2.5, 0),

            #otherwise
            steps.jet.jetPtSelector(_jet, params["thresholds"][2], 0),
            steps.jet.jetPtSelector(_jet, params["thresholds"][2], 1),
            steps.jet.jetEtaSelector(_jet,2.5,0),
            supy.steps.histos.value("%sSum%s%s"%(_jet[0], _et, _jet[1]), 50, 0, 2500, xtitle = "H_{T} (GeV) from %s%s %ss"%(_jet[0], _jet[1], _et)),
            supy.steps.filters.value(var = "%sSum%s%s"%(_jet[0], _et, _jet[1]), min = params["thresholds"][0]),
            ]
        outList += htUpper

        #one photon (zero in zMode)
        if not params["zMode"] :
            outList+=[
                #supy.steps.filters.label("photonEfficiencyPlots1"),
                #steps.gen.photonEfficiencyPlots(label = "Status1Photon", ptCut = params["thresholds"]["genPhotonPtMin"],
                #                                etaCut = 1.4, isoCut = 5.0, deltaRCut = 1.1, jets = _jet, photons = _photon),

                supy.steps.filters.pt("%sP4%s"%_photon, min = 160.0, indices = "%sIndices%s"%_photon, index = 0),
                supy.steps.filters.absEta("%sP4%s"%_photon, min = params["subdet"][0], max = params["subdet"][1], indices = "%sIndices%s"%_photon, index = 0),
                steps.filters.DeltaRGreaterSelector(jets = _jet, particles = _photon, minDeltaR = 1.0, particleIndex = 0),
                
                supy.steps.filters.multiplicity("%sIndices%s"%_photon, min = 1, max = 1),

                #steps.Other.passFilter("photonEfficiencyPlots2"),
                #steps.Gen.photonEfficiencyPlots(label = "Status1Photon", ptCut = params["thresholds"]["genPhotonPtMin"],
                #                                etaCut = 1.4, isoCut = 5.0, deltaRCut = 1.1, jets = _jet, photons = _photon),
                ]

            #outList += [steps.trigger.hltTurnOnHistogrammer("photonLeadingPtPat", (100, 70, 200), t, params["tagTriggers"]) for t in params["triggerList"]]
            #outList += [steps.trigger.hltTurnOnHistogrammer("photonLeadingPtPat", (100, 70, 200), "HLT_Photon135_v4", params["tagTriggers"]),
            #            steps.trigger.hltTurnOnHistogrammer("photonLeadingPtPat", (50, 100, 200), probe = "HLT_Photon135_v5", tags = params["tagTriggers"]),
            #            steps.trigger.hltTurnOnHistogrammer("photonLeadingPtPat", (100, 70, 200), "HLT_Photon135_v6", params["tagTriggers"]),
            #            ]
        else :
            outList+=[
                supy.steps.filters.multiplicity("%sIndices%s"%_photon, max = 0),
                ]
        
        #steps.Other.histogrammer("%sSumEt%s"%_jet,50,0,1500, title = ";H_{T} (GeV) from %s%s E_{T}s;events / bin"%_jet),

        outList += [
            #bad-jet, electron, muon, vetoes
            supy.steps.filters.multiplicity("%sIndicesOther%s"%_jet, max = 0),
            supy.steps.filters.multiplicity("%sIndices%s"%_electron, max = 0),
            supy.steps.filters.multiplicity("%sIndices%s"%_muon, max = 0),
            #supy.steps.filters.multiplicity("%sIndicesOther%s"%_muon, max = 0),
            #supy.steps.filters.multiplicity("%sIndicesUnmatched%s"%_electron, max = 0),
            #supy.steps.filters.multiplicity("%sIndicesUnmatched%s"%_photon, max = 0),
            #steps.jet.uniquelyMatchedNonisoMuons(_jet),
            supy.steps.filters.multiplicity("%sIndices%s"%_jet, min = params["nJetsMinMax"][0], max = params["nJetsMinMax"][1]),
            ]

        outList+=[
            #many plots
            supy.steps.histos.multiplicity("vertexIndices", max = 41),
            supy.steps.histos.multiplicity("%sIndices%s"%_jet, max = 10),
            
            supy.steps.filters.label("jetSumPlots"),
            steps.jet.cleanJetHtMhtHistogrammer(_jet,_etRatherThanPt),
            #supy.steps.histos.pt(_met, 100, 0.0, 500.0, xtitle=_met+" (GeV)"),
            #supy.steps.histos.pt("metP4PF", 100, 0.0, 500.0, xtitle="metP4PF (GeV)"),
            supy.steps.filters.label("kinematicPlots"),
            steps.jet.alphaHistogrammer(_jet, deltaPhiStarExtraName = "", etRatherThanPt = _etRatherThanPt),
            #supy.steps.histos.value("%sAlphaTWithPhoton1PtRatherThanMht%s"%_jet, 4, 0.0, 4*0.55,
            #                        xtitle = "#alpha_{T} using photon p_{T} rather than MHT"),
            #supy.steps.histos.histogrammer(("%sAlphaTEt%s"%_jet, "%sAlphaTWithPhoton1PtRatherThanMht%s"%_jet), (25, 25), (0.50, 0.50), (1.0, 1.0),
            #                               title = ";#alpha_{T};#alpha_{T} using photon p_{T} rather than MHT;events / bin"),
            #
            #steps.jet.photon1PtOverHtHistogrammer(jets = _jet, photons = _photon, etRatherThanPt = _etRatherThanPt),
            #
            #supy.steps.filters.value("%sMhtOverHt%s"%_jet, min = 0.4),
            
            steps.jet.alphaHistogrammer(_jet, deltaPhiStarExtraName = "", etRatherThanPt = _etRatherThanPt),            
            #steps.Other.histogrammer("%sAlphaTEt%s"%_jet, 4, 0.0, 0.55*4, title=";#alpha_{T};events / bin"),
            #supy.steps.histos.histogrammer("%sIndices%s"%_jet,10,-0.5,9.5,
            #                               title=";number of %s%s passing ID#semicolon p_{T}#semicolon #eta cuts;events / bin"%_jet,
            #                               funcString="lambda x:len(x)"),
            
            supy.steps.filters.value("%sAlphaTEt%s"%_jet, min = 0.55),
            steps.trigger.lowestUnPrescaledTriggerHistogrammer().onlyData(),
            
            #supy.steps.filters.label("purityPlots2"),
            #steps.gen.photonPurityPlots("Status1Photon", _jet, _photon).onlySim(),
            
            #supy.steps.histos.histogrammer("%sMht%sOver%s" %(_jet[0], _jet[1]+params["highPtName"], _met if params["zMode"] else _met+"Plus%s%s"%_photon),
            #                               100, 0.0, 3.0,
            #                               title = ";MHT %s%s / %s;events / bin"%(_jet[0], _jet[1], _met if params["zMode"] else _met+"Plus%s%s"%_photon)),
            #supy.steps.filters.value("%sMht%sOver%s" %(_jet[0], _jet[1]+params["highPtName"], _met if params["zMode"] else _met+"Plus%s%s"%_photon),
            #                         max = 1.25),
            steps.other.deadEcalFilter(jets = _jet, extraName = params["lowPtName"], dR = 0.3, dPhiStarCut = 0.5),

            supy.steps.histos.histogrammer("%sIndices%s"%_jet,10,-0.5,9.5,
                                           title=";number of %s%s passing ID#semicolon p_{T}#semicolon #eta cuts;events / bin"%_jet,
                                           funcString="lambda x:len(x)"),

            steps.jet.cleanJetHtMhtHistogrammer(_jet,_etRatherThanPt),
            #steps.photon.singlePhotonHistogrammer(_photon, _jet, DR = "04"),
            supy.steps.histos.value("rho", 100, 0.0, 50.0),
            steps.photon.singlePhotonHistogrammer(_photon, _jet, DR = "03"),

            #supy.steps.histos.pt("%sCorrectedP4%s"%_jet, 20, 0.0, 5*params["thresholds"][2], indices = "%sIndices%s"%_jet, index = 0,
            #                     xtitle = "jet 1 p_{T} (GeV)"),
            #supy.steps.histos.pt("%sCorrectedP4%s"%_jet, 20, 0.0, 4*params["thresholds"][2], indices = "%sIndices%s"%_jet, index = 1,
            #                     xtitle = "jet 2 p_{T} (GeV)"),
            #supy.steps.histos.pt("%sCorrectedP4%s"%_jet, 20, 0.0, 2*params["thresholds"][2], indices = "%sIndices%s"%_jet, index = 2,
            #                     xtitle = "jet 3 p_{T} (GeV)"),
            #supy.steps.histos.eta("%sCorrectedP4%s"%_jet, 6, -3.0, 3.0, indices = "%sIndices%s"%_jet, index = 2, xtitle = "jet 3"),

            #steps.Other.skimmer(),
            
            #steps.Gen.genMotherHistogrammer("genIndicesPhoton", specialPtThreshold = 100.0),
            #steps.Print.eventPrinter(),
            #steps.Print.vertexPrinter(),
            #steps.Jet.jetPrinter(_jet),
            #steps.Jet.htMhtPrinter(_jet),
            #steps.Print.particleP4Printer(_photon),
            #steps.Print.alphaTPrinter(_jet,_etRatherThanPt),
            #steps.Gen.genParticlePrinter(minPt = 10.0, minStatus = 3),
            #steps.Gen.genParticlePrinter(minPt = -1.0, minStatus = 3),
            #steps.Gen.genParticlePrinter(minPt=-10.0,minStatus=1),
            #
            #steps.Displayer.displayer(jets      = _jet,
            #                          muons     = _muon,
            #                          met       = params["objects"]["met"],
            #                          electrons = params["objects"]["electron"],
            #                          photons   = params["objects"]["photon"],                            
            #                          recHits   = params["objects"]["rechit"],recHitPtThreshold=1.0,#GeV
            #                          scale     = 400.0,#GeV
            #                          etRatherThanPt = _etRatherThanPt,
            #                          #doGenParticles = True,
            #                          deltaPhiStarExtraName = params["lowPtName"],
            #                          #deltaPhiStarExtraName = "%s%s"%("","PlusPhotons"),
            #                          mhtOverMetName = "%sMht%sOver%s"%(_jet[0], _jet[1]+params["highPtName"], _met if params["zMode"] else _met+"Plus%s%s"%_photon),
            #                          ),

            #supy.steps.histos.histogrammer("%sSumEt%s"%_jet, 40, 0, 1000, title = ";H_{T} (GeV) from %s%s E_{T}s;events / bin"%_jet),
            supy.steps.histos.multiplicity("%sIndicesBtagged2%s"%_jet),
            supy.steps.histos.histogrammer("%sSumEt%s"%_jet, 6, 375.0, 975.0, title = ";H_{T} (GeV) from %s%s E_{T}s;events / bin"%_jet),
            supy.steps.histos.generic(("%sSumEt%s"%_jet, "%sIndicesBtagged2%s"%_jet), (6, 4), (375.0, -0.5), (975.0, 3.5),
                                      title = ";H_{T} (GeV);n_{b};events / bin", funcString = "lambda x:(x[0],len(x[1]))"),
            ] + [supy.steps.filters.value("%sSumEt%s"%_jet, min = 375.0+100*iBin) for iBin in range(2,6)]
        return outList

    def listOfSampleDictionaries(self) :
        return [samples.photon17]

    def listOfSamples(self,params) :
        from supy.samples import specify

        jw2012 = calculables.other.jsonWeight("cert/Cert_190456-194479_8TeV_PromptReco_Collisions12_JSON.txt")

        data  = []
        data += specify("Photon.2012A.job171",       weights = jw2012, overrideLumi = 660.1)
        data += specify("SinglePhoton.2012B.job171", weights = jw2012, overrideLumi = 890.1)

        phw = calculables.photon.photonWeight(var = "vertexIndices", weightSet = "ZM")
        mc = specify("GJets_HT400.job174", color = r.kBlue, weights = phw)
        outList = []

        if not params["zMode"] :
            outList += data
            outList += mc
        else :
            pass
            
        return outList

    def mergeSamples(self, org) :
        org.mergeSamples(targetSpec = {"name":"2012 Data", "color":r.kBlack, "markerStyle":20},
                         sources = ["Photon.2012A.job171.jsonWeight", "SinglePhoton.2012B.job171.jsonWeight"]
                         )

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
                                   ("2012 Data","GJets_HT400.job174.photonWeight"),
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
