import supy,steps,calculables,samples,os, ROOT as r

def triggerTuple(l = [], keys = []) :
    out = []
    for item in l :
        stem = "HLT"
        for key in keys :
            stem += "_%s%s"%(key, str(item[key]).replace(".","p"))
        for version in item["v"] :
            out.append("%s_v%d"%(stem, version))
    return tuple(out)

#required to be sorted
triggers_alphaT_2012 = triggerTuple(l  = [{"HT":250, "AlphaT": 0.55, "v":range(1,10)},
                                          
                                          #{"HT":300, "AlphaT": 0.53, "v":range(1,4)},
                                          #{"HT":300, "AlphaT": 0.54, "v":range(1,10)},
                                          #
                                          #{"HT":350, "AlphaT": 0.52, "v":range(1,4)},
                                          #{"HT":350, "AlphaT": 0.53, "v":range(1,15)},
                                          #
                                          #{"HT":400, "AlphaT": 0.51, "v":range(1,15)},
                                          #{"HT":400, "AlphaT": 0.52, "v":range(1,10)},
                                          #
                                          #{"HT":450, "AlphaT": 0.51, "v":range(1,10)},
                                          ], keys = ("HT", "AlphaT"))


class hadronicLook(supy.analysis) :
    def parameters(self) :
        objects = self.vary()
        fields =                           [ "jet",                        "jetId",     "muonsInJets",           "met",
                                             "compJet",                "compJetId", "compMuonsInJets",       "compMet",
                                             "muon",                    "electron",          "photon",         "rechit"]

        objects["calo"] = dict(zip(fields, [("xcak5Jet","Pat"),       "JetIDloose",             False,   "metP4TypeIPF",
                                            ("xcak5JetPF","Pat"),     "JetIDtight",              True,        "metP4PF",
                                            ("muon","Pat"),     ("electron","Pat"),  ("photon","Pat"),           "Calo",
                                            ]))
        
        return { "objects": objects,
                 "nJetsMinMax" :      self.vary(dict([ ("ge2",(2,None)),  ("2",(2,2)),  ("ge3",(3,None)),  ("3",(3,3)) ]       [0:1] )),
                 "etRatherThanPt" : [True,False][0],
                 "lowPtThreshold" : 30.0,
                 "lowPtName" : "lowPt",
                 "highPtThreshold" : 50.0,
                 "highPtName" : "highPt",
                 "tanBeta" : [None, 3, 10, 50][0],
                 "thresholds": self.vary(dict( [("275",        (275.0, 325.0, 100.0, 50.0)),#0
                                                ("325",        (325.0, 375.0, 100.0, 50.0)),#1
                                                ("375",        (375.0, None,  100.0, 50.0)),#2
                                                ("325_scaled", (325.0, 375.0,  86.7, 43.3)),#3
                                                ("275_scaled", (275.0, 325.0,  73.3, 36.7)),#4
                                                ("675",        (675.0, None,  100.0, 50.0)),#5
                                                ][4:5] )),
                 "triggerList": triggers_alphaT_2012, 
                 }

    def calcListJet(self, obj, etRatherThanPt, ptMin, lowPtThreshold, lowPtName, highPtThreshold, highPtName, htThreshold) :
        def calcList(jet, met, photon, muon, electron, muonsInJets, jetIdFlag) :
            outList = [
                calculables.xclean.xcJet(jet,
                                         applyResidualCorrectionsToData = False,
                                         gamma = photon,
                                         gammaDR = 0.5,
                                         muon = muon,
                                         muonDR = 0.5,
                                         correctForMuons = not muonsInJets,
                                         electron = electron,
                                         electronDR = 0.5),
                calculables.jet.Indices( jet, ptMin = ptMin,           etaMax = 3.0, flagName = jetIdFlag),
                calculables.jet.Indices( jet, ptMin = lowPtThreshold,  etaMax = 3.0, flagName = jetIdFlag, extraName = lowPtName),
                calculables.jet.Indices( jet, ptMin = highPtThreshold, etaMax = 3.0, flagName = jetIdFlag, extraName = highPtName),
                calculables.jet.IndicesBtagged2(jet, tag = "CombinedSecondaryVertexBJetTags", threshold = 0.679),
                
                calculables.jet.SumP4(jet),
                calculables.jet.SumP4(jet, extraName = lowPtName ),
                calculables.jet.SumP4(jet, extraName = highPtName),
                calculables.jet.SumP4(jet, extraName = "Btagged2"),
                calculables.jet.DeltaPhiStar(jet, extraName = lowPtName),
                calculables.jet.DeltaPhiStar(jet),
                calculables.jet.MaxEmEnergyFraction(jet),
                calculables.jet.DeltaPseudoJet(jet, etRatherThanPt),
                calculables.jet.AlphaT(jet, etRatherThanPt),
                calculables.jet.AlphaTMet(jet, etRatherThanPt, met),
                calculables.jet.MhtOverMet((jet[0], highPtName+jet[1]), met),
                calculables.jet.DeadEcalDR(jet, extraName = lowPtName, minNXtals = 10),
                supy.calculables.other.fixedValue("%sFixedHtBin%s"%jet, htThreshold),
                ]
            return outList+supy.calculables.fromCollections(calculables.jet, [jet])

        outList = calcList(obj["jet"], obj["met"], obj["photon"], obj["muon"], obj["electron"], obj["muonsInJets"], obj["jetId"])
        if all([("comp"+item in obj) for item in ["Jet", "Met","MuonsInJets","JetId"]]) :
            outList += calcList(obj["compJet"], obj["compMet"], obj["photon"], obj["muon"], obj["electron"], obj["compMuonsInJets"], obj["compJetId"])
        return outList

    def calcListOther(self, obj, triggers) :
        return [
            calculables.xclean.IndicesUnmatched(collection = obj["photon"], xcjets = obj["jet"], DR = 0.5),
            calculables.xclean.IndicesUnmatched(collection = obj["electron"], xcjets = obj["jet"], DR = 0.5),

            calculables.muon.Indices( obj["muon"], ptMin = 10, isoMax = 0.20, ISO = "PfIsolationR04DeltaBCorrected", ID = "IdPog2012Tight"),
            calculables.electron.Indices( obj["electron"], ptMin = 10, flag2012 = "Veto"),
            calculables.photon.Indices(obj["photon"], ptMin = 25, flagName = "photonIDRA3Pat"),
            calculables.photon.CombinedIsoDR03RhoCorrected(obj["photon"]),

            calculables.other.RecHitSumPt(obj["rechit"]),
            calculables.other.RecHitSumP4(obj["rechit"]),
            
            calculables.vertex.ID(),
            calculables.vertex.Indices(),
            calculables.trigger.lowestUnPrescaledTrigger(triggers),
            ]
    
    def listOfCalculables(self, params) :
        obj = params["objects"]
        outList  = supy.calculables.zeroArgs(supy.calculables)
        outList += supy.calculables.zeroArgs(calculables)
        outList += supy.calculables.fromCollections(calculables.muon, [obj["muon"]])
        outList += supy.calculables.fromCollections(calculables.electron, [obj["electron"]])
        outList += supy.calculables.fromCollections(calculables.photon, [obj["photon"]])
        outList += self.calcListOther(obj, params["triggerList"])
        outList += self.calcListJet(obj, params["etRatherThanPt"], params["thresholds"][3],
                                    params["lowPtThreshold"], params["lowPtName"], params["highPtThreshold"], params["highPtName"], params["thresholds"][0])

        outList += [calculables.gen.genIndices( pdgs = [-5,5], label = "Status3b", status = [3]),
                    calculables.gen.genIndices( pdgs = [-5,5], label = "Status3bZDaughters", status = [3], motherPdgs = [23]),
                    ]
        return outList

    def stepsTrigger(self, params) :
        return [
            #steps.trigger.lowestUnPrescaledTriggerFilter().onlyData(),
            #steps.trigger.hltPrescaleHistogrammer(params["triggerList"]).onlyData(),
            steps.trigger.l1Filter("L1Tech_BPTX_plus_AND_minus.v0").onlyData(),
            steps.trigger.physicsDeclaredFilter().onlyData(),
            ]+([] if params["thresholds"][1]!=325.0 else [steps.trigger.lowestUnPrescaledTriggerFilter().onlyData()]) #apply trigger in lowest HT bin

    def stepsGenValidation(self) :
        return [supy.steps.histos.histogrammer("genpthat",200,0,2000,title=";#hat{p_{T}} (GeV);events / bin").onlySim(),
                #supy.steps.histos.histogrammer("genPartonHT",200,0,1000,title=";parton H_{T} (GeV);events / bin").onlySim(),
                ]

    def stepsEvent(self) :
        return [steps.filters.monster(),
                steps.filters.hbheNoise().onlyData(),
                #supy.steps.filters.value("beamHaloCSCTightHaloId").invert().onlyData(),
                #supy.steps.filters.value("trackingFailureFilterFlag").onlyData(),
                #supy.steps.filters.value("hcalLaserEventFilterFlag").onlyData(),
                #supy.steps.filters.value("ecalDeadCellTPFilterFlag").onlyData(),
                #supy.steps.histos.histogrammer("logErrorTooManySeeds",    2, 0.0, 1.0, title = ";logErrorTooManySeeds;events / bin"),
                #supy.steps.histos.histogrammer("logErrorTooManyClusters", 2, 0.0, 1.0, title = ";logErrorTooManyClusters;events / bin"),
                supy.steps.histos.multiplicity("vertexIndices", max = 30),
                supy.steps.filters.multiplicity("vertexIndices", min = 1),
                ]

    def stepsHtLeadingJets(self, params) :
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
               supy.steps.filters.value("%sSum%s%s"%(jet[0], et, jet[1]), min = params["thresholds"][0]),
               ]

        if params["thresholds"][1]!=None :
            out.append(supy.steps.filters.value("%sSum%s%s"%(jet[0], et, jet[1]), max = params["thresholds"][1]))
        return out

    def stepsXclean(self, params) :
        return [
            #steps.other.iterHistogrammer("ecalDeadTowerTrigPrimP4", 256, 0.0, 128.0, title=";E_{T} of ECAL TP in each dead region (GeV);TPs / bin", funcString="lambda x:x.Et()"),
            supy.steps.filters.multiplicity("%sIndices%s"%params["objects"]["muon"],     max = 0),
            supy.steps.filters.multiplicity("%sIndices%s"%params["objects"]["electron"], max = 0),
            supy.steps.filters.multiplicity("%sIndices%s"%params["objects"]["photon"],   max = 0),
            supy.steps.filters.multiplicity("%sIndicesOther%s"%params["objects"]["jet"], max = 0),
            supy.steps.filters.multiplicity("%sIndices%s"%params["objects"]["jet"], min = params["nJetsMinMax"][0], max = params["nJetsMinMax"][1]),
            #steps.jet.uniquelyMatchedNonisoMuons(_jet),
            ]

    def stepsPlotsOne(self, params) :
        _jet = params["objects"]["jet"]
        _met = params["objects"]["met"]
        _et = "Et" if params["etRatherThanPt"] else "Pt"
        return [
            supy.steps.histos.multiplicity("%sIndicesBtagged2%s"%_jet),
            supy.steps.histos.histogrammer("%sSum%s%s"%(_jet[0], _et, _jet[1]), 50, 0, 2500,
                                           title = ";H_{T} (GeV) from %s%s %ss;events / bin"%(_jet[0], _jet[1], _et)),
            supy.steps.histos.histogrammer("%sSum%s%s"%(_jet[0], _et, _jet[1]), 60, 675, 1275,
                                           title = ";H_{T} (GeV) from %s%s %ss;events / bin"%(_jet[0], _jet[1], _et)),
            supy.steps.histos.histogrammer("%sSum%s%s"%(_jet[0], _et, _jet[1]), 100, 0, 1000,
                                           title = ";H_{T} (GeV) from %s%s %ss;events / bin"%(_jet[0], _jet[1], _et)),
            
            supy.steps.histos.histogrammer("%sSumP4%s"%_jet, 50, 0, 500, title = ";MHT from %s%s (GeV);events / bin"%_jet, funcString = "lambda x:x.pt()"),
            supy.steps.histos.histogrammer("vertexIndices", 20, -0.5, 19.5, title=";N vertices;events / bin", funcString="lambda x:len(x)"),
            supy.steps.histos.histogrammer("vertexSumPt", 100, 0.0, 1.0e3, title = ";SumPt of 2nd vertex (GeV);events / bin", funcString = "lambda x:([0.0,0.0]+sorted(x))[-2]"),
            
            #steps.trigger.lowestUnPrescaledTriggerHistogrammer(),
            #supy.steps.filters.label("singleJetPlots1"),
            #steps.jet.singleJetHistogrammer(_jet),
            supy.steps.filters.label("jetSumPlots1"), 
            steps.jet.cleanJetHtMhtHistogrammer(_jet,params["etRatherThanPt"]),
            supy.steps.histos.histogrammer("%sDeltaPhiStar%s%s"%(_jet[0], params["lowPtName"], _jet[1]), 20, 0.0, r.TMath.Pi(), title = ";#Delta#phi*;events / bin", funcString = 'lambda x:x[0][0]'),
            supy.steps.histos.histogrammer("%sDeltaPhiStar%s"%(_jet[0], _jet[1]), 20, 0.0, r.TMath.Pi(), title = ";#Delta#phi*;events / bin", funcString = 'lambda x:x[0][0]'),
            supy.steps.histos.histogrammer("%sMaxEmEnergyFraction%s"%(_jet[0], _jet[1]), 20, 0.0, 1.0, title = ";MaxEmEnergyFraction;events / bin"),
            supy.steps.histos.histogrammer(_met,100,0.0,500.0,title=";"+_met+" (GeV);events / bin", funcString = "lambda x: x.pt()"),
            supy.steps.filters.label("kinematicPlots1"),

            steps.jet.alphaHistogrammer(cs = _jet, deltaPhiStarExtraName = params["lowPtName"], etRatherThanPt = params["etRatherThanPt"]),
            ]

    def stepsQcdRejection(self, params) :
        _jet = params["objects"]["jet"]
        _met = params["objects"]["met"]
        _et = "Et" if params["etRatherThanPt"] else "Pt"

        return [
            supy.steps.histos.histogrammer("%sMht%sOver%s"%(_jet[0],params["highPtName"]+_jet[1],_met), 100, 0.0, 3.0,
                                     title = ";MHT %s%s / %s;events / bin"%(_jet[0],params["highPtName"]+_jet[1],_met)),
            supy.steps.filters.value("%sMht%sOver%s"%(_jet[0],params["highPtName"]+_jet[1],_met), max = 1.25),
            steps.other.deadEcalFilter(jets = _jet, extraName = params["lowPtName"], dR = 0.3, dPhiStarCut = 0.5),
            
            steps.jet.cleanJetHtMhtHistogrammer(_jet, params["etRatherThanPt"]),
            steps.jet.alphaHistogrammer(cs = _jet, deltaPhiStarExtraName = params["lowPtName"], etRatherThanPt = params["etRatherThanPt"]),
            #steps.jet.alphaMetHistogrammer(cs = _jet, deltaPhiStarExtraName = params["lowPtName"], etRatherThanPt = _etRatherThanPt, metName = _met),

            supy.steps.histos.histogrammer("%sRecHitSumPt"%params["objects"]["rechit"], 30, 0, 300, title = ";Sum of HBHE (sev.#geq10), EB,EE (sev.#geq2) RecHit p_{T} (GeV);events / bin"),
            supy.steps.filters.value("%sRecHitSumPt"%params["objects"]["rechit"], max = 30.0),

            supy.steps.filters.value("%sAlphaT%s%s"%(_jet[0], _et, _jet[1]), min = 0.55),

            #supy.steps.histos.histogrammer("%sMaxEmEnergyFraction%s"%(_jet[0], _jet[1]), 20, 0.0, 1.0, title = ";MaxEmEnergyFraction;events / bin"),
            #supy.steps.filters.value("%sMaxEmEnergyFraction%s"%(_jet[0],_jet[1]), max = .1),
            ]

    def stepsPlotsTwo(self, params) :
        _jet = params["objects"]["jet"]
        _met = params["objects"]["met"]
        return [
            supy.steps.histos.histogrammer("vertexIndices", 20, -0.5, 19.5, title=";N vertices;events / bin", funcString="lambda x:len(x)"),
            supy.steps.histos.histogrammer("%sIndices%s"%_jet, 20, -0.5, 19.5, title=";number of %s%s passing ID#semicolon p_{T}#semicolon #eta cuts;events / bin"%_jet, funcString="lambda x:len(x)"),
            steps.jet.cleanJetHtMhtHistogrammer(_jet,params["etRatherThanPt"]),
            supy.steps.histos.histogrammer("%sDeltaPhiStar%s%s"%(_jet[0], params["lowPtName"], _jet[1]), 20, 0.0, r.TMath.Pi(), title = ";#Delta#phi*;events / bin", funcString = 'lambda x:x[0][0]'),
            supy.steps.histos.histogrammer("%sDeltaPhiStar%s"%(_jet[0], _jet[1]), 20, 0.0, r.TMath.Pi(), title = ";#Delta#phi*;events / bin", funcString = 'lambda x:x[0][0]'),
            supy.steps.histos.histogrammer("%sMht%sOver%s"%(_jet[0],params["highPtName"]+_jet[1],_met), 100, 0.0, 3.0,
                                     title = ";MHT %s%s / %s;events / bin"%(_jet[0],params["highPtName"]+_jet[1],_met)),
            ]

    def stepsHtBins(self, params) :
        _jet = params["objects"]["jet"]
        return [supy.steps.filters.value("%sSumEt%s"%_jet, min = bin) for bin in [475, 575, 675, 775, 875]]

    def stepsOptional(self, params) :
        return [
            #supy.steps.other.skimmer(),
            #steps.other.duplicateEventCheck(),
            #steps.other.pickEventSpecMaker(),
            #steps.other.cutBitHistogrammer(self.togglePfJet(_jet), self.togglePfMet(_met)),
            #steps.Print.jetPrinter(_jet),
            #steps.Print.particleP4Printer(_muon),
            #steps.Print.particleP4Printer(_photon),
            #steps.Gen.genParticlePrinter(minPt = 10.0, minStatus = 3),
            ]

    def stepsDisplayer(self, params) :
        jet = params["objects"]["jet"]
        return [
            steps.displayer.displayer(jets      = jet,
                                      muons     = params["objects"]["muon"],
                                      met       = params["objects"]["met"],
                                      electrons = params["objects"]["electron"],
                                      photons   = params["objects"]["photon"],
                                      recHits   = params["objects"]["rechit"], recHitPtThreshold = 1.0,#GeV
                                      scale = 400.0,#GeV
                                      etRatherThanPt = params["etRatherThanPt"],
                                      deltaPhiStarExtraName = params["lowPtName"],
                                      deltaPhiStarCut = 0.5,
                                      deltaPhiStarDR = 0.3,
                                      j2Factor = params["thresholds"][2]/params["thresholds"][0],
                                      mhtOverMetName = "%sMht%sOver%s"%(jet[0],params["highPtName"]+jet[1],params["objects"]["met"]),
                                      metOtherAlgo  = params["objects"]["compMet"],
                                      jetsOtherAlgo = params["objects"]["compJet"],
                                      #doGenJets = True,
                                      ),
            ]

    def stepsMbb(self, params) :
        _jet = params["objects"]["jet"]
        return [
            #supy.steps.filters.multiplicity("genIndicesStatus3b", min = 4, max = 4),
            #steps.printer.eventPrinter(),
            #steps.printer.jetPrinter(_jet),
            #steps.gen.particlePrinter(),
            supy.steps.filters.multiplicity("%sIndicesBtagged2%s"%_jet, min = 2, max = 2),
            supy.steps.histos.multiplicity("%sIndicesBtagged2%s"%_jet),
            supy.steps.histos.eta("%sCorrectedP4%s"%_jet, 24, -3.0, 3.0, indices = "%sIndicesBtagged2%s"%_jet, index = 0, xtitle = "b jet 0"),
            supy.steps.histos.eta("%sCorrectedP4%s"%_jet, 24, -3.0, 3.0, indices = "%sIndicesBtagged2%s"%_jet, index = 1, xtitle = "b jet 1"),
            supy.steps.histos.mass("%sSumP4Btagged2%s"%_jet, 24, 0.0, 1200.0, xtitle = "sum P4 {b jets}"),
            supy.steps.histos.pt("%sSumP4Btagged2%s"%_jet, 24, 0.0, 1200.0, xtitle = "sum P4 {b jets}"),
            supy.steps.filters.mass("%sSumP4Btagged2%s"%_jet, min = 450.0),
            supy.steps.histos.pt("%sSumP4Btagged2%s"%_jet, 24, 0.0, 1200.0, xtitle = "sum P4 {b jets}"),
            #steps.jet.mbbHistogrammer(_jet, drMatch = 0.2, bZDaughters = "genIndicesStatus3bZDaughters"),
            ]

    def listOfSteps(self, params) :
        _jet  = params["objects"]["jet"]
        _met  = params["objects"]["met"]
        _etRatherThanPt = params["etRatherThanPt"]
        _et = "Et" if _etRatherThanPt else "Pt"

        return ([supy.steps.printer.progressPrinter()] +
                #self.stepsGenValidation() +
                self.stepsEvent() +
                self.stepsTrigger(params) +
                self.stepsHtLeadingJets(params) +
                self.stepsXclean(params) +
                #self.stepsPlotsOne(params) +
                self.stepsQcdRejection(params) +
                self.stepsPlotsTwo(params) +
                self.stepsMbb(params) +
                #self.stepsDisplayer(params) +
                #self.stepsOptional(params) +
                #self.stepsHtBins(params) +
                [])

    def listOfSampleDictionaries(self) :
        sh = supy.samples.SampleHolder()
        sh.add("275_ge2b", '["/uscms/home/elaird/08_mbb/02_skim/2012_5fb_275_ge2b.root"]', lumi = 5.0e3)
        return [samples.ht17, samples.top17, samples.ewk17, samples.qcd17, sh]
    
    def listOfSamples(self,params) :
        from supy.samples import specify

        def data_52X() :
            jw2012 = calculables.other.jsonWeight("cert/Cert_190456-196531_8TeV_PromptReco_Collisions12_JSON.txt")

            out = []
            out += specify(names = "HT.Run2012A-PromptReco-v1.AOD.job229", weights = jw2012, overrideLumi = 707.3810)
            out += specify(names = "HTMHT.Run2012B-PromptReco-v1.AOD.job228", weights = jw2012, overrideLumi = 3354.0000)
            out += specify(names = "HTMHT.Run2012B-PromptReco-v1.AOD.job238",  weights = jw2012, overrideLumi = 923.7680)
            return out

        def data_52X_2b_skim() :
            return specify(names = "275_ge2b")

        def qcd_py6(eL) :
            low = map(lambda x:x[0],samples.__qcd17__.binsXs)[:-1]
            out = []
            for pt in low[low.index(80):] :
                out += specify("qcd_py6_pt_%d"%pt, effectiveLumi = eL)
            return out

        def qcd_b_py6(eL) :
            out = []
            out += specify("qcd_b_py6_pt_50", effectiveLumi = eL)
            out += specify("qcd_b_py6_pt_150", effectiveLumi = eL)
            return out

        def g_jets_mg(eL) :
            gM = [40,100,200]
            return specify( effectiveLumi = eL, color = r.kGreen,
                            names = [("g_jets_mg_ht_%d_%d")[:None if t[1] else -2] for t in zip(gM,gM[1:]+["inf"])] )

        def w_binned() :
            out = []
            #out += specify(names = "wj_lv_mg_ht_0_250_other_reqs", nFilesMax = 1, nEventsMax = 20000, color = r.kRed)
            out += specify(names = "wj_lv_mg_ht_250_300", color = r.kBlue)
            out += specify(names = "wj_lv_mg_ht_300_400", color = r.kGreen)
            out += specify(names = "wj_lv_mg_ht_400_inf", color = r.kCyan)
            return out

        def z_binned() :
            out = []
            out += specify("zinv_mg_ht_50_100.job214")
            out += specify("zinv_mg_ht_100_200.job234")
            out += specify("zinv_mg_ht_200_400.job233")
            out += specify("zinv_mg_ht_400_inf.job213")
            return out

        def w_inclusive() :
            return specify(names = "wj_lv_mg_ht_incl", color = r.kOrange)

        def vv() :
            out = []
            out += specify("ww_py.job188")
            out += specify("wz_py.job188")
            out += specify("zz_py.job188")
            out += specify("zinv_hbb_125_powheg.job342")
            return out

        def top() :
            out = []
            #out += specify(names = "tt_8_mg.job188")
            out += specify(names = "tt_8_mg.job315")
            out += specify(names = "ttz_8_mg.job269", nFilesMax = 1)

            out += specify("t_s_powheg.job200"    )
            #out += specify("t_t_powheg.job187"    ) #low MC stats
            out += specify("t_tw_powheg.job187"   )
            out += specify("tbar_t_powheg.job194" )
            out += specify("tbar_tw_powheg.job187")
            
            return out

        def susy(eL) :
            return specify(names = "lm6", effectiveLumi = eL, color = r.kRed)

        return (
            #data_52X() +
            #data_52X_2b_skim() +
            w_binned() +
            z_binned() +
            top() +
            vv() +
            #qcd_py6(30.0e3) +
            #qcd_b_py6(30.0e3) +
            ##w_inclusive() +
            []
            )

    def mergeSamples(self, org) :
        def md(x, y) :
            x.update(y)
            return x
        
        org.mergeSamples(targetSpec = {"name":"2012 Data", "color":r.kBlack, "markerStyle":20}, allWithPrefix = "HT")
        #org.mergeSamples(targetSpec = {"name":"2012 Data", "color":r.kBlack, "markerStyle":20}, allWithPrefix = "275_ge2b")

        mcOps = {"markerStyle":1, "lineWidth":3, "goptions":"hist"}

        qcdSources = []

        org.mergeSamples(targetSpec = md({"name":"QCD Multijet", "color":r.kGreen+3}, mcOps), allWithPrefix = "qcd_py6")
        qcdSources = ["QCD Multijet"]

        #org.mergeSamples(targetSpec = md({"name":"QCD Multijet (b-en.)", "color":r.kGreen+3}, mcOps), allWithPrefix = "qcd_b_py6")
        #qcdSources = ["QCD Multijet (b-en.)"]

        #org.mergeSamples(targetSpec = md({"name":"ttz", "color": r.kYellow}, mcOps), allWithPrefix = "ttz")
        #org.mergeSamples(targetSpec = md({"name":"tt", "color": r.kRed+1}, mcOps), allWithPrefix = "tt_")
        #org.mergeSamples(targetSpec = md({"name":"t", "color": r.kGreen}, mcOps),
        #                 sources = ["t_s_powheg.job200", "t_t_powheg.job187", "t_tw_powheg.job187", "tbar_t_powheg.job194", "tbar_tw_powheg.job187"])
        
        org.mergeSamples(targetSpec = md({"name":"tt/t/ttz", "color":r.kRed+1}, mcOps), sources = [
            "tt_8_mg.job188", "ttz_8_mg.job269",
            "t_s_powheg.job200", "t_t_powheg.job187", "t_tw_powheg.job187", "tbar_t_powheg.job194", "tbar_tw_powheg.job187"])
        org.mergeSamples(targetSpec = md({"name":"Z + jets", "color": r.kBlue}, mcOps), allWithPrefix = "zinv_mg_ht")
        org.mergeSamples(targetSpec = md({"name":"W + jets", "color": r.kOrange-3}, mcOps), allWithPrefix = "wj_lv_mg_ht_")
        org.mergeSamples(targetSpec = md({"name":"VV", "color": r.kOrange+3}, mcOps), sources = ["ww_py.job188", "wz_py.job188", "zz_py.job188"])
        org.mergeSamples(targetSpec = md({"name":"ZH", "color":r.kMagenta}, mcOps), sources = ["zinv_hbb_125_powheg.job342"])
        org.mergeSamples(targetSpec = md({"name":"LM6", "color":r.kMagenta}, mcOps), allWithPrefix = "lm6")
        ewkSources = ["tt/t/ttz", "Z + jets", "W + jets", "VV"]

        org.mergeSamples(targetSpec = md({"name":"Standard Model ", "color":r.kAzure+6}, mcOps), sources = ewkSources + qcdSources, keepSources = True)

    def conclude(self, conf) :
        org = self.organizer(conf)
        ##for skimming only
        #utils.printSkimResults(org)            

        self.mergeSamples(org)
        org.scale()
        
        self.makeStandardPlots(org)
        #self.makeIndividualPlots(org)

    def makeStandardPlots(self, org) :
        #plot
        pl = supy.plotter(org,
                          pdfFileName = self.pdfFileName(org.tag),
                          samplesForRatios = ("2012 Data","Standard Model "),
                          sampleLabelsForRatios = ("data","s.m."),
                          #samplesForRatios = ("2012 Data","tt"),
                          #sampleLabelsForRatios = ("data","tt"),
                          printRatios = True,
                          showStatBox = True,
                          rowColors = [r.kBlack, r.kViolet+4],
                          #whiteList = ["lowestUnPrescaledTrigger"],
                          #doLog = False,
                          pegMinimum = 0.1,
                          linYAfter = ("variableGreaterFilter", "xcak5JetAlphaTEtPat>=0.550 "),
                          blackList = ["lumiHisto","xsHisto","nJobsHisto"],
                          )
        pl.plotAll()
        #self.makeEfficiencyPlots(org, tag, sampleName = "LM1")

    def makeIndividualPlots(self, org) :
        ht = "H_{T}^{#color[0]{T}}"
        if "275" in org.tag :
            htLabel = "275 < %s < 325"%ht
        if "325" in org.tag :
            htLabel = "325 < %s < 375"%ht
        if "375" in org.tag :
            htLabel = "375 < %s"%ht

        #plot all
        pl = supy.plotter(org,
                          pdfFileName = self.pdfFileName(org.tag),
                          showStatBox = False,
                          doLog = False,
                          #pegMinimum = 0.1,
                          anMode = False,
                          )
        pl.individualPlots(plotSpecs = [{"plotName":"xcak5JetAlphaTRoughPat",
                                         "stepName" :"alphaHistogrammer",
                                         "stepDesc" :"xcak5JetPat",
                                         "newTitle":";#alpha_{T};events / bin",
                                         "legendCoords": (0.55, 0.60, 0.85, 0.90),
                                         "stampCoords": (0.75, 0.55)
                                         },
                                        {"plotName":"jetMultiplicity",
                                         "stepName":"singleJetHistogrammer",
                                         "stepDesc":"xcak5JetPat through index 2",
                                         "newTitle":";N_{jets};events / bin",
                                         "legendCoords": (0.7, 0.7, 0.92, 0.92),
                                         "stampCoords": (0.5, 0.28),
                                         },
                                        {"plotName":"xcak5JetHtPat",
                                         "stepName":"cleanJetHtMhtHistogrammer",
                                         "stepDesc":"xcak5JetPat",
                                         "newTitle":";H_{T} (GeV);events / bin",
                                         "legendCoords": (0.6, 0.60, 0.92, 0.92),
                                         "stampCoords": (0.45, 0.88)
                                        },
                                        ##after alphaT
                                        {"plotName":"xcak5JetDeltaPhiStarPat",
                                         "stepName":"histogrammer",
                                         "stepDesc" :"(lambda x:x[0][0])(xcak5JetDeltaPhiStarPat)",
                                         "index" : -1,
                                         "newTitle":";#Delta#phi*;events / bin",
                                         "legendCoords": (0.6, 0.6, 0.92, 0.92),
                                         "stampCoords": (0.33, 0.88),
                                         },
                                        {"plotName":"xcak5JetHtPlusMhtRoughPat",
                                         "stepName":"cleanJetHtMhtHistogrammer",
                                         "stepDesc":"xcak5JetPat",
                                         "index":-1,
                                         "newTitle":";M_{eff} (GeV);events / bin",
                                         "legendCoords": (0.7, 0.7, 0.92, 0.92),
                                         "stampCoords": (0.75, 0.4),
                                         },
                                        {"plotName":"xcak5JetIndicesPat",
                                         "stepName":"histogrammer",
                                         "stepDesc" :"(lambda x:len(x))(xcak5JetIndicesPat)",
                                         "index" : -1,
                                         "newTitle":";N_{jets};events / bin",
                                         "legendCoords": (0.6, 0.6, 0.92, 0.92),
                                         "stampCoords": (0.6, 0.38),
                                         },
                                        {"plotName":"xcak5JetMbbListPat_2b",
                                         "stepName":"mbbHistogrammer",
                                         "stepDesc" :"mbbHistogrammer",
                                         "index" : -1,
                                         "stamp": False,
                                         "newTitle":"%s;m_{bb} (GeV) [2 b-jets];events / bin / 5.0 fb^{-1}"%htLabel,
                                         "legendCoords": (0.55, 0.6, 0.8, 0.75),
                                         "reBinFactor": 5,
                                         },
                                        ][-1:],
                           newSampleNames = {"tt": "Madgraph t#bar{t}",
                                             "2012 Data": "Data"},
                           #newSampleNames = {"qcd_mg_nVtx": "Madgraph QCD",
                           #                  "g_jets_mg_nVtx": "Madgraph #gamma + jets",
                           #                  "2011 Data": "Data",
                           #                  "standard_model_nVtx": "Standard Model",
                           #                  },
                           preliminary = True,
                           tdrStyle = False,
                           )


    def makeEfficiencyPlots(self, org, tag, sampleName) :
        def sampleIndex(org, name) :
            for iSample,sample in enumerate(org.samples) :
                if sample["name"]==name : return iSample
            assert False, "could not find sample %s"%name

        def numerAndDenom(org, var) :
            d = {}
            for selection in org.selections :
                if selection.name!= "passFilter" : continue
                if   "htLabel1" in selection.title : label = "before"
                elif "htLabel2" in selection.title : label = "after"
                else : continue
                if var in selection :
                    d[label] = selection[var][sampleIndex(org, sampleName)].Clone(label)
                
            return d

        keep = []
        canvas = r.TCanvas()
        canvas.SetRightMargin(0.2)
        canvas.SetTickx()
        canvas.SetTicky()
        psFileName = "%s.ps"%tag
        canvas.Print(psFileName+"[","Lanscape")

        assert len(self.parameters()["objects"])==1
        for key,value in self.parameters()["objects"].iteritems() :
            jet = value["jet"]

        for variable in ["%sSumEt%s"%jet] :
            histos = numerAndDenom(org, variable)
            if "before" not in histos or "after" not in histos : continue
            result = histos["after"].Clone(variable)
            result.Scale(1.0/histos["before"].Integral(0, histos["before"].GetNbinsX()+1))
            result.SetMarkerStyle(20)
            result.SetStats(False)
            if result.ClassName()[2]=="1" :
                #result.GetYaxis().SetRangeUser(0.0,1.0)
                result.GetYaxis().SetTitle("efficiency")
                result.Draw()
            else :
                #result.GetZaxis().SetRangeUser(0.0,1.0)
                result.GetZaxis().SetTitle("efficiency")
                result.Draw("colz")
            canvas.Print(psFileName,"Lanscape")

        canvas.Print(psFileName+"]","Lanscape")                
        os.system("ps2pdf "+psFileName)
        os.remove(psFileName)

