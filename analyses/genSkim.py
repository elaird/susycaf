import supy,steps,calculables,samples,os,math, ROOT as r

class smsSkim(supy.analysis) :
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
                 "nJetsMinMax" :      self.vary(dict([ ("ge2",(2,None)),  ("2",(2,2)),  ("ge3",(3,None)),  ("3",(3,3)), ("e23",(2,3)), ("ge4",(4,None))][0:1] )),
                 "etRatherThanPt" : [True,False][0],
                 "lowPtThreshold" : 30.0,
                 "lowPtName" : "lowPt",
                 "highPtThreshold" : 50.0,
                 "highPtName" : "highPt",
                 "tanBeta" : [None, 3, 10, 50][0],
                 "signalScan" : False,
                 "thresholds": self.vary(dict( [("275",        (275.0, 325.0, 100.0, 50.0)),#0
                                                ("325",        (325.0, 375.0, 100.0, 50.0)),#1
                                                ("375",        (375.0, None,  100.0, 50.0)),#2
                                                ("275_scaled", (275.0, None,   10.0, 10.0)),#3*****No upper HT cut
                                                ("275_scaled", (275.0, 325.0,  73.3, 36.7)),#4
                                                ("325_scaled", (325.0, 375.0,  86.7, 43.3)),#5
                                                ("375",        (375.0, None,  100.0, 50.0)),#6
                                                ][3:4] )),
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
                                         #correctForMuons = not muonsInJets,
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

    def calcListOther(self, obj) :
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
            ]
    
    def listOfCalculables(self, params) :
        obj = params["objects"]
        outList  = supy.calculables.zeroArgs(supy.calculables)
        outList += supy.calculables.zeroArgs(calculables)
        outList += supy.calculables.fromCollections(calculables.muon, [obj["muon"]])
        outList += supy.calculables.fromCollections(calculables.electron, [obj["electron"]])
        outList += supy.calculables.fromCollections(calculables.photon, [obj["photon"]])
        outList += self.calcListJet(obj, params["etRatherThanPt"], params["thresholds"][3],
                                    params["lowPtThreshold"], params["lowPtName"], params["highPtThreshold"], params["highPtName"], params["thresholds"][0])
        ]
        return outList

    def stepsGenValidation(self,params) :
        return [supy.steps.histos.histogrammer("genPartonHT",200,0,1000,title=";parton H_{T} (GeV);events / bin").onlySim(),
                #supy.steps.histos.histogrammer("genpthat",200,0,2000,title=";#hat{p_{T}} (GeV);events / bin").onlySim(),
                ]

    def stepsPlotsOne(self, params) :
        _jet = params["objects"]["jet"]
        _met = params["objects"]["met"]
        _et = "Et" if params["etRatherThanPt"] else "Pt"
        return [
            supy.steps.histos.histogrammer("%sSum%s%s"%(_jet[0], _et, _jet[1]), 50, 0, 2500,
                                           title = ";H_{T} (GeV) from %s%s %ss;events / bin"%(_jet[0], _jet[1], _et)),
            supy.steps.histos.histogrammer("%sSum%s%s"%(_jet[0], _et, _jet[1]), 60, 675, 1275,
                                           title = ";H_{T} (GeV) from %s%s %ss;events / bin"%(_jet[0], _jet[1], _et)),
            supy.steps.histos.histogrammer("%sSum%s%s"%(_jet[0], _et, _jet[1]), 100, 0, 1000,
                                           title = ";H_{T} (GeV) from %s%s %ss;events / bin"%(_jet[0], _jet[1], _et)),
            
            supy.steps.histos.histogrammer("%sSumP4%s"%_jet, 50, 0, 500, title = ";MHT from %s%s (GeV);events / bin"%_jet, funcString = "lambda x:x.pt()"),
            supy.steps.histos.histogrammer("vertexIndices", 20, -0.5, 19.5, title=";N vertices;events / bin", funcString="lambda x:len(x)"),
            supy.steps.histos.histogrammer("vertexSumPt", 100, 0.0, 1.0e3, title = ";SumPt of 2nd vertex (GeV);events / bin", funcString = "lambda x:([0.0,0.0]+sorted(x))[-2]"),
            
            #supy.steps.filters.label("singleJetPlots1"),
            #steps.jet.singleJetHistogrammer(_jet),
            supy.steps.filters.label("jetSumPlots1"), 
            steps.jet.cleanJetHtMhtHistogrammer(_jet,params["etRatherThanPt"]),
            #supy.steps.histos.histogrammer("%sDeltaPhiStar%s%s"%(_jet[0], params["lowPtName"],
            #                                                     _jet[1]), 20, 0.0, r.TMath.Pi(), title = ";#Delta#phi*;events / bin", funcString = 'lambda x:x[0][0]'),
            #supy.steps.histos.histogrammer("%sDeltaPhiStar%s"%(_jet[0], _jet[1]), 20, 0.0, r.TMath.Pi(), title = ";#Delta#phi*;events / bin", funcString = 'lambda x:x[0][0]'),
            #supy.steps.histos.histogrammer("%sMaxEmEnergyFraction%s"%(_jet[0], _jet[1]), 20, 0.0, 1.0, title = ";MaxEmEnergyFraction;events / bin"),
            supy.steps.histos.histogrammer(_met,100,0.0,500.0,title=";"+_met+" (GeV);events / bin", funcString = "lambda x: x.pt()"),
            supy.steps.filters.label("kinematicPlots1"),

            #steps.jet.alphaHistogrammer(cs = _jet, deltaPhiStarExtraName = params["lowPtName"], etRatherThanPt = params["etRatherThanPt"]),
            ]

    def stepsOptional(self, params) :
        return [
            supy.steps.filters.value("genPartonHT", min = 10, max = 200).onlySim(),
            supy.steps.other.skimmer(),
            #steps.Print.jetPrinter(_jet),
            #steps.Print.particleP4Printer(_muon),
            #steps.Print.particleP4Printer(_photon),
            #steps.Gen.genParticlePrinter(minPt = 10.0, minStatus = 3),
            ]

    def listOfSteps(self, params) :
        _jet  = params["objects"]["jet"]
        _met  = params["objects"]["met"]
        _etRatherThanPt = params["etRatherThanPt"]
        _et = "Et" if _etRatherThanPt else "Pt"

        return ([supy.steps.printer.progressPrinter()] +
                self.stepsGenValidation(params) +
                self.stepsPlotsOne(params) +
                self.stepsOptional(params) +
                self.stepsGenValidation(params) +
                [])

    def listOfSampleDictionaries(self) :
        sh = supy.samples.SampleHolder()
        return [samples.ewk17, sh, samples.dyll17]
    
    def listOfSamples(self,params) :
        from supy.samples import specify

        def w_inclusive() :
            out = []
            #out += specify(names = "wj_lv_mg_ht_incl", color = r.kOrange, nFilesMax = 1)
            #out += specify(names = "wj_lv_mg_ht_incl_v2.job50_part1", color = r.kOrange)#, nFilesMax = 1, nEventsMax = 10000)
            #out += specify(names = "wj_lv_mg_ht_incl_v2.job50_part2", color = r.kOrange)#, nFilesMax = 1, nEventsMax = 10000)
            #out += specify(names = "wj_lv_mg_ht_incl_v2.job50_part3", color = r.kOrange)#, nFilesMax = 1, nEventsMax = 10000)
            out += specify(names = "wj_lv_mg_ht_incl_v2.job50_part4", color = r.kOrange)#, nFilesMax = 1, nEventsMax = 10000)
            return out

        def dyll() :
            out = []
            #out += specify(names = "dyll_M-50_mg_job40_part1", color = r.kOrange)#, nFilesMax = 1, nEventsMax = 10000)
            #out += specify(names = "dyll_M-50_mg_job40_part2", color = r.kOrange)#, nFilesMax = 1, nEventsMax = 10000)
            out += specify(names = "dyll_M-10To50_mg_job40", color = r.kOrange)#, nFilesMax = 1, nEventsMax = 10000)
            return out

        return (
            #w_inclusive() +
            dyll() +
            []
            )

    def mergeSamples(self, org) :
        def md(x, y) :
            x.update(y)
            return x

        mcOps = {"markerStyle":1, "lineWidth":3, "goptions":"hist"}

        org.mergeSamples(targetSpec = md({"name":"W + jets", "color": r.kOrange-3}, mcOps), allWithPrefix = "wj_lv_mg_ht_")

    def conclude(self, conf) :
        org = self.organizer(conf)
        #for skimming only
        utils.printSkimResults(org)            

        self.mergeSamples(org)
        org.scale() if not self.parameters()["signalScan"] else org.scale(100.0)
        
        self.makeStandardPlots(org)
        #self.makeIndividualPlots(org)
        #self.makeEfficiencyPlots(org)

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
                          doLog = True,
                          #pegMinimum = 0.1,
                          linYAfter = ("variableGreaterFilter", "xcak5JetAlphaTEtPat>=0.550 "),
                          blackList = ["lumiHisto","xsHisto","nJobsHisto"],
                          )
        pl.plotAll()
                                    
                                    
