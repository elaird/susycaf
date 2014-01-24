import supy,steps,calculables,samples
import ROOT as r

class aTTriggerLook(supy.analysis) :
    def useCachedFileLists(self) : return False
    
    def parameters(self) :
        objects = self.vary()
        fields =                           [ "jet",                        "jetId",     "muonsInJets",           "met",      "rechit",
                                             "compJet",                "compJetId", "compMuonsInJets",        "compMet", "compRechit",
                                             "muon",                    "electron",          "photon"]

        objects["pf"] = dict(zip(fields, [("xcak5JetPF","Pat"),       "JetIDloose",             False,   "metP4TypeIPF",       "PF",
                                            ("xcak5Jet","Pat"),     "JetIDtight",              True,   "metP4TypeIPF",         "Calo",
                                            ("muon","PF"),     ("electron","Pat"),  ("photon","Pat")]))

        objects["calo"] = dict(zip(fields, [("xcak5Jet","Pat"),       "JetIDloose",             False,   "metP4TypeIPF",       "Calo",
                                            ("xcak5JetPF","Pat"),     "JetIDtight",              True,   "metP4TypeIPF",         "PF",
                                            ("muon","Pat"),     ("electron","Pat"),  ("photon","Pat")]))

        return { "objects": objects,
                 "nJetsMinMax" :      self.vary(dict([ ("ge2",(2,None)),  ("2",(2,2)),  ("ge3",(3,None)),  ("3",(3,3)), ("e23",(2,3)), ("ge4",(4,None))][4:5] )),
                 #"nBTagJets":         self.vary(dict([ ("nbe0",(0,0)),  ("nbe1",(1,1)),  ("nbe2",(2,2)),  ("nbe3",(3,3)),  ("nbge4",(4,None)) ][0:1] )),
                 "etRatherThanPt" : True,
                 "lowPtThreshold" : 30.0,
                 "lowPtName" : "lowPt",
                 "highPtThreshold" : 50.0,
                 "highPtName" : "highPt",
                 "thresholds": self.vary(dict( [("375",        (375.0, None,  100.0, 50.0)),#0
                                                ("325_scaled", (325.0, 375.0,  86.7, 43.3)),#1
                                                ("275_scaled", (275.0, 325.0,  73.3, 36.7)),#2
                                                ("200_scaled", (200.0, 275.0,  73.3, 36.7)),#3
                                                ("875", (875.0, None,  100.0, 50.0)),#4
                                                ][0:1] )),
                 }

    def calcListJet(self, obj, etRatherThanPt, ptMin, lowPtThreshold, lowPtName, highPtThreshold, highPtName, htThreshold) :
        def calcList(jet, met, photon, muon, electron, muonsInJets, jetIdFlag) :
            print "WARNING: synchronize muon addition"
            outList = [
                calculables.xclean.xcJet(jet,
                                         gamma = photon,
                                         gammaDR = 0.5,
                                         muon = muon,
                                         muonDR = 0.5,
                                         correctForMuons = not muonsInJets,
                                         #correctForMuons = False,
                                         electron = electron,
                                         electronDR = 0.5),
                calculables.jet.Indices( jet, ptMin = ptMin,           etaMax = 3.0, flagName = jetIdFlag),
                calculables.jet.Indices( jet, ptMin = lowPtThreshold,  etaMax = 3.0, flagName = jetIdFlag, extraName = lowPtName),
                calculables.jet.Indices( jet, ptMin = highPtThreshold, etaMax = 3.0, flagName = jetIdFlag, extraName = highPtName),
                calculables.jet.IndicesBtagged2(jet, tag = "CombinedSecondaryVertexBJetTags", threshold = 0.679),
                calculables.jet.SumP4(jet),
                calculables.jet.SumP4(jet, extraName = lowPtName),
                calculables.jet.SumP4(jet, extraName = highPtName),
                calculables.jet.SumP4(jet, extraName = "Btagged2"),
                calculables.jet.DeltaPhiStar(jet, extraName = lowPtName),
                calculables.jet.DeltaPseudoJet(jet, etRatherThanPt),
                calculables.jet.MaxEmEnergyFraction(jet),
                calculables.jet.AlphaT(jet, etRatherThanPt),
                calculables.jet.AlphaTMet(jet, etRatherThanPt, met),
                calculables.jet.MhtOverMet((jet[0], highPtName+jet[1]), met),
                calculables.jet.MhtOverMet((jet[0], highPtName+jet[1]), "%sPlus%sIndices%s" % (met, muon[0],muon[1])),
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

            calculables.muon.Indices( obj["muon"], ptMin = 25, isoMax = 0.12, ISO = "PfIsolationR04DeltaBCorrected", ID = "IdPog2012Tight", absEtaMax=2.1),
            calculables.electron.Indices( obj["electron"], ptMin = 10, flag2012 = "Veto"),
            calculables.photon.Indices(obj["photon"], ptMin = 25, flagName = "photonIDRA3Pat"),
            calculables.photon.CombinedIsoDR03RhoCorrected(obj["photon"]),

            calculables.other.RecHitSumPt(obj["rechit"]),
            calculables.other.RecHitSumP4(obj["rechit"]),

            calculables.other.RecHitSumPt(obj["compRechit"]),
            calculables.other.RecHitSumP4(obj["compRechit"]),

            calculables.other.metPlusIndices(obj["met"],obj["muon"]),
            calculables.vertex.ID(),
            calculables.vertex.Indices(),
            ]
    
    def listOfCalculables(self, params) :
        obj = params["objects"]
        outList = []
        outList += supy.calculables.zeroArgs(supy.calculables)
        outList += supy.calculables.zeroArgs(calculables)
        outList += supy.calculables.fromCollections(calculables.muon, [obj["muon"]])
        outList += supy.calculables.fromCollections(calculables.electron, [obj["electron"]])
        outList += supy.calculables.fromCollections(calculables.photon, [obj["photon"]])
        outList += self.calcListOther(obj)
        outList += self.calcListJet(obj, params["etRatherThanPt"], params["thresholds"][3],
                                    params["lowPtThreshold"], params["lowPtName"], params["highPtThreshold"], params["highPtName"], params["thresholds"][0])
        return outList

    def stepsXclean(self, params):
        return [
            #steps.other.iterHistogrammer("ecalDeadTowerTrigPrimP4", 256, 0.0, 128.0, title=";E_{T} of ECAL TP in each dead region (GeV);TPs / bin", funcString="lambda x:x.Et()"),
            supy.steps.filters.multiplicity("%sIndices%s" % params["objects"]["electron"], max=0),
            supy.steps.filters.multiplicity("%sIndices%s" % params["objects"]["photon"],   max=0),
            supy.steps.filters.multiplicity("%sIndicesOther%s" % params["objects"]["jet"], max=0),
            supy.steps.filters.multiplicity("%sIndices%s" % params["objects"]["jet"], min=params["nJetsMinMax"][0], max=params["nJetsMinMax"][1]),
            #steps.jet.uniquelyMatchedNonisoMuons(_jet),
            ]
    def stepsMuonID(self, params):
        return [supy.steps.filters.multiplicity("%sIndices%s" % params["objects"]["muon"],     min=1, max=1)]

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

    def stepsPlotsOne(self, params):
        _jet = params["objects"]["jet"]
        _met = "%sPlus%sIndices%s" % (params["objects"]["met"],params["objects"]["muon"][0],params["objects"]["muon"][1])
        _et = "Et" if params["etRatherThanPt"] else "Pt"
        _muon = params["objects"]["muon"]

        return [
            supy.steps.histos.histogrammer("%sSum%s%s" % (_jet[0], _et, _jet[1]), 50, 0, 500,
                                           title=";H_{T} (GeV) from %s%s %ss;events / bin" % (_jet[0], _jet[1], _et)),
            supy.steps.histos.histogrammer("%sSumP4%s" % _jet, 50, 0, 500, title=";MHT from %s%s (GeV);events / bin" % _jet,
                                           funcString="lambda x:x.pt()"),
            supy.steps.histos.histogrammer("vertexIndices", 20, -0.5, 19.5, title=";N vertices;events / bin",
                                           funcString="lambda x:len(x)"),
            steps.jet.singleJetHistogrammer(_jet),
            steps.muon.muonHistogrammer(_muon, 0),
            supy.steps.histos.histogrammer(_met, 100, 0.0, 500.0, title=";"+_met+" (GeV);events / bin", funcString="lambda x: x.pt()"),
            supy.steps.histos.histogrammer("%sMht%sOver%s" % (_jet[0], params["highPtName"]+_jet[1], _met), 100, 0.0, 3.0,
                                           title=";MHT %s%s / %s;events / bin" % (_jet[0], params["highPtName"]+_jet[1], _met)),
            steps.jet.alphaHistogrammer(cs=_jet, deltaPhiStarExtraName=params["lowPtName"], etRatherThanPt=params["etRatherThanPt"]),
            ]

    def listOfSteps(self, params) :
        _jet = params["objects"]["jet"]
        _met = params["objects"]["met"]
        _et = "Et" if params["etRatherThanPt"] else "Pt"
        return ([supy.steps.printer.progressPrinter()]+
                 #supy.steps.filters.value("%sAlphaT%s%s" % (_jet[0], _et, _jet[1]), min=0.55)]+# max=5.0),
                #steps.trigger.hltFilterList(["HLT_IsoMu24_eta2p1_v15"])]+
                #steps.trigger.hltFilterList(["HLT_HT200_AlphaT0p57_v8"])
                #steps.trigger.hltFail(["HLT_HT200_AlphaT0p57_v8"]),
                
                #supy.steps.other.skimmer()]+
                self.stepsXclean(params)+
                self.stepsHtLeadingJets(params)+
                self.stepsMuonID(params)+
                self.stepsPlotsOne(params)+
                #supy.steps.filters.value("%sSumEt%s"%params["objects"]["jet"], min = params["thresholds"][0]),
                #supy.steps.filters.value("%sSumEt%s"%params["objects"]["jet"], max = params["thresholds"][1]),
                #supy.steps.filters.multiplicity("%sIndicesBtagged2%s"%params["objects"]["jet"], min = params["nBTagJets"][0], max = params["nBTagJets"][1]),
                #supy.steps.filters.multiplicity("%sIndices%s"%params["objects"]["jet"], min = params["nJetsMinMax"][0], max = params["nJetsMinMax"][1]),
                #supy.steps.filters.value("%sRecHitSumPt"%params["objects"]["rechit"], max = 30.0),
                #supy.steps.filters.value("%sMaxEmEnergyFraction%s"%params["objects"]["jet"], min = .1),
                [])
    
    def listOfSampleDictionaries(self) :
    
        return [samples.muon17]

    
    def listOfSamples(self,params) :
        from supy.samples import specify
        out = []
        out += specify(names="SingleMu.Run2012A-22Jan2013", nFilesMax=2)#, nEventsMax=10)
        out += specify(names="SingleMu.Run2012B-22Jan2013", nFilesMax=2)#, nEventsMax=10)
        out += specify(names="SingleMu.Run2012C-22Jan2013", nFilesMax=2)#, nEventsMax=10)
        out += specify(names="SingleMu.Run2012D-22Jan2013", nFilesMax=2)#, nEventsMax=10)
        return out 
            
    def conclude(self, conf):
        org = self.organizer(conf)
        #org.scale(toPdf=True) if not self.parameters()["signalScan"] else org.scale(100.0)
        self.makeStandardPlots(org)

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
