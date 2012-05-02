import supy,steps,calculables,samples,ROOT as r

class emuLook(supy.analysis) :
    def parameters(self) :
        objects = self.vary()
        fields =                             [ "jet",                        "jetId",     "muonsInJets",           "met",
                                               "compJet",                "compJetId", "compMuonsInJets",        "compMet",
                                               "muon",                    "electron",          "photon",         "rechit"]

        objects["calo"]   = dict(zip(fields, [("xcak5Jet","Pat"),       "JetIDloose",             False, "metP4AK5TypeII",
                                              ("xcak5JetPF","Pat"),     "JetIDtight",              True,        "metP4PF",
                                              ("muon","Pat"),     ("electron","Pat"),  ("photon","Pat"),           "Calo",
                                              ]))
        
        #objects["pf"]     = dict(zip(fields, [("xcak5JetPF","Pat"),     "JetIDtight",              True,        "metP4PF",
        #                                      ("xcak5Jet","Pat"),       "JetIDloose",             False, "metP4AK5TypeII",
        #                                      ("muon","Pat"),     ("electron","Pat"),  ("photon","Pat"),             "PF",
        #                                      ]))
        
        #objects["pf2pat"] = dict(zip(fields, [("xcak5JetPF2PAT","Pat"), "PFJetIDtight",            True,        "metP4PF",
        #                                      ("xcak5JetPF","Pat"),     "JetIDtight",              True,        "metP4PF",
        #                                      ("muon","PF"),       ("electron","PF"),  ("photon","Pat"),             "PF",
        #                                      ]))
        
        return {"objects": objects}

    def listOfCalculables(self, params) :
        obj = params["objects"]
        out  = supy.calculables.zeroArgs(supy.calculables)
        out += supy.calculables.zeroArgs(calculables)
        out += supy.calculables.fromCollections(calculables.muon, [obj["muon"]])
        out += supy.calculables.fromCollections(calculables.electron, [obj["electron"]])
        out += supy.calculables.fromCollections(calculables.photon, [obj["photon"]])
        out += supy.calculables.fromCollections(calculables.jet, [obj["jet"]])
        out += [
            calculables.muon.Indices( obj["muon"], ptMin = 10, combinedRelIsoMax = 0.15),
            calculables.electron.Indices( obj["electron"], ptMin = 10, simpleEleID = "95", useCombinedIso = True),
            calculables.photon.Indices(obj["photon"],  ptMin = 25, flagName = "photonIDTightFromTwikiPat"),
            
            calculables.vertex.ID(),
            calculables.vertex.Indices(),
            calculables.xclean.xcJet(obj["jet"],
                                     gamma = obj["photon"], gammaDR = 0.5,
                                     muon = obj["muon"], muonDR = 0.5, correctForMuons = not obj["muonsInJets"],
                                     electron = obj["electron"], electronDR = 0.5
                                     ),
            calculables.jet.Indices( obj["jet"], ptMin = 30.0, etaMax = 3.0, flagName = obj["jetId"]),
            ]
        return out
    
    def listOfSteps(self, params) :
        _jet  = params["objects"]["jet"]
        _electron = params["objects"]["electron"]
        _muon = params["objects"]["muon"]

        return [
            supy.steps.printer.progressPrinter(),

            steps.trigger.l1Filter("L1Tech_BPTX_plus_AND_minus.v0").onlyData(),
            steps.trigger.physicsDeclaredFilter().onlyData(),
            steps.filters.monster().onlyData(),
            steps.filters.hbheNoise().onlyData(),

            supy.steps.filters.multiplicity("vertexIndices",         min = 1),
            supy.steps.filters.multiplicity("%sIndices%s"%_muon,     min = 1),
            supy.steps.filters.multiplicity("%sIndices%s"%_electron, min = 1),
            supy.steps.filters.multiplicity("%sIndices%s"%_jet,      min = 2),
            #supy.steps.filters.multiplicity("%sIndicesOther%s"%_jet, max = 0),

            supy.steps.histos.multiplicity("vertexIndices"),
            
            #steps.other.skimmer(),
            #steps.Displayer.displayer(jets = _jet,
            #                          muons = _muon,
            #                          met       = params["objects"]["met"],
            #                          electrons = params["objects"]["electron"],
            #                          photons   = params["objects"]["photon"],                            
            #                          recHits   = params["objects"]["rechit"], recHitPtThreshold = 1.0,#GeV
            #                          scale = 400.0,#GeV
            #                          #doGenJets = True,
            #                          ),
            ]
    
    def listOfSampleDictionaries(self) :
        return [samples.top17]

    def listOfSamples(self,params) :
        from supy.samples import specify
        return specify(names = "ttj_8TeV_mg_job79", nFilesMax = 1, nEventsMax = 10000)
    
    def conclude(self, conf) :
        org = self.organizer(conf)
        org.scale(5.0e3)
        supy.plotter(org,
                     pdfFileName = self.pdfFileName(org.tag),
                     rowColors = [r.kBlack, r.kViolet+4],
                     #doLog = False,
                     pegMinimum = 0.1,
                     blackList = ["lumiHisto","xsHisto","nJobsHisto"],
                     ).plotAll()
