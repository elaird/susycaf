import supy,calculables,steps,samples

def nameList(t, name)  : return list(set([obj[name] for obj in dict(t).values()]))

class hadronicSkim(supy.analysis) :
    def parameters(self) :
        objects = {}
        fields =                           [ "jet",                "muon",       "muonsInJets", "jetPtMin"]
        objects["calo"] = dict(zip(fields, [("xcak5Jet","Pat"),   ("muon","Pat"),        False,      30.0 ]))
        #objects["pf"]   = dict(zip(fields, [("xcak5JetPF","Pat"), ("muon","Pat"),         True,      30.0 ]))
        return {"recoAlgos": tuple(objects.iteritems())}

    def listOfSteps(self, params) :
        stepList = [
            supy.steps.printer.progressPrinter(2,300),
            supy.steps.filters.multiplicity("genIndicesStatus1MuPlus", min = 1),
            supy.steps.filters.multiplicity("genIndicesStatus1MuMinus", min = 1),
            supy.steps.filters.multiplicity("genIndicesStatus3Z", min = 1, max = 1),
            supy.steps.histos.mass("genP4", 100, 0, 200, index = 0, indices = "genIndicesStatus3Z", xtitle="Z_{m} (GeV)"),
            supy.steps.filters.mass("genP4", index = 0, indices = "genIndicesStatus3Z", min = 65, max = 110),
            supy.steps.histos.value("xcak5JetSumEtPat",100,0,500,xtitle="H_{T} (GeV)"),
            steps.jet.htSelector(nameList(params["recoAlgos"], "jet"), 250.0),
            supy.steps.other.skimmer(),
            ]
        return stepList

    def calcListJet(self, obj) :
        outList = [
            calculables.xclean.xcJet(obj["jet"],
                                     gamma = None,
                                     gammaDR = 0.5,
                                     muon = obj["muon"],
                                     muonDR = 0.5,
                                     correctForMuons = not obj["muonsInJets"],
                                     electron = None,
                                     electronDR = 0.5),
            calculables.jet.Indices( obj["jet"], obj["jetPtMin"], etaMax = 3.0, flagName = "JetIDloose"),
            ]
        return outList+supy.calculables.fromCollections(calculables.jet, [obj["jet"]])
    
    def listOfCalculables(self, params) :
        outList = supy.calculables.zeroArgs(supy.calculables)

        for muon in nameList(params["recoAlgos"],"muon") :
            outList += supy.calculables.fromCollections(calculables.muon, [muon])
            outList += [calculables.muon.Indices(muon, ptMin = 10, combinedRelIsoMax = 0.15)]
            
        for obj in dict(params["recoAlgos"]).values() :
            outList += self.calcListJet(obj)

        outList += [calculables.gen.genIndices( pdgs = [13], label = "Status1MuPlus", status = [1]),
                    calculables.gen.genIndices( pdgs = [-13], label = "Status1MuMinus", status = [1]),
                    calculables.gen.genIndices( pdgs = [23], label = "Status3Z", status = [3]),
                    ]
        return outList
    
    def listOfSamples(self, params) :
        from supy.samples import specify
        out = []
        out += specify("DYJetsToLL_M-50.job173", markerStyle = 20)
        return out

    def listOfSampleDictionaries(self) :
        return [samples.photon17]

    def conclude(self, config) :
        org = self.organizer(config)
        supy.utils.io.printSkimResults(org)

        supy.plotter(org,
                     pdfFileName = self.pdfFileName(org.tag),
                     blackList = ["lumiHisto","xsHisto","nJobsHisto",],
                     ).plotAll()

