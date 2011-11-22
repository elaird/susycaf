import supy, calculables,steps,samples

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
        outList = supy.calculables.zeroArgs()

        for muon in nameList(params["recoAlgos"],"muon") :
            outList += supy.calculables.fromCollections(calculables.muon, [muon])
            outList += [calculables.muon.Indices(muon, ptMin = 10, combinedRelIsoMax = 0.15)]
            
        for obj in dict(params["recoAlgos"]).values() :
            outList += self.calcListJet(obj)

        return outList
    
    def listOfSamples(self, params) :
        from supy.samples import specify
        out = []
        out += specify(names = "znunu_jets_mg_ht_50_100")
        out += specify(names = "znunu_jets_mg_ht_100_200")
        out += specify(names = "znunu_jets_mg_ht_200_inf")
        return out

    def listOfSampleDictionaries(self) :
        return [samples.mc]

    def conclude(self, config) :
        utils.printSkimResults(self.organizer(config))
