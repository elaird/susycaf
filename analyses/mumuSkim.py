import supy,calculables,samples

class mumuSkim(supy.analysis) :
    def parameters(self) :
        return {"muon":("muon", "Pat")}
    
    def listOfSteps(self, params) :
        muon = params["muon"]
        return [ supy.steps.printer.progressPrinter(),
                 supy.steps.filters.multiplicity("%sIndices%s"%muon, min = 2),
                 #supy.steps.filters.pt(var = "%sP4%s"%muon, min = 17.0, indices = "%sIndices%s"%muon, index = 0),
                 supy.steps.other.skimmer(),
                 ]

    def listOfCalculables(self, params) :
        muon = params["muon"]
        return ( supy.calculables.zeroArgs(supy.calculables) +
                 supy.calculables.zeroArgs(calculables) +
                 supy.calculables.fromCollections(calculables.muon,[muon]) +
                 [calculables.muon.Indices( muon, ptMin = 10, combinedRelIsoMax = 0.50)])
    
    def listOfSamples(self, params) :
        from supy.samples import specify
        out = []

        #out += specify(names = "DoubleMu.Run2011A-05Aug2011-v1.AOD.job663"  )
        #out += specify(names = "DoubleMu.Run2011A-May10ReReco-v1.AOD.job662")
        #out += specify(names = "DoubleMu.Run2011A-PromptReco-v4.AOD.job664" )
        #out += specify(names = "DoubleMu.Run2011A-PromptReco-v6.AOD.job665" )
        #out += specify(names = "DoubleMu.Run2011B-PromptReco-v1.AOD.job666" )

        out += specify(names = "dyll_jets_mg_summer11")
        out += specify(names = "tt_jets_mg_tauola_summer11")
        
        return out
    
    def listOfSampleDictionaries(self) :
        return [samples.muon, samples.mumu, samples.mc]

    def conclude(self, pars) :
        org = self.organizer(pars)
        supy.utils.io.printSkimResults(org)
