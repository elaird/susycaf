import supy,calculables,samples

class mumuSkim(supy.analysis) :
    def parameters(self) :
        return {"muon":("muon", "Pat")}
    
    def listOfSteps(self, params) :
        return [ supy.steps.printer.progressPrinter(),
                 supy.steps.filters.multiplicity("%sIndices%s"%params["muon"], min = 2),
                 supy.steps.other.skimmer(),
                 ]

    def listOfCalculables(self, params) :
        return ( supy.calculables.zeroArgs(supy.calculables) +
                 supy.calculables.zeroArgs(calculables) +
                 supy.calculables.fromCollections(calculables.muon,[params["muon"]]) +
                 [calculables.muon.Indices( params["muon"], ptMin = 10, combinedRelIsoMax = 0.50)])
    
    def listOfSamples(self, params) :
        from supy.samples import specify
        out = []
        out += specify(names = "DoubleMu.Run2011A-05Aug2011-v1.AOD.job663"  )
        out += specify(names = "DoubleMu.Run2011A-May10ReReco-v1.AOD.job662")
        out += specify(names = "DoubleMu.Run2011A-PromptReco-v4.AOD.job664" )
        out += specify(names = "DoubleMu.Run2011A-PromptReco-v6.AOD.job665" )
        out += specify(names = "DoubleMu.Run2011B-PromptReco-v1.AOD.job666" )
        return out
    
    def listOfSampleDictionaries(self) :
        return [samples.muon, samples.mumu, samples.mc]

    def conclude(self, pars) :
        org = self.organizer(pars)
        utils.printSkimResults(org)
