import supy,calculables,samples

muon = ("muon","PF")

class muonSkim2(supy.analysis) :
    def listOfSteps(self,params) :
        stepList=[ supy.steps.printer.progressPrinter(),
                   supy.steps.filters.pt("%sP4%s"%muon, min = 24, indices = "%sIndicesAnyIso%s"%muon, index = 0),
                   supy.steps.filters.absEta("%sP4%s"%muon, max = 2.2, indices = "%sIndicesAnyIso%s"%muon, index = 0),
                   supy.steps.other.skimmer(),
                   ]
        return stepList

    def listOfCalculables(self,params) :
        return (supy.calculables.zeroArgs() +
                supy.calculables.fromCollections(calculables.Muon,[muon]) +
                [calculables.muon.Indices( muon, ptMin = 10, combinedRelIsoMax = 0.15)] )
    
    def listOfSamples(self,params) :
        return supy.samples.specify(names = ["SingleMu.Run2011A-PR-v4.FJ.Burt",
                                             "SingleMu.Run2011A-May10-v1.FJ.Burt"])

    def listOfSampleDictionaries(self) :
        return [samples.jetmet, samples.muon, samples.mc]

    def conclude(self,conf) :
        org = self.organizer(conf)
        utils.printSkimResults(org)
