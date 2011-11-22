import supy,calculables,samples

muon = ("muon","Pat")

class muonSkim(supy.analysis) :
    def listOfSteps(self,params) :
        return [ supy.steps.printer.progressPrinter(),
                 supy.steps.filters.multiplicity("%sIndices%s"%muon, min = 1),
                 #steps.Other.objectEtaSelector(muon, etaThreshold = 2.5, index = 0, p4String = "P4"),
                 supy.steps.other.skimmer(),
                 ]

    def listOfCalculables(self,params) :
        return ( supy.calculables.zeroArgs() +
                 supy.calculables.fromCollections(calculables.muon,[muon]) +
                 [calculables.muon.Indices( muon, ptMin = 10, combinedRelIsoMax = 0.50)])
    
    def listOfSamples(self,params) :
        from supy.samples import specify
        return (
            #specify(names = "Mu.Run2010A-Sep17ReReco_v2.RECO.Robin") +
            #specify(names = "Mu.Run2010B-PromptReco-v2.RECO.Arlo1") +
            #specify(names = "Mu.Run2010B-PromptReco-v2.RECO.Arlo2") +
            #specify(names = "Mu.Run2010B-PromptReco-v2.RECO.Martyn") +
            #specify(names = "v12_qcd_py6_pt30") +
            #specify(names = "v12_qcd_py6_pt80") +
            #specify(names = "v12_qcd_py6_pt170") +
            #specify(names = "v12_qcd_py6_pt300") +
            #specify(names = "Run2010B_MJ_skim3") +
            #specify(names = "Run2010B_MJ_skim4") +
            specify(names ="MultiJet.Run2010B-PromptReco-v2.RECO.RAW.Robin") +
            []
            )

    def listOfSampleDictionaries(self) :
        return [samples.jetmet, samples.muon, samples.mc]

    def conclude(self,pars) :
        org = self.organizer( pars )
        utils.printSkimResults(org)
