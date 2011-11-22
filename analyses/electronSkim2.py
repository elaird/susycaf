import supy
import samples,calculables

electron = ("electron","Pat")

class electronSkim2(supy.analysis) :
    def listOfSteps(self,params) :
        stepList=[ supy.steps.printer.progressPrinter(),
                   supy.steps.filters.pt("%sP4%s"%electron, min = 20, indices = "%sIndicesAnyIso%s"%electron, index = 0),
                   supy.steps.other.skimmer(),
                   ]
        return stepList

    def listOfCalculables(self,pars) :
        return ( supy.calculables.zeroArgs() +
                 supy.calculables.fromCollections(calculables.Electron,[electron]) +
                 [calculables.Electron.Indices( electron, ptMin = 10, simpleEleID = "95", useCombinedIso = True)])
    
    def listOfSamples(self,pars) :
        return supy.samples.specify(names = [ "EG.Run2010A-Nov4ReReco_v1.RECO.Sparrow",
                                              "Electron.Run2010B-Nov4ReReco_v1.RECO.Sparrow",
                                              ])

    def listOfSampleDictionaries(self) :
        return [samples.electron]

    def conclude(self,pars) :
        org = self.organizer(pars)
        utils.printSkimResults(org)
