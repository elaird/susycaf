import supy
import calculables,samples

class electronSkim(supy.analysis) :
    def parameters(self) :
        return { "useCombinedIso": dict([ ("combinedIso", True), ("relativeIso", False) ]  [:] ),
                 "electrons" : ("electron","Pat"),
                 }

    def listOfSteps(self,params) :
        stepList=[ supy.steps.progressPrinter(),
                   supy.steps.multiplicityFilter("%sIndices%s"%params["electrons"], nMin = 1),
                   supy.steps.skimmer(),
                   ]
        return stepList

    def listOfCalculables(self,pars) :
        return ( supy.calculables.zeroArgs() +
                 supy.calculables.fromCollections("calculablesElectron",[params["electrons"]]) +
                 [calculables.electronIndices( params["electrons"], ptMin = 20, simpleEleID = "80", useCombinedIso = params["useCombinedIso"])])
    
    def listOfSamples(self,pars) :
        return supy.samples.specify(names = ["Run2010B_MJ_skim2",
                                             "Run2010B_MJ_skim",
                                             "Run2010B_J_skim2",
                                             "Run2010B_J_skim",
                                             "Run2010A_JM_skim",
                                             "Run2010A_JMT_skim",
                                             ])

    def listOfSampleDictionaries(self) :
        return [samples.jetmet,samples.mc]

    def conclude(self, pars) :
        org = self.organizer(pars)
        utils.printSkimResults(org)
