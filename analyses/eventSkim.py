import supy,steps,samples

class eventSkim(supy.analysis) :
    def listOfSteps(self,_) :
        return [ steps.printer.progressPrinter(2,300),
                 steps.filters.runLsEvent("/home/hep/elaird1/oneEvent.txt"),
                 supy.steps.other.skimmer(),
                 ]

    def listOfCalculables(self,_) :
        return supy.calculables.zeroArgs(supy.calculables)

    def listOfSamples(self,_) :
        from supy.samples import specify
        out = []
        out += specify(names = "DoubleMu.Run2011A-05Aug2011-v1.AOD.job663_skim"  )
        out += specify(names = "DoubleMu.Run2011A-May10ReReco-v1.AOD.job662_skim")
        out += specify(names = "DoubleMu.Run2011A-PromptReco-v4.AOD.job664_skim" )
        out += specify(names = "DoubleMu.Run2011A-PromptReco-v6.AOD.job665_skim" )
        out += specify(names = "DoubleMu.Run2011B-PromptReco-v1.AOD.job666_skim" )
        return out

    def listOfSampleDictionaries(self) :
        return [samples.jetmet, samples.mumu, samples.mc]
