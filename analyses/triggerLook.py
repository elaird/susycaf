import supy,steps,samples

class triggerLook(supy.analysis) :

    def listOfCalculables(self,_) : return supy.calculables.zeroArgs()

    def listOfSampleDictionaries(self) : return [samples.ht, samples.photon, samples.mumu]

    def listOfSteps(self,params) :  return [ steps.trigger.Counts(useCache = True) ]

    def listOfSamples(self,params) :
        from supy.samples import specify
        out = []

        #out += specify(names = "Photon.Run2011A-05Aug2011-v1.AOD.job663")
        #out += specify(names = "Photon.Run2011A-May10ReReco-v1.AOD.job662")
        #out += specify(names = "Photon.Run2011A-PromptReco-v4.AOD.job664")
        #out += specify(names = "Photon.Run2011A-PromptReco-v6.AOD.job667")
        #out += specify(names = "Photon.Run2011B-PromptReco-v1.AOD.job668")

        out += specify(names = "DoubleMu.Run2011A-05Aug2011-v1.AOD.job663",  )
        out += specify(names = "DoubleMu.Run2011A-May10ReReco-v1.AOD.job662",)
        out += specify(names = "DoubleMu.Run2011A-PromptReco-v4.AOD.job664", )
        out += specify(names = "DoubleMu.Run2011A-PromptReco-v6.AOD.job665", )
        out += specify(names = "DoubleMu.Run2011B-PromptReco-v1.AOD.job666", )

        return out

    def conclude(self, pars) :
        org = self.organizer(pars)
        supy.plotter(org, psFileName = self.psFileName(org.tag), blackList = ["lumiHisto","xsHisto","nJobsHisto"]).plotAll()
