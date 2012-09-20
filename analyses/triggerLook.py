import supy,steps,samples

class triggerLook(supy.analysis) :

    #def useCachedFileLists(self) : return False
    def listOfCalculables(self,_) : return supy.calculables.zeroArgs(supy)

    def listOfSampleDictionaries(self) : return [samples.photon17]

    def listOfSteps(self,params) :  return [ steps.trigger.Counts(useCache = True) ]

    def listOfSamples(self,params) :
        from supy.samples import specify
        out  = []
        out += specify("Photon.Run2012A-PromptReco-v1.AOD.job171")
        out += specify("SinglePhoton.Run2012B-PromptReco-v1.AOD.job171")

        return out

    def conclude(self, pars) :
        org = self.organizer(pars)
        supy.plotter(org, pdfFileName = self.pdfFileName(org.tag), blackList = ["lumiHisto","xsHisto","nJobsHisto"]).plotAll()
