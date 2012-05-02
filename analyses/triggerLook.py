import supy,steps,samples

class triggerLook(supy.analysis) :

    #def useCachedFileLists(self) : return False
    def listOfCalculables(self,_) : return supy.calculables.zeroArgs(supy)

    def listOfSampleDictionaries(self) : return [samples.photon17]

    def listOfSteps(self,params) :  return [ steps.trigger.Counts(useCache = True) ]

    def listOfSamples(self,params) :
        from supy.samples import specify
        out = []

        #out += specify(names = "Photon.Run2012A-PromptReco-v1.AOD.job29")
        #out += specify(names = "Photon.Run2012A-PromptReco-v1.AOD.job44")
        #out += specify(names = "Photon.Run2012A-PromptReco-v1.AOD.job57")
        #out += specify(names = "Photon.Run2012A-PromptReco-v1.AOD.job69")
        out += specify(names = "Photon.Run2012A-PromptReco-v1.AOD.job74", nFilesMax = 1)

        return out

    def conclude(self, pars) :
        org = self.organizer(pars)
        supy.plotter(org, pdfFileName = self.pdfFileName(org.tag), blackList = ["lumiHisto","xsHisto","nJobsHisto"]).plotAll()
