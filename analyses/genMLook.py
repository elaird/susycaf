import supy,steps,samples, ROOT as r

class genMLook(supy.analysis) :

    def listOfCalculables(self,_) :  return supy.calculables.zeroArgs()

    def listOfSampleDictionaries(self) :  return [samples.mc, samples.jetmet]

    def listOfSamples(self,params) :
        return supy.samples.specify(names = "t2", nFilesMax = 1, nEventsMax = 10, color = r.kYellow-3)

    def listOfSteps(self,_) :
        outList=[
            supy.steps.printer.progressPrinter(),
            steps.Gen.genParticlePrinter(minPt = -1.0, minStatus = 3),
            #steps.Gen.genMassHistogrammer(),
            #steps.Gen.genSHatHistogrammer(),
            ]
        return outList

    def conclude(self,pars) :
        org = self.organizer(pars)
        org.scale(100.0)
        
        supy.plotter( org,
                      psFileName = self.psFileName(""),
                      ).plotAll()
