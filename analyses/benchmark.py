import supy,samples

class benchmark(supy.analysis) :

    def listOfSampleDictionaries(self) :  return [samples.jetmet]

    def listOfSamples(self,_) :  return supy.samples.specify(names = "JetMETTau.Run2010A-Nov4ReReco_v1.RECO.Burt", nFilesMax = 10)

    def listOfCalculables(self,_) :  return []
    
    def listOfSteps(self,_) :
        touch = [
            #"triggered",
            #"prescaled",
            #"L1triggered",
            #"L1prescaled"
            "ak5JetCorrectedP4Pat",
            "ak5JetJPTCorrectedP4Pat",
            "ak5JetPFCorrectedP4Pat",
            "muonP4Pat",
            "electronP4Pat",
            ]
        return [ supy.steps.other.touchstuff(touch) ]


    
